#!/bin/sh
# Checks that the library defines no global symbols other than its own API: in particular, none of the
# bundled libsais, which would clash with another copy of it linked into the same program.
#
# usage: symbols.sh <library or object file>...
#
# The list defaults to the environment variable BZIP3_TEST_LIBS. Files that do not exist are skipped.

# shellcheck disable=SC2086
[ $# -gt 0 ] || set -- $BZIP3_TEST_LIBS

NM=${NM:-nm}
command -v "$NM" > /dev/null 2>&1 || { echo "nm not available, skipping"; exit 77; }

status=0
checked=0
for lib in "$@"; do
    [ -f "$lib" ] || continue
    case $lib in
        *.so*) flags="-D -g --defined-only" ;;
        *) flags="-g --defined-only" ;;
    esac
    # shellcheck disable=SC2086
    listing=$("$NM" $flags "$lib" 2> /dev/null) || continue
    checked=$((checked + 1))
    # Names may carry a leading underscore (Mach-O, 32-bit Windows); ignore linker-defined symbols.
    stray=$(echo "$listing" | awk 'NF >= 3 && $2 ~ /^[A-Z]$/ { print $3 }' |
        grep -v -E '^_?bz3_|^_?_(init|fini|edata|end|bss_start)|^__bss_start|^_ITM_|^__gmon|^_?__' | sort -u)
    if [ -n "$stray" ]; then
        echo "FAIL: $lib exports symbols outside of the bz3_ namespace:" >&2
        echo "$stray" | sed 's/^/    /' >&2
        status=1
    fi
done

[ "$checked" -gt 0 ] || { echo "no library to inspect, skipping"; exit 77; }
[ "$status" -eq 0 ] && echo "exported symbols OK ($checked file(s))"
exit $status
