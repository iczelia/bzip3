#!/bin/sh
#  Copyright (C) 2022-2026 Kamila Szewczyk
#
#  This program is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
#  General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

dir=${BZ3_TEST_LIBS:-.libs}
syms=${BZ3_TEST_SYMS:-${srcdir:-.}/src/libbz3.sym}
set -- ${NM:-nm}
command -v "$1" > /dev/null 2>&1 || { echo "symbols.sh: no nm, skipped";  exit 77; }
test -f "$syms" || { echo "symbols.sh: no $syms" >&2;  exit 99; }
status=0
checked=0

#  The defined global symbols of a file, without the underscore Mach-O and 32-bit Windows put before them, and
#  without those the toolchain adds.
globals() {
  "$@" 2> /dev/null | awk 'NF >= 3 && $(NF - 1) ~ /^[A-TV-Z]$/ { print $NF }' | sed 's/^_//' \
    | grep -v -E '^_?(_?init|_?fini|_?edata|_?end|_bss_start|_?ITM_|_gmon|_.*)$' | sort -u
}

for lib in "$dir"/libbzip3.a; do
  test -f "$lib" || continue
  checked=`expr $checked + 1`
  stray=`globals "$@" -g "$lib" | grep -v '^bz3_'`
  if test -n "$stray"; then
    echo "symbols.sh: $lib defines symbols outside of the bz3_ namespace:"
    echo "$stray"
    status=1
  fi
done

for lib in "$dir"/libbzip3.so "$dir"/libbzip3.dylib; do
  test -f "$lib" || continue
  checked=`expr $checked + 1`
  case $lib in *.so) flags="-D -g" ;; *) flags="-g" ;; esac
  globals "$@" $flags "$lib" > symbols.sh.$$
  if ! sort -u "$syms" | cmp -s - symbols.sh.$$; then
    echo "symbols.sh: $lib does not export exactly the functions of $syms:"
    sort -u "$syms" | diff - symbols.sh.$$
    status=1
  fi
  rm -f symbols.sh.$$
done

test $checked -gt 0 || { echo "symbols.sh: no library to inspect, skipped";  exit 77; }
test $status -eq 0 && echo "symbols.sh: $checked libraries define only their own symbols"
exit $status
