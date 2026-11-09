#!/bin/sh
# Tests of the bzip3 command line tool.
#
# usage: cli.sh <path to bzip3> <path to tests/data>
#
# Both default to the environment variables BZIP3 and BZIP3_TEST_DATA.
#
# The files in tests/data were compressed by bzip3 1.5.5; being able to decompress them guards the
# backwards compatibility of the decoder.

BZIP3=${1:-${BZIP3:?usage: cli.sh <bzip3> <data directory>}}
DATA=${2:-${BZIP3_TEST_DATA:?usage: cli.sh <bzip3> <data directory>}}

case $BZIP3 in /*) ;; *) BZIP3=$PWD/$BZIP3 ;; esac
case $DATA in /*) ;; *) DATA=$PWD/$DATA ;; esac

TMP=$(mktemp -d "${TMPDIR:-/tmp}/bzip3-test.XXXXXX") || exit 99
trap 'rm -rf "$TMP"' EXIT
cd "$TMP" || exit 99

failures=0
fail() {
    echo "FAIL: $*" >&2
    failures=$((failures + 1))
}

# Reference files from an older release decode into the original data.
for packed in "$DATA"/*.bz3; do
    orig=${packed%.bz3}
    name=$(basename "$orig")
    "$BZIP3" -t "$packed" || fail "$name: -t rejects a valid file"
    "$BZIP3" -dc "$packed" | cmp -s - "$orig" || fail "$name: decoded data differs"
    "$BZIP3" -d "$packed" "$name.out" && cmp -s "$name.out" "$orig" || fail "$name: decoding to a file"
done

# Round trips, serial and (where available) threaded, through files and through pipes.
jobs=""
"$BZIP3" --help | grep -q -- "--jobs" && jobs="-j 3"
cat "$DATA"/*.bin "$DATA"/*.txt "$DATA"/*.bz3 "$DATA"/*.bin > big
for _ in 1 2 3 4 5 6; do cat big big > big2 && mv big2 big; done
for flags in "-b 1" "-b 1 $jobs" "-b 3 $jobs" ""; do
    # shellcheck disable=SC2086
    "$BZIP3" -e $flags -c big > big.bz3 || fail "encode ($flags)"
    "$BZIP3" -t big.bz3 || fail "test ($flags)"
    # shellcheck disable=SC2086
    "$BZIP3" -d $jobs -c big.bz3 | cmp -s - big || fail "round trip ($flags)"
    "$BZIP3" -dc < big.bz3 | cmp -s - big || fail "round trip through a pipe ($flags)"
done

# "-" stands for the standard streams; bz3cat and bz3grep rely on it.
"$BZIP3" -c - < big | "$BZIP3" -dc - | cmp -s - big || fail "round trip with - as the file name"
"$BZIP3" -Bcd - < big.bz3 | cmp -s - big || fail "batch mode with - as the file name"
"$BZIP3" -d big.bz3 - | cmp -s - big || fail "- as the output file name"
[ ! -e ./- ] || fail "created a file called -"

# Options are validated.
"$BZIP3" --help=x > /dev/null 2>&1 && fail "accepted an argument to --help"
"$BZIP3" --nonsense > /dev/null 2>&1 && fail "accepted an unknown option"
"$BZIP3" -b > /dev/null 2>&1 < /dev/null && fail "accepted -b without an argument"

# Empty input.
: > empty
"$BZIP3" -c empty > empty.bz3 || fail "encoding an empty file"
"$BZIP3" -dc empty.bz3 | cmp -s - empty || fail "decoding an empty file"

# Existing files are not overwritten unless asked to.
cp "$DATA/plain.txt" plain.txt
echo keep > plain.txt.bz3
"$BZIP3" plain.txt 2> /dev/null && fail "overwrote an existing file"
[ "$(cat plain.txt.bz3)" = keep ] || fail "existing output was modified"
"$BZIP3" -f plain.txt || fail "-f does not overwrite"
"$BZIP3" -dc plain.txt.bz3 | cmp -s - plain.txt || fail "-f produced bad output"

# --rm removes the input once the operation succeeded...
"$BZIP3" -f --rm plain.txt || fail "--rm failed"
[ ! -e plain.txt ] || fail "--rm kept the input"
"$BZIP3" -d --rm plain.txt.bz3 || fail "-d --rm failed"
[ -e plain.txt ] && [ ! -e plain.txt.bz3 ] || fail "-d --rm left the wrong files"
cmp -s plain.txt "$DATA/plain.txt" || fail "--rm round trip"

# ... and never after a failure, which must not leave a partial output file behind either.
"$BZIP3" plain.txt
size=$(wc -c < plain.txt.bz3)
head -c $((size - 7)) plain.txt.bz3 > truncated.bz3
"$BZIP3" -d --rm truncated.bz3 2> /dev/null && fail "accepted a truncated file"
[ -e truncated.bz3 ] || fail "--rm removed the input after a failure"
[ ! -e truncated ] || fail "left a partial output file behind"

# Damaged data is detected.
{ head -c 40 plain.txt.bz3; printf 'XXXX'; tail -c +45 plain.txt.bz3; } > damaged.bz3
"$BZIP3" -t damaged.bz3 2> /dev/null && fail "-t accepted damaged data"
"$BZIP3" -dc damaged.bz3 > /dev/null 2>&1 && fail "-d accepted damaged data"
printf 'BZ2v1xxxxxxxxxxxxxxxxx' | "$BZIP3" -dc > /dev/null 2>&1 && fail "accepted a bad signature"
printf 'BZ3v1\377\377\377\177' | "$BZIP3" -dc > /dev/null 2>&1 && fail "accepted a bad block size"

# Batch mode.
mkdir batch && cp "$DATA/plain.txt" batch/a && cp "$DATA/rle.bin" batch/b
"$BZIP3" -B batch/a batch/b || fail "batch encode"
rm batch/a batch/b
"$BZIP3" -B -d batch/a.bz3 batch/b.bz3 || fail "batch decode"
cmp -s batch/a "$DATA/plain.txt" && cmp -s batch/b "$DATA/rle.bin" || fail "batch round trip"
"$BZIP3" -B -t batch/a.bz3 batch/b.bz3 || fail "batch test"

# Invalid arguments.
for arg in "-b 0" "-b 512" "-b 4097" "-b 99999999999999999999" "-b x" "-b -1"; do
    # shellcheck disable=SC2086
    "$BZIP3" $arg -c plain.txt > /dev/null 2>&1 && fail "accepted $arg"
done

if [ "$failures" -ne 0 ]; then
    echo "$failures test(s) failed" >&2
    exit 1
fi
echo "all CLI tests passed"
