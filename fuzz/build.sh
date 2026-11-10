#!/bin/bash
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

#  Build the libFuzzer harnesses with non-recovering ASan and UBSan.  Configure first: BUILD names the directory
#  holding config.h.  FZ_SAN picks asan (the default), msan, none, or standalone for a replay tool without
#  libFuzzer; the harnesses are block, frame and rt.
#
#    fuzz/build.sh && fuzz/fuzz_block -max_len=4096 corpus/

set -eu
cd "$(dirname "$0")" || exit 2
TOP=$(cd .. && pwd)
BUILD=${BUILD:-$TOP}
OUT=${OUT:-$TOP/fuzz}
read -r -a CC_CMD <<< "${CC:-clang}"
read -r -a LD_CMD <<< "${LDFLAGS:-}"
case ${FZ_SAN:-asan} in
  asan) SAN=(-fsanitize=address,undefined -fno-sanitize-recover=all)
        OPT=-O1;  COV=-fsanitize=fuzzer-no-link;  LINK=(-fsanitize=fuzzer) ;;
  msan) SAN=(-fsanitize=memory);  OPT=-O1;  COV=-fsanitize=fuzzer-no-link;  LINK=(-fsanitize=fuzzer) ;;
  none) SAN=();  OPT=-O2;  COV=-fsanitize=fuzzer-no-link;  LINK=(-fsanitize=fuzzer) ;;
  standalone) SAN=();  OPT=-O2;  COV=;  LINK=("$TOP/fuzz/fzmain.c") ;;
  *) echo "fuzz/build.sh: FZ_SAN is asan, msan, none or standalone" >&2;  exit 2 ;;
esac
FLAGS=(-g "$OPT" -std=c99 -fno-delete-null-pointer-checks -DHAVE_CONFIG_H
       "-I$BUILD" "-I$TOP/include" "-I$TOP/src" "${SAN[@]}")
test -f "$BUILD/config.h" ||
  { echo "fuzz/build.sh: no config.h in $BUILD; configure first" >&2;  exit 2; }
HARNESSES=${*:-"block frame rt"}
OBJ=$OUT/obj
mkdir -p "$OBJ"

objs=()
build() {               # object, flags..., source
  local o=$1;  shift
  objs+=("$o")
  "${CC_CMD[@]}" "${FLAGS[@]}" $COV "$@" -c -o "$o"
}
for s in libbz3 sys/cpu sys/crc coder/cm coder/lzp coder/rle; do
  build "$OBJ/$(echo "$s" | tr / _).o" "$TOP/src/$s.c"
done
build "$OBJ/coder_bwt.o" "-I$TOP/3rdparty/libsais" "$TOP/src/coder/bwt.c"
if grep -q '^#define HAVE_SSE42' "$BUILD/config.h"; then
  build "$OBJ/crcx86.o" -msse4.2 "$TOP/src/sys/crcx86.c"
fi
if grep -q '^#define HAVE_ARM_CRC32' "$BUILD/config.h"; then
  build "$OBJ/crcarm.o" -march=armv8-a+crc "$TOP/src/sys/crcarm.c"
fi

for t in $HARNESSES; do
  "${CC_CMD[@]}" "${FLAGS[@]}" "${LINK[@]}" "-I$TOP/fuzz" \
    -o "$OUT/fuzz_$t" "$TOP/fuzz/fuzz_$t.c" "${objs[@]}" "${LD_CMD[@]}" -lpthread
done
echo "built: $HARNESSES"
