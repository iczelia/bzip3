/*  Copyright (C) 2022-2026 Kamila Szewczyk

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
    General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/*  Whatever the encoder is given, it stays within bz3_bound and the decoder restores the data, with a
    buffer of bz3_bound bytes and, where the block says it suffices, with one of the original size.  */

#include "harness.h"

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  static struct bz3_state * state;
  sz cap;
  u8 * a, * b;
  i32 packed, r;
  if (size > FZ_BLOCK) size = FZ_BLOCK;
  cap = bz3_bound(size);
  if (!state) state = bz3_new(FZ_BLOCK);
  FZ_REQUIRE(state, "no state");
  a = fz_alloc(cap);  b = fz_alloc(cap);
  memcpy(a, data, size);
  packed = bz3_encode_block(state, a, (i32) size);
  FZ_REQUIRE(packed > 0 && (sz) packed <= cap, "encoded size");
  FZ_REQUIRE(packed < 9 || (i32) rd32le(a + 4) == -1 || !(a[8] & ~2), "model");
  memcpy(b, a, (sz) packed);
  r = bz3_decode_block(state, b, cap, packed, (i32) size);
  FZ_REQUIRE(r == (i32) size && !memcmp(b, data, size), "round trip");
  if (packed >= 9 && size >= (sz) packed &&
      bz3_orig_size_sufficient_for_decode(a, (sz) packed, (i32) size) == 1) {
    u8 * c = fz_alloc(size);
    memcpy(c, a, (sz) packed);
    r = bz3_decode_block(state, c, size, packed, (i32) size);
    FZ_REQUIRE(r == (i32) size && !memcmp(c, data, size), "round trip in place");
    free(c);
  }
  free(a);  free(b);
  return 0;
}
