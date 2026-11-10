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

/*  The block decoder on untrusted data.  The input is three 32-bit numbers, the original size, the compressed
    size and the size of the buffer, and then the block.  The compressed size is that of the block unless the
    top bit of the original size asks for the number as given.  */

#include "harness.h"

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  static struct bz3_state * state;
  u32 orig, comp, cap;
  u8 * buf;
  i32 r;
  if (size < 12) return 0;
  orig = rd32le(data);  comp = rd32le(data + 4);  cap = rd32le(data + 8) % 200000;
  data += 12;  size -= 12;
  if (!(orig >> 31)) comp = (u32) size;
  orig &= 0x7fffffff;
  if (!state) state = bz3_new(FZ_BLOCK);
  FZ_REQUIRE(state, "no state");
  buf = fz_alloc(cap);
  memcpy(buf, data, MIN(size, (sz) cap));
  r = bz3_decode_block(state, buf, cap, (i32) comp, (i32) orig);
  if (r >= 0) {
    FZ_REQUIRE(r == (i32) orig && (u32) r <= cap, "decoded size");
    fz_touch(buf, (sz) r);
  }
  bz3_orig_size_sufficient_for_decode(buf, cap, (i32) orig);
  free(buf);
  return 0;
}
