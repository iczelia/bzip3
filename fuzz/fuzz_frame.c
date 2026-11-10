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

/*  The frame decoder on untrusted data.  The input is the size of the output buffer and then the frame.  */

#include "harness.h"

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size) {
  u8 * out;
  sz n, cap;
  if (size < 4) return 0;
  cap = n = rd32le(data) % 300000;
  out = fz_alloc(cap);
  if (bz3_decompress(data + 4, out, size - 4, &n) == BZ3_OK) fz_touch(out, n);
  FZ_REQUIRE(n <= cap, "decoded size");
  free(out);
  return 0;
}
