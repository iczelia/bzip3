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

#include "coder/rle.h"

/*  A bitmap of the byte values that are followed by a run length.  */
#define RLE_HEADER 32

int bz3_rle_decode(const u8 * restrict in, i32 n, u8 * restrict out,
                   i32 outlen) {
  i32 op = 0, ip = RLE_HEADER;

  if (n < RLE_HEADER) return 1;

  while (op < outlen && ip < n) {
    const u8 c = in[ip++];
    if (in[c >> 3] & (1 << (c & 7))) {
      const i32 max_run = outlen - op;
      i32 run = 1;
      u8 b;
      do {
        if (UNLIKELY(ip >= n)) return 1;
        b = in[ip++];
        if (LIKELY(run <= max_run)) run += b;
      } while (b == 255);
      if (run > max_run) run = max_run;
      memset(out + op, c, (sz) run);
      op += run;
    } else out[op++] = c;
  }

  return op != outlen;
}
