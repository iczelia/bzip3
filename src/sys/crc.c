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

#include "sys/crc.h"
#include "sys/cpu.h"

#if defined(HAVE_SSE42)
#define CRC_HW      bz3_crc32c_sse42
#define CRC_HW_NAME "sse4.2"
#elif defined(HAVE_ARM_CRC32)
#define CRC_HW      bz3_crc32c_arm
#define CRC_HW_NAME "armv8-crc"
#endif

#define POLY 0x82F63B78u

#ifdef CRC_HW

static u32 mulmodp(u32 a, u32 b) {
  u32 p = 0, m;
  for (m = (u32) 1 << 31; m; m >>= 1) {
    if (a & m) p ^= b;
    b = b >> 1 ^ (POLY & (0u - (b & 1)));
  }
  return p;
}

static u32 xpow(uint64_t n) {
  u32 p = (u32) 1 << 31, s = (u32) 1 << 30;
  for (; n; n >>= 1, s = mulmodp(s, s)) if (n & 1) p = mulmodp(p, s);
  return p;
}

static void zeros(u32 z[4][256], sz len) {
  u32 one[32];
  int i, j, k;
  one[31] = xpow((uint64_t) len * 8);
  for (i = 31; i; i--) one[i - 1] = one[i] >> 1 ^ (POLY & (0u - (one[i] & 1)));
  Fk(4,
    for (z[k][0] = 0, j = 1; j < 256; j++) {
      for (i = 0; !(j >> i & 1); i++);
      z[k][j] = z[k][j & (j - 1)] ^ one[8 * k + i];
    });
}

#endif

static void slices(u32 T[8][256]) {
  u32 i, j, c;
  Fi(256,
    c = i;
    Fj(8, c = c >> 1 ^ (POLY & (0u - (c & 1))));
    T[0][i] = c);
  Fi(256,
    for (j = 1; j < 8; j++)
      T[j][i] = T[j - 1][i] >> 8 ^ T[0][T[j - 1][i] & 0xFF]);
}

u32 bz3_crc32c_sw(const bz3_crc * t, u32 c, const void * buf, sz n) {
  const u32 (* T)[256] = t->u.slice;
  const u8 * p = buf;
  for (; n >= 8; n -= 8, p += 8) {
    u32 a = c ^ rd32le(p);
    c = T[7][a & 0xFF] ^ T[6][a >> 8 & 0xFF] ^ T[5][a >> 16 & 0xFF] ^
        T[4][a >> 24] ^ T[3][p[4]] ^ T[2][p[5]] ^ T[1][p[6]] ^ T[0][p[7]];
  }
  while (n--) c = c >> 8 ^ T[0][(c ^ *p++) & 0xFF];
  return c;
}

int bz3_crc_impl(bz3_crc * t, int i) {
  if (!i) {
    slices(t->u.slice);
    t->fn = bz3_crc32c_sw;  t->name = "portable";
    return 1;
  }
#ifdef CRC_HW
  if (i == 1 && bz3_cpu_crc()) {
    zeros(t->u.z.zlong, CRC_LONG);
    zeros(t->u.z.zshort, CRC_SHORT);
    t->fn = CRC_HW;  t->name = CRC_HW_NAME;
    return 1;
  }
#endif
  return 0;
}

void bz3_crc_init(bz3_crc * t) {
  if (!bz3_crc_impl(t, 1)) bz3_crc_impl(t, 0);
}
