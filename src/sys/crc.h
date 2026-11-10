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

#ifndef BZ3_CRC_H
#define BZ3_CRC_H

#include "sys/common.h"

typedef struct bz3_crc bz3_crc;
typedef u32 (* bz3_crc_fn)(const bz3_crc * t, u32 crc, const void * p, sz n);

#define CRC_LONG  8192
#define CRC_SHORT 256

struct bz3_crc {
  bz3_crc_fn fn;
  const char * name;
  union {
    u32 slice[8][256];                              /*  Slicing-by-8.  */
    struct { u32 zlong[4][256], zshort[4][256]; } z;  /*  Advance over CRC_LONG and CRC_SHORT zero bytes.  */
  } u;
};

/*  Prepare implementation `i': 0 is the portable one, 1 the one using the instructions of the CPU.  Returns
    zero if it was not built or the CPU can not run it.  */
int bz3_crc_impl(bz3_crc * t, int i);

/*  Prepare the fastest implementation the CPU can run.  */
void bz3_crc_init(bz3_crc * t);

static INLINE u32 bz3_crc32c(const bz3_crc * t, u32 crc, const void * p, sz n) {
  return t->fn(t, crc, p, n);
}

u32 bz3_crc32c_sw(const bz3_crc * t, u32 crc, const void * p, sz n);
u32 bz3_crc32c_sse42(const bz3_crc * t, u32 crc, const void * p, sz n);
u32 bz3_crc32c_arm(const bz3_crc * t, u32 crc, const void * p, sz n);

static INLINE u32 bz3_crc_shift(const u32 z[4][256], u32 c) {
  return z[0][c & 0xFF] ^ z[1][c >> 8 & 0xFF] ^ z[2][c >> 16 & 0xFF] ^
         z[3][c >> 24];
}

#endif
