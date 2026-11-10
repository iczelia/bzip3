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
#include <arm_acle.h>

typedef u32 crc_w;
#define CRC_W 8
static INLINE crc_w step(crc_w c, const u8 * p) {
  uint64_t w;  memcpy(&w, p, 8);  return __crc32cd(c, w);
}

#define CRC_FN bz3_crc32c_arm
#define CRC_BYTE __crc32cb
#include "sys/crchw.h"
