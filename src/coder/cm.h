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

#ifndef BZ3_CM_H
#define BZ3_CM_H

#include "sys/common.h"

typedef struct { u16 C0[256], C1[256][256], C2[512][17]; } bz3_cm;

i32 bz3_cm_encode(bz3_cm * s, const u8 * in, i32 n, u8 * out,
                  const u8 * out_end);

void bz3_cm_decode(bz3_cm * s, const u8 * in, const u8 * in_end, u8 * out,
                   i32 n);

#endif
