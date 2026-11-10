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

#ifndef BZ3_BWT_H
#define BZ3_BWT_H

#include "sys/common.h"

i32 bz3_bwt(const u8 * in, u8 * out, i32 * sa, i32 n, i32 fs);
i32 bz3_unbwt(const u8 * in, u8 * out, i32 * sa, i32 n, i32 idx);

#endif
