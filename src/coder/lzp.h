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

#ifndef BZ3_LZP_H
#define BZ3_LZP_H

#include "sys/common.h"

/*  Lempel-Ziv prediction: a hash of the last four bytes predicts where the data continues, and matches of at
    least LZP_MIN_MATCH bytes at the predicted position are replaced by their length.  */

#define LZP_DICT_BITS 18
#define LZP_DICT ((sz) 1 << LZP_DICT_BITS)

/*  Code `n' bytes of `in' into at most `n' bytes of `out'.  `lut' holds LZP_DICT entries.  Returns the size of
    the output, or -1 if the input is too short or would not shrink.  */
i32 bz3_lzp_encode(const u8 * in, u8 * out, i32 n, i32 * lut);

/*  Decode `n' bytes of `in' into at most `max' bytes of `out'.  Returns the size of the output or -1.  */
i32 bz3_lzp_decode(const u8 * in, u8 * out, i32 n, i32 max, i32 * lut);

#endif
