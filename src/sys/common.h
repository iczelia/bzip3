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

#ifndef BZ3_COMMON_H
#define BZ3_COMMON_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

typedef uint8_t  u8;  typedef int8_t  i8;
typedef uint16_t u16;  typedef int16_t i16;
typedef uint32_t u32;  typedef int32_t i32;
typedef size_t   sz;

#if CHAR_BIT != 8 || INT_MAX < 2147483647
#error "bzip3 requires 8-bit bytes and an int of at least 32 bits"
#endif
#if (-1 >> 1) != -1
#error "bzip3 requires arithmetic signed right shift"
#endif

#if defined(__GNUC__) && !defined(BZ3_NO_ATTRS)
#define INLINE      __inline__ __attribute__((__always_inline__))
#define NORETURN    __attribute__((__noreturn__))
#define LIKELY(x)   __builtin_expect(!!(x), 1)
#define UNLIKELY(x) __builtin_expect(!!(x), 0)
#define BZ3_PRINTF(f, a) __attribute__((__format__(__printf__, f, a)))
#else
#define INLINE
#define NORETURN
#define LIKELY(x)   (x)
#define UNLIKELY(x) (x)
#define BZ3_PRINTF(f, a)
#endif

static INLINE u32 rd32le(const u8 * p) {
  return p[0] | (u32) p[1] << 8 | (u32) p[2] << 16 | (u32) p[3] << 24;
}
static INLINE void wr32le(u8 * p, u32 v) {
  p[0] = (u8) v;  p[1] = (u8) (v >> 8);  p[2] = (u8) (v >> 16);  p[3] = (u8) (v >> 24);
}

#define Fi(n, ...) for (i = 0; i < (n); i++) {  __VA_ARGS__; }
#define Fj(n, ...) for (j = 0; j < (n); j++) {  __VA_ARGS__; }
#define Fk(n, ...) for (k = 0; k < (n); k++) {  __VA_ARGS__; }

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ARRAY_LEN(a) (sizeof (a) / sizeof *(a))

#define KiB(x) ((x) * 1024)
#define MiB(x) ((x) * 1024 * 1024)

/*  Limits of the block size, shared by the library and the tool.  */
#define BZ3_MIN_BLOCK KiB(65)
#define BZ3_MAX_BLOCK MiB(511)

#endif
