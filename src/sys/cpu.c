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

#include "sys/cpu.h"

#if defined(HAVE_SSE42) && (defined(__i386__) || defined(__x86_64__))
#include <cpuid.h>
#endif
#if defined(HAVE_ARM_CRC32) && !defined(__ARM_FEATURE_CRC32) \
    && defined(HAVE_GETAUXVAL) && defined(HAVE_SYS_AUXV_H)
#include <sys/auxv.h>
#ifndef HWCAP_CRC32
#define HWCAP_CRC32 (1UL << 7)
#endif
#endif

#if defined(HAVE_SSE42) && (defined(__i386__) || defined(__x86_64__))
int bz3_cpu_crc(void) {
  unsigned a, b, c, d;
  if (!__get_cpuid(1, &a, &b, &c, &d)) return 0;
  return c >> 20 & 1;
}
#elif defined(HAVE_ARM_CRC32) && defined(__ARM_FEATURE_CRC32)
int bz3_cpu_crc(void) { return 1; }
#elif defined(HAVE_ARM_CRC32) && defined(HAVE_GETAUXVAL) && \
      defined(HAVE_SYS_AUXV_H)
int bz3_cpu_crc(void) { return (getauxval(AT_HWCAP) & HWCAP_CRC32) != 0; }
#else
int bz3_cpu_crc(void) { return 0; }
#endif
