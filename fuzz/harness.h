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

#ifndef BZ3_FUZZ_HARNESS_H
#define BZ3_FUZZ_HARNESS_H

#include "sys/common.h"
#include "libbz3.h"
#include <stdio.h>

int LLVMFuzzerTestOneInput(const uint8_t * data, size_t size);

#define FZ_BLOCK KiB(65)

/*  A finding that no sanitizer would report.  */
#define FZ_REQUIRE(cond, what) \
  do { if (!(cond)) { fprintf(stderr, "fuzz: %s\n", what);  abort(); } } while (0)

static void * fz_alloc(sz n) {
  void * p = malloc(n ? n : 1);
  FZ_REQUIRE(p, "out of memory");
  return p;
}

/*  Read every byte, so that MemorySanitizer sees what was left uninitialised.  */
static void fz_touch(const u8 * p, sz n) {
  volatile u8 x = 0;
  sz i;
  Fi(n, x ^= p[i]);
  (void) x;
}

#endif
