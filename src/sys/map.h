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

#ifndef BZ3_MAP_H
#define BZ3_MAP_H

#include "sys/common.h"
#include <stdio.h>

typedef struct {
  FILE * f;
  int mapped, failed;
  uint64_t size, pos;
  u8 * win;
  uint64_t win_off;
  sz win_len;
  void * handle;
} bz3_in;

void bz3_in_open(bz3_in * in, FILE * f);
sz bz3_in_read(bz3_in * in, void * data, sz n);
void bz3_in_close(bz3_in * in);

void bz3_map_guard(const char * const * doomed);

#endif
