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

#ifndef BZ3_META_H
#define BZ3_META_H

#include "sys/common.h"
#include <stdio.h>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#endif

typedef struct {
  int known, regular;
#ifdef _WIN32
  FILETIME created, accessed, written;
  DWORD attrs, size_low, size_high;
#else
  struct stat st;
#endif
} bz3_meta;

enum {
  BZ3_META_OWNER = 1, BZ3_META_MODE = 2, BZ3_META_ACL = 4, BZ3_META_TIMES = 8, BZ3_META_ATTRS = 16
};

void bz3_meta_get(bz3_meta * m, FILE * f);

int bz3_meta_changed(const bz3_meta * m, FILE * f);

int bz3_meta_put(const bz3_meta * m, FILE * in, FILE * out, const char * in_path, const char * out_path);

#endif
