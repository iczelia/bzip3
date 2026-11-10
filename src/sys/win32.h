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

#ifndef BZ3_WIN32_H
#define BZ3_WIN32_H

#include "sys/common.h"
#include "sys/meta.h"
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <stdarg.h>
#include <stdio.h>

#ifndef BZ3_WIN_LEGACY
wchar_t * bz3_win_wide(const char * s);
char * bz3_win_utf8(const wchar_t * s);
wchar_t * bz3_win_path(const char * s);
char ** bz3_win_args(int argc, wchar_t ** wargv);
#endif

FILE * bz3_win_fopen(const char * path, const char * mode);
int bz3_win_create(const char * path, int force);
int bz3_win_remove(const char * path);

int bz3_win_is_dir(const char * path);
int bz3_win_same_file(FILE * f, const char * path);

void bz3_win_meta_get(bz3_meta * m, FILE * f);
int bz3_win_meta_put(const bz3_meta * m, const char * in, const char * out);

void bz3_win_vprintf(FILE * to, const char * fmt, va_list ap);
#endif

#endif
