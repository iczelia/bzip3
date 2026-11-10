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

#ifndef BZ3_T_HARNESS_H
#define BZ3_T_HARNESS_H

#include "sys/common.h"
#include "libbz3.h"
#include <stdio.h>

extern int xt_level;
extern int xt_tracing;
extern unsigned long xt_checks, xt_failures;
extern const char * xt_section;

extern const char * xt_data;
extern const char * xt_binary;

void xt_init(void);

void xt_section_begin(const char * name);
void xt_trace(const char * fmt, ...) BZ3_PRINTF(1, 2);
void xt_report(int ok, const char * fmt, ...) BZ3_PRINTF(2, 3);

#define CHECK(cond, ...)  xt_report((cond), __VA_ARGS__)

int xt_finish(const char * program);

void * xt_alloc(sz n);
u8 * xt_slurp(const char * path, sz * len);
void xt_spew(const char * path, const u8 * data, sz len);

const char * xt_fixture(const char * name);
const char * xt_tmp(const char * tag);
void xt_unlink(const char * path);

int xt_same_file(const char * a, const char * b);
long xt_file_size(const char * path);
int xt_file_contains(const char * path, const char * needle);

/*  Run the bzip3 binary with `args', a printf format, and send what it prints to `out'.  Returns the exit
    status.  */
int xt_run(const char * out, const char * args, ...) BZ3_PRINTF(2, 3);

typedef struct { u32 s; } xt_rng;
void xt_seed(xt_rng * r, u32 seed);
u32 xt_next(xt_rng * r, u32 n);

/*  Fill a buffer with data of a kind: 0 noise, 1 text, 2 long runs, 3 long repeats, 4 a mix of them.  */
void xt_fill(xt_rng * r, u8 * buf, sz n, int kind);

/*  The reference files: the original and what bzip3 1.5.5 made of it, with the model byte of its only block,
    -1 for a block stored verbatim.  */
typedef struct { const char * name;  int model; } xt_golden;
extern const xt_golden XT_GOLDEN[];

void xt_run_unit(void);
void xt_run_api(void);
void xt_run_files(void);
void xt_run_cli(void);

#endif
