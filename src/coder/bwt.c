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

#include "coder/bwt.h"

static void * libsais_create_ctx(void);
static void libsais_free_ctx(void * ctx);
static int32_t libsais(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_gsa(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_int(int32_t * T, int32_t * SA, int32_t n, int32_t k, int32_t fs);
static int32_t libsais_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs,
                           int32_t * freq);
static int32_t libsais_gsa_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs,
                               int32_t * freq);
static int32_t libsais_bwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_bwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs,
                               int32_t * freq, int32_t r, int32_t * I);
static int32_t libsais_bwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs,
                                  int32_t * freq);
static int32_t libsais_bwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                               int32_t fs, int32_t * freq);
static int32_t libsais_bwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                   int32_t fs, int32_t * freq, int32_t r, int32_t * I);
static int32_t libsais_bwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                      int32_t fs, int32_t * freq);
static void * libsais_unbwt_create_ctx(void);
static void libsais_unbwt_free_ctx(void * ctx);
static int32_t libsais_unbwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq,
                             int32_t i);
static int32_t libsais_unbwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                 const int32_t * freq, int32_t i);
static int32_t libsais_unbwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq,
                                 int32_t r, const int32_t * I);
static int32_t libsais_unbwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                     const int32_t * freq, int32_t r, const int32_t * I);
static int32_t libsais_unbwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                    const int32_t * freq, int32_t i);
static int32_t libsais_unbwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A,
                                        int32_t n, const int32_t * freq, int32_t i);
static int32_t libsais_plcp(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_plcp_gsa(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_plcp_int(const int32_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_lcp(const int32_t * PLCP, const int32_t * SA, int32_t * LCP, int32_t n);

#include "libsais.c"

i32 bz3_bwt(const u8 * in, u8 * out, i32 * sa, i32 n, i32 fs) {
  return libsais_bwt(in, out, sa, n, fs, NULL);
}

i32 bz3_unbwt(const u8 * in, u8 * out, i32 * sa, i32 n, i32 idx) {
  return libsais_unbwt(in, out, sa, n, NULL, idx);
}
