
/*
 * BZip3 - A spiritual successor to BZip2.
 * Copyright (C) 2022-2026 Kamila Szewczyk
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BZ3_SAIS_H
#define BZ3_SAIS_H

/* Brings the bundled, unmodified libsais (3rdparty/libsais) into the translation unit that includes this
   file, without making any of it visible to the outside. */

#include <stdint.h>

static void * libsais_create_ctx(void);
static void libsais_free_ctx(void * ctx);
static int32_t libsais(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_gsa(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_int(int32_t * T, int32_t * SA, int32_t n, int32_t k, int32_t fs);
static int32_t libsais_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_gsa_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs,
                               int32_t * freq);
static int32_t libsais_bwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_bwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq,
                               int32_t r, int32_t * I);
static int32_t libsais_bwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq);
static int32_t libsais_bwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs,
                               int32_t * freq);
static int32_t libsais_bwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs,
                                   int32_t * freq, int32_t r, int32_t * I);
static int32_t libsais_bwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                      int32_t fs, int32_t * freq);
static void * libsais_unbwt_create_ctx(void);
static void libsais_unbwt_free_ctx(void * ctx);
static int32_t libsais_unbwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i);
static int32_t libsais_unbwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                 const int32_t * freq, int32_t i);
static int32_t libsais_unbwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq,
                                 int32_t r, const int32_t * I);
static int32_t libsais_unbwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                     const int32_t * freq, int32_t r, const int32_t * I);
static int32_t libsais_unbwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq,
                                    int32_t i);
static int32_t libsais_unbwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n,
                                        const int32_t * freq, int32_t i);
static int32_t libsais_plcp(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_plcp_gsa(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_plcp_int(const int32_t * T, const int32_t * SA, int32_t * PLCP, int32_t n);
static int32_t libsais_lcp(const int32_t * PLCP, const int32_t * SA, int32_t * LCP, int32_t n);

#if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-function"
    /* GCC cannot tell that a transform with a single primary index never takes the multi-index paths. */
    #pragma GCC diagnostic ignored "-Warray-bounds"
#elif defined(_MSC_VER)
    #pragma warning(push)
    #pragma warning(disable : 4505)
#endif

#include "libsais.c"

#if defined(__GNUC__) || defined(__clang__)
    #pragma GCC diagnostic pop
#elif defined(_MSC_VER)
    #pragma warning(pop)
#endif

#endif
