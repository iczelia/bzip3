/*--

This file is a part of libsais, a library for linear time suffix array,
longest common prefix array and burrows wheeler transform construction.

   Copyright (c) 2021-2026 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 with LLVM Exceptions
   (the "License"); you may not use this file except in compliance with
   the License. You may obtain a copy of the License in the file LICENSE
   distributed with this software or at

       https://llvm.org/LICENSE.txt

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

Please see the file LICENSE for full copyright information.

--*/

#include "libsais.h"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#if defined(LIBSAIS_OPENMP)
    #include <omp.h>
#else
    #define UNUSED(_x)                  (void)(_x)
#endif

#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ >= 16)
    #pragma GCC diagnostic ignored "-Wdeprecated-openmp"
#endif

typedef int32_t                         sa_sint_t;
typedef uint32_t                        sa_uint_t;
typedef ptrdiff_t                       fast_sint_t;
typedef size_t                          fast_uint_t;

#define SAINT_BIT                       (32)
#define SAINT_MAX                       INT32_MAX
#define SAINT_MIN                       INT32_MIN

#define ALPHABET_SIZE                   (1 << CHAR_BIT)
#define UNBWT_FASTBITS                  (17)
#define UNBWT_CHAINS                    (24)
#define UNBWT_CHAINS_LARGE              (64)
#define UNBWT_CHAINS_LARGE_MIN_SIZE     (16777216)
#define UNBWT_CHAIN_MARKS               (4096)
#define UNBWT_CHAIN_MIN_SIZE            (65536)
#define UNBWT_CHAIN_THREADS_MIN_SIZE    (262144)
#define UNBWT_CHAIN_MIN_SIZE_AUX        (1048576)
#define UNBWT_CHAIN_ANY_BLOCKS_MIN_SIZE (33554432)
#define UNBWT_CHUNK_SIZE                (16320)
#define UNBWT_CHUNK_SIZE_SMALL          (4080)
#define UNBWT_DECODE_12_MIN_SIZE        (262144)
#define UNBWT_DECODE_16_MIN_SIZE        (8388608)
#define UNBWT_PACKED_MAX_SIZE           (65535)
#define UNBWT_SYMBOL_ROWS_MAX_SIZE      (8388607)
#define UNBWT_FIRST_SYMBOLS             (15)
#define UNBWT_MARK                      (((sa_uint_t)1) << (SAINT_BIT - 1))

#define SUFFIX_GROUP_BIT                (SAINT_BIT - 1)
#define SUFFIX_GROUP_MARKER             (((sa_sint_t)1) << (SUFFIX_GROUP_BIT - 1))

#define BUCKETS_INDEX2(_c, _s)          (((_c) << 1) + (_s))
#define BUCKETS_INDEX4(_c, _s)          (((_c) << 2) + (_s))

#define LIBSAIS_LOCAL_BUFFER_SIZE       (2000)
#define LIBSAIS_PER_THREAD_CACHE_SIZE   (24576)

#define LIBSAIS_NONUNIQUE_INSERTION_SORT_SIZE   (16)
#define LIBSAIS_NONUNIQUE_PREFETCH_GROUPS       (16)
#define LIBSAIS_NONUNIQUE_PREFETCH_MEMBERS      (32)
#define LIBSAIS_NONUNIQUE_DEFERRED_STEPS        (16)
#define LIBSAIS_NONUNIQUE_BUDGET_STEPS          (4)
#define LIBSAIS_NONUNIQUE_SAMPLES               (512)

#define LIBSAIS_FLAGS_NONE              (0)
#define LIBSAIS_FLAGS_BWT               (1)
#define LIBSAIS_FLAGS_GSA               (2)
#define LIBSAIS_FLAGS_CYCLIC            (4)

typedef struct LIBSAIS_THREAD_CACHE
{
        sa_sint_t                       symbol;
        sa_sint_t                       index;
} LIBSAIS_THREAD_CACHE;

typedef struct LIBSAIS_NONUNIQUE_SUFFIX
{
        sa_sint_t                       name;
        sa_sint_t                       cursor;
        sa_sint_t                       index;
} LIBSAIS_NONUNIQUE_SUFFIX;

typedef union LIBSAIS_THREAD_STATE
{
    struct
    {
        fast_sint_t                     position;
        fast_sint_t                     count;

        fast_sint_t                     m;

        sa_sint_t *                     buckets;
        LIBSAIS_THREAD_CACHE *          cache;
    } state;

    uint8_t padding[64];
} LIBSAIS_THREAD_STATE;

typedef struct LIBSAIS_CONTEXT
{
    sa_sint_t *                         buckets;
    LIBSAIS_THREAD_STATE *              thread_state;
    fast_sint_t                         threads;
} LIBSAIS_CONTEXT;

typedef struct LIBSAIS_UNBWT_CHAIN
{
    sa_uint_t                           row;
    sa_uint_t                           next;
    sa_uint_t                           lf;
    sa_uint_t                           lf_next;
    sa_uint_t                           chunk;
    sa_uint_t                           offset;
    sa_uint_t                           length;
    sa_uint_t                           successor;
    sa_uint_t                           position;
    sa_uint_t                           gap;
} LIBSAIS_UNBWT_CHAIN;

typedef struct LIBSAIS_UNBWT_SLOT
{
    fast_uint_t                         p;
    fast_uint_t                         chain;
    fast_uint_t                         chunk;
    fast_sint_t                         length;
    uint8_t *                           out;
    uint8_t *                           end;
} LIBSAIS_UNBWT_SLOT;

typedef struct LIBSAIS_UNBWT_CONTEXT
{
    sa_uint_t *                         bucket2;
    uint16_t *                          fastbits;
    sa_uint_t *                         buckets;
    sa_uint_t *                         marks;
    LIBSAIS_UNBWT_CHAIN *               chains;
    uint8_t *                           slack;
    sa_uint_t *                         links;
    fast_sint_t                         threads;
} LIBSAIS_UNBWT_CONTEXT;

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER) || defined(__INTEL_COMPILER)
    #define RESTRICT __restrict
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if defined(__GNUC__) || defined(__clang__)
    #define LIBSAIS_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
    #define LIBSAIS_NOINLINE __declspec(noinline)
#else
    #define LIBSAIS_NOINLINE
#endif

#if defined(__has_builtin)
    #if __has_builtin(__builtin_prefetch)
        #define HAS_BUILTIN_PREFETCH
    #endif
#elif defined(__GNUC__) && (((__GNUC__ == 3) && (__GNUC_MINOR__ >= 2)) || (__GNUC__ >= 4))
    #define HAS_BUILTIN_PREFETCH
#endif

#if defined(__has_builtin)
    #if __has_builtin(__builtin_bswap16)
        #define HAS_BUILTIN_BSWAP16
    #endif
#elif defined(__GNUC__) && (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 8)) || (__GNUC__ >= 5))
    #define HAS_BUILTIN_BSWAP16
#endif

#define LIBSAIS_CPU_FEATURE_NONE    0
#define LIBSAIS_CPU_FEATURE_A64     1
#define LIBSAIS_CPU_FEATURE_SSE2    2
#define LIBSAIS_CPU_FEATURE_AVX2    3

#ifndef LIBSAIS_CPU_FEATURE
    #if defined(__AVX2__)
        #define LIBSAIS_CPU_FEATURE LIBSAIS_CPU_FEATURE_AVX2
    #elif (defined(__SSE2__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2) || defined(_M_AMD64) || defined(_M_X64)) && !defined(_M_ARM64EC)
        #define LIBSAIS_CPU_FEATURE LIBSAIS_CPU_FEATURE_SSE2
    #elif (defined(__aarch64__) && !defined(__AARCH64EB__)) || defined(_M_ARM64)
        #define LIBSAIS_CPU_FEATURE LIBSAIS_CPU_FEATURE_A64
    #else
        #define LIBSAIS_CPU_FEATURE LIBSAIS_CPU_FEATURE_NONE
    #endif
#endif

#if (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_NONE) && (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_A64) && (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_SSE2) && (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_AVX2)
    #error Your compiler, configuration or platform is not supported.
#endif

#if LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2
    #if defined(_MSC_VER)
        #include <intrin.h>
    #elif defined(__GNUC__)
        #include <x86intrin.h>
    #endif
#elif LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64
    #include <arm_neon.h>
    #if defined(_MSC_VER) && !defined(__clang__)
        #pragma warning(disable: 4223)
    #endif
#endif

#if defined(HAS_BUILTIN_PREFETCH)
    #define libsais_prefetchr(address) __builtin_prefetch((const void *)(address), 0, 3)
    #define libsais_prefetchw(address) __builtin_prefetch((const void *)(address), 1, 3)
#elif (defined (_M_IX86) || defined (_M_AMD64)) && !defined (_M_ARM64EC)
    #include <intrin.h>
    #define libsais_prefetchr(address) _mm_prefetch((const void *)(address), _MM_HINT_T0)
    #define libsais_prefetchw(address) _m_prefetchw((const void *)(address))
#elif defined (_M_ARM)
    #include <intrin.h>
    #define libsais_prefetchr(address) __prefetch((const void *)(address))
    #define libsais_prefetchw(address) __prefetchw((const void *)(address))
#elif defined (_M_ARM64) || defined (_M_ARM64EC)
    #include <intrin.h>
    #define libsais_prefetchr(address) __prefetch2((const void *)(address), 0)
    #define libsais_prefetchw(address) __prefetch2((const void *)(address), 16)
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if !defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
    #if (defined(_LITTLE_ENDIAN) && !defined(_BIG_ENDIAN)) \
            || (defined(BYTE_ORDER) && defined(LITTLE_ENDIAN) && BYTE_ORDER == LITTLE_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_LITTLE_ENDIAN) && _BYTE_ORDER == _LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__LITTLE_ENDIAN) && __BYTE_ORDER == __LITTLE_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        #define __LITTLE_ENDIAN__
    #elif (defined(_BIG_ENDIAN) && !defined(_LITTLE_ENDIAN)) \
            || (defined(BYTE_ORDER) && defined(BIG_ENDIAN) && BYTE_ORDER == BIG_ENDIAN) \
            || (defined(_BYTE_ORDER) && defined(_BIG_ENDIAN) && _BYTE_ORDER == _BIG_ENDIAN) \
            || (defined(__BYTE_ORDER) && defined(__BIG_ENDIAN) && __BYTE_ORDER == __BIG_ENDIAN) \
            || (defined(__BYTE_ORDER__) && defined(__ORDER_BIG_ENDIAN__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        #define __BIG_ENDIAN__
    #elif defined(_WIN32)
        #define __LITTLE_ENDIAN__
    #else
        #error Your compiler, configuration or platform is not supported.
    #endif
#endif

#if defined(__LITTLE_ENDIAN__) && !defined(__BIG_ENDIAN__)
    #if defined(HAS_BUILTIN_BSWAP16)
        #define libsais_bswap16(x) (__builtin_bswap16(x))
    #elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
        #define libsais_bswap16(x) (_byteswap_ushort(x))
    #else
        #define libsais_bswap16(x) ((uint16_t)(x >> 8) | (uint16_t)(x << 8))
    #endif
#elif !defined(__LITTLE_ENDIAN__) && defined(__BIG_ENDIAN__)
    #define libsais_bswap16(x) (x)
#else
    #error Your compiler, configuration or platform is not supported.
#endif

#if (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_NONE) && defined(__BIG_ENDIAN__)
    #error Your compiler, configuration or platform is not supported.
#endif

static void * libsais_align_up(const void * address, size_t alignment)
{
    return (void *)((((ptrdiff_t)address) + ((ptrdiff_t)alignment) - 1) & (-((ptrdiff_t)alignment)));
}

static void * libsais_alloc_aligned(size_t size, size_t alignment)
{
    void * address = malloc(size + sizeof(short) + alignment - 1);
    if (address != NULL)
    {
        void * aligned_address = libsais_align_up((void *)((ptrdiff_t)address + (ptrdiff_t)(sizeof(short))), alignment);
        ((short *)aligned_address)[-1] = (short)((ptrdiff_t)aligned_address - (ptrdiff_t)address);

        return aligned_address;
    }

    return NULL;
}

static void libsais_free_aligned(void * aligned_address)
{
    if (aligned_address != NULL)
    {
        free((void *)((ptrdiff_t)aligned_address - ((short *)aligned_address)[-1]));
    }
}

static LIBSAIS_THREAD_STATE * libsais_alloc_thread_state(sa_sint_t threads)
{
    if ((size_t)threads > SIZE_MAX / (LIBSAIS_PER_THREAD_CACHE_SIZE * sizeof(LIBSAIS_THREAD_CACHE)) || (fast_sint_t)threads > PTRDIFF_MAX / (32 * (fast_sint_t)threads)) { return NULL; }

    LIBSAIS_THREAD_STATE *  RESTRICT thread_state    = (LIBSAIS_THREAD_STATE *)libsais_alloc_aligned((size_t)threads * sizeof(LIBSAIS_THREAD_STATE), 4096);
    sa_sint_t *             RESTRICT thread_buckets  = (sa_sint_t *)libsais_alloc_aligned((size_t)threads * 4 * ALPHABET_SIZE * sizeof(sa_sint_t), 4096);
    LIBSAIS_THREAD_CACHE *  RESTRICT thread_cache    = (LIBSAIS_THREAD_CACHE *)libsais_alloc_aligned((size_t)threads * LIBSAIS_PER_THREAD_CACHE_SIZE * sizeof(LIBSAIS_THREAD_CACHE), 4096);

    if (thread_state != NULL && thread_buckets != NULL && thread_cache != NULL)
    {
        fast_sint_t t;
        for (t = 0; t < threads; ++t)
        { 
            thread_state[t].state.buckets   = thread_buckets;   thread_buckets  += 4 * ALPHABET_SIZE;
            thread_state[t].state.cache     = thread_cache;     thread_cache    += LIBSAIS_PER_THREAD_CACHE_SIZE;
        }

        return thread_state;
    }

    libsais_free_aligned(thread_cache);
    libsais_free_aligned(thread_buckets);
    libsais_free_aligned(thread_state);
    return NULL;
}

static void libsais_free_thread_state(LIBSAIS_THREAD_STATE * thread_state)
{
    if (thread_state != NULL)
    {
        libsais_free_aligned(thread_state[0].state.cache);
        libsais_free_aligned(thread_state[0].state.buckets);
        libsais_free_aligned(thread_state);
    }
}

static LIBSAIS_CONTEXT * libsais_create_ctx_main(sa_sint_t threads)
{
    LIBSAIS_CONTEXT *       RESTRICT ctx            = (LIBSAIS_CONTEXT *)libsais_alloc_aligned(sizeof(LIBSAIS_CONTEXT), 64);
    sa_sint_t *             RESTRICT buckets        = (sa_sint_t *)libsais_alloc_aligned((size_t)8 * ALPHABET_SIZE * sizeof(sa_sint_t), 4096);
    LIBSAIS_THREAD_STATE *  RESTRICT thread_state   = threads > 1 ? libsais_alloc_thread_state(threads) : NULL;

    if (ctx != NULL && buckets != NULL && (thread_state != NULL || threads == 1))
    {
        ctx->buckets = buckets;
        ctx->threads = threads;
        ctx->thread_state = thread_state;

        return ctx;
    }

    libsais_free_thread_state(thread_state);
    libsais_free_aligned(buckets);
    libsais_free_aligned(ctx);
    return NULL;
}

static void libsais_free_ctx_main(LIBSAIS_CONTEXT * ctx)
{
    if (ctx != NULL)
    {
        libsais_free_thread_state(ctx->thread_state);
        libsais_free_aligned(ctx->buckets);
        libsais_free_aligned(ctx);
    }
}

#if defined(LIBSAIS_OPENMP)

static sa_sint_t libsais_count_negative_marked_suffixes(sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    sa_sint_t count = 0;

    fast_sint_t i; for (i = omp_block_start; i < omp_block_start + omp_block_size; ++i) { count += (SA[i] < 0); }

    return count;
}

static sa_sint_t libsais_count_zero_marked_suffixes(sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    sa_sint_t count = 0;

    fast_sint_t i; for (i = omp_block_start; i < omp_block_start + omp_block_size; ++i) { count += (SA[i] == 0); }

    return count;
}

static void libsais_place_cached_suffixes(sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SA[cache[i + prefetch_distance + 0].symbol]);
        libsais_prefetchw(&SA[cache[i + prefetch_distance + 1].symbol]);
        libsais_prefetchw(&SA[cache[i + prefetch_distance + 2].symbol]);
        libsais_prefetchw(&SA[cache[i + prefetch_distance + 3].symbol]);

        SA[cache[i + 0].symbol] = cache[i + 0].index;
        SA[cache[i + 1].symbol] = cache[i + 1].index;
        SA[cache[i + 2].symbol] = cache[i + 2].index;
        SA[cache[i + 3].symbol] = cache[i + 3].index;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[cache[i].symbol] = cache[i].index;
    }
}

static void libsais_compact_and_place_cached_suffixes(sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3, l = omp_block_start; i < j; i += 4)
    {
        libsais_prefetchw(&cache[i + prefetch_distance]);

        cache[l] = cache[i + 0]; l += cache[l].symbol >= 0;
        cache[l] = cache[i + 1]; l += cache[l].symbol >= 0;
        cache[l] = cache[i + 2]; l += cache[l].symbol >= 0;
        cache[l] = cache[i + 3]; l += cache[l].symbol >= 0;
    }

    for (j += 3; i < j; i += 1)
    {
        cache[l] = cache[i]; l += cache[l].symbol >= 0;
    }

    libsais_place_cached_suffixes(SA, cache, omp_block_start, l - omp_block_start);
}

static void libsais_accumulate_counts_s32_2(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s]; }
}

static void libsais_accumulate_counts_s32_3(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s]; }
}

static void libsais_accumulate_counts_s32_4(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s]; }
}

static void libsais_accumulate_counts_s32_5(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    sa_sint_t * RESTRICT bucket04 = bucket03 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s] + bucket04[s]; }
}

static void libsais_accumulate_counts_s32_6(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    sa_sint_t * RESTRICT bucket04 = bucket03 - bucket_stride;
    sa_sint_t * RESTRICT bucket05 = bucket04 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s] + bucket04[s] + bucket05[s]; }
}

static void libsais_accumulate_counts_s32_7(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    sa_sint_t * RESTRICT bucket04 = bucket03 - bucket_stride;
    sa_sint_t * RESTRICT bucket05 = bucket04 - bucket_stride;
    sa_sint_t * RESTRICT bucket06 = bucket05 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s] + bucket04[s] + bucket05[s] + bucket06[s]; }
}

static void libsais_accumulate_counts_s32_8(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    sa_sint_t * RESTRICT bucket04 = bucket03 - bucket_stride;
    sa_sint_t * RESTRICT bucket05 = bucket04 - bucket_stride;
    sa_sint_t * RESTRICT bucket06 = bucket05 - bucket_stride;
    sa_sint_t * RESTRICT bucket07 = bucket06 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s] + bucket04[s] + bucket05[s] + bucket06[s] + bucket07[s]; }
}

static void libsais_accumulate_counts_s32_9(sa_sint_t * RESTRICT bucket00, fast_sint_t bucket_size, fast_sint_t bucket_stride)
{
    sa_sint_t * RESTRICT bucket01 = bucket00 - bucket_stride;
    sa_sint_t * RESTRICT bucket02 = bucket01 - bucket_stride;
    sa_sint_t * RESTRICT bucket03 = bucket02 - bucket_stride;
    sa_sint_t * RESTRICT bucket04 = bucket03 - bucket_stride;
    sa_sint_t * RESTRICT bucket05 = bucket04 - bucket_stride;
    sa_sint_t * RESTRICT bucket06 = bucket05 - bucket_stride;
    sa_sint_t * RESTRICT bucket07 = bucket06 - bucket_stride;
    sa_sint_t * RESTRICT bucket08 = bucket07 - bucket_stride;
    fast_sint_t s; for (s = 0; s < bucket_size; s += 1) { bucket00[s] = bucket00[s] + bucket01[s] + bucket02[s] + bucket03[s] + bucket04[s] + bucket05[s] + bucket06[s] + bucket07[s] + bucket08[s]; }
}

static void libsais_accumulate_counts_s32(sa_sint_t * RESTRICT buckets, fast_sint_t bucket_size, fast_sint_t bucket_stride, fast_sint_t num_buckets)
{
    while (num_buckets >= 9)
    {
        libsais_accumulate_counts_s32_9(buckets - (num_buckets - 9) * bucket_stride, bucket_size, bucket_stride); num_buckets -= 8;
    }

    switch (num_buckets)
    {
        case 1: break;
        case 2: libsais_accumulate_counts_s32_2(buckets, bucket_size, bucket_stride); break;
        case 3: libsais_accumulate_counts_s32_3(buckets, bucket_size, bucket_stride); break;
        case 4: libsais_accumulate_counts_s32_4(buckets, bucket_size, bucket_stride); break;
        case 5: libsais_accumulate_counts_s32_5(buckets, bucket_size, bucket_stride); break;
        case 6: libsais_accumulate_counts_s32_6(buckets, bucket_size, bucket_stride); break;
        case 7: libsais_accumulate_counts_s32_7(buckets, bucket_size, bucket_stride); break;
        case 8: libsais_accumulate_counts_s32_8(buckets, bucket_size, bucket_stride); break;
        default: break;
    }
}

static fast_sint_t libsais_get_block_size_32s(sa_sint_t n, sa_sint_t threads)
{
    fast_sint_t block_size = (fast_sint_t)n >> 3; if (block_size < 16384) { block_size = 16384; }

    return block_size < (fast_sint_t)threads * LIBSAIS_PER_THREAD_CACHE_SIZE ? block_size : (fast_sint_t)threads * LIBSAIS_PER_THREAD_CACHE_SIZE;
}

static sa_sint_t libsais_get_block_threads_8u(fast_sint_t block_size, sa_sint_t k, sa_sint_t threads)
{
    if (block_size >= 64 * (fast_sint_t)(k > 256 ? k : 256)) { return threads; }

    fast_sint_t max_threads = threads <= 8 ? block_size / 1024 : 1; if (max_threads > threads) { max_threads = threads; }

    return max_threads > 1 ? (sa_sint_t)max_threads : 1;
}

#endif

static void libsais_flip_suffix_markers_omp(sa_sint_t * RESTRICT SA, sa_sint_t l, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && l >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (l / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : l - omp_block_start;

        fast_sint_t i; for (i = omp_block_start; i < omp_block_start + omp_block_size; ++i) { SA[i] ^= SAINT_MIN; }
    }
}

#if LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_NONE

typedef uint16_t libsais_index_t;

#endif

#if LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_AVX2

#define LIBSAIS_LANES (32)

typedef __m256i libsais_lanes_t;

static libsais_lanes_t libsais_lanes_broadcast(fast_sint_t c)
{
    return _mm256_set1_epi8((char)(c));
}

static libsais_lanes_t libsais_load_reversed(const uint8_t * RESTRICT p)
{
    const __m256i reverse_bytes = _mm256_broadcastsi128_si256(_mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15));

    __m256i v = _mm256_load_si256((const __m256i *)(const void *)p);

    return _mm256_shuffle_epi8(_mm256_permute2x128_si256(v, v, 1), reverse_bytes);
}

#define libsais_lanes_shift_in(cur, prev)   _mm256_alignr_epi8((cur), _mm256_permute2x128_si256((prev), (cur), 0x21), 15)
#define libsais_lanes_shift_out(cur, next)  _mm256_alignr_epi8(_mm256_permute2x128_si256((cur), (next), 0x21), (cur), 1)
#define libsais_lanes_type_pair(own, pred)  _mm256_sub_epi8(_mm256_and_si256((own), _mm256_set1_epi8(2)), (pred))

static uint64_t libsais_mask_greater_or_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(_mm256_max_epu8(c0, c1), c0));
}

static uint64_t libsais_mask_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(c0, c1));
}

static void libsais_lanes_from_mask(uint64_t bits, libsais_lanes_t * RESTRICT lanes)
{
    const __m256i weights = _mm256_broadcastsi128_si256(_mm_set_epi8(-128, 64, 32, 16, 8, 4, 2, 1, -128, 64, 32, 16, 8, 4, 2, 1));
    const __m256i bytes0  = _mm256_set_epi8(3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m256i bytes1  = _mm256_set_epi8(7, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4);

    __m256i v = _mm256_set1_epi64x((long long)bits);

    lanes[0] = _mm256_cmpeq_epi8(_mm256_and_si256(_mm256_shuffle_epi8(v, bytes0), weights), weights);
    lanes[1] = _mm256_cmpeq_epi8(_mm256_and_si256(_mm256_shuffle_epi8(v, bytes1), weights), weights);
}

static void libsais_store_indices(libsais_lanes_t c, libsais_lanes_t t, libsais_index_t * RESTRICT index)
{
    _mm256_storeu_si256((__m256i *)(void *)&index[0], _mm256_unpacklo_epi8(c, t));
    _mm256_storeu_si256((__m256i *)(void *)&index[16], _mm256_unpackhi_epi8(c, t));
}

#define LIBSAIS_LMS_MARGIN (36)

#if SAINT_BIT == 32

typedef __m256i libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return _mm256_set1_epi32((int)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    _mm256_storeu_si256((__m256i *)(void *)p, _mm256_sub_epi32(base, _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(const void *)lanes))));

    return _mm256_sub_epi32(base, _mm256_set1_epi32(8));
}

#else

typedef __m256i libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return _mm256_set1_epi64x((long long)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    __m128i v = _mm_loadl_epi64((const __m128i *)(const void *)lanes);

    _mm256_storeu_si256((__m256i *)(void *)&p[0], _mm256_sub_epi64(base, _mm256_cvtepu8_epi64(v)));
    _mm256_storeu_si256((__m256i *)(void *)&p[4], _mm256_sub_epi64(base, _mm256_cvtepu8_epi64(_mm_srli_si128(v, 4))));

    return _mm256_sub_epi64(base, _mm256_set1_epi64x(8));
}

#endif

#if SAINT_BIT == 32

#define LIBSAIS_LANES_32S (32)

typedef struct libsais_lanes_32s_t
{
    __m256i                             v0, v1, v2, v3;
} libsais_lanes_32s_t;

static __m256i libsais_load_reversed_8x32s(const sa_sint_t * RESTRICT p)
{
    return _mm256_shuffle_epi32(_mm256_inserti128_si256(_mm256_castsi128_si256(_mm_load_si128((const __m128i *)(const void *)&p[16])), _mm_load_si128((const __m128i *)(const void *)&p[0]), 1), _MM_SHUFFLE(0, 1, 2, 3));
}

static libsais_lanes_32s_t libsais_lanes_broadcast_32s(fast_sint_t c)
{
    libsais_lanes_32s_t r; r.v0 = r.v1 = r.v2 = r.v3 = _mm256_set1_epi32((int)(c));

    return r;
}

static libsais_lanes_32s_t libsais_load_reversed_32s(const sa_sint_t * RESTRICT p)
{
    libsais_lanes_32s_t r; r.v0 = libsais_load_reversed_8x32s(&p[12]); r.v1 = libsais_load_reversed_8x32s(&p[8]); r.v2 = libsais_load_reversed_8x32s(&p[4]); r.v3 = libsais_load_reversed_8x32s(&p[0]);

    return r;
}

static libsais_lanes_32s_t libsais_lanes_shift_in_32s(libsais_lanes_32s_t cur, libsais_lanes_32s_t prev)
{
    libsais_lanes_32s_t r;

    r.v0 = _mm256_alignr_epi8(cur.v0, _mm256_permute2x128_si256(prev.v3, cur.v3, 0x21), 12);
    r.v1 = _mm256_alignr_epi8(cur.v1, cur.v0, 12);
    r.v2 = _mm256_alignr_epi8(cur.v2, cur.v1, 12);
    r.v3 = _mm256_alignr_epi8(cur.v3, cur.v2, 12);

    return r;
}

static uint64_t libsais_mask_32s(__m256i m0, __m256i m1, __m256i m2, __m256i m3)
{
    return (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_packs_epi16(_mm256_packs_epi32(m0, m1), _mm256_packs_epi32(m2, m3)));
}

static uint64_t libsais_mask_greater_or_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(_mm256_cmpeq_epi32(_mm256_max_epi32(c0.v0, c1.v0), c0.v0), _mm256_cmpeq_epi32(_mm256_max_epi32(c0.v1, c1.v1), c0.v1), _mm256_cmpeq_epi32(_mm256_max_epi32(c0.v2, c1.v2), c0.v2), _mm256_cmpeq_epi32(_mm256_max_epi32(c0.v3, c1.v3), c0.v3));
}

static uint64_t libsais_mask_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(_mm256_cmpeq_epi32(c0.v0, c1.v0), _mm256_cmpeq_epi32(c0.v1, c1.v1), _mm256_cmpeq_epi32(c0.v2, c1.v2), _mm256_cmpeq_epi32(c0.v3, c1.v3));
}

static uint64_t libsais_mask_negative_32s(libsais_lanes_32s_t c)
{
    return libsais_mask_32s(c.v0, c.v1, c.v2, c.v3);
}

#else

#define LIBSAIS_LANES_32S (4)

typedef __m256i libsais_lanes_32s_t;

static libsais_lanes_32s_t libsais_lanes_broadcast_32s(fast_sint_t c)
{
    return _mm256_set1_epi64x((long long)(c));
}

static libsais_lanes_32s_t libsais_load_reversed_32s(const sa_sint_t * RESTRICT p)
{
    return _mm256_permute4x64_epi64(_mm256_load_si256((const __m256i *)(const void *)p), _MM_SHUFFLE(0, 1, 2, 3));
}

#define libsais_lanes_shift_in_32s(cur, prev)           _mm256_alignr_epi8((cur), _mm256_permute2x128_si256((prev), (cur), 0x21), 8)

static uint64_t libsais_mask_greater_or_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return (uint64_t)(uint32_t)(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(c1, c0))) ^ 15);
}

static uint64_t libsais_mask_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return (uint64_t)(uint32_t)_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(c0, c1)));
}

static uint64_t libsais_mask_negative_32s(libsais_lanes_32s_t c)
{
    return (uint64_t)(uint32_t)_mm256_movemask_pd(_mm256_castsi256_pd(c));
}

#endif

#if defined(LIBSAIS_OPENMP)

#if SAINT_BIT == 32

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    const __m256i zero = _mm256_setzero_si256();

    return (fast_uint_t)(uint32_t)_mm256_movemask_epi8(_mm256_or_si256(_mm256_cmpeq_epi32(_mm256_load_si256((const __m256i *)(const void *)&p[0]), zero), _mm256_cmpeq_epi32(_mm256_load_si256((const __m256i *)(const void *)&p[8]), zero)));
}

#else

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    const __m256i zero = _mm256_setzero_si256();

    __m256i z0 = _mm256_cmpeq_epi64(_mm256_load_si256((const __m256i *)(const void *)&p[0]), zero), z1 = _mm256_cmpeq_epi64(_mm256_load_si256((const __m256i *)(const void *)&p[4]), zero);
    __m256i z2 = _mm256_cmpeq_epi64(_mm256_load_si256((const __m256i *)(const void *)&p[8]), zero), z3 = _mm256_cmpeq_epi64(_mm256_load_si256((const __m256i *)(const void *)&p[12]), zero);

    return (fast_uint_t)(uint32_t)_mm256_movemask_epi8(_mm256_or_si256(_mm256_or_si256(z0, z1), _mm256_or_si256(z2, z3)));
}

#endif

#endif

#elif LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2

#define LIBSAIS_LANES (16)

typedef __m128i libsais_lanes_t;

static libsais_lanes_t libsais_lanes_broadcast(fast_sint_t c)
{
    return _mm_set1_epi8((char)(c));
}

static libsais_lanes_t libsais_load_reversed(const uint8_t * RESTRICT p)
{
    __m128i v = _mm_load_si128((const __m128i *)(const void *)p);

    v = _mm_shuffle_epi32(v, _MM_SHUFFLE(0, 1, 2, 3));
    v = _mm_shufflehi_epi16(_mm_shufflelo_epi16(v, _MM_SHUFFLE(2, 3, 0, 1)), _MM_SHUFFLE(2, 3, 0, 1));

    return _mm_or_si128(_mm_slli_epi16(v, 8), _mm_srli_epi16(v, 8));
}

#define libsais_lanes_shift_in(cur, prev)   _mm_or_si128(_mm_slli_si128((cur), 1), _mm_srli_si128((prev), 15))
#define libsais_lanes_shift_out(cur, next)  _mm_or_si128(_mm_srli_si128((cur), 1), _mm_slli_si128((next), 15))
#define libsais_lanes_type_pair(own, pred)  _mm_sub_epi8(_mm_and_si128((own), _mm_set1_epi8(2)), (pred))

static uint64_t libsais_mask_greater_or_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(_mm_max_epu8(c0, c1), c0));
}

static uint64_t libsais_mask_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(c0, c1));
}

static void libsais_lanes_from_mask(uint64_t bits, libsais_lanes_t * RESTRICT lanes)
{
    const __m128i weights = _mm_set_epi8(-128, 64, 32, 16, 8, 4, 2, 1, -128, 64, 32, 16, 8, 4, 2, 1);

    __m128i v = _mm_set_epi64x(0, (long long)bits); v = _mm_unpacklo_epi8(v, v);
    __m128i lo = _mm_unpacklo_epi16(v, v), hi = _mm_unpackhi_epi16(v, v);

    lanes[0] = _mm_cmpeq_epi8(_mm_and_si128(_mm_shuffle_epi32(lo, _MM_SHUFFLE(1, 1, 0, 0)), weights), weights);
    lanes[1] = _mm_cmpeq_epi8(_mm_and_si128(_mm_shuffle_epi32(lo, _MM_SHUFFLE(3, 3, 2, 2)), weights), weights);
    lanes[2] = _mm_cmpeq_epi8(_mm_and_si128(_mm_shuffle_epi32(hi, _MM_SHUFFLE(1, 1, 0, 0)), weights), weights);
    lanes[3] = _mm_cmpeq_epi8(_mm_and_si128(_mm_shuffle_epi32(hi, _MM_SHUFFLE(3, 3, 2, 2)), weights), weights);
}

static void libsais_store_indices(libsais_lanes_t c, libsais_lanes_t t, libsais_index_t * RESTRICT index)
{
    _mm_storeu_si128((__m128i *)(void *)&index[0], _mm_unpacklo_epi8(c, t));
    _mm_storeu_si128((__m128i *)(void *)&index[8], _mm_unpackhi_epi8(c, t));
}

#define LIBSAIS_LMS_MARGIN (36)

#if SAINT_BIT == 32

typedef __m128i libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return _mm_set1_epi32((int)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    __m128i v = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(const void *)lanes), _mm_setzero_si128());

    _mm_storeu_si128((__m128i *)(void *)&p[0], _mm_sub_epi32(base, _mm_unpacklo_epi16(v, _mm_setzero_si128())));
    _mm_storeu_si128((__m128i *)(void *)&p[4], _mm_sub_epi32(base, _mm_unpackhi_epi16(v, _mm_setzero_si128())));

    return _mm_sub_epi32(base, _mm_set1_epi32(8));
}

#else

typedef __m128i libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return _mm_set1_epi64x((long long)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    __m128i z = _mm_setzero_si128(), v = _mm_unpacklo_epi8(_mm_loadl_epi64((const __m128i *)(const void *)lanes), z);
    __m128i lo = _mm_unpacklo_epi16(v, z), hi = _mm_unpackhi_epi16(v, z);

    _mm_storeu_si128((__m128i *)(void *)&p[0], _mm_sub_epi64(base, _mm_unpacklo_epi32(lo, z)));
    _mm_storeu_si128((__m128i *)(void *)&p[2], _mm_sub_epi64(base, _mm_unpackhi_epi32(lo, z)));
    _mm_storeu_si128((__m128i *)(void *)&p[4], _mm_sub_epi64(base, _mm_unpacklo_epi32(hi, z)));
    _mm_storeu_si128((__m128i *)(void *)&p[6], _mm_sub_epi64(base, _mm_unpackhi_epi32(hi, z)));

    return _mm_sub_epi64(base, _mm_set1_epi64x(8));
}

#endif

#if SAINT_BIT == 32

#define LIBSAIS_LANES_32S (16)

typedef struct libsais_lanes_32s_t
{
    __m128i                             v0, v1, v2, v3;
} libsais_lanes_32s_t;

static __m128i libsais_load_reversed_4x32s(const sa_sint_t * RESTRICT p)
{
    return _mm_shuffle_epi32(_mm_load_si128((const __m128i *)(const void *)p), _MM_SHUFFLE(0, 1, 2, 3));
}

static libsais_lanes_32s_t libsais_lanes_broadcast_32s(fast_sint_t c)
{
    libsais_lanes_32s_t r; r.v0 = r.v1 = r.v2 = r.v3 = _mm_set1_epi32((int)(c));

    return r;
}

static libsais_lanes_32s_t libsais_load_reversed_32s(const sa_sint_t * RESTRICT p)
{
    libsais_lanes_32s_t r; r.v0 = libsais_load_reversed_4x32s(&p[12]); r.v1 = libsais_load_reversed_4x32s(&p[8]); r.v2 = libsais_load_reversed_4x32s(&p[4]); r.v3 = libsais_load_reversed_4x32s(&p[0]);

    return r;
}

static uint64_t libsais_mask_32s(__m128i m0, __m128i m1, __m128i m2, __m128i m3)
{
    return (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_packs_epi16(_mm_packs_epi32(m0, m1), _mm_packs_epi32(m2, m3)));
}

static uint64_t libsais_mask_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(_mm_cmpeq_epi32(c0.v0, c1.v0), _mm_cmpeq_epi32(c0.v1, c1.v1), _mm_cmpeq_epi32(c0.v2, c1.v2), _mm_cmpeq_epi32(c0.v3, c1.v3));
}

static uint64_t libsais_mask_negative_32s(libsais_lanes_32s_t c)
{
    return libsais_mask_32s(c.v0, c.v1, c.v2, c.v3);
}

#define libsais_lanes_funnel_32s(cur, prev)             _mm_or_si128(_mm_slli_si128((cur), 4), _mm_srli_si128((prev), 12))

static uint64_t libsais_mask_greater_or_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(_mm_cmpgt_epi32(c1.v0, c0.v0), _mm_cmpgt_epi32(c1.v1, c0.v1), _mm_cmpgt_epi32(c1.v2, c0.v2), _mm_cmpgt_epi32(c1.v3, c0.v3)) ^ 0xffff;
}

static libsais_lanes_32s_t libsais_lanes_shift_in_32s(libsais_lanes_32s_t cur, libsais_lanes_32s_t prev)
{
    libsais_lanes_32s_t r; r.v0 = libsais_lanes_funnel_32s(cur.v0, prev.v3); r.v1 = libsais_lanes_funnel_32s(cur.v1, cur.v0); r.v2 = libsais_lanes_funnel_32s(cur.v2, cur.v1); r.v3 = libsais_lanes_funnel_32s(cur.v3, cur.v2);

    return r;
}

#else

#define LIBSAIS_LANES_32S (8)

typedef struct libsais_lanes_32s_t
{
    __m128i                             v0, v1, v2, v3;
} libsais_lanes_32s_t;

static __m128i libsais_load_reversed_2x64s(const sa_sint_t * RESTRICT p)
{
    return _mm_shuffle_epi32(_mm_load_si128((const __m128i *)(const void *)p), _MM_SHUFFLE(1, 0, 3, 2));
}

static libsais_lanes_32s_t libsais_lanes_broadcast_32s(fast_sint_t c)
{
    libsais_lanes_32s_t r; r.v0 = r.v1 = r.v2 = r.v3 = _mm_set1_epi64x((long long)(c));

    return r;
}

static libsais_lanes_32s_t libsais_load_reversed_32s(const sa_sint_t * RESTRICT p)
{
    libsais_lanes_32s_t r; r.v0 = libsais_load_reversed_2x64s(&p[6]); r.v1 = libsais_load_reversed_2x64s(&p[4]); r.v2 = libsais_load_reversed_2x64s(&p[2]); r.v3 = libsais_load_reversed_2x64s(&p[0]);

    return r;
}

static __m128i libsais_lanes_high_32s(__m128i a, __m128i b)
{
    return _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(a), _mm_castsi128_ps(b), _MM_SHUFFLE(3, 1, 3, 1)));
}

static __m128i libsais_lanes_low_32s(__m128i a, __m128i b)
{
    return _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(a), _mm_castsi128_ps(b), _MM_SHUFFLE(2, 0, 2, 0)));
}

static uint64_t libsais_mask_32s(__m128i m0, __m128i m1)
{
    return (uint64_t)(uint32_t)_mm_movemask_epi8(_mm_packs_epi16(_mm_packs_epi32(m0, m1), _mm_setzero_si128())) & 0xff;
}

static __m128i libsais_lanes_greater_4x64s(__m128i h0, __m128i l0, __m128i h1, __m128i l1)
{
    const __m128i bias = _mm_set1_epi32((int)0x80000000u);

    return _mm_or_si128(_mm_cmpgt_epi32(h1, h0), _mm_and_si128(_mm_cmpeq_epi32(h1, h0), _mm_cmpgt_epi32(_mm_xor_si128(l1, bias), _mm_xor_si128(l0, bias))));
}

static uint64_t libsais_mask_greater_or_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    __m128i h0a = libsais_lanes_high_32s(c0.v0, c0.v1), h0b = libsais_lanes_high_32s(c0.v2, c0.v3), l0a = libsais_lanes_low_32s(c0.v0, c0.v1), l0b = libsais_lanes_low_32s(c0.v2, c0.v3);
    __m128i h1a = libsais_lanes_high_32s(c1.v0, c1.v1), h1b = libsais_lanes_high_32s(c1.v2, c1.v3), l1a = libsais_lanes_low_32s(c1.v0, c1.v1), l1b = libsais_lanes_low_32s(c1.v2, c1.v3);

    return libsais_mask_32s(libsais_lanes_greater_4x64s(h0a, l0a, h1a, l1a), libsais_lanes_greater_4x64s(h0b, l0b, h1b, l1b)) ^ 0xff;
}

static uint64_t libsais_mask_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    __m128i e0 = _mm_cmpeq_epi32(c0.v0, c1.v0), e1 = _mm_cmpeq_epi32(c0.v1, c1.v1), e2 = _mm_cmpeq_epi32(c0.v2, c1.v2), e3 = _mm_cmpeq_epi32(c0.v3, c1.v3);

    return libsais_mask_32s(_mm_and_si128(libsais_lanes_high_32s(e0, e1), libsais_lanes_low_32s(e0, e1)), _mm_and_si128(libsais_lanes_high_32s(e2, e3), libsais_lanes_low_32s(e2, e3)));
}

static uint64_t libsais_mask_negative_32s(libsais_lanes_32s_t c)
{
    return libsais_mask_32s(libsais_lanes_high_32s(c.v0, c.v1), libsais_lanes_high_32s(c.v2, c.v3));
}

#define libsais_lanes_funnel_32s(cur, prev)             _mm_or_si128(_mm_slli_si128((cur), 8), _mm_srli_si128((prev), 8))

static libsais_lanes_32s_t libsais_lanes_shift_in_32s(libsais_lanes_32s_t cur, libsais_lanes_32s_t prev)
{
    libsais_lanes_32s_t r; r.v0 = libsais_lanes_funnel_32s(cur.v0, prev.v3); r.v1 = libsais_lanes_funnel_32s(cur.v1, cur.v0); r.v2 = libsais_lanes_funnel_32s(cur.v2, cur.v1); r.v3 = libsais_lanes_funnel_32s(cur.v3, cur.v2);

    return r;
}

#endif

#if defined(LIBSAIS_OPENMP)

#if SAINT_BIT == 32

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    const __m128i zero = _mm_setzero_si128();

    __m128i z0 = _mm_cmpeq_epi32(_mm_load_si128((const __m128i *)(const void *)&p[0]), zero), z1 = _mm_cmpeq_epi32(_mm_load_si128((const __m128i *)(const void *)&p[4]), zero);
    __m128i z2 = _mm_cmpeq_epi32(_mm_load_si128((const __m128i *)(const void *)&p[8]), zero), z3 = _mm_cmpeq_epi32(_mm_load_si128((const __m128i *)(const void *)&p[12]), zero);

    return (fast_uint_t)(uint32_t)_mm_movemask_epi8(_mm_or_si128(_mm_or_si128(z0, z1), _mm_or_si128(z2, z3)));
}

#else

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    __m128i z = _mm_setzero_si128(); fast_sint_t j;

    for (j = 0; j < 16; j += 2)
    {
        __m128i v = _mm_load_si128((const __m128i *)(const void *)&p[j]);

        z = _mm_or_si128(z, _mm_cmpeq_epi32(_mm_or_si128(v, _mm_shuffle_epi32(v, _MM_SHUFFLE(2, 3, 0, 1))), _mm_setzero_si128()));
    }

    return (fast_uint_t)(uint32_t)_mm_movemask_epi8(z);
}

#endif

#endif

#elif LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64

#define LIBSAIS_LANES (16)

typedef uint8x16_t libsais_lanes_t;

static libsais_lanes_t libsais_lanes_broadcast(fast_sint_t c)
{
    return vdupq_n_u8((uint8_t)(c));
}

static libsais_lanes_t libsais_load_reversed(const uint8_t * RESTRICT p)
{
    uint8x16_t v = vld1q_u8(p);

    return vrev64q_u8(vextq_u8(v, v, 8));
}

#define libsais_lanes_shift_in(cur, prev)   vextq_u8((prev), (cur), 15)
#define libsais_lanes_shift_out(cur, next)  vextq_u8((cur), (next), 1)
#define libsais_lanes_type_pair(own, pred)  vsubq_u8(vandq_u8((own), vdupq_n_u8(2)), (pred))

static uint64_t libsais_mask(uint8x16_t m)
{
    const uint8x16_t weights = vcombine_u8(vcreate_u8(UINT64_C(0x8040201008040201)), vcreate_u8(UINT64_C(0x8040201008040201)));

    m = vandq_u8(m, weights);

    return (uint64_t)vaddv_u8(vget_low_u8(m)) | ((uint64_t)vaddv_u8(vget_high_u8(m)) << 8);
}

static uint64_t libsais_mask_greater_or_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return libsais_mask(vcgeq_u8(c0, c1));
}

static uint64_t libsais_mask_equal(libsais_lanes_t c0, libsais_lanes_t c1)
{
    return libsais_mask(vceqq_u8(c0, c1));
}

static void libsais_lanes_from_mask(uint64_t bits, libsais_lanes_t * RESTRICT lanes)
{
    const uint8x16_t weights = vcombine_u8(vcreate_u8(UINT64_C(0x8040201008040201)), vcreate_u8(UINT64_C(0x8040201008040201)));

    uint8x16_t v = vreinterpretq_u8_u64(vdupq_n_u64(bits)); fast_sint_t j;

    for (j = 0; j < 4; j += 1)
    {
        lanes[j] = vtstq_u8(vqtbl1q_u8(v, vcombine_u8(vdup_n_u8((uint8_t)(2 * j)), vdup_n_u8((uint8_t)(2 * j + 1)))), weights);
    }
}

static void libsais_store_indices(libsais_lanes_t c, libsais_lanes_t t, libsais_index_t * RESTRICT index)
{
    uint8x16x2_t v; v.val[0] = c; v.val[1] = t;

    vst2q_u8((uint8_t *)(void *)index, v);
}

#define LIBSAIS_LMS_MARGIN (36)

#if SAINT_BIT == 32

typedef int32x4_t libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return vdupq_n_s32((int32_t)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    uint16x8_t v = vmovl_u8(vld1_u8(lanes));

    vst1q_s32(&p[0], vsubq_s32(base, vreinterpretq_s32_u32(vmovl_u16(vget_low_u16(v)))));
    vst1q_s32(&p[4], vsubq_s32(base, vreinterpretq_s32_u32(vmovl_high_u16(v))));

    return vsubq_s32(base, vdupq_n_s32(8));
}

#else

typedef int64x2_t libsais_lms_base_t;

static libsais_lms_base_t libsais_lms_base(fast_sint_t base)
{
    return vdupq_n_s64((int64_t)base);
}

static libsais_lms_base_t libsais_store_lms_lanes(sa_sint_t * RESTRICT p, libsais_lms_base_t base, const uint8_t * RESTRICT lanes)
{
    uint16x8_t v = vmovl_u8(vld1_u8(lanes)); uint32x4_t lo = vmovl_u16(vget_low_u16(v)), hi = vmovl_high_u16(v);

    vst1q_s64(&p[0], vsubq_s64(base, vreinterpretq_s64_u64(vmovl_u32(vget_low_u32(lo)))));
    vst1q_s64(&p[2], vsubq_s64(base, vreinterpretq_s64_u64(vmovl_high_u32(lo))));
    vst1q_s64(&p[4], vsubq_s64(base, vreinterpretq_s64_u64(vmovl_u32(vget_low_u32(hi)))));
    vst1q_s64(&p[6], vsubq_s64(base, vreinterpretq_s64_u64(vmovl_high_u32(hi))));

    return vsubq_s64(base, vdupq_n_s64(8));
}

#endif

#if SAINT_BIT == 32

#define LIBSAIS_LANES_32S (4)

typedef int32x4_t libsais_lanes_32s_t;

static libsais_lanes_32s_t libsais_lanes_broadcast_32s(fast_sint_t c)
{
    return vdupq_n_s32((int32_t)(c));
}

static libsais_lanes_32s_t libsais_load_reversed_32s(const sa_sint_t * RESTRICT p)
{
    int32x4_t v = vld1q_s32(p);

    return vrev64q_s32(vextq_s32(v, v, 2));
}

#define libsais_lanes_shift_in_32s(cur, prev)           vextq_s32((prev), (cur), 3)

static uint64_t libsais_mask_32s(uint32x4_t m)
{
    static const uint32_t lanes[4] = { 1, 2, 4, 8 };

    return (uint64_t)vaddvq_u32(vandq_u32(m, vld1q_u32(lanes)));
}

static uint64_t libsais_mask_greater_or_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(vcgeq_s32(c0, c1));
}

static uint64_t libsais_mask_equal_32s(libsais_lanes_32s_t c0, libsais_lanes_32s_t c1)
{
    return libsais_mask_32s(vceqq_s32(c0, c1));
}

static uint64_t libsais_mask_negative_32s(libsais_lanes_32s_t c)
{
    return libsais_mask_32s(vcltzq_s32(c));
}

#endif

#if defined(LIBSAIS_OPENMP)

#if SAINT_BIT == 32

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    uint32x4_t z0 = vceqzq_s32(vld1q_s32(&p[0])), z1 = vceqzq_s32(vld1q_s32(&p[4])), z2 = vceqzq_s32(vld1q_s32(&p[8])), z3 = vceqzq_s32(vld1q_s32(&p[12]));

    return (fast_uint_t)vmaxvq_u32(vorrq_u32(vorrq_u32(z0, z1), vorrq_u32(z2, z3)));
}

#else

static fast_uint_t libsais_has_zero(const sa_sint_t * RESTRICT p)
{
    uint64x2_t z = vceqzq_s64(vld1q_s64(&p[0])); fast_sint_t j;

    for (j = 2; j < 16; j += 2) { z = vorrq_u64(z, vceqzq_s64(vld1q_s64(&p[j]))); }

    return (fast_uint_t)vmaxvq_u32(vreinterpretq_u32_u64(z));
}

#endif

#endif

#elif defined(LIBSAIS_OPENMP)

#define LIBSAIS_LMS_MARGIN (0)

#endif

#if defined(LIBSAIS_LANES)

#define LIBSAIS_LMS_HISTOGRAMS (4 * 4 * ALPHABET_SIZE)

static void libsais_count_indexed_suffixes_8u(const libsais_index_t * RESTRICT index, sa_sint_t * RESTRICT histograms)
{
    sa_sint_t * RESTRICT histograms1 = &histograms[1 * 4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT histograms2 = &histograms[2 * 4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT histograms3 = &histograms[3 * 4 * ALPHABET_SIZE];

    fast_sint_t j;

    for (j = 0; j < 60; j += 4)
    {
        histograms[index[j + 0]]++; histograms1[index[j + 1]]++; histograms2[index[j + 2]]++; histograms3[index[j + 3]]++;
    }

    histograms[index[60]]++; histograms1[index[61]]++; histograms2[index[62]]++;
}

static void libsais_merge_lms_histograms(sa_sint_t * RESTRICT buckets, const sa_sint_t * RESTRICT histograms)
{
    fast_sint_t j;

    for (j = 0; j < 4 * ALPHABET_SIZE; j += 1) { buckets[BUCKETS_INDEX4((fast_uint_t)j % ALPHABET_SIZE, (fast_uint_t)j / ALPHABET_SIZE)] += histograms[j] + histograms[j + 1 * 4 * ALPHABET_SIZE] + histograms[j + 2 * 4 * ALPHABET_SIZE] + histograms[j + 3 * 4 * ALPHABET_SIZE]; }
}

static const uint8_t libsais_lms_lanes[256][8] =
{
    { 0, 0, 0, 0, 0, 0, 0, 0 }, { 1, 0, 0, 0, 0, 0, 0, 1 }, { 1, 0, 0, 0, 0, 0, 0, 2 }, { 2, 0, 0, 0, 0, 0, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 3 }, { 2, 0, 0, 0, 0, 0, 3, 1 }, { 2, 0, 0, 0, 0, 0, 3, 2 }, { 3, 0, 0, 0, 0, 3, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 4 }, { 2, 0, 0, 0, 0, 0, 4, 1 }, { 2, 0, 0, 0, 0, 0, 4, 2 }, { 3, 0, 0, 0, 0, 4, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 4, 3 }, { 3, 0, 0, 0, 0, 4, 3, 1 }, { 3, 0, 0, 0, 0, 4, 3, 2 }, { 4, 0, 0, 0, 4, 3, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 5 }, { 2, 0, 0, 0, 0, 0, 5, 1 }, { 2, 0, 0, 0, 0, 0, 5, 2 }, { 3, 0, 0, 0, 0, 5, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 5, 3 }, { 3, 0, 0, 0, 0, 5, 3, 1 }, { 3, 0, 0, 0, 0, 5, 3, 2 }, { 4, 0, 0, 0, 5, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 5, 4 }, { 3, 0, 0, 0, 0, 5, 4, 1 }, { 3, 0, 0, 0, 0, 5, 4, 2 }, { 4, 0, 0, 0, 5, 4, 2, 1 },
    { 3, 0, 0, 0, 0, 5, 4, 3 }, { 4, 0, 0, 0, 5, 4, 3, 1 }, { 4, 0, 0, 0, 5, 4, 3, 2 }, { 5, 0, 0, 5, 4, 3, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 6 }, { 2, 0, 0, 0, 0, 0, 6, 1 }, { 2, 0, 0, 0, 0, 0, 6, 2 }, { 3, 0, 0, 0, 0, 6, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 6, 3 }, { 3, 0, 0, 0, 0, 6, 3, 1 }, { 3, 0, 0, 0, 0, 6, 3, 2 }, { 4, 0, 0, 0, 6, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 6, 4 }, { 3, 0, 0, 0, 0, 6, 4, 1 }, { 3, 0, 0, 0, 0, 6, 4, 2 }, { 4, 0, 0, 0, 6, 4, 2, 1 },
    { 3, 0, 0, 0, 0, 6, 4, 3 }, { 4, 0, 0, 0, 6, 4, 3, 1 }, { 4, 0, 0, 0, 6, 4, 3, 2 }, { 5, 0, 0, 6, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 6, 5 }, { 3, 0, 0, 0, 0, 6, 5, 1 }, { 3, 0, 0, 0, 0, 6, 5, 2 }, { 4, 0, 0, 0, 6, 5, 2, 1 },
    { 3, 0, 0, 0, 0, 6, 5, 3 }, { 4, 0, 0, 0, 6, 5, 3, 1 }, { 4, 0, 0, 0, 6, 5, 3, 2 }, { 5, 0, 0, 6, 5, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 6, 5, 4 }, { 4, 0, 0, 0, 6, 5, 4, 1 }, { 4, 0, 0, 0, 6, 5, 4, 2 }, { 5, 0, 0, 6, 5, 4, 2, 1 },
    { 4, 0, 0, 0, 6, 5, 4, 3 }, { 5, 0, 0, 6, 5, 4, 3, 1 }, { 5, 0, 0, 6, 5, 4, 3, 2 }, { 6, 0, 6, 5, 4, 3, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 7 }, { 2, 0, 0, 0, 0, 0, 7, 1 }, { 2, 0, 0, 0, 0, 0, 7, 2 }, { 3, 0, 0, 0, 0, 7, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 7, 3 }, { 3, 0, 0, 0, 0, 7, 3, 1 }, { 3, 0, 0, 0, 0, 7, 3, 2 }, { 4, 0, 0, 0, 7, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 7, 4 }, { 3, 0, 0, 0, 0, 7, 4, 1 }, { 3, 0, 0, 0, 0, 7, 4, 2 }, { 4, 0, 0, 0, 7, 4, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 4, 3 }, { 4, 0, 0, 0, 7, 4, 3, 1 }, { 4, 0, 0, 0, 7, 4, 3, 2 }, { 5, 0, 0, 7, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 7, 5 }, { 3, 0, 0, 0, 0, 7, 5, 1 }, { 3, 0, 0, 0, 0, 7, 5, 2 }, { 4, 0, 0, 0, 7, 5, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 5, 3 }, { 4, 0, 0, 0, 7, 5, 3, 1 }, { 4, 0, 0, 0, 7, 5, 3, 2 }, { 5, 0, 0, 7, 5, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 5, 4 }, { 4, 0, 0, 0, 7, 5, 4, 1 }, { 4, 0, 0, 0, 7, 5, 4, 2 }, { 5, 0, 0, 7, 5, 4, 2, 1 },
    { 4, 0, 0, 0, 7, 5, 4, 3 }, { 5, 0, 0, 7, 5, 4, 3, 1 }, { 5, 0, 0, 7, 5, 4, 3, 2 }, { 6, 0, 7, 5, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 7, 6 }, { 3, 0, 0, 0, 0, 7, 6, 1 }, { 3, 0, 0, 0, 0, 7, 6, 2 }, { 4, 0, 0, 0, 7, 6, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 6, 3 }, { 4, 0, 0, 0, 7, 6, 3, 1 }, { 4, 0, 0, 0, 7, 6, 3, 2 }, { 5, 0, 0, 7, 6, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 6, 4 }, { 4, 0, 0, 0, 7, 6, 4, 1 }, { 4, 0, 0, 0, 7, 6, 4, 2 }, { 5, 0, 0, 7, 6, 4, 2, 1 },
    { 4, 0, 0, 0, 7, 6, 4, 3 }, { 5, 0, 0, 7, 6, 4, 3, 1 }, { 5, 0, 0, 7, 6, 4, 3, 2 }, { 6, 0, 7, 6, 4, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 7, 6, 5 }, { 4, 0, 0, 0, 7, 6, 5, 1 }, { 4, 0, 0, 0, 7, 6, 5, 2 }, { 5, 0, 0, 7, 6, 5, 2, 1 },
    { 4, 0, 0, 0, 7, 6, 5, 3 }, { 5, 0, 0, 7, 6, 5, 3, 1 }, { 5, 0, 0, 7, 6, 5, 3, 2 }, { 6, 0, 7, 6, 5, 3, 2, 1 },
    { 4, 0, 0, 0, 7, 6, 5, 4 }, { 5, 0, 0, 7, 6, 5, 4, 1 }, { 5, 0, 0, 7, 6, 5, 4, 2 }, { 6, 0, 7, 6, 5, 4, 2, 1 },
    { 5, 0, 0, 7, 6, 5, 4, 3 }, { 6, 0, 7, 6, 5, 4, 3, 1 }, { 6, 0, 7, 6, 5, 4, 3, 2 }, { 7, 7, 6, 5, 4, 3, 2, 1 },
    { 1, 0, 0, 0, 0, 0, 0, 8 }, { 2, 0, 0, 0, 0, 0, 8, 1 }, { 2, 0, 0, 0, 0, 0, 8, 2 }, { 3, 0, 0, 0, 0, 8, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 8, 3 }, { 3, 0, 0, 0, 0, 8, 3, 1 }, { 3, 0, 0, 0, 0, 8, 3, 2 }, { 4, 0, 0, 0, 8, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 8, 4 }, { 3, 0, 0, 0, 0, 8, 4, 1 }, { 3, 0, 0, 0, 0, 8, 4, 2 }, { 4, 0, 0, 0, 8, 4, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 4, 3 }, { 4, 0, 0, 0, 8, 4, 3, 1 }, { 4, 0, 0, 0, 8, 4, 3, 2 }, { 5, 0, 0, 8, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 8, 5 }, { 3, 0, 0, 0, 0, 8, 5, 1 }, { 3, 0, 0, 0, 0, 8, 5, 2 }, { 4, 0, 0, 0, 8, 5, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 5, 3 }, { 4, 0, 0, 0, 8, 5, 3, 1 }, { 4, 0, 0, 0, 8, 5, 3, 2 }, { 5, 0, 0, 8, 5, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 5, 4 }, { 4, 0, 0, 0, 8, 5, 4, 1 }, { 4, 0, 0, 0, 8, 5, 4, 2 }, { 5, 0, 0, 8, 5, 4, 2, 1 },
    { 4, 0, 0, 0, 8, 5, 4, 3 }, { 5, 0, 0, 8, 5, 4, 3, 1 }, { 5, 0, 0, 8, 5, 4, 3, 2 }, { 6, 0, 8, 5, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 8, 6 }, { 3, 0, 0, 0, 0, 8, 6, 1 }, { 3, 0, 0, 0, 0, 8, 6, 2 }, { 4, 0, 0, 0, 8, 6, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 6, 3 }, { 4, 0, 0, 0, 8, 6, 3, 1 }, { 4, 0, 0, 0, 8, 6, 3, 2 }, { 5, 0, 0, 8, 6, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 6, 4 }, { 4, 0, 0, 0, 8, 6, 4, 1 }, { 4, 0, 0, 0, 8, 6, 4, 2 }, { 5, 0, 0, 8, 6, 4, 2, 1 },
    { 4, 0, 0, 0, 8, 6, 4, 3 }, { 5, 0, 0, 8, 6, 4, 3, 1 }, { 5, 0, 0, 8, 6, 4, 3, 2 }, { 6, 0, 8, 6, 4, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 6, 5 }, { 4, 0, 0, 0, 8, 6, 5, 1 }, { 4, 0, 0, 0, 8, 6, 5, 2 }, { 5, 0, 0, 8, 6, 5, 2, 1 },
    { 4, 0, 0, 0, 8, 6, 5, 3 }, { 5, 0, 0, 8, 6, 5, 3, 1 }, { 5, 0, 0, 8, 6, 5, 3, 2 }, { 6, 0, 8, 6, 5, 3, 2, 1 },
    { 4, 0, 0, 0, 8, 6, 5, 4 }, { 5, 0, 0, 8, 6, 5, 4, 1 }, { 5, 0, 0, 8, 6, 5, 4, 2 }, { 6, 0, 8, 6, 5, 4, 2, 1 },
    { 5, 0, 0, 8, 6, 5, 4, 3 }, { 6, 0, 8, 6, 5, 4, 3, 1 }, { 6, 0, 8, 6, 5, 4, 3, 2 }, { 7, 8, 6, 5, 4, 3, 2, 1 },
    { 2, 0, 0, 0, 0, 0, 8, 7 }, { 3, 0, 0, 0, 0, 8, 7, 1 }, { 3, 0, 0, 0, 0, 8, 7, 2 }, { 4, 0, 0, 0, 8, 7, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 7, 3 }, { 4, 0, 0, 0, 8, 7, 3, 1 }, { 4, 0, 0, 0, 8, 7, 3, 2 }, { 5, 0, 0, 8, 7, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 7, 4 }, { 4, 0, 0, 0, 8, 7, 4, 1 }, { 4, 0, 0, 0, 8, 7, 4, 2 }, { 5, 0, 0, 8, 7, 4, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 4, 3 }, { 5, 0, 0, 8, 7, 4, 3, 1 }, { 5, 0, 0, 8, 7, 4, 3, 2 }, { 6, 0, 8, 7, 4, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 7, 5 }, { 4, 0, 0, 0, 8, 7, 5, 1 }, { 4, 0, 0, 0, 8, 7, 5, 2 }, { 5, 0, 0, 8, 7, 5, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 5, 3 }, { 5, 0, 0, 8, 7, 5, 3, 1 }, { 5, 0, 0, 8, 7, 5, 3, 2 }, { 6, 0, 8, 7, 5, 3, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 5, 4 }, { 5, 0, 0, 8, 7, 5, 4, 1 }, { 5, 0, 0, 8, 7, 5, 4, 2 }, { 6, 0, 8, 7, 5, 4, 2, 1 },
    { 5, 0, 0, 8, 7, 5, 4, 3 }, { 6, 0, 8, 7, 5, 4, 3, 1 }, { 6, 0, 8, 7, 5, 4, 3, 2 }, { 7, 8, 7, 5, 4, 3, 2, 1 },
    { 3, 0, 0, 0, 0, 8, 7, 6 }, { 4, 0, 0, 0, 8, 7, 6, 1 }, { 4, 0, 0, 0, 8, 7, 6, 2 }, { 5, 0, 0, 8, 7, 6, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 6, 3 }, { 5, 0, 0, 8, 7, 6, 3, 1 }, { 5, 0, 0, 8, 7, 6, 3, 2 }, { 6, 0, 8, 7, 6, 3, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 6, 4 }, { 5, 0, 0, 8, 7, 6, 4, 1 }, { 5, 0, 0, 8, 7, 6, 4, 2 }, { 6, 0, 8, 7, 6, 4, 2, 1 },
    { 5, 0, 0, 8, 7, 6, 4, 3 }, { 6, 0, 8, 7, 6, 4, 3, 1 }, { 6, 0, 8, 7, 6, 4, 3, 2 }, { 7, 8, 7, 6, 4, 3, 2, 1 },
    { 4, 0, 0, 0, 8, 7, 6, 5 }, { 5, 0, 0, 8, 7, 6, 5, 1 }, { 5, 0, 0, 8, 7, 6, 5, 2 }, { 6, 0, 8, 7, 6, 5, 2, 1 },
    { 5, 0, 0, 8, 7, 6, 5, 3 }, { 6, 0, 8, 7, 6, 5, 3, 1 }, { 6, 0, 8, 7, 6, 5, 3, 2 }, { 7, 8, 7, 6, 5, 3, 2, 1 },
    { 5, 0, 0, 8, 7, 6, 5, 4 }, { 6, 0, 8, 7, 6, 5, 4, 1 }, { 6, 0, 8, 7, 6, 5, 4, 2 }, { 7, 8, 7, 6, 5, 4, 2, 1 },
    { 6, 0, 8, 7, 6, 5, 4, 3 }, { 7, 8, 7, 6, 5, 4, 3, 1 }, { 7, 8, 7, 6, 5, 4, 3, 2 }, { 8, 7, 6, 5, 4, 3, 2, 1 }
};

static fast_sint_t libsais_store_lms_suffixes(sa_sint_t * RESTRICT SA, fast_sint_t m, fast_sint_t i, uint64_t lms)
{
    libsais_lms_base_t base = libsais_lms_base(i + 1); fast_sint_t j;

    for (j = 0; j < 64; j += 8, lms >>= 8)
    {
        const uint8_t * RESTRICT lanes = libsais_lms_lanes[lms & 0xff]; base = libsais_store_lms_lanes(&SA[m - 7], base, lanes); m -= lanes[0];
    }

    return m;
}

static uint64_t libsais_classify_suffixes_8u(const uint8_t * RESTRICT T, fast_sint_t i, fast_uint_t f0, libsais_lanes_t prev, libsais_lanes_t * RESTRICT c)
{
    uint64_t greater_or_equal = 0, equal = 0, greater; fast_sint_t j;

    for (j = 0; j < 64 / LIBSAIS_LANES; j += 1)
    {
        libsais_lanes_t c0 = libsais_load_reversed(&T[i - LIBSAIS_LANES * j - (LIBSAIS_LANES - 1)]);
        libsais_lanes_t c1 = libsais_lanes_shift_in(c0, prev);

        greater_or_equal |= libsais_mask_greater_or_equal(c0, c1) << (LIBSAIS_LANES * j);
        equal            |= libsais_mask_equal(c0, c1) << (LIBSAIS_LANES * j);

        c[j] = c0; prev = c0;
    }

    greater = greater_or_equal ^ equal;

    return greater_or_equal ^ (equal & (greater_or_equal + greater + f0));
}

static void libsais_index_buckets_8u(const libsais_lanes_t * RESTRICT c, uint64_t types, libsais_index_t * RESTRICT index)
{
    libsais_lanes_t type_lanes[64 / LIBSAIS_LANES]; fast_sint_t j;

    libsais_lanes_from_mask(types, type_lanes);

    for (j = 0; j < 64 / LIBSAIS_LANES - 1; j += 1)
    {
        libsais_store_indices(c[j], libsais_lanes_type_pair(type_lanes[j], libsais_lanes_shift_out(type_lanes[j], type_lanes[j + 1])), &index[LIBSAIS_LANES * j]);
    }

    libsais_store_indices(c[j], libsais_lanes_type_pair(type_lanes[j], libsais_lanes_shift_out(type_lanes[j], type_lanes[j])), &index[LIBSAIS_LANES * j]);
}

#endif

#if defined(LIBSAIS_LANES_32S)

static uint64_t libsais_classify_suffixes_32s(const sa_sint_t * RESTRICT T, fast_sint_t i, fast_uint_t f0, libsais_lanes_32s_t prev, libsais_lanes_32s_t * RESTRICT c)
{
    uint64_t greater_or_equal = 0, equal = 0, greater; fast_sint_t j;

    for (j = 0; j < 64 / LIBSAIS_LANES_32S; j += 1)
    {
        libsais_lanes_32s_t c0 = libsais_load_reversed_32s(&T[i - LIBSAIS_LANES_32S * j - (LIBSAIS_LANES_32S - 1)]);
        libsais_lanes_32s_t c1 = libsais_lanes_shift_in_32s(c0, prev);

        greater_or_equal |= libsais_mask_greater_or_equal_32s(c0, c1) << (LIBSAIS_LANES_32S * j);
        equal            |= libsais_mask_equal_32s(c0, c1) << (LIBSAIS_LANES_32S * j);

        c[j] = c0; prev = c0;
    }

    greater = greater_or_equal ^ equal;

    return greater_or_equal ^ (equal & (greater_or_equal + greater + f0));
}

static uint64_t libsais_negative_suffixes_32s(const libsais_lanes_32s_t * RESTRICT c)
{
    uint64_t negative = 0; fast_sint_t j;

    for (j = 0; j < 64 / LIBSAIS_LANES_32S; j += 1) { negative |= libsais_mask_negative_32s(c[j]) << (LIBSAIS_LANES_32S * j); }

    return negative;
}

#endif

static void libsais_gather_lms_suffixes_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, fast_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 256;

        fast_sint_t i, j = omp_block_start + omp_block_size, c0 = T[omp_block_start + omp_block_size - 1], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

#if defined(LIBSAIS_LANES)
        fast_sint_t bound = ((fast_sint_t)n + 1) / 2 + LIBSAIS_LMS_MARGIN;

        for (i = omp_block_start + omp_block_size - 2; i >= omp_block_start + 63 && (((ptrdiff_t)(const void *)&T[i - 63]) & 63) != 0; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
        }

        {
            libsais_lanes_t prev = libsais_lanes_broadcast(c0);

            for (; i >= omp_block_start + 63 && m >= bound; i -= 64)
            {
                libsais_prefetchr(&T[i - prefetch_distance]);

                libsais_lanes_t c[64 / LIBSAIS_LANES]; uint64_t types = libsais_classify_suffixes_8u(T, i, f0, prev, c);

                f1 = (fast_uint_t)(types & 1); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
                m = libsais_store_lms_suffixes(SA, m, i, ~types & (types >> 1));

                prev = c[64 / LIBSAIS_LANES - 1]; c0 = T[i - 63]; f0 = (fast_uint_t)(types >> 63);
            }
        }

        for (; i >= omp_block_start; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
        }
#else
        for (i = omp_block_start + omp_block_size - 2, j = omp_block_start + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - prefetch_distance]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 0); m -= (f0 & ~f1);
            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i - 1); m -= (f1 & ~f0);
            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 2); m -= (f0 & ~f1);
        }

        for (j -= 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
        }
#endif

        SA[m] = (sa_sint_t)(i + 1);
    }
}

static sa_sint_t libsais_count_middle_lms_suffixes(const sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m)
{
    fast_sint_t c = (fast_sint_t)n - m, count = m; while (count > 0) { fast_sint_t half = count >> 1; if (SA[c + half] < (n >> 1)) { c += half + 1; count -= half + 1; } else { count = half; } }

    return (sa_sint_t)(n - c);
}

#if !defined(LIBSAIS_LANES)

static void libsais_gather_lms_suffixes_8u_split(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, fast_sint_t middle_lms_suffix, fast_sint_t middle_m)
{
    const fast_sint_t prefetch_distance = 256;

    fast_sint_t iA = middle_lms_suffix - 2, mA = (fast_sint_t)n - 1 - middle_m, cA0 = T[middle_lms_suffix - 1], cA1 = -1; fast_uint_t fA0 = 1, fA1 = 0;
    fast_sint_t iB = (fast_sint_t)n - 2, mB = (fast_sint_t)n - 1, cB0 = T[n - 1], cB1 = -1; fast_uint_t fB0 = 1, fB1 = 0;

    for (; iA >= 1 && iB >= middle_lms_suffix + 1; iA -= 2, iB -= 2)
    {
        libsais_prefetchr(&T[iA - prefetch_distance]);
        libsais_prefetchr(&T[iB - prefetch_distance]);

        cA1 = T[iA - 0]; fA1 = (fast_uint_t)(cA1 > (cA0 - (fast_sint_t)(fA0))); SA[mA] = (sa_sint_t)(iA + 1); mA -= (fA1 & ~fA0);
        cB1 = T[iB - 0]; fB1 = (fast_uint_t)(cB1 > (cB0 - (fast_sint_t)(fB0))); SA[mB] = (sa_sint_t)(iB + 1); mB -= (fB1 & ~fB0);
        cA0 = T[iA - 1]; fA0 = (fast_uint_t)(cA0 > (cA1 - (fast_sint_t)(fA1))); SA[mA] = (sa_sint_t)(iA - 0); mA -= (fA0 & ~fA1);
        cB0 = T[iB - 1]; fB0 = (fast_uint_t)(cB0 > (cB1 - (fast_sint_t)(fB1))); SA[mB] = (sa_sint_t)(iB - 0); mB -= (fB0 & ~fB1);
    }

    for (; iA >= 0; iA -= 1)
    {
        cA1 = cA0; cA0 = T[iA]; fA1 = fA0; fA0 = (fast_uint_t)(cA0 > (cA1 - (fast_sint_t)(fA1))); SA[mA] = (sa_sint_t)(iA + 1); mA -= (fA0 & ~fA1);
    }

    SA[mA] = (sa_sint_t)(iA + 1);

    for (; iB >= middle_lms_suffix; iB -= 1)
    {
        cB1 = cB0; cB0 = T[iB]; fB1 = fB0; fB0 = (fast_uint_t)(cB0 > (cB1 - (fast_sint_t)(fB1))); SA[mB] = (sa_sint_t)(iB + 1); mB -= (fB0 & ~fB1);
    }

    SA[mB] = (sa_sint_t)(iB + 1);
}

#endif

static void libsais_gather_lms_suffixes_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t middle_lms_suffix, sa_sint_t middle_m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
#if defined(LIBSAIS_LANES)
            (void)(middle_lms_suffix); (void)(middle_m);

            libsais_gather_lms_suffixes_8u(T, SA, n, (fast_sint_t)n - 1, omp_block_start, omp_block_size);
#else
            if (middle_lms_suffix > 0)
            {
                libsais_gather_lms_suffixes_8u_split(T, SA, n, middle_lms_suffix, middle_m);
            }
            else
            {
                libsais_gather_lms_suffixes_8u(T, SA, n, (fast_sint_t)n - 1, omp_block_start, omp_block_size);
            }
#endif
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < threads)
        {
            #pragma omp master
            {
                libsais_gather_lms_suffixes_8u(T, SA, n, (fast_sint_t)n - 1, 0, n);
            }
        }
        else
        {
            fast_sint_t t, m = 0; for (t = omp_num_threads - 1; t > omp_thread_num; --t) { m += thread_state[t].state.m; }

            libsais_gather_lms_suffixes_8u(T, SA, n, (fast_sint_t)n - 1 - m, omp_block_start, omp_block_size);

            #pragma omp barrier

            if (omp_thread_num < omp_num_threads - 1 && thread_state[omp_thread_num].state.m > 0)
            {
                fast_sint_t i = omp_block_start + omp_block_size - 1, j = i + 1, k = (fast_sint_t)n - 1 - m - (thread_state[omp_thread_num].state.m < LIBSAIS_LMS_MARGIN + 1 ? thread_state[omp_thread_num].state.m : LIBSAIS_LMS_MARGIN + 1);
                fast_sint_t c0 = T[i], c1 = -1; fast_uint_t f0, f1;

                while (j < n && (c1 = T[j]) == c0) { ++j; }

                for (f0 = c0 >= c1, m = (fast_sint_t)n - 1 - m; m > k; --i)
                {
                    c1 = c0; c0 = (i > 0) ? T[i - 1] : -1; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)i; m -= (f0 & ~f1);
                }
            }
        }
#endif
    }
}

static sa_sint_t libsais_gather_lms_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t             i   = n - 2;
    sa_sint_t             m   = n - 1;
    fast_uint_t           f0  = 1;
    fast_uint_t           f1  = 0;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

#if defined(LIBSAIS_LANES_32S)
    fast_sint_t bound = ((fast_sint_t)n + 1) / 2 + LIBSAIS_LMS_MARGIN;

    for (; i >= 63 && (((ptrdiff_t)(const void *)&T[i - 63]) & 63) != 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1);
    }

    {
        libsais_lanes_32s_t prev = libsais_lanes_broadcast_32s(c0);

        for (; i >= 63 && m >= bound; i -= 64)
        {
            fast_sint_t j; for (j = 0; j < 64; j += 64 / (fast_sint_t)sizeof(sa_sint_t)) { libsais_prefetchr(&T[i - prefetch_distance - j]); }

            libsais_lanes_32s_t c[64 / LIBSAIS_LANES_32S]; uint64_t types = libsais_classify_suffixes_32s(T, i, f0, prev, c);

            f1 = (fast_uint_t)(types & 1); SA[m] = i + 1; m -= (sa_sint_t)(f1 & ~f0);
            m = (sa_sint_t)libsais_store_lms_suffixes(SA, m, i, ~types & (types >> 1));

            prev = c[64 / LIBSAIS_LANES_32S - 1]; c0 = T[i - 63]; f0 = (fast_uint_t)(types >> 63);
        }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1);
    }
#else
    for (; i >= 3; i -= 4)
    {
        libsais_prefetchr(&T[i - prefetch_distance]);

        c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = i + 1; m -= (sa_sint_t)(f1 & ~f0);
        c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i - 0; m -= (sa_sint_t)(f0 & ~f1);
        c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = i - 1; m -= (sa_sint_t)(f1 & ~f0);
        c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i - 2; m -= (sa_sint_t)(f0 & ~f1);
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1);
    }
#endif

    return n - 1 - m;
}

static sa_sint_t libsais_gather_compacted_lms_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t             i   = n - 2;
    sa_sint_t             m   = n - 1;
    fast_uint_t           f0  = 1;
    fast_uint_t           f1  = 0;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

#if defined(LIBSAIS_LANES_32S)
    fast_sint_t bound = ((fast_sint_t)n + 1) / 2 + LIBSAIS_LMS_MARGIN;

    for (; i >= 63 && (((ptrdiff_t)(const void *)&T[i - 63]) & 63) != 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1 & (c1 >= 0));
    }

    {
        libsais_lanes_32s_t prev = libsais_lanes_broadcast_32s(c0);

        for (; i >= 63 && m >= bound; i -= 64)
        {
            fast_sint_t j; for (j = 0; j < 64; j += 64 / (fast_sint_t)sizeof(sa_sint_t)) { libsais_prefetchr(&T[i - prefetch_distance - j]); }

            libsais_lanes_32s_t c[64 / LIBSAIS_LANES_32S]; uint64_t types = libsais_classify_suffixes_32s(T, i, f0, prev, c);

            f1 = (fast_uint_t)(types & 1); SA[m] = i + 1; m -= (sa_sint_t)(f1 & ~f0 & (c0 >= 0));
            m = (sa_sint_t)libsais_store_lms_suffixes(SA, m, i, ~types & (types >> 1) & ~libsais_negative_suffixes_32s(c));

            prev = c[64 / LIBSAIS_LANES_32S - 1]; c0 = T[i - 63]; f0 = (fast_uint_t)(types >> 63);
        }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1 & (c1 >= 0));
    }
#else
    for (; i >= 3; i -= 4)
    {
        libsais_prefetchr(&T[i - prefetch_distance]);

        c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = i + 1; m -= (sa_sint_t)(f1 & ~f0 & (c0 >= 0));
        c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i - 0; m -= (sa_sint_t)(f0 & ~f1 & (c1 >= 0));
        c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = i - 1; m -= (sa_sint_t)(f1 & ~f0 & (c0 >= 0));
        c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i - 2; m -= (sa_sint_t)(f0 & ~f1 & (c1 >= 0));
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = i + 1; m -= (sa_sint_t)(f0 & ~f1 & (c1 >= 0));
    }
#endif

    return n - 1 - m;
}

#if defined(LIBSAIS_OPENMP)

static void libsais_count_lms_suffixes_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 64;

    memset(buckets, 0, 4 * (size_t)k * sizeof(sa_sint_t));

    sa_sint_t             i   = n - 2;
    fast_uint_t           f0  = 1;
    fast_uint_t           f1  = 0;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetchr(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 0], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 1], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 2], 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 3], 0)]);

        c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

        c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;

        c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

        c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
    }

    buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0)]++;
}

#endif

static void libsais_count_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 64;

        fast_sint_t m = omp_block_start + omp_block_size - 1;
        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;

            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        c1 = (i >= 0) ? T[i] : -1; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
        buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;
    }
}

#if defined(LIBSAIS_OPENMP)

static void libsais_count_compacted_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 64;

    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    sa_sint_t             i   = n - 2;
    fast_uint_t           f0  = 1;
    fast_uint_t           f1  = 0;
    fast_sint_t           c0  = T[n - 1];
    fast_sint_t           c1  = 0;

    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetchr(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0] & SAINT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1] & SAINT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2] & SAINT_MAX, 0)]);
        libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3] & SAINT_MAX, 0)]);

        c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
        c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

        c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;

        c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0)));
        c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

        c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
    }

    c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, 0)]++;
}

#endif

static sa_sint_t libsais_count_and_gather_lms_suffixes_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, (size_t)4 * ALPHABET_SIZE * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 256;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

#if defined(LIBSAIS_LANES)
        sa_sint_t histograms[LIBSAIS_LMS_HISTOGRAMS]; memset(histograms, 0, sizeof(histograms));

        fast_sint_t bound = omp_block_start + (omp_block_size + 1) / 2 + LIBSAIS_LMS_MARGIN;

        for (i = m - 1; i >= omp_block_start + 63 && (((ptrdiff_t)(const void *)&T[i - 63]) & 63) != 0; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }

        {
            libsais_lanes_t prev = libsais_lanes_broadcast(c0); libsais_index_t index0[64], index1[64]; libsais_index_t * RESTRICT pending = index0; libsais_index_t * RESTRICT current = index1; fast_sint_t blocks = 0;

            for (; i >= omp_block_start + 63 && m >= bound; i -= 64, blocks += 1)
            {
                libsais_prefetchr(&T[i - prefetch_distance]);

                libsais_lanes_t c[64 / LIBSAIS_LANES]; uint64_t types = libsais_classify_suffixes_8u(T, i, f0, prev, c); libsais_index_t * RESTRICT counted = pending;

                f1 = (fast_uint_t)(types & 1); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
                buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

                libsais_index_buckets_8u(c, types, current);

                if (blocks > 0) { libsais_count_indexed_suffixes_8u(pending, histograms); }

                m = libsais_store_lms_suffixes(SA, m, i, ~types & (types >> 1));

                prev = c[64 / LIBSAIS_LANES - 1]; c0 = T[i - 63]; f0 = (fast_uint_t)(types >> 63); pending = current; current = counted;
            }

            if (blocks > 0) { libsais_count_indexed_suffixes_8u(pending, histograms); }
        }

        for (; i >= omp_block_start; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }

        libsais_merge_lms_histograms(buckets, histograms);
#else
        for (i = m - 1, j = omp_block_start + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - prefetch_distance]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 0); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;

            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i - 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 2); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }
#endif

        c1 = (i >= 0) ? T[i] : -1; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t * threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(*threads) if(*threads > 1 && n >= 65536 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_lms_suffixes_8u(T, SA, n, buckets, omp_block_start, omp_block_size); *threads = 1;
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.position = omp_block_start + omp_block_size;
                thread_state[omp_thread_num].state.m = libsais_count_and_gather_lms_suffixes_8u(T, SA, n, thread_state[omp_thread_num].state.buckets, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                memset(buckets, 0, (size_t)4 * ALPHABET_SIZE * sizeof(sa_sint_t));

                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    m += (sa_sint_t)thread_state[t].state.m;

                    if (t != omp_num_threads - 1 && thread_state[t].state.m > 0)
                    {
                        memcpy(&SA[n - m], &SA[thread_state[t].state.position - thread_state[t].state.m], (size_t)thread_state[t].state.m * sizeof(sa_sint_t));
                    }

                    {
                        sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                        fast_sint_t s; for (s = 0; s < 4 * ALPHABET_SIZE; s += 1) { sa_sint_t A = buckets[s], B = temp_bucket[s]; buckets[s] = A + B; temp_bucket[s] = A; }
                    }
                }

                *threads = (sa_sint_t)omp_num_threads;
            }
        }
#endif
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 4 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 64;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 0], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 1], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 2], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX4(T[i - prefetch_distance - 3], 0)]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 0); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;

            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i - 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;

            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 2); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]++;
        }

        c1 = (i >= 0) ? T[i] : -1; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0 + f1)]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 64;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2], 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3], 0)]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 0); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;

            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i - 1); m -= (f1 & ~f0);
            buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 2); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1);
            buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        c1 = (i >= 0) ? T[i] : -1; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0);
        buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

static sa_sint_t libsais_count_and_gather_compacted_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t m = omp_block_start + omp_block_size - 1;

    if (omp_block_size > 0)
    {
        const fast_sint_t prefetch_distance = 64;

        fast_sint_t i, j = m + 1, c0 = T[m], c1 = -1;

        while (j < n && (c1 = T[j]) == c0) { ++j; }

        fast_uint_t f0 = c0 >= c1, f1 = 0;

        for (i = m - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
        {
            libsais_prefetchr(&T[i - 2 * prefetch_distance]);

            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 0] & SAINT_MAX, 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 1] & SAINT_MAX, 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 2] & SAINT_MAX, 0)]);
            libsais_prefetchw(&buckets[BUCKETS_INDEX2(T[i - prefetch_distance - 3] & SAINT_MAX, 0)]);

            c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0 & (c0 >= 0));
            c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 0); m -= (f0 & ~f1 & (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;

            c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i - 1); m -= (f1 & ~f0 & (c0 >= 0));
            c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;

            c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i - 2); m -= (f0 & ~f1 & (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        for (j -= prefetch_distance + 3; i >= j; i -= 1)
        {
            c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); SA[m] = (sa_sint_t)(i + 1); m -= (f0 & ~f1 & (c1 >= 0));
            c1 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c1, (f0 & ~f1))]++;
        }

        c1 = (i >= 0) ? T[i] : -1; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); SA[m] = (sa_sint_t)(i + 1); m -= (f1 & ~f0 & (c0 >= 0));
        c0 &= SAINT_MAX; buckets[BUCKETS_INDEX2((fast_uint_t)c0, (f1 & ~f0))]++;
    }

    return (sa_sint_t)(omp_block_start + omp_block_size - 1 - m);
}

#if defined(LIBSAIS_OPENMP)

static fast_sint_t libsais_get_bucket_stride(fast_sint_t free_space, fast_sint_t bucket_size, fast_sint_t num_buckets)
{
    fast_sint_t bucket_size_1024 = (bucket_size + 1023) & (-1024); if (free_space / (num_buckets - 1) >= bucket_size_1024) { return bucket_size_1024; }
    fast_sint_t bucket_size_16 = (bucket_size + 15) & (-16); if (free_space / (num_buckets - 1) >= bucket_size_16) { return bucket_size_16; }

    return bucket_size;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_4k_fs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(local_buckets); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_lms_suffixes_32s_4k(T, SA, n, k, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t bucket_size       = 4 * (fast_sint_t)k;
            fast_sint_t free_space        = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
            fast_sint_t bucket_stride     = libsais_get_bucket_stride(free_space, bucket_size, omp_num_threads);

            {
                thread_state[omp_thread_num].state.position = omp_block_start + omp_block_size;
                thread_state[omp_thread_num].state.count = libsais_count_and_gather_lms_suffixes_32s_4k(T, SA, n, k, buckets - (omp_thread_num * bucket_stride), omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            if (omp_thread_num == omp_num_threads - 1)
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    m += (sa_sint_t)thread_state[t].state.count;

                    if (t != omp_num_threads - 1 && thread_state[t].state.count > 0)
                    {
                        memcpy(&SA[n - m], &SA[thread_state[t].state.position - thread_state[t].state.count], (size_t)thread_state[t].state.count * sizeof(sa_sint_t));
                    }
                }
            }
            else
            {
                omp_num_threads     = omp_num_threads - 1;
                omp_block_stride    = (bucket_size / omp_num_threads) & (-16);
                omp_block_start     = omp_thread_num * omp_block_stride;
                omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : bucket_size - omp_block_start;

                libsais_accumulate_counts_s32(buckets + omp_block_start, omp_block_size, bucket_stride, omp_num_threads + 1);
            }
        }
#endif
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_2k_fs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(local_buckets); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t bucket_size       = 2 * (fast_sint_t)k;
            fast_sint_t free_space        = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
            fast_sint_t bucket_stride     = libsais_get_bucket_stride(free_space, bucket_size, omp_num_threads);

            {
                thread_state[omp_thread_num].state.position = omp_block_start + omp_block_size;
                thread_state[omp_thread_num].state.count = libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets - (omp_thread_num * bucket_stride), omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            if (omp_thread_num == omp_num_threads - 1)
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    m += (sa_sint_t)thread_state[t].state.count;

                    if (t != omp_num_threads - 1 && thread_state[t].state.count > 0)
                    {
                        memcpy(&SA[n - m], &SA[thread_state[t].state.position - thread_state[t].state.count], (size_t)thread_state[t].state.count * sizeof(sa_sint_t));
                    }
                }
            }
            else
            {
                omp_num_threads     = omp_num_threads - 1;
                omp_block_stride    = (bucket_size / omp_num_threads) & (-16);
                omp_block_start     = omp_thread_num * omp_block_stride;
                omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : bucket_size - omp_block_start;

                libsais_accumulate_counts_s32(buckets + omp_block_start, omp_block_size, bucket_stride, omp_num_threads + 1);
            }
        }
#endif
    }

    return m;
}

static void libsais_count_and_gather_compacted_lms_suffixes_32s_2k_fs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(local_buckets); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_count_and_gather_compacted_lms_suffixes_32s_2k(T, SA, n, k, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t bucket_size       = 2 * (fast_sint_t)k;
            fast_sint_t free_space        = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[(fast_sint_t)n + (fast_sint_t)n]);
            fast_sint_t bucket_stride     = libsais_get_bucket_stride(free_space, bucket_size, omp_num_threads);

            {
                thread_state[omp_thread_num].state.position = omp_block_start + omp_block_size;
                thread_state[omp_thread_num].state.count = libsais_count_and_gather_compacted_lms_suffixes_32s_2k(T, SA + n, n, k, buckets - (omp_thread_num * bucket_stride), omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, m = 0; for (t = omp_num_threads - 1; t >= omp_thread_num; --t) { m += (sa_sint_t)thread_state[t].state.count; }

                if (thread_state[omp_thread_num].state.count > 0)
                {
                    memcpy(&SA[n - m], &SA[n + thread_state[omp_thread_num].state.position - thread_state[omp_thread_num].state.count], (size_t)thread_state[omp_thread_num].state.count * sizeof(sa_sint_t));
                }
            }

            {
                omp_block_stride    = (bucket_size / omp_num_threads) & (-16);
                omp_block_start     = omp_thread_num * omp_block_stride;
                omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : bucket_size - omp_block_start;

                libsais_accumulate_counts_s32(buckets + omp_block_start, omp_block_size, bucket_stride, omp_num_threads);
            }
        }
#endif
    }
}

static void libsais_count_lms_suffixes_32s_2k_fs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(SA); UNUSED(local_buckets); UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_count_lms_suffixes_32s_2k(T, n, k, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t bucket_size       = 2 * (fast_sint_t)k;
            fast_sint_t free_space        = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
            fast_sint_t bucket_stride     = libsais_get_bucket_stride(free_space, bucket_size, omp_num_threads);

            libsais_count_lms_suffixes_32s_2k(T, n, k, buckets - (omp_thread_num * bucket_stride), omp_block_start, omp_block_size);

            #pragma omp barrier

            {
                omp_block_stride    = (bucket_size / omp_num_threads) & (-16);
                omp_block_start     = omp_thread_num * omp_block_stride;
                omp_block_size      = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : bucket_size - omp_block_start;

                libsais_accumulate_counts_s32(buckets + omp_block_start, omp_block_size, bucket_stride, omp_num_threads);
            }
        }
#endif
    }
}

#endif

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_4k_nofs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(2) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_num_threads   = 1;
#endif
        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_lms_suffixes_32s_4k(T, SA, n, k, buckets, 0, n);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_thread_num == 0)
        {
            libsais_count_lms_suffixes_32s_4k(T, n, k, buckets);
        }
        else
        {
            m = libsais_gather_lms_suffixes_32s(T, SA, n);
        }
#endif
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_2k_nofs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(2) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_num_threads   = 1;
#endif
        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets, 0, n);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_thread_num == 0)
        {
            libsais_count_lms_suffixes_32s_2k(T, n, k, buckets, 0, n);
        }
        else
        {
            m = libsais_gather_lms_suffixes_32s(T, SA, n);
        }
#endif
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads)
{
    sa_sint_t m = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(2) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_num_threads   = 1;
#endif
        if (omp_num_threads == 1)
        {
            m = libsais_count_and_gather_compacted_lms_suffixes_32s_2k(T, SA, n, k, buckets, 0, n);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_thread_num == 0)
        {
            libsais_count_compacted_lms_suffixes_32s_2k(T, n, k, buckets);
        }
        else
        {
            m = libsais_gather_compacted_lms_suffixes_32s(T, SA, n);
        }
#endif
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_4k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t m;

#if defined(LIBSAIS_OPENMP)
    fast_sint_t free_space  = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
    sa_sint_t max_threads   = (sa_sint_t)(free_space / ((4 * (fast_sint_t)k + 15) & (-16))); if (max_threads > threads) { max_threads = threads; }

    if (max_threads > 1 && n >= 65536 && n / k >= 2)
    {
        if (max_threads > n / 16 / k) { max_threads = n / 16 / k; }
        m = libsais_count_and_gather_lms_suffixes_32s_4k_fs_omp(T, SA, n, k, buckets, local_buckets, max_threads > 2 ? max_threads : 2, thread_state);
    }
    else
#else
    UNUSED(local_buckets); UNUSED(thread_state);
#endif
    {
        m = libsais_count_and_gather_lms_suffixes_32s_4k_nofs_omp(T, SA, n, k, buckets, threads);
    }

    return m;
}

static sa_sint_t libsais_count_and_gather_lms_suffixes_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t m;

#if defined(LIBSAIS_OPENMP)
    fast_sint_t free_space  = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
    sa_sint_t max_threads   = (sa_sint_t)(free_space / ((2 * (fast_sint_t)k + 15) & (-16))); if (max_threads > threads) { max_threads = threads; }

    if (max_threads > 1 && n >= 65536 && n / k >= 2)
    {
        if (max_threads > n / 8 / k) { max_threads = n / 8 / k; }
        m = libsais_count_and_gather_lms_suffixes_32s_2k_fs_omp(T, SA, n, k, buckets, local_buckets, max_threads > 2 ? max_threads : 2, thread_state);
    }
    else
#else
    UNUSED(local_buckets); UNUSED(thread_state);
#endif
    {
        m = libsais_count_and_gather_lms_suffixes_32s_2k_nofs_omp(T, SA, n, k, buckets, threads);
    }

    return m;
}

static void libsais_count_and_gather_compacted_lms_suffixes_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    fast_sint_t free_space  = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[(fast_sint_t)n + (fast_sint_t)n]);
    sa_sint_t max_threads   = (sa_sint_t)(free_space / ((2 * (fast_sint_t)k + 15) & (-16))); if (max_threads > threads) { max_threads = threads; }

    if (!local_buckets && max_threads > 1 && n >= 65536 && n / k >= 2)
    {
        if (max_threads > n / 8 / k) { max_threads = n / 8 / k; }
        libsais_count_and_gather_compacted_lms_suffixes_32s_2k_fs_omp(T, SA, n, k, buckets, local_buckets, max_threads > 2 ? max_threads : 2, thread_state);
    }
    else
#else
    UNUSED(local_buckets); UNUSED(thread_state);
#endif
    {
        libsais_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(T, SA, n, k, buckets, threads);
    }
}

static void libsais_count_lms_suffixes_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    fast_sint_t free_space  = local_buckets ? LIBSAIS_LOCAL_BUFFER_SIZE : (buckets - &SA[n]);
    sa_sint_t max_threads   = (sa_sint_t)(free_space / ((2 * (fast_sint_t)k + 15) & (-16))); if (max_threads > threads) { max_threads = threads; }

    if (max_threads > 1 && n >= 65536 && n / k >= 2)
    {
        if (max_threads > n / 8 / k) { max_threads = n / 8 / k; }
        libsais_count_lms_suffixes_32s_2k_fs_omp(T, SA, n, k, buckets, local_buckets, max_threads > 2 ? max_threads : 2);
    }
    else
#else
    UNUSED(SA); UNUSED(local_buckets); UNUSED(threads);
#endif
    {
        libsais_count_lms_suffixes_32s_2k(T, n, k, buckets, 0, n);
    }
}

static void libsais_count_suffixes_32s(const sa_sint_t * RESTRICT T, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 64;

    memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

    fast_sint_t i, j;
    for (i = 0, j = (fast_sint_t)n - 7; i < j; i += 8)
    {
        libsais_prefetchr(&T[i + prefetch_distance]);

        buckets[T[i + 0]]++;
        buckets[T[i + 1]]++;
        buckets[T[i + 2]]++;
        buckets[T[i + 3]]++;
        buckets[T[i + 4]]++;
        buckets[T[i + 5]]++;
        buckets[T[i + 6]]++;
        buckets[T[i + 7]]++;
    }

    for (j += 7; i < j; i += 1)
    {
        buckets[T[i]]++;
    }
}

static sa_sint_t libsais_initialize_buckets_start_and_end_8u(sa_sint_t * RESTRICT buckets, sa_sint_t * RESTRICT freq)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[6 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT bucket_end   = &buckets[7 * ALPHABET_SIZE];

    fast_sint_t k = -1;

    if (freq != NULL)
    {
        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        {
            sa_sint_t total = buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];

            bucket_start[j] = sum; sum += total; bucket_end[j] = sum; k = total > 0 ? j : k; freq[j] = total;
        }
    }
    else
    {
        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        {
            sa_sint_t total = buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];

            bucket_start[j] = sum; sum += total; bucket_end[j] = sum; k = total > 0 ? j : k;
        }
    }

    return (sa_sint_t)(k + 1);
}

static void libsais_initialize_buckets_start_and_end_32s_6k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[4 * (fast_sint_t)k];
    sa_sint_t * RESTRICT bucket_end   = &buckets[5 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
    {
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 2)] + buckets[i + BUCKETS_INDEX4(0, 3)];
        bucket_end[j] = sum;
    }
}

static void libsais_initialize_buckets_start_and_end_32s_4k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT bucket_end   = &buckets[3 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum = 0;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    { 
        bucket_start[j] = sum;
        sum += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        bucket_end[j] = sum;
    }
}

static void libsais_initialize_buckets_end_32s_2k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i; sa_sint_t sum0 = 0;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    { 
        sum0 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)]; buckets[i + BUCKETS_INDEX2(0, 0)] = sum0;
    }
}

static void libsais_initialize_buckets_start_and_end_32s_2k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i, j;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    {
        buckets[j] = buckets[i];
    }

    buckets[k] = 0; memcpy(&buckets[k + 1], buckets, ((size_t)k - 1) * sizeof(sa_sint_t));
}

static void libsais_initialize_buckets_start_32s_1k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i; sa_sint_t sum = 0;
    for (i = 0; i <= (fast_sint_t)k - 1; i += 1) { sa_sint_t tmp = buckets[i]; buckets[i] = sum; sum += tmp; }
}

static void libsais_initialize_buckets_end_32s_1k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    fast_sint_t i; sa_sint_t sum = 0;
    for (i = 0; i <= (fast_sint_t)k - 1; i += 1) { sum += buckets[i]; buckets[i] = sum; }
}

static sa_sint_t libsais_initialize_buckets_for_lms_suffixes_radix_sort_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    {
        fast_uint_t     f0 = 0;
        fast_uint_t     f1 = 0;
        fast_sint_t     c0 = T[first_lms_suffix];
        fast_sint_t     c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]--;
        }

        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0)]--;
    }

    {
        sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
        { 
            temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum; sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static void libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 0)]++;
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 1)]--;

    fast_sint_t i; sa_sint_t sum0 = 0, sum1 = 0;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    { 
        sum0 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        sum1 += buckets[i + BUCKETS_INDEX2(0, 1)];
        
        buckets[i + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static sa_sint_t libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    {
        fast_uint_t     f0 = 0;
        fast_uint_t     f1 = 0;
        fast_sint_t     c0 = T[first_lms_suffix];
        fast_sint_t     c1 = 0;

        for (; --first_lms_suffix >= 0; )
        {
            c1 = c0; c0 = T[first_lms_suffix]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
            buckets[BUCKETS_INDEX4((fast_uint_t)c1, f1 + f1 + f0)]--;
        }

        buckets[BUCKETS_INDEX4((fast_uint_t)c0, f0 + f0)]--;
    }

    {
        sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

        fast_sint_t i, j; sa_sint_t sum = 0;
        for (i = BUCKETS_INDEX4(0, 0), j = 0; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += 1)
        { 
            sum += buckets[i + BUCKETS_INDEX4(0, 1)] + buckets[i + BUCKETS_INDEX4(0, 3)]; temp_bucket[j] = sum;
        }

        return sum;
    }
}

static void libsais_initialize_buckets_for_radix_and_partial_sorting_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix)
{
    sa_sint_t * RESTRICT bucket_start = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT bucket_end   = &buckets[3 * (fast_sint_t)k];

    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 0)]++;
    buckets[BUCKETS_INDEX2(T[first_lms_suffix], 1)]--;

    fast_sint_t i, j; sa_sint_t sum0 = 0, sum1 = 0;
    for (i = BUCKETS_INDEX2(0, 0), j = 0; i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0), j += 1)
    { 
        bucket_start[j] = sum1;

        sum0 += buckets[i + BUCKETS_INDEX2(0, 1)];
        sum1 += buckets[i + BUCKETS_INDEX2(0, 0)] + buckets[i + BUCKETS_INDEX2(0, 1)];
        buckets[i + BUCKETS_INDEX2(0, 1)] = sum0;

        bucket_end[j] = sum1;
    }
}

static void libsais_radix_sort_lms_suffixes_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i - prefetch_distance - 0]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 2]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 3]]);

        sa_sint_t p0 = SA[i - 0]; SA[--induction_bucket[BUCKETS_INDEX2(T[p0], 0)]] = p0;
        sa_sint_t p1 = SA[i - 1]; SA[--induction_bucket[BUCKETS_INDEX2(T[p1], 0)]] = p1;
        sa_sint_t p2 = SA[i - 2]; SA[--induction_bucket[BUCKETS_INDEX2(T[p2], 0)]] = p2;
        sa_sint_t p3 = SA[i - 3]; SA[--induction_bucket[BUCKETS_INDEX2(T[p3], 0)]] = p3;
    }

    for (j -= prefetch_distance + 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[--induction_bucket[BUCKETS_INDEX2(T[p], 0)]] = p;
    }
}

static void libsais_radix_sort_lms_suffixes_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t flags, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (flags & LIBSAIS_FLAGS_GSA) { buckets[4 * ALPHABET_SIZE]--; }

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536 && m >= 65536 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_num_threads   = 1;
#endif
        if (omp_num_threads == 1)
        {
            libsais_radix_sort_lms_suffixes_8u(T, SA, &buckets[4 * ALPHABET_SIZE], (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < threads)
        {
            #pragma omp master
            {
                libsais_radix_sort_lms_suffixes_8u(T, SA, &buckets[4 * ALPHABET_SIZE], (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
            }
        }
        else
        {
            {
                sa_sint_t * RESTRICT src_bucket = &buckets[4 * ALPHABET_SIZE];
                sa_sint_t * RESTRICT dst_bucket = thread_state[omp_thread_num].state.buckets;

                fast_sint_t i, j;
                for (i = BUCKETS_INDEX2(0, 0), j = BUCKETS_INDEX4(0, 1); i <= BUCKETS_INDEX2(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX2(1, 0), j += BUCKETS_INDEX4(1, 0))
                {
                    dst_bucket[i] = src_bucket[i] - dst_bucket[j];
                }
            }

            {
                fast_sint_t t, omp_block_start = 0, omp_block_size = thread_state[omp_thread_num].state.m;
                for (t = omp_num_threads - 1; t >= omp_thread_num; --t) omp_block_start += thread_state[t].state.m;

                if (omp_block_start == (fast_sint_t)m && omp_block_size > 0)
                {
                    omp_block_start -= 1; omp_block_size -= 1;
                }

                libsais_radix_sort_lms_suffixes_8u(T, SA, thread_state[omp_thread_num].state.buckets, (fast_sint_t)n - omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

static void libsais_radix_sort_lms_suffixes_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais_prefetchr(&SA[i - 3 * prefetch_distance]);
        
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 0]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 1]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 2]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 3]]);

        libsais_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 0]]]);
        libsais_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 1]]]);
        libsais_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 2]]]);
        libsais_prefetchw(&induction_bucket[T[SA[i - prefetch_distance - 3]]]);

        sa_sint_t p0 = SA[i - 0]; SA[--induction_bucket[T[p0]]] = p0;
        sa_sint_t p1 = SA[i - 1]; SA[--induction_bucket[T[p1]]] = p1;
        sa_sint_t p2 = SA[i - 2]; SA[--induction_bucket[T[p2]]] = p2;
        sa_sint_t p3 = SA[i - 3]; SA[--induction_bucket[T[p3]]] = p3;
    }

    for (j -= 2 * prefetch_distance + 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[--induction_bucket[T[p]]] = p;
    }
}

static void libsais_radix_sort_lms_suffixes_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais_prefetchr(&SA[i - 3 * prefetch_distance]);
        
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 0]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 1]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 2]]);
        libsais_prefetchr(&T[SA[i - 2 * prefetch_distance - 3]]);

        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 0]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 1]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 2]], 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(T[SA[i - prefetch_distance - 3]], 0)]);

        sa_sint_t p0 = SA[i - 0]; SA[--induction_bucket[BUCKETS_INDEX2(T[p0], 0)]] = p0;
        sa_sint_t p1 = SA[i - 1]; SA[--induction_bucket[BUCKETS_INDEX2(T[p1], 0)]] = p1;
        sa_sint_t p2 = SA[i - 2]; SA[--induction_bucket[BUCKETS_INDEX2(T[p2], 0)]] = p2;
        sa_sint_t p3 = SA[i - 3]; SA[--induction_bucket[BUCKETS_INDEX2(T[p3], 0)]] = p3;
    }

    for (j -= 2 * prefetch_distance + 3; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; SA[--induction_bucket[BUCKETS_INDEX2(T[p], 0)]] = p;
    }
}

#if defined(LIBSAIS_OPENMP)

static void libsais_radix_sort_lms_suffixes_32s_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i + prefetch_distance + 0]]);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 1]]);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 2]]);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 3]]);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        cache[i + 0].symbol = T[cache[i + 0].index = SA[i + 0]];
        cache[i + 1].symbol = T[cache[i + 1].index = SA[i + 1]];
        cache[i + 2].symbol = T[cache[i + 2].index = SA[i + 2]];
        cache[i + 3].symbol = T[cache[i + 3].index = SA[i + 3]];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        cache[i].symbol = T[cache[i].index = SA[i]];
    }
}

static void libsais_radix_sort_lms_suffixes_32s_6k_block_sort(sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        libsais_prefetchw(&induction_bucket[cache[i - prefetch_distance - 0].symbol]);
        libsais_prefetchw(&induction_bucket[cache[i - prefetch_distance - 1].symbol]);
        libsais_prefetchw(&induction_bucket[cache[i - prefetch_distance - 2].symbol]);
        libsais_prefetchw(&induction_bucket[cache[i - prefetch_distance - 3].symbol]);

        cache[i - 0].symbol = --induction_bucket[cache[i - 0].symbol];
        cache[i - 1].symbol = --induction_bucket[cache[i - 1].symbol];
        cache[i - 2].symbol = --induction_bucket[cache[i - 2].symbol];
        cache[i - 3].symbol = --induction_bucket[cache[i - 3].symbol];
    }

    for (j -= prefetch_distance + 3; i >= j; i -= 1)
    {
        cache[i].symbol = --induction_bucket[cache[i].symbol];
    }
}

static void libsais_radix_sort_lms_suffixes_32s_2k_block_sort(sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 3; i >= j; i -= 4)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(cache[i - prefetch_distance - 0].symbol, 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(cache[i - prefetch_distance - 1].symbol, 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(cache[i - prefetch_distance - 2].symbol, 0)]);
        libsais_prefetchw(&induction_bucket[BUCKETS_INDEX2(cache[i - prefetch_distance - 3].symbol, 0)]);

        cache[i - 0].symbol = --induction_bucket[BUCKETS_INDEX2(cache[i - 0].symbol, 0)];
        cache[i - 1].symbol = --induction_bucket[BUCKETS_INDEX2(cache[i - 1].symbol, 0)];
        cache[i - 2].symbol = --induction_bucket[BUCKETS_INDEX2(cache[i - 2].symbol, 0)];
        cache[i - 3].symbol = --induction_bucket[BUCKETS_INDEX2(cache[i - 3].symbol, 0)];
    }

    for (j -= prefetch_distance + 3; i >= j; i -= 1)
    {
        cache[i].symbol = --induction_bucket[BUCKETS_INDEX2(cache[i].symbol, 0)];
    }
}

static void libsais_radix_sort_lms_suffixes_32s_6k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_radix_sort_lms_suffixes_32s_6k(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_radix_sort_lms_suffixes_32s_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                libsais_radix_sort_lms_suffixes_32s_6k_block_sort(induction_bucket, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

static void libsais_radix_sort_lms_suffixes_32s_2k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_radix_sort_lms_suffixes_32s_2k(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_radix_sort_lms_suffixes_32s_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                libsais_radix_sort_lms_suffixes_32s_2k_block_sort(induction_bucket, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

#endif

static void libsais_radix_sort_lms_suffixes_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || m < 65536)
    {
        libsais_radix_sort_lms_suffixes_32s_6k(T, SA, induction_bucket, (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < (fast_sint_t)m - 1; block_start = block_end)
        {
            block_end = block_start + (fast_sint_t)threads * LIBSAIS_PER_THREAD_CACHE_SIZE; if (block_end >= m) { block_end = (fast_sint_t)m - 1; }

            libsais_radix_sort_lms_suffixes_32s_6k_block_omp(T, SA, induction_bucket, thread_state[0].state.cache, (fast_sint_t)n - block_end, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static void libsais_radix_sort_lms_suffixes_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || m < 65536)
    {
        libsais_radix_sort_lms_suffixes_32s_2k(T, SA, induction_bucket, (fast_sint_t)n - (fast_sint_t)m + 1, (fast_sint_t)m - 1);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < (fast_sint_t)m - 1; block_start = block_end)
        {
            block_end = block_start + (fast_sint_t)threads * LIBSAIS_PER_THREAD_CACHE_SIZE; if (block_end >= m) { block_end = (fast_sint_t)m - 1; }

            libsais_radix_sort_lms_suffixes_32s_2k_block_omp(T, SA, induction_bucket, thread_state[0].state.cache, (fast_sint_t)n - block_end, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static sa_sint_t libsais_radix_sort_lms_suffixes_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t             i = n - 2;
    sa_sint_t             m = 0;
    fast_uint_t           f0 = 1;
    fast_uint_t           f1 = 0;
    fast_sint_t           c0 = T[n - 1];
    fast_sint_t           c1 = 0;
    fast_sint_t           c2 = 0;

#if defined(LIBSAIS_LANES_32S)
    for (; i >= prefetch_distance + 63 && (((ptrdiff_t)(const void *)&T[i - 63]) & 63) != 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        if (f0 & ~f1) { SA[--buckets[c2 = c1]] = i + 1; m++; }
    }

    {
        libsais_lanes_32s_t prev = libsais_lanes_broadcast_32s(c0); sa_sint_t suffixes[64];

        for (; i >= prefetch_distance + 63; i -= 64)
        {
            fast_sint_t b; for (b = 0; b < 64; b += 64 / (fast_sint_t)sizeof(sa_sint_t)) { libsais_prefetchr(&T[i - 2 * prefetch_distance - b]); }

            libsais_lanes_32s_t c[64 / LIBSAIS_LANES_32S]; uint64_t types = libsais_classify_suffixes_32s(T, i, f0, prev, c); fast_sint_t s = 63;

            f1 = (fast_uint_t)(types & 1);
            if (f1 & ~f0) { SA[--buckets[c2 = c0]] = i + 1; m++; }

            for (b = 0; b < 64; b += 1) { libsais_prefetchw(&buckets[T[i - prefetch_distance - b]]); }

            for (b = libsais_store_lms_suffixes(suffixes, s, i, ~types & (types >> 1)); s > b; s -= 1)
            {
                sa_sint_t p = suffixes[s]; SA[--buckets[c2 = T[p]]] = p; m++;
            }

            prev = c[64 / LIBSAIS_LANES_32S - 1]; c0 = T[i - 63]; f0 = (fast_uint_t)(types >> 63);
        }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        if (f0 & ~f1) { SA[--buckets[c2 = c1]] = i + 1; m++; }
    }
#else
    for (; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetchr(&T[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[T[i - prefetch_distance - 0]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 1]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 2]]);
        libsais_prefetchw(&buckets[T[i - prefetch_distance - 3]]);

        c1 = T[i - 0]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); 
        if (f1 & ~f0) { SA[--buckets[c2 = c0]] = i + 1; m++; }
        
        c0 = T[i - 1]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); 
        if (f0 & ~f1) { SA[--buckets[c2 = c1]] = i - 0; m++; }

        c1 = T[i - 2]; f1 = (fast_uint_t)(c1 > (c0 - (fast_sint_t)(f0))); 
        if (f1 & ~f0) { SA[--buckets[c2 = c0]] = i - 1; m++; }

        c0 = T[i - 3]; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1))); 
        if (f0 & ~f1) { SA[--buckets[c2 = c1]] = i - 2; m++; }
    }

    for (; i >= 0; i -= 1)
    {
        c1 = c0; c0 = T[i]; f1 = f0; f0 = (fast_uint_t)(c0 > (c1 - (fast_sint_t)(f1)));
        if (f0 & ~f1) { SA[--buckets[c2 = c1]] = i + 1; m++; }
    }
#endif

    if (m > 1)
    {
        SA[buckets[c2]] = 0;
    }

    return m;
}

static void libsais_radix_sort_set_markers_32s_6k(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&induction_bucket[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SA[induction_bucket[i + prefetch_distance + 0]]);
        libsais_prefetchw(&SA[induction_bucket[i + prefetch_distance + 1]]);
        libsais_prefetchw(&SA[induction_bucket[i + prefetch_distance + 2]]);
        libsais_prefetchw(&SA[induction_bucket[i + prefetch_distance + 3]]);

        SA[induction_bucket[i + 0]] |= SAINT_MIN;
        SA[induction_bucket[i + 1]] |= SAINT_MIN;
        SA[induction_bucket[i + 2]] |= SAINT_MIN;
        SA[induction_bucket[i + 3]] |= SAINT_MIN;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[induction_bucket[i]] |= SAINT_MIN;
    }
}

static void libsais_radix_sort_set_markers_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&induction_bucket[BUCKETS_INDEX2(i + 2 * prefetch_distance, 0)]);

        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 0, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 1, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 2, 0)]]);
        libsais_prefetchw(&SA[induction_bucket[BUCKETS_INDEX2(i + prefetch_distance + 3, 0)]]);

        SA[induction_bucket[BUCKETS_INDEX2(i + 0, 0)]] |= SUFFIX_GROUP_MARKER;
        SA[induction_bucket[BUCKETS_INDEX2(i + 1, 0)]] |= SUFFIX_GROUP_MARKER;
        SA[induction_bucket[BUCKETS_INDEX2(i + 2, 0)]] |= SUFFIX_GROUP_MARKER;
        SA[induction_bucket[BUCKETS_INDEX2(i + 3, 0)]] |= SUFFIX_GROUP_MARKER;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[induction_bucket[BUCKETS_INDEX2(i, 0)]] |= SUFFIX_GROUP_MARKER;
    }
}

static void libsais_radix_sort_set_markers_32s_6k_omp(sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && k >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (((fast_sint_t)k - 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : (fast_sint_t)k - 1 - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = (fast_sint_t)k - 1;
#endif

        libsais_radix_sort_set_markers_32s_6k(SA, induction_bucket, omp_block_start, omp_block_size);
    }
}

static void libsais_radix_sort_set_markers_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && k >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (((fast_sint_t)k - 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : (fast_sint_t)k - 1 - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = (fast_sint_t)k - 1;
#endif

        libsais_radix_sort_set_markers_32s_4k(SA, induction_bucket, omp_block_start, omp_block_size);
    }
}

static void libsais_initialize_buckets_for_partial_sorting_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

    buckets[BUCKETS_INDEX4((fast_uint_t)T[first_lms_suffix], 1)]++;

    fast_sint_t i, j; sa_sint_t sum0 = left_suffixes_count + 1, sum1 = 0;
    for (i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX4(ALPHABET_SIZE - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;

        sum0 += buckets[i + BUCKETS_INDEX4(0, 0)] + buckets[i + BUCKETS_INDEX4(0, 2)];
        sum1 += buckets[i + BUCKETS_INDEX4(0, 1)];

        buckets[j + BUCKETS_INDEX2(0, 0)] = sum0;
        buckets[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static void libsais_initialize_buckets_for_partial_sorting_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

    fast_sint_t i, j; sa_sint_t sum0 = left_suffixes_count + 1, sum1 = 0, sum2 = 0;
    for (first_lms_suffix = T[first_lms_suffix], i = BUCKETS_INDEX4(0, 0), j = BUCKETS_INDEX2(0, 0); i != BUCKETS_INDEX4((fast_sint_t)first_lms_suffix, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    {
        sa_sint_t SS = buckets[i + BUCKETS_INDEX4(0, 0)];
        sa_sint_t LS = buckets[i + BUCKETS_INDEX4(0, 1)];
        sa_sint_t SL = buckets[i + BUCKETS_INDEX4(0, 2)];
        sa_sint_t LL = buckets[i + BUCKETS_INDEX4(0, 3)];

        buckets[i + BUCKETS_INDEX4(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX4(0, 1)] = sum2;
        buckets[i + BUCKETS_INDEX4(0, 2)] = 0;
        buckets[i + BUCKETS_INDEX4(0, 3)] = 0;

        sum0 += SS + SL; sum1 += LS; sum2 += LS + LL;

        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;
        temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }

    for (sum1 += 1; i <= BUCKETS_INDEX4((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX4(1, 0), j += BUCKETS_INDEX2(1, 0))
    { 
        sa_sint_t SS = buckets[i + BUCKETS_INDEX4(0, 0)];
        sa_sint_t LS = buckets[i + BUCKETS_INDEX4(0, 1)];
        sa_sint_t SL = buckets[i + BUCKETS_INDEX4(0, 2)];
        sa_sint_t LL = buckets[i + BUCKETS_INDEX4(0, 3)];

        buckets[i + BUCKETS_INDEX4(0, 0)] = sum0;
        buckets[i + BUCKETS_INDEX4(0, 1)] = sum2;
        buckets[i + BUCKETS_INDEX4(0, 2)] = 0;
        buckets[i + BUCKETS_INDEX4(0, 3)] = 0;

        sum0 += SS + SL; sum1 += LS; sum2 += LS + LL;

        temp_bucket[j + BUCKETS_INDEX2(0, 0)] = sum0;
        temp_bucket[j + BUCKETS_INDEX2(0, 1)] = sum1;
    }
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 2);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i + 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] >= T[p0 - 1]);
        SA[induction_bucket[v0]++] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = SA[i + 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] >= T[p1 - 1]);
        SA[induction_bucket[v1]++] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]);
        SA[induction_bucket[v]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

#if defined(LIBSAIS_OPENMP)

static void libsais_partial_sorting_scan_left_to_right_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size, LIBSAIS_THREAD_STATE * RESTRICT state)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    memset(induction_bucket, 0, (size_t)2 * (size_t)k * sizeof(sa_sint_t));
    memset(distinct_names  , 0, (size_t)2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t i, j, count = 0; sa_sint_t d = 1;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 0] & SAINT_MAX] - 2);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i + prefetch_distance + 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = cache[count].index = SA[i + 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = cache[count++].symbol = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] >= T[p0 - 1]); induction_bucket[v0]++; distinct_names[v0] = d;
        sa_sint_t p1 = cache[count].index = SA[i + 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = cache[count++].symbol = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] >= T[p1 - 1]); induction_bucket[v1]++; distinct_names[v1] = d;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = cache[count].index = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = cache[count++].symbol = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]); induction_bucket[v]++; distinct_names[v] = d;
    }

    state[0].state.position   = (fast_sint_t)d - 1;
    state[0].state.count      = count;
}

static void libsais_partial_sorting_scan_left_to_right_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count, sa_sint_t d)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = 0, j = count - 1; i < j; i += 2)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        sa_sint_t p0 = cache[i + 0].index; d += (p0 < 0); sa_sint_t v0 = cache[i + 0].symbol;
        SA[induction_bucket[v0]++] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = cache[i + 1].index; d += (p1 < 0); sa_sint_t v1 = cache[i + 1].symbol;
        SA[induction_bucket[v1]++] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j += 1; i < j; i += 1)
    {
        sa_sint_t p = cache[i].index; d += (p < 0); sa_sint_t v = cache[i].symbol;
        SA[induction_bucket[v]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_left_to_right_8u(T, SA, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                d = libsais_partial_sorting_scan_left_to_right_8u(T, SA, buckets, d, block_start, block_size);
            }
        }
        else
        {
            {
                libsais_partial_sorting_scan_left_to_right_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size, &thread_state[omp_thread_num]);
            }

            #pragma omp barrier

            #pragma omp master
            {
                sa_sint_t * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
                sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

                fast_sint_t t;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    sa_sint_t * RESTRICT temp_induction_bucket    = &thread_state[t].state.buckets[0 * ALPHABET_SIZE];
                    sa_sint_t * RESTRICT temp_distinct_names      = &thread_state[t].state.buckets[2 * ALPHABET_SIZE];

                    fast_sint_t c; 
                    for (c = 0; c < 2 * k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_induction_bucket[c]; induction_bucket[c] = A + B; temp_induction_bucket[c] = A; }

                    for (d -= 1, c = 0; c < 2 * k; c += 1) { sa_sint_t A = distinct_names[c], B = temp_distinct_names[c], D = B + d; distinct_names[c] = B > 0 ? D : A; temp_distinct_names[c] = A; }
                    d += 1 + (sa_sint_t)thread_state[t].state.position; thread_state[t].state.position = (fast_sint_t)d - thread_state[t].state.position;
                }
            }

            #pragma omp barrier

            {
                libsais_partial_sorting_scan_left_to_right_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count, (sa_sint_t)thread_state[omp_thread_num].state.position);
            }
        }
#endif
    }

    return d;
}

static fast_sint_t libsais_find_zero_left_to_right(const sa_sint_t * RESTRICT SA, fast_sint_t i, fast_sint_t end)
{
#if defined(LIBSAIS_LANES)
    for (; i < end && (((ptrdiff_t)(const void *)&SA[i]) & 63) != 0; i += 1) { if (SA[i] == 0) { return i; } }
    for (; i + 16 <= end && libsais_has_zero(&SA[i]) == 0; i += 16) { }
#endif

    while (i < end && SA[i] != 0) { i++; }

    return i;
}

static fast_sint_t libsais_find_block_end_left_to_right(const sa_sint_t * RESTRICT SA, fast_sint_t block_start, fast_sint_t block_max_end)
{
    if (block_max_end > block_start + 24576) { block_max_end = block_start + 24576; }

    return libsais_find_zero_left_to_right(SA, block_start + 1, block_max_end);
}

static fast_sint_t libsais_find_block_end_left_to_right_omp(const sa_sint_t * RESTRICT SA, fast_sint_t block_end, fast_sint_t block_max_end, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t block_scan_size = (fast_sint_t)threads * 4096;

    while (block_end < block_max_end && SA[block_end] != 0)
    {
        fast_sint_t block_scan_end = block_end + block_scan_size; if (block_scan_end > block_max_end) { block_scan_end = block_max_end; }

        #pragma omp parallel num_threads(threads) if(threads > 1 && omp_get_dynamic() == 0)
        {
            fast_sint_t omp_thread_num    = omp_get_thread_num();
            fast_sint_t omp_num_threads   = omp_get_num_threads();
            fast_sint_t omp_block_stride  = ((block_scan_end - block_end) / omp_num_threads) & (-16);
            fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
            fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_scan_end - block_end - omp_block_start;

            omp_block_start += block_end;

            {
                fast_sint_t i = libsais_find_zero_left_to_right(SA, omp_block_start, omp_block_start + omp_block_size);

                thread_state[omp_thread_num].state.position = i < omp_block_start + omp_block_size ? i : block_scan_end;
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t; for (block_end = block_scan_end, t = 0; t < omp_num_threads; ++t) { if (thread_state[t].state.position < block_end) { block_end = thread_state[t].state.position; } }
            }
        }

        block_scan_size += block_scan_size;
    }

    return block_end;
}

#endif

static sa_sint_t libsais_partial_sorting_scan_left_to_right_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t left_suffixes_count, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t * RESTRICT induction_bucket = &buckets[4 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    SA[induction_bucket[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | SAINT_MIN;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    if (threads == 1 || left_suffixes_count < 65536)
    {
        d = libsais_partial_sorting_scan_left_to_right_8u(T, SA, buckets, d, 0, left_suffixes_count);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = 0; block_start < left_suffixes_count; )
        {
            if (SA[block_start] == 0)
            {
                block_start++;
            }
            else
            {
                fast_sint_t block_max_end = block_start + ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end > left_suffixes_count) { block_max_end = left_suffixes_count;}
                fast_sint_t block_end     = libsais_find_block_end_left_to_right(SA, block_start, block_max_end); if (block_end < block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_left_to_right_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_end - block_start;

                if (block_size < 32)
                {
                    for (; block_start < block_end; block_start += 1)
                    {
                        sa_sint_t p = SA[block_start]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] >= T[p - 1]);
                        SA[induction_bucket[v]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
                    }
                }
                else
                {
                    d = libsais_partial_sorting_scan_left_to_right_8u_block_omp(T, SA, k, buckets, d, block_start, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchr(&SA[i + 3 * prefetch_distance]);

        const sa_sint_t * Ts0 = &T[SA[i + 2 * prefetch_distance + 0] & SAINT_MAX]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        const sa_sint_t * Ts1 = &T[SA[i + 2 * prefetch_distance + 1] & SAINT_MAX]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i + prefetch_distance + 0] & SAINT_MAX; libsais_prefetchw(&buckets[BUCKETS_INDEX4((fast_sint_t)T[(sa_uint_t)(p0 - (p0 > 0))], 0)]);
        sa_sint_t p1 = SA[i + prefetch_distance + 1] & SAINT_MAX; libsais_prefetchw(&buckets[BUCKETS_INDEX4((fast_sint_t)T[(sa_uint_t)(p1 - (p1 > 0))], 0)]);

        sa_sint_t p2 = SA[i + 0]; d += (p2 < 0); p2 &= SAINT_MAX; sa_sint_t v2 = BUCKETS_INDEX4(T[p2 - 1], T[p2 - 2] >= T[p2 - 1]);
        sa_sint_t * bucket_v2 = &buckets[v2]; SA[bucket_v2[0]++] = (p2 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v2[2] != d) << (SAINT_BIT - 1)); bucket_v2[2] = d;

        sa_sint_t p3 = SA[i + 1]; d += (p3 < 0); p3 &= SAINT_MAX; sa_sint_t v3 = BUCKETS_INDEX4(T[p3 - 1], T[p3 - 2] >= T[p3 - 1]);
        sa_sint_t * bucket_v3 = &buckets[v3]; SA[bucket_v3[0]++] = (p3 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v3[2] != d) << (SAINT_BIT - 1)); bucket_v3[2] = d;
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX4(T[p - 1], T[p - 2] >= T[p - 1]);
        sa_sint_t * bucket_v = &buckets[v]; SA[bucket_v[0]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(bucket_v[2] != d) << (SAINT_BIT - 1)); bucket_v[2] = d;
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT distinct_names   = &buckets[0 * (fast_sint_t)k];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i + 2 * prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + 2 * prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);
        sa_sint_t s2 = SA[i + 1 * prefetch_distance + 0]; if (s2 > 0) { const fast_sint_t Ts2 = T[(s2 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts2]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts2, 0)]); }
        sa_sint_t s3 = SA[i + 1 * prefetch_distance + 1]; if (s3 > 0) { const fast_sint_t Ts3 = T[(s3 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts3]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts3, 0)]); }

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 & SAINT_MAX;
        if (p0 > 0)
        {
            SA[i + 0] = 0; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); p0 &= ~SUFFIX_GROUP_MARKER; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] < T[p0 - 1]);
            SA[induction_bucket[T[p0 - 1]]++] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] < T[p0 - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
        }

        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 & SAINT_MAX;
        if (p1 > 0)
        {
            SA[i + 1] = 0; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); p1 &= ~SUFFIX_GROUP_MARKER; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] < T[p1 - 1]);
            SA[induction_bucket[T[p1 - 1]]++] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] < T[p1 - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
        }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p & SAINT_MAX;
        if (p > 0)
        {
            SA[i] = 0; d += (p >> (SUFFIX_GROUP_BIT - 1)); p &= ~SUFFIX_GROUP_MARKER; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] < T[p - 1]);
            SA[induction_bucket[T[p - 1]]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] < T[p - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i + 2 * prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 1)]; libsais_prefetchr(Ts0 - 1);
        sa_sint_t s1 = SA[i + 2 * prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 1)]; libsais_prefetchr(Ts1 - 1);
        sa_sint_t s2 = SA[i + 1 * prefetch_distance + 0]; const sa_sint_t * Ts2 = &T[(sa_uint_t)(s2 > 0 ? s2 : 2)]; libsais_prefetchw(&induction_bucket[Ts2[-1]]); libsais_prefetchr(Ts2 - 2);
        sa_sint_t s3 = SA[i + 1 * prefetch_distance + 1]; const sa_sint_t * Ts3 = &T[(sa_uint_t)(s3 > 0 ? s3 : 2)]; libsais_prefetchw(&induction_bucket[Ts3[-1]]); libsais_prefetchr(Ts3 - 2);

        sa_sint_t p0 = SA[i + 0]; if (p0 > 0) { SA[i + 0] = 0; SA[induction_bucket[T[p0 - 1]]++] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] < T[p0 - 1]) << (SAINT_BIT - 1)); } else { SA[i + 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i + 1]; if (p1 > 0) { SA[i + 1] = 0; SA[induction_bucket[T[p1 - 1]]++] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] < T[p1 - 1]) << (SAINT_BIT - 1)); } else { SA[i + 1] = p1 & SAINT_MAX; }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { SA[i] = 0; SA[induction_bucket[T[p - 1]]++] = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] < T[p - 1]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
    }
}

#if defined(LIBSAIS_OPENMP)

static void libsais_partial_sorting_scan_left_to_right_32s_6k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        const sa_sint_t * Ts0 = &T[SA[i + prefetch_distance + 0] & SAINT_MAX]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        const sa_sint_t * Ts1 = &T[SA[i + prefetch_distance + 1] & SAINT_MAX]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t p0 = cache[i + 0].index = SA[i + 0]; sa_sint_t symbol0 = 0; p0 &= SAINT_MAX; if (p0 != 0) { symbol0 = BUCKETS_INDEX4(T[p0 - 1], T[p0 - 2] >= T[p0 - 1]); } cache[i + 0].symbol = symbol0;
        sa_sint_t p1 = cache[i + 1].index = SA[i + 1]; sa_sint_t symbol1 = 0; p1 &= SAINT_MAX; if (p1 != 0) { symbol1 = BUCKETS_INDEX4(T[p1 - 1], T[p1 - 2] >= T[p1 - 1]); } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = cache[i].index = SA[i]; sa_sint_t symbol = 0; p &= SAINT_MAX; if (p != 0) { symbol = BUCKETS_INDEX4(T[p - 1], T[p - 2] >= T[p - 1]); } cache[i].symbol = symbol;
    }
}

static void libsais_partial_sorting_scan_left_to_right_32s_4k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; if (p0 > 0) { cache[i + 0].index = p0; p0 &= ~SUFFIX_GROUP_MARKER; symbol0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] < T[p0 - 1]); p0 = 0; } cache[i + 0].symbol = symbol0; SA[i + 0] = p0 & SAINT_MAX;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; if (p1 > 0) { cache[i + 1].index = p1; p1 &= ~SUFFIX_GROUP_MARKER; symbol1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] < T[p1 - 1]); p1 = 0; } cache[i + 1].symbol = symbol1; SA[i + 1] = p1 & SAINT_MAX;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; if (p > 0) { cache[i].index = p; p &= ~SUFFIX_GROUP_MARKER; symbol = BUCKETS_INDEX2(T[p - 1], T[p - 2] < T[p - 1]); p = 0; } cache[i].symbol = symbol; SA[i] = p & SAINT_MAX;
    }
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; if (p0 > 0) { cache[i + 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] < T[p0 - 1]) << (SAINT_BIT - 1)); symbol0 = T[p0 - 1]; p0 = 0; } cache[i + 0].symbol = symbol0; SA[i + 0] = p0 & SAINT_MAX;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; if (p1 > 0) { cache[i + 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] < T[p1 - 1]) << (SAINT_BIT - 1)); symbol1 = T[p1 - 1]; p1 = 0; } cache[i + 1].symbol = symbol1; SA[i + 1] = p1 & SAINT_MAX;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; if (p > 0) { cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] < T[p - 1]) << (SAINT_BIT - 1)); symbol = T[p - 1]; p = 0; } cache[i].symbol = symbol; SA[i] = p & SAINT_MAX;
    }
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_6k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, omp_block_end = omp_block_start + omp_block_size;
    for (i = omp_block_start, j = omp_block_end - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&cache[i + 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[cache[i + prefetch_distance + 0].symbol]);
        libsais_prefetchw(&buckets[cache[i + prefetch_distance + 1].symbol]);

        sa_sint_t v0 = cache[i + 0].symbol, p0 = cache[i + 0].index; d += (p0 < 0); sa_sint_t * bucket_v0 = &buckets[v0]; cache[i + 0].symbol = bucket_v0[0]++; cache[i + 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v0[2] != d) << (SAINT_BIT - 1)); bucket_v0[2] = d;
        if (cache[i + 0].symbol < omp_block_end) { sa_sint_t s = cache[i + 0].symbol, q = (cache[s].index = cache[i + 0].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] >= T[q - 1]); }

        sa_sint_t v1 = cache[i + 1].symbol, p1 = cache[i + 1].index; d += (p1 < 0); sa_sint_t * bucket_v1 = &buckets[v1]; cache[i + 1].symbol = bucket_v1[0]++; cache[i + 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v1[2] != d) << (SAINT_BIT - 1)); bucket_v1[2] = d;
        if (cache[i + 1].symbol < omp_block_end) { sa_sint_t s = cache[i + 1].symbol, q = (cache[s].index = cache[i + 1].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] >= T[q - 1]); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t v = cache[i].symbol, p = cache[i].index; d += (p < 0); sa_sint_t * bucket_v = &buckets[v]; cache[i].symbol = bucket_v[0]++; cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)(bucket_v[2] != d) << (SAINT_BIT - 1)); bucket_v[2] = d;
        if (cache[i].symbol < omp_block_end) { sa_sint_t s = cache[i].symbol, q = (cache[s].index = cache[i].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] >= T[q - 1]); }
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_4k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT distinct_names   = &buckets[0 * (fast_sint_t)k];

    fast_sint_t i, j, omp_block_end = omp_block_start + omp_block_size;
    for (i = omp_block_start, j = omp_block_end - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&cache[i + 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i + prefetch_distance + 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 >> 1 : 0]; libsais_prefetchw(Is0); const sa_sint_t * Ds0 = &distinct_names[s0 > 0 ? s0 : 0]; libsais_prefetchw(Ds0); 
        sa_sint_t s1 = cache[i + prefetch_distance + 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 >> 1 : 0]; libsais_prefetchw(Is1); const sa_sint_t * Ds1 = &distinct_names[s1 > 0 ? s1 : 0]; libsais_prefetchw(Ds1);
        
        sa_sint_t v0 = cache[i + 0].symbol;
        if (v0 >= 0)
        {
            sa_sint_t p0 = cache[i + 0].index; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); cache[i + 0].symbol = induction_bucket[v0 >> 1]++; cache[i + 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)v0 << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
            if (cache[i + 0].symbol < omp_block_end) { sa_sint_t ni = cache[i + 0].symbol, np = cache[i + 0].index; if (np > 0) { cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] < T[np - 1]); np = 0; } cache[i + 0].index = np & SAINT_MAX; }
        }

        sa_sint_t v1 = cache[i + 1].symbol;
        if (v1 >= 0)
        {
            sa_sint_t p1 = cache[i + 1].index; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); cache[i + 1].symbol = induction_bucket[v1 >> 1]++; cache[i + 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)v1 << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
            if (cache[i + 1].symbol < omp_block_end) { sa_sint_t ni = cache[i + 1].symbol, np = cache[i + 1].index; if (np > 0) { cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] < T[np - 1]); np = 0; } cache[i + 1].index = np & SAINT_MAX; }
        }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            sa_sint_t p = cache[i].index; d += (p >> (SUFFIX_GROUP_BIT - 1)); cache[i].symbol = induction_bucket[v >> 1]++; cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)v << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
            if (cache[i].symbol < omp_block_end) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; if (np > 0) { cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] < T[np - 1]); np = 0; } cache[i].index = np & SAINT_MAX; }
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, omp_block_end = omp_block_start + omp_block_size;
    for (i = omp_block_start, j = omp_block_end - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&cache[i + 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i + prefetch_distance + 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 : 0]; libsais_prefetchw(Is0);
        sa_sint_t s1 = cache[i + prefetch_distance + 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 : 0]; libsais_prefetchw(Is1);
        
        sa_sint_t v0 = cache[i + 0].symbol;
        if (v0 >= 0)
        {
            cache[i + 0].symbol = induction_bucket[v0]++;
            if (cache[i + 0].symbol < omp_block_end) { sa_sint_t ni = cache[i + 0].symbol, np = cache[i + 0].index; if (np > 0) { cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] < T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; np = 0; } cache[i + 0].index = np & SAINT_MAX; }
        }

        sa_sint_t v1 = cache[i + 1].symbol;
        if (v1 >= 0)
        {
            cache[i + 1].symbol = induction_bucket[v1]++;
            if (cache[i + 1].symbol < omp_block_end) { sa_sint_t ni = cache[i + 1].symbol, np = cache[i + 1].index; if (np > 0) { cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] < T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; np = 0; } cache[i + 1].index = np & SAINT_MAX; }
        }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            cache[i].symbol = induction_bucket[v]++;
            if (cache[i].symbol < omp_block_end) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; if (np > 0) { cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] < T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; np = 0; } cache[i].index = np & SAINT_MAX; }
        }
    }
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_6k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_left_to_right_32s_6k(T, SA, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_left_to_right_32s_6k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                d = libsais_partial_sorting_scan_left_to_right_32s_6k_block_sort(T, buckets, d, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_4k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_left_to_right_32s_4k(T, SA, k, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_left_to_right_32s_4k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                d = libsais_partial_sorting_scan_left_to_right_32s_4k_block_sort(T, k, buckets, d, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return d;
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_partial_sorting_scan_left_to_right_32s_1k(T, SA, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_left_to_right_32s_1k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                libsais_partial_sorting_scan_left_to_right_32s_1k_block_sort(T, buckets, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

#endif

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t left_suffixes_count, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[buckets[BUCKETS_INDEX4(T[n - 1], T[n - 2] >= T[n - 1])]++] = (n - 1) | SAINT_MIN;
    buckets[2 + BUCKETS_INDEX4(T[n - 1], T[n - 2] >= T[n - 1])] = ++d;

    if (threads == 1 || left_suffixes_count < 65536)
    {
        d = libsais_partial_sorting_scan_left_to_right_32s_6k(T, SA, buckets, d, 0, left_suffixes_count);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < left_suffixes_count; block_start = block_end)
        {
            block_end = block_start + libsais_get_block_size_32s(n, threads); if (block_end > left_suffixes_count) { block_end = left_suffixes_count; }

            d = libsais_partial_sorting_scan_left_to_right_32s_6k_block_omp(T, SA, buckets, d, thread_state[0].state.cache, block_start, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_left_to_right_32s_4k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t * RESTRICT induction_bucket = &buckets[2 * (fast_sint_t)k];
    sa_sint_t * RESTRICT distinct_names   = &buckets[0 * (fast_sint_t)k];

    SA[induction_bucket[T[n - 1]]++] = (n - 1) | (sa_sint_t)((sa_uint_t)(T[n - 2] < T[n - 1]) << (SAINT_BIT - 1)) | SUFFIX_GROUP_MARKER;
    distinct_names[BUCKETS_INDEX2(T[n - 1], T[n - 2] < T[n - 1])] = ++d;

    if (threads == 1 || n < 65536)
    {
        d = libsais_partial_sorting_scan_left_to_right_32s_4k(T, SA, k, buckets, d, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < n; block_start = block_end)
        {
            block_end = block_start + libsais_get_block_size_32s(n, threads); if (block_end > n) { block_end = n; }

            d = libsais_partial_sorting_scan_left_to_right_32s_4k_block_omp(T, SA, k, buckets, d, thread_state[0].state.cache, block_start, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif

    return d;
}

static void libsais_partial_sorting_scan_left_to_right_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[buckets[T[n - 1]]++] = (n - 1) | (sa_sint_t)((sa_uint_t)(T[n - 2] < T[n - 1]) << (SAINT_BIT - 1));

    if (threads == 1 || n < 65536)
    {
       libsais_partial_sorting_scan_left_to_right_32s_1k(T, SA, buckets, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < n; block_start = block_end)
        {
            block_end = block_start + libsais_get_block_size_32s(n, threads); if (block_end > n) { block_end = n; }

            libsais_partial_sorting_scan_left_to_right_32s_1k_block_omp(T, SA, buckets, thread_state[0].state.cache, block_start, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static void libsais_partial_sorting_shift_markers_8u_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, const sa_sint_t * RESTRICT buckets, sa_sint_t threads)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT temp_bucket = &buckets[4 * ALPHABET_SIZE];

    fast_sint_t c;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel for schedule(static, 1) num_threads(threads) if(threads > 1 && n >= 65536)
#else
    UNUSED(threads); UNUSED(n);
#endif
    for (c = BUCKETS_INDEX2(ALPHABET_SIZE - 1, 0); c >= BUCKETS_INDEX2(1, 0); c -= BUCKETS_INDEX2(1, 0))
    {
        fast_sint_t i, j; sa_sint_t s = SAINT_MIN;
        for (i = (fast_sint_t)temp_bucket[c] - 1, j = (fast_sint_t)buckets[c - BUCKETS_INDEX2(1, 0)] + 3; i >= j; i -= 4)
        {
            libsais_prefetchw(&SA[i - prefetch_distance]);

            sa_sint_t p0 = SA[i - 0]; SA[i - 0] = (p0 & SAINT_MAX) | s; s = p0 & SAINT_MIN;
            sa_sint_t p1 = SA[i - 1]; SA[i - 1] = (p1 & SAINT_MAX) | s; s = p1 & SAINT_MIN;
            sa_sint_t p2 = SA[i - 2]; SA[i - 2] = (p2 & SAINT_MAX) | s; s = p2 & SAINT_MIN;
            sa_sint_t p3 = SA[i - 3]; SA[i - 3] = (p3 & SAINT_MAX) | s; s = p3 & SAINT_MIN;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            sa_sint_t p = SA[i]; SA[i] = (p & SAINT_MAX) | s; s = p & SAINT_MIN;
        }
    }
}

static void libsais_partial_sorting_shift_markers_32s_6k_omp(sa_sint_t * RESTRICT SA, sa_sint_t k, const sa_sint_t * RESTRICT buckets, sa_sint_t threads)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];
    
    fast_sint_t c;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel for schedule(static, 1) num_threads(threads) if(threads > 1 && k >= 65536)
#else
    UNUSED(threads);
#endif
    for (c = (fast_sint_t)k - 1; c >= 1; c -= 1)
    {
        fast_sint_t i, j; sa_sint_t s = SAINT_MIN;
        for (i = (fast_sint_t)buckets[BUCKETS_INDEX4(c, 0)] - 1, j = (fast_sint_t)temp_bucket[BUCKETS_INDEX2(c - 1, 0)] + 3; i >= j; i -= 4)
        {
            libsais_prefetchw(&SA[i - prefetch_distance]);

            sa_sint_t p0 = SA[i - 0]; SA[i - 0] = (p0 & SAINT_MAX) | s; s = p0 & SAINT_MIN;
            sa_sint_t p1 = SA[i - 1]; SA[i - 1] = (p1 & SAINT_MAX) | s; s = p1 & SAINT_MIN;
            sa_sint_t p2 = SA[i - 2]; SA[i - 2] = (p2 & SAINT_MAX) | s; s = p2 & SAINT_MIN;
            sa_sint_t p3 = SA[i - 3]; SA[i - 3] = (p3 & SAINT_MAX) | s; s = p3 & SAINT_MIN;
        }

        for (j -= 3; i >= j; i -= 1)
        {
            sa_sint_t p = SA[i]; SA[i] = (p & SAINT_MAX) | s; s = p & SAINT_MIN;
        }
    }
}

static void libsais_partial_sorting_shift_markers_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i; sa_sint_t s = SUFFIX_GROUP_MARKER;
    for (i = (fast_sint_t)n - 1; i >= 3; i -= 4)
    {
        libsais_prefetchw(&SA[i - prefetch_distance]);

        sa_sint_t p0 = SA[i - 0], q0 = ((p0 & SUFFIX_GROUP_MARKER) ^ s) & ((sa_sint_t)(p0 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q0; SA[i - 0] = p0 ^ q0;
        sa_sint_t p1 = SA[i - 1], q1 = ((p1 & SUFFIX_GROUP_MARKER) ^ s) & ((sa_sint_t)(p1 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q1; SA[i - 1] = p1 ^ q1;
        sa_sint_t p2 = SA[i - 2], q2 = ((p2 & SUFFIX_GROUP_MARKER) ^ s) & ((sa_sint_t)(p2 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q2; SA[i - 2] = p2 ^ q2;
        sa_sint_t p3 = SA[i - 3], q3 = ((p3 & SUFFIX_GROUP_MARKER) ^ s) & ((sa_sint_t)(p3 > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q3; SA[i - 3] = p3 ^ q3;
    }

    for (; i >= 0; i -= 1)
    {
        sa_sint_t p = SA[i], q = ((p & SUFFIX_GROUP_MARKER) ^ s) & ((sa_sint_t)(p > 0) << ((SUFFIX_GROUP_BIT - 1))); s = s ^ q; SA[i] = p ^ q;
    }
}

static void libsais_partial_sorting_shift_buckets_32s_6k(sa_sint_t k, sa_sint_t * RESTRICT buckets)
{
    sa_sint_t * RESTRICT temp_bucket = &buckets[4 * (fast_sint_t)k];

    fast_sint_t i;
    for (i = BUCKETS_INDEX2(0, 0); i <= BUCKETS_INDEX2((fast_sint_t)k - 1, 0); i += BUCKETS_INDEX2(1, 0))
    {
        buckets[2 * i + BUCKETS_INDEX4(0, 0)] = temp_bucket[i + BUCKETS_INDEX2(0, 0)];
        buckets[2 * i + BUCKETS_INDEX4(0, 1)] = temp_bucket[i + BUCKETS_INDEX2(0, 1)];
    }
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 2);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i - 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
        SA[--induction_bucket[v0]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = SA[i - 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
        SA[--induction_bucket[v1]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
        SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }

    return d;
}

static sa_sint_t libsais_partial_gsa_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 2);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = SA[i - 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
        if (v0 != 1) { SA[--induction_bucket[v0]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d; }

        sa_sint_t p1 = SA[i - 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
        if (v1 != 1) { SA[--induction_bucket[v1]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
        if (v != 1) { SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d; }
    }

    return d;
}


#if defined(LIBSAIS_OPENMP)

static void libsais_partial_sorting_scan_right_to_left_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size, LIBSAIS_THREAD_STATE * RESTRICT state)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    memset(induction_bucket, 0, (size_t)2 * (size_t)k * sizeof(sa_sint_t));
    memset(distinct_names  , 0, (size_t)2 * (size_t)k * sizeof(sa_sint_t));

    fast_sint_t i, j, count = 0; sa_sint_t d = 1;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 0] & SAINT_MAX] - 2);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 1);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1] & SAINT_MAX] - 2);

        sa_sint_t p0 = cache[count].index = SA[i - 0]; d += (p0 < 0); p0 &= SAINT_MAX; sa_sint_t v0 = cache[count++].symbol = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]); induction_bucket[v0]++; distinct_names[v0] = d;
        sa_sint_t p1 = cache[count].index = SA[i - 1]; d += (p1 < 0); p1 &= SAINT_MAX; sa_sint_t v1 = cache[count++].symbol = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]); induction_bucket[v1]++; distinct_names[v1] = d;
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = cache[count].index = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = cache[count++].symbol = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]); induction_bucket[v]++; distinct_names[v] = d;
    }

    state[0].state.position   = (fast_sint_t)d - 1;
    state[0].state.count      = count;
}

static void libsais_partial_sorting_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count, sa_sint_t d)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = 0, j = count - 1; i < j; i += 2)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        sa_sint_t p0 = cache[i + 0].index; d += (p0 < 0); sa_sint_t v0 = cache[i + 0].symbol;
        SA[--induction_bucket[v0]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d;

        sa_sint_t p1 = cache[i + 1].index; d += (p1 < 0); sa_sint_t v1 = cache[i + 1].symbol;
        SA[--induction_bucket[v1]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d;
    }

    for (j += 1; i < j; i += 1)
    {
        sa_sint_t p = cache[i].index; d += (p < 0); sa_sint_t v = cache[i].symbol;
        SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
    }
}

static void libsais_partial_gsa_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count, sa_sint_t d)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
    sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

    fast_sint_t i, j;
    for (i = 0, j = count - 1; i < j; i += 2)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        sa_sint_t p0 = cache[i + 0].index; d += (p0 < 0); sa_sint_t v0 = cache[i + 0].symbol;
        if (v0 != 1) { SA[--induction_bucket[v0]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v0] != d) << (SAINT_BIT - 1)); distinct_names[v0] = d; }

        sa_sint_t p1 = cache[i + 1].index; d += (p1 < 0); sa_sint_t v1 = cache[i + 1].symbol;
        if (v1 != 1) { SA[--induction_bucket[v1]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v1] != d) << (SAINT_BIT - 1)); distinct_names[v1] = d; }
    }

    for (j += 1; i < j; i += 1)
    {
        sa_sint_t p = cache[i].index; d += (p < 0); sa_sint_t v = cache[i].symbol;
        if (v != 1) { SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d; }
    }
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_right_to_left_8u(T, SA, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                d = libsais_partial_sorting_scan_right_to_left_8u(T, SA, buckets, d, block_start, block_size);
            }
        }
        else
        {
            {
                libsais_partial_sorting_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size, &thread_state[omp_thread_num]);
            }

            #pragma omp barrier

            #pragma omp master
            {
                sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
                sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_induction_bucket    = &thread_state[t].state.buckets[0 * ALPHABET_SIZE];
                    sa_sint_t * RESTRICT temp_distinct_names      = &thread_state[t].state.buckets[2 * ALPHABET_SIZE];

                    fast_sint_t c; 
                    for (c = 0; c < 2 * k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_induction_bucket[c]; induction_bucket[c] = A - B; temp_induction_bucket[c] = A; }

                    for (d -= 1, c = 0; c < 2 * k; c += 1) { sa_sint_t A = distinct_names[c], B = temp_distinct_names[c], D = B + d; distinct_names[c] = B > 0 ? D : A; temp_distinct_names[c] = A; }
                    d += 1 + (sa_sint_t)thread_state[t].state.position; thread_state[t].state.position = (fast_sint_t)d - thread_state[t].state.position;
                }
            }

            #pragma omp barrier

            {
                libsais_partial_sorting_scan_right_to_left_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count, (sa_sint_t)thread_state[omp_thread_num].state.position);
            }
        }
#endif
    }

    return d;
}

static sa_sint_t libsais_partial_gsa_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_gsa_scan_right_to_left_8u(T, SA, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                d = libsais_partial_gsa_scan_right_to_left_8u(T, SA, buckets, d, block_start, block_size);
            }
        }
        else
        {
            {
                libsais_partial_sorting_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size, &thread_state[omp_thread_num]);
            }

            #pragma omp barrier

            #pragma omp master
            {
                sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
                sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_induction_bucket    = &thread_state[t].state.buckets[0 * ALPHABET_SIZE];
                    sa_sint_t * RESTRICT temp_distinct_names      = &thread_state[t].state.buckets[2 * ALPHABET_SIZE];

                    fast_sint_t c; 
                    for (c = 0; c < 2 * k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_induction_bucket[c]; induction_bucket[c] = A - B; temp_induction_bucket[c] = A; }

                    for (d -= 1, c = 0; c < 2 * k; c += 1) { sa_sint_t A = distinct_names[c], B = temp_distinct_names[c], D = B + d; distinct_names[c] = B > 0 ? D : A; temp_distinct_names[c] = A; }
                    d += 1 + (sa_sint_t)thread_state[t].state.position; thread_state[t].state.position = (fast_sint_t)d - thread_state[t].state.position;
                }
            }

            #pragma omp barrier

            {
                libsais_partial_gsa_scan_right_to_left_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count, (sa_sint_t)thread_state[omp_thread_num].state.position);
            }
        }
#endif
    }

    return d;
}

static fast_sint_t libsais_find_zero_right_to_left(const sa_sint_t * RESTRICT SA, fast_sint_t i, fast_sint_t end)
{
#if defined(LIBSAIS_LANES)
    for (; i > end && (((ptrdiff_t)(const void *)&SA[i + 1]) & 63) != 0; i -= 1) { if (SA[i] == 0) { return i; } }
    for (; i - 16 >= end && libsais_has_zero(&SA[i - 15]) == 0; i -= 16) { }
#endif

    while (i > end && SA[i] != 0) { i--; }

    return i;
}

static fast_sint_t libsais_find_block_end_right_to_left(const sa_sint_t * RESTRICT SA, fast_sint_t block_start, fast_sint_t block_max_end)
{
    if (block_max_end < block_start - 24576) { block_max_end = block_start - 24576; }

    return libsais_find_zero_right_to_left(SA, block_start - 1, block_max_end);
}

static fast_sint_t libsais_find_block_end_right_to_left_omp(const sa_sint_t * RESTRICT SA, fast_sint_t block_end, fast_sint_t block_max_end, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t block_scan_size = (fast_sint_t)threads * 4096;

    while (block_end > block_max_end && SA[block_end] != 0)
    {
        fast_sint_t block_scan_end = block_end - block_scan_size; if (block_scan_end < block_max_end) { block_scan_end = block_max_end; }

        #pragma omp parallel num_threads(threads) if(threads > 1 && omp_get_dynamic() == 0)
        {
            fast_sint_t omp_thread_num    = omp_get_thread_num();
            fast_sint_t omp_num_threads   = omp_get_num_threads();
            fast_sint_t omp_block_stride  = ((block_end - block_scan_end) / omp_num_threads) & (-16);
            fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
            fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_end - block_scan_end - omp_block_start;

            omp_block_start += block_scan_end + 1;

            {
                fast_sint_t i = libsais_find_zero_right_to_left(SA, omp_block_start + omp_block_size - 1, omp_block_start - 1);

                thread_state[omp_thread_num].state.position = i >= omp_block_start ? i : block_scan_end;
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t; for (block_end = block_scan_end, t = 0; t < omp_num_threads; ++t) { if (thread_state[t].state.position > block_end) { block_end = thread_state[t].state.position; } }
            }
        }

        block_scan_size += block_scan_size;
    }

    return block_end;
}

#endif

static void libsais_partial_sorting_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t scan_start    = (fast_sint_t)left_suffixes_count + 1;
    fast_sint_t scan_end      = (fast_sint_t)n - (fast_sint_t)first_lms_suffix;

    if (threads == 1 || (scan_end - scan_start) < 65536)
    {
        libsais_partial_sorting_scan_right_to_left_8u(T, SA, buckets, d, scan_start, scan_end - scan_start);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
        sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

        fast_sint_t block_start;
        for (block_start = scan_end - 1; block_start >= scan_start; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end < scan_start) { block_max_end = scan_start - 1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
                        SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d;
                    }
                }
                else
                {
                    d = libsais_partial_sorting_scan_right_to_left_8u_block_omp(T, SA, k, buckets, d, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static void libsais_partial_gsa_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t scan_start    = (fast_sint_t)left_suffixes_count + 1;
    fast_sint_t scan_end      = (fast_sint_t)n - (fast_sint_t)first_lms_suffix;

    if (threads == 1 || (scan_end - scan_start) < 65536)
    {
        libsais_partial_gsa_scan_right_to_left_8u(T, SA, buckets, d, scan_start, scan_end - scan_start);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        sa_sint_t * RESTRICT induction_bucket = &buckets[0 * ALPHABET_SIZE];
        sa_sint_t * RESTRICT distinct_names   = &buckets[2 * ALPHABET_SIZE];

        fast_sint_t block_start;
        for (block_start = scan_end - 1; block_start >= scan_start; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end < scan_start) { block_max_end = scan_start - 1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
                        if (v != 1) { SA[--induction_bucket[v]] = (p - 1) | (sa_sint_t)((sa_uint_t)(distinct_names[v] != d) << (SAINT_BIT - 1)); distinct_names[v] = d; }
                    }
                }
                else
                {
                    d = libsais_partial_gsa_scan_right_to_left_8u_block_omp(T, SA, k, buckets, d, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchr(&SA[i - 3 * prefetch_distance]);

        const sa_sint_t * Ts0 = &T[SA[i - 2 * prefetch_distance - 0] & SAINT_MAX]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        const sa_sint_t * Ts1 = &T[SA[i - 2 * prefetch_distance - 1] & SAINT_MAX]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - prefetch_distance - 0] & SAINT_MAX; libsais_prefetchw(&buckets[BUCKETS_INDEX4((fast_sint_t)T[(sa_uint_t)(p0 - (p0 > 0))], 0)]);
        sa_sint_t p1 = SA[i - prefetch_distance - 1] & SAINT_MAX; libsais_prefetchw(&buckets[BUCKETS_INDEX4((fast_sint_t)T[(sa_uint_t)(p1 - (p1 > 0))], 0)]);

        sa_sint_t p2 = SA[i - 0]; d += (p2 < 0); p2 &= SAINT_MAX; sa_sint_t v2 = BUCKETS_INDEX4(T[p2 - 1], T[p2 - 2] > T[p2 - 1]);
        sa_sint_t * bucket_v2 = &buckets[v2]; SA[--bucket_v2[0]] = (p2 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v2[2] != d) << (SAINT_BIT - 1)); bucket_v2[2] = d;

        sa_sint_t p3 = SA[i - 1]; d += (p3 < 0); p3 &= SAINT_MAX; sa_sint_t v3 = BUCKETS_INDEX4(T[p3 - 1], T[p3 - 2] > T[p3 - 1]);
        sa_sint_t * bucket_v3 = &buckets[v3]; SA[--bucket_v3[0]] = (p3 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v3[2] != d) << (SAINT_BIT - 1)); bucket_v3[2] = d;
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; d += (p < 0); p &= SAINT_MAX; sa_sint_t v = BUCKETS_INDEX4(T[p - 1], T[p - 2] > T[p - 1]);
        sa_sint_t * bucket_v = &buckets[v]; SA[--bucket_v[0]] = (p - 1) | (sa_sint_t)((sa_uint_t)(bucket_v[2] != d) << (SAINT_BIT - 1)); bucket_v[2] = d;
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[3 * (fast_sint_t)k];
    sa_sint_t * RESTRICT distinct_names   = &buckets[0 * (fast_sint_t)k];

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i - 2 * prefetch_distance - 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - 2 * prefetch_distance - 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);
        sa_sint_t s2 = SA[i - 1 * prefetch_distance - 0]; if (s2 > 0) { const fast_sint_t Ts2 = T[(s2 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts2]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts2, 0)]); }
        sa_sint_t s3 = SA[i - 1 * prefetch_distance - 1]; if (s3 > 0) { const fast_sint_t Ts3 = T[(s3 & ~SUFFIX_GROUP_MARKER) - 1]; libsais_prefetchw(&induction_bucket[Ts3]); libsais_prefetchw(&distinct_names[BUCKETS_INDEX2(Ts3, 0)]); }

        sa_sint_t p0 = SA[i - 0];
        if (p0 > 0)
        {
            SA[i - 0] = 0; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); p0 &= ~SUFFIX_GROUP_MARKER; sa_sint_t v0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]);
            SA[--induction_bucket[T[p0 - 1]]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] > T[p0 - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
        }

        sa_sint_t p1 = SA[i - 1];
        if (p1 > 0)
        {
            SA[i - 1] = 0; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); p1 &= ~SUFFIX_GROUP_MARKER; sa_sint_t v1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]);
            SA[--induction_bucket[T[p1 - 1]]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] > T[p1 - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
        }
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i];
        if (p > 0)
        {
            SA[i] = 0; d += (p >> (SUFFIX_GROUP_BIT - 1)); p &= ~SUFFIX_GROUP_MARKER; sa_sint_t v = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]);
            SA[--induction_bucket[T[p - 1]]] = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] > T[p - 1]) << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i - 2 * prefetch_distance - 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 1)]; libsais_prefetchr(Ts0 - 1);
        sa_sint_t s1 = SA[i - 2 * prefetch_distance - 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 1)]; libsais_prefetchr(Ts1 - 1);
        sa_sint_t s2 = SA[i - 1 * prefetch_distance - 0]; const sa_sint_t * Ts2 = &T[(sa_uint_t)(s2 > 0 ? s2 : 2)]; libsais_prefetchw(&induction_bucket[Ts2[-1]]); libsais_prefetchr(Ts2 - 2);
        sa_sint_t s3 = SA[i - 1 * prefetch_distance - 1]; const sa_sint_t * Ts3 = &T[(sa_uint_t)(s3 > 0 ? s3 : 2)]; libsais_prefetchw(&induction_bucket[Ts3[-1]]); libsais_prefetchr(Ts3 - 2);

        sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { SA[i - 0] = 0; SA[--induction_bucket[T[p0 - 1]]] = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] > T[p0 - 1]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { SA[i - 1] = 0; SA[--induction_bucket[T[p1 - 1]]] = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] > T[p1 - 1]) << (SAINT_BIT - 1)); }
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { SA[i] = 0; SA[--induction_bucket[T[p - 1]]] = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] > T[p - 1]) << (SAINT_BIT - 1)); }
    }
}

#if defined(LIBSAIS_OPENMP)

static void libsais_partial_sorting_scan_right_to_left_32s_6k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        const sa_sint_t * Ts0 = &T[SA[i + prefetch_distance + 0] & SAINT_MAX]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        const sa_sint_t * Ts1 = &T[SA[i + prefetch_distance + 1] & SAINT_MAX]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t p0 = cache[i + 0].index = SA[i + 0]; sa_sint_t symbol0 = 0; p0 &= SAINT_MAX; if (p0 != 0) { symbol0 = BUCKETS_INDEX4(T[p0 - 1], T[p0 - 2] > T[p0 - 1]); } cache[i + 0].symbol = symbol0;
        sa_sint_t p1 = cache[i + 1].index = SA[i + 1]; sa_sint_t symbol1 = 0; p1 &= SAINT_MAX; if (p1 != 0) { symbol1 = BUCKETS_INDEX4(T[p1 - 1], T[p1 - 2] > T[p1 - 1]); } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = cache[i].index = SA[i]; sa_sint_t symbol = 0; p &= SAINT_MAX; if (p != 0) { symbol = BUCKETS_INDEX4(T[p - 1], T[p - 2] > T[p - 1]); } cache[i].symbol = symbol;
    }
}

static void libsais_partial_sorting_scan_right_to_left_32s_4k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 & ~SUFFIX_GROUP_MARKER : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; if (p0 > 0) { SA[i + 0] = 0; cache[i + 0].index = p0; p0 &= ~SUFFIX_GROUP_MARKER; symbol0 = BUCKETS_INDEX2(T[p0 - 1], T[p0 - 2] > T[p0 - 1]); } cache[i + 0].symbol = symbol0;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; if (p1 > 0) { SA[i + 1] = 0; cache[i + 1].index = p1; p1 &= ~SUFFIX_GROUP_MARKER; symbol1 = BUCKETS_INDEX2(T[p1 - 1], T[p1 - 2] > T[p1 - 1]); } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; if (p > 0) { SA[i] = 0; cache[i].index = p; p &= ~SUFFIX_GROUP_MARKER; symbol = BUCKETS_INDEX2(T[p - 1], T[p - 2] > T[p - 1]); } cache[i].symbol = symbol;
    }
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; if (p0 > 0) { SA[i + 0] = 0; cache[i + 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)(T[p0 - 2] > T[p0 - 1]) << (SAINT_BIT - 1)); symbol0 = T[p0 - 1]; } cache[i + 0].symbol = symbol0;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; if (p1 > 0) { SA[i + 1] = 0; cache[i + 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)(T[p1 - 2] > T[p1 - 1]) << (SAINT_BIT - 1)); symbol1 = T[p1 - 1]; } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; if (p > 0) { SA[i] = 0; cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)(T[p - 2] > T[p - 1]) << (SAINT_BIT - 1)); symbol = T[p - 1]; } cache[i].symbol = symbol;
    }
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_6k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        libsais_prefetchw(&buckets[cache[i - prefetch_distance - 0].symbol]);
        libsais_prefetchw(&buckets[cache[i - prefetch_distance - 1].symbol]);

        sa_sint_t v0 = cache[i - 0].symbol, p0 = cache[i - 0].index; d += (p0 < 0); sa_sint_t * bucket_v0 = &buckets[v0]; cache[i - 0].symbol = --bucket_v0[0]; cache[i - 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v0[2] != d) << (SAINT_BIT - 1)); bucket_v0[2] = d;
        if (cache[i - 0].symbol >= omp_block_start) { sa_sint_t s = cache[i - 0].symbol, q = (cache[s].index = cache[i - 0].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] > T[q - 1]); }

        sa_sint_t v1 = cache[i - 1].symbol, p1 = cache[i - 1].index; d += (p1 < 0); sa_sint_t * bucket_v1 = &buckets[v1]; cache[i - 1].symbol = --bucket_v1[0]; cache[i - 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)(bucket_v1[2] != d) << (SAINT_BIT - 1)); bucket_v1[2] = d;
        if (cache[i - 1].symbol >= omp_block_start) { sa_sint_t s = cache[i - 1].symbol, q = (cache[s].index = cache[i - 1].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] > T[q - 1]); }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t v = cache[i].symbol, p = cache[i].index; d += (p < 0); sa_sint_t * bucket_v = &buckets[v]; cache[i].symbol = --bucket_v[0]; cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)(bucket_v[2] != d) << (SAINT_BIT - 1)); bucket_v[2] = d;
        if (cache[i].symbol >= omp_block_start) { sa_sint_t s = cache[i].symbol, q = (cache[s].index = cache[i].index) & SAINT_MAX; cache[s].symbol = BUCKETS_INDEX4(T[q - 1], T[q - 2] > T[q - 1]); }
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_4k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT induction_bucket = &buckets[3 * (fast_sint_t)k];
    sa_sint_t * RESTRICT distinct_names   = &buckets[0 * (fast_sint_t)k];

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i - prefetch_distance - 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 >> 1 : 0]; libsais_prefetchw(Is0); const sa_sint_t * Ds0 = &distinct_names[s0 > 0 ? s0 : 0]; libsais_prefetchw(Ds0); 
        sa_sint_t s1 = cache[i - prefetch_distance - 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 >> 1 : 0]; libsais_prefetchw(Is1); const sa_sint_t * Ds1 = &distinct_names[s1 > 0 ? s1 : 0]; libsais_prefetchw(Ds1);

        sa_sint_t v0 = cache[i - 0].symbol;
        if (v0 >= 0)
        {
            sa_sint_t p0 = cache[i - 0].index; d += (p0 >> (SUFFIX_GROUP_BIT - 1)); cache[i - 0].symbol = --induction_bucket[v0 >> 1]; cache[i - 0].index = (p0 - 1) | (sa_sint_t)((sa_uint_t)v0 << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v0] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v0] = d;
            if (cache[i - 0].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 0].symbol, np = cache[i - 0].index; if (np > 0) { cache[i - 0].index = 0; cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] > T[np - 1]); } }
        }

        sa_sint_t v1 = cache[i - 1].symbol;
        if (v1 >= 0)
        {
            sa_sint_t p1 = cache[i - 1].index; d += (p1 >> (SUFFIX_GROUP_BIT - 1)); cache[i - 1].symbol = --induction_bucket[v1 >> 1]; cache[i - 1].index = (p1 - 1) | (sa_sint_t)((sa_uint_t)v1 << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v1] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v1] = d;
            if (cache[i - 1].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 1].symbol, np = cache[i - 1].index; if (np > 0) { cache[i - 1].index = 0; cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] > T[np - 1]); } }
        }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            sa_sint_t p = cache[i].index; d += (p >> (SUFFIX_GROUP_BIT - 1)); cache[i].symbol = --induction_bucket[v >> 1]; cache[i].index = (p - 1) | (sa_sint_t)((sa_uint_t)v << (SAINT_BIT - 1)) | ((sa_sint_t)(distinct_names[v] != d) << (SUFFIX_GROUP_BIT - 1)); distinct_names[v] = d;
            if (cache[i].symbol >= omp_block_start) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; if (np > 0) { cache[i].index = 0; cache[ni].index = np; np &= ~SUFFIX_GROUP_MARKER; cache[ni].symbol = BUCKETS_INDEX2(T[np - 1], T[np - 2] > T[np - 1]); } }
        }
    }

    return d;
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i - prefetch_distance - 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 : 0]; libsais_prefetchw(Is0);
        sa_sint_t s1 = cache[i - prefetch_distance - 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 : 0]; libsais_prefetchw(Is1);

        sa_sint_t v0 = cache[i - 0].symbol;
        if (v0 >= 0)
        {
            cache[i - 0].symbol = --induction_bucket[v0];
            if (cache[i - 0].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 0].symbol, np = cache[i - 0].index; if (np > 0) { cache[i - 0].index = 0; cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] > T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; } }
        }

        sa_sint_t v1 = cache[i - 1].symbol;
        if (v1 >= 0)
        {
            cache[i - 1].symbol = --induction_bucket[v1];
            if (cache[i - 1].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 1].symbol, np = cache[i - 1].index; if (np > 0) { cache[i - 1].index = 0; cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] > T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; }}
        }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            cache[i].symbol = --induction_bucket[v];
            if (cache[i].symbol >= omp_block_start) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; if (np > 0) { cache[i].index = 0; cache[ni].index = (np - 1) | (sa_sint_t)((sa_uint_t)(T[np - 2] > T[np - 1]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np - 1]; } }
        }
    }
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_6k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_right_to_left_32s_6k(T, SA, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_right_to_left_32s_6k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                d = libsais_partial_sorting_scan_right_to_left_32s_6k_block_sort(T, buckets, d, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_4k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            d = libsais_partial_sorting_scan_right_to_left_32s_4k(T, SA, k, buckets, d, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_right_to_left_32s_4k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                d = libsais_partial_sorting_scan_right_to_left_32s_4k_block_sort(T, k, buckets, d, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return d;
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(cache);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_partial_sorting_scan_right_to_left_32s_1k(T, SA, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                libsais_partial_sorting_scan_right_to_left_32s_1k_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                libsais_partial_sorting_scan_right_to_left_32s_1k_block_sort(T, buckets, cache - block_start, block_start, block_size);
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

#endif

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t scan_start    = (fast_sint_t)left_suffixes_count + 1;
    fast_sint_t scan_end      = (fast_sint_t)n - (fast_sint_t)first_lms_suffix;

    if (threads == 1 || (scan_end - scan_start) < 65536)
    {
        d = libsais_partial_sorting_scan_right_to_left_32s_6k(T, SA, buckets, d, scan_start, scan_end - scan_start);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = scan_end - 1; block_start >= scan_start; block_start = block_end)
        {
            block_end = block_start - libsais_get_block_size_32s(n, threads); if (block_end < scan_start) { block_end = scan_start - 1; }

            d = libsais_partial_sorting_scan_right_to_left_32s_6k_block_omp(T, SA, buckets, d, thread_state[0].state.cache, block_end + 1, block_start - block_end, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif

    return d;
}

static sa_sint_t libsais_partial_sorting_scan_right_to_left_32s_4k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t d, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || n < 65536)
    {
        d = libsais_partial_sorting_scan_right_to_left_32s_4k(T, SA, k, buckets, d, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; block_start = block_end)
        {
            block_end = block_start - libsais_get_block_size_32s(n, threads); if (block_end < 0) { block_end = -1; }

            d = libsais_partial_sorting_scan_right_to_left_32s_4k_block_omp(T, SA, k, buckets, d, thread_state[0].state.cache, block_end + 1, block_start - block_end, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif

    return d;
}

static void libsais_partial_sorting_scan_right_to_left_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || n < 65536)
    {
        libsais_partial_sorting_scan_right_to_left_32s_1k(T, SA, buckets, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; block_start = block_end)
        {
            block_end = block_start - libsais_get_block_size_32s(n, threads); if (block_end < 0) { block_end = -1; }

            libsais_partial_sorting_scan_right_to_left_32s_1k_block_omp(T, SA, buckets, thread_state[0].state.cache, block_end + 1, block_start - block_end, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static fast_sint_t libsais_partial_sorting_gather_lms_suffixes_32s_4k(sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3, l = omp_block_start; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + prefetch_distance]);

        sa_uint_t s0 = (sa_uint_t)SA[i + 0]; SA[l] = (sa_sint_t)((s0 - (sa_uint_t)SUFFIX_GROUP_MARKER) & (sa_uint_t)(~SUFFIX_GROUP_MARKER)); l += ((sa_sint_t)s0 < 0);
        sa_uint_t s1 = (sa_uint_t)SA[i + 1]; SA[l] = (sa_sint_t)((s1 - (sa_uint_t)SUFFIX_GROUP_MARKER) & (sa_uint_t)(~SUFFIX_GROUP_MARKER)); l += ((sa_sint_t)s1 < 0);
        sa_uint_t s2 = (sa_uint_t)SA[i + 2]; SA[l] = (sa_sint_t)((s2 - (sa_uint_t)SUFFIX_GROUP_MARKER) & (sa_uint_t)(~SUFFIX_GROUP_MARKER)); l += ((sa_sint_t)s2 < 0);
        sa_uint_t s3 = (sa_uint_t)SA[i + 3]; SA[l] = (sa_sint_t)((s3 - (sa_uint_t)SUFFIX_GROUP_MARKER) & (sa_uint_t)(~SUFFIX_GROUP_MARKER)); l += ((sa_sint_t)s3 < 0);
    }

    for (j += 3; i < j; i += 1)
    {
        sa_uint_t s = (sa_uint_t)SA[i]; SA[l] = (sa_sint_t)((s - (sa_uint_t)SUFFIX_GROUP_MARKER) & (sa_uint_t)(~SUFFIX_GROUP_MARKER)); l += ((sa_sint_t)s < 0);
    }

    return l;
}

static fast_sint_t libsais_partial_sorting_gather_lms_suffixes_32s_1k(sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3, l = omp_block_start; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + prefetch_distance]);

        sa_sint_t s0 = SA[i + 0]; SA[l] = s0 & SAINT_MAX; l += (s0 < 0);
        sa_sint_t s1 = SA[i + 1]; SA[l] = s1 & SAINT_MAX; l += (s1 < 0);
        sa_sint_t s2 = SA[i + 2]; SA[l] = s2 & SAINT_MAX; l += (s2 < 0);
        sa_sint_t s3 = SA[i + 3]; SA[l] = s3 & SAINT_MAX; l += (s3 < 0);
    }

    for (j += 3; i < j; i += 1)
    {
        sa_sint_t s = SA[i]; SA[l] = s & SAINT_MAX; l += (s < 0);
    }

    return l;
}

static void libsais_partial_sorting_gather_lms_suffixes_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_partial_sorting_gather_lms_suffixes_32s_4k(SA, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.position = omp_block_start;
                thread_state[omp_thread_num].state.count = libsais_partial_sorting_gather_lms_suffixes_32s_4k(SA, omp_block_start, omp_block_size) - omp_block_start;
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, position = 0;
                for (t = 0; t < omp_num_threads; ++t)
                { 
                    if (t > 0 && thread_state[t].state.count > 0)
                    {
                        memmove(&SA[position], &SA[thread_state[t].state.position], (size_t)thread_state[t].state.count * sizeof(sa_sint_t));
                    }

                    position += thread_state[t].state.count;
                }
            }
        }
#endif
    }
}

static void libsais_partial_sorting_gather_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_partial_sorting_gather_lms_suffixes_32s_1k(SA, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.position = omp_block_start;
                thread_state[omp_thread_num].state.count = libsais_partial_sorting_gather_lms_suffixes_32s_1k(SA, omp_block_start, omp_block_size) - omp_block_start;
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, position = 0;
                for (t = 0; t < omp_num_threads; ++t)
                { 
                    if (t > 0 && thread_state[t].state.count > 0)
                    {
                        memmove(&SA[position], &SA[thread_state[t].state.position], (size_t)thread_state[t].state.count * sizeof(sa_sint_t));
                    }

                    position += thread_state[t].state.count;
                }
            }
        }
#endif
    }
}

static void libsais_induce_partial_order_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t flags, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    memset(&buckets[2 * ALPHABET_SIZE], 0, (size_t)2 * ALPHABET_SIZE * sizeof(sa_sint_t));

    if (flags & LIBSAIS_FLAGS_GSA)
    {
        buckets[4 * ALPHABET_SIZE + BUCKETS_INDEX2(0, 1)] = buckets[4 * ALPHABET_SIZE + BUCKETS_INDEX2(1, 1)] - 1;
        libsais_flip_suffix_markers_omp(SA, buckets[4 * ALPHABET_SIZE + BUCKETS_INDEX2(0, 1)], threads);
    }

    sa_sint_t d = libsais_partial_sorting_scan_left_to_right_8u_omp(T, SA, n, k, buckets, left_suffixes_count, 0, threads, thread_state);
    libsais_partial_sorting_shift_markers_8u_omp(SA, n, buckets, threads);

    if (flags & LIBSAIS_FLAGS_GSA)
    {
        libsais_partial_gsa_scan_right_to_left_8u_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count, d, threads, thread_state);

        if (T[first_lms_suffix] == 0)
        {
            memmove(&SA[1], &SA[0], (size_t)(buckets[BUCKETS_INDEX2(1, 1)] - 1) * sizeof(sa_sint_t));
            SA[0] = first_lms_suffix | SAINT_MIN;
        }

        buckets[BUCKETS_INDEX2(0, 1)] = 0;
    }
    else
    {
        libsais_partial_sorting_scan_right_to_left_8u_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count, d, threads, thread_state);
    }
}

static void libsais_induce_partial_order_32s_6k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t first_lms_suffix, sa_sint_t left_suffixes_count, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t d = libsais_partial_sorting_scan_left_to_right_32s_6k_omp(T, SA, n, buckets, left_suffixes_count, 0, threads, thread_state);
    libsais_partial_sorting_shift_markers_32s_6k_omp(SA, k, buckets, threads);
    libsais_partial_sorting_shift_buckets_32s_6k(k, buckets);
    libsais_partial_sorting_scan_right_to_left_32s_6k_omp(T, SA, n, buckets, first_lms_suffix, left_suffixes_count, d, threads, thread_state);
}

static void libsais_induce_partial_order_32s_4k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    memset(buckets, 0, 2 * (size_t)k * sizeof(sa_sint_t));

    sa_sint_t d = libsais_partial_sorting_scan_left_to_right_32s_4k_omp(T, SA, n, k, buckets, 0, threads, thread_state);
    libsais_partial_sorting_shift_markers_32s_4k(SA, n);
    libsais_partial_sorting_scan_right_to_left_32s_4k_omp(T, SA, n, k, buckets, d, threads, thread_state);
    libsais_partial_sorting_gather_lms_suffixes_32s_4k_omp(SA, n, threads, thread_state);
}

static void libsais_induce_partial_order_32s_2k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_partial_sorting_scan_left_to_right_32s_1k_omp(T, SA, n, &buckets[1 * (fast_sint_t)k], threads, thread_state);
    libsais_partial_sorting_scan_right_to_left_32s_1k_omp(T, SA, n, &buckets[0 * (fast_sint_t)k], threads, thread_state);
    libsais_partial_sorting_gather_lms_suffixes_32s_1k_omp(SA, n, threads, thread_state);
}

static void libsais_induce_partial_order_32s_1k_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_start_32s_1k(k, buckets);
    libsais_partial_sorting_scan_left_to_right_32s_1k_omp(T, SA, n, buckets, threads, thread_state);

    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_end_32s_1k(k, buckets);
    libsais_partial_sorting_scan_right_to_left_32s_1k_omp(T, SA, n, buckets, threads, thread_state);

    libsais_partial_sorting_gather_lms_suffixes_32s_1k_omp(SA, n, threads, thread_state);
}

static sa_sint_t libsais_renumber_lms_suffixes_8u(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t name, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & SAINT_MAX) >> 1]);

        sa_sint_t p0 = SA[i + 0]; SAm[(p0 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p0 < 0;
        sa_sint_t p1 = SA[i + 1]; SAm[(p1 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p1 < 0;
        sa_sint_t p2 = SA[i + 2]; SAm[(p2 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p2 < 0;
        sa_sint_t p3 = SA[i + 3]; SAm[(p3 & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SAm[(p & SAINT_MAX) >> 1] = name | SAINT_MIN; name += p < 0;
    }

    return name;
}

static fast_sint_t libsais_gather_marked_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    l -= 1;

    fast_sint_t i, j;
    for (i = (fast_sint_t)m + omp_block_start + omp_block_size - 1, j = (fast_sint_t)m + omp_block_start + 3; i >= j; i -= 4)
    {
        libsais_prefetchr(&SA[i - prefetch_distance]);

        sa_sint_t s0 = SA[i - 0]; SA[l] = s0 & SAINT_MAX; l -= s0 < 0;
        sa_sint_t s1 = SA[i - 1]; SA[l] = s1 & SAINT_MAX; l -= s1 < 0;
        sa_sint_t s2 = SA[i - 2]; SA[l] = s2 & SAINT_MAX; l -= s2 < 0;
        sa_sint_t s3 = SA[i - 3]; SA[l] = s3 & SAINT_MAX; l -= s3 < 0;
    }

    for (j -= 3; i >= j; i -= 1)
    {
        sa_sint_t s = SA[i]; SA[l] = s & SAINT_MAX; l -= s < 0;
    }

    l += 1;

    return l;
}

static sa_sint_t libsais_renumber_lms_suffixes_8u_omp(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t name = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;

        if (omp_num_threads == 1)
        {
            name = libsais_renumber_lms_suffixes_8u(SA, m, 0, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_count_negative_marked_suffixes(SA, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, count = 0; for (t = 0; t < omp_thread_num; ++t) { count += thread_state[t].state.count; }

                if (omp_thread_num == omp_num_threads - 1)
                {
                    name = (sa_sint_t)(count + thread_state[omp_thread_num].state.count);
                }

                libsais_renumber_lms_suffixes_8u(SA, m, (sa_sint_t)count, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return name;
}

static void libsais_gather_marked_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 131072)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (((fast_sint_t)n >> 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : ((fast_sint_t)n >> 1) - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_gather_marked_lms_suffixes(SA, m, (fast_sint_t)n + (fast_sint_t)fs, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                if (omp_thread_num < omp_num_threads - 1)
                {
                    thread_state[omp_thread_num].state.position = libsais_gather_marked_lms_suffixes(SA, m, (fast_sint_t)m + omp_block_start + omp_block_size, omp_block_start, omp_block_size);
                    thread_state[omp_thread_num].state.count = (fast_sint_t)m + omp_block_start + omp_block_size - thread_state[omp_thread_num].state.position;
                }
                else
                {
                    thread_state[omp_thread_num].state.position = libsais_gather_marked_lms_suffixes(SA, m, (fast_sint_t)n + (fast_sint_t)fs, omp_block_start, omp_block_size);
                    thread_state[omp_thread_num].state.count = (fast_sint_t)n + (fast_sint_t)fs - thread_state[omp_thread_num].state.position;
                }
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, position = (fast_sint_t)n + (fast_sint_t)fs;
                    
                for (t = omp_num_threads - 1; t >= 0; --t)
                { 
                    position -= thread_state[t].state.count;
                    if (t != omp_num_threads - 1 && thread_state[t].state.count > 0)
                    {
                        memmove(&SA[position], &SA[thread_state[t].state.position], (size_t)thread_state[t].state.count * sizeof(sa_sint_t));
                    }
                }
            }
        }
#endif
    }
}

static fast_sint_t libsais_advance_nonunique_lms_suffixes(const sa_sint_t * RESTRICT SAm, fast_sint_t slots, LIBSAIS_NONUNIQUE_SUFFIX * RESTRICT suffixes, fast_sint_t count, fast_sint_t budget)
{
    fast_sint_t i;
    for (i = 0; i < count; i += 1)
    {
        fast_sint_t c = suffixes[i].cursor, s = c + 1; while (SAm[s] == 0) { s += 1; }

        suffixes[i].name = s < slots ? SAm[s] & SAINT_MAX : -1; suffixes[i].cursor = (sa_sint_t)s; budget -= 1 + ((s - c) >> 4);
    }

    return budget;
}

static fast_sint_t libsais_advance_equal_nonunique_lms_suffixes(const sa_sint_t * RESTRICT SAm, fast_sint_t slots, LIBSAIS_NONUNIQUE_SUFFIX * RESTRICT suffixes, fast_sint_t count, fast_sint_t budget, fast_sint_t limit)
{
    fast_sint_t i = count;
    while (i == count && budget >= limit)
    {
        budget = libsais_advance_nonunique_lms_suffixes(SAm, slots, suffixes, count, budget);
        for (i = 1; i < count && suffixes[i].name == suffixes[0].name; i += 1) {}
    }

    return budget;
}

static fast_sint_t libsais_sort_nonunique_lms_suffixes_pair(sa_sint_t * RESTRICT SA, const sa_sint_t * RESTRICT SAm, fast_sint_t slots, fast_sint_t group_start, fast_sint_t budget, fast_sint_t limit)
{
    sa_sint_t p0 = SA[group_start + 0] & SAINT_MAX, p1 = SA[group_start + 1] & SAINT_MAX, n0, n1;
    fast_sint_t c0 = (fast_sint_t)(p0 >> 1), c1 = (fast_sint_t)(p1 >> 1);

    do
    {
        fast_sint_t s0 = c0 + 1; while (SAm[s0] == 0) { s0 += 1; }
        fast_sint_t s1 = c1 + 1; while (SAm[s1] == 0) { s1 += 1; }

        n0 = s0 < slots ? SAm[s0] & SAINT_MAX : -1; n1 = s1 < slots ? SAm[s1] & SAINT_MAX : -1; budget -= 2 + ((s0 - c0) >> 4) + ((s1 - c1) >> 4);

        c0 = s0; c1 = s1;
    } while (n0 == n1 && budget >= limit);

    if (budget >= limit)
    {
        SA[group_start + 0] = n0 < n1 ? p0 : p1;
        SA[group_start + 1] = n0 < n1 ? p1 : p0;
    }

    return budget;
}

static fast_sint_t libsais_sort_nonunique_lms_suffixes_group(sa_sint_t * RESTRICT SA, const sa_sint_t * RESTRICT SAm, fast_sint_t slots, LIBSAIS_NONUNIQUE_SUFFIX * RESTRICT suffixes, fast_sint_t group_start, fast_sint_t group_size, fast_sint_t budget, fast_sint_t limit)
{
    if (group_size == 2) { return libsais_sort_nonunique_lms_suffixes_pair(SA, SAm, slots, group_start, budget, limit); }

    sa_sint_t * RESTRICT stack = (sa_sint_t *)(void *)&suffixes[group_size];

    fast_sint_t i, j, top = 0;
    for (i = 0; i < group_size; i += 1) { sa_sint_t p = SA[group_start + i] & SAINT_MAX; suffixes[i].cursor = p >> 1; suffixes[i].index = p; }

    budget = libsais_advance_equal_nonunique_lms_suffixes(SAm, slots, suffixes, group_size, budget, limit);

    stack[top++] = 0; stack[top++] = (sa_sint_t)group_size;
    while (top > 0 && budget >= limit)
    {
        fast_sint_t hi = stack[--top], lo = stack[--top];

        if (hi - lo == 2)
        {
            while (suffixes[lo].name == suffixes[lo + 1].name && budget >= limit) { budget = libsais_advance_nonunique_lms_suffixes(SAm, slots, &suffixes[lo], 2, budget); }
            if (suffixes[lo].name > suffixes[lo + 1].name) { LIBSAIS_NONUNIQUE_SUFFIX s = suffixes[lo]; suffixes[lo] = suffixes[lo + 1]; suffixes[lo + 1] = s; }
        }
        else if (hi - lo <= LIBSAIS_NONUNIQUE_INSERTION_SORT_SIZE)
        {
            for (i = lo + 1; i < hi; i += 1)
            {
                for (j = i; j > lo && suffixes[j - 1].name > suffixes[j].name; j -= 1) { LIBSAIS_NONUNIQUE_SUFFIX s = suffixes[j - 1]; suffixes[j - 1] = suffixes[j]; suffixes[j] = s; }
            }

            for (i = lo; i < hi; i = j)
            {
                for (j = i + 1; j < hi && suffixes[j].name == suffixes[i].name; j += 1) {}
                if (j - i > 1) { budget = libsais_advance_equal_nonunique_lms_suffixes(SAm, slots, &suffixes[i], j - i, budget, limit); stack[top++] = (sa_sint_t)i; stack[top++] = (sa_sint_t)j; }
            }
        }
        else
        {
            sa_sint_t a = suffixes[lo].name, b = suffixes[lo + ((hi - lo) >> 1)].name, c = suffixes[hi - 1].name;
            sa_sint_t x = a < b ? a : b, y = a < b ? b : a, v = c < x ? x : (c > y ? y : c);

            fast_sint_t lt = lo, gt = hi;
            for (i = lo; i < gt; )
            {
                LIBSAIS_NONUNIQUE_SUFFIX s = suffixes[i];
                if (s.name < v) { suffixes[i] = suffixes[lt]; suffixes[lt] = s; lt += 1; i += 1; }
                else if (s.name > v) { gt -= 1; suffixes[i] = suffixes[gt]; suffixes[gt] = s; }
                else { i += 1; }
            }

            if (lt - lo > 1) { budget -= lt - lo; stack[top++] = (sa_sint_t)lo; stack[top++] = (sa_sint_t)lt; }
            if (hi - gt > 1) { budget -= hi - gt; stack[top++] = (sa_sint_t)gt; stack[top++] = (sa_sint_t)hi; }
            if (gt - lt > 1) { budget = libsais_advance_equal_nonunique_lms_suffixes(SAm, slots, &suffixes[lt], gt - lt, budget, limit); stack[top++] = (sa_sint_t)lt; stack[top++] = (sa_sint_t)gt; }
        }
    }

    if (budget >= limit)
    {
        for (i = 0; i < group_size; i += 1) { SA[group_start + i] = suffixes[i].index; }
    }

    return budget;
}

static fast_sint_t libsais_sort_nonunique_lms_suffixes(sa_sint_t * RESTRICT SA, const sa_sint_t * RESTRICT SAm, fast_sint_t slots, sa_sint_t * RESTRICT records, LIBSAIS_NONUNIQUE_SUFFIX * RESTRICT suffixes, fast_sint_t capacity, fast_sint_t budget, fast_sint_t * RESTRICT deferred, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    fast_sint_t g, j, count = 0;
    for (g = omp_block_start, j = omp_block_start + omp_block_size; g < j && budget >= 0; g += 1)
    {
        if (g + LIBSAIS_NONUNIQUE_PREFETCH_GROUPS < j)
        {
            fast_sint_t p = records[2 * (g + LIBSAIS_NONUNIQUE_PREFETCH_GROUPS) + 0], q = p + records[2 * (g + LIBSAIS_NONUNIQUE_PREFETCH_GROUPS) + 1];
            for (q = q < p + LIBSAIS_NONUNIQUE_PREFETCH_MEMBERS ? q : p + LIBSAIS_NONUNIQUE_PREFETCH_MEMBERS; p < q; p += 1) { libsais_prefetchr(&SAm[((SA[p] & SAINT_MAX) >> 1) + 1]); }
        }

        fast_sint_t group_start = records[2 * g + 0], group_size = records[2 * g + 1];
        fast_sint_t limit = group_size < budget / LIBSAIS_NONUNIQUE_DEFERRED_STEPS ? budget - LIBSAIS_NONUNIQUE_DEFERRED_STEPS * group_size : 0;

        if (4 * group_size <= capacity) { budget = libsais_sort_nonunique_lms_suffixes_group(SA, SAm, slots, suffixes, group_start, group_size, budget, limit); }
        if (4 * group_size > capacity || (budget >= 0 && budget < limit)) { records[2 * g + 1] = (sa_sint_t)(-group_size); count += 1; }
    }

    *deferred = count;
    return budget;
}

static fast_sint_t libsais_gather_nonunique_lms_suffixes_8u(const sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT records, fast_sint_t capacity, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i = omp_block_start, j = omp_block_start + omp_block_size, l = 0;
    if (i > 0) { while (i < j && SA[i - 1] >= 0) { i += 1; } }

    fast_sint_t group_start = i;
    for (j -= 7; i < j && l + 10 <= capacity; i += 8)
    {
        libsais_prefetchr(&SA[i + prefetch_distance]);

        if (group_start < i || (SA[i + 0] & SA[i + 1] & SA[i + 2] & SA[i + 3] & SA[i + 4] & SA[i + 5] & SA[i + 6] & SA[i + 7]) >= 0)
        {
            sa_sint_t s0 = SA[i + 0]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 1); l += (fast_sint_t)((s0 < 0) & (i + 0 > group_start)) << 1; group_start = s0 < 0 ? i + 1 : group_start;
            sa_sint_t s1 = SA[i + 1]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 2); l += (fast_sint_t)((s1 < 0) & (i + 1 > group_start)) << 1; group_start = s1 < 0 ? i + 2 : group_start;
            sa_sint_t s2 = SA[i + 2]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 3); l += (fast_sint_t)((s2 < 0) & (i + 2 > group_start)) << 1; group_start = s2 < 0 ? i + 3 : group_start;
            sa_sint_t s3 = SA[i + 3]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 4); l += (fast_sint_t)((s3 < 0) & (i + 3 > group_start)) << 1; group_start = s3 < 0 ? i + 4 : group_start;
            sa_sint_t s4 = SA[i + 4]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 5); l += (fast_sint_t)((s4 < 0) & (i + 4 > group_start)) << 1; group_start = s4 < 0 ? i + 5 : group_start;
            sa_sint_t s5 = SA[i + 5]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 6); l += (fast_sint_t)((s5 < 0) & (i + 5 > group_start)) << 1; group_start = s5 < 0 ? i + 6 : group_start;
            sa_sint_t s6 = SA[i + 6]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 7); l += (fast_sint_t)((s6 < 0) & (i + 6 > group_start)) << 1; group_start = s6 < 0 ? i + 7 : group_start;
            sa_sint_t s7 = SA[i + 7]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 8); l += (fast_sint_t)((s7 < 0) & (i + 7 > group_start)) << 1; group_start = s7 < 0 ? i + 8 : group_start;
        }
        else
        {
            group_start = i + 8;
        }
    }

    for (j += 7; i < j && l + 2 <= capacity; i += 1)
    {
        sa_sint_t s = SA[i]; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 1); l += (fast_sint_t)((s < 0) & (i > group_start)) << 1; group_start = s < 0 ? i + 1 : group_start;
    }

    if (i < j) { return -1; }

    if (group_start < i)
    {
        while (SA[i] >= 0) { i += 1; }

        if (l + 2 > capacity) { return -1; }
        records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start + 1); l += 2;
    }

    return l >> 1;
}

static fast_sint_t libsais_gather_nonunique_lms_suffixes_32s(const sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t * RESTRICT records, fast_sint_t capacity, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i = omp_block_start, j = omp_block_start + omp_block_size, l = 0;
    sa_sint_t p = i > 0 ? SAm[SA[i - 1] >> 1] & SAINT_MAX : -1;
    if (i > 0) { while (i < j && (SAm[SA[i] >> 1] & SAINT_MAX) == p) { i += 1; } }

    fast_sint_t group_start = i;
    for (j -= prefetch_distance + 3; i < j && l + 6 <= capacity; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&SAm[SA[i + prefetch_distance + 0] >> 1]);
        libsais_prefetchr(&SAm[SA[i + prefetch_distance + 1] >> 1]);
        libsais_prefetchr(&SAm[SA[i + prefetch_distance + 2] >> 1]);
        libsais_prefetchr(&SAm[SA[i + prefetch_distance + 3] >> 1]);

        sa_sint_t q0 = SAm[SA[i + 0] >> 1] & SAINT_MAX; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i + 0 - group_start); l += (fast_sint_t)((q0 != p) & (i + 0 - group_start > 1)) << 1; group_start = q0 != p ? i + 0 : group_start; p = q0;
        sa_sint_t q1 = SAm[SA[i + 1] >> 1] & SAINT_MAX; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i + 1 - group_start); l += (fast_sint_t)((q1 != p) & (i + 1 - group_start > 1)) << 1; group_start = q1 != p ? i + 1 : group_start; p = q1;
        sa_sint_t q2 = SAm[SA[i + 2] >> 1] & SAINT_MAX; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i + 2 - group_start); l += (fast_sint_t)((q2 != p) & (i + 2 - group_start > 1)) << 1; group_start = q2 != p ? i + 2 : group_start; p = q2;
        sa_sint_t q3 = SAm[SA[i + 3] >> 1] & SAINT_MAX; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i + 3 - group_start); l += (fast_sint_t)((q3 != p) & (i + 3 - group_start > 1)) << 1; group_start = q3 != p ? i + 3 : group_start; p = q3;
    }

    for (j += prefetch_distance + 3; i < j && l + 2 <= capacity; i += 1)
    {
        sa_sint_t q = SAm[SA[i] >> 1] & SAINT_MAX; records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start); l += (fast_sint_t)((q != p) & (i - group_start > 1)) << 1; group_start = q != p ? i : group_start; p = q;
    }

    if (i < j) { return -1; }

    if (group_start < j)
    {
        while (i < m && (SAm[SA[i] >> 1] & SAINT_MAX) == p) { i += 1; }

        if (i - group_start > 1)
        {
            if (l + 2 > capacity) { return -1; }
            records[l + 0] = (sa_sint_t)group_start; records[l + 1] = (sa_sint_t)(i - group_start); l += 2;
        }
    }

    return l >> 1;
}

static fast_sint_t libsais_gather_nonunique_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t * RESTRICT records, fast_sint_t capacity, sa_sint_t marked, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t groups = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;

        if (omp_num_threads == 1)
        {
            groups = marked
                ? libsais_gather_nonunique_lms_suffixes_8u(SA, records, capacity, omp_block_start, omp_block_size)
                : libsais_gather_nonunique_lms_suffixes_32s(SA, m, records, capacity, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t omp_capacity = (capacity / omp_num_threads) & (-2);

            thread_state[omp_thread_num].state.count = marked
                ? libsais_gather_nonunique_lms_suffixes_8u(SA, &records[omp_thread_num * omp_capacity], omp_capacity, omp_block_start, omp_block_size)
                : libsais_gather_nonunique_lms_suffixes_32s(SA, m, &records[omp_thread_num * omp_capacity], omp_capacity, omp_block_start, omp_block_size);

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = 0; t < omp_num_threads && groups >= 0; ++t)
                {
                    fast_sint_t count = thread_state[t].state.count;

                    if (t > 0 && count > 0) { memmove(&records[2 * groups], &records[t * omp_capacity], (size_t)(2 * count) * sizeof(sa_sint_t)); }
                    groups = count >= 0 ? groups + count : -1;
                }
            }
        }
#endif
    }

    return groups;
}

static void libsais_rank_nonunique_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t m, const sa_sint_t * RESTRICT records, fast_sint_t groups, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t lo = 0, hi = groups;
    while (lo < hi)
    {
        fast_sint_t mid = lo + ((hi - lo) >> 1), size = records[2 * mid + 1] > 0 ? records[2 * mid + 1] : -records[2 * mid + 1];
        if ((fast_sint_t)records[2 * mid + 0] + size <= omp_block_start) { lo = mid + 1; } else { hi = mid; }
    }

    fast_sint_t g, i = omp_block_start, j = omp_block_start + omp_block_size;
    for (g = lo; g <= groups && i < j; g += 1)
    {
        if (g == groups || records[2 * g + 1] < 0)
        {
            fast_sint_t group_start = g < groups ? (fast_sint_t)records[2 * g + 0] : j, group_end = g < groups ? group_start - (fast_sint_t)records[2 * g + 1] : j;

            for (; i < (group_start < j ? group_start : j); i += 1)
            {
                if (i + prefetch_distance < j) { libsais_prefetchw(&SAm[(SA[i + prefetch_distance] & SAINT_MAX) >> 1]); }
                SAm[(SA[i] & SAINT_MAX) >> 1] = (sa_sint_t)(i + 1) | SAINT_MIN;
            }

            for (; i < (group_end < j ? group_end : j); i += 1) { SAm[(SA[i] & SAINT_MAX) >> 1] = (sa_sint_t)(group_start + 1) | SAINT_MIN; }
        }
    }
}

static void libsais_rank_nonunique_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t m, const sa_sint_t * RESTRICT records, fast_sint_t groups, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = m;
#endif
        libsais_rank_nonunique_lms_suffixes(SA, m, records, groups, omp_block_start, omp_block_size);
    }
}

static void libsais_rename_nonunique_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t m, const sa_sint_t * RESTRICT records, fast_sint_t groups, sa_sint_t marked)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, g = 0, ties = 0;
    for (i = 0; i < m; i += 1)
    {
        if (i + prefetch_distance < m) { libsais_prefetchw(&SAm[(SA[i + prefetch_distance] & SAINT_MAX) >> 1]); }

        fast_sint_t group_start = g < groups ? (fast_sint_t)records[2 * g + 0] : m, group_size = g < groups ? (fast_sint_t)records[2 * g + 1] : 1;
        if (i < group_start)
        {
            sa_sint_t name = (sa_sint_t)(i - ties); SAm[(SA[i] & SAINT_MAX) >> 1] = marked ? name | SAINT_MIN : (name + 1) | SAINT_MIN;
        }
        else
        {
            sa_sint_t name = (sa_sint_t)(group_start - ties); SAm[(SA[i] & SAINT_MAX) >> 1] = marked ? name | SAINT_MIN : name + 1;
            if (i == group_start + group_size - 1) { ties += group_size - 1; g += 1; }
        }
    }
}

static sa_sint_t libsais_sample_nonunique_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t marked, fast_sint_t limit)
{
    fast_sint_t slots = (fast_sint_t)n >> 1, samples = ((fast_sint_t)m >> 7) < LIBSAIS_NONUNIQUE_SAMPLES ? ((fast_sint_t)m >> 7) : LIBSAIS_NONUNIQUE_SAMPLES;
    fast_sint_t stride = samples > 0 ? ((fast_sint_t)m - 1) / samples : 0, steps = 0, scanned = 0, i;

    sa_sint_t * RESTRICT SAm = &SA[m]; SAm[slots] = SAINT_MIN;

    for (i = 0; i < samples * stride && scanned <= (fast_sint_t)m && 8 * steps <= limit * samples; i += stride)
    {
        sa_sint_t p0 = SA[i] & SAINT_MAX, p1 = SA[i + 1] & SAINT_MAX;

        if (marked ? SA[i] >= 0 : (SAm[p0 >> 1] & SAINT_MAX) == (SAm[p1 >> 1] & SAINT_MAX))
        {
            fast_sint_t c0 = (fast_sint_t)(p0 >> 1), c1 = (fast_sint_t)(p1 >> 1), t = 0; sa_sint_t n0, n1;

            do
            {
                fast_sint_t s0 = c0 + 1; while (SAm[s0] == 0) { s0 += 1; }
                fast_sint_t s1 = c1 + 1; while (SAm[s1] == 0) { s1 += 1; }

                n0 = s0 < slots ? SAm[s0] & SAINT_MAX : -1; n1 = s1 < slots ? SAm[s1] & SAINT_MAX : -1;
                scanned += (s0 - c0) + (s1 - c1); c0 = s0; c1 = s1; t += 1;
            } while (n0 == n1 && t < LIBSAIS_NONUNIQUE_DEFERRED_STEPS);

            steps += t;
        }
    }

    return scanned <= (fast_sint_t)m && 8 * steps <= limit * samples;
}

static sa_sint_t libsais_sort_nonunique_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t name, sa_sint_t * RESTRICT records, fast_sint_t groups, sa_sint_t marked, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t slots = (fast_sint_t)n >> 1, budget = LIBSAIS_NONUNIQUE_BUDGET_STEPS * (fast_sint_t)m, deferred = 0;
    fast_sint_t capacity = (fast_sint_t)n - (fast_sint_t)m - slots - 1 - 2 * groups;

    sa_sint_t * RESTRICT SAm = &SA[m]; SAm[slots] = SAINT_MIN;
    sa_sint_t * RESTRICT scratch = &records[2 * groups];

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && groups >= 1024)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (groups / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : groups - omp_block_start;

        if (omp_num_threads == 1)
        {
            budget = libsais_sort_nonunique_lms_suffixes(SA, SAm, slots, records, (LIBSAIS_NONUNIQUE_SUFFIX *)(void *)scratch, capacity, budget, &deferred, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            fast_sint_t omp_capacity = capacity / omp_num_threads;

            {
                fast_sint_t count = 0;
                thread_state[omp_thread_num].state.position = libsais_sort_nonunique_lms_suffixes(SA, SAm, slots, records, (LIBSAIS_NONUNIQUE_SUFFIX *)(void *)&scratch[omp_thread_num * omp_capacity], omp_capacity, budget, &count, omp_block_start, omp_block_size);
                thread_state[omp_thread_num].state.count = count;
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, used = 0, failed = 0;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    fast_sint_t left = thread_state[t].state.position; failed |= left < 0;
                    if (left >= 0) { used = budget - left < budget - used ? used + (budget - left) : budget; }
                    deferred += thread_state[t].state.count;
                }

                budget = failed || (deferred > 0 && used == budget) ? -1 : budget - used;
            }
        }
#endif
    }

    if (budget < 0) { return name; }
    if (deferred == 0) { return m; }

    libsais_rank_nonunique_lms_suffixes_omp(SA, m, records, groups, threads);

    fast_sint_t g;
    for (g = 0; g < groups && budget >= 0; g += 1)
    {
        if (records[2 * g + 1] < 0)
        {
            fast_sint_t group_start = records[2 * g + 0], group_size = -(fast_sint_t)records[2 * g + 1];
            records[2 * g + 1] = (sa_sint_t)group_size;

            budget = 4 * group_size <= capacity ? libsais_sort_nonunique_lms_suffixes_group(SA, SAm, slots, (LIBSAIS_NONUNIQUE_SUFFIX *)(void *)scratch, group_start, group_size, budget, 0) : -1;

            if (budget >= 0)
            {
                fast_sint_t i; for (i = 0; i < group_size; i += 1) { SAm[(SA[group_start + i] & SAINT_MAX) >> 1] = (sa_sint_t)(group_start + i + 1) | SAINT_MIN; }
            }
        }
    }

    for (; g < groups; g += 1) { if (records[2 * g + 1] < 0) { records[2 * g + 1] = -records[2 * g + 1]; } }

    if (budget >= 0) { return m; }

    libsais_rename_nonunique_lms_suffixes(SA, m, records, groups, marked);
    return name;
}

static sa_sint_t libsais_resolve_nonunique_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t name, sa_sint_t marked, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fast_sint_t slots = (fast_sint_t)n >> 1, ties = (fast_sint_t)m - (fast_sint_t)name, limit = threads > 1 ? ((fast_sint_t)m < 262144 ? 6 : 5) : 4;

    if (ties > 0 && (fast_sint_t)m + slots + 1 <= (fast_sint_t)n && (ties <= ((fast_sint_t)m >> 3) + ((fast_sint_t)m >> 4) || (ties <= ((fast_sint_t)m >> 3) * limit && libsais_sample_nonunique_lms_suffixes(SA, n, m, marked, limit))))
    {
        sa_sint_t * RESTRICT records = &SA[(fast_sint_t)m + slots + 1];

        fast_sint_t groups = libsais_gather_nonunique_lms_suffixes_omp(SA, m, records, (fast_sint_t)n - (fast_sint_t)m - slots - 1, marked, threads, thread_state);
        if (groups >= 0) { name = libsais_sort_nonunique_lms_suffixes_omp(SA, n, m, name, records, groups, marked, threads, thread_state); }
    }

    return name;
}

static sa_sint_t libsais_renumber_and_gather_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    memset(&SA[m], 0, ((size_t)n >> 1) * sizeof(sa_sint_t));

    sa_sint_t name = libsais_renumber_lms_suffixes_8u_omp(SA, m, threads, thread_state);
    name = libsais_resolve_nonunique_lms_suffixes_omp(SA, n, m, name, 1, threads, thread_state);
    if (name < m)
    {
        libsais_gather_marked_lms_suffixes_omp(SA, n, m, fs, threads, thread_state);
    }
    else
    {
        fast_sint_t i; for (i = 0; i < m; i += 1) { SA[i] &= SAINT_MAX; }
    }

    return name;
}

static sa_sint_t libsais_renumber_distinct_lms_suffixes_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t name, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j; sa_sint_t p0, p1, p2, p3 = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 0] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 1] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 2] & SAINT_MAX) >> 1]);
        libsais_prefetchw(&SAm[(SA[i + prefetch_distance + 3] & SAINT_MAX) >> 1]);

        p0 = SA[i + 0]; SAm[(SA[i + 0] = p0 & SAINT_MAX) >> 1] = name | (p0 & p3 & SAINT_MIN); name += p0 < 0;
        p1 = SA[i + 1]; SAm[(SA[i + 1] = p1 & SAINT_MAX) >> 1] = name | (p1 & p0 & SAINT_MIN); name += p1 < 0;
        p2 = SA[i + 2]; SAm[(SA[i + 2] = p2 & SAINT_MAX) >> 1] = name | (p2 & p1 & SAINT_MIN); name += p2 < 0;
        p3 = SA[i + 3]; SAm[(SA[i + 3] = p3 & SAINT_MAX) >> 1] = name | (p3 & p2 & SAINT_MIN); name += p3 < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        p2 = p3; p3 = SA[i]; SAm[(SA[i] = p3 & SAINT_MAX) >> 1] = name | (p3 & p2 & SAINT_MIN); name += p3 < 0;
    }

    return name;
}

static void libsais_mark_distinct_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j; sa_sint_t p0, p1, p2, p3, mask = SAINT_MAX;
    for (i = (fast_sint_t)m + omp_block_start, j = (fast_sint_t)m + omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais_prefetchw(&SA[i + prefetch_distance]);

        p0 = SA[i + 0]; SA[i + 0] = p0 & mask; mask = (p0 != 0) ? (p0 | SAINT_MAX) : mask;
        p1 = SA[i + 1]; SA[i + 1] = p1 & mask; mask = (p1 != 0) ? (p1 | SAINT_MAX) : mask;
        p2 = SA[i + 2]; SA[i + 2] = p2 & mask; mask = (p2 != 0) ? (p2 | SAINT_MAX) : mask;
        p3 = SA[i + 3]; SA[i + 3] = p3 & mask; mask = (p3 != 0) ? (p3 | SAINT_MAX) : mask;
    }

    for (j += 3; i < j; i += 1)
    {
        p3 = SA[i]; SA[i] = p3 & mask; mask = (p3 != 0) ? (p3 | SAINT_MAX) : mask;
    }
}

static void libsais_clamp_lms_suffixes_length_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais_prefetchw(&SAm[i + prefetch_distance]);

        SAm[i + 0] = (SAm[i + 0] < 0 ? SAm[i + 0] : 0) & SAINT_MAX;
        SAm[i + 1] = (SAm[i + 1] < 0 ? SAm[i + 1] : 0) & SAINT_MAX;
        SAm[i + 2] = (SAm[i + 2] < 0 ? SAm[i + 2] : 0) & SAINT_MAX;
        SAm[i + 3] = (SAm[i + 3] < 0 ? SAm[i + 3] : 0) & SAINT_MAX;
    }

    for (j += 3; i < j; i += 1)
    {
        SAm[i] = (SAm[i] < 0 ? SAm[i] : 0) & SAINT_MAX;
    }
}

static sa_sint_t libsais_renumber_distinct_lms_suffixes_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t name = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;

        if (omp_num_threads == 1)
        {
            name = libsais_renumber_distinct_lms_suffixes_32s_4k(SA, m, 1, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_count_negative_marked_suffixes(SA, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, count = 1; for (t = 0; t < omp_thread_num; ++t) { count += thread_state[t].state.count; }

                if (omp_thread_num == omp_num_threads - 1)
                {
                    name = (sa_sint_t)(count + thread_state[omp_thread_num].state.count);
                }

                libsais_renumber_distinct_lms_suffixes_32s_4k(SA, m, (sa_sint_t)count, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return name - 1;
}

static void libsais_mark_distinct_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 131072)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (((fast_sint_t)n >> 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : ((fast_sint_t)n >> 1) - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = (fast_sint_t)n >> 1;
#endif
        libsais_mark_distinct_lms_suffixes_32s(SA, m, omp_block_start, omp_block_size);
    }
}

static void libsais_clamp_lms_suffixes_length_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 131072)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (((fast_sint_t)n >> 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : ((fast_sint_t)n >> 1) - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = (fast_sint_t)n >> 1;
#endif
        libsais_clamp_lms_suffixes_length_32s(SA, m, omp_block_start, omp_block_size);
    }
}

static sa_sint_t libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    memset(&SA[m], 0, ((size_t)n >> 1) * sizeof(sa_sint_t));

    sa_sint_t name = libsais_renumber_distinct_lms_suffixes_32s_4k_omp(SA, m, threads, thread_state);
    name = libsais_resolve_nonunique_lms_suffixes_omp(SA, n, m, name, 0, threads, thread_state);
    if (name < m)
    {
        libsais_mark_distinct_lms_suffixes_32s_omp(SA, n, m, threads);
    }

    return name;
}

static sa_sint_t libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    {
        libsais_gather_lms_suffixes_32s(T, SA, n);

        memset(&SA[m], 0, ((size_t)n - (size_t)m - (size_t)m) * sizeof(sa_sint_t));

        fast_sint_t i, j;
        for (i = (fast_sint_t)n - (fast_sint_t)m, j = (fast_sint_t)n - 1 - prefetch_distance - 3; i < j; i += 4)
        {
            libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 0]) >> 1]);
            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 1]) >> 1]);
            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 2]) >> 1]);
            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 3]) >> 1]);

            SAm[((sa_uint_t)SA[i + 0]) >> 1] = SA[i + 1] - SA[i + 0] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 1]) >> 1] = SA[i + 2] - SA[i + 1] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 2]) >> 1] = SA[i + 3] - SA[i + 2] + 1 + SAINT_MIN;
            SAm[((sa_uint_t)SA[i + 3]) >> 1] = SA[i + 4] - SA[i + 3] + 1 + SAINT_MIN;
        }

        for (j += prefetch_distance + 3; i < j; i += 1)
        {
            SAm[((sa_uint_t)SA[i]) >> 1] = SA[i + 1] - SA[i] + 1 + SAINT_MIN;
        }

        SAm[((sa_uint_t)SA[n - 1]) >> 1] = 1 + SAINT_MIN;
    }

    {
        libsais_clamp_lms_suffixes_length_32s_omp(SA, n, m, threads);
    }

    sa_sint_t name = 1;

    {
        fast_sint_t i, j, p = SA[0], plen = SAm[p >> 1]; sa_sint_t pdiff = SAINT_MIN;
        for (i = 1, j = m - prefetch_distance - 1; i < j; i += 2)
        {
            libsais_prefetchr(&SA[i + 2 * prefetch_distance]);
            
            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 0]) >> 1]); libsais_prefetchr(&T[((sa_uint_t)SA[i + prefetch_distance + 0])]);
            libsais_prefetchw(&SAm[((sa_uint_t)SA[i + prefetch_distance + 1]) >> 1]); libsais_prefetchr(&T[((sa_uint_t)SA[i + prefetch_distance + 1])]);

            fast_sint_t q = SA[i + 0], qlen = SAm[q >> 1]; sa_sint_t qdiff = SAINT_MIN;
            if (plen == qlen) { fast_sint_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < qlen); qdiff = (sa_sint_t)(l - qlen) & SAINT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = SA[i + 1]; plen = SAm[p >> 1]; pdiff = SAINT_MIN;
            if (qlen == plen) { fast_sint_t l = 0; do { if (T[q + l] != T[p + l]) { break; } } while (++l < plen); pdiff = (sa_sint_t)(l - plen) & SAINT_MIN; }
            SAm[q >> 1] = name | (qdiff & pdiff); name += (pdiff < 0);
        }

        for (j += prefetch_distance + 1; i < j; i += 1)
        {
            fast_sint_t q = SA[i], qlen = SAm[q >> 1]; sa_sint_t qdiff = SAINT_MIN;
            if (plen == qlen) { fast_sint_t l = 0; do { if (T[p + l] != T[q + l]) { break; } } while (++l < plen); qdiff = (sa_sint_t)(l - plen) & SAINT_MIN; }
            SAm[p >> 1] = name | (pdiff & qdiff); name += (qdiff < 0);

            p = q; plen = qlen; pdiff = qdiff;
        }

        SAm[p >> 1] = name | pdiff; name++;
    }

    name = libsais_resolve_nonunique_lms_suffixes_omp(SA, n, m, name - 1, 0, threads, thread_state) + 1;

    if (name <= m)
    {
        libsais_mark_distinct_lms_suffixes_32s_omp(SA, n, m, threads);
    }

    return name - 1;
}

static void libsais_reconstruct_lms_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT SAnm = &SA[n - m];

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&SAnm[SA[i + prefetch_distance + 0]]);
        libsais_prefetchr(&SAnm[SA[i + prefetch_distance + 1]]);
        libsais_prefetchr(&SAnm[SA[i + prefetch_distance + 2]]);
        libsais_prefetchr(&SAnm[SA[i + prefetch_distance + 3]]);

        SA[i + 0] = SAnm[SA[i + 0]];
        SA[i + 1] = SAnm[SA[i + 1]];
        SA[i + 2] = SAnm[SA[i + 2]];
        SA[i + 3] = SAnm[SA[i + 3]];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        SA[i] = SAnm[SA[i]];
    }
}

static void libsais_reconstruct_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = m;
#endif

        libsais_reconstruct_lms_suffixes(SA, n, m, omp_block_start, omp_block_size);
    }
}

static void libsais_place_lms_suffixes_interval_8u(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t flags, sa_sint_t * RESTRICT buckets)
{
    if (flags & LIBSAIS_FLAGS_GSA) { buckets[7 * ALPHABET_SIZE]--; }

    {
        const sa_sint_t * RESTRICT bucket_end = &buckets[7 * ALPHABET_SIZE];

        fast_sint_t c, j = n;
        for (c = ALPHABET_SIZE - 2; c >= 0; --c)
        {
            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(c, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(c, 1)];
            if (l > 0)
            {
                fast_sint_t i = bucket_end[c];
                if (j - i > 0)
                {
                    memset(&SA[i], 0, (size_t)(j - i) * sizeof(sa_sint_t));
                }

                memmove(&SA[j = (i - l)], &SA[m -= (sa_sint_t)l], (size_t)l * sizeof(sa_sint_t));
            }
        }

        memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
    }

    if (flags & LIBSAIS_FLAGS_GSA) { buckets[7 * ALPHABET_SIZE]++; }
}

#if defined(__GNUC__) || defined(__clang__)
    #define libsais_ctz64(x) ((fast_sint_t)__builtin_ctzll((unsigned long long)(x)))
#elif defined(_M_X64) || defined(_M_ARM64)
    static fast_sint_t libsais_ctz64(uint64_t x) { unsigned long index; _BitScanForward64(&index, x); return (fast_sint_t)index; }
#else
    static fast_sint_t libsais_ctz64(uint64_t x) { unsigned long index; if ((uint32_t)x != 0) { _BitScanForward(&index, (unsigned long)(uint32_t)x); return (fast_sint_t)index; } _BitScanForward(&index, (unsigned long)(uint32_t)(x >> 32)); return (fast_sint_t)index + 32; }
#endif

#define LIBSAIS_PLACE_CHUNK_SIZE        (16)

static fast_sint_t libsais_place_lms_suffixes_bucket_32s(sa_sint_t * RESTRICT SA, fast_sint_t i, fast_sint_t j, fast_sint_t l, fast_sint_t m)
{
    if ((((fast_uint_t)(j - i)) | ((fast_uint_t)l)) < (fast_uint_t)LIBSAIS_PLACE_CHUNK_SIZE && i - m >= LIBSAIS_PLACE_CHUNK_SIZE && m >= LIBSAIS_PLACE_CHUNK_SIZE)
    {
        sa_sint_t chunk[LIBSAIS_PLACE_CHUNK_SIZE];

        memcpy(chunk, &SA[m - LIBSAIS_PLACE_CHUNK_SIZE], sizeof(chunk));
        memset(&SA[j - LIBSAIS_PLACE_CHUNK_SIZE], 0, sizeof(chunk));
        memcpy(&SA[i - LIBSAIS_PLACE_CHUNK_SIZE], chunk, sizeof(chunk));
    }
    else
    {
        if (j - i > 0)
        {
            memset(&SA[i], 0, (size_t)(j - i) * sizeof(sa_sint_t));
        }

        memmove(&SA[i - l], &SA[m - l], (size_t)l * sizeof(sa_sint_t));
    }

    return i - l;
}

static fast_sint_t libsais_place_lms_suffixes_gap_32s(sa_sint_t * RESTRICT SA, fast_sint_t bucket_end, fast_sint_t l, fast_sint_t i)
{
    if (l - bucket_end < LIBSAIS_PLACE_CHUNK_SIZE && l - LIBSAIS_PLACE_CHUNK_SIZE >= i)
    {
        memset(&SA[l - LIBSAIS_PLACE_CHUNK_SIZE], 0, LIBSAIS_PLACE_CHUNK_SIZE * sizeof(sa_sint_t));
    }
    else
    {
        memset(&SA[bucket_end], 0, (size_t)(l - bucket_end) * sizeof(sa_sint_t));
    }

    return bucket_end;
}

static void libsais_place_lms_suffixes_interval_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[3 * (fast_sint_t)k];

    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; c -= 64)
    {
        fast_sint_t s, group_end = c >= 63 ? c - 63 : 0; uint64_t mask = 0;
        for (s = c; s >= group_end; --s) { mask |= (uint64_t)((fast_sint_t)buckets[BUCKETS_INDEX2(s, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)] > 0) << (c - s); }

        while (mask != 0)
        {
            s = c - libsais_ctz64(mask); mask &= mask - 1;

            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)];
            j = libsais_place_lms_suffixes_bucket_32s(SA, bucket_end[s], j, l, m); m -= (sa_sint_t)l;
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais_place_lms_suffixes_interval_32s_2k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; c -= 64)
    {
        fast_sint_t s, group_end = c >= 63 ? c - 63 : 0; uint64_t mask = 0;
        for (s = c; s >= group_end; --s) { mask |= (uint64_t)((fast_sint_t)buckets[BUCKETS_INDEX2(s, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)] > 0) << (c - s); }

        while (mask != 0)
        {
            s = c - libsais_ctz64(mask); mask &= mask - 1;

            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1) + BUCKETS_INDEX2(1, 0)] - (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)];
            j = libsais_place_lms_suffixes_bucket_32s(SA, buckets[BUCKETS_INDEX2(s, 0)], j, l, m); m -= (sa_sint_t)l;
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais_place_lms_suffixes_interval_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t m, sa_sint_t * RESTRICT buckets)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t c = k - 1; fast_sint_t i, l = buckets[c];
    for (i = (fast_sint_t)m - 1; i >= prefetch_distance + 3; i -= 4)
    {
        libsais_prefetchr(&SA[i - 2 * prefetch_distance]);

        libsais_prefetchr(&T[SA[i - prefetch_distance - 0]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 1]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 2]]);
        libsais_prefetchr(&T[SA[i - prefetch_distance - 3]]);

        sa_sint_t p0 = SA[i - 0]; if (T[p0] != c) { c = T[p0]; l = libsais_place_lms_suffixes_gap_32s(SA, buckets[c], l, i - 0); } SA[--l] = p0;
        sa_sint_t p1 = SA[i - 1]; if (T[p1] != c) { c = T[p1]; l = libsais_place_lms_suffixes_gap_32s(SA, buckets[c], l, i - 1); } SA[--l] = p1;
        sa_sint_t p2 = SA[i - 2]; if (T[p2] != c) { c = T[p2]; l = libsais_place_lms_suffixes_gap_32s(SA, buckets[c], l, i - 2); } SA[--l] = p2;
        sa_sint_t p3 = SA[i - 3]; if (T[p3] != c) { c = T[p3]; l = libsais_place_lms_suffixes_gap_32s(SA, buckets[c], l, i - 3); } SA[--l] = p3;
    }

    for (; i >= 0; i -= 1)
    {
        sa_sint_t p = SA[i]; if (T[p] != c) { c = T[p]; l = libsais_place_lms_suffixes_gap_32s(SA, buckets[c], l, i); } SA[--l] = p;
    }

    memset(&SA[0], 0, (size_t)l * sizeof(sa_sint_t));
}

static void libsais_place_lms_suffixes_histogram_32s_6k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[5 * (fast_sint_t)k];

    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; c -= 64)
    {
        fast_sint_t s, group_end = c >= 63 ? c - 63 : 0; uint64_t mask = 0;
        for (s = c; s >= group_end; --s) { mask |= (uint64_t)((fast_sint_t)buckets[BUCKETS_INDEX4(s, 1)] > 0) << (c - s); }

        while (mask != 0)
        {
            s = c - libsais_ctz64(mask); mask &= mask - 1;

            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX4(s, 1)];
            j = libsais_place_lms_suffixes_bucket_32s(SA, bucket_end[s], j, l, m); m -= (sa_sint_t)l;
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais_place_lms_suffixes_histogram_32s_4k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    const sa_sint_t * RESTRICT bucket_end = &buckets[3 * (fast_sint_t)k];

    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; c -= 64)
    {
        fast_sint_t s, group_end = c >= 63 ? c - 63 : 0; uint64_t mask = 0;
        for (s = c; s >= group_end; --s) { mask |= (uint64_t)((fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)] > 0) << (c - s); }

        while (mask != 0)
        {
            s = c - libsais_ctz64(mask); mask &= mask - 1;

            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)];
            j = libsais_place_lms_suffixes_bucket_32s(SA, bucket_end[s], j, l, m); m -= (sa_sint_t)l;
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais_place_lms_suffixes_histogram_32s_2k(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, const sa_sint_t * RESTRICT buckets)
{
    fast_sint_t c, j = n;
    for (c = (fast_sint_t)k - 2; c >= 0; c -= 64)
    {
        fast_sint_t s, group_end = c >= 63 ? c - 63 : 0; uint64_t mask = 0;
        for (s = c; s >= group_end; --s) { mask |= (uint64_t)((fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)] > 0) << (c - s); }

        while (mask != 0)
        {
            s = c - libsais_ctz64(mask); mask &= mask - 1;

            fast_sint_t l = (fast_sint_t)buckets[BUCKETS_INDEX2(s, 1)];
            j = libsais_place_lms_suffixes_bucket_32s(SA, buckets[BUCKETS_INDEX2(s, 0)], j, l, m); m -= (sa_sint_t)l;
        }
    }

    memset(&SA[0], 0, (size_t)j * sizeof(sa_sint_t));
}

static void libsais_final_bwt_scan_left_to_right_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i + 0]; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; SA[induction_bucket[T[p0]]++] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); } else { SA[i + 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i + 1]; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; SA[induction_bucket[T[p1]]++] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); } else { SA[i + 1] = p1 & SAINT_MAX; }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
    }
}

static void libsais_final_bwt_aux_scan_left_to_right_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i + 0]; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; SA[induction_bucket[T[p0]]++] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); if ((p0 & rm) == 0) { I[p0 / (rm + 1)] = induction_bucket[T[p0]]; }} else { SA[i + 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i + 1]; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; SA[induction_bucket[T[p1]]++] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); if ((p1 & rm) == 0) { I[p1 / (rm + 1)] = induction_bucket[T[p1]]; }} else { SA[i + 1] = p1 & SAINT_MAX; }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]]; } } else { SA[i] = p & SAINT_MAX; }
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_left_to_right_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t origin, sa_sint_t index, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i + 0]; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; SA[induction_bucket[T[p0]]++] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); if (p0 == origin) { index = induction_bucket[T[p0]] - 1; } } else { SA[i + 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i + 1]; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; SA[induction_bucket[T[p1]]++] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); if (p1 == origin) { index = induction_bucket[T[p1]] - 1; } } else { SA[i + 1] = p1 & SAINT_MAX; }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); if (p == origin) { index = induction_bucket[T[p]] - 1; } } else { SA[i] = p & SAINT_MAX; }
    }

    return index;
}

static void libsais_final_sorting_scan_left_to_right_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); }
    }
}

static void libsais_final_sorting_scan_left_to_right_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 2 * prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i + 2 * prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 1)]; libsais_prefetchr(Ts0 - 1);
        sa_sint_t s1 = SA[i + 2 * prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 1)]; libsais_prefetchr(Ts1 - 1);
        sa_sint_t s2 = SA[i + 1 * prefetch_distance + 0]; const sa_sint_t * Ts2 = &T[(sa_uint_t)(s2 > 0 ? s2 : 2)]; libsais_prefetchw(&induction_bucket[Ts2[-1]]); libsais_prefetchr(Ts2 - 2);
        sa_sint_t s3 = SA[i + 1 * prefetch_distance + 1]; const sa_sint_t * Ts3 = &T[(sa_uint_t)(s3 > 0 ? s3 : 2)]; libsais_prefetchw(&induction_bucket[Ts3[-1]]); libsais_prefetchr(Ts3 - 2);

        sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; SA[induction_bucket[T[p0]]++] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); }
        sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; SA[induction_bucket[T[p1]]++] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); }
    }

    for (j += 2 * prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); }
    }
}

#if defined(LIBSAIS_OPENMP)

static fast_sint_t libsais_final_bwt_scan_left_to_right_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
   {
       libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i + 0]; if (p0 > 0) { p0--; SA[i + 0] = T[p0] | SAINT_MIN; buckets[cache[count].symbol = T[p0]]++; cache[count++].index = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); } else { SA[i + 0] = p0 & SAINT_MAX; }
       sa_sint_t p1 = SA[i + 1]; if (p1 > 0) { p1--; SA[i + 1] = T[p1] | SAINT_MIN; buckets[cache[count].symbol = T[p1]]++; cache[count++].index = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); } else { SA[i + 1] = p1 & SAINT_MAX; }
   }

   for (j += prefetch_distance + 1; i < j; i += 1)
   {
       sa_sint_t p = SA[i]; if (p > 0) { p--; SA[i] = T[p] | SAINT_MIN; buckets[cache[count].symbol = T[p]]++; cache[count++].index = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
   }

   return count;
}

static fast_sint_t libsais_final_sorting_scan_left_to_right_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
   {
       libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i + prefetch_distance + 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i + prefetch_distance + 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; buckets[cache[count].symbol = T[p0]]++; cache[count++].index = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); }
       sa_sint_t p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; buckets[cache[count].symbol = T[p1]]++; cache[count++].index = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); }
   }

   for (j += prefetch_distance + 1; i < j; i += 1)
   {
       sa_sint_t p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; buckets[cache[count].symbol = T[p]]++; cache[count++].index = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); }
   }

   return count;
}

static void libsais_final_order_scan_left_to_right_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = count - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[buckets[cache[i + 0].symbol]++] = cache[i + 0].index;
        SA[buckets[cache[i + 1].symbol]++] = cache[i + 1].index;
        SA[buckets[cache[i + 2].symbol]++] = cache[i + 2].index;
        SA[buckets[cache[i + 3].symbol]++] = cache[i + 3].index;
    }

    for (j += 3; i < j; i += 1)
    {
        SA[buckets[cache[i].symbol]++] = cache[i].index;
    }
}

static void libsais_final_bwt_aux_scan_left_to_right_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = count - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[buckets[cache[i + 0].symbol]++] = cache[i + 0].index; if ((cache[i + 0].index & rm) == 0) { I[(cache[i + 0].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i + 0].symbol]; }
        SA[buckets[cache[i + 1].symbol]++] = cache[i + 1].index; if ((cache[i + 1].index & rm) == 0) { I[(cache[i + 1].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i + 1].symbol]; }
        SA[buckets[cache[i + 2].symbol]++] = cache[i + 2].index; if ((cache[i + 2].index & rm) == 0) { I[(cache[i + 2].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i + 2].symbol]; }
        SA[buckets[cache[i + 3].symbol]++] = cache[i + 3].index; if ((cache[i + 3].index & rm) == 0) { I[(cache[i + 3].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i + 3].symbol]; }
    }

    for (j += 3; i < j; i += 1)
    {
        SA[buckets[cache[i].symbol]++] = cache[i].index; if ((cache[i].index & rm) == 0) { I[(cache[i].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i].symbol]; }
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_left_to_right_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t origin, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j; sa_sint_t index = -1;
    for (i = 0, j = count - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[buckets[cache[i + 0].symbol]++] = cache[i + 0].index; if ((cache[i + 0].index & SAINT_MAX) == origin) { index = buckets[cache[i + 0].symbol] - 1; }
        SA[buckets[cache[i + 1].symbol]++] = cache[i + 1].index; if ((cache[i + 1].index & SAINT_MAX) == origin) { index = buckets[cache[i + 1].symbol] - 1; }
        SA[buckets[cache[i + 2].symbol]++] = cache[i + 2].index; if ((cache[i + 2].index & SAINT_MAX) == origin) { index = buckets[cache[i + 2].symbol] - 1; }
        SA[buckets[cache[i + 3].symbol]++] = cache[i + 3].index; if ((cache[i + 3].index & SAINT_MAX) == origin) { index = buckets[cache[i + 3].symbol] - 1; }
    }

    for (j += 3; i < j; i += 1)
    {
        SA[buckets[cache[i].symbol]++] = cache[i].index; if ((cache[i].index & SAINT_MAX) == origin) { index = buckets[cache[i].symbol] - 1; }
    }

    return index;
}

static void libsais_final_sorting_scan_left_to_right_32s_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; SA[i + 0] = p0 ^ SAINT_MIN; if (p0 > 0) { p0--; cache[i + 0].index = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] < T[p0]) << (SAINT_BIT - 1)); symbol0 = T[p0]; } cache[i + 0].symbol = symbol0;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; SA[i + 1] = p1 ^ SAINT_MIN; if (p1 > 0) { p1--; cache[i + 1].index = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] < T[p1]) << (SAINT_BIT - 1)); symbol1 = T[p1]; } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; SA[i] = p ^ SAINT_MIN; if (p > 0) { p--; cache[i].index = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); symbol = T[p]; } cache[i].symbol = symbol;
    }
}

static void libsais_final_sorting_scan_left_to_right_32s_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, omp_block_end = omp_block_start + omp_block_size;
    for (i = omp_block_start, j = omp_block_end - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&cache[i + 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i + prefetch_distance + 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 : 0]; libsais_prefetchw(Is0);
        sa_sint_t s1 = cache[i + prefetch_distance + 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 : 0]; libsais_prefetchw(Is1);
        
        sa_sint_t v0 = cache[i + 0].symbol;
        if (v0 >= 0)
        {
            cache[i + 0].symbol = induction_bucket[v0]++;
            if (cache[i + 0].symbol < omp_block_end) { sa_sint_t ni = cache[i + 0].symbol, np = cache[i + 0].index; cache[i + 0].index = np ^ SAINT_MIN; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] < T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } }
        }

        sa_sint_t v1 = cache[i + 1].symbol;
        if (v1 >= 0)
        {
            cache[i + 1].symbol = induction_bucket[v1]++;
            if (cache[i + 1].symbol < omp_block_end) { sa_sint_t ni = cache[i + 1].symbol, np = cache[i + 1].index; cache[i + 1].index = np ^ SAINT_MIN; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] < T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } }
        }
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            cache[i].symbol = induction_bucket[v]++;
            if (cache[i].symbol < omp_block_end) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; cache[i].index = np ^ SAINT_MIN; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] < T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } }
        }
    }
}

#define LIBSAIS_INDUCING_WALK_RATIO     (32)


static fast_sint_t libsais_popcount64(uint64_t x)
{
    x = x - ((x >> 1) & (uint64_t)0x5555555555555555ull);
    x = (x & (uint64_t)0x3333333333333333ull) + ((x >> 2) & (uint64_t)0x3333333333333333ull);
    x = (x + (x >> 4)) & (uint64_t)0x0f0f0f0f0f0f0f0full;

    return (fast_sint_t)((x * (uint64_t)0x0101010101010101ull) >> 56);
}

static fast_sint_t libsais_final_sorting_scan_32s_block_masks(const LIBSAIS_THREAD_CACHE * RESTRICT cache, uint64_t * RESTRICT masks, fast_sint_t block_start, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    fast_sint_t i, j, switches = 0, omp_block_end = omp_block_start + omp_block_size; uint64_t last = 0;
    for (i = omp_block_start; i < omp_block_end; i += 64)
    {
        fast_sint_t group_size = omp_block_end - i < 64 ? omp_block_end - i : 64; uint64_t negative = 0;

#if (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2) && (SAINT_BIT == 32)
        if (group_size == 64)
        {
            for (j = 0; j < 64; j += 8)
            {
                __m128 s0 = _mm_loadu_ps((const float *)(const void *)&cache[i + j + 0]);
                __m128 s1 = _mm_loadu_ps((const float *)(const void *)&cache[i + j + 2]);
                __m128 s2 = _mm_loadu_ps((const float *)(const void *)&cache[i + j + 4]);
                __m128 s3 = _mm_loadu_ps((const float *)(const void *)&cache[i + j + 6]);

                negative |= (uint64_t)(_mm_movemask_ps(_mm_shuffle_ps(s0, s1, _MM_SHUFFLE(2, 0, 2, 0))) | (_mm_movemask_ps(_mm_shuffle_ps(s2, s3, _MM_SHUFFLE(2, 0, 2, 0))) << 4)) << j;
            }
        }
        else
#endif
        {
            for (j = 0; j < group_size; j += 1) { negative |= (uint64_t)((sa_uint_t)cache[i + j].symbol >> (SAINT_BIT - 1)) << j; }
        }

        {
            uint64_t valid = group_size < 64 ? (((uint64_t)1 << group_size) - 1) : ~(uint64_t)0, mask = ~negative & valid;

            switches += libsais_popcount64((mask ^ (mask >> 1)) & (valid >> 1)) + (i > omp_block_start ? (fast_sint_t)((mask ^ last) & 1) : 0);
            last = mask >> (group_size - 1);

            masks[(i - block_start) >> 6] = mask;
        }
    }

    return switches;
}

static LIBSAIS_NOINLINE void libsais_final_sorting_scan_left_to_right_32s_block_walk(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, uint64_t * RESTRICT masks, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    fast_sint_t g, omp_block_end = omp_block_start + omp_block_size, groups = (omp_block_size + 63) >> 6;
    for (g = 0; g < groups; g += 1)
    {
        fast_sint_t group_start = omp_block_start + (g << 6);

        if (group_start + 192 <= omp_block_end)
        {
            const LIBSAIS_THREAD_CACHE * RESTRICT next = &cache[group_start + 64]; fast_sint_t p;

            libsais_prefetchw(&next[64]); libsais_prefetchw(&next[72]); libsais_prefetchw(&next[80]); libsais_prefetchw(&next[88]);
            libsais_prefetchw(&next[96]); libsais_prefetchw(&next[104]); libsais_prefetchw(&next[112]); libsais_prefetchw(&next[120]);

            for (p = 0; p < 64; p += 8)
            {
                sa_sint_t s0 = next[p + 0].symbol; libsais_prefetchw(&induction_bucket[s0 > 0 ? s0 : 0]);
                sa_sint_t s1 = next[p + 1].symbol; libsais_prefetchw(&induction_bucket[s1 > 0 ? s1 : 0]);
                sa_sint_t s2 = next[p + 2].symbol; libsais_prefetchw(&induction_bucket[s2 > 0 ? s2 : 0]);
                sa_sint_t s3 = next[p + 3].symbol; libsais_prefetchw(&induction_bucket[s3 > 0 ? s3 : 0]);
                sa_sint_t s4 = next[p + 4].symbol; libsais_prefetchw(&induction_bucket[s4 > 0 ? s4 : 0]);
                sa_sint_t s5 = next[p + 5].symbol; libsais_prefetchw(&induction_bucket[s5 > 0 ? s5 : 0]);
                sa_sint_t s6 = next[p + 6].symbol; libsais_prefetchw(&induction_bucket[s6 > 0 ? s6 : 0]);
                sa_sint_t s7 = next[p + 7].symbol; libsais_prefetchw(&induction_bucket[s7 > 0 ? s7 : 0]);
            }
        }

        uint64_t mask = masks[g];
        while (mask != 0)
        {
            fast_sint_t i = group_start + libsais_ctz64(mask); mask &= mask - 1;

            sa_sint_t ni = induction_bucket[cache[i].symbol]++; cache[i].symbol = ni;
            if (ni < omp_block_end)
            {
                sa_sint_t np = cache[i].index; cache[i].index = np ^ SAINT_MIN;
                if (np > 0)
                {
                    np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] < T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np];

                    fast_sint_t d = (fast_sint_t)ni - omp_block_start; uint64_t bit = (uint64_t)1 << (d & 63);
                    masks[d >> 6] |= bit; mask |= bit & ((uint64_t)0 - (uint64_t)((d >> 6) == g));
                }
            }
        }
    }
}

static void libsais_final_bwt_scan_left_to_right_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_bwt_scan_left_to_right_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_bwt_scan_left_to_right_8u(T, SA, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_scan_left_to_right_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A + B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_order_scan_left_to_right_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static void libsais_final_bwt_aux_scan_left_to_right_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_bwt_aux_scan_left_to_right_8u(T, SA, rm, I, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_bwt_aux_scan_left_to_right_8u(T, SA, rm, I, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_scan_left_to_right_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A + B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_bwt_aux_scan_left_to_right_8u_block_place(SA, rm, I, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_left_to_right_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t origin, sa_sint_t index, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            index = libsais_final_bwt_cyclic_scan_left_to_right_8u(T, SA, origin, index, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                index = libsais_final_bwt_cyclic_scan_left_to_right_8u(T, SA, origin, index, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_scan_left_to_right_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A + B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                sa_sint_t block_index = libsais_final_bwt_cyclic_scan_left_to_right_8u_block_place(SA, origin, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
                if (block_index >= 0) { index = block_index; }
            }
        }
#endif
    }

    return index;
}

static void libsais_final_sorting_scan_left_to_right_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_sorting_scan_left_to_right_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_sorting_scan_left_to_right_8u(T, SA, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_sorting_scan_left_to_right_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = 0; t < omp_num_threads; ++t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A + B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_order_scan_left_to_right_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static void libsais_final_sorting_scan_left_to_right_32s_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_STATE * RESTRICT thread_state, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-64);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_sorting_scan_left_to_right_32s(T, SA, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            LIBSAIS_THREAD_CACHE * RESTRICT cache = thread_state[0].state.cache; uint64_t * RESTRICT masks = (uint64_t *)(void *)thread_state[0].state.buckets;

            {
                libsais_final_sorting_scan_left_to_right_32s_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
                thread_state[omp_thread_num].state.count = libsais_final_sorting_scan_32s_block_masks(cache - block_start, masks, block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, switches = 0; for (t = 0; t < omp_num_threads; ++t) { switches += thread_state[t].state.count; }

                if (LIBSAIS_INDUCING_WALK_RATIO * switches > block_size)
                {
                    libsais_final_sorting_scan_left_to_right_32s_block_walk(T, buckets, cache - block_start, masks, block_start, block_size);
                }
                else
                {
                    libsais_final_sorting_scan_left_to_right_32s_block_sort(T, buckets, cache - block_start, block_start, block_size);
                }
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

static fast_sint_t libsais_final_order_scan_left_to_right_8u_block_end(const sa_sint_t * RESTRICT SA, fast_sint_t block_start, fast_sint_t block_max_end)
{
    if (block_max_end > block_start + 32) { block_max_end = block_start + 32; }

    fast_sint_t block_end = block_start + 1; while (block_end < block_max_end && SA[block_end] != 0) { block_end++; }

    return block_end;
}

static fast_sint_t libsais_final_order_scan_left_to_right_8u_bucket_head(sa_sint_t k, const sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_end, fast_sint_t block_max_end)
{
    fast_sint_t c = 0, count = k; while (count > 0) { fast_sint_t half = count >> 1; if (induction_bucket[c + half] < block_end) { c += half + 1; count -= half + 1; } else { count = half; } }

    return c < k && induction_bucket[c] < block_max_end ? induction_bucket[c] : block_max_end;
}

static fast_sint_t libsais_final_order_scan_left_to_right_8u_block_start(const sa_sint_t * RESTRICT SA, fast_sint_t block_start, fast_sint_t n)
{
    for (; block_start + 16 <= n; block_start += 16)
    {
        const sa_sint_t * RESTRICT S = &SA[block_start];

        if ((S[0] | S[1] | S[2] | S[3] | S[4] | S[5] | S[6] | S[7] | S[8] | S[9] | S[10] | S[11] | S[12] | S[13] | S[14] | S[15]) != 0) { break; }
    }

    return block_start;
}

#endif

static void libsais_final_bwt_scan_left_to_right_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | (sa_sint_t)((sa_uint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    if (threads == 1 || n < 65536)
    {
        libsais_final_bwt_scan_left_to_right_8u(T, SA, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = 0; block_start < n; )
        {
            if (SA[block_start] == 0)
            {
                block_start = libsais_final_order_scan_left_to_right_8u_block_start(SA, block_start + 1, n);
            }
            else
            {
                fast_sint_t block_max_end = block_start + ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end > n) { block_max_end = n;}
                fast_sint_t block_end     = libsais_final_order_scan_left_to_right_8u_block_end(SA, block_start, block_max_end); if (block_end < block_max_end && SA[block_end] != 0) { block_end = libsais_final_order_scan_left_to_right_8u_bucket_head(k, induction_bucket, block_end, block_max_end); }
                fast_sint_t block_size    = block_end - block_start;

                if (block_size < 32)
                {
                    for (; block_start < block_end; block_start += 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; SA[block_start] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_bwt_scan_left_to_right_8u_block_omp(T, SA, k, induction_bucket, block_start, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static void libsais_final_bwt_aux_scan_left_to_right_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | (sa_sint_t)((sa_uint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    if ((((sa_sint_t)n - 1) & rm) == 0) { I[((sa_sint_t)n - 1) / (rm + 1)] = induction_bucket[T[(sa_sint_t)n - 1]]; }

    if (threads == 1 || n < 65536)
    {
        libsais_final_bwt_aux_scan_left_to_right_8u(T, SA, rm, I, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = 0; block_start < n; )
        {
            if (SA[block_start] == 0)
            {
                block_start = libsais_final_order_scan_left_to_right_8u_block_start(SA, block_start + 1, n);
            }
            else
            {
                fast_sint_t block_max_end = block_start + ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end > n) { block_max_end = n;}
                fast_sint_t block_end     = libsais_final_order_scan_left_to_right_8u_block_end(SA, block_start, block_max_end); if (block_end < block_max_end && SA[block_end] != 0) { block_end = libsais_final_order_scan_left_to_right_8u_bucket_head(k, induction_bucket, block_end, block_max_end); }
                fast_sint_t block_size    = block_end - block_start;

                if (block_size < 32)
                {
                    for (; block_start < block_end; block_start += 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; SA[block_start] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]]; } } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_bwt_aux_scan_left_to_right_8u_block_omp(T, SA, k, rm, I, induction_bucket, block_start, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static sa_sint_t libsais_final_bwt_cyclic_scan_left_to_right_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t origin, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t index = -1;

    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | (sa_sint_t)((sa_uint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    if (((sa_sint_t)n - 1) == origin) { index = induction_bucket[T[(sa_sint_t)n - 1]] - 1; }

    if (threads == 1 || n < 65536)
    {
        index = libsais_final_bwt_cyclic_scan_left_to_right_8u(T, SA, origin, index, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = 0; block_start < n; )
        {
            if (SA[block_start] == 0)
            {
                block_start = libsais_final_order_scan_left_to_right_8u_block_start(SA, block_start + 1, n);
            }
            else
            {
                fast_sint_t block_max_end = block_start + ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end > n) { block_max_end = n;}
                fast_sint_t block_end     = libsais_final_order_scan_left_to_right_8u_block_end(SA, block_start, block_max_end); if (block_end < block_max_end && SA[block_end] != 0) { block_end = libsais_final_order_scan_left_to_right_8u_bucket_head(k, induction_bucket, block_end, block_max_end); }
                fast_sint_t block_size    = block_end - block_start;

                if (block_size < 32)
                {
                    for (; block_start < block_end; block_start += 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; SA[block_start] = T[p] | SAINT_MIN; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); if (p == origin) { index = induction_bucket[T[p]] - 1; } } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    index = libsais_final_bwt_cyclic_scan_left_to_right_8u_block_omp(T, SA, k, origin, index, induction_bucket, block_start, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif

    return index;
}

static void libsais_final_sorting_scan_left_to_right_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[induction_bucket[T[(sa_sint_t)n - 1]]++] = ((sa_sint_t)n - 1) | (sa_sint_t)((sa_uint_t)(T[(sa_sint_t)n - 2] < T[(sa_sint_t)n - 1]) << (SAINT_BIT - 1));

    if (threads == 1 || n < 65536)
    {
        libsais_final_sorting_scan_left_to_right_8u(T, SA, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = 0; block_start < n; )
        {
            if (SA[block_start] == 0)
            {
                block_start = libsais_final_order_scan_left_to_right_8u_block_start(SA, block_start + 1, n);
            }
            else
            {
                fast_sint_t block_max_end = block_start + ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end > n) { block_max_end = n;}
                fast_sint_t block_end     = libsais_final_order_scan_left_to_right_8u_block_end(SA, block_start, block_max_end); if (block_end < block_max_end && SA[block_end] != 0) { block_end = libsais_final_order_scan_left_to_right_8u_bucket_head(k, induction_bucket, block_end, block_max_end); }
                fast_sint_t block_size    = block_end - block_start;

                if (block_size < 32)
                {
                    for (; block_start < block_end; block_start += 1)
                    {
                        sa_sint_t p = SA[block_start]; SA[block_start] = p ^ SAINT_MIN; if (p > 0) { p--; SA[induction_bucket[T[p]]++] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] < T[p]) << (SAINT_BIT - 1)); }
                    }
                }
                else
                {
                    libsais_final_sorting_scan_left_to_right_8u_block_omp(T, SA, k, induction_bucket, block_start, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static void libsais_final_sorting_scan_left_to_right_32s_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    SA[induction_bucket[T[n - 1]]++] = (n - 1) | (sa_sint_t)((sa_uint_t)(T[n - 2] < T[n - 1]) << (SAINT_BIT - 1));

    if (threads == 1 || n < 65536)
    {
        libsais_final_sorting_scan_left_to_right_32s(T, SA, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = 0; block_start < n; block_start = block_end)
        {
            block_end = block_start + libsais_get_block_size_32s(n, threads); if (block_end > n) { block_end = n; }

            libsais_final_sorting_scan_left_to_right_32s_block_omp(T, SA, induction_bucket, thread_state, block_start, block_end - block_start, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static sa_sint_t libsais_final_bwt_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, index = -1;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - 0];
        if (p0 > 0) { p0--; uint8_t c0 = T[(sa_uint_t)(p0 - (p0 > 0))], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; } else { SA[i - 0] = p0 & SAINT_MAX; index = (p0 == 0) ? (i - 0) : index; }

        sa_sint_t p1 = SA[i - 1];
        if (p1 > 0) { p1--; uint8_t c0 = T[(sa_uint_t)(p1 - (p1 > 0))], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; } else { SA[i - 1] = p1 & SAINT_MAX; index = (p1 == 0) ? (i - 1) : index; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i];
        if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; } else { SA[i] = p & SAINT_MAX; index = (p == 0) ? i : index; }
    }

    return (sa_sint_t)index;
}

static void libsais_final_bwt_aux_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - 0];
        if (p0 > 0) { p0--; uint8_t c0 = T[(sa_uint_t)(p0 - (p0 > 0))], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; if ((p0 & rm) == 0) { I[p0 / (rm + 1)] = induction_bucket[T[p0]] + 1; } } else { SA[i - 0] = p0 & SAINT_MAX; }

        sa_sint_t p1 = SA[i - 1];
        if (p1 > 0) { p1--; uint8_t c0 = T[(sa_uint_t)(p1 - (p1 > 0))], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; if ((p1 & rm) == 0) { I[p1 / (rm + 1)] = induction_bucket[T[p1]] + 1; } } else { SA[i - 1] = p1 & SAINT_MAX; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i];
        if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]] + 1; } } else { SA[i] = p & SAINT_MAX; }
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t origin, sa_sint_t index, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - 0];
        if (p0 > 0) { p0--; uint8_t c0 = T[(p0 > 0 ? p0 : n) - 1], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p0 : t; if (p0 == origin) { index = induction_bucket[T[p0]]; } } else { SA[i - 0] = p0 & SAINT_MAX; }

        sa_sint_t p1 = SA[i - 1];
        if (p1 > 0) { p1--; uint8_t c0 = T[(p1 > 0 ? p1 : n) - 1], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p1 : t; if (p1 == origin) { index = induction_bucket[T[p1]]; } } else { SA[i - 1] = p1 & SAINT_MAX; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i];
        if (p > 0) { p--; uint8_t c0 = T[(p > 0 ? p : n) - 1], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; if (p == origin) { index = induction_bucket[T[p]]; } } else { SA[i] = p & SAINT_MAX; }
    }

    return index;
}

static void libsais_final_sorting_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] > T[p0]) << (SAINT_BIT - 1)); } else { SA[i - 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] > T[p1]) << (SAINT_BIT - 1)); } else { SA[i - 1] = p1 & SAINT_MAX; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
    }
}

static void libsais_final_gsa_scan_right_to_left_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        sa_sint_t p0 = SA[i - 0]; if (p0 > 0 && T[p0 - 1] > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] > T[p0]) << (SAINT_BIT - 1)); } else { SA[i - 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i - 1]; if (p1 > 0 && T[p1 - 1] > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] > T[p1]) << (SAINT_BIT - 1)); } else { SA[i - 1] = p1 & SAINT_MAX; }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; if (p > 0 && T[p - 1] > 0) { p--; SA[--induction_bucket[T[p]]] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
    }
}

static void libsais_final_sorting_scan_right_to_left_32s(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT induction_bucket, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + 2 * prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&SA[i - 3 * prefetch_distance]);

        sa_sint_t s0 = SA[i - 2 * prefetch_distance - 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 1)]; libsais_prefetchr(Ts0 - 1);
        sa_sint_t s1 = SA[i - 2 * prefetch_distance - 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 1)]; libsais_prefetchr(Ts1 - 1);
        sa_sint_t s2 = SA[i - 1 * prefetch_distance - 0]; const sa_sint_t * Ts2 = &T[(sa_uint_t)(s2 > 0 ? s2 : 2)]; libsais_prefetchw(&induction_bucket[Ts2[-1]]); libsais_prefetchr(Ts2 - 2);
        sa_sint_t s3 = SA[i - 1 * prefetch_distance - 1]; const sa_sint_t * Ts3 = &T[(sa_uint_t)(s3 > 0 ? s3 : 2)]; libsais_prefetchw(&induction_bucket[Ts3[-1]]); libsais_prefetchr(Ts3 - 2);

        sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; SA[--induction_bucket[T[p0]]] = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] > T[p0]) << (SAINT_BIT - 1)); } else { SA[i - 0] = p0 & SAINT_MAX; }
        sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; SA[--induction_bucket[T[p1]]] = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] > T[p1]) << (SAINT_BIT - 1)); } else { SA[i - 1] = p1 & SAINT_MAX; }
    }

    for (j -= 2 * prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t p = SA[i]; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
    }
}

#if defined(LIBSAIS_OPENMP)

static fast_sint_t libsais_final_bwt_scan_right_to_left_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
   {
       libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; uint8_t c0 = T[(sa_uint_t)(p0 - (p0 > 0))], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count++].index = (c0 <= c1) ? p0 : t; } else { SA[i - 0] = p0 & SAINT_MAX; }
       sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; uint8_t c0 = T[(sa_uint_t)(p1 - (p1 > 0))], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count++].index = (c0 <= c1) ? p1 : t; } else { SA[i - 1] = p1 & SAINT_MAX; }
   }

   for (j -= prefetch_distance + 1; i >= j; i -= 1)
   {
       sa_sint_t p = SA[i]; if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count++].index = (c0 <= c1) ? p : t; } else { SA[i] = p & SAINT_MAX; }
   }

   return count;
}

static fast_sint_t libsais_final_bwt_aux_scan_right_to_left_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
   {
       libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; uint8_t c0 = T[(sa_uint_t)(p0 - (p0 > 0))], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p0 : t; cache[count + 1].index = p0; count += 2; } else { SA[i - 0] = p0 & SAINT_MAX; }
       sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; uint8_t c0 = T[(sa_uint_t)(p1 - (p1 > 0))], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p1 : t; cache[count + 1].index = p1; count += 2; } else { SA[i - 1] = p1 & SAINT_MAX; }
   }

   for (j -= prefetch_distance + 1; i >= j; i -= 1)
   {
       sa_sint_t p = SA[i]; if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p : t; cache[count + 1].index = p; count += 2; } else { SA[i] = p & SAINT_MAX; }
   }

   return count;
}

static fast_sint_t libsais_final_bwt_cyclic_scan_right_to_left_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
   {
       libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; uint8_t c0 = T[(p0 > 0 ? p0 : n) - 1], c1 = T[p0]; SA[i - 0] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p0 : t; cache[count + 1].index = p0; count += 2; } else { SA[i - 0] = p0 & SAINT_MAX; }
       sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; uint8_t c0 = T[(p1 > 0 ? p1 : n) - 1], c1 = T[p1]; SA[i - 1] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p1 : t; cache[count + 1].index = p1; count += 2; } else { SA[i - 1] = p1 & SAINT_MAX; }
   }

   for (j -= prefetch_distance + 1; i >= j; i -= 1)
   {
       sa_sint_t p = SA[i]; if (p > 0) { p--; uint8_t c0 = T[(p > 0 ? p : n) - 1], c1 = T[p]; SA[i] = c1; sa_sint_t t = c0 | SAINT_MIN; buckets[cache[count].symbol = c1]++; cache[count].index = (c0 <= c1) ? p : t; cache[count + 1].index = p; count += 2; } else { SA[i] = p & SAINT_MAX; }
   }

   return count;
}

static fast_sint_t libsais_final_sorting_scan_right_to_left_8u_block_prepare(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
   const fast_sint_t prefetch_distance = 64;

   memset(buckets, 0, (size_t)k * sizeof(sa_sint_t));

   fast_sint_t i, j, count = 0;
   for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
   {
       libsais_prefetchw(&SA[i - 2 * prefetch_distance]);

       sa_sint_t s0 = SA[i - prefetch_distance - 0]; const uint8_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
       sa_sint_t s1 = SA[i - prefetch_distance - 1]; const uint8_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

       sa_sint_t p0 = SA[i - 0]; if (p0 > 0) { p0--; buckets[cache[count].symbol = T[p0]]++; cache[count++].index = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] > T[p0]) << (SAINT_BIT - 1)); } else { SA[i - 0] = p0 & SAINT_MAX; }
       sa_sint_t p1 = SA[i - 1]; if (p1 > 0) { p1--; buckets[cache[count].symbol = T[p1]]++; cache[count++].index = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] > T[p1]) << (SAINT_BIT - 1)); } else { SA[i - 1] = p1 & SAINT_MAX; }
   }

   for (j -= prefetch_distance + 1; i >= j; i -= 1)
   {
       sa_sint_t p = SA[i]; if (p > 0) { p--; buckets[cache[count].symbol = T[p]]++; cache[count++].index = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[i] = p & SAINT_MAX; }
   }

   return count;
}

static void libsais_final_order_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = count - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[--buckets[cache[i + 0].symbol]] = cache[i + 0].index;
        SA[--buckets[cache[i + 1].symbol]] = cache[i + 1].index;
        SA[--buckets[cache[i + 2].symbol]] = cache[i + 2].index;
        SA[--buckets[cache[i + 3].symbol]] = cache[i + 3].index;
    }

    for (j += 3; i < j; i += 1)
    {
        SA[--buckets[cache[i].symbol]] = cache[i].index;
    }
}

static void libsais_final_gsa_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = count - 3; i < j; i += 4)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        if (cache[i + 0].symbol > 0) { SA[--buckets[cache[i + 0].symbol]] = cache[i + 0].index; }
        if (cache[i + 1].symbol > 0) { SA[--buckets[cache[i + 1].symbol]] = cache[i + 1].index; }
        if (cache[i + 2].symbol > 0) { SA[--buckets[cache[i + 2].symbol]] = cache[i + 2].index; }
        if (cache[i + 3].symbol > 0) { SA[--buckets[cache[i + 3].symbol]] = cache[i + 3].index; }
    }

    for (j += 3; i < j; i += 1)
    {
        if (cache[i].symbol > 0) { SA[--buckets[cache[i].symbol]] = cache[i].index; }
    }
}

static void libsais_final_bwt_aux_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = count - 6; i < j; i += 8)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[--buckets[cache[i + 0].symbol]] = cache[i + 0].index; if ((cache[i + 1].index & rm) == 0) { I[cache[i + 1].index / (rm + 1)] = buckets[cache[i + 0].symbol] + 1; }
        SA[--buckets[cache[i + 2].symbol]] = cache[i + 2].index; if ((cache[i + 3].index & rm) == 0) { I[cache[i + 3].index / (rm + 1)] = buckets[cache[i + 2].symbol] + 1; }
        SA[--buckets[cache[i + 4].symbol]] = cache[i + 4].index; if ((cache[i + 5].index & rm) == 0) { I[cache[i + 5].index / (rm + 1)] = buckets[cache[i + 4].symbol] + 1; }
        SA[--buckets[cache[i + 6].symbol]] = cache[i + 6].index; if ((cache[i + 7].index & rm) == 0) { I[cache[i + 7].index / (rm + 1)] = buckets[cache[i + 6].symbol] + 1; }
    }

    for (j += 6; i < j; i += 2)
    {
        SA[--buckets[cache[i].symbol]] = cache[i].index; if ((cache[i + 1].index & rm) == 0) { I[(cache[i + 1].index & SAINT_MAX) / (rm + 1)] = buckets[cache[i].symbol] + 1; }
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_right_to_left_8u_block_place(sa_sint_t * RESTRICT SA, sa_sint_t origin, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t count)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j; sa_sint_t index = -1;
    for (i = 0, j = count - 6; i < j; i += 8)
    {
        libsais_prefetchr(&cache[i + prefetch_distance]);

        SA[--buckets[cache[i + 0].symbol]] = cache[i + 0].index; if (cache[i + 1].index == origin) { index = buckets[cache[i + 0].symbol]; }
        SA[--buckets[cache[i + 2].symbol]] = cache[i + 2].index; if (cache[i + 3].index == origin) { index = buckets[cache[i + 2].symbol]; }
        SA[--buckets[cache[i + 4].symbol]] = cache[i + 4].index; if (cache[i + 5].index == origin) { index = buckets[cache[i + 4].symbol]; }
        SA[--buckets[cache[i + 6].symbol]] = cache[i + 6].index; if (cache[i + 7].index == origin) { index = buckets[cache[i + 6].symbol]; }
    }

    for (j += 6; i < j; i += 2)
    {
        SA[--buckets[cache[i].symbol]] = cache[i].index; if (cache[i + 1].index == origin) { index = buckets[cache[i].symbol]; }
    }

    return index;
}

static void libsais_final_sorting_scan_right_to_left_32s_block_gather(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 1; i < j; i += 2)
    {
        libsais_prefetchw(&SA[i + 2 * prefetch_distance]);

        sa_sint_t s0 = SA[i + prefetch_distance + 0]; const sa_sint_t * Ts0 = &T[(sa_uint_t)(s0 > 0 ? s0 : 2)]; libsais_prefetchr(Ts0 - 1); libsais_prefetchr(Ts0 - 2);
        sa_sint_t s1 = SA[i + prefetch_distance + 1]; const sa_sint_t * Ts1 = &T[(sa_uint_t)(s1 > 0 ? s1 : 2)]; libsais_prefetchr(Ts1 - 1); libsais_prefetchr(Ts1 - 2);

        libsais_prefetchw(&cache[i + prefetch_distance]);

        sa_sint_t symbol0 = SAINT_MIN, p0 = SA[i + 0]; if (p0 > 0) { p0--; cache[i + 0].index = p0 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p0 - (p0 > 0))] > T[p0]) << (SAINT_BIT - 1)); symbol0 = T[p0]; } else { SA[i + 0] = p0 & SAINT_MAX; } cache[i + 0].symbol = symbol0;
        sa_sint_t symbol1 = SAINT_MIN, p1 = SA[i + 1]; if (p1 > 0) { p1--; cache[i + 1].index = p1 | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p1 - (p1 > 0))] > T[p1]) << (SAINT_BIT - 1)); symbol1 = T[p1]; } else { SA[i + 1] = p1 & SAINT_MAX; } cache[i + 1].symbol = symbol1;
    }

    for (j += prefetch_distance + 1; i < j; i += 1)
    {
        sa_sint_t symbol = SAINT_MIN, p = SA[i]; if (p > 0) { p--; cache[i].index = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); symbol = T[p]; } else { SA[i] = p & SAINT_MAX; } cache[i].symbol = symbol;
    }
}

static void libsais_final_sorting_scan_right_to_left_32s_block_sort(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start + omp_block_size - 1, j = omp_block_start + prefetch_distance + 1; i >= j; i -= 2)
    {
        libsais_prefetchw(&cache[i - 2 * prefetch_distance]);

        sa_sint_t s0 = cache[i - prefetch_distance - 0].symbol; const sa_sint_t * Is0 = &induction_bucket[s0 > 0 ? s0 : 0]; libsais_prefetchw(Is0);
        sa_sint_t s1 = cache[i - prefetch_distance - 1].symbol; const sa_sint_t * Is1 = &induction_bucket[s1 > 0 ? s1 : 0]; libsais_prefetchw(Is1);

        sa_sint_t v0 = cache[i - 0].symbol;
        if (v0 >= 0)
        {
            cache[i - 0].symbol = --induction_bucket[v0];
            if (cache[i - 0].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 0].symbol, np = cache[i - 0].index; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] > T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } else { cache[i - 0].index = np & SAINT_MAX; } }
        }

        sa_sint_t v1 = cache[i - 1].symbol;
        if (v1 >= 0)
        {
            cache[i - 1].symbol = --induction_bucket[v1];
            if (cache[i - 1].symbol >= omp_block_start) { sa_sint_t ni = cache[i - 1].symbol, np = cache[i - 1].index; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] > T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } else { cache[i - 1].index = np & SAINT_MAX; } }
        }
    }

    for (j -= prefetch_distance + 1; i >= j; i -= 1)
    {
        sa_sint_t v = cache[i].symbol;
        if (v >= 0)
        {
            cache[i].symbol = --induction_bucket[v];
            if (cache[i].symbol >= omp_block_start) { sa_sint_t ni = cache[i].symbol, np = cache[i].index; if (np > 0) { np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] > T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np]; } else { cache[i].index = np & SAINT_MAX; } }
        }
    }
}

#if defined(__GNUC__) || defined(__clang__)
    #define libsais_bsr64(x) ((fast_sint_t)(63 - __builtin_clzll((unsigned long long)(x))))
#elif defined(_M_X64) || defined(_M_ARM64)
    static fast_sint_t libsais_bsr64(uint64_t x) { unsigned long index; _BitScanReverse64(&index, x); return (fast_sint_t)index; }
#else
    static fast_sint_t libsais_bsr64(uint64_t x) { unsigned long index; if ((x >> 32) != 0) { _BitScanReverse(&index, (unsigned long)(uint32_t)(x >> 32)); return (fast_sint_t)index + 32; } _BitScanReverse(&index, (unsigned long)(uint32_t)x); return (fast_sint_t)index; }
#endif

static LIBSAIS_NOINLINE void libsais_final_sorting_scan_right_to_left_32s_block_walk(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT induction_bucket, LIBSAIS_THREAD_CACHE * RESTRICT cache, uint64_t * RESTRICT masks, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    fast_sint_t g;
    for (g = ((omp_block_size + 63) >> 6) - 1; g >= 0; g -= 1)
    {
        fast_sint_t group_start = omp_block_start + (g << 6);

        if (group_start - 128 >= omp_block_start)
        {
            const LIBSAIS_THREAD_CACHE * RESTRICT next = &cache[group_start - 64]; fast_sint_t p;

            libsais_prefetchw(&next[-64]); libsais_prefetchw(&next[-56]); libsais_prefetchw(&next[-48]); libsais_prefetchw(&next[-40]);
            libsais_prefetchw(&next[-32]); libsais_prefetchw(&next[-24]); libsais_prefetchw(&next[-16]); libsais_prefetchw(&next[-8]);

            for (p = 0; p < 64; p += 8)
            {
                sa_sint_t s0 = next[p + 0].symbol; libsais_prefetchw(&induction_bucket[s0 > 0 ? s0 : 0]);
                sa_sint_t s1 = next[p + 1].symbol; libsais_prefetchw(&induction_bucket[s1 > 0 ? s1 : 0]);
                sa_sint_t s2 = next[p + 2].symbol; libsais_prefetchw(&induction_bucket[s2 > 0 ? s2 : 0]);
                sa_sint_t s3 = next[p + 3].symbol; libsais_prefetchw(&induction_bucket[s3 > 0 ? s3 : 0]);
                sa_sint_t s4 = next[p + 4].symbol; libsais_prefetchw(&induction_bucket[s4 > 0 ? s4 : 0]);
                sa_sint_t s5 = next[p + 5].symbol; libsais_prefetchw(&induction_bucket[s5 > 0 ? s5 : 0]);
                sa_sint_t s6 = next[p + 6].symbol; libsais_prefetchw(&induction_bucket[s6 > 0 ? s6 : 0]);
                sa_sint_t s7 = next[p + 7].symbol; libsais_prefetchw(&induction_bucket[s7 > 0 ? s7 : 0]);
            }
        }

        uint64_t mask = masks[g];
        while (mask != 0)
        {
            fast_sint_t j = libsais_bsr64(mask); mask ^= (uint64_t)1 << j;
            fast_sint_t i = group_start + j;

            sa_sint_t ni = --induction_bucket[cache[i].symbol]; cache[i].symbol = ni;
            if (ni >= omp_block_start)
            {
                sa_sint_t np = cache[i].index;
                if (np > 0)
                {
                    np--; cache[ni].index = np | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(np - (np > 0))] > T[np]) << (SAINT_BIT - 1)); cache[ni].symbol = T[np];

                    fast_sint_t d = (fast_sint_t)ni - omp_block_start; uint64_t bit = (uint64_t)1 << (d & 63);
                    masks[d >> 6] |= bit; mask |= bit & ((uint64_t)0 - (uint64_t)((d >> 6) == g));
                }
                else
                {
                    cache[i].index = np & SAINT_MAX;
                }
            }
        }
    }
}

static void libsais_final_bwt_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_bwt_scan_right_to_left_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_bwt_scan_right_to_left_8u(T, SA, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A - B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_order_scan_right_to_left_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static void libsais_final_bwt_aux_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_bwt_aux_scan_right_to_left_8u(T, SA, rm, I, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_bwt_aux_scan_right_to_left_8u(T, SA, rm, I, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_aux_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A - B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_bwt_aux_scan_right_to_left_8u_block_place(SA, rm, I, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static sa_sint_t libsais_final_bwt_cyclic_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t origin, sa_sint_t index, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            index = libsais_final_bwt_cyclic_scan_right_to_left_8u(T, SA, n, origin, index, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                index = libsais_final_bwt_cyclic_scan_right_to_left_8u(T, SA, n, origin, index, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_bwt_cyclic_scan_right_to_left_8u_block_prepare(T, SA, n, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A - B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                sa_sint_t block_index = libsais_final_bwt_cyclic_scan_right_to_left_8u_block_place(SA, origin, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
                if (block_index >= 0) { index = block_index; }
            }
        }
#endif
    }

    return index;
}

static void libsais_final_sorting_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_sorting_scan_right_to_left_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_sorting_scan_right_to_left_8u(T, SA, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_sorting_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A - B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_order_scan_right_to_left_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static void libsais_final_gsa_scan_right_to_left_8u_block_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    sa_sint_t block_threads = libsais_get_block_threads_8u(block_size, k, threads);
    #pragma omp parallel num_threads(block_threads) if(block_threads > 1 && omp_get_dynamic() == 0)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(k); UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_gsa_scan_right_to_left_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else if (omp_num_threads < block_threads)
        {
            #pragma omp master
            {
                libsais_final_gsa_scan_right_to_left_8u(T, SA, induction_bucket, block_start, block_size);
            }
        }
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_final_sorting_scan_right_to_left_8u_block_prepare(T, SA, k, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t;
                for (t = omp_num_threads - 1; t >= 0; --t)
                {
                    sa_sint_t * RESTRICT temp_bucket = thread_state[t].state.buckets;
                    fast_sint_t c; for (c = 0; c < k; c += 1) { sa_sint_t A = induction_bucket[c], B = temp_bucket[c]; induction_bucket[c] = A - B; temp_bucket[c] = A; }
                }
            }

            #pragma omp barrier

            {
                libsais_final_gsa_scan_right_to_left_8u_block_place(SA, thread_state[omp_thread_num].state.buckets, thread_state[omp_thread_num].state.cache, thread_state[omp_thread_num].state.count);
            }
        }
#endif
    }
}

static void libsais_final_sorting_scan_right_to_left_32s_block_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT buckets, LIBSAIS_THREAD_STATE * RESTRICT thread_state, fast_sint_t block_start, fast_sint_t block_size, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && block_size >= 16384)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (block_size / omp_num_threads) & (-64);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : block_size - omp_block_start;

        omp_block_start += block_start;

        if (omp_num_threads == 1)
        {
            libsais_final_sorting_scan_right_to_left_32s(T, SA, buckets, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            LIBSAIS_THREAD_CACHE * RESTRICT cache = thread_state[0].state.cache; uint64_t * RESTRICT masks = (uint64_t *)(void *)thread_state[0].state.buckets;

            {
                libsais_final_sorting_scan_right_to_left_32s_block_gather(T, SA, cache - block_start, omp_block_start, omp_block_size);
                thread_state[omp_thread_num].state.count = libsais_final_sorting_scan_32s_block_masks(cache - block_start, masks, block_start, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, switches = 0; for (t = 0; t < omp_num_threads; ++t) { switches += thread_state[t].state.count; }

                if (LIBSAIS_INDUCING_WALK_RATIO * switches > block_size)
                {
                    libsais_final_sorting_scan_right_to_left_32s_block_walk(T, buckets, cache - block_start, masks, block_start, block_size);
                }
                else
                {
                    libsais_final_sorting_scan_right_to_left_32s_block_sort(T, buckets, cache - block_start, block_start, block_size);
                }
            }

            #pragma omp barrier

            {
                libsais_compact_and_place_cached_suffixes(SA, cache - block_start, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

#endif

static sa_sint_t libsais_final_bwt_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t index = -1;

    if (threads == 1 || n < 65536)
    {
        index = libsais_final_bwt_scan_right_to_left_8u(T, SA, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; )
        {
            if (SA[block_start] == 0)
            {
                index = (sa_sint_t)block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end < 0) { block_max_end = -1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[block_start] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_bwt_scan_right_to_left_8u_block_omp(T, SA, k, induction_bucket, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif

    return index;
}

static void libsais_final_bwt_aux_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t rm, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || n < 65536)
    {
        libsais_final_bwt_aux_scan_right_to_left_8u(T, SA, rm, I, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * ((LIBSAIS_PER_THREAD_CACHE_SIZE - 32 * (fast_sint_t)threads) / 2); if (block_max_end < 0) { block_max_end = -1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; uint8_t c0 = T[(sa_uint_t)(p - (p > 0))], c1 = T[p]; SA[block_start] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; if ((p & rm) == 0) { I[p / (rm + 1)] = induction_bucket[T[p]] + 1; } } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_bwt_aux_scan_right_to_left_8u_block_omp(T, SA, k, rm, I, induction_bucket, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static sa_sint_t libsais_final_bwt_cyclic_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t origin, sa_sint_t index, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || n < 65536)
    {
        index = libsais_final_bwt_cyclic_scan_right_to_left_8u(T, SA, n, origin, index, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * ((LIBSAIS_PER_THREAD_CACHE_SIZE - 32 * (fast_sint_t)threads) / 2); if (block_max_end < 0) { block_max_end = -1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; uint8_t c0 = T[(p > 0 ? p : n) - 1], c1 = T[p]; SA[block_start] = c1; sa_sint_t t = c0 | SAINT_MIN; SA[--induction_bucket[c1]] = (c0 <= c1) ? p : t; if (p == origin) { index = induction_bucket[T[p]]; } } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    index = libsais_final_bwt_cyclic_scan_right_to_left_8u_block_omp(T, SA, n, k, origin, index, induction_bucket, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif

    return index;
}

static void libsais_final_sorting_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || omp_block_size < 65536)
    {
        libsais_final_sorting_scan_right_to_left_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = omp_block_start + omp_block_size - 1; block_start >= omp_block_start; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end < omp_block_start) { block_max_end = omp_block_start - 1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0) { p--; SA[--induction_bucket[T[p]]] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_sorting_scan_right_to_left_8u_block_omp(T, SA, k, induction_bucket, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static void libsais_final_gsa_scan_right_to_left_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, fast_sint_t omp_block_start, fast_sint_t omp_block_size, sa_sint_t k, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || omp_block_size < 65536)
    {
        libsais_final_gsa_scan_right_to_left_8u(T, SA, induction_bucket, omp_block_start, omp_block_size);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start;
        for (block_start = omp_block_start + omp_block_size - 1; block_start >= omp_block_start; )
        {
            if (SA[block_start] == 0)
            {
                block_start--;
            }
            else
            {
                fast_sint_t block_max_end = block_start - ((fast_sint_t)threads) * (LIBSAIS_PER_THREAD_CACHE_SIZE - 16 * (fast_sint_t)threads); if (block_max_end < omp_block_start) { block_max_end = omp_block_start - 1; }
                fast_sint_t block_end     = libsais_find_block_end_right_to_left(SA, block_start, block_max_end); if (block_end > block_max_end && SA[block_end] != 0) { block_end = libsais_find_block_end_right_to_left_omp(SA, block_end, block_max_end, threads, thread_state); }
                fast_sint_t block_size    = block_start - block_end;

                if (block_size < 32)
                {
                    for (; block_start > block_end; block_start -= 1)
                    {
                        sa_sint_t p = SA[block_start]; if (p > 0 && T[p - 1] > 0) { p--; SA[--induction_bucket[T[p]]] = p | (sa_sint_t)((sa_uint_t)(T[(sa_uint_t)(p - (p > 0))] > T[p]) << (SAINT_BIT - 1)); } else { SA[block_start] = p & SAINT_MAX; }
                    }
                }
                else
                {
                    libsais_final_gsa_scan_right_to_left_8u_block_omp(T, SA, k, induction_bucket, block_end + 1, block_size, threads, thread_state);
                    block_start = block_end;
                }
            }
        }
    }
#else
    UNUSED(k); UNUSED(thread_state);
#endif
}

static void libsais_final_sorting_scan_right_to_left_32s_omp(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t * RESTRICT induction_bucket, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (threads == 1 || n < 65536)
    {
        libsais_final_sorting_scan_right_to_left_32s(T, SA, induction_bucket, 0, n);
    }
#if defined(LIBSAIS_OPENMP)
    else
    {
        fast_sint_t block_start, block_end;
        for (block_start = (fast_sint_t)n - 1; block_start >= 0; block_start = block_end)
        {
            block_end = block_start - libsais_get_block_size_32s(n, threads); if (block_end < 0) { block_end = -1; }

            libsais_final_sorting_scan_right_to_left_32s_block_omp(T, SA, induction_bucket, thread_state, block_end + 1, block_start - block_end, threads);
        }
    }
#else
    UNUSED(thread_state);
#endif
}

static void libsais_clear_lms_suffixes_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT bucket_start, sa_sint_t * RESTRICT bucket_end, sa_sint_t threads)
{
    fast_sint_t c;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel for schedule(static, 1) num_threads(threads) if(threads > 1 && n >= 65536)
#else
    UNUSED(threads); UNUSED(n);
#endif
    for (c = 0; c < k; ++c)
    {
        if (bucket_end[c] > bucket_start[c])
        {
            memset(&SA[bucket_start[c]], 0, ((size_t)bucket_end[c] - (size_t)bucket_start[c]) * sizeof(sa_sint_t));
        }
    }
}

static sa_sint_t libsais_induce_final_order_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t flags, sa_sint_t r, sa_sint_t * RESTRICT I, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if ((flags & LIBSAIS_FLAGS_BWT) == 0)
    {
        if (flags & LIBSAIS_FLAGS_GSA) { buckets[6 * ALPHABET_SIZE] = buckets[7 * ALPHABET_SIZE] - 1; }

        libsais_final_sorting_scan_left_to_right_8u_omp(T, SA, n, k, &buckets[6 * ALPHABET_SIZE], threads, thread_state);
        if (threads > 1 && n >= 65536) { libsais_clear_lms_suffixes_omp(SA, n, ALPHABET_SIZE, &buckets[6 * ALPHABET_SIZE], &buckets[7 * ALPHABET_SIZE], threads); }

        if (flags & LIBSAIS_FLAGS_GSA)
        {
            libsais_flip_suffix_markers_omp(SA, buckets[7 * ALPHABET_SIZE], threads);
            libsais_final_gsa_scan_right_to_left_8u_omp(T, SA, buckets[7 * ALPHABET_SIZE], (fast_sint_t)n - buckets[7 * ALPHABET_SIZE], k, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
        }
        else
        {
            libsais_final_sorting_scan_right_to_left_8u_omp(T, SA, 0, n, k, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
        }

        return 0;
    }
    else if (I != NULL)
    {
        libsais_final_bwt_aux_scan_left_to_right_8u_omp(T, SA, n, k, r - 1, I, &buckets[6 * ALPHABET_SIZE], threads, thread_state);
        if (threads > 1 && n >= 65536) { libsais_clear_lms_suffixes_omp(SA, n, ALPHABET_SIZE, &buckets[6 * ALPHABET_SIZE], &buckets[7 * ALPHABET_SIZE], threads); }
        libsais_final_bwt_aux_scan_right_to_left_8u_omp(T, SA, n, k, r - 1, I, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
        return 0;
    }
    else if (flags & LIBSAIS_FLAGS_CYCLIC)
    {
        sa_sint_t index = libsais_final_bwt_cyclic_scan_left_to_right_8u_omp(T, SA, n, k, r, &buckets[6 * ALPHABET_SIZE], threads, thread_state);
        if (threads > 1 && n >= 65536) { libsais_clear_lms_suffixes_omp(SA, n, ALPHABET_SIZE, &buckets[6 * ALPHABET_SIZE], &buckets[7 * ALPHABET_SIZE], threads); }
        return libsais_final_bwt_cyclic_scan_right_to_left_8u_omp(T, SA, n, k, r, index, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
    }
    else
    {
        libsais_final_bwt_scan_left_to_right_8u_omp(T, SA, n, k, &buckets[6 * ALPHABET_SIZE], threads, thread_state);
        if (threads > 1 && n >= 65536) { libsais_clear_lms_suffixes_omp(SA, n, ALPHABET_SIZE, &buckets[6 * ALPHABET_SIZE], &buckets[7 * ALPHABET_SIZE], threads); }
        return libsais_final_bwt_scan_right_to_left_8u_omp(T, SA, n, k, &buckets[7 * ALPHABET_SIZE], threads, thread_state);
    }
}

static void libsais_induce_final_order_32s_6k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_final_sorting_scan_left_to_right_32s_omp(T, SA, n, &buckets[4 * (fast_sint_t)k], threads, thread_state);
    libsais_final_sorting_scan_right_to_left_32s_omp(T, SA, n, &buckets[5 * (fast_sint_t)k], threads, thread_state);
}

static void libsais_induce_final_order_32s_4k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_final_sorting_scan_left_to_right_32s_omp(T, SA, n, &buckets[2 * (fast_sint_t)k], threads, thread_state);
    libsais_final_sorting_scan_right_to_left_32s_omp(T, SA, n, &buckets[3 * (fast_sint_t)k], threads, thread_state);
}

static void libsais_induce_final_order_32s_2k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_final_sorting_scan_left_to_right_32s_omp(T, SA, n, &buckets[1 * (fast_sint_t)k], threads, thread_state);
    libsais_final_sorting_scan_right_to_left_32s_omp(T, SA, n, &buckets[0 * (fast_sint_t)k], threads, thread_state);
}

static void libsais_induce_final_order_32s_1k(const sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t * RESTRICT buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_start_32s_1k(k, buckets);
    libsais_final_sorting_scan_left_to_right_32s_omp(T, SA, n, buckets, threads, thread_state);

    libsais_count_suffixes_32s(T, n, k, buckets);
    libsais_initialize_buckets_end_32s_1k(k, buckets);
    libsais_final_sorting_scan_right_to_left_32s_omp(T, SA, n, buckets, threads, thread_state);
}

static sa_sint_t libsais_renumber_unique_and_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t f, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    sa_sint_t i, j;
    for (i = (sa_sint_t)omp_block_start, j = (sa_sint_t)omp_block_start + (sa_sint_t)omp_block_size - 2 * (sa_sint_t)prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 3 * prefetch_distance]);

        libsais_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 0]) >> 1]);
        libsais_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 1]) >> 1]);
        libsais_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 2]) >> 1]);
        libsais_prefetchw(&SAm[((sa_uint_t)SA[i + 2 * prefetch_distance + 3]) >> 1]);

        sa_uint_t q0 = (sa_uint_t)SA[i + prefetch_distance + 0]; sa_sint_t * RESTRICT Tq0 = &T[q0]; libsais_prefetchw(SAm[q0 >> 1] < 0 ? Tq0 : &SAm[q0 >> 1]);
        sa_uint_t q1 = (sa_uint_t)SA[i + prefetch_distance + 1]; sa_sint_t * RESTRICT Tq1 = &T[q1]; libsais_prefetchw(SAm[q1 >> 1] < 0 ? Tq1 : &SAm[q1 >> 1]);
        sa_uint_t q2 = (sa_uint_t)SA[i + prefetch_distance + 2]; sa_sint_t * RESTRICT Tq2 = &T[q2]; libsais_prefetchw(SAm[q2 >> 1] < 0 ? Tq2 : &SAm[q2 >> 1]);
        sa_uint_t q3 = (sa_uint_t)SA[i + prefetch_distance + 3]; sa_sint_t * RESTRICT Tq3 = &T[q3]; libsais_prefetchw(SAm[q3 >> 1] < 0 ? Tq3 : &SAm[q3 >> 1]);

        sa_uint_t p0 = (sa_uint_t)SA[i + 0]; sa_sint_t s0 = SAm[p0 >> 1]; if (s0 < 0) { T[p0] |= SAINT_MIN; f++; s0 = i + 0 + SAINT_MIN + f; } SAm[p0 >> 1] = s0 - f;
        sa_uint_t p1 = (sa_uint_t)SA[i + 1]; sa_sint_t s1 = SAm[p1 >> 1]; if (s1 < 0) { T[p1] |= SAINT_MIN; f++; s1 = i + 1 + SAINT_MIN + f; } SAm[p1 >> 1] = s1 - f;
        sa_uint_t p2 = (sa_uint_t)SA[i + 2]; sa_sint_t s2 = SAm[p2 >> 1]; if (s2 < 0) { T[p2] |= SAINT_MIN; f++; s2 = i + 2 + SAINT_MIN + f; } SAm[p2 >> 1] = s2 - f;
        sa_uint_t p3 = (sa_uint_t)SA[i + 3]; sa_sint_t s3 = SAm[p3 >> 1]; if (s3 < 0) { T[p3] |= SAINT_MIN; f++; s3 = i + 3 + SAINT_MIN + f; } SAm[p3 >> 1] = s3 - f;
    }

    for (j += 2 * (sa_sint_t)prefetch_distance + 3; i < j; i += 1)
    {
        sa_uint_t p = (sa_uint_t)SA[i]; sa_sint_t s = SAm[p >> 1]; if (s < 0) { T[p] |= SAINT_MIN; f++; s = i + SAINT_MIN + f; } SAm[p >> 1] = s - f;
    }

    return f;
}

static void libsais_compact_unique_and_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t * pl, fast_sint_t * pr, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_uint_t * RESTRICT SAl = (sa_uint_t *)&SA[0];
    sa_uint_t * RESTRICT SAr = (sa_uint_t *)&SA[0];

    fast_sint_t i, j, l = *pl - 1, r = *pr - 1;
    for (i = (fast_sint_t)m + omp_block_start + omp_block_size - 1, j = (fast_sint_t)m + omp_block_start + 3; i >= j; i -= 4)
    {
        libsais_prefetchr(&SA[i - prefetch_distance]);

        sa_uint_t p0 = (sa_uint_t)SA[i - 0]; SAl[l] = p0 & SAINT_MAX; l -= (sa_sint_t)p0 < 0; SAr[r] = p0 - 1; r -= (sa_sint_t)p0 > 0;
        sa_uint_t p1 = (sa_uint_t)SA[i - 1]; SAl[l] = p1 & SAINT_MAX; l -= (sa_sint_t)p1 < 0; SAr[r] = p1 - 1; r -= (sa_sint_t)p1 > 0;
        sa_uint_t p2 = (sa_uint_t)SA[i - 2]; SAl[l] = p2 & SAINT_MAX; l -= (sa_sint_t)p2 < 0; SAr[r] = p2 - 1; r -= (sa_sint_t)p2 > 0;
        sa_uint_t p3 = (sa_uint_t)SA[i - 3]; SAl[l] = p3 & SAINT_MAX; l -= (sa_sint_t)p3 < 0; SAr[r] = p3 - 1; r -= (sa_sint_t)p3 > 0;
    }

    for (j -= 3; i >= j; i -= 1)
    {
        sa_uint_t p = (sa_uint_t)SA[i]; SAl[l] = p & SAINT_MAX; l -= (sa_sint_t)p < 0; SAr[r] = p - 1; r -= (sa_sint_t)p > 0;
    }
    
    *pl = l + 1; *pr = r + 1;
}

#if defined(LIBSAIS_OPENMP)

static sa_sint_t libsais_count_unique_suffixes(sa_sint_t * RESTRICT SA, sa_sint_t m, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    sa_sint_t * RESTRICT SAm = &SA[m];

    fast_sint_t i, j; sa_sint_t f0 = 0, f1 = 0, f2 = 0, f3 = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchr(&SAm[((sa_uint_t)SA[i + prefetch_distance + 0]) >> 1]);
        libsais_prefetchr(&SAm[((sa_uint_t)SA[i + prefetch_distance + 1]) >> 1]);
        libsais_prefetchr(&SAm[((sa_uint_t)SA[i + prefetch_distance + 2]) >> 1]);
        libsais_prefetchr(&SAm[((sa_uint_t)SA[i + prefetch_distance + 3]) >> 1]);

        f0 += SAm[((sa_uint_t)SA[i + 0]) >> 1] < 0;
        f1 += SAm[((sa_uint_t)SA[i + 1]) >> 1] < 0;
        f2 += SAm[((sa_uint_t)SA[i + 2]) >> 1] < 0;
        f3 += SAm[((sa_uint_t)SA[i + 3]) >> 1] < 0;
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        f0 += SAm[((sa_uint_t)SA[i]) >> 1] < 0;
    }

    return f0 + f1 + f2 + f3;
}

#endif

static sa_sint_t libsais_renumber_unique_and_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t f = 0;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;

        if (omp_num_threads == 1)
        {
            f = libsais_renumber_unique_and_nonunique_lms_suffixes_32s(T, SA, m, 0, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_count_unique_suffixes(SA, m, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, count = 0; for (t = 0; t < omp_thread_num; ++t) { count += thread_state[t].state.count; }

                if (omp_thread_num == omp_num_threads - 1)
                {
                    f = (sa_sint_t)(count + thread_state[omp_thread_num].state.count);
                }

                libsais_renumber_unique_and_nonunique_lms_suffixes_32s(T, SA, m, (sa_sint_t)count, omp_block_start, omp_block_size);
            }
        }
#endif
    }

    return f;
}

static void libsais_compact_unique_and_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t f, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 131072 && m < fs)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (((fast_sint_t)n >> 1) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : ((fast_sint_t)n >> 1) - omp_block_start;

        if (omp_num_threads == 1)
        {
            fast_sint_t l = m, r = (fast_sint_t)n + (fast_sint_t)fs;
            libsais_compact_unique_and_nonunique_lms_suffixes_32s(SA, m, &l, &r, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.position   = (fast_sint_t)m + ((fast_sint_t)n >> 1) + omp_block_start + omp_block_size;
                thread_state[omp_thread_num].state.count      = (fast_sint_t)m + omp_block_start + omp_block_size;

                libsais_compact_unique_and_nonunique_lms_suffixes_32s(SA, m, &thread_state[omp_thread_num].state.position, &thread_state[omp_thread_num].state.count, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, position;

                for (position = m, t = omp_num_threads - 1; t >= 0; --t)
                { 
                    fast_sint_t omp_block_end     = t < omp_num_threads - 1 ? omp_block_stride * (t + 1) : ((fast_sint_t)n >> 1);
                    fast_sint_t count             = ((fast_sint_t)m + ((fast_sint_t)n >> 1) + omp_block_end - thread_state[t].state.position);

                    if (count > 0)
                    {
                        position -= count; memcpy(&SA[position], &SA[thread_state[t].state.position], (size_t)count * sizeof(sa_sint_t));
                    }
                }

                for (position = (fast_sint_t)n + (fast_sint_t)fs, t = omp_num_threads - 1; t >= 0; --t)
                {
                    fast_sint_t omp_block_end     = t < omp_num_threads - 1 ? omp_block_stride * (t + 1) : ((fast_sint_t)n >> 1);
                    fast_sint_t count             = ((fast_sint_t)m + omp_block_end - thread_state[t].state.count);

                    if (count > 0)
                    {
                        position -= count; memcpy(&SA[position], &SA[thread_state[t].state.count], (size_t)count * sizeof(sa_sint_t));
                    }
                }
            }
        }
#endif
    }

    memcpy(&SA[(fast_sint_t)n + (fast_sint_t)fs - (fast_sint_t)m], &SA[(fast_sint_t)m - (fast_sint_t)f], (size_t)f * sizeof(sa_sint_t));
}

static sa_sint_t libsais_compact_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t f = libsais_renumber_unique_and_nonunique_lms_suffixes_32s_omp(T, SA, m, threads, thread_state);
    libsais_compact_unique_and_nonunique_lms_suffixes_32s_omp(SA, n, m, fs, f, threads, thread_state);

    return f;
}

#if (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_AVX2) && (SAINT_BIT == 32)

static fast_uint_t libsais_negative_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    return (fast_uint_t)(uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_loadu_si256((const __m256i *)(const void *)p)));
}

static fast_uint_t libsais_zero_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    return (fast_uint_t)(uint32_t)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpeq_epi32(_mm256_loadu_si256((const __m256i *)(const void *)p), _mm256_setzero_si256())));
}

static void libsais_clear_marks_8x32s(sa_sint_t * RESTRICT p)
{
    _mm256_storeu_si256((__m256i *)(void *)p, _mm256_and_si256(_mm256_loadu_si256((const __m256i *)(const void *)p), _mm256_set1_epi32(SAINT_MAX)));
}

#elif (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2) && (SAINT_BIT == 32)

static fast_uint_t libsais_negative_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    return (fast_uint_t)(uint32_t)(_mm_movemask_ps(_mm_castsi128_ps(_mm_loadu_si128((const __m128i *)(const void *)&p[0]))) | (_mm_movemask_ps(_mm_castsi128_ps(_mm_loadu_si128((const __m128i *)(const void *)&p[4]))) << 4));
}

static fast_uint_t libsais_zero_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    const __m128i zero = _mm_setzero_si128();

    return (fast_uint_t)(uint32_t)(_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_loadu_si128((const __m128i *)(const void *)&p[0]), zero))) | (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(_mm_loadu_si128((const __m128i *)(const void *)&p[4]), zero))) << 4));
}

static void libsais_clear_marks_8x32s(sa_sint_t * RESTRICT p)
{
    const __m128i max = _mm_set1_epi32(SAINT_MAX);

    _mm_storeu_si128((__m128i *)(void *)&p[0], _mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&p[0]), max));
    _mm_storeu_si128((__m128i *)(void *)&p[4], _mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&p[4]), max));
}

#elif (LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64) && (SAINT_BIT == 32)

static fast_uint_t libsais_lane_mask_8x32s(uint32x4_t lo, uint32x4_t hi)
{
    static const int32_t lanes[4] = { 0, 1, 2, 3 };

    int32x4_t shift = vld1q_s32(lanes);

    return (fast_uint_t)(vaddvq_u32(vshlq_u32(vshrq_n_u32(lo, 31), shift)) | (vaddvq_u32(vshlq_u32(vshrq_n_u32(hi, 31), shift)) << 4));
}

static fast_uint_t libsais_negative_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    return libsais_lane_mask_8x32s(vreinterpretq_u32_s32(vld1q_s32(&p[0])), vreinterpretq_u32_s32(vld1q_s32(&p[4])));
}

static fast_uint_t libsais_zero_mask_8x32s(const sa_sint_t * RESTRICT p)
{
    return libsais_lane_mask_8x32s(vceqzq_s32(vld1q_s32(&p[0])), vceqzq_s32(vld1q_s32(&p[4])));
}

static void libsais_clear_marks_8x32s(sa_sint_t * RESTRICT p)
{
    const int32x4_t max = vdupq_n_s32(SAINT_MAX);

    vst1q_s32(&p[0], vandq_s32(vld1q_s32(&p[0]), max));
    vst1q_s32(&p[4], vandq_s32(vld1q_s32(&p[4]), max));
}

#endif

static void libsais_merge_unique_lms_suffixes_32s(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT SAnm = &SA[(fast_sint_t)n - (fast_sint_t)m - 1 + l];

#if (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_NONE) && (SAINT_BIT == 32)
    {
        const sa_sint_t * RESTRICT SAnm8 = SAnm; fast_sint_t g, h;
        for (g = omp_block_start, h = omp_block_start + omp_block_size - 7; g < h; g += 8)
        {
            libsais_prefetchr(&T[g + prefetch_distance]);

            fast_uint_t mask = libsais_negative_mask_8x32s(&T[g]);
            if (mask != 0)
            {
                libsais_clear_marks_8x32s(&T[g]);
                do { SA[*SAnm8++] = (sa_sint_t)(g + libsais_ctz64(mask)); mask &= mask - 1; } while (mask != 0);
            }
        }

        SAnm = SAnm8; omp_block_size -= g - omp_block_start; omp_block_start = g;
    }
#endif

    sa_sint_t i, j; fast_sint_t tmp = *SAnm++;
    for (i = (sa_sint_t)omp_block_start, j = (sa_sint_t)omp_block_start + (sa_sint_t)omp_block_size - 6; i < j; i += 4)
    {
        libsais_prefetchr(&T[i + prefetch_distance]);

        sa_sint_t c0 = T[i + 0]; if (c0 < 0) { T[i + 0] = c0 & SAINT_MAX; SA[tmp] = i + 0; i++; tmp = *SAnm++; }
        sa_sint_t c1 = T[i + 1]; if (c1 < 0) { T[i + 1] = c1 & SAINT_MAX; SA[tmp] = i + 1; i++; tmp = *SAnm++; }
        sa_sint_t c2 = T[i + 2]; if (c2 < 0) { T[i + 2] = c2 & SAINT_MAX; SA[tmp] = i + 2; i++; tmp = *SAnm++; }
        sa_sint_t c3 = T[i + 3]; if (c3 < 0) { T[i + 3] = c3 & SAINT_MAX; SA[tmp] = i + 3; i++; tmp = *SAnm++; }
    }

    for (j += 6; i < j; i += 1)
    {
        sa_sint_t c = T[i]; if (c < 0) { T[i] = c & SAINT_MAX; SA[tmp] = i; i++; tmp = *SAnm++; }
    }
}

static void libsais_merge_nonunique_lms_suffixes_32s(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, fast_sint_t l, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    const sa_sint_t * RESTRICT SAnm = &SA[(fast_sint_t)n - (fast_sint_t)m - 1 + l];

#if (LIBSAIS_CPU_FEATURE != LIBSAIS_CPU_FEATURE_NONE) && (SAINT_BIT == 32)
    {
        const sa_sint_t * RESTRICT SAnm8 = SAnm; fast_sint_t g, h;
        for (g = omp_block_start, h = omp_block_start + omp_block_size - 7; g < h; g += 8)
        {
            libsais_prefetchr(&SA[g + prefetch_distance]);

            fast_uint_t mask = libsais_zero_mask_8x32s(&SA[g]);
            while (mask != 0) { SA[g + libsais_ctz64(mask)] = *SAnm8++; mask &= mask - 1; }
        }

        SAnm = SAnm8; omp_block_size -= g - omp_block_start; omp_block_start = g;
    }
#endif

    fast_sint_t i, j; sa_sint_t tmp = *SAnm++;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + prefetch_distance]);

        if (SA[i + 0] == 0) { SA[i + 0] = tmp; tmp = *SAnm++; }
        if (SA[i + 1] == 0) { SA[i + 1] = tmp; tmp = *SAnm++; }
        if (SA[i + 2] == 0) { SA[i + 2] = tmp; tmp = *SAnm++; }
        if (SA[i + 3] == 0) { SA[i + 3] = tmp; tmp = *SAnm++; }
    }

    for (j += 3; i < j; i += 1)
    {
        if (SA[i] == 0) { SA[i] = tmp; tmp = *SAnm++; }
    }
}

static void libsais_merge_unique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_merge_unique_lms_suffixes_32s(T, SA, n, m, 0, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_count_negative_marked_suffixes(T, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, count = 0; for (t = 0; t < omp_thread_num; ++t) { count += thread_state[t].state.count; }

                libsais_merge_unique_lms_suffixes_32s(T, SA, n, m, count, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

static void libsais_merge_nonunique_lms_suffixes_32s_omp(sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t f, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && m >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads); UNUSED(thread_state);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (m / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : m - omp_block_start;

        if (omp_num_threads == 1)
        {
            libsais_merge_nonunique_lms_suffixes_32s(SA, n, m, f, omp_block_start, omp_block_size);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            {
                thread_state[omp_thread_num].state.count = libsais_count_zero_marked_suffixes(SA, omp_block_start, omp_block_size);
            }

            #pragma omp barrier

            {
                fast_sint_t t, count = f; for (t = 0; t < omp_thread_num; ++t) { count += thread_state[t].state.count; }

                libsais_merge_nonunique_lms_suffixes_32s(SA, n, m, count, omp_block_start, omp_block_size);
            }
        }
#endif
    }
}

static void libsais_merge_compacted_lms_suffixes_32s_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t f, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    libsais_merge_unique_lms_suffixes_32s_omp(T, SA, n, m, threads, thread_state);
    libsais_merge_nonunique_lms_suffixes_32s_omp(SA, n, m, f, threads, thread_state);
}

static void libsais_reconstruct_compacted_lms_suffixes_32s_2k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t m, sa_sint_t fs, sa_sint_t f, sa_sint_t * RESTRICT buckets, sa_sint_t local_buckets, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (f > 0)
    {
#if defined(LIBSAIS_LANES_32S)
        memmove(&SA[m - f], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais_count_and_gather_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, local_buckets, threads, thread_state);
        memmove(&SA[n - m - 1], &SA[m - f], (size_t)f * sizeof(sa_sint_t));
#else
        memmove(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais_count_and_gather_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, local_buckets, threads, thread_state);
#endif
        libsais_reconstruct_lms_suffixes_omp(SA, n, m - f, threads);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(sa_sint_t));
        memset(&SA[0], 0, (size_t)m * sizeof(sa_sint_t));

        libsais_merge_compacted_lms_suffixes_32s_omp(T, SA, n, m, f, threads, thread_state);
    }
    else
    {
        libsais_count_and_gather_lms_suffixes_32s_2k(T, SA, n, k, buckets, 0, n);
        libsais_reconstruct_lms_suffixes_omp(SA, n, m, threads);
    }
}

static void libsais_reconstruct_compacted_lms_suffixes_32s_1k_omp(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t m, sa_sint_t fs, sa_sint_t f, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    if (f > 0)
    {
#if defined(LIBSAIS_LANES_32S)
        memmove(&SA[m - f], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais_gather_compacted_lms_suffixes_32s(T, SA, n);
        memmove(&SA[n - m - 1], &SA[m - f], (size_t)f * sizeof(sa_sint_t));
#else
        memmove(&SA[n - m - 1], &SA[n + fs - m], (size_t)f * sizeof(sa_sint_t));

        libsais_gather_compacted_lms_suffixes_32s(T, SA, n);
#endif
        libsais_reconstruct_lms_suffixes_omp(SA, n, m - f, threads);

        memcpy(&SA[n - m - 1 + f], &SA[0], ((size_t)m - (size_t)f) * sizeof(sa_sint_t));
        memset(&SA[0], 0, (size_t)m * sizeof(sa_sint_t));

        libsais_merge_compacted_lms_suffixes_32s_omp(T, SA, n, m, f, threads, thread_state);
    }
    else
    {
        libsais_gather_lms_suffixes_32s(T, SA, n);
        libsais_reconstruct_lms_suffixes_omp(SA, n, m, threads);
    }
}

static sa_sint_t libsais_main_32s_recursion(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t fs, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state, sa_sint_t * RESTRICT local_buffer)
{
    fs = fs < (SAINT_MAX - n) ? fs : (SAINT_MAX - n);

    if (k > 0 && ((fs / k >= 6) || (LIBSAIS_LOCAL_BUFFER_SIZE / k >= 6)))
    {
        sa_sint_t alignment = (fs - 1024) / k >= 6 ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = (fs - alignment) / k >= 6 ? (sa_sint_t *)libsais_align_up(&SA[n + fs - 6 * (fast_sint_t)k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : &SA[n + fs - 6 * (fast_sint_t)k];
        buckets = (LIBSAIS_LOCAL_BUFFER_SIZE > fs) ? local_buffer : buckets;

        sa_sint_t m = libsais_count_and_gather_lms_suffixes_32s_4k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads, thread_state);
        if (m > 1)
        {
            memset(SA, 0, ((size_t)n - (size_t)m) * sizeof(sa_sint_t));

            sa_sint_t first_lms_suffix    = SA[n - m];
            sa_sint_t left_suffixes_count = libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(T, k, buckets, first_lms_suffix);

            libsais_radix_sort_lms_suffixes_32s_6k_omp(T, SA, n, m, &buckets[4 * (fast_sint_t)k], threads, thread_state);

            if ((n / 8192) < k) { libsais_radix_sort_set_markers_32s_6k_omp(SA, k, &buckets[4 * (fast_sint_t)k], threads); }
            if (threads > 1 && n >= 65536) { memset(&SA[(fast_sint_t)n - (fast_sint_t)m], 0, (size_t)m * sizeof(sa_sint_t)); }

            libsais_initialize_buckets_for_partial_sorting_32s_6k(T, k, buckets, first_lms_suffix, left_suffixes_count);
            libsais_induce_partial_order_32s_6k_omp(T, SA, n, k, buckets, first_lms_suffix, left_suffixes_count, threads, thread_state);

            sa_sint_t names = (n / 8192) < k
                ? libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(SA, n, m, threads, thread_state)
                : libsais_renumber_and_gather_lms_suffixes_omp(SA, n, m, fs, threads, thread_state);

            if (names < m)
            {
                sa_sint_t f = (n / 8192) < k
                    ? libsais_compact_lms_suffixes_32s_omp(T, SA, n, m, fs, threads, thread_state)
                    : 0;

                if (libsais_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, threads, thread_state, local_buffer) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, m, fs, f, buckets, buckets == local_buffer, threads, thread_state);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads);
            }

            libsais_initialize_buckets_start_and_end_32s_4k(k, buckets);
            libsais_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
            libsais_induce_final_order_32s_4k(T, SA, n, k, buckets, threads, thread_state);
        }
        else
        {
            SA[0] = SA[n - 1];

            libsais_initialize_buckets_start_and_end_32s_6k(k, buckets);
            libsais_place_lms_suffixes_histogram_32s_6k(SA, n, k, m, buckets);
            libsais_induce_final_order_32s_6k(T, SA, n, k, buckets, threads, thread_state);
        }

        return 0;
    }
    else if (k > 0 && (n <= SAINT_MAX / 2) && ((fs / k >= 4) || (LIBSAIS_LOCAL_BUFFER_SIZE / k >= 4)))
    {
        sa_sint_t alignment = (fs - 1024) / k >= 4 ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = (fs - alignment) / k >= 4 ? (sa_sint_t *)libsais_align_up(&SA[n + fs - 4 * (fast_sint_t)k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : &SA[n + fs - 4 * (fast_sint_t)k];
        buckets = (LIBSAIS_LOCAL_BUFFER_SIZE > fs) ? local_buffer : buckets;

        sa_sint_t m = libsais_count_and_gather_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads, thread_state);
        if (m > 1)
        {
            libsais_initialize_buckets_for_radix_and_partial_sorting_32s_4k(T, k, buckets, SA[n - m]);

            libsais_radix_sort_lms_suffixes_32s_2k_omp(T, SA, n, m, &buckets[1], threads, thread_state);
            libsais_radix_sort_set_markers_32s_4k_omp(SA, k, &buckets[1], threads);
            
            libsais_place_lms_suffixes_interval_32s_4k(SA, n, k, m - 1, buckets);
            libsais_induce_partial_order_32s_4k_omp(T, SA, n, k, buckets, threads, thread_state);

            sa_sint_t names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(SA, n, m, threads, thread_state);
            if (names < m)
            {
                sa_sint_t f = libsais_compact_lms_suffixes_32s_omp(T, SA, n, m, fs, threads, thread_state);

                if (libsais_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, threads, thread_state, local_buffer) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, m, fs, f, buckets, buckets == local_buffer, threads, thread_state);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads);
            }
        }
        else
        {
            SA[0] = SA[n - 1];
        }

        libsais_initialize_buckets_start_and_end_32s_4k(k, buckets);
        libsais_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
        libsais_induce_final_order_32s_4k(T, SA, n, k, buckets, threads, thread_state);

        return 0;
    }
    else if (k > 0 && ((fs / k >= 2) || (LIBSAIS_LOCAL_BUFFER_SIZE / k >= 2)))
    {
        sa_sint_t alignment = (fs - 1024) / k >= 2 ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = (fs - alignment) / k >= 2 ? (sa_sint_t *)libsais_align_up(&SA[n + fs - 2 * (fast_sint_t)k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : &SA[n + fs - 2 * (fast_sint_t)k];
        buckets = (LIBSAIS_LOCAL_BUFFER_SIZE > fs) ? local_buffer : buckets;

        sa_sint_t m = libsais_count_and_gather_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads, thread_state);
        if (m > 1)
        {
            libsais_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(T, k, buckets, SA[n - m]);

            libsais_radix_sort_lms_suffixes_32s_2k_omp(T, SA, n, m, &buckets[1], threads, thread_state);
            libsais_place_lms_suffixes_interval_32s_2k(SA, n, k, m - 1, buckets);

            libsais_initialize_buckets_start_and_end_32s_2k(k, buckets);
            libsais_induce_partial_order_32s_2k_omp(T, SA, n, k, buckets, threads, thread_state);

            sa_sint_t names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(T, SA, n, m, threads, thread_state);
            if (names < m)
            {
                sa_sint_t f = libsais_compact_lms_suffixes_32s_omp(T, SA, n, m, fs, threads, thread_state);

                if (libsais_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, threads, thread_state, local_buffer) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_2k_omp(T, SA, n, k, m, fs, f, buckets, buckets == local_buffer, threads, thread_state);
            }
            else
            {
                libsais_count_lms_suffixes_32s_2k_omp(T, SA, n, k, buckets, buckets == local_buffer, threads);
            }
        }
        else
        {
            SA[0] = SA[n - 1];
        }

        libsais_initialize_buckets_end_32s_2k(k, buckets);
        libsais_place_lms_suffixes_histogram_32s_2k(SA, n, k, m, buckets);

        libsais_initialize_buckets_start_and_end_32s_2k(k, buckets);
        libsais_induce_final_order_32s_2k(T, SA, n, k, buckets, threads, thread_state);

        return 0;
    }
    else
    {
        sa_sint_t * buffer = fs < k ? (sa_sint_t *)libsais_alloc_aligned((size_t)k * sizeof(sa_sint_t), 4096) : (sa_sint_t *)NULL;

        sa_sint_t alignment = fs - 1024 >= k ? (sa_sint_t)1024 : (sa_sint_t)16;
        sa_sint_t * RESTRICT buckets = fs - alignment >= k ? (sa_sint_t *)libsais_align_up(&SA[n + fs - k - alignment], (size_t)alignment * sizeof(sa_sint_t)) : fs >= k ? &SA[n + fs - k] : buffer;

        if (buckets == NULL) { return -2; }

        memset(SA, 0, (size_t)n * sizeof(sa_sint_t));

        libsais_count_suffixes_32s(T, n, k, buckets); 
        libsais_initialize_buckets_end_32s_1k(k, buckets);

        sa_sint_t m = libsais_radix_sort_lms_suffixes_32s_1k(T, SA, n, buckets);
        if (m > 1)
        {
            libsais_induce_partial_order_32s_1k_omp(T, SA, n, k, buckets, threads, thread_state);

            sa_sint_t names = libsais_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(T, SA, n, m, threads, thread_state);
            if (names < m)
            {
                if (buffer != NULL) { libsais_free_aligned(buffer); buckets = NULL; }

                sa_sint_t f = libsais_compact_lms_suffixes_32s_omp(T, SA, n, m, fs, threads, thread_state);

                if (libsais_main_32s_recursion(SA + n + fs - m + f, SA, m - f, names - f, fs + n - 2 * m + f, threads, thread_state, local_buffer) != 0)
                {
                    return -2;
                }

                libsais_reconstruct_compacted_lms_suffixes_32s_1k_omp(T, SA, n, m, fs, f, threads, thread_state);

                if (buckets == NULL) { buckets = buffer = (sa_sint_t *)libsais_alloc_aligned((size_t)k * sizeof(sa_sint_t), 4096); }
                if (buckets == NULL) { return -2; }
            }
            
            libsais_count_suffixes_32s(T, n, k, buckets);
            libsais_initialize_buckets_end_32s_1k(k, buckets);
            libsais_place_lms_suffixes_interval_32s_1k(T, SA, k, m, buckets);
        }

        libsais_induce_final_order_32s_1k(T, SA, n, k, buckets, threads, thread_state);
        libsais_free_aligned(buffer);

        return 0;
    }
}

static sa_sint_t libsais_main_32s_entry(sa_sint_t * RESTRICT T, sa_sint_t * RESTRICT SA, sa_sint_t n, sa_sint_t k, sa_sint_t fs, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    sa_sint_t local_buffer[2 * LIBSAIS_LOCAL_BUFFER_SIZE];

    return libsais_main_32s_recursion(T, SA, n, k, fs, threads, thread_state, local_buffer + LIBSAIS_LOCAL_BUFFER_SIZE);
}

static sa_sint_t libsais_main_8u(const uint8_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t * RESTRICT buckets, sa_sint_t flags, sa_sint_t r, sa_sint_t * RESTRICT I, sa_sint_t fs, sa_sint_t * freq, sa_sint_t threads, LIBSAIS_THREAD_STATE * RESTRICT thread_state)
{
    fs = fs < (SAINT_MAX - n) ? fs : (SAINT_MAX - n);

    sa_sint_t lms_threads = threads;
    sa_sint_t m = libsais_count_and_gather_lms_suffixes_8u_omp(T, SA, n, buckets, &lms_threads, thread_state);
    sa_sint_t k = libsais_initialize_buckets_start_and_end_8u(buckets, freq);

    if ((flags & LIBSAIS_FLAGS_GSA) && (buckets[0] != 0 || buckets[2] != 0 || buckets[3] != 1))
    {
        return -1;
    }

    if (m > 0)
    {
        sa_sint_t first_lms_suffix    = SA[n - m];
        sa_sint_t left_suffixes_count = libsais_initialize_buckets_for_lms_suffixes_radix_sort_8u(T, buckets, first_lms_suffix);
        sa_sint_t middle_m            = n >= 65536 ? libsais_count_middle_lms_suffixes(SA, n, m) : 0;
        sa_sint_t middle_lms_suffix   = middle_m > 0 ? SA[n - middle_m] : 0;

        if (threads > 1 && n >= 65536) { memset(SA, 0, ((size_t)n - (size_t)m) * sizeof(sa_sint_t)); }
        libsais_radix_sort_lms_suffixes_8u_omp(T, SA, n, m, flags, buckets, lms_threads, thread_state);
        if (threads > 1 && n >= 65536) { memset(&SA[(fast_sint_t)n - (fast_sint_t)m], 0, (size_t)m * sizeof(sa_sint_t)); }

        libsais_initialize_buckets_for_partial_sorting_8u(T, buckets, first_lms_suffix, left_suffixes_count);
        libsais_induce_partial_order_8u_omp(T, SA, n, k, flags, buckets, first_lms_suffix, left_suffixes_count, threads, thread_state);

        sa_sint_t names = libsais_renumber_and_gather_lms_suffixes_omp(SA, n, m, fs, threads, thread_state);
        if (names < m)
        {
            if (libsais_main_32s_entry(SA + n + fs - m, SA, m, names, fs + n - 2 * m, threads, thread_state) != 0)
            {
                return -2;
            }

            libsais_gather_lms_suffixes_8u_omp(T, SA, n, middle_lms_suffix, middle_m, lms_threads, thread_state);
            libsais_reconstruct_lms_suffixes_omp(SA, n, m, threads);
        }

        libsais_place_lms_suffixes_interval_8u(SA, n, m, flags, buckets);
    }
    else
    {
        memset(SA, 0, (size_t)n * sizeof(sa_sint_t));
    }

    return libsais_induce_final_order_8u_omp(T, SA, n, k, flags, r, I, buckets, threads, thread_state);
}

static sa_sint_t libsais_main(const uint8_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t flags, sa_sint_t r, sa_sint_t * I, sa_sint_t fs, sa_sint_t * freq, sa_sint_t threads)
{
    threads = n >= 65536 ? threads : 1;

    LIBSAIS_THREAD_STATE *  RESTRICT thread_state   = threads > 1 ? libsais_alloc_thread_state(threads) : NULL;
    sa_sint_t *             RESTRICT buckets        = (sa_sint_t *)libsais_alloc_aligned((size_t)8 * ALPHABET_SIZE * sizeof(sa_sint_t), 4096);

    sa_sint_t index = buckets != NULL && (thread_state != NULL || threads == 1)
        ? libsais_main_8u(T, SA, n, buckets, flags, r, I, fs, freq, threads, thread_state)
        : -2;

    libsais_free_aligned(buckets);
    libsais_free_thread_state(thread_state);

    return index;
}

static sa_sint_t libsais_main_int(sa_sint_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t k, sa_sint_t fs, sa_sint_t threads)
{
    threads = n >= 65536 ? threads : 1;

    LIBSAIS_THREAD_STATE * RESTRICT thread_state = threads > 1 ? libsais_alloc_thread_state(threads) : NULL;

    sa_sint_t index = thread_state != NULL || threads == 1
        ? libsais_main_32s_entry(T, SA, n, k, fs, threads, thread_state)
        : -2;

    libsais_free_thread_state(thread_state);

    return index;
}

static sa_sint_t libsais_main_ctx(const LIBSAIS_CONTEXT * ctx, const uint8_t * T, sa_sint_t * SA, sa_sint_t n, sa_sint_t flags, sa_sint_t r, sa_sint_t * I, sa_sint_t fs, sa_sint_t * freq)
{
    return ctx != NULL && (ctx->buckets != NULL && (ctx->thread_state != NULL || ctx->threads == 1))
        ? libsais_main_8u(T, SA, n, ctx->buckets, flags, r, I, fs, freq, (sa_sint_t)ctx->threads, ctx->thread_state)
        : -2;
}

static void libsais_bwt_copy_8u(uint8_t * RESTRICT U, sa_sint_t * RESTRICT A, sa_sint_t n)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = 0, j = (fast_sint_t)n - 7; i < j; i += 8)
    {
        libsais_prefetchr(&A[i + prefetch_distance]);

        U[i + 0] = (uint8_t)A[i + 0];
        U[i + 1] = (uint8_t)A[i + 1];
        U[i + 2] = (uint8_t)A[i + 2];
        U[i + 3] = (uint8_t)A[i + 3];
        U[i + 4] = (uint8_t)A[i + 4];
        U[i + 5] = (uint8_t)A[i + 5];
        U[i + 6] = (uint8_t)A[i + 6];
        U[i + 7] = (uint8_t)A[i + 7];
    }

    for (j += 7; i < j; i += 1)
    {
        U[i] = (uint8_t)A[i];
    }
}

static void libsais_bwt_copy_8u_omp(uint8_t * RESTRICT U, sa_sint_t * RESTRICT A, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
        fast_sint_t omp_block_stride  = ((fast_sint_t)n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : (fast_sint_t)n - omp_block_start;
#else
        UNUSED(threads);

        fast_sint_t omp_block_start   = 0;
        fast_sint_t omp_block_size    = (fast_sint_t)n;
#endif

        libsais_bwt_copy_8u(U + omp_block_start, A + omp_block_start, (sa_sint_t)omp_block_size);
    }
}

static fast_uint_t libsais_minimum_symbol_8u(const uint8_t * RESTRICT T, fast_sint_t n)
{
    uint8_t c = ALPHABET_SIZE - 1; fast_sint_t i; for (i = 0; i < n; ++i) { if (T[i] < c) { c = T[i]; } } return c;
}

static fast_uint_t libsais_minimum_symbol_8u_omp(const uint8_t * RESTRICT T, sa_sint_t n, sa_sint_t threads)
{
    fast_uint_t c = libsais_minimum_symbol_8u(T, n < 4096 ? n : 4096);

#if defined(LIBSAIS_OPENMP)
    if (c > 0 && threads > 1 && n >= 65536)
    {
        #pragma omp parallel num_threads(threads)
        {
            fast_sint_t omp_thread_num    = omp_get_thread_num();
            fast_sint_t omp_num_threads   = omp_get_num_threads();
            fast_sint_t omp_block_stride  = ((fast_sint_t)n / omp_num_threads) & (-16);
            fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
            fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : (fast_sint_t)n - omp_block_start;

            fast_uint_t block_c = libsais_minimum_symbol_8u(T + omp_block_start, omp_block_size);

            #pragma omp critical(libsais_minimum_symbol)
            {
                c = block_c < c ? block_c : c;
            }
        }
    }
    else
#else
    UNUSED(threads);
#endif
    if (c > 0 && n > 4096)
    {
        c = libsais_minimum_symbol_8u(T, n);
    }

    return c;
}

static fast_sint_t libsais_find_minimum_symbol_8u(const uint8_t * RESTRICT T, fast_sint_t i, fast_sint_t n, fast_uint_t c)
{
    const uint8_t * p = (const uint8_t *)memchr(&T[i], (int)c, (size_t)(n - i));

    return p != NULL ? (fast_sint_t)(p - T) : n;
}

static fast_sint_t libsais_gather_longest_runs_8u(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT A, fast_uint_t c, fast_sint_t * length, fast_sint_t count, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    fast_sint_t i, j, l = *length;
    for (i = omp_block_start + l, j = omp_block_start + omp_block_size; i < j; )
    {
        if (T[i] != c) { i = l > 4 ? i + l : libsais_find_minimum_symbol_8u(T, i, j, c); continue; }

        fast_sint_t run_start = i; while (T[run_start - 1] == c) { run_start--; }
        fast_sint_t run_end   = i; while (T[run_end] == c) { run_end++; }

        if (run_end - run_start > l) { l = run_end - run_start; count = 0; }
        if (run_end - run_start == l) { A[count++] = (sa_sint_t)run_start; }

        i = run_end + l;
    }

    *length = l;

    return count;
}

static fast_sint_t libsais_least_rotation_8u(const uint8_t * RESTRICT T, const sa_sint_t * RESTRICT A, sa_sint_t n, fast_sint_t length, fast_sint_t count)
{
    fast_sint_t i, j;
    for (i = 0, j = 1; i < count && j < count; j += (i == j))
    {
        fast_sint_t l = length, p = A[i] + l, q = A[j] + l, c0 = 0, c1 = 0;
        while (p < n && q < n && (c0 = T[p]) == (c1 = T[q])) { l++; p++; q++; }
        while (l < n && (c0 = T[p - (p >= n ? n : 0)]) == (c1 = T[q - (q >= n ? n : 0)])) { l++; p++; q++; }

        if (l == n) { break; }

        if (c0 > c1) { while (i < count && A[i] <= p) { i++; } }
        else         { while (j < count && A[j] <= q) { j++; } }
    }

    return A[i < j ? i : j];
}

static sa_sint_t libsais_least_rotation_8u_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT A, sa_sint_t n, sa_sint_t threads)
{
    fast_uint_t c = libsais_minimum_symbol_8u_omp(T, n, threads);

    fast_sint_t head = 0; while (head < n && T[head] == c) { head++; }
    if (head == n) { return 0; }

    fast_sint_t tail = 0; while (T[n - 1 - tail] == c) { tail++; }

    fast_sint_t length = head + tail, start = -1;

    threads = threads < n ? threads : n;

#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (((fast_sint_t)n - head - tail) / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = head + omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_end     = omp_thread_num < omp_num_threads - 1 ? omp_block_start + omp_block_stride : (fast_sint_t)n - tail;

        while (T[omp_block_start] == c) { omp_block_start++; }
        while (omp_block_end < (fast_sint_t)n - tail && T[omp_block_end] == c) { omp_block_end++; }

        sa_sint_t * RESTRICT block_runs = &A[omp_block_start]; fast_sint_t block_length = head + tail, block_count = 0;

        if (omp_thread_num == 0 && head > 0 && tail == 0) { block_runs -= 1; block_runs[block_count++] = 0; }

        block_count = libsais_gather_longest_runs_8u(T, block_runs, c, &block_length, block_count, omp_block_start, omp_block_end - omp_block_start);

        if (omp_thread_num == omp_num_threads - 1 && tail > 0 && head + tail == block_length) { block_runs[block_count++] = (sa_sint_t)(n - tail); }

        if (omp_num_threads == 1)
        {
            length = block_length; start = libsais_least_rotation_8u(T, block_runs, n, length, block_count);
        }
#if defined(LIBSAIS_OPENMP)
        else
        {
            #pragma omp critical(libsais_least_rotation)
            {
                if (block_count > 0 && block_length > length) { length = block_length; }
            }

            #pragma omp barrier

            fast_sint_t block_start = block_count > 0 && block_length == length ? libsais_least_rotation_8u(T, block_runs, n, length, block_count) : -1;

            #pragma omp barrier

            A[omp_thread_num] = (sa_sint_t)block_start;

            #pragma omp barrier

            #pragma omp master
            {
                fast_sint_t t, count = 0; for (t = 0; t < omp_num_threads; ++t) { if (A[t] >= 0) { A[count++] = A[t]; } }

                start = libsais_least_rotation_8u(T, A, n, length, count);
            }
        }
#endif
    }

    return (sa_sint_t)(start > 0 ? n - start : 0);
}

static void libsais_bwt_cyclic_rotate_8u(const uint8_t * T, uint8_t * U, sa_sint_t * RESTRICT A, sa_sint_t n, sa_sint_t origin)
{
    if (T != U)
    {
        memcpy(U, &T[n - origin], (size_t)origin * sizeof(uint8_t)); memcpy(&U[origin], T, (size_t)(n - origin) * sizeof(uint8_t));
    }
    else if (origin <= n - origin)
    {
        memcpy(A, &U[n - origin], (size_t)origin * sizeof(uint8_t)); memmove(&U[origin], U, (size_t)(n - origin) * sizeof(uint8_t)); memcpy(U, A, (size_t)origin * sizeof(uint8_t));
    }
    else
    {
        memcpy(A, U, (size_t)(n - origin) * sizeof(uint8_t)); memmove(U, &U[n - origin], (size_t)origin * sizeof(uint8_t)); memcpy(&U[origin], A, (size_t)(n - origin) * sizeof(uint8_t));
    }
}

void * libsais_create_ctx(void)
{
    return (void *)libsais_create_ctx_main(1);
}

void libsais_free_ctx(void * ctx)
{
    libsais_free_ctx_main((LIBSAIS_CONTEXT *)ctx);
}

int32_t libsais(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    return libsais_main(T, SA, n, LIBSAIS_FLAGS_NONE, 0, NULL, fs, freq, 1);
}

int32_t libsais_gsa(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (n > 0 && T[n - 1] != 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    return libsais_main(T, SA, n, LIBSAIS_FLAGS_GSA, 0, NULL, fs, freq, 1);
}

int32_t libsais_int(int32_t * T, int32_t * SA, int32_t n, int32_t k, int32_t fs)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (n > 1 && k <= 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { SA[0] = 0; }
        return 0;
    }

    return libsais_main_int(T, SA, n, k, fs, 1);
}

int32_t libsais_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq)
{
    if ((ctx == NULL) || (T == NULL) || (SA == NULL) || (n < 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    return libsais_main_ctx((const LIBSAIS_CONTEXT *)ctx, T, SA, n, LIBSAIS_FLAGS_NONE, 0, NULL, fs, freq);
}

int32_t libsais_gsa_ctx(const void * ctx, const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq)
{
    if ((ctx == NULL) || (T == NULL) || (SA == NULL) || (n < 0) || (n > 0 && T[n - 1] != 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    return libsais_main_ctx((const LIBSAIS_CONTEXT *)ctx, T, SA, n, LIBSAIS_FLAGS_GSA, 0, NULL, fs, freq);
}

int32_t libsais_bwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0))
    { 
        return -1; 
    }
    else if (n <= 1) 
    { 
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return n;
    }

    sa_sint_t index = libsais_main(T, A, n, LIBSAIS_FLAGS_BWT, 0, NULL, fs, freq, 1);
    if (index >= 0) 
    { 
        index++;

        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, index - 1, 1);
        libsais_bwt_copy_8u_omp(U + index, A + index, n - index, 1);
    }

    return index;
}

int32_t libsais_bwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq, int32_t r, int32_t * I)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0) || (r < 2) || ((r & (r - 1)) != 0) || (I == NULL))
    { 
        return -1; 
    }
    else if (n <= 1) 
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        I[0] = n;
        return 0;
    }

    sa_sint_t index = libsais_main(T, A, n, LIBSAIS_FLAGS_BWT, r, I, fs, freq, 1);
    if (index == 0) 
    { 
        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, I[0] - 1, 1);
        libsais_bwt_copy_8u_omp(U + I[0], A + I[0], n - I[0], 1);
    }

    return index;
}

int32_t libsais_bwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    sa_sint_t origin = libsais_least_rotation_8u_omp(T, A, n, 1);
    if (origin == 0)
    {
        sa_sint_t index = libsais_bwt(T, U, A, n, fs, freq);
        return index >= 0 ? index - 1 : index;
    }

    libsais_bwt_cyclic_rotate_8u(T, U, A, n, origin);

    sa_sint_t index = libsais_main(U, A, n, LIBSAIS_FLAGS_BWT | LIBSAIS_FLAGS_CYCLIC, origin, NULL, fs, freq, 1);
    if (index >= 0)
    {
        libsais_bwt_copy_8u_omp(U, A, n, 1);
    }
    else if (T == U)
    {
        libsais_bwt_cyclic_rotate_8u(U, U, A, n, n - origin);
    }

    return index;
}

int32_t libsais_bwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq)
{
    if ((ctx == NULL) || (T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0))
    { 
        return -1; 
    }
    else if (n <= 1) 
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return n;
    }

    sa_sint_t index = libsais_main_ctx((const LIBSAIS_CONTEXT *)ctx, T, A, n, LIBSAIS_FLAGS_BWT, 0, NULL, fs, freq);
    if (index >= 0) 
    { 
        index++;

        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, index - 1, (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
        libsais_bwt_copy_8u_omp(U + index, A + index, n - index, (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
    }

    return index;
}

int32_t libsais_bwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq, int32_t r, int32_t * I)
{
    if ((ctx == NULL) || (T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0) || (r < 2) || ((r & (r - 1)) != 0) || (I == NULL))
    { 
        return -1; 
    }
    else if (n <= 1) 
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        I[0] = n;
        return 0;
    }

    sa_sint_t index = libsais_main_ctx((const LIBSAIS_CONTEXT *)ctx, T, A, n, LIBSAIS_FLAGS_BWT, r, I, fs, freq);
    if (index == 0) 
    { 
        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, I[0] - 1, (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
        libsais_bwt_copy_8u_omp(U + I[0], A + I[0], n - I[0], (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
    }

    return index;
}

int32_t libsais_bwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq)
{
    if ((ctx == NULL) || (T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    sa_sint_t origin = libsais_least_rotation_8u_omp(T, A, n, (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
    if (origin == 0)
    {
        sa_sint_t index = libsais_bwt_ctx(ctx, T, U, A, n, fs, freq);
        return index >= 0 ? index - 1 : index;
    }

    libsais_bwt_cyclic_rotate_8u(T, U, A, n, origin);

    sa_sint_t index = libsais_main_ctx((const LIBSAIS_CONTEXT *)ctx, U, A, n, LIBSAIS_FLAGS_BWT | LIBSAIS_FLAGS_CYCLIC, origin, NULL, fs, freq);
    if (index >= 0)
    {
        libsais_bwt_copy_8u_omp(U, A, n, (sa_sint_t)((const LIBSAIS_CONTEXT *)ctx)->threads);
    }
    else if (T == U)
    {
        libsais_bwt_cyclic_rotate_8u(U, U, A, n, n - origin);
    }

    return index;
}

#if defined(LIBSAIS_OPENMP)

void * libsais_create_ctx_omp(int32_t threads)
{
    if (threads < 0) { return NULL; }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return (void *)libsais_create_ctx_main(threads);
}

int32_t libsais_omp(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (fs < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return libsais_main(T, SA, n, LIBSAIS_FLAGS_NONE, 0, NULL, fs, freq, threads);
}

int32_t libsais_gsa_omp(const uint8_t * T, int32_t * SA, int32_t n, int32_t fs, int32_t * freq, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (n > 0 && T[n - 1] != 0) || (fs < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { SA[0] = 0; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return libsais_main(T, SA, n, LIBSAIS_FLAGS_GSA, 0, NULL, fs, freq, threads);
}

int32_t libsais_int_omp(int32_t * T, int32_t * SA, int32_t n, int32_t k, int32_t fs, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (n < 0) || (n > 1 && k <= 0) || (fs < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { SA[0] = 0; }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return libsais_main_int(T, SA, n, k, fs, threads);
}

int32_t libsais_bwt_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq, int32_t threads)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return n;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    sa_sint_t index = libsais_main(T, A, n, LIBSAIS_FLAGS_BWT, 0, NULL, fs, freq, threads);
    if (index >= 0)
    {
        index++;

        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, index - 1, threads);
        libsais_bwt_copy_8u_omp(U + index, A + index, n - index, threads);
    }

    return index;
}

int32_t libsais_bwt_aux_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq, int32_t r, int32_t * I, int32_t threads)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0) || (r < 2) || ((r & (r - 1)) != 0) || (I == NULL) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        I[0] = n;
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    sa_sint_t index = libsais_main(T, A, n, LIBSAIS_FLAGS_BWT, r, I, fs, freq, threads);
    if (index == 0) 
    { 
        U[0] = T[n - 1];
        libsais_bwt_copy_8u_omp(U + 1, A, I[0] - 1, threads);
        libsais_bwt_copy_8u_omp(U + I[0], A + I[0], n - I[0], threads);
    }

    return index;
}

int32_t libsais_bwt_cyclic_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, int32_t fs, int32_t * freq, int32_t threads)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (fs < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (freq != NULL) { memset(freq, 0, ALPHABET_SIZE * sizeof(int32_t)); }
        if (n == 1) { U[0] = T[0]; if (freq != NULL) { freq[T[0]]++; } }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    sa_sint_t origin = libsais_least_rotation_8u_omp(T, A, n, threads);
    if (origin == 0)
    {
        sa_sint_t index = libsais_bwt_omp(T, U, A, n, fs, freq, threads);
        return index >= 0 ? index - 1 : index;
    }

    libsais_bwt_cyclic_rotate_8u(T, U, A, n, origin);

    sa_sint_t index = libsais_main(U, A, n, LIBSAIS_FLAGS_BWT | LIBSAIS_FLAGS_CYCLIC, origin, NULL, fs, freq, threads);
    if (index >= 0)
    {
        libsais_bwt_copy_8u_omp(U, A, n, threads);
    }
    else if (T == U)
    {
        libsais_bwt_cyclic_rotate_8u(U, U, A, n, n - origin);
    }

    return index;
}

#endif

static sa_sint_t libsais_unbwt_chains_threads(sa_sint_t n, sa_sint_t threads)
{
    sa_sint_t chains_threads = n / 262144 > 8 ? n / 262144 : 8; if (chains_threads > threads) { chains_threads = threads; } if (chains_threads > 64) { chains_threads = 64; }
    if (n < UNBWT_CHAIN_THREADS_MIN_SIZE) { chains_threads = 1; }

    return chains_threads;
}

static fast_uint_t libsais_unbwt_chunk_size(sa_sint_t n)
{
    return n < UNBWT_CHAIN_THREADS_MIN_SIZE ? UNBWT_CHUNK_SIZE_SMALL : UNBWT_CHUNK_SIZE;
}

static fast_uint_t libsais_unbwt_chains_slots(sa_sint_t n, sa_sint_t chains_threads)
{
    return n >= UNBWT_CHAINS_LARGE_MIN_SIZE && chains_threads <= 4 ? UNBWT_CHAINS_LARGE : UNBWT_CHAINS;
}

static fast_uint_t libsais_unbwt_chain_marks(sa_sint_t n, sa_sint_t chains_threads)
{
    fast_uint_t chain_marks = (fast_uint_t)32 * libsais_unbwt_chains_slots(n, chains_threads) * (fast_uint_t)chains_threads;
    fast_uint_t min_marks   = (fast_uint_t)n >> 8; if (min_marks > UNBWT_CHAIN_MARKS) { min_marks = UNBWT_CHAIN_MARKS; }

    return chain_marks > min_marks ? chain_marks : min_marks;
}

static LIBSAIS_UNBWT_CONTEXT * libsais_unbwt_create_ctx_main(sa_sint_t threads)
{
    if ((size_t)threads > SIZE_MAX / ((2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE)) * sizeof(sa_uint_t))) { return NULL; }

    sa_sint_t   chains_threads  = libsais_unbwt_chains_threads(SAINT_MAX, threads);
    fast_uint_t chain_marks     = libsais_unbwt_chain_marks(SAINT_MAX, chains_threads);

    LIBSAIS_UNBWT_CONTEXT *     RESTRICT ctx            = (LIBSAIS_UNBWT_CONTEXT *)libsais_alloc_aligned(sizeof(LIBSAIS_UNBWT_CONTEXT), 64);
    sa_uint_t *                 RESTRICT bucket2        = (sa_uint_t *)libsais_alloc_aligned(ALPHABET_SIZE * ALPHABET_SIZE * sizeof(sa_uint_t), 4096);
    uint16_t *                  RESTRICT fastbits       = (uint16_t *)libsais_alloc_aligned((1 + (1 << UNBWT_FASTBITS)) * sizeof(uint16_t), 4096);
    sa_uint_t *                 RESTRICT buckets        = threads > 1 ? (sa_uint_t *)libsais_alloc_aligned((size_t)threads * (2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE)) * sizeof(sa_uint_t), 4096) : NULL;
    sa_uint_t *                 RESTRICT marks          = (sa_uint_t *)libsais_alloc_aligned((2 * (size_t)chain_marks + 2) * sizeof(sa_uint_t), 4096);
    LIBSAIS_UNBWT_CHAIN *       RESTRICT chains         = (LIBSAIS_UNBWT_CHAIN *)libsais_alloc_aligned((2 * (size_t)chain_marks + 2) * sizeof(LIBSAIS_UNBWT_CHAIN), 4096);
    uint8_t *                   RESTRICT slack          = (uint8_t *)libsais_alloc_aligned((size_t)libsais_unbwt_chains_slots(SAINT_MAX, chains_threads) * (size_t)chains_threads * UNBWT_CHUNK_SIZE, 4096);
    sa_uint_t *                 RESTRICT links          = (sa_uint_t *)libsais_alloc_aligned(((size_t)SAINT_MAX / UNBWT_CHUNK_SIZE + (size_t)libsais_unbwt_chains_slots(SAINT_MAX, chains_threads) * (size_t)chains_threads) * sizeof(sa_uint_t), 4096);

    if (marks == NULL || chains == NULL || slack == NULL || links == NULL)
    {
        libsais_free_aligned(links); libsais_free_aligned(slack); libsais_free_aligned(chains); libsais_free_aligned(marks);
        marks = NULL; chains = NULL; slack = NULL; links = NULL;
    }

    if (ctx != NULL && bucket2 != NULL && fastbits != NULL && (buckets != NULL || threads == 1))
    {
        ctx->bucket2    = bucket2;
        ctx->fastbits   = fastbits;
        ctx->buckets    = buckets;
        ctx->marks      = marks;
        ctx->chains     = chains;
        ctx->slack      = slack;
        ctx->links      = links;
        ctx->threads    = threads;

        return ctx;
    }

    libsais_free_aligned(links);
    libsais_free_aligned(slack);
    libsais_free_aligned(chains);
    libsais_free_aligned(marks);
    libsais_free_aligned(buckets);
    libsais_free_aligned(fastbits);
    libsais_free_aligned(bucket2);
    libsais_free_aligned(ctx);

    return NULL;
}

static void libsais_unbwt_free_ctx_main(LIBSAIS_UNBWT_CONTEXT * ctx)
{
    if (ctx != NULL)
    {
        libsais_free_aligned(ctx->links);
        libsais_free_aligned(ctx->slack);
        libsais_free_aligned(ctx->chains);
        libsais_free_aligned(ctx->marks);
        libsais_free_aligned(ctx->buckets);
        libsais_free_aligned(ctx->fastbits);
        libsais_free_aligned(ctx->bucket2);
        libsais_free_aligned(ctx);
    }
}

#if (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_BLOCK_LANES (64)

static fast_uint_t libsais_unbwt_uniform_block(const uint8_t * RESTRICT T)
{
    __m128i v0 = _mm_load_si128((const __m128i *)(const void *)&T[0]), v1 = _mm_load_si128((const __m128i *)(const void *)&T[16]);
    __m128i v2 = _mm_load_si128((const __m128i *)(const void *)&T[32]), v3 = _mm_load_si128((const __m128i *)(const void *)&T[48]);
    __m128i b  = _mm_shuffle_epi32(_mm_shufflelo_epi16(_mm_unpacklo_epi8(v0, v0), 0), 0);

    return (fast_uint_t)(_mm_movemask_epi8(_mm_and_si128(_mm_and_si128(_mm_cmpeq_epi8(v0, b), _mm_cmpeq_epi8(v1, b)), _mm_and_si128(_mm_cmpeq_epi8(v2, b), _mm_cmpeq_epi8(v3, b)))) == 0xffff);
}

#elif (LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_BLOCK_LANES (64)

static fast_uint_t libsais_unbwt_uniform_block(const uint8_t * RESTRICT T)
{
    uint8x16_t v0 = vld1q_u8(&T[0]), v1 = vld1q_u8(&T[16]), v2 = vld1q_u8(&T[32]), v3 = vld1q_u8(&T[48]), b = vdupq_laneq_u8(v0, 0);

    return (fast_uint_t)(vminvq_u8(vandq_u8(vandq_u8(vceqq_u8(v0, b), vceqq_u8(v1, b)), vandq_u8(vceqq_u8(v2, b), vceqq_u8(v3, b)))) == 0xff);
}

#endif

static void libsais_unbwt_compute_histogram(const uint8_t * RESTRICT T, fast_sint_t n, sa_uint_t * RESTRICT count)
{
    const fast_sint_t prefetch_distance = 256;

    const uint8_t * RESTRICT T_p = T;

    if (n >= 1024)
    {
        sa_uint_t copy[4 * (ALPHABET_SIZE + 16)];

        memset(copy, 0, (size_t)4 * (ALPHABET_SIZE + 16) * sizeof(sa_uint_t));

        sa_uint_t * RESTRICT copy0 = copy + 0 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy1 = copy + 1 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy2 = copy + 2 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy3 = copy + 3 * (ALPHABET_SIZE + 16);

        for (; T_p < (uint8_t * )((ptrdiff_t)(T + 63) & (-64)); T_p += 1) { copy0[T_p[0]]++; }

        fast_uint_t x = ((const uint32_t *)(const void *)T_p)[0], y = ((const uint32_t *)(const void *)T_p)[1];

        for (; T_p < (uint8_t * )((ptrdiff_t)(T + n - 8) & (-64)); T_p += 64)
        { 
            libsais_prefetchr(&T_p[prefetch_distance]);

#if defined(LIBSAIS_UNBWT_BLOCK_LANES)
            if (x == y && libsais_unbwt_uniform_block(T_p)) { copy0[T_p[0]] += LIBSAIS_UNBWT_BLOCK_LANES; x = ((const uint32_t *)(const void *)T_p)[16]; y = ((const uint32_t *)(const void *)T_p)[17]; continue; }
#endif

            fast_uint_t z = ((const uint32_t *)(const void *)T_p)[2], w = ((const uint32_t *)(const void *)T_p)[3];
            copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
            copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

            x = ((const uint32_t *)(const void *)T_p)[4]; y = ((const uint32_t *)(const void *)T_p)[5];
            copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
            copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

            z = ((const uint32_t *)(const void *)T_p)[6]; w = ((const uint32_t *)(const void *)T_p)[7];
            copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
            copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

            x = ((const uint32_t *)(const void *)T_p)[8]; y = ((const uint32_t *)(const void *)T_p)[9];
            copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
            copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

            z = ((const uint32_t *)(const void *)T_p)[10]; w = ((const uint32_t *)(const void *)T_p)[11];
            copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
            copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

            x = ((const uint32_t *)(const void *)T_p)[12]; y = ((const uint32_t *)(const void *)T_p)[13];
            copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
            copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

            z = ((const uint32_t *)(const void *)T_p)[14]; w = ((const uint32_t *)(const void *)T_p)[15];
            copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
            copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

            x = ((const uint32_t *)(const void *)T_p)[16]; y = ((const uint32_t *)(const void *)T_p)[17];
            copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
            copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;
        }

        copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
        copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

        T_p += 8;

        fast_uint_t i; for (i = 0; i < ALPHABET_SIZE; i++) { count[i] += copy0[i] + copy1[i] + copy2[i] + copy3[i]; }
    }

    for (; T_p < T + n; T_p += 1) { count[T_p[0]]++; }
}

static fast_sint_t libsais_unbwt_compute_histogram_checkpoints(const uint8_t * RESTRICT T, fast_sint_t n, sa_uint_t * RESTRICT count, sa_uint_t * RESTRICT checkpoints, fast_uint_t shift)
{
    const fast_sint_t prefetch_distance = 256;

    const uint8_t * RESTRICT T_p = T; fast_sint_t checkpoints_count = 0;

    if (n >= 1024)
    {
        sa_uint_t copy[4 * (ALPHABET_SIZE + 16)];

        memset(copy, 0, (size_t)4 * (ALPHABET_SIZE + 16) * sizeof(sa_uint_t));

        sa_uint_t * RESTRICT copy0 = copy + 0 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy1 = copy + 1 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy2 = copy + 2 * (ALPHABET_SIZE + 16);
        sa_uint_t * RESTRICT copy3 = copy + 3 * (ALPHABET_SIZE + 16);

        for (; T_p < (uint8_t * )((ptrdiff_t)(T + 63) & (-64)); T_p += 1) { copy0[T_p[0]]++; }

        fast_uint_t x = ((const uint32_t *)(const void *)T_p)[0], y = ((const uint32_t *)(const void *)T_p)[1];

        const uint8_t * T_blocks_end = (uint8_t * )((ptrdiff_t)(T + n - 8) & (-64));
        while (T_p < T_blocks_end)
        {
            sa_uint_t * RESTRICT checkpoint = &checkpoints[checkpoints_count++ << 8];

            fast_uint_t c; for (c = 0; c < ALPHABET_SIZE; ++c) { checkpoint[c] = copy0[c] + copy1[c] + copy2[c] + copy3[c]; }

            const uint8_t * T_chunk_end = T_blocks_end - T_p > ((fast_sint_t)1 << shift) ? T_p + ((fast_sint_t)1 << shift) : T_blocks_end;
            for (; T_p < T_chunk_end; T_p += 64)
            {
                libsais_prefetchr(&T_p[prefetch_distance]);

#if defined(LIBSAIS_UNBWT_BLOCK_LANES)
                if (x == y && libsais_unbwt_uniform_block(T_p)) { copy0[T_p[0]] += LIBSAIS_UNBWT_BLOCK_LANES; x = ((const uint32_t *)(const void *)T_p)[16]; y = ((const uint32_t *)(const void *)T_p)[17]; continue; }
#endif

                fast_uint_t z = ((const uint32_t *)(const void *)T_p)[2], w = ((const uint32_t *)(const void *)T_p)[3];
                copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
                copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

                x = ((const uint32_t *)(const void *)T_p)[4]; y = ((const uint32_t *)(const void *)T_p)[5];
                copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
                copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

                z = ((const uint32_t *)(const void *)T_p)[6]; w = ((const uint32_t *)(const void *)T_p)[7];
                copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
                copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

                x = ((const uint32_t *)(const void *)T_p)[8]; y = ((const uint32_t *)(const void *)T_p)[9];
                copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
                copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

                z = ((const uint32_t *)(const void *)T_p)[10]; w = ((const uint32_t *)(const void *)T_p)[11];
                copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
                copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

                x = ((const uint32_t *)(const void *)T_p)[12]; y = ((const uint32_t *)(const void *)T_p)[13];
                copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
                copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;

                z = ((const uint32_t *)(const void *)T_p)[14]; w = ((const uint32_t *)(const void *)T_p)[15];
                copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
                copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

                x = ((const uint32_t *)(const void *)T_p)[16]; y = ((const uint32_t *)(const void *)T_p)[17];
                copy0[(uint8_t)z]++; z >>= 8; copy1[(uint8_t)z]++; z >>= 8; copy2[(uint8_t)z]++; z >>= 8; copy3[z]++;
                copy0[(uint8_t)w]++; w >>= 8; copy1[(uint8_t)w]++; w >>= 8; copy2[(uint8_t)w]++; w >>= 8; copy3[w]++;
            }
        }

        copy0[(uint8_t)x]++; x >>= 8; copy1[(uint8_t)x]++; x >>= 8; copy2[(uint8_t)x]++; x >>= 8; copy3[x]++;
        copy0[(uint8_t)y]++; y >>= 8; copy1[(uint8_t)y]++; y >>= 8; copy2[(uint8_t)y]++; y >>= 8; copy3[y]++;

        T_p += 8;

        fast_uint_t i; for (i = 0; i < ALPHABET_SIZE; i++) { count[i] += copy0[i] + copy1[i] + copy2[i] + copy3[i]; }
    }

    for (; T_p < T + n; T_p += 1) { count[T_p[0]]++; }

    return checkpoints_count;
}

static void libsais_unbwt_transpose_bucket2(sa_uint_t * RESTRICT bucket2, fast_uint_t symbol_start, fast_uint_t symbol_end)
{
    fast_uint_t x, y;
    for (x = symbol_start; x != symbol_end; x += 4)
    {
        {
            sa_uint_t * RESTRICT bucket2_xx = &bucket2[(x << 8) + x];

            sa_uint_t a01 = bucket2_xx[0 * 256 + 1], a02 = bucket2_xx[0 * 256 + 2], a03 = bucket2_xx[0 * 256 + 3], a12 = bucket2_xx[1 * 256 + 2], a13 = bucket2_xx[1 * 256 + 3], a23 = bucket2_xx[2 * 256 + 3];
            sa_uint_t a10 = bucket2_xx[1 * 256 + 0], a20 = bucket2_xx[2 * 256 + 0], a30 = bucket2_xx[3 * 256 + 0], a21 = bucket2_xx[2 * 256 + 1], a31 = bucket2_xx[3 * 256 + 1], a32 = bucket2_xx[3 * 256 + 2];

            bucket2_xx[0 * 256 + 1] = a10; bucket2_xx[0 * 256 + 2] = a20; bucket2_xx[0 * 256 + 3] = a30; bucket2_xx[1 * 256 + 2] = a21; bucket2_xx[1 * 256 + 3] = a31; bucket2_xx[2 * 256 + 3] = a32;
            bucket2_xx[1 * 256 + 0] = a01; bucket2_xx[2 * 256 + 0] = a02; bucket2_xx[3 * 256 + 0] = a03; bucket2_xx[2 * 256 + 1] = a12; bucket2_xx[3 * 256 + 1] = a13; bucket2_xx[3 * 256 + 2] = a23;
        }

        for (y = x + 4; y != symbol_end; y += 4)
        {
            sa_uint_t * RESTRICT bucket2_xy = &bucket2[(x << 8) + y];
            sa_uint_t * RESTRICT bucket2_yx = &bucket2[(y << 8) + x];

            sa_uint_t a00 = bucket2_xy[0 * 256 + 0], a01 = bucket2_xy[0 * 256 + 1], a02 = bucket2_xy[0 * 256 + 2], a03 = bucket2_xy[0 * 256 + 3];
            sa_uint_t a10 = bucket2_xy[1 * 256 + 0], a11 = bucket2_xy[1 * 256 + 1], a12 = bucket2_xy[1 * 256 + 2], a13 = bucket2_xy[1 * 256 + 3];
            sa_uint_t a20 = bucket2_xy[2 * 256 + 0], a21 = bucket2_xy[2 * 256 + 1], a22 = bucket2_xy[2 * 256 + 2], a23 = bucket2_xy[2 * 256 + 3];
            sa_uint_t a30 = bucket2_xy[3 * 256 + 0], a31 = bucket2_xy[3 * 256 + 1], a32 = bucket2_xy[3 * 256 + 2], a33 = bucket2_xy[3 * 256 + 3];

            sa_uint_t b00 = bucket2_yx[0 * 256 + 0], b01 = bucket2_yx[0 * 256 + 1], b02 = bucket2_yx[0 * 256 + 2], b03 = bucket2_yx[0 * 256 + 3];
            sa_uint_t b10 = bucket2_yx[1 * 256 + 0], b11 = bucket2_yx[1 * 256 + 1], b12 = bucket2_yx[1 * 256 + 2], b13 = bucket2_yx[1 * 256 + 3];
            sa_uint_t b20 = bucket2_yx[2 * 256 + 0], b21 = bucket2_yx[2 * 256 + 1], b22 = bucket2_yx[2 * 256 + 2], b23 = bucket2_yx[2 * 256 + 3];
            sa_uint_t b30 = bucket2_yx[3 * 256 + 0], b31 = bucket2_yx[3 * 256 + 1], b32 = bucket2_yx[3 * 256 + 2], b33 = bucket2_yx[3 * 256 + 3];

            bucket2_xy[0 * 256 + 0] = b00; bucket2_xy[0 * 256 + 1] = b10; bucket2_xy[0 * 256 + 2] = b20; bucket2_xy[0 * 256 + 3] = b30;
            bucket2_xy[1 * 256 + 0] = b01; bucket2_xy[1 * 256 + 1] = b11; bucket2_xy[1 * 256 + 2] = b21; bucket2_xy[1 * 256 + 3] = b31;
            bucket2_xy[2 * 256 + 0] = b02; bucket2_xy[2 * 256 + 1] = b12; bucket2_xy[2 * 256 + 2] = b22; bucket2_xy[2 * 256 + 3] = b32;
            bucket2_xy[3 * 256 + 0] = b03; bucket2_xy[3 * 256 + 1] = b13; bucket2_xy[3 * 256 + 2] = b23; bucket2_xy[3 * 256 + 3] = b33;

            bucket2_yx[0 * 256 + 0] = a00; bucket2_yx[0 * 256 + 1] = a10; bucket2_yx[0 * 256 + 2] = a20; bucket2_yx[0 * 256 + 3] = a30;
            bucket2_yx[1 * 256 + 0] = a01; bucket2_yx[1 * 256 + 1] = a11; bucket2_yx[1 * 256 + 2] = a21; bucket2_yx[1 * 256 + 3] = a31;
            bucket2_yx[2 * 256 + 0] = a02; bucket2_yx[2 * 256 + 1] = a12; bucket2_yx[2 * 256 + 2] = a22; bucket2_yx[2 * 256 + 3] = a32;
            bucket2_yx[3 * 256 + 0] = a03; bucket2_yx[3 * 256 + 1] = a13; bucket2_yx[3 * 256 + 2] = a23; bucket2_yx[3 * 256 + 3] = a33;
        }
    }
}

static void libsais_unbwt_compute_bigram_histogram(const uint8_t * RESTRICT T, const sa_uint_t * RESTRICT bucket1, const sa_uint_t * RESTRICT count1, sa_uint_t * RESTRICT bucket2, fast_uint_t index, fast_uint_t symbol_start, fast_uint_t symbol_end, const sa_uint_t * RESTRICT checkpoints, fast_sint_t checkpoints_count, fast_uint_t checkpoints_shift)
{
    fast_sint_t head = (fast_sint_t)((uint8_t * )((ptrdiff_t)(T + 63) & (-64)) - T);

    fast_uint_t c;
    for (c = 0; c < ALPHABET_SIZE; ++c)
    {
        fast_uint_t prev = bucket1[c], sum = prev + count1[c];
        if (prev != sum)
        {
            sa_uint_t * RESTRICT bucket2_p = &bucket2[c << 8];

            fast_sint_t start = (fast_sint_t)(prev - (prev > index)), end = (fast_sint_t)(sum - (sum > index));
            fast_sint_t first = (start - head + ((fast_sint_t)1 << checkpoints_shift) - 1) >> checkpoints_shift, last = (end - head) >> checkpoints_shift; if (last > checkpoints_count - 1) { last = checkpoints_count - 1; }

            if (checkpoints_count > 0 && last > first)
            {
                const sa_uint_t * RESTRICT first_p = &checkpoints[first << 8];
                const sa_uint_t * RESTRICT last_p  = &checkpoints[last << 8];

                fast_uint_t d; for (d = 0; d < ALPHABET_SIZE; ++d) { bucket2_p[d] = last_p[d] - first_p[d]; }

                libsais_unbwt_compute_histogram(&T[start], head + (first << checkpoints_shift) - start, bucket2_p);
                libsais_unbwt_compute_histogram(&T[head + (last << checkpoints_shift)], end - (head + (last << checkpoints_shift)), bucket2_p);
            }
            else
            {
                libsais_unbwt_compute_histogram(&T[start], end - start, bucket2_p);
            }
        }
    }

    libsais_unbwt_transpose_bucket2(bucket2, symbol_start, symbol_end);
}

#if (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_BIGRAM_LANES (16)

static fast_uint_t libsais_unbwt_prefix_bigrams(sa_uint_t * RESTRICT bucket2, fast_uint_t * RESTRICT sum, sa_uint_t * RESTRICT ends)
{
    __m128i c0 = _mm_load_si128((const __m128i *)(const void *)&bucket2[0]), c1 = _mm_load_si128((const __m128i *)(const void *)&bucket2[4]);
    __m128i c2 = _mm_load_si128((const __m128i *)(const void *)&bucket2[8]), c3 = _mm_load_si128((const __m128i *)(const void *)&bucket2[12]);
    __m128i s  = _mm_set1_epi32((int)*sum), z = _mm_setzero_si128();

    if (_mm_movemask_epi8(_mm_cmpeq_epi32(_mm_or_si128(_mm_or_si128(c0, c1), _mm_or_si128(c2, c3)), z)) == 0xffff)
    {
        _mm_store_si128((__m128i *)(void *)&bucket2[0], s); _mm_store_si128((__m128i *)(void *)&bucket2[4], s);
        _mm_store_si128((__m128i *)(void *)&bucket2[8], s); _mm_store_si128((__m128i *)(void *)&bucket2[12], s);

        return 0;
    }

    __m128i e0 = _mm_add_epi32(c0, _mm_slli_si128(c0, 4)); e0 = _mm_add_epi32(_mm_add_epi32(e0, _mm_slli_si128(e0, 8)), s);
    __m128i e1 = _mm_add_epi32(c1, _mm_slli_si128(c1, 4)); e1 = _mm_add_epi32(_mm_add_epi32(e1, _mm_slli_si128(e1, 8)), _mm_shuffle_epi32(e0, _MM_SHUFFLE(3, 3, 3, 3)));
    __m128i e2 = _mm_add_epi32(c2, _mm_slli_si128(c2, 4)); e2 = _mm_add_epi32(_mm_add_epi32(e2, _mm_slli_si128(e2, 8)), _mm_shuffle_epi32(e1, _MM_SHUFFLE(3, 3, 3, 3)));
    __m128i e3 = _mm_add_epi32(c3, _mm_slli_si128(c3, 4)); e3 = _mm_add_epi32(_mm_add_epi32(e3, _mm_slli_si128(e3, 8)), _mm_shuffle_epi32(e2, _MM_SHUFFLE(3, 3, 3, 3)));

    _mm_store_si128((__m128i *)(void *)&bucket2[0], _mm_sub_epi32(e0, c0)); _mm_store_si128((__m128i *)(void *)&bucket2[4], _mm_sub_epi32(e1, c1));
    _mm_store_si128((__m128i *)(void *)&bucket2[8], _mm_sub_epi32(e2, c2)); _mm_store_si128((__m128i *)(void *)&bucket2[12], _mm_sub_epi32(e3, c3));
    _mm_storeu_si128((__m128i *)(void *)&ends[0], e0); _mm_storeu_si128((__m128i *)(void *)&ends[4], e1);
    _mm_storeu_si128((__m128i *)(void *)&ends[8], e2); _mm_storeu_si128((__m128i *)(void *)&ends[12], e3);

    *sum = (fast_uint_t)(uint32_t)_mm_cvtsi128_si32(_mm_shuffle_epi32(e3, _MM_SHUFFLE(3, 3, 3, 3)));

    return 0xffff & ~((fast_uint_t)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(c0, z))) | ((fast_uint_t)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(c1, z))) << 4)
        | ((fast_uint_t)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(c2, z))) << 8) | ((fast_uint_t)_mm_movemask_ps(_mm_castsi128_ps(_mm_cmpeq_epi32(c3, z))) << 12));
}

#elif (LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_BIGRAM_LANES (16)

static fast_uint_t libsais_unbwt_prefix_bigrams(sa_uint_t * RESTRICT bucket2, fast_uint_t * RESTRICT sum, sa_uint_t * RESTRICT ends)
{
    uint32x4_t c0 = vld1q_u32(&bucket2[0]), c1 = vld1q_u32(&bucket2[4]), c2 = vld1q_u32(&bucket2[8]), c3 = vld1q_u32(&bucket2[12]);
    uint32x4_t s  = vdupq_n_u32((uint32_t)*sum), z = vdupq_n_u32(0);

    if (vmaxvq_u32(vorrq_u32(vorrq_u32(c0, c1), vorrq_u32(c2, c3))) == 0)
    {
        vst1q_u32(&bucket2[0], s); vst1q_u32(&bucket2[4], s); vst1q_u32(&bucket2[8], s); vst1q_u32(&bucket2[12], s);

        return 0;
    }

    uint32x4_t e0 = vaddq_u32(c0, vextq_u32(z, c0, 3)); e0 = vaddq_u32(vaddq_u32(e0, vextq_u32(z, e0, 2)), s);
    uint32x4_t e1 = vaddq_u32(c1, vextq_u32(z, c1, 3)); e1 = vaddq_u32(vaddq_u32(e1, vextq_u32(z, e1, 2)), vdupq_laneq_u32(e0, 3));
    uint32x4_t e2 = vaddq_u32(c2, vextq_u32(z, c2, 3)); e2 = vaddq_u32(vaddq_u32(e2, vextq_u32(z, e2, 2)), vdupq_laneq_u32(e1, 3));
    uint32x4_t e3 = vaddq_u32(c3, vextq_u32(z, c3, 3)); e3 = vaddq_u32(vaddq_u32(e3, vextq_u32(z, e3, 2)), vdupq_laneq_u32(e2, 3));

    vst1q_u32(&bucket2[0], vsubq_u32(e0, c0)); vst1q_u32(&bucket2[4], vsubq_u32(e1, c1)); vst1q_u32(&bucket2[8], vsubq_u32(e2, c2)); vst1q_u32(&bucket2[12], vsubq_u32(e3, c3));
    vst1q_u32(&ends[0], e0); vst1q_u32(&ends[4], e1); vst1q_u32(&ends[8], e2); vst1q_u32(&ends[12], e3);

    *sum = (fast_uint_t)vgetq_lane_u32(e3, 3);

    return (fast_uint_t)libsais_mask(vcombine_u8(vmovn_u16(vcombine_u16(vmovn_u32(vtstq_u32(c0, c0)), vmovn_u32(vtstq_u32(c1, c1)))), vmovn_u16(vcombine_u16(vmovn_u32(vtstq_u32(c2, c2)), vmovn_u32(vtstq_u32(c3, c3))))));
}

#endif

#if defined(LIBSAIS_UNBWT_BIGRAM_LANES)
    #if defined(__GNUC__) || defined(__clang__)
        #define libsais_unbwt_ctz(x) ((fast_uint_t)__builtin_ctz((unsigned int)(x)))
    #else
        static fast_uint_t libsais_unbwt_ctz(fast_uint_t x) { unsigned long index; _BitScanForward(&index, (unsigned long)x); return (fast_uint_t)index; }
    #endif
#endif

static void libsais_unbwt_calculate_fastbits(sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t lastc, fast_uint_t shift, fast_uint_t symbol_start, fast_uint_t symbol_end, fast_uint_t v, fast_uint_t sum)
{
    fast_uint_t w, c, d;
    for (w = symbol_start << 8, c = symbol_start; c < symbol_end; ++c)
    {
        if (c == lastc) { sum += 1; }

#if defined(LIBSAIS_UNBWT_BIGRAM_LANES)
        for (d = 0; d < ALPHABET_SIZE; d += LIBSAIS_UNBWT_BIGRAM_LANES, w += LIBSAIS_UNBWT_BIGRAM_LANES)
        {
            sa_uint_t ends[LIBSAIS_UNBWT_BIGRAM_LANES]; fast_uint_t nonempty;

            for (nonempty = libsais_unbwt_prefix_bigrams(&bucket2[w], &sum, ends); nonempty != 0; nonempty &= nonempty - 1)
            {
                fast_uint_t k = libsais_unbwt_ctz(nonempty);
                for (; v <= (((fast_uint_t)ends[k] - 1) >> shift); ++v) { fastbits[v] = (uint16_t)(w + k); }
            }
        }
#else
        for (d = 0; d < ALPHABET_SIZE; ++d, ++w)
        {
            fast_uint_t prev = sum; sum += bucket2[w]; bucket2[w] = (sa_uint_t)prev;
            if (prev != sum)
            {
                for (; v <= ((sum - 1) >> shift); ++v) { fastbits[v] = (uint16_t)w; }
            }
        }
#endif
    }
}

static void libsais_unbwt_calculate_bigram_starts(sa_uint_t * RESTRICT bucket2, fast_uint_t lastc, fast_uint_t symbol_start, fast_uint_t symbol_end, fast_uint_t sum)
{
    fast_uint_t w, c, d;
    for (w = symbol_start << 8, c = symbol_start; c < symbol_end; ++c)
    {
        if (c == lastc) { sum += 1; }

#if defined(LIBSAIS_UNBWT_BIGRAM_LANES)
        for (d = 0; d < ALPHABET_SIZE; d += LIBSAIS_UNBWT_BIGRAM_LANES, w += LIBSAIS_UNBWT_BIGRAM_LANES)
        {
            sa_uint_t ends[LIBSAIS_UNBWT_BIGRAM_LANES]; libsais_unbwt_prefix_bigrams(&bucket2[w], &sum, ends);
        }
#else
        for (d = 0; d < ALPHABET_SIZE; ++d, ++w)
        {
            fast_uint_t prev = sum; sum += bucket2[w]; bucket2[w] = (sa_uint_t)prev;
        }
#endif
    }
}

static void libsais_unbwt_calculate_first_symbols(const sa_uint_t * RESTRICT bucket1, sa_sint_t n, uint8_t * RESTRICT first_symbols)
{
    sa_uint_t * RESTRICT symbol_ends = (sa_uint_t *)(void *)(first_symbols + ((fast_uint_t)1 << UNBWT_FIRST_SYMBOLS));

    fast_uint_t first_shift = 0; while ((n >> first_shift) >= ((sa_sint_t)1 << UNBWT_FIRST_SYMBOLS)) { first_shift++; }
    fast_uint_t c, j, cells = ((fast_uint_t)n >> first_shift) + 1;

    for (c = 0; c < ALPHABET_SIZE - 1; ++c) { symbol_ends[c] = bucket1[c + 1]; }
    symbol_ends[ALPHABET_SIZE - 1] = (sa_uint_t)n + 1;

    for (c = 0; symbol_ends[c] <= 1; ++c) { }
    for (j = 0; j < cells; ++c)
    {
        fast_uint_t e = (((fast_uint_t)symbol_ends[c] - 1) >> first_shift) + 1; if (e > cells) { e = cells; }
        if (e > j) { memset(&first_symbols[j], (int)c, (size_t)(e - j)); j = e; }
    }
}

#if (LIBSAIS_CPU_FEATURE >= LIBSAIS_CPU_FEATURE_SSE2) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_RUN_LANES (16)

static fast_uint_t libsais_unbwt_uniform_symbols(const uint8_t * RESTRICT T)
{
    __m128i v = _mm_loadu_si128((const __m128i *)(const void *)T);

    return (fast_uint_t)((_mm_movemask_epi8(_mm_cmpeq_epi8(v, _mm_srli_si128(v, 1))) | 0x8000) == 0xffff);
}

static void libsais_unbwt_store_rows(sa_uint_t * RESTRICT P, fast_uint_t i)
{
    __m128i r = _mm_add_epi32(_mm_set1_epi32((int)i), _mm_setr_epi32(0, 1, 2, 3)), four = _mm_set1_epi32(4);

    _mm_storeu_si128((__m128i *)(void *)&P[0], r); r = _mm_add_epi32(r, four);
    _mm_storeu_si128((__m128i *)(void *)&P[4], r); r = _mm_add_epi32(r, four);
    _mm_storeu_si128((__m128i *)(void *)&P[8], r); r = _mm_add_epi32(r, four);
    _mm_storeu_si128((__m128i *)(void *)&P[12], r);
}

static void libsais_unbwt_store_symbol_rows(sa_uint_t * RESTRICT P, fast_uint_t i, fast_uint_t c)
{
    __m128i r = _mm_add_epi32(_mm_set1_epi32((int)((i << 8) + c)), _mm_setr_epi32(0, 256, 512, 768)), step = _mm_set1_epi32(1024);

    _mm_storeu_si128((__m128i *)(void *)&P[0], r); r = _mm_add_epi32(r, step);
    _mm_storeu_si128((__m128i *)(void *)&P[4], r); r = _mm_add_epi32(r, step);
    _mm_storeu_si128((__m128i *)(void *)&P[8], r); r = _mm_add_epi32(r, step);
    _mm_storeu_si128((__m128i *)(void *)&P[12], r);
}

static void libsais_unbwt_add_rows(sa_uint_t * P_target, const sa_uint_t * P_source)
{
    const __m128i low = _mm_set1_epi32(0xffff);

    __m128i s0 = _mm_slli_epi32(_mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&P_source[0]), low), 16);
    __m128i s1 = _mm_slli_epi32(_mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&P_source[4]), low), 16);
    __m128i s2 = _mm_slli_epi32(_mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&P_source[8]), low), 16);
    __m128i s3 = _mm_slli_epi32(_mm_and_si128(_mm_loadu_si128((const __m128i *)(const void *)&P_source[12]), low), 16);

    _mm_storeu_si128((__m128i *)(void *)&P_target[0], _mm_add_epi32(_mm_loadu_si128((const __m128i *)(const void *)&P_target[0]), s0));
    _mm_storeu_si128((__m128i *)(void *)&P_target[4], _mm_add_epi32(_mm_loadu_si128((const __m128i *)(const void *)&P_target[4]), s1));
    _mm_storeu_si128((__m128i *)(void *)&P_target[8], _mm_add_epi32(_mm_loadu_si128((const __m128i *)(const void *)&P_target[8]), s2));
    _mm_storeu_si128((__m128i *)(void *)&P_target[12], _mm_add_epi32(_mm_loadu_si128((const __m128i *)(const void *)&P_target[12]), s3));
}

#elif (LIBSAIS_CPU_FEATURE == LIBSAIS_CPU_FEATURE_A64) && (SAINT_BIT == 32)

#define LIBSAIS_UNBWT_RUN_LANES (16)

static fast_uint_t libsais_unbwt_uniform_symbols(const uint8_t * RESTRICT T)
{
    uint8x16_t v = vld1q_u8(T);

    return (fast_uint_t)(vminvq_u8(vceqq_u8(v, vextq_u8(v, v, 1))) == 0xff);
}

static void libsais_unbwt_store_rows(sa_uint_t * RESTRICT P, fast_uint_t i)
{
    const uint32_t lanes[4] = { 0, 1, 2, 3 };

    uint32x4_t r = vaddq_u32(vdupq_n_u32((uint32_t)i), vld1q_u32(lanes)), four = vdupq_n_u32(4);

    vst1q_u32(&P[0], r); r = vaddq_u32(r, four);
    vst1q_u32(&P[4], r); r = vaddq_u32(r, four);
    vst1q_u32(&P[8], r); r = vaddq_u32(r, four);
    vst1q_u32(&P[12], r);
}

static void libsais_unbwt_store_symbol_rows(sa_uint_t * RESTRICT P, fast_uint_t i, fast_uint_t c)
{
    const uint32_t lanes[4] = { 0, 256, 512, 768 };

    uint32x4_t r = vaddq_u32(vdupq_n_u32((uint32_t)((i << 8) + c)), vld1q_u32(lanes)), step = vdupq_n_u32(1024);

    vst1q_u32(&P[0], r); r = vaddq_u32(r, step);
    vst1q_u32(&P[4], r); r = vaddq_u32(r, step);
    vst1q_u32(&P[8], r); r = vaddq_u32(r, step);
    vst1q_u32(&P[12], r);
}

static void libsais_unbwt_add_rows(sa_uint_t * P_target, const sa_uint_t * P_source)
{
    const uint32x4_t low = vdupq_n_u32(0xffff);

    uint32x4_t s0 = vshlq_n_u32(vandq_u32(vld1q_u32(&P_source[0]), low), 16), s1 = vshlq_n_u32(vandq_u32(vld1q_u32(&P_source[4]), low), 16);
    uint32x4_t s2 = vshlq_n_u32(vandq_u32(vld1q_u32(&P_source[8]), low), 16), s3 = vshlq_n_u32(vandq_u32(vld1q_u32(&P_source[12]), low), 16);

    vst1q_u32(&P_target[0], vaddq_u32(vld1q_u32(&P_target[0]), s0)); vst1q_u32(&P_target[4], vaddq_u32(vld1q_u32(&P_target[4]), s1));
    vst1q_u32(&P_target[8], vaddq_u32(vld1q_u32(&P_target[8]), s2)); vst1q_u32(&P_target[12], vaddq_u32(vld1q_u32(&P_target[12]), s3));
}

#endif

static void libsais_unbwt_calculate_biPSI(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, const sa_uint_t * RESTRICT bucket1_start, sa_uint_t * RESTRICT bucket2, fast_uint_t index, fast_sint_t omp_block_start, fast_sint_t omp_block_end)
{
    sa_uint_t bucket1[ALPHABET_SIZE]; memcpy(bucket1, bucket1_start, ALPHABET_SIZE * sizeof(sa_uint_t));

    {
        fast_sint_t i = omp_block_start, j = (fast_sint_t)index; if (omp_block_end < j) { j = omp_block_end; }

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i < j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i]))
            {
                fast_uint_t c = T[i], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_rows(&P[d], (fast_uint_t)i); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e;
            for (e = i + LIBSAIS_UNBWT_RUN_LANES; i < e; i += 4)
            {
                {
                    fast_uint_t c = T[i + 0];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                }
                {
                    fast_uint_t c = T[i + 1];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                }
                {
                    fast_uint_t c = T[i + 2];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                }
                {
                    fast_uint_t c = T[i + 3];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                }
            }
        }
#endif

        for (; i < j; ++i)
        {
            fast_uint_t c = T[i];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)i;
            }
        }
    }

    {
        fast_sint_t i = (fast_sint_t)index, j = omp_block_end; if (omp_block_start > i) { i = omp_block_start; }

        i += 1;

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i <= j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i - 1]))
            {
                fast_uint_t c = T[i - 1], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_rows(&P[d], (fast_uint_t)i); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e;
            for (e = i + LIBSAIS_UNBWT_RUN_LANES; i < e; i += 4)
            {
                {
                    fast_uint_t c = T[i + 0 - 1];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                }
                {
                    fast_uint_t c = T[i + 1 - 1];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                }
                {
                    fast_uint_t c = T[i + 2 - 1];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                }
                {
                    fast_uint_t c = T[i + 3 - 1];
                    fast_uint_t p = bucket1[c]++;
                    fast_sint_t t = (fast_sint_t)(index - p);
                    if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                }
            }
        }
#endif

        for (; i <= j; ++i)
        {
            fast_uint_t c = T[i - 1];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)i;
            }
        }
    }
}

static fast_uint_t libsais_unbwt_mark_row(fast_uint_t mark, fast_uint_t marks_shift, fast_uint_t end)
{
    uint32_t hash = (uint32_t)(mark - (end >> marks_shift)) * (uint32_t)0x9E3779B9;

    return (mark << marks_shift) + ((end + (fast_uint_t)(((uint64_t)hash << marks_shift) >> 32)) & (((fast_uint_t)1 << marks_shift) - 1));
}

static void libsais_unbwt_calculate_biPSI_marks(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, const sa_uint_t * RESTRICT bucket1_start, sa_uint_t * RESTRICT bucket2, fast_uint_t index, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, fast_uint_t end, fast_sint_t omp_block_start, fast_sint_t omp_block_end)
{
    sa_uint_t bucket1[ALPHABET_SIZE]; memcpy(bucket1, bucket1_start, ALPHABET_SIZE * sizeof(sa_uint_t));

    fast_uint_t row  = (fast_uint_t)omp_block_start + ((fast_uint_t)omp_block_start >= index);
    fast_uint_t mark = row >> marks_shift, mark_row = libsais_unbwt_mark_row(mark, marks_shift, end); if (mark_row < row) { mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

    {
        fast_sint_t i = omp_block_start, j = (fast_sint_t)index; if (omp_block_end < j) { j = omp_block_end; }

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i < j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i]))
            {
                fast_uint_t c = T[i], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    while (mark_row < (fast_uint_t)i + LIBSAIS_UNBWT_RUN_LANES) { marks[mark] = (sa_uint_t)(p + (mark_row - (fast_uint_t)i)); mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_rows(&P[d], (fast_uint_t)i); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e = i + LIBSAIS_UNBWT_RUN_LANES;
            if ((fast_uint_t)e <= mark_row)
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                    }
                    {
                        fast_uint_t c = T[i + 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                    }
                    {
                        fast_uint_t c = T[i + 2];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                    }
                    {
                        fast_uint_t c = T[i + 3];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                    }
                }
            }
            else
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 0) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                    }
                    {
                        fast_uint_t c = T[i + 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 1) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                    }
                    {
                        fast_uint_t c = T[i + 2];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 2) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                    }
                    {
                        fast_uint_t c = T[i + 3];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 3) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                    }
                }
            }
        }
#endif

        for (; i < j; ++i)
        {
            fast_uint_t c = T[i];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if ((fast_uint_t)i == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)i;
            }
        }
    }

    if (mark_row == index) { mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

    {
        fast_sint_t i = (fast_sint_t)index, j = omp_block_end; if (omp_block_start > i) { i = omp_block_start; }

        i += 1;

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i <= j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i - 1]))
            {
                fast_uint_t c = T[i - 1], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    while (mark_row < (fast_uint_t)i + LIBSAIS_UNBWT_RUN_LANES) { marks[mark] = (sa_uint_t)(p + (mark_row - (fast_uint_t)i)); mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_rows(&P[d], (fast_uint_t)i); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e = i + LIBSAIS_UNBWT_RUN_LANES;
            if ((fast_uint_t)e <= mark_row)
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                    }
                    {
                        fast_uint_t c = T[i + 1 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                    }
                    {
                        fast_uint_t c = T[i + 2 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                    }
                    {
                        fast_uint_t c = T[i + 3 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                    }
                }
            }
            else
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 0) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 0); }
                    }
                    {
                        fast_uint_t c = T[i + 1 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 1) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 1); }
                    }
                    {
                        fast_uint_t c = T[i + 2 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 2) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 2); }
                    }
                    {
                        fast_uint_t c = T[i + 3 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 3) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(i + 3); }
                    }
                }
            }
        }
#endif

        for (; i <= j; ++i)
        {
            fast_uint_t c = T[i - 1];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if ((fast_uint_t)i == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)i;
            }
        }
    }
}

static void libsais_unbwt_calculate_biPSI_marks_symbol_rows(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, const sa_uint_t * RESTRICT bucket1_start, sa_uint_t * RESTRICT bucket2, fast_uint_t index, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, fast_uint_t end, fast_sint_t omp_block_start, fast_sint_t omp_block_end)
{
    sa_uint_t bucket1[ALPHABET_SIZE]; memcpy(bucket1, bucket1_start, ALPHABET_SIZE * sizeof(sa_uint_t));

    fast_uint_t row  = (fast_uint_t)omp_block_start + ((fast_uint_t)omp_block_start >= index);
    fast_uint_t mark = row >> marks_shift, mark_row = libsais_unbwt_mark_row(mark, marks_shift, end); if (mark_row < row) { mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

    {
        fast_sint_t i = omp_block_start, j = (fast_sint_t)index; if (omp_block_end < j) { j = omp_block_end; }

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i < j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i]))
            {
                fast_uint_t c = T[i], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    while (mark_row < (fast_uint_t)i + LIBSAIS_UNBWT_RUN_LANES) { marks[mark] = (sa_uint_t)(p + (mark_row - (fast_uint_t)i)); mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_symbol_rows(&P[d], (fast_uint_t)i, c); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e = i + LIBSAIS_UNBWT_RUN_LANES;
            if ((fast_uint_t)e <= mark_row)
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 0) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 1) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 2];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 2) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 3];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 3) << 8) + c); }
                    }
                }
            }
            else
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 0) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 0) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 1) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 1) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 2];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 2) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 2) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 3];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 3) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 3) << 8) + c); }
                    }
                }
            }
        }
#endif

        for (; i < j; ++i)
        {
            fast_uint_t c = T[i];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if ((fast_uint_t)i == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)i << 8) + c);
            }
        }
    }

    if (mark_row == index) { mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

    {
        fast_sint_t i = (fast_sint_t)index, j = omp_block_end; if (omp_block_start > i) { i = omp_block_start; }

        i += 1;

#if defined(LIBSAIS_UNBWT_RUN_LANES)
        for (; i <= j - (LIBSAIS_UNBWT_RUN_LANES - 1); )
        {
            if (libsais_unbwt_uniform_symbols(&T[i - 1]))
            {
                fast_uint_t c = T[i - 1], p = bucket1[c];

                if ((p + LIBSAIS_UNBWT_RUN_LANES - 1 < index || p > index) && libsais_unbwt_uniform_symbols(&T[p - (p > index)]))
                {
                    fast_uint_t w = ((fast_uint_t)T[p - (p > index)] << 8) + c, d = bucket2[w];

                    while (mark_row < (fast_uint_t)i + LIBSAIS_UNBWT_RUN_LANES) { marks[mark] = (sa_uint_t)(p + (mark_row - (fast_uint_t)i)); mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

                    bucket1[c] = (sa_uint_t)(p + LIBSAIS_UNBWT_RUN_LANES); bucket2[w] = (sa_uint_t)(d + LIBSAIS_UNBWT_RUN_LANES);
                    libsais_unbwt_store_symbol_rows(&P[d], (fast_uint_t)i, c); i += LIBSAIS_UNBWT_RUN_LANES; continue;
                }
            }

            fast_sint_t e = i + LIBSAIS_UNBWT_RUN_LANES;
            if ((fast_uint_t)e <= mark_row)
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 0) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 1 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 1) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 2 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 2) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 3 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 3) << 8) + c); }
                    }
                }
            }
            else
            {
                for (; i < e; i += 4)
                {
                    {
                        fast_uint_t c = T[i + 0 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 0) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 0) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 1 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 1) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 1) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 2 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 2) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 2) << 8) + c); }
                    }
                    {
                        fast_uint_t c = T[i + 3 - 1];
                        fast_uint_t p = bucket1[c]++;
                        fast_sint_t t = (fast_sint_t)(index - p);
                        if ((fast_uint_t)(i + 3) == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }
                        if (t != 0) { fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c; P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)(i + 3) << 8) + c); }
                    }
                }
            }
        }
#endif

        for (; i <= j; ++i)
        {
            fast_uint_t c = T[i - 1];
            fast_uint_t p = bucket1[c]++;
            fast_sint_t t = (fast_sint_t)(index - p);

            if ((fast_uint_t)i == mark_row) { marks[mark] = (sa_uint_t)p; mark_row = libsais_unbwt_mark_row(++mark, marks_shift, end); }

            if (t != 0)
            {
                fast_uint_t w = (((fast_uint_t)T[p + (fast_uint_t)(t >> ((sizeof(fast_sint_t) * 8) - 1))]) << 8) + c;
                P[bucket2[w]++] = (sa_uint_t)(((fast_uint_t)i << 8) + c);
            }
        }
    }
}

static void libsais_unbwt_init_single(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, sa_sint_t n, const sa_sint_t * freq, fast_uint_t index, fast_uint_t lastc, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, fast_uint_t end)
{
    sa_uint_t bucket1[ALPHABET_SIZE];
    sa_uint_t count1[ALPHABET_SIZE];

    fast_uint_t shift = 0; while ((n >> shift) > ((sa_sint_t)1 << UNBWT_FASTBITS)) { shift++; }
    fast_uint_t checkpoints_shift = 11; while (checkpoints_shift < 16 && ((fast_sint_t)1 << (2 * checkpoints_shift - 1)) < (fast_sint_t)n) { checkpoints_shift++; }
    fast_sint_t checkpoints_count = 0;

    if (freq != NULL)
    {
        memcpy(count1, freq, ALPHABET_SIZE * sizeof(sa_uint_t));
    }
    else
    {
        memset(count1, 0, ALPHABET_SIZE * sizeof(sa_uint_t));
        checkpoints_count = libsais_unbwt_compute_histogram_checkpoints(T, n, count1, &P[1], checkpoints_shift);
    }

    {
        fast_uint_t sum, c;
        for (sum = 1, c = 0; c < ALPHABET_SIZE; ++c) { bucket1[c] = (sa_uint_t)sum; sum += count1[c]; }
    }

    fast_uint_t symbol_start = 0;             while (symbol_start < ALPHABET_SIZE - 1 && count1[symbol_start] == 0) { symbol_start++; }
    fast_uint_t symbol_end   = ALPHABET_SIZE; while (symbol_end > symbol_start + 1 && count1[symbol_end - 1] == 0) { symbol_end--; }

    symbol_start = symbol_start & ~(fast_uint_t)15; symbol_end = (symbol_end + 15) & ~(fast_uint_t)15;

    memset(&bucket2[symbol_start << 8], 0, (size_t)((symbol_end - symbol_start) << 8) * sizeof(sa_uint_t));
    libsais_unbwt_compute_bigram_histogram(T, bucket1, count1, bucket2, index, symbol_start, symbol_end, &P[1], checkpoints_count, checkpoints_shift);

    if (lastc < ALPHABET_SIZE) { P[bucket1[lastc]] = 0; }

    if (marks != NULL && n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
    {
        libsais_unbwt_calculate_bigram_starts(bucket2, lastc, symbol_start, symbol_end, 1);
        libsais_unbwt_calculate_first_symbols(bucket1, n, (uint8_t *)(void *)fastbits);
        libsais_unbwt_calculate_biPSI_marks_symbol_rows(T, P, bucket1, bucket2, index, marks, marks_shift, end, 0, n);
    }
    else
    {
        libsais_unbwt_calculate_fastbits(bucket2, fastbits, lastc, shift, symbol_start, symbol_end, 0, 1);
        if (marks != NULL)
        {
            libsais_unbwt_calculate_biPSI_marks(T, P, bucket1, bucket2, index, marks, marks_shift, end, 0, n);
        }
        else
        {
            libsais_unbwt_calculate_biPSI(T, P, bucket1, bucket2, index, 0, n);
        }
    }
}

#if defined(LIBSAIS_OPENMP)

static void libsais_unbwt_init_parallel(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, sa_sint_t n, const sa_sint_t * freq, fast_uint_t index, fast_uint_t lastc, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, sa_uint_t * RESTRICT buckets, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, fast_uint_t end, sa_sint_t threads)
{
    sa_uint_t bucket1[ALPHABET_SIZE];

    fast_uint_t shift = 0; while ((n >> shift) > ((sa_sint_t)1 << UNBWT_FASTBITS)) { shift++; }

    fast_uint_t symbol_start = 0, symbol_end = ALPHABET_SIZE;

    memset(bucket1, 0, ALPHABET_SIZE * sizeof(sa_uint_t));

    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
    {
        fast_sint_t omp_thread_num  = omp_get_thread_num();
        fast_sint_t omp_num_threads = omp_get_num_threads();

        if (omp_num_threads == 1)
        {
            libsais_unbwt_init_single(T, P, n, freq, index, lastc, bucket2, fastbits, marks, marks_shift, end);
        }
        else
        {
            sa_uint_t * RESTRICT bucket1_local  = buckets + omp_thread_num * (2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE));
            sa_uint_t * RESTRICT count1_local   = bucket1_local + ALPHABET_SIZE;
            sa_uint_t * RESTRICT bucket2_local  = count1_local + ALPHABET_SIZE;

            fast_sint_t omp_block_stride        = (n / omp_num_threads) & (-16);
            fast_sint_t omp_block_start         = omp_thread_num * omp_block_stride;
            fast_sint_t omp_block_size          = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

            {
                memset(count1_local, 0, ALPHABET_SIZE * sizeof(sa_uint_t));
                libsais_unbwt_compute_histogram(T + omp_block_start, omp_block_size, count1_local);
            }

            #pragma omp barrier

            #pragma omp master
            {
                {
                    sa_uint_t * RESTRICT bucket1_temp = buckets;
                    sa_uint_t * RESTRICT count1_temp  = buckets + ALPHABET_SIZE;

                    fast_sint_t t;
                    for (t = 0; t < omp_num_threads; ++t, bucket1_temp += 2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE), count1_temp += 2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE))
                    {
                        fast_sint_t c; for (c = 0; c < ALPHABET_SIZE; c += 1) { sa_uint_t A = bucket1[c], B = count1_temp[c]; bucket1[c] = A + B; bucket1_temp[c] = A; }
                    }
                }

                {
                    while (symbol_start < ALPHABET_SIZE - 1 && bucket1[symbol_start] == 0) { symbol_start++; }
                    while (symbol_end > symbol_start + 1 && bucket1[symbol_end - 1] == 0) { symbol_end--; }

                    symbol_start = symbol_start & ~(fast_uint_t)15; symbol_end = (symbol_end + 15) & ~(fast_uint_t)15;
                }

                {
                    fast_uint_t sum, c;
                    for (sum = 1, c = 0; c < ALPHABET_SIZE; ++c) { fast_uint_t prev = sum; sum += bucket1[c]; bucket1[c] = (sa_uint_t)prev; }
                }

                if (lastc < ALPHABET_SIZE) { P[bucket1[lastc]] = 0; }

                if (marks != NULL && n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
                {
                    libsais_unbwt_calculate_first_symbols(bucket1, n, (uint8_t *)(void *)fastbits);
                }
            }

            #pragma omp barrier

            {
                fast_sint_t c; for (c = 0; c < ALPHABET_SIZE; c += 1) { sa_uint_t A = bucket1[c], B = bucket1_local[c]; bucket1_local[c] = A + B; }

                memset(&bucket2_local[symbol_start << 8], 0, (size_t)((symbol_end - symbol_start) << 8) * sizeof(sa_uint_t));
                libsais_unbwt_compute_bigram_histogram(T, bucket1_local, count1_local, bucket2_local, index, symbol_start, symbol_end, NULL, 0, 0);
            }

            #pragma omp barrier

            {
                fast_sint_t omp_bucket2_stride  = ((fast_sint_t)((symbol_end - symbol_start) << 8) / omp_num_threads) & (-ALPHABET_SIZE);
                fast_sint_t omp_bucket2_start   = (fast_sint_t)(symbol_start << 8) + omp_thread_num * omp_bucket2_stride;
                fast_sint_t omp_bucket2_size    = omp_thread_num < omp_num_threads - 1 ? omp_bucket2_stride : (fast_sint_t)(symbol_end << 8) - omp_bucket2_start;

                memset(&bucket2[omp_bucket2_start], 0, (size_t)omp_bucket2_size * sizeof(sa_uint_t));

                {
                    sa_uint_t * RESTRICT bucket2_temp = buckets + 2 * ALPHABET_SIZE;

                    fast_sint_t t;
                    for (t = 0; t < omp_num_threads; ++t, bucket2_temp += 2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE))
                    {
                        fast_sint_t c; for (c = omp_bucket2_start; c < omp_bucket2_start + omp_bucket2_size; c += 1) { bucket2[c] += bucket2_temp[c]; }
                    }
                }

                {
                    fast_uint_t row_start = (fast_uint_t)omp_bucket2_start >> 8, row_end = (fast_uint_t)(omp_bucket2_start + omp_bucket2_size) >> 8, v = 0, c;

                    if (marks != NULL && n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
                    {
                        libsais_unbwt_calculate_bigram_starts(bucket2, lastc, row_start, row_end, bucket1[row_start]);
                    }
                    else
                    {
                        for (c = row_start; c > 0; --c) { if ((fast_uint_t)(bucket1[c] - bucket1[c - 1]) > (fast_uint_t)(c - 1 == lastc)) { v = (((fast_uint_t)bucket1[c] - 1) >> shift) + 1; break; } }

                        libsais_unbwt_calculate_fastbits(bucket2, fastbits, lastc, shift, row_start, row_end, v, bucket1[row_start]);
                    }
                }

                {
                    sa_uint_t * RESTRICT bucket2_temp = buckets + 2 * ALPHABET_SIZE;

                    fast_sint_t t;
                    for (t = 0; t < omp_num_threads; ++t, bucket2_temp += 2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE))
                    {
                        fast_sint_t c; for (c = omp_bucket2_start; c < omp_bucket2_start + omp_bucket2_size; c += 1) { sa_uint_t A = bucket2[c], B = bucket2_temp[c]; bucket2[c] = A + B; bucket2_temp[c] = A; }
                    }
                }
            }

            #pragma omp barrier

            {
                if (marks != NULL)
                {
                    if (n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
                    {
                        libsais_unbwt_calculate_biPSI_marks_symbol_rows(T, P, bucket1_local, bucket2_local, index, marks, marks_shift, end, omp_block_start, omp_block_start + omp_block_size);
                    }
                    else
                    {
                        libsais_unbwt_calculate_biPSI_marks(T, P, bucket1_local, bucket2_local, index, marks, marks_shift, end, omp_block_start, omp_block_start + omp_block_size);
                    }
                }
                else
                {
                    libsais_unbwt_calculate_biPSI(T, P, bucket1_local, bucket2_local, index, omp_block_start, omp_block_start + omp_block_size);
                }
            }
        }
    }
}

#endif

static void libsais_unbwt_decode_1(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t * i0, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;

    fast_uint_t i, p0 = *i0;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
    }

    *i0 = p0;
}

static void libsais_unbwt_decode_2(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
    }

    *i0 = p0; *i1 = p1;
}

static void libsais_unbwt_decode_3(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
    }

    *i0 = p0; *i1 = p1; *i2 = p2;
}

static void libsais_unbwt_decode_4(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3;
}

static void libsais_unbwt_decode_5(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4;
}

static void libsais_unbwt_decode_6(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5;
}

static void libsais_unbwt_decode_7(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6;
}

static void libsais_unbwt_decode_8(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7;
}

static void libsais_unbwt_decode_9(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8;
}

static void libsais_unbwt_decode_10(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9;
}

static void libsais_unbwt_decode_11(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10;
}

static void libsais_unbwt_decode_12(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t * i11, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);
    uint16_t * RESTRICT U11 = (uint16_t *)(void *)(((uint8_t *)U10) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10, p11 = *i11;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
        uint16_t c11 = fastbits[p11 >> shift]; if (bucket2[c11] <= p11) { do { c11++; } while (bucket2[c11] <= p11); } p11 = P[p11]; U11[i] = libsais_bswap16(c11);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10; *i11 = p11;
}

static void libsais_unbwt_decode_13(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t * i11, fast_uint_t * i12, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);
    uint16_t * RESTRICT U11 = (uint16_t *)(void *)(((uint8_t *)U10) + r);
    uint16_t * RESTRICT U12 = (uint16_t *)(void *)(((uint8_t *)U11) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10, p11 = *i11, p12 = *i12;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
        uint16_t c11 = fastbits[p11 >> shift]; if (bucket2[c11] <= p11) { do { c11++; } while (bucket2[c11] <= p11); } p11 = P[p11]; U11[i] = libsais_bswap16(c11);
        uint16_t c12 = fastbits[p12 >> shift]; if (bucket2[c12] <= p12) { do { c12++; } while (bucket2[c12] <= p12); } p12 = P[p12]; U12[i] = libsais_bswap16(c12);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10; *i11 = p11; *i12 = p12;
}

static void libsais_unbwt_decode_14(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t * i11, fast_uint_t * i12, fast_uint_t * i13, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);
    uint16_t * RESTRICT U11 = (uint16_t *)(void *)(((uint8_t *)U10) + r);
    uint16_t * RESTRICT U12 = (uint16_t *)(void *)(((uint8_t *)U11) + r);
    uint16_t * RESTRICT U13 = (uint16_t *)(void *)(((uint8_t *)U12) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10, p11 = *i11, p12 = *i12, p13 = *i13;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
        uint16_t c11 = fastbits[p11 >> shift]; if (bucket2[c11] <= p11) { do { c11++; } while (bucket2[c11] <= p11); } p11 = P[p11]; U11[i] = libsais_bswap16(c11);
        uint16_t c12 = fastbits[p12 >> shift]; if (bucket2[c12] <= p12) { do { c12++; } while (bucket2[c12] <= p12); } p12 = P[p12]; U12[i] = libsais_bswap16(c12);
        uint16_t c13 = fastbits[p13 >> shift]; if (bucket2[c13] <= p13) { do { c13++; } while (bucket2[c13] <= p13); } p13 = P[p13]; U13[i] = libsais_bswap16(c13);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10; *i11 = p11; *i12 = p12; *i13 = p13;
}

static void libsais_unbwt_decode_15(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t * i11, fast_uint_t * i12, fast_uint_t * i13, fast_uint_t * i14, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);
    uint16_t * RESTRICT U11 = (uint16_t *)(void *)(((uint8_t *)U10) + r);
    uint16_t * RESTRICT U12 = (uint16_t *)(void *)(((uint8_t *)U11) + r);
    uint16_t * RESTRICT U13 = (uint16_t *)(void *)(((uint8_t *)U12) + r);
    uint16_t * RESTRICT U14 = (uint16_t *)(void *)(((uint8_t *)U13) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10, p11 = *i11, p12 = *i12, p13 = *i13, p14 = *i14;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
        uint16_t c11 = fastbits[p11 >> shift]; if (bucket2[c11] <= p11) { do { c11++; } while (bucket2[c11] <= p11); } p11 = P[p11]; U11[i] = libsais_bswap16(c11);
        uint16_t c12 = fastbits[p12 >> shift]; if (bucket2[c12] <= p12) { do { c12++; } while (bucket2[c12] <= p12); } p12 = P[p12]; U12[i] = libsais_bswap16(c12);
        uint16_t c13 = fastbits[p13 >> shift]; if (bucket2[c13] <= p13) { do { c13++; } while (bucket2[c13] <= p13); } p13 = P[p13]; U13[i] = libsais_bswap16(c13);
        uint16_t c14 = fastbits[p14 >> shift]; if (bucket2[c14] <= p14) { do { c14++; } while (bucket2[c14] <= p14); } p14 = P[p14]; U14[i] = libsais_bswap16(c14);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10; *i11 = p11; *i12 = p12; *i13 = p13; *i14 = p14;
}

static void libsais_unbwt_decode_16(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, fast_uint_t r, fast_uint_t * i0, fast_uint_t * i1, fast_uint_t * i2, fast_uint_t * i3, fast_uint_t * i4, fast_uint_t * i5, fast_uint_t * i6, fast_uint_t * i7, fast_uint_t * i8, fast_uint_t * i9, fast_uint_t * i10, fast_uint_t * i11, fast_uint_t * i12, fast_uint_t * i13, fast_uint_t * i14, fast_uint_t * i15, fast_uint_t k)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)U;
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(((uint8_t *)U0) + r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(((uint8_t *)U1) + r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(((uint8_t *)U2) + r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(((uint8_t *)U3) + r);
    uint16_t * RESTRICT U5 = (uint16_t *)(void *)(((uint8_t *)U4) + r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(((uint8_t *)U5) + r);
    uint16_t * RESTRICT U7 = (uint16_t *)(void *)(((uint8_t *)U6) + r);
    uint16_t * RESTRICT U8 = (uint16_t *)(void *)(((uint8_t *)U7) + r);
    uint16_t * RESTRICT U9 = (uint16_t *)(void *)(((uint8_t *)U8) + r);
    uint16_t * RESTRICT U10 = (uint16_t *)(void *)(((uint8_t *)U9) + r);
    uint16_t * RESTRICT U11 = (uint16_t *)(void *)(((uint8_t *)U10) + r);
    uint16_t * RESTRICT U12 = (uint16_t *)(void *)(((uint8_t *)U11) + r);
    uint16_t * RESTRICT U13 = (uint16_t *)(void *)(((uint8_t *)U12) + r);
    uint16_t * RESTRICT U14 = (uint16_t *)(void *)(((uint8_t *)U13) + r);
    uint16_t * RESTRICT U15 = (uint16_t *)(void *)(((uint8_t *)U14) + r);

    fast_uint_t i, p0 = *i0, p1 = *i1, p2 = *i2, p3 = *i3, p4 = *i4, p5 = *i5, p6 = *i6, p7 = *i7, p8 = *i8, p9 = *i9, p10 = *i10, p11 = *i11, p12 = *i12, p13 = *i13, p14 = *i14, p15 = *i15;

    for (i = 0; i != k; ++i)
    {
        uint16_t c0 = fastbits[p0 >> shift]; if (bucket2[c0] <= p0) { do { c0++; } while (bucket2[c0] <= p0); } p0 = P[p0]; U0[i] = libsais_bswap16(c0);
        uint16_t c1 = fastbits[p1 >> shift]; if (bucket2[c1] <= p1) { do { c1++; } while (bucket2[c1] <= p1); } p1 = P[p1]; U1[i] = libsais_bswap16(c1);
        uint16_t c2 = fastbits[p2 >> shift]; if (bucket2[c2] <= p2) { do { c2++; } while (bucket2[c2] <= p2); } p2 = P[p2]; U2[i] = libsais_bswap16(c2);
        uint16_t c3 = fastbits[p3 >> shift]; if (bucket2[c3] <= p3) { do { c3++; } while (bucket2[c3] <= p3); } p3 = P[p3]; U3[i] = libsais_bswap16(c3);
        uint16_t c4 = fastbits[p4 >> shift]; if (bucket2[c4] <= p4) { do { c4++; } while (bucket2[c4] <= p4); } p4 = P[p4]; U4[i] = libsais_bswap16(c4);
        uint16_t c5 = fastbits[p5 >> shift]; if (bucket2[c5] <= p5) { do { c5++; } while (bucket2[c5] <= p5); } p5 = P[p5]; U5[i] = libsais_bswap16(c5);
        uint16_t c6 = fastbits[p6 >> shift]; if (bucket2[c6] <= p6) { do { c6++; } while (bucket2[c6] <= p6); } p6 = P[p6]; U6[i] = libsais_bswap16(c6);
        uint16_t c7 = fastbits[p7 >> shift]; if (bucket2[c7] <= p7) { do { c7++; } while (bucket2[c7] <= p7); } p7 = P[p7]; U7[i] = libsais_bswap16(c7);
        uint16_t c8 = fastbits[p8 >> shift]; if (bucket2[c8] <= p8) { do { c8++; } while (bucket2[c8] <= p8); } p8 = P[p8]; U8[i] = libsais_bswap16(c8);
        uint16_t c9 = fastbits[p9 >> shift]; if (bucket2[c9] <= p9) { do { c9++; } while (bucket2[c9] <= p9); } p9 = P[p9]; U9[i] = libsais_bswap16(c9);
        uint16_t c10 = fastbits[p10 >> shift]; if (bucket2[c10] <= p10) { do { c10++; } while (bucket2[c10] <= p10); } p10 = P[p10]; U10[i] = libsais_bswap16(c10);
        uint16_t c11 = fastbits[p11 >> shift]; if (bucket2[c11] <= p11) { do { c11++; } while (bucket2[c11] <= p11); } p11 = P[p11]; U11[i] = libsais_bswap16(c11);
        uint16_t c12 = fastbits[p12 >> shift]; if (bucket2[c12] <= p12) { do { c12++; } while (bucket2[c12] <= p12); } p12 = P[p12]; U12[i] = libsais_bswap16(c12);
        uint16_t c13 = fastbits[p13 >> shift]; if (bucket2[c13] <= p13) { do { c13++; } while (bucket2[c13] <= p13); } p13 = P[p13]; U13[i] = libsais_bswap16(c13);
        uint16_t c14 = fastbits[p14 >> shift]; if (bucket2[c14] <= p14) { do { c14++; } while (bucket2[c14] <= p14); } p14 = P[p14]; U14[i] = libsais_bswap16(c14);
        uint16_t c15 = fastbits[p15 >> shift]; if (bucket2[c15] <= p15) { do { c15++; } while (bucket2[c15] <= p15); } p15 = P[p15]; U15[i] = libsais_bswap16(c15);
    }

    *i0 = p0; *i1 = p1; *i2 = p2; *i3 = p3; *i4 = p4; *i5 = p5; *i6 = p6; *i7 = p7; *i8 = p8; *i9 = p9; *i10 = p10; *i11 = p11; *i12 = p12; *i13 = p13; *i14 = p14; *i15 = p15;
}

static void libsais_unbwt_decode(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t r, const sa_uint_t * RESTRICT I, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_sint_t blocks, fast_uint_t remainder)
{
    fast_uint_t shift       = 0; while ((n >> shift) > ((sa_sint_t)1 << UNBWT_FASTBITS)) { shift++; }
    fast_uint_t offset      = 0;

    if (n >= UNBWT_DECODE_16_MIN_SIZE)
    {
        while (blocks > 16)
        {
            fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11], i12 = I[12], i13 = I[13], i14 = I[14], i15 = I[15];
            libsais_unbwt_decode_16(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, &i14, &i15, (fast_uint_t)r >> 1);
            I += 16; blocks -= 16; offset += 16 * (fast_uint_t)r;
        }
    }
    else
    {
        while (blocks > (n >= UNBWT_DECODE_12_MIN_SIZE ? 12 : 8))
        {
            fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7];
            libsais_unbwt_decode_8(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, (fast_uint_t)r >> 1);
            I += 8; blocks -= 8; offset += 8 * (fast_uint_t)r;
        }
    }

    if (blocks == 1)
    {
        fast_uint_t i0 = I[0];
        libsais_unbwt_decode_1(U + offset, P, bucket2, fastbits, shift, &i0, remainder >> 1);
    }
    else if (blocks == 2)
    {
        fast_uint_t i0 = I[0], i1 = I[1];
        libsais_unbwt_decode_2(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, remainder >> 1);
        libsais_unbwt_decode_1(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, &i0, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 3)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2];
        libsais_unbwt_decode_3(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, remainder >> 1);
        libsais_unbwt_decode_2(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 4)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3];
        libsais_unbwt_decode_4(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, remainder >> 1);
        libsais_unbwt_decode_3(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 5)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4];
        libsais_unbwt_decode_5(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, remainder >> 1);
        libsais_unbwt_decode_4(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 6)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5];
        libsais_unbwt_decode_6(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, remainder >> 1);
        libsais_unbwt_decode_5(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 7)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6];
        libsais_unbwt_decode_7(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, remainder >> 1);
        libsais_unbwt_decode_6(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 8)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7];
        libsais_unbwt_decode_8(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, remainder >> 1);
        libsais_unbwt_decode_7(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 9)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8];
        libsais_unbwt_decode_9(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, remainder >> 1);
        libsais_unbwt_decode_8(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 10)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9];
        libsais_unbwt_decode_10(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, remainder >> 1);
        libsais_unbwt_decode_9(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 11)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10];
        libsais_unbwt_decode_11(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, remainder >> 1);
        libsais_unbwt_decode_10(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 12)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11];
        libsais_unbwt_decode_12(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, remainder >> 1);
        libsais_unbwt_decode_11(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 13)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11], i12 = I[12];
        libsais_unbwt_decode_13(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, remainder >> 1);
        libsais_unbwt_decode_12(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 14)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11], i12 = I[12], i13 = I[13];
        libsais_unbwt_decode_14(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, remainder >> 1);
        libsais_unbwt_decode_13(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else if (blocks == 15)
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11], i12 = I[12], i13 = I[13], i14 = I[14];
        libsais_unbwt_decode_15(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, &i14, remainder >> 1);
        libsais_unbwt_decode_14(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
    else
    {
        fast_uint_t i0 = I[0], i1 = I[1], i2 = I[2], i3 = I[3], i4 = I[4], i5 = I[5], i6 = I[6], i7 = I[7], i8 = I[8], i9 = I[9], i10 = I[10], i11 = I[11], i12 = I[12], i13 = I[13], i14 = I[14], i15 = I[15];
        libsais_unbwt_decode_16(U + offset, P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, &i14, &i15, remainder >> 1);
        libsais_unbwt_decode_15(U + offset + 2 * (remainder >> 1), P, bucket2, fastbits, shift, (fast_uint_t)r, &i0, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9, &i10, &i11, &i12, &i13, &i14, ((fast_uint_t)r >> 1) - (remainder >> 1));
    }
}

static void libsais_unbwt_decode_omp(uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t r, const sa_uint_t * RESTRICT I, fast_uint_t lastc, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, sa_sint_t threads)
{
    fast_sint_t blocks      = 1 + (((fast_sint_t)n - 1) / (fast_sint_t)r);
    fast_uint_t remainder   = (fast_uint_t)n - ((fast_uint_t)r * ((fast_uint_t)blocks - 1));

#if defined(LIBSAIS_OPENMP)
    fast_sint_t max_threads = blocks < threads ? blocks : threads;
    #pragma omp parallel num_threads(max_threads) if(max_threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num      = omp_get_thread_num();
        fast_sint_t omp_num_threads     = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num      = 0;
        fast_sint_t omp_num_threads     = 1;
#endif

        fast_sint_t omp_block_stride    = blocks / omp_num_threads;
        fast_sint_t omp_block_remainder = blocks % omp_num_threads;
        fast_sint_t omp_block_size      = omp_block_stride + (omp_thread_num < omp_block_remainder);
        fast_sint_t omp_block_start     = omp_block_stride * omp_thread_num + (omp_thread_num < omp_block_remainder ? omp_thread_num : omp_block_remainder);

        libsais_unbwt_decode(U + r * omp_block_start, P, n, r, I + omp_block_start, bucket2, fastbits, omp_block_size, omp_thread_num < omp_num_threads - 1 ? (fast_uint_t)r : remainder);
    }

    U[n - 1] = (uint8_t)lastc;
}

static fast_uint_t libsais_unbwt_take(fast_uint_t * counter, fast_uint_t count)
{
    fast_uint_t value;

#if defined(LIBSAIS_OPENMP)
    if (omp_get_num_threads() > 1)
    {
        #pragma omp critical(libsais_unbwt_take)
        {
            value = *counter; *counter = value + count;
        }
    }
    else
#endif
    {
        value = *counter; *counter = value + count;
    }

    return value;
}

static fast_uint_t libsais_unbwt_next_chain(fast_uint_t * batch, fast_uint_t batch_size, fast_uint_t * chains_taken)
{
    if ((*batch & (batch_size - 1)) == 0) { *batch = libsais_unbwt_take(chains_taken, batch_size); }

    return (*batch)++;
}

static uint8_t * libsais_unbwt_chunk(uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, fast_uint_t chunk, fast_uint_t chunk_size)
{
    return chunk < U_chunks ? U + chunk * chunk_size : slack + (chunk - U_chunks) * chunk_size;
}

static sa_sint_t libsais_unbwt_next_chunk(LIBSAIS_UNBWT_SLOT * RESTRICT slot, uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, fast_uint_t * chunks_taken, fast_uint_t chunk_size)
{
    fast_uint_t chunk = libsais_unbwt_take(chunks_taken, 1); if (chunk >= chunks_count) { return 0; }

    links[slot->chunk] = (sa_uint_t)chunk; slot->length += chunk_size;
    slot->chunk = chunk; slot->out = libsais_unbwt_chunk(U, slack, U_chunks, chunk, chunk_size); slot->end = slot->out + chunk_size;

    return 1;
}

static void libsais_unbwt_start_chain(LIBSAIS_UNBWT_SLOT * RESTRICT slot, LIBSAIS_UNBWT_CHAIN * RESTRICT chain, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift)
{
    chain->chunk = (sa_uint_t)slot->chunk; chain->offset = (sa_uint_t)(slot->out - (slot->end - UNBWT_CHUNK_SIZE)); slot->length = -(fast_sint_t)chain->offset;

    fast_uint_t p = chain->row; uint16_t c = fastbits[p >> shift]; if (bucket2[c] <= p) { do { c++; } while (bucket2[c] <= p); }
    *(uint16_t *)(void *)slot->out = libsais_bswap16(c); slot->out += 2; slot->p = chain->next;
}

static void libsais_unbwt_start_chain_symbol_rows(LIBSAIS_UNBWT_SLOT * RESTRICT slot, LIBSAIS_UNBWT_CHAIN * RESTRICT chain, const uint8_t * RESTRICT first_symbols, fast_uint_t first_shift, const sa_uint_t * RESTRICT symbol_ends, fast_uint_t chunk_size)
{
    chain->chunk = (sa_uint_t)slot->chunk; chain->offset = (sa_uint_t)(slot->out - (slot->end - chunk_size)); slot->length = -(fast_sint_t)chain->offset;

    fast_uint_t p = chain->row, q = chain->next, c = first_symbols[p >> first_shift]; if (symbol_ends[c] <= p) { do { c++; } while (symbol_ends[c] <= p); }
    uint8_t * RESTRICT target = slot->out; target[0] = (uint8_t)c; target[1] = (uint8_t)q; slot->out = target + 2; slot->p = q >> 8;
}

static void libsais_unbwt_decode_chains(sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, LIBSAIS_UNBWT_CHAIN * RESTRICT chains, fast_uint_t chains_count, fast_uint_t batch_size, fast_uint_t * chains_taken, uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, fast_uint_t * chunks_taken)
{
    LIBSAIS_UNBWT_SLOT slots[UNBWT_CHAINS];

    fast_uint_t s, active, batch = 0;
    for (active = 0; active < UNBWT_CHAINS; ++active)
    {
        LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[active];

        if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count) { break; }
        if ((slot->chunk = libsais_unbwt_take(chunks_taken, 1)) >= chunks_count) { return; }

        slot->out = libsais_unbwt_chunk(U, slack, U_chunks, slot->chunk, UNBWT_CHUNK_SIZE); slot->end = slot->out + UNBWT_CHUNK_SIZE;
        libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
    }

    while (active == UNBWT_CHAINS)
    {
        fast_uint_t room = (fast_uint_t)(slots[0].end - slots[0].out);
        for (s = 1; s < UNBWT_CHAINS; ++s) { fast_uint_t slot_room = (fast_uint_t)(slots[s].end - slots[s].out); if (room > slot_room) { room = slot_room; } }

        fast_uint_t offset;
        for (offset = 0; offset < room; offset += 2)
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot;
            for (slot = &slots[0]; slot != &slots[UNBWT_CHAINS]; ++slot)
            {
                fast_uint_t p = slot->p, q = P[p];
                if ((q & UNBWT_MARK) == 0)
                {
                    libsais_prefetchr(&P[q]);

                    uint16_t c = fastbits[p >> shift]; if (bucket2[c] <= p) { do { c++; } while (bucket2[c] <= p); }
                    *(uint16_t *)(void *)(slot->out + offset) = libsais_bswap16(c); slot->p = q;
                }
                else
                {
                    LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                    slot->out += offset;
                    chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                    if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                    {
                        LIBSAIS_UNBWT_SLOT * RESTRICT other;
                        for (other = &slots[0]; other != &slots[UNBWT_CHAINS]; ++other) { if (other != slot) { other->out += offset + (other < slot ? 2 : 0); } }

                        active -= 1; if (slot != &slots[active]) { *slot = slots[active]; }
                        break;
                    }

                    libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
                    slot->out -= offset + 2;
                }
            }

            if (active != UNBWT_CHAINS) { break; }
        }

        if (active == UNBWT_CHAINS) { for (s = 0; s < UNBWT_CHAINS; ++s) { slots[s].out += room; } }

        for (s = 0; s < active; ++s)
        {
            if (slots[s].out == slots[s].end && !libsais_unbwt_next_chunk(&slots[s], U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }
        }
    }

    while (active > 0)
    {
        for (s = 0; s < active; )
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[s];

            fast_uint_t p = slot->p, q = P[p];
            if ((q & UNBWT_MARK) == 0)
            {
                libsais_prefetchr(&P[q]);

                uint16_t c = fastbits[p >> shift]; if (bucket2[c] <= p) { do { c++; } while (bucket2[c] <= p); }
                *(uint16_t *)(void *)slot->out = libsais_bswap16(c); slot->out += 2; slot->p = q;
            }
            else
            {
                LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                {
                    active -= 1; if (s != active) { *slot = slots[active]; } continue;
                }

                libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
            }

            if (slot->out == slot->end && !libsais_unbwt_next_chunk(slot, U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }

            s += 1;
        }
    }
}

static void libsais_unbwt_decode_chains_symbol_rows(sa_uint_t * RESTRICT P, const uint8_t * RESTRICT first_symbols, fast_uint_t first_shift, const sa_uint_t * RESTRICT symbol_ends, LIBSAIS_UNBWT_CHAIN * RESTRICT chains, fast_uint_t chains_count, fast_uint_t batch_size, fast_uint_t * chains_taken, uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, fast_uint_t * chunks_taken)
{
    LIBSAIS_UNBWT_SLOT slots[UNBWT_CHAINS];

    fast_uint_t s, active, batch = 0;
    for (active = 0; active < UNBWT_CHAINS; ++active)
    {
        LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[active];

        if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count) { break; }
        if ((slot->chunk = libsais_unbwt_take(chunks_taken, 1)) >= chunks_count) { return; }

        slot->out = libsais_unbwt_chunk(U, slack, U_chunks, slot->chunk, UNBWT_CHUNK_SIZE); slot->end = slot->out + UNBWT_CHUNK_SIZE;
        libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE);
    }

    while (active == UNBWT_CHAINS)
    {
        fast_uint_t room = (fast_uint_t)(slots[0].end - slots[0].out);
        for (s = 1; s < UNBWT_CHAINS; ++s) { fast_uint_t slot_room = (fast_uint_t)(slots[s].end - slots[s].out); if (room > slot_room) { room = slot_room; } }

        fast_uint_t offset;
        for (offset = 0; offset < room; offset += 2)
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot;
            for (slot = &slots[0]; slot != &slots[UNBWT_CHAINS]; ++slot)
            {
                fast_uint_t p = slot->p, q = P[p];
                if ((q & UNBWT_MARK) == 0)
                {
                    libsais_prefetchr(&P[q >> 8]);

                    fast_uint_t c = first_symbols[p >> first_shift]; if (symbol_ends[c] <= p) { do { c++; } while (symbol_ends[c] <= p); }
                    uint8_t * RESTRICT target = slot->out + offset; target[0] = (uint8_t)c; target[1] = (uint8_t)q; slot->p = q >> 8;
                }
                else
                {
                    LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                    slot->out += offset;
                    chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                    if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                    {
                        LIBSAIS_UNBWT_SLOT * RESTRICT other;
                        for (other = &slots[0]; other != &slots[UNBWT_CHAINS]; ++other) { if (other != slot) { other->out += offset + (other < slot ? 2 : 0); } }

                        active -= 1; if (slot != &slots[active]) { *slot = slots[active]; }
                        break;
                    }

                    libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE);
                    slot->out -= offset + 2;
                }
            }

            if (active != UNBWT_CHAINS) { break; }
        }

        if (active == UNBWT_CHAINS) { for (s = 0; s < UNBWT_CHAINS; ++s) { slots[s].out += room; } }

        for (s = 0; s < active; ++s)
        {
            if (slots[s].out == slots[s].end && !libsais_unbwt_next_chunk(&slots[s], U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }
        }
    }

    while (active > 0)
    {
        for (s = 0; s < active; )
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[s];

            fast_uint_t p = slot->p, q = P[p];
            if ((q & UNBWT_MARK) == 0)
            {
                libsais_prefetchr(&P[q >> 8]);

                fast_uint_t c = first_symbols[p >> first_shift]; if (symbol_ends[c] <= p) { do { c++; } while (symbol_ends[c] <= p); }
                uint8_t * RESTRICT target = slot->out; target[0] = (uint8_t)c; target[1] = (uint8_t)q; slot->out = target + 2; slot->p = q >> 8;
            }
            else
            {
                LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                {
                    active -= 1; if (s != active) { *slot = slots[active]; } continue;
                }

                libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE);
            }

            if (slot->out == slot->end && !libsais_unbwt_next_chunk(slot, U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }

            s += 1;
        }
    }
}

static LIBSAIS_NOINLINE void libsais_unbwt_decode_chains_symbol_rows_small(sa_uint_t * RESTRICT P, const uint8_t * RESTRICT first_symbols, fast_uint_t first_shift, const sa_uint_t * RESTRICT symbol_ends, LIBSAIS_UNBWT_CHAIN * RESTRICT chains, fast_uint_t chains_count, fast_uint_t batch_size, fast_uint_t * chains_taken, uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, fast_uint_t * chunks_taken)
{
    LIBSAIS_UNBWT_SLOT slots[UNBWT_CHAINS];

    fast_uint_t s, active, batch = 0;
    for (active = 0; active < UNBWT_CHAINS; ++active)
    {
        LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[active];

        if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count) { break; }
        if ((slot->chunk = libsais_unbwt_take(chunks_taken, 1)) >= chunks_count) { return; }

        slot->out = libsais_unbwt_chunk(U, slack, U_chunks, slot->chunk, UNBWT_CHUNK_SIZE_SMALL); slot->end = slot->out + UNBWT_CHUNK_SIZE_SMALL;
        libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE_SMALL);
    }

    while (active == UNBWT_CHAINS)
    {
        fast_uint_t room = (fast_uint_t)(slots[0].end - slots[0].out);
        for (s = 1; s < UNBWT_CHAINS; ++s) { fast_uint_t slot_room = (fast_uint_t)(slots[s].end - slots[s].out); if (room > slot_room) { room = slot_room; } }

        fast_uint_t offset;
        for (offset = 0; offset < room; offset += 2)
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot;
            for (slot = &slots[0]; slot != &slots[UNBWT_CHAINS]; ++slot)
            {
                fast_uint_t p = slot->p, q = P[p];
                if ((q & UNBWT_MARK) == 0)
                {
                    libsais_prefetchr(&P[q >> 8]);

                    fast_uint_t c = first_symbols[p >> first_shift]; if (symbol_ends[c] <= p) { do { c++; } while (symbol_ends[c] <= p); }
                    uint8_t * RESTRICT target = slot->out + offset; target[0] = (uint8_t)c; target[1] = (uint8_t)q; slot->p = q >> 8;
                }
                else
                {
                    LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                    slot->out += offset;
                    chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE_SMALL)));

                    if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                    {
                        LIBSAIS_UNBWT_SLOT * RESTRICT other;
                        for (other = &slots[0]; other != &slots[UNBWT_CHAINS]; ++other) { if (other != slot) { other->out += offset + (other < slot ? 2 : 0); } }

                        active -= 1; if (slot != &slots[active]) { *slot = slots[active]; }
                        break;
                    }

                    libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE_SMALL);
                    slot->out -= offset + 2;
                }
            }

            if (active != UNBWT_CHAINS) { break; }
        }

        if (active == UNBWT_CHAINS) { for (s = 0; s < UNBWT_CHAINS; ++s) { slots[s].out += room; } }

        for (s = 0; s < active; ++s)
        {
            if (slots[s].out == slots[s].end && !libsais_unbwt_next_chunk(&slots[s], U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE_SMALL)) { return; }
        }
    }

    while (active > 0)
    {
        for (s = 0; s < active; )
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[s];

            fast_uint_t p = slot->p, q = P[p];
            if ((q & UNBWT_MARK) == 0)
            {
                libsais_prefetchr(&P[q >> 8]);

                fast_uint_t c = first_symbols[p >> first_shift]; if (symbol_ends[c] <= p) { do { c++; } while (symbol_ends[c] <= p); }
                uint8_t * RESTRICT target = slot->out; target[0] = (uint8_t)c; target[1] = (uint8_t)q; slot->out = target + 2; slot->p = q >> 8;
            }
            else
            {
                LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE_SMALL)));

                if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                {
                    active -= 1; if (s != active) { *slot = slots[active]; } continue;
                }

                libsais_unbwt_start_chain_symbol_rows(slot, &chains[slot->chain], first_symbols, first_shift, symbol_ends, UNBWT_CHUNK_SIZE_SMALL);
            }

            if (slot->out == slot->end && !libsais_unbwt_next_chunk(slot, U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE_SMALL)) { return; }

            s += 1;
        }
    }
}

static LIBSAIS_NOINLINE void libsais_unbwt_decode_chains_large(sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, fast_uint_t shift, LIBSAIS_UNBWT_CHAIN * RESTRICT chains, fast_uint_t chains_count, fast_uint_t batch_size, fast_uint_t * chains_taken, uint8_t * RESTRICT U, uint8_t * RESTRICT slack, fast_uint_t U_chunks, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, fast_uint_t * chunks_taken)
{
    LIBSAIS_UNBWT_SLOT slots[UNBWT_CHAINS_LARGE];

    fast_uint_t s, active, batch = 0;
    for (active = 0; active < UNBWT_CHAINS_LARGE; ++active)
    {
        LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[active];

        if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count) { break; }
        if ((slot->chunk = libsais_unbwt_take(chunks_taken, 1)) >= chunks_count) { return; }

        slot->out = libsais_unbwt_chunk(U, slack, U_chunks, slot->chunk, UNBWT_CHUNK_SIZE); slot->end = slot->out + UNBWT_CHUNK_SIZE;
        libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
    }

    while (active == UNBWT_CHAINS_LARGE)
    {
        fast_uint_t room = (fast_uint_t)(slots[0].end - slots[0].out);
        for (s = 1; s < UNBWT_CHAINS_LARGE; ++s) { fast_uint_t slot_room = (fast_uint_t)(slots[s].end - slots[s].out); if (room > slot_room) { room = slot_room; } }

        fast_uint_t offset;
        for (offset = 0; offset < room; offset += 2)
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot;
            for (slot = &slots[0]; slot != &slots[UNBWT_CHAINS_LARGE]; ++slot)
            {
                fast_uint_t p = slot->p, q = P[p];
                if ((q & UNBWT_MARK) == 0)
                {
                    libsais_prefetchr(&P[q]);

                    uint16_t c = fastbits[p >> shift]; if (bucket2[c] <= p) { do { c++; } while (bucket2[c] <= p); }
                    *(uint16_t *)(void *)(slot->out + offset) = libsais_bswap16(c); slot->p = q;
                }
                else
                {
                    LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                    slot->out += offset;
                    chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                    if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                    {
                        LIBSAIS_UNBWT_SLOT * RESTRICT other;
                        for (other = &slots[0]; other != &slots[UNBWT_CHAINS_LARGE]; ++other) { if (other != slot) { other->out += offset + (other < slot ? 2 : 0); } }

                        active -= 1; if (slot != &slots[active]) { *slot = slots[active]; }
                        break;
                    }

                    libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
                    slot->out -= offset + 2;
                }
            }

            if (active != UNBWT_CHAINS_LARGE) { break; }
        }

        if (active == UNBWT_CHAINS_LARGE) { for (s = 0; s < UNBWT_CHAINS_LARGE; ++s) { slots[s].out += room; } }

        for (s = 0; s < active; ++s)
        {
            if (slots[s].out == slots[s].end && !libsais_unbwt_next_chunk(&slots[s], U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }
        }
    }

    while (active > 0)
    {
        for (s = 0; s < active; )
        {
            LIBSAIS_UNBWT_SLOT * RESTRICT slot = &slots[s];

            fast_uint_t p = slot->p, q = P[p];
            if ((q & UNBWT_MARK) == 0)
            {
                libsais_prefetchr(&P[q]);

                uint16_t c = fastbits[p >> shift]; if (bucket2[c] <= p) { do { c++; } while (bucket2[c] <= p); }
                *(uint16_t *)(void *)slot->out = libsais_bswap16(c); slot->out += 2; slot->p = q;
            }
            else
            {
                LIBSAIS_UNBWT_CHAIN * RESTRICT chain = &chains[slot->chain];

                chain->successor = (sa_uint_t)(q ^ UNBWT_MARK); chain->length = (sa_uint_t)(slot->length + (slot->out - (slot->end - UNBWT_CHUNK_SIZE)));

                if ((slot->chain = libsais_unbwt_next_chain(&batch, batch_size, chains_taken)) >= chains_count)
                {
                    active -= 1; if (s != active) { *slot = slots[active]; } continue;
                }

                libsais_unbwt_start_chain(slot, &chains[slot->chain], bucket2, fastbits, shift);
            }

            if (slot->out == slot->end && !libsais_unbwt_next_chunk(slot, U, slack, U_chunks, links, chunks_count, chunks_taken, UNBWT_CHUNK_SIZE)) { return; }

            s += 1;
        }
    }
}

static sa_sint_t libsais_unbwt_decode_split(const uint8_t * RESTRICT T, uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_sint_t n, fast_uint_t index, fast_uint_t start, fast_uint_t end, fast_uint_t gap, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, const sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, LIBSAIS_UNBWT_CHAIN * RESTRICT chains, uint8_t * RESTRICT slack, sa_uint_t * RESTRICT links, fast_uint_t chunks_count, sa_sint_t threads)
{
    fast_uint_t shift           = 0; while ((n >> shift) > ((sa_sint_t)1 << UNBWT_FASTBITS)) { shift++; }
    fast_uint_t U_chunks        = (fast_uint_t)n / libsais_unbwt_chunk_size(n);
    fast_uint_t end_lf          = marks[end >> marks_shift];
    fast_uint_t end_next        = P[end];
    fast_uint_t end_lf_next     = P[end_lf];
    fast_uint_t count           = 1;
    fast_uint_t cycle           = 0;

    {
        chains[0].row = (sa_uint_t)start; chains[0].next = P[start]; chains[0].length = 0; chains[0].successor = ~(sa_uint_t)0;

        fast_uint_t mark;
        for (mark = 0; mark <= ((fast_uint_t)n >> marks_shift); ++mark)
        {
            fast_uint_t row = libsais_unbwt_mark_row(mark, marks_shift, end);
            if (row <= (fast_uint_t)n && row != 0 && row != start && row != end_lf)
            {
                fast_uint_t lf = marks[mark];
                if (lf != start && libsais_unbwt_mark_row(lf >> marks_shift, marks_shift, end) != lf)
                {
                    chains[count].row = (sa_uint_t)row; chains[count].lf = (sa_uint_t)lf; count++;
                }
            }
        }

        fast_sint_t c;

#if defined(LIBSAIS_OPENMP)
        #pragma omp parallel for schedule(static, 64) num_threads(threads) if(threads > 1)
#endif
        for (c = 1; c < (fast_sint_t)count; ++c)
        {
            fast_uint_t row = chains[c].row, lf = chains[c].lf;

            chains[c].next = P[row]; chains[c].lf_next = P[lf]; chains[c].gap = T[row - (row > index)]; chains[c].length = 0; chains[c].successor = ~(sa_uint_t)0;
            P[row] = (sa_uint_t)(UNBWT_MARK | ((fast_uint_t)c << 1)); P[lf] = (sa_uint_t)(UNBWT_MARK | ((fast_uint_t)c << 1) | 1);
        }

        chains[count].gap = (sa_uint_t)gap;

        P[start] = UNBWT_MARK; P[end] = (sa_uint_t)(UNBWT_MARK | (count << 1)); P[end_lf] = (sa_uint_t)(UNBWT_MARK | (count << 1) | 1);
    }

    {
        fast_uint_t chains_taken = 0, chunks_taken = 0;

        const uint8_t * RESTRICT first_symbols = (const uint8_t *)(const void *)fastbits;
        const sa_uint_t * RESTRICT symbol_ends = (const sa_uint_t *)(const void *)(first_symbols + ((fast_uint_t)1 << UNBWT_FIRST_SYMBOLS));
        fast_uint_t first_shift = 0; while ((n >> first_shift) >= ((sa_sint_t)1 << UNBWT_FIRST_SYMBOLS)) { first_shift++; }

        if (n < UNBWT_CHAIN_THREADS_MIN_SIZE)
        {
            /* a team of one: outside a region omp_get_num_threads() counts the caller's team, and libsais_unbwt_take would lock */
#if defined(LIBSAIS_OPENMP)
            #pragma omp parallel num_threads(1)
#endif
            libsais_unbwt_decode_chains_symbol_rows_small(P, first_symbols, first_shift, symbol_ends, chains, count, 1, &chains_taken, U, slack, U_chunks, links, chunks_count, &chunks_taken);
        }
        else
        {
#if defined(LIBSAIS_OPENMP)
            #pragma omp parallel num_threads(threads) if(threads > 1)
#endif
            {
                if (n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
                {
                    libsais_unbwt_decode_chains_symbol_rows(P, first_symbols, first_shift, symbol_ends, chains, count, threads >= 3 ? 8 : 1, &chains_taken, U, slack, U_chunks, links, chunks_count, &chunks_taken);
                }
                else if (libsais_unbwt_chains_slots(n, threads) == UNBWT_CHAINS_LARGE)
                {
                    libsais_unbwt_decode_chains_large(P, bucket2, fastbits, shift, chains, count, threads >= 3 ? 8 : 1, &chains_taken, U, slack, U_chunks, links, chunks_count, &chunks_taken);
                }
                else
                {
                    libsais_unbwt_decode_chains(P, bucket2, fastbits, shift, chains, count, threads >= 3 ? 8 : 1, &chains_taken, U, slack, U_chunks, links, chunks_count, &chunks_taken);
                }
            }
        }
    }

    {
        fast_uint_t c, linked = 0, position = 0;
        for (c = 0; c < count; ++c) { chains[c].position = (sa_uint_t)n; }

        c = 0;
        while (c < count && linked < count && chains[c].length <= (fast_uint_t)n - position)
        {
            chains[c].position = (sa_uint_t)position; position += chains[c].length + (chains[c].successor & 1); c = chains[c].successor >> 1; linked += 1;
        }

        if (c != count || position == 0)
        {
            P[start] = (sa_uint_t)chains[0].next; P[end] = (sa_uint_t)end_next; P[end_lf] = (sa_uint_t)end_lf_next;
            for (c = 1; c < count; ++c) { P[chains[c].row] = (sa_uint_t)chains[c].next; P[chains[c].lf] = (sa_uint_t)chains[c].lf_next; }

            return 1;
        }

        cycle = position;
    }

    {
        uint8_t * RESTRICT A = (uint8_t *)(void *)P;

        fast_sint_t c;

        if (n < UNBWT_CHAIN_THREADS_MIN_SIZE)
        {
            for (c = 0; c < (fast_sint_t)count; ++c)
            {
                if (chains[c].position >= cycle) { continue; }

                uint8_t * RESTRICT target = A + chains[c].position;
                fast_uint_t chunk = chains[c].chunk, offset = chains[c].offset, length = chains[c].length;

                while (length > 0)
                {
                    fast_uint_t size = UNBWT_CHUNK_SIZE_SMALL - offset; if (size > length) { size = length; }

                    memcpy(target, libsais_unbwt_chunk(U, slack, U_chunks, chunk, UNBWT_CHUNK_SIZE_SMALL) + offset, (size_t)size);
                    target += size; length -= size; offset = 0; if (length > 0) { chunk = links[chunk]; }
                }

                if (chains[c].successor & 1) { *target = (uint8_t)chains[chains[c].successor >> 1].gap; }
            }
        }
        else
        {
#if defined(LIBSAIS_OPENMP)
            #pragma omp parallel for schedule(static, 64) num_threads(threads) if(threads > 1)
#endif
            for (c = 0; c < (fast_sint_t)count; ++c)
            {
                if (chains[c].position >= cycle) { continue; }

                uint8_t * RESTRICT target = A + chains[c].position;
                fast_uint_t chunk = chains[c].chunk, offset = chains[c].offset, length = chains[c].length;

                while (length > 0)
                {
                    fast_uint_t size = UNBWT_CHUNK_SIZE - offset; if (size > length) { size = length; }

                    memcpy(target, libsais_unbwt_chunk(U, slack, U_chunks, chunk, UNBWT_CHUNK_SIZE) + offset, (size_t)size);
                    target += size; length -= size; offset = 0; if (length > 0) { chunk = links[chunk]; }
                }

                if (chains[c].successor & 1) { *target = (uint8_t)chains[chains[c].successor >> 1].gap; }
            }
        }

        {
            fast_uint_t size; for (size = cycle; size < (fast_uint_t)n; size += size) { memcpy(A + size, A, (size_t)(size < (fast_uint_t)n - size ? size : (fast_uint_t)n - size)); }
        }

#if defined(LIBSAIS_OPENMP)
        #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
        {
#if defined(LIBSAIS_OPENMP)
            fast_sint_t omp_thread_num    = omp_get_thread_num();
            fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
            UNUSED(threads);

            fast_sint_t omp_thread_num    = 0;
            fast_sint_t omp_num_threads   = 1;
#endif
            fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
            fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
            fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

            memcpy(U + omp_block_start, A + omp_block_start, (size_t)omp_block_size);
        }
    }

    return 0;
}

static void libsais_unbwt_touch(sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (((fast_sint_t)n + 1) / omp_num_threads) & (-1024);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : (fast_sint_t)n + 1 - omp_block_start;

        fast_sint_t i; for (i = omp_block_start; i < omp_block_start + omp_block_size; i += 1024) { P[i] = 0; }
    }
}

static void libsais_unbwt_init(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, sa_sint_t n, const sa_sint_t * freq, fast_uint_t index, fast_uint_t lastc, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, sa_uint_t * buckets, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, fast_uint_t end, sa_sint_t init_threads)
{
#if defined(LIBSAIS_OPENMP)
    if (init_threads > 1 && n >= 262144)
    {
        libsais_unbwt_init_parallel(T, P, n, freq, index, lastc, bucket2, fastbits, buckets, marks, marks_shift, end, init_threads);
    }
    else
#else
    UNUSED(buckets); UNUSED(init_threads);
#endif
    {
        libsais_unbwt_init_single(T, P, n, freq, index, lastc, bucket2, fastbits, marks, marks_shift, end);
    }
}

static sa_sint_t libsais_unbwt_core(const uint8_t * RESTRICT T, uint8_t * RESTRICT U, sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t flags, const sa_sint_t * freq, sa_sint_t r, const sa_uint_t * RESTRICT I, sa_uint_t * RESTRICT bucket2, uint16_t * RESTRICT fastbits, sa_uint_t * buckets, sa_uint_t * RESTRICT marks, fast_uint_t marks_shift, LIBSAIS_UNBWT_CHAIN * chains, uint8_t * slack, sa_uint_t * links, fast_uint_t chunks_count, sa_sint_t init_threads, sa_sint_t chains_threads, sa_sint_t threads)
{
    fast_uint_t index = (flags & LIBSAIS_FLAGS_CYCLIC) ? 0             : I[0];
    fast_uint_t lastc = (flags & LIBSAIS_FLAGS_CYCLIC) ? ALPHABET_SIZE : T[0];
    fast_uint_t lastu = (flags & LIBSAIS_FLAGS_CYCLIC) ? T[I[0] - 1]   : T[0];
    fast_uint_t end   = (flags & LIBSAIS_FLAGS_CYCLIC) ? I[0]          : 0;

    libsais_unbwt_touch(P, n, init_threads);
    libsais_unbwt_init(T, P, n, freq, index, lastc, bucket2, fastbits, buckets, marks, marks_shift, end, init_threads);

    if (marks == NULL || libsais_unbwt_decode_split(T, U, P, n, index, I[0], end, lastu, bucket2, fastbits, marks, marks_shift, chains, slack, links, chunks_count, chains_threads) != 0)
    {
        if (marks != NULL && n <= UNBWT_SYMBOL_ROWS_MAX_SIZE)
        {
            libsais_unbwt_init(T, P, n, freq, index, lastc, bucket2, fastbits, buckets, NULL, 0, end, init_threads);
        }

        libsais_unbwt_decode_omp(U, P, n, r, I, lastu, bucket2, fastbits, threads);
    }

    return 0;
}

static void libsais_unbwt_scatter_packed_rows(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, sa_uint_t * RESTRICT bucket_lo, sa_uint_t * RESTRICT bucket_hi, fast_uint_t index, fast_uint_t m, fast_uint_t n, fast_uint_t second)
{
    fast_uint_t j = 0;

#if defined(LIBSAIS_UNBWT_RUN_LANES)
    for (; j + LIBSAIS_UNBWT_RUN_LANES <= m; j += LIBSAIS_UNBWT_RUN_LANES)
    {
        fast_uint_t j0 = j, j1 = m + j, r0 = j0 + (j0 >= index), r1 = j1 + (j1 >= index), k;
        fast_uint_t u0 = (j0 + LIBSAIS_UNBWT_RUN_LANES - 1 < index || j0 >= index) && libsais_unbwt_uniform_symbols(&T[j0]);
        fast_uint_t u1 = (j1 + LIBSAIS_UNBWT_RUN_LANES - 1 < index || j1 >= index) && libsais_unbwt_uniform_symbols(&T[j1]);

        if (u0) { fast_uint_t c0 = T[j0]; if (second) { libsais_unbwt_add_rows(&P[bucket_lo[c0]], &P[r0]); } else { libsais_unbwt_store_rows(&P[bucket_lo[c0]], r0); } bucket_lo[c0] += LIBSAIS_UNBWT_RUN_LANES; }
        if (u1) { fast_uint_t c1 = T[j1]; if (second) { libsais_unbwt_add_rows(&P[bucket_hi[c1]], &P[r1]); } else { libsais_unbwt_store_rows(&P[bucket_hi[c1]], r1); } bucket_hi[c1] += LIBSAIS_UNBWT_RUN_LANES; }

        if (!u0 && !u1)
        {
            for (k = 0; k < LIBSAIS_UNBWT_RUN_LANES; ++k)
            {
                fast_uint_t a = j0 + k, b = j1 + k, ra = a + (a >= index), rb = b + (b >= index), ca = T[a], cb = T[b];
                if (second) { P[bucket_lo[ca]++] += (P[ra] & 0xffff) << 16; P[bucket_hi[cb]++] += (P[rb] & 0xffff) << 16; }
                else { P[bucket_lo[ca]++] = (sa_uint_t)ra; P[bucket_hi[cb]++] = (sa_uint_t)rb; }
            }
        }
        else if (!u0)
        {
            for (k = 0; k < LIBSAIS_UNBWT_RUN_LANES; ++k) { fast_uint_t a = j0 + k, ra = a + (a >= index), ca = T[a]; if (second) { P[bucket_lo[ca]++] += (P[ra] & 0xffff) << 16; } else { P[bucket_lo[ca]++] = (sa_uint_t)ra; } }
        }
        else if (!u1)
        {
            for (k = 0; k < LIBSAIS_UNBWT_RUN_LANES; ++k) { fast_uint_t b = j1 + k, rb = b + (b >= index), cb = T[b]; if (second) { P[bucket_hi[cb]++] += (P[rb] & 0xffff) << 16; } else { P[bucket_hi[cb]++] = (sa_uint_t)rb; } }
        }
    }
#endif

    for (; j < m; ++j)
    {
        fast_uint_t a = j, b = m + j, ra = a + (a >= index), rb = b + (b >= index), ca = T[a], cb = T[b];
        if (second) { P[bucket_lo[ca]++] += (P[ra] & 0xffff) << 16; P[bucket_hi[cb]++] += (P[rb] & 0xffff) << 16; }
        else { P[bucket_lo[ca]++] = (sa_uint_t)ra; P[bucket_hi[cb]++] = (sa_uint_t)rb; }
    }

    for (j = m + m; j < n; ++j)
    {
        fast_uint_t rb = j + (j >= index), cb = T[j];
        if (second) { P[bucket_hi[cb]++] += (P[rb] & 0xffff) << 16; } else { P[bucket_hi[cb]++] = (sa_uint_t)rb; }
    }
}

static void libsais_unbwt_init_packed(const uint8_t * RESTRICT T, sa_uint_t * RESTRICT P, sa_sint_t n, const sa_sint_t * freq, fast_uint_t index)
{
    sa_uint_t count_lo[ALPHABET_SIZE], count_hi[ALPHABET_SIZE], start_lo[ALPHABET_SIZE], start_hi[ALPHABET_SIZE], bucket_lo[ALPHABET_SIZE], bucket_hi[ALPHABET_SIZE];

    fast_uint_t m = (fast_uint_t)n >> 1;

    memset(count_lo, 0, ALPHABET_SIZE * sizeof(sa_uint_t)); libsais_unbwt_compute_histogram(T, (fast_sint_t)m, count_lo);
    if (freq != NULL)
    {
        fast_uint_t c; for (c = 0; c < ALPHABET_SIZE; ++c) { count_hi[c] = (sa_uint_t)freq[c] - count_lo[c]; }
    }
    else
    {
        memset(count_hi, 0, ALPHABET_SIZE * sizeof(sa_uint_t)); libsais_unbwt_compute_histogram(T + m, (fast_sint_t)((fast_uint_t)n - m), count_hi);
    }

    {
        fast_uint_t sum, c;
        for (sum = 1, c = 0; c < ALPHABET_SIZE; ++c) { start_lo[c] = (sa_uint_t)sum; start_hi[c] = (sa_uint_t)(sum + count_lo[c]); sum += (fast_uint_t)count_lo[c] + (fast_uint_t)count_hi[c]; }
    }

    memcpy(bucket_lo, start_lo, ALPHABET_SIZE * sizeof(sa_uint_t)); memcpy(bucket_hi, start_hi, ALPHABET_SIZE * sizeof(sa_uint_t));
    libsais_unbwt_scatter_packed_rows(T, P, bucket_lo, bucket_hi, index, m, (fast_uint_t)n, 0);

    P[0] = (sa_uint_t)index;

    memcpy(bucket_lo, start_lo, ALPHABET_SIZE * sizeof(sa_uint_t)); memcpy(bucket_hi, start_hi, ALPHABET_SIZE * sizeof(sa_uint_t));
    libsais_unbwt_scatter_packed_rows(T, P, bucket_lo, bucket_hi, index, m, (fast_uint_t)n, 1);

    fast_uint_t c, x = 1;
    for (c = 0; c < ALPHABET_SIZE; ++c)
    {
        fast_uint_t e = x + (fast_uint_t)count_lo[c] + (fast_uint_t)count_hi[c];
        for (; x < e; ++x)
        {
            fast_uint_t b = P[x] >> 16;
            P[x] = (sa_uint_t)(b | ((fast_uint_t)T[b - (b >= index)] << 16) | (c << 24));
        }
    }
}

static void libsais_unbwt_decode_packed_1(uint16_t * RESTRICT U0, const sa_uint_t * RESTRICT P, fast_uint_t p0, fast_uint_t k)
{
    fast_uint_t i;
    for (i = 0; i != k; ++i)
    {
        fast_uint_t q0 = P[p0]; U0[i] = libsais_bswap16((uint16_t)(q0 >> 16)); p0 = q0 & 0xffff;
    }
}

static void libsais_unbwt_decode_packed_4(uint8_t * RESTRICT U, const sa_uint_t * RESTRICT P, fast_uint_t r, const sa_uint_t * RESTRICT I, fast_uint_t k, fast_uint_t last)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)(U + 0 * r);
    uint16_t * RESTRICT U1 = (uint16_t *)(void *)(U + 1 * r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(U + 2 * r);
    uint16_t * RESTRICT U3 = (uint16_t *)(void *)(U + 3 * r);

    fast_uint_t i, p0 = I[0], p1 = I[1], p2 = I[2], p3 = I[3];
    for (i = 0; i != last; ++i)
    {
        fast_uint_t q0 = P[p0], q1 = P[p1], q2 = P[p2], q3 = P[p3];
        U0[i] = libsais_bswap16((uint16_t)(q0 >> 16)); p0 = q0 & 0xffff;
        U1[i] = libsais_bswap16((uint16_t)(q1 >> 16)); p1 = q1 & 0xffff;
        U2[i] = libsais_bswap16((uint16_t)(q2 >> 16)); p2 = q2 & 0xffff;
        U3[i] = libsais_bswap16((uint16_t)(q3 >> 16)); p3 = q3 & 0xffff;
    }

    for (; i != k; ++i)
    {
        fast_uint_t q0 = P[p0], q1 = P[p1], q2 = P[p2];
        U0[i] = libsais_bswap16((uint16_t)(q0 >> 16)); p0 = q0 & 0xffff;
        U1[i] = libsais_bswap16((uint16_t)(q1 >> 16)); p1 = q1 & 0xffff;
        U2[i] = libsais_bswap16((uint16_t)(q2 >> 16)); p2 = q2 & 0xffff;
    }
}

static void libsais_unbwt_decode_packed_8(uint8_t * RESTRICT U, const sa_uint_t * RESTRICT P, fast_uint_t r, const sa_uint_t * RESTRICT I, fast_uint_t k, fast_uint_t last)
{
    uint16_t * RESTRICT U0 = (uint16_t *)(void *)(U + 0 * r); uint16_t * RESTRICT U1 = (uint16_t *)(void *)(U + 1 * r);
    uint16_t * RESTRICT U2 = (uint16_t *)(void *)(U + 2 * r); uint16_t * RESTRICT U3 = (uint16_t *)(void *)(U + 3 * r);
    uint16_t * RESTRICT U4 = (uint16_t *)(void *)(U + 4 * r); uint16_t * RESTRICT U5 = (uint16_t *)(void *)(U + 5 * r);
    uint16_t * RESTRICT U6 = (uint16_t *)(void *)(U + 6 * r); uint16_t * RESTRICT U7 = (uint16_t *)(void *)(U + 7 * r);

    fast_uint_t i, p0 = I[0], p1 = I[1], p2 = I[2], p3 = I[3], p4 = I[4], p5 = I[5], p6 = I[6], p7 = I[7];
    for (i = 0; i != last; ++i)
    {
        fast_uint_t q0 = P[p0], q1 = P[p1], q2 = P[p2], q3 = P[p3], q4 = P[p4], q5 = P[p5], q6 = P[p6], q7 = P[p7];
        U0[i] = libsais_bswap16((uint16_t)(q0 >> 16)); p0 = q0 & 0xffff; U1[i] = libsais_bswap16((uint16_t)(q1 >> 16)); p1 = q1 & 0xffff;
        U2[i] = libsais_bswap16((uint16_t)(q2 >> 16)); p2 = q2 & 0xffff; U3[i] = libsais_bswap16((uint16_t)(q3 >> 16)); p3 = q3 & 0xffff;
        U4[i] = libsais_bswap16((uint16_t)(q4 >> 16)); p4 = q4 & 0xffff; U5[i] = libsais_bswap16((uint16_t)(q5 >> 16)); p5 = q5 & 0xffff;
        U6[i] = libsais_bswap16((uint16_t)(q6 >> 16)); p6 = q6 & 0xffff; U7[i] = libsais_bswap16((uint16_t)(q7 >> 16)); p7 = q7 & 0xffff;
    }

    for (; i != k; ++i)
    {
        fast_uint_t q0 = P[p0], q1 = P[p1], q2 = P[p2], q3 = P[p3], q4 = P[p4], q5 = P[p5], q6 = P[p6];
        U0[i] = libsais_bswap16((uint16_t)(q0 >> 16)); p0 = q0 & 0xffff; U1[i] = libsais_bswap16((uint16_t)(q1 >> 16)); p1 = q1 & 0xffff;
        U2[i] = libsais_bswap16((uint16_t)(q2 >> 16)); p2 = q2 & 0xffff; U3[i] = libsais_bswap16((uint16_t)(q3 >> 16)); p3 = q3 & 0xffff;
        U4[i] = libsais_bswap16((uint16_t)(q4 >> 16)); p4 = q4 & 0xffff; U5[i] = libsais_bswap16((uint16_t)(q5 >> 16)); p5 = q5 & 0xffff;
        U6[i] = libsais_bswap16((uint16_t)(q6 >> 16)); p6 = q6 & 0xffff;
    }
}

static void libsais_unbwt_decode_packed(uint8_t * RESTRICT U, const sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t r, const sa_uint_t * RESTRICT I)
{
    fast_uint_t blocks = 1 + ((fast_uint_t)n - 1) / (fast_uint_t)r, remainder = (fast_uint_t)n - (fast_uint_t)r * (blocks - 1), steps = (fast_uint_t)r >> 1, t = 0;

    for (; t + 8 <= blocks; t += 8)
    {
        libsais_unbwt_decode_packed_8(U + t * (fast_uint_t)r, P, (fast_uint_t)r, I + t, steps, t + 8 == blocks ? remainder >> 1 : steps);
    }

    for (; t + 4 <= blocks; t += 4)
    {
        libsais_unbwt_decode_packed_4(U + t * (fast_uint_t)r, P, (fast_uint_t)r, I + t, steps, t + 4 == blocks ? remainder >> 1 : steps);
    }

    for (; t < blocks; ++t)
    {
        libsais_unbwt_decode_packed_1((uint16_t *)(void *)(U + t * (fast_uint_t)r), P, I[t], t + 1 == blocks ? remainder >> 1 : steps);
    }
}

static sa_sint_t libsais_unbwt_packed(const uint8_t * T, uint8_t * U, sa_uint_t * RESTRICT P, sa_sint_t n, sa_sint_t flags, const sa_sint_t * freq, sa_sint_t r, const sa_uint_t * RESTRICT I)
{
    fast_uint_t index = (flags & LIBSAIS_FLAGS_CYCLIC) ? 0           : I[0];
    fast_uint_t lastu = (flags & LIBSAIS_FLAGS_CYCLIC) ? T[I[0] - 1] : T[0];

    libsais_unbwt_init_packed(T, P, n, freq, index);
    libsais_unbwt_decode_packed(U, P, n, r, I);

    U[n - 1] = (uint8_t)lastu;

    return 0;
}

static sa_sint_t libsais_unbwt_main(const uint8_t * T, uint8_t * U, sa_uint_t * P, sa_sint_t n, sa_sint_t flags, const sa_sint_t * freq, sa_sint_t r, const sa_uint_t * I, sa_sint_t threads)
{
    if (n <= UNBWT_PACKED_MAX_SIZE) { return libsais_unbwt_packed(T, U, P, n, flags, freq, r, I); }

    fast_uint_t shift = 0; while ((n >> shift) > ((sa_sint_t)1 << UNBWT_FASTBITS)) { shift++; }

    sa_uint_t *     RESTRICT bucket2        = (sa_uint_t *)libsais_alloc_aligned(ALPHABET_SIZE * ALPHABET_SIZE * sizeof(sa_uint_t), 4096);
    uint16_t *      RESTRICT fastbits       = (uint16_t *)libsais_alloc_aligned(((size_t)1 + (size_t)(n >> shift)) * sizeof(uint16_t), 4096);
    sa_sint_t       init_threads            = n / 262144 > 8 ? n / 262144 : 8; if (init_threads > threads) { init_threads = threads; }
    size_t          buckets_size            = init_threads > 1 && n >= 262144 ? (size_t)init_threads * (2 * ALPHABET_SIZE + (ALPHABET_SIZE * ALPHABET_SIZE)) * sizeof(sa_uint_t) : 0;

    sa_uint_t *             RESTRICT marks  = NULL;
    LIBSAIS_UNBWT_CHAIN *            chains = NULL;
    uint8_t *                        slack  = NULL;
    sa_uint_t *                      links  = NULL;

    fast_uint_t marks_shift = 0, chunks_count = 0;
    size_t chains_size = 0, slack_size = 0;
    sa_sint_t chains_threads = libsais_unbwt_chains_threads(n, threads), decode_in_tables = 0;

    if (n >= (r < n ? UNBWT_CHAIN_MIN_SIZE_AUX : UNBWT_CHAIN_MIN_SIZE) && (((n - 1) / r + 1) / 8 < threads || n >= UNBWT_CHAIN_ANY_BLOCKS_MIN_SIZE || libsais_unbwt_chains_slots(n, chains_threads) == UNBWT_CHAINS_LARGE))
    {
        fast_uint_t chain_marks = libsais_unbwt_chain_marks(n, chains_threads);
        while (((fast_uint_t)n >> (marks_shift + 1)) >= chain_marks) { marks_shift++; }

        chunks_count = (fast_uint_t)n / libsais_unbwt_chunk_size(n) + libsais_unbwt_chains_slots(n, chains_threads) * (fast_uint_t)chains_threads;

        chains_size = (((size_t)(n >> marks_shift) + 2) * sizeof(LIBSAIS_UNBWT_CHAIN) + 4095) & ~(size_t)4095;
        slack_size  = (size_t)libsais_unbwt_chains_slots(n, chains_threads) * (size_t)chains_threads * libsais_unbwt_chunk_size(n);

        /* below UNBWT_CHAIN_THREADS_MIN_SIZE one thread decodes: its marks follow the first symbols and their ends in fastbits, and its
           chain records, slack chunks and chunk links lie in bucket2, which the symbol-rows decode no longer reads */
        decode_in_tables = n < UNBWT_CHAIN_THREADS_MIN_SIZE && bucket2 != NULL && fastbits != NULL
            && chains_size + slack_size + (size_t)chunks_count * sizeof(sa_uint_t) <= ALPHABET_SIZE * ALPHABET_SIZE * sizeof(sa_uint_t)
            && ((size_t)1 << UNBWT_FIRST_SYMBOLS) + ((size_t)ALPHABET_SIZE + (size_t)(n >> marks_shift) + 2) * sizeof(sa_uint_t) <= ((size_t)1 + (size_t)(n >> shift)) * sizeof(uint16_t);

        marks       = decode_in_tables
            ? (sa_uint_t *)(void *)((uint8_t *)(void *)fastbits + ((size_t)1 << UNBWT_FIRST_SYMBOLS) + ALPHABET_SIZE * sizeof(sa_uint_t))
            : (sa_uint_t *)libsais_alloc_aligned(((size_t)(n >> marks_shift) + 2) * sizeof(sa_uint_t), 4096);
    }

    /* the decode first uses its chain records, slack chunks and chunk links after the setup is done with its per-thread tables, so
       one block holds both and the decode reuses the pages the setup touched */
    size_t      decode_size = marks != NULL ? chains_size + slack_size + (size_t)chunks_count * sizeof(sa_uint_t) : 0;
    uint8_t *   work        = decode_in_tables ? (uint8_t *)(void *)bucket2 : buckets_size > 0 || decode_size > 0 ? (uint8_t *)libsais_alloc_aligned(buckets_size > decode_size ? buckets_size : decode_size, 4096) : NULL;

    if (work == NULL && buckets_size > 0 && decode_size > buckets_size) { work = (uint8_t *)libsais_alloc_aligned(buckets_size, 4096); decode_size = 0; }

    sa_uint_t * buckets     = buckets_size > 0 ? (sa_uint_t *)(void *)work : NULL;

    if (work != NULL && decode_size > 0)
    {
        chains = (LIBSAIS_UNBWT_CHAIN *)(void *)work; slack = work + chains_size; links = (sa_uint_t *)(void *)(work + chains_size + slack_size);
    }
    else
    {
        libsais_free_aligned(marks); marks = NULL;
    }

    sa_sint_t index = bucket2 != NULL && fastbits != NULL && (buckets != NULL || init_threads == 1 || n < 262144)
        ? libsais_unbwt_core(T, U, P, n, flags, freq, r, I, bucket2, fastbits, buckets, marks, marks_shift, chains, slack, links, chunks_count, init_threads, chains_threads, threads)
        : -2;

    if (!decode_in_tables)
    {
        libsais_free_aligned(work);
        libsais_free_aligned(marks);
    }

    libsais_free_aligned(fastbits);
    libsais_free_aligned(bucket2);

    return index;
}

static sa_sint_t libsais_unbwt_main_ctx(const LIBSAIS_UNBWT_CONTEXT * ctx, const uint8_t * T, uint8_t * U, sa_uint_t * P, sa_sint_t n, sa_sint_t flags, const sa_sint_t * freq, sa_sint_t r, const sa_uint_t * I)
{
    if (ctx == NULL || ctx->bucket2 == NULL || ctx->fastbits == NULL || (ctx->buckets == NULL && ctx->threads != 1))
    {
        return -2;
    }

    if (n <= UNBWT_PACKED_MAX_SIZE) { return libsais_unbwt_packed(T, U, P, n, flags, freq, r, I); }

    sa_sint_t threads = (sa_sint_t)ctx->threads, chains_threads = libsais_unbwt_chains_threads(n, threads);

    if (ctx->marks != NULL && n >= (r < n ? UNBWT_CHAIN_MIN_SIZE_AUX : UNBWT_CHAIN_MIN_SIZE) && (((n - 1) / r + 1) / 8 < threads || n >= UNBWT_CHAIN_ANY_BLOCKS_MIN_SIZE || libsais_unbwt_chains_slots(n, chains_threads) == UNBWT_CHAINS_LARGE))
    {
        fast_uint_t chain_marks = libsais_unbwt_chain_marks(n, chains_threads), marks_shift = 0;
        while (((fast_uint_t)n >> (marks_shift + 1)) >= chain_marks) { marks_shift++; }

        return libsais_unbwt_core(T, U, P, n, flags, freq, r, I, ctx->bucket2, ctx->fastbits, ctx->buckets, ctx->marks, marks_shift, ctx->chains, ctx->slack, ctx->links, (fast_uint_t)n / libsais_unbwt_chunk_size(n) + libsais_unbwt_chains_slots(n, chains_threads) * (fast_uint_t)chains_threads, threads, chains_threads, threads);
    }

    return libsais_unbwt_core(T, U, P, n, flags, freq, r, I, ctx->bucket2, ctx->fastbits, ctx->buckets, NULL, 0, NULL, NULL, NULL, 0, threads, 1, threads);
}

void * libsais_unbwt_create_ctx(void)
{
    return (void *)libsais_unbwt_create_ctx_main(1);
}

void libsais_unbwt_free_ctx(void * ctx)
{
    libsais_unbwt_free_ctx_main((LIBSAIS_UNBWT_CONTEXT *)ctx);
}

int32_t libsais_unbwt(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i)
{
    return libsais_unbwt_aux(T, U, A, n, freq, n, &i);
}

int32_t libsais_unbwt_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i)
{
    return libsais_unbwt_aux_ctx(ctx, T, U, A, n, freq, n, &i);
}

int32_t libsais_unbwt_aux(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t r, const int32_t * I)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || ((r != n) && ((r < 2) || ((r & (r - 1)) != 0))) || (I == NULL))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (I[0] != n) { return -1; }
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    fast_sint_t t; for (t = 0; t <= (n - 1) / r; ++t) { if (I[t] <= 0 || I[t] > n) { return -1; } }

    return libsais_unbwt_main(T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_NONE, freq, r, (const sa_uint_t *)I, 1);
}

int32_t libsais_unbwt_aux_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t r, const int32_t * I)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || ((r != n) && ((r < 2) || ((r & (r - 1)) != 0))) || (I == NULL))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (I[0] != n) { return -1; }
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    fast_sint_t t; for (t = 0; t <= (n - 1) / r; ++t) { if (I[t] <= 0 || I[t] > n) { return -1; } }

    return libsais_unbwt_main_ctx((const LIBSAIS_UNBWT_CONTEXT *)ctx, T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_NONE, freq, r, (const sa_uint_t *)I);
}

int32_t libsais_unbwt_cyclic(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (i < 0) || (i >= (n > 1 ? n : 1)))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    sa_uint_t index = (sa_uint_t)i + 1;

    return libsais_unbwt_main(T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_CYCLIC, freq, n, &index, 1);
}

int32_t libsais_unbwt_cyclic_ctx(const void * ctx, const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (i < 0) || (i >= (n > 1 ? n : 1)))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    sa_uint_t index = (sa_uint_t)i + 1;

    return libsais_unbwt_main_ctx((const LIBSAIS_UNBWT_CONTEXT *)ctx, T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_CYCLIC, freq, n, &index);
}

#if defined(LIBSAIS_OPENMP)

void * libsais_unbwt_create_ctx_omp(int32_t threads)
{
    if (threads < 0) { return NULL; }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return (void *)libsais_unbwt_create_ctx_main(threads);
}

int32_t libsais_unbwt_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i, int32_t threads)
{
    return libsais_unbwt_aux_omp(T, U, A, n, freq, n, &i, threads);
}

int32_t libsais_unbwt_aux_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t r, const int32_t * I, int32_t threads)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || ((r != n) && ((r < 2) || ((r & (r - 1)) != 0))) || (I == NULL) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (I[0] != n) { return -1; }
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    fast_sint_t t; for (t = 0; t <= (n - 1) / r; ++t) { if (I[t] <= 0 || I[t] > n) { return -1; } }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    return libsais_unbwt_main(T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_NONE, freq, r, (const sa_uint_t *)I, threads);
}

int32_t libsais_unbwt_cyclic_omp(const uint8_t * T, uint8_t * U, int32_t * A, int32_t n, const int32_t * freq, int32_t i, int32_t threads)
{
    if ((T == NULL) || (U == NULL) || (A == NULL) || (n < 0) || (i < 0) || (i >= (n > 1 ? n : 1)) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { U[0] = T[0]; }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    sa_uint_t index = (sa_uint_t)i + 1;

    return libsais_unbwt_main(T, U, (sa_uint_t *)A, n, LIBSAIS_FLAGS_CYCLIC, freq, n, &index, threads);
}

#endif

static void libsais_compute_phi(const sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT PLCP, sa_sint_t n, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j; sa_sint_t k = omp_block_start > 0 ? SA[omp_block_start - 1] : n;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);

        libsais_prefetchw(&PLCP[SA[i + prefetch_distance + 0]]);
        libsais_prefetchw(&PLCP[SA[i + prefetch_distance + 1]]);

        PLCP[SA[i + 0]] = k; k = SA[i + 0];
        PLCP[SA[i + 1]] = k; k = SA[i + 1];

        libsais_prefetchw(&PLCP[SA[i + prefetch_distance + 2]]);
        libsais_prefetchw(&PLCP[SA[i + prefetch_distance + 3]]);

        PLCP[SA[i + 2]] = k; k = SA[i + 2];
        PLCP[SA[i + 3]] = k; k = SA[i + 3];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        PLCP[SA[i]] = k; k = SA[i];
    }
}

static void libsais_compute_phi_omp(const sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT PLCP, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        libsais_compute_phi(SA, PLCP, n, omp_block_start, omp_block_size);
    }
}

static void libsais_compute_plcp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, fast_sint_t n, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance; i < j; i += 1)
    {
        libsais_prefetchw(&PLCP[i + 2 * prefetch_distance]);
        libsais_prefetchr(&T[PLCP[i + prefetch_distance] + l]);

        fast_sint_t k = PLCP[i], m = n - (i > k ? i : k);
        while (l < m && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }

    for (j += prefetch_distance; i < j; i += 1)
    {
        fast_sint_t k = PLCP[i], m = n - (i > k ? i : k);
        while (l < m && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }
}

static void libsais_compute_plcp_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        libsais_compute_plcp(T, PLCP, n, omp_block_start, omp_block_size);
    }
}

static void libsais_compute_plcp_gsa(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance; i < j; i += 1)
    {
        libsais_prefetchw(&PLCP[i + 2 * prefetch_distance]);
        libsais_prefetchr(&T[PLCP[i + prefetch_distance] + l]);

        fast_sint_t k = PLCP[i];
        while (T[i + l] > 0 && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }

    for (j += prefetch_distance; i < j; i += 1)
    {
        fast_sint_t k = PLCP[i];
        while (T[i + l] > 0 && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }
}

static void libsais_compute_plcp_gsa_omp(const uint8_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        libsais_compute_plcp_gsa(T, PLCP, omp_block_start, omp_block_size);
    }
}

static void libsais_compute_plcp_int(const int32_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, fast_sint_t n, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j, l = 0;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance; i < j; i += 1)
    {
        libsais_prefetchw(&PLCP[i + 2 * prefetch_distance]);
        libsais_prefetchr(&T[PLCP[i + prefetch_distance] + l]);

        fast_sint_t k = PLCP[i], m = n - (i > k ? i : k);
        while (l < m && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }

    for (j += prefetch_distance; i < j; i += 1)
    {
        fast_sint_t k = PLCP[i], m = n - (i > k ? i : k);
        while (l < m && T[i + l] == T[k + l]) { l++; }

        PLCP[i] = (sa_sint_t)l; l -= (l != 0);
    }
}

static void libsais_compute_plcp_int_omp(const int32_t * RESTRICT T, sa_sint_t * RESTRICT PLCP, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        libsais_compute_plcp_int(T, PLCP, n, omp_block_start, omp_block_size);
    }
}

static void libsais_compute_lcp(const sa_sint_t * RESTRICT PLCP, const sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT LCP, fast_sint_t omp_block_start, fast_sint_t omp_block_size)
{
    const fast_sint_t prefetch_distance = 64;

    fast_sint_t i, j;
    for (i = omp_block_start, j = omp_block_start + omp_block_size - prefetch_distance - 3; i < j; i += 4)
    {
        libsais_prefetchr(&SA[i + 2 * prefetch_distance]);
        libsais_prefetchw(&LCP[i + prefetch_distance]);

        libsais_prefetchr(&PLCP[SA[i + prefetch_distance + 0]]);
        libsais_prefetchr(&PLCP[SA[i + prefetch_distance + 1]]);

        LCP[i + 0] = PLCP[SA[i + 0]];
        LCP[i + 1] = PLCP[SA[i + 1]];

        libsais_prefetchr(&PLCP[SA[i + prefetch_distance + 2]]);
        libsais_prefetchr(&PLCP[SA[i + prefetch_distance + 3]]);

        LCP[i + 2] = PLCP[SA[i + 2]];
        LCP[i + 3] = PLCP[SA[i + 3]];
    }

    for (j += prefetch_distance + 3; i < j; i += 1)
    {
        LCP[i] = PLCP[SA[i]];
    }
}

static void libsais_compute_lcp_omp(const sa_sint_t * RESTRICT PLCP, const sa_sint_t * RESTRICT SA, sa_sint_t * RESTRICT LCP, sa_sint_t n, sa_sint_t threads)
{
#if defined(LIBSAIS_OPENMP)
    #pragma omp parallel num_threads(threads) if(threads > 1 && n >= 65536)
#endif
    {
#if defined(LIBSAIS_OPENMP)
        fast_sint_t omp_thread_num    = omp_get_thread_num();
        fast_sint_t omp_num_threads   = omp_get_num_threads();
#else
        UNUSED(threads);

        fast_sint_t omp_thread_num    = 0;
        fast_sint_t omp_num_threads   = 1;
#endif
        fast_sint_t omp_block_stride  = (n / omp_num_threads) & (-16);
        fast_sint_t omp_block_start   = omp_thread_num * omp_block_stride;
        fast_sint_t omp_block_size    = omp_thread_num < omp_num_threads - 1 ? omp_block_stride : n - omp_block_start;

        libsais_compute_lcp(PLCP, SA, LCP, omp_block_start, omp_block_size);
    }
}

int32_t libsais_plcp(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }

    libsais_compute_phi_omp(SA, PLCP, n, 1);
    libsais_compute_plcp_omp(T, PLCP, n, 1);

    return 0;
}

int32_t libsais_plcp_gsa(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0) || (n > 0 && T[n - 1] != 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }

    libsais_compute_phi_omp(SA, PLCP, n, 1);
    libsais_compute_plcp_gsa_omp(T, PLCP, n, 1);

    return 0;
}

int32_t libsais_plcp_int(const int32_t * T, const int32_t * SA, int32_t * PLCP, int32_t n)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }

    libsais_compute_phi_omp(SA, PLCP, n, 1);
    libsais_compute_plcp_int_omp(T, PLCP, n, 1);

    return 0;
}

int32_t libsais_lcp(const int32_t * PLCP, const int32_t * SA, int32_t * LCP, int32_t n)
{
    if ((PLCP == NULL) || (SA == NULL) || (LCP == NULL) || (n < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { LCP[0] = PLCP[SA[0]]; }
        return 0;
    }

    libsais_compute_lcp_omp(PLCP, SA, LCP, n, 1);

    return 0;
}

#if defined(LIBSAIS_OPENMP)

int32_t libsais_plcp_omp(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }
    
    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    libsais_compute_phi_omp(SA, PLCP, n, threads);
    libsais_compute_plcp_omp(T, PLCP, n, threads);

    return 0;
}

int32_t libsais_plcp_gsa_omp(const uint8_t * T, const int32_t * SA, int32_t * PLCP, int32_t n, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0) || (n > 0 && T[n - 1] != 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }
    
    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    libsais_compute_phi_omp(SA, PLCP, n, threads);
    libsais_compute_plcp_gsa_omp(T, PLCP, n, threads);

    return 0;
}

int32_t libsais_plcp_int_omp(const int32_t * T, const int32_t * SA, int32_t * PLCP, int32_t n, int32_t threads)
{
    if ((T == NULL) || (SA == NULL) || (PLCP == NULL) || (n < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { PLCP[0] = 0; }
        return 0;
    }
    
    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    libsais_compute_phi_omp(SA, PLCP, n, threads);
    libsais_compute_plcp_int_omp(T, PLCP, n, threads);

    return 0;
}

int32_t libsais_lcp_omp(const int32_t * PLCP, const int32_t * SA, int32_t * LCP, int32_t n, int32_t threads)
{
    if ((PLCP == NULL) || (SA == NULL) || (LCP == NULL) || (n < 0) || (threads < 0))
    {
        return -1;
    }
    else if (n <= 1)
    {
        if (n == 1) { LCP[0] = PLCP[SA[0]]; }
        return 0;
    }

    threads = threads > 0 ? threads : (omp_get_max_threads() / omp_get_num_threads());
    threads = threads > 0 ? threads : 1;

    libsais_compute_lcp_omp(PLCP, SA, LCP, n, threads);

    return 0;
}

#endif
