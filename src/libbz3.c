
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

#include "libbz3.h"

#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "sais.h"

#ifndef VERSION
    #define VERSION "unknown"
#endif

/* Block size limits and the amount of auxiliary memory allocated for a given block size. */

#define BZ3_MIN_BLOCK_SIZE KiB(65)
#define BZ3_MAX_BLOCK_SIZE MiB(511)
#define BWT_BOUND(x) (bz3_bound(x) + 128)

/* Blocks shorter than this are stored verbatim; they would not benefit from entropy coding. */
#define BZ3_STORED_THRESHOLD 64

/* Bits of the `model' byte in the block header. */
#define BZ3_MODEL_LZP 2
#define BZ3_MODEL_RLE 4

/* CRC32C (Castagnoli), with an initial value of 1 and no final xor.

   The dedicated CPU instructions are used directly if the compiler targets a CPU that has them. Otherwise, on
   x86 they are still used if the CPU the library ends up running on turns out to have them; this is decided
   once per state, so no global state is involved. Everything else gets the slicing-by-8 algorithm. */

#if defined(__SSE4_2__) && (defined(__x86_64__) || defined(__i386__))
    #include <nmmintrin.h>
    #define BZ3_CRC32C_X86
    #define BZ3_CRC32C_NATIVE
#elif defined(__ARM_FEATURE_CRC32)
    #include <arm_acle.h>
    #define BZ3_CRC32C_NATIVE
#elif (defined(__x86_64__) || defined(__i386__)) && (defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 5))
    #include <nmmintrin.h>
    #define BZ3_CRC32C_X86
    #define BZ3_CRC32C_DISPATCH __attribute__((target("sse4.2")))
#endif

#if defined(BZ3_CRC32C_X86)

    #ifndef BZ3_CRC32C_DISPATCH
        #define BZ3_CRC32C_DISPATCH
    #endif

BZ3_CRC32C_DISPATCH static u32 crc32c_x86(u32 crc, const u8 * buf, size_t size) {
    #if defined(__x86_64__)
    u64 crc64 = crc;
    for (; size >= 8; size -= 8, buf += 8) {
        u64 v;
        memcpy(&v, buf, sizeof(v));
        crc64 = _mm_crc32_u64(crc64, v);
    }
    crc = (u32)crc64;
    #endif
    for (; size; size--) crc = _mm_crc32_u8(crc, *buf++);
    return crc;
}

#endif

#if defined(BZ3_CRC32C_NATIVE)

/* Nothing to precompute. */
typedef struct {
    char unused;
} crc32c_tables;

static void crc32c_init(crc32c_tables * t) { (void)t; }

static u32 crc32sum(const crc32c_tables * t, u32 crc, const u8 * buf, size_t size) {
    (void)t;
    #if defined(BZ3_CRC32C_X86)
    return crc32c_x86(crc, buf, size);
    #else
    for (; size >= 4; size -= 4, buf += 4) crc = __crc32cw(crc, (u32)read_neutral_s32(buf));
    for (; size; size--) crc = __crc32cb(crc, *buf++);
    return crc;
    #endif
}

#else

static const u32 crc32c_table[256] = {
    0x00000000u, 0xF26B8303u, 0xE13B70F7u, 0x1350F3F4u, 0xC79A971Fu, 0x35F1141Cu, 0x26A1E7E8u, 0xD4CA64EBu, 0x8AD958CFu,
    0x78B2DBCCu, 0x6BE22838u, 0x9989AB3Bu, 0x4D43CFD0u, 0xBF284CD3u, 0xAC78BF27u, 0x5E133C24u, 0x105EC76Fu, 0xE235446Cu,
    0xF165B798u, 0x030E349Bu, 0xD7C45070u, 0x25AFD373u, 0x36FF2087u, 0xC494A384u, 0x9A879FA0u, 0x68EC1CA3u, 0x7BBCEF57u,
    0x89D76C54u, 0x5D1D08BFu, 0xAF768BBCu, 0xBC267848u, 0x4E4DFB4Bu, 0x20BD8EDEu, 0xD2D60DDDu, 0xC186FE29u, 0x33ED7D2Au,
    0xE72719C1u, 0x154C9AC2u, 0x061C6936u, 0xF477EA35u, 0xAA64D611u, 0x580F5512u, 0x4B5FA6E6u, 0xB93425E5u, 0x6DFE410Eu,
    0x9F95C20Du, 0x8CC531F9u, 0x7EAEB2FAu, 0x30E349B1u, 0xC288CAB2u, 0xD1D83946u, 0x23B3BA45u, 0xF779DEAEu, 0x05125DADu,
    0x1642AE59u, 0xE4292D5Au, 0xBA3A117Eu, 0x4851927Du, 0x5B016189u, 0xA96AE28Au, 0x7DA08661u, 0x8FCB0562u, 0x9C9BF696u,
    0x6EF07595u, 0x417B1DBCu, 0xB3109EBFu, 0xA0406D4Bu, 0x522BEE48u, 0x86E18AA3u, 0x748A09A0u, 0x67DAFA54u, 0x95B17957u,
    0xCBA24573u, 0x39C9C670u, 0x2A993584u, 0xD8F2B687u, 0x0C38D26Cu, 0xFE53516Fu, 0xED03A29Bu, 0x1F682198u, 0x5125DAD3u,
    0xA34E59D0u, 0xB01EAA24u, 0x42752927u, 0x96BF4DCCu, 0x64D4CECFu, 0x77843D3Bu, 0x85EFBE38u, 0xDBFC821Cu, 0x2997011Fu,
    0x3AC7F2EBu, 0xC8AC71E8u, 0x1C661503u, 0xEE0D9600u, 0xFD5D65F4u, 0x0F36E6F7u, 0x61C69362u, 0x93AD1061u, 0x80FDE395u,
    0x72966096u, 0xA65C047Du, 0x5437877Eu, 0x4767748Au, 0xB50CF789u, 0xEB1FCBADu, 0x197448AEu, 0x0A24BB5Au, 0xF84F3859u,
    0x2C855CB2u, 0xDEEEDFB1u, 0xCDBE2C45u, 0x3FD5AF46u, 0x7198540Du, 0x83F3D70Eu, 0x90A324FAu, 0x62C8A7F9u, 0xB602C312u,
    0x44694011u, 0x5739B3E5u, 0xA55230E6u, 0xFB410CC2u, 0x092A8FC1u, 0x1A7A7C35u, 0xE811FF36u, 0x3CDB9BDDu, 0xCEB018DEu,
    0xDDE0EB2Au, 0x2F8B6829u, 0x82F63B78u, 0x709DB87Bu, 0x63CD4B8Fu, 0x91A6C88Cu, 0x456CAC67u, 0xB7072F64u, 0xA457DC90u,
    0x563C5F93u, 0x082F63B7u, 0xFA44E0B4u, 0xE9141340u, 0x1B7F9043u, 0xCFB5F4A8u, 0x3DDE77ABu, 0x2E8E845Fu, 0xDCE5075Cu,
    0x92A8FC17u, 0x60C37F14u, 0x73938CE0u, 0x81F80FE3u, 0x55326B08u, 0xA759E80Bu, 0xB4091BFFu, 0x466298FCu, 0x1871A4D8u,
    0xEA1A27DBu, 0xF94AD42Fu, 0x0B21572Cu, 0xDFEB33C7u, 0x2D80B0C4u, 0x3ED04330u, 0xCCBBC033u, 0xA24BB5A6u, 0x502036A5u,
    0x4370C551u, 0xB11B4652u, 0x65D122B9u, 0x97BAA1BAu, 0x84EA524Eu, 0x7681D14Du, 0x2892ED69u, 0xDAF96E6Au, 0xC9A99D9Eu,
    0x3BC21E9Du, 0xEF087A76u, 0x1D63F975u, 0x0E330A81u, 0xFC588982u, 0xB21572C9u, 0x407EF1CAu, 0x532E023Eu, 0xA145813Du,
    0x758FE5D6u, 0x87E466D5u, 0x94B49521u, 0x66DF1622u, 0x38CC2A06u, 0xCAA7A905u, 0xD9F75AF1u, 0x2B9CD9F2u, 0xFF56BD19u,
    0x0D3D3E1Au, 0x1E6DCDEEu, 0xEC064EEDu, 0xC38D26C4u, 0x31E6A5C7u, 0x22B65633u, 0xD0DDD530u, 0x0417B1DBu, 0xF67C32D8u,
    0xE52CC12Cu, 0x1747422Fu, 0x49547E0Bu, 0xBB3FFD08u, 0xA86F0EFCu, 0x5A048DFFu, 0x8ECEE914u, 0x7CA56A17u, 0x6FF599E3u,
    0x9D9E1AE0u, 0xD3D3E1ABu, 0x21B862A8u, 0x32E8915Cu, 0xC083125Fu, 0x144976B4u, 0xE622F5B7u, 0xF5720643u, 0x07198540u,
    0x590AB964u, 0xAB613A67u, 0xB831C993u, 0x4A5A4A90u, 0x9E902E7Bu, 0x6CFBAD78u, 0x7FAB5E8Cu, 0x8DC0DD8Fu, 0xE330A81Au,
    0x115B2B19u, 0x020BD8EDu, 0xF0605BEEu, 0x24AA3F05u, 0xD6C1BC06u, 0xC5914FF2u, 0x37FACCF1u, 0x69E9F0D5u, 0x9B8273D6u,
    0x88D28022u, 0x7AB90321u, 0xAE7367CAu, 0x5C18E4C9u, 0x4F48173Du, 0xBD23943Eu, 0xF36E6F75u, 0x0105EC76u, 0x12551F82u,
    0xE03E9C81u, 0x34F4F86Au, 0xC69F7B69u, 0xD5CF889Du, 0x27A40B9Eu, 0x79B737BAu, 0x8BDCB4B9u, 0x988C474Du, 0x6AE7C44Eu,
    0xBE2DA0A5u, 0x4C4623A6u, 0x5F16D052u, 0xAD7D5351u
};

/* The table above, followed by the tables that advance the checksum over one to seven further zero bytes. */
typedef struct {
    u32 table[8][256];
    int use_x86; /* Whether the running CPU has the CRC32 instruction. */
} crc32c_tables;

static void crc32c_init(crc32c_tables * t) {
    t->use_x86 = 0;
    #if defined(BZ3_CRC32C_X86)
    __builtin_cpu_init();
    t->use_x86 = __builtin_cpu_supports("sse4.2") != 0;
    #endif
    memcpy(t->table[0], crc32c_table, sizeof(crc32c_table));
    for (int k = 1; k < 8; k++)
        for (int i = 0; i < 256; i++)
            t->table[k][i] = (t->table[k - 1][i] >> 8) ^ crc32c_table[t->table[k - 1][i] & 0xff];
}

static u32 crc32sum(const crc32c_tables * t, u32 crc, const u8 * buf, size_t size) {
    #if defined(BZ3_CRC32C_X86)
    if (t->use_x86) return crc32c_x86(crc, buf, size);
    #endif
    for (; size >= 8; size -= 8, buf += 8) {
        const u32 lo = crc ^ (u32)read_neutral_s32(buf), hi = (u32)read_neutral_s32(buf + 4);
        crc = t->table[7][lo & 0xff] ^ t->table[6][(lo >> 8) & 0xff] ^ t->table[5][(lo >> 16) & 0xff] ^
              t->table[4][lo >> 24] ^ t->table[3][hi & 0xff] ^ t->table[2][(hi >> 8) & 0xff] ^
              t->table[1][(hi >> 16) & 0xff] ^ t->table[0][hi >> 24];
    }
    for (; size; size--) crc = crc32c_table[(crc ^ *buf++) & 0xff] ^ (crc >> 8);
    return crc;
}

#endif

/* LZP code. */

#define LZP_DICTIONARY 18
#define LZP_DICTIONARY_SIZE ((size_t)1 << LZP_DICTIONARY)
#define LZP_MIN_MATCH 40
#define LZP_MATCH 0xf2

#define LZP_HASH(ctx) ((((ctx) >> 15) ^ (ctx) ^ ((ctx) >> 3)) & (u32)(LZP_DICTIONARY_SIZE - 1))

static inline u32 lzp_load32(const u8 * ptr) {
    u32 val;
    memcpy(&val, ptr, sizeof(val));
    return val;
}

/* The four bytes preceding `ptr', most recent one in the lowest bits. */
static inline u32 lzp_context(const u8 * ptr) {
    return ((u32)ptr[-1]) | (((u32)ptr[-2]) << 8) | (((u32)ptr[-3]) << 16) | (((u32)ptr[-4]) << 24);
}

static s32 lzp_encode_block(const u8 * RESTRICT in, const u8 * in_end, u8 * RESTRICT out, u8 * out_end,
                            s32 * RESTRICT lut) {
    const u8 * const ins = in;
    const u8 * const outs = out;
    const u8 * const out_eob = out_end - 8;
    const u8 * const in_match_end = in_end - LZP_MIN_MATCH - 32;
    const u8 * heur = in;

    for (s32 i = 0; i < 4; ++i) *out++ = *in++;

    u32 ctx = lzp_context(in);

    while (in < in_match_end && out < out_eob) {
        const u32 idx = LZP_HASH(ctx);
        const s32 val = lut[idx];
        lut[idx] = (s32)(in - ins);
        if (val > 0) {
            const u8 * ref = ins + val;
            if (lzp_load32(in + LZP_MIN_MATCH - 4) == lzp_load32(ref + LZP_MIN_MATCH - 4) &&
                lzp_load32(in) == lzp_load32(ref)) {
                if (heur > in && lzp_load32(heur) != lzp_load32(ref + (heur - in))) goto not_found;

                s32 len = 4;
                for (; in + len < in_match_end; len += sizeof(u32)) {
                    if (lzp_load32(in + len) != lzp_load32(ref + len)) break;
                }

                if (len < LZP_MIN_MATCH) {
                    if (heur < in + len) heur = in + len;
                    goto not_found;
                }

                len += in[len] == ref[len];
                len += in[len] == ref[len];
                len += in[len] == ref[len];

                in += len;
                ctx = lzp_context(in);

                *out++ = LZP_MATCH;

                len -= LZP_MIN_MATCH;
                while (len >= 254) {
                    len -= 254;
                    *out++ = 254;
                    if (out >= out_eob) break;
                }

                *out++ = (u8)len;
            } else {
            not_found:;
                const u8 next = *out++ = *in++;
                ctx = ctx << 8 | next;
                if (next == LZP_MATCH) *out++ = 255;
            }
        } else {
            ctx = (ctx << 8) | (*out++ = *in++);
        }
    }

    ctx = lzp_context(in);

    while (in < in_end && out < out_eob) {
        const u32 idx = LZP_HASH(ctx);
        const s32 val = lut[idx];
        lut[idx] = (s32)(in - ins);

        const u8 next = *out++ = *in++;
        ctx = ctx << 8 | next;
        if (next == LZP_MATCH && val > 0) *out++ = 255;
    }

    return out >= out_eob ? -1 : (s32)(out - outs);
}

static s32 lzp_decode_block(const u8 * RESTRICT in, const u8 * in_end, s32 * RESTRICT lut, u8 * RESTRICT out,
                            const u8 * out_end) {
    const u8 * const outs = out;

    for (s32 i = 0; i < 4; ++i) *out++ = *in++;

    u32 ctx = lzp_context(out);

    while (in < in_end && out < out_end) {
        const u32 idx = LZP_HASH(ctx);
        const s32 val = lut[idx];
        lut[idx] = (s32)(out - outs);
        if (*in == LZP_MATCH && val > 0) {
            in++;
            /* A match marker must be followed by a length or by an escape byte. */
            if (UNLIKELY(in == in_end)) return -1;
            if (*in != 255) {
                /* The match is clamped to the remaining output space, so stop accumulating the length
                   once it is past it. This keeps `len' from overflowing on untrusted data. */
                const s32 max_len = (s32)(out_end - out);
                s32 len = LZP_MIN_MATCH;
                while (1) {
                    if (UNLIKELY(in == in_end)) return -1;
                    if (LIKELY(len <= max_len)) len += *in;
                    if (*in++ != 254) break;
                }
                if (UNLIKELY(len > max_len)) len = max_len;

                /* `val' is the position of an earlier byte of the output, hence 0 < val < out - outs. */
                const u8 * ref = outs + val;
                if (ref + len <= out) {
                    memcpy(out, ref, (size_t)len);
                    out += len;
                } else {
                    /* The match overlaps its own output and has to be replicated byte by byte. */
                    for (const u8 * const oe = out + len; out < oe;) *out++ = *ref++;
                }

                ctx = lzp_context(out);
            } else {
                in++;
                ctx = (ctx << 8) | (*out++ = LZP_MATCH);
            }
        } else {
            ctx = (ctx << 8) | (*out++ = *in++);
        }
    }

    return (s32)(out - outs);
}

/* Returns the size of the LZP-compressed data or -1 if the data is too short or would not shrink. */
static s32 lzp_compress(const u8 * RESTRICT in, u8 * RESTRICT out, s32 n, s32 * RESTRICT lut) {
    if (n < LZP_MIN_MATCH + 32) return -1;

    memset(lut, 0, sizeof(s32) * LZP_DICTIONARY_SIZE);

    return lzp_encode_block(in, in + n, out, out + n, lut);
}

/* Decodes `n' bytes of LZP data into at most `max' bytes. Returns the amount of bytes decoded or -1. */
static s32 lzp_decompress(const u8 * RESTRICT in, u8 * RESTRICT out, s32 n, s32 max, s32 * RESTRICT lut) {
    if (n < 4 || max < 4) return -1;

    memset(lut, 0, sizeof(s32) * LZP_DICTIONARY_SIZE);

    return lzp_decode_block(in, in + n, lut, out, out + max);
}

/* RLE code. */

#define MRLE_HEADER 32

static s32 mrlec(const u8 * RESTRICT in, s32 inlen, u8 * RESTRICT out) {
    const u8 * const in_end = in + inlen;
    s32 gain[256] = { 0 };

    /* A run of n symbols occupies 2 + (n - 1) / 255 bytes when collapsed. */
    for (const u8 * ip = in; ip < in_end;) {
        const u8 c = *ip;
        const u8 * run_end = ip + 1;
        while (run_end < in_end && *run_end == c) run_end++;
        const s32 n = (s32)(run_end - ip);
        gain[c] += n - 2 - (n - 1) / 255;
        ip = run_end;
    }

    u8 * op = out;
    for (s32 i = 0; i < MRLE_HEADER; ++i) {
        u8 mask = 0;
        for (s32 j = 0; j < 8; ++j) mask |= (u8)((gain[i * 8 + j] > 0) << j);
        *op++ = mask;
    }

    for (const u8 * ip = in; ip < in_end;) {
        const u8 c = *ip;
        if (gain[c] > 0) {
            const u8 * run_end = ip + 1;
            while (run_end < in_end && *run_end == c) run_end++;
            s32 n = (s32)(run_end - ip);
            *op++ = c;
            for (; n > 255; n -= 255) *op++ = 255;
            *op++ = (u8)(n - 1);
            ip = run_end;
        } else {
            *op++ = c;
            ip++;
        }
    }

    return (s32)(op - out);
}

/* Decodes at most `inlen' bytes of RLE data into exactly `outlen' bytes. Returns nonzero on failure. */
static int mrled(const u8 * RESTRICT in, s32 inlen, u8 * RESTRICT out, s32 outlen) {
    s32 op = 0, ip = MRLE_HEADER;

    if (inlen < MRLE_HEADER) return 1;

    while (op < outlen && ip < inlen) {
        const u8 c = in[ip++];
        if (in[c >> 3] & (1 << (c & 7))) {
            const s32 max_run = outlen - op;
            s32 run = 1;
            u8 b;
            do {
                if (UNLIKELY(ip >= inlen)) return 1;
                b = in[ip++];
                if (LIKELY(run <= max_run)) run += b;
            } while (b == 255);
            if (run > max_run) run = max_run;
            memset(out + op, c, (size_t)run);
            op += run;
        } else {
            out[op++] = c;
        }
    }

    return op != outlen;
}

/* The entropy coder. */

typedef struct {
    u16 C0[256], C1[256][256], C2[512][17];
} cm_state;

/* The most bytes that coding of a single symbol can emit: up to four per bit. */
#define CM_MAX_SYMBOL_BYTES 32

static void cm_init(cm_state * s) {
    for (int i = 0; i < 256; i++) s->C0[i] = 1 << 15;
    for (int i = 0; i < 256; i++)
        for (int j = 0; j < 256; j++) s->C1[i][j] = 1 << 15;
    for (int i = 0; i < 512; i++)
        for (int k = 0; k < 17; k++) s->C2[i][k] = (u16)((k << 12) - (k == 16));
}

#define CM_PREDICT(S, ctx_, o1_, o2_, f_)                            \
    u16 * const t0##S = &s->C0[ctx_];                                \
    u16 * const t1##S = &(o1_)[ctx_];                                \
    const int p##S = ((*t0##S + *t1##S) * 7 + (o2_)[ctx_] * 2) >> 4; \
    u16 * const t2##S = &s->C2[2 * (ctx_) + (f_)][p##S >> 12];       \
    const int xl##S = t2##S[0], xh##S = t2##S[1];                    \
    const u32 pr##S = (u32)((xl##S + (((xh##S - xl##S) * (p##S & 4095)) >> 12)) * 3 + p##S)

#define CM_ADAPT(t, mask, rate) *(t) = (u16)(*(t) - (((int)*(t) - ((mask) & (65536 - (1 << (rate))))) >> (rate)))

/* Encode `size' bytes of `in'. Returns the amount of bytes written to `out', or -1 if the output would
   not fit before `out_end'. */
static s32 cm_encode(cm_state * RESTRICT s, const u8 * RESTRICT in, s32 size, u8 * RESTRICT out, const u8 * out_end) {
    const u8 * const out_start = out;
    u8 spill[CM_MAX_SYMBOL_BYTES];
    u32 high = 0xFFFFFFFF, low = 0, c1 = 0, c2 = 0, run = 0;

    for (s32 i = 0; i < size; i++) {
        const int spilling = out_end - out < CM_MAX_SYMBOL_BYTES;
        u8 * op = UNLIKELY(spilling) ? spill : out;

        run = c1 == c2 ? run + 1 : 0;
        const int f = run > 2;
        const u16 * const o2 = s->C1[c2];
        u16 * const o1 = s->C1[c1];

        const u32 c = in[i];
        u32 ctx = 1;
        for (int k = 7; k >= 0; k--) {
            CM_PREDICT(, ctx, o1, o2, f);
            const u32 mid = low + (u32)(((u64)(high - low) * pr) >> 18);
            const u32 bit = (c >> k) & 1;

            high = bit ? mid : high;
            low = bit ? low : mid + 1;

            /* Shift out the identical leading bytes. */
            while ((low ^ high) < (1 << 24)) {
                *op++ = (u8)(low >> 24);
                low <<= 8;
                high = (high << 8) + 0xFF;
            }

            const int mask = -(int)bit;
            CM_ADAPT(t0, mask, 2);
            CM_ADAPT(t1, mask, 4);
            CM_ADAPT(&t2[0], mask, 6);
            CM_ADAPT(&t2[1], mask, 6);
            ctx = ctx * 2 + bit;
        }

        c2 = c1;
        c1 = c;

        if (UNLIKELY(spilling)) {
            const size_t n = (size_t)(op - spill);
            if (n > (size_t)(out_end - out)) return -1;
            memcpy(out, spill, n);
            out += n;
        } else {
            out = op;
        }
    }

    if (out_end - out < 4) return -1;
    for (int i = 0; i < 4; i++) {
        *out++ = (u8)(low >> 24);
        low <<= 8;
    }

    return (s32)(out - out_start);
}

/* Decode `size' bytes into `out'.. */
static void cm_decode(cm_state * RESTRICT s, const u8 * RESTRICT in, const u8 * in_end, u8 * RESTRICT out, s32 size) {
    u32 high = 0xFFFFFFFF, low = 0, code = 0;
    u32 c1 = 0, run = 1;

#define CM_READ() (in < in_end ? (u32) * in++ : 0xFFFFFFFFu)

/* Decode a bit with the current prediction and adjust the counters behind it. */
#define CM_DECODE_BIT()                                          \
    const u32 mid = low + (u32)(((u64)(high - low) * pr) >> 18); \
    const u32 bit = code <= mid;                                 \
    high = bit ? mid : high;                                     \
    low = bit ? low : mid + 1;                                   \
    while ((low ^ high) < (1 << 24)) {                           \
        low <<= 8;                                               \
        high = (high << 8) + 0xFF;                               \
        code = (code << 8) + CM_READ();                          \
    }                                                            \
    const int mask = -(int)bit;                                  \
    CM_ADAPT(t0, mask, 2);                                       \
    CM_ADAPT(t1, mask, 4);                                       \
    CM_ADAPT(&t2[0], mask, 6);                                   \
    CM_ADAPT(&t2[1], mask, 6)

/* Make the prediction for the successor chosen by the bit the current one. */
#define CM_SELECT()       \
    t0 = bit ? t0B : t0A; \
    t1 = bit ? t1B : t1A; \
    t2 = bit ? t2B : t2A; \
    pr = bit ? prB : prA

    for (int i = 0; i < 4; i++) code = (code << 8) + CM_READ();

    /* The state of the model ahead of the first symbol: both preceding symbols are zero, which counts as
       the beginning of a run. */
    u16 * o1 = s->C1[0];
    const u16 * o2 = s->C1[0];
    u32 f = 0;

    /* The prediction for the first bit of the first symbol. */
    u16 *t0, *t1, *t2;
    u32 pr;
    {
        CM_PREDICT(A, 1, o1, o2, f);
        t0 = t0A, t1 = t1A, t2 = t2A, pr = prA;
    }

    for (s32 i = 0; i < size; i++) {
        u32 ctx = 1;

        for (int k = 0; k < 7; k++) {
            CM_PREDICT(A, ctx * 2, o1, o2, f);
            CM_PREDICT(B, ctx * 2 + 1, o1, o2, f);
            CM_DECODE_BIT();
            CM_SELECT();
            ctx = ctx * 2 + bit;
        }

        /* The last bit of a symbol: its successors are the first bits of the two symbols it can turn out
           to be. The symbol decoded right now becomes the second to last one for them. */
        const u32 sA = (ctx * 2) & 255, sB = sA + 1;
        const u32 runA = sA == c1 ? run + 1 : 0, runB = sB == c1 ? run + 1 : 0;
        const u32 fA = runA > 2, fB = runB > 2;
        const u16 * const o2n = s->C1[c1];
        CM_PREDICT(A, 1, s->C1[sA], o2n, fA);
        CM_PREDICT(B, 1, s->C1[sB], o2n, fB);
        CM_DECODE_BIT();
        CM_SELECT();

        c1 = bit ? sB : sA;
        run = bit ? runB : runA;
        f = bit ? fB : fA;
        o1 = s->C1[c1];
        o2 = o2n;
        out[i] = (u8)c1;
    }

#undef CM_SELECT
#undef CM_DECODE_BIT
#undef CM_READ
}

/* Public API. */

struct bz3_state {
    u8 * swap_buffer;
    s32 block_size;
    s32 *sais_array, *lzp_lut;
    cm_state * cm_state;
    s8 last_error;
    crc32c_tables crc;
};

BZIP3_API s8 bz3_last_error(struct bz3_state * state) { return state->last_error; }

BZIP3_API const char * bz3_version(void) { return VERSION; }

BZIP3_API size_t bz3_bound(size_t input_size) { return input_size + input_size / 50 + 50; }

BZIP3_API const char * bz3_strerror(struct bz3_state * state) {
    switch (state->last_error) {
        case BZ3_OK:
            return "No error";
        case BZ3_ERR_OUT_OF_BOUNDS:
            return "Data index out of bounds";
        case BZ3_ERR_BWT:
            return "Burrows-Wheeler transform failed";
        case BZ3_ERR_CRC:
            return "CRC32 check failed";
        case BZ3_ERR_MALFORMED_HEADER:
            return "Malformed header";
        case BZ3_ERR_TRUNCATED_DATA:
            return "Truncated data";
        case BZ3_ERR_DATA_TOO_BIG:
            return "Too much data";
        case BZ3_ERR_INIT:
            return "Initialization failed";
        case BZ3_ERR_DATA_SIZE_TOO_SMALL:
            return "Size of buffer `buffer_size` passed to the block decoder (bz3_decode_block) is too small. See "
                   "function docs for details.";
        default:
            return "Unknown error";
    }
}

static int bz3_valid_block_size(s32 block_size) {
    return block_size >= BZ3_MIN_BLOCK_SIZE && block_size <= BZ3_MAX_BLOCK_SIZE;
}

BZIP3_API struct bz3_state * bz3_new(s32 block_size) {
    if (!bz3_valid_block_size(block_size)) return NULL;

    struct bz3_state * state = malloc(sizeof(struct bz3_state));
    if (!state) return NULL;

    state->cm_state = malloc(sizeof(cm_state));
    state->swap_buffer = calloc(bz3_bound((size_t)block_size), 1);
    state->sais_array = malloc(BWT_BOUND((size_t)block_size) * sizeof(s32));
    state->lzp_lut = malloc(LZP_DICTIONARY_SIZE * sizeof(s32));

    if (!state->cm_state || !state->swap_buffer || !state->sais_array || !state->lzp_lut) {
        bz3_free(state);
        return NULL;
    }

    state->block_size = block_size;
    state->last_error = BZ3_OK;
    crc32c_init(&state->crc);

    return state;
}

BZIP3_API void bz3_free(struct bz3_state * state) {
    if (!state) return;
    free(state->swap_buffer);
    free(state->sais_array);
    free(state->cm_state);
    free(state->lzp_lut);
    free(state);
}

BZIP3_API s32 bz3_encode_block(struct bz3_state * state, u8 * buffer, s32 data_size) {
    if (data_size < 0 || data_size > state->block_size) {
        state->last_error = BZ3_ERR_DATA_TOO_BIG;
        return -1;
    }

    /* The caller guarantees room for this many bytes in `buffer'. */
    const size_t capacity = bz3_bound((size_t)data_size);
    const u32 crc32 = crc32sum(&state->crc, 1, buffer, (size_t)data_size);

    if (data_size < BZ3_STORED_THRESHOLD) {
        memmove(buffer + 8, buffer, (size_t)data_size);
        write_neutral_s32(buffer, (s32)crc32);
        write_neutral_s32(buffer + 4, -1);
        state->last_error = BZ3_OK;
        return data_size + 8;
    }

    /* `src' holds the data produced by the most recent stage, `dst' is free to be overwritten. */
    u8 *src = buffer, *dst = state->swap_buffer, *tmp;
    int model = 0;

    const s32 rle_size = mrlec(src, data_size, dst);
    if (rle_size < data_size) {
        tmp = src, src = dst, dst = tmp;
        data_size = rle_size;
        model |= BZ3_MODEL_RLE;
    }

    const s32 lzp_size = lzp_compress(src, dst, data_size, state->lzp_lut);
    if (lzp_size > 0 && lzp_size < data_size) {
        tmp = src, src = dst, dst = tmp;
        data_size = lzp_size;
        model |= BZ3_MODEL_LZP;
    }

    /* The suffix array construction makes use of the spare room in its work area. */
    const s32 bwt_idx = libsais_bwt(src, dst, state->sais_array, data_size,
                                    (s32)(BWT_BOUND((size_t)state->block_size) - (size_t)data_size), NULL);
    if (bwt_idx < 0) {
        state->last_error = BZ3_ERR_BWT;
        return -1;
    }

    /* Header: CRC32, BWT index, model, then the optional LZP and RLE sizes. */
    s32 header_size = 9;
    if (model & BZ3_MODEL_LZP) header_size += 4;
    if (model & BZ3_MODEL_RLE) header_size += 4;

    /* The BWT output lives in `dst', the entropy coder overwrites `src'. */
    cm_init(state->cm_state);
    const s32 coded_size = cm_encode(state->cm_state, dst, data_size, src + header_size, src + capacity);
    if (coded_size < 0) {
        state->last_error = BZ3_ERR_DATA_TOO_BIG;
        return -1;
    }

    write_neutral_s32(src, (s32)crc32);
    write_neutral_s32(src + 4, bwt_idx);
    src[8] = (u8)model;

    u8 * field = src + 9;
    if (model & BZ3_MODEL_LZP) write_neutral_s32(field, lzp_size), field += 4;
    if (model & BZ3_MODEL_RLE) write_neutral_s32(field, rle_size);

    if (src != buffer) memcpy(buffer, src, (size_t)(header_size + coded_size));

    state->last_error = BZ3_OK;
    return header_size + coded_size;
}

/* The sizes found in the header of a block that is not stored verbatim. */
struct bz3_block_header {
    s32 header_size;        /* The offset of the entropy coded data. */
    s32 lzp_size, rle_size; /* -1 if the corresponding stage was not used. */
};

/* Parse the model byte and the sizes that follow it. Returns 0 if `available' bytes are not enough. */
static int bz3_read_block_header(const u8 * block, size_t available, struct bz3_block_header * header) {
    const int model = block[8];
    header->header_size = 9;
    header->lzp_size = header->rle_size = -1;
    if (model & BZ3_MODEL_LZP) header->header_size += 4;
    if (model & BZ3_MODEL_RLE) header->header_size += 4;
    if (available < (size_t)header->header_size) return 0;

    const u8 * field = block + 9;
    if (model & BZ3_MODEL_LZP) header->lzp_size = read_neutral_s32(field), field += 4;
    if (model & BZ3_MODEL_RLE) header->rle_size = read_neutral_s32(field);
    return 1;
}

/* Data passed to the last decoding step can be one of the following:
   - original data
   - original data + LZP
   - original data + RLE
   - original data + RLE + LZP */
static int bz3_check_buffer_size(size_t buffer_size, s32 lzp_size, s32 rle_size, s32 orig_size) {
    return (lzp_size < 0 || (size_t)lzp_size <= buffer_size) && (rle_size < 0 || (size_t)rle_size <= buffer_size) &&
           (orig_size < 0 || (size_t)orig_size <= buffer_size);
}

static s32 bz3_decode_error(struct bz3_state * state, s8 error) {
    state->last_error = error;
    return -1;
}

BZIP3_API s32 bz3_decode_block(struct bz3_state * state, u8 * buffer, size_t buffer_size, s32 compressed_size,
                               s32 orig_size) {
    /* Need the minimum amount of bytes for the initial header, and the compressed data needs to fit within the
       claimed buffer size. A negative compressed size is rejected here too. */
    if (buffer_size < 9 || compressed_size < 0 || buffer_size < (size_t)compressed_size)
        return bz3_decode_error(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);

    const size_t max_size = bz3_bound((size_t)state->block_size);

    if ((size_t)compressed_size > max_size) return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    /* The original size can never exceed the block size, regardless of the block type. */
    if (orig_size > state->block_size || orig_size < 0) return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    /* Nothing past the compressed data is ever interpreted, even where the buffer extends further. */
    if (compressed_size < 8) return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    const u32 crc32 = (u32)read_neutral_s32(buffer);
    const s32 bwt_idx = read_neutral_s32(buffer + 4);

    if (bwt_idx == -1) {
        /* A block stored verbatim. */
        if (compressed_size - 8 > 64 || compressed_size - 8 != orig_size)
            return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

        memmove(buffer, buffer + 8, (size_t)orig_size);

        if (crc32sum(&state->crc, 1, buffer, (size_t)orig_size) != crc32) return bz3_decode_error(state, BZ3_ERR_CRC);

        state->last_error = BZ3_OK;
        return orig_size;
    }

    /* The model byte and the sizes that follow it must be a part of the compressed data. */
    struct bz3_block_header header;
    if (compressed_size < 9 || !bz3_read_block_header(buffer, (size_t)compressed_size, &header))
        return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    const s32 lzp_size = header.lzp_size, rle_size = header.rle_size;
    const int has_lzp = (buffer[8] & BZ3_MODEL_LZP) != 0;
    const int has_rle = (buffer[8] & BZ3_MODEL_RLE) != 0;

    if ((has_lzp && (lzp_size < 0 || (size_t)lzp_size > max_size)) ||
        (has_rle && (rle_size < 0 || (size_t)rle_size > max_size)))
        return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    /* Size that undoing BWT+CM should decompress into. */
    const s32 bwt_size = has_lzp ? lzp_size : has_rle ? rle_size : orig_size;

    if (!bz3_check_buffer_size(buffer_size, lzp_size, rle_size, orig_size))
        return bz3_decode_error(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);

    if (bwt_idx > bwt_size) return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    u8 * const swap = state->swap_buffer;

    /* Undo the entropy coding: buffer -> swap. bwt_size <= max_size, the size of the swap buffer. */
    cm_init(state->cm_state);
    cm_decode(state->cm_state, buffer + header.header_size, buffer + compressed_size, swap, bwt_size);

    /* Undo the BWT: swap -> buffer. bwt_size <= buffer_size was verified above. */
    if (bwt_size <= 1) {
        /* Degenerate transforms; any index was accepted for them historically. */
        if (bwt_size == 1) buffer[0] = swap[0];
    } else if (libsais_unbwt(swap, buffer, state->sais_array, bwt_size, NULL, bwt_idx) < 0) {
        return bz3_decode_error(state, BZ3_ERR_BWT);
    }

    /* `src' holds the data produced by the most recent stage, `dst' is free to be overwritten. */
    u8 *src = buffer, *dst = swap, *tmp;
    s32 size = bwt_size;

    if (has_lzp) {
        /* buffer -> swap, which holds max_size bytes. */
        size = lzp_decompress(src, dst, lzp_size, (s32)max_size, state->lzp_lut);
        if (size == -1) return bz3_decode_error(state, BZ3_ERR_CRC);
        /* The headers passed the validation above, but the data itself may still decode into more than fits
           into the buffer of the caller. Don't allow it to. */
        if ((size_t)size > buffer_size) return bz3_decode_error(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);
        tmp = src, src = dst, dst = tmp;
    }

    if (has_rle) {
        /* The output is capped at orig_size, which fits into either buffer: it does not exceed the block size
           and was checked against buffer_size above. */
        if (mrled(src, size, dst, orig_size)) return bz3_decode_error(state, BZ3_ERR_CRC);
        size = orig_size;
        tmp = src, src = dst, dst = tmp;
    }

    /* Callers rely on orig_size bytes of `buffer' being initialised with decoded data. */
    if (size != orig_size) return bz3_decode_error(state, BZ3_ERR_MALFORMED_HEADER);

    if (src != buffer) memcpy(buffer, src, (size_t)size);

    if (crc32 != crc32sum(&state->crc, 1, buffer, (size_t)size)) return bz3_decode_error(state, BZ3_ERR_CRC);

    state->last_error = BZ3_OK;
    return size;
}

#ifdef PTHREAD

    #include <pthread.h>

typedef struct {
    struct bz3_state * state;
    u8 * buffer;
    size_t buffer_size;
    s32 size;
    s32 orig_size;
    int decode;
    int started;
    pthread_t thread;
} bz3_job;

static void * bz3_run_job(void * arg) {
    bz3_job * job = arg;
    if (job->decode)
        bz3_decode_block(job->state, job->buffer, job->buffer_size, job->size, job->orig_size);
    else
        job->size = bz3_encode_block(job->state, job->buffer, job->size);
    return NULL;
}

/* Run all the jobs concurrently. The last one runs on the calling thread, and so does every job for which
   a thread could not be started. */
static void bz3_run_jobs(bz3_job * jobs, s32 n) {
    for (s32 i = 0; i < n - 1; i++) jobs[i].started = pthread_create(&jobs[i].thread, NULL, bz3_run_job, &jobs[i]) == 0;
    jobs[n - 1].started = 0;
    for (s32 i = n - 1; i >= 0; i--)
        if (!jobs[i].started) bz3_run_job(&jobs[i]);
    for (s32 i = 0; i < n; i++)
        if (jobs[i].started) pthread_join(jobs[i].thread, NULL);
}

BZIP3_API void bz3_encode_blocks(struct bz3_state * states[], u8 * buffers[], s32 sizes[], s32 n) {
    if (n <= 0) return;

    bz3_job * jobs = calloc((size_t)n, sizeof(bz3_job));
    if (!jobs) {
        /* Out of memory for the bookkeeping: process the blocks one by one. */
        for (s32 i = 0; i < n; i++) sizes[i] = bz3_encode_block(states[i], buffers[i], sizes[i]);
        return;
    }

    for (s32 i = 0; i < n; i++) {
        jobs[i].state = states[i];
        jobs[i].buffer = buffers[i];
        jobs[i].size = sizes[i];
    }
    bz3_run_jobs(jobs, n);
    for (s32 i = 0; i < n; i++) sizes[i] = jobs[i].size;

    free(jobs);
}

BZIP3_API void bz3_decode_blocks(struct bz3_state * states[], u8 * buffers[], size_t buffer_sizes[], s32 sizes[],
                                 s32 orig_sizes[], s32 n) {
    if (n <= 0) return;

    bz3_job * jobs = calloc((size_t)n, sizeof(bz3_job));
    if (!jobs) {
        for (s32 i = 0; i < n; i++) bz3_decode_block(states[i], buffers[i], buffer_sizes[i], sizes[i], orig_sizes[i]);
        return;
    }

    for (s32 i = 0; i < n; i++) {
        jobs[i].state = states[i];
        jobs[i].buffer = buffers[i];
        jobs[i].buffer_size = buffer_sizes[i];
        jobs[i].size = sizes[i];
        jobs[i].orig_size = orig_sizes[i];
        jobs[i].decode = 1;
    }
    bz3_run_jobs(jobs, n);

    free(jobs);
}

#endif

/* High level API implementations. */

#define BZ3_FRAME_HEADER_SIZE 13

BZIP3_API int bz3_compress(u32 block_size, const u8 * const in, u8 * out, size_t in_size, size_t * out_size) {
    if (block_size > in_size) block_size = (u32)bz3_bound(in_size);
    block_size = block_size <= BZ3_MIN_BLOCK_SIZE ? BZ3_MIN_BLOCK_SIZE : block_size;
    if (block_size > BZ3_MAX_BLOCK_SIZE) return BZ3_ERR_INIT;

    const size_t buf_max = *out_size;
    *out_size = 0;

    const size_t n_blocks = in_size / block_size + (in_size % block_size != 0);

    if (buf_max < BZ3_FRAME_HEADER_SIZE || buf_max < bz3_bound(in_size) || n_blocks > INT32_MAX)
        return BZ3_ERR_DATA_TOO_BIG;

    struct bz3_state * state = bz3_new((s32)block_size);
    u8 * compression_buf = malloc(bz3_bound(block_size));
    int result = BZ3_OK;
    if (!state || !compression_buf) {
        result = BZ3_ERR_INIT;
        goto cleanup;
    }

    memcpy(out, "BZ3v1", 5);
    write_neutral_s32(out + 5, (s32)block_size);
    write_neutral_s32(out + 9, (s32)n_blocks);
    size_t out_pos = BZ3_FRAME_HEADER_SIZE;

    /* Compress and write the blocks. */
    size_t in_offset = 0;
    for (size_t i = 0; i < n_blocks; i++) {
        const s32 size = (s32)(in_size - in_offset < block_size ? in_size - in_offset : block_size);
        memcpy(compression_buf, in + in_offset, (size_t)size);
        const s32 out_size_block = bz3_encode_block(state, compression_buf, size);
        if (out_size_block < 0) {
            result = bz3_last_error(state);
            goto cleanup;
        }
        if (buf_max - out_pos < (size_t)out_size_block + 8) {
            result = BZ3_ERR_DATA_TOO_BIG;
            goto cleanup;
        }
        write_neutral_s32(out + out_pos, out_size_block);
        write_neutral_s32(out + out_pos + 4, size);
        memcpy(out + out_pos + 8, compression_buf, (size_t)out_size_block);
        out_pos += (size_t)out_size_block + 8;
        in_offset += (size_t)size;
    }

    *out_size = out_pos;

cleanup:
    bz3_free(state);
    free(compression_buf);
    return result;
}

BZIP3_API int bz3_decompress(const uint8_t * in, uint8_t * out, size_t in_size, size_t * out_size) {
    if (in_size < BZ3_FRAME_HEADER_SIZE || memcmp(in, "BZ3v1", 5) != 0) return BZ3_ERR_MALFORMED_HEADER;

    const s32 block_size = read_neutral_s32(in + 5);
    const u32 n_blocks = (u32)read_neutral_s32(in + 9);
    in_size -= BZ3_FRAME_HEADER_SIZE;
    in += BZ3_FRAME_HEADER_SIZE;

    const size_t buf_max = *out_size;
    size_t out_pos = 0;
    *out_size = 0;

    if (!bz3_valid_block_size(block_size)) return BZ3_ERR_INIT;

    const size_t compression_buf_size = bz3_bound((size_t)block_size);
    struct bz3_state * state = bz3_new(block_size);
    u8 * compression_buf = malloc(compression_buf_size);
    int result = BZ3_OK;
    if (!state || !compression_buf) {
        result = BZ3_ERR_INIT;
        goto cleanup;
    }

    for (u32 i = 0; i < n_blocks; i++) {
        if (in_size < 8) {
            result = BZ3_ERR_MALFORMED_HEADER;
            goto cleanup;
        }
        const s32 size = read_neutral_s32(in);
        const s32 orig_size = read_neutral_s32(in + 4);
        if (size < 0 || (size_t)size > compression_buf_size) {
            result = BZ3_ERR_MALFORMED_HEADER;
            goto cleanup;
        }
        if (in_size - 8 < (size_t)size) {
            result = BZ3_ERR_TRUNCATED_DATA;
            goto cleanup;
        }
        if (orig_size < 0 || orig_size > block_size) {
            result = BZ3_ERR_MALFORMED_HEADER;
            goto cleanup;
        }
        if (buf_max - out_pos < (size_t)orig_size) {
            result = BZ3_ERR_DATA_TOO_BIG;
            goto cleanup;
        }
        memcpy(compression_buf, in + 8, (size_t)size);
        if (bz3_decode_block(state, compression_buf, compression_buf_size, size, orig_size) < 0) {
            result = bz3_last_error(state);
            goto cleanup;
        }
        memcpy(out + out_pos, compression_buf, (size_t)orig_size);
        out_pos += (size_t)orig_size;
        in += (size_t)size + 8;
        in_size -= (size_t)size + 8;
    }

cleanup:
    *out_size = out_pos;
    bz3_free(state);
    free(compression_buf);
    return result;
}

BZIP3_API size_t bz3_min_memory_needed(int32_t block_size) {
    if (!bz3_valid_block_size(block_size)) return 0;

    /* This mirrors the allocations made by bz3_new. */
    return sizeof(struct bz3_state) + sizeof(cm_state) + bz3_bound((size_t)block_size) +
           BWT_BOUND((size_t)block_size) * sizeof(s32) + LZP_DICTIONARY_SIZE * sizeof(s32);
}

BZIP3_API int bz3_orig_size_sufficient_for_decode(const u8 * block, size_t block_size, s32 orig_size) {
    /* Need at least 9 bytes for the initial header (4 bytes CRC + 4 bytes BWT index + 1 byte model). */
    if (block_size < 9) return -1;

    /* Blocks stored verbatim never need more than the original size. */
    if (read_neutral_s32(block + 4) == -1) return 1;

    struct bz3_block_header header;
    if (!bz3_read_block_header(block, block_size, &header)) return -1;

    return bz3_check_buffer_size(orig_size < 0 ? 0 : (size_t)orig_size, header.lzp_size, header.rle_size, orig_size);
}
