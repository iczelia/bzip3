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

#include "coder/lzp.h"

#define LZP_MIN_MATCH 40
#define LZP_MATCH 0xf2

#define LZP_HASH(ctx) ((((ctx) >> 15) ^ (ctx) ^ ((ctx) >> 3)) & (u32) (LZP_DICT - 1))

static INLINE u32 load32(const u8 * p) {
  u32 v;  memcpy(&v, p, sizeof v);  return v;
}

/*  The four bytes preceding `p', the most recent one in the lowest bits.  */
static INLINE u32 context(const u8 * p) {
  return p[-1] | (u32) p[-2] << 8 | (u32) p[-3] << 16 | (u32) p[-4] << 24;
}

static i32 encode(const u8 * restrict in, const u8 * in_end, u8 * restrict out,
                  u8 * out_end, i32 * restrict lut) {
  const u8 * const ins = in, * const outs = out;
  const u8 * const out_eob = out_end - 8;
  const u8 * const in_match_end = in_end - LZP_MIN_MATCH - 32;
  const u8 * heur = in;
  u32 ctx;
  int i;

  Fi(4, *out++ = *in++);
  ctx = context(in);

  while (in < in_match_end && out < out_eob) {
    const u32 idx = LZP_HASH(ctx);
    const i32 val = lut[idx];
    lut[idx] = (i32) (in - ins);
    if (val > 0) {
      const u8 * ref = ins + val;
      if (load32(in + LZP_MIN_MATCH - 4) == load32(ref + LZP_MIN_MATCH - 4) &&
          load32(in) == load32(ref)) {
        i32 len = 4;
        if (heur > in && load32(heur) != load32(ref + (heur - in)))
          goto not_found;

        for (; in + len < in_match_end; len += (i32) sizeof (u32))
          if (load32(in + len) != load32(ref + len)) break;

        if (len < LZP_MIN_MATCH) {
          if (heur < in + len) heur = in + len;
          goto not_found;
        }

        len += in[len] == ref[len];
        len += in[len] == ref[len];
        len += in[len] == ref[len];

        in += len;
        ctx = context(in);

        *out++ = LZP_MATCH;

        len -= LZP_MIN_MATCH;
        while (len >= 254) {
          len -= 254;
          *out++ = 254;
          if (out >= out_eob) break;
        }

        *out++ = (u8) len;
      } else {
        u8 next;
      not_found:
        next = *out++ = *in++;
        ctx = ctx << 8 | next;
        if (next == LZP_MATCH) *out++ = 255;
      }
    } else ctx = ctx << 8 | (*out++ = *in++);
  }

  ctx = context(in);

  while (in < in_end && out < out_eob) {
    const u32 idx = LZP_HASH(ctx);
    const i32 val = lut[idx];
    const u8 next = *out++ = *in++;
    lut[idx] = (i32) (in - 1 - ins);
    ctx = ctx << 8 | next;
    if (next == LZP_MATCH && val > 0) *out++ = 255;
  }

  return out >= out_eob ? -1 : (i32) (out - outs);
}

static i32 decode(const u8 * restrict in, const u8 * in_end, i32 * restrict lut,
                  u8 * restrict out, const u8 * out_end) {
  const u8 * const outs = out;
  u32 ctx;
  int i;

  Fi(4, *out++ = *in++);
  ctx = context(out);

  while (in < in_end && out < out_end) {
    const u32 idx = LZP_HASH(ctx);
    const i32 val = lut[idx];
    lut[idx] = (i32) (out - outs);
    if (*in == LZP_MATCH && val > 0) {
      in++;
      /*  A match marker is followed by a length or by an escape byte.  */
      if (UNLIKELY(in == in_end)) return -1;
      if (*in != 255) {
        const i32 max_len = (i32) (out_end - out);
        /*  `val' is the position of an earlier byte of the output, hence 0 < val < out - outs.  */
        const u8 * ref = outs + val;
        i32 len = LZP_MIN_MATCH;
        for (;;) {
          if (UNLIKELY(in == in_end)) return -1;
          if (LIKELY(len <= max_len)) len += *in;
          if (*in++ != 254) break;
        }
        if (UNLIKELY(len > max_len)) len = max_len;

        if (ref + len <= out) { memcpy(out, ref, (sz) len);  out += len; }
        else {
          /*  The match overlaps its own output and is replicated byte by byte.  */
          const u8 * const oe = out + len;
          while (out < oe) *out++ = *ref++;
        }

        ctx = context(out);
      } else {
        in++;
        ctx = ctx << 8 | (*out++ = LZP_MATCH);
      }
    } else ctx = ctx << 8 | (*out++ = *in++);
  }

  return (i32) (out - outs);
}

i32 bz3_lzp_encode(const u8 * in, u8 * out, i32 n, i32 * lut) {
  if (n < LZP_MIN_MATCH + 32) return -1;
  memset(lut, 0, sizeof (i32) * LZP_DICT);
  return encode(in, in + n, out, out + n, lut);
}

i32 bz3_lzp_decode(const u8 * in, u8 * out, i32 n, i32 max, i32 * lut) {
  if (n < 4 || max < 4) return -1;
  memset(lut, 0, sizeof (i32) * LZP_DICT);
  return decode(in, in + n, lut, out, out + max);
}
