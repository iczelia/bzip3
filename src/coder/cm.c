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

#include "coder/cm.h"

#define CM_MAX_SYMBOL 32

static void cm_init(bz3_cm * s) {
  int i, j, k;
  Fi(256, s->C0[i] = 1 << 15);
  Fi(256, Fj(256, s->C1[i][j] = 1 << 15));
  Fi(512, Fk(17, s->C2[i][k] = (u16) ((k << 12) - (k == 16))));
}

#define CM_PREDICT(S, ctx_, o1_, o2_, f_) \
  const int p##S = ((s->C0[ctx_] + (o1_)[ctx_]) * 7 + (o2_)[ctx_] * 2) >> 4; \
  u16 * const t2##S = &s->C2[2 * (ctx_) + (f_)][p##S >> 12]; \
  const int xl##S = t2##S[0], xh##S = t2##S[1]; \
  const u32 pr##S = (u32) ((xl##S + (((xh##S - xl##S) * (p##S & 4095)) >> 12)) * 3 + p##S)

/*  Move a counter towards the bit; `mask' is all ones for a set bit.  */
#define CM_ADAPT(t, mask, rate) \
  *(t) = (u16) (*(t) - (((int) *(t) - ((mask) & (65536 - (1 << (rate))))) >> (rate)))

static i32 encode(bz3_cm * restrict s, const u8 * restrict in, i32 n,
                  u8 * restrict out, const u8 * out_end) {
  const u8 * const out_start = out;
  u8 spill[CM_MAX_SYMBOL];
  u32 high = 0xFFFFFFFF, low = 0, c1 = 0, c2 = 0, run = 0;
  i32 i;

  for (i = 0; i < n; i++) {
    /*  Close to the end of the output a symbol is coded aside and copied over if it fits.  */
    const int spilling = out_end - out < CM_MAX_SYMBOL;
    u8 * op = UNLIKELY(spilling) ? spill : out;
    const u32 c = in[i];
    const u16 * o2;
    u16 * o1;
    u32 ctx = 1;
    int f, k;

    run = c1 == c2 ? run + 1 : 0;
    f = run > 2;
    o2 = s->C1[c2];
    o1 = s->C1[c1];

    for (k = 7; k >= 0; k--) {
      u16 * const t0 = &s->C0[ctx], * const t1 = &o1[ctx];
      CM_PREDICT(, ctx, o1, o2, f);
      const u32 mid = low + (u32) (((uint64_t) (high - low) * pr) >> 18);
      const u32 bit = (c >> k) & 1;
      const int mask = -(int) bit;

      high = bit ? mid : high;
      low = bit ? low : mid + 1;

      /*  Shift out the identical leading bytes.  */
      while ((low ^ high) < (1 << 24)) {
        *op++ = (u8) (low >> 24);
        low <<= 8;
        high = (high << 8) + 0xFF;
      }

      CM_ADAPT(t0, mask, 2);
      CM_ADAPT(t1, mask, 4);
      CM_ADAPT(&t2[0], mask, 6);
      CM_ADAPT(&t2[1], mask, 6);
      ctx = ctx * 2 + bit;
    }

    c2 = c1;
    c1 = c;

    if (UNLIKELY(spilling)) {
      const sz m = (sz) (op - spill);
      if (m > (sz) (out_end - out)) return -1;
      memcpy(out, spill, m);
      out += m;
    } else out = op;
  }

  if (out_end - out < 4) return -1;
  Fi(4, *out++ = (u8) (low >> 24);  low <<= 8);

  return (i32) (out - out_start);
}

i32 bz3_cm_encode(bz3_cm * s, const u8 * in, i32 n, u8 * out,
                  const u8 * out_end) {
  cm_init(s);
  return encode(s, in, n, out, out_end);
}

#define CM_READ() (in < in_end ? (u32) *in++ : 0xFFFFFFFFu)

/*  Decode a bit with the current prediction and adjust the counters behind it.  */
#define CM_DECODE_BIT() \
  const u32 mid = low + (u32) (((uint64_t) (high - low) * pr) >> 18); \
  const u32 bit = code <= mid; \
  high = bit ? mid : high; \
  low = bit ? low : mid + 1; \
  while ((low ^ high) < (1 << 24)) { \
    low <<= 8; \
    high = (high << 8) + 0xFF; \
    code = (code << 8) + CM_READ(); \
  } \
  CM_ADAPT(t0, -(int) bit, 2); \
  CM_ADAPT(t1, -(int) bit, 4); \
  CM_ADAPT(&t2[0], -(int) bit, 6); \
  CM_ADAPT(&t2[1], -(int) bit, 6)

#define CM_SELECT(ctx_, o1_) \
  t0 = &s->C0[ctx_]; \
  t1 = &(o1_)[ctx_]; \
  t2 = bit ? t2B : t2A; \
  pr = bit ? prB : prA

void bz3_cm_decode(bz3_cm * restrict s, const u8 * restrict in,
                   const u8 * in_end, u8 * restrict out, i32 n) {
  u32 high = 0xFFFFFFFF, low = 0, code = 0, c1 = 0, run = 1, f = 0, pr;
  /*  The model ahead of the first symbol: both preceding symbols are zero, which counts as a run.  */
  u16 * o1 = s->C1[0], * t0, * t1, * t2;
  const u16 * o2 = s->C1[0];
  i32 i;
  int k;

  cm_init(s);

  Fk(4, code = (code << 8) + CM_READ());

  /*  The prediction for the first bit of the first symbol.  */
  { CM_PREDICT(A, 1, o1, o2, f);
    t0 = &s->C0[1];  t1 = &o1[1];  t2 = t2A;  pr = prA; }

  for (i = 0; i < n; i++) {
    u32 ctx = 1;

    for (k = 0; k < 7; k++) {
      CM_PREDICT(A, ctx * 2, o1, o2, f);
      CM_PREDICT(B, ctx * 2 + 1, o1, o2, f);
      CM_DECODE_BIT();
      ctx = ctx * 2 + bit;
      CM_SELECT(ctx, o1);
    }

    /*  The last bit of a symbol: its successors are the first bits of the two symbols it can turn out to
        be.  The symbol decoded right now becomes the second to last one for them.  */
    { const u32 sA = (ctx * 2) & 255, sB = sA + 1;
      const u32 runA = sA == c1 ? run + 1 : 0, runB = sB == c1 ? run + 1 : 0;
      const u32 fA = runA > 2, fB = runB > 2;
      const u16 * const o2n = s->C1[c1];
      CM_PREDICT(A, 1, s->C1[sA], o2n, fA);
      CM_PREDICT(B, 1, s->C1[sB], o2n, fB);
      CM_DECODE_BIT();

      c1 = bit ? sB : sA;
      run = bit ? runB : runA;
      f = bit ? fB : fA;
      o1 = s->C1[c1];
      o2 = o2n;
      CM_SELECT(1, o1);
      out[i] = (u8) c1; }
  }
}
