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

#include "t_harness.h"
#include "sys/cpu.h"
#include "sys/crc.h"
#include "coder/bwt.h"
#include "coder/cm.h"
#include "coder/lzp.h"
#include "coder/rle.h"

/*  One bit at a time, straight from the definition.  */
static u32 crc_ref(u32 c, const u8 * p, sz n) {
  int k;
  while (n--) {
    c ^= *p++;
    Fk(8, c = c >> 1 ^ (0x82F63B78u & (0u - (c & 1))));
  }
  return c;
}

static void t_crc(void) {
  enum { N = 3 * CRC_LONG + 3 * CRC_SHORT + 100 };
  bz3_crc * t = xt_alloc(sizeof *t), * best = xt_alloc(sizeof *best);
  u8 * buf = xt_alloc(N + 16);
  xt_rng r;
  sz n;
  int i, k;
  xt_section_begin("crc");
  xt_seed(&r, 1);
  Fi(N + 16, buf[i] = (u8) xt_next(&r, 256));

  CHECK(bz3_crc_impl(t, 0) && !strcmp(t->name, "portable"),
        "the portable CRC is always there");
  CHECK(!bz3_crc_impl(t, 2), "there are two implementations at most");
  bz3_crc_init(best);
  xt_trace("crc: %s dispatched", best->name);
  CHECK(!strcmp(best->name, "portable") == !bz3_crc_impl(t, 1) &&
        (bz3_cpu_crc() || !strcmp(best->name, "portable")),
        "the dispatched CRC %s follows the CPU", best->name);

  Fi(2,
    if (!bz3_crc_impl(t, i)) continue;
    CHECK((bz3_crc32c(t, ~0u, "123456789", 9) ^ ~0u) == 0xE3069283u,
          "%s: the check value of CRC-32C", t->name);
    CHECK(bz3_crc32c(t, 1, buf, 0) == 1, "%s: no data leaves the register alone", t->name);
    /*  Every alignment and every tail, around the lane sizes too.  */
    Fk(9,
      for (n = 0; n < 300; n++)
        CHECK(bz3_crc32c(t, 1, buf + k, n) == crc_ref(1, buf + k, n),
              "%s: %lu bytes at offset %d", t->name, (unsigned long) n, k));
    for (n = 3 * CRC_SHORT - 9; n < N; n += n < 3 * CRC_SHORT + 40 ? 1 : 997)
      CHECK(bz3_crc32c(t, 0xDEADBEEF, buf + 3, n) == crc_ref(0xDEADBEEF, buf + 3, n),
            "%s: %lu bytes", t->name, (unsigned long) n);
    for (n = 3 * CRC_LONG - 20; n < 3 * CRC_LONG + 3 * CRC_SHORT + 20; n += 7)
      CHECK(bz3_crc32c(t, 7, buf + 5, n) == crc_ref(7, buf + 5, n),
            "%s: %lu bytes", t->name, (unsigned long) n);
    CHECK(bz3_crc32c(t, 1, buf + 1, N) == crc_ref(1, buf + 1, N),
          "%s: %d bytes", t->name, (int) N);
    /*  A checksum continues where it left off.  */
    CHECK(bz3_crc32c(t, bz3_crc32c(t, 1, buf, 12345), buf + 12345, N - 12345) ==
          bz3_crc32c(t, 1, buf, N), "%s: split in two", t->name));

  free(t);  free(best);  free(buf);
}

static void t_cm(void) {
  enum { N = 40000 };
  bz3_cm * s = xt_alloc(sizeof *s);
  u8 * in = xt_alloc(N), * out = xt_alloc(N + N / 50 + 64), * back = xt_alloc(N);
  xt_rng r;
  int kind;
  xt_section_begin("entropy coder");
  xt_seed(&r, 2);
  for (kind = 0; kind < 5; kind++) {
    i32 sizes[] = { 0, 1, 2, 100, N }, packed, cut;
    sz i;
    Fi(ARRAY_LEN(sizes),
      i32 n = sizes[i];
      xt_fill(&r, in, (sz) n, kind);
      packed = bz3_cm_encode(s, in, n, out, out + N + N / 50 + 64);
      CHECK(packed >= 4, "kind %d, %d bytes: coded", kind, (int) n);
      if (packed < 4) continue;
      memset(back, 0xAA, N);
      bz3_cm_decode(s, out, out + packed, back, n);
      CHECK(!memcmp(in, back, (sz) n), "kind %d, %d bytes: round trip", kind, (int) n);
      /*  The coder never writes past the end it was given, and says so when the data does not fit.  */
      for (cut = packed > 40 ? packed - 40 : 0; cut < packed; cut++) {
        out[cut] = 0x5A;
        CHECK(bz3_cm_encode(s, in, n, out, out + cut) == -1 && out[cut] == 0x5A,
              "kind %d, %d bytes: %d bytes of room refused", kind, (int) n, (int) cut);
      }
      /*  Truncated input decodes into something, without reading past it.  */
      packed = bz3_cm_encode(s, in, n, out, out + N + N / 50 + 64);
      bz3_cm_decode(s, out, out + packed / 2, back, n));
  }
  free(s);  free(in);  free(out);  free(back);
}

static void t_lzp(void) {
  enum { N = 200000 };
  u8 * in = xt_alloc(N), * out = xt_alloc(N), * back = xt_alloc(N);
  i32 * lut = xt_alloc(LZP_DICT * sizeof *lut);
  xt_rng r;
  i32 packed, n, cut;
  int kind;
  xt_section_begin("lzp");
  xt_seed(&r, 3);
  CHECK(bz3_lzp_encode(in, out, 71, lut) == -1, "input too short to code");
  CHECK(bz3_lzp_decode(in, out, 3, N, lut) == -1, "input too short to decode");
  for (kind = 1; kind < 5; kind++) {
    xt_fill(&r, in, N, kind);
    /*  The escape byte among the data.  */
    for (n = 500; n < N; n += 7919) in[n] = 0xf2;
    for (n = 72; ; n = MIN(N, n < 1000 ? n + 131 : n * 3)) {
      packed = bz3_lzp_encode(in, out, n, lut);
      CHECK(packed == -1 || (packed > 0 && packed < n),
            "kind %d, %d bytes: output is smaller or refused", kind, (int) n);
      if (packed < 0) { if (n == N) break;  continue; }
      CHECK(bz3_lzp_decode(out, back, packed, N, lut) == n && !memcmp(in, back, (sz) n),
            "kind %d, %d bytes: round trip", kind, (int) n);
      /*  The output space is never exceeded, whatever the data asks for.  */
      back[n / 2] = 0x33;
      CHECK(bz3_lzp_decode(out, back, packed, n / 2, lut) <= n / 2 && back[n / 2] == 0x33,
            "kind %d, %d bytes: output capped", kind, (int) n);
      for (cut = 4; cut < packed; cut += packed / 13 + 1)
        CHECK(bz3_lzp_decode(out, back, cut, N, lut) <= n,
              "kind %d, %d bytes: cut at %d", kind, (int) n, (int) cut);
      if (n == N) break;
    }
  }
  free(in);  free(out);  free(back);  free(lut);
}

static void t_rle(void) {
  u8 in[64], out[600];
  int i;
  xt_section_begin("rle");
  memset(in, 0, sizeof in);
  in['a' >> 3] |= 1 << ('a' & 7);     /*  Runs of 'a' carry a length.  */
  in[32] = 'x';
  in[33] = 'a';  in[34] = 4;          /*  1 + 4 times 'a'.  */
  in[35] = 'y';
  in[36] = 'a';  in[37] = 255;  in[38] = 10;  /*  1 + 255 + 10.  */
  CHECK(!bz3_rle_decode(in, 39, out, 273), "a valid stream decodes");
  CHECK(out[0] == 'x' && out[1] == 'a' && out[5] == 'a' && out[6] == 'y' &&
        out[7] == 'a' && out[272] == 'a', "into the right bytes");
  CHECK(bz3_rle_decode(in, 39, out, 274), "one byte more than the stream holds is refused");
  CHECK(bz3_rle_decode(in, 31, out, 0), "a stream without the header is refused");
  CHECK(bz3_rle_decode(in, 34, out, 6), "a run without its length is refused");
  CHECK(bz3_rle_decode(in, 38, out, 273), "a run cut inside its length is refused");
  /*  A run longer than the room left is clamped to it.  */
  Fi(273,
    out[i] = 0x11;
    CHECK(!bz3_rle_decode(in, 39, out, i) && out[i] == 0x11, "no write past %d bytes of output", i));
}

static void t_bwt(void) {
  enum { N = 70000 };
  u8 * in = xt_alloc(N), * out = xt_alloc(N), * back = xt_alloc(N);
  i32 * sa = xt_alloc((N + 200) * sizeof *sa);
  xt_rng r;
  i32 idx, n;
  int kind;
  xt_section_begin("bwt");
  xt_seed(&r, 4);
  for (kind = 0; kind < 5; kind++)
    for (n = 2; n <= N; n = n * 5 + 3 > N && n < N ? N : n * 5 + 3) {
      xt_fill(&r, in, (sz) n, kind);
      idx = bz3_bwt(in, out, sa, n, 200);
      CHECK(idx > 0 && idx <= n, "kind %d, %d bytes: transformed", kind, (int) n);
      CHECK(bz3_unbwt(out, back, sa, n, idx) == 0 && !memcmp(in, back, (sz) n),
            "kind %d, %d bytes: round trip", kind, (int) n);
      CHECK(bz3_unbwt(out, back, sa, n, 0) < 0 && bz3_unbwt(out, back, sa, n, n + 1) < 0,
            "kind %d, %d bytes: an index outside the data is refused", kind, (int) n);
    }
  free(in);  free(out);  free(back);  free(sa);
}

void xt_run_unit(void) {
  t_crc();
  t_cm();
  t_lzp();
  t_rle();
  t_bwt();
}
