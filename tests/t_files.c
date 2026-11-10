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

/*  Decode a file written by the command line tool: a 9 byte header, then blocks behind their two sizes.
    Returns the decoded data, or NULL.  `model' receives the model byte of the first block, -1 if stored.  */
static u8 * unpack(const u8 * in, sz n, sz * out_len, int * model) {
  struct bz3_state * state;
  u8 * out = NULL, * buf;
  sz cap, len = 0, at = 9;
  i32 block_size;
  *model = -2;
  if (n < 9 || memcmp(in, "BZ3v1", 5)) return NULL;
  block_size = (i32) rd32le(in + 5);
  state = bz3_new(block_size);
  if (!state) return NULL;
  cap = bz3_bound((sz) block_size);
  buf = xt_alloc(cap);
  while (at < n) {
    i32 size, orig;
    if (n - at < 8) goto bad;
    size = (i32) rd32le(in + at);  orig = (i32) rd32le(in + at + 4);
    at += 8;
    if (size < 8 || (sz) size > n - at || (sz) size > cap) goto bad;
    if (*model == -2)
      *model = (i32) rd32le(in + at + 4) == -1 ? -1 : size > 8 ? in[at + 8] : -2;
    memcpy(buf, in + at, (sz) size);
    if (bz3_decode_block(state, buf, cap, size, orig) != orig) goto bad;
    out = realloc(out, len + (sz) orig + 1);
    if (!out) exit(99);
    memcpy(out + len, buf, (sz) orig);
    len += (sz) orig;
    at += (sz) size;
  }
  if (!out) out = xt_alloc(1);
  *out_len = len;
  bz3_free(state);  free(buf);
  return out;
bad:
  bz3_free(state);  free(buf);  free(out);
  return NULL;
}

/*  Files compressed by bzip3 1.5.5 decode into their originals.  They hold every kind of block, those of
    the RLE transform included, which nothing else produces any more.  */
static void t_golden(void) {
  const xt_golden * g;
  char name[256];
  xt_section_begin("reference files");
  if (!xt_data) { CHECK(0, "BZ3_TEST_DATA is not set");  return; }
  for (g = XT_GOLDEN; g->name; g++) {
    sz n, pn, un = 0, cap, fn;
    u8 * orig = xt_slurp(xt_fixture(g->name), &n), * packed, * back, * frame;
    int model, err;
    snprintf(name, sizeof name, "%s.bz3", g->name);
    packed = xt_slurp(xt_fixture(name), &pn);
    CHECK(orig && packed, "%s and its compressed form are there", g->name);
    if (!orig || !packed) { free(orig);  free(packed);  continue; }

    back = unpack(packed, pn, &un, &model);
    CHECK(back && un == n && !memcmp(orig, back, n), "%s decodes into the original", name);
    CHECK(model == g->model, "%s holds a block of model %d, not %d", name, g->model, model);
    free(back);

    /*  What this version makes of the original decodes too.  */
    cap = fn = bz3_bound(n) + 64;
    frame = xt_alloc(cap);
    back = xt_alloc(n + 1);
    err = bz3_compress(MiB(1), orig, frame, n, &fn);
    un = n;
    CHECK(err == BZ3_OK && bz3_decompress(frame, back, fn, &un) == BZ3_OK &&
          un == n && !memcmp(orig, back, n), "%s round trip", g->name);
    free(frame);  free(back);  free(orig);  free(packed);
  }
}

void xt_run_files(void) {
  t_golden();
}
