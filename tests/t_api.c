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

static xt_rng rng;

static void t_basics(void) {
  struct bz3_state * state;
  xt_section_begin("api basics");
  CHECK(bz3_version() != NULL && *bz3_version(), "the library has a version");
  CHECK(bz3_bound(0) >= 13 && bz3_bound(1000) > 1000, "bz3_bound leaves room");
  CHECK(!bz3_new(0) && !bz3_new(KiB(65) - 1) && !bz3_new(-1) && !bz3_new(MiB(512)),
        "block sizes out of range are refused");
  CHECK(!bz3_min_memory_needed(KiB(64)) && bz3_min_memory_needed(KiB(65)) > (sz) KiB(65) * 5,
        "bz3_min_memory_needed counts the buffers");
  bz3_free(NULL);

  state = bz3_new(KiB(65));
  CHECK(state != NULL, "a state for the smallest block");
  if (state) {
    CHECK(bz3_last_error(state) == BZ3_OK && bz3_strerror(state) != NULL,
          "a new state has no error");
    bz3_free(state);
  }
}

/*  Compress and decompress `size' bytes with the frame API.  */
static void t_frame(sz size, int kind, u32 block_size) {
  sz cap = bz3_bound(size), out_size = cap, back_size = size;
  u8 * in = xt_alloc(size + 1), * out = xt_alloc(cap), * back = xt_alloc(size + 1);
  int err;
  xt_fill(&rng, in, size, kind);

  err = bz3_compress(block_size, in, out, size, &out_size);
  CHECK(err == BZ3_OK && out_size <= cap, "%lu bytes of kind %d compress",
        (unsigned long) size, kind);
  if (err == BZ3_OK) {
    err = bz3_decompress(out, back, out_size, &back_size);
    CHECK(err == BZ3_OK && back_size == size && !memcmp(in, back, size),
          "%lu bytes of kind %d round trip", (unsigned long) size, kind);

    /*  An output buffer that is too small is reported, not overrun.  */
    if (size > 0) {
      back_size = size - 1;
      CHECK(bz3_decompress(out, back, out_size, &back_size) == BZ3_ERR_DATA_TOO_BIG,
            "%lu bytes of kind %d: short output buffer", (unsigned long) size, kind);
    }
  }

  /*  So is one for the compressed data.  */
  out_size = cap - 1;
  CHECK(bz3_compress(block_size, in, out, size, &out_size) == BZ3_ERR_DATA_TOO_BIG,
        "%lu bytes of kind %d: short compression buffer", (unsigned long) size, kind);

  free(in);  free(out);  free(back);
}

static void t_frames(void) {
  static const sz sizes[] = {
    0, 1, 63, 64, 65, 100, 4096, KiB(65) - 1, KiB(65), KiB(65) + 1, 2 * KiB(65), 300000
  };
  sz i;
  int kind;
  xt_section_begin("api frames");
  Fi(ARRAY_LEN(sizes), for (kind = 0; kind < 5; kind++) t_frame(sizes[i], kind, KiB(65)));
  t_frame(300000, 4, MiB(1));
  t_frame(300000, 1, 1);
}

/*  Damaged frames are rejected, or at the very least decode into the original data.  */
static void t_damage(void) {
  enum { SIZE = 3000 };
  u8 in[SIZE], back[SIZE];
  sz cap = bz3_bound(SIZE), out_size = cap, back_size, i, len;
  u8 * out = xt_alloc(cap), * copy = xt_alloc(cap);
  int bit, err;
  xt_section_begin("api damage");
  xt_fill(&rng, in, SIZE, 4);
  CHECK(bz3_compress(KiB(65), in, out, SIZE, &out_size) == BZ3_OK, "compress");

  /*  The 13 byte frame header carries no checksum: another block size or block count there can still
      describe a well-formed frame.  Everything after it is covered.  */
  for (i = 13; i < out_size; i++)
    for (bit = 0; bit < 8; bit += 3) {
      memcpy(copy, out, out_size);
      copy[i] ^= (u8) (1 << bit);
      back_size = SIZE;
      err = bz3_decompress(copy, back, out_size, &back_size);
      CHECK(err != BZ3_OK || (back_size == SIZE && !memcmp(in, back, SIZE)),
            "bit %d of byte %lu flipped", bit, (unsigned long) i);
    }

  /*  Every proper prefix lacks data.  */
  for (len = 0; len < out_size; len++) {
    memcpy(copy, out, len);
    back_size = SIZE;
    CHECK(bz3_decompress(copy, back, len, &back_size) != BZ3_OK,
          "a prefix of %lu bytes", (unsigned long) len);
  }

  free(out);  free(copy);
}

static void t_blocks(void) {
  enum { BLOCK = KiB(65) };
  static const int sizes[] = { 0, 1, 63, 64, 65, 71, 72, 73, 1000, 30000, BLOCK - 1, BLOCK };
  struct bz3_state * state = bz3_new(BLOCK);
  sz cap = bz3_bound(BLOCK), i;
  u8 * in = xt_alloc(BLOCK), * buf = xt_alloc(cap);
  int kind;
  xt_section_begin("api blocks");
  CHECK(state != NULL, "a state");
  if (!state) return;

  Fi(ARRAY_LEN(sizes),
    for (kind = 0; kind < 5; kind++) {
      int size = sizes[i];
      i32 packed;
      xt_fill(&rng, in, (sz) size, kind);
      memcpy(buf, in, (sz) size);
      packed = bz3_encode_block(state, buf, size);
      CHECK(packed > 0 && (sz) packed <= bz3_bound((sz) size) &&
            bz3_last_error(state) == BZ3_OK, "%d bytes of kind %d encode", size, kind);
      if (packed <= 0) continue;

      /*  The encoder never uses the RLE transform.  */
      CHECK(packed < 9 || (i32) rd32le(buf + 4) == -1 || !(buf[8] & ~2),
            "%d bytes of kind %d: the model is plain or LZP", size, kind);

      /*  8 bytes stand for a stored empty block, which has no model byte to inspect.  */
      if (packed >= 9)
        CHECK(bz3_orig_size_sufficient_for_decode(buf, (sz) packed, size) >= 0,
              "%d bytes of kind %d: the header parses", size, kind);

      /*  Too small a buffer is refused before anything is written past it.  */
      if (packed > 9)
        CHECK(bz3_decode_block(state, buf, (sz) packed - 1, packed, size) == -1,
              "%d bytes of kind %d: short buffer", size, kind);
      CHECK(bz3_decode_block(state, buf, cap, packed, size + 1) == -1 &&
            bz3_last_error(state) != BZ3_OK, "%d bytes of kind %d: wrong size", size, kind);

      /*  The failed attempts above may have clobbered the buffer.  */
      memcpy(buf, in, (sz) size);
      packed = bz3_encode_block(state, buf, size);
      CHECK(bz3_decode_block(state, buf, cap, packed, size) == size &&
            bz3_last_error(state) == BZ3_OK && !memcmp(buf, in, (sz) size),
            "%d bytes of kind %d round trip", size, kind);
    });

  CHECK(bz3_encode_block(state, buf, BLOCK + 1) == -1 &&
        bz3_last_error(state) == BZ3_ERR_DATA_TOO_BIG, "a block too large is refused");
  CHECK(bz3_encode_block(state, buf, -1) == -1, "a negative size is refused");
  CHECK(bz3_decode_block(state, buf, cap, -1, 10) == -1, "a negative compressed size");
  CHECK(bz3_decode_block(state, buf, cap, 100, -1) == -1, "a negative original size");
  CHECK(bz3_decode_block(state, buf, 8, 8, 0) == -1, "a buffer below the header size");

  bz3_free(state);
  free(in);  free(buf);
}

/*  The blocks come out the same whether they are coded at once or one by one.  */
static void t_parallel(void) {
  enum { N = 5, BLOCK = KiB(65) };
  struct bz3_state * states[N], * one = bz3_new(BLOCK);
  u8 * bufs[N], * orig[N], * ref = xt_alloc(bz3_bound(BLOCK));
  sz caps[N];
  i32 sizes[N], orig_sizes[N], packed;
  int i;
  xt_section_begin("api parallel");
  Fi(N,
    states[i] = bz3_new(BLOCK);
    caps[i] = bz3_bound(BLOCK);
    bufs[i] = xt_alloc(caps[i]);
    orig[i] = xt_alloc(BLOCK);
    orig_sizes[i] = sizes[i] = BLOCK - i * 9000;
    xt_fill(&rng, orig[i], (sz) sizes[i], i);
    memcpy(bufs[i], orig[i], (sz) sizes[i]));
  bz3_encode_blocks(states, bufs, sizes, N);
  Fi(N,
    CHECK(sizes[i] > 0 && bz3_last_error(states[i]) == BZ3_OK, "block %d encodes", i);
    memcpy(ref, orig[i], (sz) orig_sizes[i]);
    packed = bz3_encode_block(one, ref, orig_sizes[i]);
    CHECK(packed == sizes[i] && !memcmp(ref, bufs[i], (sz) packed),
          "block %d matches the block coded alone", i));
  bz3_decode_blocks(states, bufs, caps, sizes, orig_sizes, N);
  Fi(N,
    CHECK(bz3_last_error(states[i]) == BZ3_OK &&
          !memcmp(bufs[i], orig[i], (sz) orig_sizes[i]), "block %d round trip", i);
    bz3_free(states[i]);
    free(bufs[i]);  free(orig[i]));
  bz3_free(one);
  free(ref);
}

void xt_run_api(void) {
  xt_seed(&rng, 5);
  t_basics();
  t_frames();
  t_damage();
  t_blocks();
  t_parallel();
}
