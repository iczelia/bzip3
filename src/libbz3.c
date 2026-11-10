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

#include "libbz3.h"
#include "sys/common.h"
#include "sys/crc.h"
#include "coder/bwt.h"
#include "coder/cm.h"
#include "coder/lzp.h"
#include "coder/rle.h"

#if defined(HAVE_THREADS) && defined(_WIN32)
#include <windows.h>
#elif defined(HAVE_THREADS)
#include <pthread.h>
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

/*  The entries of the suffix array for a block of x bytes.  */
#define BWT_BOUND(x) (bz3_bound(x) + 128)

/*  Blocks shorter than this are stored verbatim; entropy coding would not shrink them.  */
#define BZ3_STORED 64

/*  Bits of the model byte in a block header.  */
#define BZ3_MODEL_LZP 2
#define BZ3_MODEL_RLE 4

#define BZ3_FRAME_HEADER 13

struct bz3_state {
  u8 * swap;
  i32 block_size;
  i32 * sais, * lzp_lut;
  bz3_cm * cm;
  i8 last_error;
  bz3_crc crc;
};

BZIP3_API int8_t bz3_last_error(struct bz3_state * state) {
  return state->last_error;
}

BZIP3_API const char * bz3_version(void) { return PACKAGE_VERSION; }

BZIP3_API size_t bz3_bound(size_t input_size) {
  return input_size + input_size / 50 + 50;
}

BZIP3_API const char * bz3_strerror(struct bz3_state * state) {
  switch (state->last_error) {
  case BZ3_OK: return "No error";
  case BZ3_ERR_OUT_OF_BOUNDS: return "Data index out of bounds";
  case BZ3_ERR_BWT: return "Burrows-Wheeler transform failed";
  case BZ3_ERR_CRC: return "CRC32 check failed";
  case BZ3_ERR_MALFORMED_HEADER: return "Malformed header";
  case BZ3_ERR_TRUNCATED_DATA: return "Truncated data";
  case BZ3_ERR_DATA_TOO_BIG: return "Too much data";
  case BZ3_ERR_INIT: return "Initialization failed";
  case BZ3_ERR_DATA_SIZE_TOO_SMALL:
    return "Size of buffer `buffer_size` passed to the block decoder (bz3_decode_block) is too small. See "
           "function docs for details.";
  default: return "Unknown error";
  }
}

static int valid_block_size(i32 block_size) {
  return block_size >= BZ3_MIN_BLOCK && block_size <= BZ3_MAX_BLOCK;
}

BZIP3_API struct bz3_state * bz3_new(int32_t block_size) {
  struct bz3_state * state;
  if (!valid_block_size(block_size)) return NULL;

  state = malloc(sizeof *state);
  if (!state) return NULL;

  state->cm = malloc(sizeof *state->cm);
  state->swap = calloc(bz3_bound((sz) block_size), 1);
  state->sais = malloc(BWT_BOUND((sz) block_size) * sizeof (i32));
  state->lzp_lut = malloc(LZP_DICT * sizeof (i32));

  if (!state->cm || !state->swap || !state->sais || !state->lzp_lut) {
    bz3_free(state);
    return NULL;
  }

  state->block_size = block_size;
  state->last_error = BZ3_OK;
  bz3_crc_init(&state->crc);

  return state;
}

BZIP3_API void bz3_free(struct bz3_state * state) {
  if (!state) return;
  free(state->swap);
  free(state->sais);
  free(state->cm);
  free(state->lzp_lut);
  free(state);
}

static i32 fail(struct bz3_state * state, i8 error) {
  state->last_error = error;
  return -1;
}

BZIP3_API int32_t bz3_encode_block(struct bz3_state * state, uint8_t * buffer,
                                   int32_t data_size) {
  /*  `src' holds the output of the most recent stage, `dst' is free to be overwritten.  */
  u8 * src = buffer, * dst = state->swap, * tmp;
  sz capacity;
  u32 crc;
  i32 lzp_size, bwt_idx, header_size, coded_size;
  int model = 0;

  if (data_size < 0 || data_size > state->block_size)
    return fail(state, BZ3_ERR_DATA_TOO_BIG);

  /*  The caller guarantees room for this many bytes in `buffer'.  */
  capacity = bz3_bound((sz) data_size);
  crc = bz3_crc32c(&state->crc, 1, buffer, (sz) data_size);

  if (data_size < BZ3_STORED) {
    memmove(buffer + 8, buffer, (sz) data_size);
    wr32le(buffer, crc);
    wr32le(buffer + 4, (u32) -1);
    state->last_error = BZ3_OK;
    return data_size + 8;
  }

  lzp_size = bz3_lzp_encode(src, dst, data_size, state->lzp_lut);
  if (lzp_size > 0 && lzp_size < data_size) {
    tmp = src;  src = dst;  dst = tmp;
    data_size = lzp_size;
    model |= BZ3_MODEL_LZP;
  }

  /*  The suffix array construction makes use of the spare room in its work area.  */
  bwt_idx = bz3_bwt(src, dst, state->sais, data_size,
                    (i32) (BWT_BOUND((sz) state->block_size) - (sz) data_size));
  if (bwt_idx < 0) return fail(state, BZ3_ERR_BWT);

  /*  Header: CRC, BWT index, model, then the optional LZP size.  */
  header_size = model & BZ3_MODEL_LZP ? 13 : 9;

  /*  The BWT output lives in `dst', the entropy coder overwrites `src'.  */
  coded_size = bz3_cm_encode(state->cm, dst, data_size, src + header_size,
                             src + capacity);
  if (coded_size < 0) return fail(state, BZ3_ERR_DATA_TOO_BIG);

  wr32le(src, crc);
  wr32le(src + 4, (u32) bwt_idx);
  src[8] = (u8) model;
  if (model & BZ3_MODEL_LZP) wr32le(src + 9, (u32) lzp_size);

  if (src != buffer) memcpy(buffer, src, (sz) (header_size + coded_size));

  state->last_error = BZ3_OK;
  return header_size + coded_size;
}

/*  The sizes found in the header of a block that is not stored verbatim.  */
typedef struct {
  i32 header_size;         /*  The offset of the entropy coded data.  */
  i32 lzp_size, rle_size;  /*  -1 if the stage was not used.  */
} block_header;

/*  Parse the model byte and the sizes that follow it.  Returns 0 if `available' bytes are not enough.  */
static int read_block_header(const u8 * block, sz available, block_header * h) {
  const int model = block[8];
  const u8 * field = block + 9;
  h->header_size = 9;
  h->lzp_size = h->rle_size = -1;
  if (model & BZ3_MODEL_LZP) h->header_size += 4;
  if (model & BZ3_MODEL_RLE) h->header_size += 4;
  if (available < (sz) h->header_size) return 0;

  if (model & BZ3_MODEL_LZP) { h->lzp_size = (i32) rd32le(field);  field += 4; }
  if (model & BZ3_MODEL_RLE) h->rle_size = (i32) rd32le(field);
  return 1;
}

/*  Whether every intermediate form of a block fits a buffer: the original data, with LZP, with RLE, or with
    both applied.  */
static int sizes_fit(sz buffer_size, i32 lzp_size, i32 rle_size, i32 orig_size) {
  return (lzp_size < 0 || (sz) lzp_size <= buffer_size) &&
         (rle_size < 0 || (sz) rle_size <= buffer_size) &&
         (orig_size < 0 || (sz) orig_size <= buffer_size);
}

BZIP3_API int32_t bz3_decode_block(struct bz3_state * state, uint8_t * buffer,
                                   size_t buffer_size, int32_t compressed_size,
                                   int32_t orig_size) {
  const sz max_size = bz3_bound((sz) state->block_size);
  u8 * const swap = state->swap;
  /*  `src' holds the output of the most recent stage, `dst' is free to be overwritten.  */
  u8 * src = buffer, * dst = swap, * tmp;
  block_header h;
  u32 crc;
  i32 bwt_idx, bwt_size, size;
  int has_lzp, has_rle;

  /*  The initial header must be there, and the compressed data must fit the claimed buffer size.  */
  if (buffer_size < 9 || compressed_size < 0 ||
      buffer_size < (sz) compressed_size)
    return fail(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);

  if ((sz) compressed_size > max_size)
    return fail(state, BZ3_ERR_MALFORMED_HEADER);

  /*  The original size never exceeds the block size, whatever the block type.  */
  if (orig_size > state->block_size || orig_size < 0)
    return fail(state, BZ3_ERR_MALFORMED_HEADER);

  /*  Nothing past the compressed data is interpreted, even where the buffer extends further.  */
  if (compressed_size < 8) return fail(state, BZ3_ERR_MALFORMED_HEADER);

  crc = rd32le(buffer);
  bwt_idx = (i32) rd32le(buffer + 4);

  if (bwt_idx == -1) {
    /*  A block stored verbatim.  */
    if (compressed_size - 8 > BZ3_STORED || compressed_size - 8 != orig_size)
      return fail(state, BZ3_ERR_MALFORMED_HEADER);

    memmove(buffer, buffer + 8, (sz) orig_size);

    if (bz3_crc32c(&state->crc, 1, buffer, (sz) orig_size) != crc)
      return fail(state, BZ3_ERR_CRC);

    state->last_error = BZ3_OK;
    return orig_size;
  }

  /*  The model byte and the sizes that follow it are a part of the compressed data.  */
  if (compressed_size < 9 ||
      !read_block_header(buffer, (sz) compressed_size, &h))
    return fail(state, BZ3_ERR_MALFORMED_HEADER);

  has_lzp = (buffer[8] & BZ3_MODEL_LZP) != 0;
  has_rle = (buffer[8] & BZ3_MODEL_RLE) != 0;

  if ((has_lzp && (h.lzp_size < 0 || (sz) h.lzp_size > max_size)) ||
      (has_rle && (h.rle_size < 0 || (sz) h.rle_size > max_size)))
    return fail(state, BZ3_ERR_MALFORMED_HEADER);

  /*  The size that undoing the entropy coding and the BWT yields.  */
  bwt_size = has_lzp ? h.lzp_size : has_rle ? h.rle_size : orig_size;

  if (!sizes_fit(buffer_size, h.lzp_size, h.rle_size, orig_size))
    return fail(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);

  if (bwt_idx > bwt_size) return fail(state, BZ3_ERR_MALFORMED_HEADER);

  /*  Undo the entropy coding: buffer -> swap.  bwt_size <= max_size, the size of the swap buffer.  */
  bz3_cm_decode(state->cm, buffer + h.header_size, buffer + compressed_size,
                swap, bwt_size);

  /*  Undo the BWT: swap -> buffer.  bwt_size <= buffer_size was verified above.  */
  if (bwt_size <= 1) {
    /*  Degenerate transforms; any index was accepted for them historically.  */
    if (bwt_size == 1) buffer[0] = swap[0];
  } else if (bz3_unbwt(swap, buffer, state->sais, bwt_size, bwt_idx) < 0)
    return fail(state, BZ3_ERR_BWT);

  size = bwt_size;

  if (has_lzp) {
    /*  buffer -> swap, which holds max_size bytes.  */
    size = bz3_lzp_decode(src, dst, h.lzp_size, (i32) max_size, state->lzp_lut);
    if (size == -1) return fail(state, BZ3_ERR_CRC);
    /*  The headers passed the validation above, but the data itself may still decode into more than fits
        the buffer of the caller.  */
    if ((sz) size > buffer_size) return fail(state, BZ3_ERR_DATA_SIZE_TOO_SMALL);
    tmp = src;  src = dst;  dst = tmp;
  }

  if (has_rle) {
    /*  The output is capped at orig_size, which fits either buffer: it does not exceed the block size and
        was checked against buffer_size above.  */
    if (bz3_rle_decode(src, size, dst, orig_size)) return fail(state, BZ3_ERR_CRC);
    size = orig_size;
    tmp = src;  src = dst;  dst = tmp;
  }

  /*  Callers rely on orig_size bytes of `buffer' holding decoded data.  */
  if (size != orig_size) return fail(state, BZ3_ERR_MALFORMED_HEADER);

  if (src != buffer) memcpy(buffer, src, (sz) size);

  if (crc != bz3_crc32c(&state->crc, 1, buffer, (sz) size))
    return fail(state, BZ3_ERR_CRC);

  state->last_error = BZ3_OK;
  return size;
}

/*  A block to code, possibly on a thread of its own.  */
typedef struct {
  struct bz3_state * state;
  u8 * buffer;
  sz buffer_size;
  i32 size, orig_size;
  int decode, started;
#if defined(HAVE_THREADS) && defined(_WIN32)
  HANDLE thread;
#elif defined(HAVE_THREADS)
  pthread_t thread;
#endif
} job;

static void job_run(job * j) {
  if (j->decode)
    bz3_decode_block(j->state, j->buffer, j->buffer_size, j->size, j->orig_size);
  else j->size = bz3_encode_block(j->state, j->buffer, j->size);
}

#if defined(HAVE_THREADS) && defined(_WIN32)

static DWORD WINAPI job_thread(LPVOID arg) { job_run(arg);  return 0; }

static int job_start(job * j) {
  j->thread = CreateThread(NULL, 0, job_thread, j, 0, NULL);
  return j->thread != NULL;
}

static void job_join(job * j) {
  WaitForSingleObject(j->thread, INFINITE);
  CloseHandle(j->thread);
}

#elif defined(HAVE_THREADS)

static void * job_thread(void * arg) { job_run(arg);  return NULL; }

static int job_start(job * j) {
  return !pthread_create(&j->thread, NULL, job_thread, j);
}

static void job_join(job * j) { pthread_join(j->thread, NULL); }

#else

static int job_start(job * j) { (void) j;  return 0; }
static void job_join(job * j) { (void) j; }

#endif

/*  Run the jobs at once.  The last one runs on the calling thread, and so does every job for which no thread
    could be started.  */
static void jobs_run(job * jobs, i32 n) {
  i32 i;
  Fi(n - 1, jobs[i].started = job_start(&jobs[i]));
  jobs[n - 1].started = 0;
  for (i = n - 1; i >= 0; i--) if (!jobs[i].started) job_run(&jobs[i]);
  Fi(n, if (jobs[i].started) job_join(&jobs[i]));
}

BZIP3_API void bz3_encode_blocks(struct bz3_state * states[],
                                 uint8_t * buffers[], int32_t sizes[],
                                 int32_t n) {
  job * jobs;
  i32 i;
  if (n <= 0) return;

  jobs = calloc((sz) n, sizeof *jobs);
  if (!jobs) {
    /*  No memory for the bookkeeping: code the blocks one by one.  */
    Fi(n, sizes[i] = bz3_encode_block(states[i], buffers[i], sizes[i]));
    return;
  }

  Fi(n,
    jobs[i].state = states[i];
    jobs[i].buffer = buffers[i];
    jobs[i].size = sizes[i]);
  jobs_run(jobs, n);
  Fi(n, sizes[i] = jobs[i].size);

  free(jobs);
}

BZIP3_API void bz3_decode_blocks(struct bz3_state * states[],
                                 uint8_t * buffers[], size_t buffer_sizes[],
                                 int32_t sizes[], int32_t orig_sizes[],
                                 int32_t n) {
  job * jobs;
  i32 i;
  if (n <= 0) return;

  jobs = calloc((sz) n, sizeof *jobs);
  if (!jobs) {
    Fi(n, bz3_decode_block(states[i], buffers[i], buffer_sizes[i], sizes[i],
                           orig_sizes[i]));
    return;
  }

  Fi(n,
    jobs[i].state = states[i];
    jobs[i].buffer = buffers[i];
    jobs[i].buffer_size = buffer_sizes[i];
    jobs[i].size = sizes[i];
    jobs[i].orig_size = orig_sizes[i];
    jobs[i].decode = 1);
  jobs_run(jobs, n);

  free(jobs);
}

BZIP3_API int bz3_compress(uint32_t block_size, const uint8_t * in,
                           uint8_t * out, size_t in_size, size_t * out_size) {
  const sz buf_max = *out_size;
  struct bz3_state * state;
  u8 * buf;
  sz n_blocks, out_pos = BZ3_FRAME_HEADER, in_offset = 0, i;
  int result = BZ3_OK;

  if (block_size > in_size) block_size = (u32) bz3_bound(in_size);
  if (block_size <= BZ3_MIN_BLOCK) block_size = BZ3_MIN_BLOCK;
  if (block_size > BZ3_MAX_BLOCK) return BZ3_ERR_INIT;

  *out_size = 0;
  n_blocks = in_size / block_size + (in_size % block_size != 0);

  if (buf_max < BZ3_FRAME_HEADER || buf_max < bz3_bound(in_size) ||
      n_blocks > INT32_MAX)
    return BZ3_ERR_DATA_TOO_BIG;

  state = bz3_new((i32) block_size);
  buf = malloc(bz3_bound(block_size));
  if (!state || !buf) { result = BZ3_ERR_INIT;  goto done; }

  memcpy(out, "BZ3v1", 5);
  wr32le(out + 5, block_size);
  wr32le(out + 9, (u32) n_blocks);

  Fi(n_blocks,
    const i32 size = (i32) MIN(in_size - in_offset, (sz) block_size);
    i32 packed;
    memcpy(buf, in + in_offset, (sz) size);
    packed = bz3_encode_block(state, buf, size);
    if (packed < 0) { result = bz3_last_error(state);  goto done; }
    if (buf_max - out_pos < (sz) packed + 8) {
      result = BZ3_ERR_DATA_TOO_BIG;  goto done;
    }
    wr32le(out + out_pos, (u32) packed);
    wr32le(out + out_pos + 4, (u32) size);
    memcpy(out + out_pos + 8, buf, (sz) packed);
    out_pos += (sz) packed + 8;
    in_offset += (sz) size);

  *out_size = out_pos;

done:
  bz3_free(state);
  free(buf);
  return result;
}

BZIP3_API int bz3_decompress(const uint8_t * in, uint8_t * out, size_t in_size,
                             size_t * out_size) {
  const sz buf_max = *out_size;
  struct bz3_state * state;
  u8 * buf;
  sz buf_size, out_pos = 0;
  u32 n_blocks, i;
  i32 block_size;
  int result = BZ3_OK;

  if (in_size < BZ3_FRAME_HEADER || memcmp(in, "BZ3v1", 5))
    return BZ3_ERR_MALFORMED_HEADER;

  block_size = (i32) rd32le(in + 5);
  n_blocks = rd32le(in + 9);
  in_size -= BZ3_FRAME_HEADER;
  in += BZ3_FRAME_HEADER;

  *out_size = 0;

  if (!valid_block_size(block_size)) return BZ3_ERR_INIT;

  buf_size = bz3_bound((sz) block_size);
  state = bz3_new(block_size);
  buf = malloc(buf_size);
  if (!state || !buf) { result = BZ3_ERR_INIT;  goto done; }

  Fi(n_blocks,
    i32 size, orig_size;
    if (in_size < 8) { result = BZ3_ERR_MALFORMED_HEADER;  goto done; }
    size = (i32) rd32le(in);
    orig_size = (i32) rd32le(in + 4);
    if (size < 0 || (sz) size > buf_size) {
      result = BZ3_ERR_MALFORMED_HEADER;  goto done;
    }
    if (in_size - 8 < (sz) size) { result = BZ3_ERR_TRUNCATED_DATA;  goto done; }
    if (orig_size < 0 || orig_size > block_size) {
      result = BZ3_ERR_MALFORMED_HEADER;  goto done;
    }
    if (buf_max - out_pos < (sz) orig_size) {
      result = BZ3_ERR_DATA_TOO_BIG;  goto done;
    }
    memcpy(buf, in + 8, (sz) size);
    if (bz3_decode_block(state, buf, buf_size, size, orig_size) < 0) {
      result = bz3_last_error(state);  goto done;
    }
    memcpy(out + out_pos, buf, (sz) orig_size);
    out_pos += (sz) orig_size;
    in += (sz) size + 8;
    in_size -= (sz) size + 8);

done:
  *out_size = out_pos;
  bz3_free(state);
  free(buf);
  return result;
}

BZIP3_API size_t bz3_min_memory_needed(int32_t block_size) {
  if (!valid_block_size(block_size)) return 0;

  /*  This mirrors the allocations made by bz3_new.  */
  return sizeof (struct bz3_state) + sizeof (bz3_cm) +
         bz3_bound((sz) block_size) +
         BWT_BOUND((sz) block_size) * sizeof (i32) + LZP_DICT * sizeof (i32);
}

BZIP3_API int bz3_orig_size_sufficient_for_decode(const uint8_t * block,
                                                  size_t block_size,
                                                  int32_t orig_size) {
  block_header h;

  /*  The initial header takes 9 bytes: CRC, BWT index and the model.  */
  if (block_size < 9) return -1;

  /*  Blocks stored verbatim never need more than the original size.  */
  if ((i32) rd32le(block + 4) == -1) return 1;

  if (!read_block_header(block, block_size, &h)) return -1;

  return sizes_fit(orig_size < 0 ? 0 : (sz) orig_size, h.lzp_size, h.rle_size,
                   orig_size);
}
