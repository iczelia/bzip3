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

#include "sys/common.h"
#include "sys/map.h"
#include "sys/mem.h"
#include "sys/meta.h"
#include "sys/thread.h"
#include "sys/win32.h"
#include "sys/yarg.h"
#include "libbz3.h"

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <sys/stat.h>
/*  Compiled by hand, without configure: the headers are where each kind of system keeps them.  */
#if !defined(HAVE_CONFIG_H) && defined(_WIN32)
#define HAVE_IO_H 1
#elif !defined(HAVE_CONFIG_H)
#define HAVE_UNISTD_H 1
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_IO_H
#include <io.h>
#endif

#ifndef O_BINARY
#define O_BINARY 0
#endif
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

enum { M_DECODE, M_ENCODE, M_TEST, M_RECOVER };

#define BZ3_MAX_JOBS 64
#define BZ3_DEFAULT_BLOCK MiB(16)
#define BZ3_OUT_BUFFER MiB(1)

/*  The regular file being written, if any.  It is removed when the operation fails, so that no truncated
    output is left behind.  */
static const char * partial;
static FILE * partial_file;

#ifdef _WIN32
#define sys_fopen  bz3_win_fopen
#define sys_remove bz3_win_remove
#define sys_is_dir bz3_win_is_dir
#else
#define sys_fopen  fopen
#define sys_remove remove
static int sys_is_dir(const char * path) {
  struct stat st;
  return !stat(path, &st) && S_ISDIR(st.st_mode);
}
#endif

static void vmsg(const char * fmt, va_list ap) BZ3_PRINTF(1, 0);

static void vmsg(const char * fmt, va_list ap) {
#ifdef _WIN32
  bz3_win_vprintf(stderr, fmt, ap);
#else
  vfprintf(stderr, fmt, ap);
#endif
}

static void msg(const char * fmt, ...) BZ3_PRINTF(1, 2);

static void msg(const char * fmt, ...) {
  va_list ap;
  va_start(ap, fmt);  vmsg(fmt, ap);  va_end(ap);
}

static void partial_remove(void) {
  if (partial_file) { fclose(partial_file);  partial_file = NULL; }
  if (partial) { sys_remove(partial);  partial = NULL; }
}

static NORETURN void fatal(const char * fmt, ...) BZ3_PRINTF(1, 2);

static void fatal(const char * fmt, ...) {
  va_list ap;
  msg("bzip3: ");
  va_start(ap, fmt);  vmsg(fmt, ap);  va_end(ap);
  msg("\n");
  partial_remove();
  exit(1);
}

static void banner(FILE * to) {
  fprintf(to,
    "bzip3 " PACKAGE_VERSION " -- a better and stronger spiritual successor to bzip2\n"
    "Written by Kamila Szewczyk <k@iczelia.net>.\n"
    "License GNU Lesser GPL version 3 <https://www.gnu.org/licenses/lgpl-3.0.html>.\n");
}

static void usage(FILE * to) {
  banner(to);
  fprintf(to,
    "\n"
    "usage: bzip3 [-e|-d|-t|-r] [-b N] [-j N] [-cfkv] [FILE [OUTPUT]]\n"
    "       bzip3 -B [-e|-d|-t|-r] [-b N] [-j N] [-cfkv] FILE...\n"
    "       bzip3 -h | -V\n"
    "\n"
    "  -e, -z, --encode  compress (the default)\n"
    "  -d, --decode      decompress\n"
    "  -t, --test        check the integrity of compressed data\n"
    "  -r, --recover     decompress, writing out blocks that fail their checks\n"
    "  -b, --block=N     block size in MiB, 1 to %d; %d is the default\n"
    "  -j, --jobs=N      code N blocks at once, up to %d; 1 is the default\n"
    "  -B, --batch       treat every file as an input\n"
    "  -c, --stdout      write to the standard output\n"
    "  -f, --force       overwrite existing output files\n"
    "      --rm          remove the input files once they were processed\n"
    "  -k, --keep        keep the input files (the default)\n"
    "  -v, --verbose     print the sizes and the compression ratio\n"
    "  -h, --help        show help        -V, --version   show version\n"
    "\n"
    "Without FILE, or when it is -, the standard input is read.\n"
    "Report bugs to <https://github.com/iczelia/bzip3>.\n",
    BZ3_MAX_BLOCK / MiB(1), BZ3_DEFAULT_BLOCK / MiB(1), BZ3_MAX_JOBS);
}

static void version(void) {
  banner(stdout);
#ifdef BZ3_HOST_TRIPLE
  printf("Build host %s.\n", BZ3_HOST_TRIPLE);
#endif
  printf("Copyright (C) 2022-2026 Kamila Szewczyk.\n"
         "This is free software: you are free to change and redistribute it.\n"
         "There is NO WARRANTY, to the extent permitted by law.\n");
  printf("Library: libbzip3 %s\n", bz3_version());
}

static void * xcalloc(sz n, sz size) {
  void * p = calloc(n ? n : 1, size ? size : 1);
  if (!p) fatal("out of memory");
  return p;
}

static void close_output(FILE * f, const bz3_meta * meta, FILE * in, const char * in_path,
                         const char * out_path) {
  const int preserve = partial && meta->known && meta->regular;
  int lost = 0;
  if (!f) return;
  if (fflush(f)) fatal("write error: %s", strerror(errno));
#ifndef _WIN32
  if (preserve) lost = bz3_meta_put(meta, in, f, in_path, out_path);
#endif
#ifdef HAVE_FSYNC
  /*  Pipes, terminals and the like can not be synchronised; only failures to store the data count.  */
  while (fsync(fileno(f)) == -1) {
    if (errno == EINTR) continue;
    if (errno == EIO || errno == ENOSPC
#ifdef EDQUOT
        || errno == EDQUOT
#endif
       ) fatal("write error: %s", strerror(errno));
    break;
  }
#endif
  if (f != stdout) {
    partial_file = NULL;
    if (fclose(f)) fatal("write error: %s", strerror(errno));
  }
#ifdef _WIN32
  if (preserve) lost = bz3_meta_put(meta, in, NULL, in_path, out_path);
#endif
  partial = NULL;

  if (lost & BZ3_META_OWNER) msg("bzip3: %s: cannot preserve the owner\n", out_path);
  if (lost & BZ3_META_MODE) msg("bzip3: %s: cannot preserve the permissions\n", out_path);
  if (lost & BZ3_META_ACL) msg("bzip3: %s: cannot preserve the access control list\n", out_path);
  if (lost & BZ3_META_TIMES) msg("bzip3: %s: cannot preserve the times\n", out_path);
  if (lost & BZ3_META_ATTRS) msg("bzip3: %s: cannot preserve the file attributes\n", out_path);
}

typedef struct {
  int n;
  u8 ** buffers;
  i32 * sizes, * orig_sizes;
  int8_t * errors;
  const char ** reasons;
} batch;

typedef struct {
  int set, fatal;
  char text[256];
} failure;

typedef struct {
  bz3_in * in;
  FILE * out;
  int mode, enc, jobs, eof, stop;
  i32 block_size;
  struct bz3_state ** states;
  sz * caps;
  uint64_t n_in, n_out;
  failure read_failure, write_failure;
#ifdef HAVE_THREADS
  int threaded;
  bz3_mutex lock;
  bz3_queue * to_read, * to_code, * to_write;
#endif
} pipeline;

#define BZ3_MAX_BATCHES 4

static void fail(failure * f, int is_fatal, const char * fmt, ...) BZ3_PRINTF(3, 4);

static void fail(failure * f, int is_fatal, const char * fmt, ...) {
  va_list ap;
  if (f->set) return;
  f->set = 1;  f->fatal = is_fatal;
  va_start(ap, fmt);  vsnprintf(f->text, sizeof f->text, fmt, ap);  va_end(ap);
}

static void pipe_stop(pipeline * p) {
#ifdef HAVE_THREADS
  if (p->threaded) { bz3_mutex_lock(&p->lock);  p->stop = 1;  bz3_mutex_unlock(&p->lock);  return; }
#endif
  p->stop = 1;
}

static int pipe_stopped(pipeline * p) {
#ifdef HAVE_THREADS
  if (p->threaded) {
    int stop;
    bz3_mutex_lock(&p->lock);  stop = p->stop;  bz3_mutex_unlock(&p->lock);
    return stop;
  }
#endif
  return p->stop;
}

static long long stage_fill(pipeline * p, void * data, sz n) {
  const sz got = bz3_in_read(p->in, data, n);
  if (p->in->failed) {
    fail(&p->read_failure, 1, "read error: %s", strerror(errno));
    return -1;
  }
  return (long long) got;
}

static int stage_read(pipeline * p, batch * b) {
  u8 header[8];
  long long got;
  int n;
  b->n = 0;
  for (n = 0; n < p->jobs && !p->eof; n++) {
    if (p->enc) {
      got = stage_fill(p, b->buffers[n], (sz) p->block_size);
      if (got < p->block_size) p->eof = 1;
      if (got <= 0) break;
      b->sizes[n] = b->orig_sizes[n] = (i32) got;
      p->n_in += (uint64_t) got;
    } else {
      got = stage_fill(p, header, 8);
      if (got != 8) {
        if (got > 0) fail(&p->read_failure, 1, "the input is truncated");
        p->eof = 1;
        break;
      }
      b->sizes[n] = (i32) rd32le(header);
      b->orig_sizes[n] = (i32) rd32le(header + 4);
      if (b->sizes[n] < 0 || (sz) b->sizes[n] > p->caps[n] ||
          b->orig_sizes[n] < 0 || b->orig_sizes[n] > p->block_size) {
        fail(&p->read_failure, 0, "failed to decode a block: inconsistent headers");
        p->eof = 1;
        break;
      }
      got = stage_fill(p, b->buffers[n], (sz) b->sizes[n]);
      if (got != b->sizes[n]) {
        if (got >= 0) fail(&p->read_failure, 1, "the input is truncated");
        p->eof = 1;
        break;
      }
      p->n_in += 8 + (uint64_t) got;
    }
  }
  b->n = n;
  return n > 0;
}

static void stage_code(pipeline * p, batch * b) {
  int i;
  if (pipe_stopped(p)) return;
  if (b->n > 1) {
    if (p->enc) bz3_encode_blocks(p->states, b->buffers, b->sizes, b->n);
    else bz3_decode_blocks(p->states, b->buffers, p->caps, b->sizes, b->orig_sizes, b->n);
  } else if (b->n == 1) {
    if (p->enc) b->sizes[0] = bz3_encode_block(p->states[0], b->buffers[0], b->sizes[0]);
    else bz3_decode_block(p->states[0], b->buffers[0], p->caps[0], b->sizes[0], b->orig_sizes[0]);
  }
  Fi(b->n,
    b->errors[i] = bz3_last_error(p->states[i]);
    b->reasons[i] = bz3_strerror(p->states[i]));
}

static int stage_put(pipeline * p, const void * data, sz n) {
  if (!n || fwrite(data, 1, n, p->out) == n) return 1;
  fail(&p->write_failure, 1, "write error: %s", strerror(errno));
  return 0;
}

static void stage_write(pipeline * p, batch * b) {
  u8 header[8];
  int i;
  if (pipe_stopped(p)) return;
  Fi(b->n,
    if (b->errors[i] != BZ3_OK) {
      msg("bzip3: %s a block: %s\n",
          p->enc ? "failed to encode" : p->mode == M_RECOVER
            ? "writing an invalid" : "failed to decode",
          b->reasons[i]);
      if (p->mode != M_RECOVER) {
        fail(&p->write_failure, 0, "%s", "");
        break;
      }
    }
    if (p->enc) {
      wr32le(header, (u32) b->sizes[i]);
      wr32le(header + 4, (u32) b->orig_sizes[i]);
      if (!stage_put(p, header, 8) || !stage_put(p, b->buffers[i], (sz) b->sizes[i])) break;
      p->n_out += 8 + (uint64_t) b->sizes[i];
    } else {
      if (p->mode != M_TEST && !stage_put(p, b->buffers[i], (sz) b->orig_sizes[i])) break;
      p->n_out += (uint64_t) b->orig_sizes[i];
    });
  if (p->write_failure.set) pipe_stop(p);
}

#ifdef HAVE_THREADS

static void code_thread(void * arg) {
  pipeline * p = arg;
  batch * b;
  while ((b = bz3_queue_get(p->to_code)) != NULL) {
    stage_code(p, b);
    bz3_queue_put(p->to_write, b);
  }
  bz3_queue_put(p->to_write, NULL);
}

static void write_thread(void * arg) {
  pipeline * p = arg;
  batch * b;
  while ((b = bz3_queue_get(p->to_write)) != NULL) {
    stage_write(p, b);
    bz3_queue_put(p->to_read, b);
  }
}

static int pipe_threads(pipeline * p, batch * batches, int n) {
  bz3_thread coder, writer;
  batch * b;
  int i, ok = 0;
  p->to_read = bz3_queue_new(n + 1);
  p->to_code = bz3_queue_new(n + 1);
  p->to_write = bz3_queue_new(n + 1);
  if (p->to_read && p->to_code && p->to_write) {
    bz3_mutex_init(&p->lock);
    p->threaded = 1;
    if (bz3_thread_start(&coder, code_thread, p)) {
      if (bz3_thread_start(&writer, write_thread, p)) {
        ok = 1;
        Fi(n, bz3_queue_put(p->to_read, &batches[i]));
        for (;;) {
          b = bz3_queue_get(p->to_read);
          if (pipe_stopped(p) || !stage_read(p, b)) break;
          bz3_queue_put(p->to_code, b);
        }
      }
      bz3_queue_put(p->to_code, NULL);
      bz3_thread_join(&coder);
      if (ok) bz3_thread_join(&writer);
    }
    p->threaded = 0;
    bz3_mutex_free(&p->lock);
  }
  bz3_queue_free(p->to_read);  bz3_queue_free(p->to_code);  bz3_queue_free(p->to_write);
  return ok;
}

#endif

static void pipe_run(pipeline * p, batch * batches, int n) {
#ifdef HAVE_THREADS
  if (n > 1 && pipe_threads(p, batches, n)) return;
#else
  (void) n;
#endif
  while (!pipe_stopped(p) && stage_read(p, &batches[0])) {
    stage_code(p, &batches[0]);
    stage_write(p, &batches[0]);
  }
}

static void plan(const char * name, int enc, i32 block_size, int * jobs, int * batches) {
  const uint64_t buffer = bz3_bound((sz) block_size);
  const uint64_t job = bz3_min_memory_needed(block_size) + buffer;
  const uint64_t total = bz3_mem_total(), avail = bz3_mem_avail();
  uint64_t need, spare;
  int fit = *jobs;

  if (job > SIZE_MAX)
    fatal("%s: %s takes %" PRIu64 " MiB of memory, more than a program can have on this system",
          name, enc ? "compression" : "decompression", job / MiB(1));
  if (total && job > total)
    msg("bzip3: warning: %s: %s takes %" PRIu64 " MiB of memory, and only %" PRIu64 " MiB can be had\n",
        name, enc ? "compression" : "decompression", job / MiB(1), total / MiB(1));

  while (fit > 1 && ((total && (uint64_t) fit * job > total) || (uint64_t) fit * job > SIZE_MAX)) fit--;
  if (fit != *jobs)
    msg("bzip3: warning: %s: running %d %s and not %d, which would take %" PRIu64 " MiB of memory\n",
        name, fit, fit == 1 ? "job" : "jobs", *jobs, (uint64_t) *jobs * job / MiB(1));
  *jobs = fit;

  need = (uint64_t) fit * job;
  if (avail) spare = avail > need ? (avail - need) / 2 : 0;
  else if (total) spare = total > need ? (total - need) / 4 : 0;
  else spare = 0;
  if (sizeof (sz) < 8) spare = MIN(spare, (uint64_t) SIZE_MAX / 4);
  *batches = 1 + (int) MIN((uint64_t) (BZ3_MAX_BATCHES - 1), spare / ((uint64_t) fit * buffer));
#ifndef HAVE_THREADS
  *batches = 1;
#endif
}

static int process(FILE * fin, FILE * out, int mode, i32 block_size, int jobs,
                   int verbose, const char * name) {
  const int enc = mode == M_ENCODE;
  const char * shown = name ? name : "(stdin)";
  u8 header[9];
  bz3_in in;
  pipeline p;
  batch * batches;
  int n_batches, i, j;

  if (isatty(fileno(enc ? out : fin))) {
    msg("bzip3: compressed data is not %s a terminal\n",
            enc ? "written to" : "read from");
    return 1;
  }

  memset(&p, 0, sizeof p);
  bz3_in_open(&in, fin);
  if (in.failed) fatal("read error: %s", strerror(errno));
  p.in = &in;  p.out = out;  p.mode = mode;  p.enc = enc;

  if (enc) {
    memcpy(header, "BZ3v1", 5);
    wr32le(header + 5, (u32) block_size);
    if (fwrite(header, 1, 9, out) != 9) fatal("write error: %s", strerror(errno));
    p.n_out += 9;
  } else {
    const sz got = bz3_in_read(&in, header, 9);
    if (in.failed) fatal("read error: %s", strerror(errno));
    if (got < 5 || memcmp(header, "BZ3v1", 5)) {
      msg("bzip3: %s: not a bzip3 file\n", shown);
      bz3_in_close(&in);
      return 1;
    }
    if (got != 9) fatal("the input is truncated");
    block_size = (i32) rd32le(header + 5);
    if (block_size < BZ3_MIN_BLOCK || block_size > BZ3_MAX_BLOCK) {
      msg("bzip3: %s: invalid block size in the header\n", shown);
      if (mode != M_RECOVER) { bz3_in_close(&in);  return 1; }
      msg("bzip3: recovering with the largest block size\n");
      block_size = BZ3_MAX_BLOCK;
    }
    p.n_in += 9;
  }

  if (jobs < 1) jobs = 1;
#ifndef HAVE_THREADS
  jobs = 1;
#endif
  plan(shown, enc, block_size, &jobs, &n_batches);
  p.jobs = jobs;  p.block_size = block_size;

  p.states = xcalloc((sz) jobs, sizeof *p.states);
  p.caps = xcalloc((sz) jobs, sizeof *p.caps);
  Fi(jobs,
    p.states[i] = bz3_new(block_size);
    if (!p.states[i]) fatal("out of memory");
    p.caps[i] = bz3_bound((sz) block_size));

  batches = xcalloc((sz) n_batches, sizeof *batches);
  Fj(n_batches,
    batch * b = &batches[j];
    b->buffers = xcalloc((sz) jobs, sizeof *b->buffers);
    b->sizes = xcalloc((sz) jobs, sizeof *b->sizes);
    b->orig_sizes = xcalloc((sz) jobs, sizeof *b->orig_sizes);
    b->errors = xcalloc((sz) jobs, sizeof *b->errors);
    b->reasons = xcalloc((sz) jobs, sizeof *b->reasons);
    Fi(jobs,
      /*  Zeroed, so that the recovery mode never writes out uninitialised memory.  */
      b->buffers[i] = calloc(p.caps[i], 1);
      if (!b->buffers[i] && !j) fatal("out of memory");
      if (!b->buffers[i]) break);
    if (i < jobs) {
      Fi(jobs, free(b->buffers[i]));
      free(b->buffers);  free(b->sizes);  free(b->orig_sizes);  free(b->errors);  free(b->reasons);
      n_batches = j;
      break;
    });

  pipe_run(&p, batches, n_batches);

  Fj(n_batches,
    batch * b = &batches[j];
    Fi(jobs, free(b->buffers[i]));
    free(b->buffers);  free(b->sizes);  free(b->orig_sizes);  free(b->errors);  free(b->reasons));
  free(batches);
  Fi(jobs, bz3_free(p.states[i]));
  free(p.states);  free(p.caps);
  bz3_in_close(&in);

  { const failure * f = p.write_failure.set ? &p.write_failure : &p.read_failure;
    if (f->set && f->fatal) fatal("%s", f->text);
    if (f->set) {
      if (*f->text) msg("bzip3: %s\n", f->text);
      return 1;
    } }

  if (out && fflush(out)) fatal("write error: %s", strerror(errno));

  if (verbose) {
    /*  The size of the compressed data relative to the size of the original.  */
    const double packed = (double) (enc ? p.n_out : p.n_in);
    const double unpacked = (double) (enc ? p.n_in : p.n_out);
    if (name) msg(" %s:", name);
    msg("\t%s%" PRIu64 " -> %" PRIu64 " bytes, %.2f%%, %.2f bpb\n",
            mode == M_ENCODE || mode == M_DECODE ? "" : "OK, ", p.n_in, p.n_out,
            packed * 100.0 / unpacked, packed * 8.0 / unpacked);
  }

  return 0;
}

/*  Parse a decimal number of at most `max'.  Returns 0 on failure.  */
static int parse_number(const char * s, long max, long * value) {
  char * end;
  unsigned long n;
  if (!s || *s < '0' || *s > '9') return 0;
  errno = 0;
  n = strtoul(s, &end, 10);
  if (errno || *end || n > (unsigned long) max) return 0;
  *value = (long) n;
  return 1;
}

static FILE * open_input(const char * path) {
  FILE * f;
  if (!path) return stdin;
  if (sys_is_dir(path)) fatal("%s is a directory", path);
  f = sys_fopen(path, "rb");
  if (!f) fatal("cannot open %s: %s", path, strerror(errno));
  return f;
}

/*  Unless `force' is set, the output file is created exclusively, so that an existing file, or a symbolic
    link planted in its place, is never written to.  The new file gets the access permissions of the input,
    if known.  */
static FILE * open_output(const char * path, int force, FILE * in, const bz3_meta * meta) {
  struct stat so;
  FILE * f;
  int fd;
  if (!path) return stdout;
  if (sys_is_dir(path)) fatal("%s is a directory", path);

#ifdef _WIN32
  if (in != stdin && bz3_win_same_file(in, path))
    fatal("%s is both input and output", path);

  fd = bz3_win_create(path, force);
#else
  { struct stat si;
    const int have_si = in != stdin && !fstat(fileno(in), &si);
    const int perm = meta->known && meta->regular ? 0600 : 0666;

    if (have_si && si.st_ino && !stat(path, &so) && si.st_dev == so.st_dev &&
        si.st_ino == so.st_ino)
      fatal("%s is both input and output", path);

    fd = open(path, O_WRONLY | O_CREAT | O_BINARY | (force ? O_TRUNC : O_EXCL),
              perm); }
#endif
  if (fd == -1 && errno == EEXIST && !force)
    fatal("%s already exists; use -f to overwrite it", path);

  f = fd == -1 ? NULL : fdopen(fd, "wb");
  if (!f) fatal("cannot open %s: %s", path, strerror(errno));
  setvbuf(f, NULL, _IOFBF, BZ3_OUT_BUFFER);

  /*  Only regular files are cleaned up after a failure; never devices and the like.  */
  if (!fstat(fd, &so) && S_ISREG(so.st_mode)) { partial = path;  partial_file = f; }
  return f;
}

/*  The name of the output file for an input file.  NULL if a file to decompress lacks the .bz3 extension.  */
static char * out_name(const char * in, int mode) {
  sz n = strlen(in);
  char * s = xcalloc(n + 5, 1);
  memcpy(s, in, n);
  if (mode == M_ENCODE) memcpy(s + n, ".bz3", 5);
  else if (n > 4 && !strcmp(in + n - 4, ".bz3")) s[n - 4] = 0;
  else {
    msg("bzip3: %s: unknown extension, skipping\n", in);
    free(s);
    return NULL;
  }
  return s;
}

typedef struct {
  int mode, force, to_stdout, jobs, batch, verbose, rm;
  i32 block_size;
} options;

/*  Process one file.  A NULL input or output stands for the standard streams, and so does "-".  */
static int do_one(const options * o, const char * in, const char * out) {
  FILE * fi, * fo;
  bz3_meta meta;
  int status, changed;
  if (in && !strcmp(in, "-")) in = NULL;
  if (out && !strcmp(out, "-")) out = NULL;

  fi = open_input(in);
  memset(&meta, 0, sizeof meta);
  if (in) bz3_meta_get(&meta, fi);
  fo = o->mode != M_TEST ? open_output(out, o->force, fi, &meta) : NULL;

  status = process(fi, fo, o->mode, o->block_size, o->jobs, o->verbose, in);

  if (status) {
    if (fi != stdin) fclose(fi);
    if (fo && fo != stdout && fo != partial_file) fclose(fo);
    partial_remove();
    return status;
  }

  changed = in && meta.regular && bz3_meta_changed(&meta, fi);
  close_output(fo, &meta, fi, in, out);
  if (fi != stdin) fclose(fi);
  if (changed) {
    msg("bzip3: %s: the file changed while it was read%s\n", in,
        o->rm && out && o->mode != M_TEST ? "; it is kept" : "");
    return 1;
  }

  /*  The input is removed only once the output is complete and in a file.  */
  if (o->rm && in && out && o->mode != M_TEST && sys_remove(in))
    fatal("cannot remove %s: %s", in, strerror(errno));

  return 0;
}

static int do_batch(const options * o, const yarg_result * r) {
  int status = 0, i;
  Fi(r->pos_argc,
    const char * in = r->pos_args[i];
    char * out = NULL;
    if (o->mode != M_TEST && !o->to_stdout && strcmp(in, "-")) {
      out = out_name(in, o->mode);
      if (!out) { status = 1;  continue; }
    }
    if (do_one(o, in, out)) status = 1;
    free(out));
  return status;
}

enum { OPT_RM = 256 };

/*  Carry out what the parsed command line asks for.  Returns the exit status.  */
static int run(const yarg_result * r) {
  options o;
  const char * in, * out;
  char * derived = NULL;
  long v;
  int status, i;

  memset(&o, 0, sizeof o);
  o.mode = M_ENCODE;
  o.block_size = BZ3_DEFAULT_BLOCK;

  Fi(r->argc,
    const char * arg = r->args[i].arg;
    switch (r->args[i].opt) {
    case 'e': case 'z': o.mode = M_ENCODE;  break;
    case 'd': o.mode = M_DECODE;  break;
    case 'r': o.mode = M_RECOVER;  break;
    case 't': o.mode = M_TEST;  break;
    case 'c': o.to_stdout = 1;  break;
    case 'f': o.force = 1;  break;
    case OPT_RM: o.rm = 1;  break;
    case 'k': o.rm = 0;  break;
    case 'B': o.batch = 1;  break;
    case 'v': o.verbose = 1;  break;
    case 'h': usage(stdout);  return 0;
    case 'V': version();  return 0;
    case 'b':
      if (!parse_number(arg, BZ3_MAX_BLOCK / MiB(1), &v) || v < 1) {
        msg("bzip3: the block size is 1 to %d MiB, not '%s'\n",
                BZ3_MAX_BLOCK / MiB(1), arg);
        return 1;
      }
      o.block_size = MiB((i32) v);
      break;
    case 'j':
      if (!parse_number(arg, BZ3_MAX_JOBS, &v)) {
        msg("bzip3: the number of jobs is 0 to %d, not '%s'\n",
                BZ3_MAX_JOBS, arg);
        return 1;
      }
      o.jobs = (int) v;
      break;
    });

  bz3_map_guard(&partial);
  setvbuf(stdout, NULL, _IOFBF, BZ3_OUT_BUFFER);

#if defined(_WIN32) && defined(HAVE_IO_H)
  setmode(fileno(stdin), O_BINARY);
  setmode(fileno(stdout), O_BINARY);
#endif

  if (o.batch && r->pos_argc) status = do_batch(&o, r);
  else {
    if (r->pos_argc > 2) {
      msg("bzip3: too many files; use -B to process each of them\n");
      return 1;
    }

    /*  With one file name the output name is derived from it; the second one names the output.  */
    in = r->pos_argc > 0 ? r->pos_args[0] : NULL;
    out = r->pos_argc > 1 ? r->pos_args[1] : NULL;
    if (in && strcmp(in, "-") && !out && o.mode != M_TEST && !o.to_stdout) {
      out = derived = out_name(in, o.mode);
      if (!out) return 1;
    }

    status = do_one(&o, in, out);
    free(derived);
  }

  if (fclose(stdout)) {
    msg("bzip3: write error: %s\n", strerror(errno));
    return 1;
  }

  return status;
}

static int bz3_main(int argc, char ** argv) {
  yarg_options opt[] = {
    { 'e', no_argument, "encode" },
    { 'z', no_argument, NULL },
    { 'd', no_argument, "decode" },
    { 't', no_argument, "test" },
    { 'r', no_argument, "recover" },
    { 'c', no_argument, "stdout" },
    { 'f', no_argument, "force" },
    { OPT_RM, no_argument, "rm" },
    { 'k', no_argument, "keep" },
    { 'h', no_argument, "help" },
    { 'V', no_argument, "version" },
    { 'v', no_argument, "verbose" },
    { 'b', required_argument, "block" },
    { 'B', no_argument, "batch" },
    { 'j', required_argument, "jobs" },
    { 0, no_argument, NULL }
  };
  yarg_settings set;
  yarg_result * r;
  int status;

  set.dash_dash = true;
  set.style = YARG_STYLE_UNIX;

  r = yarg_parse(argc, argv, opt, set);
  if (!r || r->error) {
    msg("bzip3: %s\n", r ? r->error : "out of memory");
    if (r) msg("Try 'bzip3 --help' for more information.\n");
    yarg_destroy(r);
    return 1;
  }

  status = run(r);
  yarg_destroy(r);
  return status;
}

#if defined(_WIN32) && !defined(BZ3_WIN_LEGACY)
int wmain(int argc, wchar_t ** wargv);

int wmain(int argc, wchar_t ** wargv) {
  char ** argv = bz3_win_args(argc, wargv);
  if (!argv) fatal("out of memory");
  return bz3_main(argc, argv);
}
#else
int main(int argc, char ** argv) { return bz3_main(argc, argv); }
#endif
