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

/*  The regular file being written, if any.  It is removed when the operation fails, so that no truncated
    output is left behind.  */
static const char * partial;

static void partial_remove(void) {
  if (partial) { remove(partial);  partial = NULL; }
}

static NORETURN void fatal(const char * fmt, ...) BZ3_PRINTF(1, 2);

static void fatal(const char * fmt, ...) {
  va_list ap;
  fputs("bzip3: ", stderr);
  va_start(ap, fmt);  vfprintf(stderr, fmt, ap);  va_end(ap);
  fputc('\n', stderr);
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

static void xwrite(const void * data, sz n, FILE * f) {
  if (n && fwrite(data, 1, n, f) != n) fatal("write error: %s", strerror(errno));
}

/*  Read up to `n' bytes; fewer are returned only at the end of the file.  */
static sz xread(void * data, sz n, FILE * f) {
  sz got = fread(data, 1, n, f);
  if (ferror(f)) fatal("read error: %s", strerror(errno));
  return got;
}

static void xread_exact(void * data, sz n, FILE * f) {
  if (xread(data, n, f) != n) fatal("the input is truncated");
}

static void close_output(FILE * f) {
  if (!f) return;
  if (fflush(f)) fatal("write error: %s", strerror(errno));
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
  if (f != stdout && fclose(f)) fatal("write error: %s", strerror(errno));
  partial = NULL;
}

/*  The buffers and the coder states of all the jobs.  */
typedef struct {
  int n;
  struct bz3_state ** states;
  u8 ** buffers;
  sz * caps;
  i32 * sizes, * orig_sizes;
} pool;

static void pool_init(pool * p, int n, i32 block_size) {
  int i;
  p->n = n;
  p->states = xcalloc((sz) n, sizeof *p->states);
  p->buffers = xcalloc((sz) n, sizeof *p->buffers);
  p->caps = xcalloc((sz) n, sizeof *p->caps);
  p->sizes = xcalloc((sz) n, sizeof *p->sizes);
  p->orig_sizes = xcalloc((sz) n, sizeof *p->orig_sizes);
  Fi(n,
    p->states[i] = bz3_new(block_size);
    if (!p->states[i]) fatal("out of memory");
    p->caps[i] = bz3_bound((sz) block_size);
    /*  Zeroed, so that the recovery mode never writes out uninitialised memory.  */
    p->buffers[i] = xcalloc(p->caps[i], 1));
}

static void pool_free(pool * p) {
  int i;
  Fi(p->n, bz3_free(p->states[i]);  free(p->buffers[i]));
  free(p->states);  free(p->buffers);  free(p->caps);
  free(p->sizes);  free(p->orig_sizes);
}

/*  Code the first `n' blocks held by the pool.  */
static void pool_run(pool * p, int n, int enc) {
  if (n > 1) {
    if (enc) bz3_encode_blocks(p->states, p->buffers, p->sizes, n);
    else bz3_decode_blocks(p->states, p->buffers, p->caps, p->sizes,
                           p->orig_sizes, n);
  } else if (n == 1) {
    if (enc) p->sizes[0] = bz3_encode_block(p->states[0], p->buffers[0], p->sizes[0]);
    else bz3_decode_block(p->states[0], p->buffers[0], p->caps[0], p->sizes[0],
                          p->orig_sizes[0]);
  }
}

static int process(FILE * in, FILE * out, int mode, i32 block_size, int jobs,
                   int verbose, const char * name) {
  const int enc = mode == M_ENCODE;
  uint64_t n_in = 0, n_out = 0;
  u8 header[9];
  pool p;
  int status = 0, eof = 0, n, i;

  if (isatty(fileno(enc ? out : in))) {
    fprintf(stderr, "bzip3: compressed data is not %s a terminal\n",
            enc ? "written to" : "read from");
    return 1;
  }

  if (enc) {
    memcpy(header, "BZ3v1", 5);
    wr32le(header + 5, (u32) block_size);
    xwrite(header, 9, out);
    n_out += 9;
  } else {
    if (xread(header, 5, in) != 5 || memcmp(header, "BZ3v1", 5)) {
      fprintf(stderr, "bzip3: %s: not a bzip3 file\n", name ? name : "(stdin)");
      return 1;
    }
    xread_exact(header + 5, 4, in);
    block_size = (i32) rd32le(header + 5);
    if (block_size < BZ3_MIN_BLOCK || block_size > BZ3_MAX_BLOCK) {
      fprintf(stderr, "bzip3: %s: invalid block size in the header\n",
              name ? name : "(stdin)");
      if (mode != M_RECOVER) return 1;
      fprintf(stderr, "bzip3: recovering with the largest block size\n");
      block_size = BZ3_MAX_BLOCK;
    }
    n_in += 9;
  }

  if (jobs < 1) jobs = 1;
#ifndef HAVE_THREADS
  jobs = 1;
#endif
  pool_init(&p, jobs, block_size);

  while (!eof && !status) {
    /*  Fill as many buffers as there is data for.  */
    for (n = 0; n < jobs; n++) {
      if (enc) {
        sz got = xread(p.buffers[n], (sz) block_size, in);
        eof = got < (sz) block_size;
        if (!got) break;
        p.sizes[n] = p.orig_sizes[n] = (i32) got;
        n_in += got;
        if (eof) { n++;  break; }
      } else {
        sz got = xread(header, 8, in);
        if (!got) { eof = 1;  break; }
        if (got != 8) fatal("the input is truncated");
        p.sizes[n] = (i32) rd32le(header);
        p.orig_sizes[n] = (i32) rd32le(header + 4);
        if (p.sizes[n] < 0 || (sz) p.sizes[n] > bz3_bound((sz) block_size) ||
            p.orig_sizes[n] < 0 || p.orig_sizes[n] > block_size) {
          fprintf(stderr, "bzip3: failed to decode a block: inconsistent headers\n");
          status = 1;
          break;
        }
        xread_exact(p.buffers[n], (sz) p.sizes[n], in);
        n_in += 8 + (uint64_t) p.sizes[n];
      }
    }
    if (status) break;

    pool_run(&p, n, enc);

    Fi(n,
      if (bz3_last_error(p.states[i]) != BZ3_OK) {
        fprintf(stderr, "bzip3: %s a block: %s\n",
                enc ? "failed to encode" : mode == M_RECOVER
                  ? "writing an invalid" : "failed to decode",
                bz3_strerror(p.states[i]));
        if (mode != M_RECOVER) { status = 1;  break; }
      }
      if (enc) {
        wr32le(header, (u32) p.sizes[i]);
        wr32le(header + 4, (u32) p.orig_sizes[i]);
        xwrite(header, 8, out);
        xwrite(p.buffers[i], (sz) p.sizes[i], out);
        n_out += 8 + (uint64_t) p.sizes[i];
      } else {
        if (mode != M_TEST) xwrite(p.buffers[i], (sz) p.orig_sizes[i], out);
        n_out += (uint64_t) p.orig_sizes[i];
      });
  }

  pool_free(&p);
  if (status) return status;

  if (out && fflush(out)) fatal("write error: %s", strerror(errno));

  if (verbose) {
    /*  The size of the compressed data relative to the size of the original.  */
    const double packed = (double) (enc ? n_out : n_in);
    const double unpacked = (double) (enc ? n_in : n_out);
    if (name) fprintf(stderr, " %s:", name);
    fprintf(stderr, "\t%s%" PRIu64 " -> %" PRIu64 " bytes, %.2f%%, %.2f bpb\n",
            mode == M_ENCODE || mode == M_DECODE ? "" : "OK, ", n_in, n_out,
            packed * 100.0 / unpacked, packed * 8.0 / unpacked);
  }

  return 0;
}

static int is_dir(const char * path) {
  struct stat st;
  return !stat(path, &st) && S_ISDIR(st.st_mode);
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
  if (is_dir(path)) fatal("%s is a directory", path);
  f = fopen(path, "rb");
  if (!f) fatal("cannot open %s: %s", path, strerror(errno));
  return f;
}

/*  Unless `force' is set, the output file is created exclusively, so that an existing file, or a symbolic
    link planted in its place, is never written to.  The new file gets the access permissions of the input,
    if known.  */
static FILE * open_output(const char * path, int force, FILE * in) {
  struct stat si, so;
  FILE * f;
  int have_si, perm, fd;
  if (!path) return stdout;
  if (is_dir(path)) fatal("%s is a directory", path);

  have_si = in != stdin && !fstat(fileno(in), &si);
  perm = have_si && S_ISREG(si.st_mode) ? (int) (si.st_mode & 0777) : 0666;

  if (have_si && si.st_ino && !stat(path, &so) && si.st_dev == so.st_dev &&
      si.st_ino == so.st_ino)
    fatal("%s is both input and output", path);

  fd = open(path, O_WRONLY | O_CREAT | O_BINARY | (force ? O_TRUNC : O_EXCL),
            perm);
  if (fd == -1 && errno == EEXIST && !force)
    fatal("%s already exists; use -f to overwrite it", path);

  f = fd == -1 ? NULL : fdopen(fd, "wb");
  if (!f) fatal("cannot open %s: %s", path, strerror(errno));

  /*  Only regular files are cleaned up after a failure; never devices and the like.  */
  if (!fstat(fd, &so) && S_ISREG(so.st_mode)) partial = path;
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
    fprintf(stderr, "bzip3: %s: unknown extension, skipping\n", in);
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
  int status;
  if (in && !strcmp(in, "-")) in = NULL;
  if (out && !strcmp(out, "-")) out = NULL;

  fi = open_input(in);
  fo = o->mode != M_TEST ? open_output(out, o->force, fi) : NULL;

  status = process(fi, fo, o->mode, o->block_size, o->jobs, o->verbose, in);

  if (fi != stdin) fclose(fi);

  if (status) {
    if (fo && fo != stdout) fclose(fo);
    partial_remove();
    return status;
  }

  close_output(fo);

  /*  The input is removed only once the output is complete and in a file.  */
  if (o->rm && in && out && o->mode != M_TEST && remove(in))
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
        fprintf(stderr, "bzip3: the block size is 1 to %d MiB, not '%s'\n",
                BZ3_MAX_BLOCK / MiB(1), arg);
        return 1;
      }
      o.block_size = MiB((i32) v);
      break;
    case 'j':
      if (!parse_number(arg, BZ3_MAX_JOBS, &v)) {
        fprintf(stderr, "bzip3: the number of jobs is 0 to %d, not '%s'\n",
                BZ3_MAX_JOBS, arg);
        return 1;
      }
      o.jobs = (int) v;
      break;
    });

#if defined(_WIN32) && defined(HAVE_IO_H)
  setmode(fileno(stdin), O_BINARY);
  setmode(fileno(stdout), O_BINARY);
#endif

  if (o.batch && r->pos_argc) status = do_batch(&o, r);
  else {
    if (r->pos_argc > 2) {
      fprintf(stderr, "bzip3: too many files; use -B to process each of them\n");
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
    fprintf(stderr, "bzip3: write error: %s\n", strerror(errno));
    return 1;
  }

  return status;
}

int main(int argc, char ** argv) {
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
    fprintf(stderr, "bzip3: %s\n", r ? r->error : "out of memory");
    if (r) fprintf(stderr, "Try 'bzip3 --help' for more information.\n");
    yarg_destroy(r);
    return 1;
  }

  status = run(r);
  yarg_destroy(r);
  return status;
}
