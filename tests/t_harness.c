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
#include <stdarg.h>

#ifdef HAVE_SYS_WAIT_H
#include <sys/wait.h>
#endif

int xt_level = 1;
int xt_tracing;
unsigned long xt_checks, xt_failures;
const char * xt_section = "";
const char * xt_data;
const char * xt_binary;

const xt_golden XT_GOLDEN[] = {
  { "stored.bin", -1 }, { "plain.txt", 0 }, { "lzp.txt", 2 }, { "rle.bin", 4 },
  { "rle-lzp.bin", 6 }, { "binary.bin", 0 }, { NULL, 0 }
};

static const char * env(const char * name) {
  const char * v = getenv(name);
  return v && *v ? v : NULL;
}

void xt_init(void) {
  const char * v = env("BZ3_TEST_LEVEL");
  if (v) {
    if (!strcmp(v, "quick")) xt_level = 1;
    else if (!strcmp(v, "full")) xt_level = 4;
    else if (!strcmp(v, "torture")) xt_level = 16;
  }
  v = env("BZ3_TEST_TRACE");
  xt_tracing = v && strcmp(v, "0");
  xt_data = env("BZ3_TEST_DATA");
  xt_binary = env("BZ3");
}

void xt_section_begin(const char * name) {
  xt_section = name;
  if (xt_tracing) fprintf(stderr, "  == %s\n", name);
}

void xt_trace(const char * fmt, ...) {
  va_list ap;
  if (!xt_tracing) return;
  va_start(ap, fmt);
  fputs("     ", stderr);  vfprintf(stderr, fmt, ap);  fputc('\n', stderr);
  va_end(ap);
}

void xt_report(int ok, const char * fmt, ...) {
  va_list ap;
  xt_checks++;
  if (ok) return;
  xt_failures++;
  fprintf(stderr, "FAIL [%s] ", xt_section);
  va_start(ap, fmt);  vfprintf(stderr, fmt, ap);  va_end(ap);
  fputc('\n', stderr);
}

int xt_finish(const char * program) {
  fprintf(stderr, "%s: %lu checks, %lu failed\n", program, xt_checks,
          xt_failures);
  return xt_failures ? 1 : 0;
}

void * xt_alloc(sz n) {
  void * p = malloc(n ? n : 1);
  if (!p) { fprintf(stderr, "t_suite: out of memory\n");  exit(99); }
  return p;
}

u8 * xt_slurp(const char * path, sz * len) {
  FILE * f = fopen(path, "rb");
  u8 * b;
  long n;
  *len = 0;
  if (!f) return NULL;
  fseek(f, 0, SEEK_END);  n = ftell(f);  fseek(f, 0, SEEK_SET);
  b = xt_alloc((sz) (n > 0 ? n : 0));
  if (n > 0) *len = fread(b, 1, (sz) n, f);
  fclose(f);
  return b;
}

void xt_spew(const char * path, const u8 * data, sz len) {
  FILE * f = fopen(path, "wb");
  if (!f || fwrite(data, 1, len, f) != len || fclose(f)) {
    fprintf(stderr, "t_suite: cannot write %s\n", path);  exit(99);
  }
}

const char * xt_fixture(const char * name) {
  static char buf[4][4096];
  static int next;
  char * s = buf[next++ & 3];
  snprintf(s, sizeof buf[0], "%s/%s", xt_data ? xt_data : ".", name);
  return s;
}

const char * xt_tmp(const char * tag) {
  static char buf[64][128], tags[64][32];
  static int used;
  int i;
  Fi(used, if (!strcmp(tags[i], tag)) return buf[i]);
  if (used == 64 || strlen(tag) >= sizeof tags[0]) {
    fprintf(stderr, "t_suite: bad scratch name %s\n", tag);  exit(99);
  }
  strcpy(tags[used], tag);
  snprintf(buf[used], sizeof buf[used], "t_suite-%s.tmp", tag);
  return buf[used++];
}

void xt_unlink(const char * path) { remove(path); }

int xt_same_file(const char * a, const char * b) {
  sz an, bn;
  u8 * x = xt_slurp(a, &an), * y = xt_slurp(b, &bn);
  int eq = x && y && an == bn && !memcmp(x, y, an);
  free(x);  free(y);
  return eq;
}

long xt_file_size(const char * path) {
  long n;
  FILE * f = fopen(path, "rb");
  if (!f) return -1;
  fseek(f, 0, SEEK_END);  n = ftell(f);  fclose(f);
  return n;
}

int xt_file_contains(const char * path, const char * needle) {
  sz n, m = strlen(needle), i;
  u8 * b = xt_slurp(path, &n);
  if (!b) return 0;
  for (i = 0; i + m <= n; i++)
    if (!memcmp(b + i, needle, m)) { free(b);  return 1; }
  free(b);
  return 0;
}

int xt_run(const char * out, const char * args, ...) {
  char a[8192], cmd[16384];
  va_list ap;
  int st, n;
  if (!xt_binary) return -1;
  va_start(ap, args);  n = vsnprintf(a, sizeof a, args, ap);  va_end(ap);
  if (n < 0 || n >= (int) sizeof a) return -1;
#ifdef _WIN32
  if (out) n = snprintf(cmd, sizeof cmd, "\"\"%s\" %s > \"%s\" 2>&1\"",
                        xt_binary, a, out);
  else     n = snprintf(cmd, sizeof cmd, "\"\"%s\" %s\"", xt_binary, a);
#else
  if (out) n = snprintf(cmd, sizeof cmd, "\"%s\" %s > \"%s\" 2>&1",
                        xt_binary, a, out);
  else     n = snprintf(cmd, sizeof cmd, "\"%s\" %s", xt_binary, a);
#endif
  if (n < 0 || n >= (int) sizeof cmd) return -1;
  xt_trace("run: %s", cmd);
  fflush(stdout);  fflush(stderr);
  st = system(cmd);
  if (st == -1) return -1;
#if defined(WIFEXITED) && defined(WEXITSTATUS)
  if (WIFEXITED(st)) return WEXITSTATUS(st);
  return 128;
#else
  return st;
#endif
}

void xt_seed(xt_rng * r, u32 seed) { r->s = seed ? seed : 20261110UL; }

u32 xt_next(xt_rng * r, u32 n) {
  r->s = r->s * 1664525UL + 1013904223UL;
  return (r->s >> 8) % n;
}

void xt_fill(xt_rng * r, u8 * buf, sz n, int kind) {
  static const char * const words[] = {
    "block ", "sorting ", "transform ", "entropy ", "context ", "mixing\n"
  };
  sz i = 0, j;
  while (i < n) {
    int k = kind == 4 ? (int) xt_next(r, 4) : kind;
    sz len = 1 + xt_next(r, 600);
    if (len > n - i) len = n - i;
    switch (k) {
    case 0: Fj(len, buf[i + j] = (u8) xt_next(r, 256));  break;
    case 1:
      for (j = 0; j < len;) {
        const char * w = words[xt_next(r, 6)];
        for (; *w && j < len; w++) buf[i + j++] = (u8) *w;
      }
      break;
    case 2: memset(buf + i, xt_next(r, 3) ? 0 : (int) xt_next(r, 256), len);  break;
    default:
      if (i < 200) { Fj(len, buf[i + j] = (u8) ('a' + xt_next(r, 26))); }
      else {
        sz from = xt_next(r, (u32) (i - 100));
        Fj(len, buf[i + j] = buf[from + j]);
      }
    }
    i += len;
  }
}
