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
#ifndef _WIN32
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#endif

static const char * bz3_of(const char * path) {
  static char buf[4][160];
  static int next;
  char * s = buf[next++ & 3];
  snprintf(s, sizeof buf[0], "%s.bz3", path);
  return s;
}

static int exists(const char * path) { return xt_file_size(path) >= 0; }

static void t_options(void) {
  const char * log = xt_tmp("cli.log"), * in = xt_fixture("plain.txt");
  static const char * const bad[] = {
    "--help=x", "--nonsense", "-b", "--r", "-b 0", "-b 512", "-b 4097",
    "-b 99999999999999999999", "-b x", "-b -1", "-b 1x", "-b +1", "-j 65",
    "-j x", "-j -1", NULL
  };
  int i;
  xt_section_begin("cli options");
  CHECK(!xt_run(log, "--help") && xt_file_contains(log, "usage: bzip3"), "--help");
  CHECK(!xt_run(log, "-h") && xt_file_contains(log, "--batch"), "-h");
  CHECK(!xt_run(log, "--version") && xt_file_contains(log, "bzip3 ") &&
        xt_file_contains(log, bz3_version()), "--version");
  CHECK(!xt_run(log, "-V") && xt_file_contains(log, "Library"), "-V");
  for (i = 0; bad[i]; i++)
    CHECK(xt_run(log, "%s -c \"%s\"", bad[i], in) == 1, "'%s' is refused", bad[i]);
  CHECK(!xt_run(log, "--enc --std \"%s\"", in), "unambiguous prefixes of long options work");
  CHECK(!xt_run(log, "-z -c -j 0 -b 511 \"%s\"", xt_fixture("stored.bin")),
        "the limits of -b and -j are accepted");
  CHECK(xt_run(log, "a b c") == 1 && xt_file_contains(log, "too many files"),
        "three files need -B");
  CHECK(xt_run(log, "-c t_suite-no-such-file.tmp") == 1 && xt_file_contains(log, "cannot open"),
        "a missing input fails");
  CHECK(xt_run(log, "-d \"%s\"", in) == 1 && xt_file_contains(log, "unknown extension"),
        "decoding needs the .bz3 extension to name the output");
  CHECK(!xt_run(log, "-v -c \"%s\"", in) && xt_file_contains(log, " -> ") &&
        xt_file_contains(log, "bpb"), "-v reports the sizes");
  xt_unlink(log);
}

/*  The reference files, from the command line.  */
static void t_golden(void) {
  const char * log = xt_tmp("cli.log"), * out = xt_tmp("gold.out");
  const xt_golden * g;
  char packed[4096];
  xt_section_begin("cli reference files");
  for (g = XT_GOLDEN; g->name; g++) {
    const char * orig = xt_fixture(g->name);
    snprintf(packed, sizeof packed, "%s.bz3", orig);
    CHECK(!xt_run(log, "-t \"%s\"", packed), "%s: -t accepts it", g->name);
    CHECK(!xt_run(out, "-dc \"%s\"", packed) && xt_same_file(out, orig),
          "%s: decodes to the standard output", g->name);
    xt_unlink(out);
    CHECK(!xt_run(log, "-d \"%s\" \"%s\"", packed, out) && xt_same_file(out, orig),
          "%s: decodes to a file", g->name);
    xt_unlink(out);
    CHECK(!xt_run(out, "-dc -j 4 \"%s\"", packed) && xt_same_file(out, orig),
          "%s: decodes with -j", g->name);
    xt_unlink(out);
  }
  xt_unlink(log);
}

/*  A few megabytes made of the reference files, so that small blocks make many of them.  */
static void make_big(const char * path) {
  const xt_golden * g;
  FILE * f = fopen(path, "wb");
  int rep;
  if (!f) { CHECK(0, "cannot create %s", path);  return; }
  for (rep = 0; rep < 140 * xt_level; rep++)
    for (g = XT_GOLDEN; g->name; g++) {
      sz n;
      u8 * b = xt_slurp(xt_fixture(g->name), &n);
      if (b) fwrite(b, 1, n, f);
      free(b);
    }
  fclose(f);
}

static void t_round_trips(void) {
  static const char * const flags[] = { "-b 1", "-b 1 -j 3", "-b 3 -j 3", "-j 64 -b 1", "", NULL };
  const char * log = xt_tmp("cli.log"), * big = xt_tmp("big"), * out = xt_tmp("big.out");
  const char * packed = xt_tmp("big.bz3"), * one = xt_tmp("big.one");
  int i;
  xt_section_begin("cli round trips");
  make_big(big);
  CHECK(!xt_run(one, "-b 1 -c \"%s\"", big), "encode with one job");
  for (i = 0; flags[i]; i++) {
    const char * fl = flags[i];
    CHECK(!xt_run(packed, "-e %s -c \"%s\"", fl, big), "encode (%s)", fl);
    if (strstr(fl, "-b 1"))
      CHECK(xt_same_file(packed, one), "the output does not depend on the jobs (%s)", fl);
    CHECK(!xt_run(log, "-t \"%s\"", packed), "test (%s)", fl);
    CHECK(!xt_run(log, "-t -j 5 \"%s\"", packed), "test with jobs (%s)", fl);
    CHECK(!xt_run(out, "-d -j 3 -c \"%s\"", packed) && xt_same_file(out, big),
          "round trip (%s)", fl);
    CHECK(!xt_run(out, "-dc < \"%s\"", packed) && xt_same_file(out, big),
          "round trip through the standard input (%s)", fl);
  }

  /*  "-" stands for the standard streams; bz3cat and bz3grep rely on it.  */
  CHECK(!xt_run(packed, "-c - < \"%s\"", big) && !xt_run(out, "-dc - < \"%s\"", packed) &&
        xt_same_file(out, big), "round trip with - as the file name");
  CHECK(!xt_run(out, "-Bcd - < \"%s\"", packed) && xt_same_file(out, big),
        "batch mode with - as the file name");
  CHECK(!xt_run(out, "-d \"%s\" -", packed) && xt_same_file(out, big),
        "- as the output file name");
  CHECK(!exists("-"), "no file called - is created");

  xt_unlink(big);  xt_unlink(out);  xt_unlink(packed);  xt_unlink(one);  xt_unlink(log);
}

static void t_empty(void) {
  const char * log = xt_tmp("cli.log"), * in = xt_tmp("empty"), * out = xt_tmp("empty.out");
  const char * packed = xt_tmp("empty.bz3");
  xt_section_begin("cli empty input");
  xt_spew(in, (const u8 *) "", 0);
  CHECK(!xt_run(packed, "-c \"%s\"", in) && xt_file_size(packed) == 9,
        "an empty file becomes a bare header");
  CHECK(!xt_run(packed, "-j 4 -c \"%s\"", in) && xt_file_size(packed) == 9,
        "with jobs too");
  CHECK(!xt_run(out, "-dc \"%s\"", packed) && xt_file_size(out) == 0, "and decodes to nothing");
  CHECK(!xt_run(log, "-t \"%s\"", packed), "-t accepts it");
  xt_unlink(in);  xt_unlink(out);  xt_unlink(packed);  xt_unlink(log);
}

/*  Existing files, input removal and what a failure leaves behind.  */
static void t_files(void) {
  const char * log = xt_tmp("cli.log"), * in = xt_tmp("plain"), * out = xt_tmp("plain.out");
  const char * packed = bz3_of(in), * cut = bz3_of(xt_tmp("cut")), * hurt = xt_tmp("hurt.bz3");
  const char * orig = xt_fixture("plain.txt");
  sz n, pn, i;
  u8 * data = xt_slurp(orig, &n), * p;
  xt_section_begin("cli files");
  if (!data) { CHECK(0, "no plain.txt");  return; }
  xt_spew(in, data, n);
  xt_spew(packed, (const u8 *) "keep\n", 5);

  CHECK(xt_run(log, "\"%s\"", in) == 1 && xt_file_contains(log, "already exists") &&
        xt_file_size(packed) == 5, "an existing file is not overwritten");
  CHECK(!xt_run(log, "-f \"%s\"", in) && !xt_run(out, "-dc \"%s\"", packed) &&
        xt_same_file(out, orig), "-f overwrites it");
  CHECK(xt_run(log, "-f \"%s\" \"%s\"", in, in) == 1 && xt_same_file(in, orig),
        "a file is not compressed onto itself");

  /*  --rm removes the input once the operation succeeded...  */
  CHECK(!xt_run(log, "-f --rm \"%s\"", in) && !exists(in), "--rm removes the input");
  CHECK(!xt_run(log, "-d --rm \"%s\"", packed) && exists(in) && !exists(packed) &&
        xt_same_file(in, orig), "-d --rm restores the file");
  CHECK(!xt_run(log, "--rm -k \"%s\"", in) && exists(in), "-k after --rm keeps the input");
  CHECK(!xt_run(log, "--rm -t \"%s\"", packed) && exists(packed), "--rm -t removes nothing");
  CHECK(!xt_run(out, "--rm -dc \"%s\"", packed) && exists(packed),
        "--rm with the standard output removes nothing");

  /*  ... and never after a failure, which leaves no partial output behind either.  */
  p = xt_slurp(packed, &pn);
  if (p && pn > 60) {
    xt_spew(cut, p, pn - 7);
    CHECK(xt_run(log, "-d --rm \"%s\"", cut) == 1, "a truncated file is refused");
    CHECK(exists(cut), "--rm keeps the input after a failure");
    CHECK(!exists(xt_tmp("cut")), "no partial output is left behind");
    for (i = 10; i < pn; i += pn / 37 + 1) {
      xt_spew(cut, p, i);
      CHECK(xt_run(log, "-t \"%s\"", cut) == 1, "a prefix of %lu bytes is refused",
            (unsigned long) i);
    }

    /*  Damaged data is detected, and the recovery mode writes it out all the same.  */
    memcpy(p + 40, "XXXX", 4);
    xt_spew(hurt, p, pn);
    CHECK(xt_run(log, "-t \"%s\"", hurt) == 1, "-t refuses damaged data");
    CHECK(xt_run(out, "-dc \"%s\"", hurt) == 1, "-d refuses damaged data");
    xt_unlink(out);
    CHECK(!xt_run(log, "-r \"%s\" \"%s\"", hurt, out) && xt_file_size(out) == (long) n &&
          xt_file_contains(log, "invalid"), "-r writes the damaged block out");
    xt_unlink(out);
  } else CHECK(0, "no compressed file to damage");
  free(p);

  xt_spew(hurt, (const u8 *) "BZ2v1xxxxxxxxxxxxxxxxx", 22);
  CHECK(xt_run(log, "-dc \"%s\"", hurt) == 1 && xt_file_contains(log, "not a bzip3 file"),
        "a bad signature is refused");
  xt_spew(hurt, (const u8 *) "BZ3v1\377\377\377\177", 9);
  CHECK(xt_run(log, "-dc \"%s\"", hurt) == 1 && xt_file_contains(log, "block size"),
        "a bad block size is refused");
  xt_spew(hurt, (const u8 *) "BZ3v1\0\0\20\0\377\377\377\177\1\0\0\0", 17);
  CHECK(xt_run(log, "-dc \"%s\"", hurt) == 1, "a bad block header is refused");

#ifndef _WIN32
  /*  The output gets the permissions of the input, and a symbolic link in its place is not followed.  */
  { struct stat st;
    const char * target = xt_tmp("target");
    xt_unlink(packed);
    CHECK(!chmod(in, 0640) && !xt_run(log, "\"%s\"", in) && !stat(packed, &st) &&
          (st.st_mode & 0777) == 0640, "the output has the mode of the input");
    xt_unlink(packed);
    { struct timeval tv[2];
      const mode_t mask = umask(077);
      tv[0].tv_sec = 1000000000;  tv[0].tv_usec = 0;
      tv[1].tv_sec = 1200000000;  tv[1].tv_usec = 0;
      CHECK(!chmod(in, 02644) && !utimes(in, tv) && !xt_run(log, "-k \"%s\"", in) &&
            !stat(packed, &st) && (st.st_mode & 07777) == 0644 && st.st_mtime == 1200000000 &&
            st.st_atime == 1000000000, "the output has the times of the input and no special bits");
      xt_unlink(out);
      CHECK(!xt_run(log, "-d \"%s\" \"%s\"", packed, out) && !stat(out, &st) &&
            (st.st_mode & 07777) == 0644 && st.st_mtime == 1200000000,
            "and so does what is decoded from it");
      xt_unlink(out);
      CHECK(!xt_run(log, "-c \"%s\" > \"%s\"", in, out) && !stat(out, &st) &&
            st.st_mtime != 1200000000, "the standard output is left as it is");
      umask(mask);  chmod(in, 0644); }
    xt_unlink(out);  xt_unlink(packed);
    xt_spew(target, (const u8 *) "precious", 8);
    CHECK(!symlink(target, packed) && xt_run(log, "\"%s\"", in) == 1 &&
          xt_file_size(target) == 8, "a planted symbolic link is not written through");
    xt_unlink(packed);  xt_unlink(target); }
#endif

  free(data);
  xt_unlink(in);  xt_unlink(out);  xt_unlink(packed);  xt_unlink(cut);  xt_unlink(hurt);
  xt_unlink(log);
}

static void t_batch(void) {
  const char * log = xt_tmp("cli.log"), * a = xt_tmp("batch-a"), * b = xt_tmp("batch-b");
  const char * pa = bz3_of(a), * pb = bz3_of(b);
  const char * oa = xt_fixture("plain.txt"), * ob = xt_fixture("rle.bin");
  sz n;
  u8 * d;
  xt_section_begin("cli batch");
  d = xt_slurp(oa, &n);  if (d) xt_spew(a, d, n);  free(d);
  d = xt_slurp(ob, &n);  if (d) xt_spew(b, d, n);  free(d);
  CHECK(!xt_run(log, "-B --rm \"%s\" \"%s\"", a, b) && exists(pa) && exists(pb) &&
        !exists(a) && !exists(b), "batch encode");
  CHECK(!xt_run(log, "-B -t \"%s\" \"%s\"", pa, pb), "batch test");
  CHECK(!xt_run(log, "-Bd \"%s\" \"%s\"", pa, pb) && xt_same_file(a, oa) &&
        xt_same_file(b, ob), "batch decode");
  /*  A failure does not stop the batch, but it shows in the exit status, and --rm spares its input.  */
  xt_unlink(a);
  CHECK(xt_run(log, "-Bd --rm \"%s\" \"%s\"", pa, pb) == 1 && xt_same_file(a, oa) &&
        !exists(pa) && exists(pb), "a failed file fails the batch and stays");
  CHECK(xt_run(log, "-Bd \"%s\" \"%s\"", a, pb) == 1 && xt_file_contains(log, "unknown extension"),
        "a file without the extension fails the batch");
  xt_unlink(a);  xt_unlink(b);  xt_unlink(pa);  xt_unlink(pb);  xt_unlink(log);
}

void xt_run_cli(void) {
  if (!xt_binary || !xt_data) {
    xt_section_begin("cli");
    CHECK(0, "BZ3 and BZ3_TEST_DATA name the binary and the test data");
    return;
  }
  t_options();
  t_golden();
  t_round_trips();
  t_empty();
  t_files();
  t_batch();
}
