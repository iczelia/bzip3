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

/*  Without libFuzzer: run a harness on each file named on the command line, such as the inputs of earlier
    findings.  */

#include "harness.h"

int main(int argc, char ** argv) {
  int i, done = 0;
  for (i = 1; i < argc; i++) {
    FILE * f = fopen(argv[i], "rb");
    u8 * b;
    long n;
    if (!f) { fprintf(stderr, "fzmain: cannot open %s\n", argv[i]);  return 2; }
    fseek(f, 0, SEEK_END);  n = ftell(f);  fseek(f, 0, SEEK_SET);
    b = fz_alloc((sz) (n > 0 ? n : 0));
    n = n > 0 ? (long) fread(b, 1, (sz) n, f) : 0;
    fclose(f);
    LLVMFuzzerTestOneInput(b, (sz) n);
    free(b);
    done++;
  }
  fprintf(stderr, "fzmain: replayed %d inputs\n", done);
  return 0;
}
