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

/*  An example of the high level interface: compress a file in memory and decompress it again.  With bzip3
    installed, build and run it in this directory:

      cc hl-api.c -o hl-api $(pkg-config --cflags --libs bzip3) && ./hl-api  */

#include <libbz3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MB (1024 * 1024)

int main(void) {
  FILE * f = fopen("shakespeare.txt", "rb");
  uint8_t * in, * packed, * back;
  size_t size, packed_size, back_size;
  int err;

  if (!f) { perror("shakespeare.txt");  return 1; }
  fseek(f, 0, SEEK_END);  size = (size_t) ftell(f);  fseek(f, 0, SEEK_SET);

  /*  The output buffer of bz3_compress holds at least bz3_bound(size) bytes.  */
  packed_size = bz3_bound(size);
  in = malloc(size);  packed = malloc(packed_size);  back = malloc(size);
  if (!in || !packed || !back || fread(in, 1, size, f) != size) {
    fprintf(stderr, "cannot read shakespeare.txt\n");
    return 1;
  }
  fclose(f);

  err = bz3_compress(1 * MB, in, packed, size, &packed_size);
  if (err != BZ3_OK) {
    fprintf(stderr, "bz3_compress failed with the error code %d\n", err);
    return 1;
  }
  printf("%lu => %lu\n", (unsigned long) size, (unsigned long) packed_size);

  /*  Decompress into a second buffer and compare it with the original.  */
  back_size = size;
  err = bz3_decompress(packed, back, packed_size, &back_size);
  if (err != BZ3_OK) {
    fprintf(stderr, "bz3_decompress failed with the error code %d\n", err);
    return 1;
  }
  printf("%lu => %lu\n", (unsigned long) packed_size, (unsigned long) back_size);

  if (back_size != size || memcmp(in, back, size)) {
    fprintf(stderr, "the decompressed data differs from the original\n");
    return 1;
  }

  free(in);  free(packed);  free(back);
  return 0;
}
