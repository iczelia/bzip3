/* An example of the high level API: compress a file in memory and decompress it again.
 *
 * Build and run in this directory, with bzip3 installed:
 *
 *   cc hl-api.c -o hl-api $(pkg-config --cflags --libs bzip3) && ./hl-api
 */

#include <libbz3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MB (1024 * 1024)

int main(void) {
    printf("Compressing shakespeare.txt back and forth in memory.\n");

    /* Read the entire "shakespeare.txt" file to memory. */
    FILE * fp = fopen("shakespeare.txt", "rb");
    if (!fp) {
        perror("shakespeare.txt");
        return 1;
    }
    fseek(fp, 0, SEEK_END);
    size_t size = (size_t)ftell(fp);
    fseek(fp, 0, SEEK_SET);
    uint8_t * buffer = malloc(size);
    if (!buffer || fread(buffer, 1, size, fp) != size) {
        fprintf(stderr, "Failed to read shakespeare.txt\n");
        return 1;
    }
    fclose(fp);

    /* Compress the file. The output buffer must hold at least bz3_bound(size) bytes. */
    size_t out_size = bz3_bound(size);
    uint8_t * outbuf = malloc(out_size);
    if (!outbuf) return 1;
    int bzerr = bz3_compress(1 * MB, buffer, outbuf, size, &out_size);
    if (bzerr != BZ3_OK) {
        printf("bz3_compress() failed with error code %d\n", bzerr);
        return 1;
    }

    printf("%zu => %zu\n", size, out_size);

    /* Decompress the file into a second buffer and compare it with the original. */
    size_t back_size = size;
    uint8_t * back = malloc(back_size);
    if (!back) return 1;
    bzerr = bz3_decompress(outbuf, back, out_size, &back_size);
    if (bzerr != BZ3_OK) {
        printf("bz3_decompress() failed with error code %d\n", bzerr);
        return 1;
    }

    printf("%zu => %zu\n", out_size, back_size);

    if (back_size != size || memcmp(buffer, back, size) != 0) {
        printf("The decompressed data differs from the original!\n");
        return 1;
    }

    free(buffer);
    free(outbuf);
    free(back);
    return 0;
}
