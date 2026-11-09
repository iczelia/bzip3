/* Tests of the libbz3 API: round trips through both API levels and rejection of damaged input. */

#include <libbz3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KiB(x) ((x) * 1024)

static int failures = 0;

#define CHECK(cond)                                                                  \
    do {                                                                             \
        if (!(cond)) {                                                               \
            fprintf(stderr, "%s:%d: check failed: %s\n", __FILE__, __LINE__, #cond); \
            failures++;                                                              \
        }                                                                            \
    } while (0)

static uint32_t rng_state = 0x2545F491;

static uint32_t rng(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return rng_state;
}

/* Fill a buffer with data of the given kind: noise, text-like, long runs, long repeats or a mix. */
static void generate(uint8_t * buf, size_t size, int kind) {
    static const char * words[] = { "block ", "sorting ", "transform ", "entropy ", "context ", "mixing\n" };
    size_t i = 0;
    while (i < size) {
        int k = kind == 4 ? (int)(rng() % 4) : kind;
        size_t len = 1 + rng() % 600, j;
        if (len > size - i) len = size - i;
        switch (k) {
            case 0:
                for (j = 0; j < len; j++) buf[i + j] = (uint8_t)rng();
                break;
            case 1:
                for (j = 0; j < len;) {
                    const char * w = words[rng() % 6];
                    for (; *w && j < len; w++) buf[i + j++] = (uint8_t)*w;
                }
                break;
            case 2:
                memset(buf + i, (int)(rng() % 3 ? 0 : rng()), len);
                break;
            default:
                if (i < 200) {
                    for (j = 0; j < len; j++) buf[i + j] = (uint8_t)(rng() % 26 + 'a');
                } else {
                    size_t from = rng() % (i - 100);
                    for (j = 0; j < len; j++) buf[i + j] = buf[from + j];
                }
                break;
        }
        i += len;
    }
}

static void test_basics(void) {
    CHECK(bz3_version() != NULL);
    CHECK(bz3_bound(0) >= 13);
    CHECK(bz3_bound(1000) > 1000);
    CHECK(bz3_new(0) == NULL);
    CHECK(bz3_new(KiB(65) - 1) == NULL);
    CHECK(bz3_new(-1) == NULL);
    CHECK(bz3_new(512 * 1024 * 1024) == NULL);
    CHECK(bz3_min_memory_needed(KiB(64)) == 0);
    CHECK(bz3_min_memory_needed(KiB(65)) > (size_t)KiB(65) * 5);
    bz3_free(NULL);

    struct bz3_state * state = bz3_new(KiB(65));
    CHECK(state != NULL);
    if (state) {
        CHECK(bz3_last_error(state) == BZ3_OK);
        CHECK(bz3_strerror(state) != NULL);
        bz3_free(state);
    }
}

/* Compress and decompress `size' bytes with the frame API. */
static void test_frame(size_t size, int kind, uint32_t block_size) {
    uint8_t * in = malloc(size + 1);
    size_t cap = bz3_bound(size), out_size = cap, back_size = size;
    uint8_t * out = malloc(cap);
    uint8_t * back = malloc(size + 1);
    generate(in, size, kind);

    int err = bz3_compress(block_size, in, out, size, &out_size);
    CHECK(err == BZ3_OK);
    CHECK(out_size <= cap);
    if (err == BZ3_OK) {
        err = bz3_decompress(out, back, out_size, &back_size);
        CHECK(err == BZ3_OK);
        CHECK(back_size == size);
        CHECK(memcmp(in, back, size) == 0);

        /* An output buffer that is too small is reported, not overrun. */
        if (size > 0) {
            back_size = size - 1;
            CHECK(bz3_decompress(out, back, out_size, &back_size) == BZ3_ERR_DATA_TOO_BIG);
        }
    }

    /* So is one for the compressed data. */
    out_size = cap - 1;
    CHECK(bz3_compress(block_size, in, out, size, &out_size) == BZ3_ERR_DATA_TOO_BIG);

    free(in);
    free(out);
    free(back);
}

/* Damaged frames must be rejected, or at the very least decode into the original data. */
static void test_damage(void) {
    enum { SIZE = 3000 };
    uint8_t in[SIZE], back[SIZE];
    size_t cap = bz3_bound(SIZE), out_size = cap, back_size;
    uint8_t * out = malloc(cap);
    uint8_t * copy = malloc(cap);
    generate(in, SIZE, 4);
    CHECK(bz3_compress(KiB(65), in, out, SIZE, &out_size) == BZ3_OK);

    /* The 13 byte frame header carries no checksum: a different block size or block count there can
       still describe a well-formed frame. Everything after it is covered. */
    for (size_t i = 13; i < out_size; i++) {
        for (int bit = 0; bit < 8; bit += 3) {
            memcpy(copy, out, out_size);
            copy[i] ^= (uint8_t)(1 << bit);
            back_size = SIZE;
            int err = bz3_decompress(copy, back, out_size, &back_size);
            CHECK(err != BZ3_OK || (back_size == SIZE && memcmp(in, back, SIZE) == 0));
        }
    }

    /* Every proper prefix lacks data. */
    for (size_t len = 0; len < out_size; len++) {
        memcpy(copy, out, len);
        back_size = SIZE;
        CHECK(bz3_decompress(copy, back, len, &back_size) != BZ3_OK);
    }

    free(out);
    free(copy);
}

static void test_blocks(void) {
    enum { BLOCK = KiB(65) };
    struct bz3_state * state = bz3_new(BLOCK);
    size_t cap = bz3_bound(BLOCK);
    uint8_t * in = malloc(BLOCK);
    uint8_t * buf = malloc(cap);
    CHECK(state && in && buf);

    static const int sizes[] = { 0, 1, 63, 64, 65, 71, 72, 73, 1000, 30000, BLOCK - 1, BLOCK };
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
        for (int kind = 0; kind < 5; kind++) {
            int size = sizes[i];
            generate(in, (size_t)size, kind);
            memcpy(buf, in, (size_t)size);
            int32_t packed = bz3_encode_block(state, buf, size);
            CHECK(packed > 0 && (size_t)packed <= bz3_bound((size_t)size));
            CHECK(bz3_last_error(state) == BZ3_OK);
            if (packed <= 0) continue;

            /* 8 bytes stand for a stored empty block, which has no model byte to inspect. */
            if (packed >= 9) CHECK(bz3_orig_size_sufficient_for_decode(buf, (size_t)packed, size) >= 0);

            /* Too small a buffer is refused before anything is written past it. */
            if (packed > 9) CHECK(bz3_decode_block(state, buf, (size_t)packed - 1, packed, size) == -1);
            CHECK(bz3_decode_block(state, buf, cap, packed, size + 1) == -1);
            CHECK(bz3_last_error(state) != BZ3_OK);

            /* The failed attempts above may have clobbered the buffer. */
            memcpy(buf, in, (size_t)size);
            packed = bz3_encode_block(state, buf, size);
            CHECK(bz3_decode_block(state, buf, cap, packed, size) == size);
            CHECK(bz3_last_error(state) == BZ3_OK);
            CHECK(memcmp(buf, in, (size_t)size) == 0);
        }
    }

    CHECK(bz3_encode_block(state, buf, BLOCK + 1) == -1);
    CHECK(bz3_last_error(state) == BZ3_ERR_DATA_TOO_BIG);
    CHECK(bz3_encode_block(state, buf, -1) == -1);
    CHECK(bz3_decode_block(state, buf, cap, -1, 10) == -1);
    CHECK(bz3_decode_block(state, buf, cap, 100, -1) == -1);
    CHECK(bz3_decode_block(state, buf, 8, 8, 0) == -1);

    bz3_free(state);
    free(in);
    free(buf);
}

#ifdef PTHREAD
static void test_parallel(void) {
    enum { N = 5, BLOCK = KiB(65) };
    struct bz3_state * states[N];
    uint8_t *bufs[N], *orig[N];
    size_t caps[N];
    int32_t sizes[N], orig_sizes[N];
    for (int i = 0; i < N; i++) {
        states[i] = bz3_new(BLOCK);
        caps[i] = bz3_bound(BLOCK);
        bufs[i] = malloc(caps[i]);
        orig[i] = malloc(BLOCK);
        orig_sizes[i] = sizes[i] = BLOCK - i * 9000;
        generate(orig[i], (size_t)sizes[i], i);
        memcpy(bufs[i], orig[i], (size_t)sizes[i]);
    }
    bz3_encode_blocks(states, bufs, sizes, N);
    for (int i = 0; i < N; i++) CHECK(sizes[i] > 0 && bz3_last_error(states[i]) == BZ3_OK);
    bz3_decode_blocks(states, bufs, caps, sizes, orig_sizes, N);
    for (int i = 0; i < N; i++) {
        CHECK(bz3_last_error(states[i]) == BZ3_OK);
        CHECK(memcmp(bufs[i], orig[i], (size_t)orig_sizes[i]) == 0);
        bz3_free(states[i]);
        free(bufs[i]);
        free(orig[i]);
    }
}
#endif

int main(void) {
    test_basics();

    static const size_t sizes[] = {
        0, 1, 63, 64, 65, 100, 4096, KiB(65) - 1, KiB(65), KiB(65) + 1, 2 * KiB(65), 300000
    };
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++)
        for (int kind = 0; kind < 5; kind++) test_frame(sizes[i], kind, KiB(65));
    test_frame(300000, 4, 1024 * 1024);
    test_frame(300000, 1, 1);

    test_damage();
    test_blocks();
#ifdef PTHREAD
    test_parallel();
#endif

    if (failures) {
        fprintf(stderr, "%d check(s) failed\n", failures);
        return 1;
    }
    puts("all API tests passed");
    return 0;
}
