
/*
 * BZip3 - A spiritual successor to BZip2.
 * Copyright (C) 2022-2026 Kamila Szewczyk
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#if defined __MSVCRT__
    #include <io.h>
#endif

#ifndef O_BINARY
    #define O_BINARY 0
#endif

#include "common.h"
#include "libbz3.h"
#include "yarg.h"

#define MODE_DECODE 0
#define MODE_ENCODE 1
#define MODE_TEST 2
#define MODE_RECOVER 3

#define MIN_BLOCK_SIZE KiB(65)
#define MAX_BLOCK_SIZE MiB(511)
#define MAX_WORKERS 64

/* The regular file currently being written, if any. It is removed when the operation fails so that
   no truncated output is left behind. */
static const char * partial_output = NULL;

static void remove_partial_output(void) {
    if (partial_output) {
        remove(partial_output);
        partial_output = NULL;
    }
}

#if defined(__GNUC__) || defined(__clang__)
__attribute__((format(printf, 1, 2), noreturn))
#endif
static void fatal(const char * fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    remove_partial_output();
    exit(1);
}

static void version(void) {
    fprintf(stdout, "bzip3 " VERSION
                    "\n"
                    "Copyright (C) by Kamila Szewczyk, 2022-2026.\n"
                    "License: GNU Lesser GPL version 3 <https://www.gnu.org/licenses/lgpl-3.0.en.html>\n");
}

static void help(void) {
    fprintf(stdout,
            "bzip3 - better and stronger spiritual successor to bzip2.\n"
            "Usage: bzip3 [-e/-z/-d/-t/-c/-h/-V] [-b block_size] [-j jobs] files...\n"
            "Operations:\n"
            "  -e/-z, --encode   compress data (default)\n"
            "  -d, --decode      decompress data\n"
            "  -r, --recover     attempt at recovering corrupted data\n"
            "  -t, --test        verify validity of compressed data\n"
            "  -h, --help        display an usage overview\n"
            "  -f, --force       force overwriting output if it already exists\n"
            "      --rm          remove input files after successful (de)compression\n"
            "  -k, --keep        keep (don't delete) input files (default)\n"
            "  -v, --verbose     verbose mode (display more information)\n"
            "  -V, --version     display version information\n"
            "Extra flags:\n"
            "  -c, --stdout      force writing to standard output\n"
            "  -b N, --block=N   set block size in MiB {16}\n"
            "  -B, --batch       process all files specified as inputs\n"
#ifdef PTHREAD
            "  -j N, --jobs=N    set the amount of parallel threads {1}\n"
#endif
            "\n"
            "Report bugs to: https://github.com/iczelia/bzip3\n");
}

static void * xmalloc_zeroed(size_t count, size_t size) {
    void * ptr = calloc(count ? count : 1, size ? size : 1);
    if (!ptr) fatal("Failed to allocate memory.\n");
    return ptr;
}

static void xwrite(const void * data, size_t size, FILE * des) {
    if (size == 0) return;
    if (fwrite(data, 1, size, des) != size) fatal("Write error: %s\n", strerror(errno));
}

/* Read up to `size' bytes; fewer are returned only at the end of the file. */
static size_t xread(void * data, size_t size, FILE * des) {
    size_t count = fread(data, 1, size, des);
    if (ferror(des)) fatal("Read error: %s\n", strerror(errno));
    return count;
}

/* Read exactly `size' bytes. */
static void xread_exact(void * data, size_t size, FILE * des) {
    if (xread(data, size, des) != size) fatal("Error: Corrupt file\n");
}

static void close_out_file(FILE * des) {
    if (!des) return;

    if (fflush(des)) fatal("Error: Failed on fflush: %s\n", strerror(errno));

#ifdef __linux__
    while (fsync(fileno(des)) == -1) {
        if (errno == EINVAL) break; /* Not a file that can be synchronised, e.g. a pipe. */
        if (errno == EINTR) continue;
        fatal("Error: Failed on fsync: %s\n", strerror(errno));
    }
#endif

    if (des != stdout && fclose(des)) fatal("Error: Failed on fclose: %s\n", strerror(errno));

    partial_output = NULL;
}

/* The buffers and the coder states of all the workers. */
struct workers {
    int count;
    struct bz3_state ** states;
    u8 ** buffers;
    size_t * buffer_sizes;
    s32 *sizes, *orig_sizes;
};

static void workers_init(struct workers * w, int count, s32 block_size) {
    w->count = count;
    w->states = xmalloc_zeroed((size_t)count, sizeof(*w->states));
    w->buffers = xmalloc_zeroed((size_t)count, sizeof(*w->buffers));
    w->buffer_sizes = xmalloc_zeroed((size_t)count, sizeof(*w->buffer_sizes));
    w->sizes = xmalloc_zeroed((size_t)count, sizeof(*w->sizes));
    w->orig_sizes = xmalloc_zeroed((size_t)count, sizeof(*w->orig_sizes));
    for (int i = 0; i < count; i++) {
        w->states[i] = bz3_new(block_size);
        if (!w->states[i]) fatal("Failed to create a block encoder state.\n");
        w->buffer_sizes[i] = bz3_bound((size_t)block_size);
        /* Zeroed, so that the recovery mode never writes out uninitialised memory. */
        w->buffers[i] = xmalloc_zeroed(w->buffer_sizes[i], 1);
    }
}

static void workers_free(struct workers * w) {
    for (int i = 0; i < w->count; i++) {
        bz3_free(w->states[i]);
        free(w->buffers[i]);
    }
    free(w->states);
    free(w->buffers);
    free(w->buffer_sizes);
    free(w->sizes);
    free(w->orig_sizes);
}

/* Encode or decode the first `n' blocks held by the workers. */
static void workers_run(struct workers * w, int n, int encode) {
#ifdef PTHREAD
    if (n > 1) {
        if (encode)
            bz3_encode_blocks(w->states, w->buffers, w->sizes, n);
        else
            bz3_decode_blocks(w->states, w->buffers, w->buffer_sizes, w->sizes, w->orig_sizes, n);
        return;
    }
#endif
    for (int i = 0; i < n; i++) {
        if (encode)
            w->sizes[i] = bz3_encode_block(w->states[i], w->buffers[i], w->sizes[i]);
        else
            bz3_decode_block(w->states[i], w->buffers[i], w->buffer_sizes[i], w->sizes[i], w->orig_sizes[i]);
    }
}

/* Validate the compressed and original sizes stored in a block header. */
static int valid_block_header(s32 new_size, s32 old_size, s32 block_size) {
    return new_size >= 0 && (size_t)new_size <= bz3_bound((size_t)block_size) && old_size >= 0 &&
           old_size <= block_size;
}

static int process(FILE * input_des, FILE * output_des, int mode, s32 block_size, int workers, int verbose,
                   const char * file_name) {
    uint64_t bytes_read = 0, bytes_written = 0;
    u8 header[9];

    if ((mode == MODE_ENCODE && isatty(fileno(output_des))) || (mode != MODE_ENCODE && isatty(fileno(input_des)))) {
        fprintf(stderr, "Refusing to read/write binary data from/to the terminal.\n");
        return 1;
    }

    /* Reset errno after the isatty() call. */
    errno = 0;

    if (mode == MODE_ENCODE) {
        memcpy(header, "BZ3v1", 5);
        write_neutral_s32(header + 5, block_size);
        xwrite(header, 9, output_des);
        bytes_written += 9;
    } else {
        if (xread(header, 5, input_des) != 5 || memcmp(header, "BZ3v1", 5) != 0) {
            fprintf(stderr, "Invalid signature.\n");
            return 1;
        }

        xread_exact(header + 5, 4, input_des);
        block_size = read_neutral_s32(header + 5);

        if (block_size < MIN_BLOCK_SIZE || block_size > MAX_BLOCK_SIZE) {
            fprintf(stderr, "The input file is corrupted. Reason: Invalid block size in the header.\n");
            if (mode != MODE_RECOVER) return 1;
            fprintf(stderr, "Recovery mode: Proceeding.\n");
            block_size = MAX_BLOCK_SIZE;
        }

        bytes_read += 9;
    }

    if (workers < 1) workers = 1;
#ifndef PTHREAD
    workers = 1;
#endif

    struct workers w;
    workers_init(&w, workers, block_size);

    int status = 0;
    for (int eof = 0; !eof && !status;) {
        /* Fill as many buffers as there is data for. */
        int n = 0;
        for (; n < workers; n++) {
            if (mode == MODE_ENCODE) {
                size_t count = xread(w.buffers[n], (size_t)block_size, input_des);
                eof = count < (size_t)block_size;
                if (count == 0) break;
                w.sizes[n] = w.orig_sizes[n] = (s32)count;
                bytes_read += count;
                if (eof) {
                    n++;
                    break;
                }
            } else {
                size_t count = xread(header, 8, input_des);
                if (count == 0) {
                    eof = 1;
                    break;
                }
                if (count != 8) fatal("Error: Corrupt file\n");
                w.sizes[n] = read_neutral_s32(header);
                w.orig_sizes[n] = read_neutral_s32(header + 4);
                if (!valid_block_header(w.sizes[n], w.orig_sizes[n], block_size)) {
                    fprintf(stderr, "Failed to decode a block: Inconsistent headers.\n");
                    status = 1;
                    break;
                }
                xread_exact(w.buffers[n], (size_t)w.sizes[n], input_des);
                bytes_read += 8 + (uint64_t)w.sizes[n];
            }
        }
        if (status) break;

        workers_run(&w, n, mode == MODE_ENCODE);

        for (int i = 0; i < n && !status; i++) {
            if (bz3_last_error(w.states[i]) != BZ3_OK) {
                if (mode == MODE_ENCODE) {
                    fprintf(stderr, "Failed to encode a block: %s\n", bz3_strerror(w.states[i]));
                    status = 1;
                } else if (mode == MODE_RECOVER) {
                    fprintf(stderr, "Writing invalid block: %s\n", bz3_strerror(w.states[i]));
                } else {
                    fprintf(stderr, "Failed to decode a block: %s\n", bz3_strerror(w.states[i]));
                    status = 1;
                }
                if (status) break;
            }

            if (mode == MODE_ENCODE) {
                write_neutral_s32(header, w.sizes[i]);
                write_neutral_s32(header + 4, w.orig_sizes[i]);
                xwrite(header, 8, output_des);
                xwrite(w.buffers[i], (size_t)w.sizes[i], output_des);
                bytes_written += 8 + (uint64_t)w.sizes[i];
            } else {
                if (mode != MODE_TEST) xwrite(w.buffers[i], (size_t)w.orig_sizes[i], output_des);
                bytes_written += (uint64_t)w.orig_sizes[i];
            }
        }
    }

    workers_free(&w);

    if (status) return status;

    if (output_des && fflush(output_des)) fatal("Write error: %s\n", strerror(errno));

    if (verbose) {
        /* The size of the compressed data relative to the size of the original. */
        const double packed = (double)(mode == MODE_ENCODE ? bytes_written : bytes_read);
        const double unpacked = (double)(mode == MODE_ENCODE ? bytes_read : bytes_written);
        if (file_name) fprintf(stderr, " %s:", file_name);
        fprintf(stderr, "\t%s%" PRIu64 " -> %" PRIu64 " bytes, %.2f%%, %.2f bpb\n",
                mode == MODE_ENCODE || mode == MODE_DECODE ? "" : "OK, ", bytes_read, bytes_written,
                packed * 100.0 / unpacked, packed * 8.0 / unpacked);
    }

    return 0;
}

static int is_dir(const char * path) {
    struct stat sb;
    return stat(path, &sb) == 0 && S_ISDIR(sb.st_mode);
}

/* Parse a decimal number between `min' and `max'. Returns 0 on failure. */
static int parse_number(const char * str, long min, long max, long * result) {
    if (!*str) return 0;
    for (const char * p = str; *p; p++)
        if (!isdigit((unsigned char)*p)) return 0;
    errno = 0;
    long value = strtol(str, NULL, 10);
    if (errno || value < min || value > max) return 0;
    *result = value;
    return 1;
}

static FILE * open_input(const char * input) {
    if (input == NULL) return stdin;

    if (is_dir(input)) fatal("Error: input `%s' is a directory.\n", input);

    FILE * input_des = fopen(input, "rb");
    if (input_des == NULL) fatal("Error: failed to open input file `%s': %s\n", input, strerror(errno));
    return input_des;
}

/* Open the output file. Unless `force' is set, the file is created exclusively, so that an existing file
   (or a symbolic link planted in its place) is never written to. The new file gets the access permissions
   of the input, if known. */
static FILE * open_output(const char * output, int force, FILE * input_des) {
    if (output == NULL) return stdout;

    if (is_dir(output)) fatal("Error: output file `%s' is a directory.\n", output);

    struct stat in_stat, out_stat;
    int have_in_stat = input_des != stdin && fstat(fileno(input_des), &in_stat) == 0;
    int mode = have_in_stat && S_ISREG(in_stat.st_mode) ? (int)(in_stat.st_mode & 0777) : 0666;

    if (have_in_stat && in_stat.st_ino != 0 && stat(output, &out_stat) == 0 && in_stat.st_dev == out_stat.st_dev &&
        in_stat.st_ino == out_stat.st_ino)
        fatal("Error: input and output are the same file `%s'.\n", output);

    int fd = open(output, O_WRONLY | O_CREAT | O_BINARY | (force ? O_TRUNC : O_EXCL), mode);
    if (fd == -1 && errno == EEXIST && !force)
        fatal("Error: output file `%s' already exists. Use -f to force overwrite.\n", output);

    FILE * output_des = fd == -1 ? NULL : fdopen(fd, "wb");
    if (output_des == NULL) fatal("Error: failed to open output file `%s': %s\n", output, strerror(errno));

    /* Only regular files are cleaned up after a failure; never devices and the like. */
    if (fstat(fd, &out_stat) == 0 && S_ISREG(out_stat.st_mode)) partial_output = output;

    return output_des;
}

/* Derive the name of the output file from the name of the input. Returns NULL if a file to decompress
   does not carry the .bz3 extension. */
static char * output_name_for(const char * input, int mode) {
    size_t len = strlen(input);
    char * output = xmalloc_zeroed(len + 5, 1);
    memcpy(output, input, len);
    if (mode == MODE_ENCODE) {
        memcpy(output + len, ".bz3", 5);
    } else if (len > 4 && !strcmp(input + len - 4, ".bz3")) {
        output[len - 4] = 0;
    } else {
        fprintf(stderr, "Warning: file %s has an unknown extension, skipping.\n", input);
        free(output);
        return NULL;
    }
    return output;
}

/* Process a single file. A NULL input or output stands for the standard streams, and so does "-". */
static int process_file(const char * input, const char * output, int mode, s32 block_size, int workers, int verbose,
                        int force, int remove_input_file) {
    if (input && !strcmp(input, "-")) input = NULL;
    if (output && !strcmp(output, "-")) output = NULL;

    FILE * input_des = open_input(input);
    FILE * output_des = mode != MODE_TEST ? open_output(output, force, input_des) : NULL;

    int status = process(input_des, output_des, mode, block_size, workers, verbose, input);

    if (input_des != stdin) fclose(input_des);

    if (status) {
        if (output_des && output_des != stdout) fclose(output_des);
        remove_partial_output();
        return status;
    }

    close_out_file(output_des);

    /* The input is removed only once the output is known to be complete and written to a file. */
    if (remove_input_file && input && output && mode != MODE_TEST && remove(input))
        fatal("Error: failed to remove input file `%s': %s\n", input, strerror(errno));

    return 0;
}

enum { RM_OPTION = CHAR_MAX + 1 };

/* Carry out what the parsed command line asks for. Returns the exit status. */
static int run(const yarg_result * res) {
    int mode = MODE_ENCODE;
    int force = 0, force_stdstreams = 0, workers = 0, batch = 0, verbose = 0, remove_input_file = 0;
    s32 block_size = MiB(16);
    long number;

    for (int i = 0; i < res->argc; i++) {
        switch (res->args[i].opt) {
            case 'e':
            case 'z':
                mode = MODE_ENCODE;
                break;
            case 'd':
                mode = MODE_DECODE;
                break;
            case 'r':
                mode = MODE_RECOVER;
                break;
            case 't':
                mode = MODE_TEST;
                break;
            case 'c':
                force_stdstreams = 1;
                break;
            case 'f':
                force = 1;
                break;
            case RM_OPTION:
                remove_input_file = 1;
                break;
            case 'k':
                remove_input_file = 0;
                break;
            case 'h':
                help();
                return 0;
            case 'V':
                version();
                return 0;
            case 'B':
                batch = 1;
                break;
            case 'v':
                verbose = 1;
                break;
            case 'b':
                if (!parse_number(res->args[i].arg, 0, LONG_MAX, &number)) {
                    fprintf(stderr, "bzip3: invalid block size: %s\n", res->args[i].arg);
                    return 1;
                }
                if (number < 1 || number > MAX_BLOCK_SIZE / MiB(1)) {
                    fprintf(stderr, "Block size must be between 65 KiB and 511 MiB.\n");
                    return 1;
                }
                block_size = MiB((s32)number);
                break;
#ifdef PTHREAD
            case 'j':
                if (!parse_number(res->args[i].arg, 0, LONG_MAX, &number)) {
                    fprintf(stderr, "bzip3: invalid amount of jobs: %s\n", res->args[i].arg);
                    return 1;
                }
                if (number > MAX_WORKERS) {
                    fprintf(stderr, "Number of workers must be between 0 and %d.\n", MAX_WORKERS);
                    return 1;
                }
                workers = (int)number;
                break;
#endif
        }
    }

#if defined(__MSVCRT__)
    setmode(STDIN_FILENO, O_BINARY);
    setmode(STDOUT_FILENO, O_BINARY);
#endif

    int status = 0;

    if (batch && res->pos_argc) {
        /* Process each of the files on its own. */
        for (int i = 0; i < res->pos_argc; i++) {
            char * input = res->pos_args[i];
            char * output = NULL;
            if (mode != MODE_TEST && !force_stdstreams && strcmp(input, "-") != 0) {
                output = output_name_for(input, mode);
                if (!output) {
                    status = 1;
                    continue;
                }
            }
            if (process_file(input, output, mode, block_size, workers, verbose, force, remove_input_file)) status = 1;
            free(output);
        }
    } else {
        if (res->pos_argc > 2) {
            fprintf(stderr, "Error: too many files specified.\n");
            return 1;
        }

        /* With one file name, the output name is derived from it; with two, the second one is the output. */
        char * input = res->pos_argc > 0 ? res->pos_args[0] : NULL;
        char * output = res->pos_argc > 1 ? res->pos_args[1] : NULL;
        char * derived = NULL;
        if (input && strcmp(input, "-") != 0 && !output && mode != MODE_TEST && !force_stdstreams) {
            output = derived = output_name_for(input, mode);
            if (!output) return 1;
        }

        status = process_file(input, output, mode, block_size, workers, verbose, force, remove_input_file);
        free(derived);
    }

    if (fclose(stdout)) {
        fprintf(stderr, "Error: Failed on fclose(stdout): %s\n", strerror(errno));
        return 1;
    }

    return status;
}

int main(int argc, char * argv[]) {
    yarg_options opt[] = { { 'e', no_argument, "encode" },
                           { 'z', no_argument, "encode" }, /* alias */
                           { 'd', no_argument, "decode" },
                           { 't', no_argument, "test" },
                           { 'c', no_argument, "stdout" },
                           { 'f', no_argument, "force" },
                           { 'r', no_argument, "recover" },
                           { 'h', no_argument, "help" },
                           { RM_OPTION, no_argument, "rm" },
                           { 'k', no_argument, "keep" },
                           { 'V', no_argument, "version" },
                           { 'v', no_argument, "verbose" },
                           { 'b', required_argument, "block" },
                           { 'B', no_argument, "batch" },
#ifdef PTHREAD
                           { 'j', required_argument, "jobs" },
#endif
                           { 0, no_argument, NULL } };
    yarg_settings settings = {
        .dash_dash = true,
        .style = YARG_STYLE_UNIX,
    };
    yarg_result * res = yarg_parse(argc, argv, opt, settings);
    if (!res) {
        fprintf(stderr, "bzip3: out of memory.\n");
        return 1;
    }

    int status = 1;
    if (res->error) {
        fputs(res->error, stderr);
        fputs("Try 'bzip3 --help' for more information.\n", stderr);
    } else {
        status = run(res);
    }

    yarg_destroy(res);
    return status;
}
