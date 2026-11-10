# bzip3

[![Build](https://github.com/iczelia/bzip3/actions/workflows/build.yml/badge.svg)](https://github.com/iczelia/bzip3/actions/workflows/build.yml)

A better, faster and stronger spiritual successor to bzip2. It compresses
better and runs faster thanks to an order-0 context mixing entropy coder, a
fast Burrows-Wheeler transform built on suffix arrays, and a Lempel-Ziv
prediction pass based on LZ77-style string matching and PPM-style context
modeling.

bzip3 is licensed under the GNU Lesser GPL version 3 or later. See
[COPYING](COPYING) and [AUTHORS](AUTHORS). Report issues at
<https://github.com/iczelia/bzip3>.

## Quick start

```sh
bzip3 file              # writes file.bz3
bzip3 -d file.bz3       # restores file
bzip3 -B -j 4 *.txt     # process every file, four blocks at a time
tar c dir | bzip3 -b 64 > dir.tar.bz3
bz3cat dir.tar.bz3 | tar x
```

Encoding is the default setting. `-d` decompresses and `-t` checks a file. With one
file name, the output name is derived by adding or removing `.bz3`; a second
name is the output, and `-c` writes to the standard output. `-B` treats every
name as an input. Existing files are kept unless `-f` is given, and `--rm`
removes the inputs once they were processed.

`-b` sets the block size in MiB, from 1 to 511; the default is 16. Larger
blocks compress better and need about six times their size in memory, to
compress and to decompress. `-j` codes that many blocks at once.

The `bunzip3`, `bz3cat`, `bz3less`, `bz3more`, `bz3most` and `bz3grep` scripts
work like their bzip2 namesakes. See `bzip3(1)`.

## Installation

Use your package manager or download a release from GitHub Releases.

[![Packaging status](https://repology.org/badge/vertical-allrepos/bzip3.svg?columns=3)](https://repology.org/project/bzip3/versions)

To build a release tarball, run

```sh
./configure
make
make check
sudo make install
```

The project uses C99 and, for `-j`, POSIX or Windows threads. libsais is
bundled unmodified in `3rdparty/libsais`.

| Configure option | Effect |
| --- | --- |
| `--enable-sanitizers` | Enable ASan and UBSan for tests |
| `--disable-simd` | Build the portable CRC-32C only |
| `--disable-threads` | Code the blocks of `-j` and `bz3_encode_blocks` one at a time |
| `--enable-static-exe` | Link `bzip3` statically |
| `--disable-shared`, `--disable-static` | Choose the kinds of `libbzip3` to build |

Run `./bootstrap` first when building from a Git checkout. It requires
autoconf, automake and libtool. The version comes from `git describe`, so a
checkout needs its tags.

The build is portable by default: pass `CFLAGS="-O2 -march=native"` to
`configure` to tune it for the machine it is built on.

## Library

`libbzip3` and its header `libbz3.h` are installed with the tool, with a
pkg-config file:

```sh
cc example.c $(pkg-config --cflags --libs bzip3)
```

`bz3_compress` and `bz3_decompress` code a buffer in one call;
`bz3_encode_block` and `bz3_decode_block` code single blocks and are what the
tool is built on. `include/libbz3.h` documents them and `examples/hl-api.c`
uses them. The format is described in `doc/bzip3_format.md`.

Every release reads the files of the releases before it. The shared library
exports only the functions of `libbz3.h`, and the static library defines no
symbol outside of the `bz3_` namespace.

## Performance and portability

Files are portable across hosts and byte orders. The checksum is CRC-32C; on
x86 and 64-bit ARM a kernel using the CPU's CRC instructions is built next to
the portable one and chosen at run time when the CPU has them.

Windows builds use MinGW as follows:

```sh
./configure --host=x86_64-w64-mingw32 CC=x86_64-w64-mingw32-gcc --enable-static-exe
make
```

## Testing and fuzzing

`make check` runs unit tests of the coders, tests of the library interface and
of the command line tool, and decodes reference files written by bzip3 1.5.5.
`make check-full` and `make check-torture` code more data.

`fuzz/build.sh` builds libFuzzer harnesses for the block decoder, the frame
decoder and the round trip with ASan and UBSan; `FZ_SAN=msan` and
`FZ_SAN=standalone` select MemorySanitizer or a plain tool that replays
inputs.
