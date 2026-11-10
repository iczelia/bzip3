# BZip3

[![Build](https://github.com/iczelia/bzip3/actions/workflows/build.yml/badge.svg)](https://github.com/iczelia/bzip3/actions/workflows/build.yml)

A better, faster and stronger spiritual successor to BZip2. Features higher compression ratios and better performance thanks to a order-0 context mixing entropy coder, a fast Burrows-Wheeler transform code making use of suffix arrays and a Lempel Ziv+Prediction pass based on LZ77-style string matching and PPM-style context modeling.

Like its ancestor, **BZip3 excels at compressing text or code**.

## Installation

```console
# If using a git clone (not needed for source packages), first...
$ ./bootstrap.sh

# All...
$ ./configure
$ make
$ make check        # optional: run the test suite
$ sudo make install
```

Alternatively, you might be able to install bzip3 using your system's package manager:

[![Packaging status](https://repology.org/badge/vertical-allrepos/bzip3.svg?columns=3)](https://repology.org/project/bzip3/versions)

On macOS, you can use [Homebrew](https://brew.sh) to easily install:

```console
$ brew install bzip3
```

## Licensing

A breakdown of components and their licenses follows:

- (runtime) The codebase as a whole: Copyright 2022-2026, Kamila Szewczyk (k@iczelia.net); LGPL (LICENSE)
- (runtime) The Burrows-Wheeler transform (libsais, bundled unmodified in `3rdparty/libsais`) and LZP code: 2021-2026, Ilya Grebnov (ilya.grebnov@gmail.com); Apache 2.0 with LLVM exceptions (3rdparty/libsais/LICENSE)
- (compile-time) `build-aux`: Copyright 2011, Daniel Richard G (skunk@iSKUNK.ORG), 2019, Marc Stevens (marc.stevens@cwi.nl), 2008, Steven G. Johnson (stevenj@alum.mit.edu); GPL-3+ with AutoConf exception
- (compile-time) `build-aux/ax_check_compile_flag.m4`: Copyright 2008, Guido U. Draheim (guidod@gmx.de), 2011, Maarten Bosmans (mkbosmans@gmail.com); FSFAP
- (compile-time) `build-aux/git-version-gen`: Copyright 2007-2012, Free Software Foundation, Inc; GPLv3
- (runtime) `include/yarg.h` (argument parsing in the command line tool): Copyright 2026, Kamila Szewczyk (k@iczelia.net); 0BSD
- (runtime) `bz3grep`: Copyright 2003, Thomas Klausner; BSD-2-clause

`bzip3` as a whole is licensed under LGPLv3 only. It is not dual-licensed under LGPLv3 and Apache 2.0.

## Thanks

- Ilya Grebnov for his `libsais` library used for BWT construction in BZip3 and the LZP encoder which I had used as a reference implementation to improve myself.
- Caleb Maclennan for configuring autotools as a packaging-friendly build system for BZip3.
- Ilya Muravyov for his public domain BWT post-coder, a derivative of which is used in this project.
