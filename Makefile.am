##  Copyright (C) 2022-2026 Kamila Szewczyk
##
##  This program is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or (at your
##  option) any later version.
##
##  This program is distributed in the hope that it will be useful, but
##  WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
##  General Public License for more details.
##
##  You should have received a copy of the GNU Lesser General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.

ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -std=c99 $(BZ3_WARN_CFLAGS) $(BZ3_HARDEN_CFLAGS) $(BZ3_SAN_CFLAGS) $(BZ3_WIN_CFLAGS)
AM_LDFLAGS = $(BZ3_SAN_LDFLAGS) $(BZ3_WIN_LDFLAGS)
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir)/include -I$(top_srcdir)/src $(BZ3_WIN_CPPFLAGS)

lib_LTLIBRARIES = libbzip3.la
include_HEADERS = include/libbz3.h
pkgconfig_DATA = bzip3.pc
bin_PROGRAMS = bzip3
dist_bin_SCRIPTS = scripts/bz3cat scripts/bz3more scripts/bz3less scripts/bz3most scripts/bz3grep scripts/bunzip3
man1_MANS = bzip3.1 scripts/bz3grep.1 scripts/bz3less.1 scripts/bz3more.1 scripts/bz3most.1
dist_man1_MANS = scripts/bz3cat.1 scripts/bunzip3.1

##  Sources by directory.
BZ3_SYS = src/sys/cpu.c src/sys/crc.c
BZ3_CODER = src/coder/cm.c src/coder/lzp.c src/coder/rle.c
BZ3_SOURCES = src/libbz3.c $(BZ3_SYS) $(BZ3_CODER)

noinst_HEADERS = \
  src/sys/common.h src/sys/cpu.h src/sys/crc.h src/sys/crchw.h src/sys/map.h src/sys/mem.h \
  src/sys/meta.h src/sys/thread.h src/sys/win32.h src/sys/yarg.h \
  src/coder/bwt.h src/coder/cm.h src/coder/lzp.h src/coder/rle.h \
  3rdparty/libsais/libsais.h

##  The library's objects, shared by the installed library and the tests, which reach past its interface.
noinst_LTLIBRARIES = libbz3core.la
libbz3core_la_SOURCES = $(BZ3_SOURCES)
libbz3core_la_LIBADD = libbz3sais.la

##  The bundled libsais is compiled as a part of bwt.c, without the warnings this project holds its own code to.
noinst_LTLIBRARIES += libbz3sais.la
libbz3sais_la_SOURCES  = src/coder/bwt.c
libbz3sais_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/3rdparty/libsais
libbz3sais_la_CFLAGS   = -std=c99 $(BZ3_HARDEN_CFLAGS) $(BZ3_SAN_CFLAGS) $(BZ3_WIN_CFLAGS)

##  The kernels built for instruction sets beyond the baseline, which run when the CPU has them.
if BZ3_SSE42
noinst_LTLIBRARIES += libbz3_sse42.la
libbz3_sse42_la_SOURCES = src/sys/crcx86.c
libbz3_sse42_la_CFLAGS  = $(AM_CFLAGS) $(BZ3_SSE42_CFLAGS)
libbz3core_la_LIBADD += libbz3_sse42.la
endif
if BZ3_ARMCRC
noinst_LTLIBRARIES += libbz3_armcrc.la
libbz3_armcrc_la_SOURCES = src/sys/crcarm.c
libbz3_armcrc_la_CFLAGS  = $(AM_CFLAGS) $(BZ3_ARMCRC_CFLAGS)
libbz3core_la_LIBADD += libbz3_armcrc.la
endif

##  Only the functions of include/libbz3.h leave the shared library.
libbzip3_la_SOURCES =
libbzip3_la_LIBADD = libbz3core.la
libbzip3_la_LDFLAGS = $(AM_LDFLAGS) -no-undefined -version-info 1:1:0 \
                      -export-symbols $(top_srcdir)/src/libbz3.sym
EXTRA_libbzip3_la_DEPENDENCIES = $(top_srcdir)/src/libbz3.sym

bzip3_SOURCES = src/main.c src/sys/map.c src/sys/mem.c src/sys/meta.c src/sys/thread.c src/sys/win32.c
bzip3_LDADD = libbzip3.la
bzip3_LDFLAGS = $(AM_LDFLAGS) $(BZ3_WIN_EXE_LDFLAGS)
if BZ3_STATIC_EXE
bzip3_LDFLAGS += -all-static
endif

##  Test suite linked to the library's objects and run against the built binary.
check_PROGRAMS = tests/t_suite
TESTS = tests/t_suite tests/symbols.sh
TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)

tests_t_suite_SOURCES = tests/t_suite.c tests/t_harness.c tests/t_unit.c \
                        tests/t_api.c tests/t_files.c tests/t_cli.c
tests_t_suite_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/tests
tests_t_suite_LDADD = libbz3core.la

test-programs: $(check_PROGRAMS)

##  Override with `make check BZ3_TEST_TIMEOUT_SECONDS=60`.
LOG_COMPILER = $(BZ3_TIMEOUT_PROG) $(BZ3_TEST_TIMEOUT_SECONDS)

##  Convert MSYS2 paths for the native test binary.
AM_TESTS_ENVIRONMENT = \
  BZ3=`cygpath -m '$(BZ3_TEST_BINARY)' 2>/dev/null || echo '$(BZ3_TEST_BINARY)'`; \
                                                      export BZ3; \
  BZ3_TEST_DATA=`cygpath -m '$(abs_top_srcdir)/tests/data' 2>/dev/null \
                 || echo '$(abs_top_srcdir)/tests/data'`; export BZ3_TEST_DATA; \
  BZ3_TEST_LEVEL="$${BZ3_TEST_LEVEL-quick}";          export BZ3_TEST_LEVEL; \
  BZ3_TEST_TRACE="$${BZ3_TEST_TRACE-0}";              export BZ3_TEST_TRACE; \
  BZ3_TEST_LIBS='$(abs_top_builddir)/.libs';          export BZ3_TEST_LIBS; \
  BZ3_TEST_SYMS='$(abs_top_srcdir)/src/libbz3.sym';   export BZ3_TEST_SYMS; \
  NM='$(NM)';                                         export NM;

BZ3_TEST_BINARY = $(abs_top_builddir)/bzip3$(EXEEXT)

##  Files compressed by bzip3 1.5.5; decoding them guards the backwards compatibility of the decoder.
BZ3_TEST_FILES = \
  tests/data/binary.bin tests/data/binary.bin.bz3 tests/data/lzp.txt tests/data/lzp.txt.bz3 \
  tests/data/plain.txt tests/data/plain.txt.bz3 tests/data/rle.bin tests/data/rle.bin.bz3 \
  tests/data/rle-lzp.bin tests/data/rle-lzp.bin.bz3 tests/data/stored.bin tests/data/stored.bin.bz3

EXTRA_DIST = \
  bootstrap build-aux/version.sh README.md src/libbz3.sym \
  3rdparty/libsais/libsais.c 3rdparty/libsais/LICENSE \
  doc/bzip3_format.md doc/overview.md examples/hl-api.c \
  tests/t_harness.h tests/symbols.sh $(BZ3_TEST_FILES) \
  fuzz/harness.h fuzz/fuzz_block.c fuzz/fuzz_frame.c fuzz/fuzz_rt.c fuzz/fzmain.c \
  fuzz/build.sh fuzz/seeds/63_byte_file.bin fuzz/seeds/65_byte_file.bin fuzz/seeds/readme.txt \
  contrib/gitlog-to-changelog

CLEANFILES = t_suite-*.tmp t_suite-*.tmp.bz3

.PHONY: check-full check-torture print-distfiles test-programs

##  Full tests code more data.
check-full:
	$(MAKE) check BZ3_TEST_LEVEL=full

check-torture:
	$(MAKE) check BZ3_TEST_LEVEL=torture

##  Print the files included in the distribution tarball.
print-distfiles:
	@for f in $(DISTFILES) $(DIST_COMMON) $(am__DIST_COMMON) \
	          $(ACLOCAL_M4); do \
	  case $$f in \
	    $(srcdir)/*)     echo "$${f#$(srcdir)/}" ;; \
	    $(top_srcdir)/*) echo "$${f#$(top_srcdir)/}" ;; \
	    *)               echo "$$f" ;; \
	  esac; \
	done

##  Record the version and date, and write the ChangeLog from the git history or keep the one of an unpacked package.
dist-hook:
	@echo '$(VERSION)' > "$(distdir)/.tarball-version"
	@echo '$(BZ3_DATE)' > "$(distdir)/.tarball-date"
	@if test -e "$(top_srcdir)/.git"; then \
	  $(top_srcdir)/contrib/gitlog-to-changelog --srcdir="$(top_srcdir)" \
	    --format='%s%n%n%b%n' > "$(distdir)/ChangeLog.tmp" && \
	  mv "$(distdir)/ChangeLog.tmp" "$(distdir)/ChangeLog"; \
	elif test -f "$(top_srcdir)/ChangeLog"; then \
	  cp "$(top_srcdir)/ChangeLog" "$(distdir)/ChangeLog"; \
	else \
	  echo "dist-hook: no git history and no ChangeLog" >&2; exit 1; \
	fi
