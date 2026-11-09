ACLOCAL_AMFLAGS = -I build-aux

AM_CFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/src -I$(top_srcdir)/3rdparty/libsais

EXTRA_DIST = LICENSE NEWS PORTING.md README.md bootstrap.sh build-aux/git-version-gen \
	CMakeLists.txt bzip3-config.cmake.in \
	doc/bzip3_format.md doc/overview.md \
	3rdparty/libsais/LICENSE 3rdparty/libsais/libsais.c \
	$(BZ3_GOLDEN_FILES) tests/cli.sh tests/symbols.sh

pkgconfig_DATA = bzip3.pc

include_HEADERS = include/libbz3.h
noinst_HEADERS = include/common.h \
				 include/yarg.h \
				 src/sais.h \
				 3rdparty/libsais/libsais.h

lib_LTLIBRARIES = libbzip3.la
libbzip3_la_SOURCES = src/libbz3.c
libbzip3_la_LDFLAGS = -no-undefined -version-info 1:1:0

bin_PROGRAMS = bzip3
bzip3_CFLAGS = $(AM_CFLAGS)
bzip3_SOURCES = src/main.c
if ENABLE_STATIC
bzip3_SOURCES += $(libbzip3_la_SOURCES)
else
bzip3_LDADD = libbzip3.la
endif

dist_man_MANS = bzip3.1 bz3cat.1 bz3more.1 bz3less.1 bz3most.1 bz3grep.1 bunzip3.1

dist_bin_SCRIPTS = bz3cat bz3more bz3less bz3most bz3grep bunzip3

CLEANFILES = $(bin_PROGRAMS)

# End standard generic autotools stuff

# Begin special handling for autoconf VERSION being updated on commit

BUILT_SOURCES = .version
CLEANFILES += $(BUILT_SOURCES) .version-prev

src/bzip3-main.$(OBJEXT): .version

_BRANCH_REF != $(AWK) '{print ".git/" $$2}' .git/HEAD 2>/dev/null ||:

.version: $(_BRANCH_REF)
	@if [ -e "$(srcdir)/.tarball-version" ]; then \
		printf "$(VERSION)" > $@; \
	else \
		touch "$@-prev"; \
		if [ -e "$@" ]; then \
			cp "$@" "$@-prev"; \
		fi; \
		./build-aux/git-version-gen "$(srcdir)/.tarball-version" > $@; \
		cmp -s "$@" "$@-prev" || autoreconf configure.ac --force; \
	fi

dist-hook:
	printf "$(VERSION)" > "$(distdir)/.tarball-version"

# Files compressed by an older release; decoding them guards the backwards compatibility of the decoder.
BZ3_GOLDEN_FILES = tests/data/binary.bin tests/data/binary.bin.bz3 \
	tests/data/lzp.txt tests/data/lzp.txt.bz3 \
	tests/data/plain.txt tests/data/plain.txt.bz3 \
	tests/data/rle.bin tests/data/rle.bin.bz3 \
	tests/data/rle-lzp.bin tests/data/rle-lzp.bin.bz3 \
	tests/data/stored.bin tests/data/stored.bin.bz3

check_PROGRAMS = tests/api
tests_api_SOURCES = tests/api.c
tests_api_LDADD = libbzip3.la

TEST_EXTENSIONS = .sh
SH_LOG_COMPILER = $(SHELL)
AM_TESTS_ENVIRONMENT = \
	BZIP3='$(abs_builddir)/bzip3$(EXEEXT)'; export BZIP3; \
	BZIP3_TEST_DATA='$(abs_srcdir)/tests/data'; export BZIP3_TEST_DATA; \
	BZIP3_TEST_LIBS='$(abs_builddir)/.libs/libbzip3.a $(abs_builddir)/.libs/libbzip3.so $(abs_builddir)/.libs/libbzip3.dylib'; \
	export BZIP3_TEST_LIBS;
TESTS = tests/api$(EXEEXT) tests/cli.sh tests/symbols.sh

# Begin developer convenience targets

.PHONY: format
format: $(bzip3_SOURCES) $(libbzip3_la_SOURCES) include/libbz3.h include/common.h src/sais.h
	clang-format -i $^ tests/*.c examples/*.c

.PHONY: cloc
cloc: $(bzip3_SOURCES) $(libbzip3_la_SOURCES) $(include_HEADERS) $(noinst_HEADERS)
	cloc $^

CLEANFILES += LICENSE2 LICENSE.bz3
.PHONY: roundtrip

BZIP3 := bzip3$(EXEEXT)

roundtrip: $(BZIP3)
	rm -f $(builddir)/LICENSE2
	./$(BZIP3) -v -feb 6 $(srcdir)/LICENSE $(builddir)/LICENSE.bz3
	./$(BZIP3) -v -d $(builddir)/LICENSE.bz3 $(builddir)/LICENSE2
	cmp $(srcdir)/LICENSE $(builddir)/LICENSE2
	-command -v md5sum >/dev/null 2>&1 && md5sum $(builddir)/LICENSE.bz3

# Needs the examples directory, which is only a part of the git repository.
.PHONY: test
test: $(BZIP3)
	./$(BZIP3) -d < $(srcdir)/examples/shakespeare.txt.bz3 | cmp - $(srcdir)/examples/shakespeare.txt
