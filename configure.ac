dnl  Copyright (C) 2022-2026 Kamila Szewczyk
dnl
dnl  This program is free software: you can redistribute it and/or modify it
dnl  under the terms of the GNU Lesser General Public License as published by
dnl  the Free Software Foundation, either version 3 of the License, or (at your
dnl  option) any later version.
dnl
dnl  This program is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
dnl  General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with this program.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.70])
AC_INIT([bzip3], m4_esyscmd_s([build-aux/version.sh]), [https://github.com/iczelia/bzip3])
AC_CONFIG_SRCDIR([src/libbz3.c])
AC_CONFIG_HEADERS([config.h])
AC_COPYRIGHT([Copyright (C) 2022-2026 Kamila Szewczyk])
AC_LANG([C])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

AC_CANONICAL_HOST
AM_INIT_AUTOMAKE([-Wall foreign color-tests subdir-objects tar-ustar dist-bzip2 dist-xz dist-zip])
AM_SILENT_RULES([yes])

AC_PROG_CC
AC_USE_SYSTEM_EXTENSIONS
AM_PROG_AR
LT_INIT([win32-dll])

AC_DEFINE_UNQUOTED([BZ3_HOST_TRIPLE], ["$host"],
                   [Canonical host this build targets])

AC_SYS_LARGEFILE

AC_CHECK_HEADERS([unistd.h io.h sys/auxv.h sys/wait.h])
AC_CHECK_FUNCS([fsync getauxval])

dnl  The directory of bzip3.pc, without requiring pkg-config and its macros to build.
AC_ARG_WITH([pkgconfigdir],
  [AS_HELP_STRING([--with-pkgconfigdir=DIR], [install bzip3.pc in DIR (default: LIBDIR/pkgconfig)])],
  [pkgconfigdir=$withval], [pkgconfigdir='${libdir}/pkgconfig'])
AC_SUBST([pkgconfigdir])

case "$host_os" in
  mingw*) bz3_is_windows=yes ;;
  *)      bz3_is_windows=no  ;;
esac

dnl  Native Windows has threads of its own; elsewhere they are POSIX threads.
AC_ARG_ENABLE([threads],
  [AS_HELP_STRING([--disable-threads], [code the blocks of bz3_encode_blocks and -j one at a time])],
  [enable_threads=$enableval], [enable_threads=yes])
bz3_have_threads=no
if test "x$enable_threads" = "xyes"; then
  if test "x$bz3_is_windows" = "xyes"; then
    bz3_have_threads=yes
  else
    AC_CHECK_HEADER([pthread.h],
      [save_LIBS="$LIBS"
       AC_SEARCH_LIBS([pthread_create], [pthread],
         [bz3_have_threads=yes], [LIBS="$save_LIBS"])])
  fi
fi
if test "x$bz3_have_threads" = "xyes"; then
  AC_DEFINE([HAVE_THREADS], [1], [Define if thread support is available])
fi

dnl  The whole program in one file that depends on no shared library.
AC_ARG_ENABLE([static-exe],
  [AS_HELP_STRING([--enable-static-exe], [link bzip3 statically])],
  [enable_static_exe=$enableval], [enable_static_exe=no])
AM_CONDITIONAL([BZ3_STATIC_EXE], [test "x$enable_static_exe" = "xyes"])

case "$host_cpu" in
  i?86|x86_64*|amd64*) bz3_is_x86=yes ;;
  *)                   bz3_is_x86=no  ;;
esac
bz3_skip_simd=no
AC_ARG_ENABLE([simd],
  [AS_HELP_STRING([--disable-simd],[build the portable CRC only])],
  [test "x$enableval" = "xno" && bz3_skip_simd=yes])

BZ3_SSE42_CFLAGS=""
BZ3_ARMCRC_CFLAGS=""
if test "x$bz3_skip_simd" = "xno" && test "x$bz3_is_x86" = "xyes"; then
  BZ3_CHECK_ISA([SSE4.2 CRC32 intrinsics], [-msse4.2], [[
      #include <nmmintrin.h>
      #include <cpuid.h>
    ]], [[
      unsigned r[4];
      __get_cpuid(1, r, r + 1, r + 2, r + 3);
      return (int) (_mm_crc32_u32(_mm_crc32_u8(0, 1), r[2]));
    ]], [HAVE_SSE42], [Define if the SSE4.2 CRC-32C is built], [BZ3_SSE42_CFLAGS])
fi
if test "x$bz3_skip_simd" = "xno" && test "x$host_cpu" = "xaarch64"; then
  BZ3_CHECK_ISA([ARMv8 CRC32 intrinsics], [-march=armv8-a+crc], [[
      #include <arm_acle.h>
    ]], [[
      return (int) __crc32cd(__crc32cb(0, 1), 2);
    ]], [HAVE_ARM_CRC32], [Define if the ARMv8 CRC-32C is built], [BZ3_ARMCRC_CFLAGS])
fi
AC_SUBST([BZ3_SSE42_CFLAGS])
AM_CONDITIONAL([BZ3_SSE42], [test "x$BZ3_SSE42_CFLAGS" != "x"])
AC_SUBST([BZ3_ARMCRC_CFLAGS])
AM_CONDITIONAL([BZ3_ARMCRC], [test "x$BZ3_ARMCRC_CFLAGS" != "x"])

AX_APPEND_COMPILE_FLAGS([-Wall -Wextra -Wshadow -Wpointer-arith -Wcast-align \
  -Wwrite-strings \
  -Wmissing-prototypes -Wstrict-prototypes -Wredundant-decls \
  -Wno-unused-parameter -Wvla -Wformat=2 -Wundef], [BZ3_WARN_CFLAGS])
AC_SUBST([BZ3_WARN_CFLAGS])

AX_APPEND_COMPILE_FLAGS([-fno-delete-null-pointer-checks], [BZ3_HARDEN_CFLAGS])
AC_SUBST([BZ3_HARDEN_CFLAGS])

AC_ARG_ENABLE([sanitizers],
  [AS_HELP_STRING([--enable-sanitizers=LIST],
    [build with address,undefined (or yes), undefined or thread sanitizers])],
  [enable_sanitizers=$enableval], [enable_sanitizers=no])
case "$enable_sanitizers" in
  no) ;;
  yes|address|address,undefined|undefined,address)
    enable_sanitizers=address,undefined
    AX_APPEND_COMPILE_FLAGS([-fsanitize=address -fsanitize=undefined \
                             -fno-omit-frame-pointer], [BZ3_SAN_CFLAGS])
    AX_APPEND_LINK_FLAGS([-fsanitize=address -fsanitize=undefined],
                         [BZ3_SAN_LDFLAGS])
    ;;
  undefined)
    AX_APPEND_COMPILE_FLAGS([-fsanitize=undefined -fno-omit-frame-pointer],
                            [BZ3_SAN_CFLAGS])
    AX_APPEND_LINK_FLAGS([-fsanitize=undefined], [BZ3_SAN_LDFLAGS])
    ;;
  thread)
    AX_APPEND_COMPILE_FLAGS([-fsanitize=thread -fno-omit-frame-pointer],
                            [BZ3_SAN_CFLAGS])
    AX_APPEND_LINK_FLAGS([-fsanitize=thread], [BZ3_SAN_LDFLAGS])
    ;;
  *) AC_MSG_ERROR([unknown sanitizer selection: $enable_sanitizers]) ;;
esac
AC_SUBST([BZ3_SAN_CFLAGS])
AC_SUBST([BZ3_SAN_LDFLAGS])

dnl  Tests run uncapped when timeout(1) is unavailable.
AC_ARG_VAR([BZ3_TEST_TIMEOUT_SECONDS],
           [seconds the test suite may run (default 3600)])
test -n "$BZ3_TEST_TIMEOUT_SECONDS" || BZ3_TEST_TIMEOUT_SECONDS=3600
AC_CHECK_PROGS([BZ3_TIMEOUT_PROG], [timeout gtimeout], [])
if test -z "$BZ3_TIMEOUT_PROG"; then
  BZ3_TEST_TIMEOUT_SECONDS=
fi
AC_SUBST([BZ3_TIMEOUT_PROG])

dnl  The date of the source, for the manual pages.
AC_SUBST([BZ3_DATE], [m4_esyscmd_s([build-aux/version.sh --date])])

AC_CONFIG_FILES([Makefile bzip3.pc bzip3.1])
AC_CONFIG_FILES([scripts/bz3grep.1 scripts/bz3less.1 scripts/bz3more.1 scripts/bz3most.1])
AC_OUTPUT

AC_MSG_NOTICE([

  bzip3 $PACKAGE_VERSION configuration
    host          : $host
    threads       : $bz3_have_threads
    CRC-32C       : portable${BZ3_SSE42_CFLAGS:+, sse4.2}${BZ3_ARMCRC_CFLAGS:+, armv8-crc}
    static program: $enable_static_exe
    sanitizers    : $enable_sanitizers
    CFLAGS        : $CFLAGS
])
