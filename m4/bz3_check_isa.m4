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

dnl  BZ3_CHECK_ISA(WHAT, FLAGS, INCLUDES, BODY, DEFINE, DESCRIPTION, VARIABLE)
dnl  Compile BODY after INCLUDES with FLAGS added to CFLAGS; on success define DEFINE and set VARIABLE to FLAGS.
AC_DEFUN([BZ3_CHECK_ISA],
[AC_MSG_CHECKING([for $1])
bz3_isa_save_CFLAGS=$CFLAGS
CFLAGS="$CFLAGS $2"
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([$3], [$4])],
  [AC_MSG_RESULT([yes])
   AC_DEFINE([$5], [1], [$6])
   $7="$2"],
  [AC_MSG_RESULT([no])])
CFLAGS=$bz3_isa_save_CFLAGS
])
