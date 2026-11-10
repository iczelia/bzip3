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

AC_DEFUN([BZ3_TOOL_PATHS], [
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_CANONICAL_HOST])
_BZ3_TOOL_PATH([LD], [ld])
_BZ3_TOOL_PATH([NM], [nm])
])

AC_DEFUN([_BZ3_TOOL_PATH], [
if test -z "$$1"; then
  bz3_tool=`$CC -print-prog-name=$2 2>/dev/null`
  case "$bz3_tool" in
    [[\\/]]* | ?:[[\\/]]*) ;;
    *)
      AC_PATH_TOOL([bz3_cv_path_$2], [$2], [])
      bz3_tool=$bz3_cv_path_$2 ;;
  esac
  case "$bz3_tool" in
    *[[\ \	]]*)
      AC_MSG_CHECKING([for a name of $2 without whitespace])
      bz3_short=`cygpath -m -s "$bz3_tool" 2>/dev/null`
      case "$bz3_short" in
        "" | *[[\ \	]]*)
          bz3_short=
          for bz3_name in "$host_alias-$2" "$2"; do
            test "x$bz3_name" = "x-$2" && continue
            bz3_found=`command -v "$bz3_name" 2>/dev/null`
            if test -n "$bz3_found"; then bz3_short=$bz3_name;  break; fi
          done ;;
      esac
      if test -z "$bz3_short"; then
        AC_MSG_RESULT([none])
        AC_MSG_ERROR([$bz3_tool is in a directory whose name contains whitespace, which libtool can not handle; add that directory to PATH or set $1 to a name free of whitespace])
      fi
      AC_MSG_RESULT([$bz3_short])
      $1=$bz3_short ;;
  esac
fi
])
