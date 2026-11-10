#!/bin/sh
#  Copyright (C) 2022-2026 Kamila Szewczyk
#
#  This program is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
#  General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

if test "x$1" = x--date; then
  cd "${2:-.}" || exit 1
  if test -e .git && d=`git log -1 --format=%cd --date=short 2>/dev/null` && test -n "$d"; then
    printf '%s' "$d"
  elif test -f .tarball-date; then
    tr -d '\n' < .tarball-date
  else
    date -u +%Y-%m-%d | tr -d '\n'
  fi
  exit 0
fi
cd "${1:-.}" || exit 1
if test -e .git && v=`git describe --tags --match '[0-9]*' --abbrev=7 --dirty 2>/dev/null`; then
  printf '%s' "$v" | sed 's/-\([0-9][0-9]*\)-g/.\1-/'
elif test -f .tarball-version; then
  tr -d '\n' < .tarball-version
else
  printf UNKNOWN
fi
