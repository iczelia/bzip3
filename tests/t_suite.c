/*  Copyright (C) 2022-2026 Kamila Szewczyk

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
    General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

#include "t_harness.h"

int main(int argc, char ** argv) {
  int i;
  xt_init();
  if (argc > 1) {
    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "unit")) xt_run_unit();
      else if (!strcmp(argv[i], "api")) xt_run_api();
      else if (!strcmp(argv[i], "files")) xt_run_files();
      else if (!strcmp(argv[i], "cli")) xt_run_cli();
      else {
        fprintf(stderr, "t_suite: unknown section %s\n", argv[i]);  return 2;
      }
    }
    return xt_finish("t_suite");
  }
  xt_run_unit();
  xt_run_api();
  xt_run_files();
  xt_run_cli();
  return xt_finish("t_suite");
}
