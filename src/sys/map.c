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

#include "sys/map.h"

#include <errno.h>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <io.h>
#else
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_MMAP
#include <signal.h>
#include <sys/mman.h>
#endif
#endif

#define WINDOW ((uint64_t) (sizeof (void *) > 4 ? MiB(256) : MiB(32)))

#if defined(_WIN32) || defined(HAVE_MMAP)
#define HAVE_MAP 1
#endif

#ifdef HAVE_MAP

static void unmap(bz3_in * in) {
  if (!in->win) return;
#ifdef _WIN32
  UnmapViewOfFile(in->win);
#else
  munmap(in->win, in->win_len);
#endif
  in->win = NULL;
}

static int slide(bz3_in * in) {
  const uint64_t off = in->pos & ~(WINDOW - 1);
  const sz len = (sz) MIN(WINDOW, in->size - off);
  unmap(in);
#ifdef _WIN32
  in->win = MapViewOfFile((HANDLE) in->handle, FILE_MAP_READ, (DWORD) (off >> 32), (DWORD) off, len);
#else
  in->win = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fileno(in->f), (off_t) off);
  if (in->win == MAP_FAILED) in->win = NULL;
#if defined(HAVE_POSIX_MADVISE) && defined(POSIX_MADV_SEQUENTIAL)
  if (in->win) posix_madvise(in->win, len, POSIX_MADV_SEQUENTIAL);
#endif
#endif
  in->win_off = off;  in->win_len = len;
  return in->win != NULL;
}

static int unmapped(bz3_in * in) {
  unmap(in);
  in->mapped = 0;
#ifdef _WIN32
  if (in->handle) CloseHandle((HANDLE) in->handle);
  in->handle = NULL;
  return !_fseeki64(in->f, (__int64) in->pos, SEEK_SET);
#else
  return !fseeko(in->f, (off_t) in->pos, SEEK_SET);
#endif
}

#endif

void bz3_in_open(bz3_in * in, FILE * f) {
  memset(in, 0, sizeof *in);
  in->f = f;
#if defined(_WIN32)
  { HANDLE h = (HANDLE) _get_osfhandle(_fileno(f));
    DWORD lo, hi = 0;
    if (h == INVALID_HANDLE_VALUE || GetFileType(h) != FILE_TYPE_DISK) return;
    SetLastError(NO_ERROR);
    lo = GetFileSize(h, &hi);
    if (lo == INVALID_FILE_SIZE && GetLastError() != NO_ERROR) return;
    in->size = (uint64_t) hi << 32 | lo;
    if (!in->size) return;
    hi = 0;  SetLastError(NO_ERROR);
    lo = SetFilePointer(h, 0, (LONG *) &hi, FILE_CURRENT);
    if (lo == INVALID_SET_FILE_POINTER && GetLastError() != NO_ERROR) return;
    in->pos = (uint64_t) hi << 32 | lo;
    if (in->pos >= in->size) { in->pos = 0;  return; }
    in->handle = CreateFileMappingA(h, NULL, PAGE_READONLY, 0, 0, NULL);
    if (!in->handle) return;
    in->mapped = 1; }
#elif defined(HAVE_MMAP)
  { struct stat st;
    off_t at;
    if (fstat(fileno(f), &st) || !S_ISREG(st.st_mode) || st.st_size <= 0) return;
    at = lseek(fileno(f), 0, SEEK_CUR);
    if (at < 0 || at >= st.st_size) return;
    in->size = (uint64_t) st.st_size;  in->pos = (uint64_t) at;
    in->mapped = 1; }
#endif
#ifdef HAVE_MAP
  if (in->mapped && !slide(in)) {
    if (!unmapped(in)) in->failed = 1;
  }
#endif
}

sz bz3_in_read(bz3_in * in, void * data, sz n) {
  sz done = 0, got;
#ifdef HAVE_MAP
  while (in->mapped && done < n) {
    sz k;
    if (in->pos == in->size) {
      if (!unmapped(in)) { in->failed = 1;  return done; }
      break;
    }
    if (in->pos - in->win_off >= in->win_len && !slide(in)) {
      if (!unmapped(in)) { in->failed = 1;  return done; }
      break;
    }
    k = (sz) MIN((uint64_t) (n - done), in->win_off + in->win_len - in->pos);
    memcpy((u8 *) data + done, in->win + (in->pos - in->win_off), k);
    in->pos += k;  done += k;
  }
  if (done == n) return done;
#endif
  got = fread((u8 *) data + done, 1, n - done, in->f);
  if (ferror(in->f)) in->failed = 1;
  return done + got;
}

void bz3_in_close(bz3_in * in) {
#ifdef HAVE_MAP
  unmap(in);
#ifdef _WIN32
  if (in->handle) CloseHandle((HANDLE) in->handle);
  in->handle = NULL;
#endif
#endif
  in->mapped = 0;
}

#if defined(_WIN32)

static LONG WINAPI fault(EXCEPTION_POINTERS * e) {
  static const char text[] = "bzip3: read error: the input file changed or could not be read\r\n";
  DWORD done;
  if (e->ExceptionRecord->ExceptionCode != EXCEPTION_IN_PAGE_ERROR)
    return EXCEPTION_CONTINUE_SEARCH;
  WriteFile(GetStdHandle(STD_ERROR_HANDLE), text, sizeof text - 1, &done, NULL);
  ExitProcess(1);
  return EXCEPTION_EXECUTE_HANDLER;
}

void bz3_map_guard(const char * const * doomed) {
  (void) doomed;
  SetUnhandledExceptionFilter(fault);
}

#elif defined(HAVE_MMAP) && defined(SIGBUS) && defined(HAVE_SIGACTION)

static const char * const * guarded;

static void fault(int sig) {
  static const char text[] = "bzip3: read error: the input file changed or could not be read\n";
  const char * path = guarded ? *guarded : NULL;
  (void) sig;
  if (write(2, text, sizeof text - 1) < 0) { }
  if (path) unlink(path);
  _exit(1);
}

void bz3_map_guard(const char * const * doomed) {
  struct sigaction sa;
  guarded = doomed;
  memset(&sa, 0, sizeof sa);
  sa.sa_handler = fault;
  sigemptyset(&sa.sa_mask);
  sigaction(SIGBUS, &sa, NULL);
}

#else

void bz3_map_guard(const char * const * doomed) { (void) doomed; }

#endif
