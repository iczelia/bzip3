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

#include "sys/mem.h"

#include <stdio.h>
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif
#ifdef HAVE_SYS_SYSCTL_H
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#endif

#if !defined(_WIN32)

#ifdef __linux__
static uint64_t meminfo(const char * key) {
  FILE * f = fopen("/proc/meminfo", "r");
  const sz n = strlen(key);
  char line[160];
  uint64_t v = 0;
  if (!f) return 0;
  while (fgets(line, sizeof line, f))
    if (!strncmp(line, key, n) && line[n] == ':') {
      v = (uint64_t) strtoull(line + n + 1, NULL, 10) * 1024;
      break;
    }
  fclose(f);
  return v;
}

static uint64_t cgroup(const char * path) {
  FILE * f = fopen(path, "r");
  char line[64], * end;
  uint64_t v = 0;
  if (!f) return 0;
  if (fgets(line, sizeof line, f)) {
    v = (uint64_t) strtoull(line, &end, 10);
    if (end == line) v = 0;
  }
  fclose(f);
  return v;
}
#endif

static uint64_t pages(int name) {
#if defined(HAVE_SYSCONF) && defined(_SC_PAGESIZE)
  const long n = sysconf(name), size = sysconf(_SC_PAGESIZE);
  if (n > 0 && size > 0) return (uint64_t) n * (uint64_t) size;
#else
  (void) name;
#endif
  return 0;
}

static uint64_t lesser(uint64_t a, uint64_t b) {
  return !a ? b : !b ? a : MIN(a, b);
}

static uint64_t limit(void) {
#if defined(HAVE_GETRLIMIT) && defined(RLIMIT_AS)
  struct rlimit r;
  if (!getrlimit(RLIMIT_AS, &r) && r.rlim_cur != RLIM_INFINITY)
    return (uint64_t) r.rlim_cur;
#endif
  return 0;
}

#endif

uint64_t bz3_mem_total(void) {
  uint64_t v = 0;
#if defined(_WIN32) && defined(BZ3_WIN_LEGACY)
  MEMORYSTATUS ms;
  ms.dwLength = sizeof ms;
  GlobalMemoryStatus(&ms);
  v = ms.dwTotalPhys;
#elif defined(_WIN32)
  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof ms;
  if (GlobalMemoryStatusEx(&ms)) v = ms.ullTotalPhys;
#else
#if defined(_SC_PHYS_PAGES)
  v = pages(_SC_PHYS_PAGES);
#endif
#if defined(HAVE_SYSCTLBYNAME)
  if (!v) {
    uint64_t wide = 0;
    unsigned long narrow = 0;
    sz len = sizeof wide;
    if (!sysctlbyname("hw.memsize", &wide, &len, NULL, 0) && len == sizeof wide) v = wide;
    else {
      len = sizeof wide;
      if (!sysctlbyname("hw.physmem", &wide, &len, NULL, 0) && len == sizeof wide) v = wide;
      else {
        len = sizeof narrow;
        if (!sysctlbyname("hw.physmem", &narrow, &len, NULL, 0) && len == sizeof narrow) v = narrow;
      }
    }
  }
#endif
#ifdef __linux__
  v = lesser(v, cgroup("/sys/fs/cgroup/memory.max"));
#endif
  v = lesser(v, limit());
#endif
  return v;
}

uint64_t bz3_mem_avail(void) {
  uint64_t v = 0;
#if defined(_WIN32) && defined(BZ3_WIN_LEGACY)
  MEMORYSTATUS ms;
  ms.dwLength = sizeof ms;
  GlobalMemoryStatus(&ms);
  v = ms.dwAvailPhys;
#elif defined(_WIN32)
  MEMORYSTATUSEX ms;
  ms.dwLength = sizeof ms;
  if (GlobalMemoryStatusEx(&ms)) v = ms.ullAvailPhys;
#else
#ifdef __linux__
  v = meminfo("MemAvailable");
  { const uint64_t max = cgroup("/sys/fs/cgroup/memory.max");
    const uint64_t used = cgroup("/sys/fs/cgroup/memory.current");
    if (max) v = lesser(v, max > used ? max - used : 1); }
#endif
#if defined(_SC_AVPHYS_PAGES)
  if (!v) v = pages(_SC_AVPHYS_PAGES);
#endif
  v = lesser(v, limit());
#endif
  return v;
}
