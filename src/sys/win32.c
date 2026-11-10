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

#include "sys/win32.h"

#ifdef _WIN32
#include <errno.h>
#include <fcntl.h>
#include <io.h>
#include <sys/stat.h>
#include <wchar.h>
#ifndef BZ3_WIN_LEGACY
#include <aclapi.h>
#endif

extern int _dowildcard;
int _dowildcard = -1;

#ifndef BZ3_WIN_LEGACY

static wchar_t * to_wide(UINT cp, const char * s, int strict) {
  const DWORD flags = strict ? MB_ERR_INVALID_CHARS : 0;
  int n = MultiByteToWideChar(cp, flags, s, -1, NULL, 0);
  wchar_t * w;
  if (!n) { errno = EINVAL;  return NULL; }
  w = malloc((sz) n * sizeof *w);
  if (!w) { errno = ENOMEM;  return NULL; }
  if (!MultiByteToWideChar(cp, flags, s, -1, w, n)) {
    free(w);  errno = EINVAL;  return NULL;
  }
  return w;
}

static char * from_wide(UINT cp, const wchar_t * w) {
  int n = WideCharToMultiByte(cp, 0, w, -1, NULL, 0, NULL, NULL);
  char * s;
  if (!n) { errno = EINVAL;  return NULL; }
  s = malloc((sz) n);
  if (!s) { errno = ENOMEM;  return NULL; }
  if (!WideCharToMultiByte(cp, 0, w, -1, s, n, NULL, NULL)) {
    free(s);  errno = EINVAL;  return NULL;
  }
  return s;
}

wchar_t * bz3_win_wide(const char * s) { return to_wide(CP_UTF8, s, 1); }
char * bz3_win_utf8(const wchar_t * w) { return from_wide(CP_UTF8, w); }

wchar_t * bz3_win_path(const char * s) {
  wchar_t * raw = bz3_win_wide(s), * full, * out;
  DWORD need, got;
  if (!raw) return NULL;
  if (raw[0] == L'\\' && (raw[1] == L'\\' || raw[1] == L'?') &&
      (raw[2] == L'?' || raw[2] == L'.') && raw[3] == L'\\')
    return raw;
  need = GetFullPathNameW(raw, 0, NULL, NULL);
  if (!need || need > 32760) {
    free(raw);  errno = need ? ENAMETOOLONG : EINVAL;  return NULL;
  }
  full = malloc(((sz) need + 1) * sizeof *full);
  if (!full) { free(raw);  errno = ENOMEM;  return NULL; }
  got = GetFullPathNameW(raw, need + 1, full, NULL);
  free(raw);
  if (!got || got > need) {
    free(full);  errno = got ? ENAMETOOLONG : EINVAL;  return NULL;
  }
  if (got >= 4 && full[0] == L'\\' && full[1] == L'\\' && full[2] == L'.' && full[3] == L'\\')
    return full;
  out = malloc(((sz) got + 9) * sizeof *out);
  if (!out) { free(full);  errno = ENOMEM;  return NULL; }
  if (got >= 2 && full[0] == L'\\' && full[1] == L'\\') {
    memcpy(out, L"\\\\?\\UNC\\", 8 * sizeof *out);
    memcpy(out + 8, full + 2, ((sz) got - 1) * sizeof *out);
  } else {
    memcpy(out, L"\\\\?\\", 4 * sizeof *out);
    memcpy(out + 4, full, ((sz) got + 1) * sizeof *out);
  }
  free(full);
  return out;
}

char ** bz3_win_args(int argc, wchar_t ** wargv) {
  char ** argv = calloc((sz) argc + 1, sizeof *argv);
  int i;
  if (!argv) return NULL;
  Fi(argc,
    argv[i] = bz3_win_utf8(wargv[i]);
    if (!argv[i]) return NULL);
  return argv;
}

#endif

FILE * bz3_win_fopen(const char * path, const char * mode) {
#ifdef BZ3_WIN_LEGACY
  return fopen(path, mode);
#else
  wchar_t * w = bz3_win_path(path), wmode[8];
  FILE * f;
  int e, i;
  if (!w) return NULL;
  for (i = 0; mode[i] && i < 7; i++) wmode[i] = (wchar_t) mode[i];
  wmode[i] = 0;
  f = _wfopen(w, wmode);
  e = errno;  free(w);  errno = e;
  return f;
#endif
}

int bz3_win_create(const char * path, int force) {
  const int flags = _O_WRONLY | _O_CREAT | _O_BINARY | (force ? _O_TRUNC : _O_EXCL);
#ifdef BZ3_WIN_LEGACY
  return _open(path, flags, _S_IREAD | _S_IWRITE);
#else
  wchar_t * w = bz3_win_path(path);
  int fd, e;
  if (!w) return -1;
  fd = _wopen(w, flags, _S_IREAD | _S_IWRITE);
  e = errno;  free(w);  errno = e;
  return fd;
#endif
}

int bz3_win_remove(const char * path) {
#ifdef BZ3_WIN_LEGACY
  return remove(path);
#else
  wchar_t * w = bz3_win_path(path);
  int r, e;
  if (!w) return -1;
  r = _wremove(w);
  e = errno;  free(w);  errno = e;
  return r;
#endif
}

int bz3_win_is_dir(const char * path) {
  DWORD attrs;
#ifdef BZ3_WIN_LEGACY
  attrs = GetFileAttributesA(path);
#else
  wchar_t * w = bz3_win_path(path);
  if (!w) return 0;
  attrs = GetFileAttributesW(w);
  free(w);
#endif
  return attrs != INVALID_FILE_ATTRIBUTES && (attrs & FILE_ATTRIBUTE_DIRECTORY);
}

int bz3_win_same_file(FILE * f, const char * path) {
  const DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
  BY_HANDLE_FILE_INFORMATION x, y;
  HANDLE ha = (HANDLE) _get_osfhandle(_fileno(f)), hb;
  int same;
#ifdef BZ3_WIN_LEGACY
  hb = CreateFileA(path, 0, share, NULL, OPEN_EXISTING, 0, NULL);
#else
  wchar_t * w = bz3_win_path(path);
  if (!w) return 0;
  hb = CreateFileW(w, 0, share, NULL, OPEN_EXISTING, 0, NULL);
  free(w);
#endif
  if (hb == INVALID_HANDLE_VALUE) return 0;
  same = ha != INVALID_HANDLE_VALUE && GetFileType(ha) == FILE_TYPE_DISK &&
         GetFileType(hb) == FILE_TYPE_DISK &&
         GetFileInformationByHandle(ha, &x) && GetFileInformationByHandle(hb, &y) &&
         x.dwVolumeSerialNumber == y.dwVolumeSerialNumber &&
         x.nFileIndexHigh == y.nFileIndexHigh && x.nFileIndexLow == y.nFileIndexLow;
  CloseHandle(hb);
  return same;
}

void bz3_win_meta_get(bz3_meta * m, FILE * f) {
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(f));
  BY_HANDLE_FILE_INFORMATION info;
  memset(m, 0, sizeof *m);
  if (h == INVALID_HANDLE_VALUE || GetFileType(h) != FILE_TYPE_DISK) return;
  if (!GetFileInformationByHandle(h, &info)) return;
  m->known = m->regular = 1;
  m->created = info.ftCreationTime;
  m->accessed = info.ftLastAccessTime;
  m->written = info.ftLastWriteTime;
  m->attrs = info.dwFileAttributes;
  m->size_low = info.nFileSizeLow;  m->size_high = info.nFileSizeHigh;
}

#define ATTRS (FILE_ATTRIBUTE_READONLY | FILE_ATTRIBUTE_HIDDEN | FILE_ATTRIBUTE_SYSTEM | \
               FILE_ATTRIBUTE_ARCHIVE | 0x2000 )

#ifdef BZ3_WIN_LEGACY

int bz3_win_meta_put(const bz3_meta * m, const char * in, const char * out) {
  HANDLE h;
  DWORD now;
  int failed = 0;
  (void) in;
  if (!m->known) return 0;
  h = CreateFileA(out, GENERIC_WRITE, FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, NULL);
  if (h == INVALID_HANDLE_VALUE || !SetFileTime(h, &m->created, &m->accessed, &m->written))
    failed |= BZ3_META_TIMES;
  if (h != INVALID_HANDLE_VALUE) CloseHandle(h);
  now = GetFileAttributesA(out);
  if (now == INVALID_FILE_ATTRIBUTES ||
      ((now ^ m->attrs) & ATTRS && !SetFileAttributesA(out, (now & ~(DWORD) ATTRS) | (m->attrs & ATTRS))))
    failed |= BZ3_META_ATTRS;
  return failed;
}

#else

static int privilege(const wchar_t * name) {
  TOKEN_PRIVILEGES tp;
  HANDLE token;
  int ok;
  if (!OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES | TOKEN_QUERY, &token)) return 0;
  tp.PrivilegeCount = 1;
  tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
  ok = LookupPrivilegeValueW(NULL, name, &tp.Privileges[0].Luid) &&
       AdjustTokenPrivileges(token, FALSE, &tp, 0, NULL, NULL) && GetLastError() == ERROR_SUCCESS;
  CloseHandle(token);
  return ok;
}

static int can_restore, can_audit;

static void privileges(void) {
  static int done;
  if (done) return;
  done = 1;
  privilege(L"SeBackupPrivilege");
  can_restore = privilege(L"SeRestorePrivilege");
  privilege(L"SeTakeOwnershipPrivilege");
  can_audit = privilege(L"SeSecurityPrivilege");
}

static HANDLE open_meta(const char * path, DWORD access) {
  wchar_t * w = bz3_win_path(path);
  HANDLE h;
  if (!w) return INVALID_HANDLE_VALUE;
  h = CreateFileW(w, access, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, NULL,
                  OPEN_EXISTING, FILE_FLAG_BACKUP_SEMANTICS, NULL);
  free(w);
  return h;
}

static int copy_security(const char * in, const char * out, SECURITY_INFORMATION what) {
  const int audit = (what & SACL_SECURITY_INFORMATION) != 0;
  PSECURITY_DESCRIPTOR sd = NULL;
  PSID owner = NULL, group = NULL;
  PACL dacl = NULL, sacl = NULL;
  SECURITY_DESCRIPTOR_CONTROL control = 0;
  DWORD revision, access = 0;
  HANDLE from, to;
  int ok = 0;

  from = open_meta(in, READ_CONTROL | (audit ? ACCESS_SYSTEM_SECURITY : 0));
  if (from == INVALID_HANDLE_VALUE) return 0;
  if (GetSecurityInfo(from, SE_FILE_OBJECT, what, &owner, &group, &dacl, &sacl, &sd) != ERROR_SUCCESS) {
    CloseHandle(from);  return 0;
  }
  CloseHandle(from);

  GetSecurityDescriptorControl(sd, &control, &revision);
  if (what & DACL_SECURITY_INFORMATION) {
    what |= control & SE_DACL_PROTECTED ? PROTECTED_DACL_SECURITY_INFORMATION
                                        : UNPROTECTED_DACL_SECURITY_INFORMATION;
    access |= WRITE_DAC;
  }
  if (audit) {
    what |= control & SE_SACL_PROTECTED ? PROTECTED_SACL_SECURITY_INFORMATION
                                        : UNPROTECTED_SACL_SECURITY_INFORMATION;
    access |= ACCESS_SYSTEM_SECURITY;
  }
  if (what & (OWNER_SECURITY_INFORMATION | GROUP_SECURITY_INFORMATION)) access |= WRITE_OWNER;

  to = open_meta(out, access);
  if (to != INVALID_HANDLE_VALUE) {
    ok = SetSecurityInfo(to, SE_FILE_OBJECT, what, owner, group, dacl, sacl) == ERROR_SUCCESS;
    CloseHandle(to);
  }
  LocalFree(sd);
  return ok;
}

int bz3_win_meta_put(const bz3_meta * m, const char * in, const char * out) {
  const SECURITY_INFORMATION owner = OWNER_SECURITY_INFORMATION | GROUP_SECURITY_INFORMATION;
  HANDLE h;
  DWORD now;
  int failed = 0;
  wchar_t * w;
  if (!m->known) return 0;
  privileges();

  if (!(can_restore && can_audit &&
        copy_security(in, out, owner | DACL_SECURITY_INFORMATION | SACL_SECURITY_INFORMATION))) {
    if (!copy_security(in, out, DACL_SECURITY_INFORMATION)) failed |= BZ3_META_ACL;
    if (!copy_security(in, out, owner) && !copy_security(in, out, OWNER_SECURITY_INFORMATION) &&
        can_restore)
      failed |= BZ3_META_OWNER;
  }

  h = open_meta(out, FILE_WRITE_ATTRIBUTES);
  if (h == INVALID_HANDLE_VALUE || !SetFileTime(h, &m->created, &m->accessed, &m->written))
    failed |= BZ3_META_TIMES;
  if (h != INVALID_HANDLE_VALUE) CloseHandle(h);

  w = bz3_win_path(out);
  now = w ? GetFileAttributesW(w) : INVALID_FILE_ATTRIBUTES;
  if (now == INVALID_FILE_ATTRIBUTES ||
      ((now ^ m->attrs) & ATTRS && !SetFileAttributesW(w, (now & ~(DWORD) ATTRS) | (m->attrs & ATTRS))))
    failed |= BZ3_META_ATTRS;
  free(w);
  return failed;
}

#endif

static char * format(char * buf, sz cap, const char * fmt, va_list ap) {
  char * s = buf;
  for (;;) {
    va_list aq;
    int n;
    va_copy(aq, ap);  n = vsnprintf(s, cap, fmt, aq);  va_end(aq);
    if (n >= 0 && (sz) n < cap) return s;
    if (s != buf) free(s);
    cap = n >= 0 ? (sz) n + 1 : cap * 2;
    if (cap > MiB(1)) return NULL;
    s = malloc(cap);
    if (!s) return NULL;
  }
}

static int console(HANDLE h, const char * s) {
  DWORD done;
  int ok;
#ifdef BZ3_WIN_LEGACY
  char * out;
  int n = MultiByteToWideChar(CP_ACP, 0, s, -1, NULL, 0), m;
  wchar_t * w;
  if (!n) return 0;
  w = malloc((sz) n * sizeof *w);  out = malloc((sz) n * 4);
  ok = w && out && MultiByteToWideChar(CP_ACP, 0, s, -1, w, n) &&
       (m = WideCharToMultiByte(GetConsoleOutputCP(), 0, w, -1, out, n * 4, NULL, NULL)) > 0 &&
       WriteFile(h, out, (DWORD) (m - 1), &done, NULL);
  free(w);  free(out);
#else
  wchar_t * w = to_wide(CP_UTF8, s, 0), * p;
  if (!w) return 0;
  ok = 1;
  for (p = w; ok && *p; p += done) {
    const sz left = wcslen(p);
    ok = WriteConsoleW(h, p, (DWORD) MIN(left, 4096), &done, NULL) && done;
  }
  free(w);
#endif
  return ok;
}

void bz3_win_vprintf(FILE * to, const char * fmt, va_list ap) {
  char buf[512], * s = format(buf, sizeof buf, fmt, ap);
  HANDLE h = (HANDLE) _get_osfhandle(_fileno(to));
  DWORD mode;
  if (!s) { vfprintf(to, fmt, ap);  return; }
  fflush(to);
  if (h == INVALID_HANDLE_VALUE || !GetConsoleMode(h, &mode) || !console(h, s))
    fputs(s, to);
  if (s != buf) free(s);
}

#else

typedef int bz3_win32_unused;

#endif
