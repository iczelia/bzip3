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

#include "sys/meta.h"

#ifdef _WIN32

#include "sys/win32.h"

void bz3_meta_get(bz3_meta * m, FILE * f) { bz3_win_meta_get(m, f); }

int bz3_meta_changed(const bz3_meta * m, FILE * f) {
  bz3_meta now;
  if (!m->known) return 0;
  bz3_win_meta_get(&now, f);
  return !now.known || now.size_low != m->size_low || now.size_high != m->size_high ||
         CompareFileTime(&now.written, &m->written);
}

int bz3_meta_put(const bz3_meta * m, FILE * in, FILE * out, const char * in_path, const char * out_path) {
  (void) in;  (void) out;
  return bz3_win_meta_put(m, in_path, out_path);
}

#else

#include <errno.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#if defined(HAVE_COPYFILE_H) && defined(HAVE_FCOPYFILE)
#include <copyfile.h>
#elif defined(HAVE_SYS_XATTR_H) && defined(HAVE_FLISTXATTR) && defined(__linux__)
#include <sys/xattr.h>
#define HAVE_LINUX_XATTR 1
#elif defined(HAVE_SYS_ACL_H) && defined(HAVE_ACL_GET_FD) && defined(HAVE_ACL_SET_FD)
#include <sys/acl.h>
#define HAVE_POSIX_ACL 1
#endif

#if defined(HAVE_STRUCT_STAT_ST_ATIM_TV_NSEC)
#define ST_ATIME_NS(st) ((st)->st_atim.tv_nsec)
#define ST_MTIME_NS(st) ((st)->st_mtim.tv_nsec)
#elif defined(HAVE_STRUCT_STAT_ST_ATIMESPEC_TV_NSEC)
#define ST_ATIME_NS(st) ((st)->st_atimespec.tv_nsec)
#define ST_MTIME_NS(st) ((st)->st_mtimespec.tv_nsec)
#elif defined(HAVE_STRUCT_STAT_ST_ATIMENSEC)
#define ST_ATIME_NS(st) ((st)->st_atimensec)
#define ST_MTIME_NS(st) ((st)->st_mtimensec)
#else
#define ST_ATIME_NS(st) 0
#define ST_MTIME_NS(st) 0
#endif

void bz3_meta_get(bz3_meta * m, FILE * f) {
  memset(m, 0, sizeof *m);
  if (fstat(fileno(f), &m->st)) return;
  m->known = 1;
  m->regular = S_ISREG(m->st.st_mode);
}

int bz3_meta_changed(const bz3_meta * m, FILE * f) {
  struct stat now;
  if (!m->known) return 0;
  if (fstat(fileno(f), &now)) return 1;
  return now.st_size != m->st.st_size || now.st_mtime != m->st.st_mtime ||
         ST_MTIME_NS(&now) != ST_MTIME_NS(&m->st);
}

#ifdef HAVE_LINUX_XATTR
static int copy_xattrs(int from, int to) {
  char * names, * value = NULL, * name;
  ssize_t size = flistxattr(from, NULL, 0), n;
  sz cap = 0;
  int ok = 1;
  if (size <= 0) return 1;
  names = malloc((sz) size);
  if (!names) return 0;
  size = flistxattr(from, names, (sz) size);
  for (name = names; size > 0 && name < names + size; name += strlen(name) + 1) {
    const int acl = !strncmp(name, "system.posix_acl_", 17) || !strcmp(name, "system.nfs4_acl");
    n = fgetxattr(from, name, NULL, 0);
    if (n < 0) { if (acl) ok = 0;  continue; }
    if ((sz) n > cap) {
      char * grown = realloc(value, (sz) n);
      if (!grown) { ok = 0;  break; }
      value = grown;  cap = (sz) n;
    }
    n = fgetxattr(from, name, value, cap);
    if (n < 0 || fsetxattr(to, name, value, (sz) n, 0)) { if (acl) ok = 0; }
  }
  free(names);  free(value);
  return ok;
}
#endif

#ifdef HAVE_POSIX_ACL
static int copy_acl(int from, int to) {
  acl_t acl;
  int ok;
#if defined(HAVE_ACL_GET_FD_NP) && defined(ACL_TYPE_NFS4)
  acl = acl_get_fd_np(from, ACL_TYPE_NFS4);
  if (acl) {
    ok = !acl_set_fd_np(to, acl, ACL_TYPE_NFS4);
    acl_free(acl);
    return ok || errno == EOPNOTSUPP;
  }
#endif
  acl = acl_get_fd(from);
  if (!acl) return errno == EOPNOTSUPP || errno == ENOSYS || errno == EINVAL;
  ok = !acl_set_fd(to, acl) || errno == EOPNOTSUPP || errno == ENOSYS;
  acl_free(acl);
  return ok;
}
#endif

int bz3_meta_put(const bz3_meta * m, FILE * in, FILE * out, const char * in_path, const char * out_path) {
  const struct stat * st = &m->st;
  const int from = fileno(in), to = fileno(out);
  mode_t mode = st->st_mode & (S_IRWXU | S_IRWXG | S_IRWXO);
  struct stat now;
  int failed = 0, have_now;
  (void) in_path;  (void) out_path;  (void) from;
  if (!m->known) return 0;

  have_now = !fstat(to, &now);
  if (have_now) {
#ifdef HAVE_FCHOWN
    if (now.st_gid != st->st_gid && fchown(to, (uid_t) -1, st->st_gid)) {
      mode &= ~(mode_t) S_IRWXG | (mode_t) ((mode & S_IRWXO) << 3);
      if (!geteuid()) failed |= BZ3_META_OWNER;
    }
    if (now.st_uid != st->st_uid && fchown(to, st->st_uid, (gid_t) -1)) {
      if (!geteuid()) failed |= BZ3_META_OWNER;
    }
#else
    if (now.st_gid != st->st_gid) mode &= ~(mode_t) S_IRWXG | (mode_t) ((mode & S_IRWXO) << 3);
#endif
  }

#if defined(HAVE_COPYFILE_H) && defined(HAVE_FCOPYFILE)
  if (fcopyfile(from, to, NULL, COPYFILE_ACL | COPYFILE_XATTR) && errno != ENOTSUP)
    failed |= BZ3_META_ACL;
#elif defined(HAVE_LINUX_XATTR)
  if (!copy_xattrs(from, to)) failed |= BZ3_META_ACL;
#elif defined(HAVE_POSIX_ACL)
  if (!copy_acl(from, to)) failed |= BZ3_META_ACL;
#endif

#ifdef HAVE_FCHMOD
  if (fchmod(to, mode)) failed |= BZ3_META_MODE;
#else
  failed |= BZ3_META_MODE;
#endif

#if defined(HAVE_FUTIMENS)
  { struct timespec ts[2];
    ts[0].tv_sec = st->st_atime;  ts[0].tv_nsec = ST_ATIME_NS(st);
    ts[1].tv_sec = st->st_mtime;  ts[1].tv_nsec = ST_MTIME_NS(st);
    if (futimens(to, ts)) failed |= BZ3_META_TIMES; }
#elif defined(HAVE_FUTIMES)
  { struct timeval tv[2];
    tv[0].tv_sec = st->st_atime;  tv[0].tv_usec = ST_ATIME_NS(st) / 1000;
    tv[1].tv_sec = st->st_mtime;  tv[1].tv_usec = ST_MTIME_NS(st) / 1000;
    if (futimes(to, tv)) failed |= BZ3_META_TIMES; }
#else
  failed |= BZ3_META_TIMES;
#endif

#if defined(HAVE_FCHFLAGS) && defined(HAVE_STRUCT_STAT_ST_FLAGS)
  if (have_now && st->st_flags != now.st_flags && fchflags(to, st->st_flags) && !geteuid())
    failed |= BZ3_META_ATTRS;
#endif

  return failed;
}

#endif
