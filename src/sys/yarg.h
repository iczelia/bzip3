/*  Written by Kamila Szewczyk (k@iczelia.net), released to
    the public domain (0BSD).  */

#ifndef YARG_H
#define YARG_H

#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

typedef enum {
  no_argument, required_argument, optional_argument
} yarg_arg_type;

typedef struct {
  int opt;  yarg_arg_type type;  const char * long_opt;
} yarg_options;

typedef enum {
  YARG_STYLE_WINDOWS, YARG_STYLE_UNIX, YARG_STYLE_UNIX_SHORT
} yarg_style;

typedef struct {
  bool dash_dash;  yarg_style style;
} yarg_settings;

typedef struct {
  int opt;  const char * long_opt;  char * arg;
} yarg_option;

typedef struct {
  yarg_option * args;
  int argc;
  char ** pos_args;
  int pos_argc;
  char * error;
} yarg_result;

static char yarg_oom[] = "Out of memory";

#if defined(__GNUC__)
__attribute__((__format__(__printf__, 2, 0)))
#endif
static int yarg_asprintf(char ** strp, const char * fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  int len = vsnprintf(NULL, 0, fmt, ap);
  va_end(ap);
  if (len < 0) { *strp = yarg_oom;  return -1; }
  *strp = (char *) malloc((size_t) len + 1);
  if (!*strp) { *strp = yarg_oom;  return -1; }
  va_start(ap, fmt);
  len = vsnprintf(*strp, (size_t) len + 1, fmt, ap);
  va_end(ap);
  return len;
}

static yarg_options * yarg_find_long(yarg_options * opt, const char * name,
                                     size_t len, int * ambiguous) {
  yarg_options * o = NULL;
  *ambiguous = 0;
  if (!len) return NULL;
  for (int j = 0; opt[j].opt; j++) {
    if (!opt[j].long_opt || strncmp(opt[j].long_opt, name, len)) continue;
    if (!opt[j].long_opt[len]) { *ambiguous = 0;  return &opt[j]; }
    if (!o) o = &opt[j];
    else if (o->opt != opt[j].opt) *ambiguous = 1;
  }
  return *ambiguous ? NULL : o;
}

static char * yarg_strdup(const char * str) {
  char * new_str = (char *) calloc(strlen(str) + 1, 1);
  if (!new_str) return NULL;
  strcpy(new_str, str);
  return new_str;
}

static int yarg_parse_unix(int argc, char * argv[], yarg_options opt[],
                           yarg_result * res, bool dash_dash) {
  size_t no_args = 0, no_pos_args = 0;
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && argv[i][1]) {
      if (argv[i][1] == '-') {
        if (dash_dash && argv[i][2] == '\0')
          { no_pos_args += (size_t) (argc - i - 1);  break; }
        char * long_opt = argv[i] + 2;
        size_t len = 0;
    while (long_opt[len] && long_opt[len] != '=') len++;
        int amb = 0;
        yarg_options * o = yarg_find_long(opt, long_opt, len, &amb);
        if (!o) {
          if (amb)
            yarg_asprintf(&res->error, "ambiguous option '--%.*s'",
                          (int) len, long_opt);
          else
            yarg_asprintf(&res->error, "unknown option '--%.*s'",
                          (int) len, long_opt);
          return 0;
        }
        if (o->type == no_argument) {
          if (long_opt[len] == '=') {
            yarg_asprintf(&res->error, "option '--%s' takes no argument",
                          o->long_opt);
            return 0;
          }
        } else if (o->type == required_argument) {
          if (long_opt[len] == '=') ;
          else if (i + 1 < argc && argv[i + 1][0] != '-') i++;
          else {
            yarg_asprintf(&res->error, "option '--%s' needs an argument",
                          o->long_opt);
            return 0;
          }
        } else if (o->type == optional_argument) {
          if (long_opt[len] == '=') ;
          else if (i + 1 < argc && argv[i + 1][0] != '-') i++;
        }
        no_args++;
      } else {
        for (int j = 1; argv[i][j]; j++) {
          unsigned char c = (unsigned char) argv[i][j];
          yarg_options * o = NULL;
          for (int k = 0; opt[k].opt; k++)
            if (opt[k].opt == c)
              { o = &opt[k]; break; }
          if (!o) {
            yarg_asprintf(&res->error, "unknown option '-%c'", c);
            return 0;
          }
          if (o->type == required_argument) {
            if (argv[i][j + 1])  ;
            else if (i + 1 < argc && argv[i + 1][0] != '-') i++;
            else {
              yarg_asprintf(&res->error, "option '-%c' needs an argument", c);
              return 0;
            }
            no_args++;
            break;
          } else if(o->type == optional_argument) {
            if (argv[i][j + 1])
              {   no_args++;  break; }
            else if (i + 1 < argc && argv[i + 1][0] != '-')
              { i++;  no_args++;  break; }
          }
          no_args++;
        }
      }
    } else no_pos_args++;
  }
  res->args = (yarg_option *) calloc(no_args + 1, sizeof(yarg_option));
  res->pos_args = (char **) calloc(no_pos_args + 1, sizeof(char *));
  if(!res->args || !res->pos_args)
    { res->error = yarg_oom;  return 0; }
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-' && argv[i][1]) {
      if (argv[i][1] == '-') {
        if (dash_dash && argv[i][2] == '\0') {
          for (int j = i + 1; j < argc; j++)
            if(!(res->pos_args[res->pos_argc++] = yarg_strdup(argv[j])))
              { res->error = yarg_oom;  return 0; }
          break;
        }
        char * long_opt = argv[i] + 2;
        size_t len = 0;
    while (long_opt[len] && long_opt[len] != '=') len++;
        int amb = 0;
        yarg_options * o = yarg_find_long(opt, long_opt, len, &amb);
        if (!o) {
          if (amb)
            yarg_asprintf(&res->error, "ambiguous option '--%.*s'",
                          (int) len, long_opt);
          else
            yarg_asprintf(&res->error, "unknown option '--%.*s'",
                          (int) len, long_opt);
          return 0;
        }
        res->args[res->argc].opt = o->opt;
        res->args[res->argc].long_opt = o->long_opt;
        if (o->type == required_argument || o->type == optional_argument) {
          if (long_opt[len] == '=') {
            if(!(res->args[res->argc].arg = yarg_strdup(long_opt + len + 1)))
              { res->error = yarg_oom;  return 0; }
          } else if (i + 1 < argc && argv[i + 1][0] != '-') {
            if(!(res->args[res->argc].arg = yarg_strdup(argv[++i])))
              { res->error = yarg_oom;  return 0; }
          }
        }
        res->argc++;
      } else {
        for (int j = 1; argv[i][j]; j++) {
          unsigned char c = (unsigned char) argv[i][j];
          yarg_options * o = NULL;
          for (int k = 0; opt[k].opt; k++)
            if (opt[k].opt == c)
              { o = &opt[k]; break; }
          if (!o) {
            yarg_asprintf(&res->error, "unknown option '-%c'", c);
            return 0;
          }
          res->args[res->argc].opt = c;
          res->args[res->argc].long_opt = o->long_opt;
          if (o->type == required_argument || o->type == optional_argument) {
            if (argv[i][j + 1]) {
              if(!(res->args[res->argc++].arg = yarg_strdup(argv[i] + j + 1))) 
                { res->error = yarg_oom;  return 0; }
              break;
            } else if (i + 1 < argc && argv[i + 1][0] != '-') {
              if(!(res->args[res->argc++].arg = yarg_strdup(argv[++i])))
                { res->error = yarg_oom;  return 0; }
              break;
            }
          }
          res->argc++;
        }
      }
    } else if(!(res->pos_args[res->pos_argc++] = yarg_strdup(argv[i])))
      { res->error = yarg_oom;  return 0; }
  }
  return 1;
}

static int yarg_parse_unix_short(int argc, char * argv[], yarg_options opt[],
                                 yarg_result * res, bool dash_dash, char opt_char) {
  size_t no_args = 0, no_pos_args = 0;
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == opt_char && (argv[i][1] || dash_dash)) {
      if (dash_dash && argv[i][1] == '\0')
        { no_pos_args += (size_t) (argc - i - 1);  break; }
      char * long_opt = argv[i] + 1;
      size_t len = 0;
    while (long_opt[len] && long_opt[len] != '=') len++;
      int amb = 0;
      yarg_options * o = yarg_find_long(opt, long_opt, len, &amb);
      if (!o) {
        if (amb)
          yarg_asprintf(&res->error, "ambiguous option '%c%.*s'",
                        opt_char, (int) len, long_opt);
        else
          yarg_asprintf(&res->error, "unknown option '%c%.*s'",
                        opt_char, (int) len, long_opt);
        return 0;
      }
      if (o->type == no_argument) {
        if (long_opt[len] == '=') {
          yarg_asprintf(&res->error, "option '%c%s' takes no argument", opt_char, o->long_opt);
          return 0;
        }
      } else if (o->type == required_argument) {
        if (long_opt[len] == '=') ;
        else if (i + 1 < argc && argv[i + 1][0] != opt_char) i++;
        else {
          yarg_asprintf(&res->error, "option '%c%s' needs an argument", opt_char, o->long_opt);
          return 0;
        }
      } else if (o->type == optional_argument) {
        if (long_opt[len] == '=') ;
        else if (i + 1 < argc && argv[i + 1][0] != opt_char) i++;
      }
      no_args++;
    } else no_pos_args++;
  }
  res->args = (yarg_option *) calloc(no_args + 1, sizeof(yarg_option));
  res->pos_args = (char **) calloc(no_pos_args + 1, sizeof(char *));
  if (!res->args || !res->pos_args)
    { res->error = yarg_oom;  return 0; }
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == opt_char && (argv[i][1] || dash_dash)) {
      if (dash_dash && argv[i][1] == '\0') {
        for (int j = i + 1; j < argc; j++)
          if(!(res->pos_args[res->pos_argc++] = yarg_strdup(argv[j])))
            { res->error = yarg_oom;  return 0; }
        break;
      }
      char * long_opt = argv[i] + 1;
      size_t len = 0;
    while (long_opt[len] && long_opt[len] != '=') len++;
      int amb = 0;
      yarg_options * o = yarg_find_long(opt, long_opt, len, &amb);
      if (!o) {
        if (amb)
          yarg_asprintf(&res->error, "ambiguous option '%c%.*s'",
                        opt_char, (int) len, long_opt);
        else
          yarg_asprintf(&res->error, "unknown option '%c%.*s'",
                        opt_char, (int) len, long_opt);
        return 0;
      }
      res->args[res->argc].opt = o->opt;
      res->args[res->argc].long_opt = o->long_opt;
      if (o->type == required_argument || o->type == optional_argument) {
        if (long_opt[len] == '=') {
          if(!(res->args[res->argc].arg = yarg_strdup(long_opt + len + 1)))
            { res->error = yarg_oom;  return 0; }
        } else if (i + 1 < argc && argv[i + 1][0] != opt_char) {
          if(!(res->args[res->argc].arg = yarg_strdup(argv[++i])))
            { res->error = yarg_oom;  return 0; }
        }
      }
      res->argc++;
    } else if(!(res->pos_args[res->pos_argc++] = yarg_strdup(argv[i])))
      { res->error = yarg_oom;  return 0; }
  }

  return 1;
}

static inline void yarg_destroy(yarg_result * r) {
  if(r) {
    if(r->args) for (int i = 0; i < r->argc; i++) free(r->args[i].arg);
    free(r->args);
    if(r->pos_args) for (int i = 0; i < r->pos_argc; i++) free(r->pos_args[i]);
    free(r->pos_args);
    if (r->error != yarg_oom) free(r->error);
  }
  free(r);
}

static inline yarg_result * yarg_parse(int argc, char * argv[], yarg_options opt[], yarg_settings settings) {
  yarg_result * res = (yarg_result *) calloc(1, sizeof(yarg_result));
  if (!res) return NULL;
  switch (settings.style) {
    case YARG_STYLE_WINDOWS:
      yarg_parse_unix_short(argc, argv, opt, res, false, '/');
      break;
    case YARG_STYLE_UNIX:
      yarg_parse_unix(argc, argv, opt, res, settings.dash_dash);
      break;
    case YARG_STYLE_UNIX_SHORT:
      yarg_parse_unix_short(argc, argv, opt, res, settings.dash_dash, '-');
      break;
  }
  return res;
}

#endif

#ifndef YARG_VERB_H
#define YARG_VERB_H

typedef struct {
  int verb;  const char * name;
  const yarg_options * opts;  const yarg_options * inherit;
} yarg_verb;

typedef enum {
  YARG_VERB_OK, YARG_VERB_ABSENT, YARG_VERB_UNKNOWN, YARG_VERB_AMBIGUOUS
} yarg_verb_status;

typedef struct {
  yarg_result * res;
  const yarg_verb * verb;
  yarg_verb_status status;
  char * cands;
  yarg_options * merged;
} yarg_verb_result;

static yarg_verb_status yarg_find_verb(const yarg_verb * v, const char * name,
                                       const yarg_verb ** out, char ** cands) {
  size_t len = strlen(name), need = 1;
  const yarg_verb * first = NULL;
  int n = 0;
  if (!len) return YARG_VERB_UNKNOWN;
  for (int i = 0; v[i].name; i++) {
    if (strncmp(v[i].name, name, len)) continue;
    if (!v[i].name[len]) { *out = &v[i];  return YARG_VERB_OK; }
    if (!first) first = &v[i];
    n++;  need += strlen(v[i].name) + 2;
  }
  if (n == 1) { *out = first;  return YARG_VERB_OK; }
  if (!n) return YARG_VERB_UNKNOWN;
  if ((*cands = (char *) malloc(need)) != NULL) {
    char * p = *cands;  *p = '\0';
    for (int i = 0; v[i].name; i++) {
      if (strncmp(v[i].name, name, len)) continue;
      if (p != *cands) { *p++ = ',';  *p++ = ' '; }
      strcpy(p, v[i].name);  p += strlen(v[i].name);
    }
  }
  return YARG_VERB_AMBIGUOUS;
}

static size_t yarg_opt_len(const yarg_options * o) {
  size_t n = 0;  while (o && o[n].opt) n++;  return n;
}

static yarg_options * yarg_merge_opts(const yarg_options * a,
                                      const yarg_options * b,
                                      const yarg_options * c) {
  size_t na = yarg_opt_len(a), nb = yarg_opt_len(b), nc = yarg_opt_len(c);
  size_t k = 0;
  yarg_options * m = (yarg_options *) calloc(na + nb + nc + 1, sizeof *m);
  if (!m) return NULL;
  for (size_t i = 0; i < na; i++) m[k++] = a[i];
  for (size_t i = 0; i < nb; i++) m[k++] = b[i];
  for (size_t i = 0; i < nc; i++) m[k++] = c[i];
  return m;
}

static inline yarg_verb_result * yarg_parse_verb(int argc, char * argv[],
                                                 const yarg_verb * verbs,
                                                 const yarg_options * global,
                                                 yarg_settings settings) {
  yarg_verb_result * r = (yarg_verb_result *) calloc(1, sizeof *r);
  int shift = 0;
  if (!r) return NULL;
  r->status = YARG_VERB_ABSENT;
  if (argc > 1 && argv[1][0] != '-') {
    r->status = yarg_find_verb(verbs, argv[1], &r->verb, &r->cands);
    if (r->status == YARG_VERB_AMBIGUOUS) return r;
    if (r->status == YARG_VERB_OK) shift = 1;
  }
  r->merged = yarg_merge_opts(global, r->verb ? r->verb->opts : NULL,
                              r->verb ? r->verb->inherit : NULL);
  if (r->merged)
    r->res = yarg_parse(argc - shift, argv + shift, r->merged, settings);
  return r;
}

static inline void yarg_verb_destroy(yarg_verb_result * r) {
  if (!r) return;
  yarg_destroy(r->res);
  free(r->merged);  free(r->cands);  free(r);
}

#endif
