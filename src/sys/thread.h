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

#ifndef BZ3_THREAD_H
#define BZ3_THREAD_H

#include "sys/common.h"

#ifdef HAVE_THREADS

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
typedef CRITICAL_SECTION bz3_mutex;
typedef struct { HANDLE handle;  void (* fn)(void *);  void * arg; } bz3_thread;
#else
#include <pthread.h>
typedef pthread_mutex_t bz3_mutex;
typedef struct { pthread_t handle;  void (* fn)(void *);  void * arg; } bz3_thread;
#endif

int bz3_thread_start(bz3_thread * t, void (* fn)(void *), void * arg);
void bz3_thread_join(bz3_thread * t);

void bz3_mutex_init(bz3_mutex * m);
void bz3_mutex_lock(bz3_mutex * m);
void bz3_mutex_unlock(bz3_mutex * m);
void bz3_mutex_free(bz3_mutex * m);

typedef struct bz3_queue bz3_queue;
bz3_queue * bz3_queue_new(int cap);
void bz3_queue_put(bz3_queue * q, void * item);
void * bz3_queue_get(bz3_queue * q);
void bz3_queue_free(bz3_queue * q);

#endif

#endif
