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

#include "sys/thread.h"

#ifdef HAVE_THREADS

struct bz3_queue {
  void ** items;
  int cap, head, count;
  bz3_mutex lock;
#ifdef _WIN32
  HANDLE filled;
#else
  pthread_cond_t filled;
#endif
};

#ifdef _WIN32

static DWORD WINAPI thread_main(LPVOID arg) {
  bz3_thread * t = arg;
  t->fn(t->arg);
  return 0;
}

int bz3_thread_start(bz3_thread * t, void (* fn)(void *), void * arg) {
  t->fn = fn;  t->arg = arg;
  t->handle = CreateThread(NULL, 0, thread_main, t, 0, NULL);
  return t->handle != NULL;
}

void bz3_thread_join(bz3_thread * t) {
  WaitForSingleObject(t->handle, INFINITE);
  CloseHandle(t->handle);
}

void bz3_mutex_init(bz3_mutex * m) { InitializeCriticalSection(m); }
void bz3_mutex_lock(bz3_mutex * m) { EnterCriticalSection(m); }
void bz3_mutex_unlock(bz3_mutex * m) { LeaveCriticalSection(m); }
void bz3_mutex_free(bz3_mutex * m) { DeleteCriticalSection(m); }

#else

static void * thread_main(void * arg) {
  bz3_thread * t = arg;
  t->fn(t->arg);
  return NULL;
}

int bz3_thread_start(bz3_thread * t, void (* fn)(void *), void * arg) {
  t->fn = fn;  t->arg = arg;
  return !pthread_create(&t->handle, NULL, thread_main, t);
}

void bz3_thread_join(bz3_thread * t) { pthread_join(t->handle, NULL); }

void bz3_mutex_init(bz3_mutex * m) { pthread_mutex_init(m, NULL); }
void bz3_mutex_lock(bz3_mutex * m) { pthread_mutex_lock(m); }
void bz3_mutex_unlock(bz3_mutex * m) { pthread_mutex_unlock(m); }
void bz3_mutex_free(bz3_mutex * m) { pthread_mutex_destroy(m); }

#endif

bz3_queue * bz3_queue_new(int cap) {
  bz3_queue * q = calloc(1, sizeof *q);
  if (!q) return NULL;
  q->items = calloc((sz) cap, sizeof *q->items);
  if (!q->items) { free(q);  return NULL; }
  q->cap = cap;
  bz3_mutex_init(&q->lock);
#ifdef _WIN32
  q->filled = CreateSemaphoreA(NULL, 0, cap, NULL);
  if (!q->filled) { bz3_mutex_free(&q->lock);  free(q->items);  free(q);  return NULL; }
#else
  pthread_cond_init(&q->filled, NULL);
#endif
  return q;
}

void bz3_queue_put(bz3_queue * q, void * item) {
  bz3_mutex_lock(&q->lock);
  q->items[(q->head + q->count++) % q->cap] = item;
#ifdef _WIN32
  bz3_mutex_unlock(&q->lock);
  ReleaseSemaphore(q->filled, 1, NULL);
#else
  pthread_cond_signal(&q->filled);
  bz3_mutex_unlock(&q->lock);
#endif
}

void * bz3_queue_get(bz3_queue * q) {
  void * item;
#ifdef _WIN32
  WaitForSingleObject(q->filled, INFINITE);
  bz3_mutex_lock(&q->lock);
#else
  bz3_mutex_lock(&q->lock);
  while (!q->count) pthread_cond_wait(&q->filled, &q->lock);
#endif
  item = q->items[q->head];
  q->head = (q->head + 1) % q->cap;  q->count--;
  bz3_mutex_unlock(&q->lock);
  return item;
}

void bz3_queue_free(bz3_queue * q) {
  if (!q) return;
#ifdef _WIN32
  CloseHandle(q->filled);
#else
  pthread_cond_destroy(&q->filled);
#endif
  bz3_mutex_free(&q->lock);
  free(q->items);  free(q);
}

#else

typedef int bz3_thread_unused;

#endif
