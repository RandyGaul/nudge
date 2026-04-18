// See LICENSE for licensing info.
// threads.h -- minimal cross-platform thread spawn primitive.
//
// Follows nudge's opaque-handle convention: Thread is a uint64_t id, with
// the platform-specific payload (Win32 HANDLE or pthread_t) hidden behind
// it. Covers just what the worker pool needs: spawn a thread running
// fn(arg), join on shutdown. Synchronisation lives elsewhere (stdatomic.h
// for atomics, simd_pause() for spin hints, _Thread_local for TLS).
//
// Uses only libc + platform SDK. Zero external deps.
#ifndef NUDGE_THREADS_H
#define NUDGE_THREADS_H

#include <stdlib.h>
#include <stdint.h>

#ifdef _WIN32
	#ifndef WIN32_LEAN_AND_MEAN
	#define WIN32_LEAN_AND_MEAN
	#endif
	#ifndef NOMINMAX
	#define NOMINMAX
	#endif
	#include <windows.h>
#else
	#include <pthread.h>
#endif

typedef struct Thread { uint64_t id; } Thread;

// Caller-facing thread entry. Runs until it returns; caller joins via
// nudge_thread_join (or lets process exit clean up).
typedef void (*nudge_thread_fn)(void* arg);

typedef struct NudgeThreadImpl
{
#ifdef _WIN32
	HANDLE h;
#else
	pthread_t p;
#endif
	nudge_thread_fn fn;
	void* arg;
} NudgeThreadImpl;

#ifdef _WIN32
static DWORD WINAPI nudge_thread_trampoline(LPVOID p)
{
	NudgeThreadImpl* impl = (NudgeThreadImpl*)p;
	impl->fn(impl->arg);
	return 0;
}
#else
static void* nudge_thread_trampoline(void* p)
{
	NudgeThreadImpl* impl = (NudgeThreadImpl*)p;
	impl->fn(impl->arg);
	return NULL;
}
#endif

// Spawn a thread. Returns Thread{0} on failure.
static inline Thread nudge_thread_create(nudge_thread_fn fn, void* arg)
{
	NudgeThreadImpl* impl = (NudgeThreadImpl*)malloc(sizeof(NudgeThreadImpl));
	if (!impl) return (Thread){0};
	impl->fn = fn;
	impl->arg = arg;
#ifdef _WIN32
	impl->h = CreateThread(NULL, 0, nudge_thread_trampoline, impl, 0, NULL);
	if (!impl->h) { free(impl); return (Thread){0}; }
#else
	if (pthread_create(&impl->p, NULL, nudge_thread_trampoline, impl) != 0) { free(impl); return (Thread){0}; }
#endif
	return (Thread){ (uint64_t)(uintptr_t)impl };
}

// Block until the thread returns, then release resources.
static inline void nudge_thread_join(Thread t)
{
	if (!t.id) return;
	NudgeThreadImpl* impl = (NudgeThreadImpl*)(uintptr_t)t.id;
#ifdef _WIN32
	WaitForSingleObject(impl->h, INFINITE);
	CloseHandle(impl->h);
#else
	pthread_join(impl->p, NULL);
#endif
	free(impl);
}

#endif // NUDGE_THREADS_H
