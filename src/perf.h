// perf.h -- high-resolution performance counter (platform abstraction)
#ifndef PERF_H
#define PERF_H

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#undef near
#undef far
static double perf_freq;
static void perf_init() { LARGE_INTEGER f; QueryPerformanceFrequency(&f); perf_freq = (double)f.QuadPart; }
static double perf_now() { LARGE_INTEGER c; QueryPerformanceCounter(&c); return (double)c.QuadPart / perf_freq; }
#else
#include <time.h>
static void perf_init() {}
static double perf_now() { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
#endif

#endif
