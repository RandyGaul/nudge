// See LICENSE for licensing info.
// simd.h -- Cross-platform 4-wide float SIMD abstraction.
//
// Provides a thin wrapper over SSE (x86), NEON (ARM), or scalar fallback.
// All operations are on 4-wide float lanes. The abstraction is minimal:
// types and operations map 1:1 to hardware intrinsics.
//
// Usage: include this header, then use simd4f / simd4i types and simd_* ops.
// The v3 type in vmath.h is built on top of this.
// -----------------------------------------------------------------------------
#ifndef SIMD_H
#define SIMD_H

// Platform detection.
#if defined(__SSE2__) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2)
	#define SIMD_SSE 1
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
	#define SIMD_NEON 1
#else
	#define SIMD_SCALAR 1
#endif

// Includes.
#if SIMD_SSE
	#include <xmmintrin.h>
	#include <emmintrin.h>
	#include <smmintrin.h>
#elif SIMD_NEON
	#include <arm_neon.h>
#endif

// Compiler hints.
#if defined(_MSC_VER)
	#define SIMD_FORCEINLINE __forceinline
	#define SIMD_NOINLINE    __declspec(noinline)
	#define SIMD_ASSUME(x)   __assume(x)
	#define SIMD_VECTORCALL  __vectorcall
#elif defined(__GNUC__) || defined(__clang__)
	#define SIMD_FORCEINLINE __attribute__((always_inline)) inline
	#define SIMD_NOINLINE    __attribute__((noinline))
	#define SIMD_ASSUME(x)   do { if (!(x)) __builtin_unreachable(); } while(0)
	#define SIMD_VECTORCALL
#else
	#define SIMD_FORCEINLINE inline
	#define SIMD_NOINLINE
	#define SIMD_ASSUME(x)
	#define SIMD_VECTORCALL
#endif

// -----------------------------------------------------------------------------
// Types.

#if SIMD_SSE

typedef __m128  simd4f;
typedef __m128i simd4i;

#elif SIMD_NEON

typedef float32x4_t simd4f;
typedef int32x4_t   simd4i;

#else // scalar

typedef struct simd4f { float v[4]; } simd4f;
typedef struct simd4i { int   v[4]; } simd4i;

#endif

// -----------------------------------------------------------------------------
// 4-wide float operations.

// Set all lanes to the same value.
static inline simd4f simd_set1(float x)
{
#if SIMD_SSE
	return _mm_set1_ps(x);
#elif SIMD_NEON
	return vdupq_n_f32(x);
#else
	return (simd4f){{ x, x, x, x }};
#endif
}

// Set individual lanes.
static inline simd4f simd_set(float x, float y, float z, float w)
{
#if SIMD_SSE
	return _mm_set_ps(w, z, y, x);
#elif SIMD_NEON
	float v[4] = { x, y, z, w };
	return vld1q_f32(v);
#else
	return (simd4f){{ x, y, z, w }};
#endif
}

static inline simd4f simd_zero()
{
#if SIMD_SSE
	return _mm_setzero_ps();
#elif SIMD_NEON
	return vdupq_n_f32(0);
#else
	return (simd4f){{ 0, 0, 0, 0 }};
#endif
}

// Arithmetic.
static inline simd4f simd_add(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_add_ps(a, b);
#elif SIMD_NEON
	return vaddq_f32(a, b);
#else
	return (simd4f){{ a.v[0]+b.v[0], a.v[1]+b.v[1], a.v[2]+b.v[2], a.v[3]+b.v[3] }};
#endif
}

static inline simd4f simd_sub(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_sub_ps(a, b);
#elif SIMD_NEON
	return vsubq_f32(a, b);
#else
	return (simd4f){{ a.v[0]-b.v[0], a.v[1]-b.v[1], a.v[2]-b.v[2], a.v[3]-b.v[3] }};
#endif
}

static inline simd4f simd_mul(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_mul_ps(a, b);
#elif SIMD_NEON
	return vmulq_f32(a, b);
#else
	return (simd4f){{ a.v[0]*b.v[0], a.v[1]*b.v[1], a.v[2]*b.v[2], a.v[3]*b.v[3] }};
#endif
}

static inline simd4f simd_div(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_div_ps(a, b);
#elif SIMD_NEON
	return vdivq_f32(a, b);
#else
	return (simd4f){{ a.v[0]/b.v[0], a.v[1]/b.v[1], a.v[2]/b.v[2], a.v[3]/b.v[3] }};
#endif
}

static inline simd4f simd_min(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_min_ps(a, b);
#elif SIMD_NEON
	return vminq_f32(a, b);
#else
	simd4f r;
	for (int i = 0; i < 4; i++) r.v[i] = a.v[i] < b.v[i] ? a.v[i] : b.v[i];
	return r;
#endif
}

static inline simd4f simd_max(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_max_ps(a, b);
#elif SIMD_NEON
	return vmaxq_f32(a, b);
#else
	simd4f r;
	for (int i = 0; i < 4; i++) r.v[i] = a.v[i] > b.v[i] ? a.v[i] : b.v[i];
	return r;
#endif
}

// Negation.
static inline simd4f simd_neg(simd4f a) { return simd_sub(simd_zero(), a); }

// Extract first lane to scalar.
static inline float simd_get_x(simd4f a)
{
#if SIMD_SSE
	return _mm_cvtss_f32(a);
#elif SIMD_NEON
	return vgetq_lane_f32(a, 0);
#else
	return a.v[0];
#endif
}

// Broadcast lane to all positions. Lane index must be a compile-time constant.
#if SIMD_SSE
	#define simd_splat(v, lane) _mm_shuffle_ps(v, v, _MM_SHUFFLE(lane, lane, lane, lane))
#elif SIMD_NEON
	#define simd_splat(v, lane) vdupq_laneq_f32(v, lane)
#else
	#define simd_splat(v, lane) simd_set1((v).v[lane])
#endif

// Sqrt of first lane only (scalar sqrt in SIMD register).
static inline simd4f simd_sqrt_ss(simd4f a)
{
#if SIMD_SSE
	return _mm_sqrt_ss(a);
#elif SIMD_NEON
	float s = sqrtf(vgetq_lane_f32(a, 0));
	return vsetq_lane_f32(s, a, 0);
#else
	simd4f r = a; r.v[0] = sqrtf(a.v[0]); return r;
#endif
}

// Sqrt of all lanes.
static inline simd4f simd_sqrt(simd4f a)
{
#if SIMD_SSE
	return _mm_sqrt_ps(a);
#elif SIMD_NEON
	return vsqrtq_f32(a);
#else
	return (simd4f){{ sqrtf(a.v[0]), sqrtf(a.v[1]), sqrtf(a.v[2]), sqrtf(a.v[3]) }};
#endif
}

// Divide first lane only.
static inline simd4f simd_div_ss(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_div_ss(a, b);
#elif SIMD_NEON
	float d = vgetq_lane_f32(a, 0) / vgetq_lane_f32(b, 0);
	return vsetq_lane_f32(d, a, 0);
#else
	simd4f r = a; r.v[0] = a.v[0] / b.v[0]; return r;
#endif
}

// Set first lane, zero rest.
static inline simd4f simd_set_ss(float x)
{
#if SIMD_SSE
	return _mm_set_ss(x);
#elif SIMD_NEON
	return vsetq_lane_f32(x, vdupq_n_f32(0), 0);
#else
	return (simd4f){{ x, 0, 0, 0 }};
#endif
}

// -----------------------------------------------------------------------------
// Bitwise operations.

static inline simd4f simd_and(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_and_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int ai, bi, ri; memcpy(&ai, &a.v[i], 4); memcpy(&bi, &b.v[i], 4); ri = ai & bi; memcpy(&r.v[i], &ri, 4); }
	return r;
#endif
}

static inline simd4f simd_andnot(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_andnot_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(b), vreinterpretq_u32_f32(a)));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int ai, bi, ri; memcpy(&ai, &a.v[i], 4); memcpy(&bi, &b.v[i], 4); ri = ~ai & bi; memcpy(&r.v[i], &ri, 4); }
	return r;
#endif
}

static inline simd4f simd_or(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_or_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vorrq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int ai, bi, ri; memcpy(&ai, &a.v[i], 4); memcpy(&bi, &b.v[i], 4); ri = ai | bi; memcpy(&r.v[i], &ri, 4); }
	return r;
#endif
}

static inline simd4f simd_xor(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_xor_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int ai, bi, ri; memcpy(&ai, &a.v[i], 4); memcpy(&bi, &b.v[i], 4); ri = ai ^ bi; memcpy(&r.v[i], &ri, 4); }
	return r;
#endif
}

// -----------------------------------------------------------------------------
// Comparisons (return mask: all-1s or all-0s per lane).

static inline simd4f simd_cmpge(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_cmpge_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vcgeq_f32(a, b));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int m = a.v[i] >= b.v[i] ? -1 : 0; memcpy(&r.v[i], &m, 4); }
	return r;
#endif
}

static inline simd4f simd_cmpgt(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_cmpgt_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vcgtq_f32(a, b));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int m = a.v[i] > b.v[i] ? -1 : 0; memcpy(&r.v[i], &m, 4); }
	return r;
#endif
}

static inline simd4f simd_cmple(simd4f a, simd4f b)
{
#if SIMD_SSE
	return _mm_cmple_ps(a, b);
#elif SIMD_NEON
	return vreinterpretq_f32_u32(vcleq_f32(a, b));
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int m = a.v[i] <= b.v[i] ? -1 : 0; memcpy(&r.v[i], &m, 4); }
	return r;
#endif
}

// Blendv: select a where mask MSB is set, b otherwise.
static inline simd4f simd_blendv(simd4f a, simd4f b, simd4f mask)
{
#if SIMD_SSE
	return _mm_blendv_ps(a, b, mask);
#elif SIMD_NEON
	uint32x4_t m = vreinterpretq_u32_f32(mask);
	return vbslq_f32(m, b, a);
#else
	simd4f r;
	for (int i = 0; i < 4; i++) { int mi; memcpy(&mi, &mask.v[i], 4); r.v[i] = mi < 0 ? b.v[i] : a.v[i]; }
	return r;
#endif
}

// Extract sign bits of each lane (4-bit mask).
static inline int simd_movemask(simd4f a)
{
#if SIMD_SSE
	return _mm_movemask_ps(a);
#elif SIMD_NEON
	uint32x4_t u = vreinterpretq_u32_f32(a);
	uint32x4_t shift = { 0, 1, 2, 3 };
	uint32x4_t bits = vshrq_n_u32(u, 31);
	bits = vshlq_u32(bits, vreinterpretq_s32_u32(shift));
	return (int)vaddvq_u32(bits);
#else
	int r = 0;
	for (int i = 0; i < 4; i++) { int mi; memcpy(&mi, &a.v[i], 4); if (mi < 0) r |= (1 << i); }
	return r;
#endif
}

// -----------------------------------------------------------------------------
// Integer operations (for index tracking in hull scans).

static inline simd4i simd_set1_i(int x)
{
#if SIMD_SSE
	return _mm_set1_epi32(x);
#elif SIMD_NEON
	return vdupq_n_s32(x);
#else
	return (simd4i){{ x, x, x, x }};
#endif
}

static inline simd4i simd_add_i(simd4i a, simd4i b)
{
#if SIMD_SSE
	return _mm_add_epi32(a, b);
#elif SIMD_NEON
	return vaddq_s32(a, b);
#else
	return (simd4i){{ a.v[0]+b.v[0], a.v[1]+b.v[1], a.v[2]+b.v[2], a.v[3]+b.v[3] }};
#endif
}

// Reinterpret between float and int SIMD types.
static inline simd4f simd_cast_itof(simd4i a)
{
#if SIMD_SSE
	return _mm_castsi128_ps(a);
#elif SIMD_NEON
	return vreinterpretq_f32_s32(a);
#else
	simd4f r; memcpy(&r, &a, sizeof(r)); return r;
#endif
}

static inline simd4i simd_cast_ftoi(simd4f a)
{
#if SIMD_SSE
	return _mm_castps_si128(a);
#elif SIMD_NEON
	return vreinterpretq_s32_f32(a);
#else
	simd4i r; memcpy(&r, &a, sizeof(r)); return r;
#endif
}

// Load/store.
static inline simd4f simd_load(const float* p)
{
#if SIMD_SSE
	return _mm_loadu_ps(p);
#elif SIMD_NEON
	return vld1q_f32(p);
#else
	return (simd4f){{ p[0], p[1], p[2], p[3] }};
#endif
}

static inline void simd_store(float* p, simd4f v)
{
#if SIMD_SSE
	_mm_storeu_ps(p, v);
#elif SIMD_NEON
	vst1q_f32(p, v);
#else
	for (int i = 0; i < 4; i++) p[i] = v.v[i];
#endif
}

static inline void simd_store_i(int* p, simd4i v)
{
#if SIMD_SSE
	_mm_storeu_si128((__m128i*)p, v);
#elif SIMD_NEON
	vst1q_s32(p, v);
#else
	for (int i = 0; i < 4; i++) p[i] = v.v[i];
#endif
}

// Shuffle helpers (SSE-style constants).
// For cross-platform, provide the most common patterns as named functions.
#if SIMD_SSE
	#define simd_shuffle(a, b, imm) _mm_shuffle_ps(a, b, imm)
	#define simd_unpacklo(a, b)     _mm_unpacklo_ps(a, b)
	#define simd_unpackhi(a, b)     _mm_unpackhi_ps(a, b)
	#define simd_movelh(a, b)       _mm_movelh_ps(a, b)
	#define simd_movehl(a, b)       _mm_movehl_ps(a, b)
#elif SIMD_NEON
	// NEON shuffle emulation for the patterns used in v3 ops.
	static inline simd4f simd_unpacklo(simd4f a, simd4f b) { return vzip1q_f32(a, b); }
	static inline simd4f simd_unpackhi(simd4f a, simd4f b) { return vzip2q_f32(a, b); }
	static inline simd4f simd_movelh(simd4f a, simd4f b) {
		return vcombine_f32(vget_low_f32(a), vget_low_f32(b));
	}
	static inline simd4f simd_movehl(simd4f a, simd4f b) {
		return vcombine_f32(vget_high_f32(b), vget_high_f32(a));
	}
	// General shuffle: use a lookup table of common patterns.
	// For the specific shuffles used in vmath (0x00, 0x55, 0xAA for broadcasts):
	#define simd_shuffle(a, b, imm) _simd_neon_shuffle(a, b, imm)
	static inline simd4f _simd_neon_shuffle(simd4f a, simd4f b, int imm) {
		// Only handle the patterns actually used: broadcast x/y/z and cross product shuffle.
		(void)b;
		switch (imm) {
		case 0x00: return vdupq_laneq_f32(a, 0);
		case 0x55: return vdupq_laneq_f32(a, 1);
		case 0xAA: return vdupq_laneq_f32(a, 2);
		case 0xC9: { // _MM_SHUFFLE(3,0,2,1) = yzxw
			float v[4] = { vgetq_lane_f32(a,1), vgetq_lane_f32(a,2), vgetq_lane_f32(a,0), vgetq_lane_f32(a,3) };
			return vld1q_f32(v);
		}
		default: return a; // fallback
		}
	}
#else
	#define simd_shuffle(a, b, imm) _simd_scalar_shuffle(a, b, imm)
	static inline simd4f _simd_scalar_shuffle(simd4f a, simd4f b, int imm) {
		return (simd4f){{ a.v[imm & 3], a.v[(imm >> 2) & 3], b.v[(imm >> 4) & 3], b.v[(imm >> 6) & 3] }};
	}
	static inline simd4f simd_unpacklo(simd4f a, simd4f b) { return (simd4f){{ a.v[0], b.v[0], a.v[1], b.v[1] }}; }
	static inline simd4f simd_unpackhi(simd4f a, simd4f b) { return (simd4f){{ a.v[2], b.v[2], a.v[3], b.v[3] }}; }
	static inline simd4f simd_movelh(simd4f a, simd4f b) { return (simd4f){{ a.v[0], a.v[1], b.v[0], b.v[1] }}; }
	static inline simd4f simd_movehl(simd4f a, simd4f b) { return (simd4f){{ b.v[2], b.v[3], a.v[2], a.v[3] }}; }
#endif

// Sign mask constant.
static inline simd4f simd_sign_mask() { return simd_cast_itof(simd_set1_i((int)0x80000000)); }

#endif // SIMD_H
