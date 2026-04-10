// See LICENSE for licensing info.
#ifndef VMATH_H
#define VMATH_H

#include <math.h>
#include <float.h>

// -----------------------------------------------------------------------------
// Types.

// v3: 16-byte aligned vector (w=0) for native SSE operation.
#include <xmmintrin.h>
#include <smmintrin.h>

typedef union v3
{
	struct { float x, y, z, _w; };
	__m128 m;
} v3;

typedef struct quat
{ float x, y, z, w; } quat;
typedef struct mat4
{ float m[16]; } mat4; // column-major
typedef struct Transform
{ v3 position; quat rotation; } Transform;

// Constructors.
#define V3(vx, vy, vz) ((v3){ .m = _mm_set_ps(0, vz, vy, vx) })

// -----------------------------------------------------------------------------
// v3 implementation (SSE-backed).

static inline v3 v3_add(v3 a, v3 b) { return (v3){ .m = _mm_add_ps(a.m, b.m) }; }
static inline v3 v3_sub(v3 a, v3 b) { return (v3){ .m = _mm_sub_ps(a.m, b.m) }; }
static inline v3 v3_scale(v3 a, float s) { return (v3){ .m = _mm_mul_ps(a.m, _mm_set1_ps(s)) }; }
static inline float v3_dot(v3 a, v3 b) {
	__m128 m = _mm_mul_ps(a.m, b.m);
	__m128 s = _mm_add_ss(m, _mm_shuffle_ps(m, m, 1));
	return _mm_cvtss_f32(_mm_add_ss(s, _mm_shuffle_ps(m, m, 2)));
}
static inline v3 v3_cross(v3 a, v3 b) {
	__m128 a1 = _mm_shuffle_ps(a.m, a.m, _MM_SHUFFLE(3,0,2,1));
	__m128 b1 = _mm_shuffle_ps(b.m, b.m, _MM_SHUFFLE(3,0,2,1));
	__m128 r = _mm_sub_ps(_mm_mul_ps(a.m, b1), _mm_mul_ps(a1, b.m));
	return (v3){ .m = _mm_shuffle_ps(r, r, _MM_SHUFFLE(3,0,2,1)) };
}
static inline float v3_len2(v3 a) { return v3_dot(a, a); }
static inline float v3_len(v3 a) { return sqrtf(v3_len2(a)); }
static inline v3 v3_norm(v3 a) { float l = v3_len(a); return v3_scale(a, 1.0f/l); }
static inline v3 v3_neg(v3 a) { return (v3){ .m = _mm_sub_ps(_mm_setzero_ps(), a.m) }; }

// -----------------------------------------------------------------------------
// quat implementation.

static inline quat quat_identity() { return (quat){ 0, 0, 0, 1 }; }
static inline quat quat_inv(quat q) { return (quat){ -q.x, -q.y, -q.z, q.w }; }

static inline quat quat_mul(quat a, quat b)
{
	return (quat){
		a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
		a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
		a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
		a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
	};
}

static __forceinline v3 quat_rotate(quat q, v3 v)
{
	// Load quat as __m128 [x,y,z,w], mask w to get u=[x,y,z,0]
	__m128 qm = _mm_loadu_ps(&q.x);
	v3 u = { .m = _mm_and_ps(qm, _mm_castsi128_ps(_mm_set_epi32(0, ~0, ~0, ~0))) };
	float s = q.w;
	float uv2 = 2.0f * v3_dot(u, v);
	float ss_uu = s*s - v3_dot(u, u);
	v3 uxv = v3_cross(u, v);
	return (v3){ .m = _mm_add_ps(_mm_add_ps(_mm_mul_ps(u.m, _mm_set1_ps(uv2)), _mm_mul_ps(v.m, _mm_set1_ps(ss_uu))), _mm_mul_ps(uxv.m, _mm_set1_ps(2.0f * s))) };
}

// -----------------------------------------------------------------------------
// Transform implementation (rotation + translation).

static inline Transform xform_identity() { return (Transform){ V3(0,0,0), quat_identity() }; }

static inline Transform xform_mul(Transform a, Transform b)
{
	return (Transform){
		v3_add(a.position, quat_rotate(a.rotation, b.position)),
		quat_mul(a.rotation, b.rotation),
	};
}

static inline v3 xform_apply(Transform t, v3 p)
{
	return v3_add(t.position, quat_rotate(t.rotation, p));
}

static inline Transform xform_inv(Transform t)
{
	quat inv_r = quat_inv(t.rotation);
	return (Transform){ quat_rotate(inv_r, v3_neg(t.position)), inv_r };
}

static inline v3 xform_inv_apply(Transform t, v3 p)
{
	return quat_rotate(quat_inv(t.rotation), v3_sub(p, t.position));
}

// -----------------------------------------------------------------------------
// mat4 implementation (column-major, m[col*4 + row]).

static inline mat4 mat4_identity()
{
	mat4 r = {0};
	r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0f;
	return r;
}

static inline mat4 mat4_mul(mat4 a, mat4 b)
{
	mat4 r = {0};
	for (int c = 0; c < 4; c++)
		for (int row = 0; row < 4; row++)
			for (int k = 0; k < 4; k++)
				r.m[c*4 + row] += a.m[k*4 + row] * b.m[c*4 + k];
	return r;
}

static inline mat4 mat4_perspective(float fovy, float aspect, float znear, float zfar)
{
	float f = 1.0f / tanf(fovy * 0.5f);
	float nf = 1.0f / (znear - zfar);
	mat4 r = {0};
	r.m[0]  = f / aspect;
	r.m[5]  = f;
	r.m[10] = (zfar + znear) * nf;
	r.m[11] = -1.0f;
	r.m[14] = 2.0f * zfar * znear * nf;
	return r;
}

static inline mat4 mat4_ortho(float l, float r, float b, float t, float n, float f)
{
	mat4 m = {0};
	m.m[0]  = 2.0f / (r - l);
	m.m[5]  = 2.0f / (t - b);
	m.m[10] = -2.0f / (f - n);
	m.m[12] = -(r + l) / (r - l);
	m.m[13] = -(t + b) / (t - b);
	m.m[14] = -(f + n) / (f - n);
	m.m[15] = 1.0f;
	return m;
}

static inline mat4 mat4_look_at(v3 eye, v3 center, v3 up)
{
	v3 f = v3_norm(v3_sub(center, eye));
	v3 r = v3_norm(v3_cross(f, up));
	v3 u = v3_cross(r, f);
	mat4 m = {0};
	m.m[0]  = r.x;   m.m[4]  = r.y;   m.m[8]  = r.z;   m.m[12] = -v3_dot(r, eye);
	m.m[1]  = u.x;   m.m[5]  = u.y;   m.m[9]  = u.z;   m.m[13] = -v3_dot(u, eye);
	m.m[2]  = -f.x;  m.m[6]  = -f.y;  m.m[10] = -f.z;  m.m[14] = v3_dot(f, eye);
	m.m[15] = 1.0f;
	return m;
}

static inline mat4 mat4_translate(v3 t)
{
	mat4 r = mat4_identity();
	r.m[12] = t.x; r.m[13] = t.y; r.m[14] = t.z;
	return r;
}

static inline mat4 mat4_scale_v3(v3 s)
{
	mat4 r = {0};
	r.m[0] = s.x; r.m[5] = s.y; r.m[10] = s.z; r.m[15] = 1.0f;
	return r;
}

static inline mat4 mat4_from_quat(quat q)
{
	float xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
	float xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
	float wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;
	mat4 r = {0};
	r.m[0]  = 1 - 2*(yy+zz); r.m[1]  = 2*(xy+wz);     r.m[2]  = 2*(xz-wy);
	r.m[4]  = 2*(xy-wz);     r.m[5]  = 1 - 2*(xx+zz); r.m[6]  = 2*(yz+wx);
	r.m[8]  = 2*(xz+wy);     r.m[9]  = 2*(yz-wx);     r.m[10] = 1 - 2*(xx+yy);
	r.m[15] = 1.0f;
	return r;
}

static inline mat4 mat4_trs(v3 pos, quat rot, v3 s)
{
	mat4 t = mat4_translate(pos);
	mat4 r = mat4_from_quat(rot);
	mat4 sc = mat4_scale_v3(s);
	return mat4_mul(t, mat4_mul(r, sc));
}

// -----------------------------------------------------------------------------
// v3 component-wise operations.

static inline v3 v3_mul(v3 a, v3 b) { return (v3){ .m = _mm_mul_ps(a.m, b.m) }; }
static inline v3 v3_min(v3 a, v3 b) { return (v3){ .m = _mm_min_ps(a.m, b.m) }; }
static inline v3 v3_max(v3 a, v3 b) { return (v3){ .m = _mm_max_ps(a.m, b.m) }; }
static inline v3 v3_rcp(v3 a) {
	// Safe reciprocal: zero where input is zero
	__m128 zero = _mm_setzero_ps();
	__m128 mask = _mm_cmpneq_ps(a.m, zero);
	__m128 r = _mm_div_ps(_mm_set1_ps(1.0f), _mm_blendv_ps(_mm_set1_ps(1.0f), a.m, mask));
	return (v3){ .m = _mm_and_ps(r, mask) };
}

// -----------------------------------------------------------------------------
// m3x3: row-major 3x3 matrix.

typedef struct m3x3
{ float m[9]; } m3x3;

static inline m3x3 diag(v3 d)
{
	return (m3x3){{ d.x, 0, 0,  0, d.y, 0,  0, 0, d.z }};
}

static inline m3x3 skew(v3 v)
{
	return (m3x3){{ 0, -v.z, v.y,  v.z, 0, -v.x,  -v.y, v.x, 0 }};
}

static inline m3x3 m3x3_add(m3x3 a, m3x3 b)
{
	m3x3 r;
	for (int i = 0; i < 9; i++) r.m[i] = a.m[i] + b.m[i];
	return r;
}

static inline m3x3 m3x3_sub(m3x3 a, m3x3 b)
{
	m3x3 r;
	for (int i = 0; i < 9; i++) r.m[i] = a.m[i] - b.m[i];
	return r;
}

static inline m3x3 m3x3_scale(m3x3 a, float s)
{
	m3x3 r;
	for (int i = 0; i < 9; i++) r.m[i] = a.m[i] * s;
	return r;
}

// Row-major: m3x3_mul_v3 = M * v
static inline v3 m3x3_mul_v3(m3x3 a, v3 v)
{
	return (v3){
		a.m[0]*v.x + a.m[1]*v.y + a.m[2]*v.z,
		a.m[3]*v.x + a.m[4]*v.y + a.m[5]*v.z,
		a.m[6]*v.x + a.m[7]*v.y + a.m[8]*v.z };
}

static inline m3x3 m3x3_mul(m3x3 a, m3x3 b)
{
	m3x3 r;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++) {
			r.m[i*3+j] = 0;
			for (int k = 0; k < 3; k++)
				r.m[i*3+j] += a.m[i*3+k] * b.m[k*3+j];
		}
	return r;
}

static inline m3x3 m3x3_identity(void) { return (m3x3){{ 1,0,0, 0,1,0, 0,0,1 }}; }

static inline m3x3 m3x3_neg(m3x3 a)
{
	m3x3 r;
	for (int i = 0; i < 9; i++) r.m[i] = -a.m[i];
	return r;
}

static inline m3x3 m3x3_transpose(m3x3 a)
{
	return (m3x3){{ a.m[0], a.m[3], a.m[6],  a.m[1], a.m[4], a.m[7],  a.m[2], a.m[5], a.m[8] }};
}

// Solve A*x = b via Cramer's rule. Returns zero vector if singular.
static inline v3 solve(m3x3 a, v3 b)
{
	float* m = a.m;
	float det = m[0]*(m[4]*m[8] - m[5]*m[7])
	          - m[1]*(m[3]*m[8] - m[5]*m[6])
	          + m[2]*(m[3]*m[7] - m[4]*m[6]);
	if (fabsf(det) < 1e-12f) return V3(0,0,0);
	float d = 1.0f / det;
	return (v3){
		((m[4]*m[8]-m[5]*m[7])*b.x - (m[1]*m[8]-m[2]*m[7])*b.y + (m[1]*m[5]-m[2]*m[4])*b.z) * d,
		(-(m[3]*m[8]-m[5]*m[6])*b.x + (m[0]*m[8]-m[2]*m[6])*b.y - (m[0]*m[5]-m[2]*m[3])*b.z) * d,
		((m[3]*m[7]-m[4]*m[6])*b.x - (m[0]*m[7]-m[1]*m[6])*b.y + (m[0]*m[4]-m[1]*m[3])*b.z) * d };
}

// -----------------------------------------------------------------------------
// _Generic polymorphic math API.
//
// add(a, b)     -- v3+v3
// sub(a, b)     -- v3-v3
// mul(a, b)     -- mat4*mat4, Transform*Transform, Transform*v3, quat*quat, quat*v3
// scale(a, s)   -- v3*float
// dot(a, b)     -- v3.v3
// cross(a, b)   -- v3xv3
// neg(a)        -- -v3
// norm(a)       -- normalize v3
// len(a)        -- length v3
// len2(a)       -- squared length v3
// inv(a)        -- inverse quat or Transform

#define add(a, b)    _Generic((a), v3: v3_add, m3x3: m3x3_add)(a, b)
#define sub(a, b)    _Generic((a), v3: v3_sub, m3x3: m3x3_sub)(a, b)
#define scale(a, s)  _Generic((a), v3: v3_scale, m3x3: m3x3_scale)(a, s)
#define dot(a, b)    _Generic((a), v3: v3_dot)(a, b)
#define cross(a, b)  _Generic((a), v3: v3_cross)(a, b)
#define neg(a)       _Generic((a), v3: v3_neg, m3x3: m3x3_neg)(a)
#define transpose(a) _Generic((a), m3x3: m3x3_transpose)(a)
#define hmul(a, b)   _Generic((a), v3: v3_mul)(a, b)
#define rcp(a)       _Generic((a), v3: v3_rcp)(a)
#define norm(a)      _Generic((a), v3: v3_norm)(a)
#define len(a)       _Generic((a), v3: v3_len)(a)
#define len2(a)      _Generic((a), v3: v3_len2)(a)

#define inv(a)       _Generic((a), quat: quat_inv, Transform: xform_inv)(a)

// mul(a, b): same-type multiply (mat4*mat4, Transform*Transform, quat*quat).
// For mixed-type transforms use: rotate(quat, v3), xform(Transform, v3).
#define mul(a, b)    _Generic((a), mat4: mat4_mul, Transform: xform_mul, quat: quat_mul, m3x3: m3x3_mul)(a, b)
#define rotate(q, v) quat_rotate(q, v)
#define xform(t, p)  xform_apply(t, p)

// ---------------------------------------------------------------------------
// Generic NxN block LDL^T, single precision, packed lower-triangular storage.
// Max n = 6 (21 packed floats). Suitable for solving coupled multi-DOF
// constraints (joints, contact blocks, preconditioner blocks).

#define BLOCK_MAX_DOF 16 // max DOFs per block (contacts can exceed 6 per body pair)

// Packed lower-triangular index. Symmetric: BTRI(r,c) == BTRI(c,r).
#define BTRI(r, c) ((r) >= (c) ? (r)*((r)+1)/2 + (c) : (c)*((c)+1)/2 + (r))

// In-place LDL^T factorization. A is packed lower-tri (overwritten with L,
// unit diagonal implicit). D is separate diagonal. Returns 0 on success,
// -1 if a non-positive pivot is encountered. n up to BLOCK_MAX_DOF.
static inline int block_ldl_f(float* A, float* D, int n)
{
	for (int j = 0; j < n; j++) {
		float dj = A[BTRI(j,j)];
		for (int k = 0; k < j; k++)
			dj -= A[BTRI(j,k)] * A[BTRI(j,k)] * D[k];
		if (dj <= 1e-12f) return -1;
		D[j] = dj;
		float inv_dj = 1.0f / dj;
		for (int i = j + 1; i < n; i++) {
			float lij = A[BTRI(i,j)];
			for (int k = 0; k < j; k++)
				lij -= A[BTRI(i,k)] * A[BTRI(j,k)] * D[k];
			A[BTRI(i,j)] = lij * inv_dj;
		}
	}
	return 0;
}

// Solve LDL^T x = b. L and D from block_ldl_f. b is not modified.
static inline void block_solve_f(const float* L, const float* D,
                                  const float* b, float* x, int n)
{
	// Forward: L y = b
	for (int i = 0; i < n; i++) {
		x[i] = b[i];
		for (int k = 0; k < i; k++) x[i] -= L[BTRI(i,k)] * x[k];
	}
	// Diagonal: D z = y
	for (int i = 0; i < n; i++) x[i] /= D[i];
	// Backward: L^T x = z
	for (int i = n - 1; i >= 0; i--)
		for (int k = i + 1; k < n; k++)
			x[i] -= L[BTRI(k,i)] * x[k];
}

#endif
