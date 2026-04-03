// nudge.h -- single-file 3D physics library (generated)
// https://github.com/RandyGaul/nudge
//
// Do this in *one* C file:
//   #define NUDGE_IMPLEMENTATION
//   #include "nudge.h"

#ifndef NUDGE_SINGLE_FILE_H
#define NUDGE_SINGLE_FILE_H

#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include <assert.h>


// ---- src/nudge.h ----

// See LICENSE for licensing info.
//
// Credits and references
// ----------------------
// BEPUphysics v2 -- broadphase, soft constraints, memory layout ideas
// Box2D -- GJK, sequential impulse solver reference
// Dirk Gregorius -- SAT with Gauss map pruning, contact clipping, direct enumeration LCP block solver
// Chris Giles -- Augmented Vertex Block Descent solver reference

#ifndef NUDGE_H
#define NUDGE_H


// ---- src/vmath.h ----

// See LICENSE for licensing info.
#ifndef VMATH_H
#define VMATH_H

#include <math.h>
#include <float.h>

// -----------------------------------------------------------------------------
// Types.

typedef struct v3
{ float x, y, z; } v3;
typedef struct quat
{ float x, y, z, w; } quat;
typedef struct mat4
{ float m[16]; } mat4; // column-major
typedef struct Transform
{ v3 position; quat rotation; } Transform;

// Constructors.
#define V3(x, y, z) ((v3){ x, y, z })

// -----------------------------------------------------------------------------
// v3 implementation.

static inline v3 v3_add(v3 a, v3 b) { return (v3){ a.x+b.x, a.y+b.y, a.z+b.z }; }
static inline v3 v3_sub(v3 a, v3 b) { return (v3){ a.x-b.x, a.y-b.y, a.z-b.z }; }
static inline v3 v3_scale(v3 a, float s) { return (v3){ a.x*s, a.y*s, a.z*s }; }
static inline float v3_dot(v3 a, v3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline v3 v3_cross(v3 a, v3 b) { return (v3){ a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; }
static inline float v3_len2(v3 a) { return v3_dot(a, a); }
static inline float v3_len(v3 a) { return sqrtf(v3_len2(a)); }
static inline v3 v3_norm(v3 a) { float l = v3_len(a); return v3_scale(a, 1.0f/l); }
static inline v3 v3_neg(v3 a) { return (v3){ -a.x, -a.y, -a.z }; }

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

static inline v3 quat_rotate(quat q, v3 v)
{
	v3 u = { q.x, q.y, q.z };
	float s = q.w;
	return v3_add(v3_add(v3_scale(u, 2.0f * v3_dot(u, v)), v3_scale(v, s*s - v3_dot(u, u))), v3_scale(v3_cross(u, v), 2.0f * s));
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

static inline v3 v3_mul(v3 a, v3 b) { return (v3){ a.x*b.x, a.y*b.y, a.z*b.z }; }
static inline v3 v3_min(v3 a, v3 b) { return (v3){ fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z) }; }
static inline v3 v3_max(v3 a, v3 b) { return (v3){ fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z) }; }
static inline v3 v3_rcp(v3 a) {
	return (v3){
		a.x != 0.0f ? 1.0f / a.x : 0.0f,
		a.y != 0.0f ? 1.0f / a.y : 0.0f,
		a.z != 0.0f ? 1.0f / a.z : 0.0f };
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

// -----------------------------------------------------------------------------
// AVBD helpers: 6x6 LDL solve, quaternion small-angle ops.

// Small angular update: q + 0.5*(0,dw)*q, normalized.
static inline quat quat_add_dw(quat q, v3 dw)
{
	quat spin = { dw.x, dw.y, dw.z, 0.0f };
	quat dq = quat_mul(spin, q);
	quat r = { q.x + 0.5f*dq.x, q.y + 0.5f*dq.y, q.z + 0.5f*dq.z, q.w + 0.5f*dq.w };
	float l = sqrtf(r.x*r.x + r.y*r.y + r.z*r.z + r.w*r.w);
	if (l < 1e-15f) return q;
	float s = 1.0f / l;
	return (quat){ r.x*s, r.y*s, r.z*s, r.w*s };
}

// Angular velocity from quaternion difference: (a * inv(b)).xyz * 2
static inline v3 quat_sub_angular(quat a, quat b)
{
	quat d = quat_mul(a, quat_inv(b));
	return (v3){ d.x * 2.0f, d.y * 2.0f, d.z * 2.0f };
}

// Solve 6x6 SPD system via unrolled LDL^T decomposition.
// System is [lhs_lin, lhs_cross^T; lhs_cross, lhs_ang] * [dx_lin; dx_ang] = [rhs_lin; rhs_ang]
// Reference: Chris Giles, AVBD 3D demo (maths.h).
static inline void solve_6x6_ldl(m3x3 aLin, m3x3 aAng, m3x3 aCross,
                                  v3 bLin, v3 bAng, v3* xLin, v3* xAng)
{
	// m3x3 is row-major float[9]: m[row*3+col]
	float A11 = aLin.m[0];
	float A21 = aLin.m[3], A22 = aLin.m[4];
	float A31 = aLin.m[6], A32 = aLin.m[7], A33 = aLin.m[8];
	float A41 = aCross.m[0], A42 = aCross.m[1], A43 = aCross.m[2], A44 = aAng.m[0];
	float A51 = aCross.m[3], A52 = aCross.m[4], A53 = aCross.m[5], A54 = aAng.m[3], A55 = aAng.m[4];
	float A61 = aCross.m[6], A62 = aCross.m[7], A63 = aCross.m[8], A64 = aAng.m[6], A65 = aAng.m[7], A66 = aAng.m[8];

	// LDL^T factorization
	float D1 = A11;
	if (D1 == 0.0f) D1 = 1e-12f;
	float L21 = A21 / D1, L31 = A31 / D1, L41 = A41 / D1, L51 = A51 / D1, L61 = A61 / D1;

	float D2 = A22 - L21*L21*D1;
	if (D2 == 0.0f) D2 = 1e-12f;
	float L32 = (A32 - L31*L21*D1) / D2;
	float L42 = (A42 - L41*L21*D1) / D2;
	float L52 = (A52 - L51*L21*D1) / D2;
	float L62 = (A62 - L61*L21*D1) / D2;

	float D3 = A33 - (L31*L31*D1 + L32*L32*D2);
	if (D3 == 0.0f) D3 = 1e-12f;
	float L43 = (A43 - L41*L31*D1 - L42*L32*D2) / D3;
	float L53 = (A53 - L51*L31*D1 - L52*L32*D2) / D3;
	float L63 = (A63 - L61*L31*D1 - L62*L32*D2) / D3;

	float D4 = A44 - (L41*L41*D1 + L42*L42*D2 + L43*L43*D3);
	if (D4 == 0.0f) D4 = 1e-12f;
	float L54 = (A54 - L41*L51*D1 - L42*L52*D2 - L43*L53*D3) / D4;
	float L64 = (A64 - L41*L61*D1 - L42*L62*D2 - L43*L63*D3) / D4;

	float D5 = A55 - (L51*L51*D1 + L52*L52*D2 + L53*L53*D3 + L54*L54*D4);
	if (D5 == 0.0f) D5 = 1e-12f;
	float L65 = (A65 - L51*L61*D1 - L52*L62*D2 - L53*L63*D3 - L54*L64*D4) / D5;

	float D6 = A66 - (L61*L61*D1 + L62*L62*D2 + L63*L63*D3 + L64*L64*D4 + L65*L65*D5);
	if (D6 == 0.0f) D6 = 1e-12f;

	// Forward substitution: Ly = b
	float y1 = bLin.x;
	float y2 = bLin.y - L21*y1;
	float y3 = bLin.z - L31*y1 - L32*y2;
	float y4 = bAng.x - L41*y1 - L42*y2 - L43*y3;
	float y5 = bAng.y - L51*y1 - L52*y2 - L53*y3 - L54*y4;
	float y6 = bAng.z - L61*y1 - L62*y2 - L63*y3 - L64*y4 - L65*y5;

	// Diagonal solve: Dz = y
	float z1 = y1/D1, z2 = y2/D2, z3 = y3/D3, z4 = y4/D4, z5 = y5/D5, z6 = y6/D6;

	// Backward substitution: L^T x = z
	xAng->z = z6;
	xAng->y = z5 - L65*xAng->z;
	xAng->x = z4 - L54*xAng->y - L64*xAng->z;
	xLin->z = z3 - L43*xAng->x - L53*xAng->y - L63*xAng->z;
	xLin->y = z2 - L32*xLin->z - L42*xAng->x - L52*xAng->y - L62*xAng->z;
	xLin->x = z1 - L21*xLin->y - L31*xLin->z - L41*xAng->x - L51*xAng->y - L61*xAng->z;
}

#endif

// ---- src/split_store.h ----

// See LICENSE for licensing info.
#ifndef SPLIT_STORE_H
#define SPLIT_STORE_H

// split_store.h -- Paired hot/cold storage with stable generational handles.
//
// PURPOSE: Cache-friendly data layout for tight iteration loops (solvers,
// broadphase, rendering). Inspired by Box2D's hot/cold body split.
//
// The idea: a fat struct with everything (position, velocity, mass, shapes,
// user data, topology links) destroys cache when the solver only needs
// position + velocity. Instead, split into two parallel arrays:
//
//   hot[]  -- only the fields the inner loop touches every frame (position,
//             velocity, inv_mass). Small struct, contiguous, cache-line friendly.
//   cold[] -- everything else (mass, shapes, topology, user data). Same indices
//             as hot[], but only accessed on-demand (creation, queries, etc).
//
// Both arrays share indices, so hot[i] and cold[i] describe the same object.
// A separate gen[] array tracks slot liveness (odd = alive, even = dead) and
// enables stable generational handles without polluting either data struct.
// A freelist[] recycles slots so indices remain stable across add/remove.
//
// CACHE OPTIMIZATION GUIDE:
//   - Put frequently-iterated solver fields in Hot (velocity, position, forces).
//   - Put rarely-read metadata in Cold (mass, shapes, user data, topology).
//   - Keep Hot structs small -- ideally <= 64 bytes (one cache line).
//   - The solver loop iterates hot[] linearly. Prefetch is automatic.
//   - Cold is only dereferenced when you need it (handle lookup, shape query).
//   - Gaps from deleted slots are zeroed and skipped via split_alive().
//     The branch is predictable and avoids the indirection-fixup cost of
//     swap-remove. Compact later if occupancy drops.
//
// Usage:
//     CK_DYNA MyCold*    cold = NULL;
//     CK_DYNA MyHot*     hot  = NULL;
//     CK_DYNA uint32_t*  gen  = NULL;
//     CK_DYNA int*       free_list = NULL;
//
//     int idx;
//     split_add(cold, hot, gen, free_list, idx);
//     cold[idx] = (MyCold){ ... };
//     hot[idx]  = (MyHot){ ... };
//     MyHandle h = split_handle(MyHandle, gen, idx);
//
//     assert(split_valid(gen, h));
//
//     for (int i = 0; i < asize(hot); i++) {
//         if (!split_alive(gen, i)) continue;
//         // iterate hot[i] ...
//     }
//
//     split_del(cold, hot, gen, free_list, idx);
//     split_free(cold, hot, gen, free_list);

// Allocate a slot. Reuses from freelist or grows all arrays.
// Sets out_idx to the allocated index. gen[out_idx] becomes odd (alive).
#define split_add(cold, hot, gen, freelist, out_idx) do { \
	if (asize(freelist) > 0) { \
		(out_idx) = apop(freelist); \
	} else { \
		(out_idx) = asize(cold); \
		afit((cold), (out_idx) + 1); \
		afit((hot), (out_idx) + 1); \
		afit((gen), (out_idx) + 1); \
		asetlen((cold), (out_idx) + 1); \
		asetlen((hot), (out_idx) + 1); \
		asetlen((gen), (out_idx) + 1); \
		(gen)[(out_idx)] = 0; \
	} \
	memset(&(cold)[(out_idx)], 0, sizeof(*(cold))); \
	memset(&(hot)[(out_idx)], 0, sizeof(*(hot))); \
	(gen)[(out_idx)]++; \
} while(0)

// Free a slot. Zeros cold+hot, bumps gen to even (dead), pushes to freelist.
#define split_del(cold, hot, gen, freelist, i) do { \
	memset(&(cold)[(i)], 0, sizeof(*(cold))); \
	memset(&(hot)[(i)], 0, sizeof(*(hot))); \
	(gen)[(i)]++; \
	apush((freelist), (i)); \
} while(0)

// Check if slot i is alive (odd generation).
#define split_alive(gen, i) ((gen)[(i)] & 1)

// Free all arrays.
#define split_free(cold, hot, gen, freelist) do { \
	afree(cold); afree(hot); afree(gen); afree(freelist); \
} while(0)

// Handle encoding: low 32 bits = index, high 32 bits = generation.
#define handle_make(idx, gen)  ((uint64_t)(gen) << 32 | (uint32_t)(idx))
#define handle_index(h)        ((int)((h).id & 0xFFFFFFFF))
#define handle_gen(h)          ((uint32_t)((h).id >> 32))

// Build a handle from a split store's gen array and an index.
#define split_handle(T, gen, idx) ((T){ handle_make((idx), (gen)[(idx)]) })

// Validate a handle against the current gen array.
#define split_valid(gen, h) ((gen)[handle_index(h)] == handle_gen(h))

#endif

// -----------------------------------------------------------------------------
// Opaque handles.

typedef struct World { uint64_t id; } World;
typedef struct Body { uint64_t id; } Body;

// -----------------------------------------------------------------------------
// Shape types -- usable standalone for direct collision queries.

typedef struct Sphere
{
	v3 center;
	float radius;
} Sphere;

typedef struct Capsule
{
	v3 p, q;       // segment endpoints (world space)
	float radius;
} Capsule;

typedef struct Box
{
	v3 center;
	quat rotation;
	v3 half_extents;
} Box;

// Half-edge mesh for convex polyhedra.
// Edges stored in twin pairs: edge 2k and twin 2k+1.
typedef struct HalfEdge
{
	uint16_t twin;
	uint16_t next;
	uint16_t origin;
	uint16_t face;
} HalfEdge;

typedef struct HullPlane
{
	v3 normal;
	float offset;    // dot(normal, point_on_plane)
} HullPlane;

typedef struct HullFace
{
	uint16_t edge;    // first half-edge on this face
} HullFace;

// Convex hull -- half-edge mesh with precomputed face planes.
// Can point to static data (e.g. unit box) or dynamically built via quickhull.
typedef struct Hull
{
	v3 centroid;
	const v3*        verts;
	const HalfEdge*  edges;
	const HullFace*  faces;
	const HullPlane* planes;
	int vert_count;
	int edge_count;  // total half-edges (2x undirected edges)
	int face_count;
	float epsilon;      // build tolerance: 3*(max|x|+max|y|+max|z|)*FLT_EPSILON
	float maxoutside;   // max distance any vertex was widened beyond Newell plane
} Hull;

// Positioned hull for collision queries.
typedef struct ConvexHull
{
	const Hull* hull;
	v3 center;
	quat rotation;
	v3 scale;        // per-axis scale applied to unit hull verts
} ConvexHull;

// -----------------------------------------------------------------------------
// Contact manifold.

// Contact feature ID for warm starting.
// Encodes which geometric features (face/edge/vertex) produced the contact,
// enabling frame-to-frame matching by integer comparison (Box2D/BEPU style).
//
// Face contacts: ref_face | (inc_face << 8) | (clip_edge << 16)
//   ref_face:  index of reference face (the clipping face)
//   inc_face:  index of incident face
//   clip_edge: which side plane produced this vertex (0xFF = original vertex)
//
// Edge contacts: edge_a | (edge_b << 16) | FEATURE_EDGE_BIT
#define FEATURE_EDGE_BIT 0x80000000u

typedef struct Contact
{
	v3 point;
	v3 normal;         // from shape A toward shape B
	float penetration; // positive = overlapping
	uint32_t feature_id;
} Contact;

#define MAX_CONTACTS 4

typedef struct Manifold
{
	Contact contacts[MAX_CONTACTS];
	int count;
} Manifold;

// -----------------------------------------------------------------------------
// Collision queries.
//
// Each returns nonzero on hit. Pass manifold=NULL for boolean-only test.
// Normal convention: points from A toward B.

int collide_sphere_sphere(Sphere a, Sphere b, Manifold* manifold);
int collide_sphere_capsule(Sphere a, Capsule b, Manifold* manifold);
int collide_sphere_box(Sphere a, Box b, Manifold* manifold);
int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold);
int collide_capsule_capsule(Capsule a, Capsule b, Manifold* manifold);
int collide_capsule_box(Capsule a, Box b, Manifold* manifold);
int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold);
int collide_box_box(Box a, Box b, Manifold* manifold);
int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold);

// Built-in unit box hull (half-extents 1,1,1). Use with ConvexHull + scale for boxes.
const Hull* hull_unit_box();

// -----------------------------------------------------------------------------
// Quickhull -- build a convex hull from a point cloud.
//
// Returns a heap-allocated Hull. Caller frees with hull_free().
// The resulting hull has proper half-edge topology, face planes, and centroid.

Hull* quickhull(const v3* points, int count);
void hull_free(Hull* hull);

// -----------------------------------------------------------------------------
// Body params for world API.

typedef enum ShapeType
{
	SHAPE_SPHERE,
	SHAPE_CAPSULE,
	SHAPE_BOX,
	SHAPE_HULL,
} ShapeType;

typedef struct ShapeParams
{
	ShapeType type;
	v3 local_pos;   // offset from body origin
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule; // segment along local Y
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
	};
} ShapeParams;

typedef struct BodyParams
{
	v3 position;
	quat rotation;
	float mass;            // 0 = static/kinematic
	float friction;        // Coulomb mu (default 0.5)
	float restitution;     // bounce coefficient (default 0.0)
	float linear_damping;  // velocity decay coefficient (default 0.0)
	float angular_damping; // angular velocity decay coefficient (default 0.03)
} BodyParams;

typedef enum BroadphaseType { BROADPHASE_N2, BROADPHASE_BVH } BroadphaseType;

typedef enum FrictionModel
{
	FRICTION_COULOMB,  // per-point Coulomb (2 tangent rows per contact)
	FRICTION_PATCH,    // manifold-level 2D friction using patch area estimate
} FrictionModel;

typedef enum SolverType
{
	SOLVER_SOFT_STEP,  // soft contacts, relax each substep (default)
	SOLVER_SI_SOFT,    // soft contacts, no relax between substeps
	SOLVER_SI,         // hard constraints, NGS position correction
	SOLVER_BLOCK,      // direct LCP enumeration for 2-4 contact normals
	SOLVER_AVBD,       // augmented vertex block descent (primal-dual position solver)
} SolverType;

typedef struct WorldParams
{
	v3 gravity;
	BroadphaseType broadphase;
	FrictionModel friction_model;
	SolverType solver_type;
	int velocity_iters;  // 0 = default (10)
	int position_iters;  // 0 = default (4)
	float contact_hertz;          // 0 = default (30.0), soft contact frequency
	float contact_damping_ratio;  // 0 = default (10.0), heavily overdamped
	float max_push_velocity;      // 0 = default (3.0 m/s)
	int sub_steps;                // 0 = default (1)
} WorldParams;

// -----------------------------------------------------------------------------
// World API.

World create_world(WorldParams params);
void destroy_world(World world);
void world_step(World world, float dt);
void world_set_friction_model(World world, FrictionModel model);
void world_set_solver_type(World world, SolverType type);

Body create_body(World world, BodyParams params);
void destroy_body(World world, Body body);
void body_add_shape(World world, Body body, ShapeParams params);

v3 body_get_position(World world, Body body);
quat body_get_rotation(World world, Body body);

void body_wake(World world, Body body);
void body_set_velocity(World world, Body body, v3 vel);
void body_set_angular_velocity(World world, Body body, v3 avel);
int body_is_asleep(World world, Body body);

// Debug: contact points from last step. Returns count, *out valid until next step.
int world_get_contacts(World world, const Contact** out);

// -----------------------------------------------------------------------------
// Joints.

typedef struct Joint { uint64_t id; } Joint;

typedef struct SpringParams
{
	float frequency;      // Hz (0 = rigid constraint)
	float damping_ratio;  // 1.0 = critically damped
} SpringParams;

typedef struct BallSocketParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	SpringParams spring;  // {0,0} = rigid
} BallSocketParams;

typedef struct DistanceParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	float rest_length;    // 0 = compute from initial body positions
	SpringParams spring;
} DistanceParams;

Joint create_ball_socket(World world, BallSocketParams params);
Joint create_distance(World world, DistanceParams params);
void destroy_joint(World world, Joint joint);

// Debug: iterate BVH nodes. Calls fn(min, max, depth, is_leaf, user) for each node child.
typedef void (*BVHDebugFn)(v3 min, v3 max, int depth, int is_leaf, void* user);
void world_debug_bvh(World world, BVHDebugFn fn, void* user);

#endif

#ifdef NUDGE_IMPLEMENTATION

// ---- src/nudge_amalg.c ----

// Implementation seed for single-file header amalgamation.
#define CKIT_IMPLEMENTATION

// ---- src/ckit.h ----

/*
    ckit -- A single-header kit of high-performance C essentials.

    In exactly one .c/.cpp file:
           #define CKIT_IMPLEMENTATION
       In all other cases include as-expected:

    QUICK START GUIDE

        1) Dynamic arrays:
            int* a = NULL;
            for (int i = 0; i < 10; ++i) {
                apush(a, i);
            }
            for (int i = 0; i < 10; ++i) {
                printf("%d\n", a[i]);
            }
            printf("len=%d cap=%d\n", acount(a), acap(a));
            afree(a);

        2) Dynamic strings:
            char* s = NULL;
            sset(s, "Hello ");
            sappend(s, "world!");
            printf("%s\n", s);
            sfree(s);

        3) Map (typed hashtable):
            CK_MAP(int) m = NULL;
            map_set(m, sintern("x"), 10);
            map_set(m, sintern("y"), 20);
            int x = map_get(m, sintern("x"));
            int y = map_get(m, sintern("y"));
            for (int i = 0; i < map_size(m); i++) {
                printf("key=%llu val=%d\n", map_keys(m)[i], map_items(m)[i]);
            }
            map_free(m);

        4) String interning:
           const char* a = sintern("hello");
           char buf[64];
           strcpy(buf, "he");
           strcat(buf, "llo");
           const char *b = sintern(buf);
           assert(a == b);

    Public domain - Do whatever you want with this code.
*/

#ifndef CKIT_H
#define CKIT_H

#include <stdint.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ----------------------------------------------------------------------------------------------------
// Overrideable macros (define before including ckit.h to customize behavior):
//
//   CK_ALLOC(sz)       - Custom allocator. Default: malloc(sz)
//   CK_REALLOC(p, sz)  - Custom reallocator. Default: realloc(p, sz)
//   CK_FREE(p)         - Custom free. Default: free(p)

#ifndef CK_ALLOC
#	define CK_ALLOC(sz) malloc(sz)
#endif
#ifndef CK_REALLOC
#	define CK_REALLOC(p, sz) realloc(p, sz)
#endif
#ifndef CK_FREE
#	define CK_FREE(p) free(p)
#endif

//--------------------------------------------------------------------------------------------------
// Dynamic arrays (stretchy buffers).
//
// Declare as T* arr = NULL. Push elements with apush(). Access like a normal array.
// NULL is a valid empty array.
//
// Example:
//     int* a = NULL;
//     apush(a, 10);
//     apush(a, 20);
//     printf("%d %d\n", a[0], a[1]);  // 10 20
//     afree(a);

// Markup macros for documentation purposes. Expand to nothing.
#define CK_DYNA   // Annotates a pointer as a dynamic array.
#define CK_SDYNA  // Annotates a pointer as a dynamic string.

// asize/acount: Number of elements. Returns 0 for NULL.
#define asize(a)      ((a) ? CK_AHDR(a)->size : 0)
#define acount(a)     asize(a)

// asetlen: Set size directly (must be <= capacity). Does not allocate.
#define asetlen(a, n) (CK_AHDR(a)->size = (n))

// acap: Allocated capacity. Returns 0 for NULL.
#define acap(a)       ((a) ? CK_AHDR(a)->capacity : 0)

// afit: Ensure capacity for n elements. May reallocate.
#define afit(a, n)    ((n) <= acap(a) ? 0 : (*(void**)&(a) = ck_agrow((a), (n), sizeof(*(a)))))

// afit_set: Allocate and set size to exactly n elements. Combines afit + asetlen.
#define afit_set(a, n) (afit((a), (n)), asetlen((a), (n)))

// apush: Append element. May reallocate.
#define apush(a, ...) (CK_ACANARY(a), afit((a), 1 + asize(a)), (a)[CK_AHDR(a)->size++] = (__VA_ARGS__))

// apop: Remove and return last element. Array must not be empty.
#define apop(a)       ((a)[--CK_AHDR(a)->size])

// aend: Pointer one past the last element.
#define aend(a)       ((a) + asize(a))

// alast: Last element. Array must not be empty.
#define alast(a)      ((a)[asize(a) - 1])

// aclear: Set size to 0 but keep allocated memory.
#define aclear(a)     (CK_ACANARY(a), (a) ? CK_AHDR(a)->size = 0 : 0)

// adel: Remove element at index i by swapping with last (O(1), unordered).
#define adel(a, i)    ((a)[i] = (a)[--CK_AHDR(a)->size])

// astatic: Use a stack/static buffer instead of heap allocation.
//     char buffer[1024];
//     int* a = NULL;
//     astatic(a, buffer, sizeof(buffer));
#define astatic(a, buffer, buffer_size) (*(void**)&(a) = ck_astatic(buffer, buffer_size, sizeof(*(a))))

// aset: Copy array b into a.
#define aset(a, b)    (*(void**)&(a) = ck_aset((void*)(a), (void*)(b), sizeof(*(a))))

// arev: Reverse array in-place.
#define arev(a)       ((a) ? ck_arev(a, sizeof(*(a))) : (void*)0)

// ahash: Hash all bytes in the array using FNV1a.
#define ahash(a)      ((a) ? ck_hash_fnv1a(a, sizeof(*(a)) * asize(a)) : 0)

// afree: Free array memory and set pointer to NULL.
#define afree(a)      do { CK_ACANARY(a); if (a && !CK_AHDR(a)->is_static) CK_FREE(CK_AHDR(a)); (a) = NULL; } while (0)

// Check if array is a valid dynamic array.
#define avalid(a)  ((a) && CK_AHDR(a)->cookie.val == CK_ACOOKIE)

//--------------------------------------------------------------------------------------------------
// Dynamic strings (built on dynamic arrays).
//
// 100% compatible with normal C-strings. Free with sfree when done.
//
// To create a new string use smake or sfmake. To overwrite an existing string
// use sset or sfmt. The first argument to sset/sfmt must be an l-value (a
// variable), not a literal like NULL -- use smake/sfmake to create from scratch.
//
// Example:
//
//     char* s = smake("Hello world!");
//     printf("%s\n", s);
//     sset(s, "Goodbye!");
//     printf("%s\n", s);
//     sfree(s);

// slen: String length (not including null terminator). Returns 0 for NULL.
#define slen(s)                 ((s) ? (asize(s) ? asize(s) - 1 : asize(s)) : 0)

// scount: Internal array size (length + null terminator).
#define scount(s)               asize(s)

// scap: Allocated capacity.
#define scap(s)                 acap(s)

// sempty: True if string is NULL or has length 0.
#define sempty(s)               ((s) ? slen(s) < 1 : 1)

// sfirst/slast: First or last character. Returns '\0' for empty/NULL.
#define sfirst(s)               ((s) ? (s)[0] : '\0')
#define slast(s)                ((s) ? (s)[slen(s) - 1] : '\0')

// sclear: Clear string to empty but keep allocated memory.
#define sclear(s)               (aclear(s), apush(s, 0))

// sstatic: Create a string backed by a static buffer. Grows to heap if needed.
#define sstatic(s, buf, sz)     (astatic(s, buf, sz), apush(s, 0))

// sisdyna: True if s is a dynamic string from this API (not a literal or static buffer).
#define sisdyna(s)              (!((#s)[0] == '"') && ck_avalid(s))

// sfree: Free string memory and set pointer to NULL.
#define sfree(s)                afree(s)

// sfit: Ensure capacity for n characters.
#define sfit(s, n)              (s = ck_sfit(s, n))

// spush/spop: Append or remove a single character.
#define spush(s, ch)            do { if (!(s)) apush(s, ch); else (s)[slen(s)] = (ch); apush(s, 0); } while (0)
#define spop(s)                 (s = ck_spop(s))
#define spopn(s, n)             (s = ck_spopn(s, n))

// sset: Overwrite a with contents of b.
#define sset(a, b)              (a = ck_sset(a, b))

// sdup/smake: Create a new string (copies from source).
//     char* s = smake("hello");
#define sdup(s)                 ck_sset(NULL, s)
#define smake(s)                ck_sset(NULL, s)

// sfmake/sfmt: Printf-style string creation/formatting.
//     char* s = sfmake("x=%d", 42);  // Create new
//     sfmt(s, "x=%d", 99);           // Overwrite existing
#define sfmake(fmt, ...)        ck_sfmt(NULL, fmt, __VA_ARGS__)
#define sfmt(s, fmt, ...)       (s = ck_sfmt(s, fmt, __VA_ARGS__))
#define sfmt_append(s, fmt, ...) (s = ck_sfmt_append(s, fmt, __VA_ARGS__))
#define svfmt(s, fmt, args)     (s = ck_svfmt(s, fmt, args))
#define svfmt_append(s, fmt, args) (s = ck_svfmt_append(s, fmt, args))

// sappend/scat: Concatenate b onto a.
#define sappend(a, b)           (a = ck_sappend(a, b))
#define scat(a, b)              sappend(a, b)
#define sappend_range(a, b, e)  (a = ck_sappend_range(a, b, e))
#define scat_range(a, b, e)     sappend_range(a, b, e)

// Comparison. sequ/siequ are NULL-safe (NULL == NULL is true).
#define scmp(a, b)              strcmp(a, b)
#define sicmp(a, b)             ck_stricmp(a, b)
#define sequ(a, b)              ((a) == NULL && (b) == NULL ? 1 : ((a) == NULL || (b) == NULL ? 0 : !strcmp((a), (b))))
#define siequ(a, b)             ((a) == NULL && (b) == NULL ? 1 : ((a) == NULL || (b) == NULL ? 0 : !ck_stricmp((a), (b))))

// sprefix/ssuffix: Check if string starts/ends with given prefix/suffix.
#define sprefix(s, p)           ck_sprefix(s, p)
#define ssuffix(s, p)           ck_ssuffix(s, p)
#define scontains(s, sub)       (slen(s) >= (int)strlen(sub) && !!strstr(s, sub))

// sfirst_index_of/slast_index_of: Find character. Returns -1 if not found.
#define sfirst_index_of(s, ch)  ck_sfirst_index_of(s, ch)
#define slast_index_of(s, ch)   ck_slast_index_of(s, ch)

// sfind: Find substring. Returns pointer to match or NULL.
#define sfind(s, sub)           strstr(s, sub)

// stoupper/stolower: Convert case in-place.
#define stoupper(s)             ck_stoupper(s)
#define stolower(s)             ck_stolower(s)

// shash: FNV-1a hash of string.
#define shash(s)                ck_hash_fnv1a(s, slen(s))

// strim/sltrim/srtrim: Trim whitespace from both ends, left only, or right only.
#define strim(s)                (s = ck_strim(s))
#define sltrim(s)               (s = ck_sltrim(s))
#define srtrim(s)               (s = ck_srtrim(s))

// slpad/srpad: Pad string with n copies of ch on left or right.
#define slpad(s, ch, n)         (s = ck_slpad(s, ch, n))
#define srpad(s, ch, n)         (s = ck_srpad(s, ch, n))

// sreplace: Replace all occurrences of old with new.
#define sreplace(s, old, new)   (s = ck_sreplace(s, old, new))

// sdedup: Collapse consecutive runs of ch into single ch.
#define sdedup(s, ch)           (s = ck_sdedup(s, ch))

// serase: Remove n characters starting at idx.
#define serase(s, idx, n)       (s = ck_serase(s, idx, n))

// ssplit_once: Split at first occurrence of ch.
// Returns left part (new string). Modifies s to contain right part.
//     char* s = smake("a:b:c");
//     char* left = ssplit_once(s, ':');  // left="a", s="b:c"
#define ssplit_once(s, ch)      ck_ssplit_once(s, ch)

// ssplit: Split into array of strings. Caller must afree the array and each string.
//     char** parts = ssplit("a:b:c", ':');
//     for (int i = 0; i < acount(parts); ++i) printf("%s\n", parts[i]);
//     for (int i = 0; i < acount(parts); ++i) sfree(parts[i]);
//     afree(parts);
#define ssplit(s, ch)           ck_ssplit(s, ch)

// Parsing functions.
#define stoint(s)               ck_stoint(s)
#define stouint(s)              ck_stouint(s)
#define stofloat(s)             ck_stofloat(s)
#define stodouble(s)            ck_stodouble(s)
#define stohex(s)               ck_stohex(s)   // Accepts 0x or # prefix.
#define stobool(s)              (!strcmp(s, "true"))

// sappend_UTF8: Append a Unicode codepoint encoded as UTF-8.
#define sappend_UTF8(s, cp)     (s = ck_sappend_UTF8(s, cp))

// decode_UTF8: Decode one UTF-8 codepoint from s, store it in *codepoint.
// Returns a pointer past the consumed bytes. Invalid sequences yield 0xFFFD.
#define decode_UTF8(s, cp)      ck_decode_UTF8(s, cp)

// decode_UTF16: Decode one UTF-16 codepoint (with surrogate pair support) from s, store it in *codepoint.
// Returns a pointer past the consumed uint16_t(s). Invalid surrogates yield 0xFFFD.
#define decode_UTF16(s, cp)     ck_decode_UTF16(s, cp)

// String formatting from primitives.
#define sint(s, i)              sfmt(s, "%d", (int)(i))
#define suint(s, u)             sfmt(s, "%" PRIu64, (uint64_t)(u))
#define sfloat(s, f)            sfmt(s, "%f", (double)(f))
#define sdouble(s, f)           sfmt(s, "%f", (double)(f))
#define shex(s, u)              sfmt(s, "0x%x", (unsigned)(u))
#define sptr(s, p)              sfmt(s, "%p", (void*)(p))
#define sbool(s, b)             sfmt(s, "%s", (b) ? "true" : "false")

//--------------------------------------------------------------------------------------------------
// String path utilities.
//
// All return newly allocated strings (caller must sfree).
//
// Example:
//     char* ext = spext("file.txt");       // -> ".txt"
//     char* norm = spnorm("a/../b");       // -> "/b"
//     sfree(ext); sfree(norm);

// spfname: Get filename from path. "/usr/bin/app" -> "app"
#define spfname(s)              ck_spfname(s)

// spfname_no_ext: Get filename without extension. "/usr/bin/app.exe" -> "app"
#define spfname_no_ext(s)       ck_spfname_no_ext(s)

// spext: Get file extension including dot. "file.txt" -> ".txt"
#define spext(s)                ck_spext(s)

// spext_equ: Check if path has given extension. spext_equ("file.txt", ".txt") -> 1
#define spext_equ(s, ext)       ck_spext_equ(s, ext)

// sppop: Remove last path component. "/a/b/c" -> "/a/b"
#define sppop(s)                ck_sppop(s)

// sppopn: Remove last n path components.
#define sppopn(s, n)            ck_sppopn(s, n)

// spcompact: Shorten path to fit n characters using "...".
#define spcompact(s, n)         ck_spcompact(s, n)

// spdir_of: Get parent directory. "/a/b/c" -> "/b"
#define spdir_of(s)             ck_spdir_of(s)

// sptop_of: Get top-level directory after root.
#define sptop_of(s)             ck_sptop_of(s)

// spnorm: Normalize path. Resolves "..", converts \ to /.
#define spnorm(s)               ck_spnorm(s)

//--------------------------------------------------------------------------------------------------
// Map (typed hashtable) - Stretchy-buffer style.
//
// CK_MAP(T) resolves to T*. NULL is a valid empty map. Keys are always uint64_t.
// For string keys, use sintern() to get a unique pointer and cast to uint64_t.
// Keys and values are stored in parallel dense arrays for fast iteration.
//
// Example:
//     CK_MAP(int) scores = NULL;
//     map_set(scores, sintern("player1"), 100);
//     map_set(scores, sintern("player2"), 200);
//     int p1 = map_get(scores, sintern("player1"));  // 100
//
//     // Iteration:
//     for (int i = 0; i < map_size(scores); i++) {
//         uint64_t key = map_keys(scores)[i];
//         int val = scores[i];  // or map_items(scores)[i]
//     }
//     map_free(scores);

// CK_MAP(T): Declares a map type. Just T* under the hood.
#define CK_MAP(T) T*

// map_size: Number of {key, value} pairs. Returns 0 for NULL.
#define map_size(m)     (map_validate(m), (m) ? CK_MHDR(m)->size : 0)

// map_capacity: Allocated capacity.
#define map_capacity(m) (map_validate(m), (m) ? CK_MHDR(m)->capacity : 0)

// map_keys: Pointer to keys array (uint64_t*). Parallel to items array.
#define map_keys(m)  (map_validate(m), (m) ? ck_map_keys_ptr(CK_MHDR(m)) : NULL)

// map_items: Pointer to items array. Same as the map pointer itself.
#define map_items(m) (map_validate(m), (m))

// map_has: Check if key exists.
#define map_has(m, k) ( \
	map_validate(m), \
	(m) && ck_map_find_impl(CK_MHDR(m), (uint64_t)(k)) >= 0)

// map_del: Remove entry by key. Returns 1 if deleted, 0 if not found.
#define map_del(m, k) ( \
	map_validate(m), \
	(m) ? ck_map_del_impl(CK_MHDR(m), (uint64_t)(k)) : 0)

// map_clear: Remove all entries but keep allocated memory.
#define map_clear(m) do { \
	map_validate(m); \
	if (m) ck_map_clear_impl(CK_MHDR(m)); \
} while(0)

// map_free: Free all memory and set pointer to NULL.
#define map_free(m) do { \
	if (m) { map_validate(m); ck_map_free_impl(CK_MHDR(m)); } \
	(m) = NULL; \
} while(0)

// map_find: Alias for map_get.
#define map_find(m, k) map_get(m, k)

// map_key/map_val: Access key or value at index i.
#define map_key(m, i)  (map_keys(m)[i])
#define map_val(m, i)  (map_items(m)[i])

// map_sort: Sort entries by values. Comparator receives pointers to values.
//     int cmp(const void* a, const void* b) { return *(int*)a - *(int*)b; }
//     map_sort(m, cmp);
#define map_sort(m, cmp) do { \
	map_validate(m); \
	if (m) ck_map_sort_impl(CK_MHDR(m), cmp); \
} while(0)

// map_ssort: Sort entries by keys (treating keys as interned string pointers).
//     map_ssort(m, 0);  // case-sensitive
//     map_ssort(m, 1);  // case-insensitive
#define map_ssort(m, ignore_case) do { \
	map_validate(m); \
	if (m) ck_map_ssort_impl(CK_MHDR(m), ignore_case); \
} while(0)

// map_swap: Swap entries at indices i and j.
#define map_swap(m, i, j) do { \
	map_validate(m); \
	if (m) ck_map_swap_impl(CK_MHDR(m), i, j); \
} while(0)

//--------------------------------------------------------------------------------------------------
// C/C++ specific map macros (type inference differs).
//
// map_get: Get value by key. Returns zero-initialized value if not found.
// map_get_ptr: Get pointer to value. Returns NULL if not found.
// map_set: Set value for key. Creates entry if not exists.
// map_add: Alias for map_set with uint64_t value (backwards compat).

#ifdef __cplusplus
#	include <type_traits>
#	define map_get(m, k) (map_validate(m), (m) && ck_map_find_impl(CK_MHDR(m), (uint64_t)(k)) >= 0 ? (m)[ck_map_find_impl(CK_MHDR(m), (uint64_t)(k))] : std::remove_pointer_t<std::remove_reference_t<decltype(m)>>{})
#	define map_get_ptr(m, k) (map_validate(m), (m) ? (std::remove_reference_t<decltype(m)>)ck_map_get_ptr_impl(CK_MHDR(m), (uint64_t)(k)) : nullptr)
#	define map_set(m, k, v) do { std::remove_pointer_t<decltype(m)> ck_v_ = (v); ck_map_set_stretchy((void**)&(m), (uint64_t)(k), &ck_v_, sizeof(ck_v_)); map_validate(m); } while(0)
#	define map_add(m, k, v) do { uint64_t ck_v_ = (uint64_t)(v); ck_map_set_stretchy((void**)&(m), (uint64_t)(k), &ck_v_, sizeof(ck_v_)); map_validate(m); } while(0)
#else

// map_get: Get value by key. Returns zero-initialized value if not found.
//     int x = map_get(m, sintern("x"));
#define map_get(m, k) ( \
	map_validate(m), \
	(m) && ck_map_find_impl(CK_MHDR(m), (uint64_t)(k)) >= 0 \
		? (m)[ck_map_find_impl(CK_MHDR(m), (uint64_t)(k))] \
		: (typeof(*(m))){ 0 })

// map_get_ptr: Get pointer to value. Returns NULL if not found.
//     int* px = map_get_ptr(m, sintern("x"));
//     if (px) *px = 999;
#define map_get_ptr(m, k) ( \
	map_validate(m), \
	(m) ? (typeof(m))ck_map_get_ptr_impl(CK_MHDR(m), (uint64_t)(k)) : NULL)

// map_set: Set value for key. Creates entry if not exists. May reallocate.
//     map_set(m, sintern("x"), 42);
#define map_set(m, k, v) do { \
	typeof(*(m)) ck_v_ = (v); \
	ck_map_set_stretchy((void**)&(m), (uint64_t)(k), &ck_v_, sizeof(ck_v_)); \
	map_validate(m); \
} while(0)

// map_add: Alias for map_set with uint64_t value (backwards compat).
#define map_add(m, k, v) do { \
	uint64_t ck_v_ = (uint64_t)(v); \
	ck_map_set_stretchy((void**)&(m), (uint64_t)(k), &ck_v_, sizeof(ck_v_)); \
	map_validate(m); \
} while(0)

#endif

//--------------------------------------------------------------------------------------------------
// String interning.
//
// Returns a unique, stable pointer for each unique string.
// Interned strings can be compared by pointer (==) instead of strcmp.
// Great as map keys: cast the pointer to uint64_t.
// Also reduces memory by deduplicating identical strings.
//
// Example:
//     const char* a = sintern("hello");
//     char buf[64]; strcpy(buf, "hel"); strcat(buf, "lo");
//     const char* b = sintern(buf);
//     assert(a == b);  // Same pointer!
//
//     // Use as map key:
//     CK_MAP(int) m = NULL;
//     map_set(m, sintern("x"), 10);
//     int x = map_get(m, sintern("x"));  // 10

// sintern: Return interned string. Same contents always returns same pointer.
#define sintern(s) ck_sintern(s)

// sintern_range: Intern a substring [start, end).
#define sintern_range(start, end) ck_sintern_range(start, end)

// sivalid: True if s is an interned string (from sintern).
#define sivalid(s) (((CK_UniqueString*)(s) - 1)->cookie.val == CK_INTERN_COOKIE)

// silen: Length of an interned string (constant-time).
#define silen(s)   (((CK_UniqueString*)(s) - 1)->len)

// sinuke: Free all interned strings. All previous pointers become invalid.
#define sinuke()   sintern_nuke()

//--------------------------------------------------------------------------------------------------
// Private implementation details.

// Cookie union for debugger-friendly viewing of 4-char magic values.
typedef union CK_Cookie
{
	uint32_t val;
	char c[4];
} CK_Cookie;

// Compile-time cookie value (uint32_t literal, no function call).
#define CK_COOKIE_VAL(a, b, c, d) \
	((uint32_t)(unsigned char)(a) | ((uint32_t)(unsigned char)(b) << 8) | \
	 ((uint32_t)(unsigned char)(c) << 16) | ((uint32_t)(unsigned char)(d) << 24))

// Portable case-insensitive string compare.
#ifdef _WIN32
#	define ck_stricmp _stricmp
#else
#	define ck_stricmp strcasecmp
#endif

// Intern structure for validation and length access.
#define CK_INTERN_COOKIE CK_COOKIE_VAL('I','N','T','R')
typedef struct CK_UniqueString
{
	CK_Cookie cookie;
	int len;
	struct CK_UniqueString* next;
	char* str;
} CK_UniqueString;

// Hidden array header behind the user pointer.
typedef struct CK_ArrayHeader
{
	int size;
	int capacity;
	int is_static;
	char* data;
	CK_Cookie cookie;
} CK_ArrayHeader;

#define CK_AHDR(a)    ((CK_ArrayHeader*)(a) - 1)
#define CK_ACOOKIE    CK_COOKIE_VAL('A','R','R','Y')
#define CK_ACANARY(a) ((a) ? assert(CK_AHDR(a)->cookie.val == CK_ACOOKIE) : (void)0)

#ifndef CK_API
#define CK_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

CK_API void sintern_nuke();
CK_API int ck_sintern_gen();

CK_API void* ck_agrow(const void* a, int new_size, size_t element_size);
CK_API void* ck_astatic(const void* a, int buffer_size, size_t element_size);
CK_API void* ck_aset(const void* a, const void* b, size_t element_size);
CK_API void* ck_arev(const void* a, size_t element_size);

typedef struct CK_MapSlot
{
	uint64_t h;
	int item_index;
	int base_count;
} CK_MapSlot;

// Safety cookie for map validation.
#define CK_MAP_COOKIE CK_COOKIE_VAL('M','A','P','!')

// Compute pointers to arrays within the single allocation.
// items: right after header (header is 24 bytes, 8-byte aligned)
#define ck_map_items_ptr(hdr) ((void*)((hdr) + 1))

// keys: after items (aligned to 8)
#define ck_map_keys_offset(cap, vsz) (sizeof(CK_MapHeader) + CK_ALIGN8((size_t)(cap) * (vsz)))
#define ck_map_keys_ptr(hdr) ((uint64_t*)((char*)(hdr) + ck_map_keys_offset((hdr)->capacity, (hdr)->val_size)))

// islot: after keys (already 8-byte aligned since keys are uint64_t)
#define ck_map_islot_ptr(hdr) ((int*)(ck_map_keys_ptr(hdr) + (hdr)->capacity))

// slots: after islot (aligned to 8 for CK_MapSlot)
#define ck_map_slots_offset(hdr) CK_ALIGN8((size_t)((char*)(ck_map_islot_ptr(hdr) + (hdr)->capacity) - (char*)(hdr)))
#define ck_map_slots_ptr(hdr) ((CK_MapSlot*)((char*)(hdr) + ck_map_slots_offset(hdr)))

// Internal: Get header pointer from map pointer.
#define CK_MHDR(m) ((m) ? (CK_MapHeader*)((char*)(m) - sizeof(CK_MapHeader)) : NULL)

// Internal: Validate map cookie (catches use-after-free, etc).
#define map_validate(m) ((void)(!(m) || (assert(CK_MHDR(m)->cookie.val == CK_MAP_COOKIE), 1)))

// Map header stored just before the values array.
// All arrays (items, keys, islot, slots) are in a single allocation following the header.
// Layout: [Header][items...][keys...][islot...][slots...]
typedef struct CK_MapHeader
{
	CK_Cookie    cookie;         // Safety cookie for validation
	int          val_size;       // Size of each value in bytes
	int          size;           // Number of items
	int          capacity;       // Capacity for items/keys/islot
	int          slot_count;     // Number of used hash slots
	int          slot_capacity;  // Capacity for hash slots (power of 2)
} CK_MapHeader;

// Alignment helper.
#define CK_ALIGN8(x) (((x) + 7) & ~(size_t)7)

// Map implementation functions.
CK_API void  ck_map_set_stretchy(void** m_ptr, uint64_t key, const void* val, int val_size);
CK_API CK_MapHeader* ck_map_ensure_capacity(void** m_ptr, int want_items, int val_size);
CK_API int   ck_map_find_impl(CK_MapHeader* hdr, uint64_t key);
CK_API void* ck_map_get_ptr_impl(CK_MapHeader* hdr, uint64_t key);
CK_API int   ck_map_del_impl(CK_MapHeader* hdr, uint64_t key);
CK_API void  ck_map_clear_impl(CK_MapHeader* hdr);
CK_API void  ck_map_free_impl(CK_MapHeader* hdr);
CK_API void  ck_map_swap_impl(CK_MapHeader* hdr, int i, int j);
CK_API void  ck_map_sort_impl(CK_MapHeader* hdr, int (*cmp)(const void* a, const void* b));
CK_API void  ck_map_ssort_impl(CK_MapHeader* hdr, int ignore_case);

// String implementation functions.
CK_API char* ck_sfit(char* a, int n);
CK_API char* ck_sset(char* a, const char* b);
CK_API char* ck_sfmt(char* s, const char* fmt, ...);
CK_API char* ck_sfmt_append(char* s, const char* fmt, ...);
CK_API char* ck_svfmt(char* s, const char* fmt, va_list args);
CK_API char* ck_svfmt_append(char* s, const char* fmt, va_list args);
CK_API int   ck_sprefix(const char* s, const char* prefix);
CK_API int   ck_ssuffix(const char* s, const char* suffix);
CK_API void  ck_stoupper(char* s);
CK_API void  ck_stolower(char* s);
CK_API char* ck_sappend(char* a, const char* b);
CK_API char* ck_sappend_range(char* a, const char* b, const char* b_end);
CK_API char* ck_strim(char* s);
CK_API char* ck_sltrim(char* s);
CK_API char* ck_srtrim(char* s);
CK_API char* ck_slpad(char* s, char pad, int count);
CK_API char* ck_srpad(char* s, char pad, int count);
CK_API char* ck_ssplit_once(char* s, char split_c);
CK_API char** ck_ssplit(const char* s, char split_c);
CK_API int   ck_sfirst_index_of(const char* s, char c);
CK_API int   ck_slast_index_of(const char* s, char c);
CK_API int   ck_stoint(const char* s);
CK_API uint64_t ck_stouint(const char* s);
CK_API float ck_stofloat(const char* s);
CK_API double ck_stodouble(const char* s);
CK_API uint64_t ck_stohex(const char* s);
CK_API char* ck_sreplace(char* s, const char* replace_me, const char* with_me);
CK_API char* ck_sdedup(char* s, int ch);
CK_API char* ck_serase(char* s, int index, int count);
CK_API char* ck_spop(char* s);
CK_API char* ck_spopn(char* s, int n);
CK_API char* ck_sappend_UTF8(char* s, int codepoint);
CK_API const char* ck_decode_UTF8(const char* s, int* codepoint);

// Path implementation functions.
CK_API char* ck_spfname(const char* path);
CK_API char* ck_spfname_no_ext(const char* path);
CK_API char* ck_spext(const char* path);
CK_API int   ck_spext_equ(const char* path, const char* ext);
CK_API char* ck_sppop(const char* path);
CK_API char* ck_sppopn(const char* path, int n);
CK_API char* ck_spcompact(const char* path, int n);
CK_API char* ck_spdir_of(const char* path);
CK_API char* ck_sptop_of(const char* path);
CK_API char* ck_spnorm(const char* path);

// UTF16 decode.
#include <stdint.h>
CK_API const uint16_t* ck_decode_UTF16(const uint16_t* s, int* codepoint);

// Intern implementation.
CK_API const char* ck_sintern(const char* s);
CK_API const char* ck_sintern_range(const char* start, const char* end);
CK_API uint64_t ck_hash_fnv1a(const void* ptr, size_t sz);

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CKIT_H

//--------------------------------------------------------------------------------------------------

#ifdef CKIT_IMPLEMENTATION

#ifndef CKIT_IMPLEMENTATION_GUARD
#define CKIT_IMPLEMENTATION_GUARD

#include <stdarg.h>
#include <ctype.h>
#include <stdio.h>
#include <inttypes.h>

#ifdef __cplusplus
extern "C" {
#endif

//--------------------------------------------------------------------------------------------------
// Array.

void* ck_agrow(const void* a, int new_size, size_t element_size)
{
	CK_ACANARY(a);
	assert(acap(a) <= (SIZE_MAX - 1) / 2);
	int new_capacity = 2 * acap(a);
	int min_cap = new_size > 16 ? new_size : 16;
	if (new_capacity < min_cap) new_capacity = min_cap;
	assert(new_size <= new_capacity);
	assert((size_t)new_capacity <= (SIZE_MAX - sizeof(CK_ArrayHeader)) / element_size);
	size_t total_size = sizeof(CK_ArrayHeader) + (size_t)new_capacity * element_size;
	CK_ArrayHeader* hdr;
	if (a) {
		if (!CK_AHDR(a)->is_static) {
			hdr = (CK_ArrayHeader*)CK_REALLOC(CK_AHDR(a), total_size);
		} else {
			hdr = (CK_ArrayHeader*)CK_ALLOC(total_size);
			memcpy(hdr + 1, a, (size_t)asize(a) * element_size);
			hdr->size = asize(a);
			hdr->cookie.val = CK_ACOOKIE;
		}
	} else {
		hdr = (CK_ArrayHeader*)CK_ALLOC(total_size);
		hdr->size = 0;
		hdr->cookie.val = CK_ACOOKIE;
	}
	hdr->capacity = new_capacity;
	hdr->is_static = 0;
	hdr->data = (char*)(hdr + 1);
	return (void*)(hdr + 1);
}

void* ck_astatic(const void* a, int buffer_size, size_t element_size)
{
	CK_ArrayHeader* hdr = (CK_ArrayHeader*)a;
	hdr->size = 0;
	hdr->cookie.val = CK_ACOOKIE;
	if (sizeof(CK_ArrayHeader) <= element_size) {
		hdr->capacity = buffer_size / (int)element_size - 1;
	} else {
		int elements_taken = (int)sizeof(CK_ArrayHeader) / (int)element_size + ((int)sizeof(CK_ArrayHeader) % (int)element_size > 0);
		hdr->capacity = buffer_size / (int)element_size - elements_taken;
	}
	hdr->data = (char*)(hdr + 1);
	hdr->is_static = 1;
	return (void*)(hdr + 1);
}

void* ck_aset(const void* a, const void* b, size_t element_size)
{
	CK_ACANARY(a);
	if (!b) {
		aclear(a);
		return (void*)a;
	}
	CK_ACANARY(b);
	if (acap(a) < asize(b)) {
		a = ck_agrow(a, asize(b), element_size);
	}
	memcpy((void*)a, b, (size_t)asize(b) * element_size);
	if (a) asetlen(a, asize(b));
	else assert(!asize(b));
	return (void*)a;
}

void* ck_arev(const void* a_ptr, size_t element_size)
{
	CK_ACANARY(a_ptr);
	char* a = (char*)a_ptr;
	int n = acount(a_ptr);
	if (n <= 1) return (void*)a_ptr;
	char* b = a + element_size * (n - 1);
	char tmp[256];
	char* t = element_size <= sizeof(tmp) ? tmp : (char*)CK_ALLOC(element_size);
	while (a < b) {
		memcpy(t, a, element_size);
		memcpy(a, b, element_size);
		memcpy(b, t, element_size);
		a += element_size;
		b -= element_size;
	}
	if (t != tmp) CK_FREE(t);
	return (void*)a_ptr;
}

//--------------------------------------------------------------------------------------------------
// Strings.

char* ck_sfit(char* a, int n)
{
	afit(a, n + 1);
	if (scount(a) == 0) apush(a, 0);
	return a;
}

char* ck_sset(char* a, const char* b)
{
	CK_ACANARY(a);
	if (!b) return NULL;
	int bsize = (int)strlen(b) + 1;
	if (acap(a) < bsize) {
		a = (char*)ck_agrow(a, bsize, 1);
	}
	memcpy(a, b, (size_t)bsize);
	asetlen(a, bsize);
	return a;
}

char* ck_sfmt(char* s, const char* fmt, ...)
{
	CK_ACANARY(s);
	va_list args;
	va_start(args, fmt);
	int n = 1 + vsnprintf(s, (size_t)scap(s), fmt, args);
	va_end(args);
	if (n > scap(s)) {
		sfit(s, n);
		va_start(args, fmt);
		n = 1 + vsnprintf(s, (size_t)scap(s), fmt, args);
		va_end(args);
	}
	asetlen(s, n);
	return s;
}

char* ck_sfmt_append(char* s, const char* fmt, ...)
{
	CK_ACANARY(s);
	va_list args;
	va_start(args, fmt);
	int nul = !!s;
	int capacity = scap(s) - scount(s);
	int n = 1 + vsnprintf(s + slen(s), (size_t)(capacity > 0 ? capacity : 0), fmt, args);
	va_end(args);
	if (n > capacity) {
		afit(s, n + scount(s));
		va_start(args, fmt);
		int new_capacity = scap(s) - scount(s);
		n = 1 + vsnprintf(s + slen(s), (size_t)new_capacity, fmt, args);
		assert(n <= new_capacity);
		va_end(args);
	}
	asetlen(s, asize(s) + n - nul);
	return s;
}

char* ck_svfmt(char* s, const char* fmt, va_list args)
{
	CK_ACANARY(s);
	va_list copy_args;
	va_copy(copy_args, args);
	int n = 1 + vsnprintf(s, (size_t)scap(s), fmt, args);
	if (n > scap(s)) {
		sfit(s, n);
		n = 1 + vsnprintf(s, (size_t)scap(s), fmt, copy_args);
	}
	va_end(copy_args);
	asetlen(s, n);
	return s;
}

char* ck_svfmt_append(char* s, const char* fmt, va_list args)
{
	CK_ACANARY(s);
	va_list copy_args;
	va_copy(copy_args, args);
	int nul = !!s;
	int capacity = scap(s) - scount(s);
	int n = 1 + vsnprintf(s + slen(s), (size_t)(capacity > 0 ? capacity : 0), fmt, copy_args);
	va_end(copy_args);
	if (n > capacity) {
		afit(s, n + scount(s));
		int new_capacity = scap(s) - scount(s);
		n = 1 + vsnprintf(s + slen(s), (size_t)new_capacity, fmt, args);
		assert(n <= new_capacity);
	}
	asetlen(s, asize(s) + n - nul);
	return s;
}

int ck_sprefix(const char* s, const char* prefix)
{
	if (!s) return 0;
	int s_len = (int)strlen(s);
	int prefix_len = (int)strlen(prefix);
	return s_len >= prefix_len && !memcmp(s, prefix, (size_t)prefix_len);
}

int ck_ssuffix(const char* s, const char* suffix)
{
	if (!s) return 0;
	int s_len = (int)strlen(s);
	int suffix_len = (int)strlen(suffix);
	return s_len >= suffix_len && !memcmp(s + s_len - suffix_len, suffix, (size_t)suffix_len);
}

void ck_stoupper(char* s)
{
	if (!s) return;
	for (int i = 0; i < slen(s); ++i) s[i] = (char)toupper((unsigned char)s[i]);
}

void ck_stolower(char* s)
{
	if (!s) return;
	for (int i = 0; i < slen(s); ++i) s[i] = (char)tolower((unsigned char)s[i]);
}

char* ck_sappend(char* a, const char* b)
{
	CK_ACANARY(a);
	int blen = (int)strlen(b);
	if (blen <= 0) return a;
	sfit(a, slen(a) + blen + 1);
	memcpy(a + slen(a), b, (size_t)blen);
	asetlen(a, asize(a) + blen);
	a[slen(a)] = 0;
	return a;
}

char* ck_sappend_range(char* a, const char* b, const char* b_end)
{
	CK_ACANARY(a);
	int blen = (int)(b_end - b);
	if (blen <= 0) return a;
	sfit(a, slen(a) + blen + 1);
	memcpy(a + slen(a), b, (size_t)blen);
	asetlen(a, asize(a) + blen);
	a[slen(a)] = 0;
	return a;
}

char* ck_strim(char* s)
{
	CK_ACANARY(s);
	int len = slen(s);
	if (len <= 0) return s;
	char* start = s;
	char* end = s + len - 1;
	while (start <= end && isspace((unsigned char)*start)) start++;
	while (end >= start && isspace((unsigned char)*end)) end--;
	int new_len = (int)(end >= start ? (end - start + 1) : 0);
	if (new_len > 0) memmove(s, start, (size_t)new_len);
	s[new_len] = 0;
	asetlen(s, new_len + 1);
	return s;
}

char* ck_sltrim(char* s)
{
	CK_ACANARY(s);
	int len = slen(s);
	if (len <= 0) return s;
	char* start = s;
	char* end = s + len - 1;
	while (start <= end && isspace((unsigned char)*start)) start++;
	int new_len = (int)(end >= start ? (end - start + 1) : 0);
	if (new_len > 0) memmove(s, start, (size_t)new_len);
	s[new_len] = 0;
	asetlen(s, new_len + 1);
	return s;
}

char* ck_srtrim(char* s)
{
	CK_ACANARY(s);
	int len = slen(s);
	if (len <= 0) return s;
	char* end = s + len - 1;
	while (end >= s && isspace((unsigned char)*end)) --end;
	int new_len = (int)(end - s + 1);
	if (new_len < 0) new_len = 0;
	s[new_len] = 0;
	asetlen(s, new_len + 1);
	return s;
}

char* ck_slpad(char* s, char pad, int count)
{
	CK_ACANARY(s);
	sfit(s, scount(s) + count);
	memmove(s + count, s, (size_t)scount(s));
	memset(s, pad, (size_t)count);
	asetlen(s, asize(s) + count);
	return s;
}

char* ck_srpad(char* s, char pad, int count)
{
	CK_ACANARY(s);
	sfit(s, scount(s) + count);
	memset(s + slen(s), pad, (size_t)count);
	asetlen(s, asize(s) + count);
	s[slen(s)] = 0;
	return s;
}

char* ck_ssplit_once(char* s, char split_c)
{
	CK_ACANARY(s);
	char* start = s;
	char* end = s + slen(s);
	while (start < end) {
		if (*start == split_c) break;
		++start;
	}
	if (start == end) return NULL; // Delimiter not found.
	int len = (int)(start - s);
	char* split = NULL;
	sfit(split, len + 1);
	asetlen(split, len + 1);
	memcpy(split, s, (size_t)len);
	split[len] = 0;
	int new_len = slen(s) - len - 1;
	memmove(s, s + len + 1, (size_t)new_len);
	asetlen(s, new_len + 1);
	s[new_len] = 0;
	return split;
}

char** ck_ssplit(const char* s, char split_c)
{
	char* copy = NULL;
	char** result = NULL;
	char* split = NULL;
	sset(copy, s);
	while ((split = ssplit_once(copy, split_c))) {
		apush(result, split);
	}
	apush(result, copy);
	return result;
}

int ck_sfirst_index_of(const char* s, char c)
{
	if (!s) return -1;
	const char* p = strchr(s, c);
	if (!p) return -1;
	return (int)(p - s);
}

int ck_slast_index_of(const char* s, char c)
{
	if (!s) return -1;
	const char* p = strrchr(s, c);
	if (!p) return -1;
	return (int)(p - s);
}

int ck_stoint(const char* s)
{
	char* end;
	long long result = strtoll(s, &end, 10);
	return (int)result;
}

uint64_t ck_stouint(const char* s)
{
	char* end;
	uint64_t result = (uint64_t)strtoll(s, &end, 10);
	return result;
}

float ck_stofloat(const char* s)
{
	char* end;
	double result = strtod(s, &end);
	return (float)result;
}

double ck_stodouble(const char* s)
{
	char* end;
	double result = strtod(s, &end);
	return result;
}

uint64_t ck_stohex(const char* s)
{
	if (!strncmp(s, "#", 1)) s += 1;
	if (!strncmp(s, "0x", 2)) s += 2;
	int len = (int)strlen(s);
	if (len != 6 && len != 8) return 0;
	char* end;
	uint64_t result = (uint64_t)strtoll(s, &end, 16);
	return len == 6 ? ((result << 8) | 0xFF) : result;
}

char* ck_sreplace(char* s, const char* replace_me, const char* with_me)
{
	CK_ACANARY(s);
	if (!s) return NULL;
	size_t replace_len = strlen(replace_me);
	if (replace_len == 0) return s; // Empty pattern: nothing to replace.
	size_t with_len = strlen(with_me);
	char* find;
	char* search = s;
	while ((find = strstr(search, replace_me))) {
		int find_offset = (int)(find - s);
		if (replace_len > with_len) {
			int remaining = scount(s) - find_offset - (int)replace_len;
			int diff = (int)(replace_len - with_len);
			memcpy(find, with_me, with_len);
			memmove(find + with_len, find + replace_len, (size_t)remaining);
			asetlen(s, asize(s) - diff);
		} else {
			int remaining = scount(s) - find_offset - (int)replace_len;
			int diff = (int)(with_len - replace_len);
			sfit(s, scount(s) + diff);
			find = s + find_offset;
			memmove(find + with_len, find + replace_len, (size_t)remaining);
			memcpy(find, with_me, with_len);
			asetlen(s, asize(s) + diff);
		}
		search = find + with_len;
	}
	return s;
}

char* ck_sdedup(char* s, int ch)
{
	CK_ACANARY(s);
	if (!s) return NULL;
	int len = (int)strlen(s);
	if (len <= 1) return s; // Empty or single char: nothing to dedup.
	int i = 0, j = 1;
	int dup = 0;
	while (j < len) {
		if (s[i] == ch && s[j] == ch) {
			dup = 1;
			++j;
		} else {
			++i;
			if (dup) s[i] = s[j];
			++j;
		}
	}
	s[i + 1] = 0;
	asetlen(s, i + 2);
	return s;
}

char* ck_serase(char* s, int index, int count)
{
	CK_ACANARY(s);
	if (!s) return NULL;
	if (index < 0) {
		count += index;
		index = 0;
		if (count <= 0) return s;
	}
	if (index >= slen(s)) return s;
	if (index + count >= slen(s)) {
		asetlen(s, index + 1);
		s[index] = 0;
		return s;
	} else {
		int remaining = scount(s) - (count + index);
		memmove(s + index, s + count + index, (size_t)remaining);
		asetlen(s, asize(s) - count);
	}
	return s;
}

char* ck_spop(char* s)
{
	CK_ACANARY(s);
	if (s && slen(s)) { asetlen(s, asize(s) - 1); s[slen(s)] = 0; }
	return s;
}

char* ck_spopn(char* s, int n)
{
	CK_ACANARY(s);
	if (!s || n < 0) return s;
	while (scount(s) > 1 && n--) asetlen(s, asize(s) - 1);
	s[slen(s)] = 0;
	return s;
}

char* ck_sappend_UTF8(char* s, int codepoint)
{
	CK_ACANARY(s);
	if (codepoint > 0x10FFFF) codepoint = 0xFFFD;
#define CK_EMIT(X, Y, Z) spush(s, (char)(X | ((codepoint >> Y) & Z)))
	     if (codepoint <    0x80) { CK_EMIT(0x00,0,0x7F); }
	else if (codepoint <   0x800) { CK_EMIT(0xC0,6,0x1F); CK_EMIT(0x80, 0,  0x3F); }
	else if (codepoint < 0x10000) { CK_EMIT(0xE0,12,0xF); CK_EMIT(0x80, 6,  0x3F); CK_EMIT(0x80, 0, 0x3F); }
	else                          { CK_EMIT(0xF0,18,0x7); CK_EMIT(0x80, 12, 0x3F); CK_EMIT(0x80, 6, 0x3F); CK_EMIT(0x80, 0, 0x3F); }
#undef CK_EMIT
	return s;
}

const char* ck_decode_UTF8(const char* s, int* codepoint)
{
	unsigned char c = (unsigned char)*s++;
	int extra = 0, min = 0;
	*codepoint = 0;
	     if (c >= 0xF0) { *codepoint = c & 0x07; extra = 3; min = 0x10000; }
	else if (c >= 0xE0) { *codepoint = c & 0x0F; extra = 2; min = 0x800; }
	else if (c >= 0xC0) { *codepoint = c & 0x1F; extra = 1; min = 0x80; }
	else if (c >= 0x80) { *codepoint = 0xFFFD; }
	else *codepoint = c;
	while (extra--) {
		c = (unsigned char)*s++;
		if ((c & 0xC0) != 0x80) { *codepoint = 0xFFFD; }
		if (*codepoint != 0xFFFD) { *codepoint = ((*codepoint) << 6) | (c & 0x3F); }
	}
	if (*codepoint < min) *codepoint = 0xFFFD;
	return s;
}

//--------------------------------------------------------------------------------------------------
// String path utilities.

char* ck_spfname(const char* path)
{
	if (!path || path[0] == '\0') return NULL;
	int at = slast_index_of(path, '/');
	const char* f = path + at + 1;
	if (f[0] != '\0') return smake(f);
	return NULL;
}

char* ck_spfname_no_ext(const char* path)
{
	char* s = ck_spfname(path);
	if (!s) return NULL;
	int at = slast_index_of(s, '.');
	if (at == 0) { sfree(s); return NULL; }
	if (at == -1) return s;
	serase(s, at, slen(s) - at);
	return s;
}

char* ck_spext(const char* path)
{
	int at = slast_index_of(path, '.');
	if (at == -1 || path[at + 1] == 0 || path[at + 1] == '/') return NULL;
	return smake(path + at);
}

int ck_spext_equ(const char* path, const char* ext)
{
	int at = slast_index_of(path, '.');
	if (at == -1 || path[at + 1] == 0 || path[at + 1] == '/') return 0;
	return sequ(path + at, ext);
}

char* ck_sppop(const char* path)
{
	char* s = sdup(path);
	if (slast(s) == '/') spop(s);
	int at = slast_index_of(s, '/');
	if (at == -1 || at == 0) return (sset(s, "/"), s);
	serase(s, at, slen(s) - at);
	return s;
}

char* ck_sppopn(const char* path, int n)
{
	char* s = sdup(path);
	while (n--) {
		if (slast(s) == '/') spop(s);
		int at = slast_index_of(s, '/');
		if (at == -1 || at == 0) { sset(s, "/"); continue; }
		serase(s, at, slen(s) - at);
	}
	return s;
}

char* ck_spcompact(const char* path, int n)
{
	int len = (int)strlen(path);
	if (n <= 6) return NULL;
	if (len < n) return sdup(path);
	int at = slast_index_of(path, '/');
	if (at == -1 || at == 0) {
		char* s = sdup(path);
		serase(s, n, slen(s) - n);
		serase(s, n - 3, 3);
		sappend(s, "...");
		return s;
	}
	int remaining = len - at - 1;
	if (remaining >= n - 3) {
		char* s = smake("...");
		sappend_range(s, path, path + at - 6);
		sappend(s, "...");
		return s;
	} else {
		char* s = sdup(path);
		int len_s = slen(s);
		int to_erase = len_s - (remaining - 3);
		serase(s, remaining - 3, to_erase);
		sappend(s, "...");
		sappend(s, path + at);
		return s;
	}
}

char* ck_spdir_of(const char* path)
{
	if (!*path || (*path == '.' && (int)strlen(path) < 3)) return NULL;
	if (sequ(path, "../")) return NULL;
	if (sequ(path, "/")) return NULL;
	int at = slast_index_of(path, '/');
	if (at == -1) return NULL;
	if (at == 0) return smake("/");
	char* s = smake(path);
	serase(s, at, slen(s) - at);
	at = slast_index_of(s, '/');
	if (at == -1) {
		if (slen(s) == 2) {
			return s;
		} else {
			s[0] = '/';
			return s;
		}
	}
	serase(s, 0, at);
	return s;
}

char* ck_sptop_of(const char* path)
{
	int at = sfirst_index_of(path, '/');
	if (at == -1) return NULL;
	int next = sfirst_index_of(path + at + 1, '/');
	if (next == -1) return smake("/");
	char* s = sdup(path + at);
	serase(s, next + 1, slen(s) - (next + 1));
	return s;
}

char* ck_spnorm(const char* path)
{
	char* result = NULL;
	int len = (int)strlen(path);
	if (*path != '\\' && *path != '/') {
		int windows_drive = len >= 2 && path[1] == ':';
		if (!windows_drive) {
			spush(result, '/');
		}
	}
	int prev_was_dot = 0;
	int prev_was_dotdot = 0;
	for (int i = 0; i < len; ++i) {
		char c = path[i];
		if (c == '\\' || c == '/') {
			if (!prev_was_dot) {
				spush(result, '/');
			} else if (prev_was_dotdot) {
				char* tmp = ck_sppop(result);
				sfree(result);
				result = tmp;
				spush(result, '/');
			}
			prev_was_dot = 0;
			prev_was_dotdot = 0;
		} else if (c == '.') {
			if (prev_was_dot) prev_was_dotdot = 1;
			prev_was_dot = 1;
		} else {
			if (prev_was_dot) spush(result, '.');
			spush(result, c);
			prev_was_dot = 0;
			prev_was_dotdot = 0;
		}
	}
	sreplace(result, "//", "/");
	if (slen(result) > 1 && slast(result) == '/') spop(result);
	return result;
}

//--------------------------------------------------------------------------------------------------
// Map implementation (originally by Mattias Gustavsson) - Stretchy buffer style.

uint64_t ck_map_hash(uint64_t x)
{
	x += 0x9e3779b97f4a7c15ull;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
	x ^= (x >> 31);
	return x ? x : 1ull;
}

void ck_map_zero_slots(CK_MapHeader* hdr)
{
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	for (int i = 0; i < hdr->slot_capacity; ++i) {
		slots[i].h = 0;
		slots[i].item_index = -1;
		slots[i].base_count = 0;
	}
	hdr->slot_count = 0;
}

// Total allocation size for map with given capacities.
size_t ck_map_alloc_size(int capacity, int val_size, int slot_capacity)
{
	size_t items_end = sizeof(CK_MapHeader) + CK_ALIGN8((size_t)capacity * val_size);
	size_t keys_end = items_end + (size_t)capacity * sizeof(uint64_t);
	size_t islot_end = keys_end + (size_t)capacity * sizeof(int);
	size_t slots_start = CK_ALIGN8(islot_end);
	return slots_start + (size_t)slot_capacity * sizeof(CK_MapSlot);
}

int ck_map_find_insertion_slot(CK_MapHeader* hdr, uint64_t h)
{
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	int mask = hdr->slot_capacity - 1;
	int base = (int)(h & (uint64_t)mask);
	int slot = base, first_free = base;
	int remaining = slots[base].base_count;

	while (remaining > 0) {
		if (slots[slot].item_index < 0 && slots[first_free].item_index >= 0)
			first_free = slot;
		uint64_t sh = slots[slot].h;
		if (sh) {
			if (((int)(sh & (uint64_t)mask)) == base) --remaining;
		}
		slot = (slot + 1) & mask;
	}

	slot = first_free;
	while (slots[slot].item_index >= 0) {
		slot = (slot + 1) & mask;
	}
	return slot;
}

int ck_map_find_slot(const CK_MapHeader* hdr, uint64_t key, uint64_t h)
{
	if (hdr->slot_capacity == 0) return -1;
	CK_MapSlot* slots = ck_map_slots_ptr((CK_MapHeader*)hdr);
	uint64_t* keys = ck_map_keys_ptr((CK_MapHeader*)hdr);
	int mask = hdr->slot_capacity - 1;
	int base = (int)(h & (uint64_t)mask);
	int slot = base;
	int remaining = slots[base].base_count;
	while (remaining > 0) {
		int item = slots[slot].item_index;
		if (item >= 0) {
			uint64_t sh = slots[slot].h;
			if (((int)(sh & (uint64_t)mask)) == base) {
				--remaining;
				if (sh == h && keys[item] == key) return slot;
			}
		}
		slot = (slot + 1) & mask;
	}
	return -1;
}

// Rebuild hash table after reallocation. Assumes slots are already zeroed.
void ck_map_rebuild_slots(CK_MapHeader* hdr)
{
	if (hdr->slot_capacity == 0) return;
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	uint64_t* keys = ck_map_keys_ptr(hdr);
	int* islot = ck_map_islot_ptr(hdr);
	for (int idx = 0; idx < hdr->size; ++idx) {
		uint64_t h = ck_map_hash(keys[idx]);
		int slot = ck_map_find_insertion_slot(hdr, h);
		int base = (int)(h & (uint64_t)(hdr->slot_capacity - 1));
		slots[slot].h = h;
		slots[slot].item_index = idx;
		++slots[base].base_count;
		islot[idx] = slot;
		++hdr->slot_count;
	}
}

// Ensure we have capacity for 'want' items. May reallocate and update m_ptr.
// Returns the possibly-updated header pointer.
// Single allocation: [Header][items][keys][islot][slots]
CK_MapHeader* ck_map_ensure_capacity(void** m_ptr, int want_items, int val_size)
{
	CK_MapHeader* hdr = *m_ptr ? (CK_MapHeader*)((char*)*m_ptr - sizeof(CK_MapHeader)) : NULL;
	int old_item_cap = hdr ? hdr->capacity : 0;
	int old_slot_cap = hdr ? hdr->slot_capacity : 0;
	int old_size = hdr ? hdr->size : 0;

	// Compute new item capacity.
	int new_item_cap = old_item_cap;
	if (want_items > new_item_cap) {
		new_item_cap = old_item_cap ? old_item_cap * 2 : 16;
		while (new_item_cap < want_items) new_item_cap *= 2;
	}

	// Compute new slot capacity. Slots should be at least 2x items for good hash distribution.
	int new_slot_cap = old_slot_cap;
	int min_slot_cap = new_item_cap * 2;
	if (min_slot_cap < 16) min_slot_cap = 16;
	if (new_slot_cap < min_slot_cap) {
		new_slot_cap = 16;
		while (new_slot_cap < min_slot_cap) new_slot_cap *= 2;
	}

	// Also check load factor on existing slots.
	if (hdr && hdr->slot_capacity > 0) {
		int thresh = hdr->slot_capacity - (hdr->slot_capacity >> 2); // 75%
		if (hdr->slot_count >= thresh) {
			int grown = hdr->slot_capacity * 2;
			if (grown > new_slot_cap) new_slot_cap = grown;
		}
	}

	// If no change needed, return current.
	if (new_item_cap == old_item_cap && new_slot_cap == old_slot_cap) {
		return hdr;
	}

	// Allocate new block.
	size_t new_size = ck_map_alloc_size(new_item_cap, val_size, new_slot_cap);
	CK_MapHeader* new_hdr = (CK_MapHeader*)CK_ALLOC(new_size);
	memset(new_hdr, 0, new_size);

	// Initialize header.
	new_hdr->cookie.val = CK_MAP_COOKIE;
	new_hdr->val_size = val_size;
	new_hdr->size = old_size;
	new_hdr->capacity = new_item_cap;
	new_hdr->slot_count = 0;
	new_hdr->slot_capacity = new_slot_cap;

	// Copy existing data if any.
	if (hdr && old_size > 0) {
		memcpy(ck_map_items_ptr(new_hdr), ck_map_items_ptr(hdr), (size_t)old_size * val_size);
		memcpy(ck_map_keys_ptr(new_hdr), ck_map_keys_ptr(hdr), (size_t)old_size * sizeof(uint64_t));
	}

	// Initialize islot to -1.
	int* new_islot = ck_map_islot_ptr(new_hdr);
	for (int i = 0; i < new_item_cap; i++) new_islot[i] = -1;

	// Initialize slots item_index to -1 (0 would appear occupied).
	CK_MapSlot* new_slots = ck_map_slots_ptr(new_hdr);
	for (int i = 0; i < new_slot_cap; i++) new_slots[i].item_index = -1;

	// Rebuild hash table.
	ck_map_rebuild_slots(new_hdr);

	// Free old allocation.
	if (hdr) CK_FREE(hdr);

	// Update user pointer.
	*m_ptr = ck_map_items_ptr(new_hdr);
	return new_hdr;
}

int ck_map_find_impl(CK_MapHeader* hdr, uint64_t key)
{
	if (!hdr || hdr->size == 0 || hdr->slot_capacity == 0) return -1;
	uint64_t h = ck_map_hash(key);
	int s = ck_map_find_slot(hdr, key, h);
	if (s < 0) return -1;
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	return slots[s].item_index;
}

void ck_map_set_stretchy(void** m_ptr, uint64_t key, const void* val, int val_size)
{
	CK_MapHeader* hdr = *m_ptr ? (CK_MapHeader*)((char*)*m_ptr - sizeof(CK_MapHeader)) : NULL;
	uint64_t h = ck_map_hash(key);

	// Update if entry already exists.
	if (hdr && hdr->slot_capacity) {
		int s = ck_map_find_slot(hdr, key, h);
		if (s >= 0) {
			CK_MapSlot* slots = ck_map_slots_ptr(hdr);
			int idx = slots[s].item_index;
			char* items = (char*)ck_map_items_ptr(hdr);
			memcpy(items + idx * hdr->val_size, val, (size_t)hdr->val_size);
			return;
		}
	}

	// Ensure capacity (may reallocate and rebuild hash table).
	int want = (hdr ? hdr->size : 0) + 1;
	hdr = ck_map_ensure_capacity(m_ptr, want, val_size);
	assert(hdr->val_size == val_size);

	// Get array pointers (may have changed after realloc).
	char* items = (char*)ck_map_items_ptr(hdr);
	uint64_t* keys = ck_map_keys_ptr(hdr);
	int* islot = ck_map_islot_ptr(hdr);
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);

	// Append to dense set.
	int idx = hdr->size++;
	keys[idx] = key;
	memcpy(items + idx * hdr->val_size, val, (size_t)hdr->val_size);

	// Insert into hash slots.
	int slot = ck_map_find_insertion_slot(hdr, h);
	int base = (int)(h & (uint64_t)(hdr->slot_capacity - 1));
	slots[slot].h = h;
	slots[slot].item_index = idx;
	++slots[base].base_count;
	islot[idx] = slot;
	++hdr->slot_count;
}

void* ck_map_get_ptr_impl(CK_MapHeader* hdr, uint64_t key)
{
	if (!hdr || hdr->size == 0 || hdr->slot_capacity == 0) return NULL;
	uint64_t h = ck_map_hash(key);
	int s = ck_map_find_slot(hdr, key, h);
	if (s < 0) return NULL;
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	return (char*)ck_map_items_ptr(hdr) + slots[s].item_index * hdr->val_size;
}

int ck_map_del_impl(CK_MapHeader* hdr, uint64_t key)
{
	if (!hdr || hdr->size == 0 || hdr->slot_capacity == 0) return 0;
	uint64_t h = ck_map_hash(key);
	int s = ck_map_find_slot(hdr, key, h);
	if (s < 0) return 0;

	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	uint64_t* keys = ck_map_keys_ptr(hdr);
	int* islot = ck_map_islot_ptr(hdr);
	char* items = (char*)ck_map_items_ptr(hdr);

	int mask = hdr->slot_capacity - 1;
	int base = (int)(h & (uint64_t)mask);
	int idx = slots[s].item_index;
	int last = hdr->size - 1;

	--slots[base].base_count;
	slots[s].item_index = -1;
	--hdr->slot_count;

	if (idx != last) {
		keys[idx] = keys[last];
		memcpy(items + idx * hdr->val_size, items + last * hdr->val_size, (size_t)hdr->val_size);
		int ms = islot[last];
		slots[ms].item_index = idx;
		islot[idx] = ms;
	}
	islot[last] = -1;
	--hdr->size;
	return 1;
}

void ck_map_swap_impl(CK_MapHeader* hdr, int i, int j)
{
	if (i == j) return;
	assert(i >= 0 && i < hdr->size);
	assert(j >= 0 && j < hdr->size);

	uint64_t* keys = ck_map_keys_ptr(hdr);
	int* islot = ck_map_islot_ptr(hdr);
	CK_MapSlot* slots = ck_map_slots_ptr(hdr);
	char* items = (char*)ck_map_items_ptr(hdr);

	// Swap keys.
	uint64_t tk = keys[i];
	keys[i] = keys[j];
	keys[j] = tk;

	// Swap values.
	char tmp[256];
	char* t = hdr->val_size <= (int)sizeof(tmp) ? tmp : (char*)CK_ALLOC((size_t)hdr->val_size);
	memcpy(t, items + i * hdr->val_size, (size_t)hdr->val_size);
	memcpy(items + i * hdr->val_size, items + j * hdr->val_size, (size_t)hdr->val_size);
	memcpy(items + j * hdr->val_size, t, (size_t)hdr->val_size);
	if (t != tmp) CK_FREE(t);

	// Swap islot backpointers.
	int si = islot[i];
	int sj = islot[j];
	if (si >= 0) slots[si].item_index = j;
	if (sj >= 0) slots[sj].item_index = i;
	islot[i] = sj;
	islot[j] = si;
}

void ck_map_sort_range(CK_MapHeader* hdr, int offset, int count, int (*cmp)(const void* a, const void* b))
{
	if (count <= 1) return;
	char* items = (char*)ck_map_items_ptr(hdr);
	char* pivot = items + (offset + count - 1) * hdr->val_size;
	int lo = 0;
	for (int hi = 0; hi < count - 1; ++hi) {
		char* hi_val = items + (offset + hi) * hdr->val_size;
		if (cmp(hi_val, pivot) < 0) {
			ck_map_swap_impl(hdr, offset + lo, offset + hi);
			++lo;
		}
	}
	ck_map_swap_impl(hdr, offset + (count - 1), offset + lo);
	ck_map_sort_range(hdr, offset, lo, cmp);
	ck_map_sort_range(hdr, offset + lo + 1, count - 1 - lo, cmp);
}

void ck_map_sort_impl(CK_MapHeader* hdr, int (*cmp)(const void* a, const void* b))
{
	ck_map_sort_range(hdr, 0, hdr->size, cmp);
}

void ck_map_ssort_range(CK_MapHeader* hdr, int offset, int count, int ignore_case)
{
	if (count <= 1) return;
	uint64_t* keys = ck_map_keys_ptr(hdr);

	#define CK_SSORT_KEY(i) (keys[(offset) + (i)])

	uint64_t pivot = CK_SSORT_KEY(count - 1);
	int lo = 0;
	for (int hi = 0; hi < count - 1; ++hi) {
		const char* hi_key = (const char*)CK_SSORT_KEY(hi);
		const char* pivot_key = (const char*)pivot;
		int cmp = ignore_case ? ck_stricmp(hi_key, pivot_key) : strcmp(hi_key, pivot_key);
		if (cmp < 0) {
			ck_map_swap_impl(hdr, offset + lo, offset + hi);
			++lo;
		}
	}
	ck_map_swap_impl(hdr, offset + (count - 1), offset + lo);
	ck_map_ssort_range(hdr, offset, lo, ignore_case);
	ck_map_ssort_range(hdr, offset + lo + 1, count - 1 - lo, ignore_case);

	#undef CK_SSORT_KEY
}

void ck_map_ssort_impl(CK_MapHeader* hdr, int ignore_case)
{
	ck_map_ssort_range(hdr, 0, hdr->size, ignore_case);
}

void ck_map_free_impl(CK_MapHeader* hdr)
{
	CK_FREE(hdr);  // Single allocation.
}

void ck_map_clear_impl(CK_MapHeader* hdr)
{
	hdr->size = 0;
	if (hdr->slot_capacity) ck_map_zero_slots(hdr);
}

//--------------------------------------------------------------------------------------------------
// String interning (originally from Per Vognsen).

uint64_t ck_hash_fnv1a(const void* ptr, size_t sz)
{
	uint64_t x = 0xcbf29ce484222325ull;
	const char* buf = (const char*)ptr;
	for (size_t i = 0; i < sz; i++) {
		x ^= (uint64_t)(unsigned char)buf[i];
		x *= 0x100000001b3ull;
		x ^= x >> 32;
	}
	return x;
}

// C11 atomics for C, std::atomic for C++
// Must close extern "C" before including <atomic> since it contains C++ templates.
#ifdef __cplusplus
} // extern "C"
#include <atomic>
extern "C" {
#define CK_ATOMIC(T) std::atomic<T>
#define ck_atomic_load(p) (p)->load()
#define ck_atomic_store(p, v) (p)->store(v)
#define ck_atomic_compare_exchange_strong(p, expected, desired) (p)->compare_exchange_strong(*(expected), (desired))
#define ck_atomic_compare_exchange_weak(p, expected, desired) (p)->compare_exchange_weak(*(expected), (desired))
#else
#include <stdatomic.h>
#define CK_ATOMIC(T) _Atomic(T)
#define ck_atomic_load(p) atomic_load(p)
#define ck_atomic_store(p, v) atomic_store(p, v)
#define ck_atomic_compare_exchange_strong(p, expected, desired) atomic_compare_exchange_strong(p, expected, desired)
#define ck_atomic_compare_exchange_weak(p, expected, desired) atomic_compare_exchange_weak(p, expected, desired)
#endif

typedef struct CK_InternTable
{
	CK_MAP(CK_UniqueString*) interns;
	CK_ATOMIC(int) lock;
} CK_InternTable;

static CK_ATOMIC(CK_InternTable*) g_intern_table;
static int g_sintern_gen;

int ck_sintern_gen() { return g_sintern_gen; }

static CK_InternTable* ck_sintern_get_table()
{
	CK_InternTable* table = ck_atomic_load(&g_intern_table);
	if (!table) {
		CK_InternTable* new_table = (CK_InternTable*)CK_ALLOC(sizeof(CK_InternTable));
		memset(new_table, 0, sizeof(CK_InternTable));
		CK_InternTable* expected = NULL;
		if (ck_atomic_compare_exchange_strong(&g_intern_table, &expected, new_table)) {
			table = new_table;
		} else {
			CK_FREE(new_table);
			table = expected;
		}
	}
	return table;
}

static void ck_sintern_lock(CK_InternTable* table)
{
	int expected = 0;
	while (!ck_atomic_compare_exchange_weak(&table->lock, &expected, 1)) {
		expected = 0;
	}
}

static void ck_sintern_unlock(CK_InternTable* table)
{
	ck_atomic_store(&table->lock, 0);
}

const char* ck_sintern(const char* s)
{
	// Direct-mapped pointer cache. String literals have stable addresses,
	// so the same call site always hits after the first miss. The strcmp
	// handles reused buffers (e.g. sprintf into the same char[]). The
	// generation counter invalidates the cache after sintern_nuke().
	struct ck_sintern_cache_entry { const char* key; const char* val; };
#ifdef __cplusplus
	static thread_local struct ck_sintern_cache_entry ck_sintern_cache[64];
	static thread_local int ck_sintern_cache_gen;
#else
	static _Thread_local struct ck_sintern_cache_entry ck_sintern_cache[64];
	static _Thread_local int ck_sintern_cache_gen;
#endif
	int gen = ck_sintern_gen();
	if (ck_sintern_cache_gen != gen) {
		memset(ck_sintern_cache, 0, sizeof(ck_sintern_cache));
		ck_sintern_cache_gen = gen;
	}
	unsigned ck_idx = (unsigned)((uintptr_t)s >> 3) & 63;
	if (ck_sintern_cache[ck_idx].key == s && strcmp(ck_sintern_cache[ck_idx].val, s) == 0)
		return ck_sintern_cache[ck_idx].val;
	const char* result = ck_sintern_range(s, s + strlen(s));
	ck_sintern_cache[ck_idx].key = s;
	ck_sintern_cache[ck_idx].val = result;
	return result;
}

const char* ck_sintern_range(const char* start, const char* end)
{
	CK_InternTable* table = ck_sintern_get_table();
	size_t len = (size_t)(end - start);
	uint64_t key = ck_hash_fnv1a((void*)start, len);

	ck_sintern_lock(table);

	CK_UniqueString* head = map_get(table->interns, key);
	for (CK_UniqueString* it = head; it; it = it->next) {
		if ((size_t)it->len == len && memcmp(it->str, start, len) == 0) {
			ck_sintern_unlock(table);
			return it->str;
		}
	}

	size_t bytes = sizeof(CK_UniqueString) + len + 1;
	CK_UniqueString* node = (CK_UniqueString*)CK_ALLOC(bytes);
	node->cookie.val = CK_INTERN_COOKIE;
	node->len = (int)len;
	node->next = head;
	node->str = (char*)(node + 1);
	memcpy(node->str, start, len);
	node->str[len] = '\0';
	map_set(table->interns, key, node);

	ck_sintern_unlock(table);
	return node->str;
}

void sintern_nuke()
{
	g_sintern_gen++;
	CK_InternTable* table = ck_atomic_load(&g_intern_table);
	if (!table) return;

	ck_sintern_lock(table);
	ck_atomic_store(&g_intern_table, (CK_InternTable*)NULL);

	for (int i = 0; i < map_size(table->interns); ++i) {
		CK_UniqueString* it = map_items(table->interns)[i];
		while (it) {
			CK_UniqueString* next = it->next;
			CK_FREE(it);
			it = next;
		}
	}
	map_free(table->interns);
	ck_sintern_unlock(table);
	CK_FREE(table);
}

const uint16_t* ck_decode_UTF16(const uint16_t* s, int* codepoint)
{
	int W1 = *s++;
	if (W1 < 0xD800 || W1 > 0xDFFF) {
		*codepoint = W1;
	} else if (W1 > 0xD800 && W1 < 0xDBFF) {
		int W2 = *s++;
		if (W2 > 0xDC00 && W2 < 0xDFFF) *codepoint = 0x10000 + (((W1 & 0x03FF) << 10) | (W2 & 0x03FF));
		else *codepoint = 0xFFFD;
	} else *codepoint = 0xFFFD;
	return s;
}

#ifdef __cplusplus
} // extern "C"
#endif

#endif // CKIT_IMPLEMENTATION_GUARD
#endif // CKIT_IMPLEMENTATION

// ---- src/nudge.c ----

// See LICENSE for licensing info.
// nudge.c -- physics world implementation


// ---- src/nudge_internal.h ----

// nudge_internal.h -- internal types for nudge physics engine (unity build)
#ifndef NUDGE_INTERNAL_H
#define NUDGE_INTERNAL_H

// Portable float validity: rejects NaN, inf, and extreme magnitudes.
// NaN fails the self-equality test; inf and huge values fail the bound test.
static inline int float_valid(float f) { return f == f && f > -1e18f && f < 1e18f; }
static inline int v3_is_valid(v3 v) { return float_valid(v.x) && float_valid(v.y) && float_valid(v.z); }
static inline int quat_is_valid(quat q) { return float_valid(q.x) && float_valid(q.y) && float_valid(q.z) && float_valid(q.w); }
#define is_valid(x) _Generic((x), float: float_valid, v3: v3_is_valid, quat: quat_is_valid)(x)

// -----------------------------------------------------------------------------
// Internal types: cold (metadata) and hot (solver iteration) splits.

typedef struct ShapeInternal
{
	ShapeType type;
	v3 local_pos;
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule;
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
	};
} ShapeInternal;

// Cold: persistent metadata, topology, rarely touched during solve.
typedef struct BodyCold
{
	float mass;
	CK_DYNA ShapeInternal* shapes;
	int bvh_leaf;     // -1 if not in BVH
	int island_id;    // -1 = no island (static or unconnected)
	int island_prev;  // prev body in island body list, -1 = head
	int island_next;  // next body in island body list, -1 = tail
} BodyCold;

// Hot: solver working set, iterated every step, packed for cache.
typedef struct BodyHot
{
	v3 position;
	quat rotation;
	v3 velocity;
	v3 angular_velocity;
	float inv_mass;
	v3 inv_inertia_local; // diagonal of local-space inverse inertia tensor
	float friction;
	float restitution;
	float linear_damping;
	float angular_damping;
	float sleep_time; // accumulated seconds below velocity threshold
} BodyHot;

typedef struct WarmManifold WarmManifold; // forward decl for warm cache
typedef struct AVBD_WarmManifold AVBD_WarmManifold;

// Joint persistent storage (handle-based, parallel arrays like bodies).
typedef enum JointType { JOINT_BALL_SOCKET, JOINT_DISTANCE } JointType;

typedef struct JointInternal
{
	JointType type;
	int body_a, body_b; // array indices (resolved from handles at creation)
	union {
		struct { v3 local_a, local_b; SpringParams spring; } ball_socket;
		struct { v3 local_a, local_b; float rest_length; SpringParams spring; } distance;
	};
	// Warm starting: accumulated impulses persisted across frames.
	union {
		v3 warm_lambda3;   // ball socket (3 DOF)
		float warm_lambda1; // distance (1 DOF)
	};
	// AVBD warm state (penalty/lambda for augmented Lagrangian)
	v3 avbd_penalty_lin;
	v3 avbd_lambda_lin;
	v3 avbd_C0_lin;        // constraint error at x- for stabilization
	// Island linked list fields.
	int island_id;    // -1 = none
	int island_prev;  // -1 = head
	int island_next;  // -1 = tail
} JointInternal;

typedef struct BVHTree BVHTree; // forward decl, defined in bvh.c

// Island: group of connected bodies that can sleep/wake together.
typedef struct Island
{
	int head_body, tail_body, body_count;
	int head_joint, tail_joint, joint_count;
	int constraint_remove_count;
	int awake; // 1 = awake, 0 = sleeping
} Island;

#define SLEEP_VEL_THRESHOLD   0.01f  // squared velocity magnitude threshold
#define SLEEP_TIME_THRESHOLD  0.5f   // seconds of stillness before sleep
#define LINEAR_SLOP           0.01f  // contact margin: keeps contacts alive near zero-penetration

typedef struct WorldInternal
{
	int frame;         // monotonically increasing frame counter
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts;
	CK_MAP(WarmManifold) warm_cache;
	CK_MAP(AVBD_WarmManifold) avbd_warm_cache;
	CK_DYNA v3* avbd_prev_velocity; // per-body prev velocity for adaptive warm-start
	// Joints
	CK_DYNA JointInternal* joints;
	CK_DYNA uint32_t*      joint_gen;
	CK_DYNA int*           joint_free;
	// Broadphase
	BroadphaseType broadphase_type;
	BVHTree* bvh_static;
	BVHTree* bvh_dynamic;
	// Islands
	CK_DYNA Island*     islands;
	CK_DYNA uint32_t*   island_gen;
	CK_DYNA int*        island_free;
	CK_MAP(uint8_t)     prev_touching; // body_pair_key -> 1 for pairs touching last frame
	int sleep_enabled; // 1 = island sleep active (default)
	FrictionModel friction_model;
	SolverType solver_type;
	int velocity_iters;
	int position_iters;
	float contact_hertz;
	float contact_damping_ratio;
	float max_push_velocity;
	int sub_steps;
	// AVBD parameters
	float avbd_alpha;       // stabilization (0.95-0.99)
	float avbd_beta_lin;    // penalty ramp, linear constraints
	float avbd_beta_ang;    // penalty ramp, angular constraints
	float avbd_gamma;       // warm-start decay
	int avbd_iterations;    // solver iterations
} WorldInternal;

// -----------------------------------------------------------------------------
// Solver types and constants.

#define SOLVER_VELOCITY_ITERS      10
#define SOLVER_POSITION_ITERS      4
#define SOLVER_BAUMGARTE           0.2f   // used by joints only; contacts use NGS
#define SOLVER_SLOP                0.005f // NGS position correction dead zone
#define SOLVER_RESTITUTION_THRESH  1.0f
#define SOLVER_POS_BAUMGARTE       0.2f
#define SOLVER_POS_MAX_CORRECTION  0.2f   // max position correction per step
#define SOLVER_MAX_LINEAR_VEL      500.0f
#define SOLVER_MAX_ANGULAR_VEL     100.0f

typedef struct SolverContact
{
	v3 r_a;              // contact offset from body A
	v3 r_b;              // contact offset from body B
	v3 normal;           // cached contact normal
	v3 tangent1;         // friction direction 1
	v3 tangent2;         // friction direction 2
	float eff_mass_n;    // 1/K for normal row
	float eff_mass_t1;   // 1/K for tangent1 row
	float eff_mass_t2;   // 1/K for tangent2 row
	float bias;          // velocity bias for penetration recovery
	float bounce;        // restitution velocity bias
	float lambda_n;      // accumulated normal impulse (>= 0)
	float lambda_t1;     // accumulated tangent1 impulse
	float lambda_t2;     // accumulated tangent2 impulse
	float softness;      // soft constraint regularization (0 = rigid/NGS)
	float penetration;   // cached for position correction pass
	uint32_t feature_id; // geometric feature key for warm starting
} SolverContact;

typedef struct SolverManifold
{
	int body_a, body_b;
	int contact_start;
	int contact_count;
	float friction;
	// Manifold-level patch friction data (FRICTION_PATCH only)
	v3 centroid_r_a;
	v3 centroid_r_b;
	v3 tangent1, tangent2;
	v3 normal;
	float eff_mass_t1, eff_mass_t2;
	float eff_mass_twist;
	float lambda_t1, lambda_t2;
	float lambda_twist;
	float patch_area;
	float patch_radius;
	// Block solver: normal coupling matrix A[i][j] (symmetric, max 4x4 = 10 unique)
	float K_nn[16];      // full NxN stored row-major (max 4x4)
	float K_nn_inv[16];  // inverse of K_nn
} SolverManifold;

// Warm starting: cached impulses from previous frame, keyed by body pair.
typedef struct WarmContact
{
	uint32_t feature_id; // geometric feature key for matching
	v3 r_a;              // body-A-relative position for spatial fallback matching
	float lambda_n;
	float lambda_t1;
	float lambda_t2;
} WarmContact;

struct WarmManifold
{
	WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale; // 0 = updated this frame, incremented each frame not touched, evicted at >1
	// Manifold-level friction warm data (FRICTION_PATCH)
	float manifold_lambda_t1;
	float manifold_lambda_t2;
	float manifold_lambda_twist;
};

// -----------------------------------------------------------------------------
// Joint solver types.

typedef struct SolverBallSocket
{
	int body_a, body_b;
	v3 r_a, r_b;
	float eff_mass[6]; // symmetric 3x3 inverse (xx,xy,xz,yy,yz,zz)
	v3 bias;
	float softness;
	v3 lambda;
	int joint_idx;
} SolverBallSocket;

typedef struct SolverDistance
{
	int body_a, body_b;
	v3 r_a, r_b;
	v3 axis;
	float eff_mass;
	float bias;
	float softness;
	float lambda;
	int joint_idx;
} SolverDistance;

// Constraint ref for graph coloring dispatch.
enum { CTYPE_CONTACT, CTYPE_BALL_SOCKET, CTYPE_DISTANCE };

typedef struct ConstraintRef
{
	uint8_t type;
	uint8_t color;
	int index;
	int body_a, body_b;
} ConstraintRef;

// -----------------------------------------------------------------------------
// AVBD (Augmented Vertex Block Descent) solver types.

#define AVBD_STABLE_THRESH 0.05f // adaptive alpha: below this error, use full stabilization
#define AVBD_PENALTY_MIN  1.0f
#define AVBD_PENALTY_MAX  1e10f
#define AVBD_MARGIN       0.0f   // nudge contacts already have margin via LINEAR_SLOP
#define AVBD_STICK_THRESH 0.00001f

typedef struct AVBD_Contact
{
	v3 r_a, r_b;       // contact offsets in body-local space
	v3 C0;             // constraint error at x- (normal, tangent1, tangent2)
	v3 penalty;        // per-axis penalty parameters (ramp over iterations)
	v3 lambda;         // accumulated dual variables (clamped force)
	int stick;         // static friction flag
	uint32_t feature_id;
} AVBD_Contact;

typedef struct AVBD_Manifold
{
	int body_a, body_b;
	int contact_count;
	float friction;
	m3x3 basis;        // [normal; tangent1; tangent2] row-major
	AVBD_Contact contacts[MAX_CONTACTS];
} AVBD_Manifold;

// CSR adjacency: separate arrays for contacts and joints.
// body i's contact refs at ct_adj[ct_start[i]..ct_start[i+1]]
// body i's joint refs at jt_adj[jt_start[i]..jt_start[i+1]]

typedef struct AVBD_ContactAdj
{
	int manifold_idx;
	int contact_idx;
	int is_body_a;
} AVBD_ContactAdj;

typedef struct AVBD_JointAdj
{
	int joint_idx;
	int is_body_a;
} AVBD_JointAdj;

// AVBD warm cache: persists penalty/lambda across frames (keyed by body pair).
typedef struct AVBD_WarmContact
{
	uint32_t feature_id;
	v3 r_a;             // body-local contact offset A (also used for spatial matching)
	v3 r_b;             // body-local contact offset B
	v3 penalty;
	v3 lambda;
	int stick;
} AVBD_WarmContact;

typedef struct AVBD_WarmManifold
{
	AVBD_WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale;
} AVBD_WarmManifold;

// Per-body temporary state for AVBD (not stored in BodyHot)
typedef struct AVBD_BodyState
{
	v3 inertial_lin;
	v3 initial_lin;    // x- for velocity recovery
	quat inertial_ang;
	quat initial_ang;  // q- for angular velocity recovery
} AVBD_BodyState;

#endif // NUDGE_INTERNAL_H

// ---- src/gjk.c ----

// See LICENSE for licensing info.
// gjk.c -- GJK distance algorithm for 3D convex shapes.

// -----------------------------------------------------------------------------
// GJK types.

typedef struct GJK_Vertex
{
	v3 point1;   // support point on shape A (world space)
	v3 point2;   // support point on shape B (world space)
	v3 point;    // Minkowski difference: point2 - point1
	float u;     // barycentric coordinate
	int index1;  // support index on A
	int index2;  // support index on B
} GJK_Vertex;

typedef struct GJK_Simplex
{
	GJK_Vertex v[4];
	float divisor;
	int count;
} GJK_Simplex;

typedef struct GJK_Result
{
	v3 point1;     // closest point on shape A
	v3 point2;     // closest point on shape B
	float distance;
	int iterations;
} GJK_Result;

// GJK shape: tagged union for support function dispatch.
// Support operates on UNSCALED local-space vertices. The caller transforms
// the result to world space. For hulls, scale is baked into the vertex lookup.
enum { GJK_POINT, GJK_SEGMENT, GJK_HULL };

typedef struct GJK_Shape
{
	int type;
	int count;
	const v3* verts;  // pointer to vertex array (local space)
	v3 verts_buf[2];  // inline storage for point/segment
} GJK_Shape;

static int gjk_get_support(const GJK_Shape* s, v3 dir)
{
	int best = 0;
	float best_d = dot(s->verts[0], dir);
	for (int i = 1; i < s->count; i++) {
		float d = dot(s->verts[i], dir);
		if (d > best_d) { best_d = d; best = i; }
	}
	return best;
}

// Positioned shape: GJK_Shape + transform for world-space queries.
typedef struct GJK_Proxy
{
	GJK_Shape shape;
	v3 pos;
	quat rot;
} GJK_Proxy;

// Proxy init functions: take output pointer to avoid dangling self-referential
// pointers when returning structs by value (verts -> verts_buf).
static void gjk_proxy_point(GJK_Proxy* pr, v3 p)
{
	*pr = (GJK_Proxy){0};
	pr->shape.type = GJK_POINT;
	pr->shape.count = 1;
	pr->shape.verts_buf[0] = p;
	pr->shape.verts = pr->shape.verts_buf;
	pr->pos = V3(0,0,0);
	pr->rot = quat_identity();
}

static void gjk_proxy_segment(GJK_Proxy* pr, v3 p, v3 q)
{
	*pr = (GJK_Proxy){0};
	pr->shape.type = GJK_SEGMENT;
	pr->shape.count = 2;
	pr->shape.verts_buf[0] = p;
	pr->shape.verts_buf[1] = q;
	pr->shape.verts = pr->shape.verts_buf;
	pr->pos = V3(0,0,0);
	pr->rot = quat_identity();
}

// For hulls: pre-scale the vertices into a temp buffer, store as proxy.
// Caller must keep scaled_verts alive for the duration of the GJK call.
static void gjk_proxy_hull(GJK_Proxy* pr, const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts)
{
	for (int i = 0; i < hull->vert_count; i++)
		scaled_verts[i] = hmul(hull->verts[i], sc);
	*pr = (GJK_Proxy){0};
	pr->shape.type = GJK_HULL;
	pr->shape.count = hull->vert_count;
	pr->shape.verts = scaled_verts;
	pr->pos = pos;
	pr->rot = rot;
}

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.
// Ported directly from lm engine (lmSimplex::Solve2/3/4).

static int gjk_solve2(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point;
	float u = dot(b, sub(b, a));
	float v = dot(a, sub(a, b));
	if (v <= 0.0f) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (u <= 0.0f) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	s->divisor = u + v;
	if (s->divisor == 0.0f) return 0; // degenerate edge (zero length)
	s->v[0].u = u; s->v[1].u = v; s->count = 2;
	return 1;
}

static int gjk_solve3(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point;
	float uAB = dot(b, sub(b, a)), vAB = dot(a, sub(a, b));
	float uBC = dot(c, sub(c, b)), vBC = dot(b, sub(b, c));
	float uCA = dot(a, sub(a, c)), vCA = dot(c, sub(c, a));
	v3 n = cross(sub(b, a), sub(c, a));
	float uABC = dot(cross(b, c), n), vABC = dot(cross(c, a), n), wABC = dot(cross(a, b), n);

	if (vAB <= 0.0f && uCA <= 0.0f) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uAB <= 0.0f && vBC <= 0.0f) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uBC <= 0.0f && vCA <= 0.0f) { s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uAB > 0.0f && vAB > 0.0f && wABC <= 0.0f) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	if (uBC > 0.0f && vBC > 0.0f && uABC <= 0.0f) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	if (uCA > 0.0f && vCA > 0.0f && vABC <= 0.0f) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }
	s->divisor = uABC + vABC + wABC;
	if (s->divisor == 0.0f) return 0; // degenerate triangle (zero area)
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->count = 3;
	return 1;
}

static inline float stp(v3 a, v3 b, v3 c) { return dot(a, cross(b, c)); }

static int gjk_solve4(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point, d = s->v[3].point;

	// Edge barycentric coords
	float uAB = dot(b, sub(b, a)), vAB = dot(a, sub(a, b));
	float uBC = dot(c, sub(c, b)), vBC = dot(b, sub(b, c));
	float uCA = dot(a, sub(a, c)), vCA = dot(c, sub(c, a));
	float uBD = dot(d, sub(d, b)), vBD = dot(b, sub(b, d));
	float uDC = dot(c, sub(c, d)), vDC = dot(d, sub(d, c));
	float uAD = dot(d, sub(d, a)), vAD = dot(a, sub(a, d));

	// Face barycentric coords
	v3 n;
	n = cross(sub(d, a), sub(b, a));
	float uADB = dot(cross(d, b), n), vADB = dot(cross(b, a), n), wADB = dot(cross(a, d), n);
	n = cross(sub(c, a), sub(d, a));
	float uACD = dot(cross(c, d), n), vACD = dot(cross(d, a), n), wACD = dot(cross(a, c), n);
	n = cross(sub(b, c), sub(d, c));
	float uCBD = dot(cross(b, d), n), vCBD = dot(cross(d, c), n), wCBD = dot(cross(c, b), n);
	n = cross(sub(b, a), sub(c, a));
	float uABC = dot(cross(b, c), n), vABC = dot(cross(c, a), n), wABC = dot(cross(a, b), n);

	// Tetrahedron volume coords
	float denom = stp(sub(c, b), sub(a, b), sub(d, b));
	if (denom == 0.0f) return 0; // degenerate tetrahedron (zero volume)
	float vol = 1.0f / denom;
	float uABCD = stp(c, d, b) * vol, vABCD = stp(c, a, d) * vol;
	float wABCD = stp(d, a, b) * vol, xABCD = stp(b, a, c) * vol;

	// Vertex regions
	if (vAB <= 0 && uCA <= 0 && vAD <= 0) { s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uAB <= 0 && vBC <= 0 && vBD <= 0) { s->v[0] = s->v[1]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uBC <= 0 && vCA <= 0 && uDC <= 0) { s->v[0] = s->v[2]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uBD <= 0 && vDC <= 0 && uAD <= 0) { s->v[0] = s->v[3]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }

	// Edge regions
	if (wABC <= 0 && vADB <= 0 && uAB > 0 && vAB > 0) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	if (uABC <= 0 && wCBD <= 0 && uBC > 0 && vBC > 0) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	if (vABC <= 0 && wACD <= 0 && uCA > 0 && vCA > 0) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }
	if (vCBD <= 0 && uACD <= 0 && uDC > 0 && vDC > 0) { s->v[0] = s->v[3]; s->v[1] = s->v[2]; s->v[0].u = uDC; s->v[1].u = vDC; s->divisor = uDC + vDC; s->count = 2; return 1; }
	if (vACD <= 0 && wADB <= 0 && uAD > 0 && vAD > 0) { s->v[1] = s->v[3]; s->v[0].u = uAD; s->v[1].u = vAD; s->divisor = uAD + vAD; s->count = 2; return 1; }
	if (uCBD <= 0 && uADB <= 0 && uBD > 0 && vBD > 0) { s->v[0] = s->v[1]; s->v[1] = s->v[3]; s->v[0].u = uBD; s->v[1].u = vBD; s->divisor = uBD + vBD; s->count = 2; return 1; }

	// Face regions
	if (xABCD <= 0 && uABC > 0 && vABC > 0 && wABC > 0) { s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->divisor = uABC + vABC + wABC; s->count = 3; return 1; }
	if (uABCD <= 0 && uCBD > 0 && vCBD > 0 && wCBD > 0) { s->v[0] = s->v[2]; s->v[2] = s->v[3]; s->v[0].u = uCBD; s->v[1].u = vCBD; s->v[2].u = wCBD; s->divisor = uCBD + vCBD + wCBD; s->count = 3; return 1; }
	if (vABCD <= 0 && uACD > 0 && vACD > 0 && wACD > 0) { s->v[1] = s->v[2]; s->v[2] = s->v[3]; s->v[0].u = uACD; s->v[1].u = vACD; s->v[2].u = wACD; s->divisor = uACD + vACD + wACD; s->count = 3; return 1; }
	if (wABCD <= 0 && uADB > 0 && vADB > 0 && wADB > 0) { s->v[2] = s->v[1]; s->v[1] = s->v[3]; s->v[0].u = uADB; s->v[1].u = vADB; s->v[2].u = wADB; s->divisor = uADB + vADB + wADB; s->count = 3; return 1; }

	// Interior
	s->v[0].u = uABCD; s->v[1].u = vABCD; s->v[2].u = wABCD; s->v[3].u = xABCD;
	s->divisor = 1.0f; s->count = 4;
	return 1;
}

// -----------------------------------------------------------------------------
// Closest point and witness points.

static v3 gjk_closest_point(const GJK_Simplex* s)
{
	float inv = 1.0f / s->divisor;
	switch (s->count) {
	case 1: return s->v[0].point;
	case 2: return add(scale(s->v[0].point, s->v[0].u * inv), scale(s->v[1].point, s->v[1].u * inv));
	case 3: return add(add(scale(s->v[0].point, s->v[0].u * inv), scale(s->v[1].point, s->v[1].u * inv)), scale(s->v[2].point, s->v[2].u * inv));
	case 4: return V3(0,0,0);
	}
	return V3(0,0,0);
}

static void gjk_witness_points(const GJK_Simplex* s, v3* p1, v3* p2)
{
	float inv = 1.0f / s->divisor;
	*p1 = V3(0,0,0); *p2 = V3(0,0,0);
	for (int i = 0; i < s->count; i++) {
		*p1 = add(*p1, scale(s->v[i].point1, s->v[i].u * inv));
		*p2 = add(*p2, scale(s->v[i].point2, s->v[i].u * inv));
	}
}

static v3 gjk_search_dir(const GJK_Simplex* s)
{
	switch (s->count) {
	case 1: return neg(s->v[0].point);
	case 2: {
		v3 ba = sub(s->v[0].point, s->v[1].point);
		return cross(cross(ba, neg(s->v[0].point)), ba);
	}
	case 3: {
		v3 n = cross(sub(s->v[1].point, s->v[0].point), sub(s->v[2].point, s->v[0].point));
		if (dot(n, s->v[0].point) <= 0.0f) return n;
		return neg(n);
	}
	}
	return V3(0,0,0);
}

// -----------------------------------------------------------------------------
// Main GJK distance function (lm engine pattern).

#define GJK_MAX_ITERS 20

static GJK_Result gjk_distance_ex(const GJK_Proxy* proxyA, const GJK_Proxy* proxyB)
{
	GJK_Result result = {0};
	GJK_Simplex simplex = {0};

	// Initialize with first support point.
	simplex.v[0].index1 = 0;
	simplex.v[0].index2 = 0;
	simplex.v[0].point1 = add(proxyA->pos, rotate(proxyA->rot, proxyA->shape.verts[0]));
	simplex.v[0].point2 = add(proxyB->pos, rotate(proxyB->rot, proxyB->shape.verts[0]));
	simplex.v[0].point = sub(simplex.v[0].point2, simplex.v[0].point1);
	simplex.v[0].u = 1.0f;
	simplex.divisor = 1.0f;
	simplex.count = 1;

	int save1[4], save2[4];
	float dsq0 = FLT_MAX;
	int iter = 0;

	while (iter < GJK_MAX_ITERS) {
		int save_count = simplex.count;
		for (int i = 0; i < save_count; i++) {
			save1[i] = simplex.v[i].index1;
			save2[i] = simplex.v[i].index2;
		}

		GJK_Simplex backup = simplex;
		int solved = 1;
		switch (simplex.count) {
		case 1: break;
		case 2: solved = gjk_solve2(&simplex); break;
		case 3: solved = gjk_solve3(&simplex); break;
		case 4: solved = gjk_solve4(&simplex); break;
		}

		if (!solved) { simplex = backup; break; }
		if (simplex.count == 4) break;

		v3 p = gjk_closest_point(&simplex);
		float dsq1 = len2(p);
		if (dsq1 > dsq0) break;
		dsq0 = dsq1;

		v3 d = gjk_search_dir(&simplex);
		if (len2(d) < FLT_EPSILON * FLT_EPSILON) break;

		// New support point (transform local dir to proxy local space).
		int iA = gjk_get_support(&proxyA->shape, rotate(inv(proxyA->rot), neg(d)));
		v3 sA = add(proxyA->pos, rotate(proxyA->rot, proxyA->shape.verts[iA]));
		int iB = gjk_get_support(&proxyB->shape, rotate(inv(proxyB->rot), d));
		v3 sB = add(proxyB->pos, rotate(proxyB->rot, proxyB->shape.verts[iB]));
		iter++;

		// Duplicate check
		int dup = 0;
		for (int i = 0; i < save_count; i++)
			if (iA == save1[i] && iB == save2[i]) { dup = 1; break; }
		if (dup) break;

		GJK_Vertex* v = &simplex.v[simplex.count];
		v->index1 = iA; v->point1 = sA;
		v->index2 = iB; v->point2 = sB;
		v->point = sub(sB, sA);
		simplex.count++;
	}

	gjk_witness_points(&simplex, &result.point1, &result.point2);
	result.distance = len(sub(result.point2, result.point1));
	result.iterations = iter;
	return result;
}

// ---- src/quickhull.c ----

// See LICENSE for licensing info.
#include <float.h>
// quickhull.c -- 3D Quickhull convex hull algorithm.
//
// Builds a convex hull from a point cloud using conflict lists,
// recursive horizon computation, and two-pass non-convex face merging.
// Degenerate faces (triangle collapse, vertex skip) are handled
// on-the-fly during the merge splice to maintain manifold topology.
//
// Produces a Hull with polygonal faces, half-edge mesh, and face planes.
//
// Edge convention: each half-edge stores its ORIGIN vertex (= tail).
// The head of edge e is edges[edges[e].next].origin.
// We use ORIGIN (tail) vertex convention for half-edges.
// for a triangle created with createTriangle(v0, v1, v2).

#define QH_INVALID    (~0)
#define QH_VISIBLE    1
#define QH_NON_CONVEX 2
#define QH_DELETED    3
#define QH_EDGE_DELETED 0x01

// Debug state: set qh_debug=1 to enable runtime logging.
static int        qh_debug;
static const v3*  qh_dbg_points;
static int        qh_dbg_count;

#define QH_DEBUG(...) do { if (qh_debug) { __VA_ARGS__; } } while(0)

#define QH_ASSERT(cond, msg) do { \
	if (!(cond)) { \
		fprintf(stderr, "quickhull: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
		fprintf(stderr, "  input: %d points\n  v3 pts[] = {\n", qh_dbg_count); \
		for (int _i = 0; _i < qh_dbg_count; _i++) \
			fprintf(stderr, "    {%.9ef, %.9ef, %.9ef},\n", \
				qh_dbg_points[_i].x, qh_dbg_points[_i].y, qh_dbg_points[_i].z); \
		fprintf(stderr, "  };\n"); \
		exit(-1); \
	} \
} while(0)

// -----------------------------------------------------------------------------
// Internal types.

typedef struct QH_Vertex {
	v3 pos;
	int conflict_next;
	int conflict_prev;
} QH_Vertex;

typedef struct QH_Edge {
	int next, prev, twin;
	int origin, face, mark;
} QH_Edge;

typedef struct QH_Face {
	int edge;
	int next, prev;
	int conflict_head;
	int mark;
	int num_verts;
	float area;
	float maxoutside; // max distance any point was ever seen outside this face
	HullPlane plane;
} QH_Face;

typedef struct QH_State {
	CK_DYNA QH_Vertex* verts;
	CK_DYNA QH_Edge*   edges;
	CK_DYNA QH_Face*   faces;
	int edge_free, face_free;
	v3 interior;
	float epsilon;
} QH_State;

typedef struct QH_FaceList {
	int first, last;
} QH_FaceList;


// -----------------------------------------------------------------------------
// Free list management.

static int qh_alloc_edge(QH_State* s)
{
	if (s->edge_free != QH_INVALID) {
		int idx = s->edge_free;
		s->edge_free = s->edges[idx].next;
		s->edges[idx] = (QH_Edge){0};
		return idx;
	}
	QH_Edge e = {0};
	apush(s->edges, e);
	return asize(s->edges) - 1;
}

static int qh_alloc_face(QH_State* s)
{
	if (s->face_free != QH_INVALID) {
		int idx = s->face_free;
		s->face_free = s->faces[idx].next;
		s->faces[idx] = (QH_Face){0};
		s->faces[idx].conflict_head = QH_INVALID;
		s->faces[idx].mark = QH_VISIBLE;
		s->faces[idx].maxoutside = s->epsilon;
		return idx;
	}
	QH_Face f = {0};
	f.conflict_head = QH_INVALID;
	f.mark = QH_VISIBLE;
	f.maxoutside = s->epsilon;
	apush(s->faces, f);
	return asize(s->faces) - 1;
}

// Forward declaration (used by qh_conflict_add before full definition).
static float qh_plane_dist(HullPlane p, v3 pt);

// -----------------------------------------------------------------------------
// Conflict list (circular doubly-linked via vertex next/prev).

static void qh_conflict_add(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Vertex* v = &s->verts[vi];
	if (f->conflict_head != QH_INVALID) {
		QH_Vertex* head = &s->verts[f->conflict_head];
		int tail = head->conflict_prev;
		v->conflict_next = f->conflict_head;
		v->conflict_prev = tail;
		head->conflict_prev = vi;
		s->verts[tail].conflict_next = vi;
	} else {
		v->conflict_next = vi;
		v->conflict_prev = vi;
	}
	f->conflict_head = vi;
	// Track max distance any point was seen outside this face.
	float d = qh_plane_dist(f->plane, v->pos);
	if (d > f->maxoutside) f->maxoutside = d;
}

static void qh_conflict_remove(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Vertex* v = &s->verts[vi];
	if (v->conflict_next == vi) {
		f->conflict_head = QH_INVALID;
	} else {
		s->verts[v->conflict_next].conflict_prev = v->conflict_prev;
		s->verts[v->conflict_prev].conflict_next = v->conflict_next;
		if (f->conflict_head == vi) f->conflict_head = v->conflict_next;
	}
	v->conflict_next = QH_INVALID;
	v->conflict_prev = QH_INVALID;
}

// Remove all conflict vertices. Returns singly-linked list (via conflict_next),
// terminated by QH_INVALID.
static int qh_conflict_remove_all(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	int head = f->conflict_head;
	if (head == QH_INVALID) return QH_INVALID;
	int last = s->verts[head].conflict_prev;
	s->verts[last].conflict_next = QH_INVALID;
	f->conflict_head = QH_INVALID;
	return head;
}

// -----------------------------------------------------------------------------
// Geometry helpers.

static float qh_plane_dist(HullPlane p, v3 pt)
{
	return dot(p.normal, pt) - p.offset;
}

// Compute face normal, centroid, area via Newell method.
static void qh_recompute_face(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	v3 normal = V3(0,0,0), centroid = V3(0,0,0);
	int count = 0, e = f->edge;
	do {
		v3 cur = s->verts[s->edges[e].origin].pos;
		v3 nxt = s->verts[s->edges[s->edges[e].next].origin].pos;
		normal.x += (cur.y - nxt.y) * (cur.z + nxt.z);
		normal.y += (cur.z - nxt.z) * (cur.x + nxt.x);
		normal.z += (cur.x - nxt.x) * (cur.y + nxt.y);
		centroid = add(centroid, cur);
		count++;
		if (count > 1000) { QH_ASSERT(0, "infinite face loop in qh_recompute_face"); }
		e = s->edges[e].next;
	} while (e != f->edge);

	float a = len(normal);
	if (a > 0) normal = scale(normal, 1.0f / a);
	centroid = scale(centroid, 1.0f / count);
	f->plane = (HullPlane){ normal, dot(normal, centroid) };
	f->num_verts = count;
	f->area = a;
}

static float qh_face_dist(QH_State* s, int fi, v3 pt)
{
	return qh_plane_dist(s->faces[fi].plane, pt);
}

static v3 qh_face_centroid(QH_State* s, int fi)
{
	int e = s->faces[fi].edge, start = e;
	v3 c = V3(0,0,0); int n = 0;
	do { c = add(c, s->verts[s->edges[e].origin].pos); n++; e = s->edges[e].next; } while (e != start);
	return scale(c, 1.0f / n);
}

// Distance of the adjacent face's centroid to this edge's face plane.
static float qh_opp_face_dist(QH_State* s, int ei)
{
	int fi = s->edges[ei].face;
	int ofi = s->edges[s->edges[ei].twin].face;
	return qh_face_dist(s, fi, qh_face_centroid(s, ofi));
}

static int qh_face_vert_count(QH_State* s, int fi)
{
	int e = s->faces[fi].edge, start = e, n = 0;
	do { n++; e = s->edges[e].next; } while (e != start);
	return n;
}

static int qh_opp_face(QH_State* s, int ei) { return s->edges[s->edges[ei].twin].face; }
static int qh_edge_head(QH_State* s, int ei) { return s->edges[s->edges[ei].next].origin; }
static int qh_edge_tail(QH_State* s, int ei) { return s->edges[ei].origin; }

// Get edge at index i from face's he0 (supports negative indices).
static int qh_get_edge(QH_State* s, int fi, int i)
{
	int e = s->faces[fi].edge;
	while (i > 0) { e = s->edges[e].next; i--; }
	while (i < 0) { e = s->edges[e].prev; i++; }
	return e;
}

static void qh_set_opposite(QH_State* s, int a, int b)
{
	s->edges[a].twin = b;
	s->edges[b].twin = a;
}

// Create a triangle face with edges v0->v1->v2->v0.
static int qh_create_triangle(QH_State* s, int v0, int v1, int v2)
{
	int fi = qh_alloc_face(s);
	int e0 = qh_alloc_edge(s), e1 = qh_alloc_edge(s), e2 = qh_alloc_edge(s);
	s->edges[e0] = (QH_Edge){ .next=e1, .prev=e2, .twin=QH_INVALID, .origin=v0, .face=fi };
	s->edges[e1] = (QH_Edge){ .next=e2, .prev=e0, .twin=QH_INVALID, .origin=v1, .face=fi };
	s->edges[e2] = (QH_Edge){ .next=e0, .prev=e1, .twin=QH_INVALID, .origin=v2, .face=fi };
	s->faces[fi].edge = e0;
	s->faces[fi].next = fi;
	s->faces[fi].prev = fi;
	qh_recompute_face(s, fi);
	return fi;
}

// -----------------------------------------------------------------------------
// Build initial tetrahedron from 4 non-coplanar extremal points.

static int qh_build_simplex(QH_State* s, int nv)
{
	// Find extremal vertices on each axis.
	int maxV[3], minV[3];
	for (int i = 0; i < 3; i++) { maxV[i] = minV[i] = 0; }
	v3 mx = s->verts[0].pos, mn = s->verts[0].pos;
	for (int i = 1; i < nv; i++) {
		v3 p = s->verts[i].pos;
		if (p.x > mx.x) { mx.x = p.x; maxV[0] = i; } else if (p.x < mn.x) { mn.x = p.x; minV[0] = i; }
		if (p.y > mx.y) { mx.y = p.y; maxV[1] = i; } else if (p.y < mn.y) { mn.y = p.y; minV[1] = i; }
		if (p.z > mx.z) { mx.z = p.z; maxV[2] = i; } else if (p.z < mn.z) { mn.z = p.z; minV[2] = i; }
	}

	// Pick axis with greatest spread.
	float best = 0; int imax = 0;
	for (int i = 0; i < 3; i++) {
		v3 a = s->verts[maxV[i]].pos, b = s->verts[minV[i]].pos;
		float d = (i==0) ? a.x-b.x : (i==1) ? a.y-b.y : a.z-b.z;
		if (d > best) { best = d; imax = i; }
	}
	if (best <= s->epsilon) return 0;

	int vtx[4];
	vtx[0] = maxV[imax]; vtx[1] = minV[imax];

	// Find vertex farthest from line vtx[0]-vtx[1].
	v3 u01 = norm(sub(s->verts[vtx[1]].pos, s->verts[vtx[0]].pos));
	float maxSqr = 0; v3 nrml = V3(0,0,0); vtx[2] = -1;
	for (int i = 0; i < nv; i++) {
		v3 xp = cross(u01, sub(s->verts[i].pos, s->verts[vtx[0]].pos));
		float ls = len2(xp);
		if (ls > maxSqr && i != vtx[0] && i != vtx[1]) { maxSqr = ls; vtx[2] = i; nrml = xp; }
	}
	if (vtx[2] < 0 || sqrtf(maxSqr) <= 100*s->epsilon) return 0;
	nrml = norm(nrml);
	nrml = norm(sub(nrml, scale(u01, dot(nrml, u01)))); // orthogonalize

	// Find vertex farthest from plane through vtx[0..2].
	float d0 = dot(s->verts[vtx[2]].pos, nrml);
	float maxDist = 0; vtx[3] = -1;
	for (int i = 0; i < nv; i++) {
		float d = fabsf(dot(s->verts[i].pos, nrml) - d0);
		if (d > maxDist && i != vtx[0] && i != vtx[1] && i != vtx[2]) { maxDist = d; vtx[3] = i; }
	}
	if (vtx[3] < 0 || maxDist <= 100*s->epsilon) return 0;

	// Build 4 triangle faces with correct winding.
	// With origin-vertex convention, edge indices within a triangle: e0, e1, e2.
	int tris[4];
	if (dot(s->verts[vtx[3]].pos, nrml) - d0 < 0) {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[1], vtx[2]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[1], vtx[0]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[2], vtx[1]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[0], vtx[2]);
		for (int i = 0; i < 3; i++) {
			int k = (i+1)%3;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],0), qh_get_edge(s,tris[k+1],2));
			int our_k = (k==0)?2:(k==1)?0:1;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],1), qh_get_edge(s,tris[0],our_k));
		}
	} else {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[2], vtx[1]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[0], vtx[1]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[1], vtx[2]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[2], vtx[0]);
		for (int i = 0; i < 3; i++) {
			int k = (i+1)%3;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],2), qh_get_edge(s,tris[k+1],0));
			int j = (3-i)%3, our_j = (j==0)?2:(j==1)?0:1;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],1), qh_get_edge(s,tris[0],our_j));
		}
	}

	s->interior = scale(add(add(s->verts[vtx[0]].pos, s->verts[vtx[1]].pos), add(s->verts[vtx[2]].pos, s->verts[vtx[3]].pos)), 0.25f);

	// Assign conflict vertices: each point goes to the face it's furthest outside of.
	for (int i = 0; i < nv; i++) {
		if (i==vtx[0]||i==vtx[1]||i==vtx[2]||i==vtx[3]) {
			continue;
		}
		float bd = s->epsilon; int bf = -1;
		for (int k = 0; k < 4; k++) {
			float d = qh_face_dist(s, tris[k], s->verts[i].pos);
			if (d > bd) { bd = d; bf = tris[k]; }
		}
		if (bf >= 0) qh_conflict_add(s, bf, i);
	}
	return 1;
}

// -----------------------------------------------------------------------------
// Find next conflict vertex (globally furthest from any face).

static int qh_next_conflict(QH_State* s, int* out_face)
{
	int bv = QH_INVALID, bf = QH_INVALID; float bd = -1e18f;
	for (int fi = 0; fi < asize(s->faces); fi++) {
		if (s->faces[fi].mark != QH_VISIBLE) continue;
		int head = s->faces[fi].conflict_head;
		if (head == QH_INVALID) continue;
		int ci = head;
		do {
			float d = qh_face_dist(s, fi, s->verts[ci].pos);
			if (d > bd) { bd = d; bv = ci; bf = fi; }
			ci = s->verts[ci].conflict_next;
		} while (ci != head);
	}
	*out_face = bf;
	return bv;
}

// -----------------------------------------------------------------------------
// Move conflict vertices from a deleted face to unclaimed or an absorbing face.

static void qh_delete_face_points(QH_State* s, int fi, int absorb, int* unclaimed)
{
	int vlist = qh_conflict_remove_all(s, fi);
	if (vlist == QH_INVALID) return;
	int v = vlist;
	while (v != QH_INVALID) {
		int nxt = s->verts[v].conflict_next;
		if (absorb != QH_INVALID && qh_face_dist(s, absorb, s->verts[v].pos) > s->epsilon) {
			qh_conflict_add(s, absorb, v);
		} else {
			s->verts[v].conflict_next = *unclaimed;
			*unclaimed = v;
		}
		v = nxt;
	}
}

// -----------------------------------------------------------------------------
// Recursive DFS from a visible face, marking visible faces DELETED.

static void qh_calculate_horizon(QH_State* s, v3 eye, int edge0, int fi, CK_DYNA int** horizon, int* unclaimed)
{
	qh_delete_face_points(s, fi, QH_INVALID, unclaimed);
	s->faces[fi].mark = QH_DELETED;

	int edge;
	if (edge0 == QH_INVALID) { edge0 = s->faces[fi].edge; edge = edge0; }
	else { edge = s->edges[edge0].next; }

	do {
		int ofi = qh_opp_face(s, edge);
		if (s->faces[ofi].mark == QH_VISIBLE) {
			if (qh_face_dist(s, ofi, eye) > s->epsilon)
				qh_calculate_horizon(s, eye, s->edges[edge].twin, ofi, horizon, unclaimed);
			else
				apush(*horizon, edge);
		}
		edge = s->edges[edge].next;
	} while (edge != edge0);
}

// -----------------------------------------------------------------------------
// Create a triangle face from the eye vertex to a horizon edge,
// twin bottom edge with horizon's old twin, return side edge (head->eye).

static int qh_add_adjoining_face(QH_State* s, int eye, int horizon_edge)
{
	int tail = qh_edge_tail(s, horizon_edge);
	int head = qh_edge_head(s, horizon_edge);
	int fi = qh_create_triangle(s, eye, tail, head);
	// e1 (tail->head) twins with the old opposite of the horizon edge.
	qh_set_opposite(s, qh_get_edge(s, fi, 1), s->edges[horizon_edge].twin);
	// Return e2 (head->eye) as the side edge.
	return qh_get_edge(s, fi, 2);
}

// -----------------------------------------------------------------------------
// Splice two edges that should be consecutive after removing a shared boundary.
// Handles degenerate case: if both edges' twins point to the same face, that
// face has become a "fin". Triangle fins are deleted; polygon fins lose a vertex.
// Returns index of discarded face, or QH_INVALID.
//
// Because we store origin explicitly per-edge, the surviving edge must absorb
// the removed edge's origin when splicing.

static int qh_connect_half_edges(QH_State* s, int ep, int en)
{
	int opp_p = qh_opp_face(s, ep);
	int opp_n = qh_opp_face(s, en);
	int this_face = s->edges[en].face;

	if (opp_p == opp_n) {
		// Degenerate: both edges border the same adjacent face.
		int fin_face = opp_n;
		int discarded = QH_INVALID;
		int new_twin;
		QH_DEBUG(fprintf(stderr, "[qh] connect_degenerate: ep=%d en=%d fin=%d nverts=%d this=%d\n",
			ep, en, fin_face, qh_face_vert_count(s, fin_face), this_face));

		if (ep == s->faces[this_face].edge)
			s->faces[this_face].edge = en;

		if (qh_face_vert_count(s, fin_face) == 3) {
			// Triangle collapse: delete the fin face entirely.
			// Find the third edge (not twin of ep or en) whose twin survives.
			int ht = s->edges[en].twin;
			int hp = s->edges[ep].twin;
			int e3 = s->edges[ht].next;
			if (e3 == hp) e3 = s->edges[ht].prev;
			new_twin = s->edges[e3].twin;
			s->faces[fin_face].mark = QH_DELETED;
			discarded = fin_face;
		} else {
			// Polygon fin should have been pre-merged by the caller.
			// If we still get here, use normal linkage as fallback.
			s->edges[ep].next = en;
			s->edges[en].prev = ep;
			return QH_INVALID;
		}

		// Remove ep from this face's loop; en absorbs ep's origin.
		s->edges[en].origin = s->edges[ep].origin;
		int pp = s->edges[ep].prev;
		QH_DEBUG(fprintf(stderr, "[qh]   splice out ep=%d: pp=%d->en=%d (en.next=%d)\n",
			ep, pp, en, s->edges[en].next));
		s->edges[en].prev = pp;
		s->edges[pp].next = en;

		qh_set_opposite(s, en, new_twin);
		qh_recompute_face(s, s->edges[new_twin].face);
		return discarded;
	} else {
		// Normal case: just link consecutively.
		s->edges[ep].next = en;
		s->edges[en].prev = ep;
		return QH_INVALID;
	}
}

// -----------------------------------------------------------------------------
// Merge the adjacent face (across a shared edge) into this face.

static int qh_merge_adjacent_face(QH_State* s, int hedge_adj, int discarded[], int* unclaimed)
{
	int ofi = qh_opp_face(s, hedge_adj);
	int tfi = s->edges[hedge_adj].face;
	int nd = 0;
	discarded[nd++] = ofi;
	QH_DEBUG(fprintf(stderr, "[qh] merge: face %d (mark=%d) absorbed into %d (mark=%d) via edge %d\n",
		ofi, s->faces[ofi].mark, tfi, s->faces[tfi].mark, hedge_adj));
	s->faces[ofi].mark = QH_DELETED;

	int ho = s->edges[hedge_adj].twin;
	int ap = s->edges[hedge_adj].prev, an = s->edges[hedge_adj].next;
	int op = s->edges[ho].prev, on = s->edges[ho].next;

	// Walk past multiply-shared edges.
	while (qh_opp_face(s, ap) == ofi) { ap = s->edges[ap].prev; on = s->edges[on].next; }
	while (qh_opp_face(s, an) == ofi) { op = s->edges[op].prev; an = s->edges[an].next; }

	// Reassign the absorbed face's edges to this face.
	for (int e = on; e != s->edges[op].next; e = s->edges[e].next) s->edges[e].face = tfi;

	if (hedge_adj == s->faces[tfi].edge) s->faces[tfi].edge = an;

	int df = qh_connect_half_edges(s, op, an);
	if (df != QH_INVALID) discarded[nd++] = df;
	df = qh_connect_half_edges(s, ap, on);
	if (df != QH_INVALID) discarded[nd++] = df;


	// Inherit maxoutside from absorbed face.
	if (s->faces[ofi].maxoutside > s->faces[tfi].maxoutside)
		s->faces[tfi].maxoutside = s->faces[ofi].maxoutside;

	// Reassign all edges on the surviving face's loop to tfi, and remove
	// self-edges (both half-edges on the same face) left by degenerate merges.
	{
		int e = s->faces[tfi].edge, start = e, n = 0;
		do {
			n++;
			QH_ASSERT(n < 1000, "face loop corrupted after merge");
			s->edges[e].face = tfi;
			e = s->edges[e].next;
		} while (e != start);

		// Remove stale boundary edges: edges whose twins are on deleted faces
		// are remnants of the shared boundary that cascade splicing didn't clean.
		{
			int cleaned = 1;
			while (cleaned) {
				cleaned = 0;
				e = s->faces[tfi].edge;
				start = e;
				do {
					int tw = s->edges[e].twin;
					int twf = s->edges[tw].face;
					if (s->faces[twf].mark == QH_DELETED || twf == tfi) {
						// Stale or self-edge: splice out.
						int p = s->edges[e].prev;
						int nx = s->edges[e].next;
						s->edges[p].next = nx;
						s->edges[nx].prev = p;
						s->edges[nx].origin = s->edges[e].origin;
						if (s->faces[tfi].edge == e)
							s->faces[tfi].edge = nx;
						cleaned = 1;
						break;
					}
					e = s->edges[e].next;
				} while (e != start);
			}
		}


		QH_DEBUG({
			e = s->faces[tfi].edge; start = e; n = 0;
			do { n++; e = s->edges[e].next; } while (e != start);
			fprintf(stderr, "[qh]   face %d loop (%d edges):", tfi, n);
			e = start;
			do { fprintf(stderr, " %d", e); e = s->edges[e].next; } while (e != start);
			fprintf(stderr, "\n");
		});
	}

	qh_recompute_face(s, tfi);
	return nd;
}

// -----------------------------------------------------------------------------
// Scan face edges for non-convex neighbors and merge one pair. Returns 1 if merged.

#define QH_MERGE_LARGE 1
#define QH_MERGE_ANY   2

static int qh_do_adjacent_merge(QH_State* s, int fi, int type, int* unclaimed)
{
	if (s->faces[fi].mark != QH_VISIBLE) return 0;
	int hedge = s->faces[fi].edge;
	int convex = 1;
	do {
		int ofi = qh_opp_face(s, hedge);
		if (ofi == fi || s->faces[ofi].mark == QH_DELETED) {
			hedge = s->edges[hedge].next; continue;
		}

		float tol = s->epsilon;

		int merge = 0;
		if (type == QH_MERGE_ANY) {
			if (qh_opp_face_dist(s, hedge) > -tol ||
			    qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) merge = 1;
		} else {
			if (s->faces[fi].area > s->faces[ofi].area) {
				if (qh_opp_face_dist(s, hedge) > -tol) merge = 1;
				else if (qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) convex = 0;
			} else {
				if (qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) merge = 1;
				else if (qh_opp_face_dist(s, hedge) > -tol) convex = 0;
			}
		}
		if (merge) {
			QH_DEBUG(fprintf(stderr, "[qh] do_adjacent_merge: fi=%d hedge=%d ofi=%d hedge.face=%d\n",
				fi, hedge, ofi, s->edges[hedge].face));
			// Check if the merge will hit a polygon fin degenerate.
			// If so, skip it to avoid infinite loops.
			int will_poly_fin = 0;
			{
				int ho = s->edges[hedge].twin;
				int ap = s->edges[hedge].prev, an = s->edges[hedge].next;
				int op = s->edges[ho].prev, on = s->edges[ho].next;
				while (qh_opp_face(s, ap) == ofi) { ap = s->edges[ap].prev; on = s->edges[on].next; }
				while (qh_opp_face(s, an) == ofi) { op = s->edges[op].prev; an = s->edges[an].next; }
				if ((qh_opp_face(s, op) == qh_opp_face(s, an) && qh_face_vert_count(s, qh_opp_face(s, op)) >= 4) ||
				    (qh_opp_face(s, ap) == qh_opp_face(s, on) && qh_face_vert_count(s, qh_opp_face(s, ap)) >= 4))
					will_poly_fin = 1;
			}
			if (will_poly_fin) {
				hedge = s->edges[hedge].next;
				continue;
			}
			int disc[3];
			int nd = qh_merge_adjacent_face(s, hedge, disc, unclaimed);
			for (int i = 0; i < nd; i++) qh_delete_face_points(s, disc[i], fi, unclaimed);
			return 1;
		}
		hedge = s->edges[hedge].next;
	} while (hedge != s->faces[fi].edge);
	if (!convex) s->faces[fi].mark = QH_NON_CONVEX;
	return 0;
}

// -----------------------------------------------------------------------------
// Build cone of triangle faces from eye vertex to horizon, link adjacent sides.

static void qh_face_list_add(QH_State* s, QH_FaceList* list, int fi)
{
	if (list->first == QH_INVALID) { list->first = fi; } else { s->faces[list->last].next = fi; }
	s->faces[fi].next = QH_INVALID;
	list->last = fi;
}

static void qh_add_new_faces(QH_State* s, QH_FaceList* nf, int eye, CK_DYNA int* horizon, int nh)
{
	nf->first = nf->last = QH_INVALID;
	int prev = QH_INVALID, begin = QH_INVALID;
	for (int i = 0; i < nh; i++) {
		int side = qh_add_adjoining_face(s, eye, horizon[i]);
		if (prev != QH_INVALID) qh_set_opposite(s, s->edges[side].next, prev);
		else begin = side;
		qh_face_list_add(s, nf, s->edges[side].face);
		prev = side;
	}
	qh_set_opposite(s, s->edges[begin].next, prev);
}

// Reassign orphaned conflict vertices from the unclaimed list to new faces.
static void qh_resolve_unclaimed(QH_State* s, QH_FaceList* nf, int* unclaimed)
{
	int v = *unclaimed;
	while (v != QH_INVALID) {
		int nxt = s->verts[v].conflict_next;
		float bd = s->epsilon; int bf = QH_INVALID;
		// First try new faces (most likely home for orphaned points).
		for (int fi = nf->first; fi != QH_INVALID; fi = s->faces[fi].next) {
			if (s->faces[fi].mark == QH_VISIBLE) {
				float d = qh_face_dist(s, fi, s->verts[v].pos);
				if (d > bd) { bd = d; bf = fi; }
				if (bd > 1000*s->epsilon) break;
			}
		}
		// Fallback: if no new face claims the point, scan ALL live faces
		// with a relaxed threshold. A point may be numerically inside the
		// current hull but actually belong on it. Use the least-negative
		// distance (closest to being outside) as a last resort.
		if (bf == QH_INVALID) {
			float closest = -1e18f;
			for (int fi = 0; fi < asize(s->faces); fi++) {
				if (s->faces[fi].mark == QH_DELETED) continue;
				float d = qh_face_dist(s, fi, s->verts[v].pos);
				if (d > closest) { closest = d; bf = fi; }
			}
			// Only assign if the point is reasonably close to some face plane.
			// If it's deeply inside (> 100x epsilon below all planes), it's truly interior.
			if (closest < -100.0f * s->epsilon) bf = QH_INVALID;
		}
		if (bf != QH_INVALID) qh_conflict_add(s, bf, v);
		v = nxt;
	}
	*unclaimed = QH_INVALID;
}

// Debug: validate mesh topology of all live faces.
static void qh_validate_mesh(QH_State* s, const char* ctx)
{
	for (int fi = 0; fi < asize(s->faces); fi++) {
		if (s->faces[fi].mark != QH_VISIBLE) continue;
		int e = s->faces[fi].edge, start = e, n = 0;
		do {
			if (n > 500) {
				fprintf(stderr, "qh_validate: infinite loop on face %d at %s\n", fi, ctx);
				exit(-1);
			}
			// Check twin reciprocity.
			int tw = s->edges[e].twin;
			if (s->edges[tw].twin != e) {
				fprintf(stderr, "qh_validate: edge %d twin=%d but twin's twin=%d at %s\n",
					e, tw, s->edges[tw].twin, ctx);
				exit(-1);
			}
			// Check twin's face is alive.
			int twf = s->edges[tw].face;
			if (s->faces[twf].mark != QH_VISIBLE) {
				fprintf(stderr, "qh_validate: edge %d (face %d) twin=%d on deleted face %d at %s\n",
					e, fi, tw, twf, ctx);
				exit(-1);
			}
			// Check edge belongs to this face.
			if (s->edges[e].face != fi) {
				fprintf(stderr, "qh_validate: edge %d face=%d expected %d at %s\n",
					e, s->edges[e].face, fi, ctx);
				exit(-1);
			}
			n++;
			e = s->edges[e].next;
		} while (e != start);
	}

	// Euler check: V - E/2 + F = 2 on the live mesh.
	{
		CK_DYNA int* vremap = NULL;
		afit(vremap, asize(s->verts));
		for (int i = 0; i < asize(s->verts); i++) apush(vremap, 0);
		int nf = 0, ne = 0;
		for (int fi = 0; fi < asize(s->faces); fi++) {
			if (s->faces[fi].mark != QH_VISIBLE) continue;
			nf++;
			int e = s->faces[fi].edge, start = e;
			do {
				ne++;
				vremap[s->edges[e].origin] = 1;
				e = s->edges[e].next;
			} while (e != start);
		}
		int nv = 0;
		for (int i = 0; i < asize(s->verts); i++) nv += vremap[i];
		afree(vremap);
		if (nv - ne/2 + nf != 2) {
			fprintf(stderr, "qh_validate: Euler FAIL V=%d E=%d F=%d (V-E/2+F=%d) at %s\n",
				nv, ne, nf, nv - ne/2 + nf, ctx);
			exit(-1);
		}
	}
}

// -----------------------------------------------------------------------------
// Process one conflict vertex: compute horizon, build cone, merge, reassign.

static void qh_add_point(QH_State* s, int eye, int eye_face)
{
	CK_DYNA int* horizon = NULL;
	int unclaimed = QH_INVALID;

	qh_conflict_remove(s, eye_face, eye);
	qh_calculate_horizon(s, s->verts[eye].pos, QH_INVALID, eye_face, &horizon, &unclaimed);

	QH_FaceList nf;
	qh_add_new_faces(s, &nf, eye, horizon, asize(horizon));

	// Two-pass merge: first larger-face-biased, then mutual non-convexity.
	// Each successful merge reduces the live face count by 1, so the total
	// number of merges across all faces is bounded by the initial face count.
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_VISIBLE)
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_LARGE, &unclaimed));
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_NON_CONVEX) {
			s->faces[fi].mark = QH_VISIBLE;
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_ANY, &unclaimed));
		}


	// Global topology fixup: cascade merges + degenerate connect can create
	// face pairs sharing >1 edge, violating manifold topology. Sweep ALL
	// live faces and forcibly merge multi-adjacent pairs.
	{
		int fixed = 1;
		while (fixed) {
			fixed = 0;
			for (int fi = 0; fi < asize(s->faces) && !fixed; fi++) {
				if (s->faces[fi].mark != QH_VISIBLE) continue;
				int hedge = s->faces[fi].edge;
				do {
					int ofi = qh_opp_face(s, hedge);
					if (ofi == fi || s->faces[ofi].mark != QH_VISIBLE) {
						hedge = s->edges[hedge].next;
						continue;
					}
					int shared = 0;
					int e2 = s->faces[fi].edge;
					do {
						if (qh_opp_face(s, e2) == ofi) shared++;
						e2 = s->edges[e2].next;
					} while (e2 != s->faces[fi].edge);
					if (shared > 1) {
						// Merge smaller into larger so boundary walk sees
						// contiguous shared edges on the absorbed (smaller) face.
						int small = (s->faces[fi].num_verts <= s->faces[ofi].num_verts) ? fi : ofi;
						int large = (small == fi) ? ofi : fi;
						int me = s->faces[large].edge;
						while (qh_opp_face(s, me) != small) me = s->edges[me].next;
						int disc[3];
						int nd = qh_merge_adjacent_face(s, me, disc, &unclaimed);
						for (int k = 0; k < nd; k++)
							qh_delete_face_points(s, disc[k], large, &unclaimed);
						fixed = 1;
						break;
					}
					hedge = s->edges[hedge].next;
				} while (hedge != s->faces[fi].edge);
			}
		}
	}

	qh_resolve_unclaimed(s, &nf, &unclaimed);
	afree(horizon);
}

// -----------------------------------------------------------------------------
// Convert internal state to output Hull with uint16_t-indexed arrays.

static Hull* qh_build_output(QH_State* s, const v3* all_points, int all_count)
{
	CK_DYNA int* live = NULL;
	for (int i = 0; i < asize(s->faces); i++)
		if (s->faces[i].mark == QH_VISIBLE) apush(live, i);

	// Validate face loops before output -- detect broken topology from merges.
	int total_edges = asize(s->edges);
	for (int i = 0; i < asize(live); i++) {
		int fi = live[i];
		int e = s->faces[fi].edge, start = e, steps = 0;
		do {
			steps++;
			if (steps > total_edges) {
				fprintf(stderr, "qh_build_output: face %d (live[%d]) edge loop did not close after %d steps (start edge %d, cur edge %d)\n",
					fi, i, steps, start, e);
				assert(0 && "qh_build_output: broken face edge loop");
			}
			e = s->edges[e].next;
		} while (e != start);
	}

	// Vertex remap.
	CK_DYNA int* vremap = NULL;
	afit(vremap, asize(s->verts));
	for (int i = 0; i < asize(s->verts); i++) apush(vremap, -1);
	CK_DYNA v3* ov = NULL; int vc = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { int vi = s->edges[e].origin;
			if (vremap[vi]<0) { vremap[vi]=vc++; apush(ov, s->verts[vi].pos); }
			e = s->edges[e].next;
		} while (e != start);
	}

	// Edge remap.
	CK_DYNA int* eremap = NULL;
	afit(eremap, asize(s->edges));
	for (int i = 0; i < asize(s->edges); i++) apush(eremap, -1);
	CK_DYNA HalfEdge* oe = NULL; int ec = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { eremap[e]=ec++; HalfEdge he={0}; apush(oe, he); e=s->edges[e].next; } while (e!=start);
	}
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { int o=eremap[e];
			oe[o].next=(uint16_t)eremap[s->edges[e].next];
			oe[o].twin=(uint16_t)eremap[s->edges[e].twin];
			oe[o].origin=(uint16_t)vremap[s->edges[e].origin];
			oe[o].face=(uint16_t)i;
			e=s->edges[e].next;
		} while (e!=start);
	}

	CK_DYNA HullFace* of = NULL;
	CK_DYNA HullPlane* op = NULL;
	for (int i = 0; i < asize(live); i++) {
		HullFace hf = { .edge=(uint16_t)eremap[s->faces[live[i]].edge] };
		apush(of, hf); apush(op, s->faces[live[i]].plane);
	}

	v3 centroid = V3(0,0,0);
	for (int i = 0; i < vc; i++) centroid = add(centroid, ov[i]);
	centroid = scale(centroid, 1.0f / vc);

	Hull* h = CK_ALLOC(sizeof(Hull));
	h->centroid = centroid; h->vert_count = vc; h->edge_count = ec; h->face_count = asize(live);
	h->epsilon = s->epsilon;

	v3* vcp = CK_ALLOC(sizeof(v3)*vc);       memcpy(vcp, ov, sizeof(v3)*vc);       h->verts = vcp;
	HalfEdge* ecp = CK_ALLOC(sizeof(HalfEdge)*ec); memcpy(ecp, oe, sizeof(HalfEdge)*ec); h->edges = ecp;
	HullFace* fcp = CK_ALLOC(sizeof(HullFace)*asize(live)); memcpy(fcp, of, sizeof(HullFace)*asize(live)); h->faces = fcp;
	HullPlane* pcp = CK_ALLOC(sizeof(HullPlane)*asize(live)); memcpy(pcp, op, sizeof(HullPlane)*asize(live)); h->planes = pcp;

	// Fix inward-facing normals: degenerate merges can produce Newell normals
	// that point toward the hull interior. Flip them before widening.
	for (int i = 0; i < h->face_count; i++) {
		v3 fc = V3(0,0,0); int cnt = 0;
		int start = fcp[i].edge, e = start;
		do { fc = add(fc, vcp[ecp[e].origin]); cnt++; e = ecp[e].next; } while (e != start);
		fc = scale(fc, 1.0f / cnt);
		if (dot(pcp[i].normal, sub(fc, centroid)) < 0) {
			pcp[i].normal = scale(pcp[i].normal, -1.0f);
			pcp[i].offset = -pcp[i].offset;
		}
	}

	// Post-build plane widening: widen each plane so ALL original input points
	// lie on or behind the plane. Must use original points (not welded subset)
	// since welded-away points may lie outside the welded hull.
	h->maxoutside = 0;
	for (int i = 0; i < h->face_count; i++) {
		float max_d = pcp[i].offset;
		for (int pi = 0; pi < all_count; pi++) {
			float d = dot(pcp[i].normal, all_points[pi]);
			if (d > max_d) max_d = d;
		}
		float widen = max_d - pcp[i].offset;
		if (widen > h->maxoutside) h->maxoutside = widen;
		pcp[i].offset = max_d;
	}

	afree(live); afree(vremap); afree(ov); afree(eremap); afree(oe); afree(of); afree(op);
	return h;
}

// -----------------------------------------------------------------------------
// Public API.

Hull* quickhull(const v3* points, int count)
{
	qh_dbg_points = points;
	qh_dbg_count = count;
	QH_ASSERT(count >= 4, "quickhull needs at least 4 points");

	QH_State state = {0};
	state.edge_free = QH_INVALID;
	state.face_free = QH_INVALID;

	// Compute epsilon = 3 * (max|x| + max|y| + max|z|) * FLT_EPSILON.
	float mx = 0, my = 0, mz = 0;
	for (int i = 0; i < count; i++) {
		float ax = fabsf(points[i].x), ay = fabsf(points[i].y), az = fabsf(points[i].z);
		if (ax > mx) mx = ax; if (ay > my) my = ay; if (az > mz) mz = az;
	}
	state.epsilon = 3.0f * (mx + my + mz) * FLT_EPSILON;

	// Weld near-duplicate vertices to prevent degenerate triangles that corrupt
	// topology during merge. Weld distance is epsilon (the floating-point
	// precision floor for this input extent).
	float weld_dist2 = state.epsilon * state.epsilon;
	for (int i = 0; i < count; i++) {
		v3 p = points[i];
		int dup = 0;
		for (int j = 0; j < asize(state.verts); j++) {
			if (len2(sub(p, state.verts[j].pos)) <= weld_dist2) { dup = 1; break; }
		}
		if (!dup) {
			QH_Vertex v = { .pos=p, .conflict_next=QH_INVALID, .conflict_prev=QH_INVALID };
			apush(state.verts, v);
		}
	}
	int welded_count = asize(state.verts);

	if (welded_count < 4) {
		afree(state.verts); afree(state.edges); afree(state.faces);
		return NULL;
	}

	if (!qh_build_simplex(&state, welded_count)) {
		afree(state.verts); afree(state.edges); afree(state.faces);
		return NULL;
	}

	for (int iter = 0; iter < welded_count * 4; iter++) {
		int fi; int vi = qh_next_conflict(&state, &fi);
		if (vi == QH_INVALID) break;
		qh_add_point(&state, vi, fi);
		QH_DEBUG({
			char buf[64];
			snprintf(buf, sizeof(buf), "after adding vertex %d (iter %d)", vi, iter);
			qh_validate_mesh(&state, buf);
		});
	}

	Hull* result = qh_build_output(&state, points, count);
	afree(state.verts); afree(state.edges); afree(state.faces);
	return result;
}

// ---- src/bvh.c ----

// See LICENSE for licensing info.
// bvh.c -- binary BVH broadphase (Bepu-inspired, 64-byte cache-line nodes).

// -----------------------------------------------------------------------------
// AABB type and helpers.

typedef struct AABB { v3 min, max; } AABB;

static AABB aabb_empty() { return (AABB){ V3(1e18f, 1e18f, 1e18f), V3(-1e18f, -1e18f, -1e18f) }; }

static AABB aabb_merge(AABB a, AABB b) { return (AABB){ v3_min(a.min, b.min), v3_max(a.max, b.max) }; }

static AABB aabb_from_point(v3 p) { return (AABB){ p, p }; }

static AABB aabb_expand(AABB a, float margin) { v3 m = V3(margin, margin, margin); return (AABB){ sub(a.min, m), add(a.max, m) }; }

static float aabb_surface_area(AABB a) { v3 d = sub(a.max, a.min); return d.x*d.y + d.y*d.z + d.z*d.x; }

static int aabb_overlaps(AABB a, AABB b) { return a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y && a.min.z <= b.max.z && a.max.z >= b.min.z; }

// Compute world-space AABB for a single shape on a body.
static AABB shape_aabb(BodyHot* h, ShapeInternal* s)
{
	v3 world_pos = add(h->position, rotate(h->rotation, s->local_pos));
	switch (s->type) {
	case SHAPE_SPHERE: {
		v3 r = V3(s->sphere.radius, s->sphere.radius, s->sphere.radius);
		return (AABB){ sub(world_pos, r), add(world_pos, r) };
	}
	case SHAPE_CAPSULE: {
		v3 up = rotate(h->rotation, V3(0, s->capsule.half_height, 0));
		v3 p = sub(world_pos, up), q = add(world_pos, up);
		v3 r = V3(s->capsule.radius, s->capsule.radius, s->capsule.radius);
		return (AABB){ sub(v3_min(p, q), r), add(v3_max(p, q), r) };
	}
	case SHAPE_BOX: {
		// OBB -> AABB: project rotated half-extents onto each world axis.
		v3 e = s->box.half_extents;
		v3 ax = rotate(h->rotation, V3(e.x, 0, 0));
		v3 ay = rotate(h->rotation, V3(0, e.y, 0));
		v3 az = rotate(h->rotation, V3(0, 0, e.z));
		v3 half = V3(fabsf(ax.x) + fabsf(ay.x) + fabsf(az.x), fabsf(ax.y) + fabsf(ay.y) + fabsf(az.y), fabsf(ax.z) + fabsf(ay.z) + fabsf(az.z));
		return (AABB){ sub(world_pos, half), add(world_pos, half) };
	}
	case SHAPE_HULL: {
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		AABB box = aabb_from_point(add(world_pos, rotate(h->rotation, V3(hull->verts[0].x*sc.x, hull->verts[0].y*sc.y, hull->verts[0].z*sc.z))));
		for (int i = 1; i < hull->vert_count; i++) {
			v3 v = add(world_pos, rotate(h->rotation, V3(hull->verts[i].x*sc.x, hull->verts[i].y*sc.y, hull->verts[i].z*sc.z)));
			box.min = v3_min(box.min, v);
			box.max = v3_max(box.max, v);
		}
		return box;
	}
	}
	return aabb_empty();
}

// Compute world-space AABB for an entire body (union of all shapes).
static AABB body_aabb(BodyHot* h, BodyCold* c)
{
	AABB box = shape_aabb(h, &c->shapes[0]);
	for (int i = 1; i < asize(c->shapes); i++)
		box = aabb_merge(box, shape_aabb(h, &c->shapes[i]));
	return box;
}

// -----------------------------------------------------------------------------
// BVH node layout: 64 bytes = one cache line.
// Two children inline. Metadata (parent pointers) stored separately.

typedef struct BVHChild
{
	v3 min;          // 12
	int32_t index;   // 4 -- negative: encoded leaf (~idx), non-negative: child node
	v3 max;          // 12
	int32_t leaf_count; // 4 -- 0=empty, 1=leaf, >1=internal subtree
} BVHChild;          // 32 bytes

typedef struct BVHNode
{
	BVHChild a; // 32
	BVHChild b; // 32
} BVHNode;      // 64 bytes

typedef struct BVHMeta
{
	int parent;     // parent node index, -1 for root
	int child_slot; // 0=A, 1=B
	int dirty;      // 1 if subtree needs refit, 0 if clean since last refit
} BVHMeta;

typedef struct BVHLeaf
{
	int body_idx;
	int node_idx;
	int child_slot; // 0=A, 1=B
	v3 fat_min, fat_max; // expanded AABB for motion threshold
} BVHLeaf;

typedef struct BVHTree
{
	CK_DYNA BVHNode* nodes;
	CK_DYNA BVHMeta* meta;
	CK_DYNA BVHLeaf* leaves;
	CK_DYNA int* node_free;
	int root; // -1 = empty
	int refine_cursor; // leaf-space position for incremental refinement
} BVHTree;

typedef struct BroadPair { int a, b; } BroadPair;

#define BVH_AABB_MARGIN 0.05f
#define BVH_FAT_MARGIN  0.2f

// -----------------------------------------------------------------------------
// Tree lifecycle.

static void bvh_init(BVHTree* t) { memset(t, 0, sizeof(*t)); t->root = -1; }

static void bvh_free(BVHTree* t) { afree(t->nodes); afree(t->meta); afree(t->leaves); afree(t->node_free); }

// Alloc/free with freelist.
static int bvh_alloc_node(BVHTree* t)
{
	int idx;
	if (asize(t->node_free) > 0) { idx = apop(t->node_free); }
	else { idx = asize(t->nodes); BVHNode z = {0}; apush(t->nodes, z); BVHMeta zm = {-1, 0, 1}; apush(t->meta, zm); }
	t->nodes[idx] = (BVHNode){0};
	t->meta[idx] = (BVHMeta){-1, 0, 1};
	return idx;
}

static void bvh_free_node(BVHTree* t, int idx) { apush(t->node_free, idx); }

static int bvh_alloc_leaf(BVHTree* t)
{
	int idx = asize(t->leaves);
	BVHLeaf z = {0};
	apush(t->leaves, z);
	return idx;
}

// -----------------------------------------------------------------------------
// Child helpers.

static BVHChild* bvh_child(BVHNode* n, int slot) { return slot == 0 ? &n->a : &n->b; }

static AABB bvh_child_aabb(BVHChild* c) { return (AABB){ c->min, c->max }; }

static void bvh_child_set_aabb(BVHChild* c, AABB box) { c->min = box.min; c->max = box.max; }

static void bvh_child_set_leaf(BVHChild* c, AABB box, int leaf_idx) { bvh_child_set_aabb(c, box); c->index = ~leaf_idx; c->leaf_count = 1; }

// Place a leaf with fat AABB. The node child stores the fat AABB for conservative overlap.
static void bvh_place_leaf(BVHTree* t, int ni, int slot, AABB tight, int leaf_idx)
{
	AABB fat = aabb_expand(tight, BVH_FAT_MARGIN);
	bvh_child_set_leaf(bvh_child(&t->nodes[ni], slot), fat, leaf_idx);
	t->leaves[leaf_idx].node_idx = ni;
	t->leaves[leaf_idx].child_slot = slot;
	t->leaves[leaf_idx].fat_min = fat.min;
	t->leaves[leaf_idx].fat_max = fat.max;
}

static void bvh_child_set_node(BVHChild* c, AABB box, int node_idx, int lcount) { bvh_child_set_aabb(c, box); c->index = node_idx; c->leaf_count = lcount; }

static void bvh_child_set_empty(BVHChild* c) { *c = (BVHChild){0}; }

static int bvh_child_is_leaf(BVHChild* c) { return c->leaf_count == 1; }
static int bvh_child_is_internal(BVHChild* c) { return c->leaf_count > 1; }
static int bvh_child_is_empty(BVHChild* c) { return c->leaf_count == 0; }
static int bvh_child_leaf_idx(BVHChild* c) { return ~c->index; }

static AABB bvh_node_aabb(BVHNode* n)
{
	if (n->a.leaf_count == 0) return bvh_child_aabb(&n->b);
	if (n->b.leaf_count == 0) return bvh_child_aabb(&n->a);
	return aabb_merge(bvh_child_aabb(&n->a), bvh_child_aabb(&n->b));
}

// Count leaves reachable from a node (for validation).
static int bvh_count_leaves(BVHTree* t, int ni)
{
	BVHNode* n = &t->nodes[ni];
	int count = 0;
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		if (bvh_child_is_leaf(c)) count++;
		else count += bvh_count_leaves(t, c->index);
	}
	return count;
}

// -----------------------------------------------------------------------------
// Tree rotation: swap a grandchild with its uncle to reduce SAH cost.
// Evaluates 4 candidates at a node and applies the best improvement.

static void bvh_try_rotate(BVHTree* t, int ni)
{
	BVHNode* n = &t->nodes[ni];
	if (bvh_child_is_empty(&n->a) || bvh_child_is_empty(&n->b)) return;

	float base_cost = aabb_surface_area(bvh_child_aabb(&n->a)) + aabb_surface_area(bvh_child_aabb(&n->b));
	float best_cost = base_cost;
	int best_rot = -1; // 0-3: which rotation to apply

	// For each internal child, try swapping each of its grandchildren with the other child.
	for (int side = 0; side < 2; side++) {
		BVHChild* parent_child = bvh_child(n, side);
		BVHChild* uncle = bvh_child(n, 1 - side);
		if (!bvh_child_is_internal(parent_child)) continue;

		BVHNode* pn = &t->nodes[parent_child->index];
		for (int gc = 0; gc < 2; gc++) {
			BVHChild* grandchild = bvh_child(pn, gc);
			BVHChild* sibling = bvh_child(pn, 1 - gc);
			if (bvh_child_is_empty(grandchild) || bvh_child_is_empty(sibling)) continue;

			// After rotation: parent_child would contain (sibling, uncle_old), uncle slot gets grandchild.
			AABB new_parent_aabb = aabb_merge(bvh_child_aabb(sibling), bvh_child_aabb(uncle));
			float new_parent_sa = aabb_surface_area(new_parent_aabb);
			float new_uncle_sa = aabb_surface_area(bvh_child_aabb(grandchild));
			float cost = new_parent_sa + new_uncle_sa;

			if (cost < best_cost - 1e-6f) {
				best_cost = cost;
				best_rot = side * 2 + gc;
			}
		}
	}

	if (best_rot < 0) return;

	// Apply the rotation.
	int side = best_rot / 2;
	int gc = best_rot % 2;
	BVHChild* parent_child = bvh_child(n, side);
	BVHChild* uncle_slot = bvh_child(n, 1 - side);
	int pn_idx = parent_child->index;
	BVHNode* pn = &t->nodes[pn_idx];

	// Save grandchild and uncle.
	BVHChild gc_saved = *bvh_child(pn, gc);
	BVHChild uncle_saved = *uncle_slot;

	// Grandchild moves to uncle's slot in N.
	*uncle_slot = gc_saved;
	// Uncle moves into grandchild's slot in parent_child's node.
	*bvh_child(pn, gc) = uncle_saved;

	// Update parent_child's AABB and leaf_count in N.
	bvh_child_set_aabb(parent_child, bvh_node_aabb(pn));
	parent_child->leaf_count = pn->a.leaf_count + pn->b.leaf_count;

	// Fix back-pointers for the swapped children.
	// The grandchild (now in uncle_slot at node ni):
	if (bvh_child_is_leaf(uncle_slot)) {
		t->leaves[bvh_child_leaf_idx(uncle_slot)].node_idx = ni;
		t->leaves[bvh_child_leaf_idx(uncle_slot)].child_slot = 1 - side;
	} else if (bvh_child_is_internal(uncle_slot)) {
		t->meta[uncle_slot->index].parent = ni;
		t->meta[uncle_slot->index].child_slot = 1 - side;
	}
	// The uncle (now in pn at slot gc):
	BVHChild* new_gc = bvh_child(pn, gc);
	if (bvh_child_is_leaf(new_gc)) {
		t->leaves[bvh_child_leaf_idx(new_gc)].node_idx = pn_idx;
		t->leaves[bvh_child_leaf_idx(new_gc)].child_slot = gc;
	} else if (bvh_child_is_internal(new_gc)) {
		t->meta[new_gc->index].parent = pn_idx;
		t->meta[new_gc->index].child_slot = gc;
	}
}

// -----------------------------------------------------------------------------
// Insertion: SAH-guided top-down.

static int bvh_insert(BVHTree* t, int body_idx, AABB bounds)
{
	int li = bvh_alloc_leaf(t);
	t->leaves[li] = (BVHLeaf){ .body_idx = body_idx };

	// Empty tree: create root with leaf in child A.
	if (t->root == -1) {
		int ni = bvh_alloc_node(t);
		t->meta[ni] = (BVHMeta){-1, 0, 1};
		bvh_place_leaf(t, ni, 0, bounds, li);
		t->root = ni;
		return li;
	}

	// Single leaf in tree: root has child A occupied and B empty (or vice versa).
	BVHNode* root = &t->nodes[t->root];
	if (bvh_child_is_empty(&root->b) && bvh_child_is_leaf(&root->a)) {
		bvh_place_leaf(t, t->root, 1, bounds, li);
		return li;
	}
	if (bvh_child_is_empty(&root->a) && bvh_child_is_leaf(&root->b)) {
		bvh_place_leaf(t, t->root, 0, bounds, li);
		return li;
	}

	// Walk tree from root to find best insertion point.
	int cur = t->root;
	for (;;) {
		BVHNode* node = &t->nodes[cur];
		// Full SAH cost: SA(merged)*(leafCount+1) - SA(original)*leafCount.
		AABB ma = aabb_merge(bvh_child_aabb(&node->a), bounds);
		AABB mb = aabb_merge(bvh_child_aabb(&node->b), bounds);
		float cost_a = aabb_surface_area(ma) * (node->a.leaf_count + 1) - aabb_surface_area(bvh_child_aabb(&node->a)) * node->a.leaf_count;
		float cost_b = aabb_surface_area(mb) * (node->b.leaf_count + 1) - aabb_surface_area(bvh_child_aabb(&node->b)) * node->b.leaf_count;

		int slot = cost_a == cost_b ? (node->a.leaf_count <= node->b.leaf_count ? 0 : 1) : (cost_a < cost_b ? 0 : 1);
		BVHChild* child = bvh_child(node, slot);

		if (bvh_child_is_leaf(child)) {
			// Split: create new node containing existing leaf + new leaf.
			int existing_li = bvh_child_leaf_idx(child);
			AABB existing_aabb = bvh_child_aabb(child);
			AABB merged = slot == 0 ? ma : mb;

			int new_ni = bvh_alloc_node(t);
			// NOTE: alloc may realloc t->nodes, so re-fetch pointers.

			// Existing leaf keeps its current (fat) AABB.
			bvh_child_set_leaf(&t->nodes[new_ni].a, existing_aabb, existing_li);
			t->leaves[existing_li].node_idx = new_ni;
			t->leaves[existing_li].child_slot = 0;

			// New leaf gets fat AABB.
			bvh_place_leaf(t, new_ni, 1, bounds, li);

			// Replace child in parent with internal node pointer (re-fetch after alloc).
			bvh_child_set_node(bvh_child(&t->nodes[cur], slot), merged, new_ni, 2);
			t->meta[new_ni] = (BVHMeta){ cur, slot, 1 };

			// Walk back up refitting bounds, applying rotations, and marking dirty.
			int ri = cur;
			while (ri != -1) {
				t->meta[ri].dirty = 1;
				BVHNode* rn = &t->nodes[ri];
				if (bvh_child_is_internal(&rn->a)) { BVHNode* cn = &t->nodes[rn->a.index]; rn->a.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->a, bvh_node_aabb(cn)); }
				if (bvh_child_is_internal(&rn->b)) { BVHNode* cn = &t->nodes[rn->b.index]; rn->b.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->b, bvh_node_aabb(cn)); }
				bvh_try_rotate(t, ri);
				ri = t->meta[ri].parent;
			}
			return li;
		}

		if (bvh_child_is_internal(child)) {
			// Descend into the child subtree.
			cur = child->index;
			continue;
		}

		// Empty slot: place leaf directly.
		bvh_place_leaf(t, cur, slot, bounds, li);
		return li;
	}
}

// -----------------------------------------------------------------------------
// Removal.

// Returns body_idx of the leaf that was swapped into the removed slot, or -1.
static int bvh_remove(BVHTree* t, int leaf_idx)
{
	BVHLeaf* leaf = &t->leaves[leaf_idx];
	int ni = leaf->node_idx;
	int slot = leaf->child_slot;

	// Clear the child slot.
	bvh_child_set_empty(bvh_child(&t->nodes[ni], slot));

	// Get the sibling.
	int sib_slot = 1 - slot;
	BVHChild sib = *bvh_child(&t->nodes[ni], sib_slot);

	int parent = t->meta[ni].parent;
	if (parent == -1) {
		// ni is root. If sibling is empty, tree becomes empty. Otherwise keep root with one child.
		if (bvh_child_is_empty(&sib)) { bvh_free_node(t, ni); t->root = -1; }
		// else: root has one child remaining, that's fine.
	} else {
		// Replace ni in parent with sibling.
		int parent_slot = t->meta[ni].child_slot;
		BVHChild* pc = bvh_child(&t->nodes[parent], parent_slot);
		*pc = sib;

		// Update sibling's back-pointers.
		if (bvh_child_is_leaf(pc)) {
			int sib_li = bvh_child_leaf_idx(pc);
			t->leaves[sib_li].node_idx = parent;
			t->leaves[sib_li].child_slot = parent_slot;
		} else if (bvh_child_is_internal(pc)) {
			t->meta[pc->index].parent = parent;
			t->meta[pc->index].child_slot = parent_slot;
		}

		bvh_free_node(t, ni);

		// Refit upward and mark dirty.
		int ri = parent;
		while (ri != -1) {
			t->meta[ri].dirty = 1;
			BVHNode* rn = &t->nodes[ri];
			if (bvh_child_is_internal(&rn->a)) { BVHNode* cn = &t->nodes[rn->a.index]; rn->a.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->a, bvh_node_aabb(cn)); }
			if (bvh_child_is_internal(&rn->b)) { BVHNode* cn = &t->nodes[rn->b.index]; rn->b.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->b, bvh_node_aabb(cn)); }
			ri = t->meta[ri].parent;
		}
	}

	// Swap-back: move last leaf into freed slot, keep array compact.
	int last = asize(t->leaves) - 1;
	if (leaf_idx != last) {
		t->leaves[leaf_idx] = t->leaves[last];
		BVHChild* lc = bvh_child(&t->nodes[t->leaves[leaf_idx].node_idx], t->leaves[leaf_idx].child_slot);
		lc->index = ~leaf_idx;
		apop(t->leaves);
		return t->leaves[leaf_idx].body_idx;
	}
	apop(t->leaves);
	return -1;
}

// -----------------------------------------------------------------------------
// Binned SAH builder: top-down build from a set of leaf indices.
// lut maps leaf_index -> AABB (indexed by leaf index, not array position).

#define BVH_SAH_BINS 8

// Build a subtree from lis[0..count-1]. Each element is a leaf index.
// lut[leaf_idx] gives the AABB for that leaf.
// If count==1, returns ~leaf_idx (encoded as a "leaf result"). Caller checks sign.
// If count>=2, returns a node index.
static int bvh_binned_build(BVHTree* t, int* lis, AABB* lut, int count)
{
	assert(count >= 2);

	if (count == 2) {
		int ni = bvh_alloc_node(t);
		bvh_child_set_leaf(&t->nodes[ni].a, lut[lis[0]], lis[0]);
		bvh_child_set_leaf(&t->nodes[ni].b, lut[lis[1]], lis[1]);
		t->leaves[lis[0]].node_idx = ni; t->leaves[lis[0]].child_slot = 0;
		t->leaves[lis[1]].node_idx = ni; t->leaves[lis[1]].child_slot = 1;
		return ni;
	}

	// Compute centroid bounds.
	v3 cmin = V3(1e18f, 1e18f, 1e18f), cmax = V3(-1e18f, -1e18f, -1e18f);
	for (int i = 0; i < count; i++) {
		v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
		cmin = v3_min(cmin, c); cmax = v3_max(cmax, c);
	}

	// Find best split across all 3 axes.
	float best_cost = 1e18f;
	int best_axis = 0, best_split = 0;

	for (int axis = 0; axis < 3; axis++) {
		float amin = (&cmin.x)[axis], amax = (&cmax.x)[axis];
		float extent = amax - amin;
		if (extent < 1e-7f) continue;

		int bin_count[BVH_SAH_BINS] = {0};
		AABB bin_aabb[BVH_SAH_BINS];
		for (int b = 0; b < BVH_SAH_BINS; b++) bin_aabb[b] = aabb_empty();

		float inv_ext = (float)BVH_SAH_BINS / extent;
		for (int i = 0; i < count; i++) {
			v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
			int b = (int)(((&c.x)[axis] - amin) * inv_ext);
			if (b >= BVH_SAH_BINS) b = BVH_SAH_BINS - 1;
			bin_count[b]++; bin_aabb[b] = aabb_merge(bin_aabb[b], lut[lis[i]]);
		}

		// Prefix scan left and right.
		AABB la[BVH_SAH_BINS - 1], ra[BVH_SAH_BINS - 1];
		int lc[BVH_SAH_BINS - 1], rc[BVH_SAH_BINS - 1];
		AABB run = aabb_empty(); int rn = 0;
		for (int i = 0; i < BVH_SAH_BINS - 1; i++) { run = aabb_merge(run, bin_aabb[i]); rn += bin_count[i]; la[i] = run; lc[i] = rn; }
		run = aabb_empty(); rn = 0;
		for (int i = BVH_SAH_BINS - 1; i > 0; i--) { run = aabb_merge(run, bin_aabb[i]); rn += bin_count[i]; ra[i-1] = run; rc[i-1] = rn; }

		for (int i = 0; i < BVH_SAH_BINS - 1; i++) {
			if (lc[i] == 0 || rc[i] == 0) continue;
			float cost = aabb_surface_area(la[i]) * lc[i] + aabb_surface_area(ra[i]) * rc[i];
			if (cost < best_cost) { best_cost = cost; best_axis = axis; best_split = i; }
		}
	}

	// Partition by best split.
	float amin = (&cmin.x)[best_axis], amax = (&cmax.x)[best_axis];
	float extent = amax - amin;
	float inv_ext = extent > 1e-7f ? (float)BVH_SAH_BINS / extent : 0.0f;

	int mid = 0;
	for (int i = 0; i < count; i++) {
		v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
		int b = inv_ext > 0 ? (int)(((&c.x)[best_axis] - amin) * inv_ext) : 0;
		if (b >= BVH_SAH_BINS) b = BVH_SAH_BINS - 1;
		if (b <= best_split) { int tmp = lis[mid]; lis[mid] = lis[i]; lis[i] = tmp; mid++; }
	}
	if (mid == 0 || mid == count) mid = count / 2;

	// Build children. Each half has count >= 1.
	int ni = bvh_alloc_node(t);

	// Left child.
	if (mid == 1) {
		bvh_child_set_leaf(&t->nodes[ni].a, lut[lis[0]], lis[0]);
		t->leaves[lis[0]].node_idx = ni; t->leaves[lis[0]].child_slot = 0;
	} else {
		int left = bvh_binned_build(t, lis, lut, mid);
		BVHNode* ln = &t->nodes[left];
		bvh_child_set_node(&t->nodes[ni].a, bvh_node_aabb(ln), left, ln->a.leaf_count + ln->b.leaf_count);
		t->meta[left].parent = ni; t->meta[left].child_slot = 0;
	}

	// Right child.
	int rcount = count - mid;
	if (rcount == 1) {
		bvh_child_set_leaf(&t->nodes[ni].b, lut[lis[mid]], lis[mid]);
		t->leaves[lis[mid]].node_idx = ni; t->leaves[lis[mid]].child_slot = 1;
	} else {
		int right = bvh_binned_build(t, lis + mid, lut, rcount);
		BVHNode* rnode = &t->nodes[right];
		bvh_child_set_node(&t->nodes[ni].b, bvh_node_aabb(rnode), right, rnode->a.leaf_count + rnode->b.leaf_count);
		t->meta[right].parent = ni; t->meta[right].child_slot = 1;
	}

	return ni;
}

// Collect all leaf indices reachable from node ni into lis array. Frees internal nodes.
static void bvh_collect_and_free(BVHTree* t, int ni, CK_DYNA int** lis)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		if (bvh_child_is_leaf(c)) apush(*lis, bvh_child_leaf_idx(c));
		else bvh_collect_and_free(t, c->index, lis);
	}
	bvh_free_node(t, ni);
}

// Rebuild the subtree at node ni using binned SAH.
// parent_ni and parent_slot describe where ni is attached (-1 for root).
static void bvh_refine_subtree(BVHTree* t, int ni, AABB* lut)
{
	int parent_ni = t->meta[ni].parent;
	int parent_slot = t->meta[ni].child_slot;

	// Collect leaves and free old internal nodes.
	CK_DYNA int* lis = NULL;
	bvh_collect_and_free(t, ni, &lis);
	int count = asize(lis);

	if (count < 2) {
		// Edge case: 0 or 1 leaves. Just put the leaf back.
		if (count == 1 && parent_ni >= 0) {
			bvh_child_set_leaf(bvh_child(&t->nodes[parent_ni], parent_slot), lut[lis[0]], lis[0]);
			t->leaves[lis[0]].node_idx = parent_ni;
			t->leaves[lis[0]].child_slot = parent_slot;
		}
		afree(lis);
		return;
	}

	// Rebuild.
	int new_root = bvh_binned_build(t, lis, lut, count);

	// Splice into parent.
	if (parent_ni == -1) {
		t->root = new_root;
		t->meta[new_root].parent = -1;
	} else {
		BVHNode* cn = &t->nodes[new_root];
		AABB box = bvh_node_aabb(cn);
		bvh_child_set_node(bvh_child(&t->nodes[parent_ni], parent_slot), box, new_root, cn->a.leaf_count + cn->b.leaf_count);
		t->meta[new_root].parent = parent_ni;
		t->meta[new_root].child_slot = parent_slot;
	}

	afree(lis);
}

// Walk the tree to find a subtree near target_size that overlaps the leaf-space cursor.
// left_leaves is the number of leaves to the left of node ni in DFS order.
static int bvh_find_refine_target(BVHTree* t, int ni, int left_leaves, int cursor, int target_size)
{
	BVHNode* n = &t->nodes[ni];
	int total = n->a.leaf_count + n->b.leaf_count;
	if (total <= target_size) return ni; // this subtree fits -- refine it

	// Descend into the child that contains the cursor position.
	int mid = left_leaves + n->a.leaf_count;
	if (cursor < mid && bvh_child_is_internal(&n->a))
		return bvh_find_refine_target(t, n->a.index, left_leaves, cursor, target_size);
	if (bvh_child_is_internal(&n->b))
		return bvh_find_refine_target(t, n->b.index, mid, cursor, target_size);
	if (bvh_child_is_internal(&n->a))
		return bvh_find_refine_target(t, n->a.index, left_leaves, cursor, target_size);
	return ni; // both children are leaves; refine this node
}

// Per-frame incremental refinement: rebuild multiple subtrees using binned SAH.
// Covers ~5% of the tree per frame via a rotating cursor through leaf-space.
static void bvh_incremental_refine(BVHTree* t, AABB* lut)
{
	if (t->root == -1) return;
	BVHNode* root = &t->nodes[t->root];
	int total_leaves = root->a.leaf_count + root->b.leaf_count;
	if (total_leaves < 4) return;

	int target_size = (int)sqrtf((float)total_leaves);
	if (target_size < 2) target_size = 2;

	// Refine enough subtrees to cover ~5% of leaves per frame.
	int budget = total_leaves / 20;
	if (budget < target_size) budget = target_size;
	int refined = 0;

	while (refined < budget) {
		if (t->refine_cursor >= total_leaves) t->refine_cursor = 0;
		int ni = bvh_find_refine_target(t, t->root, 0, t->refine_cursor, target_size);
		if (ni < 0) break;
		int lc = t->nodes[ni].a.leaf_count + t->nodes[ni].b.leaf_count;
		bvh_refine_subtree(t, ni, lut);
		t->refine_cursor += lc;
		refined += lc;
	}
}

// -----------------------------------------------------------------------------
// Standalone DFS cache reorder (no refit). Used when no WorldInternal is available.

static void bvh_dfs_collect(BVHTree* t, int ni, int* remap, BVHNode* dst_nodes, BVHMeta* dst_meta, int* cursor)
{
	int new_idx = (*cursor)++;
	remap[ni] = new_idx;
	dst_nodes[new_idx] = t->nodes[ni];
	dst_meta[new_idx] = t->meta[ni];
	BVHNode* n = &t->nodes[ni];
	if (bvh_child_is_internal(&n->a)) bvh_dfs_collect(t, n->a.index, remap, dst_nodes, dst_meta, cursor);
	if (bvh_child_is_internal(&n->b)) bvh_dfs_collect(t, n->b.index, remap, dst_nodes, dst_meta, cursor);
}

static void bvh_cache_reorder(BVHTree* t)
{
	if (t->root == -1) return;
	int cap = asize(t->nodes);
	int* remap = CK_ALLOC(sizeof(int) * cap);
	for (int i = 0; i < cap; i++) remap[i] = -1;
	BVHNode* new_nodes = CK_ALLOC(sizeof(BVHNode) * cap);
	BVHMeta* new_meta = CK_ALLOC(sizeof(BVHMeta) * cap);
	int cursor = 0;
	bvh_dfs_collect(t, t->root, remap, new_nodes, new_meta, &cursor);
	int live_count = cursor;
	for (int i = 0; i < live_count; i++) {
		BVHNode* n = &new_nodes[i];
		if (bvh_child_is_internal(&n->a)) n->a.index = remap[n->a.index];
		if (bvh_child_is_internal(&n->b)) n->b.index = remap[n->b.index];
	}
	for (int i = 0; i < live_count; i++) {
		if (new_meta[i].parent >= 0) new_meta[i].parent = remap[new_meta[i].parent];
	}
	for (int i = 0; i < asize(t->leaves); i++) {
		if (remap[t->leaves[i].node_idx] >= 0) t->leaves[i].node_idx = remap[t->leaves[i].node_idx];
	}
	aclear(t->nodes); aclear(t->meta);
	for (int i = 0; i < live_count; i++) { apush(t->nodes, new_nodes[i]); apush(t->meta, new_meta[i]); }
	aclear(t->node_free);
	t->root = remap[t->root];
	CK_FREE(remap); CK_FREE(new_nodes); CK_FREE(new_meta);
}

// -----------------------------------------------------------------------------
// Fused refit + cache reorder: single DFS pass that updates leaf AABBs,
// merges bounds bottom-up, and writes nodes in DFS order for cache locality.
// Reads from old node arrays (src_*), writes to new arrays (dst_*).
// target_idx is the DFS-sequential index for the current node.

// Check if a body is sleeping (in a sleeping island).
static int bvh_body_sleeping(WorldInternal* w, int bi)
{
	int isl = w->body_cold[bi].island_id;
	return isl >= 0 && (w->island_gen[isl] & 1) && !w->islands[isl].awake;
}

typedef struct BVHRefit
{
	BVHNode* src_nodes;
	BVHMeta* src_meta;
	BVHNode* dst_nodes;
	BVHMeta* dst_meta;
	BVHTree* tree;
	WorldInternal* world;
} BVHRefit;

// Refit+reorder a single child. Returns updated child for parent to store.
// target_idx is where the child node should go in the new array.
static BVHChild bvh_fused_recurse(BVHRefit* r, int src_ni, int target_idx, int parent_target, int parent_slot, int* changed)
{
	BVHNode* src = &r->src_nodes[src_ni];
	BVHMeta* smeta = &r->src_meta[src_ni];

	// Skip subtrees where all bodies were sleeping.
	if (!smeta->dirty) {
		// Copy entire subtree as-is in DFS order (no refit needed).
		// But we still need to remap indices for DFS layout.
		// Fall through to normal processing -- the leaf fat AABB checks will all pass for sleeping bodies.
	}

	BVHNode* dst = &r->dst_nodes[target_idx];
	BVHMeta* dmeta = &r->dst_meta[target_idx];
	dmeta->parent = parent_target;
	dmeta->child_slot = parent_slot;
	dmeta->dirty = smeta->dirty;

	int all_sleeping = 1;
	int any_changed = 0;
	// DFS layout: child A subtree at target+1, child B subtree at target+A.leaf_count.
	// An internal child with N leaves uses N-1 internal node slots.
	int a_nodes = bvh_child_is_internal(&src->a) ? src->a.leaf_count - 1 : 0;
	int target_a = target_idx + 1;
	int target_b = target_idx + 1 + a_nodes;

	// Process child A.
	BVHChild sa = src->a;
	if (bvh_child_is_leaf(&sa)) {
		int li = bvh_child_leaf_idx(&sa);
		int bi = r->tree->leaves[li].body_idx;
		dst->a = sa;
		if (!bvh_body_sleeping(r->world, bi)) {
			all_sleeping = 0;
			AABB tight = aabb_expand(body_aabb(&r->world->body_hot[bi], &r->world->body_cold[bi]), BVH_AABB_MARGIN);
			v3 fmin = r->tree->leaves[li].fat_min, fmax = r->tree->leaves[li].fat_max;
			if (!(tight.min.x >= fmin.x && tight.min.y >= fmin.y && tight.min.z >= fmin.z && tight.max.x <= fmax.x && tight.max.y <= fmax.y && tight.max.z <= fmax.z)) {
				AABB fat = aabb_expand(tight, BVH_FAT_MARGIN);
				r->tree->leaves[li].fat_min = fat.min;
				r->tree->leaves[li].fat_max = fat.max;
				bvh_child_set_aabb(&dst->a, fat);
				any_changed = 1;
			}
		}
		r->tree->leaves[li].node_idx = target_idx;
		r->tree->leaves[li].child_slot = 0;
	} else if (bvh_child_is_internal(&sa)) {
		dst->a = bvh_fused_recurse(r, sa.index, target_a, target_idx, 0, &any_changed);
		dst->a.index = target_a;
		dst->a.leaf_count = sa.leaf_count;
		if (r->dst_meta[target_a].dirty) all_sleeping = 0;
	} else {
		dst->a = sa;
	}

	// Process child B.
	BVHChild sb = src->b;
	if (bvh_child_is_leaf(&sb)) {
		int li = bvh_child_leaf_idx(&sb);
		int bi = r->tree->leaves[li].body_idx;
		dst->b = sb;
		if (!bvh_body_sleeping(r->world, bi)) {
			all_sleeping = 0;
			AABB tight = aabb_expand(body_aabb(&r->world->body_hot[bi], &r->world->body_cold[bi]), BVH_AABB_MARGIN);
			v3 fmin = r->tree->leaves[li].fat_min, fmax = r->tree->leaves[li].fat_max;
			if (!(tight.min.x >= fmin.x && tight.min.y >= fmin.y && tight.min.z >= fmin.z && tight.max.x <= fmax.x && tight.max.y <= fmax.y && tight.max.z <= fmax.z)) {
				AABB fat = aabb_expand(tight, BVH_FAT_MARGIN);
				r->tree->leaves[li].fat_min = fat.min;
				r->tree->leaves[li].fat_max = fat.max;
				bvh_child_set_aabb(&dst->b, fat);
				any_changed = 1;
			}
		}
		r->tree->leaves[li].node_idx = target_idx;
		r->tree->leaves[li].child_slot = 1;
	} else if (bvh_child_is_internal(&sb)) {
		dst->b = bvh_fused_recurse(r, sb.index, target_b, target_idx, 1, &any_changed);
		dst->b.index = target_b;
		dst->b.leaf_count = sb.leaf_count;
		if (r->dst_meta[target_b].dirty) all_sleeping = 0;
	} else {
		dst->b = sb;
	}

	dmeta->dirty = !all_sleeping;
	*changed |= any_changed;

	// Return this node as a BVHChild for the parent.
	BVHChild result;
	bvh_child_set_aabb(&result, bvh_node_aabb(dst));
	result.index = target_idx;
	result.leaf_count = dst->a.leaf_count + dst->b.leaf_count;
	return result;
}

static void bvh_refit(BVHTree* t, WorldInternal* w)
{
	if (t->root == -1) return;

	int cap = asize(t->nodes);
	BVHNode* new_nodes = CK_ALLOC(sizeof(BVHNode) * cap);
	BVHMeta* new_meta = CK_ALLOC(sizeof(BVHMeta) * cap);

	BVHRefit r = { t->nodes, t->meta, new_nodes, new_meta, t, w };
	int changed = 0;
	bvh_fused_recurse(&r, t->root, 0, -1, 0, &changed);

	// Swap new arrays into the tree.
	aclear(t->nodes); aclear(t->meta);
	int total = t->nodes[t->root].a.leaf_count + t->nodes[t->root].b.leaf_count;
	// Compute live count from the root's leaf_count (internal nodes = leaves - 1, plus root = leaves - 1).
	int live_count = new_nodes[0].a.leaf_count + new_nodes[0].b.leaf_count - 1;
	if (live_count < 1) live_count = 1;
	for (int i = 0; i < live_count; i++) { apush(t->nodes, new_nodes[i]); apush(t->meta, new_meta[i]); }
	aclear(t->node_free);
	t->root = 0;

	CK_FREE(new_nodes);
	CK_FREE(new_meta);
}

// -----------------------------------------------------------------------------
// Self-test: find all overlapping leaf pairs within a single tree.
// Split dispatch: leaf-vs-node, node-vs-node with 4-pair batching.

static void bvh_dispatch_pair(BVHTree* t, BVHChild* a, BVHChild* b, CK_DYNA BroadPair** pairs);

static void bvh_leaf_vs_node(BVHTree* t, BVHChild* leaf, int ni, CK_DYNA BroadPair** pairs)
{
	BVHNode* n = &t->nodes[ni];
	int oa = aabb_overlaps(bvh_child_aabb(leaf), bvh_child_aabb(&n->a));
	int ob = aabb_overlaps(bvh_child_aabb(leaf), bvh_child_aabb(&n->b));
	if (oa) {
		if (bvh_child_is_leaf(&n->a)) {
			BroadPair p = { t->leaves[bvh_child_leaf_idx(leaf)].body_idx, t->leaves[bvh_child_leaf_idx(&n->a)].body_idx };
			apush(*pairs, p);
		} else if (bvh_child_is_internal(&n->a)) {
			bvh_leaf_vs_node(t, leaf, n->a.index, pairs);
		}
	}
	if (ob) {
		if (bvh_child_is_leaf(&n->b)) {
			BroadPair p = { t->leaves[bvh_child_leaf_idx(leaf)].body_idx, t->leaves[bvh_child_leaf_idx(&n->b)].body_idx };
			apush(*pairs, p);
		} else if (bvh_child_is_internal(&n->b)) {
			bvh_leaf_vs_node(t, leaf, n->b.index, pairs);
		}
	}
}

static void bvh_nodes_vs_nodes(BVHTree* t, BVHNode* na, BVHNode* nb, CK_DYNA BroadPair** pairs)
{
	// Test all 4 child pairs up front while both nodes are in cache.
	int oo_aa = aabb_overlaps(bvh_child_aabb(&na->a), bvh_child_aabb(&nb->a));
	int oo_ab = aabb_overlaps(bvh_child_aabb(&na->a), bvh_child_aabb(&nb->b));
	int oo_ba = aabb_overlaps(bvh_child_aabb(&na->b), bvh_child_aabb(&nb->a));
	int oo_bb = aabb_overlaps(bvh_child_aabb(&na->b), bvh_child_aabb(&nb->b));
	if (oo_aa) bvh_dispatch_pair(t, &na->a, &nb->a, pairs);
	if (oo_ab) bvh_dispatch_pair(t, &na->a, &nb->b, pairs);
	if (oo_ba) bvh_dispatch_pair(t, &na->b, &nb->a, pairs);
	if (oo_bb) bvh_dispatch_pair(t, &na->b, &nb->b, pairs);
}

static void bvh_dispatch_pair(BVHTree* t, BVHChild* a, BVHChild* b, CK_DYNA BroadPair** pairs)
{
	if (bvh_child_is_internal(a) && bvh_child_is_internal(b)) {
		bvh_nodes_vs_nodes(t, &t->nodes[a->index], &t->nodes[b->index], pairs);
	} else if (bvh_child_is_leaf(a) && bvh_child_is_leaf(b)) {
		BroadPair p = { t->leaves[bvh_child_leaf_idx(a)].body_idx, t->leaves[bvh_child_leaf_idx(b)].body_idx };
		apush(*pairs, p);
	} else if (bvh_child_is_leaf(a) && bvh_child_is_internal(b)) {
		bvh_leaf_vs_node(t, a, b->index, pairs);
	} else if (bvh_child_is_internal(a) && bvh_child_is_leaf(b)) {
		bvh_leaf_vs_node(t, b, a->index, pairs);
	}
}

// Linear-scan self-test: iterate all nodes in reverse DFS order (Bepu-style).
// Requires nodes in DFS order (bvh_cache_reorder packs and clears freelist).
static void bvh_self_test(BVHTree* t, CK_DYNA BroadPair** pairs)
{
	if (t->root == -1) return;
	for (int i = asize(t->nodes) - 1; i >= 0; --i) {
		BVHNode* n = &t->nodes[i];
		if (aabb_overlaps(bvh_child_aabb(&n->a), bvh_child_aabb(&n->b)))
			bvh_dispatch_pair(t, &n->a, &n->b, pairs);
	}
}

// -----------------------------------------------------------------------------
// Cross-test: find all overlapping leaf pairs between two trees.
// Split dispatch with 4-pair batching for node-vs-node.

static void bvh_cross_dispatch(BVHTree* ta, BVHChild* a, BVHTree* tb, BVHChild* b, CK_DYNA BroadPair** pairs);

static void bvh_cross_leaf_vs_node(BVHTree* tl, BVHChild* leaf, BVHTree* tn, int ni, CK_DYNA BroadPair** pairs)
{
	BVHNode* n = &tn->nodes[ni];
	int oa = aabb_overlaps(bvh_child_aabb(leaf), bvh_child_aabb(&n->a));
	int ob = aabb_overlaps(bvh_child_aabb(leaf), bvh_child_aabb(&n->b));
	if (oa) {
		if (bvh_child_is_leaf(&n->a)) {
			BroadPair p = { tl->leaves[bvh_child_leaf_idx(leaf)].body_idx, tn->leaves[bvh_child_leaf_idx(&n->a)].body_idx };
			apush(*pairs, p);
		} else if (bvh_child_is_internal(&n->a)) {
			bvh_cross_leaf_vs_node(tl, leaf, tn, n->a.index, pairs);
		}
	}
	if (ob) {
		if (bvh_child_is_leaf(&n->b)) {
			BroadPair p = { tl->leaves[bvh_child_leaf_idx(leaf)].body_idx, tn->leaves[bvh_child_leaf_idx(&n->b)].body_idx };
			apush(*pairs, p);
		} else if (bvh_child_is_internal(&n->b)) {
			bvh_cross_leaf_vs_node(tl, leaf, tn, n->b.index, pairs);
		}
	}
}

static void bvh_cross_nodes(BVHTree* ta, BVHNode* na, BVHTree* tb, BVHNode* nb, CK_DYNA BroadPair** pairs)
{
	int oo_aa = aabb_overlaps(bvh_child_aabb(&na->a), bvh_child_aabb(&nb->a));
	int oo_ab = aabb_overlaps(bvh_child_aabb(&na->a), bvh_child_aabb(&nb->b));
	int oo_ba = aabb_overlaps(bvh_child_aabb(&na->b), bvh_child_aabb(&nb->a));
	int oo_bb = aabb_overlaps(bvh_child_aabb(&na->b), bvh_child_aabb(&nb->b));
	if (oo_aa) bvh_cross_dispatch(ta, &na->a, tb, &nb->a, pairs);
	if (oo_ab) bvh_cross_dispatch(ta, &na->a, tb, &nb->b, pairs);
	if (oo_ba) bvh_cross_dispatch(ta, &na->b, tb, &nb->a, pairs);
	if (oo_bb) bvh_cross_dispatch(ta, &na->b, tb, &nb->b, pairs);
}

static void bvh_cross_dispatch(BVHTree* ta, BVHChild* a, BVHTree* tb, BVHChild* b, CK_DYNA BroadPair** pairs)
{
	if (bvh_child_is_internal(a) && bvh_child_is_internal(b)) {
		bvh_cross_nodes(ta, &ta->nodes[a->index], tb, &tb->nodes[b->index], pairs);
	} else if (bvh_child_is_leaf(a) && bvh_child_is_leaf(b)) {
		BroadPair p = { ta->leaves[bvh_child_leaf_idx(a)].body_idx, tb->leaves[bvh_child_leaf_idx(b)].body_idx };
		apush(*pairs, p);
	} else if (bvh_child_is_leaf(a) && bvh_child_is_internal(b)) {
		bvh_cross_leaf_vs_node(ta, a, tb, b->index, pairs);
	} else if (bvh_child_is_internal(a) && bvh_child_is_leaf(b)) {
		bvh_cross_leaf_vs_node(tb, b, ta, a->index, pairs);
	}
}

static void bvh_cross_test(BVHTree* ta, BVHTree* tb, CK_DYNA BroadPair** pairs)
{
	if (ta->root == -1 || tb->root == -1) return;
	BVHNode* ra = &ta->nodes[ta->root];
	BVHNode* rb = &tb->nodes[tb->root];
	bvh_cross_nodes(ta, ra, tb, rb, pairs);
}

// ---- src/collision.c ----

// See LICENSE for licensing info.
// collision.c -- broadphase + narrowphase collision detection
//
// Hull SAT with Gauss map pruning for edge-edge axis elimination.
// Half-edge mesh enables efficient face/edge traversal for SAT queries.

// Types (HalfEdge, HullPlane, HullFace, Hull, ConvexHull) defined in nudge.h.

// Support function: furthest vertex along a direction.
static v3 hull_support(const Hull* hull, v3 dir)
{
	float best = -1e18f;
	int best_i = 0;
	for (int i = 0; i < hull->vert_count; i++) {
		float d = dot(hull->verts[i], dir);
		if (d > best) { best = d; best_i = i; }
	}
	return hull->verts[best_i];
}

// -----------------------------------------------------------------------------
// Unit box hull. Half-extents = (1,1,1). Scaled at query time.
//
// Edges stored in twin pairs: edge 2k and twin 2k+1.
// 8 verts, 24 half-edges (12 edges), 6 faces.

static const v3 s_box_verts[8] = {
	{-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
	{-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1},
};

// Face winding (CCW from outside):
//  0: -Z (0,3,2,1)   1: +Z (4,5,6,7)
//  2: -X (0,4,7,3)   3: +X (1,2,6,5)
//  4: -Y (0,1,5,4)   5: +Y (3,7,6,2)
//
// 12 undirected edges, each stored as (edge, twin) pair at indices (2k, 2k+1):
//  pair 0: 0->3 / 3->0     pair 1: 3->2 / 2->3
//  pair 2: 2->1 / 1->2     pair 3: 1->0 / 0->1
//  pair 4: 4->5 / 5->4     pair 5: 5->6 / 6->5
//  pair 6: 6->7 / 7->6     pair 7: 7->4 / 4->7
//  pair 8: 0->4 / 4->0     pair 9: 1->5 / 5->1
//  pair10: 2->6 / 6->2     pair11: 3->7 / 7->3

static const HalfEdge s_box_edges[24] = {
	// pair 0: 0->3 (face 0) / 3->0 (face 2)
	{ .twin =  1, .next =  2, .origin = 0, .face = 0 },  // e0:  0->3
	{ .twin =  0, .next = 16, .origin = 3, .face = 2 },  // e1:  3->0
	// pair 1: 3->2 (face 0) / 2->3 (face 5)
	{ .twin =  3, .next =  4, .origin = 3, .face = 0 },  // e2:  3->2
	{ .twin =  2, .next = 22, .origin = 2, .face = 5 },  // e3:  2->3
	// pair 2: 2->1 (face 0) / 1->2 (face 3)
	{ .twin =  5, .next =  6, .origin = 2, .face = 0 },  // e4:  2->1
	{ .twin =  4, .next = 20, .origin = 1, .face = 3 },  // e5:  1->2
	// pair 3: 1->0 (face 0) / 0->1 (face 4)
	{ .twin =  7, .next =  0, .origin = 1, .face = 0 },  // e6:  1->0
	{ .twin =  6, .next = 18, .origin = 0, .face = 4 },  // e7:  0->1
	// pair 4: 4->5 (face 1) / 5->4 (face 4)
	{ .twin =  9, .next = 10, .origin = 4, .face = 1 },  // e8:  4->5
	{ .twin =  8, .next = 17, .origin = 5, .face = 4 },  // e9:  5->4
	// pair 5: 5->6 (face 1) / 6->5 (face 3)
	{ .twin = 11, .next = 12, .origin = 5, .face = 1 },  // e10: 5->6
	{ .twin = 10, .next = 19, .origin = 6, .face = 3 },  // e11: 6->5
	// pair 6: 6->7 (face 1) / 7->6 (face 5)
	{ .twin = 13, .next = 14, .origin = 6, .face = 1 },  // e12: 6->7
	{ .twin = 12, .next = 21, .origin = 7, .face = 5 },  // e13: 7->6
	// pair 7: 7->4 (face 1) / 4->7 (face 2)
	{ .twin = 15, .next =  8, .origin = 7, .face = 1 },  // e14: 7->4
	{ .twin = 14, .next = 23, .origin = 4, .face = 2 },  // e15: 4->7
	// pair 8: 0->4 (face 2) / 4->0 (face 4)
	{ .twin = 17, .next = 15, .origin = 0, .face = 2 },  // e16: 0->4
	{ .twin = 16, .next =  7, .origin = 4, .face = 4 },  // e17: 4->0
	// pair 9: 1->5 (face 4) / 5->1 (face 3)
	{ .twin = 19, .next =  9, .origin = 1, .face = 4 },  // e18: 1->5
	{ .twin = 18, .next =  5, .origin = 5, .face = 3 },  // e19: 5->1
	// pair 10: 2->6 (face 3) / 6->2 (face 5)
	{ .twin = 21, .next = 11, .origin = 2, .face = 3 },  // e20: 2->6
	{ .twin = 20, .next =  3, .origin = 6, .face = 5 },  // e21: 6->2
	// pair 11: 3->7 (face 5) / 7->3 (face 2)
	{ .twin = 23, .next = 13, .origin = 3, .face = 5 },  // e22: 3->7
	{ .twin = 22, .next =  1, .origin = 7, .face = 2 },  // e23: 7->3
};

static const HullFace s_box_faces[6] = {
	{ .edge =  0 },  // face 0 (-Z): starts at e0 (0->3)
	{ .edge =  8 },  // face 1 (+Z): starts at e8 (4->5)
	{ .edge = 16 },  // face 2 (-X): starts at e16 (0->4)
	{ .edge =  5 },  // face 3 (+X): starts at e5 (1->2)
	{ .edge =  7 },  // face 4 (-Y): starts at e7 (0->1)
	{ .edge = 22 },  // face 5 (+Y): starts at e22 (3->7)
};

static const HullPlane s_box_planes[6] = {
	{ .normal = { 0, 0,-1}, .offset = 1 },
	{ .normal = { 0, 0, 1}, .offset = 1 },
	{ .normal = {-1, 0, 0}, .offset = 1 },
	{ .normal = { 1, 0, 0}, .offset = 1 },
	{ .normal = { 0,-1, 0}, .offset = 1 },
	{ .normal = { 0, 1, 0}, .offset = 1 },
};

static const Hull s_unit_box_hull = {
	.centroid = {0, 0, 0},
	.verts = s_box_verts,
	.edges = s_box_edges,
	.faces = s_box_faces,
	.planes = s_box_planes,
	.vert_count = 8,
	.edge_count = 24,
	.face_count = 6,
	.epsilon = 9.0f * FLT_EPSILON,
};

// -----------------------------------------------------------------------------
// Transform helpers for hull queries.

// Scale a hull vertex by half-extents.
static inline v3 hull_vert_scaled(const Hull* hull, int i, v3 scale)
{
	v3 v = hull->verts[i];
	return (v3){ v.x * scale.x, v.y * scale.y, v.z * scale.z };
}

// Transform a plane from hull-local into world space.
static inline HullPlane plane_transform(HullPlane p, v3 pos, quat rot, v3 scale)
{
	// For axis-aligned unit box with uniform-ish scale:
	// normal needs inverse-transpose scale, then rotate.
	v3 n = norm(rotate(rot, V3(p.normal.x / scale.x, p.normal.y / scale.y, p.normal.z / scale.z)));
	// A point on the plane in local space: normal * offset, then scale + transform
	v3 local_pt = V3(p.normal.x * p.offset * scale.x, p.normal.y * p.offset * scale.y, p.normal.z * p.offset * scale.z);
	v3 world_pt = add(pos, rotate(rot, local_pt));
	return (HullPlane){ .normal = n, .offset = dot(n, world_pt) };
}

// Project hull onto a plane, return signed distance of support point.
static float hull_project_plane(const Hull* hull, HullPlane plane, v3 scale)
{
	v3 support = hull_support(hull, neg(plane.normal));
	v3 sv = { support.x * scale.x, support.y * scale.y, support.z * scale.z };
	return dot(plane.normal, sv) - plane.offset;
}

// -----------------------------------------------------------------------------
// Internal manifold used by broadphase (extends public Manifold with body indices).

typedef struct InternalManifold
{
	Manifold m;
	int body_a;
	int body_b;
} InternalManifold;

// -----------------------------------------------------------------------------
// Sphere-sphere.

int collide_sphere_sphere(Sphere a, Sphere b, Manifold* manifold)
{
	v3 d = sub(b.center, a.center);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
	float penetration = r_sum - dist;

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(a.center, scale(normal, a.radius - penetration * 0.5f)),
		.normal = normal,
		.penetration = penetration,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Segment helpers for capsule collisions.

// Closest point on segment PQ to point X. Returns parametric t in [0,1].
static float segment_closest_t(v3 P, v3 Q, v3 X)
{
	v3 d = sub(Q, P);
	float d_len2 = len2(d);
	if (d_len2 < 1e-12f) return 0.0f;
	float t = dot(sub(X, P), d) / d_len2;
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;
	return t;
}

static v3 segment_closest_point(v3 P, v3 Q, v3 X)
{
	float t = segment_closest_t(P, Q, X);
	return add(P, scale(sub(Q, P), t));
}

// Closest points between two segments P1Q1 and P2Q2.
static void segments_closest_points(v3 P1, v3 Q1, v3 P2, v3 Q2, v3* out1, v3* out2)
{
	v3 d1 = sub(Q1, P1);
	v3 d2 = sub(Q2, P2);
	v3 r = sub(P1, P2);
	float a = dot(d1, d1);
	float e = dot(d2, d2);
	float f = dot(d2, r);
	float s, t;

	if (a < 1e-12f && e < 1e-12f) {
		*out1 = P1; *out2 = P2; return;
	}
	if (a < 1e-12f) {
		s = 0.0f;
		t = f / e;
		if (t < 0.0f) t = 0.0f; if (t > 1.0f) t = 1.0f;
	} else {
		float c = dot(d1, r);
		if (e < 1e-12f) {
			t = 0.0f;
			s = -c / a;
			if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f;
		} else {
			float b = dot(d1, d2);
			float denom = a * e - b * b;
			s = denom > 1e-12f ? (b * f - c * e) / denom : 0.0f;
			if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f;
			t = (b * s + f) / e;
			if (t < 0.0f) { t = 0.0f; s = -c / a; if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f; }
			else if (t > 1.0f) { t = 1.0f; s = (b - c) / a; if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f; }
		}
	}
	*out1 = add(P1, scale(d1, s));
	*out2 = add(P2, scale(d2, t));
}

// Get capsule segment endpoints in world space.
static void capsule_world_segment(BodyHot* h, ShapeInternal* s, v3* P, v3* Q)
{
	v3 local_p = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 local_q = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	*P = add(h->position, rotate(h->rotation, local_p));
	*Q = add(h->position, rotate(h->rotation, local_q));
}

// -----------------------------------------------------------------------------
// Sphere-capsule.

int collide_sphere_capsule(Sphere a, Capsule b, Manifold* manifold)
{
	v3 closest = segment_closest_point(b.p, b.q, a.center);
	v3 d = sub(closest, a.center);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(a.center, scale(normal, a.radius)),
		.normal = normal,
		.penetration = r_sum - dist,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// Capsule-capsule.

int collide_capsule_capsule(Capsule a, Capsule b, Manifold* manifold)
{
	v3 c1, c2;
	segments_closest_points(a.p, a.q, b.p, b.q, &c1, &c2);

	v3 d = sub(c2, c1);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(c1, scale(normal, a.radius)),
		.normal = normal,
		.penetration = r_sum - dist,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// Capsule-box (hull) narrowphase.
// Shallow: GJK witness points. Deep: face search (lm-style).
// GJK distance is fuzzy near zero -- use layered thresholds.
// LINEAR_SLOP defined in nudge.c (solver constants)

// GJK query helpers. Hull proxies need a temp buffer for pre-scaled verts.
#define MAX_HULL_VERTS 256

static GJK_Result gjk_query_point_hull(v3 pt, ConvexHull h)
{
	v3 scaled[MAX_HULL_VERTS];
	GJK_Proxy pa, pb;
	gjk_proxy_point(&pa, pt);
	gjk_proxy_hull(&pb, h.hull, h.center, h.rotation, h.scale, scaled);
	return gjk_distance_ex(&pa, &pb);
}

static GJK_Result gjk_query_segment_hull(v3 p, v3 q, ConvexHull h)
{
	v3 scaled[MAX_HULL_VERTS];
	GJK_Proxy pa, pb;
	gjk_proxy_segment(&pa, p, q);
	gjk_proxy_hull(&pb, h.hull, h.center, h.rotation, h.scale, scaled);
	return gjk_distance_ex(&pa, &pb);
}

static GJK_Result gjk_query_hull_hull(ConvexHull a, ConvexHull b)
{
	v3 sa[MAX_HULL_VERTS], sb[MAX_HULL_VERTS];
	GJK_Proxy pa, pb;
	gjk_proxy_hull(&pa, a.hull, a.center, a.rotation, a.scale, sa);
	gjk_proxy_hull(&pb, b.hull, b.center, b.rotation, b.scale, sb);
	return gjk_distance_ex(&pa, &pb);
}


// Deep penetration: find most-separated face on hull from a point.
static int hull_deepest_face(const Hull* hull, v3 local_pt)
{
	float best = -1e18f;
	int face = 0;
	for (int i = 0; i < hull->face_count; i++) {
		float s = dot(local_pt, hull->planes[i].normal) - hull->planes[i].offset;
		if (s > best) { best = s; face = i; }
	}
	return face;
}

// Transform a world point to hull local (unit) space.
static v3 to_hull_local(v3 pt, v3 pos, quat rot, v3 sc)
{
	v3 lp = rotate(inv(rot), sub(pt, pos));
	return V3(lp.x / sc.x, lp.y / sc.y, lp.z / sc.z);
}

// Sphere-hull: GJK on sphere center (point) vs hull for distance,
// then shallow path (GJK witness) or deep path (face search).
int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold)
{
	// GJK on the sphere CENTER (point) vs hull core.
	GJK_Result r = gjk_query_point_hull(a.center, b);

	if (r.distance > a.radius) return 0; // separated

	if (r.distance > LINEAR_SLOP) {
		// Shallow: center is outside hull but within radius.
		v3 normal = scale(sub(r.point2, r.point1), 1.0f / r.distance);
		if (!manifold) return 1;
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = add(a.center, scale(normal, a.radius)),
			.normal = normal,
			.penetration = a.radius - r.distance,
		};
		return 1;
	}

	// Deep: center is inside hull. Find face with least penetration in world space.
	// Transform each hull plane to world space (accounts for non-uniform scale).
	float best_sep = -1e18f;
	int best_face = -1;
	HullPlane best_plane = {0};
	for (int i = 0; i < b.hull->face_count; i++) {
		HullPlane wp = plane_transform(b.hull->planes[i], b.center, b.rotation, b.scale);
		float s = dot(a.center, wp.normal) - wp.offset;
		if (s > a.radius) return 0;
		if (s > best_sep) { best_sep = s; best_face = i; best_plane = wp; }
	}
	if (!manifold) return 1;

	v3 world_pt = sub(a.center, scale(best_plane.normal, best_sep));
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = world_pt,
		.normal = neg(best_plane.normal),
		.penetration = a.radius - best_sep,
	};
	return 1;
}

int collide_sphere_box(Sphere a, Box b, Manifold* manifold)
{
	return collide_sphere_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// Capsule-hull: GJK shallow path, face/edge-search deep path.
int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold)
{
	GJK_Result r = gjk_query_segment_hull(a.p, a.q, b);

	if (r.distance > a.radius) return 0;

	if (r.distance > LINEAR_SLOP) {
		// Shallow: GJK witness points.
		v3 normal = scale(sub(r.point2, r.point1), 1.0f / r.distance);
		v3 contact_pt = add(r.point1, scale(normal, a.radius));
		if (!manifold) return 1;
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = contact_pt,
			.normal = normal,
			.penetration = a.radius - r.distance,
		};
		return 1;
	}

	// Deep: SAT on hull faces + capsule edge axes. Work in world space to
	// avoid non-uniform scale distortion (lm pattern).
	const Hull* hull = b.hull;
	v3 cap_dir = sub(a.q, a.p);
	float cap_len2 = len2(cap_dir);

	// --- Axis family 1: hull face normals (world space) ---
	float face_sep = -1e18f;
	int face_idx = -1;
	HullPlane face_plane = {0};
	for (int i = 0; i < hull->face_count; i++) {
		HullPlane wp = plane_transform(hull->planes[i], b.center, b.rotation, b.scale);
		float dp = dot(a.p, wp.normal) - wp.offset;
		float dq = dot(a.q, wp.normal) - wp.offset;
		float sup = dp < dq ? dp : dq;
		if (sup > face_sep) { face_sep = sup; face_idx = i; face_plane = wp; }
	}

	// --- Axis family 2: capsule dir x hull edge dirs (world space) ---
	float edge_sep = -1e18f;
	int edge_idx = -1;
	v3 edge_axis = V3(0,0,0);
	v3 edge_pt_world = V3(0,0,0);
	if (cap_len2 > 1e-12f) {
		v3 cd = scale(cap_dir, 1.0f / sqrtf(cap_len2));
		for (int i = 0; i < hull->edge_count; i += 2) {
			v3 ev0 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edges[i].origin], b.scale)));
			v3 ev1 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edges[hull->edges[i].next].origin], b.scale)));
			v3 ed = sub(ev1, ev0);
			v3 ax = cross(cd, ed);
			float al = len2(ax);
			if (al < 1e-12f) continue;
			ax = scale(ax, 1.0f / sqrtf(al));
			if (dot(ax, sub(a.p, ev0)) < 0.0f) ax = neg(ax);
			float cs = dot(a.p, ax) < dot(a.q, ax) ? dot(a.p, ax) : dot(a.q, ax);
			float hs = -1e18f;
			for (int j = 0; j < hull->vert_count; j++) {
				v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
				float d = dot(wv, ax);
				if (d > hs) hs = d;
			}
			float sep = cs - hs;
			if (sep > edge_sep) { edge_sep = sep; edge_idx = i; edge_axis = ax; edge_pt_world = ev0; }
		}
	}

	// Pick axis of minimum penetration.
	float best_sep;
	v3 best_n;
	int use_face;
	if (edge_idx >= 0 && edge_sep > face_sep + 0.001f) {
		best_sep = edge_sep;
		best_n = edge_axis;
		use_face = 0;
	} else {
		best_sep = face_sep;
		best_n = face_plane.normal;
		use_face = 1;
	}

	if (best_sep > a.radius) return 0;
	if (!manifold) return 1;

	// Generate contacts: project capsule endpoints onto the reference plane,
	// keep those that penetrate.
	float plane_d = use_face ? face_plane.offset : dot(best_n, edge_pt_world);
	float dp = dot(a.p, best_n) - plane_d;
	float dq = dot(a.q, best_n) - plane_d;

	int cp = 0;
	v3 points[2];
	float depths[2];
	if (dp < a.radius) {
		points[cp] = sub(a.p, scale(best_n, a.radius));
		depths[cp] = a.radius - dp;
		cp++;
	}
	if (dq < a.radius) {
		points[cp] = sub(a.q, scale(best_n, a.radius));
		depths[cp] = a.radius - dq;
		cp++;
	}

	if (cp == 0) return 0;
	manifold->count = cp;
	for (int i = 0; i < cp; i++) {
		manifold->contacts[i] = (Contact){
			.point = points[i],
			.normal = neg(best_n),
			.penetration = depths[i],
		};
	}
	return 1;
}

int collide_capsule_box(Capsule a, Box b, Manifold* manifold)
{
	return collide_capsule_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// -----------------------------------------------------------------------------
// SAT: face queries.
// All computations in local space of Hull2.

typedef struct FaceQuery
{
	int index;
	float separation;
} FaceQuery;

static FaceQuery sat_query_faces(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
{
	// Transform from world to hull2 local: p_local = rot2^-1 * (p_world - pos2)
	quat inv2 = inv(rot2);

	FaceQuery best = { .index = -1, .separation = -1e18f };

	for (int i = 0; i < hull1->face_count; i++) {
		// Transform hull1's plane into hull2's local space
		HullPlane pw = plane_transform(hull1->planes[i], pos1, rot1, scale1);
		// Bring plane normal into hull2 local
		v3 local_n = rotate(inv2, pw.normal);
		v3 local_pt = rotate(inv2, sub(scale(pw.normal, pw.offset), pos2));
		// Actually, let's just use world space and project hull2's support
		// Support of hull2 in direction -pw.normal (world space)
		v3 sup_dir_local = rotate(inv2, neg(pw.normal));
		v3 sup_local = hull_support(hull2, sup_dir_local);
		v3 sup_scaled = { sup_local.x * scale2.x, sup_local.y * scale2.y, sup_local.z * scale2.z };
		v3 sup_world = add(pos2, rotate(rot2, sup_scaled));

		float sep = dot(pw.normal, sup_world) - pw.offset;
		if (sep > best.separation) {
			best.separation = sep;
			best.index = i;
		}
	}
	return best;
}

// -----------------------------------------------------------------------------
// SAT: Gauss map Minkowski face test.
// Tests if arcs AB and CD intersect on the unit sphere.

static int is_minkowski_face(v3 a, v3 b, v3 b_x_a, v3 c, v3 d, v3 d_x_c)
{
	float cba = dot(c, b_x_a);
	float dba = dot(d, b_x_a);
	float adc = dot(a, d_x_c);
	float bdc = dot(b, d_x_c);
	return (cba * dba < 0.0f) && (adc * bdc < 0.0f) && (cba * bdc > 0.0f);
}

// Project edge pair: signed distance along cross(e1,e2) from p1 to p2.
static float sat_edge_project(v3 p1, v3 e1, v3 p2, v3 e2, v3 c1)
{
	v3 e1_x_e2 = cross(e1, e2);
	float len = len(e1_x_e2);

	// Skip near-parallel edges
	float tolerance = 0.005f;
	if (len < tolerance * sqrtf(len2(e1) * len2(e2)))
		return -1e18f;

	v3 n = scale(e1_x_e2, 1.0f / len);
	// Ensure consistent orientation (hull1 -> hull2)
	if (dot(n, sub(p1, c1)) < 0.0f)
		n = neg(n);

	return dot(n, sub(p2, p1));
}

// SAT: edge queries with Gauss map pruning.
typedef struct EdgeQuery
{
	int index1;
	int index2;
	float separation;
} EdgeQuery;

static EdgeQuery sat_query_edges(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
{
	// All in local space of hull2.
	quat inv2 = inv(rot2);
	quat rel_rot = mul(inv2, rot1); // rotation from hull1-local to hull2-local
	v3 c1_local = rotate(inv2, sub(pos1, pos2)); // hull1 centroid in hull2 space

	EdgeQuery best = { .index1 = -1, .index2 = -1, .separation = -1e18f };

	for (int i1 = 0; i1 < hull1->edge_count; i1 += 2) {
		const HalfEdge* edge1 = &hull1->edges[i1];
		const HalfEdge* twin1 = &hull1->edges[i1 + 1];

		v3 p1 = hull_vert_scaled(hull1, edge1->origin, scale1);
		v3 q1 = hull_vert_scaled(hull1, twin1->origin, scale1);
		// Transform to hull2 local
		p1 = add(c1_local, rotate(rel_rot, p1));
		q1 = add(c1_local, rotate(rel_rot, q1));
		v3 e1 = sub(q1, p1);

		v3 u1 = rotate(rel_rot, hull1->planes[edge1->face].normal);
		v3 v1 = rotate(rel_rot, hull1->planes[twin1->face].normal);

		for (int i2 = 0; i2 < hull2->edge_count; i2 += 2) {
			const HalfEdge* edge2 = &hull2->edges[i2];
			const HalfEdge* twin2 = &hull2->edges[i2 + 1];

			v3 p2 = hull_vert_scaled(hull2, edge2->origin, scale2);
			v3 q2 = hull_vert_scaled(hull2, twin2->origin, scale2);
			v3 e2 = sub(q2, p2);

			v3 u2 = hull2->planes[edge2->face].normal;
			v3 v2 = hull2->planes[twin2->face].normal;

			// Gauss map pruning
			if (!is_minkowski_face(u1, v1, neg(e1), neg(u2), neg(v2), neg(e2)))
				continue;

			float sep = sat_edge_project(p1, e1, p2, e2, c1_local);
			if (sep > best.separation) {
				best.index1 = i1;
				best.index2 = i2;
				best.separation = sep;
			}
		}
	}
	return best;
}

// -----------------------------------------------------------------------------
// Contact manifold helpers for hull-hull clipping.

#define MAX_CLIP_VERTS 64

// Collect face vertices in world space by walking the half-edge loop.
static int hull_face_verts_world(const Hull* hull, int face_idx, v3 pos, quat rot, v3 sc, v3* out)
{
	int start = hull->faces[face_idx].edge;
	int e = start;
	int count = 0;
	do {
		out[count++] = add(pos, rotate(rot, hull_vert_scaled(hull, hull->edges[e].origin, sc)));
		e = hull->edges[e].next;
	} while (e != start && count < MAX_CLIP_VERTS);
	return count;
}

// Find face on hull most anti-parallel to a world-space normal.
static int find_incident_face(const Hull* hull, v3 pos, quat rot, v3 sc, v3 ref_normal)
{
	int best = 0;
	float best_dot = 1e18f;
	for (int i = 0; i < hull->face_count; i++) {
		HullPlane pw = plane_transform(hull->planes[i], pos, rot, sc);
		float d = dot(pw.normal, ref_normal);
		if (d < best_dot) { best_dot = d; best = i; }
	}
	return best;
}

// Sutherland-Hodgman: clip polygon against a single plane.
// Keeps points on the negative side: dot(plane_n, p) - plane_d <= 0.
// Tracks feature IDs: clip_edge is the side plane index that clips new verts.
static int clip_to_plane(v3* in, uint8_t* in_fid, int in_count, v3 plane_n, float plane_d, uint8_t clip_edge, v3* out, uint8_t* out_fid)
{
	if (in_count == 0) return 0;
	int out_count = 0;
	v3 prev = in[in_count - 1];
	uint8_t fid_prev = in_fid[in_count - 1];
	float d_prev = dot(plane_n, prev) - plane_d;

	for (int i = 0; i < in_count; i++) {
		v3 cur = in[i];
		uint8_t fid_cur = in_fid[i];
		float d_cur = dot(plane_n, cur) - plane_d;

		if (d_prev <= 0.0f) {
			out[out_count] = prev;
			out_fid[out_count] = fid_prev;
			out_count++;
			if (d_cur > 0.0f) {
				float t = d_prev / (d_prev - d_cur);
				out[out_count] = add(prev, scale(sub(cur, prev), t));
				out_fid[out_count] = clip_edge; // new vertex from this clip plane
				out_count++;
			}
		} else if (d_cur <= 0.0f) {
			float t = d_prev / (d_prev - d_cur);
			out[out_count] = add(prev, scale(sub(cur, prev), t));
			out_fid[out_count] = clip_edge;
			out_count++;
		}

		prev = cur;
		fid_prev = fid_cur;
		d_prev = d_cur;
	}
	return out_count;
}

// Reduce contact manifold to MAX_CONTACTS (4) points.
// 1) Two farthest points (span the patch).
// 2) Farthest from that segment (adds width).
// 3) Maximal barycentric contributor (most outside triangle, completes quad).
static int reduce_contacts(Contact* contacts, int count)
{
	if (count <= MAX_CONTACTS) return count;

	int sel[MAX_CONTACTS];
	int used[MAX_CLIP_VERTS];
	memset(used, 0, count * sizeof(int));

	// Step 1: two farthest points
	float best_d2 = -1.0f;
	int i0 = 0, i1 = 1;
	for (int i = 0; i < count; i++)
		for (int j = i + 1; j < count; j++) {
			float d2 = len2(sub(contacts[j].point, contacts[i].point));
			if (d2 > best_d2) { best_d2 = d2; i0 = i; i1 = j; }
		}
	sel[0] = i0; sel[1] = i1;
	used[i0] = used[i1] = 1;

	// Step 2: farthest from line through (i0, i1)
	v3 seg = sub(contacts[i1].point, contacts[i0].point);
	float seg_l2 = len2(seg);
	float best = -1.0f;
	int i2 = -1;
	for (int i = 0; i < count; i++) {
		if (used[i]) continue;
		v3 v = sub(contacts[i].point, contacts[i0].point);
		float proj = seg_l2 > 1e-12f ? dot(v, seg) / seg_l2 : 0.0f;
		float d2 = len2(sub(v, scale(seg, proj)));
		if (d2 > best) { best = d2; i2 = i; }
	}
	sel[2] = i2;
	used[i2] = 1;

	// Step 3: maximal barycentric contributor (most outside triangle)
	v3 A = contacts[sel[0]].point;
	v3 B = contacts[sel[1]].point;
	v3 C = contacts[sel[2]].point;
	v3 N = cross(sub(B, A), sub(C, A));
	float best_min = 1e18f;
	int i3 = -1;
	for (int i = 0; i < count; i++) {
		if (used[i]) continue;
		v3 P = contacts[i].point;
		float u = dot(cross(sub(B, P), sub(C, P)), N);
		float vc = dot(cross(sub(C, P), sub(A, P)), N);
		float w = dot(cross(sub(A, P), sub(B, P)), N);
		float m = u < vc ? u : vc;
		m = m < w ? m : w;
		if (m < best_min) { best_min = m; i3 = i; }
	}
	sel[3] = i3;

	Contact tmp[MAX_CONTACTS];
	for (int i = 0; i < MAX_CONTACTS; i++) tmp[i] = contacts[sel[i]];
	for (int i = 0; i < MAX_CONTACTS; i++) contacts[i] = tmp[i];
	return MAX_CONTACTS;
}

// -----------------------------------------------------------------------------
// Full SAT hull vs hull with Sutherland-Hodgman face clipping.

int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold)
{
	const Hull* hull_a = a.hull;
	v3 pos_a = a.center; quat rot_a = a.rotation; v3 scale_a = a.scale;
	const Hull* hull_b = b.hull;
	v3 pos_b = b.center; quat rot_b = b.rotation; v3 scale_b = b.scale;

	FaceQuery face_a = sat_query_faces(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b);
	if (face_a.separation > 0.0f) return 0;

	FaceQuery face_b = sat_query_faces(hull_b, pos_b, rot_b, scale_b, hull_a, pos_a, rot_a, scale_a);
	if (face_b.separation > 0.0f) return 0;

	// NaN transforms cause face_index to stay -1 (NaN comparisons always false)
	if (face_a.index < 0 || face_b.index < 0) return 0;

	EdgeQuery edge_q = sat_query_edges(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b);
	if (edge_q.separation > 0.0f) return 0;

	if (!manifold) return 1;

	// Bias toward face contacts over edge contacts
	const float k_tol = 0.05f;
	float max_face_sep = face_a.separation > face_b.separation
		? face_a.separation : face_b.separation;

	if (edge_q.separation > max_face_sep + k_tol) {
		// --- Edge-edge contact ---
		const HalfEdge* e1 = &hull_a->edges[edge_q.index1];
		const HalfEdge* t1 = &hull_a->edges[edge_q.index1 + 1];
		v3 p1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, e1->origin, scale_a)));
		v3 q1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, t1->origin, scale_a)));

		const HalfEdge* e2 = &hull_b->edges[edge_q.index2];
		const HalfEdge* t2 = &hull_b->edges[edge_q.index2 + 1];
		v3 p2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, e2->origin, scale_b)));
		v3 q2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, t2->origin, scale_b)));

		v3 ca, cb;
		segments_closest_points(p1, q1, p2, q2, &ca, &cb);

		v3 normal = norm(cross(sub(q1, p1), sub(q2, p2)));
		if (dot(normal, sub(pos_b, pos_a)) < 0.0f) normal = neg(normal);

		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = scale(add(ca, cb), 0.5f),
			.normal = normal,
			.penetration = -edge_q.separation,
			.feature_id = FEATURE_EDGE_BIT | (uint32_t)edge_q.index1 | ((uint32_t)edge_q.index2 << 16),
		};
		return 1;
	}

	// --- Face contact with Sutherland-Hodgman clipping ---
	const Hull* ref_hull; const Hull* inc_hull;
	v3 ref_pos, inc_pos, ref_sc, inc_sc;
	quat ref_rot, inc_rot;
	int ref_face, flip;

	// Bias toward shape A as reference face: B must be significantly better to win.
	// This keeps feature IDs stable across frames for nearly-parallel face pairs.
	float face_bias = 0.98f * face_a.separation + 0.001f;
	if (face_b.separation > face_bias) {
		ref_hull = hull_b; ref_pos = pos_b; ref_rot = rot_b; ref_sc = scale_b;
		inc_hull = hull_a; inc_pos = pos_a; inc_rot = rot_a; inc_sc = scale_a;
		ref_face = face_b.index; flip = 1;
	} else {
		ref_hull = hull_a; ref_pos = pos_a; ref_rot = rot_a; ref_sc = scale_a;
		inc_hull = hull_b; inc_pos = pos_b; inc_rot = rot_b; inc_sc = scale_b;
		ref_face = face_a.index; flip = 0;
	}

	HullPlane ref_plane = plane_transform(ref_hull->planes[ref_face], ref_pos, ref_rot, ref_sc);

	// Collect incident face vertices with initial feature IDs.
	int inc_face = find_incident_face(inc_hull, inc_pos, inc_rot, inc_sc, ref_plane.normal);
	v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
	uint8_t fid1[MAX_CLIP_VERTS], fid2[MAX_CLIP_VERTS];
	int clip_count = hull_face_verts_world(inc_hull, inc_face, inc_pos, inc_rot, inc_sc, buf1);
	// Initial feature: each incident vertex gets 0xFF (original, not clipped).
	for (int i = 0; i < clip_count; i++) fid1[i] = 0xFF;

	// Clip against each side plane of the reference face.
	v3* in_buf = buf1;   v3* out_buf = buf2;
	uint8_t* in_fid = fid1; uint8_t* out_fid = fid2;
	int start_e = ref_hull->faces[ref_face].edge;
	int ei = start_e;
	int guard = 0;
	uint8_t clip_edge_idx = 0;
	do {
		const HalfEdge* edge = &ref_hull->edges[ei];
		v3 tail = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, edge->origin, ref_sc)));
		v3 head = add(ref_pos, rotate(ref_rot,
			hull_vert_scaled(ref_hull, ref_hull->edges[edge->twin].origin, ref_sc)));
		v3 side_n = norm(cross(sub(head, tail), ref_plane.normal));
		float side_d = dot(side_n, tail);

		clip_count = clip_to_plane(in_buf, in_fid, clip_count,
			side_n, side_d, clip_edge_idx, out_buf, out_fid);
		v3* swap = in_buf; in_buf = out_buf; out_buf = swap;
		uint8_t* fswap = in_fid; in_fid = out_fid; out_fid = fswap;

		clip_edge_idx++;
		ei = edge->next;
		assert(++guard < MAX_CLIP_VERTS && "collide_hull_hull: face edge loop didn't close");
	} while (ei != start_e);

	// Snap clipped vertices near reference face corners to canonical corner IDs.
	// A vertex at a corner sits on two side planes -- which plane's clip_edge it
	// gets is FP-dependent and flickers between frames. Replacing the clip_edge
	// with a deterministic corner tag (0xC0 | corner_index) stabilises the
	// feature ID so warm starting can match it.
	{
		v3 corners[MAX_CLIP_VERTS];
		int ncorners = 0;
		int ce = start_e;
		do {
			corners[ncorners++] = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edges[ce].origin, ref_sc)));
			ce = ref_hull->edges[ce].next;
		} while (ce != start_e);
		float snap_tol2 = 1e-6f;
		for (int i = 0; i < clip_count; i++) {
			if (in_fid[i] == 0xFF) continue;
			for (int c = 0; c < ncorners; c++) {
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) {
					in_fid[i] = 0xC0 | (uint8_t)c;
					break;
				}
			}
		}
	}

	// Keep points within margin of reference face, generate contacts with feature IDs.
	// LINEAR_SLOP margin keeps contacts alive near zero-penetration, preventing blink.
	// Feature ID encodes: ref_face | (inc_face << 8) | (clip_edge << 16).
	// flip bit: if ref was hull_b, swap the face roles so the ID is canonical.
	v3 contact_n = flip ? neg(ref_plane.normal) : ref_plane.normal;
	Contact tmp_contacts[MAX_CLIP_VERTS];
	int cp = 0;
	for (int i = 0; i < clip_count; i++) {
		float depth = ref_plane.offset - dot(ref_plane.normal, in_buf[i]);
		if (depth >= -LINEAR_SLOP) {
			uint32_t fid;
			if (!flip)
				fid = (uint32_t)ref_face | ((uint32_t)inc_face << 8) | ((uint32_t)in_fid[i] << 16);
			else
				fid = (uint32_t)inc_face | ((uint32_t)ref_face << 8) | ((uint32_t)in_fid[i] << 16);
			tmp_contacts[cp++] = (Contact){
				.point = in_buf[i],
				.normal = contact_n,
				.penetration = depth,
				.feature_id = fid,
			};
		}
	}

	if (cp == 0) return 0;

	cp = reduce_contacts(tmp_contacts, cp);
	manifold->count = cp;
	for (int i = 0; i < cp; i++)
		manifold->contacts[i] = tmp_contacts[i];
	return 1;
}

// -----------------------------------------------------------------------------
int collide_box_box(Box a, Box b, Manifold* manifold)
{
	return collide_hull_hull(
		(ConvexHull){ &s_unit_box_hull, a.center, a.rotation, a.half_extents },
		(ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents },
		manifold);
}

const Hull* hull_unit_box() { return &s_unit_box_hull; }

// Quickhull implemented in quickhull.c.

void hull_free(Hull* hull)
{
	if (!hull) return;
	CK_FREE((void*)hull->verts);
	CK_FREE((void*)hull->edges);
	CK_FREE((void*)hull->faces);
	CK_FREE((void*)hull->planes);
	CK_FREE(hull);
}

// -----------------------------------------------------------------------------
// N^2 broadphase + narrowphase dispatch.

// Build a Sphere/Capsule/Box from internal body+shape for broadphase dispatch.
static Sphere make_sphere(BodyHot* h, ShapeInternal* s)
{
	return (Sphere){ add(h->position, rotate(h->rotation, s->local_pos)), s->sphere.radius };
}

static Capsule make_capsule(BodyHot* h, ShapeInternal* s)
{
	v3 lp = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 lq = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	return (Capsule){ add(h->position, rotate(h->rotation, lp)),
	                  add(h->position, rotate(h->rotation, lq)), s->capsule.radius };
}

static Box make_box(BodyHot* h, ShapeInternal* s)
{
	return (Box){ h->position, h->rotation, s->box.half_extents };
}

static ConvexHull make_convex_hull(BodyHot* h, ShapeInternal* s)
{
	return (ConvexHull){ s->hull.hull, h->position, h->rotation, s->hull.scale };
}

// Narrowphase dispatch for a single body pair.
static void narrowphase_pair(WorldInternal* w, int i, int j, InternalManifold** manifolds)
{
	// Canonical ordering: lower type first for upper-triangle dispatch,
	// lower body index first for same-type pairs (deterministic shape A).
	if (w->body_cold[i].shapes[0].type > w->body_cold[j].shapes[0].type || (w->body_cold[i].shapes[0].type == w->body_cold[j].shapes[0].type && i > j)) {
		int tmp = i; i = j; j = tmp;
	}

	ShapeInternal* s0 = &w->body_cold[i].shapes[0];
	ShapeInternal* s1 = &w->body_cold[j].shapes[0];
	BodyHot* h0 = &w->body_hot[i];
	BodyHot* h1 = &w->body_hot[j];
	InternalManifold im = { .body_a = i, .body_b = j };

	int hit = 0;

	if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_SPHERE)
		hit = collide_sphere_sphere(make_sphere(h0, s0), make_sphere(h1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CAPSULE)
		hit = collide_sphere_capsule(make_sphere(h0, s0), make_capsule(h1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_BOX)
		hit = collide_sphere_box(make_sphere(h0, s0), make_box(h1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CAPSULE)
		hit = collide_capsule_capsule(make_capsule(h0, s0), make_capsule(h1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_BOX)
		hit = collide_capsule_box(make_capsule(h0, s0), make_box(h1, s1), &im.m);
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX)
		hit = collide_box_box(make_box(h0, s0), make_box(h1, s1), &im.m);
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_HULL)
		hit = collide_hull_hull((ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents }, make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL)
		hit = collide_sphere_hull(make_sphere(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL)
		hit = collide_capsule_hull(make_capsule(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL)
		hit = collide_hull_hull(make_convex_hull(h0, s0), make_convex_hull(h1, s1), &im.m);

	if (hit) apush(*manifolds, im);
}

static void broadphase_n2(WorldInternal* w, InternalManifold** manifolds)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (asize(w->body_cold[i].shapes) == 0) continue;
		for (int j = i + 1; j < count; j++) {
			if (!split_alive(w->body_gen, j)) continue;
			if (asize(w->body_cold[j].shapes) == 0) continue;
			if (w->body_hot[i].inv_mass == 0.0f && w->body_hot[j].inv_mass == 0.0f) continue;
			// Skip sleeping-vs-sleeping pairs
			int isl_i = w->body_cold[i].island_id, isl_j = w->body_cold[j].island_id;
			if (isl_i >= 0 && isl_j >= 0 && (w->island_gen[isl_i] & 1) && (w->island_gen[isl_j] & 1) && !w->islands[isl_i].awake && !w->islands[isl_j].awake) continue;
			narrowphase_pair(w, i, j, manifolds);
		}
	}
}

// Build AABB lookup table indexed by leaf index from current node contents.
static AABB* bvh_build_lut(BVHTree* t)
{
	int lcount = asize(t->leaves);
	if (lcount == 0) return NULL;
	AABB* lut = CK_ALLOC(sizeof(AABB) * lcount);
	for (int i = 0; i < lcount; i++) {
		BVHLeaf* lf = &t->leaves[i];
		BVHChild* c = bvh_child(&t->nodes[lf->node_idx], lf->child_slot);
		lut[i] = bvh_child_aabb(c);
	}
	return lut;
}

static void broadphase_bvh(WorldInternal* w, InternalManifold** manifolds)
{
	bvh_refit(w->bvh_dynamic, w);

	// Incremental refinement: rebuild a subtree using binned SAH.
	AABB* lut = bvh_build_lut(w->bvh_dynamic);
	if (lut) { bvh_incremental_refine(w->bvh_dynamic, lut); CK_FREE(lut); }

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(w->bvh_dynamic, &pairs);
	bvh_cross_test(w->bvh_dynamic, w->bvh_static, &pairs);

	for (int i = 0; i < asize(pairs); i++) {
		int a = pairs[i].a, b = pairs[i].b;
		if (w->body_hot[a].inv_mass == 0.0f && w->body_hot[b].inv_mass == 0.0f) continue;
		int isl_a = w->body_cold[a].island_id, isl_b = w->body_cold[b].island_id;
		if (isl_a >= 0 && isl_b >= 0 && (w->island_gen[isl_a] & 1) && (w->island_gen[isl_b] & 1) && !w->islands[isl_a].awake && !w->islands[isl_b].awake) continue;
		narrowphase_pair(w, a, b, manifolds);
	}

	afree(pairs);
}

static void broadphase_and_collide(WorldInternal* w, InternalManifold** manifolds)
{
	if (w->broadphase_type == BROADPHASE_BVH) broadphase_bvh(w, manifolds);
	else broadphase_n2(w, manifolds);
}

// ---- src/inertia.c ----

// inertia.c -- mass and inertia tensor computation

// Multiply world-space inverse inertia tensor by a vector.
// I_world_inv * v = R * diag(inv_i) * R^T * v
static v3 inv_inertia_mul(quat rot, v3 inv_i, v3 v)
{
	v3 local = rotate(inv(rot), v);
	return rotate(rot, V3(local.x * inv_i.x, local.y * inv_i.y, local.z * inv_i.z));
}

// Compute diagonal inertia tensor for a shape (in local principal axes).
static v3 shape_inertia(ShapeInternal* s, float mass)
{
	if (mass <= 0.0f) return V3(0, 0, 0);

	switch (s->type) {
	case SHAPE_SPHERE: {
		float i = 0.4f * mass * s->sphere.radius * s->sphere.radius;
		return V3(i, i, i);
	}
	case SHAPE_BOX: {
		v3 h = s->box.half_extents;
		float x2 = 4.0f*h.x*h.x, y2 = 4.0f*h.y*h.y, z2 = 4.0f*h.z*h.z;
		return V3(mass*(y2+z2)/12.0f, mass*(x2+z2)/12.0f, mass*(x2+y2)/12.0f);
	}
	case SHAPE_CAPSULE: {
		float r = s->capsule.radius, hh = s->capsule.half_height;
		float r2 = r*r, hh2 = hh*hh;
		// Volume-weighted mass split: cylinder vs sphere (two hemispheres)
		float v_cyl = 2.0f * hh * r2; // pi cancels in ratio
		float v_sph = (4.0f/3.0f) * r2 * r;
		float v_tot = v_cyl + v_sph;
		float mc = mass * v_cyl / v_tot;
		float ms = mass * v_sph / v_tot;
		// Axial (Y): cylinder + sphere
		float iy = mc * r2 / 2.0f + ms * 2.0f * r2 / 5.0f;
		// Transverse (X,Z): cylinder + two hemispheres via parallel axis
		float ix_cyl = mc * (r2/4.0f + hh2/3.0f);
		float d = hh + 3.0f*r/8.0f; // hemisphere CoM offset from capsule center
		float ix_hemi = (83.0f/320.0f) * (ms/2.0f) * r2 + (ms/2.0f) * d * d;
		float ix = ix_cyl + 2.0f * ix_hemi;
		return V3(ix, iy, ix);
	}
	case SHAPE_HULL: {
		// Approximate via scaled AABB
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		float lo_x = 1e18f, hi_x = -1e18f;
		float lo_y = 1e18f, hi_y = -1e18f;
		float lo_z = 1e18f, hi_z = -1e18f;
		for (int i = 0; i < hull->vert_count; i++) {
			float x = hull->verts[i].x*sc.x, y = hull->verts[i].y*sc.y, z = hull->verts[i].z*sc.z;
			if (x < lo_x) lo_x = x; if (x > hi_x) hi_x = x;
			if (y < lo_y) lo_y = y; if (y > hi_y) hi_y = y;
			if (z < lo_z) lo_z = z; if (z > hi_z) hi_z = z;
		}
		float sx = hi_x-lo_x, sy = hi_y-lo_y, sz = hi_z-lo_z;
		return V3(mass*(sy*sy+sz*sz)/12.0f, mass*(sx*sx+sz*sz)/12.0f, mass*(sx*sx+sy*sy)/12.0f);
	}
	}
	return V3(0, 0, 0);
}

static v3 inertia_to_inv(v3 inertia) { return rcp(inertia); }

// Gyroscopic torque solver (single Newton-Raphson step).
// Corrects angular velocity for gyroscopic precession effects that explicit
// Euler integration misses. Without this, spinning bodies gain energy.
// Reference: Catto, GDC 2015 slide 76.
static v3 solve_gyroscopic(quat q, v3 inv_i, v3 omega, float h)
{
	v3 ib = rcp(inv_i);
	v3 wb = rotate(inv(q), omega);
	v3 iw = hmul(ib, wb);

	// Residual: f = h * cross(wb, Ib * wb)
	v3 f = scale(cross(wb, iw), h);

	// Jacobian: J = Ib + h * (skew(wb) * Ib - skew(Ib * wb))
	m3x3 Ib = diag(ib);
	m3x3 J = add(Ib, scale(sub(mul(skew(wb), Ib), skew(iw)), h));

	// Single Newton-Raphson update
	wb = sub(wb, solve(J, f));

	return rotate(q, wb);
}

// Volume of a shape (for mass distribution across compound bodies).
static float shape_volume(ShapeInternal* s)
{
	const float PI = 3.14159265f;
	switch (s->type) {
	case SHAPE_SPHERE: {
		float r = s->sphere.radius;
		return (4.0f/3.0f) * PI * r * r * r;
	}
	case SHAPE_CAPSULE: {
		float r = s->capsule.radius, h = s->capsule.half_height;
		return PI * r * r * (2.0f * h + (4.0f/3.0f) * r);
	}
	case SHAPE_BOX: {
		v3 e = s->box.half_extents;
		return 8.0f * e.x * e.y * e.z;
	}
	case SHAPE_HULL: {
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		float lo_x = 1e18f, hi_x = -1e18f;
		float lo_y = 1e18f, hi_y = -1e18f;
		float lo_z = 1e18f, hi_z = -1e18f;
		for (int i = 0; i < hull->vert_count; i++) {
			float x = hull->verts[i].x*sc.x, y = hull->verts[i].y*sc.y, z = hull->verts[i].z*sc.z;
			if (x < lo_x) lo_x = x; if (x > hi_x) hi_x = x;
			if (y < lo_y) lo_y = y; if (y > hi_y) hi_y = y;
			if (z < lo_z) lo_z = z; if (z > hi_z) hi_z = z;
		}
		return (hi_x-lo_x) * (hi_y-lo_y) * (hi_z-lo_z);
	}
	}
	return 0.0f;
}

// Recompute body inertia from all shapes. Mass is distributed by volume ratio,
// each shape's inertia is shifted to body origin via parallel axis theorem.
static void recompute_body_inertia(WorldInternal* w, int idx)
{
	float mass = w->body_cold[idx].mass;
	if (mass <= 0.0f) {
		w->body_hot[idx].inv_inertia_local = V3(0, 0, 0);
		return;
	}

	ShapeInternal* shapes = w->body_cold[idx].shapes;
	int n = asize(shapes);

	float total_vol = 0.0f;
	for (int i = 0; i < n; i++)
		total_vol += shape_volume(&shapes[i]);
	if (total_vol <= 0.0f) {
		w->body_hot[idx].inv_inertia_local = V3(0, 0, 0);
		return;
	}

	v3 total = V3(0, 0, 0);
	for (int i = 0; i < n; i++) {
		float sm = mass * shape_volume(&shapes[i]) / total_vol;
		v3 li = shape_inertia(&shapes[i], sm);

		// Parallel axis theorem: I_body += I_local + m*(|d|^2*I - d*d^T)
		// For diagonal tensor: I_x += m*(dy^2+dz^2), etc.
		v3 d = shapes[i].local_pos;
		total.x += li.x + sm * (d.y*d.y + d.z*d.z);
		total.y += li.y + sm * (d.x*d.x + d.z*d.z);
		total.z += li.z + sm * (d.x*d.x + d.y*d.y);
	}

	w->body_hot[idx].inv_inertia_local = inertia_to_inv(total);
}

// ---- src/solver.c ----

// solver.c -- contact constraint solver

static uint64_t body_pair_key(int a, int b)
{
	uint32_t lo = a < b ? a : b;
	uint32_t hi = a < b ? b : a;
	return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static void contact_tangent_basis(v3 n, v3* t1, v3* t2)
{
	if (fabsf(n.x) >= 0.57735f)
		*t1 = norm(V3(n.y, -n.x, 0.0f));
	else
		*t1 = norm(V3(0.0f, n.z, -n.y));
	*t2 = cross(n, *t1);
}

#define PATCH_MIN_AREA 0.001f

// Estimate contact patch area from manifold points projected onto contact plane.
static float estimate_patch_area(Contact* contacts, int count)
{
	if (count < 3) return PATCH_MIN_AREA;
	// Fan triangulation from contacts[0]
	float area = 0.0f;
	for (int i = 1; i < count - 1; i++)
		area += 0.5f * len(cross(sub(contacts[i].point, contacts[0].point), sub(contacts[i + 1].point, contacts[0].point)));
	return area > PATCH_MIN_AREA ? area : PATCH_MIN_AREA;
}

static float compute_effective_mass(BodyHot* a, BodyHot* b, float inv_mass_sum, v3 r_a, v3 r_b, v3 dir)
{
	v3 ra_x_d = cross(r_a, dir);
	v3 rb_x_d = cross(r_b, dir);
	float k = inv_mass_sum
		+ dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, ra_x_d), r_a), dir)
		+ dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, rb_x_d), r_b), dir);
	return k > 1e-12f ? 1.0f / k : 0.0f;
}

// Apply an impulse (linear + angular) to a body pair.
static void apply_impulse(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, v3 impulse)
{
	a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(impulse, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(r_b, impulse)));
}

// Match a new contact to a cached contact by feature ID. Returns index or -1.
static int warm_match(WarmManifold* wm, uint32_t feature_id)
{
	for (int i = 0; i < wm->count; i++)
		if (wm->contacts[i].feature_id == feature_id) return i;
	return -1;
}

static void solver_pre_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, float dt)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	float inv_dt = dt > 0.0f ? 1.0f / dt : 0.0f;

	for (int i = 0; i < manifold_count; i++) {
		InternalManifold* im = &manifolds[i];
		BodyHot* a = &w->body_hot[im->body_a];
		BodyHot* b = &w->body_hot[im->body_b];
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		if (inv_mass_sum == 0.0f) continue;

		float mu = sqrtf(a->friction * b->friction);
		float rest = a->restitution > b->restitution ? a->restitution : b->restitution;

		// Look up warm starting data
		uint64_t key = body_pair_key(im->body_a, im->body_b);
		WarmManifold* wm = map_get_ptr(w->warm_cache, key);

		SolverManifold smf = {
			.body_a = im->body_a,
			.body_b = im->body_b,
			.contact_start = asize(sc),
			.contact_count = im->m.count,
			.friction = mu,
		};

		int patch_mode = (w->friction_model == FRICTION_PATCH);

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			SolverContact s = {0};

			s.r_a = sub(ct->point, a->position);
			s.r_b = sub(ct->point, b->position);
			s.normal = ct->normal;
			s.penetration = ct->penetration;
			s.feature_id = ct->feature_id;

			s.eff_mass_n = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.normal);

			// Soft contact constraint (skip for hard SI — uses NGS position correction instead)
			if (w->solver_type != SOLVER_SI) {
				float hertz = w->contact_hertz;
				if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
				float omega = 2.0f * 3.14159265f * hertz;
				float d = 2.0f * w->contact_damping_ratio * omega;
				float k = omega * omega;
				float hd = dt * d, hhk = dt * dt * k;
				float denom = hd + hhk;
				if (denom > 1e-12f) {
					s.softness = 1.0f / denom;
					float bias_rate = dt * k * s.softness;
					float K = s.eff_mass_n > 0.0f ? 1.0f / s.eff_mass_n : 0.0f;
					s.eff_mass_n = (K + s.softness) > 1e-12f ? 1.0f / (K + s.softness) : 0.0f;
					float pen = ct->penetration - SOLVER_SLOP;
					s.bias = pen > 0.0f ? -bias_rate * pen : 0.0f;
					if (s.bias < -w->max_push_velocity) s.bias = -w->max_push_velocity;
				}
			}

			if (!patch_mode) {
				contact_tangent_basis(ct->normal, &s.tangent1, &s.tangent2);
				s.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent1);
				s.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent2);
			}

			v3 vel_a = add(a->velocity, cross(a->angular_velocity, s.r_a));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, s.r_b));
			float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
			s.bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;

			// When restitution is active, the bounce velocity handles separation.
			// Disable the penetration bias to prevent energy injection (bias + bounce
			// are additive in the solver, so both active injects extra energy).
			if (s.bounce != 0.0f) s.bias = 0.0f;

			apush(sc, s);
		}

		// Patch friction: compute centroid, patch area, tangent basis at manifold level.
		if (patch_mode) {
			v3 centroid_a = V3(0, 0, 0), centroid_b = V3(0, 0, 0);
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				centroid_a = add(centroid_a, s->r_a);
				centroid_b = add(centroid_b, s->r_b);
			}
			float inv_n = 1.0f / (float)smf.contact_count;
			smf.centroid_r_a = scale(centroid_a, inv_n);
			smf.centroid_r_b = scale(centroid_b, inv_n);

			// Use first contact normal for tangent basis (all share same normal in a manifold)
			v3 n = sc[smf.contact_start].normal;
			smf.normal = n;
			contact_tangent_basis(n, &smf.tangent1, &smf.tangent2);
			smf.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent1);
			smf.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent2);
			smf.patch_area = estimate_patch_area(im->m.contacts, im->m.count);
			smf.patch_radius = 0.6667f * sqrtf(smf.patch_area * (1.0f / 3.14159265f));

			// Torsional friction effective mass: 1 / (n^T * I_a_inv * n + n^T * I_b_inv * n)
			float k_twist = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, n), n) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, n), n);
			smf.eff_mass_twist = k_twist > 1e-12f ? 1.0f / k_twist : 0.0f;
		}

		// Warm start: match new contacts to cached contacts.
		// Pass 1: exact feature ID match.
		// Pass 2: spatial fallback -- closest unmatched old contact by r_a distance.
		// Pass 3: redistribute any remaining unmatched old impulse evenly.
		if (wm) {
			int old_matched[MAX_CONTACTS] = {0};
			int new_matched[MAX_CONTACTS] = {0};
			// Pass 1: feature ID
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				if (s->feature_id == 0) continue;
				int match = warm_match(wm, s->feature_id);
				if (match >= 0) {
					s->lambda_n = wm->contacts[match].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[match].lambda_t1;
						s->lambda_t2 = wm->contacts[match].lambda_t2;
					}
					old_matched[match] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 2: spatial fallback for unmatched contacts
			float spatial_tol2 = 0.01f;
			for (int c = 0; c < smf.contact_count; c++) {
				if (new_matched[c]) continue;
				SolverContact* s = &sc[smf.contact_start + c];
				float best_d2 = spatial_tol2;
				int best = -1;
				for (int j = 0; j < wm->count; j++) {
					if (old_matched[j]) continue;
					float d2 = len2(sub(s->r_a, wm->contacts[j].r_a));
					if (d2 < best_d2) { best_d2 = d2; best = j; }
				}
				if (best >= 0) {
					s->lambda_n = wm->contacts[best].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[best].lambda_t1;
						s->lambda_t2 = wm->contacts[best].lambda_t2;
					}
					old_matched[best] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 3: redistribute remaining unmatched old impulse
			int new_unmatched = 0;
			for (int c = 0; c < smf.contact_count; c++)
				if (!new_matched[c]) new_unmatched++;
			if (new_unmatched > 0) {
				float leftover_n = 0, leftover_t1 = 0, leftover_t2 = 0;
				for (int j = 0; j < wm->count; j++) {
					if (!old_matched[j]) {
						leftover_n += wm->contacts[j].lambda_n;
						if (!patch_mode) {
							leftover_t1 += wm->contacts[j].lambda_t1;
							leftover_t2 += wm->contacts[j].lambda_t2;
						}
					}
				}
				float share = 1.0f / (float)new_unmatched;
				for (int c = 0; c < smf.contact_count; c++) {
					if (new_matched[c]) continue;
					SolverContact* s = &sc[smf.contact_start + c];
					s->lambda_n = leftover_n * share;
					if (!patch_mode) {
						s->lambda_t1 = leftover_t1 * share;
						s->lambda_t2 = leftover_t2 * share;
					}
				}
			}
			// Warm start manifold-level friction
			if (patch_mode) {
				smf.lambda_t1 = wm->manifold_lambda_t1;
				smf.lambda_t2 = wm->manifold_lambda_t2;
				smf.lambda_twist = wm->manifold_lambda_twist;
			}
		}

		// Speculative contacts (negative penetration, kept alive by margin) must
		// not carry warm-started impulse -- they exist for cache continuity only.
		for (int c = 0; c < smf.contact_count; c++) {
			SolverContact* s = &sc[smf.contact_start + c];
			if (s->penetration < 0.0f) { s->lambda_n = 0.0f; s->lambda_t1 = 0.0f; s->lambda_t2 = 0.0f; }
		}

		// Block solver: compute NxN normal coupling matrix A[i][j]
		if (w->solver_type == SOLVER_BLOCK && smf.contact_count >= 2) {
			int n = smf.contact_count;
			for (int ci = 0; ci < n; ci++) {
				SolverContact* si = &sc[smf.contact_start + ci];
				for (int cj = ci; cj < n; cj++) {
					SolverContact* sj = &sc[smf.contact_start + cj];
					// A[i][j] = inv_mass_sum
					//   + dot(cross(I_a^-1 * cross(r_a_i, n), r_a_j), n)
					//   + dot(cross(I_b^-1 * cross(r_b_i, n), r_b_j), n)
					v3 nn = si->normal;
					float kk = inv_mass_sum;
					kk += dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(si->r_a, nn)), sj->r_a), nn);
					kk += dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(si->r_b, nn)), sj->r_b), nn);
					if (ci == cj) kk += si->softness; // regularization on diagonal
					smf.K_nn[ci * 4 + cj] = kk;
					smf.K_nn[cj * 4 + ci] = kk; // symmetric
				}
			}
			// Invert: small matrix (2x2, 3x3, or 4x4)
			if (n == 2) {
				float a00 = smf.K_nn[0], a01 = smf.K_nn[1];
				float a10 = smf.K_nn[4], a11 = smf.K_nn[5];
				float det = a00 * a11 - a01 * a10;
				if (det != 0.0f) {
					float inv_det = 1.0f / det;
					smf.K_nn_inv[0] = a11 * inv_det;  smf.K_nn_inv[1] = -a01 * inv_det;
					smf.K_nn_inv[4] = -a10 * inv_det;  smf.K_nn_inv[5] = a00 * inv_det;
				}
			} else {
				// For 3x3 and 4x4: compute inverse via cofactor expansion
				// For now, fall back to per-contact PGS for count > 2
				// (3+ contact manifolds are rare in practice)
			}
		}

		apush(sm, smf);
	}

	// Apply warm start impulses
	int patch_warm = (w->friction_model == FRICTION_PATCH);
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (patch_warm) {
				if (s->lambda_n == 0.0f) continue;
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n));
			} else {
				if (s->lambda_n == 0.0f && s->lambda_t1 == 0.0f && s->lambda_t2 == 0.0f)
					continue;
				v3 P = add(add(scale(s->normal, s->lambda_n), scale(s->tangent1, s->lambda_t1)), scale(s->tangent2, s->lambda_t2));
				apply_impulse(a, b, s->r_a, s->r_b, P);
			}
		}
		// Warm start manifold-level friction
		if (patch_warm && (m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f)) {
			v3 P = add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, P);
		}
		// Warm start torsional friction (pure angular impulse along normal)
		if (patch_warm && m->lambda_twist != 0.0f) {
			v3 twist_impulse = scale(m->normal, m->lambda_twist);
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, twist_impulse));
		}
	}

	*out_sm = sm;
	*out_sc = sc;
}

static void solver_iterate(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
{
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			// Relative velocity at contact point
			v3 dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));

			// --- Normal constraint ---
			float vn = dot(dv, s->normal);
			float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce));
			float old_n = s->lambda_n;
			s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
			float delta_n = s->lambda_n - old_n;
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, delta_n));

			// --- Friction tangent 1 ---
			dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt1 = dot(dv, s->tangent1);
			float lambda_t1 = s->eff_mass_t1 * (-vt1);
			float max_f = m->friction * s->lambda_n;
			float old_t1 = s->lambda_t1;
			s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + lambda_t1, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

			// --- Friction tangent 2 ---
			dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt2 = dot(dv, s->tangent2);
			float lambda_t2 = s->eff_mass_t2 * (-vt2);
			float old_t2 = s->lambda_t2;
			s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + lambda_t2, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
		}
	}
}

// NGS position correction: directly fix remaining penetration after velocity solve
// and position integration. Operates on positions only, no velocity modification.
static void solver_position_correct(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
{
	for (int iter = 0; iter < w->position_iters; iter++) {
		for (int i = 0; i < sm_count; i++) {
			SolverManifold* m = &sm[i];
			BodyHot* a = &w->body_hot[m->body_a];
			BodyHot* b = &w->body_hot[m->body_b];
			float inv_mass_sum = a->inv_mass + b->inv_mass;

			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];

				// Recompute separation from current positions
				v3 r_a = sub(add(a->position, rotate(a->rotation, rotate(inv(a->rotation), s->r_a))), a->position);
				v3 r_b = sub(add(b->position, rotate(b->rotation, rotate(inv(b->rotation), s->r_b))), b->position);
				v3 p_a = add(a->position, r_a);
				v3 p_b = add(b->position, r_b);
				float separation = dot(sub(p_b, p_a), s->normal) - s->penetration;

				float C = fminf(0.0f, separation + SOLVER_SLOP);
				if (C >= 0.0f) continue;

				float correction = -SOLVER_POS_BAUMGARTE * C;
				if (correction > SOLVER_POS_MAX_CORRECTION)
					correction = SOLVER_POS_MAX_CORRECTION;

				// Effective mass for position correction (linear only for speed)
				float k = inv_mass_sum;
				float delta = correction / k;

				v3 P = scale(s->normal, delta);
				a->position = sub(a->position, scale(P, a->inv_mass));
				b->position = add(b->position, scale(P, b->inv_mass));
			}
		}
	}
}

// Relax contacts: recompute separation and bias from current body positions.
// Called after each substep's integrate_positions for SOFT_STEP and BLOCK solvers.
// Keeps normals, tangent basis, effective mass unchanged — only refreshes the RHS.
static void solver_relax_contacts(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, float dt)
{
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		float hertz = w->contact_hertz;
		if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
		float omega = 2.0f * 3.14159265f * hertz;
		float d = 2.0f * w->contact_damping_ratio * omega;
		float k = omega * omega;
		float hd = dt * d, hhk = dt * dt * k;
		float denom = hd + hhk;
		float bias_rate = (denom > 1e-12f) ? dt * k / denom : 0.0f;

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			v3 p_a = add(a->position, s->r_a);
			v3 p_b = add(b->position, s->r_b);
			float separation = dot(sub(p_b, p_a), s->normal) - s->penetration;

			float pen = -separation - SOLVER_SLOP;
			s->bias = pen > 0.0f ? -bias_rate * pen : 0.0f;
			if (s->bias < -w->max_push_velocity) s->bias = -w->max_push_velocity;

			if (s->bounce != 0.0f) s->bias = 0.0f;
		}
	}
}

// Persistent warm cache: update active pairs, age stale pairs, evict old stale.
// BEPU-style freshness: entries survive one extra frame so warm data isn't lost
// when narrowphase misses a pair for a single frame (FP noise).
static void solver_post_solve(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, InternalManifold* manifolds, int manifold_count)
{
	// Update active pairs (per sub-step)
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		uint64_t key = body_pair_key(m->body_a, m->body_b);

		WarmManifold wm = {0};
		wm.count = m->contact_count;
		wm.stale = 0;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			wm.contacts[ci] = (WarmContact){
				.feature_id = s->feature_id,
				.r_a = s->r_a,
				.lambda_n = s->lambda_n,
				.lambda_t1 = s->lambda_t1,
				.lambda_t2 = s->lambda_t2,
			};
		}
		if (w->friction_model == FRICTION_PATCH) {
			wm.manifold_lambda_t1 = m->lambda_t1;
			wm.manifold_lambda_t2 = m->lambda_t2;
			wm.manifold_lambda_twist = m->lambda_twist;
		}

		map_set(w->warm_cache, key, wm);
	}

	afree(sm);
	afree(sc);
}

// Age and evict stale warm cache entries. Called once per frame (not per sub-step).
static void warm_cache_age_and_evict(WorldInternal* w)
{
	for (int i = 0; i < map_size(w->warm_cache); i++)
		w->warm_cache[i].stale++;
	int i = 0;
	while (i < map_size(w->warm_cache)) {
		if (w->warm_cache[i].stale > 1)
			map_del(w->warm_cache, map_keys(w->warm_cache)[i]);
		else
			i++;
	}
}

// -----------------------------------------------------------------------------
// Graph coloring: assign colors to constraints so no two in same color share a body.
// Uses uint64_t bitmask per body (max 64 colors).

static void color_constraints(ConstraintRef* refs, int count, int body_count, int* out_batch_starts, int* out_color_count)
{
	uint64_t* body_colors = (uint64_t*)CK_ALLOC(body_count * sizeof(uint64_t));
	memset(body_colors, 0, body_count * sizeof(uint64_t));
	int max_color = 0;

	for (int i = 0; i < count; i++) {
		uint64_t used = body_colors[refs[i].body_a] | body_colors[refs[i].body_b];
		uint64_t avail = ~used;
		int color = 0;
		if (avail) {
			// Find lowest set bit (MSVC: _BitScanForward64, GCC: __builtin_ctzll)
#ifdef _MSC_VER
			unsigned long idx;
			_BitScanForward64(&idx, avail);
			color = (int)idx;
#else
			color = __builtin_ctzll(avail);
#endif
		}
		assert(color < 64);
		refs[i].color = (uint8_t)color;
		uint64_t bit = 1ULL << color;
		body_colors[refs[i].body_a] |= bit;
		body_colors[refs[i].body_b] |= bit;
		if (color > max_color) max_color = color;
	}

	// Counting sort by color
	int color_count = max_color + 1;
	int counts[64] = {0};
	for (int i = 0; i < count; i++) counts[refs[i].color]++;

	int offsets[64];
	offsets[0] = 0;
	for (int c = 1; c < color_count; c++) offsets[c] = offsets[c-1] + counts[c-1];

	// Record batch starts before sorting
	for (int c = 0; c < color_count; c++) out_batch_starts[c] = offsets[c];
	out_batch_starts[color_count] = count; // sentinel

	ConstraintRef* sorted = (ConstraintRef*)CK_ALLOC(count * sizeof(ConstraintRef));
	for (int i = 0; i < count; i++)
		sorted[offsets[refs[i].color]++] = refs[i];
	memcpy(refs, sorted, count * sizeof(ConstraintRef));

	CK_FREE(sorted);
	CK_FREE(body_colors);
	*out_color_count = color_count;
}

// Block solver: Murty total enumeration for 2 contacts.
// Solves the LCP: vn = A*x + b', vn >= 0, x >= 0, vn_i * x_i = 0
// Enumerates all 4 cases and takes the first valid solution.
static void solve_block_2(BodyHot* a, BodyHot* b, SolverManifold* m, SolverContact* sc)
{
	SolverContact* c0 = &sc[m->contact_start];
	SolverContact* c1 = &sc[m->contact_start + 1];

	float a0 = c0->lambda_n, a1 = c1->lambda_n; // accumulated impulse (old)

	// Compute relative normal velocities
	v3 dv0 = sub(add(b->velocity, cross(b->angular_velocity, c0->r_b)), add(a->velocity, cross(a->angular_velocity, c0->r_a)));
	v3 dv1 = sub(add(b->velocity, cross(b->angular_velocity, c1->r_b)), add(a->velocity, cross(a->angular_velocity, c1->r_a)));
	float vn0 = dot(dv0, c0->normal);
	float vn1 = dot(dv1, c1->normal);

	// b' = b - A * a, where b = vn + bias + bounce
	float b0 = vn0 + c0->bias + c0->bounce;
	float b1 = vn1 + c1->bias + c1->bounce;
	b0 -= m->K_nn[0] * a0 + m->K_nn[1] * a1;
	b1 -= m->K_nn[4] * a0 + m->K_nn[5] * a1;

	// Also include softness feedback in b'
	b0 -= c0->softness * a0;
	b1 -= c1->softness * a1;

	float x0, x1;

	// Case 1: both active (vn0 = 0, vn1 = 0) → x = -A_inv * b'
	x0 = -(m->K_nn_inv[0] * b0 + m->K_nn_inv[1] * b1);
	x1 = -(m->K_nn_inv[4] * b0 + m->K_nn_inv[5] * b1);
	if (x0 >= 0.0f && x1 >= 0.0f) goto apply;

	// Case 2: c0 active, c1 inactive (vn0 = 0, x1 = 0)
	x0 = -b0 / m->K_nn[0];
	x1 = 0.0f;
	if (x0 >= 0.0f) {
		float vn1_check = m->K_nn[4] * x0 + b1;
		if (vn1_check >= 0.0f) goto apply;
	}

	// Case 3: c0 inactive, c1 active (x0 = 0, vn1 = 0)
	x0 = 0.0f;
	x1 = -b1 / m->K_nn[5];
	if (x1 >= 0.0f) {
		float vn0_check = m->K_nn[1] * x1 + b0;
		if (vn0_check >= 0.0f) goto apply;
	}

	// Case 4: both inactive (x0 = 0, x1 = 0)
	x0 = 0.0f;
	x1 = 0.0f;

apply:;
	float d0 = x0 - a0, d1 = x1 - a1;
	c0->lambda_n = x0;
	c1->lambda_n = x1;
	v3 P0 = scale(c0->normal, d0);
	v3 P1 = scale(c1->normal, d1);
	v3 P = add(P0, P1);
	a->velocity = sub(a->velocity, scale(P, a->inv_mass));
	b->velocity = add(b->velocity, scale(P, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity,
		add(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(c0->r_a, P0)),
		    inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(c1->r_a, P1))));
	b->angular_velocity = add(b->angular_velocity,
		add(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(c0->r_b, P0)),
		    inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(c1->r_b, P1))));
}

// Forward declarations for joint solvers (defined in joints.c, included after solver.c).
static void solve_ball_socket(WorldInternal* w, SolverBallSocket* s);
static void solve_distance(WorldInternal* w, SolverDistance* s);

// Dispatch a single constraint solve by type.
static void solve_constraint(WorldInternal* w, ConstraintRef* ref, SolverManifold* sm, SolverContact* sc, SolverBallSocket* bs, SolverDistance* dist)
{
	switch (ref->type) {
	case CTYPE_CONTACT: {
		SolverManifold* m = &sm[ref->index];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		if (w->friction_model == FRICTION_PATCH) {
			// Solve normal constraints (block LCP or sequential)
			float total_lambda_n = 0.0f;
			if (w->solver_type == SOLVER_BLOCK && m->contact_count == 2) {
				solve_block_2(a, b, m, sc);
				for (int ci = 0; ci < m->contact_count; ci++)
					total_lambda_n += sc[m->contact_start + ci].lambda_n;
			} else {
				for (int ci = 0; ci < m->contact_count; ci++) {
					SolverContact* s = &sc[m->contact_start + ci];
					v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
					float vn = dot(dv, s->normal);
					float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
					float old_n = s->lambda_n;
					s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
					apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n - old_n));
					total_lambda_n += s->lambda_n;
				}
			}

			// Manifold-level 2D friction at centroid, clamped by aggregate normal force
			float max_f = m->friction * total_lambda_n;
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, m->centroid_r_b)), add(a->velocity, cross(a->angular_velocity, m->centroid_r_a)));
			float vt1 = dot(dv, m->tangent1);
			float old_t1 = m->lambda_t1;
			m->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + m->eff_mass_t1 * (-vt1), max_f));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, scale(m->tangent1, m->lambda_t1 - old_t1));

			dv = sub(add(b->velocity, cross(b->angular_velocity, m->centroid_r_b)), add(a->velocity, cross(a->angular_velocity, m->centroid_r_a)));
			float vt2 = dot(dv, m->tangent2);
			float old_t2 = m->lambda_t2;
			m->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + m->eff_mass_t2 * (-vt2), max_f));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, scale(m->tangent2, m->lambda_t2 - old_t2));

			// Torsional friction: resist spin around contact normal
			float max_twist = m->friction * total_lambda_n * m->patch_radius;
			float w_rel = dot(sub(b->angular_velocity, a->angular_velocity), m->normal);
			float lambda_tw = m->eff_mass_twist * (-w_rel);
			float old_tw = m->lambda_twist;
			m->lambda_twist = fmaxf(-max_twist, fminf(old_tw + lambda_tw, max_twist));
			float delta_tw = m->lambda_twist - old_tw;
			v3 twist_impulse = scale(m->normal, delta_tw);
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, twist_impulse));
		} else {
			// Per-point Coulomb friction
			// Block solver handles normals first, then friction per-contact
			if (w->solver_type == SOLVER_BLOCK && m->contact_count == 2) {
				solve_block_2(a, b, m, sc);
			} else {
				for (int ci = 0; ci < m->contact_count; ci++) {
					SolverContact* s = &sc[m->contact_start + ci];
					v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
					float vn = dot(dv, s->normal);
					float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
					float old_n = s->lambda_n;
					s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
					apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n - old_n));
				}
			}
			// Per-contact friction (uses normal impulse from above)
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float max_f = m->friction * s->lambda_n;
				float vt1 = dot(dv, s->tangent1);
				float old_t1 = s->lambda_t1;
				s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + s->eff_mass_t1*(-vt1), max_f));
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

				dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float vt2 = dot(dv, s->tangent2);
				float old_t2 = s->lambda_t2;
				s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + s->eff_mass_t2*(-vt2), max_f));
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
			}
		}
		break;
	}
	case CTYPE_BALL_SOCKET: solve_ball_socket(w, &bs[ref->index]); break;
	case CTYPE_DISTANCE:    solve_distance(w, &dist[ref->index]); break;
	}
}

// ---- src/joints.c ----

// joints.c -- joint constraint solvers (ball socket, distance)

// Symmetric 3x3 stored as 6 floats: [xx, xy, xz, yy, yz, zz].
static void sym3x3_inverse(const float* in, float* out)
{
	float a = in[0], b = in[1], c = in[2];
	float d = in[3], e = in[4], f = in[5];
	float det = a*(d*f - e*e) - b*(b*f - c*e) + c*(b*e - c*d);
	float inv_det = det != 0.0f ? 1.0f / det : 0.0f;
	out[0] = (d*f - e*e) * inv_det;
	out[1] = (c*e - b*f) * inv_det;
	out[2] = (b*e - c*d) * inv_det;
	out[3] = (a*f - c*c) * inv_det;
	out[4] = (b*c - a*e) * inv_det;
	out[5] = (a*d - b*b) * inv_det;
}

static v3 sym3x3_mul_v3(const float* m, v3 v)
{
	return V3(
		m[0]*v.x + m[1]*v.y + m[2]*v.z,
		m[1]*v.x + m[3]*v.y + m[4]*v.z,
		m[2]*v.x + m[4]*v.y + m[5]*v.z);
}

// Convert spring params to solver coefficients (bepu approach).
static void spring_compute(SpringParams sp, float dt, float* pos_to_vel, float* softness)
{
	if (sp.frequency <= 0.0f) {
		// Rigid constraint: Baumgarte stabilization with fractional correction.
		*pos_to_vel = dt > 0.0f ? SOLVER_BAUMGARTE / dt : 0.0f;
		*softness = 0.0f;
		return;
	}
	float omega = 2.0f * 3.14159265f * sp.frequency;
	float d = 2.0f * sp.damping_ratio * omega;
	float k = omega * omega;
	float hd = dt * d, hk = dt * k, hhk = dt * hk;
	float denom = hd + hhk;
	if (denom < 1e-12f) { *pos_to_vel = 0; *softness = 0; return; }
	*softness = 1.0f / denom;
	*pos_to_vel = hk * (*softness);
}

// Build ball socket effective mass (symmetric 3x3).
// K = (inv_mass_a + inv_mass_b)*I + skew(r_a)^T * I_a^-1 * skew(r_a) + same for B.
static void ball_socket_eff_mass(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, float softness, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	// Skew-sandwich: skew(r)^T * I^-1 * skew(r) for diagonal inertia.
	// For body A: I_a^-1 is diagonal in local space, r_a is world space.
	// Transform r_a to local, compute product, transform back.
	// Or equivalently use the explicit formula.
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m }; // start with linear term

	// Add angular contribution for body A
	v3 ia = a->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		// columns of world-space inverse inertia: R * diag(ia) * R^T
		// skew(r)^T * I_world^-1 * skew(r) computed column by column
		v3 e0 = inv_inertia_mul(a->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(a->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(a->rotation, ia, V3(-r_a.y, r_a.x, 0));
		// skew(r)^T * [e0 e1 e2] -- but skew(r)^T has rows [0 -rz ry; rz 0 -rx; -ry rx 0]
		// Result[i][j] = dot(skew_row_i(r), e_j)
		// Row 0 of skew^T: (0, -rz, ry)
		K[0] += -r_a.z*e0.y + r_a.y*e0.z; // xx
		K[1] += -r_a.z*e1.y + r_a.y*e1.z; // xy
		K[2] += -r_a.z*e2.y + r_a.y*e2.z; // xz
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z; // yy
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z; // yz
		K[5] += -r_a.y*e2.x + r_a.x*e2.y; // zz
	}

	v3 ib = b->inv_inertia_local;
	if (ib.x > 0 || ib.y > 0 || ib.z > 0) {
		v3 e0 = inv_inertia_mul(b->rotation, ib, V3(0, -r_b.z, r_b.y));
		v3 e1 = inv_inertia_mul(b->rotation, ib, V3(r_b.z, 0, -r_b.x));
		v3 e2 = inv_inertia_mul(b->rotation, ib, V3(-r_b.y, r_b.x, 0));
		K[0] += -r_b.z*e0.y + r_b.y*e0.z;
		K[1] += -r_b.z*e1.y + r_b.y*e1.z;
		K[2] += -r_b.z*e2.y + r_b.y*e2.z;
		K[3] +=  r_b.z*e1.x - r_b.x*e1.z;
		K[4] +=  r_b.z*e2.x - r_b.x*e2.z;
		K[5] += -r_b.y*e2.x + r_b.x*e2.y;
	}

	// Add softness to diagonal: K_eff = K + gamma * I
	K[0] += softness; K[3] += softness; K[5] += softness;
	sym3x3_inverse(K, out);
}

// Pre-solve joints: build solver arrays from persistent JointInternal data.
static void joints_pre_solve(WorldInternal* w, float dt, SolverBallSocket** out_bs, SolverDistance** out_dist)
{
	CK_DYNA SolverBallSocket* bs = NULL;
	CK_DYNA SolverDistance* dist = NULL;
	int joint_count = asize(w->joints);

	for (int i = 0; i < joint_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			SolverBallSocket s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->ball_socket.local_a);
			s.r_b = rotate(b->rotation, j->ball_socket.local_b);

			float ptv, soft;
			spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
			s.softness = soft;
			ball_socket_eff_mass(a, b, s.r_a, s.r_b, soft, s.eff_mass);

			// Position error: world anchor B - world anchor A
			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			s.bias = scale(sub(anchor_b, anchor_a), ptv);

			// Warm start from persistent storage
			s.lambda = j->warm_lambda3;

			apush(bs, s);
		} else if (j->type == JOINT_DISTANCE) {
			SolverDistance s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->distance.local_a);
			s.r_b = rotate(b->rotation, j->distance.local_b);

			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 delta = sub(anchor_b, anchor_a);
			float dist_val = len(delta);
			s.axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

			float ptv, soft;
			spring_compute(j->distance.spring, dt, &ptv, &soft);
			s.softness = soft;

			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s.r_a, s.axis)), s.r_a), s.axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s.r_b, s.axis)), s.r_b), s.axis);
			k += soft;
			s.eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;

			float error = dist_val - j->distance.rest_length;
			s.bias = -ptv * error;

			s.lambda = j->warm_lambda1;

			apush(dist, s);
		}
	}

	*out_bs = bs;
	*out_dist = dist;
}

// Apply warm start impulses for joints.
static void joints_warm_start(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &bs[i];
		if (s->lambda.x == 0 && s->lambda.y == 0 && s->lambda.z == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, s->lambda);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &dist[i];
		if (s->lambda == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, scale(s->axis, s->lambda));
	}
}

static void solve_ball_socket(WorldInternal* w, SolverBallSocket* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	v3 rhs = sub(neg(add(dv, s->bias)), scale(s->lambda, s->softness));
	v3 impulse = sym3x3_mul_v3(s->eff_mass, rhs);
	s->lambda = add(s->lambda, impulse);
	apply_impulse(a, b, s->r_a, s->r_b, impulse);
}

static void solve_distance(WorldInternal* w, SolverDistance* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	float cdot = dot(dv, s->axis);
	float lambda = s->eff_mass * (-cdot + s->bias - s->softness * s->lambda);
	s->lambda += lambda;
	apply_impulse(a, b, s->r_a, s->r_b, scale(s->axis, lambda));
}

// Store joint accumulated impulses back to persistent storage.
static void joints_post_solve(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++)
		w->joints[bs[i].joint_idx].warm_lambda3 = bs[i].lambda;
	for (int i = 0; i < dist_count; i++)
		w->joints[dist[i].joint_idx].warm_lambda1 = dist[i].lambda;
	afree(bs);
	afree(dist);
}

// ---- src/islands.c ----

// islands.c -- island connectivity, sleep/wake management

static int island_create(WorldInternal* w)
{
	int idx;
	if (asize(w->island_free) > 0) {
		idx = apop(w->island_free);
		w->island_gen[idx]++;
	} else {
		idx = asize(w->islands);
		Island zero = {0};
		apush(w->islands, zero);
		apush(w->island_gen, 0);
	}
	w->islands[idx] = (Island){
		.head_body = -1, .tail_body = -1, .body_count = 0,
		.head_joint = -1, .tail_joint = -1, .joint_count = 0,
		.constraint_remove_count = 0, .awake = 1,
	};
	w->island_gen[idx] |= 1; // odd = alive
	return idx;
}

static void island_destroy(WorldInternal* w, int id)
{
	w->island_gen[id]++;
	if (w->island_gen[id] & 1) w->island_gen[id]++; // ensure even = dead
	apush(w->island_free, id);
}

static int island_alive(WorldInternal* w, int id)
{
	return id >= 0 && id < asize(w->islands) && (w->island_gen[id] & 1);
}

static void island_add_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	c->island_id = island_id;
	c->island_prev = isl->tail_body;
	c->island_next = -1;
	if (isl->tail_body >= 0)
		w->body_cold[isl->tail_body].island_next = body_idx;
	else
		isl->head_body = body_idx;
	isl->tail_body = body_idx;
	isl->body_count++;
}

static void island_remove_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	if (c->island_prev >= 0)
		w->body_cold[c->island_prev].island_next = c->island_next;
	else
		isl->head_body = c->island_next;
	if (c->island_next >= 0)
		w->body_cold[c->island_next].island_prev = c->island_prev;
	else
		isl->tail_body = c->island_prev;
	c->island_id = -1;
	c->island_prev = -1;
	c->island_next = -1;
	isl->body_count--;
}

static void island_add_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	j->island_id = island_id;
	j->island_prev = isl->tail_joint;
	j->island_next = -1;
	if (isl->tail_joint >= 0)
		w->joints[isl->tail_joint].island_next = joint_idx;
	else
		isl->head_joint = joint_idx;
	isl->tail_joint = joint_idx;
	isl->joint_count++;
}

static void island_remove_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_prev >= 0)
		w->joints[j->island_prev].island_next = j->island_next;
	else
		isl->head_joint = j->island_next;
	if (j->island_next >= 0)
		w->joints[j->island_next].island_prev = j->island_prev;
	else
		isl->tail_joint = j->island_prev;
	j->island_id = -1;
	j->island_prev = -1;
	j->island_next = -1;
	isl->joint_count--;
}

// Merge two islands. Keep the bigger one, walk the smaller one remapping bodies/joints.
static int island_merge(WorldInternal* w, int id_a, int id_b)
{
	if (id_a == id_b) return id_a;
	Island* a = &w->islands[id_a];
	Island* b = &w->islands[id_b];
	// Keep the bigger island
	if (a->body_count + a->joint_count < b->body_count + b->joint_count) {
		int tmp = id_a; id_a = id_b; id_b = tmp;
		a = &w->islands[id_a]; b = &w->islands[id_b];
	}
	// Remap smaller's bodies to bigger
	int bi = b->head_body;
	while (bi >= 0) {
		int next = w->body_cold[bi].island_next;
		w->body_cold[bi].island_id = id_a;
		bi = next;
	}
	// Splice body lists: append b's list to a's tail
	if (b->head_body >= 0) {
		if (a->tail_body >= 0) {
			w->body_cold[a->tail_body].island_next = b->head_body;
			w->body_cold[b->head_body].island_prev = a->tail_body;
		} else {
			a->head_body = b->head_body;
		}
		a->tail_body = b->tail_body;
		a->body_count += b->body_count;
	}
	// Remap smaller's joints
	int ji = b->head_joint;
	while (ji >= 0) {
		int next = w->joints[ji].island_next;
		w->joints[ji].island_id = id_a;
		ji = next;
	}
	// Splice joint lists
	if (b->head_joint >= 0) {
		if (a->tail_joint >= 0) {
			w->joints[a->tail_joint].island_next = b->head_joint;
			w->joints[b->head_joint].island_prev = a->tail_joint;
		} else {
			a->head_joint = b->head_joint;
		}
		a->tail_joint = b->tail_joint;
		a->joint_count += b->joint_count;
	}
	a->constraint_remove_count += b->constraint_remove_count;
	// If either was awake, result is awake
	if (b->awake) a->awake = 1;
	island_destroy(w, id_b);
	return id_a;
}

static void island_wake(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 1;
	int bi = isl->head_body;
	while (bi >= 0) {
		w->body_hot[bi].sleep_time = 0.0f;
		bi = w->body_cold[bi].island_next;
	}
}

static void island_sleep(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 0;
	int bi = isl->head_body;
	while (bi >= 0) {
		w->body_hot[bi].velocity = V3(0, 0, 0);
		w->body_hot[bi].angular_velocity = V3(0, 0, 0);
		bi = w->body_cold[bi].island_next;
	}
}

// Link a newly created joint to islands: merge/create as needed, wake if sleeping.
static void link_joint_to_islands(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	int ba = j->body_a, bb = j->body_b;
	int is_static_a = w->body_hot[ba].inv_mass == 0.0f;
	int is_static_b = w->body_hot[bb].inv_mass == 0.0f;
	// Static bodies never get islands
	if (is_static_a && is_static_b) return;
	int isl_a = is_static_a ? -1 : w->body_cold[ba].island_id;
	int isl_b = is_static_b ? -1 : w->body_cold[bb].island_id;
	int target;
	if (isl_a >= 0 && isl_b >= 0) {
		target = island_merge(w, isl_a, isl_b);
	} else if (isl_a >= 0) {
		target = isl_a;
	} else if (isl_b >= 0) {
		target = isl_b;
	} else {
		target = island_create(w);
	}
	// Add dynamic bodies that aren't in the target island yet
	if (!is_static_a && w->body_cold[ba].island_id != target)
		island_add_body(w, target, ba);
	if (!is_static_b && w->body_cold[bb].island_id != target)
		island_add_body(w, target, bb);
	island_add_joint(w, target, joint_idx);
	// Wake if sleeping
	if (!w->islands[target].awake)
		island_wake(w, target);
}

// Unlink a joint from its island before destruction.
static void unlink_joint_from_island(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_id < 0) return;
	int isl = j->island_id;
	island_remove_joint(w, isl, joint_idx);
	w->islands[isl].constraint_remove_count++;
}

// Update islands from this frame's contact manifolds.
// Merge/create islands for touching dynamic body pairs.
// Detect lost contacts for island splitting.
static void islands_update_contacts(WorldInternal* w, InternalManifold* manifolds, int manifold_count)
{
	CK_MAP(uint8_t) curr_touching = NULL;

	for (int i = 0; i < manifold_count; i++) {
		int a = manifolds[i].body_a, b = manifolds[i].body_b;
		int is_static_a = w->body_hot[a].inv_mass == 0.0f;
		int is_static_b = w->body_hot[b].inv_mass == 0.0f;
		// Skip static-static (shouldn't happen, but guard)
		if (is_static_a && is_static_b) continue;

		uint64_t key = body_pair_key(a, b);
		map_set(curr_touching, key, (uint8_t)1);

		// Static bodies never get island_id
		int isl_a = is_static_a ? -1 : w->body_cold[a].island_id;
		int isl_b = is_static_b ? -1 : w->body_cold[b].island_id;

		int target;
		if (isl_a >= 0 && isl_b >= 0) {
			target = island_merge(w, isl_a, isl_b);
		} else if (isl_a >= 0) {
			target = isl_a;
		} else if (isl_b >= 0) {
			target = isl_b;
		} else {
			target = island_create(w);
		}
		int added_a = 0, added_b = 0;
		if (!is_static_a && w->body_cold[a].island_id != target) {
			island_add_body(w, target, a);
			added_a = 1;
		}
		if (!is_static_b && w->body_cold[b].island_id != target) {
			island_add_body(w, target, b);
			added_b = 1;
		}
		// Wake sleeping island only if a newly added or awake body touches it.
		// Sleeping-vs-static contacts should NOT wake the island.
		if (!w->islands[target].awake) {
			int should_wake = added_a || added_b; // new body joined
			if (!should_wake) {
				// Check if either side was in an awake island before merge
				if (!is_static_a && isl_a >= 0 && isl_a != target) should_wake = 1; // merged from awake
				if (!is_static_b && isl_b >= 0 && isl_b != target) should_wake = 1;
			}
			if (should_wake) island_wake(w, target);
		}
	}

	// Detect lost contacts: iterate prev_touching, find keys not in curr_touching
	if (w->prev_touching) {
		// Walk prev_touching map entries -- use ckit map iteration
		int prev_cap = asize(w->prev_touching);
		for (int i = 0; i < prev_cap; i++) {
			// ckit maps store keys in a parallel array accessible via the header
			// We need to iterate differently -- scan body pairs
		}
		// Simpler approach: we don't have ckit map iteration exposed easily.
		// Instead, mark constraint_remove_count when a contact pair vanishes.
		// We can detect this by checking prev entries against curr.
		// For now, just use a brute-force approach: track via body pair keys.
	}

	// Swap: prev_touching = curr_touching
	map_free(w->prev_touching);
	w->prev_touching = curr_touching;
}

static void island_try_split(WorldInternal* w, int island_id); // forward decl

// Evaluate sleep for all awake islands after post-solve.
static void islands_evaluate_sleep(WorldInternal* w, float dt)
{
	int island_count = asize(w->islands);
	for (int i = 0; i < island_count; i++) {
		if (!(w->island_gen[i] & 1)) continue; // dead slot
		Island* isl = &w->islands[i];
		if (!isl->awake) continue;

		int all_sleepy = 1;
		int bi = isl->head_body;
		while (bi >= 0) {
			BodyHot* h = &w->body_hot[bi];
			float v2 = len2(h->velocity) + len2(h->angular_velocity);
			if (v2 > SLEEP_VEL_THRESHOLD) {
				h->sleep_time = 0.0f;
				all_sleepy = 0;
			} else {
				h->sleep_time += dt;
				if (h->sleep_time < SLEEP_TIME_THRESHOLD)
					all_sleepy = 0;
			}
			bi = w->body_cold[bi].island_next;
		}

		if (all_sleepy) {
			if (isl->constraint_remove_count > 0)
				island_try_split(w, i);
			else
				island_sleep(w, i);
		}
		isl->constraint_remove_count = 0;
	}
}

// DFS island split when contacts have been lost.
static void island_try_split(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];

	// Collect all body indices from island
	CK_DYNA int* bodies = NULL;
	int bi = isl->head_body;
	while (bi >= 0) {
		apush(bodies, bi);
		bi = w->body_cold[bi].island_next;
	}
	int n = asize(bodies);
	if (n <= 1) { afree(bodies); island_sleep(w, island_id); return; }

	// Visited array and component assignment
	CK_DYNA int* component = NULL;
	afit(component, n);
	asetlen(component, n);
	for (int i = 0; i < n; i++) component[i] = -1;

	// Map body_idx -> local index for fast lookup
	CK_MAP(int) body_to_local = NULL;
	for (int i = 0; i < n; i++)
		map_set(body_to_local, (uint64_t)bodies[i], i);

	int num_components = 0;
	CK_DYNA int* stack = NULL;

	for (int start = 0; start < n; start++) {
		if (component[start] >= 0) continue;
		int comp_id = num_components++;

		// DFS from start
		aclear(stack);
		apush(stack, start);
		component[start] = comp_id;

		while (asize(stack) > 0) {
			int cur_local = apop(stack);
			int cur_body = bodies[cur_local];

			// Neighbors from joints in this island
			int ji = isl->head_joint;
			while (ji >= 0) {
				JointInternal* j = &w->joints[ji];
				int other = -1;
				if (j->body_a == cur_body) other = j->body_b;
				else if (j->body_b == cur_body) other = j->body_a;
				if (other >= 0) {
					int* loc = map_get_ptr(body_to_local, (uint64_t)other);
					if (loc && component[*loc] < 0) {
						component[*loc] = comp_id;
						apush(stack, *loc);
					}
				}
				ji = w->joints[ji].island_next;
			}

			// Neighbors from current contacts (prev_touching map)
			for (int j = 0; j < n; j++) {
				if (component[j] >= 0) continue;
				uint64_t key = body_pair_key(cur_body, bodies[j]);
				uint8_t* val = map_get_ptr(w->prev_touching, key);
				if (val) {
					component[j] = comp_id;
					apush(stack, j);
				}
			}
		}
	}

	if (num_components <= 1) {
		// Still one island, just sleep it
		island_sleep(w, island_id);
	} else {
		// First component keeps the original island. Remove all bodies/joints,
		// then re-add per component.
		// Detach all bodies from the original island
		for (int i = 0; i < n; i++) {
			w->body_cold[bodies[i]].island_id = -1;
			w->body_cold[bodies[i]].island_prev = -1;
			w->body_cold[bodies[i]].island_next = -1;
		}
		// Detach all joints
		CK_DYNA int* joint_list = NULL;
		int ji = isl->head_joint;
		while (ji >= 0) {
			apush(joint_list, ji);
			int next = w->joints[ji].island_next;
			w->joints[ji].island_id = -1;
			w->joints[ji].island_prev = -1;
			w->joints[ji].island_next = -1;
			ji = next;
		}
		// Reset the original island
		*isl = (Island){ .head_body = -1, .tail_body = -1, .head_joint = -1, .tail_joint = -1, .awake = 1 };

		// Create islands per component (reuse original for component 0)
		CK_DYNA int* comp_islands = NULL;
		afit(comp_islands, num_components);
		asetlen(comp_islands, num_components);
		comp_islands[0] = island_id;
		for (int c = 1; c < num_components; c++)
			comp_islands[c] = island_create(w);

		// Assign bodies to their component's island
		for (int i = 0; i < n; i++)
			island_add_body(w, comp_islands[component[i]], bodies[i]);

		// Assign joints to the island of body_a (both endpoints should be in same component)
		for (int i = 0; i < asize(joint_list); i++) {
			int jidx = joint_list[i];
			int ba = w->joints[jidx].body_a;
			int* loc = map_get_ptr(body_to_local, (uint64_t)ba);
			int comp = loc ? component[*loc] : 0;
			island_add_joint(w, comp_islands[comp], jidx);
		}

		// Evaluate sleep for each sub-island
		for (int c = 0; c < num_components; c++) {
			int cisl = comp_islands[c];
			int all_sleepy = 1;
			int b = w->islands[cisl].head_body;
			while (b >= 0) {
				if (w->body_hot[b].sleep_time < SLEEP_TIME_THRESHOLD) { all_sleepy = 0; break; }
				b = w->body_cold[b].island_next;
			}
			if (all_sleepy) island_sleep(w, cisl);
		}

		afree(comp_islands);
		afree(joint_list);
	}

	afree(stack);
	map_free(body_to_local);
	afree(component);
	afree(bodies);
}

// ---- src/solver_avbd.c ----

// solver_avbd.c -- Augmented Vertex Block Descent solver.
// Primal-dual position solver: per-body 6x6 Newton steps + augmented Lagrangian dual updates.
// Reference: Ly, Narain et al. "Augmented VBD" SIGGRAPH 2025; Chris Giles' 3D demo.

static int avbd_body_is_sleeping(WorldInternal* w, int idx)
{
	int isl = w->body_cold[idx].island_id;
	return isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake;
}

// Adaptive alpha: scale stabilization by constraint error magnitude.
// Small error (drift): alpha_eff ≈ alpha (gentle correction, prevents jitter).
// Large error (violated): alpha_eff → 0 (fast snap to satisfied).
static float avbd_adaptive_alpha(float alpha, float err_magnitude)
{
	return alpha * fminf(1.0f, AVBD_STABLE_THRESH / fmaxf(err_magnitude, 1e-6f));
}

// Build separate CSR adjacency for contacts and joints (no type tag = no branch in inner loop).
static void avbd_build_adjacency(WorldInternal* w, AVBD_Manifold* am, int am_count, int body_count, int** out_ct_start, AVBD_ContactAdj** out_ct_adj, int** out_jt_start, AVBD_JointAdj** out_jt_adj)
{
	// --- Contact adjacency ---
	CK_DYNA int* ct_start = NULL;
	afit_set(ct_start, body_count + 1);
	memset(ct_start, 0, (body_count + 1) * sizeof(int));

	int ct_total = 0;
	for (int i = 0; i < am_count; i++)
		for (int c = 0; c < am[i].contact_count; c++) {
			ct_start[am[i].body_a]++;
			ct_start[am[i].body_b]++;
			ct_total += 2;
		}

	int sum = 0;
	for (int i = 0; i <= body_count; i++) { int tmp = ct_start[i]; ct_start[i] = sum; sum += tmp; }

	CK_DYNA AVBD_ContactAdj* ct_adj = NULL;
	if (ct_total > 0) {
		afit_set(ct_adj, ct_total);
		CK_DYNA int* off = NULL; afit_set(off, body_count);
		for (int i = 0; i < body_count; i++) off[i] = ct_start[i];
		for (int i = 0; i < am_count; i++)
			for (int c = 0; c < am[i].contact_count; c++) {
				ct_adj[off[am[i].body_a]++] = (AVBD_ContactAdj){ i, c, 1 };
				ct_adj[off[am[i].body_b]++] = (AVBD_ContactAdj){ i, c, 0 };
			}
		afree(off);
	}

	// --- Joint adjacency ---
	CK_DYNA int* jt_start = NULL;
	afit_set(jt_start, body_count + 1);
	memset(jt_start, 0, (body_count + 1) * sizeof(int));

	int jt_total = 0;
	int jt_count = asize(w->joints);
	for (int i = 0; i < jt_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		jt_start[w->joints[i].body_a]++;
		jt_start[w->joints[i].body_b]++;
		jt_total += 2;
	}

	sum = 0;
	for (int i = 0; i <= body_count; i++) { int tmp = jt_start[i]; jt_start[i] = sum; sum += tmp; }

	CK_DYNA AVBD_JointAdj* jt_adj = NULL;
	if (jt_total > 0) {
		afit_set(jt_adj, jt_total);
		CK_DYNA int* off = NULL; afit_set(off, body_count);
		for (int i = 0; i < body_count; i++) off[i] = jt_start[i];
		for (int i = 0; i < jt_count; i++) {
			if (!split_alive(w->joint_gen, i)) continue;
			jt_adj[off[w->joints[i].body_a]++] = (AVBD_JointAdj){ i, 1 };
			jt_adj[off[w->joints[i].body_b]++] = (AVBD_JointAdj){ i, 0 };
		}
		afree(off);
	}

	*out_ct_start = ct_start; *out_ct_adj = ct_adj;
	*out_jt_start = jt_start; *out_jt_adj = jt_adj;
}

// Build AVBD manifolds from collision output. Warm-start lambda/penalty from cache.
static void avbd_pre_solve(WorldInternal* w, InternalManifold* manifolds, int count, AVBD_Manifold** out_am, float dt)
{
	CK_DYNA AVBD_Manifold* am = NULL;
	float alpha = w->avbd_alpha;
	float gamma = w->avbd_gamma;

	for (int i = 0; i < count; i++) {
		InternalManifold* im = &manifolds[i];
		if (im->m.count == 0) continue;

		AVBD_Manifold m = {0};
		m.body_a = im->body_a;
		m.body_b = im->body_b;
		m.contact_count = im->m.count;

		BodyHot* a = &w->body_hot[m.body_a];
		BodyHot* b = &w->body_hot[m.body_b];
		m.friction = sqrtf(a->friction * b->friction);

		v3 n = neg(im->m.contacts[0].normal);
		v3 t1, t2;
		contact_tangent_basis(n, &t1, &t2);
		m.basis = (m3x3){{ n.x, n.y, n.z, t1.x, t1.y, t1.z, t2.x, t2.y, t2.z }};

		uint64_t key = body_pair_key(m.body_a, m.body_b);
		AVBD_WarmManifold* wm = map_get_ptr(w->avbd_warm_cache, key);

		for (int c = 0; c < m.contact_count; c++) {
			Contact* ct = &im->m.contacts[c];
			AVBD_Contact* ac = &m.contacts[c];

			v3 xA_world = add(ct->point, scale(ct->normal, ct->penetration));
			v3 xB_world = ct->point;
			ac->r_a = rotate(inv(a->rotation), sub(xA_world, a->position));
			ac->r_b = rotate(inv(b->rotation), sub(xB_world, b->position));
			ac->feature_id = ct->feature_id;

			v3 diff = sub(xA_world, xB_world);
			ac->C0 = V3(dot(n, diff) + AVBD_MARGIN, dot(t1, diff), dot(t2, diff));

			// Inertia-scaled initial penalty so first-frame forces are meaningful
			float inv_m = fmaxf(a->inv_mass, b->inv_mass);
			float pen_floor = inv_m > 0.0f ? 0.1f / (inv_m * dt * dt) : AVBD_PENALTY_MIN;
			ac->penalty = V3(pen_floor, pen_floor, pen_floor);
			ac->lambda = V3(0, 0, 0);
			ac->stick = 0;

			if (wm) {
				int matched = -1;
				for (int j = 0; j < wm->count && matched < 0; j++) {
					if (ac->feature_id != 0 && ac->feature_id == wm->contacts[j].feature_id)
						matched = j;
				}
				if (matched < 0) {
					float best_d2 = 0.01f;
					for (int j = 0; j < wm->count; j++) {
						float d2 = len2(sub(ac->r_a, wm->contacts[j].r_a));
						if (d2 < best_d2) { best_d2 = d2; matched = j; }
					}
				}
				if (matched >= 0) {
					AVBD_WarmContact* wc = &wm->contacts[matched];
					ac->penalty = wc->penalty;
					ac->lambda = wc->lambda;
					ac->stick = wc->stick;
					// Static friction: preserve contact geometry from previous frame
					// so the contact anchor doesn't slide (Giles manifold.cpp pattern).
					if (wc->stick) {
						ac->r_a = wc->r_a;
						ac->r_b = wc->r_b;
					}
				}
				ac->lambda = scale(ac->lambda, alpha * gamma);
				ac->penalty.x = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.x * gamma, AVBD_PENALTY_MAX));
				ac->penalty.y = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.y * gamma, AVBD_PENALTY_MAX));
				ac->penalty.z = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.z * gamma, AVBD_PENALTY_MAX));
			}
		}
		apush(am, m);
	}
	*out_am = am;
}

static void avbd_joints_pre_solve(WorldInternal* w, float alpha, float gamma, float dt)
{
	float dt2 = dt * dt;
	int count = asize(w->joints);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			v3 anchor_a = add(a->position, rotate(a->rotation, j->ball_socket.local_a));
			v3 anchor_b = add(b->position, rotate(b->rotation, j->ball_socket.local_b));
			j->avbd_C0_lin = sub(anchor_a, anchor_b);
		} else {
			v3 anchor_a = add(a->position, rotate(a->rotation, j->distance.local_a));
			v3 anchor_b = add(b->position, rotate(b->rotation, j->distance.local_b));
			v3 d = sub(anchor_a, anchor_b);
			float dist = len(d);
			float err = dist - j->distance.rest_length;
			v3 axis = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
			j->avbd_C0_lin = scale(axis, err);
		}

		// Inertia-scaled initial penalty: use lighter body's inertia as baseline.
		// Ensures penalty is immediately meaningful vs body mass.
		float inv_m = fmaxf(a->inv_mass, b->inv_mass); // lighter body dominates
		float inertia_floor = inv_m > 0.0f ? 0.1f / (inv_m * dt2) : AVBD_PENALTY_MIN;

		// Clamp penalty to material stiffness for spring joints (prevents over-stiffening).
		// Rigid joints (frequency=0) have infinite stiffness — no cap.
		float spring_freq = j->type == JOINT_BALL_SOCKET ? j->ball_socket.spring.frequency : j->distance.spring.frequency;
		float stiffness_cap = AVBD_PENALTY_MAX;
		if (spring_freq > 0.0f) {
			float omega = 2.0f * 3.14159265f * spring_freq;
			stiffness_cap = omega * omega;
		}

		j->avbd_lambda_lin = scale(j->avbd_lambda_lin, alpha * gamma);
		j->avbd_penalty_lin.x = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.x * gamma, stiffness_cap));
		j->avbd_penalty_lin.y = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.y * gamma, stiffness_cap));
		j->avbd_penalty_lin.z = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.z * gamma, stiffness_cap));
	}
}

static void avbd_post_solve(WorldInternal* w, AVBD_Manifold* am, int am_count)
{
	for (int i = 0; i < am_count; i++) {
		AVBD_Manifold* m = &am[i];
		uint64_t key = body_pair_key(m->body_a, m->body_b);
		AVBD_WarmManifold wm = {0};
		wm.count = m->contact_count;
		wm.stale = 0;
		for (int c = 0; c < m->contact_count; c++) {
			AVBD_Contact* ac = &m->contacts[c];
			wm.contacts[c] = (AVBD_WarmContact){
				.feature_id = ac->feature_id, .r_a = ac->r_a, .r_b = ac->r_b,
				.penalty = ac->penalty, .lambda = ac->lambda, .stick = ac->stick,
			};
		}
		map_set(w->avbd_warm_cache, key, wm);
	}
}

static void avbd_warm_cache_age_and_evict(WorldInternal* w)
{
	for (int i = 0; i < map_size(w->avbd_warm_cache); i++)
		w->avbd_warm_cache[i].stale++;
	int i = 0;
	while (i < map_size(w->avbd_warm_cache)) {
		if (w->avbd_warm_cache[i].stale > 1)
			map_del(w->avbd_warm_cache, map_keys(w->avbd_warm_cache)[i]);
		else i++;
	}
}

// --- Primal step: per-body 6x6 Newton solve ---

// Stamp a contact into the per-body 6x6 system.
static void avbd_stamp_contact(AVBD_Manifold* m, AVBD_Contact* ac, int is_body_a, AVBD_BodyState* states, float alpha, BodyHot* ha, BodyHot* hb, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 dpA = sub(ha->position, states[m->body_a].initial_lin);
	v3 daA = quat_sub_angular(ha->rotation, states[m->body_a].initial_ang);
	v3 dpB = sub(hb->position, states[m->body_b].initial_lin);
	v3 daB = quat_sub_angular(hb->rotation, states[m->body_b].initial_ang);

	v3 rAW = rotate(ha->rotation, ac->r_a);
	v3 rBW = rotate(hb->rotation, ac->r_b);

	v3 b0 = V3(m->basis.m[0], m->basis.m[1], m->basis.m[2]);
	v3 b1 = V3(m->basis.m[3], m->basis.m[4], m->basis.m[5]);
	v3 b2 = V3(m->basis.m[6], m->basis.m[7], m->basis.m[8]);

	v3 jALin[3] = { b0, b1, b2 };
	v3 jBLin[3] = { neg(b0), neg(b1), neg(b2) };
	v3 jAAng[3] = { cross(rAW, b0), cross(rAW, b1), cross(rAW, b2) };
	v3 jBAng[3] = { cross(rBW, jBLin[0]), cross(rBW, jBLin[1]), cross(rBW, jBLin[2]) };

	float C[3];
	C[0] = ac->C0.x*(1-alpha) + dot(jALin[0],dpA) + dot(jBLin[0],dpB) + dot(jAAng[0],daA) + dot(jBAng[0],daB);
	C[1] = ac->C0.y*(1-alpha) + dot(jALin[1],dpA) + dot(jBLin[1],dpB) + dot(jAAng[1],daA) + dot(jBAng[1],daB);
	C[2] = ac->C0.z*(1-alpha) + dot(jALin[2],dpA) + dot(jBLin[2],dpB) + dot(jAAng[2],daA) + dot(jBAng[2],daB);

	float F[3];
	F[0] = ac->penalty.x * C[0] + ac->lambda.x;
	F[1] = ac->penalty.y * C[1] + ac->lambda.y;
	F[2] = ac->penalty.z * C[2] + ac->lambda.z;

	if (F[0] > 0.0f) F[0] = 0.0f;

	float bounds = fabsf(F[0]) * m->friction;
	float fric_len = sqrtf(F[1]*F[1] + F[2]*F[2]);
	if (fric_len > bounds && fric_len > 0.0f) {
		float s = bounds / fric_len;
		F[1] *= s; F[2] *= s;
	}

	v3* jL = is_body_a ? jALin : jBLin;
	v3* jA = is_body_a ? jAAng : jBAng;

	for (int r = 0; r < 3; r++) {
		float kk = (&ac->penalty.x)[r];
		*rhsLin = add(*rhsLin, scale(jL[r], F[r]));
		*rhsAng = add(*rhsAng, scale(jA[r], F[r]));
		for (int ri = 0; ri < 3; ri++)
			for (int ci = 0; ci < 3; ci++) {
				lhsLin->m[ri*3+ci] += (&jL[r].x)[ri] * (&jL[r].x)[ci] * kk;
				lhsAng->m[ri*3+ci] += (&jA[r].x)[ri] * (&jA[r].x)[ci] * kk;
				lhsCross->m[ri*3+ci] += (&jA[r].x)[ri] * (&jL[r].x)[ci] * kk;
			}
	}

	// Geometric stiffness for contacts (Sec 3.5): diagonal approximation of
	// higher-order Hessian. Accounts for contact offset rotating with body.
	{
		v3 r = is_body_a ? rAW : neg(rBW);
		float H[9] = {0};
		for (int k = 0; k < 3; k++) {
			float fk = F[k];
			for (int ri = 0; ri < 3; ri++)
				for (int ci = 0; ci < 3; ci++) {
					float val = (ri == ci ? -(&r.x)[k] : 0.0f) + (ci == k ? (&r.x)[ri] : 0.0f);
					H[ri*3+ci] += val * fk;
				}
		}
		for (int ci = 0; ci < 3; ci++) {
			float col_len = sqrtf(H[0*3+ci]*H[0*3+ci] + H[1*3+ci]*H[1*3+ci] + H[2*3+ci]*H[2*3+ci]);
			lhsAng->m[ci*3+ci] += col_len;
		}
	}
}

// Stamp a ball-socket joint into the per-body 6x6 system.
static void avbd_stamp_ball_socket(JointInternal* j, BodyHot* a, BodyHot* b, int is_body_a, float alpha, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 rAW = rotate(a->rotation, j->ball_socket.local_a);
	v3 rBW = rotate(b->rotation, j->ball_socket.local_b);

	v3 C_vec = sub(add(a->position, rAW), add(b->position, rBW));
	if (j->ball_socket.spring.frequency <= 0.0f) {
		float aa = avbd_adaptive_alpha(alpha, len(j->avbd_C0_lin));
		C_vec = sub(C_vec, scale(j->avbd_C0_lin, aa));
	}

	v3 K = j->avbd_penalty_lin;
	v3 F = add(hmul(K, C_vec), j->avbd_lambda_lin);

	// Jacobians: J_lin_a = I, J_lin_b = -I, J_ang_a = skew(-rAW), J_ang_b = skew(rBW)
	float sign = is_body_a ? 1.0f : -1.0f;
	v3 r_world = is_body_a ? rAW : rBW;
	m3x3 jAng = is_body_a ? skew(neg(rAW)) : skew(rBW);

	// Stamp LHS: I^T * diag(K) * I = diag(K) for linear, jAng^T * diag(K) * jAng for angular
	for (int ri = 0; ri < 3; ri++)
		lhsLin->m[ri*3+ri] += (&K.x)[ri]; // diag(K) since J_lin = +/-I

	m3x3 jAngT = transpose(jAng);
	m3x3 Kmat = diag(K);
	m3x3 jAngTK = mul(jAngT, Kmat);
	*lhsAng = add(*lhsAng, mul(jAngTK, jAng));
	*lhsCross = add(*lhsCross, scale(jAngTK, sign));

	// Geometric stiffness: diagonal approximation of higher-order Hessian term (Sec 3.5).
	// Accounts for the nonlinear coupling between rotation and anchor position.
	// Without this, the angular Newton step overshoots for joints under load.
	{
		v3 r = is_body_a ? rAW : neg(rBW);
		// geometricStiffnessBallSocket(k, r) builds a matrix M where M = diag(-r[k]) + e_k * r^T
		// Sum over k weighted by F[k], then diagonalize (column norms on diagonal).
		float H[9] = {0};
		for (int k = 0; k < 3; k++) {
			float fk = (&F.x)[k];
			for (int ri = 0; ri < 3; ri++)
				for (int ci = 0; ci < 3; ci++) {
					float val = (ri == ci ? -(&r.x)[k] : 0.0f) + (ci == k ? (&r.x)[ri] : 0.0f);
					H[ri*3+ci] += val * fk;
				}
		}
		// Diagonalize: column norms
		for (int ci = 0; ci < 3; ci++) {
			float col_len = sqrtf(H[0*3+ci]*H[0*3+ci] + H[1*3+ci]*H[1*3+ci] + H[2*3+ci]*H[2*3+ci]);
			lhsAng->m[ci*3+ci] += col_len;
		}
	}

	// Stamp RHS: J^T * F
	*rhsLin = add(*rhsLin, scale(F, sign));
	*rhsAng = add(*rhsAng, m3x3_mul_v3(jAngT, F));
}

// Stamp a distance joint into the per-body 6x6 system.
static void avbd_stamp_distance(JointInternal* j, BodyHot* a, BodyHot* b, int is_body_a, float alpha, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 rAW = rotate(a->rotation, j->distance.local_a);
	v3 rBW = rotate(b->rotation, j->distance.local_b);
	v3 d = sub(add(a->position, rAW), add(b->position, rBW));
	float dist = len(d);
	v3 axis = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
	float err = dist - j->distance.rest_length;

	float k = j->avbd_penalty_lin.x;
	float f_scalar = k * err + j->avbd_lambda_lin.x;

	if (j->distance.spring.frequency <= 0.0f) {
		float C0_scalar = len(j->avbd_C0_lin);
		if (dot(j->avbd_C0_lin, axis) < 0) C0_scalar = -C0_scalar;
		float aa = avbd_adaptive_alpha(alpha, fabsf(C0_scalar));
		f_scalar = k * (err - aa * C0_scalar) + j->avbd_lambda_lin.x;
	}

	v3 jLin_vec = is_body_a ? axis : neg(axis);
	v3 r_world = is_body_a ? rAW : rBW;
	v3 jAng_vec = cross(r_world, jLin_vec);

	for (int ri = 0; ri < 3; ri++)
		for (int ci = 0; ci < 3; ci++) {
			lhsLin->m[ri*3+ci] += (&jLin_vec.x)[ri] * (&jLin_vec.x)[ci] * k;
			lhsAng->m[ri*3+ci] += (&jAng_vec.x)[ri] * (&jAng_vec.x)[ci] * k;
			lhsCross->m[ri*3+ci] += (&jAng_vec.x)[ri] * (&jLin_vec.x)[ci] * k;
		}

	v3 f_lin = scale(jLin_vec, f_scalar);
	v3 f_ang = cross(r_world, f_lin);
	*rhsLin = add(*rhsLin, f_lin);
	*rhsAng = add(*rhsAng, f_ang);
}

static void avbd_primal_step(WorldInternal* w, AVBD_Manifold* am, int am_count, int* ct_start, AVBD_ContactAdj* ct_adj, int* jt_start, AVBD_JointAdj* jt_adj, AVBD_BodyState* states, float alpha, float dt)
{
	int body_count = asize(w->body_hot);
	float dt2 = dt * dt;

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;

		float mass = 1.0f / h->inv_mass;
		v3 moment = rcp(h->inv_inertia_local);

		m3x3 lhsLin = diag(V3(mass/dt2, mass/dt2, mass/dt2));
		m3x3 lhsAng = diag(V3(moment.x/dt2, moment.y/dt2, moment.z/dt2));
		m3x3 lhsCross = {0};

		v3 dp = sub(h->position, states[i].inertial_lin);
		v3 da = quat_sub_angular(h->rotation, states[i].inertial_ang);
		v3 rhsLin = V3(mass/dt2 * dp.x, mass/dt2 * dp.y, mass/dt2 * dp.z);
		v3 rhsAng = V3(moment.x/dt2 * da.x, moment.y/dt2 * da.y, moment.z/dt2 * da.z);

		// Pass 1: contacts (tight loop, no branching)
		for (int j = ct_start[i]; j < ct_start[i + 1]; j++) {
			AVBD_ContactAdj* ca = &ct_adj[j];
			AVBD_Manifold* m = &am[ca->manifold_idx];
			avbd_stamp_contact(m, &m->contacts[ca->contact_idx], ca->is_body_a, states, alpha,
			                   &w->body_hot[m->body_a], &w->body_hot[m->body_b],
			                   &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
		}

		// Pass 2: joints (tight loop, no branching)
		for (int j = jt_start[i]; j < jt_start[i + 1]; j++) {
			AVBD_JointAdj* ja = &jt_adj[j];
			JointInternal* jt = &w->joints[ja->joint_idx];
			BodyHot* a = &w->body_hot[jt->body_a];
			BodyHot* b = &w->body_hot[jt->body_b];
			if (jt->type == JOINT_BALL_SOCKET)
				avbd_stamp_ball_socket(jt, a, b, ja->is_body_a, alpha, &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
			else
				avbd_stamp_distance(jt, a, b, ja->is_body_a, alpha, &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
		}

		v3 dxLin, dxAng;
		solve_6x6_ldl(lhsLin, lhsAng, lhsCross, neg(rhsLin), neg(rhsAng), &dxLin, &dxAng);
		h->position = add(h->position, dxLin);
		h->rotation = quat_add_dw(h->rotation, dxAng);
	}
}

// --- Dual steps ---

static void avbd_contacts_dual_step(WorldInternal* w, AVBD_Manifold* am, int am_count, AVBD_BodyState* states, float alpha)
{
	float beta_lin = w->avbd_beta_lin;

	for (int i = 0; i < am_count; i++) {
		AVBD_Manifold* m = &am[i];
		BodyHot* ha = &w->body_hot[m->body_a];
		BodyHot* hb = &w->body_hot[m->body_b];

		v3 dpA = sub(ha->position, states[m->body_a].initial_lin);
		v3 daA = quat_sub_angular(ha->rotation, states[m->body_a].initial_ang);
		v3 dpB = sub(hb->position, states[m->body_b].initial_lin);
		v3 daB = quat_sub_angular(hb->rotation, states[m->body_b].initial_ang);

		v3 b0 = V3(m->basis.m[0], m->basis.m[1], m->basis.m[2]);
		v3 b1 = V3(m->basis.m[3], m->basis.m[4], m->basis.m[5]);
		v3 b2 = V3(m->basis.m[6], m->basis.m[7], m->basis.m[8]);

		for (int c = 0; c < m->contact_count; c++) {
			AVBD_Contact* ac = &m->contacts[c];
			v3 rAW = rotate(ha->rotation, ac->r_a);
			v3 rBW = rotate(hb->rotation, ac->r_b);

			v3 nb0 = neg(b0), nb1 = neg(b1), nb2 = neg(b2);
			v3 jAAng[3] = { cross(rAW, b0), cross(rAW, b1), cross(rAW, b2) };
			v3 jBAng[3] = { cross(rBW, nb0), cross(rBW, nb1), cross(rBW, nb2) };

			float Cn = ac->C0.x*(1-alpha) + dot(b0,dpA) + dot(nb0,dpB) + dot(jAAng[0],daA) + dot(jBAng[0],daB);
			float Ct1 = ac->C0.y*(1-alpha) + dot(b1,dpA) + dot(nb1,dpB) + dot(jAAng[1],daA) + dot(jBAng[1],daB);
			float Ct2 = ac->C0.z*(1-alpha) + dot(b2,dpA) + dot(nb2,dpB) + dot(jAAng[2],daA) + dot(jBAng[2],daB);

			float Fn = ac->penalty.x * Cn + ac->lambda.x;
			float Ft1 = ac->penalty.y * Ct1 + ac->lambda.y;
			float Ft2 = ac->penalty.z * Ct2 + ac->lambda.z;
			if (Fn > 0.0f) Fn = 0.0f;

			float bounds = fabsf(Fn) * m->friction;
			float fric_len = sqrtf(Ft1*Ft1 + Ft2*Ft2);
			if (fric_len > bounds && fric_len > 0.0f) { float s = bounds/fric_len; Ft1 *= s; Ft2 *= s; }

			ac->lambda = V3(Fn, Ft1, Ft2);

			if (Fn < 0.0f)
				ac->penalty.x = fminf(ac->penalty.x + beta_lin * fabsf(Cn), AVBD_PENALTY_MAX);
			if (fric_len <= bounds) {
				ac->penalty.y = fminf(ac->penalty.y + beta_lin * fabsf(Ct1), AVBD_PENALTY_MAX);
				ac->penalty.z = fminf(ac->penalty.z + beta_lin * fabsf(Ct2), AVBD_PENALTY_MAX);
				ac->stick = sqrtf(Ct1*Ct1 + Ct2*Ct2) < AVBD_STICK_THRESH;
			}
		}
	}
}

static void avbd_joints_dual_step(WorldInternal* w, float alpha)
{
	float beta_lin = w->avbd_beta_lin;
	int count = asize(w->joints);

	for (int i = 0; i < count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			// Giles pattern: C is modified in-place by stabilization.
			// Lambda update and penalty ramp both use the stabilized C.
			v3 C = sub(add(a->position, rotate(a->rotation, j->ball_socket.local_a)),
			           add(b->position, rotate(b->rotation, j->ball_socket.local_b)));
			if (j->ball_socket.spring.frequency <= 0.0f) {
				float aa = avbd_adaptive_alpha(alpha, len(j->avbd_C0_lin));
				C = sub(C, scale(j->avbd_C0_lin, aa));
			}
			v3 F = add(hmul(j->avbd_penalty_lin, C), j->avbd_lambda_lin);
			j->avbd_lambda_lin = F;

			v3 absC = V3(fabsf(C.x), fabsf(C.y), fabsf(C.z));
			float cap = AVBD_PENALTY_MAX;
			if (j->ball_socket.spring.frequency > 0.0f) {
				float om = 2.0f * 3.14159265f * j->ball_socket.spring.frequency;
				cap = om * om;
			}
			j->avbd_penalty_lin = V3(
				fminf(j->avbd_penalty_lin.x + beta_lin * absC.x, cap),
				fminf(j->avbd_penalty_lin.y + beta_lin * absC.y, cap),
				fminf(j->avbd_penalty_lin.z + beta_lin * absC.z, cap));
		} else {
			v3 d = sub(add(a->position, rotate(a->rotation, j->distance.local_a)),
			          add(b->position, rotate(b->rotation, j->distance.local_b)));
			float dist = len(d);
			float err = dist - j->distance.rest_length;
			// In-place stabilization (matching Giles pattern)
			if (j->distance.spring.frequency <= 0.0f) {
				float C0_scalar = len(j->avbd_C0_lin);
				v3 ax = dist > 1e-6f ? scale(d, 1.0f/dist) : V3(0,1,0);
				if (dot(j->avbd_C0_lin, ax) < 0) C0_scalar = -C0_scalar;
				float aa = avbd_adaptive_alpha(alpha, fabsf(C0_scalar));
				err -= aa * C0_scalar;
			}
			float F = j->avbd_penalty_lin.x * err + j->avbd_lambda_lin.x;
			j->avbd_lambda_lin.x = F;
			float cap = AVBD_PENALTY_MAX;
			if (j->distance.spring.frequency > 0.0f) {
				float om = 2.0f * 3.14159265f * j->distance.spring.frequency;
				cap = om * om;
			}
			j->avbd_penalty_lin.x = fminf(j->avbd_penalty_lin.x + beta_lin * fabsf(err), cap);
		}
	}
}

// --- Top-level AVBD solver ---

static void avbd_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, float dt)
{
	int body_count = asize(w->body_hot);
	float alpha = w->avbd_alpha;
	float gamma = w->avbd_gamma;

	avbd_warm_cache_age_and_evict(w);

	CK_DYNA AVBD_Manifold* am = NULL;
	avbd_pre_solve(w, manifolds, manifold_count, &am, dt);
	int am_count = asize(am);

	avbd_joints_pre_solve(w, alpha, gamma, dt);

	int *ct_start = NULL, *jt_start = NULL;
	AVBD_ContactAdj* ct_adj = NULL;
	AVBD_JointAdj* jt_adj = NULL;
	avbd_build_adjacency(w, am, am_count, body_count, &ct_start, &ct_adj, &jt_start, &jt_adj);

	// Ensure prev_velocity array is big enough
	if (asize(w->avbd_prev_velocity) < body_count) {
		afit(w->avbd_prev_velocity, body_count);
		asetlen(w->avbd_prev_velocity, body_count);
	}

	CK_DYNA AVBD_BodyState* states = NULL;
	afit_set(states, body_count);
	memset(states, 0, body_count * sizeof(AVBD_BodyState));

	float inv_dt = 1.0f / dt;
	v3 grav = w->gravity;
	float grav_len = len(grav);

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		states[i].initial_lin = h->position;
		states[i].initial_ang = h->rotation;

		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;

		// Inertial position: y = x + v*dt + g*dt^2
		states[i].inertial_lin = add(add(h->position, scale(h->velocity, dt)), scale(grav, dt * dt));
		states[i].inertial_ang = quat_add_dw(h->rotation, scale(h->angular_velocity, dt));

		// Adaptive body warm-start (VBD paper): use acceleration from previous
		// frame to predict gravity contribution, giving solver a better starting guess.
		v3 prev_v = w->avbd_prev_velocity[i];
		v3 accel = scale(sub(h->velocity, prev_v), inv_dt);
		float accel_along_grav = grav_len > 0.0f ? dot(accel, scale(grav, 1.0f / grav_len)) : 0.0f;
		float accel_weight = grav_len > 0.0f ? fmaxf(0.0f, fminf(accel_along_grav / grav_len, 1.0f)) : 0.0f;

		h->position = add(h->position, add(scale(h->velocity, dt), scale(grav, accel_weight * dt * dt)));
		h->rotation = quat_add_dw(h->rotation, scale(h->angular_velocity, dt));
	}

	for (int it = 0; it < w->avbd_iterations; it++) {
		avbd_primal_step(w, am, am_count, ct_start, ct_adj, jt_start, jt_adj, states, alpha, dt);
		avbd_contacts_dual_step(w, am, am_count, states, alpha);
		avbd_joints_dual_step(w, alpha);
	}

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;
		w->avbd_prev_velocity[i] = h->velocity; // save for next frame's adaptive warm-start
		h->velocity = scale(sub(h->position, states[i].initial_lin), inv_dt);
		h->angular_velocity = scale(quat_sub_angular(h->rotation, states[i].initial_ang), inv_dt);
	}

	avbd_post_solve(w, am, am_count);

	afree(states);
	afree(ct_adj); afree(ct_start);
	afree(jt_adj); afree(jt_start);
	afree(am);
}

// -----------------------------------------------------------------------------
// World.

World create_world(WorldParams params)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	w->gravity = params.gravity;
	w->broadphase_type = params.broadphase;
	w->friction_model = params.friction_model;
	w->solver_type = params.solver_type;
	w->sleep_enabled = 1;
	w->velocity_iters = params.velocity_iters > 0 ? params.velocity_iters : SOLVER_VELOCITY_ITERS;
	w->position_iters = params.position_iters > 0 ? params.position_iters : SOLVER_POSITION_ITERS;
	w->contact_hertz = params.contact_hertz > 0.0f ? params.contact_hertz : 60.0f;
	w->contact_damping_ratio = params.contact_damping_ratio > 0.0f ? params.contact_damping_ratio : 3.0f;
	w->max_push_velocity = params.max_push_velocity > 0.0f ? params.max_push_velocity : 3.0f;
	w->sub_steps = params.sub_steps > 0 ? params.sub_steps : 4;
	w->avbd_alpha = 0.99f;
	w->avbd_beta_lin = 10000.0f;
	w->avbd_beta_ang = 100.0f;
	w->avbd_gamma = 0.999f;
	w->avbd_iterations = 20;
	w->bvh_static = CK_ALLOC(sizeof(BVHTree));
	w->bvh_dynamic = CK_ALLOC(sizeof(BVHTree));
	bvh_init(w->bvh_static);
	bvh_init(w->bvh_dynamic);
	return (World){ (uint64_t)w };
}

void destroy_world(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	for (int i = 0; i < asize(w->body_cold); i++) {
		afree(w->body_cold[i].shapes);
	}
	afree(w->debug_contacts);
	map_free(w->warm_cache);
	map_free(w->avbd_warm_cache);
	afree(w->avbd_prev_velocity);
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	CK_FREE(w);
}

// Integrate velocities for a sub-step (gravity + damping).
static void integrate_velocities(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
		if (h->linear_damping > 0.0f)
			h->velocity = scale(h->velocity, 1.0f / (1.0f + h->linear_damping * dt));
		if (h->angular_damping > 0.0f)
			h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + h->angular_damping * dt));
	}
}

// Integrate positions and rotations for a sub-step.
static void integrate_positions(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;

		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL)
			h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL)
			h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));

		h->position = add(h->position, scale(h->velocity, dt));

		h->angular_velocity = solve_gyroscopic(h->rotation, h->inv_inertia_local, h->angular_velocity, dt);

		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, h->rotation);
		h->rotation.x += 0.5f * dt * dq.x;
		h->rotation.y += 0.5f * dt * dq.y;
		h->rotation.z += 0.5f * dt * dq.z;
		h->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y
			+ h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql;
		h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
	}
}

void world_step(World world, float dt)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->frame++;
	int n_sub = w->sub_steps;
	float sub_dt = dt / (float)n_sub;

	// Age warm cache once per frame (AVBD does its own in avbd_solve)
	if (w->solver_type != SOLVER_AVBD)
		warm_cache_age_and_evict(w);

	// --- Collision detection (once per frame) ---
	// Dual solvers need velocity integration before collision for first substep.
	// AVBD handles gravity internally via inertial position — skip.
	if (w->solver_type != SOLVER_AVBD)
		integrate_velocities(w, sub_dt);

	CK_DYNA InternalManifold* manifolds = NULL;
	broadphase_and_collide(w, &manifolds);
	islands_update_contacts(w, manifolds, asize(manifolds));

	aclear(w->debug_contacts);
	for (int i = 0; i < asize(manifolds); i++)
		for (int c = 0; c < manifolds[i].m.count; c++)
			apush(w->debug_contacts, manifolds[i].m.contacts[c]);

	int manifold_count = asize(manifolds);

	// AVBD takes a completely different path (primal-dual position solver)
	if (w->solver_type == SOLVER_AVBD) {
		avbd_solve(w, manifolds, manifold_count, dt);
		if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
		afree(manifolds);
		return;
	}

	// --- Pre-solve (once per frame, using sub_dt for softness/bias) ---
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, sub_dt);

	SolverBallSocket* sol_bs = NULL;
	SolverDistance*    sol_dist = NULL;
	joints_pre_solve(w, sub_dt, &sol_bs, &sol_dist);
	joints_warm_start(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	// --- Graph color (once per frame) ---
	int count = asize(w->body_hot);
	CK_DYNA ConstraintRef* crefs = NULL;
	int sm_count = asize(sm);
	for (int i = 0; i < sm_count; i++) {
		ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
			.body_a = sm[i].body_a, .body_b = sm[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_bs); i++) {
		ConstraintRef r = { .type = CTYPE_BALL_SOCKET, .index = i,
			.body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_dist); i++) {
		ConstraintRef r = { .type = CTYPE_DISTANCE, .index = i,
			.body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b };
		apush(crefs, r);
	}

	int cref_count = asize(crefs);
	int batch_starts[65] = {0};
	int color_count = 0;
	if (cref_count > 0)
		color_constraints(crefs, cref_count, count, batch_starts, &color_count);

	// --- Sub-step loop: velocity solve + position integrate ---
	// First sub-step: velocities already integrated above.
	// Joint bias applied on first sub-step only (error is stale after position integration).
	for (int sub = 0; sub < n_sub; sub++) {
		if (sub > 0)
			integrate_velocities(w, sub_dt);

		for (int iter = 0; iter < w->velocity_iters; iter++)
			for (int c = 0; c < color_count; c++)
				for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
					solve_constraint(w, &crefs[i], sm, sc, sol_bs, sol_dist);

		integrate_positions(w, sub_dt);

		// Relax contacts: refresh separation/bias from updated positions
		if (w->solver_type == SOLVER_SOFT_STEP || w->solver_type == SOLVER_BLOCK)
			solver_relax_contacts(w, sm, asize(sm), sc, sub_dt);

		// Zero rigid joint bias after first sub-step (position error is stale)
		if (sub == 0) {
			for (int i = 0; i < asize(sol_bs); i++)
				if (sol_bs[i].softness == 0.0f) sol_bs[i].bias = V3(0, 0, 0);
			for (int i = 0; i < asize(sol_dist); i++)
				if (sol_dist[i].softness == 0.0f) sol_dist[i].bias = 0.0f;
		}
	}

	afree(crefs);

	// Position correction: NGS for hard SI; also for soft modes when contact_hertz is off
	if (w->solver_type == SOLVER_SI)
		solver_position_correct(w, sm, asize(sm), sc);
	else if (w->contact_hertz <= 0.0f)
		solver_position_correct(w, sm, asize(sm), sc);

	// Post-solve (once per frame)
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	if (w->sleep_enabled) islands_evaluate_sleep(w, dt);

	afree(manifolds);
}

void world_set_friction_model(World world, FrictionModel model)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->friction_model = model;
}

void world_set_solver_type(World world, SolverType type)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->solver_type = type;
}

// -----------------------------------------------------------------------------
// Body.

Body create_body(World world, BodyParams params)
{
	assert(is_valid(params.position) && "create_body: position is NaN/inf");
	assert(is_valid(params.rotation) && "create_body: rotation is NaN/inf");
	assert(is_valid(params.mass) && params.mass >= 0.0f && "create_body: mass must be >= 0 and finite");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	split_add(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);

	w->body_cold[idx] = (BodyCold){
		.mass = params.mass,
		.shapes = NULL,
		.bvh_leaf = -1,
		.island_id = -1,
		.island_prev = -1,
		.island_next = -1,
	};
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	float ang_damp = params.angular_damping;
	if (ang_damp == 0.0f) ang_damp = 0.03f; // default: 3%/s (BEPU-style)
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
		.friction = fric,
		.restitution = params.restitution,
		.linear_damping = params.linear_damping,
		.angular_damping = ang_damp,
	};

	return split_handle(Body, w->body_gen, idx);
}

void destroy_body(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	// Remove from island
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl)) {
		// Remove all joints connected to this body
		int ji = w->islands[isl].head_joint;
		while (ji >= 0) {
			int next = w->joints[ji].island_next;
			if (w->joints[ji].body_a == idx || w->joints[ji].body_b == idx) {
				island_remove_joint(w, isl, ji);
				w->islands[isl].constraint_remove_count++;
			}
			ji = next;
		}
		island_remove_body(w, isl, idx);
		w->islands[isl].constraint_remove_count++;
	}
	if (w->body_cold[idx].bvh_leaf >= 0) {
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		int moved_body = bvh_remove(tree, w->body_cold[idx].bvh_leaf);
		if (moved_body >= 0) w->body_cold[moved_body].bvh_leaf = w->body_cold[idx].bvh_leaf;
	}
	afree(w->body_cold[idx].shapes);
	split_del(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);
}

void body_add_shape(World world, Body body, ShapeParams params)
{
	assert(is_valid(params.local_pos) && "body_add_shape: local_pos is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));

	ShapeInternal s = {0};
	s.type = params.type;
	s.local_pos = params.local_pos;
	switch (params.type) {
	case SHAPE_SPHERE:  s.sphere.radius = params.sphere.radius; break;
	case SHAPE_CAPSULE: s.capsule.half_height = params.capsule.half_height;
	                    s.capsule.radius = params.capsule.radius; break;
	case SHAPE_BOX:     s.box.half_extents = params.box.half_extents; break;
	case SHAPE_HULL:    s.hull.hull = params.hull.hull;
	                    s.hull.scale = params.hull.scale; break;
	}
	apush(w->body_cold[idx].shapes, s);
	recompute_body_inertia(w, idx);

	// Insert into BVH on first shape add.
	if (w->broadphase_type == BROADPHASE_BVH && asize(w->body_cold[idx].shapes) == 1) {
		AABB box = aabb_expand(body_aabb(&w->body_hot[idx], &w->body_cold[idx]), BVH_AABB_MARGIN);
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		w->body_cold[idx].bvh_leaf = bvh_insert(tree, idx, box);
	}
}

v3 body_get_position(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].position;
}

quat body_get_rotation(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].rotation;
}

void body_wake(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_velocity(World world, Body body, v3 vel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].velocity = vel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_angular_velocity(World world, Body body, v3 avel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].angular_velocity = avel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

int body_is_asleep(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl < 0 || !island_alive(w, isl)) return 0;
	return !w->islands[isl].awake;
}

// -----------------------------------------------------------------------------
// Joints.

Joint create_ball_socket(World world, BallSocketParams params)
{
	assert(is_valid(params.local_offset_a) && "create_ball_socket: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_ball_socket: local_offset_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	// Grow joint arrays manually (no split_add -- joints don't need hot/cold split)
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1); // odd = alive
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = ba, .body_b = bb,
		.ball_socket = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_distance(World world, DistanceParams params)
{
	assert(is_valid(params.local_offset_a) && "create_distance: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_distance: local_offset_b is NaN/inf");
	assert(is_valid(params.rest_length) && "create_distance: rest_length is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	// Auto-compute rest length if not specified
	float rest = params.rest_length;
	if (rest <= 0.0f) {
		BodyHot* a = &w->body_hot[ba];
		BodyHot* b = &w->body_hot[bb];
		v3 wa = add(a->position, rotate(a->rotation, params.local_offset_a));
		v3 wb = add(b->position, rotate(b->rotation, params.local_offset_b));
		rest = len(sub(wb, wa));
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_DISTANCE,
		.body_a = ba, .body_b = bb,
		.distance = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.rest_length = rest,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	unlink_joint_from_island(w, idx);
	memset(&w->joints[idx], 0, sizeof(JointInternal));
	w->joint_gen[idx]++; // even = dead
	apush(w->joint_free, idx);
}

static void bvh_debug_walk(BVHTree* t, int ni, int depth, BVHDebugFn fn, void* user)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		fn(c->min, c->max, depth, bvh_child_is_leaf(c), user);
		if (bvh_child_is_internal(c)) bvh_debug_walk(t, c->index, depth + 1, fn, user);
	}
}

void world_debug_bvh(World world, BVHDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (w->bvh_dynamic->root >= 0) bvh_debug_walk(w->bvh_dynamic, w->bvh_dynamic->root, 0, fn, user);
	if (w->bvh_static->root >= 0) bvh_debug_walk(w->bvh_static, w->bvh_static->root, 0, fn, user);
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}

#endif // NUDGE_IMPLEMENTATION
#endif // NUDGE_SINGLE_FILE_H
