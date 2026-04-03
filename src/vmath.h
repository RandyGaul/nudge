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
