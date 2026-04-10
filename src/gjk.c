// See LICENSE for licensing info.
// gjk.c -- GJK distance algorithm for 3D convex shapes.
//
// Per-shape-type support dispatch:
//   - Sphere/capsule: world-space core point/segment, radius applied post-hoc
//   - Box: analytical sign-based corner selection
//   - Hull: local-space vertex scan with rotation
//   - Cylinder: implicit surface support (no tessellation)
// Termination: 3-condition (containment, monotonic progress, support progress).

// -----------------------------------------------------------------------------
// SSE v3: a v3 packed into __m128 (w=0). Used internally for SIMD math.

#include <xmmintrin.h>
#include <smmintrin.h> // SSE4.1 for _mm_dp_ps

typedef __m128 s3; // SSE v3: [x, y, z, 0]

#define s3_load(a)       _mm_set_ps(0, (a).z, (a).y, (a).x)
#define s3_add(a, b)     _mm_add_ps(a, b)
#define s3_sub(a, b)     _mm_sub_ps(a, b)
#define s3_neg(a)        _mm_sub_ps(_mm_setzero_ps(), a)
#define s3_scale(a, s)   _mm_mul_ps(a, _mm_set1_ps(s))
#define s3_dot(a, b)     _mm_cvtss_f32(_mm_dp_ps(a, b, 0x71))
#define s3_len2(a)       s3_dot(a, a)
#define s3_stp(a, b, c)  s3_dot(a, s3_cross(b, c))

static inline v3 s3_store(s3 a) { float t[4]; _mm_storeu_ps(t, a); return V3(t[0], t[1], t[2]); }

static inline s3 s3_cross(s3 a, s3 b) {
	s3 a1 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1));
	s3 b1 = _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1));
	s3 r = _mm_sub_ps(_mm_mul_ps(a, b1), _mm_mul_ps(a1, b));
	return _mm_shuffle_ps(r, r, _MM_SHUFFLE(3,0,2,1));
}

static inline s3 s3_quat_rotate(quat q, s3 v) {
	s3 u = _mm_set_ps(0, q.z, q.y, q.x);
	float s = q.w;
	float uv = s3_dot(u, v);
	float uu = s3_dot(u, u);
	s3 uXv = s3_cross(u, v);
	return s3_add(s3_add(s3_scale(u, 2.0f * uv), s3_scale(v, s*s - uu)), s3_scale(uXv, 2.0f * s));
}

// -----------------------------------------------------------------------------
// Simplex types.

typedef struct GJK_Vertex
{
	s3 point1;  // support point on shape A (SSE __m128)
	s3 point2;  // support point on shape B (SSE __m128)
	s3 point;   // Minkowski difference: point2 - point1
	float u;    // barycentric coordinate
	int feat1;  // feature ID on shape A
	int feat2;  // feature ID on shape B
} GJK_Vertex;

typedef struct GJK_Simplex
{
	GJK_Vertex v[4];
	float divisor;
	int count;
} GJK_Simplex;

// Byte size of a simplex with n active vertices (vertices + divisor + count).
#define gjk_simplex_size(n) (sizeof(GJK_Vertex) * (n) + sizeof(float) + sizeof(int))

typedef struct GJK_Result
{
	v3 point1;  // closest point on shape A
	v3 point2;  // closest point on shape B
	float distance;
	int iterations;
	int feat1;  // feature ID on shape A (from contributing simplex vertex)
	int feat2;  // feature ID on shape B
} GJK_Result;

// -----------------------------------------------------------------------------
// Shape types and constructors.

enum { GJK_POINT, GJK_SEGMENT, GJK_BOX, GJK_HULL, GJK_CYLINDER };

typedef struct GJK_Shape
{
	int type;
	float radius; // applied post-hoc (sphere/capsule); 0 for surface shapes
	union {
		struct { v3 center; } point;
		struct { v3 p, q; } segment;
		struct { v3 center; quat rot; quat inv_rot; v3 half_extents; } box;
		struct { v3 center; quat rot; quat inv_rot; const v3* verts; int count; } hull;
		struct { v3 p, q; float radius; v3 axis; float inv_axis_len; } cylinder;
	};
} GJK_Shape;

static GJK_Shape gjk_sphere(v3 center, float radius) { return (GJK_Shape){ .type = GJK_POINT, .radius = radius, .point.center = center }; }
static GJK_Shape gjk_capsule(v3 p, v3 q, float radius) { return (GJK_Shape){ .type = GJK_SEGMENT, .radius = radius, .segment.p = p, .segment.q = q }; }
static GJK_Shape gjk_box(v3 center, quat rot, v3 half_extents) { return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.rot = rot, .box.inv_rot = inv(rot), .box.half_extents = half_extents }; }
static GJK_Shape gjk_hull(v3 center, quat rot, const v3* verts, int count) { return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.rot = rot, .hull.inv_rot = inv(rot), .hull.verts = verts, .hull.count = count }; }
static GJK_Shape gjk_cylinder(v3 p, v3 q, float radius) {
	v3 axis = sub(q, p);
	float al = len(axis);
	float inv_al = al > FLT_EPSILON ? 1.0f / al : 0.0f;
	return (GJK_Shape){ .type = GJK_CYLINDER, .cylinder.p = p, .cylinder.q = q, .cylinder.radius = radius, .cylinder.axis = scale(axis, inv_al), .cylinder.inv_axis_len = inv_al };
}

// Hull convenience: pre-scale vertices and build shape.
// Caller must keep scaled_verts alive for the duration of the GJK call.
static GJK_Shape gjk_hull_scaled(const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts)
{
	for (int i = 0; i < hull->vert_count; i++)
		scaled_verts[i] = hmul(hull->verts[i], sc);
	return gjk_hull(pos, rot, scaled_verts, hull->vert_count);
}

// Support function: returns world-space support point + feature ID.
#define gjk_support(shape, dir, out_feat, out_point) do {                                                        \
	const GJK_Shape* sp = (shape); v3 sd = (dir);                                                                \
	switch (sp->type) {                                                                                          \
	case GJK_POINT: *(out_feat) = 0; (out_point) = sp->point.center; break;                                      \
	case GJK_SEGMENT: {                                                                                          \
		int sf = dot(sd, sub(sp->segment.q, sp->segment.p)) >= 0.0f;                                             \
		*(out_feat) = sf; (out_point) = sf ? sp->segment.q : sp->segment.p; break;                               \
	}                                                                                                            \
	case GJK_BOX: {                                                                                              \
		v3 ld = rotate(sp->box.inv_rot, sd); v3 he = sp->box.half_extents;                                       \
		v3 lc = V3(ld.x >= 0 ? he.x : -he.x, ld.y >= 0 ? he.y : -he.y, ld.z >= 0 ? he.z : -he.z);                \
		*(out_feat) = (ld.x >= 0.0f) | ((ld.y >= 0.0f) << 1) | ((ld.z >= 0.0f) << 2);                            \
		(out_point) = add(sp->box.center, rotate(sp->box.rot, lc)); break;                                       \
	}                                                                                                            \
	case GJK_HULL: {                                                                                             \
		v3 ld = rotate(sp->hull.inv_rot, sd);                                                                    \
		float hbest = -1e18f; int hbi = 0;                                                                       \
		for (int hi = 0; hi < sp->hull.count; hi++) {                                                            \
			float hd = dot(sp->hull.verts[hi], ld); if (hd > hbest) { hbest = hd; hbi = hi; }                    \
		}                                                                                                        \
		*(out_feat) = hbi; (out_point) = add(sp->hull.center, rotate(sp->hull.rot, sp->hull.verts[hbi])); break; \
	}                                                                                                            \
	case GJK_CYLINDER: {                                                                                         \
		v3 cu = sp->cylinder.axis;                                                                               \
		if (sp->cylinder.inv_axis_len == 0.0f) { *(out_feat) = 0; (out_point) = sp->cylinder.p; break; }         \
		float cda = dot(sd, cu); int cf = cda >= 0.0f; *(out_feat) = cf;                                         \
		v3 cbase = cf ? sp->cylinder.q : sp->cylinder.p;                                                         \
		v3 cdp = sub(sd, scale(cu, cda)); float cpl = len(cdp);                                                  \
		(out_point) = (cpl > FLT_EPSILON) ? add(cbase, scale(cdp, sp->cylinder.radius / cpl)) : cbase; break;    \
	}                                                                                                            \
	default: *(out_feat) = 0; (out_point) = V3(0,0,0); break;                                                    \
	}                                                                                                            \
} while(0)

static v3 gjk_center(const GJK_Shape* s)
{
	switch (s->type) {
	case GJK_POINT:    return s->point.center;
	case GJK_SEGMENT:  return scale(add(s->segment.p, s->segment.q), 0.5f);
	case GJK_BOX:      return s->box.center;
	case GJK_HULL:     return s->hull.center;
	case GJK_CYLINDER: return scale(add(s->cylinder.p, s->cylinder.q), 0.5f);
	}
	return V3(0,0,0);
}

#define gjk_closest_point(simplex, out) do {                                                                                        \
	const GJK_Simplex* cs = (simplex);                                                                                              \
	float cinv = 1.0f / cs->divisor;                                                                                                \
	switch (cs->count) {                                                                                                            \
	case 1: (out) = cs->v[0].point; break;                                                                                          \
	case 2: (out) = s3_add(s3_scale(cs->v[0].point, cs->v[0].u * cinv), s3_scale(cs->v[1].point, cs->v[1].u * cinv)); break;        \
	case 3: (out) = s3_add(s3_add(s3_scale(cs->v[0].point, cs->v[0].u * cinv), s3_scale(cs->v[1].point, cs->v[1].u * cinv)),        \
	                        s3_scale(cs->v[2].point, cs->v[2].u * cinv)); break;                                                    \
	default: (out) = _mm_setzero_ps(); break;                                                                                       \
	}                                                                                                                               \
} while(0)

static void gjk_witness_points(const GJK_Simplex* s, v3* p1, v3* p2, int* f1, int* f2)
{
	float inv = 1.0f / s->divisor;
	s3 sp1 = _mm_setzero_ps(), sp2 = _mm_setzero_ps();
	float best_u = -1.0f;
	*f1 = 0; *f2 = 0;
	for (int i = 0; i < s->count; i++) {
		float w = s->v[i].u * inv;
		sp1 = s3_add(sp1, s3_scale(s->v[i].point1, w));
		sp2 = s3_add(sp2, s3_scale(s->v[i].point2, w));
		if (s->v[i].u > best_u) { best_u = s->v[i].u; *f1 = s->v[i].feat1; *f2 = s->v[i].feat2; }
	}
	*p1 = s3_store(sp1); *p2 = s3_store(sp2);
}

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.

static int gjk_solve2(GJK_Simplex* s)
{
	s3 a = s->v[0].point, b = s->v[1].point;
	s3 ba = s3_sub(b, a);
	float u = s3_dot(b, ba);
	float v = -s3_dot(a, ba);
	if (v <= 0.0f) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (u <= 0.0f) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	s->divisor = u + v;
	if (s->divisor == 0.0f) return 0;
	s->v[0].u = u; s->v[1].u = v; s->count = 2;
	return 1;
}

static int gjk_solve3(GJK_Simplex* s)
{
	s3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point;
	s3 ba = s3_sub(b, a), cb = s3_sub(c, b), ac = s3_sub(a, c);
	float uAB = s3_dot(b, ba), vAB = -s3_dot(a, ba);
	float uBC = s3_dot(c, cb), vBC = -s3_dot(b, cb);
	float uCA = s3_dot(a, ac), vCA = -s3_dot(c, ac);
	s3 n = s3_cross(ba, s3_sub(c, a));
	float uABC = s3_dot(s3_cross(b, c), n), vABC = s3_dot(s3_cross(c, a), n), wABC = s3_dot(s3_cross(a, b), n);

	if (vAB <= 0.0f && uCA <= 0.0f) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uAB <= 0.0f && vBC <= 0.0f) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uBC <= 0.0f && vCA <= 0.0f) { s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (uAB > 0.0f && vAB > 0.0f && wABC <= 0.0f) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	if (uBC > 0.0f && vBC > 0.0f && uABC <= 0.0f) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	if (uCA > 0.0f && vCA > 0.0f && vABC <= 0.0f) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }
	s->divisor = uABC + vABC + wABC;
	if (s->divisor == 0.0f) return 0;
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->count = 3;
	return 1;
}

#define s3_stp(a, b, c) s3_dot(a, s3_cross(b, c))

static int gjk_solve4(GJK_Simplex* s)
{
	s3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point, d = s->v[3].point;

	s3 ba = s3_sub(b, a), cb = s3_sub(c, b), ac = s3_sub(a, c);
	s3 db = s3_sub(d, b), bd = s3_sub(b, d), dc = s3_sub(d, c), cd = s3_sub(c, d);
	s3 da = s3_sub(d, a), ad = s3_sub(a, d);
	float uAB = s3_dot(b, ba), vAB = -s3_dot(a, ba);
	float uBC = s3_dot(c, cb), vBC = -s3_dot(b, cb);
	float uCA = s3_dot(a, ac), vCA = -s3_dot(c, ac);
	float uBD = s3_dot(d, db), vBD = -s3_dot(b, db);
	float uDC = s3_dot(c, cd), vDC = -s3_dot(d, cd);
	float uAD = s3_dot(d, da), vAD = -s3_dot(a, da);

	s3 n;
	n = s3_cross(da, ba);
	float uADB = s3_dot(s3_cross(d, b), n), vADB = s3_dot(s3_cross(b, a), n), wADB = s3_dot(s3_cross(a, d), n);
	n = s3_cross(s3_sub(c, a), da);
	float uACD = s3_dot(s3_cross(c, d), n), vACD = s3_dot(s3_cross(d, a), n), wACD = s3_dot(s3_cross(a, c), n);
	n = s3_cross(s3_sub(b, c), dc);
	float uCBD = s3_dot(s3_cross(b, d), n), vCBD = s3_dot(s3_cross(d, c), n), wCBD = s3_dot(s3_cross(c, b), n);
	n = s3_cross(ba, s3_sub(c, a));
	float uABC = s3_dot(s3_cross(b, c), n), vABC = s3_dot(s3_cross(c, a), n), wABC = s3_dot(s3_cross(a, b), n);

	float denom = s3_stp(cb, s3_sub(a, b), db);
	if (denom == 0.0f) return 0;
	float vol = 1.0f / denom;
	float uABCD = s3_stp(c, d, b) * vol, vABCD = s3_stp(c, a, d) * vol;
	float wABCD = s3_stp(d, a, b) * vol, xABCD = s3_stp(b, a, c) * vol;

	if (vAB <= 0 && uCA <= 0 && vAD <= 0) { s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uAB <= 0 && vBC <= 0 && vBD <= 0) { s->v[0] = s->v[1]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uBC <= 0 && vCA <= 0 && uDC <= 0) { s->v[0] = s->v[2]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }
	if (uBD <= 0 && vDC <= 0 && uAD <= 0) { s->v[0] = s->v[3]; s->v[0].u = 1; s->divisor = 1; s->count = 1; return 1; }

	if (wABC <= 0 && vADB <= 0 && uAB > 0 && vAB > 0) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	if (uABC <= 0 && wCBD <= 0 && uBC > 0 && vBC > 0) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	if (vABC <= 0 && wACD <= 0 && uCA > 0 && vCA > 0) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }
	if (vCBD <= 0 && uACD <= 0 && uDC > 0 && vDC > 0) { s->v[0] = s->v[3]; s->v[1] = s->v[2]; s->v[0].u = uDC; s->v[1].u = vDC; s->divisor = uDC + vDC; s->count = 2; return 1; }
	if (vACD <= 0 && wADB <= 0 && uAD > 0 && vAD > 0) { s->v[1] = s->v[3]; s->v[0].u = uAD; s->v[1].u = vAD; s->divisor = uAD + vAD; s->count = 2; return 1; }
	if (uCBD <= 0 && uADB <= 0 && uBD > 0 && vBD > 0) { s->v[0] = s->v[1]; s->v[1] = s->v[3]; s->v[0].u = uBD; s->v[1].u = vBD; s->divisor = uBD + vBD; s->count = 2; return 1; }

	if (xABCD <= 0 && uABC > 0 && vABC > 0 && wABC > 0) { s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->divisor = uABC + vABC + wABC; s->count = 3; return 1; }
	if (uABCD <= 0 && uCBD > 0 && vCBD > 0 && wCBD > 0) { s->v[0] = s->v[2]; s->v[2] = s->v[3]; s->v[0].u = uCBD; s->v[1].u = vCBD; s->v[2].u = wCBD; s->divisor = uCBD + vCBD + wCBD; s->count = 3; return 1; }
	if (vABCD <= 0 && uACD > 0 && vACD > 0 && wACD > 0) { s->v[1] = s->v[2]; s->v[2] = s->v[3]; s->v[0].u = uACD; s->v[1].u = vACD; s->v[2].u = wACD; s->divisor = uACD + vACD + wACD; s->count = 3; return 1; }
	if (wABCD <= 0 && uADB > 0 && vADB > 0 && wADB > 0) { s->v[2] = s->v[1]; s->v[1] = s->v[3]; s->v[0].u = uADB; s->v[1].u = vADB; s->v[2].u = wADB; s->divisor = uADB + vADB + wADB; s->count = 3; return 1; }

	s->v[0].u = uABCD; s->v[1].u = vABCD; s->v[2].u = wABCD; s->v[3].u = xABCD;
	s->divisor = 1.0f; s->count = 4;
	return 1;
}

// -----------------------------------------------------------------------------
// Main GJK distance function.
// 3-condition termination + post-hoc radius for sphere/capsule.

#define GJK_MAX_ITERS         64
#define GJK_CONTAINMENT_EPS2  1e-8f
#define GJK_PROGRESS_EPS      1e-7f

static GJK_Result gjk_distance(GJK_Shape shapeA, GJK_Shape shapeB)
{
	GJK_Result result = {0};
	GJK_Simplex simplex = {0};

	v3 init_d = sub(gjk_center(&shapeB), gjk_center(&shapeA));
	if (len2(init_d) < FLT_EPSILON) init_d = V3(1, 0, 0);

	int fA, fB;
	v3 sA; gjk_support(&shapeA, init_d, &fA, sA);
	v3 sB; gjk_support(&shapeB, neg(init_d), &fB, sB);
	simplex.v[0].point1 = s3_load(sA);
	simplex.v[0].point2 = s3_load(sB);
	simplex.v[0].point = s3_sub(s3_load(sB), s3_load(sA));
	simplex.v[0].feat1 = fA;
	simplex.v[0].feat2 = fB;
	simplex.v[0].u = 1.0f;
	simplex.divisor = 1.0f;
	simplex.count = 1;

	float dsq_prev = FLT_MAX;
	int iter = 0;

	while (iter < GJK_MAX_ITERS) {
		int solved = 1;
		if (simplex.count > 1) {
			GJK_Simplex backup;
			memcpy(&backup, &simplex, gjk_simplex_size(simplex.count));
			switch (simplex.count) {
			case 2: solved = gjk_solve2(&simplex); break;
			case 3: solved = gjk_solve3(&simplex); break;
			case 4: solved = gjk_solve4(&simplex); break;
			}
			if (!solved) { memcpy(&simplex, &backup, gjk_simplex_size(backup.count)); break; }
		}
		if (simplex.count == 4) break;

		s3 closest; gjk_closest_point(&simplex, closest);
		float dsq = s3_len2(closest);

		if (dsq <= GJK_CONTAINMENT_EPS2) break;
		if (dsq >= dsq_prev) break;
		dsq_prev = dsq;

		v3 cv = s3_store(closest);
		gjk_support(&shapeA, cv, &fA, sA);
		gjk_support(&shapeB, neg(cv), &fB, sB);
		s3 sw = s3_sub(s3_load(sB), s3_load(sA));

		float max_vert2 = 0.0f;
		for (int i = 0; i < simplex.count; i++) {
			float v2 = s3_len2(simplex.v[i].point);
			if (v2 > max_vert2) max_vert2 = v2;
		}
		float progress = dsq - s3_dot(sw, closest);
		if (progress <= max_vert2 * GJK_PROGRESS_EPS) break;

		iter++;
		GJK_Vertex* vert = &simplex.v[simplex.count];
		vert->point1 = s3_load(sA);
		vert->point2 = s3_load(sB);
		vert->point = sw;
		vert->feat1 = fA;
		vert->feat2 = fB;
		simplex.count++;
	}

	gjk_witness_points(&simplex, &result.point1, &result.point2, &result.feat1, &result.feat2);
	result.distance = len(sub(result.point2, result.point1));
	result.iterations = iter;

	// Post-hoc radius for sphere/capsule core shapes.
	float rA = shapeA.radius;
	float rB = shapeB.radius;
	if (rA != 0.0f || rB != 0.0f) {
		if (result.distance > FLT_EPSILON) {
			v3 normal = scale(sub(result.point2, result.point1), 1.0f / result.distance);
			result.point1 = add(result.point1, scale(normal, rA));
			result.point2 = sub(result.point2, scale(normal, rB));
			result.distance = fmaxf(0.0f, result.distance - rA - rB);
		} else {
			v3 diff = sub(result.point2, result.point1);
			float d2 = len2(diff);
			v3 normal = d2 > FLT_EPSILON * FLT_EPSILON ? scale(diff, 1.0f / sqrtf(d2)) : V3(1, 0, 0);
			result.point1 = add(result.point1, scale(normal, rA));
			result.point2 = sub(result.point2, scale(normal, rB));
			result.distance = 0.0f;
		}
	}

	return result;
}
