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
// Simplex types.
typedef struct GJK_Vertex
{
	v3 point1;  // support point on shape A (world space)
	v3 point2;  // support point on shape B (world space)
	v3 point;   // Minkowski difference: point2 - point1
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
typedef struct GJK_Cache
{
	v3 dir;       // separating direction from previous call
	int hintA;    // hull vertex hint for shape A
	int hintB;    // hull vertex hint for shape B
	int count;    // cached simplex vertex count (0 = direction-only warm start)
	int feat1[2]; // feature IDs on shape A (max 2 cached vertices)
	int feat2[2]; // feature IDs on shape B
	v3 point1[2]; // cached world-space support points on A
	v3 point2[2]; // cached world-space support points on B
} GJK_Cache;
// -----------------------------------------------------------------------------
// Shape types and constructors.
enum { GJK_POINT, GJK_SEGMENT, GJK_BOX, GJK_HULL, GJK_CYLINDER, GJK_TRIANGLE };
typedef struct GJK_Shape
{
	int type;
	float radius; // applied post-hoc (sphere/capsule); 0 for surface shapes
	union {
		struct { v3 center; } point;
		struct { v3 p, q; } segment;
		struct { v3 center; v3 col0, col1, col2; v3 half_extents; } box;
		struct { v3 center; v3 col0, col1, col2; v3 scale; const v3* verts; const float* soa; const HalfEdge* edges; const int* vert_edge; int count; int hint; } hull;
		struct { v3 mid; v3 half_axis; float radius; v3 axis; float inv_axis_len; } cylinder;
		struct { v3 a, b, c; } tri;
	};
} GJK_Shape;
static GJK_Shape gjk_sphere(v3 center, float radius) { return (GJK_Shape){ .type = GJK_POINT, .radius = radius, .point.center = center }; }
static GJK_Shape gjk_capsule(v3 p, v3 q, float radius) { return (GJK_Shape){ .type = GJK_SEGMENT, .radius = radius, .segment.p = p, .segment.q = q }; }
static GJK_Shape gjk_triangle(v3 a, v3 b, v3 c) { return (GJK_Shape){ .type = GJK_TRIANGLE, .tri.a = a, .tri.b = b, .tri.c = c }; }
// Mat3x3 transpose rotate: R^T * v (column-wise multiply-add). Used for inverse rotation.
// 3 broadcasts + 3 muls + 2 adds = 8 SSE ops (vs 18 for row-dot version).
static inline v3 gjk_mat_rotate_t(v3 r0, v3 r1, v3 r2, v3 v) {
	return (v3){ .m = _mm_add_ps(_mm_add_ps(
		_mm_mul_ps(r0.m, _mm_shuffle_ps(v.m, v.m, 0x00)),
		_mm_mul_ps(r1.m, _mm_shuffle_ps(v.m, v.m, 0x55))),
		_mm_mul_ps(r2.m, _mm_shuffle_ps(v.m, v.m, 0xAA))) };
}

// Mat3x3 forward rotate: R * v. Transpose rows to columns, then broadcast-mul.
static inline v3 gjk_mat_rotate(v3 r0, v3 r1, v3 r2, v3 v) {
	__m128 t01lo = _mm_unpacklo_ps(r0.m, r1.m);
	__m128 t01hi = _mm_unpackhi_ps(r0.m, r1.m);
	__m128 t2lo  = _mm_unpacklo_ps(r2.m, _mm_setzero_ps());
	__m128 t2hi  = _mm_unpackhi_ps(r2.m, _mm_setzero_ps());
	__m128 col0 = _mm_movelh_ps(t01lo, t2lo);
	__m128 col1 = _mm_movehl_ps(t2lo, t01lo);
	__m128 col2 = _mm_movelh_ps(t01hi, t2hi);
	return (v3){ .m = _mm_add_ps(_mm_add_ps(
		_mm_mul_ps(col0, _mm_shuffle_ps(v.m, v.m, 0x00)),
		_mm_mul_ps(col1, _mm_shuffle_ps(v.m, v.m, 0x55))),
		_mm_mul_ps(col2, _mm_shuffle_ps(v.m, v.m, 0xAA))) };
}
static GJK_Shape gjk_box(v3 center, quat rot, v3 half_extents) {
	// Columns of rotation matrix = rotated basis vectors
	v3 c0 = quat_rotate(rot, V3(1,0,0)), c1 = quat_rotate(rot, V3(0,1,0)), c2 = quat_rotate(rot, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.col0 = c0, .box.col1 = c1, .box.col2 = c2, .box.half_extents = half_extents };
}
static GJK_Shape gjk_box_m(v3 center, v3 col0, v3 col1, v3 col2, v3 half_extents) {
	return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.col0 = col0, .box.col1 = col1, .box.col2 = col2, .box.half_extents = half_extents };
}
static GJK_Shape gjk_hull(v3 center, quat rot, v3 sc, const v3* verts, int count, const float* soa, const HalfEdge* edges, const int* vert_edge) {
	v3 c0 = quat_rotate(rot, V3(1,0,0)), c1 = quat_rotate(rot, V3(0,1,0)), c2 = quat_rotate(rot, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.col0 = c0, .hull.col1 = c1, .hull.col2 = c2, .hull.scale = sc, .hull.verts = verts, .hull.soa = soa, .hull.edges = edges, .hull.vert_edge = vert_edge, .hull.count = count };
}
static GJK_Shape gjk_hull_m(v3 center, v3 col0, v3 col1, v3 col2, v3 sc, const v3* verts, int count, const float* soa, const HalfEdge* edges, const int* vert_edge) {
	return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.col0 = col0, .hull.col1 = col1, .hull.col2 = col2, .hull.scale = sc, .hull.verts = verts, .hull.soa = soa, .hull.edges = edges, .hull.vert_edge = vert_edge, .hull.count = count };
}
static GJK_Shape gjk_cylinder(v3 p, v3 q, float radius) {
	v3 axis = sub(q, p);
	float al = len(axis);
	float inv_al = al > FLT_EPSILON ? 1.0f / al : 0.0f;
	v3 mid = scale(add(p, q), 0.5f);
	v3 half_axis = scale(axis, 0.5f);
	return (GJK_Shape){ .type = GJK_CYLINDER, .cylinder.mid = mid, .cylinder.half_axis = half_axis, .cylinder.radius = radius, .cylinder.axis = scale(axis, inv_al), .cylinder.inv_axis_len = inv_al };
}
// Hull convenience: pre-scale vertices and build shape.
// Caller must keep scaled_verts alive for the duration of the GJK call.
// For hulls > 32 verts, also builds SoA layout in soa_buf for fast SIMD scan.
// soa_buf must hold 3*count floats (or NULL for small hulls).
// vert_edge_buf must hold hull->vert_count ints. Builds per-vertex first-edge lookup.
static GJK_Shape gjk_hull_scaled(const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts, float* soa_buf)
{
	int n = hull->vert_count;
	// Store raw vertices (not scaled) — scale applied in support function via dot(v*sc, ld) = dot(v, sc*ld)
	const v3* raw_verts = hull->verts;
	(void)scaled_verts; // unused now — keeping parameter for API compat
	const float* soa = NULL;
	// Cache per-vertex first-edge lookup: topology-only, rebuild when hull pointer changes.
	#define VE_CACHE_SLOTS 4
	static const Hull* ve_cache_hull[VE_CACHE_SLOTS] = {0};
	static int ve_cache_buf[VE_CACHE_SLOTS][1024];
	if (hull->edges && n <= 1024) {
		int slot = -1;
		for (int s = 0; s < VE_CACHE_SLOTS; s++) if (ve_cache_hull[s] == hull) { slot = s; break; }
		if (slot < 0) {
			// Evict oldest (round-robin)
			static int ve_next_slot = 0;
			slot = ve_next_slot; ve_next_slot = (ve_next_slot + 1) % VE_CACHE_SLOTS;
			for (int i = 0; i < n; i++) ve_cache_buf[slot][i] = -1;
			for (int i = 0; i < hull->edge_count; i++)
				if (ve_cache_buf[slot][hull->edges[i].origin] < 0) ve_cache_buf[slot][hull->edges[i].origin] = i;
			ve_cache_hull[slot] = hull;
		}
		return gjk_hull(pos, rot, sc, raw_verts, n, soa, hull->edges, ve_cache_buf[slot]);
	}
	return gjk_hull(pos, rot, sc, raw_verts, n, soa, hull->edges, NULL);
}
// Hull support scan: extracted to a function to reduce macro expansion size,
// allowing the compiler to inline quat_rotate in the caller.
// Hill-climbing support: walk vertex adjacency graph from a start vertex.
// O(sqrt(n)) expected for convex hulls vs O(n) for linear scan.
// Hill-climbing support: walk vertex adjacency graph from vertex 0.
// Ring walk: next outgoing edge from v = edges[current ^ 1].next
static int gjk_hull_support_climb(const v3* __restrict verts, const HalfEdge* __restrict edges, const int* __restrict vert_edge, v3 ld, int start)
{
	int best = start;
	float best_d = dot(verts[start], ld);
	for (;;) {
		int e = vert_edge[best];
		if (e < 0) break;
		int start = e;
		int improved = 0;
		do {
			int neighbor = edges[e ^ 1].origin;
			float nd = dot(verts[neighbor], ld);
			if (nd > best_d) { best_d = nd; best = neighbor; improved = 1; }
			e = edges[e ^ 1].next; // walk ring: twin's next gives next outgoing from same vertex
		} while (e != start);
		if (!improved) break;
	}
	return best;
}

static int gjk_hull_support_scan(const v3* __restrict verts, int count, const float* __restrict soa, v3 ld)
{
	__m128 ldx = _mm_shuffle_ps(ld.m, ld.m, 0x00);
	__m128 ldy = _mm_shuffle_ps(ld.m, ld.m, 0x55);
	__m128 ldz = _mm_shuffle_ps(ld.m, ld.m, 0xAA);
	__m128 vbest = _mm_set1_ps(-1e18f);
	__m128i ibest = _mm_setzero_si128();
	int hi = 0;
	if (soa) {
		const float* sx = soa, *sy = sx + count, *sz = sy + count;
		__m128 vbest2 = _mm_set1_ps(-1e18f);
		__m128i ibest2 = _mm_setzero_si128();
		__m128i idx0 = _mm_set_epi32(3,2,1,0), idx1 = _mm_set_epi32(7,6,5,4);
		__m128i eight = _mm_set1_epi32(8);
		for (; hi + 7 < count; hi += 8) {
			__m128 d0 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_loadu_ps(sx+hi), ldx), _mm_mul_ps(_mm_loadu_ps(sy+hi), ldy)), _mm_mul_ps(_mm_loadu_ps(sz+hi), ldz));
			__m128 d1 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_loadu_ps(sx+hi+4), ldx), _mm_mul_ps(_mm_loadu_ps(sy+hi+4), ldy)), _mm_mul_ps(_mm_loadu_ps(sz+hi+4), ldz));
			__m128 m0 = _mm_cmpgt_ps(d0, vbest), m1 = _mm_cmpgt_ps(d1, vbest2);
			vbest = _mm_blendv_ps(vbest, d0, m0);
			ibest = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(ibest), _mm_castsi128_ps(idx0), m0));
			vbest2 = _mm_blendv_ps(vbest2, d1, m1);
			ibest2 = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(ibest2), _mm_castsi128_ps(idx1), m1));
			idx0 = _mm_add_epi32(idx0, eight); idx1 = _mm_add_epi32(idx1, eight);
		}
		__m128 mg = _mm_cmpgt_ps(vbest2, vbest);
		vbest = _mm_blendv_ps(vbest, vbest2, mg);
		ibest = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(ibest), _mm_castsi128_ps(ibest2), mg));
		for (; hi + 3 < count; hi += 4) {
			__m128 dots = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_loadu_ps(sx+hi), ldx), _mm_mul_ps(_mm_loadu_ps(sy+hi), ldy)), _mm_mul_ps(_mm_loadu_ps(sz+hi), ldz));
			__m128i idx = _mm_set_epi32(hi+3,hi+2,hi+1,hi);
			__m128 mask = _mm_cmpgt_ps(dots, vbest);
			vbest = _mm_blendv_ps(vbest, dots, mask);
			ibest = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(ibest), _mm_castsi128_ps(idx), mask));
		}
	} else {
		for (; hi + 3 < count; hi += 4) {
			__m128 v0 = verts[hi].m, v1 = verts[hi+1].m, v2 = verts[hi+2].m, v3r = verts[hi+3].m;
			__m128 t0 = _mm_unpacklo_ps(v0, v1), t1 = _mm_unpacklo_ps(v2, v3r);
			__m128 t2 = _mm_unpackhi_ps(v0, v1), t3 = _mm_unpackhi_ps(v2, v3r);
			__m128 dots = _mm_add_ps(_mm_add_ps(_mm_mul_ps(_mm_movelh_ps(t0, t1), ldx), _mm_mul_ps(_mm_movehl_ps(t1, t0), ldy)), _mm_mul_ps(_mm_movelh_ps(t2, t3), ldz));
			__m128i idx = _mm_set_epi32(hi+3, hi+2, hi+1, hi);
			__m128 mask = _mm_cmpgt_ps(dots, vbest);
			vbest = _mm_blendv_ps(vbest, dots, mask);
			ibest = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(ibest), _mm_castsi128_ps(idx), mask));
		}
	}
	float hbests[4]; int hbidxs[4];
	_mm_storeu_ps(hbests, vbest); _mm_storeu_si128((__m128i*)hbidxs, ibest);
	float hbest = hbests[0]; int hbi = hbidxs[0];
	for (int k = 1; k < 4; k++) { if (hbests[k] > hbest) { hbest = hbests[k]; hbi = hbidxs[k]; } }
	for (; hi < count; hi++) { float hd = dot(verts[hi], ld); if (hd > hbest) { hbest = hd; hbi = hi; } }
	return hbi;
}
// Box/cylinder support: separate functions to keep gjk_support macro small for inlining budget.
static v3 gjk_box_support(const GJK_Shape* __restrict sp, v3 sd, int* __restrict feat)
{
	v3 ld = gjk_mat_rotate(sp->box.col0, sp->box.col1, sp->box.col2, sd);
	// SSE sign select: copysign(he, ld) — flip he components where ld < 0
	__m128 sign_bits = _mm_and_ps(ld.m, _mm_castsi128_ps(_mm_set1_epi32((int)0x80000000)));
	v3 lc = { .m = _mm_xor_ps(sp->box.half_extents.m, sign_bits) };
	*feat = _mm_movemask_ps(_mm_cmpge_ps(ld.m, _mm_setzero_ps())) & 7;
	return add(sp->box.center, gjk_mat_rotate_t(sp->box.col0, sp->box.col1, sp->box.col2, lc));
}
static v3 gjk_cylinder_support(const GJK_Shape* __restrict sp, v3 sd, int* __restrict feat)
{
	v3 cu = sp->cylinder.axis;
	if (sp->cylinder.inv_axis_len == 0.0f) { *feat = 0; return sp->cylinder.mid; }
	float cda = dot(sd, cu);
	int cap = cda >= 0.0f;
	// Branchless base: mid + copysign(half_axis, da)
	__m128 sign = _mm_and_ps(_mm_set1_ps(cda), _mm_castsi128_ps(_mm_set1_epi32((int)0x80000000)));
	v3 cbase = add(sp->cylinder.mid, (v3){ .m = _mm_xor_ps(sp->cylinder.half_axis.m, sign) });
	v3 cdp = sub(sd, scale(cu, cda));
	__m128 cdp2 = v3_dot_m(cdp, cdp);
	__m128 cpl = _mm_sqrt_ss(cdp2);
	float cplf = _mm_cvtss_f32(cpl);
	if (cplf <= FLT_EPSILON) { *feat = cap; return cbase; }
	float inv_cpl = 1.0f / cplf;
	// Pack normalized perp direction as 3x10-bit signed + 1 cap bit = 31 bits
	int ix = (int)(cdp.x * inv_cpl * 511.0f), iy = (int)(cdp.y * inv_cpl * 511.0f), iz = (int)(cdp.z * inv_cpl * 511.0f);
	*feat = cap | ((ix & 0x3FF) << 1) | ((iy & 0x3FF) << 11) | ((iz & 0x3FF) << 21);
	return add(cbase, v3_scale_m(cdp, _mm_div_ss(_mm_set_ss(sp->cylinder.radius), cpl)));
}
// Support macro: dispatches per shape type. Box/cylinder/hull-scan are functions to reduce code size.
#define gjk_support(shape, dir, out_feat, out_point) do {                                                                 \
	GJK_Shape* sp = (shape); v3 sd = (dir);                                                                               \
	switch (sp->type) {                                                                                                   \
	case GJK_POINT: *(out_feat) = 0; (out_point) = sp->point.center; break;                                               \
	case GJK_SEGMENT: {                                                                                                   \
		int sf = dot(sd, sub(sp->segment.q, sp->segment.p)) >= 0.0f;                                                      \
		*(out_feat) = sf; (out_point) = sf ? sp->segment.q : sp->segment.p; break;                                         \
	}                                                                                                                     \
	case GJK_BOX: (out_point) = gjk_box_support(sp, sd, (out_feat)); break;                                                \
	case GJK_HULL: {                                                                                                      \
		v3 ld = gjk_mat_rotate(sp->hull.col0, sp->hull.col1, sp->hull.col2, sd);                                          \
		v3 sld = hmul(ld, sp->hull.scale); /* scale-weighted direction for dot(v, sld) = dot(v*sc, ld) */                  \
		int hbi = (sp->hull.vert_edge && sp->hull.count > 0)                                                              \
			? gjk_hull_support_climb(sp->hull.verts, sp->hull.edges, sp->hull.vert_edge, sld, sp->hull.hint)               \
			: gjk_hull_support_scan(sp->hull.verts, sp->hull.count, sp->hull.soa, sld);                                   \
		sp->hull.hint = hbi;                                                                                               \
		*(out_feat) = hbi;                                                                                                \
		(out_point) = add(sp->hull.center, gjk_mat_rotate_t(sp->hull.col0, sp->hull.col1, sp->hull.col2,                  \
		                  hmul(sp->hull.verts[hbi], sp->hull.scale)));                                                     \
		break;                                                                                                            \
	}                                                                                                                     \
	case GJK_CYLINDER: (out_point) = gjk_cylinder_support(sp, sd, (out_feat)); break;                                      \
	case GJK_TRIANGLE: {                                                                                                  \
		/* 3 dots in parallel via AoS→SoA transpose */                                                                     \
		__m128 v0 = sp->tri.a.m, v1 = sp->tri.b.m, v2 = sp->tri.c.m;                                                     \
		__m128 t01lo = _mm_unpacklo_ps(v0, v1), t01hi = _mm_unpackhi_ps(v0, v1);                                          \
		__m128 t2lo = _mm_unpacklo_ps(v2, _mm_setzero_ps());                                                              \
		__m128 t2hi = _mm_unpackhi_ps(v2, _mm_setzero_ps());                                                              \
		__m128 xs = _mm_movelh_ps(t01lo, t2lo), ys = _mm_movehl_ps(t2lo, t01lo), zs = _mm_movelh_ps(t01hi, t2hi);         \
		__m128 dx = _mm_shuffle_ps(sd.m, sd.m, 0x00), dy = _mm_shuffle_ps(sd.m, sd.m, 0x55);                              \
		__m128 dz = _mm_shuffle_ps(sd.m, sd.m, 0xAA);                                                                     \
		__m128 dots = _mm_add_ps(_mm_add_ps(_mm_mul_ps(xs, dx), _mm_mul_ps(ys, dy)), _mm_mul_ps(zs, dz));                 \
		/* dots = {da, db, dc, 0} — find max */                                                                            \
		float td[4]; _mm_storeu_ps(td, dots);                                                                              \
		if (td[0] >= td[1] && td[0] >= td[2]) { *(out_feat) = 0; (out_point) = sp->tri.a; }                               \
		else if (td[1] >= td[2]) { *(out_feat) = 1; (out_point) = sp->tri.b; }                                            \
		else { *(out_feat) = 2; (out_point) = sp->tri.c; }                                                                \
		break;                                                                                                            \
	}                                                                                                                     \
	default: *(out_feat) = 0; (out_point) = V3(0,0,0); break;                                                             \
	}                                                                                                                     \
} while(0)

// Reconstruct a support point from a cached feature index (no direction search needed).
static inline v3 gjk_support_feature(const GJK_Shape* sp, int feat)
{
	switch (sp->type) {
	case GJK_POINT: return sp->point.center;
	case GJK_SEGMENT: return feat ? sp->segment.q : sp->segment.p;
	case GJK_BOX: {
		v3 he = sp->box.half_extents;
		v3 lc = V3((feat & 1) ? he.x : -he.x, (feat & 2) ? he.y : -he.y, (feat & 4) ? he.z : -he.z);
		return add(sp->box.center, gjk_mat_rotate_t(sp->box.col0, sp->box.col1, sp->box.col2, lc));
	}
	case GJK_HULL:
		return add(sp->hull.center, gjk_mat_rotate_t(sp->hull.col0, sp->hull.col1, sp->hull.col2, hmul(sp->hull.verts[feat], sp->hull.scale)));
	case GJK_CYLINDER: {
		int cap = feat & 1;
		__m128 sign = cap ? _mm_setzero_ps() : _mm_castsi128_ps(_mm_set1_epi32((int)0x80000000));
		v3 cbase = add(sp->cylinder.mid, (v3){ .m = _mm_xor_ps(sp->cylinder.half_axis.m, sign) });
		// Decode 3x10-bit signed perpendicular direction
		int ix = (feat >> 1) & 0x3FF, iy = (feat >> 11) & 0x3FF, iz = (feat >> 21) & 0x3FF;
		if (ix >= 0x200) ix -= 0x400; // sign extend 10-bit
		if (iy >= 0x200) iy -= 0x400;
		if (iz >= 0x200) iz -= 0x400;
		if ((ix | iy | iz) == 0) return cbase;
		// Direction is approximately unit-length after quantization; skip renormalization.
		return add(cbase, scale(V3((float)ix, (float)iy, (float)iz), sp->cylinder.radius / 511.0f));
	}
	case GJK_TRIANGLE:
		return feat == 0 ? sp->tri.a : (feat == 1 ? sp->tri.b : sp->tri.c);
	}
	return V3(0,0,0);
}

static v3 gjk_center(const GJK_Shape* s)
{
	switch (s->type) {
	case GJK_POINT:    return s->point.center;
	case GJK_SEGMENT:  return scale(add(s->segment.p, s->segment.q), 0.5f);
	case GJK_BOX:      return s->box.center;
	case GJK_HULL:     return s->hull.center;
	case GJK_CYLINDER: return s->cylinder.mid;
	case GJK_TRIANGLE: return scale(add(add(s->tri.a, s->tri.b), s->tri.c), 1.0f/3.0f);
	}
	return V3(0,0,0);
}
#define gjk_closest_point(simplex, out) do {                                                                        \
	const GJK_Simplex* cs = (simplex);                                                                              \
	float cinv = 1.0f / cs->divisor;                                                                                \
	switch (cs->count) {                                                                                            \
	case 1: (out) = cs->v[0].point; break;                                                                          \
	case 2: (out) = add(scale(cs->v[0].point, cs->v[0].u * cinv), scale(cs->v[1].point, cs->v[1].u * cinv)); break; \
	case 3: (out) = add(add(scale(cs->v[0].point, cs->v[0].u * cinv), scale(cs->v[1].point, cs->v[1].u * cinv)),    \
	                     scale(cs->v[2].point, cs->v[2].u * cinv)); break;                                          \
	default: (out) = V3(0,0,0); break;                                                                              \
	}                                                                                                               \
} while(0)
#define gjk_witness_points(simplex, out_p1, out_p2, out_f1, out_f2) do {                                              \
	const GJK_Simplex* ws = (simplex);                                                                                \
	float winv = 1.0f / ws->divisor;                                                                                  \
	switch (ws->count) {                                                                                              \
	case 1:                                                                                                           \
		(out_p1) = ws->v[0].point1; (out_p2) = ws->v[0].point2;                                                       \
		*(out_f1) = ws->v[0].feat1; *(out_f2) = ws->v[0].feat2; break;                                                \
	case 2: {                                                                                                         \
		float w0 = ws->v[0].u * winv, w1 = ws->v[1].u * winv;                                                         \
		(out_p1) = add(scale(ws->v[0].point1, w0), scale(ws->v[1].point1, w1));                                        \
		(out_p2) = add(scale(ws->v[0].point2, w0), scale(ws->v[1].point2, w1));                                        \
		int wi = ws->v[1].u > ws->v[0].u;                                                                             \
		*(out_f1) = ws->v[wi].feat1; *(out_f2) = ws->v[wi].feat2; break;                                              \
	}                                                                                                                 \
	case 3: {                                                                                                         \
		float w0 = ws->v[0].u * winv, w1 = ws->v[1].u * winv, w2 = ws->v[2].u * winv;                                \
		(out_p1) = add(add(scale(ws->v[0].point1, w0), scale(ws->v[1].point1, w1)), scale(ws->v[2].point1, w2));       \
		(out_p2) = add(add(scale(ws->v[0].point2, w0), scale(ws->v[1].point2, w1)), scale(ws->v[2].point2, w2));       \
		int wi = 0; if (ws->v[1].u > ws->v[wi].u) wi = 1; if (ws->v[2].u > ws->v[wi].u) wi = 2;                      \
		*(out_f1) = ws->v[wi].feat1; *(out_f2) = ws->v[wi].feat2; break;                                              \
	}                                                                                                                 \
	default: (out_p1) = V3(0,0,0); (out_p2) = V3(0,0,0); *(out_f1) = 0; *(out_f2) = 0; break;                        \
	}                                                                                                                 \
} while(0)
// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.
static __declspec(noinline) int gjk_solve2(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point;
	v3 ba = sub(b, a);
	float u = dot(b, ba);
	float v = -dot(a, ba);
	// Branchless: compute edge case always, then conditionally override for vertex cases
	float div = u + v;
	if (div == 0.0f) return 0;
	int va = v <= 0.0f, vb = u <= 0.0f;
	if (va | vb) {
		// One of the vertex regions. If vb, swap v[0]=v[1] first.
		if (vb) s->v[0] = s->v[1];
		s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1;
		return 1;
	}
	s->v[0].u = u; s->v[1].u = v; s->divisor = div; s->count = 2;
	return 1;
}
static __declspec(noinline) int gjk_solve3(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point;
	v3 ba = sub(b, a), cb = sub(c, b), ac = sub(a, c);
	float uAB = dot(b, ba), vAB = -dot(a, ba);
	float uBC = dot(c, cb), vBC = -dot(b, cb);
	float uCA = dot(a, ac), vCA = -dot(c, ac);
	v3 n = cross(ba, sub(c, a));
	float uABC = dot(cross(b, c), n), vABC = dot(cross(c, a), n), wABC = dot(cross(a, b), n);
	// Pack sign bits: bit=1 means value > 0
	int signs = (uAB > 0) | ((vAB > 0) << 1) | ((uBC > 0) << 2) | ((vBC > 0) << 3) |
	            ((uCA > 0) << 4) | ((vCA > 0) << 5) | ((uABC > 0) << 6) | ((vABC > 0) << 7) | ((wABC > 0) << 8);
	// Vertex A: vAB<=0 && uCA<=0 → bits 1,4 both 0
	if ((signs & 0x12) == 0) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	// Vertex B: uAB<=0 && vBC<=0 → bits 0,3 both 0
	if ((signs & 0x09) == 0) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	// Vertex C: uBC<=0 && vCA<=0 → bits 2,5 both 0
	if ((signs & 0x24) == 0) { s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	// Edge AB: uAB>0 && vAB>0 && wABC<=0 → bits 0,1 set, bit 8 clear
	if ((signs & 0x103) == 0x03) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	// Edge BC: uBC>0 && vBC>0 && uABC<=0 → bits 2,3 set, bit 6 clear
	if ((signs & 0x4C) == 0x0C) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	// Edge CA: uCA>0 && vCA>0 && vABC<=0 → bits 4,5 set, bit 7 clear
	if ((signs & 0xB0) == 0x30) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }
	// Face ABC
	float divABC = uABC + vABC + wABC;
	if (divABC == 0.0f) return 0;
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->divisor = divABC; s->count = 3;
	return 1;
}
#define stp(a, b, c) dot(a, cross(b, c))
static __declspec(noinline) int gjk_solve4(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point, d = s->v[3].point;
	float uAB = dot(b, sub(b, a)), vAB = dot(a, sub(a, b));
	float uBC = dot(c, sub(c, b)), vBC = dot(b, sub(b, c));
	float uCA = dot(a, sub(a, c)), vCA = dot(c, sub(c, a));
	float uBD = dot(d, sub(d, b)), vBD = dot(b, sub(b, d));
	float uDC = dot(c, sub(c, d)), vDC = dot(d, sub(d, c));
	float uAD = dot(d, sub(d, a)), vAD = dot(a, sub(a, d));
	v3 n;
	n = cross(sub(d, a), sub(b, a));
	float uADB = dot(cross(d, b), n), vADB = dot(cross(b, a), n), wADB = dot(cross(a, d), n);
	n = cross(sub(c, a), sub(d, a));
	float uACD = dot(cross(c, d), n), vACD = dot(cross(d, a), n), wACD = dot(cross(a, c), n);
	n = cross(sub(b, c), sub(d, c));
	float uCBD = dot(cross(b, d), n), vCBD = dot(cross(d, c), n), wCBD = dot(cross(c, b), n);
	n = cross(sub(b, a), sub(c, a));
	float uABC = dot(cross(b, c), n), vABC = dot(cross(c, a), n), wABC = dot(cross(a, b), n);
	float denom = stp(sub(c, b), sub(a, b), sub(d, b));
	if (denom == 0.0f) return 0;
	float vol = 1.0f / denom;
	float uABCD = stp(c, d, b) * vol, vABCD = stp(c, a, d) * vol;
	float wABCD = stp(d, a, b) * vol, xABCD = stp(b, a, c) * vol;
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
#define GJK_PROGRESS_EPS      1e-4f
// cache: optional warm-start from previous frame (separating direction + hull hints).
static GJK_Result gjk_distance(GJK_Shape* __restrict shapeA, GJK_Shape* __restrict shapeB, GJK_Cache* cache)
{
	GJK_Result result;
	GJK_Simplex simplex;
	if (cache) {
		if (shapeA->type == GJK_HULL) shapeA->hull.hint = cache->hintA;
		if (shapeB->type == GJK_HULL) shapeB->hull.hint = cache->hintB;
	}
	int fA, fB;
	v3 sA, sB;
	// Reconstruct simplex from cache: use world-space points (fast) or feature replay (accurate for cylinders).
	int use_simplex_cache = cache && cache->count >= 1 && cache->count <= 2;
	if (use_simplex_cache) {
		int has_cyl = shapeA->type == GJK_CYLINDER || shapeB->type == GJK_CYLINDER;
		for (int i = 0; i < cache->count; i++) {
			v3 p1, p2;
			if (has_cyl) {
				p1 = gjk_support_feature(shapeA, cache->feat1[i]);
				p2 = gjk_support_feature(shapeB, cache->feat2[i]);
			} else {
				p1 = cache->point1[i];
				p2 = cache->point2[i];
			}
			simplex.v[i].point1 = p1;
			simplex.v[i].point2 = p2;
			simplex.v[i].point = sub(p2, p1);
			simplex.v[i].feat1 = cache->feat1[i];
			simplex.v[i].feat2 = cache->feat2[i];
			simplex.v[i].u = 1.0f;
		}
		simplex.divisor = (float)cache->count;
		simplex.count = cache->count;
	} else {
		v3 init_d;
		if (cache && len2(cache->dir) > FLT_EPSILON) {
			init_d = cache->dir;
		} else {
			init_d = sub(gjk_center(shapeB), gjk_center(shapeA));
			if (len2(init_d) < FLT_EPSILON) init_d = V3(1, 0, 0);
		}
		gjk_support(shapeA, init_d, &fA, sA);
		gjk_support(shapeB, neg(init_d), &fB, sB);
		simplex.v[0].point1 = sA;
		simplex.v[0].point2 = sB;
		simplex.v[0].point = sub(sB, sA);
		simplex.v[0].feat1 = fA;
		simplex.v[0].feat2 = fB;
		simplex.v[0].u = 1.0f;
		simplex.divisor = 1.0f;
		simplex.count = 1;
	}
	float dsq_prev = FLT_MAX;
	int use_index_term = shapeA->type != GJK_CYLINDER && shapeB->type != GJK_CYLINDER;
	int iter = 0;
	while (iter < GJK_MAX_ITERS) {
		if (simplex.count > 1) {
			int solved;
			switch (simplex.count) {
			case 2: solved = gjk_solve2(&simplex); break;
			case 3: solved = gjk_solve3(&simplex); break;
			case 4: solved = gjk_solve4(&simplex); break;
			default: solved = 0; break;
			}
			if (!solved) break;
		}
		if (simplex.count == 4) break;
		v3 closest; gjk_closest_point(&simplex, closest);
		float dsq = len2(closest);
		if (dsq <= GJK_CONTAINMENT_EPS2) break;
		if (dsq >= dsq_prev) break;
		dsq_prev = dsq;
		gjk_support(shapeA, closest, &fA, sA);
		gjk_support(shapeB, neg(closest), &fB, sB);
		v3 w = sub(sB, sA);
		if (use_index_term) {
			int dup = 0;
			for (int i = 0; i < simplex.count; i++) if (simplex.v[i].feat1 == fA && simplex.v[i].feat2 == fB) { dup = 1; break; }
			if (dup) break;
		}
		float max_vert2 = 0.0f;
		for (int i = 0; i < simplex.count; i++) {
			float v2 = len2(simplex.v[i].point);
			if (v2 > max_vert2) max_vert2 = v2;
		}
		float progress = dsq - dot(w, closest);
		if (progress <= max_vert2 * GJK_PROGRESS_EPS) break;
		iter++;
		GJK_Vertex* vert = &simplex.v[simplex.count];
		vert->point1 = sA;
		vert->point2 = sB;
		vert->point = w;
		vert->feat1 = fA;
		vert->feat2 = fB;
		simplex.count++;
	}
	gjk_witness_points(&simplex, result.point1, result.point2, &result.feat1, &result.feat2);
	v3 sep = sub(result.point2, result.point1);
	result.distance = len(sep);
	result.iterations = iter;
	if (cache) {
		cache->dir = sep;
		if (shapeA->type == GJK_HULL) cache->hintA = shapeA->hull.hint;
		if (shapeB->type == GJK_HULL) cache->hintB = shapeB->hull.hint;
		// Cache simplex features for next frame warm-start
		int cn = simplex.count < 2 ? simplex.count : 2;
		cache->count = cn;
		for (int i = 0; i < cn; i++) { cache->feat1[i] = simplex.v[i].feat1; cache->feat2[i] = simplex.v[i].feat2; cache->point1[i] = simplex.v[i].point1; cache->point2[i] = simplex.v[i].point2; }
	}
	// Post-hoc radius for sphere/capsule core shapes.
	float rA = shapeA->radius;
	float rB = shapeB->radius;
	if (rA != 0.0f || rB != 0.0f) {
		if (result.distance > FLT_EPSILON) {
			v3 normal = scale(sep, 1.0f / result.distance);
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
// By-value convenience wrapper for callers that construct shapes inline.
// The v3* cache variant is for callers that don't need hull hint caching.
static GJK_Result gjk_distance_v(GJK_Shape a, GJK_Shape b, v3* dir_cache)
{
	GJK_Cache c = {0};
	GJK_Cache* cp = NULL;
	if (dir_cache) { c.dir = *dir_cache; cp = &c; }
	GJK_Result r = gjk_distance(&a, &b, cp);
	if (dir_cache) *dir_cache = c.dir;
	return r;
}
