// See LICENSE for licensing info.
// gjk.c -- GJK distance algorithm for 3D convex shapes.
//
// Per-shape-type support dispatch:
//   - Sphere/capsule: world-space core point/segment, radius applied post-hoc
//   - Box: analytical sign-based corner selection
//   - Hull: local-space vertex scan with rotation
//   - Cylinder: implicit surface support (no tessellation)
//   - Triangle: 3-vertex max-dot
// Termination: containment, monotonic progress, support progress, duplicate vertex.

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
	int feat1[4]; // feature IDs on shape A
	int feat2[4]; // feature IDs on shape B
	v3 point1[4]; // cached world-space support points on A
	v3 point2[4]; // cached world-space support points on B
} GJK_Cache;

// (GJK_Simplex_Out removed: EPA seeds from warm-cache Minkowski directions
// stored in EpaManifold instead, avoiding the CORE-vs-inflated support mismatch
// for sphere/capsule shapes.)

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
		struct { v3 center; v3 col0, col1, col2; v3 scale; const v3* verts; const float* soa; const int* edge_twin; const int* edge_next; const int* edge_origin; const int* vert_edge; int count; int hint; } hull;
		struct { v3 mid; v3 half_axis; float radius; v3 axis; float inv_axis_len; } cylinder;
		struct { v3 a, b, c; } tri;
	};
} GJK_Shape;

static GJK_Shape gjk_sphere(v3 center, float radius) { return (GJK_Shape){ .type = GJK_POINT, .radius = radius, .point.center = center }; }
static GJK_Shape gjk_capsule(v3 p, v3 q, float radius) { return (GJK_Shape){ .type = GJK_SEGMENT, .radius = radius, .segment.p = p, .segment.q = q }; }
static GJK_Shape gjk_triangle(v3 a, v3 b, v3 c) { return (GJK_Shape){ .type = GJK_TRIANGLE, .tri.a = a, .tri.b = b, .tri.c = c }; }

static GJK_Shape gjk_box(v3 center, quat rot, v3 half_extents)
{
	v3 c0 = quat_rotate(rot, V3(1,0,0)), c1 = quat_rotate(rot, V3(0,1,0)), c2 = quat_rotate(rot, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.col0 = c0, .box.col1 = c1, .box.col2 = c2, .box.half_extents = half_extents };
}

static GJK_Shape gjk_box_m(v3 center, v3 col0, v3 col1, v3 col2, v3 half_extents)
{
	return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.col0 = col0, .box.col1 = col1, .box.col2 = col2, .box.half_extents = half_extents };
}

static GJK_Shape gjk_hull(v3 center, quat rot, v3 sc, const v3* verts, int count, const float* soa, const int* edge_twin, const int* edge_next, const int* edge_origin, const int* vert_edge)
{
	v3 c0 = quat_rotate(rot, V3(1,0,0)), c1 = quat_rotate(rot, V3(0,1,0)), c2 = quat_rotate(rot, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.col0 = c0, .hull.col1 = c1, .hull.col2 = c2, .hull.scale = sc, .hull.verts = verts, .hull.soa = soa, .hull.edge_twin = edge_twin, .hull.edge_next = edge_next, .hull.edge_origin = edge_origin, .hull.vert_edge = vert_edge, .hull.count = count };
}

static GJK_Shape gjk_hull_m(v3 center, v3 col0, v3 col1, v3 col2, v3 sc, const v3* verts, int count, const float* soa, const int* edge_twin, const int* edge_next, const int* edge_origin, const int* vert_edge)
{
	return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.col0 = col0, .hull.col1 = col1, .hull.col2 = col2, .hull.scale = sc, .hull.verts = verts, .hull.soa = soa, .hull.edge_twin = edge_twin, .hull.edge_next = edge_next, .hull.edge_origin = edge_origin, .hull.vert_edge = vert_edge, .hull.count = count };
}

static GJK_Shape gjk_cylinder(v3 p, v3 q, float radius)
{
	v3 axis = sub(q, p);
	float al = len(axis);
	float inv_al = al > FLT_EPSILON ? 1.0f / al : 0.0f;
	v3 mid = scale(add(p, q), 0.5f);
	v3 half_axis = scale(axis, 0.5f);
	return (GJK_Shape){ .type = GJK_CYLINDER, .cylinder.mid = mid, .cylinder.half_axis = half_axis, .cylinder.radius = radius, .cylinder.axis = scale(axis, inv_al), .cylinder.inv_axis_len = inv_al };
}

// -----------------------------------------------------------------------------
// Rotation helpers.

#define mat3_tmul_v mat3_tmul_v
#define mat3_mul_v   mat3_mul_v

// -----------------------------------------------------------------------------
// Hull convenience constructor with vert-edge table caching.

static GJK_Shape gjk_hull_scaled(const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts, float* soa_buf)
{
	int n = hull->vert_count;
	const v3* raw_verts = hull->verts;
	(void)scaled_verts;
	(void)soa_buf;
	const float* soa = hull->soa_verts;

	// Cache per-vertex first-edge lookup: topology-only, rebuild when hull pointer changes.
	#define VE_CACHE_SLOTS 4
	static const Hull* ve_cache_hull[VE_CACHE_SLOTS] = {0};
	static int ve_cache_buf[VE_CACHE_SLOTS][1024];
	if (hull->edge_twin && n <= 1024) {
		int slot = -1;
		for (int s = 0; s < VE_CACHE_SLOTS; s++) if (ve_cache_hull[s] == hull) { slot = s; break; }
		if (slot < 0) {
			static int ve_next_slot = 0;
			slot = ve_next_slot; ve_next_slot = (ve_next_slot + 1) % VE_CACHE_SLOTS;
			for (int i = 0; i < n; i++) ve_cache_buf[slot][i] = -1;
			for (int i = 0; i < hull->edge_count; i++) {
				if (ve_cache_buf[slot][hull->edge_origin[i]] < 0) ve_cache_buf[slot][hull->edge_origin[i]] = i;
			}
			ve_cache_hull[slot] = hull;
		}
		return gjk_hull(pos, rot, sc, raw_verts, n, soa, hull->edge_twin, hull->edge_next, hull->edge_origin, ve_cache_buf[slot]);
	}
	return gjk_hull(pos, rot, sc, raw_verts, n, soa, NULL, NULL, NULL, NULL);
}

// -----------------------------------------------------------------------------
// Hull support functions.

// Hill-climbing support: walk vertex adjacency graph from a start vertex.
// O(sqrt(n)) expected for convex hulls vs O(n) for linear scan.
static int gjk_hull_support_climb(const v3* __restrict verts, const int* __restrict edge_twin, const int* __restrict edge_next, const int* __restrict edge_origin, const int* __restrict vert_edge, v3 ld, int start)
{
	int best = start;
	float best_d = dot(verts[start], ld);
	for (;;) {
		int e = vert_edge[best];
		if (e < 0) break;
		int start = e;
		int improved = 0;
		do {
			int twin = edge_twin[e];
			int neighbor = edge_origin[twin];
			float nd = dot(verts[neighbor], ld);
			if (nd > best_d) { best_d = nd; best = neighbor; improved = 1; }
			e = edge_next[twin];
		} while (e != start);
		if (!improved) break;
	}
	return best;
}

// Linear scan with 4/8-wide SIMD dot products.
static int gjk_hull_support_scan(const v3* __restrict verts, int count, const float* __restrict soa, v3 ld)
{
	simd4f ldx = simd_splat(ld.m, 0);
	simd4f ldy = simd_splat(ld.m, 1);
	simd4f ldz = simd_splat(ld.m, 2);
	simd4f vbest = simd_set1(-1e18f);
	simd4i ibest = simd_set1_i(0);
	int hi = 0;

	if (soa) {
		// SoA path: 8-wide with dual accumulators for ILP.
		const float* sx = soa, *sy = sx + count, *sz = sy + count;
		simd4f vbest2 = simd_set1(-1e18f);
		simd4i ibest2 = simd_set1_i(0);
		simd4i idx0 = simd_seti(0, 1, 2, 3), idx1 = simd_seti(4, 5, 6, 7);
		simd4i eight = simd_set1_i(8);
		for (; hi + 7 < count; hi += 8) {
			simd4f d0 = simd_add(simd_add(simd_mul(simd_load(sx+hi), ldx), simd_mul(simd_load(sy+hi), ldy)), simd_mul(simd_load(sz+hi), ldz));
			simd4f d1 = simd_add(simd_add(simd_mul(simd_load(sx+hi+4), ldx), simd_mul(simd_load(sy+hi+4), ldy)), simd_mul(simd_load(sz+hi+4), ldz));
			simd4f m0 = simd_cmpgt(d0, vbest), m1 = simd_cmpgt(d1, vbest2);
			vbest = simd_blendv(vbest, d0, m0);
			ibest = simd_cast_ftoi(simd_blendv(simd_cast_itof(ibest), simd_cast_itof(idx0), m0));
			vbest2 = simd_blendv(vbest2, d1, m1);
			ibest2 = simd_cast_ftoi(simd_blendv(simd_cast_itof(ibest2), simd_cast_itof(idx1), m1));
			idx0 = simd_add_i(idx0, eight); idx1 = simd_add_i(idx1, eight);
		}

		simd4f mg = simd_cmpgt(vbest2, vbest);
		vbest = simd_blendv(vbest, vbest2, mg);
		ibest = simd_cast_ftoi(simd_blendv(simd_cast_itof(ibest), simd_cast_itof(ibest2), mg));
		for (; hi + 3 < count; hi += 4) {
			simd4f dots = simd_add(simd_add(simd_mul(simd_load(sx+hi), ldx), simd_mul(simd_load(sy+hi), ldy)), simd_mul(simd_load(sz+hi), ldz));
			simd4i idx = simd_seti(hi, hi+1, hi+2, hi+3);
			simd4f mask = simd_cmpgt(dots, vbest);
			vbest = simd_blendv(vbest, dots, mask);
			ibest = simd_cast_ftoi(simd_blendv(simd_cast_itof(ibest), simd_cast_itof(idx), mask));
		}
	} else {
		// AoS path: transpose 4 vertices to SoA per iteration.
		for (; hi + 3 < count; hi += 4) {
			simd4f v0 = verts[hi].m, v1 = verts[hi+1].m, v2 = verts[hi+2].m, v3r = verts[hi+3].m;
			simd4f t0 = simd_unpacklo(v0, v1), t1 = simd_unpacklo(v2, v3r);
			simd4f t2 = simd_unpackhi(v0, v1), t3 = simd_unpackhi(v2, v3r);
			simd4f dots = simd_add(simd_add(simd_mul(simd_movelh(t0, t1), ldx), simd_mul(simd_movehl(t1, t0), ldy)), simd_mul(simd_movelh(t2, t3), ldz));
			simd4i idx = simd_seti(hi, hi+1, hi+2, hi+3);
			simd4f mask = simd_cmpgt(dots, vbest);
			vbest = simd_blendv(vbest, dots, mask);
			ibest = simd_cast_ftoi(simd_blendv(simd_cast_itof(ibest), simd_cast_itof(idx), mask));
		}
	}

	// Horizontal reduction.
	float hbests[4]; int hbidxs[4];
	simd_store(hbests, vbest); simd_store_i(hbidxs, ibest);
	float hbest = hbests[0]; int hbi = hbidxs[0];
	for (int k = 1; k < 4; k++) { if (hbests[k] > hbest) { hbest = hbests[k]; hbi = hbidxs[k]; } }
	for (; hi < count; hi++) { float hd = dot(verts[hi], ld); if (hd > hbest) { hbest = hd; hbi = hi; } }
	return hbi;
}

// Min/max support scan: find both extreme vertices in one pass over SoA data.
// Returns max and min vertex indices.  Used for cold-start 2-vertex simplex.
static void gjk_hull_minmax_scan(const v3* __restrict verts, int count, const float* __restrict soa, v3 ld, int* out_max, int* out_min)
{
	simd4f ldx = simd_splat(ld.m, 0), ldy = simd_splat(ld.m, 1), ldz = simd_splat(ld.m, 2);
	simd4f vmax = simd_set1(-1e18f), vmin = simd_set1(1e18f);
	simd4i imax = simd_set1_i(0), imin = simd_set1_i(0);
	int hi = 0;

	if (soa) {
		const float* sx = soa, *sy = sx + count, *sz = sy + count;
		simd4i idx = simd_seti(0, 1, 2, 3);
		simd4i four = simd_set1_i(4);
		for (; hi + 3 < count; hi += 4) {
			simd4f d = simd_add(simd_add(simd_mul(simd_load(sx+hi), ldx), simd_mul(simd_load(sy+hi), ldy)), simd_mul(simd_load(sz+hi), ldz));
			simd4f mx = simd_cmpgt(d, vmax);
			vmax = simd_blendv(vmax, d, mx);
			imax = simd_cast_ftoi(simd_blendv(simd_cast_itof(imax), simd_cast_itof(idx), mx));
			simd4f mn = simd_cmpgt(vmin, d);
			vmin = simd_blendv(vmin, d, mn);
			imin = simd_cast_ftoi(simd_blendv(simd_cast_itof(imin), simd_cast_itof(idx), mn));
			idx = simd_add_i(idx, four);
		}
	} else {
		for (; hi + 3 < count; hi += 4) {
			simd4f v0 = verts[hi].m, v1 = verts[hi+1].m, v2 = verts[hi+2].m, v3r = verts[hi+3].m;
			simd4f t0 = simd_unpacklo(v0, v1), t1 = simd_unpacklo(v2, v3r);
			simd4f t2 = simd_unpackhi(v0, v1), t3 = simd_unpackhi(v2, v3r);
			simd4f d = simd_add(simd_add(simd_mul(simd_movelh(t0, t1), ldx), simd_mul(simd_movehl(t1, t0), ldy)), simd_mul(simd_movelh(t2, t3), ldz));
			simd4i idx = simd_seti(hi, hi+1, hi+2, hi+3);
			simd4f mx = simd_cmpgt(d, vmax);
			vmax = simd_blendv(vmax, d, mx);
			imax = simd_cast_ftoi(simd_blendv(simd_cast_itof(imax), simd_cast_itof(idx), mx));
			simd4f mn = simd_cmpgt(vmin, d);
			vmin = simd_blendv(vmin, d, mn);
			imin = simd_cast_ftoi(simd_blendv(simd_cast_itof(imin), simd_cast_itof(idx), mn));
		}
	}

	// Horizontal reduction for max.
	float hmaxv[4]; int hmaxi[4];
	simd_store(hmaxv, vmax); simd_store_i(hmaxi, imax);
	float bmax = hmaxv[0]; int bi_max = hmaxi[0];
	for (int k = 1; k < 4; k++) { if (hmaxv[k] > bmax) { bmax = hmaxv[k]; bi_max = hmaxi[k]; } }
	// Horizontal reduction for min.
	float hminv[4]; int hmini[4];
	simd_store(hminv, vmin); simd_store_i(hmini, imin);
	float bmin = hminv[0]; int bi_min = hmini[0];
	for (int k = 1; k < 4; k++) { if (hminv[k] < bmin) { bmin = hminv[k]; bi_min = hmini[k]; } }
	// Scalar tail.
	for (; hi < count; hi++) {
		float hd = dot(verts[hi], ld);
		if (hd > bmax) { bmax = hd; bi_max = hi; }
		if (hd < bmin) { bmin = hd; bi_min = hi; }
	}
	*out_max = bi_max;
	*out_min = bi_min;
}

// -----------------------------------------------------------------------------
// Support functions (per-shape-type).

static v3 gjk_box_support(const GJK_Shape* __restrict sp, v3 sd, int* __restrict feat)
{
	v3 ld = mat3_mul_v(sp->box.col0, sp->box.col1, sp->box.col2, sd);
	simd4f sign_bits = simd_and(ld.m, simd_sign_mask());
	v3 lc = { .m = simd_xor(sp->box.half_extents.m, sign_bits) };
	*feat = simd_movemask(simd_cmpge(ld.m, simd_zero())) & 7;
	return add(sp->box.center, mat3_tmul_v(sp->box.col0, sp->box.col1, sp->box.col2, lc));
}

static v3 gjk_cylinder_support(const GJK_Shape* __restrict sp, v3 sd, int* __restrict feat)
{
	v3 cu = sp->cylinder.axis;
	if (sp->cylinder.inv_axis_len == 0.0f) { *feat = 0; return sp->cylinder.mid; }
	float cda = dot(sd, cu);
	int cap = cda >= 0.0f;

	// Branchless base: mid + copysign(half_axis, da)
	simd4f sign = simd_and(simd_set1(cda), simd_sign_mask());
	v3 cbase = add(sp->cylinder.mid, (v3){ .m = simd_xor(sp->cylinder.half_axis.m, sign) });
	v3 cdp = sub(sd, scale(cu, cda));
	simd4f cdp2 = v3_dot_m(cdp, cdp);
	simd4f cpl = simd_sqrt_ss(cdp2);
	float cplf = simd_get_x(cpl);
	if (cplf <= FLT_EPSILON) { *feat = cap; return cbase; }

	// Feature ID not needed for cylinder (warm-start uses cached world-space
	// points instead of feature reconstruction).  Skip the expensive 3x10-bit
	// packing entirely.
	*feat = cap;
	float ratio = sp->cylinder.radius / cplf;
	return add(cbase, scale(cdp, ratio));
}

// Support dispatch macro. Box/cylinder/hull use separate functions to control inlining.
#define gjk_support(shape, dir, out_feat, out_point) do {                                                                 \
	GJK_Shape* sp = (shape); v3 sd = (dir);                                                                               \
	SIMD_ASSUME(sp->type >= GJK_POINT && sp->type <= GJK_TRIANGLE);                                                          \
	switch (sp->type) {                                                                                                   \
	case GJK_POINT: *(out_feat) = 0; (out_point) = sp->point.center; break;                                               \
	case GJK_SEGMENT: {                                                                                                   \
		int sf = dot(sd, sub(sp->segment.q, sp->segment.p)) >= 0.0f;                                                      \
		*(out_feat) = sf; (out_point) = sf ? sp->segment.q : sp->segment.p; break;                                         \
	}                                                                                                                     \
	case GJK_BOX: (out_point) = gjk_box_support(sp, sd, (out_feat)); break;                                                \
	case GJK_HULL: {                                                                                                      \
		v3 ld = mat3_mul_v(sp->hull.col0, sp->hull.col1, sp->hull.col2, sd);                                          \
		v3 sld = hmul(ld, sp->hull.scale);                                                                                \
		int hbi = (sp->hull.vert_edge && sp->hull.count > 0)                                                              \
			? gjk_hull_support_climb(sp->hull.verts, sp->hull.edge_twin, sp->hull.edge_next, sp->hull.edge_origin, sp->hull.vert_edge, sld, sp->hull.hint)               \
			: gjk_hull_support_scan(sp->hull.verts, sp->hull.count, sp->hull.soa, sld);                                   \
		sp->hull.hint = hbi;                                                                                               \
		*(out_feat) = hbi;                                                                                                \
		(out_point) = add(sp->hull.center, mat3_tmul_v(sp->hull.col0, sp->hull.col1, sp->hull.col2,                  \
		                  hmul(sp->hull.verts[hbi], sp->hull.scale)));                                                     \
		break;                                                                                                            \
	}                                                                                                                     \
	case GJK_CYLINDER: (out_point) = gjk_cylinder_support(sp, sd, (out_feat)); break;                                      \
	case GJK_TRIANGLE: {                                                                                                  \
		/* 3 dots in parallel via AoS->SoA transpose */                                                                    \
		simd4f v0 = sp->tri.a.m, v1 = sp->tri.b.m, v2 = sp->tri.c.m;                                                     \
		simd4f t01lo = simd_unpacklo(v0, v1), t01hi = simd_unpackhi(v0, v1);                                              \
		simd4f t2lo = simd_unpacklo(v2, simd_zero()), t2hi = simd_unpackhi(v2, simd_zero());                              \
		simd4f xs = simd_movelh(t01lo, t2lo), ys = simd_movehl(t2lo, t01lo), zs = simd_movelh(t01hi, t2hi);              \
		simd4f dx = simd_splat(sd.m, 0), dy = simd_splat(sd.m, 1), dz = simd_splat(sd.m, 2);                             \
		simd4f dots = simd_add(simd_add(simd_mul(xs, dx), simd_mul(ys, dy)), simd_mul(zs, dz));                           \
		float td[4]; simd_store(td, dots);                                                                                 \
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
		return add(sp->box.center, mat3_tmul_v(sp->box.col0, sp->box.col1, sp->box.col2, lc));
	}
	case GJK_HULL:
		return add(sp->hull.center, mat3_tmul_v(sp->hull.col0, sp->hull.col1, sp->hull.col2, hmul(sp->hull.verts[feat], sp->hull.scale)));
	case GJK_CYLINDER: {
		// Cylinder warm-start uses cached world-space points, not feature
		// reconstruction.  Feature ID is just the cap index (0 or 1).
		simd4f sign = (feat & 1) ? simd_zero() : simd_sign_mask();
		return add(sp->cylinder.mid, (v3){ .m = simd_xor(sp->cylinder.half_axis.m, sign) });
	}
	case GJK_TRIANGLE:
		return feat == 0 ? sp->tri.a : (feat == 1 ? sp->tri.b : sp->tri.c);
	}
	return V3(0,0,0);
}

// Shape centroid (used for initial search direction when no cache is available).
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

// -----------------------------------------------------------------------------
// Closest point and witness point extraction macros.

#define gjk_closest_point(simplex, out) do {                                                                        \
	const GJK_Simplex* cs = (simplex);                                                                              \
	float cinv = 1.0f / cs->divisor;                                                                                \
	SIMD_ASSUME(cs->count >= 1 && cs->count <= 3);                                                                  \
	switch (cs->count) {                                                                                            \
	case 1: (out) = cs->v[0].point; break;                                                                          \
	case 2: (out) = add(scale(cs->v[0].point, cs->v[0].u * cinv), scale(cs->v[1].point, cs->v[1].u * cinv)); break; \
	case 3: (out) = add(add(scale(cs->v[0].point, cs->v[0].u * cinv), scale(cs->v[1].point, cs->v[1].u * cinv)),    \
	                     scale(cs->v[2].point, cs->v[2].u * cinv)); break;                                          \
	default: (out) = V3(0,0,0); break;                                                                              \
	}                                                                                                               \
} while(0)

#define gjk_witness_points(simplex, out_p1, out_p2, out_f1, out_f2) do {                                         \
	const GJK_Simplex* ws = (simplex);                                                                           \
	float winv = 1.0f / ws->divisor;                                                                             \
	switch (ws->count) {                                                                                         \
	case 1:                                                                                                      \
		(out_p1) = ws->v[0].point1; (out_p2) = ws->v[0].point2;                                                  \
		*(out_f1) = ws->v[0].feat1; *(out_f2) = ws->v[0].feat2; break;                                           \
	case 2: {                                                                                                    \
		float w0 = ws->v[0].u * winv, w1 = ws->v[1].u * winv;                                                    \
		(out_p1) = add(scale(ws->v[0].point1, w0), scale(ws->v[1].point1, w1));                                  \
		(out_p2) = add(scale(ws->v[0].point2, w0), scale(ws->v[1].point2, w1));                                  \
		int wi = ws->v[1].u > ws->v[0].u;                                                                        \
		*(out_f1) = ws->v[wi].feat1; *(out_f2) = ws->v[wi].feat2; break;                                         \
	}                                                                                                            \
	case 3: {                                                                                                    \
		float w0 = ws->v[0].u * winv, w1 = ws->v[1].u * winv, w2 = ws->v[2].u * winv;                            \
		(out_p1) = add(add(scale(ws->v[0].point1, w0), scale(ws->v[1].point1, w1)), scale(ws->v[2].point1, w2)); \
		(out_p2) = add(add(scale(ws->v[0].point2, w0), scale(ws->v[1].point2, w1)), scale(ws->v[2].point2, w2)); \
		int wi = 0; if (ws->v[1].u > ws->v[wi].u) wi = 1; if (ws->v[2].u > ws->v[wi].u) wi = 2;                  \
		*(out_f1) = ws->v[wi].feat1; *(out_f2) = ws->v[wi].feat2; break;                                         \
	}                                                                                                            \
	default: (out_p1) = V3(0,0,0); (out_p2) = V3(0,0,0); *(out_f1) = 0; *(out_f2) = 0; break;                    \
	}                                                                                                            \
} while(0)

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.

static SIMD_NOINLINE int gjk_solve2(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point;
	v3 ba = sub(b, a);
	float u = dot(b, ba);
	float v = -dot(a, ba);
	float div = u + v;
	if (div == 0.0f) return 0;

	int va = v <= 0.0f, vb = u <= 0.0f;
	if (va | vb) {
		if (vb) s->v[0] = s->v[1];
		s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1;
		return 1;
	}
	s->v[0].u = u; s->v[1].u = v; s->divisor = div; s->count = 2;
	return 1;
}

static SIMD_NOINLINE int gjk_solve3(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point, c = s->v[2].point;
	v3 ba = sub(b, a), cb = sub(c, b), ac = sub(a, c);
	float uAB = dot(b, ba), vAB = -dot(a, ba);
	float uBC = dot(c, cb), vBC = -dot(b, cb);
	float uCA = dot(a, ac), vCA = -dot(c, ac);
	v3 n = cross(ba, sub(c, a));
	float uABC = dot(cross(b, c), n), vABC = dot(cross(c, a), n), wABC = dot(cross(a, b), n);

	// Pack sign bits: bit=1 means value > 0
	int signs = (uAB > 0) | ((vAB > 0) << 1) | ((uBC > 0) << 2) | ((vBC > 0) << 3) | ((uCA > 0) << 4) | ((vCA > 0) << 5) | ((uABC > 0) << 6) | ((vABC > 0) << 7) | ((wABC > 0) << 8);

	if ((signs & 0x12) == 0) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if ((signs & 0x09) == 0) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if ((signs & 0x24) == 0) { s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if ((signs & 0x103) == 0x03) { s->v[0].u = uAB; s->v[1].u = vAB; s->divisor = uAB + vAB; s->count = 2; return 1; }
	if ((signs & 0x4C) == 0x0C) { s->v[0] = s->v[1]; s->v[1] = s->v[2]; s->v[0].u = uBC; s->v[1].u = vBC; s->divisor = uBC + vBC; s->count = 2; return 1; }
	if ((signs & 0xB0) == 0x30) { s->v[1] = s->v[0]; s->v[0] = s->v[2]; s->v[0].u = uCA; s->v[1].u = vCA; s->divisor = uCA + vCA; s->count = 2; return 1; }

	float divABC = uABC + vABC + wABC;
	if (divABC == 0.0f) return 0;
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->divisor = divABC; s->count = 3;
	return 1;
}

static SIMD_NOINLINE int gjk_solve4(GJK_Simplex* s)
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
// Main GJK distance function.
// Termination: containment + monotonic progress + support progress + duplicate vertex.
// Simplex caching: world-space points for non-cylinders, feature replay for cylinders.

#define GJK_MAX_ITERS         64
#define GJK_CONTAINMENT_EPS2  1e-8f
#define GJK_PROGRESS_EPS      1e-4f

static GJK_Result gjk_distance(GJK_Shape* __restrict shapeA, GJK_Shape* __restrict shapeB, GJK_Cache* cache)
{
	GJK_Result result;

	// Restore hull hints from cache.
	GJK_Simplex simplex;
	if (cache) {
		if (shapeA->type == GJK_HULL) shapeA->hull.hint = cache->hintA;
		if (shapeB->type == GJK_HULL) shapeB->hull.hint = cache->hintB;
	}

	// Simplex initialization: reconstruct from cache or cold-start.
	int fA, fB;
	v3 sA, sB;
	int use_simplex_cache = cache && cache->count >= 1 && cache->count <= 4;
	if (use_simplex_cache) {
		// Unrolled cache restore: avoid loop overhead for the hot path.
		#define RESTORE_VERT(i) simplex.v[i].point1 = cache->point1[i]; simplex.v[i].point2 = cache->point2[i]; simplex.v[i].point = sub(cache->point2[i], cache->point1[i]); simplex.v[i].feat1 = cache->feat1[i]; simplex.v[i].feat2 = cache->feat2[i]; simplex.v[i].u = 1.0f;
		RESTORE_VERT(0);
		if (cache->count >= 2) { RESTORE_VERT(1); }
		if (cache->count >= 3) { RESTORE_VERT(2); }
		if (cache->count >= 4) { RESTORE_VERT(3); }
		#undef RESTORE_VERT
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

		// SIMD 2-vertex cold start for hull shapes: min/max scan produces
		// both extreme vertices in one pass, giving GJK a line segment to
		// start from instead of a single point.  Saves ~1 iteration.
		int did_2v = 0;
		if (shapeA->type == GJK_HULL && shapeB->type == GJK_HULL) {
			v3 ldA = hmul(mat3_mul_v(shapeA->hull.col0, shapeA->hull.col1, shapeA->hull.col2, init_d), shapeA->hull.scale);
			v3 ldB = hmul(mat3_mul_v(shapeB->hull.col0, shapeB->hull.col1, shapeB->hull.col2, init_d), shapeB->hull.scale);
			int amx, amn, bmx, bmn;
			gjk_hull_minmax_scan(shapeA->hull.verts, shapeA->hull.count, shapeA->hull.soa, ldA, &amx, &amn);
			gjk_hull_minmax_scan(shapeB->hull.verts, shapeB->hull.count, shapeB->hull.soa, ldB, &bmx, &bmn);
			// Vertex 0: support(A, +d) vs support(B, -d)
			v3 a0 = add(shapeA->hull.center, mat3_tmul_v(shapeA->hull.col0, shapeA->hull.col1, shapeA->hull.col2, hmul(shapeA->hull.verts[amx], shapeA->hull.scale)));
			v3 b0 = add(shapeB->hull.center, mat3_tmul_v(shapeB->hull.col0, shapeB->hull.col1, shapeB->hull.col2, hmul(shapeB->hull.verts[bmn], shapeB->hull.scale)));
			// Vertex 1: support(A, -d) vs support(B, +d)
			v3 a1 = add(shapeA->hull.center, mat3_tmul_v(shapeA->hull.col0, shapeA->hull.col1, shapeA->hull.col2, hmul(shapeA->hull.verts[amn], shapeA->hull.scale)));
			v3 b1 = add(shapeB->hull.center, mat3_tmul_v(shapeB->hull.col0, shapeB->hull.col1, shapeB->hull.col2, hmul(shapeB->hull.verts[bmx], shapeB->hull.scale)));
			v3 w0 = sub(b0, a0), w1 = sub(b1, a1);
			if (len2(sub(w0, w1)) > FLT_EPSILON) {
				simplex.v[0].point1 = a0; simplex.v[0].point2 = b0; simplex.v[0].point = w0; simplex.v[0].feat1 = amx; simplex.v[0].feat2 = bmn; simplex.v[0].u = 1.0f;
				simplex.v[1].point1 = a1; simplex.v[1].point2 = b1; simplex.v[1].point = w1; simplex.v[1].feat1 = amn; simplex.v[1].feat2 = bmx; simplex.v[1].u = 1.0f;
				simplex.divisor = 2.0f; simplex.count = 2;
				did_2v = 1;
			}
		}

		if (!did_2v) {
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
	}

	// Main loop.
	float dsq_prev = FLT_MAX;
	int use_index_term = shapeA->type != GJK_CYLINDER && shapeB->type != GJK_CYLINDER;
	int iter = 0;
	while (iter < GJK_MAX_ITERS) {
		SIMD_ASSUME(simplex.count >= 1 && simplex.count <= 4);
		if (simplex.count > 1) {
			int solved;
			switch (simplex.count) {
			case 2: solved = gjk_solve2(&simplex); break;
			case 3: solved = gjk_solve3(&simplex); break;
			case 4: solved = gjk_solve4(&simplex); break;
			default: SIMD_ASSUME(0);
			}
			if (!solved) break;
		}
		if (simplex.count == 4) break;

		v3 closest; gjk_closest_point(&simplex, closest);
		float dsq = len2(closest);
		if (dsq <= GJK_CONTAINMENT_EPS2 || dsq >= dsq_prev) break;
		dsq_prev = dsq;

		gjk_support(shapeA, closest, &fA, sA);
		gjk_support(shapeB, neg(closest), &fB, sB);
		v3 w = sub(sB, sA);

		if (use_index_term) {
			// Duplicate vertex termination (non-cylinder shapes only).
			int dup = 0;
			for (int i = 0; i < simplex.count; i++) {
				if (simplex.v[i].feat1 == fA && simplex.v[i].feat2 == fB) { dup = 1; break; }
			}
			if (dup) break;
		} else {
			// Relative progress termination (cylinder only).
			float max_vert2 = 0.0f;
			for (int i = 0; i < simplex.count; i++) {
				float v2 = len2(simplex.v[i].point);
				if (v2 > max_vert2) max_vert2 = v2;
			}
			float progress = dsq - dot(w, closest);
			if (progress <= max_vert2 * GJK_PROGRESS_EPS) break;
		}


		iter++;
		GJK_Vertex* vert = &simplex.v[simplex.count];
		vert->point1 = sA;
		vert->point2 = sB;
		vert->point = w;
		vert->feat1 = fA;
		vert->feat2 = fB;
		simplex.count++;
	}

	// Extract witness points and distance.
	gjk_witness_points(&simplex, result.point1, result.point2, &result.feat1, &result.feat2);
	v3 sep = sub(result.point2, result.point1);
	result.distance = len(sep);
	result.iterations = iter;

	// Save cache for next frame.
	if (cache) {
		cache->dir = sep;
		if (shapeA->type == GJK_HULL) cache->hintA = shapeA->hull.hint;
		if (shapeB->type == GJK_HULL) cache->hintB = shapeB->hull.hint;
		cache->count = simplex.count;
		for (int i = 0; i < simplex.count; i++) { cache->feat1[i] = simplex.v[i].feat1; cache->feat2[i] = simplex.v[i].feat2; cache->point1[i] = simplex.v[i].point1; cache->point2[i] = simplex.v[i].point2; }
	}

	// Post-hoc radius for sphere/capsule core shapes.
	float rA = shapeA->radius;
	float rB = shapeB->radius;
	if (rA != 0.0f || rB != 0.0f) {
		v3 normal = scale(sep, 1.0f / result.distance);
		result.point1 = add(result.point1, scale(normal, rA));
		result.point2 = sub(result.point2, scale(normal, rB));
		result.distance = fmaxf(0.0f, result.distance - rA - rB);
	}
	return result;
}

// By-value convenience wrapper for callers that construct shapes inline.
static GJK_Result gjk_distance_v(GJK_Shape a, GJK_Shape b, v3* dir_cache)
{
	GJK_Cache c = {0};
	GJK_Cache* cp = NULL;
	if (dir_cache) { c.dir = *dir_cache; cp = &c; }
	GJK_Result r = gjk_distance(&a, &b, cp);
	if (dir_cache) *dir_cache = c.dir;
	return r;
}
