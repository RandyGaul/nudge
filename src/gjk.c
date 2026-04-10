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
		struct { v3 center; v3 rot_row0, rot_row1, rot_row2; v3 inv_row0, inv_row1, inv_row2; v3 half_extents; } box;
		struct { v3 center; v3 rot_row0, rot_row1, rot_row2; v3 inv_row0, inv_row1, inv_row2; const v3* verts; const float* soa; const HalfEdge* edges; const int* vert_edge; int count; int hint; } hull;
		struct { v3 p, q; float radius; v3 axis; float inv_axis_len; } cylinder;
	};
} GJK_Shape;
static GJK_Shape gjk_sphere(v3 center, float radius) { return (GJK_Shape){ .type = GJK_POINT, .radius = radius, .point.center = center }; }
static GJK_Shape gjk_capsule(v3 p, v3 q, float radius) { return (GJK_Shape){ .type = GJK_SEGMENT, .radius = radius, .segment.p = p, .segment.q = q }; }
// Mat3x3 rotate using SSE: multiply-add across 3 rows without scalar extraction.
static inline v3 gjk_mat_rotate(v3 r0, v3 r1, v3 r2, v3 v) {
	__m128 m0 = _mm_mul_ps(r0.m, v.m);
	__m128 m1 = _mm_mul_ps(r1.m, v.m);
	__m128 m2 = _mm_mul_ps(r2.m, v.m);
	// Horizontal sum each row: [x0+y0+z0, x1+y1+z1, x2+y2+z2, 0]
	// Transpose pairs then add
	__m128 t01lo = _mm_unpacklo_ps(m0, m1);  // x0 x1 y0 y1
	__m128 t01hi = _mm_unpackhi_ps(m0, m1);  // z0 z1 w0 w1
	__m128 t2lo  = _mm_unpacklo_ps(m2, _mm_setzero_ps()); // x2 0 y2 0
	__m128 t2hi  = _mm_unpackhi_ps(m2, _mm_setzero_ps()); // z2 0 w2 0
	__m128 xy = _mm_movelh_ps(t01lo, t2lo);  // x0 x1 x2 0
	__m128 yz = _mm_movehl_ps(t2lo, t01lo);  // y0 y1 y2 0
	__m128 zw = _mm_movelh_ps(t01hi, t2hi);  // z0 z1 z2 0
	return (v3){ .m = _mm_add_ps(_mm_add_ps(xy, yz), zw) };
}
static GJK_Shape gjk_box(v3 center, quat rot, v3 half_extents) {
	// Build rotation matrix rows from quat
	v3 r0 = quat_rotate(rot, V3(1,0,0)), r1 = quat_rotate(rot, V3(0,1,0)), r2 = quat_rotate(rot, V3(0,0,1));
	quat ir = inv(rot);
	v3 i0 = quat_rotate(ir, V3(1,0,0)), i1 = quat_rotate(ir, V3(0,1,0)), i2 = quat_rotate(ir, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_BOX, .box.center = center, .box.rot_row0 = r0, .box.rot_row1 = r1, .box.rot_row2 = r2, .box.inv_row0 = i0, .box.inv_row1 = i1, .box.inv_row2 = i2, .box.half_extents = half_extents };
}
static GJK_Shape gjk_hull(v3 center, quat rot, const v3* verts, int count, const float* soa, const HalfEdge* edges, const int* vert_edge) {
	v3 r0 = quat_rotate(rot, V3(1,0,0)), r1 = quat_rotate(rot, V3(0,1,0)), r2 = quat_rotate(rot, V3(0,0,1));
	quat ir = inv(rot);
	v3 i0 = quat_rotate(ir, V3(1,0,0)), i1 = quat_rotate(ir, V3(0,1,0)), i2 = quat_rotate(ir, V3(0,0,1));
	return (GJK_Shape){ .type = GJK_HULL, .hull.center = center, .hull.rot_row0 = r0, .hull.rot_row1 = r1, .hull.rot_row2 = r2, .hull.inv_row0 = i0, .hull.inv_row1 = i1, .hull.inv_row2 = i2, .hull.verts = verts, .hull.soa = soa, .hull.edges = edges, .hull.vert_edge = vert_edge, .hull.count = count };
}
static GJK_Shape gjk_cylinder(v3 p, v3 q, float radius) {
	v3 axis = sub(q, p);
	float al = len(axis);
	float inv_al = al > FLT_EPSILON ? 1.0f / al : 0.0f;
	return (GJK_Shape){ .type = GJK_CYLINDER, .cylinder.p = p, .cylinder.q = q, .cylinder.radius = radius, .cylinder.axis = scale(axis, inv_al), .cylinder.inv_axis_len = inv_al };
}
// Hull convenience: pre-scale vertices and build shape.
// Caller must keep scaled_verts alive for the duration of the GJK call.
// For hulls > 32 verts, also builds SoA layout in soa_buf for fast SIMD scan.
// soa_buf must hold 3*count floats (or NULL for small hulls).
// vert_edge_buf must hold hull->vert_count ints. Builds per-vertex first-edge lookup.
static GJK_Shape gjk_hull_scaled(const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts, float* soa_buf)
{
	int n = hull->vert_count;
	for (int i = 0; i < n; i++)
		scaled_verts[i] = hmul(hull->verts[i], sc);
	const float* soa = NULL;
	if (soa_buf && n > 32) {
		float* sx = soa_buf, *sy = soa_buf + n, *sz = soa_buf + n * 2;
		for (int i = 0; i < n; i++) { sx[i] = scaled_verts[i].x; sy[i] = scaled_verts[i].y; sz[i] = scaled_verts[i].z; }
		soa = soa_buf;
	}
	// Build per-vertex first-edge lookup for hill climbing (one-time O(E) setup)
	static int vert_edge_buf[1024]; // max verts
	if (hull->edges && n <= 1024) {
		for (int i = 0; i < n; i++) vert_edge_buf[i] = -1;
		for (int i = 0; i < hull->edge_count; i++)
			if (vert_edge_buf[hull->edges[i].origin] < 0) vert_edge_buf[hull->edges[i].origin] = i;
		return gjk_hull(pos, rot, scaled_verts, n, soa, hull->edges, vert_edge_buf);
	}
	return gjk_hull(pos, rot, scaled_verts, n, soa, hull->edges, NULL);
}
// Hull support scan: extracted to a function to reduce macro expansion size,
// allowing the compiler to inline quat_rotate in the caller.
// Hill-climbing support: walk vertex adjacency graph from a start vertex.
// O(sqrt(n)) expected for convex hulls vs O(n) for linear scan.
// Hill-climbing support: walk vertex adjacency graph from vertex 0.
// Ring walk: next outgoing edge from v = edges[current ^ 1].next
static int gjk_hull_support_climb(const v3* verts, const HalfEdge* edges, const int* vert_edge, v3 ld, int start)
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

static int gjk_hull_support_scan(const v3* verts, int count, const float* soa, v3 ld)
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
static v3 gjk_box_support(const GJK_Shape* sp, v3 sd, int* feat)
{
	v3 ld = gjk_mat_rotate(sp->box.inv_row0, sp->box.inv_row1, sp->box.inv_row2, sd);
	v3 he = sp->box.half_extents;
	v3 lc = V3(ld.x >= 0 ? he.x : -he.x, ld.y >= 0 ? he.y : -he.y, ld.z >= 0 ? he.z : -he.z);
	*feat = (ld.x >= 0.0f) | ((ld.y >= 0.0f) << 1) | ((ld.z >= 0.0f) << 2);
	return add(sp->box.center, gjk_mat_rotate(sp->box.rot_row0, sp->box.rot_row1, sp->box.rot_row2, lc));
}
static v3 gjk_cylinder_support(const GJK_Shape* sp, v3 sd, int* feat)
{
	v3 cu = sp->cylinder.axis;
	if (sp->cylinder.inv_axis_len == 0.0f) { *feat = 0; return sp->cylinder.p; }
	float cda = dot(sd, cu); int cf = cda >= 0.0f; *feat = cf;
	v3 cbase = cf ? sp->cylinder.q : sp->cylinder.p;
	v3 cdp = sub(sd, scale(cu, cda)); float cpl = len(cdp);
	return (cpl > FLT_EPSILON) ? add(cbase, scale(cdp, sp->cylinder.radius / cpl)) : cbase;
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
		v3 ld = gjk_mat_rotate(sp->hull.inv_row0, sp->hull.inv_row1, sp->hull.inv_row2, sd);                               \
		int hbi = (sp->hull.vert_edge && sp->hull.count > 0)                                                              \
			? gjk_hull_support_climb(sp->hull.verts, sp->hull.edges, sp->hull.vert_edge, ld, sp->hull.hint)                 \
			: gjk_hull_support_scan(sp->hull.verts, sp->hull.count, sp->hull.soa, ld);                                     \
		sp->hull.hint = hbi;                                                                                               \
		*(out_feat) = hbi;                                                                                                \
		(out_point) = add(sp->hull.center, gjk_mat_rotate(sp->hull.rot_row0, sp->hull.rot_row1, sp->hull.rot_row2,         \
		                  sp->hull.verts[hbi]));                                                                           \
		break;                                                                                                            \
	}                                                                                                                     \
	case GJK_CYLINDER: (out_point) = gjk_cylinder_support(sp, sd, (out_feat)); break;                                      \
	default: *(out_feat) = 0; (out_point) = V3(0,0,0); break;                                                             \
	}                                                                                                                     \
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
static void gjk_witness_points(const GJK_Simplex* s, v3* p1, v3* p2, int* f1, int* f2)
{
	float inv = 1.0f / s->divisor;
	*p1 = V3(0,0,0); *p2 = V3(0,0,0);
	float best_u = -1.0f;
	*f1 = 0; *f2 = 0;
	for (int i = 0; i < s->count; i++) {
		float w = s->v[i].u * inv;
		*p1 = add(*p1, scale(s->v[i].point1, w));
		*p2 = add(*p2, scale(s->v[i].point2, w));
		if (s->v[i].u > best_u) { best_u = s->v[i].u; *f1 = s->v[i].feat1; *f2 = s->v[i].feat2; }
	}
}
// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.
static int gjk_solve2(GJK_Simplex* s)
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
static int gjk_solve3(GJK_Simplex* s)
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
	s->divisor = uABC + vABC + wABC;
	if (s->divisor == 0.0f) return 0;
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC; s->count = 3;
	return 1;
}
#define stp(a, b, c) dot(a, cross(b, c))
static int gjk_solve4(GJK_Simplex* s)
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
static GJK_Result gjk_distance(GJK_Shape shapeA, GJK_Shape shapeB)
{
	GJK_Result result = {0};
	GJK_Simplex simplex = {0};
	v3 init_d = sub(gjk_center(&shapeB), gjk_center(&shapeA));
	if (len2(init_d) < FLT_EPSILON) init_d = V3(1, 0, 0);
	int fA, fB;
	v3 sA; gjk_support(&shapeA, init_d, &fA, sA);
	v3 sB; gjk_support(&shapeB, neg(init_d), &fB, sB);
	simplex.v[0].point1 = sA;
	simplex.v[0].point2 = sB;
	simplex.v[0].point = sub(sB, sA);
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
		v3 closest; gjk_closest_point(&simplex, closest);
		float dsq = len2(closest);
		if (dsq <= GJK_CONTAINMENT_EPS2) break;
		if (dsq >= dsq_prev) break;
		dsq_prev = dsq;
		gjk_support(&shapeA, closest, &fA, sA);
		gjk_support(&shapeB, neg(closest), &fB, sB);
		v3 w = sub(sB, sA);
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
