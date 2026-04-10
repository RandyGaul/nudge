// See LICENSE for licensing info.
// tests_gjk_perf.c -- GJK correctness and performance tests.
//
// Correctness: known-distance cases, brute-force reference comparison,
// near-contact tracking.
// Performance: timing across shape pairs.

#include <stdio.h>
#include <math.h>
#include <float.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
static double qpc_freq;
static void qpc_init() { LARGE_INTEGER f; QueryPerformanceFrequency(&f); qpc_freq = (double)f.QuadPart; }
static double qpc_now() { LARGE_INTEGER c; QueryPerformanceCounter(&c); return (double)c.QuadPart / qpc_freq; }
#else
#include <time.h>
static void qpc_init() {}
static double qpc_now() { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
#endif

// =============================================================================
// Double-precision math (extends dv3/dquat from solver_ldl.c).

static inline double dv3_len2(dv3 a) { return dv3_dot(a, a); }
static inline dv3 dv3_neg(dv3 a) { return DV3(-a.x, -a.y, -a.z); }
static inline double dv3_stp(dv3 a, dv3 b, dv3 c) { return dv3_dot(a, dv3_cross(b, c)); }

static inline dquat dquat_identity() { return (dquat){0,0,0,1}; }
static inline dquat dquat_inv(dquat q) { return (dquat){-q.x, -q.y, -q.z, q.w}; }
static inline dquat dquat_from_quat(quat q) { return (dquat){q.x, q.y, q.z, q.w}; }
static inline dv3 dquat_rotate(dquat q, dv3 v) {
	dv3 u = {q.x, q.y, q.z};
	double s = q.w;
	return dv3_add(dv3_add(dv3_scale(u, 2.0 * dv3_dot(u, v)), dv3_scale(v, s*s - dv3_dot(u, u))), dv3_scale(dv3_cross(u, v), 2.0 * s));
}

// =============================================================================
// Brute-force double-precision reference for hull-involving pairs.

static dv3 bf_closest_on_segment(dv3 p, dv3 a, dv3 b)
{
	dv3 ab = dv3_sub(b, a);
	double d2 = dv3_dot(ab, ab);
	if (d2 < 1e-30) return a;
	double t = dv3_dot(dv3_sub(p, a), ab) / d2;
	if (t < 0.0) t = 0.0;
	if (t > 1.0) t = 1.0;
	return dv3_add(a, dv3_scale(ab, t));
}

static void bf_closest_segments(dv3 a1, dv3 a2, dv3 b1, dv3 b2, dv3* ca, dv3* cb)
{
	dv3 da = dv3_sub(a2, a1), db = dv3_sub(b2, b1), r = dv3_sub(a1, b1);
	double a = dv3_dot(da, da), e = dv3_dot(db, db), f = dv3_dot(db, r);
	double s, t;

	if (a <= 1e-30 && e <= 1e-30) { *ca = a1; *cb = b1; return; }
	if (a <= 1e-30) {
		s = 0.0; t = f / e;
		if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
	} else {
		double c = dv3_dot(da, r);
		if (e <= 1e-30) {
			t = 0.0; s = -c / a;
			if (s < 0.0) s = 0.0; else if (s > 1.0) s = 1.0;
		} else {
			double b = dv3_dot(da, db);
			double denom = a * e - b * b;
			s = denom > 1e-30 ? (b * f - c * e) / denom : 0.0;
			if (s < 0.0) s = 0.0; else if (s > 1.0) s = 1.0;
			t = (b * s + f) / e;
			if (t < 0.0) { t = 0.0; s = -c / a; if (s < 0.0) s = 0.0; else if (s > 1.0) s = 1.0; }
			else if (t > 1.0) { t = 1.0; s = (b - c) / a; if (s < 0.0) s = 0.0; else if (s > 1.0) s = 1.0; }
		}
	}
	*ca = dv3_add(a1, dv3_scale(da, s));
	*cb = dv3_add(b1, dv3_scale(db, t));
}

static dv3 bf_closest_on_hull(dv3 point, const dv3* verts, const Hull* hull)
{
	double best2 = 1e30;
	dv3 best = DV3(0,0,0);

	for (int i = 0; i < hull->vert_count; i++) {
		double d2 = dv3_len2(dv3_sub(point, verts[i]));
		if (d2 < best2) { best2 = d2; best = verts[i]; }
	}

	for (int i = 0; i < hull->edge_count; i += 2) {
		dv3 a = verts[hull->edges[i].origin];
		dv3 b = verts[hull->edges[i ^ 1].origin];
		dv3 cp = bf_closest_on_segment(point, a, b);
		double d2 = dv3_len2(dv3_sub(point, cp));
		if (d2 < best2) { best2 = d2; best = cp; }
	}

	for (int fi = 0; fi < hull->face_count; fi++) {
		dv3 fv[32];
		int fc = 0;
		int edge = hull->faces[fi].edge;
		int first = edge;
		do {
			fv[fc++] = verts[hull->edges[edge].origin];
			edge = hull->edges[edge].next;
		} while (edge != first && fc < 32);
		if (fc < 3) continue;

		dv3 fn = dv3_cross(dv3_sub(fv[1], fv[0]), dv3_sub(fv[2], fv[0]));
		double fn_len = dv3_len(fn);
		if (fn_len < 1e-30) continue;
		fn = dv3_scale(fn, 1.0 / fn_len);

		dv3 proj = dv3_sub(point, dv3_scale(fn, dv3_dot(dv3_sub(point, fv[0]), fn)));

		// Inside test with relative epsilon: scale by max edge length squared
		// to handle irregular quickhull faces where absolute 1e-12 is too tight.
		double max_edge2 = 0;
		for (int i = 0; i < fc; i++) {
			double e2 = dv3_len2(dv3_sub(fv[(i+1)%fc], fv[i]));
			if (e2 > max_edge2) max_edge2 = e2;
		}
		double face_eps = max_edge2 * 1e-10;

		int inside = 1;
		for (int i = 0; i < fc; i++) {
			dv3 edir = dv3_sub(fv[(i + 1) % fc], fv[i]);
			dv3 outward = dv3_cross(edir, fn);
			if (dv3_dot(dv3_sub(proj, fv[i]), outward) > face_eps) { inside = 0; break; }
		}
		if (inside) {
			double d2 = dv3_len2(dv3_sub(point, proj));
			if (d2 < best2) { best2 = d2; best = proj; }
		}
	}

	return best;
}

typedef struct BruteResult { dv3 point1, point2; double distance; } BruteResult;

static void bf_hull_world_verts(dv3* out, const Hull* hull, v3 sc, dv3 pos, dquat rot)
{
	for (int i = 0; i < hull->vert_count; i++) {
		dv3 local = DV3((double)hull->verts[i].x * sc.x, (double)hull->verts[i].y * sc.y, (double)hull->verts[i].z * sc.z);
		out[i] = dv3_add(pos, dquat_rotate(rot, local));
	}
}

// Closest point on a segment to a convex face (polygon in 3D).
// Tests: project segment's closest point to face plane, check if inside polygon.
// Also clips the segment to the face plane and tests the intersection point.
static void bf_closest_segment_face(dv3 seg_a, dv3 seg_b, const dv3* fv, int fc, dv3 fn, double face_eps, double* best2, dv3* best_seg, dv3* best_face)
{
	// For each point along the segment, find its projection onto the face plane.
	// The closest segment point to the face is the one whose plane-projection is inside the face.
	// Test both endpoints and the perpendicular foot.
	dv3 seg_dir = dv3_sub(seg_b, seg_a);
	double seg_len2 = dv3_dot(seg_dir, seg_dir);

	// Test several candidate points on the segment
	for (int ep = 0; ep < 3; ep++) {
		dv3 sp;
		if (ep == 0) sp = seg_a;
		else if (ep == 1) sp = seg_b;
		else {
			// Closest point on the infinite line to the face plane (perpendicular foot)
			double denom = dv3_dot(seg_dir, fn);
			if (fabs(denom) < 1e-30) continue;
			double t = dv3_dot(dv3_sub(fv[0], seg_a), fn) / denom;
			if (t < 0.0) t = 0.0;
			if (t > 1.0) t = 1.0;
			sp = dv3_add(seg_a, dv3_scale(seg_dir, t));
		}
		// Project sp onto face plane
		dv3 proj = dv3_sub(sp, dv3_scale(fn, dv3_dot(dv3_sub(sp, fv[0]), fn)));
		// Inside check
		int inside = 1;
		for (int i = 0; i < fc; i++) {
			dv3 edir = dv3_sub(fv[(i+1)%fc], fv[i]);
			dv3 outward = dv3_cross(edir, fn);
			if (dv3_dot(dv3_sub(proj, fv[i]), outward) > face_eps) { inside = 0; break; }
		}
		if (inside) {
			double d2 = dv3_len2(dv3_sub(sp, proj));
			if (d2 < *best2) { *best2 = d2; *best_seg = sp; *best_face = proj; }
		}
	}
}

static BruteResult bf_hull_hull(const dv3* va, const Hull* ha, const dv3* vb, const Hull* hb)
{
	double best2 = 1e30;
	dv3 best_a = DV3(0,0,0), best_b = DV3(0,0,0);

	// Vertex A -> hull B
	for (int i = 0; i < ha->vert_count; i++) {
		dv3 cp = bf_closest_on_hull(va[i], vb, hb);
		double d2 = dv3_len2(dv3_sub(va[i], cp));
		if (d2 < best2) { best2 = d2; best_a = va[i]; best_b = cp; }
	}
	// Vertex B -> hull A
	for (int i = 0; i < hb->vert_count; i++) {
		dv3 cp = bf_closest_on_hull(vb[i], va, ha);
		double d2 = dv3_len2(dv3_sub(vb[i], cp));
		if (d2 < best2) { best2 = d2; best_a = cp; best_b = vb[i]; }
	}
	// Edge A x Edge B
	for (int i = 0; i < ha->edge_count; i += 2) {
		dv3 a1 = va[ha->edges[i].origin], a2 = va[ha->edges[i ^ 1].origin];
		for (int j = 0; j < hb->edge_count; j += 2) {
			dv3 b1 = vb[hb->edges[j].origin], b2 = vb[hb->edges[j ^ 1].origin];
			dv3 ca, cb;
			bf_closest_segments(a1, a2, b1, b2, &ca, &cb);
			double d2 = dv3_len2(dv3_sub(ca, cb));
			if (d2 < best2) { best2 = d2; best_a = ca; best_b = cb; }
		}
	}
	// Edge A -> face B
	for (int i = 0; i < ha->edge_count; i += 2) {
		dv3 a1 = va[ha->edges[i].origin], a2 = va[ha->edges[i ^ 1].origin];
		for (int fi = 0; fi < hb->face_count; fi++) {
			dv3 fv[32]; int fc = 0;
			int edge = hb->faces[fi].edge; int first = edge;
			do { fv[fc++] = vb[hb->edges[edge].origin]; edge = hb->edges[edge].next; } while (edge != first && fc < 32);
			if (fc < 3) continue;
			dv3 fn = dv3_cross(dv3_sub(fv[1], fv[0]), dv3_sub(fv[2], fv[0]));
			double fn_len = dv3_len(fn); if (fn_len < 1e-30) continue;
			fn = dv3_scale(fn, 1.0 / fn_len);
			double max_e2 = 0; for (int k = 0; k < fc; k++) { double e2 = dv3_len2(dv3_sub(fv[(k+1)%fc], fv[k])); if (e2 > max_e2) max_e2 = e2; }
			dv3 bs = DV3(0,0,0), bf = DV3(0,0,0);
			bf_closest_segment_face(a1, a2, fv, fc, fn, max_e2 * 1e-10, &best2, &best_a, &best_b);
		}
	}
	// Edge B -> face A
	for (int i = 0; i < hb->edge_count; i += 2) {
		dv3 b1 = vb[hb->edges[i].origin], b2 = vb[hb->edges[i ^ 1].origin];
		for (int fi = 0; fi < ha->face_count; fi++) {
			dv3 fv[32]; int fc = 0;
			int edge = ha->faces[fi].edge; int first = edge;
			do { fv[fc++] = va[ha->edges[edge].origin]; edge = ha->edges[edge].next; } while (edge != first && fc < 32);
			if (fc < 3) continue;
			dv3 fn = dv3_cross(dv3_sub(fv[1], fv[0]), dv3_sub(fv[2], fv[0]));
			double fn_len = dv3_len(fn); if (fn_len < 1e-30) continue;
			fn = dv3_scale(fn, 1.0 / fn_len);
			double max_e2 = 0; for (int k = 0; k < fc; k++) { double e2 = dv3_len2(dv3_sub(fv[(k+1)%fc], fv[k])); if (e2 > max_e2) max_e2 = e2; }
			bf_closest_segment_face(b1, b2, fv, fc, fn, max_e2 * 1e-10, &best2, &best_b, &best_a);
		}
	}
	return (BruteResult){ best_a, best_b, sqrt(best2) };
}

static BruteResult bf_segment_hull(dv3 seg_a, dv3 seg_b, const dv3* verts, const Hull* hull)
{
	double best2 = 1e30;
	dv3 best_seg = DV3(0,0,0), best_hull = DV3(0,0,0);

	for (int ep = 0; ep < 2; ep++) {
		dv3 p = ep == 0 ? seg_a : seg_b;
		dv3 cp = bf_closest_on_hull(p, verts, hull);
		double d2 = dv3_len2(dv3_sub(p, cp));
		if (d2 < best2) { best2 = d2; best_seg = p; best_hull = cp; }
	}
	for (int i = 0; i < hull->vert_count; i++) {
		dv3 cp = bf_closest_on_segment(verts[i], seg_a, seg_b);
		double d2 = dv3_len2(dv3_sub(cp, verts[i]));
		if (d2 < best2) { best2 = d2; best_seg = cp; best_hull = verts[i]; }
	}
	for (int i = 0; i < hull->edge_count; i += 2) {
		dv3 b1 = verts[hull->edges[i].origin], b2 = verts[hull->edges[i ^ 1].origin];
		dv3 ca, cb;
		bf_closest_segments(seg_a, seg_b, b1, b2, &ca, &cb);
		double d2 = dv3_len2(dv3_sub(ca, cb));
		if (d2 < best2) { best2 = d2; best_seg = ca; best_hull = cb; }
	}
	return (BruteResult){ best_seg, best_hull, sqrt(best2) };
}

// =============================================================================
// RNG.

static uint32_t gjk_perf_rng = 12345;

static float gjk_perf_randf()
{
	gjk_perf_rng = gjk_perf_rng * 1664525u + 1013904223u;
	return (float)(gjk_perf_rng >> 8) / (float)(1 << 24);
}

static quat gjk_perf_random_quat()
{
	float u1 = gjk_perf_randf(), u2 = gjk_perf_randf() * 6.2831853f, u3 = gjk_perf_randf() * 6.2831853f;
	float s1 = sqrtf(1.0f - u1), s2 = sqrtf(u1);
	return (quat){ s1 * sinf(u2), s1 * cosf(u2), s2 * sinf(u3), s2 * cosf(u3) };
}

static Hull* gjk_perf_make_cylinder_hull(float half_height, float radius, int segments)
{
	v3 verts[128];
	int n = 0;
	for (int cap = 0; cap < 2; cap++) {
		float y = cap == 0 ? -half_height : half_height;
		for (int i = 0; i < segments; i++) {
			float angle = 6.2831853f * (float)i / (float)segments;
			verts[n++] = V3(radius * cosf(angle), y, radius * sinf(angle));
		}
	}
	return quickhull(verts, n);
}

// Generate a random convex hull with approximately n_target vertices.
// Points are sampled uniformly on a sphere of given radius, then run through
// quickhull. The actual hull vertex count may be less than n_target.
static Hull* gjk_perf_make_random_hull(int n_target, float radius)
{
	v3* pts = CK_ALLOC(sizeof(v3) * n_target);
	for (int i = 0; i < n_target; i++) {
		// Random direction via rejection sampling on unit cube
		float x, y, z, d2;
		do { x = gjk_perf_randf()*2-1; y = gjk_perf_randf()*2-1; z = gjk_perf_randf()*2-1; d2 = x*x+y*y+z*z; } while (d2 < 0.01f || d2 > 1.0f);
		float inv = radius / sqrtf(d2);
		// Add slight radial jitter so quickhull doesn't degenerate to a sphere shell
		float r = radius * (0.7f + 0.3f * gjk_perf_randf());
		inv = r / sqrtf(d2);
		pts[i] = V3(x * inv, y * inv, z * inv);
	}
	Hull* h = quickhull(pts, n_target);
	CK_FREE(pts);
	return h;
}

// =============================================================================
// Correctness: known distances for every shape pair.

static void test_gjk_known_distances()
{
	GJK_Result r;
	quat id = quat_identity();

	// --- Sphere ---
	TEST_BEGIN("sphere-sphere separated");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_sphere(V3(3,0,0), 1), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.0f, 0.01f);

	TEST_BEGIN("sphere-sphere touching");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_sphere(V3(2,0,0), 1), NULL);
	TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

	TEST_BEGIN("sphere-sphere overlap");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_sphere(V3(1,0,0), 1), NULL);
	TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

	TEST_BEGIN("sphere-box face");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_box(V3(4,0,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);

	TEST_BEGIN("sphere-box edge");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 0.5f), gjk_box(V3(2,2,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, sqrtf(2.0f) - 0.5f, 0.02f);

	TEST_BEGIN("sphere-box vertex");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 0.5f), gjk_box(V3(3,3,3), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, sqrtf(12.0f) - 0.5f, 0.02f);

	TEST_BEGIN("sphere-capsule separated");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_capsule(V3(4,-1,0), V3(4,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.5f, 0.01f);

	TEST_BEGIN("sphere-cylinder separated");
	r = gjk_distance_v(gjk_sphere(V3(0,0,0), 1), gjk_cylinder(V3(4,-1,0), V3(4,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.5f, 0.02f);

	// --- Capsule ---
	TEST_BEGIN("capsule-capsule parallel separated");
	r = gjk_distance_v(gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_capsule(V3(3,-1,0), V3(3,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);

	TEST_BEGIN("capsule-capsule crossed");
	r = gjk_distance_v(gjk_capsule(V3(0,-2,0), V3(0,2,0), 0.3f), gjk_capsule(V3(2,0,-2), V3(2,0,2), 0.3f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.4f, 0.01f);

	TEST_BEGIN("capsule-capsule touching");
	r = gjk_distance_v(gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_capsule(V3(1,-1,0), V3(1,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

	TEST_BEGIN("capsule-box face");
	r = gjk_distance_v(gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_box(V3(3,0,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.5f, 0.01f);

	TEST_BEGIN("capsule-cylinder separated");
	r = gjk_distance_v(gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_cylinder(V3(3,-1,0), V3(3,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.02f);

	// --- Box ---
	TEST_BEGIN("box-box face separated");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(V3(4,0,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);

	TEST_BEGIN("box-box face touching");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(V3(2,0,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

	TEST_BEGIN("box-box edge-edge");
	{
		quat rot45z = { 0, 0, 0.3826834f, 0.9238795f };
		float expected = 1.0f; // gap between face at x=1 and rotated edge
		r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(V3(1 + sqrtf(2.0f) + expected, 0, 0), rot45z, V3(1,1,1)), NULL);
		TEST_ASSERT_FLOAT(r.distance, expected, 0.02f);
	}

	TEST_BEGIN("box-box different sizes");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(0.5f,0.5f,0.5f)), gjk_box(V3(3,0,0), id, V3(2,2,2)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 0.5f, 0.01f);

	TEST_BEGIN("box-box Y axis");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(V3(0,5,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 3.0f, 0.01f);

	TEST_BEGIN("box-box Z axis");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(V3(0,0,4), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);

	// --- Cylinder ---
	TEST_BEGIN("cylinder-box face");
	r = gjk_distance_v(gjk_cylinder(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_box(V3(3,0,0), id, V3(1,1,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.5f, 0.02f);

	TEST_BEGIN("cylinder-cylinder parallel");
	r = gjk_distance_v(gjk_cylinder(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_cylinder(V3(3,-1,0), V3(3,1,0), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.02f);

	TEST_BEGIN("cylinder-cylinder crossed");
	r = gjk_distance_v(gjk_cylinder(V3(0,-2,0), V3(0,2,0), 0.5f), gjk_cylinder(V3(3,0,-2), V3(3,0,2), 0.5f), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.02f);

	// --- Triangle ---
	TEST_BEGIN("sphere-triangle face");
	r = gjk_distance_v(gjk_sphere(V3(0,3,0), 0.5f), gjk_triangle(V3(-1,0,-1), V3(1,0,-1), V3(0,0,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.5f, 0.01f);

	TEST_BEGIN("sphere-triangle edge");
	r = gjk_distance_v(gjk_sphere(V3(3,0,0), 0.5f), gjk_triangle(V3(0,0,0), V3(1,0,0), V3(0.5f,1,0)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.5f, 0.01f);

	TEST_BEGIN("sphere-triangle vertex");
	r = gjk_distance_v(gjk_sphere(V3(5,0,0), 1.0f), gjk_triangle(V3(0,0,0), V3(1,0,0), V3(0.5f,1,0)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 3.0f, 0.01f);

	TEST_BEGIN("box-triangle separated");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_triangle(V3(3,0,0), V3(4,1,0), V3(3.5f,0,1)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.05f);

	TEST_BEGIN("capsule-triangle");
	r = gjk_distance_v(gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f), gjk_triangle(V3(3,0,0), V3(4,1,0), V3(3.5f,0,1)), NULL);
	TEST_ASSERT(r.distance > 0.0f); // just verify non-negative, exact value depends on geometry

	// --- Hull ---
	TEST_BEGIN("hull-hull (box via hull)");
	{
		const Hull* ub = hull_unit_box();
		v3 s1[8], s2[8]; float so1[24], so2[24];
		GJK_Shape ga = gjk_hull_scaled(ub, V3(0,0,0), id, V3(1,1,1), s1, so1);
		GJK_Shape gb = gjk_hull_scaled(ub, V3(4,0,0), id, V3(1,1,1), s2, so2);
		r = gjk_distance_v(ga, gb, NULL);
		TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);
	}

	TEST_BEGIN("point-hull vertex region");
	{
		const Hull* ub = hull_unit_box();
		v3 s[8]; float so[24];
		GJK_Shape pt = gjk_sphere(V3(2,2,2), 0);
		GJK_Shape bx = gjk_hull_scaled(ub, V3(0,0,0), id, V3(1,1,1), s, so);
		r = gjk_distance_v(pt, bx, NULL);
		TEST_ASSERT_FLOAT(r.distance, sqrtf(3.0f), 0.02f);
	}

	TEST_BEGIN("point-hull edge region");
	{
		const Hull* ub = hull_unit_box();
		v3 s[8]; float so[24];
		GJK_Shape pt = gjk_sphere(V3(1.5f,1.5f,0), 0);
		GJK_Shape bx = gjk_hull_scaled(ub, V3(0,0,0), id, V3(1,1,1), s, so);
		r = gjk_distance_v(pt, bx, NULL);
		TEST_ASSERT_FLOAT(r.distance, sqrtf(0.5f), 0.02f);
	}

	TEST_BEGIN("point-hull face region");
	{
		const Hull* ub = hull_unit_box();
		v3 s[8]; float so[24];
		GJK_Shape pt = gjk_sphere(V3(3,0,0), 0);
		GJK_Shape bx = gjk_hull_scaled(ub, V3(0,0,0), id, V3(1,1,1), s, so);
		r = gjk_distance_v(pt, bx, NULL);
		TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);
	}

	// --- Asymmetric scales ---
	TEST_BEGIN("box-box asymmetric half_extents");
	r = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(3,0.5f,0.5f)), gjk_box(V3(5,0,0), id, V3(0.5f,3,3)), NULL);
	TEST_ASSERT_FLOAT(r.distance, 1.5f, 0.01f);

	// --- Rotated box ---
	TEST_BEGIN("rotated box-box");
	{
		// 45 deg around Z, box corner at sqrt(2) from center
		quat rot45z = { 0, 0, 0.3826834f, 0.9238795f };
		GJK_Shape a = gjk_box(V3(0,0,0), rot45z, V3(1,1,1));
		GJK_Shape b = gjk_box(V3(5,0,0), id, V3(1,1,1));
		r = gjk_distance_v(a, b, NULL);
		float expected = 5.0f - sqrtf(2.0f) - 1.0f;
		TEST_ASSERT_FLOAT(r.distance, expected, 0.02f);
	}
}

// =============================================================================
// Correctness: brute-force comparison, all hull-involving pairs.
// 500 random configs per pair, random rotations/scales/separations.

#define BF_NCASES 500
#define BF_TOL    1e-3

static void test_gjk_bf_box_box()
{
	const Hull* ub = hull_unit_box();
	gjk_perf_rng = 777;
	double max_err_hull = 0, max_err_box = 0;

	for (int c = 0; c < BF_NCASES; c++) {
		float sep = 0.1f + gjk_perf_randf() * 10.0f;
		float sc = 0.3f + gjk_perf_randf() * 5.0f;
		quat rA = gjk_perf_random_quat(), rB = gjk_perf_random_quat();
		v3 heA = V3(sc, sc * (0.5f + gjk_perf_randf()), sc * (0.3f + gjk_perf_randf()));
		v3 heB = V3(sc * (0.4f + gjk_perf_randf()), sc, sc * (0.5f + gjk_perf_randf()));
		float extA = sqrtf(heA.x*heA.x + heA.y*heA.y + heA.z*heA.z);
		float extB = sqrtf(heB.x*heB.x + heB.y*heB.y + heB.z*heB.z);
		v3 posA = V3(0, 0, 0), posB = V3(sep + extA + extB, gjk_perf_randf() - 0.5f, gjk_perf_randf() - 0.5f);

		dv3 wa[8], wb[8];
		bf_hull_world_verts(wa, ub, heA, dv3_from_v3(posA), dquat_from_quat(rA));
		bf_hull_world_verts(wb, ub, heB, dv3_from_v3(posB), dquat_from_quat(rB));
		BruteResult ref = bf_hull_hull(wa, ub, wb, ub);

		// gjk_hull path
		v3 s1[8], s2[8]; float so1[24], so2[24];
		GJK_Shape ga = gjk_hull_scaled(ub, posA, rA, heA, s1, so1);
		GJK_Shape gb = gjk_hull_scaled(ub, posB, rB, heB, s2, so2);
		double err = fabs((double)gjk_distance_v(ga, gb, NULL).distance - ref.distance);
		if (err > max_err_hull) max_err_hull = err;

		// gjk_box path
		double err2 = fabs((double)gjk_distance_v(gjk_box(posA, rA, heA), gjk_box(posB, rB, heB), NULL).distance - ref.distance);
		if (err2 > max_err_box) max_err_box = err2;
	}
	printf("  bf box-box (%d): hull max_err=%.3e  box max_err=%.3e\n", BF_NCASES, max_err_hull, max_err_box);
	TEST_BEGIN("bf: box-box hull < tol");  TEST_ASSERT(max_err_hull < BF_TOL);
	TEST_BEGIN("bf: box-box box < tol");   TEST_ASSERT(max_err_box < BF_TOL);
}

static void test_gjk_bf_segment_hull()
{
	const Hull* ub = hull_unit_box();
	gjk_perf_rng = 999;
	double max_err = 0;

	for (int c = 0; c < BF_NCASES; c++) {
		float sep = 0.2f + gjk_perf_randf() * 8.0f;
		float sc = 0.3f + gjk_perf_randf() * 4.0f;
		quat rB = gjk_perf_random_quat();
		v3 segP = V3(gjk_perf_randf() - 0.5f, -sc, gjk_perf_randf() - 0.5f);
		v3 segQ = V3(gjk_perf_randf() - 0.5f, sc, gjk_perf_randf() - 0.5f);
		v3 boxHE = V3(sc, sc * (0.5f + gjk_perf_randf()), sc * (0.5f + gjk_perf_randf()));
		float ext = sqrtf(boxHE.x*boxHE.x + boxHE.y*boxHE.y + boxHE.z*boxHE.z);
		v3 boxPos = V3(sep + ext + sc, gjk_perf_randf() - 0.5f, 0);

		dv3 wb[8];
		bf_hull_world_verts(wb, ub, boxHE, dv3_from_v3(boxPos), dquat_from_quat(rB));
		BruteResult ref = bf_segment_hull(dv3_from_v3(segP), dv3_from_v3(segQ), wb, ub);

		v3 sb[8]; float sosb[24];
		GJK_Shape sa = gjk_capsule(segP, segQ, 0);
		GJK_Shape sbx = gjk_hull_scaled(ub, boxPos, rB, boxHE, sb, sosb);
		double err = fabs((double)gjk_distance_v(sa, sbx, NULL).distance - ref.distance);
		if (err > max_err) max_err = err;
	}
	printf("  bf seg-hull (%d): max_err=%.3e\n", BF_NCASES, max_err);
	TEST_BEGIN("bf: seg-hull < tol");  TEST_ASSERT(max_err < BF_TOL);
}

static void test_gjk_bf_capsule_box()
{
	const Hull* ub = hull_unit_box();
	gjk_perf_rng = 1234;
	double max_err = 0;

	for (int c = 0; c < BF_NCASES; c++) {
		float sep = 0.2f + gjk_perf_randf() * 8.0f;
		float sc = 0.3f + gjk_perf_randf() * 4.0f;
		float capR = sc * (0.1f + gjk_perf_randf() * 0.4f);
		quat rB = gjk_perf_random_quat();
		v3 capP = V3(0, -sc, 0), capQ = V3(0, sc, 0);
		v3 boxHE = V3(sc, sc * (0.5f + gjk_perf_randf()), sc * (0.5f + gjk_perf_randf()));
		float ext = sqrtf(boxHE.x*boxHE.x + boxHE.y*boxHE.y + boxHE.z*boxHE.z);
		v3 boxPos = V3(sep + ext + capR, gjk_perf_randf() - 0.5f, 0);

		// bf gives segment-hull distance; capsule distance = that - capR
		dv3 wb[8];
		bf_hull_world_verts(wb, ub, boxHE, dv3_from_v3(boxPos), dquat_from_quat(rB));
		BruteResult ref = bf_segment_hull(dv3_from_v3(capP), dv3_from_v3(capQ), wb, ub);
		double bf_dist = ref.distance - (double)capR;

		v3 sb[8]; float sosb[24];
		GJK_Shape sa = gjk_capsule(capP, capQ, capR);
		GJK_Shape sbx = gjk_hull_scaled(ub, boxPos, rB, boxHE, sb, sosb);
		double err = fabs((double)gjk_distance_v(sa, sbx, NULL).distance - bf_dist);
		if (err > max_err) max_err = err;

		// Also test capsule vs gjk_box
		double err2 = fabs((double)gjk_distance_v(sa, gjk_box(boxPos, rB, boxHE), NULL).distance - bf_dist);
		if (err2 > max_err) max_err = err2;
	}
	printf("  bf capsule-box (%d): max_err=%.3e\n", BF_NCASES, max_err);
	TEST_BEGIN("bf: capsule-box < tol");  TEST_ASSERT(max_err < BF_TOL);
}

static void test_gjk_bf_sphere_box()
{
	const Hull* ub = hull_unit_box();
	gjk_perf_rng = 5678;
	double max_err = 0;

	for (int c = 0; c < BF_NCASES; c++) {
		float sep = 0.2f + gjk_perf_randf() * 8.0f;
		float sc = 0.3f + gjk_perf_randf() * 4.0f;
		float sphR = sc * (0.1f + gjk_perf_randf() * 0.5f);
		quat rB = gjk_perf_random_quat();
		v3 sphC = V3(gjk_perf_randf() - 0.5f, gjk_perf_randf() - 0.5f, gjk_perf_randf() - 0.5f);
		v3 boxHE = V3(sc, sc * (0.5f + gjk_perf_randf()), sc * (0.5f + gjk_perf_randf()));
		float ext = sqrtf(boxHE.x*boxHE.x + boxHE.y*boxHE.y + boxHE.z*boxHE.z);
		v3 boxPos = V3(sep + ext + sphR + 1.0f, gjk_perf_randf() - 0.5f, gjk_perf_randf() - 0.5f);

		// bf: point-hull distance - sphR
		dv3 wb[8];
		bf_hull_world_verts(wb, ub, boxHE, dv3_from_v3(boxPos), dquat_from_quat(rB));
		dv3 cp = bf_closest_on_hull(dv3_from_v3(sphC), wb, ub);
		double bf_dist = dv3_len(dv3_sub(dv3_from_v3(sphC), cp)) - (double)sphR;

		GJK_Shape sa = gjk_sphere(sphC, sphR);
		v3 sb[8]; float sosb[24];
		GJK_Shape sbh = gjk_hull_scaled(ub, boxPos, rB, boxHE, sb, sosb);
		double err = fabs((double)gjk_distance_v(sa, sbh, NULL).distance - bf_dist);
		if (err > max_err) max_err = err;

		GJK_Shape sbb = gjk_box(boxPos, rB, boxHE);
		double err2 = fabs((double)gjk_distance_v(sa, sbb, NULL).distance - bf_dist);
		if (err2 > max_err) max_err = err2;
	}
	printf("  bf sphere-box (%d): max_err=%.3e\n", BF_NCASES, max_err);
	TEST_BEGIN("bf: sphere-box < tol");  TEST_ASSERT(max_err < BF_TOL);
}

// Box vs box at near-contact: sweep gap to verify gjk matches bf zero-crossing.
static void test_gjk_near_contact()
{
	const Hull* ub = hull_unit_box();
	quat id = quat_identity();
	int total_ok = 0;

	// Face-face
	{
		float gap = 0.1f;
		int agree = 1;
		while (gap > 1e-10f) {
			v3 posB = V3(2.0f + gap, 0, 0);
			dv3 wa[8], wb[8];
			bf_hull_world_verts(wa, ub, V3(1,1,1), DV3(0,0,0), dquat_from_quat(id));
			bf_hull_world_verts(wb, ub, V3(1,1,1), dv3_from_v3(posB), dquat_from_quat(id));
			BruteResult ref = bf_hull_hull(wa, ub, wb, ub);
			GJK_Result rg = gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(posB, id, V3(1,1,1)), NULL);
			if ((ref.distance == 0) != (rg.distance == 0)) { agree = 0; break; }
			if (ref.distance == 0 && rg.distance == 0) break;
			gap *= 0.8f;
		}
		TEST_BEGIN("near-contact: box face"); TEST_ASSERT(agree); total_ok += agree;
	}

	// Edge-edge
	{
		quat rot45z = { 0, 0, 0.3826834f, 0.9238795f };
		float contact = 1.0f + sqrtf(2.0f);
		float gap = 0.1f;
		int agree = 1;
		while (gap > 1e-10f) {
			v3 posB = V3(contact + gap, 0, 0);
			dv3 wa[8], wb[8];
			bf_hull_world_verts(wa, ub, V3(1,1,1), DV3(0,0,0), dquat_from_quat(id));
			bf_hull_world_verts(wb, ub, V3(1,1,1), dv3_from_v3(posB), dquat_from_quat(rot45z));
			BruteResult ref = bf_hull_hull(wa, ub, wb, ub);
			GJK_Result rg =gjk_distance_v(gjk_box(V3(0,0,0), id, V3(1,1,1)), gjk_box(posB, rot45z, V3(1,1,1)), NULL);
			if ((ref.distance == 0) != (rg.distance == 0)) { agree = 0; break; }
			if (ref.distance == 0 && rg.distance == 0) break;
			gap *= 0.8f;
		}
		TEST_BEGIN("near-contact: box edge"); TEST_ASSERT(agree); total_ok += agree;
	}
}

// =============================================================================
// Correctness: symmetry and consistency checks.

static void test_gjk_properties()
{
	gjk_perf_rng = 31415;
	int ncases = 200;

	// Symmetry: distance(A,B) == distance(B,A)
	for (int c = 0; c < ncases; c++) {
		float sc = 0.5f + gjk_perf_randf() * 5.0f;
		float sep = 1.0f + gjk_perf_randf() * 10.0f;
		quat rA = gjk_perf_random_quat(), rB = gjk_perf_random_quat();
		GJK_Shape a = gjk_box(V3(0,0,0), rA, V3(sc, sc*0.7f, sc*0.5f));
		GJK_Shape b = gjk_box(V3(sep, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f), rB, V3(sc*0.6f, sc, sc*0.8f));
		GJK_Result rab = gjk_distance_v(a, b, NULL);
		GJK_Result rba = gjk_distance_v(b, a, NULL);
		TEST_BEGIN("symmetry: box-box");
		TEST_ASSERT(fabsf(rab.distance - rba.distance) < 1e-4f);
	}

	// Symmetry: capsule
	gjk_perf_rng = 27182;
	for (int c = 0; c < ncases; c++) {
		float sc = 0.5f + gjk_perf_randf() * 3.0f;
		float sep = 1.0f + gjk_perf_randf() * 8.0f;
		GJK_Shape a = gjk_capsule(V3(0,-sc,0), V3(0,sc,0), sc*0.3f);
		GJK_Shape b = gjk_box(V3(sep, 0, 0), gjk_perf_random_quat(), V3(sc, sc*0.8f, sc*0.6f));
		GJK_Result rab = gjk_distance_v(a, b, NULL);
		GJK_Result rba = gjk_distance_v(b, a, NULL);
		TEST_BEGIN("symmetry: capsule-box");
		TEST_ASSERT(fabsf(rab.distance - rba.distance) < 1e-4f);
	}

	// Non-negativity
	gjk_perf_rng = 14142;
	for (int c = 0; c < ncases; c++) {
		float sc = 0.3f + gjk_perf_randf() * 5.0f;
		GJK_Shape a = gjk_box(V3(0,0,0), gjk_perf_random_quat(), V3(sc, sc*0.7f, sc*0.5f));
		GJK_Shape b = gjk_box(V3(gjk_perf_randf()*20-10, gjk_perf_randf()*20-10, gjk_perf_randf()*20-10), gjk_perf_random_quat(), V3(sc*0.6f, sc, sc*0.8f));
		GJK_Result r = gjk_distance_v(a, b, NULL);
		TEST_BEGIN("non-negative distance");
		TEST_ASSERT(r.distance >= 0.0f);
	}

	// Self-distance = 0 (overlapping with self)
	TEST_BEGIN("self-distance box");
	{
		GJK_Shape a = gjk_box(V3(1,2,3), gjk_perf_random_quat(), V3(1,1,1));
		GJK_Result r = gjk_distance_v(a, a, NULL);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);
	}

	TEST_BEGIN("self-distance capsule");
	{
		GJK_Shape a = gjk_capsule(V3(0,-1,0), V3(0,1,0), 0.5f);
		GJK_Result r = gjk_distance_v(a, a, NULL);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);
	}

	TEST_BEGIN("self-distance sphere");
	{
		GJK_Shape a = gjk_sphere(V3(3,4,5), 2.0f);
		GJK_Result r = gjk_distance_v(a, a, NULL);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);
	}

	// Witness points lie on/near shape surfaces
	gjk_perf_rng = 99999;
	for (int c = 0; c < ncases; c++) {
		float sc = 0.5f + gjk_perf_randf() * 3.0f;
		float sep = sc * 2.0f + 0.5f + gjk_perf_randf() * 5.0f;
		GJK_Shape a = gjk_box(V3(0,0,0), gjk_perf_random_quat(), V3(sc, sc*0.7f, sc*0.5f));
		GJK_Shape b = gjk_box(V3(sep, 0, 0), gjk_perf_random_quat(), V3(sc*0.6f, sc, sc*0.8f));
		GJK_Result r = gjk_distance_v(a, b, NULL);
		if (r.distance > 0.01f) {
			// witness distance should match reported distance
			float wd = len(sub(r.point2, r.point1));
			TEST_BEGIN("witness dist == reported dist");
			TEST_ASSERT(fabsf(wd - r.distance) < 1e-4f);
		}
	}
}

// =============================================================================
// Performance: shapes rotate and translate each iteration to simulate rigid body motion.

#define PERF_N       100000
#define PERF_CONFIGS 32
#define PERF_DT      (1.0f / 60.0f / 100.0f)  // sub-frame timestep

typedef struct PerfRow { double ns; float iters; } PerfRow;

// Integrate rotation: rot' = normalize(rot + 0.5 * quat(omega*dt, 0) * rot)
static quat perf_integrate_rot(quat rot, v3 omega, float dt)
{
	v3 half_w = scale(omega, 0.5f * dt);
	quat dq = { half_w.x * rot.w + half_w.y * rot.z - half_w.z * rot.y,
	            -half_w.x * rot.z + half_w.y * rot.w + half_w.z * rot.x,
	             half_w.x * rot.y - half_w.y * rot.x + half_w.z * rot.w,
	            -half_w.x * rot.x - half_w.y * rot.y - half_w.z * rot.z };
	quat r = { rot.x + dq.x, rot.y + dq.y, rot.z + dq.z, rot.w + dq.w };
	float inv_len = 1.0f / sqrtf(r.x*r.x + r.y*r.y + r.z*r.z + r.w*r.w);
	return (quat){ r.x*inv_len, r.y*inv_len, r.z*inv_len, r.w*inv_len };
}

// Per-config motion state
typedef struct PerfMotion
{
	v3 posA, posB;       // current positions
	v3 velA, velB;       // linear velocity
	quat rotA, rotB;     // current rotations
	v3 omegaA, omegaB;   // angular velocity (rad/s)
	v3 heA, heB;         // half-extents / scale
	float radiusA, radiusB; // for capsule/cylinder
} PerfMotion;

static PerfRow perf_box_box()
{
	gjk_perf_rng = 42;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.1f + (float)(c % 4) * 3.0f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.5f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc*0.7f, sc*0.5f); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			GJK_Shape a = gjk_box(mc->posA, mc->rotA, mc->heA), b = gjk_box(mc->posB, mc->rotB, mc->heB);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

static PerfRow perf_capsule_box()
{
	gjk_perf_rng = 137;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc, sc*0.8f, sc*0.6f);
		m[c].radiusA = sc * 0.3f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			v3 cap_dir = quat_rotate(mc->rotA, V3(0,1,0));
			v3 cap_p = sub(mc->posA, scale(cap_dir, mc->heA.y));
			v3 cap_q = add(mc->posA, scale(cap_dir, mc->heA.y));
			GJK_Shape a = gjk_capsule(cap_p, cap_q, mc->radiusA), b = gjk_box(mc->posB, mc->rotB, mc->heB);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

static PerfRow perf_capsule_capsule()
{
	gjk_perf_rng = 271;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc*0.8f, sc*0.8f, sc*0.8f);
		m[c].radiusA = sc * 0.3f; m[c].radiusB = sc * 0.4f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			v3 da = quat_rotate(mc->rotA, V3(0,1,0)), db = quat_rotate(mc->rotB, V3(0,1,0));
			GJK_Shape a = gjk_capsule(sub(mc->posA, scale(da, mc->heA.y)), add(mc->posA, scale(da, mc->heA.y)), mc->radiusA);
			GJK_Shape b = gjk_capsule(sub(mc->posB, scale(db, mc->heB.y)), add(mc->posB, scale(db, mc->heB.y)), mc->radiusB);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

static PerfRow perf_cylinder_box()
{
	gjk_perf_rng = 401;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc, sc*0.8f, sc*0.6f);
		m[c].radiusA = sc * 0.5f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			v3 ca = quat_rotate(mc->rotA, V3(0,1,0));
			v3 cp = sub(mc->posA, scale(ca, mc->heA.y));
			v3 cq = add(mc->posA, scale(ca, mc->heA.y));
			GJK_Shape a = gjk_cylinder(cp, cq, mc->radiusA), b = gjk_box(mc->posB, mc->rotB, mc->heB);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

static PerfRow perf_cylinder_cylinder()
{
	gjk_perf_rng = 523;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc*0.7f, sc*0.7f, sc*0.7f);
		m[c].radiusA = sc * 0.5f; m[c].radiusB = sc * 0.4f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			v3 da = quat_rotate(mc->rotA, V3(0,1,0)), db = quat_rotate(mc->rotB, V3(0,1,0));
			GJK_Shape a = gjk_cylinder(sub(mc->posA, scale(da, mc->heA.y)), add(mc->posA, scale(da, mc->heA.y)), mc->radiusA);
			GJK_Shape b = gjk_cylinder(sub(mc->posB, scale(db, mc->heB.y)), add(mc->posB, scale(db, mc->heB.y)), mc->radiusB);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// =============================================================================
// Entry point.

// Hull-hull perf: two random convex hulls of n_target verts each, rotating + translating.
static PerfRow perf_hull_hull(int n_target)
{
	gjk_perf_rng = 7000 + n_target;
	Hull* ha = gjk_perf_make_random_hull(n_target, 1.0f);
	Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);

	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*3, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc*0.7f, sc*0.5f); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int n_iters = PERF_N;
	if (n_target >= 200) n_iters = PERF_N / 10;
	if (n_target >= 1000) n_iters = PERF_N / 100;
	int total = n_iters * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < n_iters; i++) {
			GJK_Shape ga = gjk_hull_scaled(ha, mc->posA, mc->rotA, mc->heA, NULL, NULL);
			GJK_Shape gb = gjk_hull_scaled(hb, mc->posB, mc->rotB, mc->heB, NULL, NULL);
			GJK_Result r = gjk_distance(&ga, &gb, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	hull_free(ha); hull_free(hb);
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Hull-hull bf accuracy: random hulls of given size vs brute force.
static void test_gjk_bf_hull_hull(int n_target)
{
	gjk_perf_rng = 8000 + n_target;
	int ncases = 100;
	double max_err = 0;

	for (int c = 0; c < ncases; c++) {
		Hull* ha = gjk_perf_make_random_hull(n_target, 1.0f);
		Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);
		float sep = 0.5f + gjk_perf_randf() * 8.0f;
		float sc = 0.5f + gjk_perf_randf() * 3.0f;
		quat rA = gjk_perf_random_quat(), rB = gjk_perf_random_quat();
		v3 heA = V3(sc, sc*0.7f, sc*0.5f), heB = V3(sc*0.8f, sc, sc*0.6f);
		float extA = sqrtf(heA.x*heA.x + heA.y*heA.y + heA.z*heA.z);
		float extB = sqrtf(heB.x*heB.x + heB.y*heB.y + heB.z*heB.z);
		v3 posA = V3(0,0,0), posB = V3(sep + extA + extB, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);

		dv3 wa[1024], wb[1024];
		bf_hull_world_verts(wa, ha, heA, dv3_from_v3(posA), dquat_from_quat(rA));
		bf_hull_world_verts(wb, hb, heB, dv3_from_v3(posB), dquat_from_quat(rB));
		BruteResult ref = bf_hull_hull(wa, ha, wb, hb);

		v3 sa[1024], sb[1024]; float soa[3072], sob[3072];
		GJK_Shape ga = gjk_hull_scaled(ha, posA, rA, heA, sa, soa);
		GJK_Shape gb = gjk_hull_scaled(hb, posB, rB, heB, sb, sob);
		double err = fabs((double)gjk_distance_v(ga, gb, NULL).distance - ref.distance);
		if (err > max_err) max_err = err;

		hull_free(ha); hull_free(hb);
	}
	printf("  bf hull-%dv (%d): max_err=%.3e\n", n_target, ncases, max_err);
	char name[64]; snprintf(name, sizeof(name), "bf: hull-%dv < tol", n_target);
	TEST_BEGIN(name); TEST_ASSERT(max_err < BF_TOL);
}

// Sphere-hull perf: sphere vs random convex hull.
static PerfRow perf_sphere_hull(int n_target)
{
	gjk_perf_rng = 9000 + n_target;
	Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*3, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
		m[c].radiusA = sc * 0.5f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int n_iters = PERF_N;
	if (n_target >= 200) n_iters = PERF_N / 10;
	if (n_target >= 1000) n_iters = PERF_N / 100;
	int total = n_iters * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < n_iters; i++) {
			GJK_Shape gb = gjk_hull_scaled(hb, mc->posB, mc->rotB, mc->heB, NULL, NULL);
			GJK_Shape a = gjk_sphere(mc->posA, mc->radiusA);
			GJK_Result r = gjk_distance(&a, &gb, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	hull_free(hb);
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Capsule-hull perf.
static PerfRow perf_capsule_hull(int n_target)
{
	gjk_perf_rng = 9500 + n_target;
	Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*3, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
		m[c].radiusA = sc * 0.3f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int n_iters = PERF_N;
	if (n_target >= 200) n_iters = PERF_N / 10;
	if (n_target >= 1000) n_iters = PERF_N / 100;
	int total = n_iters * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < n_iters; i++) {
			v3 cap_dir = quat_rotate(mc->rotA, V3(0,1,0));
			v3 cap_p = sub(mc->posA, scale(cap_dir, mc->heA.y));
			v3 cap_q = add(mc->posA, scale(cap_dir, mc->heA.y));
			GJK_Shape gb = gjk_hull_scaled(hb, mc->posB, mc->rotB, mc->heB, NULL, NULL);
			GJK_Shape a = gjk_capsule(cap_p, cap_q, mc->radiusA);
			GJK_Result r = gjk_distance(&a, &gb, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	hull_free(hb);
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Box-hull perf.
static PerfRow perf_box_hull(int n_target)
{
	gjk_perf_rng = 10000 + n_target;
	Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*3, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc*0.7f, sc*0.5f); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int n_iters = PERF_N;
	if (n_target >= 200) n_iters = PERF_N / 10;
	if (n_target >= 1000) n_iters = PERF_N / 100;
	int total = n_iters * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < n_iters; i++) {
			GJK_Shape gb = gjk_hull_scaled(hb, mc->posB, mc->rotB, mc->heB, NULL, NULL);
			GJK_Shape a = gjk_box(mc->posA, mc->rotA, mc->heA);
			GJK_Result r = gjk_distance(&a, &gb, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	hull_free(hb);
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Cylinder-hull perf.
static PerfRow perf_cylinder_hull(int n_target)
{
	gjk_perf_rng = 10500 + n_target;
	Hull* hb = gjk_perf_make_random_hull(n_target, 1.0f);
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*3, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc, sc); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
		m[c].radiusA = sc * 0.5f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int n_iters = PERF_N;
	if (n_target >= 200) n_iters = PERF_N / 10;
	if (n_target >= 1000) n_iters = PERF_N / 100;
	int total = n_iters * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < n_iters; i++) {
			v3 ca = quat_rotate(mc->rotA, V3(0,1,0));
			v3 cp = sub(mc->posA, scale(ca, mc->heA.y));
			v3 cq = add(mc->posA, scale(ca, mc->heA.y));
			GJK_Shape gb = gjk_hull_scaled(hb, mc->posB, mc->rotB, mc->heB, NULL, NULL);
			GJK_Shape a = gjk_cylinder(cp, cq, mc->radiusA);
			GJK_Result r = gjk_distance(&a, &gb, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	hull_free(hb);
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Box-triangle perf: box vs rotating triangle.
static PerfRow perf_box_triangle()
{
	gjk_perf_rng = 11111;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotA = gjk_perf_random_quat(); m[c].rotB = gjk_perf_random_quat();
		m[c].omegaA = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heA = V3(sc, sc*0.7f, sc*0.5f); m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			GJK_Shape a = gjk_box(mc->posA, mc->rotA, mc->heA);
			// Triangle: 3 vertices rotated around posB
			v3 tv0 = add(mc->posB, quat_rotate(mc->rotB, V3(-mc->heB.x, 0, -mc->heB.z)));
			v3 tv1 = add(mc->posB, quat_rotate(mc->rotB, V3(mc->heB.x, 0, -mc->heB.z)));
			v3 tv2 = add(mc->posB, quat_rotate(mc->rotB, V3(0, 0, mc->heB.z)));
			GJK_Shape b = gjk_triangle(tv0, tv1, tv2);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotA = perf_integrate_rot(mc->rotA, mc->omegaA, PERF_DT);
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

// Sphere-triangle perf.
static PerfRow perf_sphere_triangle()
{
	gjk_perf_rng = 11222;
	PerfMotion m[PERF_CONFIGS];
	for (int c = 0; c < PERF_CONFIGS; c++) {
		float sep = 0.5f + (float)(c % 4) * 2.5f;
		float sc = 0.5f + (float)(c / 4 % 4) * 2.0f;
		m[c].posA = V3(0,0,0); m[c].posB = V3(sep + sc*2, 0, 0);
		m[c].velA = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].velB = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
		m[c].rotB = gjk_perf_random_quat();
		m[c].omegaB = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
		m[c].heB = V3(sc*0.8f, sc, sc*0.6f);
		m[c].radiusA = sc * 0.5f;
	}
	GJK_Cache cache[PERF_CONFIGS] = {0};
	int total = PERF_N * PERF_CONFIGS;
	float sum = 0; int iters = 0;
	double t0 = qpc_now();
	for (int c = 0; c < PERF_CONFIGS; c++) {
		PerfMotion* mc = &m[c];
		for (int i = 0; i < PERF_N; i++) {
			v3 tv0 = add(mc->posB, quat_rotate(mc->rotB, V3(-mc->heB.x, 0, -mc->heB.z)));
			v3 tv1 = add(mc->posB, quat_rotate(mc->rotB, V3(mc->heB.x, 0, -mc->heB.z)));
			v3 tv2 = add(mc->posB, quat_rotate(mc->rotB, V3(0, 0, mc->heB.z)));
			GJK_Shape a = gjk_sphere(mc->posA, mc->radiusA), b = gjk_triangle(tv0, tv1, tv2);
			GJK_Result r = gjk_distance(&a, &b, &cache[c]);
			sum += r.distance; iters += r.iterations;
			mc->posA = add(mc->posA, scale(mc->velA, PERF_DT));
			mc->posB = add(mc->posB, scale(mc->velB, PERF_DT));
			mc->rotB = perf_integrate_rot(mc->rotB, mc->omegaB, PERF_DT);
		}
	}
	double t1 = qpc_now();
	volatile float sink = sum; (void)sink;
	return (PerfRow){ (t1 - t0) * 1e9 / total, (float)iters / total };
}

static void run_gjk_perf_tests()
{
	qpc_init();

	// Correctness: known distances (all shape pairs)
	test_gjk_known_distances();

	// Brute-force accuracy (slow -- uncomment when needed)
	// printf("\n  --- brute-force accuracy ---\n");
	// test_gjk_bf_box_box();
	// test_gjk_bf_segment_hull();
	// test_gjk_bf_capsule_box();
	// test_gjk_bf_sphere_box();

	// Correctness: near-contact zero-crossing
	test_gjk_near_contact();

	// Correctness: symmetry, non-negativity, witness points
	test_gjk_properties();

	// Performance
	int hull_sizes[] = { 20, 50, 200, 1000 };
	int nh = 4;
	enum { N_PERF = 5 + 5*4 + 2 }; // 5 prim-prim + 5 shape types * 4 hull sizes + 2 triangle
	const char* names[N_PERF];
	PerfRow p[N_PERF];
	char name_bufs[N_PERF][20];
	int pi = 0;
	names[pi] = "box-box";       p[pi++] = perf_box_box();
	names[pi] = "capsule-box";   p[pi++] = perf_capsule_box();
	names[pi] = "capsule-cap";   p[pi++] = perf_capsule_capsule();
	names[pi] = "cylinder-box";  p[pi++] = perf_cylinder_box();
	names[pi] = "cylinder-cyl";  p[pi++] = perf_cylinder_cylinder();
	for (int h = 0; h < nh; h++) {
		int n = hull_sizes[h];
		sprintf(name_bufs[pi], "sph-h%d", n);    names[pi] = name_bufs[pi]; p[pi] = perf_sphere_hull(n);   pi++;
		sprintf(name_bufs[pi], "cap-h%d", n);    names[pi] = name_bufs[pi]; p[pi] = perf_capsule_hull(n);  pi++;
		sprintf(name_bufs[pi], "box-h%d", n);    names[pi] = name_bufs[pi]; p[pi] = perf_box_hull(n);      pi++;
		sprintf(name_bufs[pi], "cyl-h%d", n);    names[pi] = name_bufs[pi]; p[pi] = perf_cylinder_hull(n); pi++;
		sprintf(name_bufs[pi], "hull-h%d", n);   names[pi] = name_bufs[pi]; p[pi] = perf_hull_hull(n);     pi++;
	}

	// Triangle tests
	names[pi] = "sph-tri";       p[pi++] = perf_sphere_triangle();
	names[pi] = "box-tri";       p[pi++] = perf_box_triangle();

	printf("\n=== GJK Performance (%d iters x %d configs) ===\n\n", PERF_N, PERF_CONFIGS);
	printf("  shape pair        ns/call  iters\n");
	printf("  --------------- --------  -----\n");
	for (int i = 0; i < pi; i++)
		printf("  %-15s %8.1f  %5.1f\n", names[i], p[i].ns, p[i].iters);

	// Batch box-box: 4 pairs at a time
	{
		gjk_perf_rng = 77777;
		int batch_n = PERF_N;
		int batch_configs = PERF_CONFIGS;
		GJK_Shape shapes_a[32], shapes_b[32];
		GJK_Cache caches[32] = {0};
		v3 posA[32], posB[32], velA[32], velB[32];
		quat rotA[32], rotB[32];
		v3 omA[32], omB[32], heA[32], heB[32];
		for (int c = 0; c < batch_configs; c++) {
			float sep = 0.1f + (float)(c % 4) * 3.0f;
			float sc = 0.5f + (float)(c / 4 % 4) * 2.5f;
			posA[c] = V3(0,0,0); posB[c] = V3(sep + sc*2, 0, 0);
			velA[c] = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
			velB[c] = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
			rotA[c] = gjk_perf_random_quat(); rotB[c] = gjk_perf_random_quat();
			omA[c] = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
			omB[c] = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
			heA[c] = V3(sc, sc*0.7f, sc*0.5f); heB[c] = V3(sc*0.8f, sc, sc*0.6f);
		}
		float sum = 0;
		double t0 = qpc_now();
		for (int i = 0; i < batch_n; i++) {
			// Build shapes for all configs
			for (int c = 0; c < batch_configs; c++) {
				shapes_a[c] = gjk_box(posA[c], rotA[c], heA[c]);
				shapes_b[c] = gjk_box(posB[c], rotB[c], heB[c]);
			}
			// Process in batches of 4
			for (int c = 0; c < batch_configs; c += 4) {
				const GJK_Shape* pa[4] = { &shapes_a[c], &shapes_a[c+1], &shapes_a[c+2], &shapes_a[c+3] };
				const GJK_Shape* pb[4] = { &shapes_b[c], &shapes_b[c+1], &shapes_b[c+2], &shapes_b[c+3] };
				GJK_Cache* pc[4] = { &caches[c], &caches[c+1], &caches[c+2], &caches[c+3] };
				float dist[4];
				gjk_distance_batch_box(pa, pb, pc, dist);
				sum += dist[0] + dist[1] + dist[2] + dist[3];
			}
			// Integrate
			for (int c = 0; c < batch_configs; c++) {
				posA[c] = add(posA[c], scale(velA[c], PERF_DT));
				posB[c] = add(posB[c], scale(velB[c], PERF_DT));
				rotA[c] = perf_integrate_rot(rotA[c], omA[c], PERF_DT);
				rotB[c] = perf_integrate_rot(rotB[c], omB[c], PERF_DT);
			}
		}
		double t1 = qpc_now();
		volatile float sink = sum; (void)sink;
		int total = batch_n * batch_configs;
		double batch_ns = (t1 - t0) * 1e9 / total;

		// Compare with scalar
		for (int c = 0; c < batch_configs; c++) caches[c] = (GJK_Cache){0};
		gjk_perf_rng = 77777;
		for (int c = 0; c < batch_configs; c++) {
			float sep = 0.1f + (float)(c % 4) * 3.0f;
			float sc = 0.5f + (float)(c / 4 % 4) * 2.5f;
			posA[c] = V3(0,0,0); posB[c] = V3(sep + sc*2, 0, 0);
			velA[c] = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
			velB[c] = V3(gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f, gjk_perf_randf()-0.5f);
			rotA[c] = gjk_perf_random_quat(); rotB[c] = gjk_perf_random_quat();
			omA[c] = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
			omB[c] = V3((gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4, (gjk_perf_randf()-0.5f)*4);
			heA[c] = V3(sc, sc*0.7f, sc*0.5f); heB[c] = V3(sc*0.8f, sc, sc*0.6f);
		}
		sum = 0;
		t0 = qpc_now();
		for (int i = 0; i < batch_n; i++) {
			for (int c = 0; c < batch_configs; c++) {
				GJK_Shape a = gjk_box(posA[c], rotA[c], heA[c]), b = gjk_box(posB[c], rotB[c], heB[c]);
				GJK_Result r = gjk_distance(&a, &b, &caches[c]);
				sum += r.distance;
			}
			for (int c = 0; c < batch_configs; c++) {
				posA[c] = add(posA[c], scale(velA[c], PERF_DT));
				posB[c] = add(posB[c], scale(velB[c], PERF_DT));
				rotA[c] = perf_integrate_rot(rotA[c], omA[c], PERF_DT);
				rotB[c] = perf_integrate_rot(rotB[c], omB[c], PERF_DT);
			}
		}
		t1 = qpc_now();
		volatile float sink2 = sum; (void)sink2;
		double scalar_ns = (t1 - t0) * 1e9 / total;
		printf("\n  --- batch box-box (4-wide SIMD) ---\n");
		printf("  scalar:  %6.1f ns/pair\n", scalar_ns);
		printf("  batch:   %6.1f ns/pair (%.1fx)\n", batch_ns, scalar_ns / batch_ns);
	}
	printf("\n");
}
