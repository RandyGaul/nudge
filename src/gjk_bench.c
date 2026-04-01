// See LICENSE for licensing info.
// gjk_bench.c -- side-by-side comparison of index-based vs epsilon-based GJK.
//
// Sweeps shape pairs from separated to near-touching with gentle rotation,
// recording iteration counts and distances from both implementations.
// Output is CSV for easy plotting.

#define CKIT_IMPLEMENTATION
#include "ckit.h"

#include "nudge.h"
#include "nudge.c"

#include <stdio.h>

// -----------------------------------------------------------------------------
// Config.

#define BENCH_PI 3.14159265f
#define SWEEP_STEPS 200
#define TOTAL_ROTATION 0.3f // radians over full sweep -- gentle tumble

// -----------------------------------------------------------------------------
// Helpers.

static quat quat_axis_angle(v3 axis, float angle)
{
	float half = angle * 0.5f;
	float s = sinf(half);
	return (quat){ axis.x * s, axis.y * s, axis.z * s, cosf(half) };
}

// -----------------------------------------------------------------------------
// Shape wrapper for both GJK adapters.

typedef struct BenchHull
{
	const Hull* hull;
	v3 pos;
	quat rot;
} BenchHull;

// Index-based support (for gjk.c). Must inline the loop to capture vertex index.
static v3 bench_support_index(const void* shape, v3 dir, int* out_index)
{
	const BenchHull* bh = shape;
	v3 local_dir = rotate(inv(bh->rot), dir);
	float best = -1e18f;
	int best_i = 0;
	for (int i = 0; i < bh->hull->vert_count; i++) {
		float d = dot(bh->hull->verts[i], local_dir);
		if (d > best) { best = d; best_i = i; }
	}
	*out_index = best_i;
	return add(bh->pos, rotate(bh->rot, bh->hull->verts[best_i]));
}

// Epsilon-based support (for gjk_gino.c). No index needed.
static v3 bench_support_gino(const void* shape, v3 dir)
{
	const BenchHull* bh = shape;
	v3 local_dir = rotate(inv(bh->rot), dir);
	v3 local_sup = hull_support(bh->hull, local_dir);
	return add(bh->pos, rotate(bh->rot, local_sup));
}

// -----------------------------------------------------------------------------
// Shape generators. All return quickhull-allocated Hulls.

static Hull* make_box_hull(v3 he)
{
	v3 pts[8] = {
		{-he.x, -he.y, -he.z}, { he.x, -he.y, -he.z},
		{ he.x,  he.y, -he.z}, {-he.x,  he.y, -he.z},
		{-he.x, -he.y,  he.z}, { he.x, -he.y,  he.z},
		{ he.x,  he.y,  he.z}, {-he.x,  he.y,  he.z},
	};
	return quickhull(pts, 8);
}

static Hull* make_cylinder_hull(float radius, float half_height, int segments)
{
	int count = segments * 2;
	v3* pts = malloc(count * sizeof(v3));
	for (int i = 0; i < segments; i++) {
		float angle = 2.0f * BENCH_PI * (float)i / (float)segments;
		float cx = radius * cosf(angle);
		float cz = radius * sinf(angle);
		pts[i] = V3(cx, half_height, cz);
		pts[i + segments] = V3(cx, -half_height, cz);
	}
	Hull* h = quickhull(pts, count);
	free(pts);
	return h;
}

// Choose cylinder mesh fidelity based on radius.
static int cylinder_segments(float radius)
{
	if (radius < 0.5f) return 12;
	if (radius < 2.0f) return 24;
	return 32;
}

static Hull* make_icosphere_hull(float radius)
{
	// Icosahedron base vertices (12 points).
	float phi = (1.0f + sqrtf(5.0f)) * 0.5f;
	float a = radius / sqrtf(1.0f + phi * phi);
	float b = a * phi;
	v3 pts[42]; // 12 base + up to 30 edge midpoints
	pts[0]  = V3(-a,  b,  0); pts[1]  = V3( a,  b,  0);
	pts[2]  = V3(-a, -b,  0); pts[3]  = V3( a, -b,  0);
	pts[4]  = V3( 0, -a,  b); pts[5]  = V3( 0,  a,  b);
	pts[6]  = V3( 0, -a, -b); pts[7]  = V3( 0,  a, -b);
	pts[8]  = V3( b,  0, -a); pts[9]  = V3( b,  0,  a);
	pts[10] = V3(-b,  0, -a); pts[11] = V3(-b,  0,  a);

	// Add edge midpoints projected onto sphere for more vertices.
	// Icosahedron edges (30 undirected edges).
	static const int edges[][2] = {
		{0,1},{0,5},{0,7},{0,10},{0,11},{1,5},{1,7},{1,8},{1,9},
		{2,3},{2,4},{2,6},{2,10},{2,11},{3,4},{3,6},{3,8},{3,9},
		{4,5},{4,9},{4,11},{5,9},{5,11},{6,7},{6,8},{6,10},
		{7,8},{7,10},{8,9},{9,11},// should be {10,11} for last
	};
	int np = 12;
	for (int i = 0; i < 30 && np < 42; i++) {
		v3 mid = scale(add(pts[edges[i][0]], pts[edges[i][1]]), 0.5f);
		float l = len(mid);
		if (l > 1e-6f) mid = scale(mid, radius / l);
		pts[np++] = mid;
	}

	return quickhull(pts, np);
}

static Hull* make_capsule_hull(float half_height, float radius, int slices, int rings)
{
	// Generate hemisphere points at each end + equator rings.
	int cap_pts = slices * rings;
	int count = cap_pts * 2;
	v3* pts = malloc(count * sizeof(v3));
	int n = 0;
	for (int h = 0; h < 2; h++) {
		float sign = h == 0 ? 1.0f : -1.0f;
		float y_off = sign * half_height;
		for (int r = 0; r < rings; r++) {
			float lat = (BENCH_PI * 0.5f) * (float)(r + 1) / (float)rings;
			float ry = radius * cosf(lat);
			float rr = radius * sinf(lat);
			for (int s = 0; s < slices; s++) {
				float lon = 2.0f * BENCH_PI * (float)s / (float)slices;
				pts[n++] = V3(rr * cosf(lon), y_off + sign * ry, rr * sinf(lon));
			}
		}
	}
	Hull* hull = quickhull(pts, n);
	free(pts);
	return hull;
}

// -----------------------------------------------------------------------------
// Sweep runner.
//
// Sweeps two hulls along `dir` from `start_sep` to `end_sep` (center-to-center
// distance). Shape B gets a gentle rotation around `rot_axis` over the sweep.
// Negative end_sep means overlap.

static void run_sweep(const char* label, const Hull* hull_a, const Hull* hull_b, v3 dir, float start_sep, float end_sep, v3 rot_axis)
{
	v3 nd = norm(dir);
	v3 nrot = len2(rot_axis) > 1e-6f ? norm(rot_axis) : V3(0, 1, 0);

	for (int step = 0; step < SWEEP_STEPS; step++) {
		float t = (float)step / (float)(SWEEP_STEPS - 1);
		float sep = start_sep + (end_sep - start_sep) * t;
		float angle = TOTAL_ROTATION * t;

		v3 pos_a = scale(nd, -sep * 0.5f);
		v3 pos_b = scale(nd, sep * 0.5f);
		quat rot_b = quat_axis_angle(nrot, angle);

		BenchHull bh_a = { hull_a, pos_a, quat_identity() };
		BenchHull bh_b = { hull_b, pos_b, rot_b };

		GjkResult ri = gjk_distance(bench_support_index, &bh_a, bench_support_index, &bh_b);
		GinoResult rg = gino_gjk_distance(bench_support_gino, &bh_a, bench_support_gino, &bh_b);

		printf("%s,%d,%.6f,%d,%.6f,%d\n", label, step, ri.distance, ri.iterations, rg.distance, rg.iterations);
	}
}

// -----------------------------------------------------------------------------
// Benchmark scenarios.

static void bench_box_box()
{
	Hull* a = make_box_hull(V3(1, 1, 1));
	Hull* b = make_box_hull(V3(1, 1, 1));
	run_sweep("box-box", a, b, V3(1, 0, 0), 6.0f, 1.8f, V3(0, 1, 0));
	hull_free(a);
	hull_free(b);
}

static void bench_sphere_sphere()
{
	Hull* a = make_icosphere_hull(1.0f);
	Hull* b = make_icosphere_hull(1.0f);
	run_sweep("sphere-sphere", a, b, V3(1, 0.2f, 0), 6.0f, 1.8f, V3(0.3f, 1, 0.1f));
	hull_free(a);
	hull_free(b);
}

static void bench_capsule_capsule()
{
	Hull* a = make_capsule_hull(0.5f, 0.3f, 12, 4);
	Hull* b = make_capsule_hull(0.5f, 0.3f, 12, 4);
	run_sweep("capsule-capsule", a, b, V3(1, 0, 0.1f), 5.0f, 0.4f, V3(0, 0, 1));
	hull_free(a);
	hull_free(b);
}

static void bench_cylinder_cylinder()
{
	float r = 1.0f, hh = 1.0f;
	int segs = cylinder_segments(r);
	Hull* a = make_cylinder_hull(r, hh, segs);
	Hull* b = make_cylinder_hull(r, hh, segs);
	run_sweep("cylinder-cylinder", a, b, V3(1, 0, 0), 6.0f, 1.8f, V3(0, 1, 0));
	hull_free(a);
	hull_free(b);
}

static void bench_cylinder_box()
{
	float r = 1.0f, hh = 1.0f;
	int segs = cylinder_segments(r);
	Hull* a = make_cylinder_hull(r, hh, segs);
	Hull* b = make_box_hull(V3(1, 1, 1));
	run_sweep("cylinder-box", a, b, V3(1, 0.1f, 0), 6.0f, 1.5f, V3(0.2f, 1, 0));
	hull_free(a);
	hull_free(b);
}

static void bench_cylinder_sphere()
{
	float r = 1.0f, hh = 1.0f;
	int segs = cylinder_segments(r);
	Hull* a = make_cylinder_hull(r, hh, segs);
	Hull* b = make_icosphere_hull(1.0f);
	run_sweep("cylinder-sphere", a, b, V3(1, 0, 0.15f), 6.0f, 1.5f, V3(0, 1, 0.3f));
	hull_free(a);
	hull_free(b);
}

// -----------------------------------------------------------------------------

int main()
{
	printf("pair,step,dist_index,iters_index,dist_gino,iters_gino\n");
	bench_box_box();
	bench_sphere_sphere();
	bench_capsule_capsule();
	bench_cylinder_cylinder();
	bench_cylinder_box();
	bench_cylinder_sphere();
	return 0;
}
