// See LICENSE for licensing info.
// gjk_gino.c -- Gino van den Bergen style GJK with epsilon-based termination.
//
// For comparison against the index-based GJK in gjk.c.
// This version terminates when the support point progress falls below a
// relative epsilon threshold, which is the correct approach for shapes
// with continuous support functions (cylinders, spheres, capsules).
//
// Reference: "A Fast and Robust GJK Implementation for Collision Detection
// of Convex Objects" (Gino van den Bergen, 1999).

// -----------------------------------------------------------------------------
// Types.

typedef struct GinoVertex
{
	v3 point1;   // support point on shape A (world space)
	v3 point2;   // support point on shape B (world space)
	v3 point;    // Minkowski difference: point2 - point1
	float u;     // barycentric coordinate
} GinoVertex;

typedef struct GinoSimplex
{
	GinoVertex v[4];
	float divisor;
	int count;
} GinoSimplex;

typedef struct GinoResult
{
	v3 point1;     // closest point on shape A
	v3 point2;     // closest point on shape B
	float distance;
	int iterations;
} GinoResult;

// Support function signature: no index output (not needed for epsilon termination).
typedef v3 (*GinoSupportFn)(const void* shape, v3 dir);

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.
// Same Voronoi region approach as gjk.c, just using GinoVertex.

static void gino_solve2(GinoSimplex* s)
{
	v3 A = s->v[0].point;
	v3 B = s->v[1].point;

	float u = dot(sub(V3(0,0,0), B), sub(A, B));
	float v = dot(sub(V3(0,0,0), A), sub(B, A));

	if (v <= 0.0f) {
		s->v[0].u = 1.0f;
		s->divisor = 1.0f;
		s->count = 1;
		return;
	}
	if (u <= 0.0f) {
		s->v[0] = s->v[1];
		s->v[0].u = 1.0f;
		s->divisor = 1.0f;
		s->count = 1;
		return;
	}
	s->v[0].u = u;
	s->v[1].u = v;
	v3 e = sub(B, A);
	s->divisor = dot(e, e);
	s->count = 2;
}

static void gino_solve3(GinoSimplex* s)
{
	v3 A = s->v[0].point;
	v3 B = s->v[1].point;
	v3 C = s->v[2].point;
	v3 Q = V3(0,0,0);

	float uAB = dot(sub(Q, B), sub(A, B));
	float vAB = dot(sub(Q, A), sub(B, A));
	float uBC = dot(sub(Q, C), sub(B, C));
	float vBC = dot(sub(Q, B), sub(C, B));
	float uCA = dot(sub(Q, A), sub(C, A));
	float vCA = dot(sub(Q, C), sub(A, C));

	if (vAB <= 0.0f && uCA <= 0.0f) {
		s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}
	if (uAB <= 0.0f && vBC <= 0.0f) {
		s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}
	if (uBC <= 0.0f && vCA <= 0.0f) {
		s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}

	v3 n = cross(sub(B, A), sub(C, A));
	float uABC = dot(cross(sub(B, Q), sub(C, Q)), n);
	float vABC = dot(cross(sub(C, Q), sub(A, Q)), n);
	float wABC = dot(cross(sub(A, Q), sub(B, Q)), n);

	if (uAB > 0.0f && vAB > 0.0f && wABC <= 0.0f) {
		s->v[0].u = uAB; s->v[1].u = vAB;
		s->divisor = dot(sub(B, A), sub(B, A));
		s->count = 2; return;
	}
	if (uBC > 0.0f && vBC > 0.0f && uABC <= 0.0f) {
		s->v[0] = s->v[1]; s->v[1] = s->v[2];
		s->v[0].u = uBC; s->v[1].u = vBC;
		s->divisor = dot(sub(C, B), sub(C, B));
		s->count = 2; return;
	}
	if (uCA > 0.0f && vCA > 0.0f && vABC <= 0.0f) {
		s->v[1] = s->v[0]; s->v[0] = s->v[2];
		s->v[0].u = uCA; s->v[1].u = vCA;
		s->divisor = dot(sub(A, C), sub(A, C));
		s->count = 2; return;
	}

	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC;
	s->divisor = dot(n, n);
	s->count = 3;
}

static void gino_solve4(GinoSimplex* s)
{
	v3 A = s->v[0].point;
	v3 B = s->v[1].point;
	v3 C = s->v[2].point;
	v3 D = s->v[3].point;
	v3 Q = V3(0,0,0);

	float uABCD = triple(sub(B, Q), sub(C, Q), sub(D, Q));
	float vABCD = triple(sub(A, Q), sub(D, Q), sub(C, Q));
	float wABCD = triple(sub(A, Q), sub(B, Q), sub(D, Q));
	float xABCD = triple(sub(A, Q), sub(C, Q), sub(B, Q));

	if (uABCD > 0.0f && vABCD > 0.0f && wABCD > 0.0f && xABCD > 0.0f) {
		s->v[0].u = uABCD; s->v[1].u = vABCD;
		s->v[2].u = wABCD; s->v[3].u = xABCD;
		s->divisor = uABCD + vABCD + wABCD + xABCD;
		s->count = 4;
		return;
	}

	float best_dist = 1e18f;
	GinoSimplex best = *s;

	if (uABCD <= 0.0f) {
		GinoSimplex test = *s;
		test.v[0] = s->v[1]; test.v[1] = s->v[2]; test.v[2] = s->v[3];
		test.count = 3;
		gino_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	if (vABCD <= 0.0f) {
		GinoSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[2]; test.v[2] = s->v[3];
		test.count = 3;
		gino_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	if (wABCD <= 0.0f) {
		GinoSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[1]; test.v[2] = s->v[3];
		test.count = 3;
		gino_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	if (xABCD <= 0.0f) {
		GinoSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[1]; test.v[2] = s->v[2];
		test.count = 3;
		gino_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	*s = best;
}

// -----------------------------------------------------------------------------
// Closest point on the current simplex.

static v3 gino_closest_point(const GinoSimplex* s)
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

static void gino_witness_points(const GinoSimplex* s, v3* p1, v3* p2)
{
	float inv = 1.0f / s->divisor;
	*p1 = V3(0,0,0); *p2 = V3(0,0,0);
	for (int i = 0; i < s->count; i++) {
		*p1 = add(*p1, scale(s->v[i].point1, s->v[i].u * inv));
		*p2 = add(*p2, scale(s->v[i].point2, s->v[i].u * inv));
	}
}

// -----------------------------------------------------------------------------
// Main GJK distance function -- Gino van den Bergen style.
//
// Termination: the support point progress test.
// Let v = closest point on simplex, w = new support point.
// Terminate when: dot(v,v) - dot(v,w) <= eps_rel^2 * dot(v,v)
// This measures relative improvement -- when progress is negligible, stop.

#define GINO_MAX_ITERS 32
#define GINO_EPS_REL   1e-4f  // relative tolerance
#define GINO_EPS_ABS2  1e-14f // absolute tolerance squared (near-overlap)

static GinoResult gino_gjk_distance(GinoSupportFn support_a, const void* shape_a, GinoSupportFn support_b, const void* shape_b)
{
	GinoResult result = {0};

	GinoSimplex simplex = {0};
	v3 d = V3(1, 0, 0);

	// First support point
	simplex.v[0].point1 = support_a(shape_a, neg(d));
	simplex.v[0].point2 = support_b(shape_b, d);
	simplex.v[0].point = sub(simplex.v[0].point2, simplex.v[0].point1);
	simplex.v[0].u = 1.0f;
	simplex.divisor = 1.0f;
	simplex.count = 1;

	int iter = 0;
	float eps_rel2 = GINO_EPS_REL * GINO_EPS_REL;

	while (iter < GINO_MAX_ITERS) {
		// Solve for closest point on simplex to origin
		switch (simplex.count) {
		case 1: break;
		case 2: gino_solve2(&simplex); break;
		case 3: gino_solve3(&simplex); break;
		case 4: gino_solve4(&simplex); break;
		}

		// Origin inside tetrahedron = shapes overlap
		if (simplex.count == 4) break;

		// v = closest point on simplex to origin
		v3 v = gino_closest_point(&simplex);
		float v_dot_v = len2(v);

		// Near-zero distance: shapes are (nearly) touching/overlapping
		if (v_dot_v < GINO_EPS_ABS2) break;

		// Search direction: toward origin from closest point
		d = neg(v);

		// New support point w in direction -v
		GinoVertex* vertex = &simplex.v[simplex.count];
		vertex->point1 = support_a(shape_a, neg(d));
		vertex->point2 = support_b(shape_b, d);
		vertex->point = sub(vertex->point2, vertex->point1);
		iter++;

		// Gino's convergence test:
		// The maximum possible progress is dot(v, v) - dot(v, w).
		// If this is small relative to dot(v, v), we've converged.
		float v_dot_w = dot(v, vertex->point);
		float progress = v_dot_v - v_dot_w;
		if (progress <= eps_rel2 * v_dot_v) break;

		simplex.count++;
	}

	gino_witness_points(&simplex, &result.point1, &result.point2);
	result.distance = len(sub(result.point2, result.point1));
	result.iterations = iter;
	return result;
}
