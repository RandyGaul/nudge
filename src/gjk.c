// See LICENSE for licensing info.
// gjk.c -- GJK distance algorithm for 3D convex shapes.
//
// 3D GJK distance algorithm with full tetrahedron Voronoi regions.
//
// Computes closest distance between two convex shapes defined by
// support functions. Returns witness points on each shape.

// -----------------------------------------------------------------------------
// GJK types.

typedef struct GjkVertex
{
	v3 point1;   // support point on shape A (world space)
	v3 point2;   // support point on shape B (world space)
	v3 point;    // Minkowski difference: point2 - point1
	float u;     // barycentric coordinate
	int index1;  // support index on A
	int index2;  // support index on B
} GjkVertex;

typedef struct GjkSimplex
{
	GjkVertex v[4];
	float divisor;
	int count;
} GjkSimplex;

typedef struct GjkResult
{
	v3 point1;     // closest point on shape A
	v3 point2;     // closest point on shape B
	float distance;
	int iterations;
} GjkResult;

// Support function signature: returns support point in world space
// and writes the vertex index to *out_index.
typedef v3 (*GjkSupportFn)(const void* shape, v3 dir, int* out_index);

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.

// Line segment: Voronoi regions A, B, AB
static void gjk_solve2(GjkSimplex* s)
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

// Triangle: Voronoi regions A, B, C, AB, BC, CA, ABC
static void gjk_solve3(GjkSimplex* s)
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

	// Vertex regions
	if (vAB <= 0.0f && uCA <= 0.0f) {
		s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}
	if (uAB <= 0.0f && vBC <= 0.0f) {
		s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}
	if (uBC <= 0.0f && vCA <= 0.0f) {
		s->v[0] = s->v[2]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return;
	}

	// Triangle normal and area-based barycentric coordinates
	v3 n = cross(sub(B, A), sub(C, A));
	float uABC = dot(cross(sub(B, Q), sub(C, Q)), n);
	float vABC = dot(cross(sub(C, Q), sub(A, Q)), n);
	float wABC = dot(cross(sub(A, Q), sub(B, Q)), n);

	// Edge regions
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

	// Interior
	s->v[0].u = uABC; s->v[1].u = vABC; s->v[2].u = wABC;
	s->divisor = dot(n, n);
	s->count = 3;
}

// Tetrahedron: Voronoi regions for 4 verts, 6 edges, 4 faces, interior.
// Uses scalar triple product for volume-based barycentric coordinates.
static inline float triple(v3 a, v3 b, v3 c) { return dot(a, cross(b, c)); }

static void gjk_solve4(GjkSimplex* s)
{
	v3 A = s->v[0].point;
	v3 B = s->v[1].point;
	v3 C = s->v[2].point;
	v3 D = s->v[3].point;
	v3 Q = V3(0,0,0);

	// Face normals (outward from each face)
	// Compute signed volumes for the 4 sub-tetrahedra
	float uABCD = triple(sub(B, Q), sub(C, Q), sub(D, Q));
	float vABCD = triple(sub(A, Q), sub(D, Q), sub(C, Q));
	float wABCD = triple(sub(A, Q), sub(B, Q), sub(D, Q));
	float xABCD = triple(sub(A, Q), sub(C, Q), sub(B, Q));

	// If origin is inside the tetrahedron
	if (uABCD > 0.0f && vABCD > 0.0f && wABCD > 0.0f && xABCD > 0.0f) {
		s->v[0].u = uABCD; s->v[1].u = vABCD;
		s->v[2].u = wABCD; s->v[3].u = xABCD;
		s->divisor = uABCD + vABCD + wABCD + xABCD;
		s->count = 4;
		return;
	}

	// Otherwise, find closest feature among the 4 triangular faces.
	// Try each face (triangle), keep the one closest to origin.
	float best_dist = 1e18f;
	GjkSimplex best = *s;

	// Face BCD (opposite A, when uABCD <= 0)
	if (uABCD <= 0.0f) {
		GjkSimplex test = *s;
		test.v[0] = s->v[1]; test.v[1] = s->v[2]; test.v[2] = s->v[3];
		test.count = 3;
		gjk_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	// Face ACD (opposite B, when vABCD <= 0)
	if (vABCD <= 0.0f) {
		GjkSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[2]; test.v[2] = s->v[3];
		test.count = 3;
		gjk_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	// Face ABD (opposite C, when wABCD <= 0)
	if (wABCD <= 0.0f) {
		GjkSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[1]; test.v[2] = s->v[3];
		test.count = 3;
		gjk_solve3(&test);
		float inv = 1.0f / test.divisor;
		v3 cp = V3(0,0,0);
		for (int i = 0; i < test.count; i++)
			cp = add(cp, scale(test.v[i].point, test.v[i].u * inv));
		float d = len2(cp);
		if (d < best_dist) { best_dist = d; best = test; }
	}

	// Face ABC (opposite D, when xABCD <= 0)
	if (xABCD <= 0.0f) {
		GjkSimplex test = *s;
		test.v[0] = s->v[0]; test.v[1] = s->v[1]; test.v[2] = s->v[2];
		test.count = 3;
		gjk_solve3(&test);
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

static v3 gjk_closest_point(const GjkSimplex* s)
{
	float inv = 1.0f / s->divisor;
	switch (s->count) {
	case 1: return s->v[0].point;
	case 2: return add(scale(s->v[0].point, s->v[0].u * inv),
	                      scale(s->v[1].point, s->v[1].u * inv));
	case 3: return add(add(
	                      scale(s->v[0].point, s->v[0].u * inv),
	                      scale(s->v[1].point, s->v[1].u * inv)),
	                      scale(s->v[2].point, s->v[2].u * inv));
	case 4: return V3(0,0,0); // origin is inside
	}
	return V3(0,0,0);
}

// Witness points on each shape from barycentric coords.
static void gjk_witness_points(const GjkSimplex* s, v3* p1, v3* p2)
{
	float inv = 1.0f / s->divisor;
	*p1 = V3(0,0,0); *p2 = V3(0,0,0);
	for (int i = 0; i < s->count; i++) {
		*p1 = add(*p1, scale(s->v[i].point1, s->v[i].u * inv));
		*p2 = add(*p2, scale(s->v[i].point2, s->v[i].u * inv));
	}
}

// Search direction: negative of closest point for 1-simplex,
// perpendicular to edge/face for higher simplices.
static v3 gjk_search_dir(const GjkSimplex* s)
{
	switch (s->count) {
	case 1: return neg(s->v[0].point);
	case 2: {
		v3 ab = sub(s->v[1].point, s->v[0].point);
		v3 ao = neg(s->v[0].point);
		// Triple cross: ab x ao x ab (perpendicular to ab toward origin)
		v3 n = cross(cross(ab, ao), ab);
		return len2(n) > 1e-12f ? n : neg(gjk_closest_point(s));
	}
	case 3: {
		v3 ab = sub(s->v[1].point, s->v[0].point);
		v3 ac = sub(s->v[2].point, s->v[0].point);
		v3 n = cross(ab, ac);
		// Orient toward origin
		if (dot(n, s->v[0].point) > 0.0f) n = neg(n);
		return n;
	}
	}
	return V3(0,0,0);
}

// -----------------------------------------------------------------------------
// Main GJK distance function.

#define GJK_MAX_ITERS 20

static GjkResult gjk_distance(
	GjkSupportFn support_a, const void* shape_a,
	GjkSupportFn support_b, const void* shape_b)
{
	GjkResult result = {0};

	// Initialize simplex with first support point
	GjkSimplex simplex = {0};
	int idx_a, idx_b;
	v3 d = V3(1, 0, 0); // initial search direction

	simplex.v[0].point1 = support_a(shape_a, neg(d), &idx_a);
	simplex.v[0].point2 = support_b(shape_b, d, &idx_b);
	simplex.v[0].point = sub(simplex.v[0].point2, simplex.v[0].point1);
	simplex.v[0].index1 = idx_a;
	simplex.v[0].index2 = idx_b;
	simplex.v[0].u = 1.0f;
	simplex.divisor = 1.0f;
	simplex.count = 1;

	int save1[4], save2[4];
	int iter = 0;

	while (iter < GJK_MAX_ITERS) {
		// Save indices for duplicate check
		int save_count = simplex.count;
		for (int i = 0; i < save_count; i++) {
			save1[i] = simplex.v[i].index1;
			save2[i] = simplex.v[i].index2;
		}

		// Solve for closest point on simplex to origin
		switch (simplex.count) {
		case 1: break;
		case 2: gjk_solve2(&simplex); break;
		case 3: gjk_solve3(&simplex); break;
		case 4: gjk_solve4(&simplex); break;
		}

		// Origin inside tetrahedron = shapes overlap
		if (simplex.count == 4) break;

		d = gjk_search_dir(&simplex);
		if (len2(d) < 1e-14f) break;

		// New support point
		GjkVertex* vertex = &simplex.v[simplex.count];
		v3 nd = norm(d);
		vertex->point1 = support_a(shape_a, neg(nd), &idx_a);
		vertex->point2 = support_b(shape_b, nd, &idx_b);
		vertex->point = sub(vertex->point2, vertex->point1);
		vertex->index1 = idx_a;
		vertex->index2 = idx_b;
		iter++;

		// Duplicate check
		int duplicate = 0;
		for (int i = 0; i < save_count; i++) {
			if (vertex->index1 == save1[i] && vertex->index2 == save2[i]) {
				duplicate = 1; break;
			}
		}
		if (duplicate) break;

		simplex.count++;
	}

	gjk_witness_points(&simplex, &result.point1, &result.point2);
	result.distance = len(sub(result.point2, result.point1));
	result.iterations = iter;
	return result;
}
