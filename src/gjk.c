// See LICENSE for licensing info.
// gjk.c -- GJK distance algorithm for 3D convex shapes.

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

// GJK shape: tagged union for support function dispatch.
// Support operates on UNSCALED local-space vertices. The caller transforms
// the result to world space. For hulls, scale is baked into the vertex lookup.
enum { GJK_POINT, GJK_SEGMENT, GJK_HULL };

typedef struct GjkShape
{
	int type;
	int count;
	const v3* verts;  // pointer to vertex array (local space)
	v3 verts_buf[2];  // inline storage for point/segment
} GjkShape;

static int gjk_get_support(const GjkShape* s, v3 dir)
{
	int best = 0;
	float best_d = dot(s->verts[0], dir);
	for (int i = 1; i < s->count; i++) {
		float d = dot(s->verts[i], dir);
		if (d > best_d) { best_d = d; best = i; }
	}
	return best;
}

// Positioned shape: GjkShape + transform for world-space queries.
typedef struct GjkProxy
{
	GjkShape shape;
	v3 pos;
	quat rot;
} GjkProxy;

// Proxy init functions: take output pointer to avoid dangling self-referential
// pointers when returning structs by value (verts -> verts_buf).
static void gjk_proxy_point(GjkProxy* pr, v3 p)
{
	*pr = (GjkProxy){0};
	pr->shape.type = GJK_POINT;
	pr->shape.count = 1;
	pr->shape.verts_buf[0] = p;
	pr->shape.verts = pr->shape.verts_buf;
	pr->pos = V3(0,0,0);
	pr->rot = quat_identity();
}

static void gjk_proxy_segment(GjkProxy* pr, v3 p, v3 q)
{
	*pr = (GjkProxy){0};
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
static void gjk_proxy_hull(GjkProxy* pr, const Hull* hull, v3 pos, quat rot, v3 sc, v3* scaled_verts)
{
	for (int i = 0; i < hull->vert_count; i++)
		scaled_verts[i] = hmul(hull->verts[i], sc);
	*pr = (GjkProxy){0};
	pr->shape.type = GJK_HULL;
	pr->shape.count = hull->vert_count;
	pr->shape.verts = scaled_verts;
	pr->pos = pos;
	pr->rot = rot;
}

// Legacy function pointer typedef kept for gjk_bench compatibility.
typedef v3 (*GjkSupportFn)(const void* shape, v3 dir, int* out_index);

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.
// Ported directly from lm engine (lmSimplex::Solve2/3/4).

static int gjk_solve2(GjkSimplex* s)
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

static int gjk_solve3(GjkSimplex* s)
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

static int gjk_solve4(GjkSimplex* s)
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

static v3 gjk_closest_point(const GjkSimplex* s)
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

static void gjk_witness_points(const GjkSimplex* s, v3* p1, v3* p2)
{
	float inv = 1.0f / s->divisor;
	*p1 = V3(0,0,0); *p2 = V3(0,0,0);
	for (int i = 0; i < s->count; i++) {
		*p1 = add(*p1, scale(s->v[i].point1, s->v[i].u * inv));
		*p2 = add(*p2, scale(s->v[i].point2, s->v[i].u * inv));
	}
}

static v3 gjk_search_dir(const GjkSimplex* s)
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

static GjkResult gjk_distance_ex(const GjkProxy* proxyA, const GjkProxy* proxyB)
{
	GjkResult result = {0};
	GjkSimplex simplex = {0};

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

		GjkSimplex backup = simplex;
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

		GjkVertex* v = &simplex.v[simplex.count];
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
