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

// -----------------------------------------------------------------------------
// Support function: returns world-space support point + feature ID.
// Point/segment operate in world space with no rotation.
//
// Feature IDs per shape type:
//   Point:    always 0
//   Segment:  0 = endpoint p, 1 = endpoint q
//   Box:      3 sign bits (0-7), one per axis
//   Hull:     vertex index
//   Cylinder: 0 = endpoint p, 1 = endpoint q

static __forceinline v3 gjk_support(const GJK_Shape* s, v3 d, int* feat)
{
	switch (s->type) {
	case GJK_POINT: *feat = 0; return s->point.center;
	case GJK_SEGMENT: {
		int f = dot(d, sub(s->segment.q, s->segment.p)) >= 0.0f;
		*feat = f;
		return f ? s->segment.q : s->segment.p;
	}
	case GJK_BOX: {
		v3 ld = rotate(s->box.inv_rot, d);
		v3 he = s->box.half_extents;
		v3 local = V3(ld.x >= 0 ? he.x : -he.x, ld.y >= 0 ? he.y : -he.y, ld.z >= 0 ? he.z : -he.z);
		*feat = (ld.x >= 0.0f) | ((ld.y >= 0.0f) << 1) | ((ld.z >= 0.0f) << 2);
		return add(s->box.center, rotate(s->box.rot, local));
	}
	case GJK_HULL: {
		v3 ld = rotate(s->hull.inv_rot, d);
		float best = -1e18f;
		int bi = 0;
		for (int i = 0; i < s->hull.count; i++) {
			float dd = dot(s->hull.verts[i], ld);
			if (dd > best) { best = dd; bi = i; }
		}
		*feat = bi;
		return add(s->hull.center, rotate(s->hull.rot, s->hull.verts[bi]));
	}
	case GJK_CYLINDER: {
		v3 u = s->cylinder.axis;
		if (s->cylinder.inv_axis_len == 0.0f) { *feat = 0; return s->cylinder.p; }
		float da = dot(d, u);
		int f = da >= 0.0f;
		*feat = f;
		v3 base = f ? s->cylinder.q : s->cylinder.p;
		v3 dp = sub(d, scale(u, da));
		float pl = len(dp);
		if (pl > FLT_EPSILON) return add(base, scale(dp, s->cylinder.radius / pl));
		return base;
	}
	}
	*feat = 0;
	return V3(0,0,0);
}

static __forceinline v3 gjk_center(const GJK_Shape* s)
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

// -----------------------------------------------------------------------------
// Simplex solvers: find closest point on simplex to origin.

static int gjk_solve2(GJK_Simplex* s)
{
	v3 a = s->v[0].point, b = s->v[1].point;
	float u = dot(b, sub(b, a));
	float v = dot(a, sub(a, b));
	if (v <= 0.0f) { s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	if (u <= 0.0f) { s->v[0] = s->v[1]; s->v[0].u = 1.0f; s->divisor = 1.0f; s->count = 1; return 1; }
	s->divisor = u + v;
	if (s->divisor == 0.0f) return 0;
	s->v[0].u = u; s->v[1].u = v; s->count = 2;
	return 1;
}

static int gjk_solve3(GJK_Simplex* s)
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
// Closest point and witness points.

static __forceinline v3 gjk_closest_point(const GJK_Simplex* s)
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

static __forceinline void gjk_witness_points(const GJK_Simplex* s, v3* p1, v3* p2, int* f1, int* f2)
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
	v3 sA = gjk_support(&shapeA, init_d, &fA);
	v3 sB = gjk_support(&shapeB, neg(init_d), &fB);
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
		GJK_Simplex backup = simplex;
		int solved = 1;
		switch (simplex.count) {
		case 1: break;
		case 2: solved = gjk_solve2(&simplex); break;
		case 3: solved = gjk_solve3(&simplex); break;
		case 4: solved = gjk_solve4(&simplex); break;
		}
		if (!solved) { simplex = backup; break; }
		if (simplex.count == 4) break;

		v3 closest = gjk_closest_point(&simplex);
		float dsq = len2(closest);

		if (dsq <= GJK_CONTAINMENT_EPS2) break;
		if (dsq >= dsq_prev) break;
		dsq_prev = dsq;

		v3 d = neg(closest);
		sA = gjk_support(&shapeA, neg(d), &fA);
		sB = gjk_support(&shapeB, d, &fB);
		v3 w = sub(sB, sA);

		float max_vert2 = 0.0f;
		for (int i = 0; i < simplex.count; i++) {
			float v2 = len2(simplex.v[i].point);
			if (v2 > max_vert2) max_vert2 = v2;
		}
		float progress = dot(sub(w, closest), neg(closest));
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
