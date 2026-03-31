// See LICENSE for licensing info.
// collision.c -- broadphase + narrowphase collision detection
//
// Hull SAT with Gauss map pruning for edge-edge axis elimination.
// Half-edge mesh enables efficient face/edge traversal for SAT queries.

// Types (HalfEdge, HullPlane, HullFace, Hull, ConvexHull) defined in nudge.h.

// Support function: furthest vertex along a direction.
static v3 hull_support(const Hull* hull, v3 dir)
{
	float best = -1e18f;
	int best_i = 0;
	for (int i = 0; i < hull->vert_count; i++) {
		float d = dot(hull->verts[i], dir);
		if (d > best) { best = d; best_i = i; }
	}
	return hull->verts[best_i];
}

// -----------------------------------------------------------------------------
// Unit box hull. Half-extents = (1,1,1). Scaled at query time.
//
// Edges stored in twin pairs: edge 2k and twin 2k+1.
// 8 verts, 24 half-edges (12 edges), 6 faces.

static const v3 s_box_verts[8] = {
	{-1,-1,-1}, { 1,-1,-1}, { 1, 1,-1}, {-1, 1,-1},
	{-1,-1, 1}, { 1,-1, 1}, { 1, 1, 1}, {-1, 1, 1},
};

// Face winding (CCW from outside):
//  0: -Z (0,3,2,1)   1: +Z (4,5,6,7)
//  2: -X (0,4,7,3)   3: +X (1,2,6,5)
//  4: -Y (0,1,5,4)   5: +Y (3,7,6,2)
//
// 12 undirected edges, each stored as (edge, twin) pair at indices (2k, 2k+1):
//  pair 0: 0->3 / 3->0     pair 1: 3->2 / 2->3
//  pair 2: 2->1 / 1->2     pair 3: 1->0 / 0->1
//  pair 4: 4->5 / 5->4     pair 5: 5->6 / 6->5
//  pair 6: 6->7 / 7->6     pair 7: 7->4 / 4->7
//  pair 8: 0->4 / 4->0     pair 9: 1->5 / 5->1
//  pair10: 2->6 / 6->2     pair11: 3->7 / 7->3

static const HalfEdge s_box_edges[24] = {
	// pair 0: 0->3 (face 0) / 3->0 (face 2)
	{ .twin =  1, .next =  2, .origin = 0, .face = 0 },  // e0:  0->3
	{ .twin =  0, .next = 16, .origin = 3, .face = 2 },  // e1:  3->0
	// pair 1: 3->2 (face 0) / 2->3 (face 5)
	{ .twin =  3, .next =  4, .origin = 3, .face = 0 },  // e2:  3->2
	{ .twin =  2, .next = 22, .origin = 2, .face = 5 },  // e3:  2->3
	// pair 2: 2->1 (face 0) / 1->2 (face 3)
	{ .twin =  5, .next =  6, .origin = 2, .face = 0 },  // e4:  2->1
	{ .twin =  4, .next = 20, .origin = 1, .face = 3 },  // e5:  1->2
	// pair 3: 1->0 (face 0) / 0->1 (face 4)
	{ .twin =  7, .next =  0, .origin = 1, .face = 0 },  // e6:  1->0
	{ .twin =  6, .next = 18, .origin = 0, .face = 4 },  // e7:  0->1
	// pair 4: 4->5 (face 1) / 5->4 (face 4)
	{ .twin =  9, .next = 10, .origin = 4, .face = 1 },  // e8:  4->5
	{ .twin =  8, .next = 17, .origin = 5, .face = 4 },  // e9:  5->4
	// pair 5: 5->6 (face 1) / 6->5 (face 3)
	{ .twin = 11, .next = 12, .origin = 5, .face = 1 },  // e10: 5->6
	{ .twin = 10, .next = 19, .origin = 6, .face = 3 },  // e11: 6->5
	// pair 6: 6->7 (face 1) / 7->6 (face 5)
	{ .twin = 13, .next = 14, .origin = 6, .face = 1 },  // e12: 6->7
	{ .twin = 12, .next = 21, .origin = 7, .face = 5 },  // e13: 7->6
	// pair 7: 7->4 (face 1) / 4->7 (face 2)
	{ .twin = 15, .next =  8, .origin = 7, .face = 1 },  // e14: 7->4
	{ .twin = 14, .next = 23, .origin = 4, .face = 2 },  // e15: 4->7
	// pair 8: 0->4 (face 2) / 4->0 (face 4)
	{ .twin = 17, .next = 15, .origin = 0, .face = 2 },  // e16: 0->4
	{ .twin = 16, .next =  7, .origin = 4, .face = 4 },  // e17: 4->0
	// pair 9: 1->5 (face 4) / 5->1 (face 3)
	{ .twin = 19, .next =  9, .origin = 1, .face = 4 },  // e18: 1->5
	{ .twin = 18, .next =  5, .origin = 5, .face = 3 },  // e19: 5->1
	// pair 10: 2->6 (face 3) / 6->2 (face 5)
	{ .twin = 21, .next = 11, .origin = 2, .face = 3 },  // e20: 2->6
	{ .twin = 20, .next =  3, .origin = 6, .face = 5 },  // e21: 6->2
	// pair 11: 3->7 (face 5) / 7->3 (face 2)
	{ .twin = 23, .next = 13, .origin = 3, .face = 5 },  // e22: 3->7
	{ .twin = 22, .next =  1, .origin = 7, .face = 2 },  // e23: 7->3
};

static const HullFace s_box_faces[6] = {
	{ .edge =  0 },  // face 0 (-Z): starts at e0 (0->3)
	{ .edge =  8 },  // face 1 (+Z): starts at e8 (4->5)
	{ .edge = 16 },  // face 2 (-X): starts at e16 (0->4)
	{ .edge =  5 },  // face 3 (+X): starts at e5 (1->2)
	{ .edge =  7 },  // face 4 (-Y): starts at e7 (0->1)
	{ .edge = 22 },  // face 5 (+Y): starts at e22 (3->7)
};

static const HullPlane s_box_planes[6] = {
	{ .normal = { 0, 0,-1}, .offset = 1 },
	{ .normal = { 0, 0, 1}, .offset = 1 },
	{ .normal = {-1, 0, 0}, .offset = 1 },
	{ .normal = { 1, 0, 0}, .offset = 1 },
	{ .normal = { 0,-1, 0}, .offset = 1 },
	{ .normal = { 0, 1, 0}, .offset = 1 },
};

static const Hull s_unit_box_hull = {
	.centroid = {0, 0, 0},
	.verts = s_box_verts,
	.edges = s_box_edges,
	.faces = s_box_faces,
	.planes = s_box_planes,
	.vert_count = 8,
	.edge_count = 24,
	.face_count = 6,
	.epsilon = 9.0f * FLT_EPSILON,
};

// -----------------------------------------------------------------------------
// Transform helpers for hull queries.

// Scale a hull vertex by half-extents.
static inline v3 hull_vert_scaled(const Hull* hull, int i, v3 scale)
{
	v3 v = hull->verts[i];
	return (v3){ v.x * scale.x, v.y * scale.y, v.z * scale.z };
}

// Transform a plane from hull-local into world space.
static inline HullPlane plane_transform(HullPlane p, v3 pos, quat rot, v3 scale)
{
	// For axis-aligned unit box with uniform-ish scale:
	// normal needs inverse-transpose scale, then rotate.
	v3 n = norm(rotate(rot, V3(
		p.normal.x / scale.x, p.normal.y / scale.y, p.normal.z / scale.z)));
	// A point on the plane in local space: normal * offset, then scale + transform
	v3 local_pt = V3(p.normal.x * p.offset * scale.x,
	                      p.normal.y * p.offset * scale.y,
	                      p.normal.z * p.offset * scale.z);
	v3 world_pt = add(pos, rotate(rot, local_pt));
	return (HullPlane){ .normal = n, .offset = dot(n, world_pt) };
}

// Project hull onto a plane, return signed distance of support point.
static float hull_project_plane(const Hull* hull, HullPlane plane, v3 scale)
{
	v3 support = hull_support(hull, neg(plane.normal));
	v3 sv = { support.x * scale.x, support.y * scale.y, support.z * scale.z };
	return dot(plane.normal, sv) - plane.offset;
}

// -----------------------------------------------------------------------------
// Internal manifold used by broadphase (extends public Manifold with body indices).

typedef struct InternalManifold
{
	Manifold m;
	int body_a;
	int body_b;
} InternalManifold;

// -----------------------------------------------------------------------------
// Sphere-sphere.

int collide_sphere_sphere(Sphere a, Sphere b, Manifold* manifold)
{
	v3 d = sub(b.center, a.center);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
	float penetration = r_sum - dist;

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(a.center, scale(normal, a.radius - penetration * 0.5f)),
		.normal = normal,
		.penetration = penetration,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Segment helpers for capsule collisions.

// Closest point on segment PQ to point X. Returns parametric t in [0,1].
static float segment_closest_t(v3 P, v3 Q, v3 X)
{
	v3 d = sub(Q, P);
	float d_len2 = len2(d);
	if (d_len2 < 1e-12f) return 0.0f;
	float t = dot(sub(X, P), d) / d_len2;
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;
	return t;
}

static v3 segment_closest_point(v3 P, v3 Q, v3 X)
{
	float t = segment_closest_t(P, Q, X);
	return add(P, scale(sub(Q, P), t));
}

// Closest points between two segments P1Q1 and P2Q2.
static void segments_closest_points(v3 P1, v3 Q1, v3 P2, v3 Q2, v3* out1, v3* out2)
{
	v3 d1 = sub(Q1, P1);
	v3 d2 = sub(Q2, P2);
	v3 r = sub(P1, P2);
	float a = dot(d1, d1);
	float e = dot(d2, d2);
	float f = dot(d2, r);
	float s, t;

	if (a < 1e-12f && e < 1e-12f) {
		*out1 = P1; *out2 = P2; return;
	}
	if (a < 1e-12f) {
		s = 0.0f;
		t = f / e;
		if (t < 0.0f) t = 0.0f; if (t > 1.0f) t = 1.0f;
	} else {
		float c = dot(d1, r);
		if (e < 1e-12f) {
			t = 0.0f;
			s = -c / a;
			if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f;
		} else {
			float b = dot(d1, d2);
			float denom = a * e - b * b;
			s = denom > 1e-12f ? (b * f - c * e) / denom : 0.0f;
			if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f;
			t = (b * s + f) / e;
			if (t < 0.0f) { t = 0.0f; s = -c / a; if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f; }
			else if (t > 1.0f) { t = 1.0f; s = (b - c) / a; if (s < 0.0f) s = 0.0f; if (s > 1.0f) s = 1.0f; }
		}
	}
	*out1 = add(P1, scale(d1, s));
	*out2 = add(P2, scale(d2, t));
}

// Get capsule segment endpoints in world space.
static void capsule_world_segment(BodyHot* h, ShapeInternal* s, v3* P, v3* Q)
{
	v3 local_p = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 local_q = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	*P = add(h->position, rotate(h->rotation, local_p));
	*Q = add(h->position, rotate(h->rotation, local_q));
}

// -----------------------------------------------------------------------------
// Sphere-capsule.

int collide_sphere_capsule(Sphere a, Capsule b, Manifold* manifold)
{
	v3 closest = segment_closest_point(b.p, b.q, a.center);
	v3 d = sub(closest, a.center);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(a.center, scale(normal, a.radius)),
		.normal = normal,
		.penetration = r_sum - dist,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// Capsule-capsule.

int collide_capsule_capsule(Capsule a, Capsule b, Manifold* manifold)
{
	v3 c1, c2;
	segments_closest_points(a.p, a.q, b.p, b.q, &c1, &c2);

	v3 d = sub(c2, c1);
	float dist2 = len2(d);
	float r_sum = a.radius + b.radius;

	if (dist2 > r_sum * r_sum) return 0;
	if (!manifold) return 1;

	float dist = sqrtf(dist2);
	v3 normal = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(c1, scale(normal, a.radius)),
		.normal = normal,
		.penetration = r_sum - dist,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// Capsule-box (hull) narrowphase.
// Shallow: GJK witness points. Deep: face search (lm-style).

// GJK support for a capsule (core segment only, radius handled externally).
static v3 capsule_gjk_support(v3 P, v3 Q, v3 dir)
{
	return dot(P, dir) >= dot(Q, dir) ? P : Q;
}

// GJK support for a scaled hull.
static v3 hull_gjk_support(const Hull* hull, v3 pos, quat rot, v3 scale, v3 dir)
{
	v3 local_dir = rotate(inv(rot), dir);
	v3 local_sup = hull_support(hull, local_dir);
	v3 scaled = { local_sup.x * scale.x, local_sup.y * scale.y, local_sup.z * scale.z };
	return add(pos, rotate(rot, scaled));
}

int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold)
{
	v3 cap_p = a.p, cap_q = a.q;
	float cap_radius = a.radius;
	const Hull* hull = b.hull;
	v3 hull_pos = b.center;
	quat hull_rot = b.rotation;
	v3 hull_scale = b.scale;
	float hull_radius = 0.0f;
	// Find closest points between capsule core segment and hull
	// using support-function iteration (simplified GJK-like approach).
	// Project capsule onto hull faces to find separation.
	quat inv_rot = inv(hull_rot);

	// Transform capsule into hull local space
	v3 lp = rotate(inv_rot, sub(cap_p, hull_pos));
	v3 lq = rotate(inv_rot, sub(cap_q, hull_pos));
	// Scale to unit hull space
	lp = V3(lp.x / hull_scale.x, lp.y / hull_scale.y, lp.z / hull_scale.z);
	lq = V3(lq.x / hull_scale.x, lq.y / hull_scale.y, lq.z / hull_scale.z);

	float r_sum = cap_radius + hull_radius;

	// Find best separating face (SAT-style, like lm's sphere-to-hull deep path)
	float best_sep = -1e18f;
	int best_face = -1;

	for (int i = 0; i < hull->face_count; i++) {
		v3 n = hull->planes[i].normal;
		float d = hull->planes[i].offset;
		// Support of capsule segment along -n (closest point to face)
		float dp = dot(lp, n);
		float dq = dot(lq, n);
		float sup = dp < dq ? dp : dq; // min projection = closest to face
		float sep = sup - d; // negative = penetrating
		// Account for radii in scaled space (approximate)
		if (sep > best_sep) {
			best_sep = sep;
			best_face = i;
		}
	}

	if (best_sep > r_sum) return 0; // separated

	// For shallow hits, use closest point on capsule to hull face
	v3 face_n = hull->planes[best_face].normal;
	v3 world_n = rotate(hull_rot, face_n);

	// Clip capsule endpoints to reference face side planes
	// (simplified: just use the endpoint(s) that penetrate the face)
	float dp = dot(lp, face_n) - hull->planes[best_face].offset;
	float dq = dot(lq, face_n) - hull->planes[best_face].offset;

	int cp = 0;
	v3 points[2];
	float depths[2];

	if (dp - r_sum < 0.0f) {
		v3 wp = add(hull_pos, rotate(hull_rot, V3(
			lp.x * hull_scale.x, lp.y * hull_scale.y, lp.z * hull_scale.z)));
		points[cp] = add(wp, scale(world_n, cap_radius));
		depths[cp] = r_sum - dp;
		cp++;
	}
	if (dq - r_sum < 0.0f) {
		v3 wq = add(hull_pos, rotate(hull_rot, V3(
			lq.x * hull_scale.x, lq.y * hull_scale.y, lq.z * hull_scale.z)));
		points[cp] = add(wq, scale(world_n, cap_radius));
		depths[cp] = r_sum - dq;
		cp++;
	}

	if (cp == 0) return 0;
	if (!manifold) return 1;

	manifold->count = cp;
	for (int i = 0; i < cp; i++) {
		manifold->contacts[i] = (Contact){
			.point = points[i],
			.normal = neg(world_n), // from A (capsule) toward B (hull)
			.penetration = depths[i],
		};
	}
	return 1;
}

// -----------------------------------------------------------------------------
int collide_capsule_box(Capsule a, Box b, Manifold* manifold)
{
	return collide_capsule_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// Sphere-hull.

int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold)
{
	const Hull* hull = b.hull;
	quat inv_rot = inv(b.rotation);
	v3 local_pos = rotate(inv_rot, sub(a.center, b.center));
	v3 lp = V3(local_pos.x / b.scale.x, local_pos.y / b.scale.y, local_pos.z / b.scale.z);

	float best_sep = -1e18f;
	int best_face = -1;

	for (int i = 0; i < hull->face_count; i++) {
		float sep = dot(lp, hull->planes[i].normal) - hull->planes[i].offset;
		if (sep > a.radius) return 0;
		if (sep > best_sep) { best_sep = sep; best_face = i; }
	}

	if (!manifold) return 1;

	v3 world_n = rotate(b.rotation, hull->planes[best_face].normal);

	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = add(a.center, scale(world_n, -a.radius)),
		.normal = neg(world_n),
		.penetration = a.radius - best_sep,
	};
	return 1;
}

int collide_sphere_box(Sphere a, Box b, Manifold* manifold)
{
	return collide_sphere_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// -----------------------------------------------------------------------------
// SAT: face queries.
// All computations in local space of Hull2.

typedef struct FaceQuery
{
	int index;
	float separation;
} FaceQuery;

static FaceQuery sat_query_faces(
	const Hull* hull1, v3 pos1, quat rot1, v3 scale1,
	const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
{
	// Transform from world to hull2 local: p_local = rot2^-1 * (p_world - pos2)
	quat inv2 = inv(rot2);

	FaceQuery best = { .index = -1, .separation = -1e18f };

	for (int i = 0; i < hull1->face_count; i++) {
		// Transform hull1's plane into hull2's local space
		HullPlane pw = plane_transform(hull1->planes[i], pos1, rot1, scale1);
		// Bring plane normal into hull2 local
		v3 local_n = rotate(inv2, pw.normal);
		v3 local_pt = rotate(inv2, sub(scale(pw.normal, pw.offset), pos2));
		// Actually, let's just use world space and project hull2's support
		// Support of hull2 in direction -pw.normal (world space)
		v3 sup_dir_local = rotate(inv2, neg(pw.normal));
		v3 sup_local = hull_support(hull2, sup_dir_local);
		v3 sup_scaled = { sup_local.x * scale2.x, sup_local.y * scale2.y, sup_local.z * scale2.z };
		v3 sup_world = add(pos2, rotate(rot2, sup_scaled));

		float sep = dot(pw.normal, sup_world) - pw.offset;
		if (sep > best.separation) {
			best.separation = sep;
			best.index = i;
		}
	}
	return best;
}

// -----------------------------------------------------------------------------
// SAT: Gauss map Minkowski face test.
// Tests if arcs AB and CD intersect on the unit sphere.

static int is_minkowski_face(v3 a, v3 b, v3 b_x_a, v3 c, v3 d, v3 d_x_c)
{
	float cba = dot(c, b_x_a);
	float dba = dot(d, b_x_a);
	float adc = dot(a, d_x_c);
	float bdc = dot(b, d_x_c);
	return (cba * dba < 0.0f) && (adc * bdc < 0.0f) && (cba * bdc > 0.0f);
}

// Project edge pair: signed distance along cross(e1,e2) from p1 to p2.
static float sat_edge_project(v3 p1, v3 e1, v3 p2, v3 e2, v3 c1)
{
	v3 e1_x_e2 = cross(e1, e2);
	float len = len(e1_x_e2);

	// Skip near-parallel edges
	float tolerance = 0.005f;
	if (len < tolerance * sqrtf(len2(e1) * len2(e2)))
		return -1e18f;

	v3 n = scale(e1_x_e2, 1.0f / len);
	// Ensure consistent orientation (hull1 -> hull2)
	if (dot(n, sub(p1, c1)) < 0.0f)
		n = neg(n);

	return dot(n, sub(p2, p1));
}

// SAT: edge queries with Gauss map pruning.
typedef struct EdgeQuery
{
	int index1;
	int index2;
	float separation;
} EdgeQuery;

static EdgeQuery sat_query_edges(
	const Hull* hull1, v3 pos1, quat rot1, v3 scale1,
	const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
{
	// All in local space of hull2.
	quat inv2 = inv(rot2);
	quat rel_rot = mul(inv2, rot1); // rotation from hull1-local to hull2-local
	v3 c1_local = rotate(inv2, sub(pos1, pos2)); // hull1 centroid in hull2 space

	EdgeQuery best = { .index1 = -1, .index2 = -1, .separation = -1e18f };

	for (int i1 = 0; i1 < hull1->edge_count; i1 += 2) {
		const HalfEdge* edge1 = &hull1->edges[i1];
		const HalfEdge* twin1 = &hull1->edges[i1 + 1];

		v3 p1 = hull_vert_scaled(hull1, edge1->origin, scale1);
		v3 q1 = hull_vert_scaled(hull1, twin1->origin, scale1);
		// Transform to hull2 local
		p1 = add(c1_local, rotate(rel_rot, p1));
		q1 = add(c1_local, rotate(rel_rot, q1));
		v3 e1 = sub(q1, p1);

		v3 u1 = rotate(rel_rot, hull1->planes[edge1->face].normal);
		v3 v1 = rotate(rel_rot, hull1->planes[twin1->face].normal);

		for (int i2 = 0; i2 < hull2->edge_count; i2 += 2) {
			const HalfEdge* edge2 = &hull2->edges[i2];
			const HalfEdge* twin2 = &hull2->edges[i2 + 1];

			v3 p2 = hull_vert_scaled(hull2, edge2->origin, scale2);
			v3 q2 = hull_vert_scaled(hull2, twin2->origin, scale2);
			v3 e2 = sub(q2, p2);

			v3 u2 = hull2->planes[edge2->face].normal;
			v3 v2 = hull2->planes[twin2->face].normal;

			// Gauss map pruning
			if (!is_minkowski_face(u1, v1, neg(e1), neg(u2), neg(v2), neg(e2)))
				continue;

			float sep = sat_edge_project(p1, e1, p2, e2, c1_local);
			if (sep > best.separation) {
				best.index1 = i1;
				best.index2 = i2;
				best.separation = sep;
			}
		}
	}
	return best;
}

// -----------------------------------------------------------------------------
// Contact manifold helpers for hull-hull clipping.

#define MAX_CLIP_VERTS 64

// Collect face vertices in world space by walking the half-edge loop.
static int hull_face_verts_world(
	const Hull* hull, int face_idx, v3 pos, quat rot, v3 sc, v3* out)
{
	int start = hull->faces[face_idx].edge;
	int e = start;
	int count = 0;
	do {
		out[count++] = add(pos, rotate(rot, hull_vert_scaled(hull, hull->edges[e].origin, sc)));
		e = hull->edges[e].next;
	} while (e != start && count < MAX_CLIP_VERTS);
	return count;
}

// Find face on hull most anti-parallel to a world-space normal.
static int find_incident_face(
	const Hull* hull, v3 pos, quat rot, v3 sc, v3 ref_normal)
{
	int best = 0;
	float best_dot = 1e18f;
	for (int i = 0; i < hull->face_count; i++) {
		HullPlane pw = plane_transform(hull->planes[i], pos, rot, sc);
		float d = dot(pw.normal, ref_normal);
		if (d < best_dot) { best_dot = d; best = i; }
	}
	return best;
}

// Sutherland-Hodgman: clip polygon against a single plane.
// Keeps points on the negative side: dot(plane_n, p) - plane_d <= 0.
static int clip_to_plane(
	v3* in, int in_count, v3 plane_n, float plane_d, v3* out)
{
	if (in_count == 0) return 0;
	int out_count = 0;
	v3 prev = in[in_count - 1];
	float d_prev = dot(plane_n, prev) - plane_d;

	for (int i = 0; i < in_count; i++) {
		v3 cur = in[i];
		float d_cur = dot(plane_n, cur) - plane_d;

		if (d_prev <= 0.0f) {
			out[out_count++] = prev;
			if (d_cur > 0.0f) {
				float t = d_prev / (d_prev - d_cur);
				out[out_count++] = add(prev, scale(sub(cur, prev), t));
			}
		} else if (d_cur <= 0.0f) {
			float t = d_prev / (d_prev - d_cur);
			out[out_count++] = add(prev, scale(sub(cur, prev), t));
		}

		prev = cur;
		d_prev = d_cur;
	}
	return out_count;
}

// Reduce contact manifold to MAX_CONTACTS points.
// 1) Two farthest points.
// 2) Farthest from that segment.
// 3) Maximal barycentric contributor (most outside triangle).
// 4) Deepest remaining point.
static int reduce_contacts(Contact* contacts, int count)
{
	if (count <= MAX_CONTACTS) return count;

	int sel[MAX_CONTACTS];
	int used[MAX_CLIP_VERTS];
	memset(used, 0, count * sizeof(int));

	// Step 1: two farthest points
	float best_d2 = -1.0f;
	int i0 = 0, i1 = 1;
	for (int i = 0; i < count; i++)
		for (int j = i + 1; j < count; j++) {
			float d2 = len2(sub(contacts[j].point, contacts[i].point));
			if (d2 > best_d2) { best_d2 = d2; i0 = i; i1 = j; }
		}
	sel[0] = i0; sel[1] = i1;
	used[i0] = used[i1] = 1;

	// Step 2: farthest from line through (i0, i1)
	v3 seg = sub(contacts[i1].point, contacts[i0].point);
	float seg_l2 = len2(seg);
	float best = -1.0f;
	int i2 = -1;
	for (int i = 0; i < count; i++) {
		if (used[i]) continue;
		v3 v = sub(contacts[i].point, contacts[i0].point);
		float proj = seg_l2 > 1e-12f ? dot(v, seg) / seg_l2 : 0.0f;
		float d2 = len2(sub(v, scale(seg, proj)));
		if (d2 > best) { best = d2; i2 = i; }
	}
	sel[2] = i2;
	used[i2] = 1;

	// Step 3: maximal barycentric contributor (most outside triangle)
	v3 A = contacts[sel[0]].point;
	v3 B = contacts[sel[1]].point;
	v3 C = contacts[sel[2]].point;
	v3 N = cross(sub(B, A), sub(C, A));
	float best_min = 1e18f;
	int i3 = -1;
	for (int i = 0; i < count; i++) {
		if (used[i]) continue;
		v3 P = contacts[i].point;
		float u = dot(cross(sub(B, P), sub(C, P)), N);
		float vc = dot(cross(sub(C, P), sub(A, P)), N);
		float w = dot(cross(sub(A, P), sub(B, P)), N);
		float m = u < vc ? u : vc;
		m = m < w ? m : w;
		if (m < best_min) { best_min = m; i3 = i; }
	}
	sel[3] = i3;
	used[i3] = 1;

	// Step 4: deepest remaining point
	float best_depth = -1e18f;
	int i4 = -1;
	for (int i = 0; i < count; i++) {
		if (used[i]) continue;
		if (contacts[i].penetration > best_depth) {
			best_depth = contacts[i].penetration;
			i4 = i;
		}
	}
	int result = 4;
	if (i4 >= 0) { sel[4] = i4; result = 5; }

	Contact tmp[MAX_CONTACTS];
	for (int i = 0; i < result; i++) tmp[i] = contacts[sel[i]];
	for (int i = 0; i < result; i++) contacts[i] = tmp[i];
	return result;
}

// -----------------------------------------------------------------------------
// Full SAT hull vs hull with Sutherland-Hodgman face clipping.

int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold)
{
	const Hull* hull_a = a.hull;
	v3 pos_a = a.center; quat rot_a = a.rotation; v3 scale_a = a.scale;
	const Hull* hull_b = b.hull;
	v3 pos_b = b.center; quat rot_b = b.rotation; v3 scale_b = b.scale;

	FaceQuery face_a = sat_query_faces(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b);
	if (face_a.separation > 0.0f) return 0;

	FaceQuery face_b = sat_query_faces(hull_b, pos_b, rot_b, scale_b, hull_a, pos_a, rot_a, scale_a);
	if (face_b.separation > 0.0f) return 0;

	EdgeQuery edge_q = sat_query_edges(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b);
	if (edge_q.separation > 0.0f) return 0;

	if (!manifold) return 1;

	// Bias toward face contacts over edge contacts
	const float k_tol = 0.05f;
	float max_face_sep = face_a.separation > face_b.separation
		? face_a.separation : face_b.separation;

	if (edge_q.separation > max_face_sep + k_tol) {
		// --- Edge-edge contact ---
		const HalfEdge* e1 = &hull_a->edges[edge_q.index1];
		const HalfEdge* t1 = &hull_a->edges[edge_q.index1 + 1];
		v3 p1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, e1->origin, scale_a)));
		v3 q1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, t1->origin, scale_a)));

		const HalfEdge* e2 = &hull_b->edges[edge_q.index2];
		const HalfEdge* t2 = &hull_b->edges[edge_q.index2 + 1];
		v3 p2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, e2->origin, scale_b)));
		v3 q2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, t2->origin, scale_b)));

		v3 ca, cb;
		segments_closest_points(p1, q1, p2, q2, &ca, &cb);

		v3 normal = norm(cross(sub(q1, p1), sub(q2, p2)));
		if (dot(normal, sub(pos_b, pos_a)) < 0.0f) normal = neg(normal);

		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = scale(add(ca, cb), 0.5f),
			.normal = normal,
			.penetration = -edge_q.separation,
		};
		return 1;
	}

	// --- Face contact with Sutherland-Hodgman clipping ---
	const Hull* ref_hull; const Hull* inc_hull;
	v3 ref_pos, inc_pos, ref_sc, inc_sc;
	quat ref_rot, inc_rot;
	int ref_face, flip;

	if (face_a.separation > face_b.separation + k_tol) {
		ref_hull = hull_a; ref_pos = pos_a; ref_rot = rot_a; ref_sc = scale_a;
		inc_hull = hull_b; inc_pos = pos_b; inc_rot = rot_b; inc_sc = scale_b;
		ref_face = face_a.index; flip = 0;
	} else {
		ref_hull = hull_b; ref_pos = pos_b; ref_rot = rot_b; ref_sc = scale_b;
		inc_hull = hull_a; inc_pos = pos_a; inc_rot = rot_a; inc_sc = scale_a;
		ref_face = face_b.index; flip = 1;
	}

	HullPlane ref_plane = plane_transform(ref_hull->planes[ref_face], ref_pos, ref_rot, ref_sc);

	// Collect incident face vertices
	int inc_face = find_incident_face(inc_hull, inc_pos, inc_rot, inc_sc, ref_plane.normal);
	v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
	int clip_count = hull_face_verts_world(inc_hull, inc_face, inc_pos, inc_rot, inc_sc, buf1);

	// Clip against each side plane of the reference face
	v3* in_buf = buf1;
	v3* out_buf = buf2;
	int start_e = ref_hull->faces[ref_face].edge;
	int ei = start_e;
	do {
		const HalfEdge* edge = &ref_hull->edges[ei];
		v3 tail = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, edge->origin, ref_sc)));
		v3 head = add(ref_pos, rotate(ref_rot,
			hull_vert_scaled(ref_hull, ref_hull->edges[edge->twin].origin, ref_sc)));
		v3 side_n = norm(cross(sub(head, tail), ref_plane.normal));
		float side_d = dot(side_n, tail);

		clip_count = clip_to_plane(in_buf, clip_count, side_n, side_d, out_buf);
		v3* swap = in_buf; in_buf = out_buf; out_buf = swap;

		ei = edge->next;
	} while (ei != start_e);

	// Keep points below reference face, generate contacts
	v3 contact_n = flip ? neg(ref_plane.normal) : ref_plane.normal;
	Contact tmp_contacts[MAX_CLIP_VERTS];
	int cp = 0;
	for (int i = 0; i < clip_count; i++) {
		float depth = ref_plane.offset - dot(ref_plane.normal, in_buf[i]);
		if (depth >= 0.0f) {
			tmp_contacts[cp++] = (Contact){
				.point = in_buf[i],
				.normal = contact_n,
				.penetration = depth,
			};
		}
	}

	if (cp == 0) return 0;

	cp = reduce_contacts(tmp_contacts, cp);
	manifold->count = cp;
	for (int i = 0; i < cp; i++)
		manifold->contacts[i] = tmp_contacts[i];
	return 1;
}

// -----------------------------------------------------------------------------
int collide_box_box(Box a, Box b, Manifold* manifold)
{
	return collide_hull_hull(
		(ConvexHull){ &s_unit_box_hull, a.center, a.rotation, a.half_extents },
		(ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents },
		manifold);
}

const Hull* hull_unit_box() { return &s_unit_box_hull; }

// Quickhull implemented in quickhull.c.

void hull_free(Hull* hull)
{
	if (!hull) return;
	CK_FREE((void*)hull->verts);
	CK_FREE((void*)hull->edges);
	CK_FREE((void*)hull->faces);
	CK_FREE((void*)hull->planes);
	CK_FREE(hull);
}

// -----------------------------------------------------------------------------
// N^2 broadphase + narrowphase dispatch.

// Build a Sphere/Capsule/Box from internal body+shape for broadphase dispatch.
static Sphere make_sphere(BodyHot* h, ShapeInternal* s)
{
	return (Sphere){ add(h->position, rotate(h->rotation, s->local_pos)), s->sphere.radius };
}

static Capsule make_capsule(BodyHot* h, ShapeInternal* s)
{
	v3 lp = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 lq = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	return (Capsule){ add(h->position, rotate(h->rotation, lp)),
	                  add(h->position, rotate(h->rotation, lq)), s->capsule.radius };
}

static Box make_box(BodyHot* h, ShapeInternal* s)
{
	return (Box){ h->position, h->rotation, s->box.half_extents };
}

static ConvexHull make_convex_hull(BodyHot* h, ShapeInternal* s)
{
	return (ConvexHull){ s->hull.hull, h->position, h->rotation, s->hull.scale };
}

static void broadphase_and_collide(WorldInternal* w, InternalManifold** manifolds)
{
	int count = asize(w->body_hot);

	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (asize(w->body_cold[i].shapes) == 0) continue;

		for (int j = i + 1; j < count; j++) {
			if (!split_alive(w->body_gen, j)) continue;
			if (asize(w->body_cold[j].shapes) == 0) continue;

			if (w->body_hot[i].inv_mass == 0.0f && w->body_hot[j].inv_mass == 0.0f) continue;

			ShapeInternal* sa = &w->body_cold[i].shapes[0];
			ShapeInternal* sb = &w->body_cold[j].shapes[0];
			BodyHot* ha = &w->body_hot[i];
			BodyHot* hb = &w->body_hot[j];

			InternalManifold im = { .body_a = i, .body_b = j };

			// Ensure s0->type <= s1->type for upper-triangle dispatch.
			ShapeInternal* s0 = sa; ShapeInternal* s1 = sb;
			BodyHot* h0 = ha; BodyHot* h1 = hb;
			if (s0->type > s1->type) {
				ShapeInternal* tmp_s = s0; s0 = s1; s1 = tmp_s;
				BodyHot* tmp_h = h0; h0 = h1; h1 = tmp_h;
				im.body_a = j; im.body_b = i;
			}

			int hit = 0;

			if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_SPHERE)
				hit = collide_sphere_sphere(make_sphere(h0, s0), make_sphere(h1, s1), &im.m);
			else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CAPSULE)
				hit = collide_sphere_capsule(make_sphere(h0, s0), make_capsule(h1, s1), &im.m);
			else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_BOX)
				hit = collide_sphere_box(make_sphere(h0, s0), make_box(h1, s1), &im.m);
			else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CAPSULE)
				hit = collide_capsule_capsule(make_capsule(h0, s0), make_capsule(h1, s1), &im.m);
			else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_BOX)
				hit = collide_capsule_box(make_capsule(h0, s0), make_box(h1, s1), &im.m);
			else if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX)
				hit = collide_box_box(make_box(h0, s0), make_box(h1, s1), &im.m);
			else if (s0->type == SHAPE_BOX && s1->type == SHAPE_HULL)
				hit = collide_hull_hull(
					(ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents },
					make_convex_hull(h1, s1), &im.m);
			else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL)
				hit = collide_sphere_hull(make_sphere(h0, s0), make_convex_hull(h1, s1), &im.m);
			else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL)
				hit = collide_capsule_hull(make_capsule(h0, s0), make_convex_hull(h1, s1), &im.m);
			else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL)
				hit = collide_hull_hull(make_convex_hull(h0, s0), make_convex_hull(h1, s1), &im.m);

			if (hit) apush(*manifolds, im);
		}
	}
}
