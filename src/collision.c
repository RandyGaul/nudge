// See LICENSE for licensing info.
// collision.c -- broadphase + narrowphase collision detection
//
// Hull SAT with Gauss map pruning for edge-edge axis elimination.
// Half-edge mesh enables efficient face/edge traversal for SAT queries.

int g_hull_trace;


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
	v3 n = norm(rotate(rot, V3(p.normal.x / scale.x, p.normal.y / scale.y, p.normal.z / scale.z)));
	// A point on the plane in local space: normal * offset, then scale + transform
	v3 local_pt = V3(p.normal.x * p.offset * scale.x, p.normal.y * p.offset * scale.y, p.normal.z * p.offset * scale.z);
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
// GJK distance is fuzzy near zero -- use layered thresholds.
// LINEAR_SLOP defined in nudge.c (solver constants)

// GJK query helpers.
#define MAX_HULL_VERTS 256

static GJK_Result gjk_query_point_hull(v3 pt, ConvexHull h)
{
	v3 scaled[MAX_HULL_VERTS];
	GJK_Shape a = gjk_sphere(pt, 0);
	GJK_Shape b = gjk_hull_scaled(h.hull, h.center, h.rotation, h.scale, scaled);
	return gjk_distance(a, b);
}

static GJK_Result gjk_query_segment_hull(v3 p, v3 q, ConvexHull h)
{
	v3 scaled[MAX_HULL_VERTS];
	GJK_Shape a = gjk_capsule(p, q, 0);
	GJK_Shape b = gjk_hull_scaled(h.hull, h.center, h.rotation, h.scale, scaled);
	return gjk_distance(a, b);
}

static GJK_Result gjk_query_hull_hull(ConvexHull a, ConvexHull b)
{
	v3 sa[MAX_HULL_VERTS], sb[MAX_HULL_VERTS];
	GJK_Shape ga = gjk_hull_scaled(a.hull, a.center, a.rotation, a.scale, sa);
	GJK_Shape gb = gjk_hull_scaled(b.hull, b.center, b.rotation, b.scale, sb);
	return gjk_distance(ga, gb);
}


// Deep penetration: find most-separated face on hull from a point.
static int hull_deepest_face(const Hull* hull, v3 local_pt)
{
	float best = -1e18f;
	int face = 0;
	for (int i = 0; i < hull->face_count; i++) {
		float s = dot(local_pt, hull->planes[i].normal) - hull->planes[i].offset;
		if (s > best) { best = s; face = i; }
	}
	return face;
}

// Transform a world point to hull local (unit) space.
static v3 to_hull_local(v3 pt, v3 pos, quat rot, v3 sc)
{
	v3 lp = rotate(inv(rot), sub(pt, pos));
	return V3(lp.x / sc.x, lp.y / sc.y, lp.z / sc.z);
}

// Sphere-hull: GJK on sphere center (point) vs hull for distance,
// then shallow path (GJK witness) or deep path (face search).
int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold)
{
	// GJK on the sphere CENTER (point) vs hull core.
	GJK_Result r = gjk_query_point_hull(a.center, b);

	if (r.distance > a.radius) return 0; // separated

	if (r.distance > LINEAR_SLOP) {
		// Shallow: center is outside hull but within radius.
		v3 normal = scale(sub(r.point2, r.point1), 1.0f / r.distance);
		if (!manifold) return 1;
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = add(a.center, scale(normal, a.radius)),
			.normal = normal,
			.penetration = a.radius - r.distance,
		};
		return 1;
	}

	// Deep: center is inside hull. Find face with least penetration in world space.
	// Transform each hull plane to world space (accounts for non-uniform scale).
	float best_sep = -1e18f;
	int best_face = -1;
	HullPlane best_plane = {0};
	for (int i = 0; i < b.hull->face_count; i++) {
		HullPlane wp = plane_transform(b.hull->planes[i], b.center, b.rotation, b.scale);
		float s = dot(a.center, wp.normal) - wp.offset;
		if (s > a.radius) return 0;
		if (s > best_sep) { best_sep = s; best_face = i; best_plane = wp; }
	}
	if (!manifold) return 1;

	v3 world_pt = sub(a.center, scale(best_plane.normal, best_sep));
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = world_pt,
		.normal = neg(best_plane.normal),
		.penetration = a.radius - best_sep,
	};
	return 1;
}

int collide_sphere_box(Sphere a, Box b, Manifold* manifold)
{
	return collide_sphere_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// Capsule-hull: GJK shallow path, face/edge-search deep path.
int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold)
{
	GJK_Result r = gjk_query_segment_hull(a.p, a.q, b);

	if (r.distance > a.radius) return 0;

	if (r.distance > LINEAR_SLOP) {
		// Shallow: GJK witness points.
		v3 normal = scale(sub(r.point2, r.point1), 1.0f / r.distance);
		v3 contact_pt = add(r.point1, scale(normal, a.radius));
		if (!manifold) return 1;
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = contact_pt,
			.normal = normal,
			.penetration = a.radius - r.distance,
		};
		return 1;
	}

	// Deep: SAT on hull faces + capsule edge axes. Work in world space to
	// avoid non-uniform scale distortion (lm pattern).
	const Hull* hull = b.hull;
	v3 cap_dir = sub(a.q, a.p);
	float cap_len2 = len2(cap_dir);

	// --- Axis family 1: hull face normals (world space) ---
	float face_sep = -1e18f;
	int face_idx = -1;
	HullPlane face_plane = {0};
	for (int i = 0; i < hull->face_count; i++) {
		HullPlane wp = plane_transform(hull->planes[i], b.center, b.rotation, b.scale);
		float dp = dot(a.p, wp.normal) - wp.offset;
		float dq = dot(a.q, wp.normal) - wp.offset;
		float sup = dp < dq ? dp : dq;
		if (sup > face_sep) { face_sep = sup; face_idx = i; face_plane = wp; }
	}

	// --- Axis family 2: capsule dir x hull edge dirs (world space) ---
	float edge_sep = -1e18f;
	int edge_idx = -1;
	v3 edge_axis = V3(0,0,0);
	v3 edge_pt_world = V3(0,0,0);
	if (cap_len2 > 1e-12f) {
		v3 cd = scale(cap_dir, 1.0f / sqrtf(cap_len2));
		for (int i = 0; i < hull->edge_count; i += 2) {
			v3 ev0 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edges[i].origin], b.scale)));
			v3 ev1 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edges[hull->edges[i].next].origin], b.scale)));
			v3 ed = sub(ev1, ev0);
			v3 ax = cross(cd, ed);
			float al = len2(ax);
			if (al < 1e-12f) continue;
			ax = scale(ax, 1.0f / sqrtf(al));
			if (dot(ax, sub(a.p, ev0)) < 0.0f) ax = neg(ax);
			float cs = dot(a.p, ax) < dot(a.q, ax) ? dot(a.p, ax) : dot(a.q, ax);
			float hs = -1e18f;
			for (int j = 0; j < hull->vert_count; j++) {
				v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
				float d = dot(wv, ax);
				if (d > hs) hs = d;
			}
			float sep = cs - hs;
			if (sep > edge_sep) { edge_sep = sep; edge_idx = i; edge_axis = ax; edge_pt_world = ev0; }
		}
	}

	// Pick axis of minimum penetration.
	float best_sep;
	v3 best_n;
	int use_face;
	if (edge_idx >= 0 && edge_sep > face_sep + 0.001f) {
		best_sep = edge_sep;
		best_n = edge_axis;
		use_face = 0;
	} else {
		best_sep = face_sep;
		best_n = face_plane.normal;
		use_face = 1;
	}

	if (best_sep > a.radius) return 0;
	if (!manifold) return 1;

	// Generate contacts: project capsule endpoints onto the reference plane,
	// keep those that penetrate.
	float plane_d = use_face ? face_plane.offset : dot(best_n, edge_pt_world);
	float dp = dot(a.p, best_n) - plane_d;
	float dq = dot(a.q, best_n) - plane_d;

	int cp = 0;
	v3 points[2];
	float depths[2];
	if (dp < a.radius) {
		points[cp] = sub(a.p, scale(best_n, a.radius));
		depths[cp] = a.radius - dp;
		cp++;
	}
	if (dq < a.radius) {
		points[cp] = sub(a.q, scale(best_n, a.radius));
		depths[cp] = a.radius - dq;
		cp++;
	}

	if (cp == 0) return 0;
	manifold->count = cp;
	for (int i = 0; i < cp; i++) {
		manifold->contacts[i] = (Contact){
			.point = points[i],
			.normal = neg(best_n),
			.penetration = depths[i],
		};
	}
	return 1;
}

int collide_capsule_box(Capsule a, Box b, Manifold* manifold)
{
	return collide_capsule_hull(a, (ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents }, manifold);
}

// -----------------------------------------------------------------------------
// SAT: face queries.
// All computations in local space of Hull2.

typedef struct FaceQuery
{
	int index;
	float separation;
} FaceQuery;

static FaceQuery sat_query_faces(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
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
static float sat_edge_project_full(v3 e1, v3 e2, v3 c1,
	const Hull* hull1, quat rel_rot, v3 scale1,
	const Hull* hull2, v3 scale2)
{
	v3 e1_x_e2 = cross(e1, e2);
	float l = len(e1_x_e2);

	// Skip near-parallel edges
	float tolerance = 0.005f;
	if (l < tolerance * sqrtf(len2(e1) * len2(e2)))
		return -1e18f;

	v3 n = scale(e1_x_e2, 1.0f / l);
	// Orient n from hull1 toward hull2 by projecting ALL vertices of both hulls.
	// This is O(V) per edge pair but robust for any scale/aspect ratio.
	float max1 = -1e18f, min2 = 1e18f;
	for (int i = 0; i < hull1->vert_count; i++) {
		v3 v = add(c1, rotate(rel_rot, hull_vert_scaled(hull1, i, scale1)));
		float d = dot(n, v);
		if (d > max1) max1 = d;
	}
	for (int i = 0; i < hull2->vert_count; i++) {
		v3 v = hull_vert_scaled(hull2, i, scale2);
		float d = dot(n, v);
		if (d < min2) min2 = d;
	}
	return min2 - max1;
}

// SAT: edge queries with Gauss map pruning.
typedef struct EdgeQuery
{
	int index1;
	int index2;
	float separation;
} EdgeQuery;

static EdgeQuery sat_query_edges(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2)
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

			float sep = sat_edge_project_full(e1, e2, c1_local, hull1, rel_rot, scale1, hull2, scale2);
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
static int hull_face_verts_world(const Hull* hull, int face_idx, v3 pos, quat rot, v3 sc, v3* out)
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
static int find_incident_face(const Hull* hull, v3 pos, quat rot, v3 sc, v3 ref_normal)
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
// Tracks feature IDs: clip_edge is the side plane index that clips new verts.
static int clip_to_plane(v3* in, uint8_t* in_fid, int in_count, v3 plane_n, float plane_d, uint8_t clip_edge, v3* out, uint8_t* out_fid)
{
	if (in_count == 0) return 0;
	int out_count = 0;
	v3 prev = in[in_count - 1];
	uint8_t fid_prev = in_fid[in_count - 1];
	float d_prev = dot(plane_n, prev) - plane_d;

	for (int i = 0; i < in_count; i++) {
		v3 cur = in[i];
		uint8_t fid_cur = in_fid[i];
		float d_cur = dot(plane_n, cur) - plane_d;

		if (d_prev <= 0.0f) {
			out[out_count] = prev;
			out_fid[out_count] = fid_prev;
			out_count++;
			if (d_cur > 0.0f) {
				// Exiting clip plane: tag with 0x40 bit to distinguish from entry
				float t = d_prev / (d_prev - d_cur);
				out[out_count] = add(prev, scale(sub(cur, prev), t));
				out_fid[out_count] = clip_edge | 0x40;
				out_count++;
			}
		} else if (d_cur <= 0.0f) {
			// Entering clip plane: plain clip_edge
			float t = d_prev / (d_prev - d_cur);
			out[out_count] = add(prev, scale(sub(cur, prev), t));
			out_fid[out_count] = clip_edge;
			out_count++;
		}

		prev = cur;
		fid_prev = fid_cur;
		d_prev = d_cur;
	}
	return out_count;
}

// Reduce contact manifold to MAX_CONTACTS (4) points.
// 1) Two farthest points (span the patch).
// 2) Farthest from that segment (adds width).
// 3) Maximal barycentric contributor (most outside triangle, completes quad).
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

	Contact tmp[MAX_CONTACTS];
	for (int i = 0; i < MAX_CONTACTS; i++) tmp[i] = contacts[sel[i]];
	for (int i = 0; i < MAX_CONTACTS; i++) contacts[i] = tmp[i];
	return MAX_CONTACTS;
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

	// NaN transforms cause face_index to stay -1 (NaN comparisons always false)
	if (face_a.index < 0 || face_b.index < 0) return 0;

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
			.feature_id = FEATURE_EDGE_BIT | (uint32_t)edge_q.index1 | ((uint32_t)edge_q.index2 << 16),
		};
		return 1;
	}

	// --- Face contact with Sutherland-Hodgman clipping ---
	const Hull* ref_hull; const Hull* inc_hull;
	v3 ref_pos, inc_pos, ref_sc, inc_sc;
	quat ref_rot, inc_rot;
	int ref_face, flip;

	// Bias toward shape A as reference face: B must be significantly better to win.
	// This keeps feature IDs stable across frames for nearly-parallel face pairs.
	float face_bias = 0.98f * face_a.separation + 0.001f;
	if (face_b.separation > face_bias) {
		ref_hull = hull_b; ref_pos = pos_b; ref_rot = rot_b; ref_sc = scale_b;
		inc_hull = hull_a; inc_pos = pos_a; inc_rot = rot_a; inc_sc = scale_a;
		ref_face = face_b.index; flip = 1;
	} else {
		ref_hull = hull_a; ref_pos = pos_a; ref_rot = rot_a; ref_sc = scale_a;
		inc_hull = hull_b; inc_pos = pos_b; inc_rot = rot_b; inc_sc = scale_b;
		ref_face = face_a.index; flip = 0;
	}

	HullPlane ref_plane = plane_transform(ref_hull->planes[ref_face], ref_pos, ref_rot, ref_sc);

	// Collect incident face vertices with initial feature IDs.
	int inc_face = find_incident_face(inc_hull, inc_pos, inc_rot, inc_sc, ref_plane.normal);
	v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
	uint8_t fid1[MAX_CLIP_VERTS], fid2[MAX_CLIP_VERTS];
	int clip_count = hull_face_verts_world(inc_hull, inc_face, inc_pos, inc_rot, inc_sc, buf1);
	// Initial feature: each incident vertex gets a unique ID based on its
	// position in the face winding. This ensures face-face contacts where no
	// clipping or corner-snapping occurs still get distinct feature IDs for
	// warm-start matching. Range 0x80..0xBF avoids collision with clip-edge
	// indices (0..N) and corner-snap tags (0xC0..0xFF).
	for (int i = 0; i < clip_count; i++) fid1[i] = 0x80 | (uint8_t)i;

	// Clip against each side plane of the reference face.
	v3* in_buf = buf1;   v3* out_buf = buf2;
	uint8_t* in_fid = fid1; uint8_t* out_fid = fid2;
	int start_e = ref_hull->faces[ref_face].edge;
	int ei = start_e;
	int guard = 0;
	uint8_t clip_edge_idx = 0;
	do {
		const HalfEdge* edge = &ref_hull->edges[ei];
		v3 tail = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, edge->origin, ref_sc)));
		v3 head = add(ref_pos, rotate(ref_rot,
			hull_vert_scaled(ref_hull, ref_hull->edges[edge->twin].origin, ref_sc)));
		v3 side_n = norm(cross(sub(head, tail), ref_plane.normal));
		float side_d = dot(side_n, tail);

		clip_count = clip_to_plane(in_buf, in_fid, clip_count,
			side_n, side_d, clip_edge_idx, out_buf, out_fid);
		v3* swap = in_buf; in_buf = out_buf; out_buf = swap;
		uint8_t* fswap = in_fid; in_fid = out_fid; out_fid = fswap;

		clip_edge_idx++;
		ei = edge->next;
		assert(++guard < MAX_CLIP_VERTS && "collide_hull_hull: face edge loop didn't close");
	} while (ei != start_e);

	// Snap clipped vertices near reference face corners to canonical corner IDs.
	// A vertex at a corner sits on two side planes -- which plane's clip_edge it
	// gets is FP-dependent and flickers between frames. Replacing the clip_edge
	// with a deterministic corner tag (0xC0 | corner_index) stabilises the
	// feature ID so warm starting can match it.
	{
		v3 corners[MAX_CLIP_VERTS];
		int ncorners = 0;
		int ce = start_e;
		do {
			corners[ncorners++] = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edges[ce].origin, ref_sc)));
			ce = ref_hull->edges[ce].next;
		} while (ce != start_e);
		float snap_tol2 = 1e-6f;
		for (int i = 0; i < clip_count; i++) {
			for (int c = 0; c < ncorners; c++) {
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) {
					in_fid[i] = 0xC0 | (uint8_t)c;
					break;
				}
			}
		}
	}

	// Keep points within margin of reference face, generate contacts with feature IDs.
	// LINEAR_SLOP margin keeps contacts alive near zero-penetration, preventing blink.
	// Feature ID encodes: ref_face | (inc_face << 8) | (clip_edge << 16).
	// flip bit: if ref was hull_b, swap the face roles so the ID is canonical.
	v3 contact_n = flip ? neg(ref_plane.normal) : ref_plane.normal;
	Contact tmp_contacts[MAX_CLIP_VERTS];
	int cp = 0;
	for (int i = 0; i < clip_count; i++) {
		float depth = ref_plane.offset - dot(ref_plane.normal, in_buf[i]);
		if (depth >= -LINEAR_SLOP) {
			uint32_t fid;
			if (!flip)
				fid = (uint32_t)ref_face | ((uint32_t)inc_face << 8) | ((uint32_t)in_fid[i] << 16);
			else
				fid = (uint32_t)inc_face | ((uint32_t)ref_face << 8) | ((uint32_t)in_fid[i] << 16);
			tmp_contacts[cp++] = (Contact){
				.point = in_buf[i],
				.normal = contact_n,
				.penetration = depth,
				.feature_id = fid,
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

// Narrowphase dispatch for a single body pair.
static void narrowphase_pair(WorldInternal* w, int i, int j, InternalManifold** manifolds)
{
	// Canonical ordering: lower type first for upper-triangle dispatch,
	// lower body index first for same-type pairs (deterministic shape A).
	if (w->body_cold[i].shapes[0].type > w->body_cold[j].shapes[0].type || (w->body_cold[i].shapes[0].type == w->body_cold[j].shapes[0].type && i > j)) {
		int tmp = i; i = j; j = tmp;
	}

	ShapeInternal* s0 = &w->body_cold[i].shapes[0];
	ShapeInternal* s1 = &w->body_cold[j].shapes[0];
	BodyHot* h0 = &w->body_hot[i];
	BodyHot* h1 = &w->body_hot[j];
	InternalManifold im = { .body_a = i, .body_b = j };

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
		hit = collide_hull_hull((ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents }, make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL)
		hit = collide_sphere_hull(make_sphere(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL)
		hit = collide_capsule_hull(make_capsule(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL)
		hit = collide_hull_hull(make_convex_hull(h0, s0), make_convex_hull(h1, s1), &im.m);

	if (hit) apush(*manifolds, im);
}

static uint64_t body_pair_key(int a, int b)
{
	uint32_t lo = a < b ? a : b;
	uint32_t hi = a < b ? b : a;
	return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int jointed_pair_skip(CK_MAP(uint8_t) joint_pairs, int a, int b)
{
	if (!joint_pairs) return 0;
	uint64_t key = body_pair_key(a, b);
	return map_get_ptr(joint_pairs, key) != NULL;
}

static void broadphase_n2(WorldInternal* w, InternalManifold** manifolds)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (asize(w->body_cold[i].shapes) == 0) continue;
		for (int j = i + 1; j < count; j++) {
			if (!split_alive(w->body_gen, j)) continue;
			if (asize(w->body_cold[j].shapes) == 0) continue;
			if (w->body_hot[i].inv_mass == 0.0f && w->body_hot[j].inv_mass == 0.0f) continue;
			// Skip sleeping-vs-sleeping pairs
			int isl_i = w->body_cold[i].island_id, isl_j = w->body_cold[j].island_id;
			if (isl_i >= 0 && isl_j >= 0 && (w->island_gen[isl_i] & 1) && (w->island_gen[isl_j] & 1) && !w->islands[isl_i].awake && !w->islands[isl_j].awake) continue;
			if (jointed_pair_skip(w->joint_pairs, i, j)) continue;
			narrowphase_pair(w, i, j, manifolds);
		}
	}
}

// Build AABB lookup table indexed by leaf index from current node contents.
static AABB* bvh_build_lut(BVHTree* t)
{
	int lcount = asize(t->leaves);
	if (lcount == 0) return NULL;
	AABB* lut = CK_ALLOC(sizeof(AABB) * lcount);
	for (int i = 0; i < lcount; i++) {
		BVHLeaf* lf = &t->leaves[i];
		BVHChild* c = bvh_child(&t->nodes[lf->node_idx], lf->child_slot);
		lut[i] = bvh_child_aabb(c);
	}
	return lut;
}

static void broadphase_bvh(WorldInternal* w, InternalManifold** manifolds)
{
	bvh_refit(w->bvh_dynamic, w);

	// Incremental refinement: rebuild a subtree using binned SAH.
	AABB* lut = bvh_build_lut(w->bvh_dynamic);
	if (lut) { bvh_incremental_refine(w->bvh_dynamic, lut); CK_FREE(lut); }

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(w->bvh_dynamic, &pairs);
	bvh_cross_test(w->bvh_dynamic, w->bvh_static, &pairs);

	for (int i = 0; i < asize(pairs); i++) {
		int a = pairs[i].a, b = pairs[i].b;
		if (w->body_hot[a].inv_mass == 0.0f && w->body_hot[b].inv_mass == 0.0f) continue;
		int isl_a = w->body_cold[a].island_id, isl_b = w->body_cold[b].island_id;
		if (isl_a >= 0 && isl_b >= 0 && (w->island_gen[isl_a] & 1) && (w->island_gen[isl_b] & 1) && !w->islands[isl_a].awake && !w->islands[isl_b].awake) continue;
		if (jointed_pair_skip(w->joint_pairs, a, b)) continue;
		narrowphase_pair(w, a, b, manifolds);
	}

	afree(pairs);
}

// Rebuild joint_pairs map when joint topology changes (create/destroy joint).
// O(joint_count) build, O(1) lookup per broadphase pair.
static void joint_pairs_refresh(WorldInternal* w)
{
	if (w->joint_pairs_version == w->ldl_topo_version) return;
	map_free(w->joint_pairs);
	w->joint_pairs = NULL;
	int jcount = asize(w->joints);
	for (int i = 0; i < jcount; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		uint64_t key = body_pair_key(j->body_a, j->body_b);
		map_set(w->joint_pairs, key, (uint8_t)1);
	}
	w->joint_pairs_version = w->ldl_topo_version;
}

static void broadphase_and_collide(WorldInternal* w, InternalManifold** manifolds)
{
	joint_pairs_refresh(w);
	if (w->broadphase_type == BROADPHASE_BVH) broadphase_bvh(w, manifolds);
	else broadphase_n2(w, manifolds);
}

// -----------------------------------------------------------------------------
// Ray-shape intersection. Each returns 1 on hit, sets t_out and n_out.
// Direction must be normalized. Rejects hits behind the origin (t < 0).

static int ray_sphere(v3 ro, v3 rd, Sphere s, float max_t, float* t_out, v3* n_out)
{
	v3 oc = sub(ro, s.center);
	float b = dot(oc, rd);
	float c = dot(oc, oc) - s.radius * s.radius;
	float disc = b * b - c;
	if (disc < 0.0f) return 0;
	float sq = sqrtf(disc);
	float t = -b - sq;
	if (t < 0.0f) t = -b + sq;
	if (t < 0.0f || t > max_t) return 0;
	*t_out = t;
	*n_out = norm(sub(add(ro, scale(rd, t)), s.center));
	return 1;
}

static int ray_capsule(v3 ro, v3 rd, Capsule cap, float max_t, float* t_out, v3* n_out)
{
	v3 ba = sub(cap.q, cap.p);
	v3 oa = sub(ro, cap.p);
	float baba = dot(ba, ba);
	float bard = dot(ba, rd);
	float baoa = dot(ba, oa);
	float rdoa = dot(rd, oa);
	float oaoa = dot(oa, oa);
	float r2 = cap.radius * cap.radius;

	// Infinite cylinder test.
	float a = baba - bard * bard;
	float b = baba * rdoa - baoa * bard;
	float c = baba * oaoa - baoa * baoa - r2 * baba;
	if (fabsf(a) > 1e-8f) {
		float h = b * b - a * c;
		if (h >= 0.0f) {
			float t = (-b - sqrtf(h)) / a;
			if (t >= 0.0f && t <= max_t) {
				float y = baoa + t * bard;
				if (y > 0.0f && y < baba) {
					*t_out = t;
					v3 hp = add(ro, scale(rd, t));
					*n_out = norm(sub(hp, add(cap.p, scale(ba, y / baba))));
					return 1;
				}
			}
		}
	}

	// Hemisphere caps.
	float best = max_t + 1.0f;
	v3 best_n = {0};
	int hit = 0;
	for (int ci = 0; ci < 2; ci++) {
		v3 center = ci == 0 ? cap.p : cap.q;
		v3 oc2 = sub(ro, center);
		float b2 = dot(oc2, rd);
		float c2 = dot(oc2, oc2) - r2;
		float d2 = b2 * b2 - c2;
		if (d2 < 0.0f) continue;
		float t = -b2 - sqrtf(d2);
		if (t < 0.0f) t = -b2 + sqrtf(d2);
		if (t < 0.0f || t >= best) continue;
		float y = baoa + t * bard;
		if (ci == 0 && y > 0.0f) continue;
		if (ci == 1 && y < baba) continue;
		best = t; hit = 1;
		best_n = norm(sub(add(ro, scale(rd, t)), center));
	}
	if (hit && best <= max_t) { *t_out = best; *n_out = best_n; return 1; }
	return 0;
}

static int ray_box(v3 ro, v3 rd, Box box, float max_t, float* t_out, v3* n_out)
{
	quat iq = inv(box.rotation);
	v3 lo = rotate(iq, sub(ro, box.center));
	v3 ld = rotate(iq, rd);
	v3 e = box.half_extents;
	float tmin = -1e30f, tmax = max_t;
	v3 ln = {0};
	float* lp = &lo.x, *dp = &ld.x, *ep = &e.x;
	for (int ax = 0; ax < 3; ax++) {
		if (fabsf(dp[ax]) < 1e-8f) {
			if (lp[ax] < -ep[ax] || lp[ax] > ep[ax]) return 0;
		} else {
			float id = 1.0f / dp[ax];
			float t1 = (-ep[ax] - lp[ax]) * id;
			float t2 = ( ep[ax] - lp[ax]) * id;
			int sw = t1 > t2;
			if (sw) { float tmp = t1; t1 = t2; t2 = tmp; }
			if (t1 > tmin) {
				tmin = t1;
				ln = V3(ax == 0 ? (sw ? 1.0f : -1.0f) : 0.0f,
				        ax == 1 ? (sw ? 1.0f : -1.0f) : 0.0f,
				        ax == 2 ? (sw ? 1.0f : -1.0f) : 0.0f);
			}
			tmax = fminf(tmax, t2);
			if (tmin > tmax) return 0;
		}
	}
	if (tmin < 0.0f) return 0;
	*t_out = tmin;
	*n_out = rotate(box.rotation, ln);
	return 1;
}

static int ray_hull(v3 ro, v3 rd, ConvexHull ch, float max_t, float* t_out, v3* n_out)
{
	quat iq = inv(ch.rotation);
	v3 inv_sc = rcp(ch.scale);
	v3 lo = hmul(rotate(iq, sub(ro, ch.center)), inv_sc);
	v3 ld = hmul(rotate(iq, rd), inv_sc);
	const Hull* hull = ch.hull;
	float t_enter = -1e30f, t_exit = max_t;
	v3 enter_normal = {0};
	for (int i = 0; i < hull->face_count; i++) {
		v3 n = hull->planes[i].normal;
		float d = hull->planes[i].offset;
		float denom = dot(n, ld);
		float num = d - dot(n, lo);
		if (fabsf(denom) < 1e-8f) {
			if (num < 0.0f) return 0;
			continue;
		}
		float t = num / denom;
		if (denom < 0.0f) {
			if (t > t_enter) { t_enter = t; enter_normal = n; }
		} else {
			if (t < t_exit) t_exit = t;
		}
		if (t_enter > t_exit) return 0;
	}
	if (t_enter < 0.0f || t_enter > max_t) return 0;
	*t_out = t_enter;
	*n_out = norm(rotate(ch.rotation, hmul(enter_normal, inv_sc)));
	return 1;
}

// Test a ray against all shapes on a body. Returns closest hit.
static int ray_body(WorldInternal* w, int body_idx, v3 origin, v3 dir, float max_t, float* t_out, v3* n_out)
{
	BodyHot* bh = &w->body_hot[body_idx];
	BodyCold* bc = &w->body_cold[body_idx];
	float best_t = max_t;
	v3 best_n = {0};
	int found = 0;
	for (int i = 0; i < asize(bc->shapes); i++) {
		ShapeInternal* s = &bc->shapes[i];
		float t; v3 n;
		int hit = 0;
		switch (s->type) {
		case SHAPE_SPHERE:  hit = ray_sphere(origin, dir, make_sphere(bh, s), best_t, &t, &n); break;
		case SHAPE_CAPSULE: hit = ray_capsule(origin, dir, make_capsule(bh, s), best_t, &t, &n); break;
		case SHAPE_BOX:     hit = ray_box(origin, dir, make_box(bh, s), best_t, &t, &n); break;
		case SHAPE_HULL:    hit = ray_hull(origin, dir, make_convex_hull(bh, s), best_t, &t, &n); break;
		}
		if (hit && t < best_t) { best_t = t; best_n = n; found = 1; }
	}
	if (found) { *t_out = best_t; *n_out = best_n; }
	return found;
}
