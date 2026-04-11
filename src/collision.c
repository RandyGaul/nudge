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
	// Box fast path: sign selection (3 comparisons vs 8 dot products).
	if (hull->vert_count == 8 && hull->face_count == 6) return V3(dir.x >= 0.0f ? 1.0f : -1.0f, dir.y >= 0.0f ? 1.0f : -1.0f, dir.z >= 0.0f ? 1.0f : -1.0f);
	int n = hull->vert_count;
	if (hull->soa_verts) {
		// SIMD 4-wide support scan using SoA vertex data.
		const float* sx = hull->soa_verts, *sy = sx + ((n + 3) & ~3), *sz = sy + ((n + 3) & ~3);
		simd4f vdx = simd_set1(dir.x), vdy = simd_set1(dir.y), vdz = simd_set1(dir.z);
		simd4f best4 = simd_set1(-1e18f);
		simd4i besti4 = simd_set1_i(0);
		simd4i idx4 = _mm_set_epi32(3, 2, 1, 0);
		simd4i inc4 = simd_set1_i(4);
		int n4 = n & ~3;
		for (int i = 0; i < n4; i += 4) {
			simd4f d = simd_add(simd_add(simd_mul(simd_load(sx + i), vdx), simd_mul(simd_load(sy + i), vdy)), simd_mul(simd_load(sz + i), vdz));
			simd4f mask = simd_cmpgt(d, best4);
			best4 = simd_max(best4, d);
			besti4 = simd_cast_ftoi(simd_blendv(simd_cast_itof(besti4), simd_cast_itof(idx4), mask));
			idx4 = simd_add_i(idx4, inc4);
		}
		// Horizontal reduction.
		float ds[4]; simd_store(ds, best4);
		int idxs[4]; simd_store_i(idxs, besti4);
		int best_i = idxs[0]; float best = ds[0];
		for (int k = 1; k < 4; k++) if (ds[k] > best) { best = ds[k]; best_i = idxs[k]; }
		for (int i = n4; i < n; i++) {
			float d = sx[i]*dir.x + sy[i]*dir.y + sz[i]*dir.z;
			if (d > best) { best = d; best_i = i; }
		}
		return hull->verts[best_i];
	}
	float best = -1e18f;
	int best_i = 0;
	for (int i = 0; i < n; i++) {
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

static const uint16_t s_box_edge_twin[24] = { 1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14, 17,16, 19,18, 21,20, 23,22 };
static const uint16_t s_box_edge_next[24] = { 2,16, 4,22, 6,20, 0,18, 10,17, 12,19, 14,21, 8,23, 15,7, 9,5, 11,3, 13,1 };
static const uint16_t s_box_edge_origin[24] = { 0,3, 3,2, 2,1, 1,0, 4,5, 5,6, 6,7, 7,4, 0,4, 1,5, 2,6, 3,7 };
static const uint16_t s_box_edge_face[24] = { 0,2, 0,5, 0,3, 0,4, 1,4, 1,3, 1,5, 1,2, 2,4, 4,3, 3,5, 5,2 };

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
	.edge_twin = s_box_edge_twin,
	.edge_next = s_box_edge_next,
	.edge_origin = s_box_edge_origin,
	.edge_face = s_box_edge_face,
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
	// Fast path: uniform scale skips inverse-transpose normalization (rotation preserves unit length).
	if (scale.x == scale.y && scale.y == scale.z) {
		v3 n = rotate(rot, p.normal);
		float off = dot(n, pos) + p.offset * scale.x;
		return (HullPlane){ .normal = n, .offset = off };
	}
	v3 n = norm(rotate(rot, V3(p.normal.x / scale.x, p.normal.y / scale.y, p.normal.z / scale.z)));
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
	v3 scaled[MAX_HULL_VERTS]; float soa[MAX_HULL_VERTS*3];
	GJK_Shape a = gjk_sphere(pt, 0);
	GJK_Shape b = gjk_hull_scaled(h.hull, h.center, h.rotation, h.scale, scaled, soa);
	return gjk_distance(&a, &b, NULL);
}

static GJK_Result gjk_query_segment_hull(v3 p, v3 q, ConvexHull h)
{
	v3 scaled[MAX_HULL_VERTS]; float soa[MAX_HULL_VERTS*3];
	GJK_Shape a = gjk_capsule(p, q, 0);
	GJK_Shape b = gjk_hull_scaled(h.hull, h.center, h.rotation, h.scale, scaled, soa);
	return gjk_distance(&a, &b, NULL);
}

static GJK_Result gjk_query_hull_hull(ConvexHull a, ConvexHull b)
{
	v3 sa[MAX_HULL_VERTS], sb[MAX_HULL_VERTS]; float soa_a[MAX_HULL_VERTS*3], soa_b[MAX_HULL_VERTS*3];
	GJK_Shape ga = gjk_hull_scaled(a.hull, a.center, a.rotation, a.scale, sa, soa_a);
	GJK_Shape gb = gjk_hull_scaled(b.hull, b.center, b.rotation, b.scale, sb, soa_b);
	return gjk_distance(&ga, &gb, NULL);
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

	// Deep: center is inside hull. Find face with least penetration via hill-climb.
	// For large hulls (>12 faces), hill-climb from face 0 using topological adjacency.
	float best_sep = -1e18f;
	int best_face = -1;
	HullPlane best_plane = {0};
	if (b.hull->face_count > 12) {
		int cur = 0;
		HullPlane wp = plane_transform(b.hull->planes[0], b.center, b.rotation, b.scale);
		best_sep = dot(a.center, wp.normal) - wp.offset;
		best_face = 0; best_plane = wp;
		for (int iter = 0; iter < b.hull->face_count; iter++) {
			int improved = 0;
			int start_e = b.hull->faces[cur].edge, ei = start_e;
			do {
				int adj = b.hull->edge_face[b.hull->edge_twin[ei]];
				if (adj != cur) {
					HullPlane awp = plane_transform(b.hull->planes[adj], b.center, b.rotation, b.scale);
					float s = dot(a.center, awp.normal) - awp.offset;
					if (s > a.radius) return 0;
					if (s > best_sep) { best_sep = s; best_face = adj; best_plane = awp; improved = 1; }
				}
				ei = b.hull->edge_next[ei];
			} while (ei != start_e);
			if (!improved) break;
			cur = best_face;
		}
	} else {
		for (int i = 0; i < b.hull->face_count; i++) {
			HullPlane wp = plane_transform(b.hull->planes[i], b.center, b.rotation, b.scale);
			float s = dot(a.center, wp.normal) - wp.offset;
			if (s > a.radius) return 0;
			if (s > best_sep) { best_sep = s; best_face = i; best_plane = wp; }
		}
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

// Analytical sphere-box: project sphere center to box local space, clamp to extents, compute distance.
// Much faster than sphere_hull (which uses GJK + plane search).
int collide_sphere_box(Sphere a, Box b, Manifold* manifold)
{
	// Transform sphere center to box local space.
	quat inv_rot = inv(b.rotation);
	v3 local_center = rotate(inv_rot, sub(a.center, b.center));

	// Clamp to box extents = closest point on box surface.
	v3 he = b.half_extents;
	v3 clamped = V3(local_center.x < -he.x ? -he.x : (local_center.x > he.x ? he.x : local_center.x), local_center.y < -he.y ? -he.y : (local_center.y > he.y ? he.y : local_center.y), local_center.z < -he.z ? -he.z : (local_center.z > he.z ? he.z : local_center.z));

	v3 diff = sub(local_center, clamped);
	float dist2 = len2(diff);

	if (dist2 > a.radius * a.radius && dist2 > 1e-12f) return 0;

	if (dist2 > 1e-12f) {
		// Sphere center outside box — contact on box surface.
		float dist = sqrtf(dist2);
		v3 local_n = scale(diff, 1.0f / dist);
		v3 world_n = rotate(b.rotation, local_n);
		v3 world_pt = add(b.center, rotate(b.rotation, clamped));
		manifold->count = 1;
		manifold->contacts[0] = (Contact){ .point = world_pt, .normal = world_n, .penetration = a.radius - dist, .feature_id = 0 };
		return 1;
	}

	// Sphere center inside box — find closest face.
	float best_depth = he.x - fabsf(local_center.x);
	int best_face = local_center.x > 0 ? 3 : 2; // +X:3, -X:2
	float dy = he.y - fabsf(local_center.y);
	if (dy < best_depth) { best_depth = dy; best_face = local_center.y > 0 ? 5 : 4; }
	float dz = he.z - fabsf(local_center.z);
	if (dz < best_depth) { best_depth = dz; best_face = local_center.z > 0 ? 1 : 0; }

	v3 local_n = V3(0, 0, 0);
	static const int face_axis[6] = {2, 2, 0, 0, 1, 1};
	static const float face_sign[6] = {-1, 1, -1, 1, -1, 1};
	(&local_n.x)[face_axis[best_face]] = face_sign[best_face];
	v3 world_n = rotate(b.rotation, local_n);
	v3 world_pt = sub(a.center, scale(world_n, a.radius));
	manifold->count = 1;
	manifold->contacts[0] = (Contact){ .point = world_pt, .normal = world_n, .penetration = a.radius + best_depth, .feature_id = 0 };
	return 1;
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
	// Hill-climb for large hulls, full scan for small.
	float face_sep = -1e18f;
	int face_idx = -1;
	HullPlane face_plane = {0};
	if (hull->face_count > 12) {
		int cur = 0;
		HullPlane wp = plane_transform(hull->planes[0], b.center, b.rotation, b.scale);
		float dp = dot(a.p, wp.normal) - wp.offset, dq = dot(a.q, wp.normal) - wp.offset;
		face_sep = dp < dq ? dp : dq; face_idx = 0; face_plane = wp;
		for (int iter = 0; iter < hull->face_count; iter++) {
			int improved = 0;
			int start_e = hull->faces[cur].edge, ei = start_e;
			do {
				int adj = hull->edge_face[hull->edge_twin[ei]];
				if (adj != cur) {
					HullPlane awp = plane_transform(hull->planes[adj], b.center, b.rotation, b.scale);
					float adp = dot(a.p, awp.normal) - awp.offset, adq = dot(a.q, awp.normal) - awp.offset;
					float asup = adp < adq ? adp : adq;
					if (asup > face_sep) { face_sep = asup; face_idx = adj; face_plane = awp; improved = 1; }
				}
				ei = hull->edge_next[ei];
			} while (ei != start_e);
			if (!improved) break;
			cur = face_idx;
		}
	} else {
		for (int i = 0; i < hull->face_count; i++) {
			HullPlane wp = plane_transform(hull->planes[i], b.center, b.rotation, b.scale);
			float dp = dot(a.p, wp.normal) - wp.offset;
			float dq = dot(a.q, wp.normal) - wp.offset;
			float sup = dp < dq ? dp : dq;
			if (sup > face_sep) { face_sep = sup; face_idx = i; face_plane = wp; }
		}
	}

	// --- Axis family 2: capsule dir x hull edge dirs (world space) ---
	float edge_sep = -1e18f;
	int edge_idx = -1;
	v3 edge_axis = V3(0,0,0);
	v3 edge_pt_world = V3(0,0,0);
	if (cap_len2 > 1e-12f) {
		v3 cd = scale(cap_dir, 1.0f / sqrtf(cap_len2));
		for (int i = 0; i < hull->edge_count; i += 2) {
			v3 ev0 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[i]], b.scale)));
			v3 ev1 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[hull->edge_next[i]]], b.scale)));
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

// Analytical capsule-box: transform capsule to box local space, clamp segment to box,
// compute closest point, then project back. Much faster than capsule_hull (GJK + SAT).
int collide_capsule_box(Capsule a, Box b, Manifold* manifold)
{
	quat inv_rot = inv(b.rotation);
	v3 lp = rotate(inv_rot, sub(a.p, b.center));
	v3 lq = rotate(inv_rot, sub(a.q, b.center));
	v3 he = b.half_extents;

	// Find closest point on segment to box (in box local space).
	// Clamp segment endpoints to box extents to get approximate closest point on box.
	v3 seg_dir = sub(lq, lp);
	float seg_len2 = len2(seg_dir);

	// Find closest point on segment to box surface using iterative clamping.
	// Start with segment closest point to box center, then clamp to box.
	float t = seg_len2 > 1e-12f ? -dot(lp, seg_dir) / seg_len2 : 0.5f;
	t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);
	v3 seg_pt = add(lp, scale(seg_dir, t));

	// Clamp to box = closest point on box to seg_pt.
	v3 box_pt = V3(seg_pt.x < -he.x ? -he.x : (seg_pt.x > he.x ? he.x : seg_pt.x), seg_pt.y < -he.y ? -he.y : (seg_pt.y > he.y ? he.y : seg_pt.y), seg_pt.z < -he.z ? -he.z : (seg_pt.z > he.z ? he.z : seg_pt.z));

	// Re-project: find closest point on segment to box_pt.
	float t2 = seg_len2 > 1e-12f ? dot(sub(box_pt, lp), seg_dir) / seg_len2 : 0.5f;
	t2 = t2 < 0.0f ? 0.0f : (t2 > 1.0f ? 1.0f : t2);
	v3 seg_pt2 = add(lp, scale(seg_dir, t2));

	// Re-clamp box point to segment's new closest point.
	v3 box_pt2 = V3(seg_pt2.x < -he.x ? -he.x : (seg_pt2.x > he.x ? he.x : seg_pt2.x), seg_pt2.y < -he.y ? -he.y : (seg_pt2.y > he.y ? he.y : seg_pt2.y), seg_pt2.z < -he.z ? -he.z : (seg_pt2.z > he.z ? he.z : seg_pt2.z));

	v3 diff = sub(seg_pt2, box_pt2);
	float dist2 = len2(diff);

	if (dist2 > a.radius * a.radius && dist2 > 1e-12f) return 0;

	if (dist2 > 1e-12f) {
		float dist = sqrtf(dist2);
		v3 local_n = scale(diff, 1.0f / dist);
		v3 world_n = rotate(b.rotation, local_n);
		v3 world_pt = add(b.center, rotate(b.rotation, box_pt2));
		manifold->count = 1;
		manifold->contacts[0] = (Contact){ .point = world_pt, .normal = world_n, .penetration = a.radius - dist, .feature_id = 0 };
		return 1;
	}

	// Deep penetration: capsule segment inside box. Use face search like sphere-box.
	float best_depth = he.x - fabsf(seg_pt2.x);
	int best_face = seg_pt2.x > 0 ? 3 : 2;
	float dy = he.y - fabsf(seg_pt2.y);
	if (dy < best_depth) { best_depth = dy; best_face = seg_pt2.y > 0 ? 5 : 4; }
	float dz = he.z - fabsf(seg_pt2.z);
	if (dz < best_depth) { best_depth = dz; best_face = seg_pt2.z > 0 ? 1 : 0; }

	v3 local_n = V3(0, 0, 0);
	static const int face_axis[6] = {2, 2, 0, 0, 1, 1};
	static const float face_sign[6] = {-1, 1, -1, 1, -1, 1};
	(&local_n.x)[face_axis[best_face]] = face_sign[best_face];
	v3 world_n = rotate(b.rotation, local_n);
	v3 world_seg = add(b.center, rotate(b.rotation, seg_pt2));
	v3 world_pt = sub(world_seg, scale(world_n, a.radius));
	manifold->count = 1;
	manifold->contacts[0] = (Contact){ .point = world_pt, .normal = world_n, .penetration = a.radius + best_depth, .feature_id = 0 };
	return 1;
}

// -----------------------------------------------------------------------------
// SAT: face queries.
// All computations in local space of Hull2.

typedef struct FaceQuery
{
	int index;
	float separation;
} FaceQuery;

// Evaluate a single face of hull1 against hull2. Returns separation.
static float sat_eval_face_ex(const Hull* hull1, int face_idx, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, quat inv2, v3 scale2)
{
	HullPlane pw = plane_transform(hull1->planes[face_idx], pos1, rot1, scale1);
	v3 sup_dir_local = rotate(inv2, neg(pw.normal));
	v3 sup_local = hull_support(hull2, sup_dir_local);
	v3 sup_scaled = V3(sup_local.x * scale2.x, sup_local.y * scale2.y, sup_local.z * scale2.z);
	v3 sup_world = add(pos2, rotate(rot2, sup_scaled));
	return dot(pw.normal, sup_world) - pw.offset;
}

// face_hint: if >= 0, hill-climb from this face (for temporal coherence). -1 = full scan.
static FaceQuery sat_query_faces_hint(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2, int face_hint)
{
	FaceQuery best = { .index = -1, .separation = -1e18f };

	// Fast path: box vs any hull. Box has 6 faces with axis-aligned normals.
	if (hull1 == &s_unit_box_hull) {
		v3 cols[3] = { rotate(rot1, V3(1, 0, 0)), rotate(rot1, V3(0, 1, 0)), rotate(rot1, V3(0, 0, 1)) };
		quat inv2 = inv(rot2);
		for (int i = 0; i < 6; i++) {
			HullPlane lp = hull1->planes[i];
			int axis = lp.normal.x != 0.0f ? 0 : (lp.normal.y != 0.0f ? 1 : 2);
			float sign = (&lp.normal.x)[axis];
			v3 face_n = sign > 0 ? cols[axis] : neg(cols[axis]);
			float face_off = dot(face_n, pos1) + (&scale1.x)[axis];
			v3 sup_dir_local = rotate(inv2, neg(face_n));
			v3 sup_local = hull_support(hull2, sup_dir_local);
			v3 sup_scaled = V3(sup_local.x * scale2.x, sup_local.y * scale2.y, sup_local.z * scale2.z);
			v3 sup_world = add(pos2, rotate(rot2, sup_scaled));
			float sep = dot(face_n, sup_world) - face_off;
			if (sep > best.separation) { best.separation = sep; best.index = i; }
		}
		return best;
	}

	// Hill-climb from cached face: walk to topological neighbors with better separation.
	// Precompute inv(rot2) once for all face evaluations.
	quat inv2_pre = inv(rot2);
	if (face_hint >= 0 && face_hint < hull1->face_count && hull1->face_count > 8) {
		int cur = face_hint;
		float cur_sep = sat_eval_face_ex(hull1, cur, pos1, rot1, scale1, hull2, pos2, rot2, inv2_pre, scale2);
		best.index = cur; best.separation = cur_sep;
		for (int iter = 0; iter < hull1->face_count; iter++) {
			int improved = 0;
			int start_e = hull1->faces[cur].edge;
			int ei = start_e;
			do {
				int twin = hull1->edge_twin[ei];
				int adj_face = hull1->edge_face[twin];
				if (adj_face != cur) {
					float adj_sep = sat_eval_face_ex(hull1, adj_face, pos1, rot1, scale1, hull2, pos2, rot2, inv2_pre, scale2);
					if (adj_sep > best.separation) { best.separation = adj_sep; best.index = adj_face; improved = 1; }
				}
				ei = hull1->edge_next[ei];
			} while (ei != start_e);
			if (!improved) return best;
			cur = best.index;
		}
		return best;
	}

	for (int i = 0; i < hull1->face_count; i++) {
		HullPlane pw = plane_transform(hull1->planes[i], pos1, rot1, scale1);
		v3 sup_dir_local = rotate(inv2_pre, neg(pw.normal));
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

static FaceQuery sat_query_faces(const Hull* hull1, v3 pos1, quat rot1, v3 scale1, const Hull* hull2, v3 pos2, quat rot2, v3 scale2) { return sat_query_faces_hint(hull1, pos1, rot1, scale1, hull2, pos2, rot2, scale2, -1); }

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
// Uses support functions instead of full vertex scan — O(1) for boxes.
static float sat_edge_project_full(v3 e1, v3 e2, v3 c1,
	const Hull* hull1, quat rel_rot, v3 scale1,
	const Hull* hull2, v3 scale2)
{
	v3 e1_x_e2 = cross(e1, e2);
	float l2 = len2(e1_x_e2);

	// Skip near-parallel edges (squared comparison avoids 2 sqrt calls).
	float tol2 = 0.005f * 0.005f;
	if (l2 < tol2 * len2(e1) * len2(e2))
		return -1e18f;

	// Fast approximate normalize: rsqrt (1 Newton-Raphson iteration, ~5 cycles vs ~20 for sqrt+div).
	float inv_l;
#if SIMD_SSE
	{ __m128 v = _mm_set_ss(l2); v = _mm_rsqrt_ss(v); _mm_store_ss(&inv_l, v); }
	// One Newton-Raphson refinement: inv_l = inv_l * (1.5 - 0.5 * l2 * inv_l * inv_l)
	inv_l = inv_l * (1.5f - 0.5f * l2 * inv_l * inv_l);
#else
	inv_l = 1.0f / sqrtf(l2);
#endif
	v3 n = scale(e1_x_e2, inv_l);

	// For box hulls, use O(1) support function. For general hulls, keep O(V) vertex scan
	// (avoids extra inv_rot overhead per edge pair).
	if (hull1->vert_count == 8 && hull1->face_count == 6 && hull2->vert_count == 8 && hull2->face_count == 6) {
		quat inv_rel = inv(rel_rot);
		v3 sup1 = hull_support(hull1, rotate(inv_rel, n));
		float max1 = dot(n, add(c1, rotate(rel_rot, V3(sup1.x*scale1.x, sup1.y*scale1.y, sup1.z*scale1.z))));
		v3 sup2 = hull_support(hull2, neg(n));
		float min2 = dot(n, V3(sup2.x*scale2.x, sup2.y*scale2.y, sup2.z*scale2.z));
		return min2 - max1;
	}

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

	// Precompute both hulls' edge data into contiguous arrays.
	int n1 = hull1->edge_count / 2, n2 = hull2->edge_count / 2;
	assert(n1 <= 128 && n2 <= 128);
	v3 e1_arr[128], u1_arr[128], v1_arr[128], ne1_arr[128];
	v3 e2_arr[128], nu2_arr[128], nv2_arr[128], ne2_arr[128];
	for (int k = 0; k < n1; k++) {
		int i = k * 2;
		v3 p = add(c1_local, rotate(rel_rot, hull_vert_scaled(hull1, hull1->edge_origin[i], scale1)));
		v3 q = add(c1_local, rotate(rel_rot, hull_vert_scaled(hull1, hull1->edge_origin[i+1], scale1)));
		e1_arr[k] = sub(q, p);
		u1_arr[k] = rotate(rel_rot, hull1->planes[hull1->edge_face[i]].normal);
		v1_arr[k] = rotate(rel_rot, hull1->planes[hull1->edge_face[i+1]].normal);
		ne1_arr[k] = neg(e1_arr[k]);
	}
	for (int k = 0; k < n2; k++) {
		int i = k * 2;
		v3 p = hull_vert_scaled(hull2, hull2->edge_origin[i], scale2);
		v3 q = hull_vert_scaled(hull2, hull2->edge_origin[i+1], scale2);
		e2_arr[k] = sub(q, p);
		nu2_arr[k] = neg(hull2->planes[hull2->edge_face[i]].normal);
		nv2_arr[k] = neg(hull2->planes[hull2->edge_face[i+1]].normal);
		ne2_arr[k] = neg(e2_arr[k]);
	}

	// Transpose hull2 edge normals into SoA for SIMD Gauss map pruning.
	float nu2x[128], nu2y[128], nu2z[128], nv2x[128], nv2y[128], nv2z[128];
	float ne2x[128], ne2y[128], ne2z[128];
	for (int k = 0; k < n2; k++) {
		nu2x[k] = nu2_arr[k].x; nu2y[k] = nu2_arr[k].y; nu2z[k] = nu2_arr[k].z;
		nv2x[k] = nv2_arr[k].x; nv2y[k] = nv2_arr[k].y; nv2z[k] = nv2_arr[k].z;
		ne2x[k] = ne2_arr[k].x; ne2y[k] = ne2_arr[k].y; ne2z[k] = ne2_arr[k].z;
	}
	// Pad to multiple of 4.
	for (int k = n2; k < ((n2 + 3) & ~3); k++) {
		nu2x[k]=nu2y[k]=nu2z[k]=nv2x[k]=nv2y[k]=nv2z[k]=ne2x[k]=ne2y[k]=ne2z[k]=0;
	}

	for (int k1 = 0; k1 < n1; k1++) {
		v3 e1 = e1_arr[k1], u1 = u1_arr[k1], v1 = v1_arr[k1], b_x_a = ne1_arr[k1];
		// Broadcast k1 data for SIMD.
		simd4f bxa_x = simd_set1(b_x_a.x), bxa_y = simd_set1(b_x_a.y), bxa_z = simd_set1(b_x_a.z);
		simd4f u1x = simd_set1(u1.x), u1y = simd_set1(u1.y), u1z = simd_set1(u1.z);
		simd4f v1x = simd_set1(v1.x), v1y = simd_set1(v1.y), v1z = simd_set1(v1.z);
		simd4f zero = simd_zero();

		int n2_4 = n2 & ~3;
		for (int k2 = 0; k2 < n2_4; k2 += 4) {
			// 4-wide Gauss map pruning.
			simd4f cba = simd_add(simd_add(simd_mul(simd_load(nu2x+k2), bxa_x), simd_mul(simd_load(nu2y+k2), bxa_y)), simd_mul(simd_load(nu2z+k2), bxa_z));
			simd4f dba = simd_add(simd_add(simd_mul(simd_load(nv2x+k2), bxa_x), simd_mul(simd_load(nv2y+k2), bxa_y)), simd_mul(simd_load(nv2z+k2), bxa_z));
			simd4f adc = simd_add(simd_add(simd_mul(u1x, simd_load(ne2x+k2)), simd_mul(u1y, simd_load(ne2y+k2))), simd_mul(u1z, simd_load(ne2z+k2)));
			simd4f bdc = simd_add(simd_add(simd_mul(v1x, simd_load(ne2x+k2)), simd_mul(v1y, simd_load(ne2y+k2))), simd_mul(v1z, simd_load(ne2z+k2)));
			// Check: (cba*dba < 0) && (adc*bdc < 0) && (cba*bdc > 0)
			simd4f cd_neg = simd_cmpgt(zero, simd_mul(cba, dba));
			simd4f ab_neg = simd_cmpgt(zero, simd_mul(adc, bdc));
			simd4f cb_pos = simd_cmpgt(simd_mul(cba, bdc), zero);
			int mask = simd_movemask(simd_and(simd_and(cd_neg, ab_neg), cb_pos));
			if (!mask) continue;
			for (int lane = 0; lane < 4; lane++) {
				if (!(mask & (1 << lane))) continue;
				int k2i = k2 + lane;
				float sep = sat_edge_project_full(e1, e2_arr[k2i], c1_local, hull1, rel_rot, scale1, hull2, scale2);
				if (sep > best.separation) { best.index1 = k1*2; best.index2 = k2i*2; best.separation = sep; }
			}
		}
		for (int k2 = n2_4; k2 < n2; k2++) {
			v3 ne2 = ne2_arr[k2];
			float cba = dot(nu2_arr[k2], b_x_a), dba = dot(nv2_arr[k2], b_x_a);
			float adc = dot(u1, ne2), bdc = dot(v1, ne2);
			if (!((cba * dba < 0.0f) && (adc * bdc < 0.0f) && (cba * bdc > 0.0f))) continue;
			float sep = sat_edge_project_full(e1, e2_arr[k2], c1_local, hull1, rel_rot, scale1, hull2, scale2);
			if (sep > best.separation) { best.index1 = k1*2; best.index2 = k2*2; best.separation = sep; }
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
	// Fast path: box face vertices from rotation columns (avoids 4 quaternion rotations).
	if (hull->verts == s_box_verts) {
		v3 cx = scale(rotate(rot, V3(1, 0, 0)), sc.x), cy = scale(rotate(rot, V3(0, 1, 0)), sc.y), cz = scale(rotate(rot, V3(0, 0, 1)), sc.z);
		// Walk the face winding to get vertices in correct order.
		int start = hull->faces[face_idx].edge;
		int e = start;
		int count = 0;
		do {
			v3 lv = hull->verts[hull->edge_origin[e]];
			out[count++] = add(pos, add(add(scale(cx, lv.x), scale(cy, lv.y)), scale(cz, lv.z)));
			e = hull->edge_next[e];
		} while (e != start);
		return count;
	}
	int start = hull->faces[face_idx].edge;
	int e = start;
	int count = 0;
	do {
		out[count++] = add(pos, rotate(rot, hull_vert_scaled(hull, hull->edge_origin[e], sc)));
		e = hull->edge_next[e];
	} while (e != start && count < MAX_CLIP_VERTS);
	return count;
}

// Find face on hull most anti-parallel to a world-space normal.
static int find_incident_face(const Hull* hull, v3 pos, quat rot, v3 sc, v3 ref_normal)
{
	// Fast path: box faces are rotation columns. 3 dot products vs 6 plane_transforms.
	if (hull->verts == s_box_verts) {
		v3 cols[3] = { rotate(rot, V3(1, 0, 0)), rotate(rot, V3(0, 1, 0)), rotate(rot, V3(0, 0, 1)) };
		// Box faces: -Z=0,+Z=1, -X=2,+X=3, -Y=4,+Y=5. Normals: +-col[2], +-col[0], +-col[1].
		static const int axis_to_neg_face[3] = {2, 4, 0}; // -X=2, -Y=4, -Z=0
		int best = 0; float best_dot = 1e18f;
		for (int i = 0; i < 3; i++) {
			float d = dot(cols[i], ref_normal);
			// Positive face: d itself. Negative face: -d.
			if (d < best_dot) { best_dot = d; best = axis_to_neg_face[i] + 1; } // positive face
			if (-d < best_dot) { best_dot = -d; best = axis_to_neg_face[i]; }   // negative face
		}
		return best;
	}
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

// Generate face-face contact between two convex hulls using Sutherland-Hodgman clipping.
// ref_face is the separating face on (ref_hull, ref_pos, ref_rot, ref_sc).
// flip=1 means ref is actually hull B (so normal is negated for A->B convention).
static int generate_face_contact(const Hull* ref_hull, v3 ref_pos, quat ref_rot, v3 ref_sc, const Hull* inc_hull, v3 inc_pos, quat inc_rot, v3 inc_sc, int ref_face, int flip, Manifold* manifold)
{
	HullPlane ref_plane = plane_transform(ref_hull->planes[ref_face], ref_pos, ref_rot, ref_sc);
	int inc_face = find_incident_face(inc_hull, inc_pos, inc_rot, inc_sc, ref_plane.normal);
	v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
	uint8_t fid1[MAX_CLIP_VERTS], fid2[MAX_CLIP_VERTS];
	int clip_count = hull_face_verts_world(inc_hull, inc_face, inc_pos, inc_rot, inc_sc, buf1);
	for (int i = 0; i < clip_count; i++) fid1[i] = 0x80 | (uint8_t)i;

	v3* in_buf = buf1; v3* out_buf = buf2;
	uint8_t* in_fid = fid1; uint8_t* out_fid = fid2;
	int start_e = ref_hull->faces[ref_face].edge;
	int ei = start_e;
	int guard = 0;
	uint8_t clip_edge_idx = 0;
	do {
		v3 tail = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ei], ref_sc)));
		v3 head = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ref_hull->edge_twin[ei]], ref_sc)));
		v3 side_n = norm(cross(sub(head, tail), ref_plane.normal));
		float side_d = dot(side_n, tail);
		clip_count = clip_to_plane(in_buf, in_fid, clip_count, side_n, side_d, clip_edge_idx, out_buf, out_fid);
		v3* swap = in_buf; in_buf = out_buf; out_buf = swap;
		uint8_t* fswap = in_fid; in_fid = out_fid; out_fid = fswap;
		clip_edge_idx++;
		ei = ref_hull->edge_next[ei];
		assert(++guard < MAX_CLIP_VERTS && "generate_face_contact: face edge loop didn't close");
	} while (ei != start_e);

	{
		v3 corners[MAX_CLIP_VERTS];
		int ncorners = 0;
		int ce = start_e;
		do { corners[ncorners++] = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ce], ref_sc))); ce = ref_hull->edge_next[ce]; } while (ce != start_e);
		float snap_tol2 = 1e-6f;
		for (int i = 0; i < clip_count; i++)
			for (int c = 0; c < ncorners; c++)
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) { in_fid[i] = 0xC0 | (uint8_t)c; break; }
	}

	v3 contact_n = flip ? neg(ref_plane.normal) : ref_plane.normal;
	Contact tmp_contacts[MAX_CLIP_VERTS];
	int cp = 0;
	for (int i = 0; i < clip_count; i++) {
		float depth = ref_plane.offset - dot(ref_plane.normal, in_buf[i]);
		if (depth >= -LINEAR_SLOP) {
			uint32_t fid;
			if (!flip) fid = (uint32_t)ref_face | ((uint32_t)inc_face << 8) | ((uint32_t)in_fid[i] << 16);
			else fid = (uint32_t)inc_face | ((uint32_t)ref_face << 8) | ((uint32_t)in_fid[i] << 16);
			tmp_contacts[cp++] = (Contact){ .point = in_buf[i], .normal = contact_n, .penetration = depth, .feature_id = fid };
		}
	}
	if (cp == 0) return 0;
	cp = reduce_contacts(tmp_contacts, cp);
	manifold->count = cp;
	for (int i = 0; i < cp; i++) manifold->contacts[i] = tmp_contacts[i];
	return 1;
}

// -----------------------------------------------------------------------------
// Full SAT hull vs hull with Sutherland-Hodgman face clipping.

int collide_hull_hull_ex(ConvexHull a, ConvexHull b, Manifold* manifold, int* sat_hint)
{
	const Hull* hull_a = a.hull;
	v3 pos_a = a.center; quat rot_a = a.rotation; v3 scale_a = a.scale;
	const Hull* hull_b = b.hull;
	v3 pos_b = b.center; quat rot_b = b.rotation; v3 scale_b = b.scale;

	// Pass cached face hints to the face queries for hill-climb warm-start.
	int face_hint_a = -1, face_hint_b = -1;
	if (sat_hint && *sat_hint >= 0) {
		int h = *sat_hint;
		if (h < hull_a->face_count) face_hint_a = h;
		else { int fi = h - hull_a->face_count; if (fi < hull_b->face_count) face_hint_b = fi; }
	}

	FaceQuery face_a = sat_query_faces_hint(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b, face_hint_a);
	if (face_a.separation > 0.0f) { if (sat_hint) *sat_hint = face_a.index; return 0; }

	FaceQuery face_b = sat_query_faces_hint(hull_b, pos_b, rot_b, scale_b, hull_a, pos_a, rot_a, scale_a, face_hint_b);
	if (face_b.separation > 0.0f) { if (sat_hint) *sat_hint = hull_a->face_count + face_b.index; return 0; }

	// NaN transforms cause face_index to stay -1 (NaN comparisons always false)
	if (face_a.index < 0 || face_b.index < 0) return 0;

	EdgeQuery edge_q = sat_query_edges(hull_a, pos_a, rot_a, scale_a, hull_b, pos_b, rot_b, scale_b);
	if (edge_q.separation > 0.0f) { if (sat_hint) *sat_hint = -1; return 0; }

	if (!manifold) return 1;

	// Bias toward face contacts over edge contacts
	const float k_tol = 0.05f;
	float max_face_sep = face_a.separation > face_b.separation
		? face_a.separation : face_b.separation;

	if (edge_q.separation > max_face_sep + k_tol) {
		// --- Edge-edge contact ---
		v3 p1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, hull_a->edge_origin[edge_q.index1], scale_a)));
		v3 q1 = add(pos_a, rotate(rot_a, hull_vert_scaled(hull_a, hull_a->edge_origin[edge_q.index1 + 1], scale_a)));

		v3 p2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, hull_b->edge_origin[edge_q.index2], scale_b)));
		v3 q2 = add(pos_b, rotate(rot_b, hull_vert_scaled(hull_b, hull_b->edge_origin[edge_q.index2 + 1], scale_b)));

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
		v3 tail = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ei], ref_sc)));
		v3 head = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ref_hull->edge_twin[ei]], ref_sc)));
		v3 side_n = norm(cross(sub(head, tail), ref_plane.normal));
		float side_d = dot(side_n, tail);

		clip_count = clip_to_plane(in_buf, in_fid, clip_count,
			side_n, side_d, clip_edge_idx, out_buf, out_fid);
		v3* swap = in_buf; in_buf = out_buf; out_buf = swap;
		uint8_t* fswap = in_fid; in_fid = out_fid; out_fid = fswap;

		clip_edge_idx++;
		ei = ref_hull->edge_next[ei];
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
			corners[ncorners++] = add(ref_pos, rotate(ref_rot, hull_vert_scaled(ref_hull, ref_hull->edge_origin[ce], ref_sc)));
			ce = ref_hull->edge_next[ce];
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
	// Cache the winning axis for next-frame warm-start.
	if (sat_hint) *sat_hint = flip ? hull_a->face_count + ref_face : ref_face;
	return 1;
}

int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold)
{
	return collide_hull_hull_ex(a, b, manifold, NULL);
}

// -----------------------------------------------------------------------------
// Dedicated box-box SAT: Gottschalk OBB test with direct rotation column arithmetic.
// Tests all 15 separating axes without generic hull machinery (plane_transform, hull_support loops, Gauss map).
// Falls through to collide_hull_hull for contact generation when penetrating.
// sat_hint: if non-NULL, *sat_hint is the cached axis from last frame (-1 = no cache).
// On return, *sat_hint is updated to the winning axis for next frame.
static int collide_box_box_ex(Box a, Box b, Manifold* manifold, int* sat_hint)
{
	// Rotation columns for each box.
	v3 ax = rotate(a.rotation, V3(1, 0, 0)), ay = rotate(a.rotation, V3(0, 1, 0)), az = rotate(a.rotation, V3(0, 0, 1));
	v3 bx = rotate(b.rotation, V3(1, 0, 0)), by = rotate(b.rotation, V3(0, 1, 0)), bz = rotate(b.rotation, V3(0, 0, 1));

	v3 d = sub(b.center, a.center);
	float ea = a.half_extents.x, eb = a.half_extents.y, ec = a.half_extents.z;
	float fa = b.half_extents.x, fb = b.half_extents.y, fc = b.half_extents.z;

	// R[i][j] = dot(a_axis_i, b_axis_j), absR = |R| + epsilon for parallel edge robustness.
	float R00 = dot(ax, bx), R01 = dot(ax, by), R02 = dot(ax, bz);
	float R10 = dot(ay, bx), R11 = dot(ay, by), R12 = dot(ay, bz);
	float R20 = dot(az, bx), R21 = dot(az, by), R22 = dot(az, bz);
	float eps = 1e-6f;
	float A00 = fabsf(R00)+eps, A01 = fabsf(R01)+eps, A02 = fabsf(R02)+eps;
	float A10 = fabsf(R10)+eps, A11 = fabsf(R11)+eps, A12 = fabsf(R12)+eps;
	float A20 = fabsf(R20)+eps, A21 = fabsf(R21)+eps, A22 = fabsf(R22)+eps;

	// Translation in A's frame.
	float ta = dot(d, ax), tb = dot(d, ay), tc = dot(d, az);

	// Test 15 separating axes. Track minimum penetration.
	float ra, rb, sep, pen, best_pen = 1e18f;
	int best_axis = -1;

	// A's face normals (axes 0-2).
	ra = ea; rb = fa*A00 + fb*A01 + fc*A02; sep = fabsf(ta); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 0; }
	ra = eb; rb = fa*A10 + fb*A11 + fc*A12; sep = fabsf(tb); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 1; }
	ra = ec; rb = fa*A20 + fb*A21 + fc*A22; sep = fabsf(tc); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 2; }

	// B's face normals (axes 3-5).
	ra = ea*A00 + eb*A10 + ec*A20; rb = fa; sep = fabsf(ta*R00 + tb*R10 + tc*R20); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 3; }
	ra = ea*A01 + eb*A11 + ec*A21; rb = fb; sep = fabsf(ta*R01 + tb*R11 + tc*R21); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 4; }
	ra = ea*A02 + eb*A12 + ec*A22; rb = fc; sep = fabsf(ta*R02 + tb*R12 + tc*R22); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 5; }

	// Edge cross products (axes 6-14). Skip near-parallel edges (cross product near zero).
	// ax x bx
	ra = eb*A20 + ec*A10; rb = fb*A02 + fc*A01; sep = fabsf(tc*R10 - tb*R20); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 6; }
	// ax x by
	ra = eb*A21 + ec*A11; rb = fa*A02 + fc*A00; sep = fabsf(tc*R11 - tb*R21); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 7; }
	// ax x bz
	ra = eb*A22 + ec*A12; rb = fa*A01 + fb*A00; sep = fabsf(tc*R12 - tb*R22); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 8; }
	// ay x bx
	ra = ea*A20 + ec*A00; rb = fb*A12 + fc*A11; sep = fabsf(ta*R20 - tc*R00); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 9; }
	// ay x by
	ra = ea*A21 + ec*A01; rb = fa*A12 + fc*A10; sep = fabsf(ta*R21 - tc*R01); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 10; }
	// ay x bz
	ra = ea*A22 + ec*A02; rb = fa*A11 + fb*A10; sep = fabsf(ta*R22 - tc*R02); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 11; }
	// az x bx
	ra = ea*A10 + eb*A00; rb = fb*A22 + fc*A21; sep = fabsf(tb*R00 - ta*R10); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 12; }
	// az x by
	ra = ea*A11 + eb*A01; rb = fa*A22 + fc*A20; sep = fabsf(tb*R01 - ta*R11); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 13; }
	// az x bz
	ra = ea*A12 + eb*A02; rb = fa*A21 + fb*A20; sep = fabsf(tb*R02 - ta*R12); pen = ra + rb - sep; if (pen < 0) return 0; if (pen < best_pen) { best_pen = pen; best_axis = 14; }

	if (!manifold) return 1;

	// Face contact: fully inlined box-box contact gen using rotation columns.
	// No hull infrastructure (plane_transform, half-edge walk, hull_support).
	if (best_axis < 6) {
		v3 ref_cols[3], inc_cols[3], ref_pos, inc_pos, ref_he, inc_he;
		int flip;
		if (best_axis < 3) {
			ref_cols[0] = ax; ref_cols[1] = ay; ref_cols[2] = az; ref_pos = a.center; ref_he = a.half_extents;
			inc_cols[0] = bx; inc_cols[1] = by; inc_cols[2] = bz; inc_pos = b.center; inc_he = b.half_extents;
			flip = 0;
		} else {
			ref_cols[0] = bx; ref_cols[1] = by; ref_cols[2] = bz; ref_pos = b.center; ref_he = b.half_extents;
			inc_cols[0] = ax; inc_cols[1] = ay; inc_cols[2] = az; inc_pos = a.center; inc_he = a.half_extents;
			flip = 1;
		}
		int la = best_axis < 3 ? best_axis : best_axis - 3;
		float proj_sign = best_axis < 3 ? (&((v3){ta, tb, tc}).x)[la] : (&((v3){ta*R00+tb*R10+tc*R20, ta*R01+tb*R11+tc*R21, ta*R02+tb*R12+tc*R22}).x)[la];
		float nsign = proj_sign >= 0.0f ? 1.0f : -1.0f;
		v3 ref_n = scale(ref_cols[la], nsign);
		float ref_off = dot(ref_n, ref_pos) + (&ref_he.x)[la];
		int u = (la + 1) % 3, v_ax = (la + 2) % 3;

		// Incident face: find which inc column is most anti-parallel to ref_n.
		int inc_la = 0; float inc_best = 1e18f;
		for (int i = 0; i < 3; i++) { float d = dot(inc_cols[i], ref_n); if (d < inc_best) { inc_best = d; inc_la = i; } if (-d < inc_best) { inc_best = -d; inc_la = i; } }
		// Determine sign of incident face normal.
		float inc_nsign = dot(inc_cols[inc_la], ref_n) > 0 ? -1.0f : 1.0f;
		int inc_u = (inc_la + 1) % 3, inc_v = (inc_la + 2) % 3;

		// Incident face vertices (4 corners of the incident face).
		v3 inc_center = add(inc_pos, scale(inc_cols[inc_la], inc_nsign * (&inc_he.x)[inc_la]));
		v3 inc_eu = scale(inc_cols[inc_u], (&inc_he.x)[inc_u]), inc_ev = scale(inc_cols[inc_v], (&inc_he.x)[inc_v]);
		v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
		uint8_t fid1[MAX_CLIP_VERTS], fid2[MAX_CLIP_VERTS];
		buf1[0] = add(inc_center, add(inc_eu, inc_ev));
		buf1[1] = add(inc_center, sub(neg(inc_eu), neg(inc_ev)));
		buf1[2] = add(inc_center, sub(neg(inc_eu), inc_ev));
		buf1[3] = add(inc_center, sub(inc_eu, inc_ev));
		// Fix winding: ensure CCW when viewed from ref_n direction.
		if (dot(cross(sub(buf1[1], buf1[0]), sub(buf1[2], buf1[0])), ref_n) > 0) { v3 tmp = buf1[1]; buf1[1] = buf1[3]; buf1[3] = tmp; }
		for (int i = 0; i < 4; i++) fid1[i] = 0x80 | (uint8_t)i;
		int clip_count = 4;

		// 4 side planes from rotation columns (no cross product, no normalize).
		v3 side_n[4] = { ref_cols[u], neg(ref_cols[u]), ref_cols[v_ax], neg(ref_cols[v_ax]) };
		float side_d[4] = { dot(ref_cols[u], ref_pos) + (&ref_he.x)[u], dot(neg(ref_cols[u]), ref_pos) + (&ref_he.x)[u], dot(ref_cols[v_ax], ref_pos) + (&ref_he.x)[v_ax], dot(neg(ref_cols[v_ax]), ref_pos) + (&ref_he.x)[v_ax] };

		v3* in_buf = buf1; v3* out_buf = buf2;
		uint8_t* in_fid = fid1; uint8_t* out_fid = fid2;
		for (int i = 0; i < 4; i++) {
			clip_count = clip_to_plane(in_buf, in_fid, clip_count, side_n[i], side_d[i], (uint8_t)i, out_buf, out_fid);
			v3* sw = in_buf; in_buf = out_buf; out_buf = sw;
			uint8_t* fs = in_fid; in_fid = out_fid; out_fid = fs;
		}

		// Corner snap: 4 reference face corners from rotation columns.
		v3 ref_center = add(ref_pos, scale(ref_cols[la], nsign * (&ref_he.x)[la]));
		v3 ref_eu = scale(ref_cols[u], (&ref_he.x)[u]), ref_ev = scale(ref_cols[v_ax], (&ref_he.x)[v_ax]);
		v3 corners[4] = { add(ref_center, add(ref_eu, ref_ev)), add(ref_center, sub(ref_eu, ref_ev)), add(ref_center, sub(neg(ref_eu), neg(ref_ev))), add(ref_center, sub(neg(ref_eu), ref_ev)) };
		float snap_tol2 = 1e-6f;
		for (int i = 0; i < clip_count; i++)
			for (int c = 0; c < 4; c++)
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) { in_fid[i] = 0xC0 | (uint8_t)c; break; }

		// Build face indices for feature IDs. Map (la, nsign) to hull face index.
		static const int face_map[3][2] = { {2, 3}, {4, 5}, {0, 1} };
		int ref_face = face_map[la][nsign > 0 ? 1 : 0];
		int inc_face = face_map[inc_la][inc_nsign > 0 ? 1 : 0];

		v3 contact_n = flip ? neg(ref_n) : ref_n;
		Contact tmp_contacts[MAX_CLIP_VERTS];
		int cp = 0;
		for (int i = 0; i < clip_count; i++) {
			float depth = ref_off - dot(ref_n, in_buf[i]);
			if (depth >= -LINEAR_SLOP) {
				uint32_t fid = flip ? ((uint32_t)inc_face | ((uint32_t)ref_face << 8) | ((uint32_t)in_fid[i] << 16)) : ((uint32_t)ref_face | ((uint32_t)inc_face << 8) | ((uint32_t)in_fid[i] << 16));
				tmp_contacts[cp++] = (Contact){ .point = in_buf[i], .normal = contact_n, .penetration = depth, .feature_id = fid };
			}
		}
		if (cp == 0) return 0;
		cp = reduce_contacts(tmp_contacts, cp);
		manifold->count = cp;
		for (int i = 0; i < cp; i++) manifold->contacts[i] = tmp_contacts[i];
		if (sat_hint) *sat_hint = best_axis;
		return 1;
	}

	// Edge-edge contact: fall through to hull-hull (edge contacts are rare in box piles).
	if (sat_hint) *sat_hint = best_axis;
	return collide_hull_hull(
		(ConvexHull){ &s_unit_box_hull, a.center, a.rotation, a.half_extents },
		(ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents },
		manifold);
}

int collide_box_box(Box a, Box b, Manifold* manifold) { return collide_box_box_ex(a, b, manifold, NULL); }

const Hull* hull_unit_box() { return &s_unit_box_hull; }

// Quickhull implemented in quickhull.c.

void hull_free(Hull* hull)
{
	if (!hull) return;
	CK_FREE((void*)hull->verts);
	if (hull->soa_verts) CK_FREE_ALIGNED((void*)hull->soa_verts);
	CK_FREE((void*)hull->edge_twin);
	CK_FREE((void*)hull->edge_next);
	CK_FREE((void*)hull->edge_origin);
	CK_FREE((void*)hull->edge_face);
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
// Narrowphase timing accumulators (indexed by shape pair type).
// Pair encoding: type_a * 5 + type_b (upper triangle, type_a <= type_b).
#define NP_PAIR_COUNT 25
static double np_time_acc[NP_PAIR_COUNT];
static int np_call_acc[NP_PAIR_COUNT];
static int np_frame_count;

static int np_pair_idx(int ta, int tb) { return ta * 5 + tb; }

void narrowphase_reset_timers() { memset(np_time_acc, 0, sizeof(np_time_acc)); memset(np_call_acc, 0, sizeof(np_call_acc)); np_frame_count = 0; }
void narrowphase_end_frame() { np_frame_count++; }
void narrowphase_print_timers()
{
	if (np_frame_count == 0) return;
	double n = (double)np_frame_count;
	const char* names[] = {"sphere", "capsule", "box", "hull", "?"};
	printf("  --- narrowphase breakdown ---\n");
	double total = 0;
	int total_calls = 0;
	for (int a = 0; a < 4; a++) for (int b = a; b < 4; b++) {
		int idx = np_pair_idx(a, b);
		if (np_call_acc[idx] == 0) continue;
		double avg_us = np_time_acc[idx] / n * 1e6;
		int avg_calls = (int)((double)np_call_acc[idx] / n + 0.5);
		printf("  np.%-8s-%-8s %7.1f us  (%d calls, %.2f us/call)\n", names[a], names[b], avg_us, avg_calls, avg_us / avg_calls);
		total += np_time_acc[idx] / n * 1000.0;
		total_calls += avg_calls;
	}
	printf("  np.total:          %7.3f ms  (%d calls)\n", total, total_calls);
}

static uint64_t body_pair_key(int a, int b);

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
	double t0 = perf_now();

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
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX) {
		uint64_t pkey = body_pair_key(i, j);
		WarmManifold* wm = map_get_ptr(w->warm_cache, pkey);
		int hint = wm ? wm->sat_axis : -1;
		hit = collide_box_box_ex(make_box(h0, s0), make_box(h1, s1), &im.m, &hint);
		if (wm) wm->sat_axis = hint;
	}
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_HULL) {
		uint64_t pkey = body_pair_key(i, j);
		WarmManifold* wm = map_get_ptr(w->warm_cache, pkey);
		int hint = wm ? wm->sat_axis : -1;
		hit = collide_hull_hull_ex((ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents }, make_convex_hull(h1, s1), &im.m, &hint);
		if (wm) wm->sat_axis = hint;
	}
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL)
		hit = collide_sphere_hull(make_sphere(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL)
		hit = collide_capsule_hull(make_capsule(h0, s0), make_convex_hull(h1, s1), &im.m);
	else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL) {
		uint64_t pkey = body_pair_key(i, j);
		WarmManifold* wm = map_get_ptr(w->warm_cache, pkey);
		int hint = wm ? wm->sat_axis : -1;
		hit = collide_hull_hull_ex(make_convex_hull(h0, s0), make_convex_hull(h1, s1), &im.m, &hint);
		if (wm) wm->sat_axis = hint;
	}

	double dt = perf_now() - t0;
	int idx = np_pair_idx(s0->type, s1->type);
	np_time_acc[idx] += dt;
	np_call_acc[idx]++;

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

// Sweep-and-prune entry for axis-sorted broadphase.
typedef struct SAPEntry { int body_idx; float min_x, max_x; } SAPEntry;
static int sap_cmp(const void* a, const void* b) { float d = ((SAPEntry*)a)->min_x - ((SAPEntry*)b)->min_x; return (d > 0) - (d < 0); }

// Broadphase sub-phase timing accumulators (seconds, summed across frames).
double bp_refit_acc, bp_precomp_acc, bp_sweep_acc, bp_cross_acc;
int bp_frame_count;

static void broadphase_bvh(WorldInternal* w, InternalManifold** manifolds)
{
	double t0 = perf_now();
	bvh_refit(w->bvh_dynamic, w);
	bp_refit_acc += perf_now() - t0;

	// Build tight AABBs + sweep-and-prune entries for all dynamic bodies.
	double t1 = perf_now();
	int body_count = asize(w->body_hot);
	AABB* tight = CK_ALLOC(sizeof(AABB) * body_count);
	CK_DYNA SAPEntry* sap = NULL;
	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i) || asize(w->body_cold[i].shapes) == 0) { tight[i] = aabb_empty(); continue; }
		tight[i] = body_aabb(&w->body_hot[i], &w->body_cold[i]);
		if (w->body_hot[i].inv_mass > 0.0f) apush(sap, ((SAPEntry){ i, tight[i].min.x, tight[i].max.x }));
	}

	// Sort dynamic bodies by x-axis AABB min.
	int sap_count = asize(sap);
	if (sap_count > 1) qsort(sap, sap_count, sizeof(SAPEntry), sap_cmp);
	bp_precomp_acc += perf_now() - t1;

	// Sweep: test overlapping pairs along x-axis, then full 3D AABB overlap (SIMD branchless).
	// Phase 1: collect broadphase pairs (no narrowphase yet — just AABB overlap test).
	double t2 = perf_now();
	CK_DYNA BroadPair* dd_pairs = NULL;
	for (int i = 0; i < sap_count; i++) {
		float max_x = sap[i].max_x;
		int a = sap[i].body_idx;
		AABB ta = tight[a];
		int isl_a = w->body_cold[a].island_id;
		for (int j = i + 1; j < sap_count && sap[j].min_x <= max_x; j++) {
			int b = sap[j].body_idx;
			if (!aabb_overlaps(ta, tight[b])) continue;
			int isl_b = w->body_cold[b].island_id;
			if (isl_a >= 0 && isl_b >= 0 && (w->island_gen[isl_a] & 1) && (w->island_gen[isl_b] & 1) && !w->islands[isl_a].awake && !w->islands[isl_b].awake) continue;
			if (jointed_pair_skip(w->joint_pairs, a, b)) continue;
			apush(dd_pairs, ((BroadPair){ a, b }));
		}
	}

	// Collect d-s pairs too.
	bp_sweep_acc += perf_now() - t2;
	double t3 = perf_now();
	CK_DYNA BroadPair* pairs = NULL;
	bvh_cross_test(w->bvh_dynamic, w->bvh_static, &pairs);
	for (int i = 0; i < asize(pairs); i++) {
		int a = pairs[i].a, b = pairs[i].b;
		if (!aabb_overlaps(tight[a], tight[b])) continue;
		apush(dd_pairs, ((BroadPair){ a, b }));
	}
	afree(pairs);
	bp_cross_acc += perf_now() - t3;

	// Narrowphase on all collected pairs.
	// If parallel dispatch is available (thread_count > 1), output pairs for external dispatch.
	// Otherwise run narrowphase inline.
	if (w->thread_count > 1 && w->np_pairs_out) {
		// Output pairs for parallel narrowphase in nudge.c
		CK_DYNA BroadPair** out = (CK_DYNA BroadPair**)w->np_pairs_out;
		for (int i = 0; i < asize(dd_pairs); i++) apush(*out, dd_pairs[i]);
	} else {
		int dd_count = asize(dd_pairs);
		for (int i = 0; i < dd_count; i++)
			narrowphase_pair(w, dd_pairs[i].a, dd_pairs[i].b, manifolds);
	}
	afree(dd_pairs);
	bp_frame_count++;

	CK_FREE(tight);
	afree(sap);
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
