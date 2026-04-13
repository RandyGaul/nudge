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

	// SIMD path: 4-wide support scan using SoA vertex data.
	if (hull->soa_verts) {
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

		// Scalar tail for remaining vertices.
		for (int i = n4; i < n; i++) {
			float d = sx[i]*dir.x + sy[i]*dir.y + sz[i]*dir.z;
			if (d > best) { best = d; best_i = i; }
		}
		return hull->verts[best_i];
	}

	// Scalar fallback: linear scan over AoS vertices.
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
	WarmManifold* warm; // cached warm cache pointer (avoids duplicate hash lookup in pre_solve)
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

// Segment helpers (segment_closest_t, segment_closest_point, segments_closest_points)
// now in vmath.h.

// Get capsule segment endpoints in world space.
static void capsule_world_segment(BodyState* bs, ShapeInternal* s, v3* P, v3* Q)
{
	v3 local_p = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 local_q = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	*P = add(bs->position, rotate(bs->rotation, local_p));
	*Q = add(bs->position, rotate(bs->rotation, local_q));
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
	v3 da = sub(a.q, a.p), db = sub(b.q, b.p);
	float la2 = len2(da), lb2 = len2(db);
	float r_sum = a.radius + b.radius;

	// Parallel / near-parallel capsules: generate up to 2 contacts along
	// the overlapping segment range for rotational stability.
	if (la2 > 1e-8f && lb2 > 1e-8f) {
		float la = sqrtf(la2), lb = sqrtf(lb2);
		v3 ua = scale(da, 1.0f / la), ub = scale(db, 1.0f / lb);
		float par = fabsf(dot(ua, ub));
		if (par > 0.99f) {
			// Perpendicular distance between the two infinite lines.
			v3 w = sub(a.p, b.p);
			v3 perp = sub(w, scale(ua, dot(w, ua)));
			float perp_len = len(perp);
			if (perp_len > r_sum) return 0;
			if (!manifold) return 1;
			v3 normal = perp_len > 1e-6f ? scale(perp, -1.0f / perp_len) : V3(0, 1, 0);

			// Project B's endpoints onto A's axis and clip to A's range [0, la].
			float tb_p = dot(sub(b.p, a.p), ua);
			float tb_q = dot(sub(b.q, a.p), ua);
			if (tb_p > tb_q) { float tmp = tb_p; tb_p = tb_q; tb_q = tmp; }
			float t0 = tb_p > 0.0f ? tb_p : 0.0f;
			float t1 = tb_q < la ? tb_q : la;
			if (t0 >= t1 - 1e-6f) {
				// Overlap collapsed to a point — treat as single contact.
				float tm = (t0 + t1) * 0.5f;
				v3 pa = add(a.p, scale(ua, tm));
				v3 pb = segment_closest_point(b.p, b.q, pa);
				v3 d = sub(pb, pa);
				float dist = len(d);
				if (dist > r_sum) return 0;
				v3 n = dist > 1e-6f ? scale(d, 1.0f / dist) : normal;
				manifold->count = 1;
				manifold->contacts[0] = (Contact){ .point = add(pa, scale(n, a.radius)), .normal = n, .penetration = r_sum - dist };
				return 1;
			}

			// Two contacts at overlap endpoints.
			int cp = 0;
			float ts[2] = { t0, t1 };
			for (int i = 0; i < 2; i++) {
				v3 pa = add(a.p, scale(ua, ts[i]));
				v3 pb = segment_closest_point(b.p, b.q, pa);
				v3 d = sub(pb, pa);
				float dist = len(d);
				if (dist > r_sum) continue;
				v3 n = dist > 1e-6f ? scale(d, 1.0f / dist) : normal;
				manifold->contacts[cp++] = (Contact){ .point = add(pa, scale(n, a.radius)), .normal = n, .penetration = r_sum - dist, .feature_id = (uint32_t)i };
			}
			if (cp > 0) { manifold->count = cp; return 1; }
			return 0;
		}
	}

	// General case: single contact at segment-segment closest points.
	v3 c1, c2;
	segments_closest_points(a.p, a.q, b.p, b.q, &c1, &c2);

	v3 d = sub(c2, c1);
	float dist2 = len2(d);

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
		// diff points from box toward sphere (B->A). Negate for A->B convention.
		float dist = sqrtf(dist2);
		v3 local_n = scale(diff, -1.0f / dist);
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

	// face_sign gives outward box normal. Negate for A->B (sphere->box) convention.
	v3 local_n = V3(0, 0, 0);
	static const int face_axis[6] = {2, 2, 0, 0, 1, 1};
	static const float face_sign[6] = {-1, 1, -1, 1, -1, 1};
	(&local_n.x)[face_axis[best_face]] = -face_sign[best_face];
	v3 world_n = rotate(b.rotation, local_n);
	v3 world_pt = sub(a.center, scale(world_n, a.radius));
	manifold->count = 1;
	manifold->contacts[0] = (Contact){ .point = world_pt, .normal = world_n, .penetration = a.radius + best_depth, .feature_id = 0 };
	return 1;
}

// Capsule-hull narrowphase.
// Shallow (GJK dist > LINEAR_SLOP): single reference plane from GJK, project
// both endpoints.  Two contacts for face-parallel stability.
// Deep (GJK dist <= LINEAR_SLOP): per-endpoint face search.  Each endpoint
// independently finds its own minimum-penetration face, so contacts stay valid
// when endpoints are in different Voronoi regions (e.g. near a box edge/corner).
// Both paths place contacts on the hull surface for continuity at the boundary.
int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold)
{
	GJK_Result r = gjk_query_segment_hull(a.p, a.q, b);

	if (r.distance > a.radius) return 0;
	if (!manifold) return 1;

	if (r.distance > LINEAR_SLOP) {
		// --- Shallow path ---
		v3 gjk_n = scale(sub(r.point2, r.point1), 1.0f / r.distance);
		v3 outward = neg(gjk_n);  // hull outward direction at contact

		// Try to match the GJK normal to a hull face.  Only face contacts
		// can produce valid multi-contact manifolds via endpoint projection.
		// Edge/vertex contacts give a diagonal reference plane that makes
		// projected contact points fly off the face -- use single GJK
		// contact for those.
		const Hull* hull = b.hull;
		int face_match = -1;
		HullPlane matched_plane = {0};
		for (int i = 0; i < hull->face_count; i++) {
			HullPlane wp = plane_transform(hull->planes[i], b.center, b.rotation, b.scale);
			if (dot(outward, wp.normal) > 0.99f) {
				face_match = i;
				matched_plane = wp;
				break;
			}
		}

		if (face_match >= 0) {
			// Face contact: clip capsule segment against face side planes,
			// then test clipped endpoints for penetration.  This produces
			// 0-2 contacts with proper rotational stability.
			v3 fn = matched_plane.normal;
			float fd = matched_plane.offset;

			// Clip segment to face polygon (Sutherland-Hodgman on a segment).
			v3 cp_p = a.p, cp_q = a.q;
			int clipped_out = 0;
			{
				int start_e = hull->faces[face_match].edge, ei = start_e;
				do {
					v3 tail = add(b.center, rotate(b.rotation, hull_vert_scaled(hull, hull->edge_origin[ei], b.scale)));
					int twin = hull->edge_twin[ei];
					v3 head = add(b.center, rotate(b.rotation, hull_vert_scaled(hull, hull->edge_origin[twin], b.scale)));
					v3 side_n = cross(sub(head, tail), fn);
					float side_d = dot(side_n, tail);
					float dp_s = dot(side_n, cp_p) - side_d;
					float dq_s = dot(side_n, cp_q) - side_d;
					if (dp_s < 0 && dq_s < 0) { clipped_out = 1; break; }
					if (dp_s * dq_s < 0) {
						float t = dp_s / (dp_s - dq_s);
						if (dp_s < 0) cp_p = add(cp_p, scale(sub(cp_q, cp_p), t));
						else          cp_q = add(cp_p, scale(sub(cp_q, cp_p), t));
					}
					ei = hull->edge_next[ei];
				} while (ei != start_e);
			}

			if (!clipped_out) {
				float dp = dot(cp_p, fn) - fd;
				float dq = dot(cp_q, fn) - fd;
				int ncp = 0;
				if (a.radius - dp >= -LINEAR_SLOP) {
					float pen = a.radius - dp;
					if (pen < 0) pen = 0;
					manifold->contacts[ncp++] = (Contact){
						.point = sub(cp_p, scale(fn, dp)),
						.normal = neg(fn),
						.penetration = pen,
						.feature_id = 0,
					};
				}
				if (len2(sub(cp_p, cp_q)) > 1e-8f && a.radius - dq >= -LINEAR_SLOP) {
					float pen = a.radius - dq;
					if (pen < 0) pen = 0;
					manifold->contacts[ncp++] = (Contact){
						.point = sub(cp_q, scale(fn, dq)),
						.normal = neg(fn),
						.penetration = pen,
						.feature_id = 1,
					};
				}
				if (ncp > 0) {
					manifold->count = ncp;
					return 1;
				}
			}
			// Segment entirely clipped away -- fall through to single GJK contact.
		}

		// Edge/vertex contact (or face fallback): single contact from GJK.
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = add(r.point1, scale(gjk_n, a.radius)),
			.normal = gjk_n,
			.penetration = a.radius - r.distance,
			.feature_id = 2,
		};
		return 1;
	}

	// --- Deep path: per-endpoint face search. ---
	// Each endpoint independently finds the hull face of minimum penetration.
	// This handles the case where endpoints are in different Voronoi regions
	// (one near top face, the other near a side face after tumbling in).
	const Hull* hull = b.hull;
	int cp = 0;
	v3 endpoints[2] = { a.p, a.q };
	for (int ei = 0; ei < 2; ei++) {
		v3 pt = endpoints[ei];
		float best_sep = -1e18f;
		v3 best_n = V3(0, 1, 0);
		float best_d = 0;
		for (int fi = 0; fi < hull->face_count; fi++) {
			HullPlane wp = plane_transform(hull->planes[fi], b.center, b.rotation, b.scale);
			float s = dot(pt, wp.normal) - wp.offset;
			if (s > a.radius) goto next_ep;  // this endpoint is separated on this face
			if (s > best_sep) { best_sep = s; best_n = wp.normal; best_d = wp.offset; }
		}
		{
			float pen = a.radius - best_sep;
			if (pen >= -LINEAR_SLOP) {
				if (pen < 0) pen = 0;
				manifold->contacts[cp++] = (Contact){
					.point = sub(pt, scale(best_n, best_sep)),  // project onto face
					.normal = neg(best_n),
					.penetration = pen,
					.feature_id = (uint32_t)ei,
				};
			}
		}
		next_ep:;
	}

	// Fallback: cylindrical portion contacts hull but neither hemisphere does.
	if (cp == 0) {
		v3 gjk_n = r.distance > 1e-6f ? scale(sub(r.point2, r.point1), 1.0f / r.distance) : V3(0, 1, 0);
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = r.point2,
			.normal = gjk_n,
			.penetration = a.radius - r.distance,
			.feature_id = 2,
		};
		return 1;
	}

	manifold->count = cp;
	return 1;
}

// Capsule-box: route through capsule-hull with the unit box hull.
// The previous analytical approach (iterative segment-box clamping) failed for
// tilted capsules in deep penetration — the 2-iteration clamp didn't converge,
// giving wildly underreported penetration and flickering contacts. The hull path
// (GJK + SAT) handles all orientations robustly.
int collide_capsule_box(Capsule a, Box b, Manifold* manifold)
{
	ConvexHull bh = { &s_unit_box_hull, b.center, b.rotation, b.half_extents };
	return collide_capsule_hull(a, bh, manifold);
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
static int g_sat_hillclimb_enabled = 1; // toggled from WorldInternal before narrowphase

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
	if (g_sat_hillclimb_enabled && face_hint >= 0 && face_hint < hull1->face_count && hull1->face_count > 8) {
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
// SAT: Gauss map Minkowski face test (arc intersection on the unit sphere).
// 4-wide SIMD: returns movemask of lanes where arcs intersect.
#define gauss_map_prune(cba, dba, adc, bdc) simd_movemask(simd_and(simd_and(simd_cmpgt(simd_zero(), simd_mul(cba, dba)), simd_cmpgt(simd_zero(), simd_mul(adc, bdc))), simd_cmpgt(simd_mul(cba, bdc), simd_zero())))

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

	float inv_l = 1.0f / sqrtf(l2);
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

	// Reformulate: dot(n, c1 + rot(rel, v*sc)) = dot(n,c1) + dot(rot(inv_rel,n)*sc, v)
	// Precompute scaled direction in each hull's local space for SoA SIMD scan.
	float bias1 = dot(n, c1);
	quat inv_rel = inv(rel_rot);
	v3 nd1 = rotate(inv_rel, n);
	float d1x = nd1.x * scale1.x, d1y = nd1.y * scale1.y, d1z = nd1.z * scale1.z;
	float d2x = -n.x * scale2.x, d2y = -n.y * scale2.y, d2z = -n.z * scale2.z; // neg for min

	int nv1 = hull1->vert_count, nv2 = hull2->vert_count;
	if (hull1->soa_verts && hull2->soa_verts) {
		int p1 = (nv1 + 3) & ~3, p2 = (nv2 + 3) & ~3;
		const float* s1x = hull1->soa_verts, *s1y = s1x + p1, *s1z = s1y + p1;
		const float* s2x = hull2->soa_verts, *s2y = s2x + p2, *s2z = s2y + p2;
		simd4f vd1x = simd_set1(d1x), vd1y = simd_set1(d1y), vd1z = simd_set1(d1z);
		simd4f max4 = simd_set1(-1e18f);
		for (int i = 0; i < p1; i += 4) {
			simd4f d = simd_add(simd_add(simd_mul(vd1x, simd_load(s1x+i)), simd_mul(vd1y, simd_load(s1y+i))), simd_mul(vd1z, simd_load(s1z+i)));
			max4 = simd_max(max4, d);
		}
		float ds[4]; simd_store(ds, max4);
		float max1 = ds[0]; for (int k = 1; k < 4; k++) if (ds[k] > max1) max1 = ds[k];
		max1 += bias1;

		simd4f vd2x = simd_set1(d2x), vd2y = simd_set1(d2y), vd2z = simd_set1(d2z);
		simd4f max4b = simd_set1(-1e18f);
		for (int i = 0; i < p2; i += 4) {
			simd4f d = simd_add(simd_add(simd_mul(vd2x, simd_load(s2x+i)), simd_mul(vd2y, simd_load(s2y+i))), simd_mul(vd2z, simd_load(s2z+i)));
			max4b = simd_max(max4b, d);
		}
		simd_store(ds, max4b);
		float max_neg2 = ds[0]; for (int k = 1; k < 4; k++) if (ds[k] > max_neg2) max_neg2 = ds[k];
		return -max_neg2 - max1;
	}

	float max1 = -1e18f, min2 = 1e18f;
	for (int i = 0; i < nv1; i++) {
		v3 v = add(c1, rotate(rel_rot, hull_vert_scaled(hull1, i, scale1)));
		float d = dot(n, v);
		if (d > max1) max1 = d;
	}
	for (int i = 0; i < nv2; i++) {
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

	// Transpose hull2 edge normals from AoS to SoA for SIMD Gauss map pruning.
	float nu2x[128], nu2y[128], nu2z[128], nv2x[128], nv2y[128], nv2z[128];
	float ne2x[128], ne2y[128], ne2z[128];
	for (int k = 0; k < n2; k++) {
		nu2x[k] = nu2_arr[k].x; nu2y[k] = nu2_arr[k].y; nu2z[k] = nu2_arr[k].z;
		nv2x[k] = nv2_arr[k].x; nv2y[k] = nv2_arr[k].y; nv2z[k] = nv2_arr[k].z;
		ne2x[k] = ne2_arr[k].x; ne2y[k] = ne2_arr[k].y; ne2z[k] = ne2_arr[k].z;
	}
	for (int k = n2; k < ((n2 + 3) & ~3); k++) {
		nu2x[k] = nu2y[k] = nu2z[k] = nv2x[k] = nv2y[k] = nv2z[k] = ne2x[k] = ne2y[k] = ne2z[k] = 0;
	}

	// For each hull1 edge, test all hull2 edges with 4-wide SIMD Gauss map pruning.
	for (int k1 = 0; k1 < n1; k1++) {
		v3 e1 = e1_arr[k1], u1 = u1_arr[k1], v1 = v1_arr[k1], b_x_a = ne1_arr[k1];

		// Broadcast hull1 edge data for SIMD inner loop.
		simd4f bxa_x = simd_set1(b_x_a.x), bxa_y = simd_set1(b_x_a.y), bxa_z = simd_set1(b_x_a.z);
		simd4f u1x = simd_set1(u1.x), u1y = simd_set1(u1.y), u1z = simd_set1(u1.z);
		simd4f v1x = simd_set1(v1.x), v1y = simd_set1(v1.y), v1z = simd_set1(v1.z);

		// SIMD 4-wide inner loop: Gauss map arc intersection test.
		int n2_4 = n2 & ~3;
		for (int k2 = 0; k2 < n2_4; k2 += 4) {
			simd4f cba = simd_add(simd_add(simd_mul(simd_load(nu2x+k2), bxa_x), simd_mul(simd_load(nu2y+k2), bxa_y)), simd_mul(simd_load(nu2z+k2), bxa_z));
			simd4f dba = simd_add(simd_add(simd_mul(simd_load(nv2x+k2), bxa_x), simd_mul(simd_load(nv2y+k2), bxa_y)), simd_mul(simd_load(nv2z+k2), bxa_z));
			simd4f adc = simd_add(simd_add(simd_mul(u1x, simd_load(ne2x+k2)), simd_mul(u1y, simd_load(ne2y+k2))), simd_mul(u1z, simd_load(ne2z+k2)));
			simd4f bdc = simd_add(simd_add(simd_mul(v1x, simd_load(ne2x+k2)), simd_mul(v1y, simd_load(ne2y+k2))), simd_mul(v1z, simd_load(ne2z+k2)));

			int mask = gauss_map_prune(cba, dba, adc, bdc);
			if (!mask) continue;

			for (int lane = 0; lane < 4; lane++) {
				if (!(mask & (1 << lane))) continue;
				int k2i = k2 + lane;
				float sep = sat_edge_project_full(e1, e2_arr[k2i], c1_local, hull1, rel_rot, scale1, hull2, scale2);
				if (sep > best.separation) { best.index1 = k1*2; best.index2 = k2i*2; best.separation = sep; }
			}
		}

		// Scalar tail for remaining hull2 edges.
		for (int k2 = n2_4; k2 < n2; k2++) {
			v3 ne2 = ne2_arr[k2];
			if (!gauss_map_prune(simd_set1(dot(nu2_arr[k2], b_x_a)), simd_set1(dot(nv2_arr[k2], b_x_a)), simd_set1(dot(u1, ne2)), simd_set1(dot(v1, ne2)))) continue;
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
static SIMD_FORCEINLINE int clip_to_plane(v3* in, uint8_t* in_fid, int in_count, v3 plane_n, float plane_d, uint8_t clip_edge, v3* out, uint8_t* out_fid)
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
static SIMD_FORCEINLINE int reduce_contacts(Contact* contacts, int count)
{
	if (count <= MAX_CONTACTS) return count;

	int sel[MAX_CONTACTS];
	int used[MAX_CLIP_VERTS];
	memset(used, 0, count * sizeof(int));

	// Step 1: two farthest points
	float best_d2 = -1.0f;
	int i0 = 0, i1 = 1;
	for (int i = 0; i < count; i++) {
		for (int j = i + 1; j < count; j++) {
			float d2 = len2(sub(contacts[j].point, contacts[i].point));
			if (d2 > best_d2) { best_d2 = d2; i0 = i; i1 = j; }
		}
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
static SIMD_FORCEINLINE int generate_face_contact(const Hull* ref_hull, v3 ref_pos, quat ref_rot, v3 ref_sc, const Hull* inc_hull, v3 inc_pos, quat inc_rot, v3 inc_sc, int ref_face, int flip, Manifold* manifold)
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
		for (int i = 0; i < clip_count; i++) {
			for (int c = 0; c < ncorners; c++) {
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) { in_fid[i] = 0xC0 | (uint8_t)c; break; }
			}
		}
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

int collide_hull_hull_ex(ConvexHull a, ConvexHull b, Manifold* manifold, int* sat_hint, CachedFeaturePair* out_pair)
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
	float max_face_sep = face_a.separation > face_b.separation ? face_a.separation : face_b.separation;
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
		if (out_pair) *out_pair = (CachedFeaturePair){.type = 2, .edge_a = (int16_t)edge_q.index1, .edge_b = (int16_t)edge_q.index2};
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
	for (int i = 0; i < cp; i++) {
		manifold->contacts[i] = tmp_contacts[i];
	}
	// Cache the winning axis for next-frame warm-start.
	if (sat_hint) *sat_hint = flip ? hull_a->face_count + ref_face : ref_face;
	if (out_pair) *out_pair = (CachedFeaturePair){.type = 1, .ref_body = (int16_t)flip, .face_a = (int16_t)(flip ? inc_face : ref_face), .face_b = (int16_t)(flip ? ref_face : inc_face)};
	return 1;
}

int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold)
{
	return collide_hull_hull_ex(a, b, manifold, NULL, NULL);
}

// -----------------------------------------------------------------------------
// Dedicated box-box SAT: Gottschalk OBB test with direct rotation column arithmetic.
// Tests all 15 separating axes without generic hull machinery (plane_transform, hull_support loops, Gauss map).
// Falls through to collide_hull_hull for contact generation when penetrating.
// sat_hint: if non-NULL, *sat_hint is the cached axis from last frame (-1 = no cache).
// On return, *sat_hint is updated to the winning axis for next frame.
static int collide_box_box_ex(Box a, Box b, Manifold* manifold, int* sat_hint, CachedFeaturePair* out_pair)
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
		// nsign selects which face of the ref box is the contact face.
		// For A-ref (flip=0): positive proj means B is in + direction, A's + face faces B → nsign=+1.
		// For B-ref (flip=1): positive proj means B is in + direction of its own axis,
		// so B's - face faces A → nsign must be -1. Negate for flip.
		float nsign = proj_sign >= 0.0f ? 1.0f : -1.0f;
		if (flip) nsign = -nsign;
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
		for (int i = 0; i < clip_count; i++) {
			for (int c = 0; c < 4; c++) {
				if (len2(sub(in_buf[i], corners[c])) < snap_tol2) { in_fid[i] = 0xC0 | (uint8_t)c; break; }
			}
		}

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
		if (out_pair) *out_pair = (CachedFeaturePair){.type = 1, .ref_body = (int16_t)flip, .face_a = (int16_t)(flip ? inc_face : ref_face), .face_b = (int16_t)(flip ? ref_face : inc_face)};
		return 1;
	}

	// Edge-edge contact: fall through to hull-hull (edge contacts are rare in box piles).
	if (sat_hint) *sat_hint = best_axis;
	return collide_hull_hull_ex(
		(ConvexHull){ &s_unit_box_hull, a.center, a.rotation, a.half_extents },
		(ConvexHull){ &s_unit_box_hull, b.center, b.rotation, b.half_extents },
		manifold, NULL, out_pair);
}

int collide_box_box(Box a, Box b, Manifold* manifold) { return collide_box_box_ex(a, b, manifold, NULL, NULL); }

const Hull* hull_unit_box() { return &s_unit_box_hull; }

// -----------------------------------------------------------------------------
// Cylinder Voronoi region classification.
//
// Used by native cylinder narrowphase routines. Takes a world-space point
// (typically the witness on the other shape from GJK on the cylinder's axis
// segment) and returns the true nearest point on the cylinder surface, plus
// its outward normal and signed distance.
//
// The four regions partition space outside the cylinder into distinct feature
// owners, plus INSIDE for the interior:
//   SIDE:   |axial| <= hh && rad >= r   -- curved wall owns the witness
//   CAP:    |axial| >  hh && rad <= r   -- flat top/bottom owns the witness
//   RIM:    |axial| >  hh && rad >  r   -- the circular edge owns the witness
//   INSIDE: |axial| <  hh && rad <  r   -- escape via nearest side or cap face
//
// Distance is positive outside, negative inside (value = -escape_depth).

typedef enum
{
	CYL_REGION_SIDE,
	CYL_REGION_CAP,
	CYL_REGION_RIM,
	CYL_REGION_INSIDE,
} CylRegion;

typedef struct CylFeature
{
	CylRegion region;
	float distance; // signed: >0 outside, <0 inside (-escape_depth)
	v3 surface_pt;  // world-space point on cylinder surface
	v3 normal;      // unit outward cylinder normal at surface_pt (world)
} CylFeature;

static CylFeature cyl_classify_point(v3 x_world, v3 cyl_pos, quat cyl_rot, float half_height, float radius)
{
	// Transform witness into cylinder-local space.
	v3 lp = rotate(inv(cyl_rot), sub(x_world, cyl_pos));
	float rad2 = lp.x * lp.x + lp.z * lp.z;
	float rad = sqrtf(rad2);
	float axial = lp.y;
	float abs_axial = fabsf(axial);
	float sign_y = axial >= 0.0f ? 1.0f : -1.0f;

	CylFeature feat;
	v3 sp_local, n_local;

	if (abs_axial <= half_height) {
		if (rad >= radius) {
			// SIDE: perpendicular distance to the curved wall.
			feat.region = CYL_REGION_SIDE;
			feat.distance = rad - radius;
			// rad is >= radius > 0 so inv_rad is finite.
			float inv_rad = 1.0f / rad;
			float rx = lp.x * inv_rad, rz = lp.z * inv_rad;
			sp_local = V3(rx * radius, axial, rz * radius);
			n_local = V3(rx, 0.0f, rz);
		} else {
			// INSIDE: pick the nearest escape face (side vs cap).
			feat.region = CYL_REGION_INSIDE;
			float side_esc = radius - rad;
			float cap_esc = half_height - abs_axial;
			if (cap_esc < side_esc) {
				feat.distance = -cap_esc;
				sp_local = V3(lp.x, sign_y * half_height, lp.z);
				n_local = V3(0.0f, sign_y, 0.0f);
			} else {
				feat.distance = -side_esc;
				if (rad > 1e-12f) {
					float inv_rad = 1.0f / rad;
					float rx = lp.x * inv_rad, rz = lp.z * inv_rad;
					sp_local = V3(rx * radius, axial, rz * radius);
					n_local = V3(rx, 0.0f, rz);
				} else {
					// Exactly on the axis -- arbitrary radial (+X).
					sp_local = V3(radius, axial, 0.0f);
					n_local = V3(1.0f, 0.0f, 0.0f);
				}
			}
		}
	} else {
		if (rad <= radius) {
			// CAP: perpendicular distance to the flat disk.
			feat.region = CYL_REGION_CAP;
			feat.distance = abs_axial - half_height;
			sp_local = V3(lp.x, sign_y * half_height, lp.z);
			n_local = V3(0.0f, sign_y, 0.0f);
		} else {
			// RIM: closest feature is the circular edge at (rad=r, axial=+/-hh).
			feat.region = CYL_REGION_RIM;
			float dr = rad - radius;
			float da = abs_axial - half_height;
			feat.distance = sqrtf(dr * dr + da * da);
			float inv_rad = 1.0f / rad; // rad > radius > 0
			float rx = lp.x * inv_rad, rz = lp.z * inv_rad;
			sp_local = V3(rx * radius, sign_y * half_height, rz * radius);
			// Outward normal along the direction from the rim point to the witness.
			float inv_d = feat.distance > 1e-12f ? 1.0f / feat.distance : 0.0f;
			n_local = V3(rx * dr * inv_d, sign_y * da * inv_d, rz * dr * inv_d);
		}
	}

	feat.surface_pt = add(cyl_pos, rotate(cyl_rot, sp_local));
	feat.normal = rotate(cyl_rot, n_local);
	return feat;
}

// hull_build_csr, compact hull converters, face extension, hull_from_compact
// all moved to quickhull.c -- live with hull construction code.

// -----------------------------------------------------------------------------
// Incremental narrowphase: validate + refresh cached feature pairs.

// Box face validation: check if the cached local axis is still the dominant
// projection axis for the separation direction. 3 dot products + 3 fabsf.
static int validate_cached_face_box(v3 cols[3], v3 d, int cached_la)
{
	float dots[3] = { fabsf(dot(d, cols[0])), fabsf(dot(d, cols[1])), fabsf(dot(d, cols[2])) };
	for (int i = 0; i < 3; i++) {
		if (i != cached_la && dots[i] > dots[cached_la])
			return 0;
	}
	return 1;
}

// Hull face validation: check if any neighbor face of the cached face is more
// aligned with the separation direction. Uses half-edge neighbor traversal.
static int validate_cached_face_hull(const Hull* hull, v3 pos, quat rot, v3 sc, int cached_face, v3 sep_dir)
{
	HullPlane cached_plane = plane_transform(hull->planes[cached_face], pos, rot, sc);
	float cached_dot = dot(cached_plane.normal, sep_dir);
	int ei = hull->faces[cached_face].edge;
	int start = ei;
	do {
		int adj = hull->edge_face[hull->edge_twin[ei]];
		HullPlane adj_plane = plane_transform(hull->planes[adj], pos, rot, sc);
		if (dot(adj_plane.normal, sep_dir) > cached_dot + 1e-4f) return 0;
		ei = hull->edge_next[ei];
	} while (ei != start);
	return 1;
}

// Refresh box-box face contact using cached feature pair (skip full SAT).
// Returns 1 if refresh succeeded, 0 if invalidated (caller falls through to full SAT).
static int refresh_box_box_face(Box a, Box b, Manifold* manifold, CachedFeaturePair* cp)
{
	v3 ax = rotate(a.rotation, V3(1, 0, 0)), ay = rotate(a.rotation, V3(0, 1, 0)), az = rotate(a.rotation, V3(0, 0, 1));
	v3 bx = rotate(b.rotation, V3(1, 0, 0)), by = rotate(b.rotation, V3(0, 1, 0)), bz = rotate(b.rotation, V3(0, 0, 1));
	v3 d = sub(b.center, a.center);

	// face_map: local axis + sign -> hull face index. Matches collide_box_box_ex.
	// face_map_la[face] = local axis, face_map_sign[face] = +1 or -1.
	static const int face_map_la[6] = {2, 2, 0, 0, 1, 1};   // -Z=0,+Z=1,-X=2,+X=3,-Y=4,+Y=5
	static const float face_map_sign[6] = {-1, 1, -1, 1, -1, 1};

	v3 a_cols[3] = {ax, ay, az}, b_cols[3] = {bx, by, bz};
	v3 ref_cols[3], inc_cols[3], ref_pos, inc_pos, ref_he, inc_he;
	int ref_fi, inc_fi;
	int flip = cp->ref_body;
	if (!flip) { ref_fi = cp->face_a; inc_fi = cp->face_b; for (int i = 0; i < 3; i++) { ref_cols[i] = a_cols[i]; inc_cols[i] = b_cols[i]; } ref_pos = a.center; inc_pos = b.center; ref_he = a.half_extents; inc_he = b.half_extents; }
	else { ref_fi = cp->face_b; inc_fi = cp->face_a; for (int i = 0; i < 3; i++) { ref_cols[i] = b_cols[i]; inc_cols[i] = a_cols[i]; } ref_pos = b.center; inc_pos = a.center; ref_he = b.half_extents; inc_he = a.half_extents; }

	int la = face_map_la[ref_fi];
	float nsign = face_map_sign[ref_fi];
	int inc_la = face_map_la[inc_fi];

	// Validate reference side: is cached axis still dominant?
	v3 sep = flip ? neg(d) : d;
	if (!validate_cached_face_box(ref_cols, sep, la)) return 0;
	v3 ref_n = scale(ref_cols[la], nsign);

	// Recompute incident face from current transforms (3 dots, always correct).
	inc_la = 0; float inc_best = 1e18f;
	for (int i = 0; i < 3; i++) { float dd = dot(inc_cols[i], ref_n); if (dd < inc_best) { inc_best = dd; inc_la = i; } if (-dd < inc_best) { inc_best = -dd; inc_la = i; } }
	float inc_nsign = dot(inc_cols[inc_la], ref_n) > 0 ? -1.0f : 1.0f;
	static const int face_map[3][2] = { {2, 3}, {4, 5}, {0, 1} };
	inc_fi = face_map[inc_la][inc_nsign > 0 ? 1 : 0];

	// Re-run face clip with known ref / recomputed inc faces.
	float ref_off = dot(ref_n, ref_pos) + (&ref_he.x)[la];
	int u = (la + 1) % 3, v_ax = (la + 2) % 3;
	int inc_u = (inc_la + 1) % 3, inc_v = (inc_la + 2) % 3;

	v3 inc_center = add(inc_pos, scale(inc_cols[inc_la], inc_nsign * (&inc_he.x)[inc_la]));
	v3 inc_eu = scale(inc_cols[inc_u], (&inc_he.x)[inc_u]), inc_ev = scale(inc_cols[inc_v], (&inc_he.x)[inc_v]);
	v3 buf1[MAX_CLIP_VERTS], buf2[MAX_CLIP_VERTS];
	uint8_t fid1[MAX_CLIP_VERTS], fid2[MAX_CLIP_VERTS];
	buf1[0] = add(inc_center, add(inc_eu, inc_ev));
	buf1[1] = add(inc_center, sub(neg(inc_eu), neg(inc_ev)));
	buf1[2] = add(inc_center, sub(neg(inc_eu), inc_ev));
	buf1[3] = add(inc_center, sub(inc_eu, inc_ev));
	if (dot(cross(sub(buf1[1], buf1[0]), sub(buf1[2], buf1[0])), ref_n) > 0) { v3 tmp = buf1[1]; buf1[1] = buf1[3]; buf1[3] = tmp; }
	for (int i = 0; i < 4; i++) fid1[i] = 0x80 | (uint8_t)i;
	int clip_count = 4;

	v3 side_n[4] = { ref_cols[u], neg(ref_cols[u]), ref_cols[v_ax], neg(ref_cols[v_ax]) };
	float side_d[4] = { dot(ref_cols[u], ref_pos) + (&ref_he.x)[u], dot(neg(ref_cols[u]), ref_pos) + (&ref_he.x)[u], dot(ref_cols[v_ax], ref_pos) + (&ref_he.x)[v_ax], dot(neg(ref_cols[v_ax]), ref_pos) + (&ref_he.x)[v_ax] };

	v3* in_buf = buf1; v3* out_buf = buf2;
	uint8_t* in_fid = fid1; uint8_t* out_fid = fid2;
	for (int i = 0; i < 4; i++) {
		clip_count = clip_to_plane(in_buf, in_fid, clip_count, side_n[i], side_d[i], (uint8_t)i, out_buf, out_fid);
		v3* sw = in_buf; in_buf = out_buf; out_buf = sw;
		uint8_t* fs = in_fid; in_fid = out_fid; out_fid = fs;
	}

	v3 ref_center = add(ref_pos, scale(ref_cols[la], nsign * (&ref_he.x)[la]));
	v3 ref_eu = scale(ref_cols[u], (&ref_he.x)[u]), ref_ev = scale(ref_cols[v_ax], (&ref_he.x)[v_ax]);
	v3 corners[4] = { add(ref_center, add(ref_eu, ref_ev)), add(ref_center, sub(ref_eu, ref_ev)), add(ref_center, sub(neg(ref_eu), neg(ref_ev))), add(ref_center, sub(neg(ref_eu), ref_ev)) };
	float snap_tol2 = 1e-6f;
	for (int i = 0; i < clip_count; i++) {
		for (int c = 0; c < 4; c++) {
			if (len2(sub(in_buf[i], corners[c])) < snap_tol2) { in_fid[i] = 0xC0 | (uint8_t)c; break; }
		}
	}

	v3 contact_n = flip ? neg(ref_n) : ref_n;
	Contact tmp_contacts[MAX_CLIP_VERTS];
	int cpc = 0;
	for (int i = 0; i < clip_count; i++) {
		float depth = ref_off - dot(ref_n, in_buf[i]);
		if (depth >= -LINEAR_SLOP) {
			uint32_t fid = flip ? ((uint32_t)inc_fi | ((uint32_t)ref_fi << 8) | ((uint32_t)in_fid[i] << 16)) : ((uint32_t)ref_fi | ((uint32_t)inc_fi << 8) | ((uint32_t)in_fid[i] << 16));
			tmp_contacts[cpc++] = (Contact){ .point = in_buf[i], .normal = contact_n, .penetration = depth, .feature_id = fid };
		}
	}
	if (cpc == 0) return 0;
	cpc = reduce_contacts(tmp_contacts, cpc);
	manifold->count = cpc;
	for (int i = 0; i < cpc; i++) manifold->contacts[i] = tmp_contacts[i];
	return 1;
}

// Refresh hull-hull face contact using cached feature pair.
// Validates both ref and incident faces via neighbor check, then re-clips.
static int refresh_hull_hull_face(ConvexHull a, ConvexHull b, Manifold* manifold, CachedFeaturePair* cp)
{
	const Hull* hull_a = a.hull; v3 pos_a = a.center; quat rot_a = a.rotation; v3 sc_a = a.scale;
	const Hull* hull_b = b.hull; v3 pos_b = b.center; quat rot_b = b.rotation; v3 sc_b = b.scale;

	const Hull* ref_hull; const Hull* inc_hull;
	v3 ref_pos, inc_pos, ref_sc, inc_sc;
	quat ref_rot, inc_rot;
	int ref_face, flip;
	if (cp->ref_body == 0) {
		ref_hull = hull_a; ref_pos = pos_a; ref_rot = rot_a; ref_sc = sc_a;
		inc_hull = hull_b; inc_pos = pos_b; inc_rot = rot_b; inc_sc = sc_b;
		ref_face = cp->face_a; flip = 0;
	} else {
		ref_hull = hull_b; ref_pos = pos_b; ref_rot = rot_b; ref_sc = sc_b;
		inc_hull = hull_a; inc_pos = pos_a; inc_rot = rot_a; inc_sc = sc_a;
		ref_face = cp->face_b; flip = 1;
	}

	// Validate ref face: is it still the best separating face on its side?
	v3 sep_dir = norm(sub(inc_pos, ref_pos));
	if (!validate_cached_face_hull(ref_hull, ref_pos, ref_rot, ref_sc, ref_face, sep_dir)) return 0;

	// Re-clip using existing generate_face_contact (re-discovers incident face internally).
	return generate_face_contact(ref_hull, ref_pos, ref_rot, ref_sc, inc_hull, inc_pos, inc_rot, inc_sc, ref_face, flip, manifold);
}

// -----------------------------------------------------------------------------
// N^2 broadphase + narrowphase dispatch.

// Build a Sphere/Capsule/Box from internal body+shape for broadphase dispatch.
static Sphere make_sphere(BodyState* bs, ShapeInternal* s)
{
	return (Sphere){ add(bs->position, rotate(bs->rotation, s->local_pos)), s->sphere.radius };
}

static Capsule make_capsule(BodyState* bs, ShapeInternal* s)
{
	v3 lp = add(s->local_pos, V3(0, -s->capsule.half_height, 0));
	v3 lq = add(s->local_pos, V3(0,  s->capsule.half_height, 0));
	return (Capsule){ add(bs->position, rotate(bs->rotation, lp)), add(bs->position, rotate(bs->rotation, lq)), s->capsule.radius };
}

static Box make_box(BodyState* bs, ShapeInternal* s)
{
	return (Box){ bs->position, bs->rotation, s->box.half_extents };
}

static ConvexHull make_convex_hull(BodyState* bs, ShapeInternal* s)
{
	return (ConvexHull){ s->hull.hull, bs->position, bs->rotation, s->hull.scale };
}

// Build a Cylinder from an internal body+shape (mirror of make_box/make_sphere).
static Cylinder make_cylinder(BodyState* bs, ShapeInternal* s)
{
	return (Cylinder){ bs->position, bs->rotation, s->cylinder.half_height, s->cylinder.radius };
}

// -----------------------------------------------------------------------------
// Cylinder narrowphase.
//
// Analytical implementations using cyl_classify_point for surface contacts.

// Analytical cyl-sphere: classify the sphere center as a point against the
// cylinder, compare against sphere radius. No GJK iteration -- one call to
// cyl_classify_point handles all four Voronoi regions.
int collide_cylinder_sphere(Cylinder a, Sphere b, Manifold* manifold)
{
	CylFeature feat = cyl_classify_point(b.center, a.center, a.rotation, a.half_height, a.radius);

	float gap = feat.distance - b.radius;
	if (gap > 0.0f) return 0;
	if (!manifold) return 1;

	// Normal points from cyl surface toward sphere center.
	// External regions (SIDE/CAP/RIM): feat.normal already points outward and
	// the sphere center lies in that direction.
	// INSIDE: the sphere center lies OPPOSITE to feat.normal (inside the cyl),
	// so flip to get the A->B direction that the solver will use to push them
	// apart through the nearest escape face.
	v3 normal = feat.region == CYL_REGION_INSIDE ? neg(feat.normal) : feat.normal;
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = feat.surface_pt,
		.normal = normal,
		.penetration = b.radius - feat.distance,
		.feature_id = 0,
	};
	return 1;
}

// Cyl axis in world space as a segment (helper used by several pairs).
static void cylinder_axis_segment(Cylinder c, v3* p, v3* q)
{
	v3 axis = rotate(c.rotation, V3(0, 1, 0));
	*p = sub(c.center, scale(axis, c.half_height));
	*q = add(c.center, scale(axis, c.half_height));
}

// Analytical cyl-capsule: closest-points between cyl axis segment and capsule segment,
// classify the capsule-side witness, inflate by capsule radius. Detects the parallel
// SIDE-SIDE / CAP-CAP case by classifying both capsule endpoints -- when both land in
// the same region with matching normals, emit a 2-point manifold.
int collide_cylinder_capsule(Cylinder a, Capsule b, Manifold* manifold)
{
	v3 cyl_p, cyl_q;
	cylinder_axis_segment(a, &cyl_p, &cyl_q);

	v3 cpa, cpb;
	segments_closest_points(cyl_p, cyl_q, b.p, b.q, &cpa, &cpb);

	CylFeature feat = cyl_classify_point(cpb, a.center, a.rotation, a.half_height, a.radius);
	float gap = feat.distance - b.radius;
	if (gap > 0.0f) return 0;
	if (!manifold) return 1;

	// Check for parallel 2-point manifold: both capsule endpoints land in the same
	// SIDE or CAP region with matching normals (dot > 0.99 ~ 8 degrees).
	CylFeature feat_p = cyl_classify_point(b.p, a.center, a.rotation, a.half_height, a.radius);
	CylFeature feat_q = cyl_classify_point(b.q, a.center, a.rotation, a.half_height, a.radius);
	int parallel =
		feat_p.region == feat_q.region
		&& (feat_p.region == CYL_REGION_SIDE || feat_p.region == CYL_REGION_CAP)
		&& feat_p.distance - b.radius <= 0.0f
		&& feat_q.distance - b.radius <= 0.0f
		&& dot(feat_p.normal, feat_q.normal) > 0.99f;

	if (parallel) {
		manifold->count = 2;
		manifold->contacts[0] = (Contact){
			.point = feat_p.surface_pt,
			.normal = feat_p.normal,
			.penetration = b.radius - feat_p.distance,
		};
		manifold->contacts[1] = (Contact){
			.point = feat_q.surface_pt,
			.normal = feat_q.normal,
			.penetration = b.radius - feat_q.distance,
		};
		return 1;
	}

	// Single contact at the closest-approach witness.
	v3 normal = feat.region == CYL_REGION_INSIDE ? neg(feat.normal) : feat.normal;
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = feat.surface_pt,
		.normal = normal,
		.penetration = b.radius - feat.distance,
	};
	return 1;
}

// Cyl-box: delegate to cyl-hull using the unit box hull. The unit box has only
// 8 verts / 6 faces / 24 edges, so the SAT is cheap.
int collide_cylinder_box(Cylinder a, Box b, Manifold* manifold)
{
	ConvexHull cb = { &s_unit_box_hull, b.center, b.rotation, b.half_extents };
	return collide_cylinder_hull(a, cb, manifold);
}

// Cylinder support projection onto an arbitrary axis (tight, analytical).
// Returns the half-extent: max of dot(surface_point, n) over all cylinder surface points, minus center.
static float cyl_half_support(v3 cyl_axis, float hh, float radius, v3 n)
{
	float an = dot(cyl_axis, n);
	float axial = fabsf(an) * hh;
	float sin2 = 1.0f - an * an;
	float radial = sin2 > 0.0f ? sqrtf(sin2) * radius : 0.0f;
	return axial + radial;
}

// Cyl-hull: GJK(axis_segment, hull) shallow + SAT deep.
// Deep SAT tests hull face normals, cylinder axis, and cyl_axis x hull_edge.
int collide_cylinder_hull(Cylinder a, ConvexHull b, Manifold* manifold)
{
	v3 cyl_p, cyl_q;
	cylinder_axis_segment(a, &cyl_p, &cyl_q);
	v3 cyl_axis = rotate(a.rotation, V3(0, 1, 0));
	const Hull* hull = b.hull;

	// Quick bounding-sphere rejection before running full SAT.
	{
		GJK_Result r = gjk_query_segment_hull(cyl_p, cyl_q, b);
		float bound = sqrtf(a.radius * a.radius + a.half_height * a.half_height);
		if (r.distance > bound) return 0;
	}

	// --- Deep SAT ---
	// Track best separating axis across three families.
	float best_sep = -1e18f;
	v3 best_n = V3(0,0,0);
	int best_type = -1; // 0 = hull face, 1 = cap, 2 = edge

	// Family 1: hull face normals.
	HullPlane best_face_plane = {0};
	for (int i = 0; i < hull->face_count; i++) {
		HullPlane wp = plane_transform(hull->planes[i], b.center, b.rotation, b.scale);
		float cyl_min = dot(a.center, wp.normal) - cyl_half_support(cyl_axis, a.half_height, a.radius, wp.normal);
		float sep = cyl_min - wp.offset;
		if (sep > 0.0f) return 0;
		if (sep > best_sep) { best_sep = sep; best_n = wp.normal; best_type = 0; best_face_plane = wp; }
	}

	// Family 2: cylinder axis (cap face contacts).
	{
		float hull_min = 1e18f, hull_max = -1e18f;
		for (int j = 0; j < hull->vert_count; j++) {
			v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
			float p = dot(wv, cyl_axis);
			if (p < hull_min) hull_min = p;
			if (p > hull_max) hull_max = p;
		}
		float cyl_min_ax = dot(a.center, cyl_axis) - a.half_height;
		float cyl_max_ax = dot(a.center, cyl_axis) + a.half_height;
		if (hull_max < cyl_min_ax || cyl_max_ax < hull_min) return 0;
		// Two overlap measures: hull poking through bottom cap, hull poking through top cap.
		float pen_bot = hull_max - cyl_min_ax; // hull above cyl bottom
		float pen_top = cyl_max_ax - hull_min; // cyl above hull bottom
		float pen = pen_bot < pen_top ? pen_bot : pen_top;
		float sep = -pen;
		int sign = pen_bot < pen_top ? -1 : 1; // -1 = bottom cap reference, +1 = top cap
		if (sep > best_sep + 0.001f) { best_sep = sep; best_n = scale(cyl_axis, (float)sign); best_type = 1; }
	}

	// Family 3: cyl_axis x hull_edge (edge-edge contacts).
	int best_edge = -1;
	for (int i = 0; i < hull->edge_count; i += 2) {
		v3 ev0 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[i]], b.scale)));
		v3 ev1 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[hull->edge_next[i]]], b.scale)));
		v3 ed = sub(ev1, ev0);
		v3 ax = cross(cyl_axis, ed);
		float al = len2(ax);
		if (al < 1e-12f) continue;
		ax = scale(ax, 1.0f / sqrtf(al));
		if (dot(ax, sub(a.center, ev0)) < 0.0f) ax = neg(ax);
		// cyl_axis · ax ≈ 0 (ax is perpendicular to cyl_axis by construction of cross product)
		// so cylinder projects as: center · ax ± radius
		float cyl_min_e = dot(a.center, ax) - a.radius;
		float hull_max_e = -1e18f;
		for (int j = 0; j < hull->vert_count; j++) {
			v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
			float d = dot(wv, ax);
			if (d > hull_max_e) hull_max_e = d;
		}
		float sep = cyl_min_e - hull_max_e;
		if (sep > 0.0f) return 0;
		if (sep > best_sep + 0.001f) { best_sep = sep; best_n = ax; best_type = 2; best_edge = i; }
	}

	if (best_sep > 0.0f || best_type < 0) return 0;
	if (!manifold) return 1;

	// --- Manifold generation by axis type ---

	if (best_type == 0) {
		// Hull face is reference. Two sub-paths depending on the angle between
		// the face normal and the cylinder axis:
		//
		// CAP path (|dot| > 0.7): face is roughly parallel to a cylinder cap.
		// A single rim-point contact would cause pivoting/instability. Instead,
		// project hull verts onto the cap plane and clip against the cap disk,
		// giving up to MAX_CONTACTS multi-point support (same as Family 2).
		//
		// SIDE path (|dot| < 0.7): face is roughly perpendicular to the axis.
		// Use the capsule-style approach: offset axis endpoints by full radius
		// toward the face, giving a 2-point contact strip.
		v3 neg_n = neg(best_n);
		float axis_alignment = fabsf(dot(best_n, cyl_axis));

		if (axis_alignment > 0.7f) {
			// CAP path: pick the cap whose outward normal faces the hull face.
			// cap_n points outward from the chosen cap (toward the hull).
			v3 cap_n = dot(best_n, cyl_axis) > 0.0f ? neg(cyl_axis) : cyl_axis;
			v3 cap_center = add(a.center, scale(cap_n, a.half_height));
			float cap_d = dot(cap_center, cap_n); // cap plane offset

			// Try projecting hull verts onto cap plane, keep those inside the disk.
			int cp = 0;
			v3 points[MAX_CONTACTS]; float depths[MAX_CONTACTS];
			for (int j = 0; j < hull->vert_count && cp < MAX_CONTACTS; j++) {
				v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
				float d = dot(wv, cap_n) - cap_d;
				if (d > 0.0f) continue; // outside cap plane
				v3 on_cap = sub(wv, scale(cap_n, d));
				v3 radial = sub(on_cap, cap_center);
				if (len2(radial) > a.radius * a.radius) continue;
				points[cp] = wv; depths[cp] = -d; cp++;
			}
			// Fallback for large hulls (no verts inside cap disk): project cap
			// center onto the hull face and use that as the contact point.
			if (cp == 0) {
				float face_d = dot(cap_center, best_n) - best_face_plane.offset;
				if (face_d < 0.0f) {
					// Cap center is behind the face — contact at cap center projected onto face.
					points[0] = sub(cap_center, scale(best_n, face_d));
					depths[0] = -face_d;
				} else {
					// Cap center is above the face — use cap center directly.
					points[0] = cap_center;
					depths[0] = -best_sep;
				}
				cp = 1;
			}
			manifold->count = cp;
			for (int i = 0; i < cp; i++)
				manifold->contacts[i] = (Contact){ .point = points[i], .normal = neg_n, .penetration = depths[i] };
			return 1;
		}

		// SIDE path: project axis endpoints offset by full radius toward face.
		float axial_proj = dot(neg_n, cyl_axis);
		v3 perp = sub(neg_n, scale(cyl_axis, axial_proj));
		float perp_len2 = len2(perp);
		v3 surface_offset;
		if (perp_len2 > 1e-8f) surface_offset = scale(perp, a.radius / sqrtf(perp_len2));
		else surface_offset = V3(0, 0, 0);

		float plane_d = best_face_plane.offset;
		v3 sp = add(cyl_p, surface_offset);
		v3 sq = add(cyl_q, surface_offset);
		float dp = dot(sp, best_n);
		float dq = dot(sq, best_n);
		int cp = 0;
		v3 points[2]; float depths[2];
		if (dp < plane_d) { points[cp] = sp; depths[cp] = plane_d - dp; cp++; }
		if (dq < plane_d) { points[cp] = sq; depths[cp] = plane_d - dq; cp++; }
		if (cp == 0) { cp = 1; points[0] = dp < dq ? sp : sq; depths[0] = -best_sep; }
		manifold->count = cp;
		for (int i = 0; i < cp; i++)
			manifold->contacts[i] = (Contact){ .point = points[i], .normal = neg_n, .penetration = depths[i] };
		return 1;
	}

	if (best_type == 1) {
		// Cylinder cap is reference. Project hull vertices onto cap plane,
		// keep up to 4 that are within cap disk radius and penetrate.
		v3 cap_n = best_n; // outward cap normal (+/- cyl_axis)
		float cap_d = dot(a.center, cap_n) + a.half_height; // cap plane offset
		int cp = 0;
		v3 points[MAX_CONTACTS]; float depths[MAX_CONTACTS];
		for (int j = 0; j < hull->vert_count && cp < MAX_CONTACTS; j++) {
			v3 wv = add(b.center, rotate(b.rotation, hmul(hull->verts[j], b.scale)));
			float d = dot(wv, cap_n) - cap_d;
			if (d > 0.0f) continue; // above cap plane (outside)
			// Check if vertex is within cap disk radius.
			v3 on_cap = sub(wv, scale(cap_n, d)); // project onto cap plane
			v3 radial = sub(on_cap, add(a.center, scale(cap_n, a.half_height)));
			if (len2(radial) > a.radius * a.radius) continue;
			points[cp] = wv;
			depths[cp] = -d;
			cp++;
		}
		if (cp == 0) { manifold->count = 1; manifold->contacts[0] = (Contact){ .point = add(a.center, scale(best_n, a.half_height)), .normal = neg(best_n), .penetration = -best_sep }; return 1; }
		manifold->count = cp;
		for (int i = 0; i < cp; i++)
			manifold->contacts[i] = (Contact){ .point = points[i], .normal = neg(best_n), .penetration = depths[i] };
		return 1;
	}

	// best_type == 2: edge-edge contact (single point).
	// Use closest-points between cyl axis and the winning hull edge.
	{
		manifold->count = 1;
		// Recompute closest pair for the winning edge.
		int ei = best_edge;
		v3 ev0 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[ei]], b.scale)));
		v3 ev1 = add(b.center, rotate(b.rotation, hmul(hull->verts[hull->edge_origin[hull->edge_next[ei]]], b.scale)));
		v3 cpa, cpb;
		segments_closest_points(cyl_p, cyl_q, ev0, ev1, &cpa, &cpb);
		CylFeature ef = cyl_classify_point(cpb, a.center, a.rotation, a.half_height, a.radius);
		manifold->contacts[0] = (Contact){
			.point = ef.surface_pt,
			.normal = neg(best_n),
			.penetration = -best_sep,
		};
		return 1;
	}
}

// Cyl-cyl: SAT with both cylinder axes + their cross product.
// Each cylinder projects analytically via cyl_half_support.
int collide_cylinder_cylinder(Cylinder a, Cylinder b, Manifold* manifold)
{
	v3 axis_a = rotate(a.rotation, V3(0, 1, 0));
	v3 axis_b = rotate(b.rotation, V3(0, 1, 0));
	v3 a_p, a_q, b_p, b_q;
	cylinder_axis_segment(a, &a_p, &a_q);
	cylinder_axis_segment(b, &b_p, &b_q);

	// Quick bounding-sphere rejection.
	float bound_a = sqrtf(a.radius * a.radius + a.half_height * a.half_height);
	float bound_b = sqrtf(b.radius * b.radius + b.half_height * b.half_height);
	float dist2_centers = len2(sub(a.center, b.center));
	float bound_sum = bound_a + bound_b;
	if (dist2_centers > bound_sum * bound_sum) return 0;

	float best_sep = -1e18f;
	v3 best_n = V3(0,0,0);
	int best_type = -1; // 0 = axis_a, 1 = axis_b, 2 = cross
	int best_sign = 1;

	// Helper: test a candidate SAT axis. Returns 1 if separating (no collision).
	#define CYL_CYL_TEST_AXIS(n_arg, type_arg, sign_arg) do { \
		v3 n_ = (n_arg); \
		float a_half = cyl_half_support(axis_a, a.half_height, a.radius, n_); \
		float b_half = cyl_half_support(axis_b, b.half_height, b.radius, n_); \
		float a_center = dot(a.center, n_); \
		float b_center = dot(b.center, n_); \
		float a_min = a_center - a_half, a_max = a_center + a_half; \
		float b_min = b_center - b_half, b_max = b_center + b_half; \
		if (a_max < b_min || b_max < a_min) return 0; /* separating axis */ \
		float pen1 = a_max - b_min, pen2 = b_max - a_min; \
		float pen = pen1 < pen2 ? pen1 : pen2; \
		float sep = -pen; \
		int s_ = pen1 < pen2 ? 1 : -1; \
		if (sep > best_sep + 0.001f) { best_sep = sep; best_n = n_; best_type = (type_arg); best_sign = s_; } \
	} while(0)

	// Family 1: cylinder A axis.
	CYL_CYL_TEST_AXIS(axis_a, 0, 0);

	// Family 2: cylinder B axis.
	CYL_CYL_TEST_AXIS(axis_b, 1, 0);

	// Family 3: cross(axis_a, axis_b). Only if axes aren't parallel.
	{
		v3 cross_ab = cross(axis_a, axis_b);
		float cross_len2 = len2(cross_ab);
		if (cross_len2 > 1e-8f) {
			v3 ax = scale(cross_ab, 1.0f / sqrtf(cross_len2));
			if (dot(ax, sub(a.center, b.center)) < 0.0f) ax = neg(ax);
			CYL_CYL_TEST_AXIS(ax, 2, 0);
		}
	}

	// Family 4: "radial" direction — from closest points between the two axis segments.
	// When axes cross at a point (radial = 0), fall back to center-to-center direction.
	// This is the physically-preferred axis for side-side contacts, so it wins TIES
	// (uses >= instead of the biased > comparison that other families use).
	{
		v3 cpa, cpb;
		segments_closest_points(a_p, a_q, b_p, b_q, &cpa, &cpb);
		v3 radial = sub(cpb, cpa);
		float rl2 = len2(radial);
		if (rl2 < 1e-12f) { radial = sub(b.center, a.center); rl2 = len2(radial); }
		if (rl2 > 1e-12f) {
			v3 ax = scale(radial, 1.0f / sqrtf(rl2));
			if (dot(ax, sub(a.center, b.center)) < 0.0f) ax = neg(ax);
			float a_half = cyl_half_support(axis_a, a.half_height, a.radius, ax);
			float b_half = cyl_half_support(axis_b, b.half_height, b.radius, ax);
			float a_center_proj = dot(a.center, ax);
			float b_center_proj = dot(b.center, ax);
			float a_min_v = a_center_proj - a_half, a_max_v = a_center_proj + a_half;
			float b_min_v = b_center_proj - b_half, b_max_v = b_center_proj + b_half;
			if (a_max_v < b_min_v || b_max_v < a_min_v) return 0;
			float pen1 = a_max_v - b_min_v, pen2 = b_max_v - a_min_v;
			float pen = pen1 < pen2 ? pen1 : pen2;
			float sep = -pen;
			if (sep >= best_sep) { best_sep = sep; best_n = ax; best_type = 3; }
		}
	}

	#undef CYL_CYL_TEST_AXIS

	if (best_sep > 0.0f || best_type < 0) return 0;
	if (!manifold) return 1;

	// --- Contact generation ---
	if (best_type == 0) {
		// A's cap wins as reference. Project B's cylinder surface onto A's cap plane.
		// B's surface offset uses the same perp-to-axis formula as cyl-hull face gen.
		v3 cap_n = scale(axis_a, (float)best_sign);
		float cap_d = dot(a.center, cap_n) + a.half_height;
		v3 neg_cap = neg(cap_n);
		float b_axial = dot(neg_cap, axis_b);
		v3 b_perp = sub(neg_cap, scale(axis_b, b_axial));
		float bp2 = len2(b_perp);
		v3 b_off = bp2 > 1e-8f ? scale(b_perp, b.radius / sqrtf(bp2)) : V3(0,0,0);
		v3 bsp = add(b_p, b_off);
		v3 bsq = add(b_q, b_off);
		float d_bp = dot(bsp, cap_n);
		float d_bq = dot(bsq, cap_n);
		int cp = 0;
		v3 points[2]; float depths[2];
		if (d_bp < cap_d) {
			v3 radial = sub(bsp, add(a.center, scale(cap_n, a.half_height)));
			radial = sub(radial, scale(cap_n, dot(radial, cap_n)));
			if (len2(radial) <= a.radius * a.radius) { points[cp] = bsp; depths[cp] = cap_d - d_bp; cp++; }
		}
		if (d_bq < cap_d) {
			v3 radial = sub(bsq, add(a.center, scale(cap_n, a.half_height)));
			radial = sub(radial, scale(cap_n, dot(radial, cap_n)));
			if (len2(radial) <= a.radius * a.radius) { points[cp] = bsq; depths[cp] = cap_d - d_bq; cp++; }
		}
		if (cp == 0) { cp = 1; points[0] = add(a.center, scale(cap_n, a.half_height)); depths[0] = -best_sep; }
		manifold->count = cp;
		for (int i = 0; i < cp; i++)
			manifold->contacts[i] = (Contact){ .point = points[i], .normal = cap_n, .penetration = depths[i] };
		return 1;
	}

	if (best_type == 1) {
		// B's cap wins as reference. Same as type 0 with A/B roles swapped.
		v3 cap_n = scale(axis_b, (float)best_sign);
		float cap_d = dot(b.center, cap_n) + b.half_height;
		v3 neg_cap = neg(cap_n);
		float a_axial = dot(neg_cap, axis_a);
		v3 a_perp = sub(neg_cap, scale(axis_a, a_axial));
		float ap2 = len2(a_perp);
		v3 a_off = ap2 > 1e-8f ? scale(a_perp, a.radius / sqrtf(ap2)) : V3(0,0,0);
		v3 asp = add(a_p, a_off);
		v3 asq = add(a_q, a_off);
		float d_ap = dot(asp, cap_n);
		float d_aq = dot(asq, cap_n);
		int cp = 0;
		v3 points[2]; float depths[2];
		if (d_ap < cap_d) {
			v3 radial = sub(asp, add(b.center, scale(cap_n, b.half_height)));
			radial = sub(radial, scale(cap_n, dot(radial, cap_n)));
			if (len2(radial) <= b.radius * b.radius) { points[cp] = asp; depths[cp] = cap_d - d_ap; cp++; }
		}
		if (d_aq < cap_d) {
			v3 radial = sub(asq, add(b.center, scale(cap_n, b.half_height)));
			radial = sub(radial, scale(cap_n, dot(radial, cap_n)));
			if (len2(radial) <= b.radius * b.radius) { points[cp] = asq; depths[cp] = cap_d - d_aq; cp++; }
		}
		if (cp == 0) { cp = 1; points[0] = add(b.center, scale(cap_n, b.half_height)); depths[0] = -best_sep; }
		manifold->count = cp;
		for (int i = 0; i < cp; i++)
			manifold->contacts[i] = (Contact){ .point = points[i], .normal = cap_n, .penetration = depths[i] };
		return 1;
	}

	// best_type == 2 or 3: edge/radial contact — single point from classify.
	{
		v3 cpa, cpb;
		segments_closest_points(a_p, a_q, b_p, b_q, &cpa, &cpb);
		CylFeature feat = cyl_classify_point(cpb, a.center, a.rotation, a.half_height, a.radius);
		manifold->count = 1;
		manifold->contacts[0] = (Contact){
			.point = feat.surface_pt,
			.normal = neg(best_n),
			.penetration = -best_sep,
		};
		return 1;
	}
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
	const char* names[] = {"sphere", "capsule", "box", "hull", "cylinder"};
	printf("  --- narrowphase breakdown ---\n");
	double total = 0;
	int total_calls = 0;
	for (int a = 0; a < 5; a++) for (int b = a; b < 5; b++) {
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
	g_sat_hillclimb_enabled = w->sat_hillclimb_enabled;
	// Canonical ordering: lower type first for upper-triangle dispatch,
	// lower body index first for same-type pairs (deterministic shape A).
	if (w->body_cold[i].shapes[0].type > w->body_cold[j].shapes[0].type || (w->body_cold[i].shapes[0].type == w->body_cold[j].shapes[0].type && i > j)) {
		int tmp = i; i = j; j = tmp;
	}

	ShapeInternal* s0 = &w->body_cold[i].shapes[0];
	ShapeInternal* s1 = &w->body_cold[j].shapes[0];
	BodyState* bs0 = &w->body_state[i];
	BodyState* bs1 = &w->body_state[j];
	InternalManifold im = { .body_a = i, .body_b = j };

	// Warm cache lookup: used for SAT hints, geometry caching, and passed to pre_solve.
	uint64_t pkey = body_pair_key(i, j);
	WarmManifold* wm = w->warm_start_enabled ? map_get_ptr(w->warm_cache, pkey) : NULL;

	// SAT hint from warm cache (skip for box-box — 15 axes is cheap).
	int* hp = NULL;
	int hint = -1;
	int uses_sat = (s0->type >= SHAPE_BOX && s0->type != SHAPE_CYLINDER && s1->type >= SHAPE_BOX && s1->type != SHAPE_CYLINDER);
	int uses_hint = uses_sat && !(s0->type == SHAPE_BOX && s1->type == SHAPE_BOX && !w->box_use_hull);
	if (uses_hint && wm && w->sat_hint_enabled) { hint = wm->sat_axis; hp = &hint; }

	// Incremental narrowphase fast path: validate cached feature pair, re-clip without SAT.
	if (wm && wm->cached_pair.type == 1 && uses_sat && w->incremental_np_enabled) {
		int refreshed = 0;
		if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX && !w->box_use_hull)
			refreshed = refresh_box_box_face(make_box(bs0, s0), make_box(bs1, s1), &im.m, &wm->cached_pair);
		else
			refreshed = refresh_hull_hull_face((ConvexHull){ s0->type == SHAPE_BOX ? &s_unit_box_hull : s0->hull.hull, bs0->position, bs0->rotation, s0->type == SHAPE_BOX ? s0->box.half_extents : s0->hull.scale }, (ConvexHull){ s1->type == SHAPE_BOX ? &s_unit_box_hull : s1->hull.hull, bs1->position, bs1->rotation, s1->type == SHAPE_BOX ? s1->box.half_extents : s1->hull.scale }, &im.m, &wm->cached_pair);
		if (refreshed) {
			im.warm = wm;
			wm->stale = 0;
			np_call_acc[np_pair_idx(s0->type, s1->type)]++;
			apush(*manifolds, im);
			return;
		}
		wm->cached_pair.type = 0; // invalidated, fall through to full SAT
	}
	// Edge-edge cached pairs (type==2): always invalidate for V1
	if (wm && wm->cached_pair.type == 2) wm->cached_pair.type = 0;

	// Upper-triangle dispatch: simple pairs first, then SAT-based pairs.
	int hit = 0;
	CachedFeaturePair out_pair = {0};

	if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_SPHERE)
		hit = collide_sphere_sphere(make_sphere(bs0, s0), make_sphere(bs1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CAPSULE)
		hit = collide_sphere_capsule(make_sphere(bs0, s0), make_capsule(bs1, s1), &im.m);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_BOX)
		hit = collide_sphere_box(make_sphere(bs0, s0), make_box(bs1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CAPSULE)
		hit = collide_capsule_capsule(make_capsule(bs0, s0), make_capsule(bs1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_BOX)
		hit = collide_capsule_box(make_capsule(bs0, s0), make_box(bs1, s1), &im.m);
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX) {
		if (w->box_use_hull) hit = collide_hull_hull_ex((ConvexHull){ &s_unit_box_hull, bs0->position, bs0->rotation, s0->box.half_extents }, (ConvexHull){ &s_unit_box_hull, bs1->position, bs1->rotation, s1->box.half_extents }, &im.m, hp, &out_pair);
		else hit = collide_box_box_ex(make_box(bs0, s0), make_box(bs1, s1), &im.m, hp, &out_pair);
	}
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_HULL)
		hit = collide_hull_hull_ex((ConvexHull){ &s_unit_box_hull, bs0->position, bs0->rotation, s0->box.half_extents }, make_convex_hull(bs1, s1), &im.m, hp, &out_pair);
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL)
		hit = collide_sphere_hull(make_sphere(bs0, s0), make_convex_hull(bs1, s1), &im.m);
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL)
		hit = collide_capsule_hull(make_capsule(bs0, s0), make_convex_hull(bs1, s1), &im.m);
	else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL)
		hit = collide_hull_hull_ex(make_convex_hull(bs0, s0), make_convex_hull(bs1, s1), &im.m, hp, &out_pair);

	// Cylinder pairs: native narrowphase with cylinder as shape A.
	// Dispatch puts cyl as s1 (higher enum). Call with cyl as first arg, flip normals.
	else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CYLINDER) {
		hit = collide_cylinder_sphere(make_cylinder(bs1, s1), make_sphere(bs0, s0), &im.m);
		for (int c = 0; c < im.m.count; c++) im.m.contacts[c].normal = neg(im.m.contacts[c].normal);
	}
	else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CYLINDER) {
		hit = collide_cylinder_capsule(make_cylinder(bs1, s1), make_capsule(bs0, s0), &im.m);
		for (int c = 0; c < im.m.count; c++) im.m.contacts[c].normal = neg(im.m.contacts[c].normal);
	}
	else if (s0->type == SHAPE_BOX && s1->type == SHAPE_CYLINDER) {
		hit = collide_cylinder_box(make_cylinder(bs1, s1), make_box(bs0, s0), &im.m);
		for (int c = 0; c < im.m.count; c++) im.m.contacts[c].normal = neg(im.m.contacts[c].normal);
	}
	else if (s0->type == SHAPE_HULL && s1->type == SHAPE_CYLINDER) {
		hit = collide_cylinder_hull(make_cylinder(bs1, s1), make_convex_hull(bs0, s0), &im.m);
		for (int c = 0; c < im.m.count; c++) im.m.contacts[c].normal = neg(im.m.contacts[c].normal);
	}
	else if (s0->type == SHAPE_CYLINDER && s1->type == SHAPE_CYLINDER) {
		hit = collide_cylinder_cylinder(make_cylinder(bs0, s0), make_cylinder(bs1, s1), &im.m);
	}

	// Store SAT hint back to warm cache
	if (hp && wm) wm->sat_axis = hint;

	int idx = np_pair_idx(s0->type, s1->type);
	np_call_acc[idx]++;
	// Cache geometry and feature pair for next frame reuse
	if (hit) {
		if (!wm && w->warm_start_enabled) { wm = map_get_ptr(w->warm_cache, pkey); }
		if (wm && out_pair.type) wm->cached_pair = out_pair;
		im.warm = wm;
		apush(*manifolds, im);
	}
}

// Broadphase moved to broadphase.c.

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
static int ray_cylinder(v3 ro, v3 rd, Cylinder cyl, float max_t, float* t_out, v3* n_out)
{
	v3 axis = rotate(cyl.rotation, V3(0, 1, 0));
	v3 oc = sub(ro, cyl.center);
	float rd_a = dot(rd, axis);
	float oc_a = dot(oc, axis);
	float r2 = cyl.radius * cyl.radius;

	float best = max_t + 1.0f;
	v3 best_n = {0};
	int hit = 0;

	// Infinite cylinder: |oc + t*rd - (oc_a + t*rd_a)*axis|^2 = r^2
	float a = dot(rd, rd) - rd_a * rd_a;
	float b = dot(oc, rd) - oc_a * rd_a;
	float c = dot(oc, oc) - oc_a * oc_a - r2;
	if (fabsf(a) > 1e-8f) {
		float disc = b * b - a * c;
		if (disc >= 0.0f) {
			float sq = sqrtf(disc);
			float t = (-b - sq) / a;
			if (t < 0.0f) t = (-b + sq) / a;
			if (t >= 0.0f && t < best) {
				float y = oc_a + t * rd_a;
				if (y >= -cyl.half_height && y <= cyl.half_height) {
					best = t; hit = 1;
					v3 hp = add(ro, scale(rd, t));
					best_n = norm(sub(hp, add(cyl.center, scale(axis, y))));
				}
			}
		}
	}

	// Flat caps: plane at center +/- half_height * axis.
	for (int ci = 0; ci < 2; ci++) {
		float sign = ci == 0 ? -1.0f : 1.0f;
		float denom = rd_a * sign;
		if (fabsf(denom) < 1e-8f) continue;
		float t = (sign * cyl.half_height - oc_a) / (rd_a);
		if (t < 0.0f || t >= best) continue;
		v3 hp = add(ro, scale(rd, t));
		v3 d = sub(hp, add(cyl.center, scale(axis, sign * cyl.half_height)));
		if (dot(d, d) <= r2) {
			best = t; hit = 1;
			best_n = scale(axis, sign);
		}
	}

	if (hit && best <= max_t) { *t_out = best; *n_out = best_n; return 1; }
	return 0;
}

static int ray_body(WorldInternal* w, int body_idx, v3 origin, v3 dir, float max_t, float* t_out, v3* n_out)
{
	BodyState* bs = &w->body_state[body_idx];
	BodyCold* bc = &w->body_cold[body_idx];
	float best_t = max_t;
	v3 best_n = {0};
	int found = 0;
	for (int i = 0; i < asize(bc->shapes); i++) {
		ShapeInternal* s = &bc->shapes[i];
		float t; v3 n;
		int hit = 0;
		switch (s->type) {
		case SHAPE_SPHERE:   hit = ray_sphere(origin, dir, make_sphere(bs, s), best_t, &t, &n); break;
		case SHAPE_CAPSULE:  hit = ray_capsule(origin, dir, make_capsule(bs, s), best_t, &t, &n); break;
		case SHAPE_BOX:      hit = ray_box(origin, dir, make_box(bs, s), best_t, &t, &n); break;
		case SHAPE_HULL:     hit = ray_hull(origin, dir, make_convex_hull(bs, s), best_t, &t, &n); break;
		case SHAPE_CYLINDER: hit = ray_cylinder(origin, dir, make_cylinder(bs, s), best_t, &t, &n); break;
		}
		if (hit && t < best_t) { best_t = t; best_n = n; found = 1; }
	}
	if (found) { *t_out = best_t; *n_out = best_n; }
	return found;
}
