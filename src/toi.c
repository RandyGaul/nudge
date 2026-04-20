// toi.c -- Time-of-impact advancement against static world geometry.
//
// Used by the CCD post-solve pass (wired in step 6): each body enrolled in
// WorldInternal.fast_body_list by the speed test (inertia.c) is swept from
// t=0 to t=1 and clamped to the earliest static-world impact fraction.
//
// Algorithm: bilateral advancement (Erin Catto, GDC 2010 "Continuous Collision
// Detection", refined in Box2D v3 `b2TimeOfImpact` / distance.c:1138). We use
// GJK to compute the distance at the current t, extract the closest-feature
// pair from the simplex, and build a separating-axis function that gives a
// monotonic 1D root to find. Bisection+secant locates the TOI on that axis;
// when the active witness features change, the outer loop reseeds by
// re-running GJK at the new t. See Catto GDC 2010 for the derivation.
//
// Static-only simplification: shape B's pose is fixed, so we only integrate
// body A. Dynamic-vs-dynamic is handled by speculative contacts per plan.
//
// Impact response: defer to next frame. We clamp the pose at t* - slop, leave
// velocity untouched, and let the next frame's speculative pass apply
// restitution via the replace-rule (see pre_solve_manifold.inc).

// Perf telemetry: counts sep-fn classifications + Gauss decisions. Printed
// by toi_print_stats() (called from bench harnesses).
static int toi_stat_seed_edge_edge = 0;
static int toi_stat_seed_edge_reject_gauss = 0;
static int toi_stat_seed_face_a = 0;
static int toi_stat_seed_face_b = 0;
static int toi_stat_seed_points = 0;
static int toi_stat_hi_gauss_reject = 0;

static void toi_print_stats(void)
{
	int total = toi_stat_seed_edge_edge + toi_stat_seed_edge_reject_gauss
		+ toi_stat_seed_face_a + toi_stat_seed_face_b + toi_stat_seed_points;
	if (total == 0) { printf("  toi-sepfn: no classifications recorded\n"); return; }
	printf("  toi-sepfn: total=%d  POINTS=%d  FACE_A=%d  FACE_B=%d  EDGE_EDGE=%d  (gauss-reject seed=%d hi=%d)\n",
		total, toi_stat_seed_points, toi_stat_seed_face_a, toi_stat_seed_face_b,
		toi_stat_seed_edge_edge, toi_stat_seed_edge_reject_gauss, toi_stat_hi_gauss_reject);
}
static void toi_reset_stats(void)
{
	toi_stat_seed_edge_edge = 0;
	toi_stat_seed_edge_reject_gauss = 0;
	toi_stat_seed_face_a = 0;
	toi_stat_seed_face_b = 0;
	toi_stat_seed_points = 0;
	toi_stat_hi_gauss_reject = 0;
}

#define TOI_OUTER_ITER_MAX       20
#define TOI_ROOT_ITER_MAX        50
#define TOI_PUSHBACK_ITER_MAX    8
#define TOI_SUBSTEP_MAX          4
#define TOI_TARGET_SEPARATION    (0.5f * LINEAR_SLOP)
#define TOI_CONVERGE_TOLERANCE   (0.25f * LINEAR_SLOP)
#define TOI_SLOP_FRACTION        0.001f   // pull back from impact instant by this fraction of the step

// Integrate a body pose forward by fraction t of the step (linear + angular).
static void toi_pose_at(BodyState* bs, BodyHot* bh, float dt, float t, v3* out_pos, quat* out_rot)
{
	float h = t * dt;
	*out_pos = add(bs->position, scale(bh->velocity, h));
	v3 ww = bh->angular_velocity;
	quat q = bs->rotation;
	quat spin = { ww.x, ww.y, ww.z, 0.0f };
	quat dq = mul(spin, q);
	quat qn = { q.x + 0.5f * h * dq.x, q.y + 0.5f * h * dq.y, q.z + 0.5f * h * dq.z, q.w + 0.5f * h * dq.w };
	*out_rot = quat_norm(qn);
}

// Build a GJK shape for a single ShapeInternal at a given body pose.
// Duplicates the world-frame setup from collision.c's make_* helpers but
// parameterized on arbitrary pose (not the body's current state).
static GJK_Shape toi_shape_at_pose(ShapeInternal* s, v3 body_pos, quat body_rot)
{
	v3 world_pos = add(body_pos, rotate(body_rot, s->local_pos));
	quat world_rot = quat_mul(body_rot, s->local_rot);
	switch (s->type) {
	case SHAPE_SPHERE:
		return gjk_sphere(world_pos, s->sphere.radius);
	case SHAPE_CAPSULE: {
		v3 axis = rotate(world_rot, V3(0, s->capsule.half_height, 0));
		return gjk_capsule(sub(world_pos, axis), add(world_pos, axis), s->capsule.radius);
	}
	case SHAPE_BOX:
		// Route through the unit-box hull so boxes participate in hull-aware
		// sep-fn paths (edge-edge Gauss validation). Semantically identical
		// to gjk_box for support calls, but now exposes half-edge + face-plane
		// data for Gregorius edge-edge axis validation.
		return gjk_hull_scaled(hull_unit_box(), world_pos, world_rot, s->box.half_extents);
	case SHAPE_HULL: {
		const Hull* h = s->hull.hull;
		return gjk_hull_scaled(h, world_pos, world_rot, s->hull.scale);
	}
	case SHAPE_MESH:
	case SHAPE_HEIGHTFIELD:
	default:
		// Static-only. Dynamic bodies never carry these.
		return gjk_sphere(world_pos, 0.0f);
	}
}

// Declare test-only hook points so tests.c can exercise the feature-
// extraction and separation-function paths without going through full world
// setup. These are static helpers; tests include toi.c via main.c's unity
// build and can call them directly.

// -----------------------------------------------------------------------------
// Separation function (Erin Catto, GDC 2010 "Continuous Collision Detection";
// Box2D v3 `b2SeparationFunction` / distance.c:932). Generalized to 3D:
// POINTS / FACE_A / FACE_B (edge-edge falls through to FACE_A with the
// GJK-derived axis, which is correct when edge-edge witnesses lie on a
// common face tangent — acceptable loss of generality for this MVP; full
// 3D edge-edge support is future work).

// Separation function (Box2D v3 b2SeparationFunction / distance.c:932, adapted
// to 3D). Three types classified from GJK simplex unique-feature counts kA/kB:
//   POINTS (kA=1, kB=1)   : vertex-vs-vertex. Axis stored in world frame.
//                           sep(t) = dot(pB(idxB) - pA_t(idxA), axis_world).
//   FACE_A (kA>=2)        : face/edge on A. Axis stored in A's local frame;
//                           rotates with A so it tracks A's face normal.
//                           sep(t) = dot(pB(idxB) - pA_t(idxA_ref), rotate(rotA(t), axis_local)).
//   FACE_B (kA=1, kB>=2)  : face/edge on B. Since B is static, axis stored
//                           in world frame (equivalent to B-local * static rot).
// All evaluation goes through gjk_support_feature() at FIXED vertex indices —
// stable across pose sweeps, no live support calls during root find.
//
// Key differences from the previous (now-deleted) 4-type design:
//   - No separate EDGE_EDGE type. True 3D edge-vs-edge gets classified as
//     kA=2, kB=2 → FACE_A with the GJK closest-direction axis; this is
//     geometrically correct (the GJK axis is perpendicular to both edges when
//     they are skew, and lies in the edge plane when they are coplanar).
//     The old EDGE_EDGE used cross(eA, eB) which misfires whenever GJK pulls
//     multiple vertices on each side for a face-vs-edge contact (e.g.
//     capsule-vs-box face).
//   - POINTS axis stored in world, not A-local. For translational sweeps
//     this gives linear sep(t); the previous A-local storage rotated the
//     axis with A every sample and produced false zero-crossings.
//   - Evaluate uses fixed vertex indices via gjk_support_feature rather than
//     world-space witnesses re-rotated through A's pose. Matches box2d's
//     proxyA->points[cache->indexA] approach and is simpler.
// EDGE_EDGE (Gregorius GDC 2013 "The Separating Axis Test"):
// Both shapes are hulls and the GJK simplex identifies an edge-vs-edge
// closest-feature pair (kA=2, kB=2 from distinct vertices that form an
// edge on each side). Axis is the live cross product cross(eA(t), eB(t))
// normalized — tracks both bodies' rotations exactly. The pair is only a
// valid separating-axis candidate when the four arc endpoints formed by
// each edge's two adjacent face normals intersect on the unit sphere
// (Gauss-map overlap test). If the arcs stop intersecting mid-sweep, the
// edge-edge axis no longer bounds the true min-separation, so the outer
// TOI loop reseeds on the next iteration.
enum { TOI_POINTS = 0, TOI_FACE_A = 1, TOI_FACE_B = 2, TOI_EDGE_EDGE = 3 };

typedef struct ToiSepFn
{
	int type;
	int idxA;   // A vertex index (reference feature, stable across t)
	int idxB;   // B vertex index (reference feature, stable across t)
	v3 axis;    // POINTS/FACE_B: world frame; FACE_A: A-local frame
	// Edge-edge only. idxA2/idxB2 are the second endpoint of each edge.
	// Adjacent face normals stored in body-local frames so each body's
	// rotation propagates to the Gauss arc at eval time.
	int idxA2, idxB2;
	v3 nA0_local, nA1_local;  // A's two adjacent face normals (A-local)
	v3 nB0_local, nB1_local;  // B's two adjacent face normals (B-local)
	int axis_sign;             // +/-1, orienting cross(eA, eB) so +sep = B above A's edge plane
} ToiSepFn;

// Gregorius Gauss-arc overlap predicate. Arc (a0, a1) is A's edge-adjacent
// face normals; arc (b0, b1) is NEGATED B-face normals. Returns true iff
// the two arcs cross on the unit sphere — the necessary condition for
// cross(eA, eB) to be the minimum-separating axis. Branchless, ~20 flops.
static inline int toi_gauss_arcs_overlap(v3 a0, v3 a1, v3 b0, v3 b1)
{
	v3 ba = cross(a0, a1);
	v3 dc = cross(b0, b1);
	float cba = dot(b0, ba), dba = dot(b1, ba);
	float adc = dot(a0, dc), bdc = dot(a1, dc);
	return (cba * dba < 0.0f) && (adc * bdc < 0.0f) && (cba * bdc > 0.0f);
}

// Find the half-edge on a GJK hull shape whose origin = v0 and destination
// = v1. Returns half-edge index, or -1. GJK_Shape carries edge_origin[] and
// edge_twin[]; destination = origin[twin[e]].
static int toi_hull_find_halfedge(const GJK_Shape* s, int v0, int v1)
{
	if (s->type != GJK_HULL || !s->hull.edge_twin || !s->hull.edge_origin) return -1;
	int ec = s->hull.edge_count;
	for (int e = 0; e < ec; e++) {
		int o = s->hull.edge_origin[e];
		int d = s->hull.edge_origin[s->hull.edge_twin[e]];
		if (o == v0 && d == v1) return e;
	}
	return -1;
}

// Fetch adjacent face plane normals for a half-edge (in LOCAL frame of the
// hull — shared across faces either side). Returns 0 on failure, 1 on
// success. edge_face[e] gives the left-face of half-edge e; edge_face of
// its twin gives the right-face.
static int toi_hull_edge_face_normals(const GJK_Shape* s, int he, v3* out_n0, v3* out_n1)
{
	if (s->type != GJK_HULL || !s->hull.edge_face || !s->hull.planes || !s->hull.edge_twin) return 0;
	int f0 = s->hull.edge_face[he];
	int f1 = s->hull.edge_face[s->hull.edge_twin[he]];
	*out_n0 = s->hull.planes[f0].normal;
	*out_n1 = s->hull.planes[f1].normal;
	return 1;
}

// Extract (kA, kB) = counts of unique supporting features on each side of
// the Minkowski simplex. `out_feats_a[0..kA-1]` and `out_feats_b[0..kB-1]`
// return the unique support-index lists. Static helper used by the
// separation-function builder and by tests.
static void toi_simplex_feature_counts(const GJK_Simplex* s, int* out_kA, int* out_kB, int out_feats_a[4], int out_feats_b[4])
{
	int nA = 0, nB = 0;
	for (int i = 0; i < s->count; i++) {
		int fA = s->v[i].feat1, fB = s->v[i].feat2;
		int found = 0;
		for (int j = 0; j < nA; j++) if (out_feats_a[j] == fA) { found = 1; break; }
		if (!found) out_feats_a[nA++] = fA;
		found = 0;
		for (int j = 0; j < nB; j++) if (out_feats_b[j] == fB) { found = 1; break; }
		if (!found) out_feats_b[nB++] = fB;
	}
	*out_kA = nA;
	*out_kB = nB;
}

// Compute world-frame axis for this sep-fn at the given A pose. FACE_A
// rotates the stored A-local axis into world; FACE_B and POINTS already
// store world-frame axes. EDGE_EDGE isn't handled here — it recomputes the
// cross-product axis live inside its evaluate paths.
static v3 toi_axis_world(const ToiSepFn* f, quat rot_A)
{
	if (f->type == TOI_FACE_A) return rotate(rot_A, f->axis);
	return f->axis;
}


// Build a separation function from the GJK simplex at seed pose.
//   kA=1, kB=1    -> POINTS (axis world-frame = normalize(pB - pA))
//   kA>=2         -> FACE_A (axis stored A-local, rotates with A)
//   kA=1, kB>=2   -> FACE_B (axis stored world-frame; B static)
// Falls back to POINTS if GJK returned a zero-direction (degenerate).
static ToiSepFn toi_make_sep_fn(const GJK_Simplex* simplex, const GJK_Result* r0, quat rot_A, GJK_Shape* shapeA_at_seed, GJK_Shape* shapeB)
{
	ToiSepFn f = {0};
	int kA = 0, kB = 0;
	int fa[4] = {0}, fb[4] = {0};
	toi_simplex_feature_counts(simplex, &kA, &kB, fa, fb);

	v3 diff = sub(r0->point2, r0->point1);
	float dlen = v3_len(diff);
	v3 axis_world = dlen > 1e-8f ? scale(diff, 1.0f / dlen) : V3(0, 1, 0);

	// Pick representative feature indices (first entry from each side's
	// unique-supports list).
	f.idxA = fa[0];
	f.idxB = fb[0];

	// Classification. If GJK direction is degenerate (shapes overlap at seed)
	// the outer TOI loop has already bailed before we get here, so dlen is
	// effectively always > tol here; fallback is defensive.

	// Edge-edge (Gregorius): exactly 2 distinct supports on each side AND
	// both shapes carry hull half-edge + face-plane data AND the Gauss arcs
	// overlap. If the predicate fails at seed, cross(eA, eB) is not the
	// minimum-separating direction; fall through to FACE_A.
	// Disable with TOI_DISABLE_EDGE_EDGE=1 for benchmarking.
	static int ee_disabled = -1;
	if (ee_disabled < 0) ee_disabled = getenv("TOI_DISABLE_EDGE_EDGE") ? 1 : 0;
	if (!ee_disabled && kA == 2 && kB == 2 && dlen > 1e-8f
	 && shapeA_at_seed->hull.edge_twin && shapeA_at_seed->hull.planes
	 && shapeB->hull.edge_twin && shapeB->hull.planes) {
		int eA_he = toi_hull_find_halfedge(shapeA_at_seed, fa[0], fa[1]);
		int eB_he = toi_hull_find_halfedge(shapeB, fb[0], fb[1]);
		if (eA_he >= 0 && eB_he >= 0) {
			// Adjacent face-plane normals in each hull's LOCAL frame —
			// Hull.planes[].normal is stored local. That's what we want
			// to cache: future t poses rotate them into world via each
			// body's rotation.
			v3 nA0_l, nA1_l, nB0_l, nB1_l;
			if (toi_hull_edge_face_normals(shapeA_at_seed, eA_he, &nA0_l, &nA1_l)
			 && toi_hull_edge_face_normals(shapeB, eB_he, &nB0_l, &nB1_l)) {
				// Gauss arc test at seed: rotate local normals to world and
				// check overlap. shapeA is at the seed pose; shapeB's basis
				// columns encode B's static rotation.
				v3 nA0_w = rotate(rot_A, nA0_l);
				v3 nA1_w = rotate(rot_A, nA1_l);
				#define MB(col0,col1,col2,v) V3( \
					(col0).x*(v).x + (col1).x*(v).y + (col2).x*(v).z, \
					(col0).y*(v).x + (col1).y*(v).y + (col2).y*(v).z, \
					(col0).z*(v).x + (col1).z*(v).y + (col2).z*(v).z)
				v3 nB0_w = MB(shapeB->hull.col0, shapeB->hull.col1, shapeB->hull.col2, nB0_l);
				v3 nB1_w = MB(shapeB->hull.col0, shapeB->hull.col1, shapeB->hull.col2, nB1_l);
				#undef MB
				if (!toi_gauss_arcs_overlap(nA0_w, nA1_w, neg(nB0_w), neg(nB1_w))) {
					toi_stat_seed_edge_reject_gauss++;
				} else {
					v3 a0w = gjk_support_feature(shapeA_at_seed, fa[0]);
					v3 a1w = gjk_support_feature(shapeA_at_seed, fa[1]);
					v3 b0w = gjk_support_feature(shapeB, fb[0]);
					v3 b1w = gjk_support_feature(shapeB, fb[1]);
					v3 ee = cross(sub(a1w, a0w), sub(b1w, b0w));
					float L = v3_len(ee);
					if (L > 1e-8f) {
						v3 n_w = scale(ee, 1.0f / L);
						int sign = dot(sub(b0w, a0w), n_w) >= 0.0f ? 1 : -1;
						f.type = TOI_EDGE_EDGE;
						f.idxA = fa[0]; f.idxA2 = fa[1];
						f.idxB = fb[0]; f.idxB2 = fb[1];
						f.axis_sign = sign;
						f.nA0_local = nA0_l; f.nA1_local = nA1_l;
						f.nB0_local = nB0_l; f.nB1_local = nB1_l;
						toi_stat_seed_edge_edge++;
						return f;
					}
				}
			}
		}
	}

	if (kA >= 2 && dlen > 1e-8f) {
		// A-side face/edge. Axis rotates with A.
		f.type = TOI_FACE_A;
		f.axis = rotate(quat_inv(rot_A), axis_world);
		// Re-evaluate support along the stored axis to pin idxA to A's
		// deepest vertex on that axis — gives stable feature across sweep.
		v3 pa; int fi_a;
		gjk_support(shapeA_at_seed, axis_world, &fi_a, pa);
		f.idxA = fi_a;
		// And B's deepest vertex in -axis (toward A).
		v3 pb; int fi_b;
		gjk_support(shapeB, neg(axis_world), &fi_b, pb);
		f.idxB = fi_b;
		toi_stat_seed_face_a++;
		return f;
	}
	if (kB >= 2 && dlen > 1e-8f) {
		// B-side face/edge. B static, so world-frame axis is equivalent to
		// B-local axis rotated by B's fixed rotation.
		f.type = TOI_FACE_B;
		f.axis = axis_world;
		v3 pa; int fi_a;
		gjk_support(shapeA_at_seed, axis_world, &fi_a, pa);
		f.idxA = fi_a;
		v3 pb; int fi_b;
		gjk_support(shapeB, neg(axis_world), &fi_b, pb);
		f.idxB = fi_b;
		toi_stat_seed_face_b++;
		return f;
	}
	// Vertex-vertex (POINTS). Axis stored in world (fixed direction from
	// A's seed vertex to B's seed vertex). Translational sweeps produce a
	// linear sep(t); any re-frame (A-local / B-local) would inject rotation
	// noise when no face feature is involved.
	f.type = TOI_POINTS;
	f.axis = axis_world;
	toi_stat_seed_points++;
	return f;
}

// Edge-edge separation at time t: live cross(eA, eB) axis, signed projection
// of one edge endpoint pair onto it. Both edge endpoints on each side come
// from stored vertex indices — shapeA_at_t rotates A; shapeB static. If the
// axis length is degenerate (edges parallel mid-sweep) returns 0.
// `out_gauss_valid` receives the Gregorius Gauss-arc check at the current
// poses (rotate stored local face normals into world). Callers in the
// "declare separation" path (FindMinSep) use this to force a root-find
// instead of returning "separated" under an axis that's no longer the
// minimum-separating direction.
static float toi_edge_edge_sep(const ToiSepFn* f, GJK_Shape* shapeA_at_t, GJK_Shape* shapeB, quat rot_A, int* out_gauss_valid)
{
	v3 a0 = gjk_support_feature(shapeA_at_t, f->idxA);
	v3 a1 = gjk_support_feature(shapeA_at_t, f->idxA2);
	v3 b0 = gjk_support_feature(shapeB, f->idxB);
	v3 b1 = gjk_support_feature(shapeB, f->idxB2);
	v3 n = cross(sub(a1, a0), sub(b1, b0));
	float L = v3_len(n);
	if (L <= 1e-8f) { if (out_gauss_valid) *out_gauss_valid = 0; return 0.0f; }
	n = scale(n, (float)f->axis_sign / L);
	if (out_gauss_valid) {
		// Rotate A's face normals (A-local) by current rot_A; B's normals
		// (B-local) by B's static rot stored in shapeB's basis columns.
		v3 a0n = rotate(rot_A, f->nA0_local);
		v3 a1n = rotate(rot_A, f->nA1_local);
		#define MB(v) V3( \
			shapeB->hull.col0.x*(v).x + shapeB->hull.col1.x*(v).y + shapeB->hull.col2.x*(v).z, \
			shapeB->hull.col0.y*(v).x + shapeB->hull.col1.y*(v).y + shapeB->hull.col2.y*(v).z, \
			shapeB->hull.col0.z*(v).x + shapeB->hull.col1.z*(v).y + shapeB->hull.col2.z*(v).z)
		v3 b0n = neg(MB(f->nB0_local));
		v3 b1n = neg(MB(f->nB1_local));
		#undef MB
		*out_gauss_valid = toi_gauss_arcs_overlap(a0n, a1n, b0n, b1n);
	}
	float rA = shapeA_at_t->radius, rB = shapeB->radius;
	v3 pa_eff = add(a0, scale(n, rA));
	v3 pb_eff = sub(b0, scale(n, rB));
	return dot(sub(pb_eff, pa_eff), n);
}

// FindMinSeparation (box2d b2FindMinSeparation): find deepest feature along
// the current world-axis at time t and return the separation. Updates the
// stored indices to the newly-found features. For EDGE_EDGE it keeps the
// edge endpoints (idxA/A2, idxB/B2) fixed — live cross-product evaluates
// in closed form. When Gauss arcs fail at the current pose, the cross
// product is no longer the minimum-separating direction; clamp the returned
// sep to the convergence target so the caller does NOT declare "separated"
// and instead enters the root-find / outer reseed path. Under-reporting sep
// is the safe direction — worst case we do extra bisection; alternative
// would be a false "separated through step" and a missed TOI.
static float toi_find_min_sep(ToiSepFn* f, GJK_Shape* shapeA_at_t, GJK_Shape* shapeB, quat rot_A)
{
	if (f->type == TOI_EDGE_EDGE) {
		int gauss_valid = 1;
		float s = toi_edge_edge_sep(f, shapeA_at_t, shapeB, rot_A, &gauss_valid);
		if (!gauss_valid) {
			toi_stat_hi_gauss_reject++;
			if (s > TOI_TARGET_SEPARATION) s = TOI_TARGET_SEPARATION;
		}
		return s;
	}
	v3 axis_world = toi_axis_world(f, rot_A);
	int fi_a, fi_b;
	v3 pa, pb;
	gjk_support(shapeA_at_t, axis_world, &fi_a, pa);
	gjk_support(shapeB, neg(axis_world), &fi_b, pb);
	f->idxA = fi_a;
	f->idxB = fi_b;
	float rA = shapeA_at_t->radius, rB = shapeB->radius;
	v3 pa_eff = add(pa, scale(axis_world, rA));
	v3 pb_eff = sub(pb, scale(axis_world, rB));
	return dot(sub(pb_eff, pa_eff), axis_world);
}

// EvaluateSeparation (box2d b2EvaluateSeparation): compute separation at
// time t with feature indices held FIXED (from a prior FindMinSep). Does
// no support search — just re-fetches the stored vertex indices at the
// current pose via gjk_support_feature. This is what the root finder
// iterates on. Gauss validity isn't checked here — root-find operates on
// a continuous sep-fn inside an already-bracketed sign change, so local
// invalidity doesn't prevent convergence. The outer loop reseeds at the
// advanced t1.
static float toi_eval_sep(const ToiSepFn* f, GJK_Shape* shapeA_at_t, GJK_Shape* shapeB, quat rot_A)
{
	if (f->type == TOI_EDGE_EDGE) return toi_edge_edge_sep(f, shapeA_at_t, shapeB, rot_A, NULL);
	v3 axis_world = toi_axis_world(f, rot_A);
	v3 pa = gjk_support_feature(shapeA_at_t, f->idxA);
	v3 pb = gjk_support_feature(shapeB, f->idxB);
	float rA = shapeA_at_t->radius, rB = shapeB->radius;
	v3 pa_eff = add(pa, scale(axis_world, rA));
	v3 pb_eff = sub(pb, scale(axis_world, rB));
	return dot(sub(pb_eff, pa_eff), axis_world);
}

// Bilateral advancement (Catto GDC 2010; mirrors Box2D v3 b2TimeOfImpact).
// Returns the fraction t in [0, 1] of first contact, or 1.0 if separated
// through the whole step. The BodyState carries A's PRE-STEP pose.
static float toi_pair_advance(ShapeInternal* sa, BodyState* bs, BodyHot* bh, float r_max, float dt, GJK_Shape* static_b)
{
	(void)r_max;
	float t1 = 0.0f;
	float t_max = 1.0f;
	float target = TOI_TARGET_SEPARATION;
	float tol = TOI_CONVERGE_TOLERANCE;

	for (int outer = 0; outer < TOI_OUTER_ITER_MAX; outer++) {
		// Distance query at t1, with simplex extraction for feature analysis.
		v3 pos_A_t1; quat rot_A_t1;
		toi_pose_at(bs, bh, dt, t1, &pos_A_t1, &rot_A_t1);
		GJK_Shape ga_t1 = toi_shape_at_pose(sa, pos_A_t1, rot_A_t1);
		GJK_Simplex simplex;
		GJK_Result r_t1 = gjk_distance_ex(&ga_t1, static_b, NULL, &simplex);

		// "Touching or overlapping at t1": only counts as a hit if the body
		// is still APPROACHING the static surface at this instant. When
		// body A's velocity along the contact normal already points away
		// (e.g. the previous frame's solver applied a bounce impulse),
		// re-clamping to t1 would erase the post-impact motion and freeze
		// the body at the surface. Check the closing speed before bailing.
		// Already overlapping at t1: TOI can't help — clamping to pre_pos
		// just freezes the body inside the static. Let the solver's
		// position correction + speculative bias handle it.
		if (r_t1.distance <= 0.0f) return 1.0f;
		if (r_t1.distance < target + tol) {
			v3 diff = sub(r_t1.point2, r_t1.point1);
			float L = v3_len(diff);
			if (L > 1e-6f) {
				v3 n_ab = scale(diff, 1.0f / L);
				v3 r_a = sub(r_t1.point1, pos_A_t1);
				v3 vA_surf = add(bh->velocity, cross(bh->angular_velocity, r_a));
				float v_approach = dot(vA_surf, n_ab);
				if (v_approach <= 0.0f) return 1.0f;
			}
			return t1;
		}

		ToiSepFn fcn = toi_make_sep_fn(&simplex, &r_t1, rot_A_t1, &ga_t1, static_b);

		// Inner (pushback) loop: advance along the separation axis.
		float t2 = t_max;

		for (int push = 0; push < TOI_PUSHBACK_ITER_MAX; push++) {
			v3 pos_A_t2; quat rot_A_t2;
			toi_pose_at(bs, bh, dt, t2, &pos_A_t2, &rot_A_t2);
			GJK_Shape ga_t2 = toi_shape_at_pose(sa, pos_A_t2, rot_A_t2);

			// FindMinSep at t2: picks deepest feature at current axis and
			// updates fcn.idxA/idxB in-place. After this, Evaluate uses
			// the same fixed indices.
			float s2 = toi_find_min_sep(&fcn, &ga_t2, static_b, rot_A_t2);
			if (s2 > target + tol) return 1.0f;
			if (s2 > target - tol) { t1 = t2; break; }

			// Evaluate sep at t1 with the features FindMinSep just pinned.
			GJK_Shape ga_t1_ref = toi_shape_at_pose(sa, pos_A_t1, rot_A_t1);
			float s1 = toi_eval_sep(&fcn, &ga_t1_ref, static_b, rot_A_t1);
			if (s1 < target - tol) return t1;
			if (s1 <= target + tol) return t1;

			// Bisection + secant root find on eval_sep(...) - target in [t1, t2].
			float a1 = t1, a2 = t2;
			float sa1 = s1, sa2 = s2;
			float t_star = a1;
			for (int root = 0; root < TOI_ROOT_ITER_MAX; root++) {
				float t;
				if (root & 1) {
					float denom = sa2 - sa1;
					t = (fabsf(denom) > 1e-12f) ? a1 + (target - sa1) * (a2 - a1) / denom : 0.5f * (a1 + a2);
				} else {
					t = 0.5f * (a1 + a2);
				}
				v3 pos_A_t; quat rot_A_t;
				toi_pose_at(bs, bh, dt, t, &pos_A_t, &rot_A_t);
				GJK_Shape ga_t = toi_shape_at_pose(sa, pos_A_t, rot_A_t);
				float s = toi_eval_sep(&fcn, &ga_t, static_b, rot_A_t);
				if (fabsf(s - target) < tol) { t_star = t; break; }
				if (s > target) { a1 = t; sa1 = s; } else { a2 = t; sa2 = s; }
				t_star = t;
			}
			t2 = t_star;
			// Outer loop reseeds sep-fn if features changed.
		}
	}

	return t1 < 1.0f ? t1 : 1.0f;
}

// Body-vs-world: sweep every shape of body `bi` through the step against the
// static BVH, returning the earliest TOI fraction. Falls out early at t=0.
// Returns 1.0 if no hit. `pre_pos`/`pre_rot` must be the body's PRE-STEP pose.
static float toi_body_vs_static(WorldInternal* w, int bi, v3 pre_pos, quat pre_rot, float dt)
{
	BodyCold* bc = &w->body_cold[bi];
	BodyHot* bh = &w->body_hot[bi];
	BodyState pre_bs = (BodyState){ .position = pre_pos, .rotation = pre_rot };
	BodyState* bs = &pre_bs;
	float r_max = bc->r_max_body;
	int nshapes = asize(bc->shapes);
	if (nshapes == 0) return 1.0f;

	// Swept AABB: union of tight AABB at t=0 and at t=1.
	AABB start_aabb = body_aabb(bs, bc);
	v3 trans = scale(bh->velocity, dt);
	AABB end_aabb = (AABB){ add(start_aabb.min, trans), add(start_aabb.max, trans) };
	// Angular motion: inflate by r_max * |omega| * dt on all sides (worst case).
	float ang_pad = r_max * v3_len(bh->angular_velocity) * dt;
	v3 pad = V3(ang_pad, ang_pad, ang_pad);
	AABB swept = (AABB){
		.min = sub(V3(fminf(start_aabb.min.x, end_aabb.min.x), fminf(start_aabb.min.y, end_aabb.min.y), fminf(start_aabb.min.z, end_aabb.min.z)), pad),
		.max = add(V3(fmaxf(start_aabb.max.x, end_aabb.max.x), fmaxf(start_aabb.max.y, end_aabb.max.y), fmaxf(start_aabb.max.z, end_aabb.max.z)), pad),
	};

	// Query static BVH for candidate static bodies overlapping the swept AABB.
	CK_DYNA int* candidates = NULL;
	bvh_query_aabb(w->bvh_static, (AABB){ swept.min, swept.max }, &candidates);

	float earliest = 1.0f;
	int n_cand = asize(candidates);
	for (int c = 0; c < n_cand; c++) {
		int sbi = candidates[c];
		if (sbi == bi) continue; // shouldn't happen (static list), defensive
		if (!split_alive(w->body_gen, sbi)) continue;
		BodyCold* sbc = &w->body_cold[sbi];
		BodyState* sbs = &w->body_state[sbi];
		int sbshapes = asize(sbc->shapes);
		for (int sshi = 0; sshi < sbshapes; sshi++) {
			ShapeInternal* ss = &sbc->shapes[sshi];
			if (ss->type == SHAPE_MESH || ss->type == SHAPE_HEIGHTFIELD) {
				// Broad triangle sweep not implemented in this MVP — emit a
				// coarse BVH-AABB check and skip. Handled by speculative
				// contacts once within margin. Future: iterate candidate
				// triangles from the mesh/heightfield leaf.
				continue;
			}
			GJK_Shape gb = toi_shape_at_pose(ss, sbs->position, sbs->rotation);
			for (int ashi = 0; ashi < nshapes; ashi++) {
				ShapeInternal* sa = &bc->shapes[ashi];
				if (sa->type == SHAPE_MESH || sa->type == SHAPE_HEIGHTFIELD) continue;
				float t = toi_pair_advance(sa, bs, bh, r_max, dt, &gb);
				if (t < earliest) earliest = t;
				if (earliest <= 0.0f) { afree(candidates); return 0.0f; }
			}
		}
	}
	afree(candidates);
	return earliest;
}

// Public driver used by world_step (wired in CCD plan step 6). Advances
// every body in fast_body_list to its earliest static-world TOI (or through
// the whole step), writing back only pose. Velocity is left alone — the next
// frame's speculative pass applies the restitution replace-rule.
//
// Multi-pass: after the first clamp, re-sweep on the remaining (1-t*) budget
// up to TOI_SUBSTEP_MAX passes. Stops early when a pass consumes the whole
// remainder (t*=1) or hits the cap.
// Build a ShapeInternal from a public ShapeParams, as body_add_shape does,
// minus side effects (BVH insert etc.). Used by the cast API.
static ShapeInternal toi_shape_from_params(ShapeParams p)
{
	ShapeInternal s = {0};
	s.type = p.type;
	s.local_pos = V3(0, 0, 0);
	s.local_rot = quat_identity();
	s.rounding_radius = (p.type == SHAPE_BOX || p.type == SHAPE_HULL) && p.rounding_radius > 0.0f ? p.rounding_radius : 0.0f;
	switch (p.type) {
	case SHAPE_SPHERE:  s.sphere.radius = p.sphere.radius; break;
	case SHAPE_CAPSULE: s.capsule.half_height = p.capsule.half_height; s.capsule.radius = p.capsule.radius; break;
	case SHAPE_BOX:     s.box.half_extents = p.box.half_extents; break;
	case SHAPE_HULL:    s.hull.hull = p.hull.hull; s.hull.scale = p.hull.scale; break;
	default: break;
	}
	shape_compute_bounds(&s);
	return s;
}

// Shape-cast driver. Sweeps `sa` from `start_pos`/`start_rot` through
// `translation` (and optional angular velocity `omega`) against the static
// BVH. Returns the earliest TOI fraction in [0, 1]; 1.0 means no hit.
static float toi_shape_cast_world(WorldInternal* w, ShapeInternal* sa, v3 start_pos, quat start_rot, v3 translation, v3 omega, Body ignore_body, v3* out_point, v3* out_normal, int* out_body_idx, int* out_shape_idx)
{
	// Treat the sweep as one "step" of nominal dt = 1.0. Velocity = translation.
	float dt = 1.0f;
	BodyState pre = (BodyState){ .position = start_pos, .rotation = start_rot };
	BodyHot ghost = {0};
	ghost.velocity = translation;
	ghost.angular_velocity = omega;
	ghost.inv_mass = 1.0f;

	// Swept AABB for BVH query.
	// Shape AABB at start pose:
	v3 pmin = V3(1e18f, 1e18f, 1e18f), pmax = V3(-1e18f, -1e18f, -1e18f);
	// Cheap bound: sphere of radius r_max at start_pos.
	float r_max = sa->r_max;
	pmin = V3(start_pos.x - r_max, start_pos.y - r_max, start_pos.z - r_max);
	pmax = V3(start_pos.x + r_max, start_pos.y + r_max, start_pos.z + r_max);
	v3 end_pos = add(start_pos, translation);
	v3 emin = V3(end_pos.x - r_max, end_pos.y - r_max, end_pos.z - r_max);
	v3 emax = V3(end_pos.x + r_max, end_pos.y + r_max, end_pos.z + r_max);
	float ang_pad = r_max * v3_len(omega) * dt;
	v3 pad = V3(ang_pad, ang_pad, ang_pad);
	AABB swept = (AABB){
		.min = sub(V3(fminf(pmin.x, emin.x), fminf(pmin.y, emin.y), fminf(pmin.z, emin.z)), pad),
		.max = add(V3(fmaxf(pmax.x, emax.x), fmaxf(pmax.y, emax.y), fmaxf(pmax.z, emax.z)), pad),
	};

	CK_DYNA int* candidates = NULL;
	bvh_query_aabb(w->bvh_static, swept, &candidates);

	float earliest = 1.0f;
	int best_body = -1, best_shape = -1;
	v3 best_pt = V3(0,0,0), best_n = V3(0,1,0);
	int n_cand = asize(candidates);
	for (int c = 0; c < n_cand; c++) {
		int sbi = candidates[c];
		if (!split_alive(w->body_gen, sbi)) continue;
		if (ignore_body.id != 0 && split_valid(w->body_gen, ignore_body) && handle_index(ignore_body) == sbi) continue;
		BodyCold* sbc = &w->body_cold[sbi];
		BodyState* sbs = &w->body_state[sbi];
		int sbshapes = asize(sbc->shapes);
		for (int sshi = 0; sshi < sbshapes; sshi++) {
			ShapeInternal* ss = &sbc->shapes[sshi];
			if (ss->type == SHAPE_MESH || ss->type == SHAPE_HEIGHTFIELD) continue;
			GJK_Shape gb = toi_shape_at_pose(ss, sbs->position, sbs->rotation);
			float t = toi_pair_advance(sa, &pre, &ghost, sa->r_max, dt, &gb);
			if (t < earliest) {
				earliest = t;
				best_body = sbi;
				best_shape = sshi;
				// Contact point/normal: re-run GJK at the hit t to get closest-point pair.
				v3 hp; quat hq;
				toi_pose_at(&pre, &ghost, dt, t, &hp, &hq);
				GJK_Shape ga = toi_shape_at_pose(sa, hp, hq);
				GJK_Result r = gjk_distance(&ga, &gb, NULL);
				best_pt = r.point2; // on static side
				v3 diff = sub(r.point1, r.point2);
				float dlen = v3_len(diff);
				best_n = dlen > 1e-6f ? scale(diff, 1.0f / dlen) : V3(0, 1, 0);
			}
		}
	}
	afree(candidates);
	if (out_point)    *out_point    = best_pt;
	if (out_normal)   *out_normal   = best_n;
	if (out_body_idx) *out_body_idx = best_body;
	if (out_shape_idx)*out_shape_idx= best_shape;
	return earliest;
}

// Helper: lerp body pose from (pa, qa) to (pb, qb) by fraction t, with
// shortest-arc quaternion interp + renormalize.
static void toi_lerp_pose(v3 pa, quat qa, v3 pb, quat qb, float t, v3* out_pos, quat* out_rot)
{
	out_pos->x = pa.x + t * (pb.x - pa.x);
	out_pos->y = pa.y + t * (pb.y - pa.y);
	out_pos->z = pa.z + t * (pb.z - pa.z);
	float d = qa.x*qb.x + qa.y*qb.y + qa.z*qb.z + qa.w*qb.w;
	quat qbs = qb;
	if (d < 0.0f) { qbs.x = -qbs.x; qbs.y = -qbs.y; qbs.z = -qbs.z; qbs.w = -qbs.w; }
	quat q = { qa.x + t * (qbs.x - qa.x), qa.y + t * (qbs.y - qa.y), qa.z + t * (qbs.z - qa.z), qa.w + t * (qbs.w - qa.w) };
	*out_rot = quat_norm(q);
}

// Helper: is body `bi`'s COM at `pos` inside any static body's AABB?
// Used as a tunnel indicator: only "body center inside static volume"
// triggers a clamp, not just surface-tangent gravity overlap on a resting
// body. The rot parameter is unused (center check is pose-invariant) but
// kept for API symmetry with pose-carrying callers.
static int toi_center_in_static(WorldInternal* w, int bi, v3 pos, quat rot, AABB query)
{
	(void)rot;
	CK_DYNA int* cands = NULL;
	bvh_query_aabb(w->bvh_static, query, &cands);
	int inside = 0;
	for (int c = 0; c < asize(cands) && !inside; c++) {
		int sbi = cands[c];
		if (sbi == bi || !split_alive(w->body_gen, sbi)) continue;
		AABB sbbox = body_aabb(&w->body_state[sbi], &w->body_cold[sbi]);
		if (pos.x >= sbbox.min.x && pos.x <= sbbox.max.x
		 && pos.y >= sbbox.min.y && pos.y <= sbbox.max.y
		 && pos.z >= sbbox.min.z && pos.z <= sbbox.max.z) {
			inside = 1;
		}
	}
	afree(cands);
	return inside;
}

// Per-body TOI step. Pure function of (pre-pose, post-pose, static BVH).
// No shared writes except to this body's own pose, so calls across bodies
// are trivially parallel.
static void toi_advance_one_body(WorldInternal* w, int k, float dt)
{
	int bi = w->fast_body_list[k];
	if (!split_alive(w->body_gen, bi)) return;
	if (w->body_hot[bi].inv_mass == 0.0f) return;

	BodyState* bs = &w->body_state[bi];
	BodyHot*   bh = &w->body_hot[bi];
	v3   p0 = w->fast_body_pre_pos[k];
	quat q0 = w->fast_body_pre_rot[k];

	float earliest = toi_body_vs_static(w, bi, p0, q0, dt);
	if (earliest >= 1.0f) {
		// TOI sweep reported no hit. The sweep uses straight-line integration
		// from pre_pos with post-solve velocity. Multi-substep integration can
		// produce a curved path whose endpoint lands INSIDE a static body even
		// though the straight-line sweep missed. Sanity check: GJK the actual
		// post-solve pose against static candidates; if any overlap, clamp to
		// the last safe pose along the lerp from pre to post.
		//
		BodyCold* bc = &w->body_cold[bi];
		int nshapes = asize(bc->shapes); (void)nshapes;
		AABB body_aabb_now = body_aabb(bs, bc);
		CK_DYNA int* cands = NULL;
		bvh_query_aabb(w->bvh_static, body_aabb_now, &cands);
		// Tunnel detection: body CENTER either ends inside a static AABB, or
		// the pre→post CENTER path crosses fully through a static AABB (body
		// tunneled past a thin static). Sample 8 lerp points along the
		// center trajectory and check center-in-AABB at each; also test
		// whether the segment pre_pos→post_pos passes through any candidate
		// static's AABB (Liang-Barsky slab clip).
		int overlapping = 0;
		v3 post_pos_tmp = bs->position;
		for (int c = 0; c < asize(cands) && !overlapping; c++) {
			int sbi = cands[c];
			if (sbi == bi || !split_alive(w->body_gen, sbi)) continue;
			AABB sbbox = body_aabb(&w->body_state[sbi], &w->body_cold[sbi]);
			// Liang-Barsky: segment p0 + t*(post-p0), find t-range inside AABB.
			v3 d = sub(post_pos_tmp, p0);
			float tmin = 0.0f, tmax = 1.0f;
			float p_origin[3] = { p0.x, p0.y, p0.z };
			float dd[3] = { d.x, d.y, d.z };
			float bmin[3] = { sbbox.min.x, sbbox.min.y, sbbox.min.z };
			float bmax[3] = { sbbox.max.x, sbbox.max.y, sbbox.max.z };
			int hit = 1;
			for (int a = 0; a < 3 && hit; a++) {
				if (fabsf(dd[a]) < 1e-12f) {
					if (p_origin[a] < bmin[a] || p_origin[a] > bmax[a]) hit = 0;
				} else {
					float t1 = (bmin[a] - p_origin[a]) / dd[a];
					float t2 = (bmax[a] - p_origin[a]) / dd[a];
					if (t1 > t2) { float tmp = t1; t1 = t2; t2 = tmp; }
					if (t1 > tmin) tmin = t1;
					if (t2 < tmax) tmax = t2;
					if (tmin > tmax) hit = 0;
				}
			}
			if (hit) overlapping = 1;
		}
		afree(cands);
		if (overlapping) {
			// Post-solve pose overlaps static. If pre-step pose is also
			// overlapping (body already embedded), keep solver's attempt.
			// Otherwise, binary search between pre and post to find the last
			// non-overlapping lerp fraction and clamp there — this progresses
			// the body further than pre_pos (avoiding the "frozen against
			// surface" oscillation where shallow gravity overlap reverts every
			// frame) while preventing deep tunnel.
			v3 post_pos = bs->position;
			quat post_rot = bs->rotation;
			int pre_overlaps = toi_center_in_static(w, bi, p0, q0, body_aabb_now);
			if (!pre_overlaps) {
				float lo = 0.0f, hi = 1.0f;
				for (int iter = 0; iter < 8; iter++) {
					float mid = 0.5f * (lo + hi);
					v3 tp; quat tr;
					toi_lerp_pose(p0, q0, post_pos, post_rot, mid, &tp, &tr);
					if (toi_center_in_static(w, bi, tp, tr, body_aabb_now)) hi = mid;
					else lo = mid;
				}
				toi_lerp_pose(p0, q0, post_pos, post_rot, lo, &bs->position, &bs->rotation);
			}
		}
		afree(cands);
		return;
	}

	// Pull back a hair so we leave a small gap for next-frame speculative.
	float pullback = TOI_SLOP_FRACTION;
	if (earliest > pullback) earliest -= pullback;
	else earliest = 0.0f;

	// Clamp pose to fraction `earliest` of the step from the pre-step pose.
	// Velocity is intentionally left alone: next frame's speculative pass
	// applies restitution via the replace-rule.
	BodyState pre = (BodyState){ .position = p0, .rotation = q0 };
	toi_pose_at(&pre, bh, dt, earliest, &bs->position, &bs->rotation);

	// Post-clamp safety: if the clamped pose still overlaps any static,
	// binary-search along the lerp from pre_pos → clamped_pose for the
	// largest non-overlapping fraction. Catches compound-body TOI where one
	// shape's first-contact fraction pushes another already-near-contact
	// shape deeper into a static. If the lerp has no safe fraction (pre_pos
	// also overlaps; body is embedded and rotating into/out of contact),
	// keep the solver's best-effort clamped pose.
	int post_clamp_fired = 0;
	{
		BodyCold* bc = &w->body_cold[bi];
		AABB body_aabb_clamped = body_aabb(bs, bc);
		int clamped_in = toi_center_in_static(w, bi, bs->position, bs->rotation, body_aabb_clamped);
		if (clamped_in) {
			v3 clamped_pos = bs->position;
			quat clamped_rot = bs->rotation;
			float lo = -1.0f, hi = 1.0f;
			// Sample 5 points along the lerp to find a safe seed.
			for (int i = 0; i <= 4; i++) {
				float u = (float)i / 4.0f;
				v3 tp; quat tr;
				toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, u, &tp, &tr);
				if (!toi_center_in_static(w, bi, tp, tr, body_aabb_clamped)) { lo = u; break; }
			}
			if (lo >= 0.0f) {
				// Binary search between safe lo and the first overlap above it.
				for (int iter = 0; iter < 8; iter++) {
					float mid = 0.5f * (lo + hi);
					v3 tp; quat tr;
					toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, mid, &tp, &tr);
					if (toi_center_in_static(w, bi, tp, tr, body_aabb_clamped)) hi = mid;
					else lo = mid;
				}
				toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, lo, &bs->position, &bs->rotation);
				post_clamp_fired = 1;
			}
		}
	}

	// When earliest is effectively 0 (body already at boundary at pre-step)
	// AND the body would tunnel more than its inscribed sphere in one step,
	// narrowphase may have failed to emit a boundary contact at the clamped
	// pose (hull-vs-hull SAT can reject near-exactly-touching pairs). Pull
	// the body back along the sweep direction until GJK shows a visible gap
	// matching the speculative margin, so next frame's narrowphase sees the
	// contact. Pulling along -velocity (not along GJK normal) preserves
	// lateral motion at corners/edges where a single normal is ambiguous.
	float v_norm = v3_len(bh->velocity);
	float r_min = w->body_cold[bi].r_min_body;
	int body_is_fast = r_min > 0.0f && v_norm * dt > 2.0f * r_min;
	if (earliest <= TOI_SLOP_FRACTION + 1e-6f && body_is_fast) {
		BodyCold* bc = &w->body_cold[bi];
		AABB body_aabb_now = body_aabb(bs, bc);
		CK_DYNA int* cands = NULL;
		bvh_query_aabb(w->bvh_static, aabb_expand(body_aabb_now, w->speculative_margin), &cands);
		float min_dist = FLT_MAX;
		for (int c = 0; c < asize(cands); c++) {
			int sbi = cands[c];
			if (sbi == bi || !split_alive(w->body_gen, sbi)) continue;
			BodyCold* sbc = &w->body_cold[sbi];
			BodyState* sbs = &w->body_state[sbi];
			for (int sshi = 0; sshi < asize(sbc->shapes); sshi++) {
				ShapeInternal* ss = &sbc->shapes[sshi];
				if (ss->type == SHAPE_MESH || ss->type == SHAPE_HEIGHTFIELD) continue;
				GJK_Shape gb = toi_shape_at_pose(ss, sbs->position, sbs->rotation);
				for (int ashi = 0; ashi < asize(bc->shapes); ashi++) {
					ShapeInternal* sa = &bc->shapes[ashi];
					if (sa->type == SHAPE_MESH || sa->type == SHAPE_HEIGHTFIELD) continue;
					GJK_Shape ga = toi_shape_at_pose(sa, bs->position, bs->rotation);
					GJK_Result r = gjk_distance(&ga, &gb, NULL);
					if (r.distance < min_dist) min_dist = r.distance;
				}
			}
		}
		afree(cands);
		if (min_dist < w->speculative_margin) {
			float backstep = (w->speculative_margin - min_dist);
			if (v_norm > 1e-6f) bs->position = sub(bs->position, scale(bh->velocity, backstep / v_norm));
		}
	}

	// Multi-pass: after clamping, re-sweep the remaining (1 - earliest) of
	// the step so the body can slide along the struck surface or hit a
	// second primitive. Capped by TOI_SUBSTEP_MAX. Skip when the post-clamp
	// safety binary-searched us back — additional forward sweeps would
	// undo that safety by re-integrating INTO the static we just escaped.
	if (post_clamp_fired) return;
	float t_consumed = earliest;
	for (int pass = 1; pass < TOI_SUBSTEP_MAX; pass++) {
		float remaining = 1.0f - t_consumed;
		if (remaining <= 1e-4f) break;
		v3   p_now = bs->position;
		quat q_now = bs->rotation;
		float t_rel = toi_body_vs_static(w, bi, p_now, q_now, dt * remaining);
		if (t_rel >= 1.0f) break;
		if (t_rel > pullback) t_rel -= pullback;
		else t_rel = 0.0f;
		BodyState pre2 = (BodyState){ .position = p_now, .rotation = q_now };
		toi_pose_at(&pre2, bh, dt * remaining, t_rel, &bs->position, &bs->rotation);
		t_consumed += t_rel * remaining;
		if (t_rel <= 0.0f) break;
	}
}

// Parallel-for worker: advance a contiguous range of fast_body_list entries.
typedef struct ToiCtx { WorldInternal* w; float dt; } ToiCtx;
static void toi_work_fn(void* ctx, int start, int count)
{
	ToiCtx* tc = (ToiCtx*)ctx;
	for (int k = start; k < start + count; k++) toi_advance_one_body(tc->w, k, tc->dt);
}

// Public driver used by world_step (wired in CCD plan step 6). Dispatches
// per-body TOI across the thread pool. Each task reads the static BVH
// (read-only) and writes only its own body's pose, so the whole pass is
// embarrassingly parallel.
static void toi_advance_fast_bodies(WorldInternal* w, float dt)
{
	int n = asize(w->fast_body_list);
	if (n == 0) return;
	int n_workers = w->thread_count > 0 ? w->thread_count : 1;
	if (n_workers > 1 && n >= 2) {
		ToiCtx tc = { .w = w, .dt = dt };
		pool_dispatch(toi_work_fn, &tc, n, 4, n_workers);
	} else {
		for (int k = 0; k < n; k++) toi_advance_one_body(w, k, dt);
	}
}
