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
		return gjk_box(world_pos, world_rot, s->box.half_extents);
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

enum { TOI_POINTS = 0, TOI_FACE_A = 1, TOI_FACE_B = 2, TOI_EDGE_EDGE = 3 };

typedef struct ToiSepFn
{
	int type;
	// `axis_local` is a unit direction.
	//   POINTS / FACE_A: stored in A's local frame (rotates with A over t).
	//   FACE_B:          stored in world (B static, axis doesn't move).
	//   EDGE_EDGE:       stored in A's local frame. Derived from seed-time
	//                    normalize(eA_world x eB_world) and sign-locked so
	//                    positive means B is on the + side of A's edge
	//                    plane. Rotated to world via rot_A at each sample
	//                    — never recomputed via cross product, which is
	//                    what produced the sign flip at parallel edges.
	v3 axis_local;
	// Reference support points. POINTS: A-local witness on A, B-world witness
	// on B. FACE_A: A-local reference point on the face. FACE_B: B-world
	// reference point on the face. EDGE_EDGE: A-local edge endpoints in
	// edge_a_loc0/edge_a_loc1; B-world endpoints in edge_b_w0/edge_b_w1.
	v3 local_point_a;
	v3 world_point_b;
	// Edge-edge extras.
	v3 edge_a_loc0, edge_a_loc1;
	v3 edge_b_w0,   edge_b_w1;
	int axis_sign; // +1 or -1 for edge-edge axis sign lock; 0 otherwise
} ToiSepFn;

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

// Build a separation function from the GJK simplex at the current pose.
// Feature type is chosen from the unique-support counts on each side:
//   kA=1, kB=1:              POINTS
//   kA>=2, kB==1:            FACE_A
//   kA==1, kB>=2:            FACE_B
//   kA==2, kB==2 (3D):       EDGE_EDGE (with parallel-edge fallback to FACE_A)
//   kA>=3 or kB>=3:          FACE on the side with more supports
// The axis is derived from the GJK closest-point direction (coincides with
// the relevant face normal when a face feature dominates, and with the
// edge-cross-product direction for edge-edge when non-parallel).
static ToiSepFn toi_make_sep_fn(const GJK_Simplex* simplex, const GJK_Result* r0, v3 pos_A, quat rot_A)
{
	ToiSepFn f = {0};
	int kA = 0, kB = 0;
	int fa[4] = {0}, fb[4] = {0};
	toi_simplex_feature_counts(simplex, &kA, &kB, fa, fb);

	v3 diff = sub(r0->point2, r0->point1);
	float dlen = v3_len(diff);
	v3 axis_world = dlen > 1e-8f ? scale(diff, 1.0f / dlen) : V3(0, 1, 0);

	quat rot_A_inv = quat_inv(rot_A);
	v3 local_pA = rotate(rot_A_inv, sub(r0->point1, pos_A));

	// Edge-edge: two distinct supports on each side. Axis from (eA x eB).
	// Guard: when the GJK closest-point direction (dlen > tol) is a credible
	// separating axis, prefer it — cross-product axis is only correct for true
	// edge-vs-edge contact. Capsule-vs-box-face routinely returns kA=2 kB=2
	// (two capsule endpoints + two face corners) but the correct axis is the
	// face normal, not cross(capsule_axis, face_edge). dlen > tol means GJK
	// found a well-defined direction; trust it over the cross product.
	if (kA == 2 && kB == 2 && dlen <= 1e-8f) {
		// Recover the simplex vertices that produced each side's edge.
		v3 a0_w = V3(0,0,0), a1_w = V3(0,0,0), b0_w = V3(0,0,0), b1_w = V3(0,0,0);
		int nA = 0, nB = 0;
		for (int i = 0; i < simplex->count && (nA < 2 || nB < 2); i++) {
			if (nA < 2) {
				int dup = 0;
				for (int j = 0; j < nA; j++) {
					v3 w = j == 0 ? a0_w : a1_w;
					if (len2(sub(w, simplex->v[i].point1)) < 1e-12f) { dup = 1; break; }
				}
				if (!dup) {
					if (nA == 0) a0_w = simplex->v[i].point1;
					else         a1_w = simplex->v[i].point1;
					nA++;
				}
			}
			if (nB < 2) {
				int dup = 0;
				for (int j = 0; j < nB; j++) {
					v3 w = j == 0 ? b0_w : b1_w;
					if (len2(sub(w, simplex->v[i].point2)) < 1e-12f) { dup = 1; break; }
				}
				if (!dup) {
					if (nB == 0) b0_w = simplex->v[i].point2;
					else         b1_w = simplex->v[i].point2;
					nB++;
				}
			}
		}
		v3 eA_w = sub(a1_w, a0_w);
		v3 eB_w = sub(b1_w, b0_w);
		v3 n_w = cross(eA_w, eB_w);
		float nl2 = len2(n_w);
		if (nl2 > 1e-10f) {
			// Axis captured at SEED from a known non-parallel (eA × eB) and
			// then STORED IN A'S LOCAL FRAME so it rotates with A. The
			// scalar-triple-product bug came from *recomputing* the cross
			// product at every sample — crossing through the parallel locus
			// sends |cross| → 0 with an indeterminate direction, which is
			// what produced the spurious pos↔neg flips mid-step. Recording
			// the seed axis in A-local sidesteps that: axis_world(t) =
			// rotate(rot_A(t), axis_local) is always unit-length, always
			// smooth, and tracks A's edge orientation (one half of the
			// cross product) correctly through rotation. Small rotations
			// within a step match "cross product with A's current edge"
			// almost exactly; large rotations don't matter because the
			// feature will be reseeded by the outer loop before then.
			f.type = TOI_EDGE_EDGE;
			v3 axis_world = scale(n_w, 1.0f / sqrtf(nl2));
			// Sign-lock so positive = B on the + side of A's edge plane at
			// seed. A sign change in the sep function now corresponds to an
			// actual edge-plane crossing, not to an axis direction flip.
			v3 pA_seed = r0->point1, pB_seed = r0->point2;
			if (dot(sub(pB_seed, pA_seed), axis_world) < 0.0f) axis_world = neg(axis_world);
			f.axis_local = rotate(rot_A_inv, axis_world); // A-local
			f.axis_sign = 0;            // unused; kept for struct layout
			f.edge_a_loc0 = rotate(rot_A_inv, sub(a0_w, pos_A));
			f.edge_a_loc1 = rotate(rot_A_inv, sub(a1_w, pos_A));
			f.edge_b_w0 = b0_w;
			f.edge_b_w1 = b1_w;
			f.local_point_a = local_pA;
			f.world_point_b = pB_seed;
			return f;
		}
		// Fall through if edges are parallel: FACE_A handles it.
	}

	if (kB >= 2 && kA == 1) {
		f.type = TOI_FACE_B;
		f.axis_local = axis_world;
		f.world_point_b = r0->point2;
		f.local_point_a = local_pA;
		return f;
	}
	if (kA >= 2) {
		f.type = TOI_FACE_A;
		f.axis_local = rotate(rot_A_inv, axis_world);
		f.local_point_a = local_pA;
		f.world_point_b = r0->point2;
		return f;
	}
	f.type = TOI_POINTS;
	f.local_point_a = local_pA;
	f.world_point_b = r0->point2;
	// Axis must be SEED-anchored, not re-derived from live witness points.
	// If we recomputed axis = normalize(pB - pA(t)) each sample, the
	// direction flips when A crosses past B's seed side, making the root
	// finder think the shapes have already separated (positive axis-dot
	// on both sides of the crossing) and return a false miss. Store the
	// seed-time direction in A's local frame so it rotates correctly with
	// A's pose but never flips purely from translation.
	f.axis_local = rotate(rot_A_inv, axis_world);
	return f;
}

// Compute the world-frame axis of the separation function at time t, given
// A's integrated pose.
static v3 toi_axis_world(const ToiSepFn* f, v3 pos_A, quat rot_A)
{
	switch (f->type) {
	case TOI_FACE_A: return rotate(rot_A, f->axis_local);
	case TOI_FACE_B: return f->axis_local; // world-frame, B static
	case TOI_POINTS: return rotate(rot_A, f->axis_local); // seed direction in A's frame
	case TOI_EDGE_EDGE:
		// Seed axis stored in A's local frame — rotates with A, never
		// recomputed via cross product at sample time. See make_sep_fn.
		return rotate(rot_A, f->axis_local);
	default:
		// Unreachable: all ToiSepType values are covered above. The POINTS
		// case is handled as `case TOI_POINTS` earlier.
		return V3(0, 1, 0);
	}
}

// Find the witness pair giving minimum separation at time t. Returns the
// minimum separation value; writes the supporting world points (on A and B)
// into *out_pa and *out_pb. Box2D's b2FindMinSeparation analog — but since
// we operate on GJK_Shape (which internally handles any shape type), we
// store world support points directly rather than vertex indices.
// Compute the edge-edge scalar-triple-product separation at a given A-pose.
// s = ((eA(t) x eB) . (pB - pA(t))) / |eA(t) x eB|  — the signed distance
// between the two infinite lines. Smooth through the parallel-edge locus
// because both numerator and denominator vanish together (the singularity
// Edge-edge separation function. Signed projection of the current closest-
// point pair onto a SEED-LOCKED world axis — the cross product eA × eB
// measured at make_sep_fn time and sign-adjusted so positive means B is on
// the + side of A's edge plane.
//
// Design goals (requested feature):
//   * Smooth: the closest-points pair (qA, qB) is a continuous function of
//     A's pose (pose → clamped segment-segment solve); projecting onto a
//     fixed axis stays smooth.
//   * Convex-ish: for linear A-motion with fixed A-rotation, the projection
//     is affine in t (qA moves linearly, axis is fixed). Rotation makes it
//     locally smooth rather than affine, but for one TOI step the rotation
//     is small.
//   * No parallel-axis singularity: axis is recorded ONCE at seed from a
//     non-parallel cross product; we never recompute it at sample time, so
//     near-parallel sample poses don't divide by zero.
//   * Pos-to-neg sign change ONLY at actual edge-plane crossings: sign
//     flips when (qB - qA) . n̂_seed changes sign, which happens when A's
//     edge crosses the plane containing B's edge perpendicular to n̂_seed.
//     At that crossing the segments are in the same plane — geometrically
//     the "intersection event" the TOI root finder is meant to locate.
//   * Root finder gets a real sign change to bracket, not a U-shape.
//
// The old scalar-triple-product form was the same projection but with the
// axis RE-DERIVED at each sample, which reintroduced the parallel singularity
// and produced spurious sign flips when A swept across the plane mid-step.
static float toi_edge_edge_sep(const ToiSepFn* f, v3 pos_A, quat rot_A)
{
	v3 a0 = add(pos_A, rotate(rot_A, f->edge_a_loc0));
	v3 a1 = add(pos_A, rotate(rot_A, f->edge_a_loc1));
	v3 b0 = f->edge_b_w0;
	v3 b1 = f->edge_b_w1;
	v3 dA = sub(a1, a0);
	v3 dB = sub(b1, b0);
	v3 r  = sub(a0, b0);
	float a = dot(dA, dA);
	float c = dot(dB, dB);
	float b = dot(dA, dB);
	float e = dot(dA, r);
	float fv = dot(dB, r);
	const float EPS = 1e-10f;
	float s, t;
	if (a < EPS && c < EPS) { s = 0.0f; t = 0.0f; }
	else if (a < EPS) { s = 0.0f; t = fmaxf(0.0f, fminf(1.0f, fv / c)); }
	else if (c < EPS) { t = 0.0f; s = fmaxf(0.0f, fminf(1.0f, -e / a)); }
	else {
		float denom = a * c - b * b;
		s = denom > EPS ? fmaxf(0.0f, fminf(1.0f, (b * fv - c * e) / denom)) : 0.0f;
		t = (b * s + fv) / c;
		if (t < 0.0f) { t = 0.0f; s = fmaxf(0.0f, fminf(1.0f, -e / a)); }
		else if (t > 1.0f) { t = 1.0f; s = fmaxf(0.0f, fminf(1.0f, (b - e) / a)); }
	}
	v3 qA = add(a0, scale(dA, s));
	v3 qB = add(b0, scale(dB, t));
	v3 axis_world = rotate(rot_A, f->axis_local);
	return dot(sub(qB, qA), axis_world);
}

static float toi_find_min_sep(const ToiSepFn* f, ShapeInternal* sa, v3 pos_A, quat rot_A, GJK_Shape* gjk_b, v3* out_pa, v3* out_pb)
{
	GJK_Shape ga = toi_shape_at_pose(sa, pos_A, rot_A);
	v3 axis_world = toi_axis_world(f, pos_A, rot_A);

	// Support on A along +axis (closest to B), on B along -axis (closest
	// to A). Recorded for the root-finder's EvaluateSeparation path.
	int featA;
	v3 pa;
	gjk_support(&ga, axis_world, &featA, pa);
	int featB;
	v3 pb;
	gjk_support(gjk_b, neg(axis_world), &featB, pb);
	*out_pa = pa;
	*out_pb = pb;

	// For edge-edge, replace the axis-dot separation with the scalar
	// triple product (divided by the cross-product magnitude). This is
	// the signed line-to-line distance and stays smooth across the
	// parallel-edges singularity — the failure mode that plagues a
	// cross-product + sign-lock axis.
	if (f->type == TOI_EDGE_EDGE) return toi_edge_edge_sep(f, pos_A, rot_A);

	float rA = ga.radius, rB = gjk_b->radius;
	v3 pa_eff = add(pa, scale(axis_world, rA));
	v3 pb_eff = sub(pb, scale(axis_world, rB));
	return dot(sub(pb_eff, pa_eff), axis_world);
}

// Evaluate separation at time t using the witness points cached from a prior
// toi_find_min_sep call. For POINTS/FACE_A this uses local-frame points
// rotated to world by A's pose at t; for FACE_B it uses the fixed world
// point on B. Analog of b2EvaluateSeparation.
static float toi_eval_sep(const ToiSepFn* f, v3 pa_world_ref_local_A, v3 pb_world, float rA, float rB, v3 pos_A, quat rot_A)
{
	// Edge-edge uses the scalar triple product directly (depends only on
	// the stored edge anchors + A's pose, not on the cached witness
	// points). Match the FindMinSep formula so the root finder is
	// consistent.
	if (f->type == TOI_EDGE_EDGE) return toi_edge_edge_sep(f, pos_A, rot_A);
	v3 pa_world = add(pos_A, rotate(rot_A, pa_world_ref_local_A));
	v3 axis_world = toi_axis_world(f, pos_A, rot_A);
	v3 pa_eff = add(pa_world, scale(axis_world, rA));
	v3 pb_eff = sub(pb_world, scale(axis_world, rB));
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

		ToiSepFn fcn = toi_make_sep_fn(&simplex, &r_t1, pos_A_t1, rot_A_t1);

		// Inner (pushback) loop: advance along the separation axis.
		int done = 0;
		float t2 = t_max;
		float rA = ga_t1.radius;
		float rB = static_b->radius;

		for (int push = 0; push < TOI_PUSHBACK_ITER_MAX; push++) {
			v3 pos_A_t2; quat rot_A_t2;
			toi_pose_at(bs, bh, dt, t2, &pos_A_t2, &rot_A_t2);
			v3 pa2, pb2;
			float s2 = toi_find_min_sep(&fcn, sa, pos_A_t2, rot_A_t2, static_b, &pa2, &pb2);
	if (s2 > target + tol) {
				// Sep-fn at t2 claims "separated through the step". For the
				// edge-edge scalar-triple-product, the sign convention flips
				// when the body sweeps across the plane spanned by the two
				// edges — a large moving body at a glancing angle can produce
				// a positive s2 even though it tunnels deep through the
				// feature. If we're already close enough at t1 that a real
				// contact is plausible, trust GJK's distance over sep-fn's
				// sign and accept t1 as the TOI.
				if (fcn.type == TOI_EDGE_EDGE && r_t1.distance < LINEAR_SLOP * 4.0f) return t1;
				return 1.0f;
			}
			if (s2 > target - tol) { t1 = t2; break; }
			// Convert world support points back into A's local frame for
			// root-finder-stable re-evaluation. B stays in world (static).
			v3 pa_local_A = rotate(quat_inv(rot_A_t2), sub(pa2, pos_A_t2));

			float s1 = toi_eval_sep(&fcn, pa_local_A, pb2, rA, rB, pos_A_t1, rot_A_t1);
			if (s1 < target - tol) {
				// Initial overlap at t1 with this witness pair — accept t1 as TOI.
				return t1;
			}
			if (s1 <= target + tol) {
				return t1; // touching at t1
			}

			// Bisection + secant root find on eval_sep(...) - target in [t1, t2].
			float a1 = t1, a2 = t2;
			float sa1 = s1, sa2 = s2;
			float t_star = a1;
			for (int root = 0; root < TOI_ROOT_ITER_MAX; root++) {
				float t;
				if (root & 1) {
					// Secant rule
					float denom = sa2 - sa1;
					t = (fabsf(denom) > 1e-12f) ? a1 + (target - sa1) * (a2 - a1) / denom : 0.5f * (a1 + a2);
				} else {
					// Bisection guarantees progress
					t = 0.5f * (a1 + a2);
				}
				v3 pos_A_t; quat rot_A_t;
				toi_pose_at(bs, bh, dt, t, &pos_A_t, &rot_A_t);
				float s = toi_eval_sep(&fcn, pa_local_A, pb2, rA, rB, pos_A_t, rot_A_t);
				if (fabsf(s - target) < tol) { t_star = t; break; }
				if (s > target) { a1 = t; sa1 = s; } else { a2 = t; sa2 = s; }
				t_star = t;
			}
			t2 = t_star;
			// Loop: the shrunk [t1, t2] now brackets the root under the fixed
			// witness pair; re-test if the feature pair changed.
		}

		if (done) break;
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

// Helper: does body `bi` at pose (pos, rot) overlap ANY static body reachable
// from the given query AABB? Used by the post-TOI safety clamp to decide
// whether to binary-search for a safe pose.
static int toi_aabb_overlaps_any_static(WorldInternal* w, int bi, v3 pos, quat rot, AABB query)
{
	BodyCold* bc = &w->body_cold[bi];
	int nshapes = asize(bc->shapes);
	CK_DYNA int* cands = NULL;
	bvh_query_aabb(w->bvh_static, query, &cands);
	int overlapping = 0;
	for (int c = 0; c < asize(cands) && !overlapping; c++) {
		int sbi = cands[c];
		if (sbi == bi || !split_alive(w->body_gen, sbi)) continue;
		BodyCold* sbc = &w->body_cold[sbi];
		BodyState* sbs = &w->body_state[sbi];
		for (int sshi = 0; sshi < asize(sbc->shapes) && !overlapping; sshi++) {
			ShapeInternal* ss = &sbc->shapes[sshi];
			if (ss->type == SHAPE_MESH || ss->type == SHAPE_HEIGHTFIELD) continue;
			GJK_Shape gb = toi_shape_at_pose(ss, sbs->position, sbs->rotation);
			for (int ashi = 0; ashi < nshapes && !overlapping; ashi++) {
				ShapeInternal* sa = &bc->shapes[ashi];
				if (sa->type == SHAPE_MESH || sa->type == SHAPE_HEIGHTFIELD) continue;
				GJK_Shape ga = toi_shape_at_pose(sa, pos, rot);
				if (gjk_distance(&ga, &gb, NULL).distance <= 0.0f) overlapping = 1;
			}
		}
	}
	afree(cands);
	return overlapping;
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
		// Only run this check for bodies moving fast enough that the curved
		// substep path can genuinely diverge from the linear sweep (> 5x the
		// inscribed sphere in one frame). Resting/slow bodies have vanishingly
		// small curvature in their path and the check would otherwise fight
		// with the solver's speculative-contact penetration handling (shallow
		// gravity-induced overlap at each step gets erroneously reverted,
		// freezing the body against the resting surface).
		BodyCold* bc = &w->body_cold[bi];
		float v_norm_pre = v3_len(w->body_hot[bi].velocity);
		float r_min_pre = bc->r_min_body;
		// r_min == 0 means compound with no inscribed sphere around COM —
		// always run overlap check (the body has no self-contained safety zone
		// and can tunnel at any speed). For positive r_min, gate on speed.
		if (r_min_pre > 0.0f && v_norm_pre * dt < 5.0f * r_min_pre) return;
		int nshapes = asize(bc->shapes);
		AABB body_aabb_now = body_aabb(bs, bc);
		CK_DYNA int* cands = NULL;
		bvh_query_aabb(w->bvh_static, body_aabb_now, &cands);
		int overlapping = 0;
		for (int c = 0; c < asize(cands) && !overlapping; c++) {
			int sbi = cands[c];
			if (sbi == bi || !split_alive(w->body_gen, sbi)) continue;
			BodyCold* sbc = &w->body_cold[sbi];
			BodyState* sbs = &w->body_state[sbi];
			for (int sshi = 0; sshi < asize(sbc->shapes) && !overlapping; sshi++) {
				ShapeInternal* ss = &sbc->shapes[sshi];
				if (ss->type == SHAPE_MESH || ss->type == SHAPE_HEIGHTFIELD) continue;
				GJK_Shape gb = toi_shape_at_pose(ss, sbs->position, sbs->rotation);
				for (int ashi = 0; ashi < nshapes && !overlapping; ashi++) {
					ShapeInternal* sa = &bc->shapes[ashi];
					if (sa->type == SHAPE_MESH || sa->type == SHAPE_HEIGHTFIELD) continue;
					GJK_Shape ga = toi_shape_at_pose(sa, bs->position, bs->rotation);
					GJK_Result r = gjk_distance(&ga, &gb, NULL);
					if (r.distance <= 0.0f) overlapping = 1;
				}
			}
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
			int pre_overlaps = toi_aabb_overlaps_any_static(w, bi, p0, q0, body_aabb_now);
			if (!pre_overlaps) {
				float lo = 0.0f, hi = 1.0f;
				for (int iter = 0; iter < 8; iter++) {
					float mid = 0.5f * (lo + hi);
					v3 tp; quat tr;
					toi_lerp_pose(p0, q0, post_pos, post_rot, mid, &tp, &tr);
					if (toi_aabb_overlaps_any_static(w, bi, tp, tr, body_aabb_now)) hi = mid;
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
	{
		BodyCold* bc = &w->body_cold[bi];
		AABB body_aabb_clamped = body_aabb(bs, bc);
		if (toi_aabb_overlaps_any_static(w, bi, bs->position, bs->rotation, body_aabb_clamped)) {
			v3 clamped_pos = bs->position;
			quat clamped_rot = bs->rotation;
			float lo = -1.0f, hi = 1.0f;
			// Sample 5 points along the lerp to find a safe seed.
			for (int i = 0; i <= 4; i++) {
				float u = (float)i / 4.0f;
				v3 tp; quat tr;
				toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, u, &tp, &tr);
				if (!toi_aabb_overlaps_any_static(w, bi, tp, tr, body_aabb_clamped)) { lo = u; break; }
			}
			if (lo >= 0.0f) {
				// Binary search between safe lo and the first overlap above it.
				for (int iter = 0; iter < 8; iter++) {
					float mid = 0.5f * (lo + hi);
					v3 tp; quat tr;
					toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, mid, &tp, &tr);
					if (toi_aabb_overlaps_any_static(w, bi, tp, tr, body_aabb_clamped)) hi = mid;
					else lo = mid;
				}
				toi_lerp_pose(p0, q0, clamped_pos, clamped_rot, lo, &bs->position, &bs->rotation);
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
	// second primitive. Capped by TOI_SUBSTEP_MAX.
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
