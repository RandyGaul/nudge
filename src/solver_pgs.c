// solver.c -- contact constraint solver

// Analysis: Extremely well-studied and robust solver. Very simple, converges predictably. Has
// trouble propogating down longer chains, and dealing with large mass ratios. However, it's just
// so modular and well-understood it often becomes useful as a supporting algorithm in other more
// advanced solvers, or hybrid approaches.

// body_pair_key defined in broadphase.c (included before this file)

static void contact_tangent_basis(v3 n, v3* t1, v3* t2)
{
	if (fabsf(n.x) >= 0.57735f)
		*t1 = norm(V3(n.y, -n.x, 0.0f));
	else
		*t1 = norm(V3(0.0f, n.z, -n.y));
	*t2 = cross(n, *t1);
}

#define PATCH_MIN_AREA 0.001f

// Estimate contact patch area from manifold points projected onto contact plane.
// Uses dot(n, cross(a,b)) instead of len(cross(a,b)) — avoids sqrtf per triangle.
static float estimate_patch_area(Contact* contacts, int count)
{
	if (count < 3) return PATCH_MIN_AREA;
	v3 n = contacts[0].normal;
	float area = 0.0f;
	for (int i = 1; i < count - 1; i++) {
		area += 0.5f * fabsf(dot(n, cross(sub(contacts[i].point, contacts[0].point), sub(contacts[i + 1].point, contacts[0].point))));
	}
	return area > PATCH_MIN_AREA ? area : PATCH_MIN_AREA;
}


// Apply an impulse (linear + angular) to a body pair.
// Uses precomputed world-space inverse inertia (iw_diag/iw_off) for speed.
static SIMD_FORCEINLINE void apply_impulse(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, v3 impulse)
{
	a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(impulse, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_world_mul(a, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_world_mul(b, cross(r_b, impulse)));
}

// Apply a scalar impulse along a precomputed direction.
// w_a/w_b = precomputed I_w * cross(r, direction), so angular part is just scale.
static SIMD_FORCEINLINE void apply_impulse_row(BodyHot* a, BodyHot* b, v3 direction, v3 w_a, v3 w_b, float delta)
{
	v3 P = scale(direction, delta);
	a->velocity = sub(a->velocity, scale(P, a->inv_mass));
	b->velocity = add(b->velocity, scale(P, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, scale(w_a, delta));
	b->angular_velocity = add(b->angular_velocity, scale(w_b, delta));
}

// Match a new contact to a cached contact by feature ID. Returns index or -1.
static int warm_match(WarmManifold* wm, uint32_t feature_id)
{
	for (int i = 0; i < wm->count; i++) {
		if (wm->contacts[i].feature_id == feature_id) return i;
	}
	return -1;
}

#include "pre_solve_manifold.inc"

static void solver_pre_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, CK_DYNA PatchContact** out_pc, float dt)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	CK_DYNA PatchContact*   pc = NULL;
	// Fixed-stride: each manifold gets MAX_CONTACTS slots. Enables parallel dispatch.
	if (manifold_count == 0) { *out_sm = sm; *out_sc = sc; if (out_pc) *out_pc = pc; else afree(pc); return; }
	afit(sm, manifold_count); asetlen(sm, manifold_count);
	int total_contacts = manifold_count * MAX_CONTACTS;
	afit(sc, total_contacts); asetlen(sc, total_contacts);
	afit(pc, total_contacts); asetlen(pc, total_contacts);
	// Only contact_count must be zero for early-return manifolds; the full struct
	// is written by pre_solve_manifold for all active manifolds.
	for (int i = 0; i < manifold_count; i++) sm[i].contact_count = 0;

	// Precompute softness for dynamic-dynamic and dynamic-static pairs.
	float soft_dd = 0, bias_dd = 0, soft_ds = 0, bias_ds = 0;
	if (w->solver_type != SOLVER_SI) {
		float h1 = w->contact_hertz, h2 = h1 * 2.0f;
		float dr = w->contact_damping_ratio;
		float o1 = 6.28318530718f * h1, o2 = 6.28318530718f * h2;
		float d1 = 2*dr*o1, k1 = o1*o1, den1 = dt*d1 + dt*dt*k1;
		float d2 = 2*dr*o2, k2 = o2*o2, den2 = dt*d2 + dt*dt*k2;
		if (den1 > 1e-12f) { soft_dd = 1.0f / den1; bias_dd = dt * k1 * soft_dd; }
		if (den2 > 1e-12f) { soft_ds = 1.0f / den2; bias_ds = dt * k2 * soft_ds; }
	}

	for (int i = 0; i < manifold_count; i++) {
		pre_solve_manifold(w, &manifolds[i], i, sm, sc, pc, dt, soft_dd, bias_dd, soft_ds, bias_ds);
	}

	// Apply warm start impulses
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue; // skip empty (static-static filtered)
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (s->lambda_n == 0.0f) continue;
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n));
		}
		// Warm start manifold-level friction
		if ((m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f)) {
			v3 P = add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, P);
		}
		// Warm start torsional friction (pure angular impulse along normal)
		if (m->lambda_twist != 0.0f) {
			v3 twist_impulse = scale(m->normal, m->lambda_twist);
			BodyState* sa = &w->body_state[m->body_a]; BodyState* sb = &w->body_state[m->body_b];
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, twist_impulse));
		}
	}

	*out_sm = sm;
	*out_sc = sc;
	if (out_pc) *out_pc = pc; else afree(pc);
}

// NGS position correction: directly fix remaining penetration after velocity solve
// and position integration. Operates on positions only, no velocity modification.
static void solver_position_correct(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
{
	for (int iter = 0; iter < w->position_iters; iter++) {
		for (int i = 0; i < sm_count; i++) {
			SolverManifold* m = &sm[i];
			if (m->contact_count == 0) continue;
			BodyHot* a = &w->body_hot[m->body_a];
			BodyHot* b = &w->body_hot[m->body_b];
			BodyState* sa = &w->body_state[m->body_a];
			BodyState* sb = &w->body_state[m->body_b];
			float inv_mass_sum = a->inv_mass + b->inv_mass;

			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];

				// Recompute separation from current positions (r_a/r_b are close enough
				// after small position corrections — matches solver_relax_contacts approach).
				v3 p_a = add(sa->position, s->r_a);
				v3 p_b = add(sb->position, s->r_b);
				float separation = dot(sub(p_b, p_a), s->normal) - s->penetration;

				float C = fminf(0.0f, separation + SOLVER_SLOP);
				if (C >= 0.0f) continue;

				float correction = -SOLVER_POS_BAUMGARTE * C;
				if (correction > SOLVER_POS_MAX_CORRECTION)
					correction = SOLVER_POS_MAX_CORRECTION;

				// Effective mass for position correction (linear only for speed)
				float k = inv_mass_sum;
				float delta = correction / k;

				v3 P = scale(s->normal, delta);
				sa->position = sub(sa->position, scale(P, a->inv_mass));
				sb->position = add(sb->position, scale(P, b->inv_mass));
			}
		}
	}
}

// Relax contacts: recompute separation and bias from current body positions.
// Called after each substep's integrate_positions for SOFT_STEP and BLOCK solvers.
// Keeps normals, tangent basis, effective mass unchanged — only refreshes the RHS.
static void solver_relax_contacts(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, float dt)
{
	// Precompute bias rates: dynamic-dynamic and dynamic-static variants.
	float bias_rate_dd = 0.0f, bias_rate_ds = 0.0f;
	{
		float hertz = w->contact_hertz;
		float omega = 2.0f * 3.14159265f * hertz;
		float d = 2.0f * w->contact_damping_ratio * omega;
		float k = omega * omega;
		float hd = dt * d, hhk = dt * dt * k;
		float denom = hd + hhk;
		bias_rate_dd = (denom > 1e-12f) ? dt * k / denom : 0.0f;
		float omega2 = 2.0f * 3.14159265f * hertz * 2.0f;
		float d2 = 2.0f * w->contact_damping_ratio * omega2;
		float k2 = omega2 * omega2;
		float hd2 = dt * d2, hhk2 = dt * dt * k2;
		float denom2 = hd2 + hhk2;
		bias_rate_ds = (denom2 > 1e-12f) ? dt * k2 / denom2 : 0.0f;
	}
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue;
		BodyState* sa = &w->body_state[m->body_a];
		BodyState* sb = &w->body_state[m->body_b];
		float bias_rate = (m->inv_mass_a == 0.0f || m->inv_mass_b == 0.0f) ? bias_rate_ds : bias_rate_dd;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			v3 p_a = add(sa->position, s->r_a);
			v3 p_b = add(sb->position, s->r_b);
			float separation = dot(sub(p_b, p_a), s->normal) - s->penetration;

			float pen = -separation - SOLVER_SLOP;
			s->bias = pen > 0.0f ? -bias_rate * pen : 0.0f;
			if (s->bias < -w->max_push_velocity) s->bias = -w->max_push_velocity;

			if (s->bounce != 0.0f) s->bias = 0.0f;
		}
	}
}

// Persistent warm cache: update active pairs, age stale pairs, evict old stale.
// BEPU-style freshness: entries survive one extra frame so warm data isn't lost
// when narrowphase misses a pair for a single frame (FP noise).
static void solver_post_solve(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, InternalManifold* manifolds, int manifold_count)
{
	// Update active pairs in-place. Avoids double hash lookup + full copy for existing entries.
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue;
		uint64_t key = body_pair_key(m->body_a, m->body_b);

		WarmManifold* wm = map_get_ptr(w->warm_cache, key);
		if (!wm) { map_set(w->warm_cache, key, (WarmManifold){0}); wm = map_get_ptr(w->warm_cache, key); }
		wm->count = m->contact_count;
		wm->stale = 0;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			wm->contacts[ci] = (WarmContact){
				.feature_id = s->feature_id,
				.r_a = s->r_a,
				.lambda_n = s->lambda_n,
			};
		}
		wm->manifold_lambda_t1 = m->lambda_t1;
		wm->manifold_lambda_t2 = m->lambda_t2;
		wm->manifold_lambda_twist = m->lambda_twist;
	}

	afree(sm);
	afree(sc);
}

// Age and evict stale warm cache entries. Called once per frame (not per sub-step).
// Single pass: increment stale and evict in one traversal.
static void warm_cache_age_and_evict(WorldInternal* w)
{
	int i = 0;
	while (i < map_size(w->warm_cache)) {
		if (++w->warm_cache[i].stale > 1)
			map_del(w->warm_cache, map_key(w->warm_cache, i));
		else
			i++;
	}
}

// -----------------------------------------------------------------------------
// Graph coloring: assign colors to constraints so no two in same color share a body.
// Uses uint64_t bitmask per body (max 64 colors).

static void color_constraints(ConstraintRef* refs, int count, int body_count, int* out_batch_starts, int* out_color_count)
{
	uint64_t* body_colors = (uint64_t*)CK_ALLOC(body_count * sizeof(uint64_t));
	memset(body_colors, 0, body_count * sizeof(uint64_t));
	int max_color = 0;

	for (int i = 0; i < count; i++) {
		uint64_t used = body_colors[refs[i].body_a] | body_colors[refs[i].body_b];
		uint64_t avail = ~used;
		int color = 0;
		if (avail) {
			// Find lowest set bit (MSVC: _BitScanForward64, GCC: __builtin_ctzll)
#ifdef _MSC_VER
			unsigned long idx;
			_BitScanForward64(&idx, avail);
			color = (int)idx;
#else
			color = __builtin_ctzll(avail);
#endif
		}
		assert(color < 64);
		refs[i].color = (uint8_t)color;
		uint64_t bit = 1ULL << color;
		body_colors[refs[i].body_a] |= bit;
		body_colors[refs[i].body_b] |= bit;
		if (color > max_color) max_color = color;
	}

	// Counting sort by color
	int color_count = max_color + 1;
	int counts[64] = {0};
	for (int i = 0; i < count; i++) counts[refs[i].color]++;

	int offsets[64];
	offsets[0] = 0;
	for (int c = 1; c < color_count; c++) offsets[c] = offsets[c-1] + counts[c-1];

	// Record batch starts before sorting
	for (int c = 0; c < color_count; c++) out_batch_starts[c] = offsets[c];
	out_batch_starts[color_count] = count; // sentinel

	ConstraintRef* sorted = (ConstraintRef*)CK_ALLOC(count * sizeof(ConstraintRef));
	for (int i = 0; i < count; i++) {
		sorted[offsets[refs[i].color]++] = refs[i];
	}
	memcpy(refs, sorted, count * sizeof(ConstraintRef));

	CK_FREE(sorted);
	CK_FREE(body_colors);
	*out_color_count = color_count;
}

// Forward declaration for generic joint solver (defined in joints.c, included after solver.c).
static void solve_joint(WorldInternal* w, SolverJoint* s);

// Dispatch a single constraint solve by type.
static void solve_constraint(WorldInternal* w, ConstraintRef* ref, SolverManifold* sm, SolverContact* sc, SolverJoint* joints)
{
	switch (ref->type) {
	case CTYPE_CONTACT: {
		SolverManifold* m = &sm[ref->index];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		{
			float total_lambda_n = 0.0f;
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				float vn = dot(sub(b->velocity, a->velocity), s->normal) + dot(b->angular_velocity, s->rn_b) - dot(a->angular_velocity, s->rn_a);
				float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
				float old_n = s->lambda_n;
				s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
				apply_impulse_row(a, b, s->normal, s->w_n_a, s->w_n_b, s->lambda_n - old_n);
				total_lambda_n += s->lambda_n;
			}

			// Manifold-level 2D friction at centroid, clamped by aggregate normal force
			float max_f = m->friction * total_lambda_n;
			float vt1 = dot(sub(b->velocity, a->velocity), m->tangent1) + dot(b->angular_velocity, m->rct1_b) - dot(a->angular_velocity, m->rct1_a);
			float old_t1 = m->lambda_t1;
			m->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + m->eff_mass_t1 * (-vt1), max_f));
			apply_impulse_row(a, b, m->tangent1, m->w_t1_a, m->w_t1_b, m->lambda_t1 - old_t1);

			float vt2 = dot(sub(b->velocity, a->velocity), m->tangent2) + dot(b->angular_velocity, m->rct2_b) - dot(a->angular_velocity, m->rct2_a);
			float old_t2 = m->lambda_t2;
			m->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + m->eff_mass_t2 * (-vt2), max_f));
			apply_impulse_row(a, b, m->tangent2, m->w_t2_a, m->w_t2_b, m->lambda_t2 - old_t2);

			// Torsional friction: resist spin around contact normal
			float max_twist = m->friction * total_lambda_n * m->patch_radius;
			float w_rel = dot(sub(b->angular_velocity, a->angular_velocity), m->normal);
			float lambda_tw = m->eff_mass_twist * (-w_rel);
			float old_tw = m->lambda_twist;
			m->lambda_twist = fmaxf(-max_twist, fminf(old_tw + lambda_tw, max_twist));
			float delta_tw = m->lambda_twist - old_tw;
			a->angular_velocity = sub(a->angular_velocity, scale(m->w_tw_a, delta_tw));
			b->angular_velocity = add(b->angular_velocity, scale(m->w_tw_b, delta_tw));
		}
		break;
	}
	case CTYPE_JOINT: {
		solve_joint(w, &joints[ref->index]);
		break;
	}
	}
}

// --- BodyHot fast path: lean 80-byte body state for PGS iteration ---

// World-space inverse inertia multiply for BodyHot (macro to guarantee inlining).
#define sv_inertia_mul(h, v) V3((h)->iw_diag.x*(v).x + (h)->iw_off.x*(v).y + (h)->iw_off.y*(v).z, (h)->iw_off.x*(v).x + (h)->iw_diag.y*(v).y + (h)->iw_off.z*(v).z, (h)->iw_off.y*(v).x + (h)->iw_off.z*(v).y + (h)->iw_diag.z*(v).z)

// Apply impulse row to velocity-only body state. inv_mass comes from the manifold (cached).
static SIMD_FORCEINLINE void apply_impulse_row_sv(BodyHot* a, BodyHot* b, float ima, float imb, v3 direction, v3 w_a, v3 w_b, float delta)
{
	v3 P = scale(direction, delta);
	a->velocity = sub(a->velocity, scale(P, ima));
	b->velocity = add(b->velocity, scale(P, imb));
	a->angular_velocity = sub(a->angular_velocity, scale(w_a, delta));
	b->angular_velocity = add(b->angular_velocity, scale(w_b, delta));
}

// Solve one patch-friction manifold using BodyHot arrays directly.
// Body inv_mass is read from the manifold (cached in pre_solve), never from body arrays.
// BEPU-style: recompute cross+inertia inline instead of reading precomputed data.
// Trades ALU (cheap in Release) for bandwidth (expensive at 10K+ bodies).
// SolverContact only needs: r_a, r_b, eff_mass_n, bias, bounce, softness, lambda_n.
static SIMD_FORCEINLINE void solve_contact_patch_sv(BodyHot* bodies, SolverManifold* m, PatchContact* pc)
{
	BodyHot* a = &bodies[m->body_a];
	BodyHot* b = &bodies[m->body_b];
	float ima = m->inv_mass_a, imb = m->inv_mass_b;

	v3 normal = m->normal;
	float inv_mass_sum = ima + imb;
	float linear_vn = dot(sub(b->velocity, a->velocity), normal);
	float total_lambda_n = 0.0f;
	for (int ci = 0; ci < m->contact_count; ci++) {
		PatchContact* s = &pc[m->contact_start + ci];
		v3 rn_a = s->rn_a;
		v3 rn_b = s->rn_b;
		float vn = linear_vn + dot(b->angular_velocity, rn_b) - dot(a->angular_velocity, rn_a);
		float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
		float old_n = s->lambda_n;
		s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
		float delta = s->lambda_n - old_n;
		// Apply: fused linear impulse + precomputed angular direction.
		a->velocity = sub(a->velocity, scale(normal, delta * ima));
		b->velocity = add(b->velocity, scale(normal, delta * imb));
		a->angular_velocity = sub(a->angular_velocity, scale(s->w_n_a, delta));
		b->angular_velocity = add(b->angular_velocity, scale(s->w_n_b, delta));
		linear_vn += delta * inv_mass_sum;
		total_lambda_n += s->lambda_n;
	}

	float max_f = m->friction * total_lambda_n;

	// Jacobi-style friction: compute both tangent + torsional from same velocity snapshot,
	// then apply all deltas in one combined pass. Fewer velocity read/write round-trips.
	v3 dv = sub(b->velocity, a->velocity);
	v3 wa_snap = a->angular_velocity, wb_snap = b->angular_velocity;
	float vt1 = dot(dv, m->tangent1) + dot(wb_snap, m->rct1_b) - dot(wa_snap, m->rct1_a);
	float old_t1 = m->lambda_t1;
	m->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + m->eff_mass_t1 * (-vt1), max_f));
	float dt1 = m->lambda_t1 - old_t1;

	float vt2 = dot(dv, m->tangent2) + dot(wb_snap, m->rct2_b) - dot(wa_snap, m->rct2_a);
	float old_t2 = m->lambda_t2;
	m->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + m->eff_mass_t2 * (-vt2), max_f));
	float dt2 = m->lambda_t2 - old_t2;

	float max_twist = m->friction * total_lambda_n * m->patch_radius;
	float w_rel = dot(sub(wb_snap, wa_snap), m->normal);
	float old_tw = m->lambda_twist;
	m->lambda_twist = fmaxf(-max_twist, fminf(old_tw + m->eff_mass_twist * (-w_rel), max_twist));
	float dtw = m->lambda_twist - old_tw;

	// Single combined apply for all friction rows.
	v3 lin_impulse = add(scale(m->tangent1, dt1), scale(m->tangent2, dt2));
	a->velocity = sub(a->velocity, scale(lin_impulse, ima));
	b->velocity = add(b->velocity, scale(lin_impulse, imb));
	v3 ang_a = add(add(scale(m->w_t1_a, dt1), scale(m->w_t2_a, dt2)), scale(m->w_tw_a, dtw));
	v3 ang_b = add(add(scale(m->w_t1_b, dt1), scale(m->w_t2_b, dt2)), scale(m->w_tw_b, dtw));
	a->angular_velocity = sub(a->angular_velocity, ang_a);
	b->angular_velocity = add(b->angular_velocity, ang_b);
}

// --- SIMD 4-wide batch solver (BEPU-style: store minimal, recompute inline) ---
// Process 4 manifolds simultaneously using SSE. Each lane = one manifold.
// Contacts processed in lockstep: contact[0] across all 4, then [1], etc.
// Store only raw offsets (r_a, r_b, centroid_r); recompute cross products and
// inertia terms inline each iteration. Trades cheap ALU for expensive bandwidth.
#if SIMD_SSE

// Lean contact layer: raw offsets + scalar prestep. Cross/inertia recomputed inline.
typedef struct PGS_ContactLayer4
{
	v3w r_a, r_b;
	simd4f eff_mass_n, bias, bounce, softness;
	simd4f lambda_n;
} PGS_ContactLayer4;

// Persistent SoA batch: built once per frame, reused across iterations + substeps.
// Lean layout: directions + offsets stored as v3w; all angular terms recomputed inline.
typedef struct PGS_Batch4
{
	int body_a[4], body_b[4];
	int max_contacts;
	simd4f inv_mass_a, inv_mass_b;
	PGS_ContactLayer4 cp[MAX_CONTACTS];
	v3w tangent1, tangent2, normal;
	v3w centroid_r_a, centroid_r_b;
	simd4f eff_mass_t1, eff_mass_t2, eff_mass_twist;
	simd4f lambda_t1, lambda_t2, lambda_twist;
	simd4f friction, patch_radius;
	int manifold_idx[4];
	int lane_count;
} PGS_Batch4;

// --- Generic SIMD gather/scatter helpers ---
// These work with any struct array + field name. Portable across constraint types.

// Gather a v3 field from 4 struct entries into v3w via hardware transpose.
// Requires the v3 field to have a .m member (simd4f). Fast: 4 loads + 1 shuffle.
#define GATHER_V3(out, arr, i0, i1, i2, i3, field) do { \
	simd4f _r0=(arr)[i0].field.m, _r1=(arr)[i1].field.m, _r2=(arr)[i2].field.m, _r3=(arr)[i3].field.m; \
	simd_transpose4(&_r0,&_r1,&_r2,&_r3); (out).x=_r0; (out).y=_r1; (out).z=_r2; } while(0)

// Scatter v3w back to 4 struct entries via reverse transpose.
#define SCATTER_V3(arr, i0, i1, i2, i3, field, src) do { \
	simd4f _r0=(src).x, _r1=(src).y, _r2=(src).z, _r3=simd_zero(); \
	simd_transpose4(&_r0,&_r1,&_r2,&_r3); \
	(arr)[i0].field.m=_r0; (arr)[i1].field.m=_r1; (arr)[i2].field.m=_r2; (arr)[i3].field.m=_r3; } while(0)

// Gather one v3 field from up to 4 indexed structs into v3w. Zeroes inactive lanes.
// Slower than GATHER_V3 (scalar path), but works with any array + index list.
#define GATHER_V3_INDEXED(dst, arr, indices, count, field) do { \
	v3 _v[4] = {{0}}; \
	for (int _j = 0; _j < (count); _j++) { _v[_j] = (arr)[(indices)[_j]].field; } \
	dst = v3w_load4(_v[3], _v[2], _v[1], _v[0]); \
} while(0)

// Gather one float field from up to 4 indexed structs into simd4f. Zeroes inactive lanes.
#define GATHER_F1(dst, arr, indices, count, field) do { \
	float _buf[4] = {0}; \
	for (int _j = 0; _j < (count); _j++) { _buf[_j] = (arr)[(indices)[_j]].field; } \
	dst = simd_load(_buf); \
} while(0)

// Scatter one simd4f back to a float field on up to 4 indexed structs.
#define SCATTER_F1(arr, indices, count, field, src) do { \
	float _buf[4]; simd_store(_buf, src); \
	for (int _j = 0; _j < (count); _j++) { (arr)[(indices)[_j]].field = _buf[_j]; } \
} while(0)

// --- BodyWide: 4-wide SIMD body state for constraint solvers ---
// Bundles all solver-needed fields for 4 bodies. One gather call replaces 4-8 GATHER_V3 calls.
// Reusable across constraint types (contacts, joints, springs -- anything that touches bodies).
typedef struct BodyWide
{
	v3w vel, angvel;
	v3w iw_diag, iw_off;
	simd4f inv_mass;
} BodyWide;

static inline BodyWide body_wide_gather(BodyHot* arr, int i0, int i1, int i2, int i3)
{
	BodyWide bw;
	GATHER_V3(bw.vel, arr, i0, i1, i2, i3, velocity);
	GATHER_V3(bw.angvel, arr, i0, i1, i2, i3, angular_velocity);
	GATHER_V3(bw.iw_diag, arr, i0, i1, i2, i3, iw_diag);
	GATHER_V3(bw.iw_off, arr, i0, i1, i2, i3, iw_off);
	// inv_mass: scalar field, no .m member -- gather manually
	float _im[4] = { arr[i0].inv_mass, arr[i1].inv_mass, arr[i2].inv_mass, arr[i3].inv_mass };
	bw.inv_mass = simd_load(_im);
	return bw;
}

static inline void body_wide_scatter(BodyHot* arr, int i0, int i1, int i2, int i3, BodyWide* bw)
{
	SCATTER_V3(arr, i0, i1, i2, i3, velocity, bw->vel);
	SCATTER_V3(arr, i0, i1, i2, i3, angular_velocity, bw->angvel);
	// iw_diag, iw_off, inv_mass are read-only during solve -- don't scatter back
}

// Symmetric 3x3 inverse-inertia multiply: I_w * v, where I_w = diag(xx,yy,zz) + off(xy,xz,yz).
static inline v3w iw_mul(v3w iw_d, v3w iw_o, v3w v)
{
	return (v3w){
		simd_add(simd_add(simd_mul(iw_d.x, v.x), simd_mul(iw_o.x, v.y)), simd_mul(iw_o.y, v.z)),
		simd_add(simd_add(simd_mul(iw_o.x, v.x), simd_mul(iw_d.y, v.y)), simd_mul(iw_o.z, v.z)),
		simd_add(simd_add(simd_mul(iw_o.y, v.x), simd_mul(iw_o.z, v.y)), simd_mul(iw_d.z, v.z))
	};
}

static void pgs_batch4_prepare(PGS_Batch4* bt, SolverManifold* sm, int* indices, int count, SolverContact* sc, PatchContact* pc)
{
	float buf[4];
	#define GATHER1(dst, field) for (int j = 0; j < 4; j++) { buf[j] = (j < count) ? sm[indices[j]].field : 0; } dst = simd_load(buf)
	bt->lane_count = count;
	for (int j = 0; j < 4; j++) { bt->manifold_idx[j] = (j < count) ? indices[j] : 0; bt->body_a[j] = (j < count) ? sm[indices[j]].body_a : 0; bt->body_b[j] = (j < count) ? sm[indices[j]].body_b : 0; }
	GATHER1(bt->inv_mass_a, inv_mass_a); GATHER1(bt->inv_mass_b, inv_mass_b);
	GATHER1(bt->friction, friction); GATHER1(bt->patch_radius, patch_radius);
	GATHER1(bt->eff_mass_t1, eff_mass_t1); GATHER1(bt->eff_mass_t2, eff_mass_t2); GATHER1(bt->eff_mass_twist, eff_mass_twist);
	GATHER1(bt->lambda_t1, lambda_t1); GATHER1(bt->lambda_t2, lambda_t2); GATHER1(bt->lambda_twist, lambda_twist);
	#undef GATHER1
	// Pack v3 manifold fields via hardware transpose (GATHER_V3 on SolverManifold array).
	int mi0 = indices[0], mi1 = count > 1 ? indices[1] : mi0, mi2 = count > 2 ? indices[2] : mi0, mi3 = count > 3 ? indices[3] : mi0;
	GATHER_V3(bt->normal, sm, mi0, mi1, mi2, mi3, normal);
	GATHER_V3(bt->tangent1, sm, mi0, mi1, mi2, mi3, tangent1);
	GATHER_V3(bt->tangent2, sm, mi0, mi1, mi2, mi3, tangent2);
	GATHER_V3(bt->centroid_r_a, sm, mi0, mi1, mi2, mi3, centroid_r_a);
	GATHER_V3(bt->centroid_r_b, sm, mi0, mi1, mi2, mi3, centroid_r_b);

	// Pack contact layers: raw r_a/r_b from SolverContact + scalar prestep from PatchContact.
	bt->max_contacts = 0;
	for (int j = 0; j < count; j++) { if (sm[indices[j]].contact_count > bt->max_contacts) { bt->max_contacts = sm[indices[j]].contact_count; } }
	for (int cp_idx = 0; cp_idx < bt->max_contacts; cp_idx++) {
		PGS_ContactLayer4* cl = &bt->cp[cp_idx];
		v3 ra[4] = {{0}}, rb[4] = {{0}};
		float emn[4]={0}, bi[4]={0}, bnc[4]={0}, sft[4]={0}, lam[4]={0};
		for (int j = 0; j < count; j++) {
			if (cp_idx >= sm[indices[j]].contact_count) { continue; }
			int ci = sm[indices[j]].contact_start + cp_idx;
			ra[j] = sc[ci].r_a; rb[j] = sc[ci].r_b;
			PatchContact* s = &pc[ci];
			emn[j]=s->eff_mass_n; bi[j]=s->bias; bnc[j]=s->bounce; sft[j]=s->softness; lam[j]=s->lambda_n;
		}
		cl->r_a = v3w_load4(ra[3], ra[2], ra[1], ra[0]);
		cl->r_b = v3w_load4(rb[3], rb[2], rb[1], rb[0]);
		cl->eff_mass_n = simd_load(emn); cl->bias = simd_load(bi); cl->bounce = simd_load(bnc); cl->softness = simd_load(sft); cl->lambda_n = simd_load(lam);
	}
}

// Lightweight refresh: only update bias, bounce, and lambda from PatchContact/SolverManifold.
// Called on substep 2+ when structural data (normals, offsets, etc.) hasn't changed.
static void pgs_batch4_refresh(PGS_Batch4* bt, SolverManifold* sm, PatchContact* pc)
{
	int count = bt->lane_count;
	float buf[4];
	for (int j = 0; j < 4; j++) { buf[j] = (j < count) ? sm[bt->manifold_idx[j]].lambda_t1 : 0; } bt->lambda_t1 = simd_load(buf);
	for (int j = 0; j < 4; j++) { buf[j] = (j < count) ? sm[bt->manifold_idx[j]].lambda_t2 : 0; } bt->lambda_t2 = simd_load(buf);
	for (int j = 0; j < 4; j++) { buf[j] = (j < count) ? sm[bt->manifold_idx[j]].lambda_twist : 0; } bt->lambda_twist = simd_load(buf);
	for (int cp_idx = 0; cp_idx < bt->max_contacts; cp_idx++) {
		float bi[4]={0}, lam[4]={0};
		for (int j = 0; j < count; j++) {
			int mi = bt->manifold_idx[j];
			if (cp_idx >= sm[mi].contact_count) { continue; }
			PatchContact* s = &pc[sm[mi].contact_start + cp_idx];
			bi[j] = s->bias; lam[j] = s->lambda_n;
		}
		bt->cp[cp_idx].bias = simd_load(bi);
		bt->cp[cp_idx].lambda_n = simd_load(lam);
	}
}

static void solve_contact_batch4_sv(BodyHot* bodies, PGS_Batch4* b)
{
	int i0a=b->body_a[0], i1a=b->body_a[1], i2a=b->body_a[2], i3a=b->body_a[3];
	int i0b=b->body_b[0], i1b=b->body_b[1], i2b=b->body_b[2], i3b=b->body_b[3];

	// Gather body state as flat v3w locals (not struct -- avoids register spill).
	v3w va, wa, vb, wb, iw_d_a, iw_o_a, iw_d_b, iw_o_b;
	GATHER_V3(va, bodies, i0a, i1a, i2a, i3a, velocity);
	GATHER_V3(wa, bodies, i0a, i1a, i2a, i3a, angular_velocity);
	GATHER_V3(iw_d_a, bodies, i0a, i1a, i2a, i3a, iw_diag);
	GATHER_V3(iw_o_a, bodies, i0a, i1a, i2a, i3a, iw_off);
	GATHER_V3(vb, bodies, i0b, i1b, i2b, i3b, velocity);
	GATHER_V3(wb, bodies, i0b, i1b, i2b, i3b, angular_velocity);
	GATHER_V3(iw_d_b, bodies, i0b, i1b, i2b, i3b, iw_diag);
	GATHER_V3(iw_o_b, bodies, i0b, i1b, i2b, i3b, iw_off);

	simd4f zero = simd_zero();
	simd4f inv_mass_sum = simd_add(b->inv_mass_a, b->inv_mass_b);
	v3w dv = v3w_sub(vb, va);
	simd4f linear_vn = v3w_dot(dv, b->normal);
	simd4f total_lambda_n = zero;

	// Normal contacts: recompute cross(r, n) and iw*cross(r, n) inline each iteration.
	for (int ci = 0; ci < b->max_contacts; ci++) {
		PGS_ContactLayer4* cl = &b->cp[ci];
		v3w rn_a = v3w_cross(cl->r_a, b->normal);
		v3w rn_b = v3w_cross(cl->r_b, b->normal);
		v3w wn_a = iw_mul(iw_d_a, iw_o_a, rn_a);
		v3w wn_b = iw_mul(iw_d_b, iw_o_b, rn_b);
		simd4f vn = simd_add(linear_vn, simd_sub(v3w_dot(wb, rn_b), v3w_dot(wa, rn_a)));
		simd4f lambda_n = simd_mul(cl->eff_mass_n, simd_sub(simd_sub(zero, simd_add(vn, simd_add(cl->bias, cl->bounce))), simd_mul(cl->softness, cl->lambda_n)));
		simd4f new_lam = simd_max(simd_add(cl->lambda_n, lambda_n), zero);
		simd4f delta = simd_sub(new_lam, cl->lambda_n);
		cl->lambda_n = new_lam;
		simd4f da = simd_mul(delta, b->inv_mass_a), db = simd_mul(delta, b->inv_mass_b);
		va = v3w_sub(va, v3w_scale(b->normal, da));
		vb = v3w_add(vb, v3w_scale(b->normal, db));
		wa = v3w_sub(wa, v3w_scale(wn_a, delta));
		wb = v3w_add(wb, v3w_scale(wn_b, delta));
		linear_vn = simd_add(linear_vn, simd_mul(delta, inv_mass_sum));
		total_lambda_n = simd_add(total_lambda_n, new_lam);
	}

	// Jacobi friction: compute tangent + torsional from same velocity snapshot, apply once.
	simd4f max_f = simd_mul(b->friction, total_lambda_n);
	v3w dv2 = v3w_sub(vb, va);
	v3w wa_s = wa, wb_s = wb;

	// Tangent 1: recompute cross(centroid_r, tangent1) and iw*cross inline.
	v3w rct1_a = v3w_cross(b->centroid_r_a, b->tangent1), rct1_b = v3w_cross(b->centroid_r_b, b->tangent1);
	v3w wt1_a = iw_mul(iw_d_a, iw_o_a, rct1_a), wt1_b = iw_mul(iw_d_b, iw_o_b, rct1_b);
	simd4f vt1 = simd_add(v3w_dot(dv2, b->tangent1), simd_sub(v3w_dot(wb_s, rct1_b), v3w_dot(wa_s, rct1_a)));
	simd4f nt1 = simd_max(simd_neg(max_f), simd_min(simd_add(b->lambda_t1, simd_mul(b->eff_mass_t1, simd_neg(vt1))), max_f));
	simd4f dt1 = simd_sub(nt1, b->lambda_t1); b->lambda_t1 = nt1;

	// Tangent 2: recompute cross(centroid_r, tangent2) and iw*cross inline.
	v3w rct2_a = v3w_cross(b->centroid_r_a, b->tangent2), rct2_b = v3w_cross(b->centroid_r_b, b->tangent2);
	v3w wt2_a = iw_mul(iw_d_a, iw_o_a, rct2_a), wt2_b = iw_mul(iw_d_b, iw_o_b, rct2_b);
	simd4f vt2 = simd_add(v3w_dot(dv2, b->tangent2), simd_sub(v3w_dot(wb_s, rct2_b), v3w_dot(wa_s, rct2_a)));
	simd4f nt2 = simd_max(simd_neg(max_f), simd_min(simd_add(b->lambda_t2, simd_mul(b->eff_mass_t2, simd_neg(vt2))), max_f));
	simd4f dt2 = simd_sub(nt2, b->lambda_t2); b->lambda_t2 = nt2;

	// Torsional: recompute iw*normal inline.
	v3w wtw_a = iw_mul(iw_d_a, iw_o_a, b->normal), wtw_b = iw_mul(iw_d_b, iw_o_b, b->normal);
	simd4f max_tw = simd_mul(max_f, b->patch_radius);
	simd4f wrel = v3w_dot(v3w_sub(wb_s, wa_s), b->normal);
	simd4f ntw = simd_max(simd_neg(max_tw), simd_min(simd_add(b->lambda_twist, simd_mul(b->eff_mass_twist, simd_neg(wrel))), max_tw));
	simd4f dtw = simd_sub(ntw, b->lambda_twist); b->lambda_twist = ntw;

	// Single combined apply for all friction rows.
	v3w lin = v3w_add(v3w_scale(b->tangent1, dt1), v3w_scale(b->tangent2, dt2));
	va = v3w_sub(va, v3w_scale(lin, b->inv_mass_a));
	vb = v3w_add(vb, v3w_scale(lin, b->inv_mass_b));
	v3w ang_a = v3w_add(v3w_add(v3w_scale(wt1_a, dt1), v3w_scale(wt2_a, dt2)), v3w_scale(wtw_a, dtw));
	v3w ang_b = v3w_add(v3w_add(v3w_scale(wt1_b, dt1), v3w_scale(wt2_b, dt2)), v3w_scale(wtw_b, dtw));
	wa = v3w_sub(wa, ang_a);
	wb = v3w_add(wb, ang_b);

	// Scatter velocities back via reverse transpose.
	SCATTER_V3(bodies, i0a, i1a, i2a, i3a, velocity, va);
	SCATTER_V3(bodies, i0a, i1a, i2a, i3a, angular_velocity, wa);
	SCATTER_V3(bodies, i0b, i1b, i2b, i3b, velocity, vb);
	SCATTER_V3(bodies, i0b, i1b, i2b, i3b, angular_velocity, wb);
}

#endif // SIMD_SSE
