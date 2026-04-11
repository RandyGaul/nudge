// solver.c -- contact constraint solver

// Analysis: Extremely well-studied and robust solver. Very simple, converges predictably. Has
// trouble propogating down longer chains, and dealing with large mass ratios. However, it's just
// so modular and well-understood it often becomes useful as a supporting algorithm in other more
// advanced solvers, or hybrid approaches.

// body_pair_key defined in collision.c (included before this file)

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
static float estimate_patch_area(Contact* contacts, int count)
{
	if (count < 3) return PATCH_MIN_AREA;
	// Fan triangulation from contacts[0]
	float area = 0.0f;
	for (int i = 1; i < count - 1; i++)
		area += 0.5f * len(cross(sub(contacts[i].point, contacts[0].point), sub(contacts[i + 1].point, contacts[0].point)));
	return area > PATCH_MIN_AREA ? area : PATCH_MIN_AREA;
}

static float compute_effective_mass(BodyHot* a, BodyHot* b, float inv_mass_sum, v3 r_a, v3 r_b, v3 dir)
{
	v3 ra_x_d = cross(r_a, dir);
	v3 rb_x_d = cross(r_b, dir);
	float k = inv_mass_sum + dot(cross(inv_inertia_world_mul(a, ra_x_d), r_a), dir) + dot(cross(inv_inertia_world_mul(b, rb_x_d), r_b), dir);
	return k > 1e-12f ? 1.0f / k : 0.0f;
}

// Apply an impulse (linear + angular) to a body pair.
// Uses precomputed world-space inverse inertia (iw_diag/iw_off) for speed.
static void apply_impulse(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, v3 impulse)
{
	a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(impulse, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_world_mul(a, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_world_mul(b, cross(r_b, impulse)));
}

// Apply a scalar impulse along a precomputed direction.
// w_a/w_b = precomputed I_w * cross(r, direction), so angular part is just scale.
static void apply_impulse_row(BodyHot* a, BodyHot* b, v3 direction, v3 w_a, v3 w_b, float delta)
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
	for (int i = 0; i < wm->count; i++)
		if (wm->contacts[i].feature_id == feature_id) return i;
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
	memset(sm, 0, manifold_count * sizeof(SolverManifold));
	memset(sc, 0, total_contacts * sizeof(SolverContact));
	memset(pc, 0, total_contacts * sizeof(PatchContact));

	for (int i = 0; i < manifold_count; i++)
		pre_solve_manifold(w, &manifolds[i], i, sm, sc, pc, dt);

	// Apply warm start impulses
	int patch_warm = (w->friction_model == FRICTION_PATCH);
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue; // skip empty (static-static filtered)
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (patch_warm) {
				if (s->lambda_n == 0.0f) continue;
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n));
			} else {
				if (s->lambda_n == 0.0f && s->lambda_t1 == 0.0f && s->lambda_t2 == 0.0f)
					continue;
				v3 P = add(add(scale(s->normal, s->lambda_n), scale(s->tangent1, s->lambda_t1)), scale(s->tangent2, s->lambda_t2));
				apply_impulse(a, b, s->r_a, s->r_b, P);
			}
		}
		// Warm start manifold-level friction
		if (patch_warm && (m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f)) {
			v3 P = add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, P);
		}
		// Warm start torsional friction (pure angular impulse along normal)
		if (patch_warm && m->lambda_twist != 0.0f) {
			v3 twist_impulse = scale(m->normal, m->lambda_twist);
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, twist_impulse));
		}
	}

	*out_sm = sm;
	*out_sc = sc;
	if (out_pc) *out_pc = pc; else afree(pc);
}

static void solver_iterate(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
{
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			// Relative velocity at contact point
			v3 dv = sub( add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));

			// --- Normal constraint ---
			float vn = dot(dv, s->normal);
			float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce));
			float old_n = s->lambda_n;
			s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
			float delta_n = s->lambda_n - old_n;
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, delta_n));

			// --- Friction tangent 1 ---
			dv = sub( add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt1 = dot(dv, s->tangent1);
			float lambda_t1 = s->eff_mass_t1 * (-vt1);
			float max_f = m->friction * s->lambda_n;
			float old_t1 = s->lambda_t1;
			s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + lambda_t1, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

			// --- Friction tangent 2 ---
			dv = sub( add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt2 = dot(dv, s->tangent2);
			float lambda_t2 = s->eff_mass_t2 * (-vt2);
			float old_t2 = s->lambda_t2;
			s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + lambda_t2, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
		}
	}
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
			float inv_mass_sum = a->inv_mass + b->inv_mass;

			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];

				// Recompute separation from current positions
				v3 r_a = sub(add(a->position, rotate(a->rotation, rotate(inv(a->rotation), s->r_a))), a->position);
				v3 r_b = sub(add(b->position, rotate(b->rotation, rotate(inv(b->rotation), s->r_b))), b->position);
				v3 p_a = add(a->position, r_a);
				v3 p_b = add(b->position, r_b);
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
				a->position = sub(a->position, scale(P, a->inv_mass));
				b->position = add(b->position, scale(P, b->inv_mass));
			}
		}
	}
}

// Relax contacts: recompute separation and bias from current body positions.
// Called after each substep's integrate_positions for SOFT_STEP and BLOCK solvers.
// Keeps normals, tangent basis, effective mass unchanged — only refreshes the RHS.
static void solver_relax_contacts(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, float dt)
{
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue;
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		float hertz = w->contact_hertz;
		if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
		float omega = 2.0f * 3.14159265f * hertz;
		float d = 2.0f * w->contact_damping_ratio * omega;
		float k = omega * omega;
		float hd = dt * d, hhk = dt * dt * k;
		float denom = hd + hhk;
		float bias_rate = (denom > 1e-12f) ? dt * k / denom : 0.0f;

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			v3 p_a = add(a->position, s->r_a);
			v3 p_b = add(b->position, s->r_b);
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
	// Update active pairs (per sub-step)
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue; // skip empty (fixed-stride padding)
		uint64_t key = body_pair_key(m->body_a, m->body_b);

		WarmManifold wm = {0};
		wm.count = m->contact_count;
		wm.stale = 0;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			wm.contacts[ci] = (WarmContact){
				.feature_id = s->feature_id,
				.r_a = s->r_a,
				.lambda_n = s->lambda_n,
				.lambda_t1 = s->lambda_t1,
				.lambda_t2 = s->lambda_t2,
			};
		}
		if (w->friction_model == FRICTION_PATCH) {
			wm.manifold_lambda_t1 = m->lambda_t1;
			wm.manifold_lambda_t2 = m->lambda_t2;
			wm.manifold_lambda_twist = m->lambda_twist;
		}

		map_set(w->warm_cache, key, wm);
	}

	afree(sm);
	afree(sc);
}

// Age and evict stale warm cache entries. Called once per frame (not per sub-step).
static void warm_cache_age_and_evict(WorldInternal* w)
{
	map_each(w->warm_cache, i) w->warm_cache[i].stale++;
	// Evict entries stale for >1 frame. map_del swaps last entry into slot i,
	// so don't advance i when deleting (re-check the swapped entry).
	int i = 0;
	while (i < map_size(w->warm_cache)) {
		if (w->warm_cache[i].stale > 1)
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
	for (int i = 0; i < count; i++)
		sorted[offsets[refs[i].color]++] = refs[i];
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

		if (w->friction_model == FRICTION_PATCH) {
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
		} else {
			// Per-point Coulomb friction
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float vn = dot(dv, s->normal);
				float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
				float old_n = s->lambda_n;
				s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
				apply_impulse_row(a, b, s->normal, s->w_n_a, s->w_n_b, s->lambda_n - old_n);
			}
			// Per-contact friction (uses normal impulse from above)
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float max_f = m->friction * s->lambda_n;
				float vt1 = dot(dv, s->tangent1);
				float old_t1 = s->lambda_t1;
				s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + s->eff_mass_t1*(-vt1), max_f));
				apply_impulse_row(a, b, s->tangent1, s->w_t1_a, s->w_t1_b, s->lambda_t1 - old_t1);

				dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float vt2 = dot(dv, s->tangent2);
				float old_t2 = s->lambda_t2;
				s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + s->eff_mass_t2*(-vt2), max_f));
				apply_impulse_row(a, b, s->tangent2, s->w_t2_a, s->w_t2_b, s->lambda_t2 - old_t2);
			}
		}
		break;
	}
	case CTYPE_JOINT: {
		solve_joint(w, &joints[ref->index]);
		break;
	}
	}
}

// --- SolverBodyVel fast path: compact 32-byte body state for PGS iteration ---

// Sync body_hot velocity → body_vel (before PGS).
static void solver_sync_vel_in(WorldInternal* w)
{
	int count = asize(w->body_hot);
	split_ensure(w->body_vel, count - 1);
	for (int i = 0; i < count; i++) {
		w->body_vel[i].velocity = w->body_hot[i].velocity;
		w->body_vel[i].angular_velocity = w->body_hot[i].angular_velocity;
		w->body_vel[i].inv_mass = w->body_hot[i].inv_mass;
		w->body_vel[i].iw_diag = w->body_hot[i].iw_diag;
		w->body_vel[i].iw_off = w->body_hot[i].iw_off;
	}
}

// Sync body_vel → body_hot velocity (after PGS).
static void solver_sync_vel_out(WorldInternal* w)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		w->body_hot[i].velocity = w->body_vel[i].velocity;
		w->body_hot[i].angular_velocity = w->body_vel[i].angular_velocity;
	}
}

// World-space inverse inertia multiply for SolverBodyVel (macro to guarantee inlining).
#define sv_inertia_mul(h, v) V3((h)->iw_diag.x*(v).x + (h)->iw_off.x*(v).y + (h)->iw_off.y*(v).z, (h)->iw_off.x*(v).x + (h)->iw_diag.y*(v).y + (h)->iw_off.z*(v).z, (h)->iw_off.y*(v).x + (h)->iw_off.z*(v).y + (h)->iw_diag.z*(v).z)

// Apply impulse row to velocity-only body state. inv_mass comes from the manifold (cached).
static void apply_impulse_row_sv(SolverBodyVel* a, SolverBodyVel* b, float ima, float imb, v3 direction, v3 w_a, v3 w_b, float delta)
{
	v3 P = scale(direction, delta);
	a->velocity = sub(a->velocity, scale(P, ima));
	b->velocity = add(b->velocity, scale(P, imb));
	a->angular_velocity = sub(a->angular_velocity, scale(w_a, delta));
	b->angular_velocity = add(b->angular_velocity, scale(w_b, delta));
}

// Solve one patch-friction manifold using compact SolverBodyVel arrays.
// Body inv_mass is read from the manifold (cached in pre_solve), never from body arrays.
// BEPU-style: recompute cross+inertia inline instead of reading precomputed data.
// Trades ALU (cheap in Release) for bandwidth (expensive at 10K+ bodies).
// SolverContact only needs: r_a, r_b, eff_mass_n, bias, bounce, softness, lambda_n.
static SIMD_FORCEINLINE void solve_contact_patch_sv(SolverBodyVel* bodies, SolverManifold* m, PatchContact* pc)
{
	SolverBodyVel* a = &bodies[m->body_a];
	SolverBodyVel* b = &bodies[m->body_b];
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

// --- SIMD 4-wide batch solver (Box2D / BEPU style) ---
// Process 4 manifolds simultaneously using SSE. Each lane = one manifold.
// Contacts processed in lockstep: contact[0] across all 4, then [1], etc.
#if SIMD_SSE

// Pre-built SoA contact layer: 4 contacts (one per manifold in batch) for a given contact index.
typedef struct PGS_ContactLayer4
{
	__m128 rn_a_x, rn_a_y, rn_a_z, rn_b_x, rn_b_y, rn_b_z;
	__m128 wn_a_x, wn_a_y, wn_a_z, wn_b_x, wn_b_y, wn_b_z;
	__m128 eff_mass_n, bias, bounce, softness;
	__m128 lambda_n; // read/write — persists between iterations
} PGS_ContactLayer4;

// Persistent SoA batch: built ONCE per frame, reused across all iterations + substeps.
// Eliminates per-iteration manifold→SoA conversion (was 40x redundant).
typedef struct PGS_Batch4
{
	int body_a[4], body_b[4];
	int max_contacts;
	__m128 inv_mass_a, inv_mass_b, inv_mass_sum, friction, patch_radius;
	__m128 normal_x, normal_y, normal_z;
	__m128 tangent1_x, tangent1_y, tangent1_z, tangent2_x, tangent2_y, tangent2_z;
	__m128 eff_mass_t1, eff_mass_t2, eff_mass_twist;
	__m128 rct1_a_x, rct1_a_y, rct1_a_z, rct1_b_x, rct1_b_y, rct1_b_z;
	__m128 rct2_a_x, rct2_a_y, rct2_a_z, rct2_b_x, rct2_b_y, rct2_b_z;
	__m128 w_t1_a_x, w_t1_a_y, w_t1_a_z, w_t1_b_x, w_t1_b_y, w_t1_b_z;
	__m128 w_t2_a_x, w_t2_a_y, w_t2_a_z, w_t2_b_x, w_t2_b_y, w_t2_b_z;
	__m128 w_tw_a_x, w_tw_a_y, w_tw_a_z, w_tw_b_x, w_tw_b_y, w_tw_b_z;
	__m128 lambda_t1, lambda_t2, lambda_twist;
	PGS_ContactLayer4 cp[MAX_CONTACTS]; // pre-built contact layers
} PGS_Batch4;

#define SOA_DOT3(ax,ay,az,bx,by,bz) _mm_add_ps(_mm_add_ps(_mm_mul_ps(ax,bx), _mm_mul_ps(ay,by)), _mm_mul_ps(az,bz))

static void pgs_batch4_prepare(PGS_Batch4* bt, SolverManifold* sm, int* indices, int count, PatchContact* pc)
{
	float buf[4];
	#define GATHER1(dst, field) for (int j = 0; j < 4; j++) buf[j] = (j < count) ? sm[indices[j]].field : 0; dst = _mm_loadu_ps(buf)
	#define GATHER3(dx,dy,dz, field) for (int j = 0; j < 4; j++) { v3 v = (j < count) ? sm[indices[j]].field : V3(0,0,0); ((float*)&dx)[j]=v.x; ((float*)&dy)[j]=v.y; ((float*)&dz)[j]=v.z; }
	for (int j = 0; j < 4; j++) { bt->body_a[j] = (j < count) ? sm[indices[j]].body_a : 0; bt->body_b[j] = (j < count) ? sm[indices[j]].body_b : 0; }
	GATHER1(bt->inv_mass_a, inv_mass_a); GATHER1(bt->inv_mass_b, inv_mass_b);
	bt->inv_mass_sum = _mm_add_ps(bt->inv_mass_a, bt->inv_mass_b);
	GATHER1(bt->friction, friction); GATHER1(bt->patch_radius, patch_radius);
	GATHER1(bt->eff_mass_t1, eff_mass_t1); GATHER1(bt->eff_mass_t2, eff_mass_t2); GATHER1(bt->eff_mass_twist, eff_mass_twist);
	GATHER1(bt->lambda_t1, lambda_t1); GATHER1(bt->lambda_t2, lambda_t2); GATHER1(bt->lambda_twist, lambda_twist);
	GATHER3(bt->normal_x, bt->normal_y, bt->normal_z, normal);
	GATHER3(bt->tangent1_x, bt->tangent1_y, bt->tangent1_z, tangent1);
	GATHER3(bt->tangent2_x, bt->tangent2_y, bt->tangent2_z, tangent2);
	GATHER3(bt->rct1_a_x, bt->rct1_a_y, bt->rct1_a_z, rct1_a); GATHER3(bt->rct1_b_x, bt->rct1_b_y, bt->rct1_b_z, rct1_b);
	GATHER3(bt->rct2_a_x, bt->rct2_a_y, bt->rct2_a_z, rct2_a); GATHER3(bt->rct2_b_x, bt->rct2_b_y, bt->rct2_b_z, rct2_b);
	GATHER3(bt->w_t1_a_x, bt->w_t1_a_y, bt->w_t1_a_z, w_t1_a); GATHER3(bt->w_t1_b_x, bt->w_t1_b_y, bt->w_t1_b_z, w_t1_b);
	GATHER3(bt->w_t2_a_x, bt->w_t2_a_y, bt->w_t2_a_z, w_t2_a); GATHER3(bt->w_t2_b_x, bt->w_t2_b_y, bt->w_t2_b_z, w_t2_b);
	GATHER3(bt->w_tw_a_x, bt->w_tw_a_y, bt->w_tw_a_z, w_tw_a); GATHER3(bt->w_tw_b_x, bt->w_tw_b_y, bt->w_tw_b_z, w_tw_b);
	#undef GATHER1
	#undef GATHER3

	// Pre-build contact layers (constant data — reused across all iterations).
	bt->max_contacts = 0;
	for (int j = 0; j < count; j++) if (sm[indices[j]].contact_count > bt->max_contacts) bt->max_contacts = sm[indices[j]].contact_count;
	for (int cp_idx = 0; cp_idx < bt->max_contacts; cp_idx++) {
		PGS_ContactLayer4* cl = &bt->cp[cp_idx];
		float rnax[4]={0},rnay[4]={0},rnaz[4]={0},rnbx[4]={0},rnby[4]={0},rnbz[4]={0};
		float wnax[4]={0},wnay[4]={0},wnaz[4]={0},wnbx[4]={0},wnby[4]={0},wnbz[4]={0};
		float emn[4]={0},bias[4]={0},bnc[4]={0},sft[4]={0},lam[4]={0};
		for (int j = 0; j < count; j++) {
			if (cp_idx >= sm[indices[j]].contact_count) continue;
			PatchContact* s = &pc[sm[indices[j]].contact_start + cp_idx];
			rnax[j]=s->rn_a.x;rnay[j]=s->rn_a.y;rnaz[j]=s->rn_a.z;
			rnbx[j]=s->rn_b.x;rnby[j]=s->rn_b.y;rnbz[j]=s->rn_b.z;
			wnax[j]=s->w_n_a.x;wnay[j]=s->w_n_a.y;wnaz[j]=s->w_n_a.z;
			wnbx[j]=s->w_n_b.x;wnby[j]=s->w_n_b.y;wnbz[j]=s->w_n_b.z;
			emn[j]=s->eff_mass_n;bias[j]=s->bias;bnc[j]=s->bounce;sft[j]=s->softness;lam[j]=s->lambda_n;
		}
		cl->rn_a_x=_mm_loadu_ps(rnax);cl->rn_a_y=_mm_loadu_ps(rnay);cl->rn_a_z=_mm_loadu_ps(rnaz);
		cl->rn_b_x=_mm_loadu_ps(rnbx);cl->rn_b_y=_mm_loadu_ps(rnby);cl->rn_b_z=_mm_loadu_ps(rnbz);
		cl->wn_a_x=_mm_loadu_ps(wnax);cl->wn_a_y=_mm_loadu_ps(wnay);cl->wn_a_z=_mm_loadu_ps(wnaz);
		cl->wn_b_x=_mm_loadu_ps(wnbx);cl->wn_b_y=_mm_loadu_ps(wnby);cl->wn_b_z=_mm_loadu_ps(wnbz);
		cl->eff_mass_n=_mm_loadu_ps(emn);cl->bias=_mm_loadu_ps(bias);cl->bounce=_mm_loadu_ps(bnc);cl->softness=_mm_loadu_ps(sft);cl->lambda_n=_mm_loadu_ps(lam);
	}
}

static void solve_contact_batch4_sv(SolverBodyVel* bodies, PGS_Batch4* b)
{
	// Gather body velocities into SoA (only thing that changes per iteration)
	__m128 va_x, va_y, va_z, wa_x, wa_y, wa_z, vb_x, vb_y, vb_z, wb_x, wb_y, wb_z;
	{ float t[4];
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].velocity.x; va_x=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].velocity.y; va_y=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].velocity.z; va_z=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].angular_velocity.x; wa_x=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].angular_velocity.y; wa_y=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_a[j]].angular_velocity.z; wa_z=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].velocity.x; vb_x=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].velocity.y; vb_y=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].velocity.z; vb_z=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].angular_velocity.x; wb_x=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].angular_velocity.y; wb_y=_mm_loadu_ps(t);
	  for (int j=0;j<4;j++) t[j]=bodies[b->body_b[j]].angular_velocity.z; wb_z=_mm_loadu_ps(t);
	}
	__m128 zero = _mm_setzero_ps();
	__m128 linear_vn = SOA_DOT3(_mm_sub_ps(vb_x,va_x),_mm_sub_ps(vb_y,va_y),_mm_sub_ps(vb_z,va_z), b->normal_x,b->normal_y,b->normal_z);
	__m128 total_lambda_n = zero;

	// Normal contacts: read directly from pre-built contact layers (zero per-iteration gather).
	for (int cp = 0; cp < b->max_contacts; cp++) {
		PGS_ContactLayer4* cl = &b->cp[cp];
		__m128 vn = _mm_add_ps(linear_vn, _mm_sub_ps(SOA_DOT3(wb_x,wb_y,wb_z,cl->rn_b_x,cl->rn_b_y,cl->rn_b_z), SOA_DOT3(wa_x,wa_y,wa_z,cl->rn_a_x,cl->rn_a_y,cl->rn_a_z)));
		__m128 lambda_n = _mm_mul_ps(cl->eff_mass_n, _mm_sub_ps(_mm_sub_ps(zero, _mm_add_ps(vn, _mm_add_ps(cl->bias, cl->bounce))), _mm_mul_ps(cl->softness, cl->lambda_n)));
		__m128 new_lam = _mm_max_ps(_mm_add_ps(cl->lambda_n, lambda_n), zero);
		__m128 delta = _mm_sub_ps(new_lam, cl->lambda_n);
		cl->lambda_n = new_lam; // persists for next iteration
		__m128 da=_mm_mul_ps(delta,b->inv_mass_a), db=_mm_mul_ps(delta,b->inv_mass_b);
		va_x=_mm_sub_ps(va_x,_mm_mul_ps(b->normal_x,da)); va_y=_mm_sub_ps(va_y,_mm_mul_ps(b->normal_y,da)); va_z=_mm_sub_ps(va_z,_mm_mul_ps(b->normal_z,da));
		vb_x=_mm_add_ps(vb_x,_mm_mul_ps(b->normal_x,db)); vb_y=_mm_add_ps(vb_y,_mm_mul_ps(b->normal_y,db)); vb_z=_mm_add_ps(vb_z,_mm_mul_ps(b->normal_z,db));
		wa_x=_mm_sub_ps(wa_x,_mm_mul_ps(cl->wn_a_x,delta)); wa_y=_mm_sub_ps(wa_y,_mm_mul_ps(cl->wn_a_y,delta)); wa_z=_mm_sub_ps(wa_z,_mm_mul_ps(cl->wn_a_z,delta));
		wb_x=_mm_add_ps(wb_x,_mm_mul_ps(cl->wn_b_x,delta)); wb_y=_mm_add_ps(wb_y,_mm_mul_ps(cl->wn_b_y,delta)); wb_z=_mm_add_ps(wb_z,_mm_mul_ps(cl->wn_b_z,delta));
		linear_vn = _mm_add_ps(linear_vn, _mm_mul_ps(delta, b->inv_mass_sum));
		total_lambda_n = _mm_add_ps(total_lambda_n, new_lam);
	}

	// Jacobi friction: compute all three rows from same velocity snapshot, apply once.
	__m128 max_f = _mm_mul_ps(b->friction, total_lambda_n);
	__m128 dvx=_mm_sub_ps(vb_x,va_x),dvy=_mm_sub_ps(vb_y,va_y),dvz=_mm_sub_ps(vb_z,va_z);
	// Snapshot angular velocities for all friction rows
	__m128 wa_sx=wa_x,wa_sy=wa_y,wa_sz=wa_z, wb_sx=wb_x,wb_sy=wb_y,wb_sz=wb_z;

	__m128 vt1=_mm_add_ps(SOA_DOT3(dvx,dvy,dvz,b->tangent1_x,b->tangent1_y,b->tangent1_z),_mm_sub_ps(SOA_DOT3(wb_sx,wb_sy,wb_sz,b->rct1_b_x,b->rct1_b_y,b->rct1_b_z),SOA_DOT3(wa_sx,wa_sy,wa_sz,b->rct1_a_x,b->rct1_a_y,b->rct1_a_z)));
	__m128 nt1=_mm_max_ps(_mm_sub_ps(zero,max_f),_mm_min_ps(_mm_add_ps(b->lambda_t1,_mm_mul_ps(b->eff_mass_t1,_mm_sub_ps(zero,vt1))),max_f));
	__m128 dt1=_mm_sub_ps(nt1,b->lambda_t1); b->lambda_t1=nt1;

	__m128 vt2=_mm_add_ps(SOA_DOT3(dvx,dvy,dvz,b->tangent2_x,b->tangent2_y,b->tangent2_z),_mm_sub_ps(SOA_DOT3(wb_sx,wb_sy,wb_sz,b->rct2_b_x,b->rct2_b_y,b->rct2_b_z),SOA_DOT3(wa_sx,wa_sy,wa_sz,b->rct2_a_x,b->rct2_a_y,b->rct2_a_z)));
	__m128 nt2=_mm_max_ps(_mm_sub_ps(zero,max_f),_mm_min_ps(_mm_add_ps(b->lambda_t2,_mm_mul_ps(b->eff_mass_t2,_mm_sub_ps(zero,vt2))),max_f));
	__m128 dt2=_mm_sub_ps(nt2,b->lambda_t2); b->lambda_t2=nt2;

	__m128 max_tw=_mm_mul_ps(_mm_mul_ps(b->friction,total_lambda_n),b->patch_radius);
	__m128 wrel=SOA_DOT3(_mm_sub_ps(wb_sx,wa_sx),_mm_sub_ps(wb_sy,wa_sy),_mm_sub_ps(wb_sz,wa_sz),b->normal_x,b->normal_y,b->normal_z);
	__m128 ntw=_mm_max_ps(_mm_sub_ps(zero,max_tw),_mm_min_ps(_mm_add_ps(b->lambda_twist,_mm_mul_ps(b->eff_mass_twist,_mm_sub_ps(zero,wrel))),max_tw));
	__m128 dtw=_mm_sub_ps(ntw,b->lambda_twist); b->lambda_twist=ntw;

	// Single combined apply for all friction
	__m128 lin_x=_mm_add_ps(_mm_mul_ps(b->tangent1_x,dt1),_mm_mul_ps(b->tangent2_x,dt2));
	__m128 lin_y=_mm_add_ps(_mm_mul_ps(b->tangent1_y,dt1),_mm_mul_ps(b->tangent2_y,dt2));
	__m128 lin_z=_mm_add_ps(_mm_mul_ps(b->tangent1_z,dt1),_mm_mul_ps(b->tangent2_z,dt2));
	va_x=_mm_sub_ps(va_x,_mm_mul_ps(lin_x,b->inv_mass_a)); va_y=_mm_sub_ps(va_y,_mm_mul_ps(lin_y,b->inv_mass_a)); va_z=_mm_sub_ps(va_z,_mm_mul_ps(lin_z,b->inv_mass_a));
	vb_x=_mm_add_ps(vb_x,_mm_mul_ps(lin_x,b->inv_mass_b)); vb_y=_mm_add_ps(vb_y,_mm_mul_ps(lin_y,b->inv_mass_b)); vb_z=_mm_add_ps(vb_z,_mm_mul_ps(lin_z,b->inv_mass_b));
	wa_x=_mm_sub_ps(wa_x,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_a_x,dt1),_mm_mul_ps(b->w_t2_a_x,dt2)),_mm_mul_ps(b->w_tw_a_x,dtw)));
	wa_y=_mm_sub_ps(wa_y,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_a_y,dt1),_mm_mul_ps(b->w_t2_a_y,dt2)),_mm_mul_ps(b->w_tw_a_y,dtw)));
	wa_z=_mm_sub_ps(wa_z,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_a_z,dt1),_mm_mul_ps(b->w_t2_a_z,dt2)),_mm_mul_ps(b->w_tw_a_z,dtw)));
	wb_x=_mm_add_ps(wb_x,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_b_x,dt1),_mm_mul_ps(b->w_t2_b_x,dt2)),_mm_mul_ps(b->w_tw_b_x,dtw)));
	wb_y=_mm_add_ps(wb_y,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_b_y,dt1),_mm_mul_ps(b->w_t2_b_y,dt2)),_mm_mul_ps(b->w_tw_b_y,dtw)));
	wb_z=_mm_add_ps(wb_z,_mm_add_ps(_mm_add_ps(_mm_mul_ps(b->w_t1_b_z,dt1),_mm_mul_ps(b->w_t2_b_z,dt2)),_mm_mul_ps(b->w_tw_b_z,dtw)));

	// Scatter velocities back
	{ float t[4];
	  _mm_storeu_ps(t, va_x); for(int j=0;j<4;j++) bodies[b->body_a[j]].velocity.x=t[j];
	  _mm_storeu_ps(t, va_y); for(int j=0;j<4;j++) bodies[b->body_a[j]].velocity.y=t[j];
	  _mm_storeu_ps(t, va_z); for(int j=0;j<4;j++) bodies[b->body_a[j]].velocity.z=t[j];
	  _mm_storeu_ps(t, wa_x); for(int j=0;j<4;j++) bodies[b->body_a[j]].angular_velocity.x=t[j];
	  _mm_storeu_ps(t, wa_y); for(int j=0;j<4;j++) bodies[b->body_a[j]].angular_velocity.y=t[j];
	  _mm_storeu_ps(t, wa_z); for(int j=0;j<4;j++) bodies[b->body_a[j]].angular_velocity.z=t[j];
	  _mm_storeu_ps(t, vb_x); for(int j=0;j<4;j++) bodies[b->body_b[j]].velocity.x=t[j];
	  _mm_storeu_ps(t, vb_y); for(int j=0;j<4;j++) bodies[b->body_b[j]].velocity.y=t[j];
	  _mm_storeu_ps(t, vb_z); for(int j=0;j<4;j++) bodies[b->body_b[j]].velocity.z=t[j];
	  _mm_storeu_ps(t, wb_x); for(int j=0;j<4;j++) bodies[b->body_b[j]].angular_velocity.x=t[j];
	  _mm_storeu_ps(t, wb_y); for(int j=0;j<4;j++) bodies[b->body_b[j]].angular_velocity.y=t[j];
	  _mm_storeu_ps(t, wb_z); for(int j=0;j<4;j++) bodies[b->body_b[j]].angular_velocity.z=t[j];
	}
}

#endif // SIMD_SSE
