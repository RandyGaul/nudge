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
	float k = inv_mass_sum
		+ dot(cross(inv_inertia_world_mul(a, ra_x_d), r_a), dir)
		+ dot(cross(inv_inertia_world_mul(b, rb_x_d), r_b), dir);
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

static void solver_pre_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, float dt)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	float inv_dt = dt > 0.0f ? 1.0f / dt : 0.0f;

	for (int i = 0; i < manifold_count; i++) {
		InternalManifold* im = &manifolds[i];
		BodyHot* a = &w->body_hot[im->body_a];
		BodyHot* b = &w->body_hot[im->body_b];
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		if (inv_mass_sum == 0.0f) continue;

		float mu = sqrtf(a->friction * b->friction);
		float rest = a->restitution > b->restitution ? a->restitution : b->restitution;

		// Look up warm starting data
		uint64_t key = body_pair_key(im->body_a, im->body_b);
		WarmManifold* wm = map_get_ptr(w->warm_cache, key);

		SolverManifold smf = {
			.body_a = im->body_a,
			.body_b = im->body_b,
			.contact_start = asize(sc),
			.contact_count = im->m.count,
			.friction = mu,
			.inv_mass_a = a->inv_mass,
			.inv_mass_b = b->inv_mass,
		};

		int patch_mode = (w->friction_model == FRICTION_PATCH);

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			SolverContact s = {0};

			s.r_a = sub(ct->point, a->position);
			s.r_b = sub(ct->point, b->position);
			s.normal = ct->normal;
			s.penetration = ct->penetration;
			s.feature_id = ct->feature_id;

			s.eff_mass_n = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.normal);

			// Soft contact constraint (skip for hard SI — uses NGS position correction instead)
			if (w->solver_type != SOLVER_SI) {
				float hertz = w->contact_hertz;
				if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
				float omega = 2.0f * 3.14159265f * hertz;
				float d = 2.0f * w->contact_damping_ratio * omega;
				float k = omega * omega;
				float hd = dt * d, hhk = dt * dt * k;
				float denom = hd + hhk;
				if (denom > 1e-12f) {
					s.softness = 1.0f / denom;
					float bias_rate = dt * k * s.softness;
					float K = s.eff_mass_n > 0.0f ? 1.0f / s.eff_mass_n : 0.0f;
					s.eff_mass_n = (K + s.softness) > 1e-12f ? 1.0f / (K + s.softness) : 0.0f;
					float pen = ct->penetration - SOLVER_SLOP;
					s.bias = pen > 0.0f ? -bias_rate * pen : 0.0f;
					if (s.bias < -w->max_push_velocity) s.bias = -w->max_push_velocity;
				}
			}

			// Precompute cross(r, normal) and I_w * cross(r, normal) — used by both modes.
			s.rn_a = cross(s.r_a, s.normal);
			s.rn_b = cross(s.r_b, s.normal);
			s.w_n_a = inv_inertia_world_mul(a, s.rn_a);
			s.w_n_b = inv_inertia_world_mul(b, s.rn_b);

			if (!patch_mode) {
				// Coulomb mode: precompute additional tangent angular data
				contact_tangent_basis(ct->normal, &s.tangent1, &s.tangent2);
				s.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent1);
				s.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent2);
				s.w_t1_a = inv_inertia_world_mul(a, cross(s.r_a, s.tangent1));
				s.w_t1_b = inv_inertia_world_mul(b, cross(s.r_b, s.tangent1));
				s.w_t2_a = inv_inertia_world_mul(a, cross(s.r_a, s.tangent2));
				s.w_t2_b = inv_inertia_world_mul(b, cross(s.r_b, s.tangent2));
			}

			v3 vel_a = add(a->velocity, cross(a->angular_velocity, s.r_a));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, s.r_b));
			float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
			s.bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;

			// When restitution is active, the bounce velocity handles separation.
			// Disable the penetration bias to prevent energy injection (bias + bounce
			// are additive in the solver, so both active injects extra energy).
			if (s.bounce != 0.0f) s.bias = 0.0f;

			apush(sc, s);
		}

		// Patch friction: compute centroid, patch area, tangent basis at manifold level.
		if (patch_mode) {
			v3 centroid_a = V3(0, 0, 0), centroid_b = V3(0, 0, 0);
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				centroid_a = add(centroid_a, s->r_a);
				centroid_b = add(centroid_b, s->r_b);
			}
			float inv_n = 1.0f / (float)smf.contact_count;
			smf.centroid_r_a = scale(centroid_a, inv_n);
			smf.centroid_r_b = scale(centroid_b, inv_n);

			// Use first contact normal for tangent basis (all share same normal in a manifold)
			v3 n = sc[smf.contact_start].normal;
			smf.normal = n;
			contact_tangent_basis(n, &smf.tangent1, &smf.tangent2);
			smf.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent1);
			smf.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent2);
			smf.patch_area = estimate_patch_area(im->m.contacts, im->m.count);
			smf.patch_radius = 0.6667f * sqrtf(smf.patch_area * (1.0f / 3.14159265f));

			// Torsional friction effective mass: 1 / (n^T * I_a_inv * n + n^T * I_b_inv * n)
			float k_twist = dot(inv_inertia_world_mul(a, n), n) + dot(inv_inertia_world_mul(b, n), n);
			smf.eff_mass_twist = k_twist > 1e-12f ? 1.0f / k_twist : 0.0f;

			// Precompute cross(centroid_r, tangent) for dv dot product + angular impulse
			smf.rct1_a = cross(smf.centroid_r_a, smf.tangent1);
			smf.rct1_b = cross(smf.centroid_r_b, smf.tangent1);
			smf.rct2_a = cross(smf.centroid_r_a, smf.tangent2);
			smf.rct2_b = cross(smf.centroid_r_b, smf.tangent2);
			smf.w_t1_a = inv_inertia_world_mul(a, smf.rct1_a);
			smf.w_t1_b = inv_inertia_world_mul(b, smf.rct1_b);
			smf.w_t2_a = inv_inertia_world_mul(a, smf.rct2_a);
			smf.w_t2_b = inv_inertia_world_mul(b, smf.rct2_b);
			smf.w_tw_a = inv_inertia_world_mul(a, n);
			smf.w_tw_b = inv_inertia_world_mul(b, n);
		}

		// Warm start: match new contacts to cached contacts.
		// Pass 1: exact feature ID match.
		// Pass 2: spatial fallback -- closest unmatched old contact by r_a distance.
		// Pass 3: redistribute any remaining unmatched old impulse evenly.
		if (wm) {
			int old_matched[MAX_CONTACTS] = {0};
			int new_matched[MAX_CONTACTS] = {0};
			// Pass 1: feature ID
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				if (s->feature_id == 0) continue;
				int match = warm_match(wm, s->feature_id);
				if (match >= 0) {
					s->lambda_n = wm->contacts[match].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[match].lambda_t1;
						s->lambda_t2 = wm->contacts[match].lambda_t2;
					}
					old_matched[match] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 2: spatial fallback for unmatched contacts
			float spatial_tol2 = 0.01f;
			for (int c = 0; c < smf.contact_count; c++) {
				if (new_matched[c]) continue;
				SolverContact* s = &sc[smf.contact_start + c];
				float best_d2 = spatial_tol2;
				int best = -1;
				for (int j = 0; j < wm->count; j++) {
					if (old_matched[j]) continue;
					float d2 = len2(sub(s->r_a, wm->contacts[j].r_a));
					if (d2 < best_d2) { best_d2 = d2; best = j; }
				}
				if (best >= 0) {
					s->lambda_n = wm->contacts[best].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[best].lambda_t1;
						s->lambda_t2 = wm->contacts[best].lambda_t2;
					}
					old_matched[best] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 3: redistribute remaining unmatched old impulse
			int new_unmatched = 0;
			for (int c = 0; c < smf.contact_count; c++)
				if (!new_matched[c]) new_unmatched++;
			if (new_unmatched > 0) {
				float leftover_n = 0, leftover_t1 = 0, leftover_t2 = 0;
				for (int j = 0; j < wm->count; j++) {
					if (!old_matched[j]) {
						leftover_n += wm->contacts[j].lambda_n;
						if (!patch_mode) {
							leftover_t1 += wm->contacts[j].lambda_t1;
							leftover_t2 += wm->contacts[j].lambda_t2;
						}
					}
				}
				float share = 1.0f / (float)new_unmatched;
				for (int c = 0; c < smf.contact_count; c++) {
					if (new_matched[c]) continue;
					SolverContact* s = &sc[smf.contact_start + c];
					s->lambda_n = leftover_n * share;
					if (!patch_mode) {
						s->lambda_t1 = leftover_t1 * share;
						s->lambda_t2 = leftover_t2 * share;
					}
				}
			}
			// Warm start manifold-level friction
			if (patch_mode) {
				smf.lambda_t1 = wm->manifold_lambda_t1;
				smf.lambda_t2 = wm->manifold_lambda_t2;
				smf.lambda_twist = wm->manifold_lambda_twist;
			}
		}

		// Speculative contacts (negative penetration, kept alive by margin) must
		// not carry warm-started impulse -- they exist for cache continuity only.
		for (int c = 0; c < smf.contact_count; c++) {
			SolverContact* s = &sc[smf.contact_start + c];
			if (s->penetration < 0.0f) { s->lambda_n = 0.0f; s->lambda_t1 = 0.0f; s->lambda_t2 = 0.0f; }
		}

		apush(sm, smf);
	}

	// Apply warm start impulses
	int patch_warm = (w->friction_model == FRICTION_PATCH);
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
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
	for (int i = 0; i < map_size(w->warm_cache); i++)
		w->warm_cache[i].stale++;
	int i = 0;
	while (i < map_size(w->warm_cache)) {
		if (w->warm_cache[i].stale > 1)
			map_del(w->warm_cache, map_keys(w->warm_cache)[i]);
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
static void solve_contact_patch_sv(SolverBodyVel* bodies, SolverManifold* m, SolverContact* sc)
{
	SolverBodyVel* a = &bodies[m->body_a];
	SolverBodyVel* b = &bodies[m->body_b];
	float ima = m->inv_mass_a, imb = m->inv_mass_b;

	v3 normal = m->normal;
	float inv_mass_sum = ima + imb;
	float linear_vn = dot(sub(b->velocity, a->velocity), normal);
	float total_lambda_n = 0.0f;
	for (int ci = 0; ci < m->contact_count; ci++) {
		SolverContact* s = &sc[m->contact_start + ci];
		// Use precomputed cross(r, normal) from pre_solve.
		v3 rn_a = s->rn_a;
		v3 rn_b = s->rn_b;
		float vn = linear_vn + dot(b->angular_velocity, rn_b) - dot(a->angular_velocity, rn_a);
		float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
		float old_n = s->lambda_n;
		s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
		float delta = s->lambda_n - old_n;
		// Apply: use precomputed angular impulse direction.
		v3 P = scale(normal, delta);
		a->velocity = sub(a->velocity, scale(P, ima));
		b->velocity = add(b->velocity, scale(P, imb));
		a->angular_velocity = sub(a->angular_velocity, scale(s->w_n_a, delta));
		b->angular_velocity = add(b->angular_velocity, scale(s->w_n_b, delta));
		linear_vn += delta * inv_mass_sum;
		total_lambda_n += s->lambda_n;
	}

	float max_f = m->friction * total_lambda_n;

	// Manifold-level tangent friction: use precomputed cross + angular impulse directions.
	v3 dv = sub(b->velocity, a->velocity);
	float vt1 = dot(dv, m->tangent1) + dot(b->angular_velocity, m->rct1_b) - dot(a->angular_velocity, m->rct1_a);
	float old_t1 = m->lambda_t1;
	m->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + m->eff_mass_t1 * (-vt1), max_f));
	float dt1 = m->lambda_t1 - old_t1;
	v3 P1 = scale(m->tangent1, dt1);
	a->velocity = sub(a->velocity, scale(P1, ima));
	b->velocity = add(b->velocity, scale(P1, imb));
	a->angular_velocity = sub(a->angular_velocity, scale(m->w_t1_a, dt1));
	b->angular_velocity = add(b->angular_velocity, scale(m->w_t1_b, dt1));

	dv = sub(b->velocity, a->velocity);
	float vt2 = dot(dv, m->tangent2) + dot(b->angular_velocity, m->rct2_b) - dot(a->angular_velocity, m->rct2_a);
	float old_t2 = m->lambda_t2;
	m->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + m->eff_mass_t2 * (-vt2), max_f));
	float dt2 = m->lambda_t2 - old_t2;
	v3 P2 = scale(m->tangent2, dt2);
	a->velocity = sub(a->velocity, scale(P2, ima));
	b->velocity = add(b->velocity, scale(P2, imb));
	a->angular_velocity = sub(a->angular_velocity, scale(m->w_t2_a, dt2));
	b->angular_velocity = add(b->angular_velocity, scale(m->w_t2_b, dt2));

	// Torsional friction: use precomputed angular impulse direction.
	float max_twist = m->friction * total_lambda_n * m->patch_radius;
	float w_rel = dot(sub(b->angular_velocity, a->angular_velocity), m->normal);
	float old_tw = m->lambda_twist;
	m->lambda_twist = fmaxf(-max_twist, fminf(old_tw + m->eff_mass_twist * (-w_rel), max_twist));
	float dtw = m->lambda_twist - old_tw;
	a->angular_velocity = sub(a->angular_velocity, scale(m->w_tw_a, dtw));
	b->angular_velocity = add(b->angular_velocity, scale(m->w_tw_b, dtw));
}
