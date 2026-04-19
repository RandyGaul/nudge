// pre_solve_manifold.inc -- extracted per-manifold pre_solve for parallel dispatch.
// Included from solver_pgs.c. Each call writes to fixed-stride slots in sm/sc arrays.
// Thread-safe: warm_cache is read-only, each manifold writes to its own index range.

static void pre_solve_manifold(WorldInternal* w, InternalManifold* im, int idx, SolverManifold* sm, SolverContact* sc, float dt, float soft_dd, float bias_dd, float soft_ds, float bias_ds)
{
	BodyHot* a = &w->body_hot[im->body_a];
	BodyHot* b = &w->body_hot[im->body_b];
	BodyState* sa = &w->body_state[im->body_a];
	BodyState* sb = &w->body_state[im->body_b];
	float inv_mass_sum = a->inv_mass + b->inv_mass;
	if (inv_mass_sum == 0.0f) return; // sm[idx] stays zeroed (contact_count=0)

	float mu = sqrtf(sa->friction * sb->friction);
	float mu_roll = sqrtf(sa->rolling_friction * sb->rolling_friction);
	float rest = sa->restitution > sb->restitution ? sa->restitution : sb->restitution;
	WarmManifold* wm = im->warm; // cached from narrowphase (avoids duplicate hash lookup)
	int contact_base = idx * MAX_CONTACTS;

	SolverManifold smf = { .body_a = im->body_a, .body_b = im->body_b, .sub_id = im->sub_id, .contact_start = contact_base, .contact_count = im->m.count, .friction = mu, .rolling_friction = mu_roll, .inv_mass_a = a->inv_mass, .inv_mass_b = b->inv_mass };

	// Skip per-contact restitution velocity when both bodies are slow (settled stacks).
	float rest_thresh2 = SOLVER_RESTITUTION_THRESH * SOLVER_RESTITUTION_THRESH;
	int check_bounce = rest > 0.0f && (len2(a->velocity) + len2(a->angular_velocity) + len2(b->velocity) + len2(b->angular_velocity)) > rest_thresh2;

	// Select precomputed softness parameters based on static/dynamic pair type.
	int has_static = (a->inv_mass == 0.0f || b->inv_mass == 0.0f);
	float soft_val = has_static ? soft_ds : soft_dd;
	float bias_rate_val = has_static ? bias_ds : bias_dd;

	// Per-contact setup. Only initialize fields used by patch-mode PGS and post-solve.
	for (int c = 0; c < im->m.count; c++) {
		Contact* ct = &im->m.contacts[c];
		v3 r_a = sub(ct->point, sa->position);
		v3 r_b = sub(ct->point, sb->position);
		v3 rn_a = cross(r_a, ct->normal);
		v3 rn_b = cross(r_b, ct->normal);
		v3 w_n_a = inv_inertia_world_mul(a, rn_a);
		v3 w_n_b = inv_inertia_world_mul(b, rn_b);
		float K_n = inv_mass_sum + dot(w_n_a, rn_a) + dot(w_n_b, rn_b);
		float eff_mass_n = K_n > 1e-12f ? 1.0f / K_n : 0.0f;
		float softness = 0.0f, bias = 0.0f;
		if (soft_val > 0.0f) {
			float K = eff_mass_n > 0.0f ? 1.0f / eff_mass_n : 0.0f;
			// Bepu-style softened constraint, mass-aware. With extraDt = soft_val
			// (= 1/den) the constraint has:
			//   softness_impulse_scale_bepu = K / (K + den)     (dimensionless, <= 1)
			//   eff_mass_softened           = den / (K*(K+den))
			// Our PGS update has the form
			//     delta_lam = eff_mass * (-(vn+bias) - softness * lam_old)
			// i.e. the softness*lam_old term is multiplied by eff_mass. Matching
			// Bepu's decay factor requires eff_mass * softness == softness_bepu,
			// so softness = softness_bepu / eff_mass = K * Ks (Ks = K*soft_val).
			// Product eff_mass*softness = Ks/(1+Ks) <= 1 — stable for all masses.
			// Prior formulation used softness = Ks/(1+Ks) directly, making the
			// product Ks * Ks/(1+Ks)^2 / K = soft_val/(K*(1+soft_val)^2). For heavy
			// bodies this blew up to ~5.5 and left deep-pen contacts oscillating
			// between full and zero lambda every PGS iteration.
			float Ks = K * soft_val;
			float denom = 1.0f + Ks;
			softness = K * Ks / denom;
			eff_mass_n = (K * denom) > 1e-12f ? 1.0f / (K * denom) : 0.0f;
			float pen = ct->penetration - SOLVER_SLOP;
			bias = pen > 0.0f ? -bias_rate_val * pen : 0.0f;
			if (bias < -w->max_push_velocity) bias = -w->max_push_velocity;
		}
		float bounce = 0.0f;
		if (check_bounce) {
			v3 vel_a = add(a->velocity, cross(a->angular_velocity, r_a));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, r_b));
			float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
			bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;
			if (bounce != 0.0f) bias = 0.0f;
		}
		SolverContact* s = &sc[contact_base + c];
		s->r_a = r_a; s->r_b = r_b; s->normal = ct->normal; s->penetration = ct->penetration;
		s->feature_id = ct->feature_id; s->eff_mass_n = eff_mass_n;
		s->softness = softness; s->bias = bias; s->bounce = bounce;
		s->lambda_n = 0;
	}

	// Patch friction: centroid, tangent basis, effective mass
	v3 centroid_a = V3(0,0,0), centroid_b = V3(0,0,0);
	for (int c = 0; c < smf.contact_count; c++) { centroid_a = add(centroid_a, sc[contact_base+c].r_a); centroid_b = add(centroid_b, sc[contact_base+c].r_b); }
	float inv_n = 1.0f / (float)smf.contact_count;
	smf.centroid_r_a = scale(centroid_a, inv_n);
	smf.centroid_r_b = scale(centroid_b, inv_n);
	v3 n = sc[contact_base].normal;
	smf.normal = n;
	contact_tangent_basis(n, &smf.tangent1, &smf.tangent2);
	smf.patch_area = estimate_patch_area(im->m.contacts, im->m.count);
	smf.patch_radius = 0.6667f * sqrtf(smf.patch_area * (1.0f / 3.14159265f));
	// Compute cross(r, dir) and I_w * cross(r, dir) for effective mass computation.
	// Not stored in SolverManifold -- both scalar and SIMD solve paths recompute inline.
	v3 rct1_a = cross(smf.centroid_r_a, smf.tangent1), rct1_b = cross(smf.centroid_r_b, smf.tangent1);
	v3 rct2_a = cross(smf.centroid_r_a, smf.tangent2), rct2_b = cross(smf.centroid_r_b, smf.tangent2);
	v3 wt1_a = inv_inertia_world_mul(a, rct1_a), wt1_b = inv_inertia_world_mul(b, rct1_b);
	v3 wt2_a = inv_inertia_world_mul(a, rct2_a), wt2_b = inv_inertia_world_mul(b, rct2_b);
	v3 wtw_a = inv_inertia_world_mul(a, n), wtw_b = inv_inertia_world_mul(b, n);
	float k_t1 = inv_mass_sum + dot(wt1_a, rct1_a) + dot(wt1_b, rct1_b);
	smf.eff_mass_t1 = k_t1 > 1e-12f ? 1.0f / k_t1 : 0.0f;
	float k_t2 = inv_mass_sum + dot(wt2_a, rct2_a) + dot(wt2_b, rct2_b);
	smf.eff_mass_t2 = k_t2 > 1e-12f ? 1.0f / k_t2 : 0.0f;
	float k_twist = dot(wtw_a, n) + dot(wtw_b, n);
	smf.eff_mass_twist = k_twist > 1e-12f ? 1.0f / k_twist : 0.0f;

	// Rolling friction: one scalar row opposing ω_tangent (angular velocity in
	// the contact plane — the component that corresponds to rolling). Axis is
	// snapshotted here and held for the whole step; PGS treats it as a fixed
	// motor driving ω·axis toward 0, capped by mu_roll * lambda_n_total.
	// Skip cheaply when mu_roll == 0 or ω_tangent is effectively zero.
	smf.rolling_axis = V3(0, 0, 0);
	smf.eff_mass_roll = 0.0f;
	if (mu_roll > 0.0f) {
		v3 omega_rel = sub(b->angular_velocity, a->angular_velocity);
		v3 omega_tan = sub(omega_rel, scale(n, dot(omega_rel, n)));
		float l2 = dot(omega_tan, omega_tan);
		if (l2 > 1e-8f) {
			v3 axis = scale(omega_tan, 1.0f / sqrtf(l2));
			v3 iw_axis_a = inv_inertia_world_mul(a, axis);
			v3 iw_axis_b = inv_inertia_world_mul(b, axis);
			float k_roll = dot(iw_axis_a, axis) + dot(iw_axis_b, axis);
			if (k_roll > 1e-12f) {
				smf.rolling_axis = axis;
				smf.eff_mass_roll = 1.0f / k_roll;
			}
		}
	}

	// Warm start matching
	if (wm && w->warm_start_enabled) {
		int old_matched[MAX_CONTACTS] = {0}, new_matched[MAX_CONTACTS] = {0};
		for (int c = 0; c < smf.contact_count; c++) { SolverContact* s = &sc[contact_base+c]; if (s->feature_id == 0) continue; int match = warm_match(wm, s->feature_id); if (match >= 0) { s->lambda_n = wm->contacts[match].lambda_n; old_matched[match] = 1; new_matched[c] = 1; } }
		for (int c = 0; c < smf.contact_count; c++) { if (new_matched[c]) continue; SolverContact* s = &sc[contact_base+c]; float best_d2 = 0.01f; int best = -1; for (int j = 0; j < wm->count; j++) { if (old_matched[j]) continue; float d2 = len2(sub(s->r_a, wm->contacts[j].r_a)); if (d2 < best_d2) { best_d2 = d2; best = j; } } if (best >= 0) { s->lambda_n = wm->contacts[best].lambda_n; old_matched[best] = 1; new_matched[c] = 1; } }
		int new_unmatched = 0; for (int c = 0; c < smf.contact_count; c++) if (!new_matched[c]) new_unmatched++;
		if (new_unmatched > 0) { float lo_n=0; for (int j=0;j<wm->count;j++) if (!old_matched[j]) lo_n+=wm->contacts[j].lambda_n; float share=1.0f/(float)new_unmatched; for (int c=0;c<smf.contact_count;c++) { if(new_matched[c])continue; sc[contact_base+c].lambda_n=lo_n*share; } }
		smf.lambda_t1 = wm->manifold_lambda_t1; smf.lambda_t2 = wm->manifold_lambda_t2; smf.lambda_twist = wm->manifold_lambda_twist;
		// Rolling lambda only survives warm-start if the cached axis still aligns
		// with this frame's snapshot (ω_tangent direction hasn't flipped). A
		// mostly-aligned axis preserves the impulse with a projection; near-flip
		// zeroes it so we don't inject backwards energy.
		if (smf.eff_mass_roll > 0.0f) {
			float axis_dot = dot(wm->manifold_rolling_axis, smf.rolling_axis);
			smf.lambda_roll = (axis_dot > 0.5f) ? wm->manifold_lambda_roll * axis_dot : 0.0f;
		}
	}

	// Speculative/shallow contact lambda clamp: zero for negative penetration,
	// and cap warm-started lambda for near-zero penetration to prevent overshoot.
	// Without this, a stale warm lambda from a deep-pen frame gets applied to
	// a shallow contact, injecting energy and causing erratic bouncing.
	for (int c = 0; c < smf.contact_count; c++) {
		SolverContact* s = &sc[contact_base+c];
		if (s->penetration <= 0.0f) { s->lambda_n = 0; continue; }
		// Cap warm lambda: a shallow contact (pen near slop) shouldn't carry
		// more impulse than what a steady resting contact would produce.
		// Steady-state lambda ~= mass * gravity * dt / n_contacts ~= 0.08.
		// Scale: allow full lambda when pen >= 2*slop, taper to 0 at pen=0.
		float pen_scale = s->penetration / (2.0f * SOLVER_SLOP);
		if (pen_scale > 1.0f) pen_scale = 1.0f;
		s->lambda_n *= pen_scale;
	}

	sm[idx] = smf;
}
