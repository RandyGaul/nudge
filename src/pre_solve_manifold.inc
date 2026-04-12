// pre_solve_manifold.inc -- extracted per-manifold pre_solve for parallel dispatch.
// Included from solver_pgs.c. Each call writes to fixed-stride slots in sm/sc/pc arrays.
// Thread-safe: warm_cache is read-only, each manifold writes to its own index range.

static void pre_solve_manifold(WorldInternal* w, InternalManifold* im, int idx, SolverManifold* sm, SolverContact* sc, PatchContact* pc, float dt)
{
	BodyHot* a = &w->body_hot[im->body_a];
	BodyHot* b = &w->body_hot[im->body_b];
	float inv_mass_sum = a->inv_mass + b->inv_mass;
	if (inv_mass_sum == 0.0f) return; // sm[idx] stays zeroed (contact_count=0)

	float mu = sqrtf(a->friction * b->friction);
	float rest = a->restitution > b->restitution ? a->restitution : b->restitution;
	WarmManifold* wm = im->warm; // cached from narrowphase (avoids duplicate hash lookup)
	int contact_base = idx * MAX_CONTACTS;

	SolverManifold smf = { .body_a = im->body_a, .body_b = im->body_b, .contact_start = contact_base, .contact_count = im->m.count, .friction = mu, .inv_mass_a = a->inv_mass, .inv_mass_b = b->inv_mass };

	// Softness (same for all contacts in this manifold)
	float soft_val = 0.0f, bias_rate_val = 0.0f;
	if (w->solver_type != SOLVER_SI) {
		float hertz = w->contact_hertz;
		if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
		float omega = 2.0f * 3.14159265f * hertz;
		float d = 2.0f * w->contact_damping_ratio * omega;
		float kk = omega * omega;
		float hd = dt * d, hhk = dt * dt * kk;
		float denom = hd + hhk;
		if (denom > 1e-12f) { soft_val = 1.0f / denom; bias_rate_val = dt * kk * soft_val; }
	}

	// Per-contact setup
	for (int c = 0; c < im->m.count; c++) {
		Contact* ct = &im->m.contacts[c];
		SolverContact s = {0};
		s.r_a = sub(ct->point, a->position);
		s.r_b = sub(ct->point, b->position);
		s.normal = ct->normal;
		s.penetration = ct->penetration;
		s.feature_id = ct->feature_id;
		s.rn_a = cross(s.r_a, s.normal);
		s.rn_b = cross(s.r_b, s.normal);
		s.w_n_a = inv_inertia_world_mul(a, s.rn_a);
		s.w_n_b = inv_inertia_world_mul(b, s.rn_b);
		float K_n = inv_mass_sum + dot(cross(s.w_n_a, s.r_a), s.normal) + dot(cross(s.w_n_b, s.r_b), s.normal);
		s.eff_mass_n = K_n > 1e-12f ? 1.0f / K_n : 0.0f;
		if (soft_val > 0.0f) {
			float K = s.eff_mass_n > 0.0f ? 1.0f / s.eff_mass_n : 0.0f;
			s.softness = soft_val;
			s.eff_mass_n = (K + soft_val) > 1e-12f ? 1.0f / (K + soft_val) : 0.0f;
			float pen = ct->penetration - SOLVER_SLOP;
			s.bias = pen > 0.0f ? -bias_rate_val * pen : 0.0f;
			if (s.bias < -w->max_push_velocity) s.bias = -w->max_push_velocity;
		}
		v3 vel_a = add(a->velocity, cross(a->angular_velocity, s.r_a));
		v3 vel_b = add(b->velocity, cross(b->angular_velocity, s.r_b));
		float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
		s.bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;
		if (s.bounce != 0.0f) s.bias = 0.0f;
		sc[contact_base + c] = s;
		pc[contact_base + c] = (PatchContact){ .rn_a = s.rn_a, .rn_b = s.rn_b, .w_n_a = s.w_n_a, .w_n_b = s.w_n_b, .eff_mass_n = s.eff_mass_n, .bias = s.bias, .bounce = s.bounce, .softness = s.softness, .lambda_n = s.lambda_n };
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
	smf.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent1);
	smf.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent2);
	smf.patch_area = estimate_patch_area(im->m.contacts, im->m.count);
	smf.patch_radius = 0.6667f * sqrtf(smf.patch_area * (1.0f / 3.14159265f));
	float k_twist = dot(inv_inertia_world_mul(a, n), n) + dot(inv_inertia_world_mul(b, n), n);
	smf.eff_mass_twist = k_twist > 1e-12f ? 1.0f / k_twist : 0.0f;
	smf.rct1_a = cross(smf.centroid_r_a, smf.tangent1); smf.rct1_b = cross(smf.centroid_r_b, smf.tangent1);
	smf.rct2_a = cross(smf.centroid_r_a, smf.tangent2); smf.rct2_b = cross(smf.centroid_r_b, smf.tangent2);
	smf.w_t1_a = inv_inertia_world_mul(a, smf.rct1_a); smf.w_t1_b = inv_inertia_world_mul(b, smf.rct1_b);
	smf.w_t2_a = inv_inertia_world_mul(a, smf.rct2_a); smf.w_t2_b = inv_inertia_world_mul(b, smf.rct2_b);
	smf.w_tw_a = inv_inertia_world_mul(a, n); smf.w_tw_b = inv_inertia_world_mul(b, n);

	// Warm start matching
	if (wm && w->warm_start_enabled) {
		int old_matched[MAX_CONTACTS] = {0}, new_matched[MAX_CONTACTS] = {0};
		for (int c = 0; c < smf.contact_count; c++) { SolverContact* s = &sc[contact_base+c]; if (s->feature_id == 0) continue; int match = warm_match(wm, s->feature_id); if (match >= 0) { s->lambda_n = wm->contacts[match].lambda_n; old_matched[match] = 1; new_matched[c] = 1; } }
		for (int c = 0; c < smf.contact_count; c++) { if (new_matched[c]) continue; SolverContact* s = &sc[contact_base+c]; float best_d2 = 0.01f; int best = -1; for (int j = 0; j < wm->count; j++) { if (old_matched[j]) continue; float d2 = len2(sub(s->r_a, wm->contacts[j].r_a)); if (d2 < best_d2) { best_d2 = d2; best = j; } } if (best >= 0) { s->lambda_n = wm->contacts[best].lambda_n; old_matched[best] = 1; new_matched[c] = 1; } }
		int new_unmatched = 0; for (int c = 0; c < smf.contact_count; c++) if (!new_matched[c]) new_unmatched++;
		if (new_unmatched > 0) { float lo_n=0; for (int j=0;j<wm->count;j++) if (!old_matched[j]) lo_n+=wm->contacts[j].lambda_n; float share=1.0f/(float)new_unmatched; for (int c=0;c<smf.contact_count;c++) { if(new_matched[c])continue; sc[contact_base+c].lambda_n=lo_n*share; } }
		smf.lambda_t1 = wm->manifold_lambda_t1; smf.lambda_t2 = wm->manifold_lambda_t2; smf.lambda_twist = wm->manifold_lambda_twist;
	}

	// Speculative contacts: zero warm-started impulse for negative-penetration contacts
	for (int c = 0; c < smf.contact_count; c++) { SolverContact* s = &sc[contact_base+c]; if (s->penetration < 0.0f) { s->lambda_n = 0; s->lambda_t1 = 0; s->lambda_t2 = 0; } }
	// Sync warm-started lambda_n to PatchContact
	for (int c = 0; c < smf.contact_count; c++) pc[contact_base+c].lambda_n = sc[contact_base+c].lambda_n;

	sm[idx] = smf;
}
