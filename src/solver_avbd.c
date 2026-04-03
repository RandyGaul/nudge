// solver_avbd.c -- Augmented Vertex Block Descent solver.
// Primal-dual position solver: per-body 6x6 Newton steps + augmented Lagrangian dual updates.
// Reference: Ly, Narain et al. "Augmented VBD" SIGGRAPH 2025; Chris Giles' 3D demo.

static int avbd_body_is_sleeping(WorldInternal* w, int idx)
{
	int isl = w->body_cold[idx].island_id;
	return isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake;
}

// Adaptive alpha: scale stabilization by constraint error magnitude.
// Small error (drift): alpha_eff ≈ alpha (gentle correction, prevents jitter).
// Large error (violated): alpha_eff → 0 (fast snap to satisfied).
static float avbd_adaptive_alpha(float alpha, float err_magnitude)
{
	return alpha * fminf(1.0f, AVBD_STABLE_THRESH / fmaxf(err_magnitude, 1e-6f));
}

// Build separate CSR adjacency for contacts and joints (no type tag = no branch in inner loop).
static void avbd_build_adjacency(WorldInternal* w, AVBD_Manifold* am, int am_count, int body_count, int** out_ct_start, AVBD_ContactAdj** out_ct_adj, int** out_jt_start, AVBD_JointAdj** out_jt_adj)
{
	// --- Contact adjacency ---
	CK_DYNA int* ct_start = NULL;
	afit_set(ct_start, body_count + 1);
	memset(ct_start, 0, (body_count + 1) * sizeof(int));

	int ct_total = 0;
	for (int i = 0; i < am_count; i++)
		for (int c = 0; c < am[i].contact_count; c++) {
			ct_start[am[i].body_a]++;
			ct_start[am[i].body_b]++;
			ct_total += 2;
		}

	int sum = 0;
	for (int i = 0; i <= body_count; i++) { int tmp = ct_start[i]; ct_start[i] = sum; sum += tmp; }

	CK_DYNA AVBD_ContactAdj* ct_adj = NULL;
	if (ct_total > 0) {
		afit_set(ct_adj, ct_total);
		CK_DYNA int* off = NULL; afit_set(off, body_count);
		for (int i = 0; i < body_count; i++) off[i] = ct_start[i];
		for (int i = 0; i < am_count; i++)
			for (int c = 0; c < am[i].contact_count; c++) {
				ct_adj[off[am[i].body_a]++] = (AVBD_ContactAdj){ i, c, 1 };
				ct_adj[off[am[i].body_b]++] = (AVBD_ContactAdj){ i, c, 0 };
			}
		afree(off);
	}

	// --- Joint adjacency ---
	CK_DYNA int* jt_start = NULL;
	afit_set(jt_start, body_count + 1);
	memset(jt_start, 0, (body_count + 1) * sizeof(int));

	int jt_total = 0;
	int jt_count = asize(w->joints);
	for (int i = 0; i < jt_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		jt_start[w->joints[i].body_a]++;
		jt_start[w->joints[i].body_b]++;
		jt_total += 2;
	}

	sum = 0;
	for (int i = 0; i <= body_count; i++) { int tmp = jt_start[i]; jt_start[i] = sum; sum += tmp; }

	CK_DYNA AVBD_JointAdj* jt_adj = NULL;
	if (jt_total > 0) {
		afit_set(jt_adj, jt_total);
		CK_DYNA int* off = NULL; afit_set(off, body_count);
		for (int i = 0; i < body_count; i++) off[i] = jt_start[i];
		for (int i = 0; i < jt_count; i++) {
			if (!split_alive(w->joint_gen, i)) continue;
			jt_adj[off[w->joints[i].body_a]++] = (AVBD_JointAdj){ i, 1 };
			jt_adj[off[w->joints[i].body_b]++] = (AVBD_JointAdj){ i, 0 };
		}
		afree(off);
	}

	*out_ct_start = ct_start; *out_ct_adj = ct_adj;
	*out_jt_start = jt_start; *out_jt_adj = jt_adj;
}

// Build AVBD manifolds from collision output. Warm-start lambda/penalty from cache.
static void avbd_pre_solve(WorldInternal* w, InternalManifold* manifolds, int count, AVBD_Manifold** out_am, float dt)
{
	CK_DYNA AVBD_Manifold* am = NULL;
	float alpha = w->avbd_alpha;
	float gamma = w->avbd_gamma;

	for (int i = 0; i < count; i++) {
		InternalManifold* im = &manifolds[i];
		if (im->m.count == 0) continue;

		AVBD_Manifold m = {0};
		m.body_a = im->body_a;
		m.body_b = im->body_b;
		m.contact_count = im->m.count;

		BodyHot* a = &w->body_hot[m.body_a];
		BodyHot* b = &w->body_hot[m.body_b];
		m.friction = sqrtf(a->friction * b->friction);

		v3 n = neg(im->m.contacts[0].normal);
		v3 t1, t2;
		contact_tangent_basis(n, &t1, &t2);
		m.basis = (m3x3){{ n.x, n.y, n.z, t1.x, t1.y, t1.z, t2.x, t2.y, t2.z }};

		uint64_t key = body_pair_key(m.body_a, m.body_b);
		AVBD_WarmManifold* wm = map_get_ptr(w->avbd_warm_cache, key);

		for (int c = 0; c < m.contact_count; c++) {
			Contact* ct = &im->m.contacts[c];
			AVBD_Contact* ac = &m.contacts[c];

			v3 xA_world = add(ct->point, scale(ct->normal, ct->penetration));
			v3 xB_world = ct->point;
			ac->r_a = rotate(inv(a->rotation), sub(xA_world, a->position));
			ac->r_b = rotate(inv(b->rotation), sub(xB_world, b->position));
			ac->feature_id = ct->feature_id;

			v3 diff = sub(xA_world, xB_world);
			ac->C0 = V3(dot(n, diff) + AVBD_MARGIN, dot(t1, diff), dot(t2, diff));

			// Inertia-scaled initial penalty so first-frame forces are meaningful
			float inv_m = fmaxf(a->inv_mass, b->inv_mass);
			float pen_floor = inv_m > 0.0f ? 0.1f / (inv_m * dt * dt) : AVBD_PENALTY_MIN;
			ac->penalty = V3(pen_floor, pen_floor, pen_floor);
			ac->lambda = V3(0, 0, 0);
			ac->stick = 0;

			if (wm) {
				int matched = -1;
				for (int j = 0; j < wm->count && matched < 0; j++) {
					if (ac->feature_id != 0 && ac->feature_id == wm->contacts[j].feature_id)
						matched = j;
				}
				if (matched < 0) {
					float best_d2 = 0.01f;
					for (int j = 0; j < wm->count; j++) {
						float d2 = len2(sub(ac->r_a, wm->contacts[j].r_a));
						if (d2 < best_d2) { best_d2 = d2; matched = j; }
					}
				}
				if (matched >= 0) {
					AVBD_WarmContact* wc = &wm->contacts[matched];
					ac->penalty = wc->penalty;
					ac->lambda = wc->lambda;
					ac->stick = wc->stick;
					// Static friction: preserve contact geometry from previous frame
					// so the contact anchor doesn't slide (Giles manifold.cpp pattern).
					if (wc->stick) {
						ac->r_a = wc->r_a;
						ac->r_b = wc->r_b;
					}
				}
				ac->lambda = scale(ac->lambda, alpha * gamma);
				ac->penalty.x = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.x * gamma, AVBD_PENALTY_MAX));
				ac->penalty.y = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.y * gamma, AVBD_PENALTY_MAX));
				ac->penalty.z = fmaxf(AVBD_PENALTY_MIN, fminf(ac->penalty.z * gamma, AVBD_PENALTY_MAX));
			}
		}
		apush(am, m);
	}
	*out_am = am;
}

static void avbd_joints_pre_solve(WorldInternal* w, float alpha, float gamma, float dt)
{
	float dt2 = dt * dt;
	int count = asize(w->joints);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			v3 anchor_a = add(a->position, rotate(a->rotation, j->ball_socket.local_a));
			v3 anchor_b = add(b->position, rotate(b->rotation, j->ball_socket.local_b));
			j->avbd_C0_lin = sub(anchor_a, anchor_b);
		} else {
			v3 anchor_a = add(a->position, rotate(a->rotation, j->distance.local_a));
			v3 anchor_b = add(b->position, rotate(b->rotation, j->distance.local_b));
			v3 d = sub(anchor_a, anchor_b);
			float dist = len(d);
			float err = dist - j->distance.rest_length;
			v3 axis = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
			j->avbd_C0_lin = scale(axis, err);
		}

		// Inertia-scaled initial penalty: use lighter body's inertia as baseline.
		// Ensures penalty is immediately meaningful vs body mass.
		float inv_m = fmaxf(a->inv_mass, b->inv_mass); // lighter body dominates
		float inertia_floor = inv_m > 0.0f ? 0.1f / (inv_m * dt2) : AVBD_PENALTY_MIN;

		// Clamp penalty to material stiffness for spring joints (prevents over-stiffening).
		// Rigid joints (frequency=0) have infinite stiffness — no cap.
		float spring_freq = j->type == JOINT_BALL_SOCKET ? j->ball_socket.spring.frequency : j->distance.spring.frequency;
		float stiffness_cap = AVBD_PENALTY_MAX;
		if (spring_freq > 0.0f) {
			float omega = 2.0f * 3.14159265f * spring_freq;
			stiffness_cap = omega * omega;
		}

		j->avbd_lambda_lin = scale(j->avbd_lambda_lin, alpha * gamma);
		j->avbd_penalty_lin.x = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.x * gamma, stiffness_cap));
		j->avbd_penalty_lin.y = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.y * gamma, stiffness_cap));
		j->avbd_penalty_lin.z = fmaxf(inertia_floor, fminf(j->avbd_penalty_lin.z * gamma, stiffness_cap));
	}
}

static void avbd_post_solve(WorldInternal* w, AVBD_Manifold* am, int am_count)
{
	for (int i = 0; i < am_count; i++) {
		AVBD_Manifold* m = &am[i];
		uint64_t key = body_pair_key(m->body_a, m->body_b);
		AVBD_WarmManifold wm = {0};
		wm.count = m->contact_count;
		wm.stale = 0;
		for (int c = 0; c < m->contact_count; c++) {
			AVBD_Contact* ac = &m->contacts[c];
			wm.contacts[c] = (AVBD_WarmContact){
				.feature_id = ac->feature_id, .r_a = ac->r_a, .r_b = ac->r_b,
				.penalty = ac->penalty, .lambda = ac->lambda, .stick = ac->stick,
			};
		}
		map_set(w->avbd_warm_cache, key, wm);
	}
}

static void avbd_warm_cache_age_and_evict(WorldInternal* w)
{
	for (int i = 0; i < map_size(w->avbd_warm_cache); i++)
		w->avbd_warm_cache[i].stale++;
	int i = 0;
	while (i < map_size(w->avbd_warm_cache)) {
		if (w->avbd_warm_cache[i].stale > 1)
			map_del(w->avbd_warm_cache, map_keys(w->avbd_warm_cache)[i]);
		else i++;
	}
}

// --- Primal step: per-body 6x6 Newton solve ---

// Stamp a contact into the per-body 6x6 system.
static void avbd_stamp_contact(AVBD_Manifold* m, AVBD_Contact* ac, int is_body_a, AVBD_BodyState* states, float alpha, BodyHot* ha, BodyHot* hb, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 dpA = sub(ha->position, states[m->body_a].initial_lin);
	v3 daA = quat_sub_angular(ha->rotation, states[m->body_a].initial_ang);
	v3 dpB = sub(hb->position, states[m->body_b].initial_lin);
	v3 daB = quat_sub_angular(hb->rotation, states[m->body_b].initial_ang);

	v3 rAW = rotate(ha->rotation, ac->r_a);
	v3 rBW = rotate(hb->rotation, ac->r_b);

	v3 b0 = V3(m->basis.m[0], m->basis.m[1], m->basis.m[2]);
	v3 b1 = V3(m->basis.m[3], m->basis.m[4], m->basis.m[5]);
	v3 b2 = V3(m->basis.m[6], m->basis.m[7], m->basis.m[8]);

	v3 jALin[3] = { b0, b1, b2 };
	v3 jBLin[3] = { neg(b0), neg(b1), neg(b2) };
	v3 jAAng[3] = { cross(rAW, b0), cross(rAW, b1), cross(rAW, b2) };
	v3 jBAng[3] = { cross(rBW, jBLin[0]), cross(rBW, jBLin[1]), cross(rBW, jBLin[2]) };

	float C[3];
	C[0] = ac->C0.x*(1-alpha) + dot(jALin[0],dpA) + dot(jBLin[0],dpB) + dot(jAAng[0],daA) + dot(jBAng[0],daB);
	C[1] = ac->C0.y*(1-alpha) + dot(jALin[1],dpA) + dot(jBLin[1],dpB) + dot(jAAng[1],daA) + dot(jBAng[1],daB);
	C[2] = ac->C0.z*(1-alpha) + dot(jALin[2],dpA) + dot(jBLin[2],dpB) + dot(jAAng[2],daA) + dot(jBAng[2],daB);

	float F[3];
	F[0] = ac->penalty.x * C[0] + ac->lambda.x;
	F[1] = ac->penalty.y * C[1] + ac->lambda.y;
	F[2] = ac->penalty.z * C[2] + ac->lambda.z;

	if (F[0] > 0.0f) F[0] = 0.0f;

	float bounds = fabsf(F[0]) * m->friction;
	float fric_len = sqrtf(F[1]*F[1] + F[2]*F[2]);
	if (fric_len > bounds && fric_len > 0.0f) {
		float s = bounds / fric_len;
		F[1] *= s; F[2] *= s;
	}

	v3* jL = is_body_a ? jALin : jBLin;
	v3* jA = is_body_a ? jAAng : jBAng;

	for (int r = 0; r < 3; r++) {
		float kk = (&ac->penalty.x)[r];
		*rhsLin = add(*rhsLin, scale(jL[r], F[r]));
		*rhsAng = add(*rhsAng, scale(jA[r], F[r]));
		for (int ri = 0; ri < 3; ri++)
			for (int ci = 0; ci < 3; ci++) {
				lhsLin->m[ri*3+ci] += (&jL[r].x)[ri] * (&jL[r].x)[ci] * kk;
				lhsAng->m[ri*3+ci] += (&jA[r].x)[ri] * (&jA[r].x)[ci] * kk;
				lhsCross->m[ri*3+ci] += (&jA[r].x)[ri] * (&jL[r].x)[ci] * kk;
			}
	}

	// Geometric stiffness for contacts (Sec 3.5): diagonal approximation of
	// higher-order Hessian. Accounts for contact offset rotating with body.
	{
		v3 r = is_body_a ? rAW : neg(rBW);
		float H[9] = {0};
		for (int k = 0; k < 3; k++) {
			float fk = F[k];
			for (int ri = 0; ri < 3; ri++)
				for (int ci = 0; ci < 3; ci++) {
					float val = (ri == ci ? -(&r.x)[k] : 0.0f) + (ci == k ? (&r.x)[ri] : 0.0f);
					H[ri*3+ci] += val * fk;
				}
		}
		for (int ci = 0; ci < 3; ci++) {
			float col_len = sqrtf(H[0*3+ci]*H[0*3+ci] + H[1*3+ci]*H[1*3+ci] + H[2*3+ci]*H[2*3+ci]);
			lhsAng->m[ci*3+ci] += col_len;
		}
	}
}

// Stamp a ball-socket joint into the per-body 6x6 system.
static void avbd_stamp_ball_socket(JointInternal* j, BodyHot* a, BodyHot* b, int is_body_a, float alpha, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 rAW = rotate(a->rotation, j->ball_socket.local_a);
	v3 rBW = rotate(b->rotation, j->ball_socket.local_b);

	v3 C_vec = sub(add(a->position, rAW), add(b->position, rBW));
	if (j->ball_socket.spring.frequency <= 0.0f) {
		float aa = avbd_adaptive_alpha(alpha, len(j->avbd_C0_lin));
		C_vec = sub(C_vec, scale(j->avbd_C0_lin, aa));
	}

	v3 K = j->avbd_penalty_lin;
	v3 F = add(hmul(K, C_vec), j->avbd_lambda_lin);

	// Jacobians: J_lin_a = I, J_lin_b = -I, J_ang_a = skew(-rAW), J_ang_b = skew(rBW)
	float sign = is_body_a ? 1.0f : -1.0f;
	v3 r_world = is_body_a ? rAW : rBW;
	m3x3 jAng = is_body_a ? skew(neg(rAW)) : skew(rBW);

	// Stamp LHS: I^T * diag(K) * I = diag(K) for linear, jAng^T * diag(K) * jAng for angular
	for (int ri = 0; ri < 3; ri++)
		lhsLin->m[ri*3+ri] += (&K.x)[ri]; // diag(K) since J_lin = +/-I

	m3x3 jAngT = transpose(jAng);
	m3x3 Kmat = diag(K);
	m3x3 jAngTK = mul(jAngT, Kmat);
	*lhsAng = add(*lhsAng, mul(jAngTK, jAng));
	*lhsCross = add(*lhsCross, scale(jAngTK, sign));

	// Geometric stiffness: diagonal approximation of higher-order Hessian term (Sec 3.5).
	// Accounts for the nonlinear coupling between rotation and anchor position.
	// Without this, the angular Newton step overshoots for joints under load.
	{
		v3 r = is_body_a ? rAW : neg(rBW);
		// geometricStiffnessBallSocket(k, r) builds a matrix M where M = diag(-r[k]) + e_k * r^T
		// Sum over k weighted by F[k], then diagonalize (column norms on diagonal).
		float H[9] = {0};
		for (int k = 0; k < 3; k++) {
			float fk = (&F.x)[k];
			for (int ri = 0; ri < 3; ri++)
				for (int ci = 0; ci < 3; ci++) {
					float val = (ri == ci ? -(&r.x)[k] : 0.0f) + (ci == k ? (&r.x)[ri] : 0.0f);
					H[ri*3+ci] += val * fk;
				}
		}
		// Diagonalize: column norms
		for (int ci = 0; ci < 3; ci++) {
			float col_len = sqrtf(H[0*3+ci]*H[0*3+ci] + H[1*3+ci]*H[1*3+ci] + H[2*3+ci]*H[2*3+ci]);
			lhsAng->m[ci*3+ci] += col_len;
		}
	}

	// Stamp RHS: J^T * F
	*rhsLin = add(*rhsLin, scale(F, sign));
	*rhsAng = add(*rhsAng, m3x3_mul_v3(jAngT, F));
}

// Stamp a distance joint into the per-body 6x6 system.
static void avbd_stamp_distance(JointInternal* j, BodyHot* a, BodyHot* b, int is_body_a, float alpha, m3x3* lhsLin, m3x3* lhsAng, m3x3* lhsCross, v3* rhsLin, v3* rhsAng)
{
	v3 rAW = rotate(a->rotation, j->distance.local_a);
	v3 rBW = rotate(b->rotation, j->distance.local_b);
	v3 d = sub(add(a->position, rAW), add(b->position, rBW));
	float dist = len(d);
	v3 axis = dist > 1e-6f ? scale(d, 1.0f / dist) : V3(0, 1, 0);
	float err = dist - j->distance.rest_length;

	float k = j->avbd_penalty_lin.x;
	float f_scalar = k * err + j->avbd_lambda_lin.x;

	if (j->distance.spring.frequency <= 0.0f) {
		float C0_scalar = len(j->avbd_C0_lin);
		if (dot(j->avbd_C0_lin, axis) < 0) C0_scalar = -C0_scalar;
		float aa = avbd_adaptive_alpha(alpha, fabsf(C0_scalar));
		f_scalar = k * (err - aa * C0_scalar) + j->avbd_lambda_lin.x;
	}

	v3 jLin_vec = is_body_a ? axis : neg(axis);
	v3 r_world = is_body_a ? rAW : rBW;
	v3 jAng_vec = cross(r_world, jLin_vec);

	for (int ri = 0; ri < 3; ri++)
		for (int ci = 0; ci < 3; ci++) {
			lhsLin->m[ri*3+ci] += (&jLin_vec.x)[ri] * (&jLin_vec.x)[ci] * k;
			lhsAng->m[ri*3+ci] += (&jAng_vec.x)[ri] * (&jAng_vec.x)[ci] * k;
			lhsCross->m[ri*3+ci] += (&jAng_vec.x)[ri] * (&jLin_vec.x)[ci] * k;
		}

	v3 f_lin = scale(jLin_vec, f_scalar);
	v3 f_ang = cross(r_world, f_lin);
	*rhsLin = add(*rhsLin, f_lin);
	*rhsAng = add(*rhsAng, f_ang);
}

static void avbd_primal_step(WorldInternal* w, AVBD_Manifold* am, int am_count, int* ct_start, AVBD_ContactAdj* ct_adj, int* jt_start, AVBD_JointAdj* jt_adj, AVBD_BodyState* states, float alpha, float dt)
{
	int body_count = asize(w->body_hot);
	float dt2 = dt * dt;

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;

		float mass = 1.0f / h->inv_mass;
		v3 moment = rcp(h->inv_inertia_local);

		m3x3 lhsLin = diag(V3(mass/dt2, mass/dt2, mass/dt2));
		m3x3 lhsAng = diag(V3(moment.x/dt2, moment.y/dt2, moment.z/dt2));
		m3x3 lhsCross = {0};

		v3 dp = sub(h->position, states[i].inertial_lin);
		v3 da = quat_sub_angular(h->rotation, states[i].inertial_ang);
		v3 rhsLin = V3(mass/dt2 * dp.x, mass/dt2 * dp.y, mass/dt2 * dp.z);
		v3 rhsAng = V3(moment.x/dt2 * da.x, moment.y/dt2 * da.y, moment.z/dt2 * da.z);

		// Pass 1: contacts (tight loop, no branching)
		for (int j = ct_start[i]; j < ct_start[i + 1]; j++) {
			AVBD_ContactAdj* ca = &ct_adj[j];
			AVBD_Manifold* m = &am[ca->manifold_idx];
			avbd_stamp_contact(m, &m->contacts[ca->contact_idx], ca->is_body_a, states, alpha,
			                   &w->body_hot[m->body_a], &w->body_hot[m->body_b],
			                   &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
		}

		// Pass 2: joints (tight loop, no branching)
		for (int j = jt_start[i]; j < jt_start[i + 1]; j++) {
			AVBD_JointAdj* ja = &jt_adj[j];
			JointInternal* jt = &w->joints[ja->joint_idx];
			BodyHot* a = &w->body_hot[jt->body_a];
			BodyHot* b = &w->body_hot[jt->body_b];
			if (jt->type == JOINT_BALL_SOCKET)
				avbd_stamp_ball_socket(jt, a, b, ja->is_body_a, alpha, &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
			else
				avbd_stamp_distance(jt, a, b, ja->is_body_a, alpha, &lhsLin, &lhsAng, &lhsCross, &rhsLin, &rhsAng);
		}

		v3 dxLin, dxAng;
		solve_6x6_ldl(lhsLin, lhsAng, lhsCross, neg(rhsLin), neg(rhsAng), &dxLin, &dxAng);
		h->position = add(h->position, dxLin);
		h->rotation = quat_add_dw(h->rotation, dxAng);
	}
}

// --- Dual steps ---

static void avbd_contacts_dual_step(WorldInternal* w, AVBD_Manifold* am, int am_count, AVBD_BodyState* states, float alpha)
{
	float beta_lin = w->avbd_beta_lin;

	for (int i = 0; i < am_count; i++) {
		AVBD_Manifold* m = &am[i];
		BodyHot* ha = &w->body_hot[m->body_a];
		BodyHot* hb = &w->body_hot[m->body_b];

		v3 dpA = sub(ha->position, states[m->body_a].initial_lin);
		v3 daA = quat_sub_angular(ha->rotation, states[m->body_a].initial_ang);
		v3 dpB = sub(hb->position, states[m->body_b].initial_lin);
		v3 daB = quat_sub_angular(hb->rotation, states[m->body_b].initial_ang);

		v3 b0 = V3(m->basis.m[0], m->basis.m[1], m->basis.m[2]);
		v3 b1 = V3(m->basis.m[3], m->basis.m[4], m->basis.m[5]);
		v3 b2 = V3(m->basis.m[6], m->basis.m[7], m->basis.m[8]);

		for (int c = 0; c < m->contact_count; c++) {
			AVBD_Contact* ac = &m->contacts[c];
			v3 rAW = rotate(ha->rotation, ac->r_a);
			v3 rBW = rotate(hb->rotation, ac->r_b);

			v3 nb0 = neg(b0), nb1 = neg(b1), nb2 = neg(b2);
			v3 jAAng[3] = { cross(rAW, b0), cross(rAW, b1), cross(rAW, b2) };
			v3 jBAng[3] = { cross(rBW, nb0), cross(rBW, nb1), cross(rBW, nb2) };

			float Cn = ac->C0.x*(1-alpha) + dot(b0,dpA) + dot(nb0,dpB) + dot(jAAng[0],daA) + dot(jBAng[0],daB);
			float Ct1 = ac->C0.y*(1-alpha) + dot(b1,dpA) + dot(nb1,dpB) + dot(jAAng[1],daA) + dot(jBAng[1],daB);
			float Ct2 = ac->C0.z*(1-alpha) + dot(b2,dpA) + dot(nb2,dpB) + dot(jAAng[2],daA) + dot(jBAng[2],daB);

			float Fn = ac->penalty.x * Cn + ac->lambda.x;
			float Ft1 = ac->penalty.y * Ct1 + ac->lambda.y;
			float Ft2 = ac->penalty.z * Ct2 + ac->lambda.z;
			if (Fn > 0.0f) Fn = 0.0f;

			float bounds = fabsf(Fn) * m->friction;
			float fric_len = sqrtf(Ft1*Ft1 + Ft2*Ft2);
			if (fric_len > bounds && fric_len > 0.0f) { float s = bounds/fric_len; Ft1 *= s; Ft2 *= s; }

			ac->lambda = V3(Fn, Ft1, Ft2);

			if (Fn < 0.0f)
				ac->penalty.x = fminf(ac->penalty.x + beta_lin * fabsf(Cn), AVBD_PENALTY_MAX);
			if (fric_len <= bounds) {
				ac->penalty.y = fminf(ac->penalty.y + beta_lin * fabsf(Ct1), AVBD_PENALTY_MAX);
				ac->penalty.z = fminf(ac->penalty.z + beta_lin * fabsf(Ct2), AVBD_PENALTY_MAX);
				ac->stick = sqrtf(Ct1*Ct1 + Ct2*Ct2) < AVBD_STICK_THRESH;
			}
		}
	}
}

static void avbd_joints_dual_step(WorldInternal* w, float alpha)
{
	float beta_lin = w->avbd_beta_lin;
	int count = asize(w->joints);

	for (int i = 0; i < count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			// Giles pattern: C is modified in-place by stabilization.
			// Lambda update and penalty ramp both use the stabilized C.
			v3 C = sub(add(a->position, rotate(a->rotation, j->ball_socket.local_a)),
			           add(b->position, rotate(b->rotation, j->ball_socket.local_b)));
			if (j->ball_socket.spring.frequency <= 0.0f) {
				float aa = avbd_adaptive_alpha(alpha, len(j->avbd_C0_lin));
				C = sub(C, scale(j->avbd_C0_lin, aa));
			}
			v3 F = add(hmul(j->avbd_penalty_lin, C), j->avbd_lambda_lin);
			j->avbd_lambda_lin = F;

			v3 absC = V3(fabsf(C.x), fabsf(C.y), fabsf(C.z));
			float cap = AVBD_PENALTY_MAX;
			if (j->ball_socket.spring.frequency > 0.0f) {
				float om = 2.0f * 3.14159265f * j->ball_socket.spring.frequency;
				cap = om * om;
			}
			j->avbd_penalty_lin = V3(
				fminf(j->avbd_penalty_lin.x + beta_lin * absC.x, cap),
				fminf(j->avbd_penalty_lin.y + beta_lin * absC.y, cap),
				fminf(j->avbd_penalty_lin.z + beta_lin * absC.z, cap));
		} else {
			v3 d = sub(add(a->position, rotate(a->rotation, j->distance.local_a)),
			          add(b->position, rotate(b->rotation, j->distance.local_b)));
			float dist = len(d);
			float err = dist - j->distance.rest_length;
			// In-place stabilization (matching Giles pattern)
			if (j->distance.spring.frequency <= 0.0f) {
				float C0_scalar = len(j->avbd_C0_lin);
				v3 ax = dist > 1e-6f ? scale(d, 1.0f/dist) : V3(0,1,0);
				if (dot(j->avbd_C0_lin, ax) < 0) C0_scalar = -C0_scalar;
				float aa = avbd_adaptive_alpha(alpha, fabsf(C0_scalar));
				err -= aa * C0_scalar;
			}
			float F = j->avbd_penalty_lin.x * err + j->avbd_lambda_lin.x;
			j->avbd_lambda_lin.x = F;
			float cap = AVBD_PENALTY_MAX;
			if (j->distance.spring.frequency > 0.0f) {
				float om = 2.0f * 3.14159265f * j->distance.spring.frequency;
				cap = om * om;
			}
			j->avbd_penalty_lin.x = fminf(j->avbd_penalty_lin.x + beta_lin * fabsf(err), cap);
		}
	}
}

// --- Top-level AVBD solver ---

static void avbd_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, float dt)
{
	int body_count = asize(w->body_hot);
	float alpha = w->avbd_alpha;
	float gamma = w->avbd_gamma;

	avbd_warm_cache_age_and_evict(w);

	CK_DYNA AVBD_Manifold* am = NULL;
	avbd_pre_solve(w, manifolds, manifold_count, &am, dt);
	int am_count = asize(am);

	avbd_joints_pre_solve(w, alpha, gamma, dt);

	int *ct_start = NULL, *jt_start = NULL;
	AVBD_ContactAdj* ct_adj = NULL;
	AVBD_JointAdj* jt_adj = NULL;
	avbd_build_adjacency(w, am, am_count, body_count, &ct_start, &ct_adj, &jt_start, &jt_adj);

	// Ensure prev_velocity array is big enough
	if (asize(w->avbd_prev_velocity) < body_count) {
		afit(w->avbd_prev_velocity, body_count);
		asetlen(w->avbd_prev_velocity, body_count);
	}

	CK_DYNA AVBD_BodyState* states = NULL;
	afit_set(states, body_count);
	memset(states, 0, body_count * sizeof(AVBD_BodyState));

	float inv_dt = 1.0f / dt;
	v3 grav = w->gravity;
	float grav_len = len(grav);

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		states[i].initial_lin = h->position;
		states[i].initial_ang = h->rotation;

		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;

		// Inertial position: y = x + v*dt + g*dt^2
		states[i].inertial_lin = add(add(h->position, scale(h->velocity, dt)), scale(grav, dt * dt));
		states[i].inertial_ang = quat_add_dw(h->rotation, scale(h->angular_velocity, dt));

		// Adaptive body warm-start (VBD paper): use acceleration from previous
		// frame to predict gravity contribution, giving solver a better starting guess.
		v3 prev_v = w->avbd_prev_velocity[i];
		v3 accel = scale(sub(h->velocity, prev_v), inv_dt);
		float accel_along_grav = grav_len > 0.0f ? dot(accel, scale(grav, 1.0f / grav_len)) : 0.0f;
		float accel_weight = grav_len > 0.0f ? fmaxf(0.0f, fminf(accel_along_grav / grav_len, 1.0f)) : 0.0f;

		h->position = add(h->position, add(scale(h->velocity, dt), scale(grav, accel_weight * dt * dt)));
		h->rotation = quat_add_dw(h->rotation, scale(h->angular_velocity, dt));
	}

	for (int it = 0; it < w->avbd_iterations; it++) {
		avbd_primal_step(w, am, am_count, ct_start, ct_adj, jt_start, jt_adj, states, alpha, dt);
		avbd_contacts_dual_step(w, am, am_count, states, alpha);
		avbd_joints_dual_step(w, alpha);
	}

	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		if (avbd_body_is_sleeping(w, i)) continue;
		w->avbd_prev_velocity[i] = h->velocity; // save for next frame's adaptive warm-start
		h->velocity = scale(sub(h->position, states[i].initial_lin), inv_dt);
		h->angular_velocity = scale(quat_sub_angular(h->rotation, states[i].initial_ang), inv_dt);
	}

	avbd_post_solve(w, am, am_count);

	afree(states);
	afree(ct_adj); afree(ct_start);
	afree(jt_adj); afree(jt_start);
	afree(am);
}
