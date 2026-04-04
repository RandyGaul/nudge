// joints.c -- joint constraint solvers (ball socket, distance)

// Symmetric 3x3 stored as 6 floats: [xx, xy, xz, yy, yz, zz].
static void sym3x3_inverse(const float* in, float* out)
{
	float a = in[0], b = in[1], c = in[2];
	float d = in[3], e = in[4], f = in[5];
	float det = a*(d*f - e*e) - b*(b*f - c*e) + c*(b*e - c*d);
	float inv_det = det != 0.0f ? 1.0f / det : 0.0f;
	out[0] = (d*f - e*e) * inv_det;
	out[1] = (c*e - b*f) * inv_det;
	out[2] = (b*e - c*d) * inv_det;
	out[3] = (a*f - c*c) * inv_det;
	out[4] = (b*c - a*e) * inv_det;
	out[5] = (a*d - b*b) * inv_det;
}

static v3 sym3x3_mul_v3(const float* m, v3 v)
{
	return V3(
		m[0]*v.x + m[1]*v.y + m[2]*v.z,
		m[1]*v.x + m[3]*v.y + m[4]*v.z,
		m[2]*v.x + m[4]*v.y + m[5]*v.z);
}

// Convert spring params to solver coefficients (bepu approach).
static void spring_compute(SpringParams sp, float dt, float* pos_to_vel, float* softness)
{
	if (sp.frequency <= 0.0f) {
		// Rigid constraint: Baumgarte stabilization with fractional correction.
		*pos_to_vel = dt > 0.0f ? SOLVER_BAUMGARTE / dt : 0.0f;
		*softness = 0.0f;
		return;
	}
	float omega = 2.0f * 3.14159265f * sp.frequency;
	float d = 2.0f * sp.damping_ratio * omega;
	float k = omega * omega;
	float hd = dt * d, hk = dt * k, hhk = dt * hk;
	float denom = hd + hhk;
	if (denom < 1e-12f) { *pos_to_vel = 0; *softness = 0; return; }
	*softness = 1.0f / denom;
	*pos_to_vel = hk * (*softness);
}

// Build ball socket effective mass (symmetric 3x3).
// K = (inv_mass_a + inv_mass_b)*I + skew(r_a)^T * I_a^-1 * skew(r_a) + same for B.
static void ball_socket_eff_mass(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, float softness, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	// Skew-sandwich: skew(r)^T * I^-1 * skew(r) for diagonal inertia.
	// For body A: I_a^-1 is diagonal in local space, r_a is world space.
	// Transform r_a to local, compute product, transform back.
	// Or equivalently use the explicit formula.
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m }; // start with linear term

	// Add angular contribution for body A
	v3 ia = a->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		// columns of world-space inverse inertia: R * diag(ia) * R^T
		// skew(r)^T * I_world^-1 * skew(r) computed column by column
		v3 e0 = inv_inertia_mul(a->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(a->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(a->rotation, ia, V3(-r_a.y, r_a.x, 0));
		// skew(r)^T * [e0 e1 e2] -- but skew(r)^T has rows [0 -rz ry; rz 0 -rx; -ry rx 0]
		// Result[i][j] = dot(skew_row_i(r), e_j)
		// Row 0 of skew^T: (0, -rz, ry)
		K[0] += -r_a.z*e0.y + r_a.y*e0.z; // xx
		K[1] += -r_a.z*e1.y + r_a.y*e1.z; // xy
		K[2] += -r_a.z*e2.y + r_a.y*e2.z; // xz
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z; // yy
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z; // yz
		K[5] += -r_a.y*e2.x + r_a.x*e2.y; // zz
	}

	v3 ib = b->inv_inertia_local;
	if (ib.x > 0 || ib.y > 0 || ib.z > 0) {
		v3 e0 = inv_inertia_mul(b->rotation, ib, V3(0, -r_b.z, r_b.y));
		v3 e1 = inv_inertia_mul(b->rotation, ib, V3(r_b.z, 0, -r_b.x));
		v3 e2 = inv_inertia_mul(b->rotation, ib, V3(-r_b.y, r_b.x, 0));
		K[0] += -r_b.z*e0.y + r_b.y*e0.z;
		K[1] += -r_b.z*e1.y + r_b.y*e1.z;
		K[2] += -r_b.z*e2.y + r_b.y*e2.z;
		K[3] +=  r_b.z*e1.x - r_b.x*e1.z;
		K[4] +=  r_b.z*e2.x - r_b.x*e2.z;
		K[5] += -r_b.y*e2.x + r_b.x*e2.y;
	}

	// Add softness to diagonal: K_eff = K + gamma * I
	K[0] += softness; K[3] += softness; K[5] += softness;
	sym3x3_inverse(K, out);
}

// Pre-solve joints: build solver arrays from persistent JointInternal data.
static void joints_pre_solve(WorldInternal* w, float dt, SolverBallSocket** out_bs, SolverDistance** out_dist)
{
	CK_DYNA SolverBallSocket* bs = NULL;
	CK_DYNA SolverDistance* dist = NULL;
	int joint_count = asize(w->joints);

	for (int i = 0; i < joint_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			SolverBallSocket s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->ball_socket.local_a);
			s.r_b = rotate(b->rotation, j->ball_socket.local_b);

			float ptv, soft;
			spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
			s.softness = soft;
			ball_socket_eff_mass(a, b, s.r_a, s.r_b, soft, s.eff_mass);

			// Position error: world anchor B - world anchor A
			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			s.bias = scale(sub(anchor_b, anchor_a), ptv);

			// Warm start from persistent storage
			s.lambda = j->warm_lambda3;

			apush(bs, s);
		} else if (j->type == JOINT_DISTANCE) {
			SolverDistance s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->distance.local_a);
			s.r_b = rotate(b->rotation, j->distance.local_b);

			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 delta = sub(anchor_b, anchor_a);
			float dist_val = len(delta);
			s.axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

			float ptv, soft;
			spring_compute(j->distance.spring, dt, &ptv, &soft);
			s.softness = soft;

			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s.r_a, s.axis)), s.r_a), s.axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s.r_b, s.axis)), s.r_b), s.axis);
			k += soft;
			s.eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;

			float error = dist_val - j->distance.rest_length;
			s.bias = -ptv * error;

			s.lambda = j->warm_lambda1;

			apush(dist, s);
		}
	}

	*out_bs = bs;
	*out_dist = dist;
}

// Apply warm start impulses for joints.
static void joints_warm_start(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &bs[i];
		if (s->lambda.x == 0 && s->lambda.y == 0 && s->lambda.z == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, s->lambda);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &dist[i];
		if (s->lambda == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, scale(s->axis, s->lambda));
	}
}

static void solve_ball_socket(WorldInternal* w, SolverBallSocket* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	v3 rhs = sub(neg(add(dv, s->bias)), scale(s->lambda, s->softness));
	v3 impulse = sym3x3_mul_v3(s->eff_mass, rhs);
	s->lambda = add(s->lambda, impulse);
	apply_impulse(a, b, s->r_a, s->r_b, impulse);
}

static void solve_distance(WorldInternal* w, SolverDistance* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	float cdot = dot(dv, s->axis);
	float lambda = s->eff_mass * (-cdot + s->bias - s->softness * s->lambda);
	s->lambda += lambda;
	apply_impulse(a, b, s->r_a, s->r_b, scale(s->axis, lambda));
}

// Split-impulse position correction for joints.
// Runs a pseudo-velocity PGS solve targeting position error, integrates positions
// with the corrective velocities, then restores real velocities. Prevents energy
// injection because pseudo-velocities never enter the real velocity state.
static void joints_split_impulse(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count, float dt, int iters)
{
	if (iters <= 0 || dt <= 0) return;
	int body_count = asize(w->body_hot);

	// Save real velocities and zero them for pseudo-velocity accumulation
	v3* saved_vel = CK_ALLOC(body_count * sizeof(v3));
	v3* saved_ang = CK_ALLOC(body_count * sizeof(v3));
	for (int i = 0; i < body_count; i++) {
		saved_vel[i] = w->body_hot[i].velocity;
		saved_ang[i] = w->body_hot[i].angular_velocity;
		w->body_hot[i].velocity = V3(0, 0, 0);
		w->body_hot[i].angular_velocity = V3(0, 0, 0);
	}

	// Save lambda and zero for fresh pseudo-velocity accumulation
	v3* saved_bs_lambda = CK_ALLOC(bs_count * sizeof(v3));
	float* saved_dist_lambda = CK_ALLOC(dist_count * sizeof(float));
	for (int i = 0; i < bs_count; i++) { saved_bs_lambda[i] = bs[i].lambda; bs[i].lambda = V3(0,0,0); }
	for (int i = 0; i < dist_count; i++) { saved_dist_lambda[i] = dist[i].lambda; dist[i].lambda = 0; }

	// Set bias from current position error and refresh effective mass
	float ptv = 0.1f / dt; // SI gain: 0.1 is stable for shattering; 0.15+ destabilizes hub_12
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &bs[i];
		if (s->softness != 0.0f) continue;
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		s->r_a = rotate(a->rotation, w->joints[s->joint_idx].ball_socket.local_a);
		s->r_b = rotate(b->rotation, w->joints[s->joint_idx].ball_socket.local_b);
		s->bias = scale(sub(add(b->position, s->r_b), add(a->position, s->r_a)), ptv);
		ball_socket_eff_mass(a, b, s->r_a, s->r_b, 0.0f, s->eff_mass);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &dist[i];
		if (s->softness != 0.0f) continue;
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		s->r_a = rotate(a->rotation, w->joints[s->joint_idx].distance.local_a);
		s->r_b = rotate(b->rotation, w->joints[s->joint_idx].distance.local_b);
		v3 anchor_a = add(a->position, s->r_a);
		v3 anchor_b = add(b->position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		float dist_val = len(delta);
		s->axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);
		s->bias = -ptv * (dist_val - w->joints[s->joint_idx].distance.rest_length);
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s->r_a, s->axis)), s->r_a), s->axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s->r_b, s->axis)), s->r_b), s->axis);
		s->eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;
	}

	// PGS iterations on joints (pseudo-velocity solve)
	for (int iter = 0; iter < iters; iter++) {
		for (int i = 0; i < bs_count; i++) if (bs[i].softness == 0.0f) solve_ball_socket(w, &bs[i]);
		for (int i = 0; i < dist_count; i++) if (dist[i].softness == 0.0f) solve_distance(w, &dist[i]);
	}

	// Integrate positions with pseudo-velocities
	for (int i = 0; i < body_count; i++) {
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f || !split_alive(w->body_gen, i)) continue;
		h->position = add(h->position, scale(h->velocity, dt));
		v3 ww = h->angular_velocity;
		if (ww.x != 0 || ww.y != 0 || ww.z != 0) {
			quat spin = { ww.x, ww.y, ww.z, 0.0f };
			quat dq = mul(spin, h->rotation);
			h->rotation.x += 0.5f * dt * dq.x; h->rotation.y += 0.5f * dt * dq.y;
			h->rotation.z += 0.5f * dt * dq.z; h->rotation.w += 0.5f * dt * dq.w;
			float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y + h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
			float inv_ql = 1.0f / (ql > 1e-15f ? ql : 1.0f);
			h->rotation.x *= inv_ql; h->rotation.y *= inv_ql; h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
		}
	}

	// Restore real velocities and lambda
	for (int i = 0; i < body_count; i++) {
		w->body_hot[i].velocity = saved_vel[i];
		w->body_hot[i].angular_velocity = saved_ang[i];
	}
	for (int i = 0; i < bs_count; i++) bs[i].lambda = saved_bs_lambda[i];
	for (int i = 0; i < dist_count; i++) dist[i].lambda = saved_dist_lambda[i];

	// Zero bias for main velocity solve
	for (int i = 0; i < bs_count; i++) if (bs[i].softness == 0.0f) bs[i].bias = V3(0, 0, 0);
	for (int i = 0; i < dist_count; i++) if (dist[i].softness == 0.0f) dist[i].bias = 0;

	CK_FREE(saved_vel);
	CK_FREE(saved_ang);
	CK_FREE(saved_bs_lambda);
	CK_FREE(saved_dist_lambda);
}

// NGS position correction for joints. Operates on positions only, no velocity modification.
// Mirrors solver_position_correct for contacts.
static void joints_position_correct(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count, int iters)
{
	for (int iter = 0; iter < iters; iter++) {
		for (int i = 0; i < bs_count; i++) {
			SolverBallSocket* s = &bs[i];
			if (s->softness != 0.0f) continue;
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 r_a = rotate(a->rotation, w->joints[s->joint_idx].ball_socket.local_a);
			v3 r_b = rotate(b->rotation, w->joints[s->joint_idx].ball_socket.local_b);
			v3 error = sub(add(b->position, r_b), add(a->position, r_a));
			float err_len = len(error);
			if (err_len < 1e-6f) continue;
			float inv_mass_sum = a->inv_mass + b->inv_mass;
			if (inv_mass_sum <= 0) continue;
			float correction = SOLVER_POS_BAUMGARTE * err_len;
			if (correction > SOLVER_POS_MAX_CORRECTION) correction = SOLVER_POS_MAX_CORRECTION;
			v3 n = scale(error, 1.0f / err_len);
			float P = correction / inv_mass_sum;
			a->position = add(a->position, scale(n, P * a->inv_mass));
			b->position = sub(b->position, scale(n, P * b->inv_mass));
		}
		for (int i = 0; i < dist_count; i++) {
			SolverDistance* s = &dist[i];
			if (s->softness != 0.0f) continue;
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 r_a = rotate(a->rotation, w->joints[s->joint_idx].distance.local_a);
			v3 r_b = rotate(b->rotation, w->joints[s->joint_idx].distance.local_b);
			v3 delta = sub(add(b->position, r_b), add(a->position, r_a));
			float dist_val = len(delta);
			if (dist_val < 1e-6f) continue;
			v3 axis = scale(delta, 1.0f / dist_val);
			float error = dist_val - w->joints[s->joint_idx].distance.rest_length;
			if (fabsf(error) < 1e-6f) continue;
			float inv_mass_sum = a->inv_mass + b->inv_mass;
			if (inv_mass_sum <= 0) continue;
			float correction = SOLVER_POS_BAUMGARTE * error;
			if (fabsf(correction) > SOLVER_POS_MAX_CORRECTION) correction = correction > 0 ? SOLVER_POS_MAX_CORRECTION : -SOLVER_POS_MAX_CORRECTION;
			float P = correction / inv_mass_sum;
			a->position = add(a->position, scale(axis, P * a->inv_mass));
			b->position = sub(b->position, scale(axis, P * b->inv_mass));
		}
	}
}

// Store joint accumulated impulses back to persistent storage.
static void joints_post_solve(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++)
		w->joints[bs[i].joint_idx].warm_lambda3 = bs[i].lambda;
	for (int i = 0; i < dist_count; i++)
		w->joints[dist[i].joint_idx].warm_lambda1 = dist[i].lambda;
	afree(bs);
	afree(dist);
}
