// joints.c -- joint constraint solvers (ball socket, distance, hinge)

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
		// Rigid constraint: no velocity bias. Position correction via NGS or LDL.
		*pos_to_vel = 0.0f;
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

static void hinge_tangent_basis(v3 axis, v3* t1, v3* t2)
{
	v3 ref = fabsf(axis.x) < 0.9f ? V3(1, 0, 0) : V3(0, 1, 0);
	*t1 = norm(cross(axis, ref));
	*t2 = cross(axis, *t1);
}

// Initialize lo/hi bounds for a SolverJoint (bilateral: unclamped).
static void solver_joint_init_bounds(SolverJoint* s)
{
	for (int d = 0; d < JOINT_MAX_DOF; d++) { s->lo[d] = -FLT_MAX; s->hi[d] = FLT_MAX; }
}

// Pre-solve joints: build unified SolverJoint array from persistent JointInternal data.
static void joints_pre_solve(WorldInternal* w, float dt, SolverJoint** out_joints)
{
	CK_DYNA SolverJoint* joints = NULL;
	int joint_count = asize(w->joints);

	for (int i = 0; i < joint_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			SolverJoint s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.type = JOINT_BALL_SOCKET;
			s.dof = 3;
			solver_joint_init_bounds(&s);
			s.r_a = rotate(a->rotation, j->ball_socket.local_a);
			s.r_b = rotate(b->rotation, j->ball_socket.local_b);

			float ptv, soft;
			spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
			s.softness = soft;
			ball_socket_eff_mass(a, b, s.r_a, s.r_b, soft, s.bs.eff_mass);

			// Position error: world anchor B - world anchor A
			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 bias_v = scale(sub(anchor_b, anchor_a), ptv);
			s.bias[0] = bias_v.x; s.bias[1] = bias_v.y; s.bias[2] = bias_v.z;

			// Warm start from persistent storage
			s.lambda[0] = j->warm_lambda3.x; s.lambda[1] = j->warm_lambda3.y; s.lambda[2] = j->warm_lambda3.z;

			apush(joints, s);
		} else if (j->type == JOINT_DISTANCE) {
			SolverJoint s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.type = JOINT_DISTANCE;
			s.dof = 1;
			solver_joint_init_bounds(&s);
			s.r_a = rotate(a->rotation, j->distance.local_a);
			s.r_b = rotate(b->rotation, j->distance.local_b);

			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 delta = sub(anchor_b, anchor_a);
			float dist_val = len(delta);
			s.dist.axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

			float ptv, soft;
			spring_compute(j->distance.spring, dt, &ptv, &soft);
			s.softness = soft;

			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s.r_a, s.dist.axis)), s.r_a), s.dist.axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s.r_b, s.dist.axis)), s.r_b), s.dist.axis);
			k += soft;
			s.dist.eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;

			float error = dist_val - j->distance.rest_length;
			s.bias[0] = -ptv * error;

			s.lambda[0] = j->warm_lambda1;

			apush(joints, s);
		} else if (j->type == JOINT_HINGE) {
			SolverJoint s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.type = JOINT_HINGE;
			s.dof = 5;
			solver_joint_init_bounds(&s);
			s.r_a = rotate(a->rotation, j->hinge.local_a);
			s.r_b = rotate(b->rotation, j->hinge.local_b);

			float ptv, soft;
			spring_compute(j->hinge.spring, dt, &ptv, &soft);
			s.softness = soft;
			ball_socket_eff_mass(a, b, s.r_a, s.r_b, soft, s.hinge.lin_eff_mass);

			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 lin_bias = scale(sub(anchor_b, anchor_a), ptv);
			s.bias[0] = lin_bias.x; s.bias[1] = lin_bias.y; s.bias[2] = lin_bias.z;

			v3 axis_a = norm(rotate(a->rotation, j->hinge.local_axis_a));
			s.hinge.axis_b = norm(rotate(b->rotation, j->hinge.local_axis_b));
			hinge_tangent_basis(axis_a, &s.hinge.t1, &s.hinge.t2);
			s.hinge.u1 = cross(s.hinge.t1, s.hinge.axis_b);
			s.hinge.u2 = cross(s.hinge.t2, s.hinge.axis_b);

			float k1 = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, s.hinge.u1), s.hinge.u1) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, s.hinge.u1), s.hinge.u1) + soft;
			s.hinge.ang_eff_mass[0] = k1 > 1e-12f ? 1.0f / k1 : 0.0f;
			float k2 = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, s.hinge.u2), s.hinge.u2) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, s.hinge.u2), s.hinge.u2) + soft;
			s.hinge.ang_eff_mass[1] = k2 > 1e-12f ? 1.0f / k2 : 0.0f;

			s.bias[3] = ptv * dot(s.hinge.t1, s.hinge.axis_b);
			s.bias[4] = ptv * dot(s.hinge.t2, s.hinge.axis_b);

			s.lambda[0] = j->warm_hinge.lin.x; s.lambda[1] = j->warm_hinge.lin.y; s.lambda[2] = j->warm_hinge.lin.z;
			s.lambda[3] = j->warm_hinge.ang[0];
			s.lambda[4] = j->warm_hinge.ang[1];

			apush(joints, s);
		}
	}

	*out_joints = joints;
}

// Apply warm start impulses for joints.
static void joints_warm_start(WorldInternal* w, SolverJoint* joints, int count)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		if (s->type == JOINT_BALL_SOCKET) {
			v3 lam = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
			if (lam.x == 0 && lam.y == 0 && lam.z == 0) continue;
			apply_impulse(a, b, s->r_a, s->r_b, lam);
		} else if (s->type == JOINT_DISTANCE) {
			if (s->lambda[0] == 0) continue;
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->dist.axis, s->lambda[0]));
		} else if (s->type == JOINT_HINGE) {
			v3 lin_lam = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
			if (lin_lam.x != 0 || lin_lam.y != 0 || lin_lam.z != 0)
				apply_impulse(a, b, s->r_a, s->r_b, lin_lam);
			for (int d = 0; d < 2; d++) {
				if (s->lambda[3 + d] == 0) continue;
				v3 u = d == 0 ? s->hinge.u1 : s->hinge.u2;
				v3 ang = scale(u, s->lambda[3 + d]);
				a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, ang));
				b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, ang));
			}
		}
	}
}

static void solve_ball_socket(WorldInternal* w, SolverJoint* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	v3 bias_v = V3(s->bias[0], s->bias[1], s->bias[2]);
	v3 lam = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
	v3 rhs = sub(neg(add(dv, bias_v)), scale(lam, s->softness));
	v3 impulse = sym3x3_mul_v3(s->bs.eff_mass, rhs);
	s->lambda[0] += impulse.x; s->lambda[1] += impulse.y; s->lambda[2] += impulse.z;
	apply_impulse(a, b, s->r_a, s->r_b, impulse);
}

static void solve_distance(WorldInternal* w, SolverJoint* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	float cdot = dot(dv, s->dist.axis);
	float lambda = s->dist.eff_mass * (-cdot + s->bias[0] - s->softness * s->lambda[0]);
	s->lambda[0] += lambda;
	apply_impulse(a, b, s->r_a, s->r_b, scale(s->dist.axis, lambda));
}

static void solve_hinge(WorldInternal* w, SolverJoint* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];

	// Linear part (same as ball socket)
	v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
	v3 lin_bias = V3(s->bias[0], s->bias[1], s->bias[2]);
	v3 lin_lam = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
	v3 rhs = sub(neg(add(dv, lin_bias)), scale(lin_lam, s->softness));
	v3 impulse = sym3x3_mul_v3(s->hinge.lin_eff_mass, rhs);
	s->lambda[0] += impulse.x; s->lambda[1] += impulse.y; s->lambda[2] += impulse.z;
	apply_impulse(a, b, s->r_a, s->r_b, impulse);

	// Angular DOF 1
	v3 dw = sub(a->angular_velocity, b->angular_velocity);
	float cdot1 = dot(dw, s->hinge.u1);
	float lambda1 = s->hinge.ang_eff_mass[0] * (-(cdot1 + s->bias[3]) - s->softness * s->lambda[3]);
	s->lambda[3] += lambda1;
	v3 ang1 = scale(s->hinge.u1, lambda1);
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, ang1));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, ang1));

	// Angular DOF 2 (re-read after DOF 1 impulse)
	dw = sub(a->angular_velocity, b->angular_velocity);
	float cdot2 = dot(dw, s->hinge.u2);
	float lambda2 = s->hinge.ang_eff_mass[1] * (-(cdot2 + s->bias[4]) - s->softness * s->lambda[4]);
	s->lambda[4] += lambda2;
	v3 ang2 = scale(s->hinge.u2, lambda2);
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, ang2));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, ang2));
}

// Split-impulse position correction for joints.
// Runs a pseudo-velocity PGS solve targeting position error, integrates positions
// with the corrective velocities, then restores real velocities. Prevents energy
// injection because pseudo-velocities never enter the real velocity state.
static void joints_split_impulse(WorldInternal* w, SolverJoint* joints, int joint_count, float dt, int iters)
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
	float* saved_lambda = CK_ALLOC(joint_count * JOINT_MAX_DOF * sizeof(float));
	for (int i = 0; i < joint_count; i++) {
		for (int d = 0; d < JOINT_MAX_DOF; d++) { saved_lambda[i * JOINT_MAX_DOF + d] = joints[i].lambda[d]; joints[i].lambda[d] = 0; }
	}

	// Set bias from current position error and refresh effective mass
	float ptv = 0.1f / dt; // SI gain: 0.1 is stable for shattering; 0.15+ destabilizes hub_12
	for (int i = 0; i < joint_count; i++) {
		SolverJoint* s = &joints[i];
		if (s->softness != 0.0f) continue;
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		if (s->type == JOINT_BALL_SOCKET) {
			s->r_a = rotate(a->rotation, w->joints[s->joint_idx].ball_socket.local_a);
			s->r_b = rotate(b->rotation, w->joints[s->joint_idx].ball_socket.local_b);
			v3 bias_v = scale(sub(add(b->position, s->r_b), add(a->position, s->r_a)), ptv);
			s->bias[0] = bias_v.x; s->bias[1] = bias_v.y; s->bias[2] = bias_v.z;
			ball_socket_eff_mass(a, b, s->r_a, s->r_b, 0.0f, s->bs.eff_mass);
		} else if (s->type == JOINT_DISTANCE) {
			s->r_a = rotate(a->rotation, w->joints[s->joint_idx].distance.local_a);
			s->r_b = rotate(b->rotation, w->joints[s->joint_idx].distance.local_b);
			v3 anchor_a = add(a->position, s->r_a);
			v3 anchor_b = add(b->position, s->r_b);
			v3 delta = sub(anchor_b, anchor_a);
			float dist_val = len(delta);
			s->dist.axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);
			s->bias[0] = -ptv * (dist_val - w->joints[s->joint_idx].distance.rest_length);
			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s->r_a, s->dist.axis)), s->r_a), s->dist.axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s->r_b, s->dist.axis)), s->r_b), s->dist.axis);
			s->dist.eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;
		} else if (s->type == JOINT_HINGE) {
			s->r_a = rotate(a->rotation, w->joints[s->joint_idx].hinge.local_a);
			s->r_b = rotate(b->rotation, w->joints[s->joint_idx].hinge.local_b);
			v3 lin_bias = scale(sub(add(b->position, s->r_b), add(a->position, s->r_a)), ptv);
			s->bias[0] = lin_bias.x; s->bias[1] = lin_bias.y; s->bias[2] = lin_bias.z;
			ball_socket_eff_mass(a, b, s->r_a, s->r_b, 0.0f, s->hinge.lin_eff_mass);
			v3 axis_a = norm(rotate(a->rotation, w->joints[s->joint_idx].hinge.local_axis_a));
			s->hinge.axis_b = norm(rotate(b->rotation, w->joints[s->joint_idx].hinge.local_axis_b));
			hinge_tangent_basis(axis_a, &s->hinge.t1, &s->hinge.t2);
			s->hinge.u1 = cross(s->hinge.t1, s->hinge.axis_b);
			s->hinge.u2 = cross(s->hinge.t2, s->hinge.axis_b);
			s->bias[3] = ptv * dot(s->hinge.t1, s->hinge.axis_b);
			s->bias[4] = ptv * dot(s->hinge.t2, s->hinge.axis_b);
			float k1 = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, s->hinge.u1), s->hinge.u1) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, s->hinge.u1), s->hinge.u1);
			s->hinge.ang_eff_mass[0] = k1 > 1e-12f ? 1.0f / k1 : 0.0f;
			float k2 = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, s->hinge.u2), s->hinge.u2) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, s->hinge.u2), s->hinge.u2);
			s->hinge.ang_eff_mass[1] = k2 > 1e-12f ? 1.0f / k2 : 0.0f;
		}
	}

	// PGS iterations on joints (pseudo-velocity solve)
	for (int iter = 0; iter < iters; iter++) {
		for (int i = 0; i < joint_count; i++) {
			if (joints[i].softness != 0.0f) continue;
			if (joints[i].type == JOINT_BALL_SOCKET) solve_ball_socket(w, &joints[i]);
			else if (joints[i].type == JOINT_DISTANCE) solve_distance(w, &joints[i]);
			else if (joints[i].type == JOINT_HINGE) solve_hinge(w, &joints[i]);
		}
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
	for (int i = 0; i < joint_count; i++) {
		for (int d = 0; d < JOINT_MAX_DOF; d++) joints[i].lambda[d] = saved_lambda[i * JOINT_MAX_DOF + d];
	}

	// Zero bias for main velocity solve
	for (int i = 0; i < joint_count; i++) {
		if (joints[i].softness != 0.0f) continue;
		for (int d = 0; d < joints[i].dof; d++) joints[i].bias[d] = 0;
	}

	CK_FREE(saved_vel);
	CK_FREE(saved_ang);
	CK_FREE(saved_lambda);
}

// NGS position correction for joints. Operates on positions only, no velocity modification.
// Mirrors solver_position_correct for contacts.
static void joints_position_correct(WorldInternal* w, SolverJoint* joints, int joint_count, int iters)
{
	for (int iter = 0; iter < iters; iter++) {
		for (int i = 0; i < joint_count; i++) {
			SolverJoint* s = &joints[i];
			if (s->softness != 0.0f) continue;
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			if (s->type == JOINT_BALL_SOCKET || s->type == JOINT_HINGE) {
				v3 local_a, local_b;
				if (s->type == JOINT_BALL_SOCKET) { local_a = w->joints[s->joint_idx].ball_socket.local_a; local_b = w->joints[s->joint_idx].ball_socket.local_b; }
				else { local_a = w->joints[s->joint_idx].hinge.local_a; local_b = w->joints[s->joint_idx].hinge.local_b; }
				v3 r_a = rotate(a->rotation, local_a);
				v3 r_b = rotate(b->rotation, local_b);
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
			} else if (s->type == JOINT_DISTANCE) {
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
}

// Store joint accumulated impulses back to persistent storage.
static void joints_post_solve(WorldInternal* w, SolverJoint* joints, int count)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		if (s->type == JOINT_BALL_SOCKET) {
			w->joints[s->joint_idx].warm_lambda3 = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
		} else if (s->type == JOINT_DISTANCE) {
			w->joints[s->joint_idx].warm_lambda1 = s->lambda[0];
		} else if (s->type == JOINT_HINGE) {
			w->joints[s->joint_idx].warm_hinge.lin = V3(s->lambda[0], s->lambda[1], s->lambda[2]);
			w->joints[s->joint_idx].warm_hinge.ang[0] = s->lambda[3];
			w->joints[s->joint_idx].warm_hinge.ang[1] = s->lambda[4];
		}
	}
	afree(joints);
}
