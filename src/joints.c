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
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m };

	v3 ia = a->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		v3 e0 = inv_inertia_mul(a->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(a->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(a->rotation, ia, V3(-r_a.y, r_a.x, 0));
		K[0] += -r_a.z*e0.y + r_a.y*e0.z;
		K[1] += -r_a.z*e1.y + r_a.y*e1.z;
		K[2] += -r_a.z*e2.y + r_a.y*e2.z;
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z;
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z;
		K[5] += -r_a.y*e2.x + r_a.x*e2.y;
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

// -----------------------------------------------------------------------------
// Generic Jacobian helpers: velocity readback, impulse application, effective mass.

// Compute constraint velocity for one DOF row.
static float jac_velocity_f(JacobianRow* row, BodyHot* a, BodyHot* b)
{
	return row->J_a[0]*a->velocity.x + row->J_a[1]*a->velocity.y + row->J_a[2]*a->velocity.z + row->J_a[3]*a->angular_velocity.x + row->J_a[4]*a->angular_velocity.y + row->J_a[5]*a->angular_velocity.z + row->J_b[0]*b->velocity.x + row->J_b[1]*b->velocity.y + row->J_b[2]*b->velocity.z + row->J_b[3]*b->angular_velocity.x + row->J_b[4]*b->angular_velocity.y + row->J_b[5]*b->angular_velocity.z;
}

// Apply impulse for one DOF row: v += M^-1 * J^T * lambda.
static void jac_apply(JacobianRow* row, float lambda, BodyHot* a, BodyHot* b)
{
	float ima = a->inv_mass, imb = b->inv_mass;
	a->velocity.x += ima * row->J_a[0] * lambda;
	a->velocity.y += ima * row->J_a[1] * lambda;
	a->velocity.z += ima * row->J_a[2] * lambda;
	v3 ang_a = inv_inertia_mul(a->rotation, a->inv_inertia_local, V3(row->J_a[3]*lambda, row->J_a[4]*lambda, row->J_a[5]*lambda));
	a->angular_velocity = add(a->angular_velocity, ang_a);
	b->velocity.x += imb * row->J_b[0] * lambda;
	b->velocity.y += imb * row->J_b[1] * lambda;
	b->velocity.z += imb * row->J_b[2] * lambda;
	v3 ang_b = inv_inertia_mul(b->rotation, b->inv_inertia_local, V3(row->J_b[3]*lambda, row->J_b[4]*lambda, row->J_b[5]*lambda));
	b->angular_velocity = add(b->angular_velocity, ang_b);
}

// Compute scalar effective mass for one DOF row: 1 / (J * M^-1 * J^T + softness).
static float jac_eff_mass(JacobianRow* row, BodyHot* a, BodyHot* b, float softness)
{
	float ima = a->inv_mass, imb = b->inv_mass;
	float k = ima * (row->J_a[0]*row->J_a[0] + row->J_a[1]*row->J_a[1] + row->J_a[2]*row->J_a[2]) + imb * (row->J_b[0]*row->J_b[0] + row->J_b[1]*row->J_b[1] + row->J_b[2]*row->J_b[2]);
	v3 wa = inv_inertia_mul(a->rotation, a->inv_inertia_local, V3(row->J_a[3], row->J_a[4], row->J_a[5]));
	k += row->J_a[3]*wa.x + row->J_a[4]*wa.y + row->J_a[5]*wa.z;
	v3 wb = inv_inertia_mul(b->rotation, b->inv_inertia_local, V3(row->J_b[3], row->J_b[4], row->J_b[5]));
	k += row->J_b[3]*wb.x + row->J_b[4]*wb.y + row->J_b[5]*wb.z;
	k += softness;
	return k > 1e-12f ? 1.0f / k : 0.0f;
}

// -----------------------------------------------------------------------------
// Centralized Jacobian fill: the ONE function that knows about joint types.
// Fills s->r_a, s->r_b, s->rows[], s->bias[], s->rows[].eff_mass from current body state.

static void joint_fill_rows(SolverJoint* s, BodyHot* a, BodyHot* b, WorldInternal* w, float dt)
{
	float ptv, soft;
	JointInternal* j = &w->joints[s->joint_idx];

	for (int d = 0; d < s->dof; d++) { memset(s->rows[d].J_a, 0, 6 * sizeof(float)); memset(s->rows[d].J_b, 0, 6 * sizeof(float)); }

	if (s->type == JOINT_BALL_SOCKET) {
		s->r_a = rotate(a->rotation, j->ball_socket.local_a);
		s->r_b = rotate(b->rotation, j->ball_socket.local_b);
		spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 ra = s->r_a, rb = s->r_b;
		// J_a = [-I, -skew(r_a)], J_b = [I, skew(r_b)]
		s->rows[0].J_a[0] = -1; s->rows[0].J_a[4] = -ra.z; s->rows[0].J_a[5] =  ra.y;
		s->rows[0].J_b[0] =  1; s->rows[0].J_b[4] =  rb.z; s->rows[0].J_b[5] = -rb.y;
		s->rows[1].J_a[1] = -1; s->rows[1].J_a[3] =  ra.z; s->rows[1].J_a[5] = -ra.x;
		s->rows[1].J_b[1] =  1; s->rows[1].J_b[3] = -rb.z; s->rows[1].J_b[5] =  rb.x;
		s->rows[2].J_a[2] = -1; s->rows[2].J_a[3] = -ra.y; s->rows[2].J_a[4] =  ra.x;
		s->rows[2].J_b[2] =  1; s->rows[2].J_b[3] =  rb.y; s->rows[2].J_b[4] = -rb.x;

		v3 anchor_a = add(a->position, s->r_a);
		v3 anchor_b = add(b->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;
	} else if (s->type == JOINT_DISTANCE) {
		s->r_a = rotate(a->rotation, j->distance.local_a);
		s->r_b = rotate(b->rotation, j->distance.local_b);
		spring_compute(j->distance.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 anchor_a = add(a->position, s->r_a);
		v3 anchor_b = add(b->position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		float dist_val = len(delta);
		v3 axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

		v3 rxa = cross(s->r_a, axis), rxb = cross(s->r_b, axis);
		s->rows[0].J_a[0] = -axis.x; s->rows[0].J_a[1] = -axis.y; s->rows[0].J_a[2] = -axis.z;
		s->rows[0].J_a[3] = -rxa.x;  s->rows[0].J_a[4] = -rxa.y;  s->rows[0].J_a[5] = -rxa.z;
		s->rows[0].J_b[0] =  axis.x; s->rows[0].J_b[1] =  axis.y; s->rows[0].J_b[2] =  axis.z;
		s->rows[0].J_b[3] =  rxb.x;  s->rows[0].J_b[4] =  rxb.y;  s->rows[0].J_b[5] =  rxb.z;

		s->pos_error[0] = dist_val - j->distance.rest_length;
	} else if (s->type == JOINT_HINGE) {
		s->r_a = rotate(a->rotation, j->hinge.local_a);
		s->r_b = rotate(b->rotation, j->hinge.local_b);
		spring_compute(j->hinge.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 ra = s->r_a, rb = s->r_b;
		// Linear rows 0-2: same as ball socket (J_a = [-I, -skew(r_a)], J_b = [I, skew(r_b)])
		s->rows[0].J_a[0] = -1; s->rows[0].J_a[4] = -ra.z; s->rows[0].J_a[5] =  ra.y;
		s->rows[0].J_b[0] =  1; s->rows[0].J_b[4] =  rb.z; s->rows[0].J_b[5] = -rb.y;
		s->rows[1].J_a[1] = -1; s->rows[1].J_a[3] =  ra.z; s->rows[1].J_a[5] = -ra.x;
		s->rows[1].J_b[1] =  1; s->rows[1].J_b[3] = -rb.z; s->rows[1].J_b[5] =  rb.x;
		s->rows[2].J_a[2] = -1; s->rows[2].J_a[3] = -ra.y; s->rows[2].J_a[4] =  ra.x;
		s->rows[2].J_b[2] =  1; s->rows[2].J_b[3] =  rb.y; s->rows[2].J_b[4] = -rb.x;

		// Angular rows 3-4
		v3 axis_a = norm(rotate(a->rotation, j->hinge.local_axis_a));
		v3 axis_b = norm(rotate(b->rotation, j->hinge.local_axis_b));
		v3 t1, t2;
		hinge_tangent_basis(axis_a, &t1, &t2);
		v3 u1 = cross(t1, axis_b);
		v3 u2 = cross(t2, axis_b);

		s->rows[3].J_a[3] =  u1.x; s->rows[3].J_a[4] =  u1.y; s->rows[3].J_a[5] =  u1.z;
		s->rows[3].J_b[3] = -u1.x; s->rows[3].J_b[4] = -u1.y; s->rows[3].J_b[5] = -u1.z;
		s->rows[4].J_a[3] =  u2.x; s->rows[4].J_a[4] =  u2.y; s->rows[4].J_a[5] =  u2.z;
		s->rows[4].J_b[3] = -u2.x; s->rows[4].J_b[4] = -u2.y; s->rows[4].J_b[5] = -u2.z;

		v3 anchor_a = add(a->position, s->r_a);
		v3 anchor_b = add(b->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;
		s->pos_error[3] = dot(t1, axis_b);
		s->pos_error[4] = dot(t2, axis_b);
	}

	// Generic: compute eff_mass and bias from pos_error for all DOFs
	for (int d = 0; d < s->dof; d++) {
		s->rows[d].eff_mass = jac_eff_mass(&s->rows[d], a, b, soft);
		s->bias[d] = ptv * s->pos_error[d];
	}
}

// -----------------------------------------------------------------------------
// Pre-solve: build unified SolverJoint array from persistent JointInternal data.

static void joints_pre_solve(WorldInternal* w, float dt, SolverJoint** out_joints)
{
	CK_DYNA SolverJoint* joints = NULL;
	int joint_count = asize(w->joints);

	for (int i = 0; i < joint_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		SolverJoint s = {0};
		s.body_a = j->body_a;
		s.body_b = j->body_b;
		s.joint_idx = i;
		s.type = j->type;
		s.dof = j->type == JOINT_BALL_SOCKET ? 3 : j->type == JOINT_DISTANCE ? 1 : 5;
		solver_joint_init_bounds(&s);

		joint_fill_rows(&s, a, b, w, dt);

		// Warm start from persistent storage
		for (int d = 0; d < s.dof; d++) s.lambda[d] = j->warm_lambda[d];

		apush(joints, s);
	}

	*out_joints = joints;
}

// Generic warm start: iterate DOFs, apply each row's warm lambda.
static void joints_warm_start(WorldInternal* w, SolverJoint* joints, int count)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		for (int d = 0; d < s->dof; d++) {
			if (s->lambda[d] == 0) continue;
			jac_apply(&s->rows[d], s->lambda[d], a, b);
		}
	}
}

// Generic PGS joint solve: no type switch.
static void solve_joint(WorldInternal* w, SolverJoint* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	for (int d = 0; d < s->dof; d++) {
		float vel_err = jac_velocity_f(&s->rows[d], a, b);
		float rhs = -vel_err - s->bias[d] - s->softness * s->lambda[d];
		float delta = s->rows[d].eff_mass * rhs;
		s->lambda[d] += delta;
		jac_apply(&s->rows[d], delta, a, b);
	}
}

// Split-impulse position correction for joints.
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

	// Refresh Jacobians and position errors from current rotations, override bias with SI gain
	float ptv = 0.1f / dt; // SI gain: 0.1 is stable for shattering; 0.15+ destabilizes hub_12
	for (int i = 0; i < joint_count; i++) {
		SolverJoint* s = &joints[i];
		if (s->softness != 0.0f) continue;
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		joint_fill_rows(s, a, b, w, dt);
		// Override bias: joint_fill_rows uses spring ptv (0 for rigid), we want SI gain
		for (int d = 0; d < s->dof; d++) s->bias[d] = ptv * s->pos_error[d];
	}

	// PGS iterations on joints (pseudo-velocity solve)
	for (int iter = 0; iter < iters; iter++) {
		for (int i = 0; i < joint_count; i++) {
			if (joints[i].softness != 0.0f) continue;
			solve_joint(w, &joints[i]);
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
// TODO: could be made generic with pos_error[] field, but position error is inherently geometric.
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
		for (int d = 0; d < s->dof; d++) w->joints[s->joint_idx].warm_lambda[d] = s->lambda[d];
	}
	afree(joints);
}
