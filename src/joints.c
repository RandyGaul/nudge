// joints.c -- joint constraint solvers (ball socket, distance, hinge, etc.)

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
static void ball_socket_eff_mass(BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, v3 r_a, v3 r_b, float softness, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m };

	v3 ia = sa->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		v3 e0 = inv_inertia_mul(sa->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(sa->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(sa->rotation, ia, V3(-r_a.y, r_a.x, 0));
		K[0] += -r_a.z*e0.y + r_a.y*e0.z;
		K[1] += -r_a.z*e1.y + r_a.y*e1.z;
		K[2] += -r_a.z*e2.y + r_a.y*e2.z;
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z;
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z;
		K[5] += -r_a.y*e2.x + r_a.x*e2.y;
	}

	v3 ib = sb->inv_inertia_local;
	if (ib.x > 0 || ib.y > 0 || ib.z > 0) {
		v3 e0 = inv_inertia_mul(sb->rotation, ib, V3(0, -r_b.z, r_b.y));
		v3 e1 = inv_inertia_mul(sb->rotation, ib, V3(r_b.z, 0, -r_b.x));
		v3 e2 = inv_inertia_mul(sb->rotation, ib, V3(-r_b.y, r_b.x, 0));
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
static void jac_apply(JacobianRow* row, float lambda, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	float ima = a->inv_mass, imb = b->inv_mass;
	a->velocity.x += ima * row->J_a[0] * lambda;
	a->velocity.y += ima * row->J_a[1] * lambda;
	a->velocity.z += ima * row->J_a[2] * lambda;
	v3 ang_a = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, V3(row->J_a[3]*lambda, row->J_a[4]*lambda, row->J_a[5]*lambda));
	a->angular_velocity = add(a->angular_velocity, ang_a);
	b->velocity.x += imb * row->J_b[0] * lambda;
	b->velocity.y += imb * row->J_b[1] * lambda;
	b->velocity.z += imb * row->J_b[2] * lambda;
	v3 ang_b = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, V3(row->J_b[3]*lambda, row->J_b[4]*lambda, row->J_b[5]*lambda));
	b->angular_velocity = add(b->angular_velocity, ang_b);
}

// Compute scalar effective mass for one DOF row: 1 / (J * M^-1 * J^T + softness).
static float jac_eff_mass(JacobianRow* row, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, float softness)
{
	float ima = a->inv_mass, imb = b->inv_mass;
	float k = ima * (row->J_a[0]*row->J_a[0] + row->J_a[1]*row->J_a[1] + row->J_a[2]*row->J_a[2]) + imb * (row->J_b[0]*row->J_b[0] + row->J_b[1]*row->J_b[1] + row->J_b[2]*row->J_b[2]);
	v3 wa = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, V3(row->J_a[3], row->J_a[4], row->J_a[5]));
	k += row->J_a[3]*wa.x + row->J_a[4]*wa.y + row->J_a[5]*wa.z;
	v3 wb = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, V3(row->J_b[3], row->J_b[4], row->J_b[5]));
	k += row->J_b[3]*wb.x + row->J_b[4]*wb.y + row->J_b[5]*wb.z;
	k += softness;
	return k > 1e-12f ? 1.0f / k : 0.0f;
}

// Shared 3-DOF linear block step (ball socket, hinge linear, fixed linear).
// Uses precomputed sym3x3 inverse lin_inv_eff_mass. Jv and impulse apply inline
// from r_a, r_b. lambda_base is the starting DOF index (0 for all current uses).
static inline void solve_point_block(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, int lambda_base)
{
	v3 ra = s->r_a, rb = s->r_b;
	v3 jv = sub(add(b->velocity, cross(b->angular_velocity, rb)), add(a->velocity, cross(a->angular_velocity, ra)));
	float* lam = &s->lambda[lambda_base];
	float* bias = &s->bias[lambda_base];
	v3 rhs = V3(-jv.x - bias[0] - s->softness * lam[0], -jv.y - bias[1] - s->softness * lam[1], -jv.z - bias[2] - s->softness * lam[2]);
	v3 delta = sym3x3_mul_v3(s->lin_inv_eff_mass, rhs);
	lam[0] += delta.x; lam[1] += delta.y; lam[2] += delta.z;
	a->velocity = sub(a->velocity, scale(delta, a->inv_mass));
	b->velocity = add(b->velocity, scale(delta, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, cross(ra, delta)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, cross(rb, delta)));
}

// 1-DOF bounded angular solve (hinge limit/motor, angular motor, cone, twist).
// Jv = (w_a - w_b) . axis. Impulse direction = axis (world). Pure angular: no
// velocity update on linear components.
static inline void solve_bounded_angular(BoundedAxis* br, float* lambda, float bias, float lo, float hi, float softness, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	if (br->eff_mass == 0.0f) return;
	v3 axis = br->axis;
	float jv = dot(axis, sub(a->angular_velocity, b->angular_velocity));
	float r = -jv - bias - softness * (*lambda);
	float dl = br->eff_mass * r;
	float old = *lambda;
	*lambda += dl;
	if (*lambda < lo) *lambda = lo;
	if (*lambda > hi) *lambda = hi;
	float delta = *lambda - old;
	if (delta == 0.0f) return;
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, scale(axis, delta)));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, scale(axis, delta)));
}

// 1-DOF bounded linear solve (distance, prismatic motor). Jv combines linear
// axis component + angular cross terms (precomputed r_cross_a, r_cross_b).
static inline void solve_bounded_linear(BoundedAxis* br, float* lambda, float bias, float lo, float hi, float softness, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	if (br->eff_mass == 0.0f) return;
	v3 axis = br->axis, rxa = br->r_cross_a, rxb = br->r_cross_b;
	float jv = dot(axis, sub(b->velocity, a->velocity)) + dot(rxb, b->angular_velocity) - dot(rxa, a->angular_velocity);
	float r = -jv - bias - softness * (*lambda);
	float dl = br->eff_mass * r;
	float old = *lambda;
	*lambda += dl;
	if (*lambda < lo) *lambda = lo;
	if (*lambda > hi) *lambda = hi;
	float delta = *lambda - old;
	if (delta == 0.0f) return;
	a->velocity = sub(a->velocity, scale(axis, delta * a->inv_mass));
	b->velocity = add(b->velocity, scale(axis, delta * b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, scale(rxa, delta)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, scale(rxb, delta)));
}

static void solve_ball_socket(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	solve_point_block(s, a, b, sa, sb, 0);
}

// Hinge 5- or 6-DOF PGS step. Linear (DOF 0-2) uses sym3x3 inverse; angular
// alignment (DOF 3-4) uses sym2x2 inverse with stored u1/u2 axes. Bounded
// limit/motor (DOF 5, optional) reads rows[5] / eff_mass as before.
static void solve_hinge(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	solve_point_block(s, a, b, sa, sb, 0);

	// Angular 2-DOF: Jv_i = dot(u_i, w_a - w_b). Impulse along u1*d0 + u2*d1 (world frame).
	v3 u1 = s->hinge_u1, u2 = s->hinge_u2;
	v3 dw = sub(a->angular_velocity, b->angular_velocity);
	float jv0 = dot(u1, dw), jv1 = dot(u2, dw);
	float rhs0 = -jv0 - s->bias[3] - s->softness * s->lambda[3];
	float rhs1 = -jv1 - s->bias[4] - s->softness * s->lambda[4];
	float m00 = s->hinge_ang_inv_eff_mass[0], m01 = s->hinge_ang_inv_eff_mass[1], m11 = s->hinge_ang_inv_eff_mass[2];
	float d0 = m00 * rhs0 + m01 * rhs1;
	float d1 = m01 * rhs0 + m11 * rhs1;
	s->lambda[3] += d0;
	s->lambda[4] += d1;
	v3 ang_impulse = add(scale(u1, d0), scale(u2, d1));
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, ang_impulse));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, ang_impulse));

	// Optional DOF 5 (limit/motor): bounded, uses bounded[0] (axis + eff_mass).
	if (s->dof > 5)
		solve_bounded_angular(&s->bounded[0], &s->lambda[5], s->bias[5], s->lo[5], s->hi[5], s->softness, a, b, sa, sb);
}

// Shared 3-DOF angular block step (fixed rotation lock, prismatic rotation lock).
// Jv = w_a - w_b (identity Jacobian). Impulse = delta (world frame).
// Pre-computed ang3_inv_eff_mass = (I_a_world^-1 + I_b_world^-1 + softness*I)^-1.
static inline void solve_ang3_block(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, int lambda_base)
{
	v3 dw = sub(a->angular_velocity, b->angular_velocity);
	float* lam = &s->lambda[lambda_base];
	float* bias = &s->bias[lambda_base];
	v3 rhs = V3(-dw.x - bias[0] - s->softness * lam[0], -dw.y - bias[1] - s->softness * lam[1], -dw.z - bias[2] - s->softness * lam[2]);
	v3 delta = sym3x3_mul_v3(s->ang3_inv_eff_mass, rhs);
	lam[0] += delta.x; lam[1] += delta.y; lam[2] += delta.z;
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, delta));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, delta));
}

// Fixed joint 6-DOF PGS step: linear 3-DOF (point block) + angular 3-DOF (ang3 block).
static void solve_fixed(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	solve_point_block(s, a, b, sa, sb, 0);
	solve_ang3_block(s, a, b, sa, sb, 3);
}

// Prismatic 5- or 6-DOF PGS step: lateral 2-DOF linear + angular 3-DOF + optional motor.
static void solve_prismatic(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	// Lateral 2-DOF linear: Jv_i = t_i . (v_b - v_a) + (r_b x t_i) . w_b - (r_a x t_i) . w_a.
	v3 t1 = s->prism_t1, t2 = s->prism_t2;
	v3 ra = s->r_a, rb = s->r_b;
	v3 dv = sub(b->velocity, a->velocity);
	v3 rat1 = cross(ra, t1), rat2 = cross(ra, t2);
	v3 rbt1 = cross(rb, t1), rbt2 = cross(rb, t2);
	float jv0 = dot(t1, dv) + dot(rbt1, b->angular_velocity) - dot(rat1, a->angular_velocity);
	float jv1 = dot(t2, dv) + dot(rbt2, b->angular_velocity) - dot(rat2, a->angular_velocity);
	float rhs0 = -jv0 - s->bias[0] - s->softness * s->lambda[0];
	float rhs1 = -jv1 - s->bias[1] - s->softness * s->lambda[1];
	float m00 = s->prism_lateral_inv_eff_mass[0], m01 = s->prism_lateral_inv_eff_mass[1], m11 = s->prism_lateral_inv_eff_mass[2];
	float d0 = m00 * rhs0 + m01 * rhs1;
	float d1 = m01 * rhs0 + m11 * rhs1;
	s->lambda[0] += d0; s->lambda[1] += d1;
	v3 lin_impulse = add(scale(t1, d0), scale(t2, d1));
	a->velocity = sub(a->velocity, scale(lin_impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(lin_impulse, b->inv_mass));
	v3 ang_a = add(scale(rat1, d0), scale(rat2, d1));
	v3 ang_b = add(scale(rbt1, d0), scale(rbt2, d1));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, ang_a));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, ang_b));

	// Angular 3-DOF: identical to fixed joint rotation lock.
	solve_ang3_block(s, a, b, sa, sb, 2);

	// Optional DOF 5 (motor along slide axis): bounded, uses bounded[0] (linear).
	if (s->dof > 5)
		solve_bounded_linear(&s->bounded[0], &s->lambda[5], s->bias[5], s->lo[5], s->hi[5], s->softness, a, b, sa, sb);
}

// Helper: build world-frame angular K (= I_a_world^-1 + I_b_world^-1 + soft*I)
// and invert it (sym3x3). BodyHot has precomputed world-space inverse inertia
// at iw_diag/iw_off; we just add the two bodies' contributions plus softness.
static void ang3_build_inv_eff_mass(BodyHot* a, BodyHot* b, float softness, float* out)
{
	float K[6] = {
		a->iw_diag.x + b->iw_diag.x + softness,  // xx
		a->iw_off.x  + b->iw_off.x,               // xy
		a->iw_off.y  + b->iw_off.y,               // xz
		a->iw_diag.y + b->iw_diag.y + softness,  // yy
		a->iw_off.z  + b->iw_off.z,               // yz
		a->iw_diag.z + b->iw_diag.z + softness,  // zz
	};
	sym3x3_inverse(K, out);
}

// Helper: pre-compute sym2x2 inverse for prismatic lateral (2x2 K block).
// K[0,0] = (im_a+im_b) + (r_a x t1)^T I_a^-1 (r_a x t1) + (r_b x t1)^T I_b^-1 (r_b x t1) + soft
// K[1,1] = ... with t2 ... + soft
// K[0,1] = (r_a x t1)^T I_a^-1 (r_a x t2) + same for B (linear cross zero since t1 perp t2)
static void prism_build_lateral_inv_eff_mass(BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, v3 ra, v3 rb, v3 t1, v3 t2, float softness, float* out)
{
	float im_sum = a->inv_mass + b->inv_mass;
	v3 rat1 = cross(ra, t1), rat2 = cross(ra, t2);
	v3 rbt1 = cross(rb, t1), rbt2 = cross(rb, t2);
	v3 ia_rat1 = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, rat1);
	v3 ia_rat2 = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, rat2);
	v3 ib_rbt1 = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, rbt1);
	v3 ib_rbt2 = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, rbt2);
	float k00 = im_sum + dot(rat1, ia_rat1) + dot(rbt1, ib_rbt1) + softness;
	float k11 = im_sum + dot(rat2, ia_rat2) + dot(rbt2, ib_rbt2) + softness;
	float k01 = dot(rat1, ia_rat2) + dot(rbt1, ib_rbt2); // linear part is 0 since t1 perp t2
	float det = k00 * k11 - k01 * k01;
	float inv_det = det > 1e-20f ? 1.0f / det : 0.0f;
	out[0] =  k11 * inv_det;
	out[1] = -k01 * inv_det;
	out[2] =  k00 * inv_det;
}

// -----------------------------------------------------------------------------
// Centralized Jacobian fill: the ONE function that knows about joint types.

static void joint_fill_rows(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb, WorldInternal* w, float dt)
{
	float ptv, soft;
	JointInternal* j = &w->joints[s->joint_idx];


	if (s->type == JOINT_BALL_SOCKET) {
		s->r_a = rotate(sa->rotation, j->ball_socket.local_a);
		s->r_b = rotate(sb->rotation, j->ball_socket.local_b);
		spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 ra = s->r_a, rb = s->r_b;
		// J_a = [-I, skew(r_a)], J_b = [I, -skew(r_b)]

		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;

		// Precompute sym3x3 inverse for the inline block solve in solve_ball_socket.
		ball_socket_eff_mass(a, b, sa, sb, ra, rb, soft, s->lin_inv_eff_mass);
	} else if (s->type == JOINT_DISTANCE) {
		s->r_a = rotate(sa->rotation, j->distance.local_a);
		s->r_b = rotate(sb->rotation, j->distance.local_b);
		spring_compute(j->distance.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		float dist_val = len(delta);
		v3 axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

		v3 rxa = cross(s->r_a, axis), rxb = cross(s->r_b, axis);
		s->dist_axis = axis; // LDL reads this instead of reconstructing from positions
		// bounded[0]: 1-DOF linear along axis (works for both rigid and limit cases).
		s->bounded[0].axis = axis;
		s->bounded[0].r_cross_a = rxa;
		s->bounded[0].r_cross_b = rxb;
		float im_dist = a->inv_mass + b->inv_mass;
		v3 ia_rxa = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, rxa);
		v3 ib_rxb = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, rxb);
		float k_dist = im_dist + dot(rxa, ia_rxa) + dot(rxb, ib_rxb) + soft;
		s->bounded[0].eff_mass = k_dist > 1e-12f ? 1.0f / k_dist : 0.0f;

		s->pos_error[0] = dist_val - j->distance.rest_length;

		// Distance limits: override pos_error and set unilateral bounds
		float dmin = j->distance.limit_min, dmax = j->distance.limit_max;
		if (dmin > 0 || dmax > 0) {
			if (dmin > 0 && dist_val < dmin) {
				s->pos_error[0] = dist_val - dmin;
				s->lo[0] = 0.0f; s->hi[0] = FLT_MAX;
			} else if (dmax > 0 && dist_val > dmax) {
				s->pos_error[0] = dist_val - dmax;
				s->lo[0] = -FLT_MAX; s->hi[0] = 0.0f;
			} else {
				// Within bounds: constraint inactive
				s->pos_error[0] = 0;
			}
		}
	} else if (s->type == JOINT_HINGE) {
		s->r_a = rotate(sa->rotation, j->hinge.local_a);
		s->r_b = rotate(sb->rotation, j->hinge.local_b);
		spring_compute(j->hinge.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 ra = s->r_a, rb = s->r_b;
		// Linear rows 0-2: same as ball socket (J_a = [-I, -skew(r_a)], J_b = [I, skew(r_b)])

		// Angular rows 3-4
		v3 axis_a = norm(rotate(sa->rotation, j->hinge.local_axis_a));
		v3 axis_b = norm(rotate(sb->rotation, j->hinge.local_axis_b));
		v3 t1, t2;
		hinge_tangent_basis(axis_a, &t1, &t2);
		v3 u1 = cross(t1, axis_b);
		v3 u2 = cross(t2, axis_b);


		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;
		s->pos_error[3] = dot(t1, axis_b);
		s->pos_error[4] = dot(t2, axis_b);

		// Precompute inline-solve data: sym3x3 linear inverse (same as ball socket)
		// and sym2x2 angular inverse. K_ang[i,j] = u_i . (I_a^-1 + I_b^-1) u_j.
		ball_socket_eff_mass(a, b, sa, sb, ra, rb, soft, s->lin_inv_eff_mass);
		s->hinge_u1 = u1; s->hinge_u2 = u2;
		v3 ia_u1 = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, u1);
		v3 ia_u2 = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, u2);
		v3 ib_u1 = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, u1);
		v3 ib_u2 = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, u2);
		float k00 = dot(u1, add(ia_u1, ib_u1)) + soft;
		float k11 = dot(u2, add(ia_u2, ib_u2)) + soft;
		float k01 = dot(u1, add(ia_u2, ib_u2));
		float det = k00 * k11 - k01 * k01;
		float inv_det = det > 1e-20f ? 1.0f / det : 0.0f;
		s->hinge_ang_inv_eff_mass[0] =  k11 * inv_det;
		s->hinge_ang_inv_eff_mass[1] = -k01 * inv_det;
		s->hinge_ang_inv_eff_mass[2] =  k00 * inv_det;

		// Hinge DOF 5: motor and/or limit along the hinge axis.
		float hmin = j->hinge.limit_min, hmax = j->hinge.limit_max;
		float motor_max = j->hinge.motor_max_impulse;
		int has_limit = (hmin != 0.0f || hmax != 0.0f);
		int has_motor = (motor_max > 0.0f);
		if (has_limit || has_motor) {
			s->dof = 6;
			// Jacobian: relative angular velocity along hinge axis
			s->pos_error[5] = 0;
			// Default bounds: bilateral (motor or limit will narrow)
			s->lo[5] = -FLT_MAX;
			s->hi[5] = FLT_MAX;

			// Motor: bounded force, velocity target set after generic loop
			if (has_motor) {
				s->lo[5] = -motor_max;
				s->hi[5] = motor_max;
			}

			// Limit: clamp bounds to prevent exceeding angle range.
			if (has_limit) {
				v3 ref_a_w = rotate(sa->rotation, j->hinge.local_ref_a);
				v3 ref_b_w = rotate(sb->rotation, j->hinge.local_ref_b);
				float angle = atan2f(dot(cross(ref_a_w, ref_b_w), axis_a), dot(ref_a_w, ref_b_w));
				if (hmin != 0.0f && angle <= hmin) {
					s->pos_error[5] = hmin - angle;
					if (s->hi[5] > 0.0f) s->hi[5] = 0.0f;
				} else if (hmax != 0.0f && angle >= hmax) {
					s->pos_error[5] = hmax - angle;
					if (s->lo[5] < 0.0f) s->lo[5] = 0.0f;
				} else if (!has_motor) {
					// Limit inactive and no motor: zero eff_mass so solve_bounded_angular short-circuits.
					s->bounded[0].eff_mass = 0.0f;
					goto hinge_dof5_done;
				}
			}
			// Populate bounded[0] for solve_hinge (1-DOF angular along axis_a).
			s->bounded[0].axis = axis_a;
			s->bounded[0].r_cross_a = V3(0,0,0);
			s->bounded[0].r_cross_b = V3(0,0,0);
			v3 ia_axis = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, axis_a);
			v3 ib_axis = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, axis_a);
			float k = dot(axis_a, add(ia_axis, ib_axis)) + soft;
			s->bounded[0].eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;
			hinge_dof5_done:;
		}
	} else if (s->type == JOINT_FIXED) {
		s->r_a = rotate(sa->rotation, j->fixed.local_a);
		s->r_b = rotate(sb->rotation, j->fixed.local_b);
		spring_compute(j->fixed.spring, dt, &ptv, &soft);
		s->softness = soft;

		v3 ra = s->r_a, rb = s->r_b;
		// Linear rows 0-2: identical to ball socket

		// Angular rows 3-5: lock all relative rotation

		// Linear position error
		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;

		// Angular error: world-frame quaternion.
		// diff = R_b * R_b0^-1 * R_a0 * R_a^-1 = delta_b * delta_a^-1 (world frame).
		// Matches identity Jacobian [I, -I] which operates on world-frame angular velocity.
		quat q_err = mul(mul(sb->rotation, inv(j->fixed.local_rel_quat)), inv(sa->rotation));
		float sign_w = q_err.w >= 0.0f ? 1.0f : -1.0f;
		s->pos_error[3] = 2.0f * q_err.x * sign_w;
		s->pos_error[4] = 2.0f * q_err.y * sign_w;
		s->pos_error[5] = 2.0f * q_err.z * sign_w;

		// Precompute inline-solve data for solve_fixed.
		ball_socket_eff_mass(a, b, sa, sb, ra, rb, soft, s->lin_inv_eff_mass);
		ang3_build_inv_eff_mass(a, b, soft, s->ang3_inv_eff_mass);
	} else if (s->type == JOINT_PRISMATIC) {
		s->r_a = rotate(sa->rotation, j->prismatic.local_a);
		s->r_b = rotate(sb->rotation, j->prismatic.local_b);
		spring_compute(j->prismatic.spring, dt, &ptv, &soft);
		s->softness = soft;

		// Slide axis in world space, tangent plane perpendicular to it
		v3 slide_axis = norm(rotate(sa->rotation, j->prismatic.local_axis_a));
		v3 t1, t2;
		hinge_tangent_basis(slide_axis, &t1, &t2);

		v3 ra = s->r_a, rb = s->r_b;
		v3 rxa_t1 = cross(ra, t1), rxb_t1 = cross(rb, t1);
		v3 rxa_t2 = cross(ra, t2), rxb_t2 = cross(rb, t2);

		// Linear rows 0-1: lateral constraints (perpendicular to slide axis)


		// Angular rows 2-4: lock all relative rotation

		// Lateral position error
		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		s->pos_error[0] = dot(delta, t1);
		s->pos_error[1] = dot(delta, t2);

		// Angular error: world-frame quaternion (matches identity Jacobian rows 2-4).
		quat q_err = mul(mul(sb->rotation, inv(j->prismatic.local_rel_quat)), inv(sa->rotation));
		float sign_w = q_err.w >= 0.0f ? 1.0f : -1.0f;
		s->pos_error[2] = 2.0f * q_err.x * sign_w;
		s->pos_error[3] = 2.0f * q_err.y * sign_w;
		s->pos_error[4] = 2.0f * q_err.z * sign_w;

		// Precompute inline-solve data for solve_prismatic.
		s->prism_t1 = t1; s->prism_t2 = t2;
		prism_build_lateral_inv_eff_mass(a, b, sa, sb, ra, rb, t1, t2, soft, s->prism_lateral_inv_eff_mass);
		ang3_build_inv_eff_mass(a, b, soft, s->ang3_inv_eff_mass);

		// Prismatic motor: DOF 5 along slide axis
		if (j->prismatic.motor_max_impulse > 0.0f) {
			s->dof = 6;
			v3 rxa_s = cross(ra, slide_axis), rxb_s = cross(rb, slide_axis);
			s->pos_error[5] = 0;
			s->lo[5] = -j->prismatic.motor_max_impulse;
			s->hi[5] = j->prismatic.motor_max_impulse;
			// Populate bounded[0] for solve_prismatic: 1-DOF linear along slide axis.
			s->bounded[0].axis = slide_axis;
			s->bounded[0].r_cross_a = rxa_s;
			s->bounded[0].r_cross_b = rxb_s;
			float im_sum = a->inv_mass + b->inv_mass;
			v3 ia_rxa = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, rxa_s);
			v3 ib_rxb = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, rxb_s);
			float k = im_sum + dot(rxa_s, ia_rxa) + dot(rxb_s, ib_rxb) + soft;
			s->bounded[0].eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;
		}
	} else if (s->type == JOINT_ANGULAR_MOTOR) {
		s->softness = 0.0f; // motor is velocity-only; no spring
		v3 axis_a_w = norm(rotate(sa->rotation, j->angular_motor.local_axis_a));
		v3 axis_b_w = norm(rotate(sb->rotation, j->angular_motor.local_axis_b));
		v3 axis = norm(add(axis_a_w, axis_b_w)); // symmetric average for stability
		s->pos_error[0] = 0;
		s->lo[0] = -j->angular_motor.max_impulse;
		s->hi[0] = j->angular_motor.max_impulse;
		// Bias convention matches hinge motor: bias = target_speed drives d(angle)/dt = +speed.
		ptv = 0.0f;
		s->bias[0] = j->angular_motor.target_speed;
		// bounded[0]: pure-angular 1-DOF along axis.
		s->bounded[0].axis = axis;
		s->bounded[0].r_cross_a = V3(0,0,0); s->bounded[0].r_cross_b = V3(0,0,0);
		v3 ia_ax = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, axis);
		v3 ib_ax = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, axis);
		float k_am = dot(axis, add(ia_ax, ib_ax));
		s->bounded[0].eff_mass = k_am > 1e-12f ? 1.0f / k_am : 0.0f;
	} else if (s->type == JOINT_CONE_LIMIT) {
		spring_compute(j->cone_limit.spring, dt, &ptv, &soft);
		s->softness = soft;
		v3 axis_a_w = norm(rotate(sa->rotation, j->cone_limit.local_axis_a));
		v3 axis_b_w = norm(rotate(sb->rotation, j->cone_limit.local_axis_b));
		float cos_theta = dot(axis_a_w, axis_b_w);
		float cos_limit = cosf(j->cone_limit.half_angle);
		if (cos_theta < cos_limit) {
			// Violated: swing axis = axis_a x axis_b. Jv = (w_a - w_b) . swing =
			// d(cos_theta)/dt. Positive lambda drives cos_theta up (angle down).
			v3 swing = cross(axis_a_w, axis_b_w);
			float slen = len(swing);
			if (slen > 1e-6f) swing = scale(swing, 1.0f / slen);
			else { v3 t2; hinge_tangent_basis(axis_a_w, &swing, &t2); }
			s->pos_error[0] = acosf(cos_theta > 1.0f ? 1.0f : (cos_theta < -1.0f ? -1.0f : cos_theta)) - j->cone_limit.half_angle;
			s->lo[0] = 0.0f; s->hi[0] = FLT_MAX;
			// bounded[0]: pure-angular along swing axis.
			s->bounded[0].axis = swing;
			s->bounded[0].r_cross_a = V3(0,0,0); s->bounded[0].r_cross_b = V3(0,0,0);
			v3 ia_sw = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, swing);
			v3 ib_sw = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, swing);
			float k_cone = dot(swing, add(ia_sw, ib_sw)) + soft;
			s->bounded[0].eff_mass = k_cone > 1e-12f ? 1.0f / k_cone : 0.0f;
		} else {
			// Inactive: zero row, zero bounds = no-op.
			s->pos_error[0] = 0;
			s->lo[0] = 0; s->hi[0] = 0;
			s->bounded[0].eff_mass = 0.0f;
		}
	} else if (s->type == JOINT_TWIST_LIMIT) {
		spring_compute(j->twist_limit.spring, dt, &ptv, &soft);
		s->softness = soft;
		// Twist axis in world (symmetric average for stability when axes nearly aligned).
		v3 axis_a_w = norm(rotate(sa->rotation, j->twist_limit.local_axis_a));
		v3 axis_b_w = norm(rotate(sb->rotation, j->twist_limit.local_axis_b));
		v3 twist_axis = norm(add(axis_a_w, axis_b_w));
		// Twist angle via quaternion swing-twist decomposition: project q_rel onto twist axis.
		quat q_rel = mul(inv(sa->rotation), sb->rotation);
		v3 qv = V3(q_rel.x, q_rel.y, q_rel.z);
		v3 local_axis = j->twist_limit.local_axis_a;
		float proj = dot(qv, local_axis);
		quat q_twist = { proj * local_axis.x, proj * local_axis.y, proj * local_axis.z, q_rel.w };
		float tl = sqrtf(q_twist.x*q_twist.x + q_twist.y*q_twist.y + q_twist.z*q_twist.z + q_twist.w*q_twist.w);
		if (tl > 1e-12f) { q_twist.x /= tl; q_twist.y /= tl; q_twist.z /= tl; q_twist.w /= tl; }
		else q_twist = (quat){ 0, 0, 0, 1 };
		float sign = q_twist.w >= 0 ? 1.0f : -1.0f;
		float twist_angle = 2.0f * atan2f(sign * proj, sign * q_twist.w);
		float tmin = j->twist_limit.limit_min, tmax = j->twist_limit.limit_max;
		// Jv = (w_a - w_b) . axis = -d(twist)/dt (same convention as hinge limit).
		// twist > tmax: want d(twist)/dt<0, Jv>0, positive lambda -> lo=0.
		// twist < tmin: want d(twist)/dt>0, Jv<0, negative lambda -> hi=0.
		int twist_active = 0;
		if (twist_angle > tmax) {
			s->pos_error[0] = twist_angle - tmax;
			s->lo[0] = 0.0f; s->hi[0] = FLT_MAX;
			twist_active = 1;
		} else if (twist_angle < tmin) {
			s->pos_error[0] = twist_angle - tmin;
			s->lo[0] = -FLT_MAX; s->hi[0] = 0.0f;
			twist_active = 1;
		} else {
			s->pos_error[0] = 0;
			s->lo[0] = 0; s->hi[0] = 0;
		}
		if (twist_active) {
			s->bounded[0].axis = twist_axis;
			s->bounded[0].r_cross_a = V3(0,0,0); s->bounded[0].r_cross_b = V3(0,0,0);
			v3 ia_tw = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, twist_axis);
			v3 ib_tw = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, twist_axis);
			float k_tw = dot(twist_axis, add(ia_tw, ib_tw)) + soft;
			s->bounded[0].eff_mass = k_tw > 1e-12f ? 1.0f / k_tw : 0.0f;
		} else {
			s->bounded[0].eff_mass = 0.0f;
		}
	} else if (s->type == JOINT_SWING_TWIST) {
		// Linear 3-DOF ball socket at DOFs 0-2, cone at DOF 3, twist at DOF 4.
		s->r_a = rotate(sa->rotation, j->swing_twist.local_a);
		s->r_b = rotate(sb->rotation, j->swing_twist.local_b);
		spring_compute(j->swing_twist.spring, dt, &ptv, &soft);
		s->softness = soft;
		s->dof = 5;
		v3 ra = s->r_a, rb = s->r_b;
		v3 anchor_a = add(sa->position, s->r_a);
		v3 anchor_b = add(sb->position, s->r_b);
		v3 err = sub(anchor_b, anchor_a);
		s->pos_error[0] = err.x; s->pos_error[1] = err.y; s->pos_error[2] = err.z;
		ball_socket_eff_mass(a, b, sa, sb, ra, rb, soft, s->lin_inv_eff_mass);

		// Cone DOF 3.
		v3 axis_a_w = norm(rotate(sa->rotation, j->swing_twist.local_axis_a));
		v3 axis_b_w = norm(rotate(sb->rotation, j->swing_twist.local_axis_b));
		float cos_theta = dot(axis_a_w, axis_b_w);
		float cos_limit = cosf(j->swing_twist.cone_half_angle);
		if (cos_theta < cos_limit) {
			v3 swing = cross(axis_a_w, axis_b_w);
			float slen = len(swing);
			if (slen > 1e-6f) swing = scale(swing, 1.0f / slen);
			else { v3 t2; hinge_tangent_basis(axis_a_w, &swing, &t2); }
			s->pos_error[3] = acosf(cos_theta > 1.0f ? 1.0f : (cos_theta < -1.0f ? -1.0f : cos_theta)) - j->swing_twist.cone_half_angle;
			s->lo[3] = 0.0f; s->hi[3] = FLT_MAX;
			s->bounded[0].axis = swing;
			s->bounded[0].r_cross_a = V3(0,0,0); s->bounded[0].r_cross_b = V3(0,0,0);
			v3 ia_sw = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, swing);
			v3 ib_sw = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, swing);
			float k_cone = dot(swing, add(ia_sw, ib_sw)) + soft;
			s->bounded[0].eff_mass = k_cone > 1e-12f ? 1.0f / k_cone : 0.0f;
		} else {
			s->pos_error[3] = 0; s->lo[3] = 0; s->hi[3] = 0;
			s->bounded[0].eff_mass = 0.0f;
		}

		// Twist DOF 4.
		v3 twist_axis = norm(add(axis_a_w, axis_b_w));
		quat q_rel = mul(inv(sa->rotation), sb->rotation);
		v3 qv = V3(q_rel.x, q_rel.y, q_rel.z);
		v3 local_axis = j->swing_twist.local_axis_a;
		float proj = dot(qv, local_axis);
		quat q_twist = { proj * local_axis.x, proj * local_axis.y, proj * local_axis.z, q_rel.w };
		float tl = sqrtf(q_twist.x*q_twist.x + q_twist.y*q_twist.y + q_twist.z*q_twist.z + q_twist.w*q_twist.w);
		if (tl > 1e-12f) { q_twist.x /= tl; q_twist.y /= tl; q_twist.z /= tl; q_twist.w /= tl; }
		else q_twist = (quat){ 0, 0, 0, 1 };
		float sign = q_twist.w >= 0 ? 1.0f : -1.0f;
		float twist_angle = 2.0f * atan2f(sign * proj, sign * q_twist.w);
		float tmin = j->swing_twist.twist_min, tmax = j->swing_twist.twist_max;
		int twist_active = 0;
		if (twist_angle > tmax) { s->pos_error[4] = twist_angle - tmax; s->lo[4] = 0.0f; s->hi[4] = FLT_MAX; twist_active = 1; }
		else if (twist_angle < tmin) { s->pos_error[4] = twist_angle - tmin; s->lo[4] = -FLT_MAX; s->hi[4] = 0.0f; twist_active = 1; }
		else { s->pos_error[4] = 0; s->lo[4] = 0; s->hi[4] = 0; }
		if (twist_active) {
			s->bounded[1].axis = twist_axis;
			s->bounded[1].r_cross_a = V3(0,0,0); s->bounded[1].r_cross_b = V3(0,0,0);
			v3 ia_tw = inv_inertia_mul(sa->rotation, sa->inv_inertia_local, twist_axis);
			v3 ib_tw = inv_inertia_mul(sb->rotation, sb->inv_inertia_local, twist_axis);
			float k_tw = dot(twist_axis, add(ia_tw, ib_tw)) + soft;
			s->bounded[1].eff_mass = k_tw > 1e-12f ? 1.0f / k_tw : 0.0f;
		} else {
			s->bounded[1].eff_mass = 0.0f;
		}
	}

	// Generic: compute eff_mass and bias from pos_error for all DOFs
	for (int d = 0; d < s->dof; d++) {
		s->bias[d] = ptv * s->pos_error[d];
	}

	// Motor/limit bias for DOF 5.
	// When a motor is present, it dominates the bias and the lo/hi clamps handle
	// limit enforcement. Falling through to Baumgarte on any nonzero pos_error
	// would disable the motor whenever the arm sits even a float-epsilon past a
	// limit, so the motor can never drive the arm back off the limit.
	if (s->type == JOINT_HINGE && s->dof == 6) {
		float hmin = j->hinge.limit_min, hmax = j->hinge.limit_max;
		if (j->hinge.motor_max_impulse > 0.0f) {
			float speed = j->hinge.motor_speed;
			// Clamp motor speed so it can't overshoot limits in one substep.
			// Jv = omega_a.axis - omega_b.axis = -d(angle)/dt. Convergence: d(angle)/dt = bias.
			// So bias = speed drives angle at +speed. Clamp to keep within [hmin, hmax].
			if (hmin != 0.0f || hmax != 0.0f) {
				v3 axis_a_w = norm(rotate(sa->rotation, j->hinge.local_axis_a));
				v3 ref_a_w = rotate(sa->rotation, j->hinge.local_ref_a);
				v3 ref_b_w = rotate(sb->rotation, j->hinge.local_ref_b);
				float angle = atan2f(dot(cross(ref_a_w, ref_b_w), axis_a_w), dot(ref_a_w, ref_b_w));
				if (hmax != 0.0f) { float max_speed = fmaxf((hmax - angle) / dt, 0.0f); if (speed > max_speed) speed = max_speed; }
				if (hmin != 0.0f) { float min_speed = fminf((hmin - angle) / dt, 0.0f); if (speed < min_speed) speed = min_speed; }
			}
			s->bias[5] = speed;
		} else if (s->pos_error[5] != 0.0f) {
			// No motor, limit active: strong Baumgarte correction to push back.
			s->bias[5] = (0.2f / dt) * s->pos_error[5];
		}
	}
	if (s->type == JOINT_PRISMATIC && j->prismatic.motor_max_impulse > 0.0f && s->dof == 6) {
		// Prismatic Jv = (v_b - v_a).axis = d(disp)/dt. Convergence: Jv = -bias.
		// Negate so positive motor_speed drives positive displacement.
		if (s->pos_error[5] == 0.0f) s->bias[5] = -j->prismatic.motor_speed;
	}
	// Angular motor: bias = target_speed (hinge motor convention). Override the
	// generic bias=ptv*pos_error pass above since pos_error is 0 for a pure motor.
	if (s->type == JOINT_ANGULAR_MOTOR) {
		s->bias[0] = j->angular_motor.target_speed;
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
		BodyState* sa = &w->body_state[j->body_a];
		BodyState* sb = &w->body_state[j->body_b];

		SolverJoint s = {0};
		s.body_a = j->body_a;
		s.body_b = j->body_b;
		s.joint_idx = i;
		s.type = j->type;
		int base_dof = j->type == JOINT_BALL_SOCKET ? 3 : j->type == JOINT_DISTANCE ? 1 : j->type == JOINT_FIXED ? 6 : j->type == JOINT_SWING_TWIST ? 5 : (j->type == JOINT_ANGULAR_MOTOR || j->type == JOINT_CONE_LIMIT || j->type == JOINT_TWIST_LIMIT) ? 1 : 5;
		if (j->type == JOINT_HINGE && (j->hinge.limit_min != 0.0f || j->hinge.limit_max != 0.0f || j->hinge.motor_max_impulse > 0.0f)) base_dof = 6;
		if (j->type == JOINT_PRISMATIC && j->prismatic.motor_max_impulse > 0.0f) base_dof = 6;
		s.dof = base_dof;
		solver_joint_init_bounds(&s);

		joint_fill_rows(&s, a, b, sa, sb, w, dt);

		// Warm start from persistent storage
		for (int d = 0; d < s.dof; d++) s.lambda[d] = w->joint_hot[i].warm_lambda[d];

		apush(joints, s);
	}

	*out_joints = joints;
}

// Refresh joint Jacobians, lever arms, limit activation, and biases from current
// body state. Called at the start of each substep (after integrate_positions/velocities)
// so that limit DOFs activate based on current rotations, not stale frame-start state.
static void joints_refresh_substep(WorldInternal* w, SolverJoint* joints, int count, float dt)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		JointInternal* j = &w->joints[s->joint_idx];
		// Only refresh joints that go through PGS (not LDL).
		// LDL joints have their own lever arm refresh in ldl_refresh_lever_arms.
		// Also skip joints without limits -- they don't need per-substep refresh.
		int has_dof5 = (j->type == JOINT_HINGE && (j->hinge.limit_min != 0.0f || j->hinge.limit_max != 0.0f || j->hinge.motor_max_impulse > 0.0f)) || (j->type == JOINT_DISTANCE && (j->distance.limit_min > 0 || j->distance.limit_max > 0)) || (j->type == JOINT_PRISMATIC && j->prismatic.motor_max_impulse > 0.0f);
		if (!has_dof5) continue;
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		BodyState* sa = &w->body_state[s->body_a];
		BodyState* sb = &w->body_state[s->body_b];
		int old_dof = s->dof;
		// Reset DOF to base (joint_fill_rows will set final dof including limits/motors)
		s->dof = j->type == JOINT_DISTANCE ? 1 : 6; // hinge with limits or prismatic with motor = 6
		solver_joint_init_bounds(s);
		float saved_lambda[JOINT_MAX_DOF];
		memcpy(saved_lambda, s->lambda, sizeof(saved_lambda));
		joint_fill_rows(s, a, b, sa, sb, w, dt);
		memcpy(s->lambda, saved_lambda, sizeof(saved_lambda));
		if (s->dof < old_dof) {
			for (int d = s->dof; d < old_dof; d++) s->lambda[d] = 0;
		}
		if (s->dof > old_dof) {
			for (int d = old_dof; d < s->dof; d++) s->lambda[d] = 0;
		}
	}
}

// Apply a 3-DOF linear (point) impulse: lambda[0..2] is a world-space impulse
// vector at lever arms r_a, r_b. Used by ball socket / hinge / fixed / swing
// twist linear blocks during warm start. Mirrors solve_point_block's apply.
static inline void apply_point_impulse(v3 r_a, v3 r_b, v3 impulse, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	if (impulse.x == 0 && impulse.y == 0 && impulse.z == 0) return;
	a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(impulse, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, cross(r_b, impulse)));
}

// Apply a 3-DOF angular identity impulse (fixed-style): impulse vector applied
// directly to angular velocities, no linear component.
static inline void apply_ang3_impulse(v3 impulse, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	if (impulse.x == 0 && impulse.y == 0 && impulse.z == 0) return;
	a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, impulse));
	b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, impulse));
}

// Apply a bounded-1-DOF impulse (lambda * axis). is_angular controls whether
// linear velocities are touched. Mirrors solve_bounded_{angular,linear} apply.
static inline void apply_bounded_impulse(BoundedAxis* br, float lambda, int is_angular, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	if (lambda == 0.0f) return;
	v3 axis = br->axis;
	if (is_angular) {
		a->angular_velocity = add(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, scale(axis, lambda)));
		b->angular_velocity = sub(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, scale(axis, lambda)));
	} else {
		a->velocity = sub(a->velocity, scale(axis, lambda * a->inv_mass));
		b->velocity = add(b->velocity, scale(axis, lambda * b->inv_mass));
		a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, scale(br->r_cross_a, lambda)));
		b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, scale(br->r_cross_b, lambda)));
	}
}

// Type-specialized warm start: apply cached warm lambdas as constraint impulses
// using the same geometric data the per-type solvers use. No rows[] reads.
static void joints_warm_start(WorldInternal* w, SolverJoint* joints, int count)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		BodyState* sa = &w->body_state[s->body_a];
		BodyState* sb = &w->body_state[s->body_b];
		switch (s->type) {
		case JOINT_BALL_SOCKET:
			apply_point_impulse(s->r_a, s->r_b, V3(s->lambda[0], s->lambda[1], s->lambda[2]), a, b, sa, sb);
			break;
		case JOINT_HINGE:
			apply_point_impulse(s->r_a, s->r_b, V3(s->lambda[0], s->lambda[1], s->lambda[2]), a, b, sa, sb);
			apply_ang3_impulse(add(scale(s->hinge_u1, s->lambda[3]), scale(s->hinge_u2, s->lambda[4])), a, b, sa, sb);
			if (s->dof > 5) apply_bounded_impulse(&s->bounded[0], s->lambda[5], 1, a, b, sa, sb);
			break;
		case JOINT_FIXED:
			apply_point_impulse(s->r_a, s->r_b, V3(s->lambda[0], s->lambda[1], s->lambda[2]), a, b, sa, sb);
			apply_ang3_impulse(V3(s->lambda[3], s->lambda[4], s->lambda[5]), a, b, sa, sb);
			break;
		case JOINT_PRISMATIC: {
			// Lateral linear: impulse = t1*lambda[0] + t2*lambda[1] with cross terms.
			v3 t1 = s->prism_t1, t2 = s->prism_t2;
			float l0 = s->lambda[0], l1 = s->lambda[1];
			v3 lin_imp = add(scale(t1, l0), scale(t2, l1));
			if (lin_imp.x != 0 || lin_imp.y != 0 || lin_imp.z != 0) {
				a->velocity = sub(a->velocity, scale(lin_imp, a->inv_mass));
				b->velocity = add(b->velocity, scale(lin_imp, b->inv_mass));
				v3 ra_cross = add(scale(cross(s->r_a, t1), l0), scale(cross(s->r_a, t2), l1));
				v3 rb_cross = add(scale(cross(s->r_b, t1), l0), scale(cross(s->r_b, t2), l1));
				a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, ra_cross));
				b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, rb_cross));
			}
			apply_ang3_impulse(V3(s->lambda[2], s->lambda[3], s->lambda[4]), a, b, sa, sb);
			if (s->dof > 5) apply_bounded_impulse(&s->bounded[0], s->lambda[5], 0, a, b, sa, sb);
			break;
		}
		case JOINT_DISTANCE:
			apply_bounded_impulse(&s->bounded[0], s->lambda[0], 0, a, b, sa, sb);
			break;
		case JOINT_ANGULAR_MOTOR:
		case JOINT_CONE_LIMIT:
		case JOINT_TWIST_LIMIT:
			apply_bounded_impulse(&s->bounded[0], s->lambda[0], 1, a, b, sa, sb);
			break;
		case JOINT_SWING_TWIST:
			apply_point_impulse(s->r_a, s->r_b, V3(s->lambda[0], s->lambda[1], s->lambda[2]), a, b, sa, sb);
			apply_bounded_impulse(&s->bounded[0], s->lambda[3], 1, a, b, sa, sb);
			apply_bounded_impulse(&s->bounded[1], s->lambda[4], 1, a, b, sa, sb);
			break;
		}
	}
}

static int solver_joint_has_limits(SolverJoint* s)
{
	for (int d = 0; d < s->dof; d++) {
		if (s->lo[d] > -1e18f || s->hi[d] < 1e18f) return 1;
	}
	return 0;
}

// Swing-twist 5-DOF PGS step: linear 3-DOF (point block) + cone (DOF 3, bounded
// angular) + twist (DOF 4, bounded angular). Each part either fires (when active)
// or short-circuits via eff_mass==0 (set by joint_fill_rows when within range).
static void solve_swing_twist(SolverJoint* s, BodyHot* a, BodyHot* b, BodyState* sa, BodyState* sb)
{
	solve_point_block(s, a, b, sa, sb, 0);
	solve_bounded_angular(&s->bounded[0], &s->lambda[3], s->bias[3], s->lo[3], s->hi[3], s->softness, a, b, sa, sb);
	solve_bounded_angular(&s->bounded[1], &s->lambda[4], s->bias[4], s->lo[4], s->hi[4], s->softness, a, b, sa, sb);
}

// PGS joint solve: dispatches to a specialized inline solver per type. No
// generic block-LDL path -- every joint type has a specialized solver that
// uses precomputed effective-mass blocks (lin_inv_eff_mass, ang3_inv_eff_mass,
// hinge_ang_inv_eff_mass, prism_lateral_inv_eff_mass) or per-DOF bounded[].
static void solve_joint(WorldInternal* w, SolverJoint* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	BodyState* sa = &w->body_state[s->body_a];
	BodyState* sb = &w->body_state[s->body_b];

	switch (s->type) {
	case JOINT_BALL_SOCKET: solve_ball_socket(s, a, b, sa, sb); break;
	case JOINT_HINGE:       solve_hinge(s, a, b, sa, sb); break;
	case JOINT_FIXED:       solve_fixed(s, a, b, sa, sb); break;
	case JOINT_PRISMATIC:   solve_prismatic(s, a, b, sa, sb); break;
	case JOINT_DISTANCE:    solve_bounded_linear(&s->bounded[0], &s->lambda[0], s->bias[0], s->lo[0], s->hi[0], s->softness, a, b, sa, sb); break;
	case JOINT_ANGULAR_MOTOR:
	case JOINT_CONE_LIMIT:
	case JOINT_TWIST_LIMIT: solve_bounded_angular(&s->bounded[0], &s->lambda[0], s->bias[0], s->lo[0], s->hi[0], s->softness, a, b, sa, sb); break;
	case JOINT_SWING_TWIST: solve_swing_twist(s, a, b, sa, sb); break;
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
		saved_vel[i] = body_vel(w, i);
		saved_ang[i] = body_angvel(w, i);
		body_vel(w, i) = V3(0, 0, 0);
		body_angvel(w, i) = V3(0, 0, 0);
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
		BodyState* sa = &w->body_state[s->body_a];
		BodyState* sb = &w->body_state[s->body_b];
		joint_fill_rows(s, a, b, sa, sb, w, dt);
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
		BodyState* s = &w->body_state[i];
		s->position = add(s->position, scale(h->velocity, dt));
		v3 ww = h->angular_velocity;
		if (ww.x != 0 || ww.y != 0 || ww.z != 0) {
			quat spin = { ww.x, ww.y, ww.z, 0.0f };
			quat dq = mul(spin, s->rotation);
			s->rotation.x += 0.5f * dt * dq.x; s->rotation.y += 0.5f * dt * dq.y;
			s->rotation.z += 0.5f * dt * dq.z; s->rotation.w += 0.5f * dt * dq.w;
			float ql = sqrtf(s->rotation.x*s->rotation.x + s->rotation.y*s->rotation.y + s->rotation.z*s->rotation.z + s->rotation.w*s->rotation.w);
			float inv_ql = 1.0f / (ql > 1e-15f ? ql : 1.0f);
			s->rotation.x *= inv_ql; s->rotation.y *= inv_ql; s->rotation.z *= inv_ql; s->rotation.w *= inv_ql;
		}
	}

	// Restore real velocities and lambda
	for (int i = 0; i < body_count; i++) {
		body_vel(w, i) = saved_vel[i];
		body_angvel(w, i) = saved_ang[i];
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
			BodyState* sa = &w->body_state[s->body_a];
			BodyState* sb = &w->body_state[s->body_b];
			if (s->type == JOINT_BALL_SOCKET || s->type == JOINT_HINGE || s->type == JOINT_FIXED) {
				v3 local_a, local_b;
				if (s->type == JOINT_BALL_SOCKET) { local_a = w->joints[s->joint_idx].ball_socket.local_a; local_b = w->joints[s->joint_idx].ball_socket.local_b; }
				else if (s->type == JOINT_HINGE) { local_a = w->joints[s->joint_idx].hinge.local_a; local_b = w->joints[s->joint_idx].hinge.local_b; }
				else { local_a = w->joints[s->joint_idx].fixed.local_a; local_b = w->joints[s->joint_idx].fixed.local_b; }
				v3 r_a = rotate(sa->rotation, local_a);
				v3 r_b = rotate(sb->rotation, local_b);
				v3 error = sub(add(sb->position, r_b), add(sa->position, r_a));
				float err_len = len(error);
				if (err_len < 1e-6f) continue;
				float inv_mass_sum = a->inv_mass + b->inv_mass;
				if (inv_mass_sum <= 0) continue;
				float correction = SOLVER_POS_BAUMGARTE * err_len;
				if (correction > SOLVER_POS_MAX_CORRECTION) correction = SOLVER_POS_MAX_CORRECTION;
				v3 n = scale(error, 1.0f / err_len);
				float P = correction / inv_mass_sum;
				sa->position = add(sa->position, scale(n, P * a->inv_mass));
				sb->position = sub(sb->position, scale(n, P * b->inv_mass));
			} else if (s->type == JOINT_PRISMATIC) {
				v3 r_a = rotate(sa->rotation, w->joints[s->joint_idx].prismatic.local_a);
				v3 r_b = rotate(sb->rotation, w->joints[s->joint_idx].prismatic.local_b);
				v3 delta = sub(add(sb->position, r_b), add(sa->position, r_a));
				v3 slide_axis = norm(rotate(sa->rotation, w->joints[s->joint_idx].prismatic.local_axis_a));
				v3 lateral_error = sub(delta, scale(slide_axis, dot(delta, slide_axis)));
				float err_len = len(lateral_error);
				if (err_len < 1e-6f) continue;
				float inv_mass_sum = a->inv_mass + b->inv_mass;
				if (inv_mass_sum <= 0) continue;
				float correction = SOLVER_POS_BAUMGARTE * err_len;
				if (correction > SOLVER_POS_MAX_CORRECTION) correction = SOLVER_POS_MAX_CORRECTION;
				v3 n = scale(lateral_error, 1.0f / err_len);
				float P = correction / inv_mass_sum;
				sa->position = add(sa->position, scale(n, P * a->inv_mass));
				sb->position = sub(sb->position, scale(n, P * b->inv_mass));
			} else if (s->type == JOINT_DISTANCE) {
				v3 r_a = rotate(sa->rotation, w->joints[s->joint_idx].distance.local_a);
				v3 r_b = rotate(sb->rotation, w->joints[s->joint_idx].distance.local_b);
				v3 delta = sub(add(sb->position, r_b), add(sa->position, r_a));
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
				sa->position = add(sa->position, scale(axis, P * a->inv_mass));
				sb->position = sub(sb->position, scale(axis, P * b->inv_mass));
			}
		}
	}
}

// Store joint accumulated impulses back to persistent storage.
static void joints_post_solve(WorldInternal* w, SolverJoint* joints, int count)
{
	for (int i = 0; i < count; i++) {
		SolverJoint* s = &joints[i];
		for (int d = 0; d < s->dof; d++) w->joint_hot[s->joint_idx].warm_lambda[d] = s->lambda[d];
	}
	// joints NOT freed -- kept in WorldInternal.dbg_solver_joints for remote viewer.
}
