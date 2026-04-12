// tests_impulse_unit.c -- unit tests for ldl_apply_jacobian_impulse and ldl_constraint_velocity.
// These are the "apply" side of the solver: impulse application and constraint velocity readback.

#define IMP_EPS 1e-4f

// Helper: compute reference delta_v = inv_mass * J_lin^T * lambda (linear part).
static v3 ref_linear_impulse(LDL_JacobianRow* jac, int dof, int side, double* lambda, float inv_mass)
{
	double dx = 0, dy = 0, dz = 0;
	for (int d = 0; d < dof; d++) {
		double* J = side ? jac[d].J_b : jac[d].J_a;
		dx += inv_mass * J[0] * lambda[d];
		dy += inv_mass * J[1] * lambda[d];
		dz += inv_mass * J[2] * lambda[d];
	}
	return V3((float)dx, (float)dy, (float)dz);
}

// Helper: compute reference delta_w = I_world_inv * J_ang^T * lambda (angular part).
static v3 ref_angular_impulse(LDL_JacobianRow* jac, int dof, int side, double* lambda, BodyHot* body)
{
	v3 total = V3(0, 0, 0);
	for (int d = 0; d < dof; d++) {
		double* J = side ? jac[d].J_b : jac[d].J_a;
		v3 j_ang = V3((float)(J[3] * lambda[d]), (float)(J[4] * lambda[d]), (float)(J[5] * lambda[d]));
		v3 dw = inv_inertia_mul(body->rotation, body->inv_inertia_local, j_ang);
		total = add(total, dw);
	}
	return total;
}

// ============================================================================
// ldl_apply_jacobian_impulse tests

static void test_impulse_ball_socket_unit()
{
	TEST_BEGIN("impulse_ball_socket_unit");
	// Unit mass, anchor at origin. lambda = (1,0,0).
	// J_a = [-I, 0]. delta_v = inv_mass * J_a^T * lambda = 1 * [-1, 0, 0]^T * 1 = (-1, 0, 0).
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	TEST_ASSERT_FLOAT(body.velocity.x, -1.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, 0.0f, IMP_EPS);
	// No angular part (r = 0)
	TEST_ASSERT_FLOAT(body.angular_velocity.x, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, 0.0f, IMP_EPS);
}

static void test_impulse_ball_socket_side_b()
{
	TEST_BEGIN("impulse_ball_socket_side_b");
	// Side B: J_b = [+I, -skew(r_b)]. lambda = (1,0,0).
	// With r_b = 0: delta_v = inv_mass * [1, 0, 0]^T * 1 = (1, 0, 0).
	// Opposite sign from side A.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 1);

	TEST_ASSERT_FLOAT(body.velocity.x, 1.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, 0.0f, IMP_EPS);
}

static void test_impulse_ball_socket_with_lever()
{
	TEST_BEGIN("impulse_ball_socket_with_lever");
	// Lever arm r_a = (0, 0, 1). lambda = (1, 0, 0).
	// Linear: delta_v = inv_mass * J_a_lin^T * lambda = (-1, 0, 0).
	// Angular: J_a row 0 angular = [0, -r_a.z, r_a.y] = [0, -1, 0].
	// delta_w = I_inv * [0, -1, 0] * 1.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0, 0, 1), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	v3 ref_v = ref_linear_impulse(jac, 3, 0, lambda, body.inv_mass);
	v3 ref_w = ref_angular_impulse(jac, 3, 0, lambda, &body);

	TEST_ASSERT_FLOAT(body.velocity.x, ref_v.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, ref_v.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, ref_v.z, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, ref_w.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, ref_w.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, ref_w.z, IMP_EPS);
}

static void test_impulse_accumulates()
{
	TEST_BEGIN("impulse_accumulates");
	// Velocity starts at (1, 2, 3). Impulse should ADD to existing velocity.
	BodyHot body = make_body(1, 1);
	body.velocity = V3(1, 2, 3);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	// J_a = [-I, 0], so delta_v = (-1, 0, 0). Final = (0, 2, 3).
	TEST_ASSERT_FLOAT(body.velocity.x, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, 2.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, 3.0f, IMP_EPS);
}

static void test_impulse_heavy_body()
{
	TEST_BEGIN("impulse_heavy_body");
	// Heavy body: mass = 1000. Same impulse produces 1000x smaller velocity change.
	BodyHot body = make_body(1000, 1000);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 1, 1 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	TEST_ASSERT_FLOAT(body.velocity.x, -0.001f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, -0.001f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, -0.001f, IMP_EPS);
}

static void test_impulse_static_body()
{
	TEST_BEGIN("impulse_static_body");
	// Static body: inv_mass = 0, inv_inertia = 0. No velocity change.
	BodyHot body = {0};
	body.rotation = quat_identity();
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,2,3), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 100, 200, 300 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	TEST_ASSERT_FLOAT(body.velocity.x, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, 0.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, 0.0f, IMP_EPS);
}

static void test_impulse_zero_lambda()
{
	TEST_BEGIN("impulse_zero_lambda");
	// Zero impulse: no velocity change.
	BodyHot body = make_body(1, 1);
	body.velocity = V3(5, -3, 7);
	body.angular_velocity = V3(1, 2, 3);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 0, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	TEST_ASSERT_FLOAT(body.velocity.x, 5.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, -3.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, 7.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, 1.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, 2.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, 3.0f, IMP_EPS);
}

static void test_impulse_distance()
{
	TEST_BEGIN("impulse_distance");
	// Distance constraint: 1 DOF. Impulse along the axis.
	BodyHot body = make_body(2, 4);
	v3 ax = norm(V3(1, 1, 0));
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(1, 0, 0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, ax);
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[1] = { 3.0 };
	ldl_apply_jacobian_impulse(jac, 1, lambda, &body, 0);

	v3 ref_v = ref_linear_impulse(jac, 1, 0, lambda, body.inv_mass);
	v3 ref_w = ref_angular_impulse(jac, 1, 0, lambda, &body);

	TEST_ASSERT_FLOAT(body.velocity.x, ref_v.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, ref_v.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, ref_v.z, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, ref_w.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, ref_w.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, ref_w.z, IMP_EPS);
}

static void test_impulse_hinge()
{
	TEST_BEGIN("impulse_hinge");
	// Hinge: 5 DOF (3 linear + 2 angular). Full lambda vector.
	BodyHot body = make_body(3, 6);
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	test_fill_hinge_rows(&sol, norm(V3(1, 0, 0)), norm(V3(0, 1, 0)));
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[5] = { 1, -0.5, 2, 0.3, -0.7 };
	ldl_apply_jacobian_impulse(jac, 5, lambda, &body, 0);

	v3 ref_v = ref_linear_impulse(jac, 5, 0, lambda, body.inv_mass);
	v3 ref_w = ref_angular_impulse(jac, 5, 0, lambda, &body);

	TEST_ASSERT_FLOAT(body.velocity.x, ref_v.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, ref_v.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, ref_v.z, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, ref_w.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, ref_w.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, ref_w.z, IMP_EPS);
}

static void test_impulse_rotated_asymmetric()
{
	TEST_BEGIN("impulse_rotated_asymmetric");
	// Rotated body with asymmetric inertia. The angular impulse goes through
	// dinv_inertia_mul, which we fixed earlier. Verify end-to-end.
	quat rot = quat_axis_angle(norm(V3(1, 2, -1)), 0.8f);
	BodyHot body = make_body_full(2, V3(1, 5, 20), rot);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(3, -1, 2), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 2.0, -1.5, 0.7 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	v3 ref_v = ref_linear_impulse(jac, 3, 0, lambda, body.inv_mass);
	v3 ref_w = ref_angular_impulse(jac, 3, 0, lambda, &body);

	TEST_ASSERT_FLOAT(body.velocity.x, ref_v.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.y, ref_v.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.velocity.z, ref_v.z, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, ref_w.x, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, ref_w.y, IMP_EPS);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, ref_w.z, IMP_EPS);
}

static void test_impulse_both_bodies()
{
	TEST_BEGIN("impulse_both_bodies");
	// Apply to both bodies of a ball-socket. Body A gets side 0, body B gets side 1.
	// After applying the same lambda to both, the constraint velocity should change
	// by lambda / effective_mass (the K matrix diagonal, roughly).
	BodyHot a = make_body(2, 4);
	BodyHot b = make_body(5, 10);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 1, 0, 0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &a, 0);
	ldl_apply_jacobian_impulse(jac, 3, lambda, &b, 1);

	// Verify against reference for each body independently
	BodyHot a_ref = make_body(2, 4);
	BodyHot b_ref = make_body(5, 10);
	v3 a_dv = ref_linear_impulse(jac, 3, 0, lambda, a_ref.inv_mass);
	v3 a_dw = ref_angular_impulse(jac, 3, 0, lambda, &a_ref);
	v3 b_dv = ref_linear_impulse(jac, 3, 1, lambda, b_ref.inv_mass);
	v3 b_dw = ref_angular_impulse(jac, 3, 1, lambda, &b_ref);

	TEST_ASSERT_FLOAT(a.velocity.x, a_dv.x, IMP_EPS);
	TEST_ASSERT_FLOAT(a.velocity.y, a_dv.y, IMP_EPS);
	TEST_ASSERT_FLOAT(a.velocity.z, a_dv.z, IMP_EPS);
	TEST_ASSERT_FLOAT(a.angular_velocity.x, a_dw.x, IMP_EPS);
	TEST_ASSERT_FLOAT(a.angular_velocity.y, a_dw.y, IMP_EPS);
	TEST_ASSERT_FLOAT(a.angular_velocity.z, a_dw.z, IMP_EPS);
	TEST_ASSERT_FLOAT(b.velocity.x, b_dv.x, IMP_EPS);
	TEST_ASSERT_FLOAT(b.velocity.y, b_dv.y, IMP_EPS);
	TEST_ASSERT_FLOAT(b.velocity.z, b_dv.z, IMP_EPS);
	TEST_ASSERT_FLOAT(b.angular_velocity.x, b_dw.x, IMP_EPS);
	TEST_ASSERT_FLOAT(b.angular_velocity.y, b_dw.y, IMP_EPS);
	TEST_ASSERT_FLOAT(b.angular_velocity.z, b_dw.z, IMP_EPS);
}

// ============================================================================
// ldl_constraint_velocity tests

static void test_cvel_stationary()
{
	TEST_BEGIN("cvel_stationary");
	// Both bodies stationary: constraint velocity = 0.
	BodyHot a = make_body(1, 1);
	BodyHot b = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,1,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	for (int d = 0; d < 3; d++) {
		double cv = ldl_constraint_velocity(&jac[d], &a, &b);
		TEST_ASSERT(fabs(cv) < 1e-10);
	}
}

static void test_cvel_linear_separation()
{
	TEST_BEGIN("cvel_linear_separation");
	// Body B moving away from A along X. Anchors at origin.
	// Constraint velocity row 0 (X) should be positive.
	BodyHot a = make_body(1, 1);
	BodyHot b = make_body(1, 1);
	b.velocity = V3(5, 0, 0);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double cv0 = ldl_constraint_velocity(&jac[0], &a, &b);
	double cv1 = ldl_constraint_velocity(&jac[1], &a, &b);
	double cv2 = ldl_constraint_velocity(&jac[2], &a, &b);
	TEST_ASSERT(fabs(cv0 - 5.0) < 1e-10);
	TEST_ASSERT(fabs(cv1) < 1e-10);
	TEST_ASSERT(fabs(cv2) < 1e-10);
}

static void test_cvel_angular()
{
	TEST_BEGIN("cvel_angular");
	// Body A rotating about Z at w=1, lever arm r_a = (1, 0, 0).
	// Anchor velocity = w x r = (0,0,1) x (1,0,0) = (0, 1, 0).
	// Constraint velocity Y row should pick this up as -1 (J_a convention).
	BodyHot a = make_body(1, 1);
	a.angular_velocity = V3(0, 0, 1);
	BodyHot b = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double cv1 = ldl_constraint_velocity(&jac[1], &a, &b);
	// J_a row 1: [-1 in Y, angular=(r_a.z, 0, -r_a.x)=(0,0,-1)].
	// cv = J_a * v_a = 0*0 + (-1)*0 + 0*0 + 0*0 + 0*0 + (-1)*1 = -1.
	TEST_ASSERT(fabs(cv1 + 1.0) < 1e-6);
}

static void test_impulse_large_lambda()
{
	TEST_BEGIN("impulse_large_lambda");
	// Stretched constraint: lambda = 5000. Light body (mass=1) gets a huge kick.
	// Verify no precision loss in the velocity delta.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 5000, -3000, 2000 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	v3 ref_v = ref_linear_impulse(jac, 3, 0, lambda, body.inv_mass);
	v3 ref_w = ref_angular_impulse(jac, 3, 0, lambda, &body);

	TEST_ASSERT_FLOAT(body.velocity.x, ref_v.x, 1.0f); // absolute tolerance scales with magnitude
	TEST_ASSERT_FLOAT(body.velocity.y, ref_v.y, 1.0f);
	TEST_ASSERT_FLOAT(body.velocity.z, ref_v.z, 1.0f);
	TEST_ASSERT_FLOAT(body.angular_velocity.x, ref_w.x, 1.0f);
	TEST_ASSERT_FLOAT(body.angular_velocity.y, ref_w.y, 1.0f);
	TEST_ASSERT_FLOAT(body.angular_velocity.z, ref_w.z, 1.0f);
}

static void test_impulse_large_lambda_on_existing_velocity()
{
	TEST_BEGIN("impulse_large_lambda_on_existing_velocity");
	// Large impulse added to already-large velocity. Tests float32 addition
	// precision: adding 5000 to 10000 should not lose the small digits.
	BodyHot body = make_body(1, 1);
	body.velocity = V3(10000, -10000, 5000);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 5000, -3000, 1000 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &body, 0);

	// J_a = [-I, 0], delta_v = (-5000, 3000, -1000)
	TEST_ASSERT_FLOAT(body.velocity.x, 5000.0f, 1.0f);
	TEST_ASSERT_FLOAT(body.velocity.y, -7000.0f, 1.0f);
	TEST_ASSERT_FLOAT(body.velocity.z, 4000.0f, 1.0f);
}

static void test_impulse_extreme_mass_ratio()
{
	TEST_BEGIN("impulse_extreme_mass_ratio");
	// 10000:1 mass ratio. Same lambda applied to both bodies.
	// Light body: huge velocity change. Heavy body: tiny change.
	BodyHot light = make_body(1, 1);
	BodyHot heavy = make_body(10000, 10000);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double lambda[3] = { 10, 10, 10 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &light, 0);
	ldl_apply_jacobian_impulse(jac, 3, lambda, &heavy, 1);

	// Light: delta_v = inv_mass * J_a_lin^T * lambda = 1 * (-10, -10, -10)
	TEST_ASSERT_FLOAT(light.velocity.x, -10.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(light.velocity.y, -10.0f, IMP_EPS);
	TEST_ASSERT_FLOAT(light.velocity.z, -10.0f, IMP_EPS);

	// Heavy: delta_v = 0.0001 * (10, 10, 10) = (0.001, 0.001, 0.001)
	TEST_ASSERT_FLOAT(heavy.velocity.x, 0.001f, IMP_EPS);
	TEST_ASSERT_FLOAT(heavy.velocity.y, 0.001f, IMP_EPS);
	TEST_ASSERT_FLOAT(heavy.velocity.z, 0.001f, IMP_EPS);
}

static void test_cvel_large_velocities()
{
	TEST_BEGIN("cvel_large_velocities");
	// Both bodies moving fast (post-impulse from stretched constraint).
	// Large velocity + large angular Jacobian entries: precision test.
	BodyHot a = make_body(1, 1);
	BodyHot b = make_body(1, 1);
	a.velocity = V3(-5000, 3000, -1000);
	a.angular_velocity = V3(2000, -1500, 800);
	b.velocity = V3(4000, -2000, 3000);
	b.angular_velocity = V3(-1000, 500, -200);

	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(2, -1, 3), .r_b = V3(-1, 2, 0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	for (int d = 0; d < 3; d++) {
		double cv = ldl_constraint_velocity(&jac[d], &a, &b);
		// Verify against manual computation
		double* Ja = jac[d].J_a;
		double* Jb = jac[d].J_b;
		double expected = Ja[0]*a.velocity.x + Ja[1]*a.velocity.y + Ja[2]*a.velocity.z + Ja[3]*a.angular_velocity.x + Ja[4]*a.angular_velocity.y + Ja[5]*a.angular_velocity.z + Jb[0]*b.velocity.x + Jb[1]*b.velocity.y + Jb[2]*b.velocity.z + Jb[3]*b.angular_velocity.x + Jb[4]*b.angular_velocity.y + Jb[5]*b.angular_velocity.z;
		TEST_ASSERT(fabs(cv - expected) < 1.0); // large values, absolute tolerance
	}
}

static void test_cvel_roundtrip_extreme_mass_ratio()
{
	TEST_BEGIN("cvel_roundtrip_extreme_mass_ratio");
	// Roundtrip with 10000:1 mass ratio: apply lambda, verify cv = K * lambda.
	// Heavy body contributes ~0.0001 to cv, light contributes ~1.
	// The constraint velocity is the sum of tiny + large -- cancellation test.
	BodyHot a = make_body(1, 1);       // light
	BodyHot b = make_body(10000, 10000); // heavy
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Build K
	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &b, 1.0, K);

	// Apply lambda
	double lambda[3] = { 1.0, -0.5, 2.0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &a, 0);
	ldl_apply_jacobian_impulse(jac, 3, lambda, &b, 1);

	// Read constraint velocity
	double cv[3];
	for (int d = 0; d < 3; d++) {
		cv[d] = ldl_constraint_velocity(&jac[d], &a, &b);
	}

	// Expected: K * lambda
	double expected[3] = {0};
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			expected[r] += K[LDL_TRI(r, c)] * lambda[c];
		}
	}

	TEST_ASSERT(fabs(cv[0] - expected[0]) < 0.01);
	TEST_ASSERT(fabs(cv[1] - expected[1]) < 0.01);
	TEST_ASSERT(fabs(cv[2] - expected[2]) < 0.01);
}

static void test_cvel_roundtrip_with_lever()
{
	TEST_BEGIN("cvel_roundtrip_with_lever");
	// Roundtrip with large lever arms on both sides. Angular contributions
	// dominate, and the K matrix has large off-diagonals.
	BodyHot a = make_body(2, 3);
	BodyHot b = make_body(5, 8);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(10, -5, 3), .r_b = V3(-3, 8, -2) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &b, 1.0, K);

	double lambda[3] = { 3.0, -1.0, 0.5 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &a, 0);
	ldl_apply_jacobian_impulse(jac, 3, lambda, &b, 1);

	double cv[3];
	for (int d = 0; d < 3; d++) {
		cv[d] = ldl_constraint_velocity(&jac[d], &a, &b);
	}

	double expected[3] = {0};
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			expected[r] += K[LDL_TRI(r, c)] * lambda[c];
		}
	}

	// Looser tolerance: float32 accumulation with large lever arms
	for (int i = 0; i < 3; i++) {
		TEST_ASSERT(fabs(cv[i] - expected[i]) / (fabs(expected[i]) + 0.01) < 0.01);
	}
}

static void test_cvel_impulse_roundtrip()
{
	TEST_BEGIN("cvel_impulse_roundtrip");
	// Key property: applying lambda through K should change constraint velocity
	// by K * lambda. This is the fundamental consistency check between
	// K assembly, impulse application, and constraint velocity readback.
	//
	// Start with zero velocity. Apply lambda. Read constraint velocity.
	// It should equal K * lambda.
	BodyHot a = make_body(2, 4);
	BodyHot b = make_body(5, 10);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Build K
	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &b, 1.0, K);

	// Apply lambda
	double lambda[3] = { 1.0, -0.5, 2.0 };
	ldl_apply_jacobian_impulse(jac, 3, lambda, &a, 0);
	ldl_apply_jacobian_impulse(jac, 3, lambda, &b, 1);

	// Read constraint velocity
	double cv[3];
	for (int d = 0; d < 3; d++) {
		cv[d] = ldl_constraint_velocity(&jac[d], &a, &b);
	}

	// Expected: K * lambda
	double expected[3] = {0};
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			expected[r] += K[LDL_TRI(r, c)] * lambda[c];
		}
	}

	TEST_ASSERT(fabs(cv[0] - expected[0]) < 0.01);
	TEST_ASSERT(fabs(cv[1] - expected[1]) < 0.01);
	TEST_ASSERT(fabs(cv[2] - expected[2]) < 0.01);
}

// ============================================================================
// Runner

static void run_impulse_unit_tests()
{
	printf("--- impulse unit tests ---\n");

	// Apply impulse
	test_impulse_ball_socket_unit();
	test_impulse_ball_socket_side_b();
	test_impulse_ball_socket_with_lever();
	test_impulse_accumulates();
	test_impulse_heavy_body();
	test_impulse_static_body();
	test_impulse_zero_lambda();
	test_impulse_distance();
	test_impulse_hinge();
	test_impulse_rotated_asymmetric();
	test_impulse_both_bodies();
	test_impulse_large_lambda();
	test_impulse_large_lambda_on_existing_velocity();
	test_impulse_extreme_mass_ratio();

	// Constraint velocity
	test_cvel_stationary();
	test_cvel_linear_separation();
	test_cvel_angular();
	test_cvel_large_velocities();
	test_cvel_roundtrip_extreme_mass_ratio();
	test_cvel_roundtrip_with_lever();
	test_cvel_impulse_roundtrip();
}
