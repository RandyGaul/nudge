// tests_jacobian_unit.c -- unit tests for ldl_fill_jacobian.
// Verifies Jacobian rows for each constraint type: ball-socket, distance, hinge, weld.
// Tests known geometries, axis-aligned cases, and degenerate configurations.

#define JAC_EPS 1e-6

// Helper: fill SolverJoint.rows[] for ball socket from (r_a, r_b).
// J_a = [-I, -skew(r_a)], J_b = [I, skew(r_b)] where skew(r)*w = r x w.
static void test_fill_bs_rows(SolverJoint* s)
{
	v3 ra = s->r_a, rb = s->r_b;
	for (int d = 0; d < 3; d++) { memset(s->rows[d].J_a, 0, 6 * sizeof(float)); memset(s->rows[d].J_b, 0, 6 * sizeof(float)); }
	s->rows[0].J_a[0] = -1; s->rows[0].J_a[4] = -ra.z; s->rows[0].J_a[5] =  ra.y;
	s->rows[0].J_b[0] =  1; s->rows[0].J_b[4] =  rb.z; s->rows[0].J_b[5] = -rb.y;
	s->rows[1].J_a[1] = -1; s->rows[1].J_a[3] =  ra.z; s->rows[1].J_a[5] = -ra.x;
	s->rows[1].J_b[1] =  1; s->rows[1].J_b[3] = -rb.z; s->rows[1].J_b[5] =  rb.x;
	s->rows[2].J_a[2] = -1; s->rows[2].J_a[3] = -ra.y; s->rows[2].J_a[4] =  ra.x;
	s->rows[2].J_b[2] =  1; s->rows[2].J_b[3] =  rb.y; s->rows[2].J_b[4] = -rb.x;
}

// Helper: fill SolverJoint.rows[] for distance from (r_a, r_b, axis).
static void test_fill_dist_rows(SolverJoint* s, v3 axis)
{
	memset(s->rows[0].J_a, 0, 6 * sizeof(float)); memset(s->rows[0].J_b, 0, 6 * sizeof(float));
	v3 rxa = cross(s->r_a, axis), rxb = cross(s->r_b, axis);
	s->rows[0].J_a[0] = -axis.x; s->rows[0].J_a[1] = -axis.y; s->rows[0].J_a[2] = -axis.z;
	s->rows[0].J_a[3] = -rxa.x;  s->rows[0].J_a[4] = -rxa.y;  s->rows[0].J_a[5] = -rxa.z;
	s->rows[0].J_b[0] =  axis.x; s->rows[0].J_b[1] =  axis.y; s->rows[0].J_b[2] =  axis.z;
	s->rows[0].J_b[3] =  rxb.x;  s->rows[0].J_b[4] =  rxb.y;  s->rows[0].J_b[5] =  rxb.z;
}

// Helper: fill SolverJoint.rows[] for hinge from (r_a, r_b, u1, u2).
static void test_fill_hinge_rows(SolverJoint* s, v3 u1, v3 u2)
{
	// Linear rows 0-2: same as ball socket
	test_fill_bs_rows(s);
	// Angular rows 3-4
	for (int d = 3; d < 5; d++) { memset(s->rows[d].J_a, 0, 6 * sizeof(float)); memset(s->rows[d].J_b, 0, 6 * sizeof(float)); }
	s->rows[3].J_a[3] =  u1.x; s->rows[3].J_a[4] =  u1.y; s->rows[3].J_a[5] =  u1.z;
	s->rows[3].J_b[3] = -u1.x; s->rows[3].J_b[4] = -u1.y; s->rows[3].J_b[5] = -u1.z;
	s->rows[4].J_a[3] =  u2.x; s->rows[4].J_a[4] =  u2.y; s->rows[4].J_a[5] =  u2.z;
	s->rows[4].J_b[3] = -u2.x; s->rows[4].J_b[4] = -u2.y; s->rows[4].J_b[5] = -u2.z;
}

// Helper: verify a Jacobian row satisfies J*v = constraint_velocity for given body velocities.
// This is the fundamental property: the Jacobian maps body velocities to constraint-space velocity.
// v_a = (lin_a, ang_a), v_b = (lin_b, ang_b). Returns J_a*v_a + J_b*v_b.
static double jac_velocity(LDL_JacobianRow* jac, v3 lin_a, v3 ang_a, v3 lin_b, v3 ang_b)
{
	return jac->J_a[0]*lin_a.x + jac->J_a[1]*lin_a.y + jac->J_a[2]*lin_a.z + jac->J_a[3]*ang_a.x + jac->J_a[4]*ang_a.y + jac->J_a[5]*ang_a.z + jac->J_b[0]*lin_b.x + jac->J_b[1]*lin_b.y + jac->J_b[2]*lin_b.z + jac->J_b[3]*ang_b.x + jac->J_b[4]*ang_b.y + jac->J_b[5]*ang_b.z;
}

// Helper: verify J_a and J_b have opposite sign structure for linear DOFs.
// For ball-socket and hinge linear rows: J_a_lin = -I, J_b_lin = +I.
static int jac_linear_antisymmetric(LDL_JacobianRow* jac, int row)
{
	return (fabs(jac[row].J_a[row] - (-1.0)) < JAC_EPS) && (fabs(jac[row].J_b[row] - 1.0) < JAC_EPS);
}

// ============================================================================
// Ball-socket Jacobian tests
//
// J_a = [-I_3, skew(r_a)],  J_b = [I_3, -skew(r_b)]
// where skew(r) maps r to the matrix such that skew(r)*w = r x w.
//
// Physical meaning: constrains the world-space anchor points to coincide.
// The constraint velocity is: v_b + w_b x r_b - v_a - w_a x r_a = 0.

static void test_jac_ball_socket_origin()
{
	TEST_BEGIN("jac_ball_socket_origin");
	// Both anchors at body origins (r_a = r_b = 0).
	// J_a = [-I, 0], J_b = [I, 0]. Pure linear, no angular coupling.
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	for (int d = 0; d < 3; d++) {
		TEST_ASSERT(jac_linear_antisymmetric(jac, d));
		// Angular columns must be zero
		TEST_ASSERT(fabs(jac[d].J_a[3]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[d].J_a[4]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[d].J_a[5]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[d].J_b[3]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[d].J_b[4]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[d].J_b[5]) < JAC_EPS);
	}
}

static void test_jac_ball_socket_offset()
{
	TEST_BEGIN("jac_ball_socket_offset");
	// r_a = (1, 0, 0), r_b = (0, 1, 0).
	// J_a angular part = skew(r_a) = skew(1,0,0).
	// skew(1,0,0) = [0 0 0; 0 0 -1; 0 1 0]  (column-wise: skew*e_x=0, skew*e_y=(0,0,1), skew*e_z=(0,-1,0))
	// Wait -- the convention is J_a = [-I, skew(r_a)] where skew(r)*w = r x w.
	// Row 0 angular: [0, -r_a.z, r_a.y] = [0, 0, 0]
	// Row 1 angular: [r_a.z, 0, -r_a.x] = [0, 0, -1]
	// Row 2 angular: [-r_a.y, r_a.x, 0] = [0, 1, 0]
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,1,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Row 0: J_a = [-1, 0, 0,   0, -0, 0] = [-1, 0, 0, 0, 0, 0]
	TEST_ASSERT(fabs(jac[0].J_a[0] + 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[4]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[5]) < JAC_EPS);

	// Row 1: J_a = [0, -1, 0,   0, 0, -1]
	TEST_ASSERT(fabs(jac[1].J_a[1] + 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[1].J_a[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[1].J_a[4]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[1].J_a[5] + 1.0) < JAC_EPS);

	// Row 2: J_a = [0, 0, -1,   0, 1, 0]
	TEST_ASSERT(fabs(jac[2].J_a[2] + 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[2].J_a[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[2].J_a[4] - 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[2].J_a[5]) < JAC_EPS);

	// J_b angular: -skew(r_b) = -skew(0,1,0)
	// Row 0: [r_b.z, -r_b.y] -> [0, 0, 0,  0, 0, -1]
	TEST_ASSERT(fabs(jac[0].J_b[4]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_b[5] + 1.0) < JAC_EPS);

	// Row 1: [-r_b.z, r_b.x] -> [0, 0, 0,  0, 0, 0]
	TEST_ASSERT(fabs(jac[1].J_b[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[1].J_b[5]) < JAC_EPS);

	// Row 2: [r_b.y, -r_b.x] -> [0, 0, 0,  1, 0, 0]
	TEST_ASSERT(fabs(jac[2].J_b[3] - 1.0) < JAC_EPS);
}

static void test_jac_ball_socket_velocity()
{
	TEST_BEGIN("jac_ball_socket_velocity");
	// Physical test: if body A has angular velocity w_a about Z,
	// with anchor at r_a = (1, 0, 0), then w_a x r_a = (0, 0, w) x (1, 0, 0) = (0, w, 0).
	// The constraint velocity row 1 (Y) should pick this up.
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Body A: stationary linear, rotating about Z at w=1.
	// Body B: stationary.
	// Constraint velocity = v_b + w_b x r_b - v_a - w_a x r_a
	//                     = 0 - 0 - 0 - (0,0,1) x (1,0,0)
	//                     = -(0, 1, 0) = (0, -1, 0)
	double cv0 = jac_velocity(&jac[0], V3(0,0,0), V3(0,0,1), V3(0,0,0), V3(0,0,0));
	double cv1 = jac_velocity(&jac[1], V3(0,0,0), V3(0,0,1), V3(0,0,0), V3(0,0,0));
	double cv2 = jac_velocity(&jac[2], V3(0,0,0), V3(0,0,1), V3(0,0,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv0) < JAC_EPS);
	TEST_ASSERT(fabs(cv1 + 1.0) < JAC_EPS); // -1 in Y
	TEST_ASSERT(fabs(cv2) < JAC_EPS);
}

static void test_jac_ball_socket_arbitrary()
{
	TEST_BEGIN("jac_ball_socket_arbitrary");
	// Arbitrary lever arms: verify via velocity test.
	// Both bodies translating and rotating, check constraint velocity formula.
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 2, 3), .r_b = V3(-1, 0.5f, 2) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	v3 va = V3(1, -1, 2), wa = V3(0.5f, -0.3f, 0.7f);
	v3 vb = V3(-0.5f, 1, 0), wb = V3(0.1f, 0.2f, -0.4f);

	// Expected: (vb + wb x rb) - (va + wa x ra) for each component
	v3 anchor_vel_a = add(va, cross(wa, sol.r_a));
	v3 anchor_vel_b = add(vb, cross(wb, sol.r_b));
	v3 expected = sub(anchor_vel_b, anchor_vel_a);

	double cv0 = jac_velocity(&jac[0], va, wa, vb, wb);
	double cv1 = jac_velocity(&jac[1], va, wa, vb, wb);
	double cv2 = jac_velocity(&jac[2], va, wa, vb, wb);
	TEST_ASSERT(fabs(cv0 - expected.x) < 1e-4);
	TEST_ASSERT(fabs(cv1 - expected.y) < 1e-4);
	TEST_ASSERT(fabs(cv2 - expected.z) < 1e-4);
}

static void test_jac_ball_socket_coincident()
{
	TEST_BEGIN("jac_ball_socket_coincident");
	// Degenerate: bodies at the same position, anchors at the same world point.
	// Jacobian should still be well-formed. With r_a = r_b, angular parts match.
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,1,1), .r_b = V3(1,1,1) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// With identical r, J_a angular = skew(r) and J_b angular = -skew(r).
	// Pure equal/opposite angular coupling.
	for (int d = 0; d < 3; d++) {
		for (int a = 3; a < 6; a++) {
			TEST_ASSERT(fabs(jac[d].J_a[a] + jac[d].J_b[a]) < JAC_EPS);
		}
	}
}

static void test_jac_ball_socket_large_lever()
{
	TEST_BEGIN("jac_ball_socket_large_lever");
	// Large lever arms (r = 100). Angular columns should have entries ~100.
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(100, 0, 0), .r_b = V3(0, 100, 0) };
	test_fill_bs_rows(&sol);
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Row 2, J_a angular: [-r_a.y, r_a.x, 0] = [0, 100, 0]
	TEST_ASSERT(fabs(jac[2].J_a[4] - 100.0) < JAC_EPS);
	// Row 0, J_b angular: [0, r_b.z, -r_b.y] = [0, 0, -100]
	TEST_ASSERT(fabs(jac[0].J_b[5] + 100.0) < JAC_EPS);

	// Velocity test still works at this scale
	v3 wa = V3(0, 0, 0.01f);
	v3 expected_a = cross(wa, sol.r_a); // (0,0,0.01) x (100,0,0) = (0, 1, 0)
	double cv1 = jac_velocity(&jac[1], V3(0,0,0), wa, V3(0,0,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv1 + expected_a.y) < 1e-4);
}

// ============================================================================
// Distance Jacobian tests
//
// J_a = [-axis^T, -(r_a x axis)^T],  J_b = [axis^T, (r_b x axis)^T]
// 1 DOF: constrains the distance between anchor points along axis.

static void test_jac_distance_x_axis()
{
	TEST_BEGIN("jac_distance_x_axis");
	// Axis along X, anchors at body origins.
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, V3(1,0,0));
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	// J_a = [-1, 0, 0, 0, 0, 0], J_b = [1, 0, 0, 0, 0, 0]
	TEST_ASSERT(fabs(jac[0].J_a[0] + 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[1]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[2]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_b[0] - 1.0) < JAC_EPS);
	// Angular parts zero (r = 0)
	for (int a = 3; a < 6; a++) {
		TEST_ASSERT(fabs(jac[0].J_a[a]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[0].J_b[a]) < JAC_EPS);
	}
}

static void test_jac_distance_diagonal_axis()
{
	TEST_BEGIN("jac_distance_diagonal_axis");
	// Axis along (1,1,1)/sqrt(3). Linear part should be [-axis, axis].
	v3 ax = norm(V3(1, 1, 1));
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, ax);
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	TEST_ASSERT(fabs(jac[0].J_a[0] + ax.x) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[1] + ax.y) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[2] + ax.z) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_b[0] - ax.x) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_b[1] - ax.y) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_b[2] - ax.z) < JAC_EPS);
}

static void test_jac_distance_with_lever()
{
	TEST_BEGIN("jac_distance_with_lever");
	// Axis along X, r_a = (0, 1, 0). Angular part = -(r_a x axis) = -(0,1,0) x (1,0,0) = -(0,0,-1) = (0,0,1).
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,1,0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, V3(1,0,0));
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	// Angular A: -(r_a x axis) = -((0,1,0) x (1,0,0)) = -(0*0-0*0, 0*1-1*0, 1*0-0*1) = -(0, 0, -1) = (0, 0, 1)
	TEST_ASSERT(fabs(jac[0].J_a[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[4]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[0].J_a[5] - 1.0) < JAC_EPS);
}

static void test_jac_distance_velocity()
{
	TEST_BEGIN("jac_distance_velocity");
	// Physical test: body B moving along the axis should produce positive constraint velocity.
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, V3(1,0,0));
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	// B moving +X at speed 5: constraint velocity = axis . (vb - va) = 1*5 = 5
	double cv = jac_velocity(&jac[0], V3(0,0,0), V3(0,0,0), V3(5,0,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv - 5.0) < JAC_EPS);

	// B moving perpendicular (Y): constraint velocity = 0
	cv = jac_velocity(&jac[0], V3(0,0,0), V3(0,0,0), V3(0,5,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv) < JAC_EPS);
}

static void test_jac_distance_arbitrary()
{
	TEST_BEGIN("jac_distance_arbitrary");
	// Arbitrary geometry, verify via velocity formula.
	v3 ra = V3(1, 2, -1), rb = V3(-0.5f, 1, 3);
	v3 ax = norm(V3(2, -1, 1));
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = ra, .r_b = rb };
	test_fill_dist_rows(&sol, ax);
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	v3 va = V3(1, -1, 2), wa = V3(0.5f, -0.3f, 0.7f);
	v3 vb = V3(-0.5f, 1, 0), wb = V3(0.1f, 0.2f, -0.4f);

	// Expected: axis . ((vb + wb x rb) - (va + wa x ra))
	v3 anchor_vel_a = add(va, cross(wa, ra));
	v3 anchor_vel_b = add(vb, cross(wb, rb));
	v3 diff = sub(anchor_vel_b, anchor_vel_a);
	double expected = dot(ax, diff);

	double cv = jac_velocity(&jac[0], va, wa, vb, wb);
	TEST_ASSERT(fabs(cv - expected) < 1e-4);
}

static void test_jac_distance_zero_axis()
{
	TEST_BEGIN("jac_distance_zero_axis");
	// Degenerate: axis = (0,0,0). This happens when anchors coincide.
	// Jacobian should be all zeros (no constraint direction).
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_dist_rows(&sol, V3(0,0,0));
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	for (int i = 0; i < 6; i++) {
		TEST_ASSERT(fabs(jac[0].J_a[i]) < JAC_EPS);
		TEST_ASSERT(fabs(jac[0].J_b[i]) < JAC_EPS);
	}
}

// ============================================================================
// Hinge Jacobian tests
//
// Rows 0-2: linear, same as ball-socket (3 DOF).
// Rows 3-4: angular, J_a = [0, u_d], J_b = [0, -u_d] where u1, u2 are
// angular constraint axes (perpendicular to hinge axis).

static void test_jac_hinge_linear_matches_ball_socket()
{
	TEST_BEGIN("jac_hinge_linear_matches_ball_socket");
	// Linear rows of hinge must exactly match ball-socket Jacobian.
	v3 ra = V3(1, 2, 3), rb = V3(-1, 0.5f, 2);
	SolverJoint sol_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = ra, .r_b = rb };
	test_fill_bs_rows(&sol_bs);
	SolverJoint sol_h = { .type = JOINT_HINGE, .dof = 5, .r_a = ra, .r_b = rb };
	test_fill_hinge_rows(&sol_h, V3(1,0,0), V3(0,1,0));
	LDL_Constraint con_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_Constraint con_h = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac_bs[3], jac_h[5];
	ldl_fill_jacobian(&con_bs, &sol_bs, jac_bs);
	ldl_fill_jacobian(&con_h, &sol_h, jac_h);

	for (int d = 0; d < 3; d++) {
		for (int i = 0; i < 6; i++) {
			TEST_ASSERT(fabs(jac_h[d].J_a[i] - jac_bs[d].J_a[i]) < JAC_EPS);
			TEST_ASSERT(fabs(jac_h[d].J_b[i] - jac_bs[d].J_b[i]) < JAC_EPS);
		}
	}
}

static void test_jac_hinge_angular_axis_aligned()
{
	TEST_BEGIN("jac_hinge_angular_axis_aligned");
	// Hinge axis along Z. u1, u2 are the angular constraint axes
	// (perpendicular to the hinge axis, constraining rotation off-axis).
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_hinge_rows(&sol, V3(1, 0, 0), V3(0, 1, 0));
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	// Row 3: J_a = [0,0,0, u1] = [0,0,0, 1,0,0], J_b = [0,0,0, -u1] = [0,0,0, -1,0,0]
	TEST_ASSERT(fabs(jac[3].J_a[0]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[3].J_a[3] - 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[3].J_a[4]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[3].J_a[5]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[3].J_b[3] + 1.0) < JAC_EPS);

	// Row 4: J_a = [0,0,0, u2] = [0,0,0, 0,1,0], J_b = [0,0,0, -u2]
	TEST_ASSERT(fabs(jac[4].J_a[3]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[4].J_a[4] - 1.0) < JAC_EPS);
	TEST_ASSERT(fabs(jac[4].J_a[5]) < JAC_EPS);
	TEST_ASSERT(fabs(jac[4].J_b[4] + 1.0) < JAC_EPS);
}

static void test_jac_hinge_angular_velocity()
{
	TEST_BEGIN("jac_hinge_angular_velocity");
	// Physical test: rotation about the hinge axis (Z) should produce zero
	// constraint velocity on angular rows. Rotation about X or Y should not.
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_hinge_rows(&sol, V3(1, 0, 0), V3(0, 1, 0));
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	// Body A rotating about Z (hinge axis): angular rows should be zero.
	double cv3 = jac_velocity(&jac[3], V3(0,0,0), V3(0,0,1), V3(0,0,0), V3(0,0,0));
	double cv4 = jac_velocity(&jac[4], V3(0,0,0), V3(0,0,1), V3(0,0,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv3) < JAC_EPS);
	TEST_ASSERT(fabs(cv4) < JAC_EPS);

	// Body A rotating about X (off-axis): angular row 3 (u1 = X) should pick it up.
	cv3 = jac_velocity(&jac[3], V3(0,0,0), V3(1,0,0), V3(0,0,0), V3(0,0,0));
	TEST_ASSERT(fabs(cv3 - 1.0) < JAC_EPS); // u1 . w_a = 1

	// Both rotating about X equally: angular constraint velocity = 0 (relative = 0).
	cv3 = jac_velocity(&jac[3], V3(0,0,0), V3(1,0,0), V3(0,0,0), V3(1,0,0));
	TEST_ASSERT(fabs(cv3) < JAC_EPS);
}

static void test_jac_hinge_arbitrary()
{
	TEST_BEGIN("jac_hinge_arbitrary");
	// Arbitrary geometry, verify angular rows via velocity formula.
	v3 u1 = norm(V3(1, 1, 0)), u2 = norm(V3(-1, 1, 0));
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5, .r_a = V3(2, -1, 0.5f), .r_b = V3(-1, 1, 3) };
	test_fill_hinge_rows(&sol, u1, u2);
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	v3 wa = V3(0.5f, -0.3f, 0.7f), wb = V3(0.1f, 0.2f, -0.4f);
	// Angular constraint velocity for row d: u_d . (w_a - w_b)
	v3 dw = sub(wa, wb);
	double expected3 = dot(u1, dw);
	double expected4 = dot(u2, dw);

	double cv3 = jac_velocity(&jac[3], V3(0,0,0), wa, V3(0,0,0), wb);
	double cv4 = jac_velocity(&jac[4], V3(0,0,0), wa, V3(0,0,0), wb);
	TEST_ASSERT(fabs(cv3 - expected3) < 1e-4);
	TEST_ASSERT(fabs(cv4 - expected4) < 1e-4);
}

static void test_jac_hinge_zero_lever()
{
	TEST_BEGIN("jac_hinge_zero_lever");
	// Anchors at body origins: linear rows are pure [-I, I], angular rows are pure u.
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	test_fill_hinge_rows(&sol, V3(0, 0, 1), V3(1, 0, 0));
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	// Linear rows: zero angular
	for (int d = 0; d < 3; d++)
		for (int a = 3; a < 6; a++) {
			TEST_ASSERT(fabs(jac[d].J_a[a]) < JAC_EPS);
			TEST_ASSERT(fabs(jac[d].J_b[a]) < JAC_EPS);
		}

	// Angular rows: zero linear
	for (int d = 3; d < 5; d++)
		for (int l = 0; l < 3; l++) {
			TEST_ASSERT(fabs(jac[d].J_a[l]) < JAC_EPS);
			TEST_ASSERT(fabs(jac[d].J_b[l]) < JAC_EPS);
		}
}

// ============================================================================
// Synthetic weld Jacobian tests

static void test_jac_weld()
{
	TEST_BEGIN("jac_weld");
	// Synthetic weld: J_a = -I_6, J_b = +I_6. 6 DOF, all independent.
	LDL_Constraint con = { .type = -1, .dof = 6, .is_synthetic = 1 };
	LDL_JacobianRow jac[6];
	ldl_fill_jacobian(&con, NULL, jac);

	for (int d = 0; d < 6; d++) {
		for (int i = 0; i < 6; i++) {
			double expected_a = (i == d) ? -1.0 : 0.0;
			double expected_b = (i == d) ? 1.0 : 0.0;
			TEST_ASSERT(fabs(jac[d].J_a[i] - expected_a) < JAC_EPS);
			TEST_ASSERT(fabs(jac[d].J_b[i] - expected_b) < JAC_EPS);
		}
	}
}

static void test_jac_weld_velocity()
{
	TEST_BEGIN("jac_weld_velocity");
	// Weld: constraint velocity = v_b - v_a (all 6 DOFs).
	// If both have same velocity, constraint velocity is zero.
	LDL_Constraint con = { .type = -1, .dof = 6, .is_synthetic = 1 };
	LDL_JacobianRow jac[6];
	ldl_fill_jacobian(&con, NULL, jac);

	v3 v = V3(3, -2, 7), w = V3(1, -1, 0.5f);
	for (int d = 0; d < 6; d++) {
		double cv = jac_velocity(&jac[d], v, w, v, w);
		TEST_ASSERT(fabs(cv) < JAC_EPS);
	}

	// Different velocities: constraint velocity = difference
	v3 va = V3(1, 0, 0), wa = V3(0, 0, 0);
	v3 vb = V3(0, 0, 0), wb = V3(0, 0, 0);
	// Row 0 (lin X): J_a[0]=-1, J_b[0]=1. cv = -1*1 + 1*0 = -1
	double cv0 = jac_velocity(&jac[0], va, wa, vb, wb);
	TEST_ASSERT(fabs(cv0 + 1.0) < JAC_EPS);
}

// ============================================================================
// Runner

static void run_jacobian_unit_tests()
{
	printf("--- jacobian unit tests ---\n");

	// Ball-socket
	test_jac_ball_socket_origin();
	test_jac_ball_socket_offset();
	test_jac_ball_socket_velocity();
	test_jac_ball_socket_arbitrary();
	test_jac_ball_socket_coincident();
	test_jac_ball_socket_large_lever();

	// Distance
	test_jac_distance_x_axis();
	test_jac_distance_diagonal_axis();
	test_jac_distance_with_lever();
	test_jac_distance_velocity();
	test_jac_distance_arbitrary();
	test_jac_distance_zero_axis();

	// Hinge
	test_jac_hinge_linear_matches_ball_socket();
	test_jac_hinge_angular_axis_aligned();
	test_jac_hinge_angular_velocity();
	test_jac_hinge_arbitrary();
	test_jac_hinge_zero_lever();

	// Synthetic weld
	test_jac_weld();
	test_jac_weld_velocity();
}
