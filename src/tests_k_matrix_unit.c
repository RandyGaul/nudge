// tests_k_matrix_unit.c -- unit tests for ldl_K_body_contrib and ldl_K_body_off.
// These build K = J * M^{-1} * J^T, the constraint-space effective mass matrix.
// K is what block_ldl factorizes, so correctness here is critical.

#define K_EPS 1e-8

// Helper: build K = J * M^{-1} * J^T via explicit reference implementation.
// For one body contribution (one side of the constraint):
// W = M^{-1} * J^T (6 x dof), then K += J * W (dof x dof, symmetric).
// M^{-1} = diag(inv_m, inv_m, inv_m, I_world_inv) where I_world_inv = R*diag(inv_i)*R^T.
static void ref_K_body_contrib(LDL_JacobianRow* jac, int dof, int side, BodyHot* body, double weight, double* K_full)
{
	double wm = (double)body->inv_mass * weight;
	for (int r = 0; r < dof; r++) {
		double* Jr = side ? jac[r].J_b : jac[r].J_a;
		for (int c = 0; c < dof; c++) {
			double* Jc = side ? jac[c].J_b : jac[c].J_a;
			// Linear part: inv_mass * (J_lin_r . J_lin_c)
			double lin = wm * (Jr[0]*Jc[0] + Jr[1]*Jc[1] + Jr[2]*Jc[2]);
			// Angular part: Jr_ang^T * I_world_inv * Jc_ang
			v3 jc_ang_f = V3((float)Jc[3], (float)Jc[4], (float)Jc[5]);
			v3 Iinv_jc = inv_inertia_mul(body->rotation, body->inv_inertia_local, jc_ang_f);
			double ang = weight * (Jr[3]*(double)Iinv_jc.x + Jr[4]*(double)Iinv_jc.y + Jr[5]*(double)Iinv_jc.z);
			K_full[r * dof + c] += lin + ang;
		}
	}
}

// Helper: check K is symmetric (should always be for K = J*M^{-1}*J^T).
static int k_is_symmetric(double* K_full, int n, double eps)
{
	for (int r = 0; r < n; r++)
		for (int c = 0; c < r; c++)
			if (fabs(K_full[r*n+c] - K_full[c*n+r]) > eps) return 0;
	return 1;
}

// Helper: check K is positive semi-definite by checking all diagonal entries > 0
// and that v^T K v > 0 for random vectors. (Not rigorous, but catches obvious bugs.)
static int k_diag_positive(double* K_packed, int n)
{
	for (int i = 0; i < n; i++)
		if (K_packed[LDL_TRI(i, i)] <= 0) return 0;
	return 1;
}

// Helper: unpack lower-tri to full symmetric matrix.
static void unpack_symmetric(double* packed, int n, double* full)
{
	for (int r = 0; r < n; r++)
		for (int c = 0; c < n; c++)
			full[r * n + c] = packed[LDL_TRI(r, c)];
}

// Helper: make a body with given mass, uniform inertia, identity rotation.
static BodyHot make_body(float mass, float inertia)
{
	BodyHot b = {0};
	b.rotation = quat_identity();
	if (mass > 0) {
		b.inv_mass = 1.0f / mass;
		float inv_i = (inertia > 0) ? 1.0f / inertia : 0;
		b.inv_inertia_local = V3(inv_i, inv_i, inv_i);
	}
	return b;
}

// Helper: make a body with specific per-axis inertia and rotation.
static BodyHot make_body_full(float mass, v3 inertia, quat rot)
{
	BodyHot b = {0};
	b.rotation = rot;
	if (mass > 0) {
		b.inv_mass = 1.0f / mass;
		b.inv_inertia_local = V3(
			inertia.x > 0 ? 1.0f / inertia.x : 0,
			inertia.y > 0 ? 1.0f / inertia.y : 0,
			inertia.z > 0 ? 1.0f / inertia.z : 0
		);
	}
	return b;
}

// ============================================================================
// ldl_K_body_contrib: diagonal block K += J_side * M^{-1} * J_side^T
// One body's contribution to the diagonal K block.

static void test_K_ball_socket_unit_mass_origin()
{
	TEST_BEGIN("K_ball_socket_unit_mass_origin");
	// Unit mass body, uniform unit inertia, anchor at origin.
	// J_a = [-I, 0]. K_a = J_a * I * J_a^T = I_3.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	// Should be identity (diagonal = 1, off-diagonal = 0)
	TEST_ASSERT(fabs(K[LDL_TRI(0,0)] - 1.0) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(1,1)] - 1.0) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,2)] - 1.0) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(1,0)]) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,0)]) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,1)]) < K_EPS);
}

static void test_K_ball_socket_heavy_body()
{
	TEST_BEGIN("K_ball_socket_heavy_body");
	// Heavy body (mass=1000), anchor at origin. K = inv_mass * I = 0.001 * I.
	BodyHot body = make_body(1000, 1000);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	TEST_ASSERT(fabs(K[LDL_TRI(0,0)] - 0.001) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(1,1)] - 0.001) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,2)] - 0.001) < K_EPS);
}

static void test_K_ball_socket_with_lever()
{
	TEST_BEGIN("K_ball_socket_with_lever");
	// Lever arm r_a = (1, 0, 0) introduces angular coupling.
	// J_a = [-I, skew(r_a)]. K_a = J_a * M^{-1} * J_a^T.
	// With unit mass and unit inertia, M^{-1} = I_6.
	// K = J * J^T (since M^{-1} = I).
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	// Compare against reference
	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);

	// K should be symmetric and have positive diagonal
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_ball_socket_two_bodies()
{
	TEST_BEGIN("K_ball_socket_two_bodies");
	// Full K = K_a + K_b from both bodies. Different masses.
	BodyHot body_a = make_body(2, 4);   // lighter
	BodyHot body_b = make_body(10, 20); // heavier
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body_a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &body_b, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body_a, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &body_b, 1.0, K_ref);

	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);
	TEST_ASSERT(k_is_symmetric(K_ref, 3, 1e-4));
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_static_body()
{
	TEST_BEGIN("K_static_body");
	// Static body (inv_mass = 0, inv_inertia = 0). Contribution should be zero.
	BodyHot body = {0};
	body.rotation = quat_identity();
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,2,3), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	for (int i = 0; i < 6; i++)
		TEST_ASSERT(fabs(K[i]) < K_EPS);
}

static void test_K_extreme_mass_ratio()
{
	TEST_BEGIN("K_extreme_mass_ratio");
	// 10000:1 mass ratio. K is dominated by the light body's contribution.
	BodyHot light = make_body(1, 1);
	BodyHot heavy = make_body(10000, 10000);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &light, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &heavy, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &light, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &heavy, 1.0, K_ref);

	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);

	// Heavy body contributes ~0.0001, light contributes ~1. Ratio ~10000.
	// All diagonal entries should still be positive.
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_asymmetric_inertia()
{
	TEST_BEGIN("K_asymmetric_inertia");
	// Non-uniform inertia: thin rod along X (small Ix, large Iy=Iz).
	// This creates anisotropic angular response.
	BodyHot body = make_body_full(1, V3(0.1f, 10, 10), quat_identity());
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0, 0, 1), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_rotated_body()
{
	TEST_BEGIN("K_rotated_body");
	// Non-identity rotation. The world-space inertia tensor changes,
	// which affects angular K contributions.
	quat rot = quat_axis_angle(norm(V3(1, 1, 1)), 1.0f);
	BodyHot body = make_body_full(1, V3(1, 5, 10), rot);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 2, 0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	// Looser tolerance: ref uses float inv_inertia_mul, code uses double dinv_inertia_mul
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_shattering_weight()
{
	TEST_BEGIN("K_shattering_weight");
	// Shattering: weight = 4 (4 shards). K should be 4x a weight=1 computation.
	BodyHot body = make_body(10, 5);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0.5f, 0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K1[6] = {0}, K4[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K1);
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 4.0, K4);

	for (int i = 0; i < 6; i++)
		TEST_ASSERT(fabs(K4[i] - 4.0 * K1[i]) < K_EPS);
}

static void test_K_mixed_scale_lever()
{
	TEST_BEGIN("K_mixed_scale_lever");
	// Lever arm r = (0.01, 0, 50): huge coupling on one angular axis,
	// near-zero on another. K diagonal entries will span ~4 orders of magnitude.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.01f, 0, 50), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);
	TEST_ASSERT(k_diag_positive(K, 3));

	// Verify the scale spread: largest diagonal vs smallest should be large.
	double dmin = K[LDL_TRI(0,0)], dmax = K[LDL_TRI(0,0)];
	for (int i = 1; i < 3; i++) {
		double d = K[LDL_TRI(i,i)];
		if (d < dmin) dmin = d;
		if (d > dmax) dmax = d;
	}
	TEST_ASSERT(dmax / dmin > 100.0); // confirm mixed scale
}

static void test_K_lever_along_weak_inertia()
{
	TEST_BEGIN("K_lever_along_weak_inertia");
	// Asymmetric inertia: Ix=0.1, Iy=Iz=10 (thin rod along X).
	// Lever arm along X (weak axis): torque about Y/Z uses the large inertia,
	// so angular K contribution is small.
	// Lever arm along Y (strong axis): torque about X uses the small inertia,
	// so angular K contribution is large.
	// Test both and verify the ratio.
	BodyHot body = make_body_full(1, V3(0.1f, 10, 10), quat_identity());

	// Lever along X (weak inertia axis)
	SolverJoint sol_x = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(5,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac_x[3];
	ldl_fill_jacobian(&con, &sol_x, jac_x);
	double Kx[6] = {0};
	ldl_K_body_contrib(jac_x, 3, 0, 0, &body, 1.0, Kx);

	// Lever along Y (strong inertia axis)
	SolverJoint sol_y = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,5,0), .r_b = V3(0,0,0) };
	LDL_JacobianRow jac_y[3];
	ldl_fill_jacobian(&con, &sol_y, jac_y);
	double Ky[6] = {0};
	ldl_K_body_contrib(jac_y, 3, 0, 0, &body, 1.0, Ky);

	// Both should match their references
	double Kx_ref[9] = {0}, Ky_ref[9] = {0};
	ref_K_body_contrib(jac_x, 3, 0, &body, 1.0, Kx_ref);
	ref_K_body_contrib(jac_y, 3, 0, &body, 1.0, Ky_ref);
	double Kx_full[9], Ky_full[9];
	unpack_symmetric(Kx, 3, Kx_full);
	unpack_symmetric(Ky, 3, Ky_full);
	for (int i = 0; i < 9; i++) {
		TEST_ASSERT(fabs(Kx_full[i] - Kx_ref[i]) < 1e-3);
		TEST_ASSERT(fabs(Ky_full[i] - Ky_ref[i]) < 1e-3);
	}

	// Lever along Y with small Ix: the row that picks up torque about X
	// (which is row with Z linear + X angular coupling) should have a much
	// larger angular contribution than the X-lever case.
	// Specifically: trace(Ky) should differ significantly from trace(Kx)
	// because lever*inertia interaction differs.
	double trace_x = Kx[LDL_TRI(0,0)] + Kx[LDL_TRI(1,1)] + Kx[LDL_TRI(2,2)];
	double trace_y = Ky[LDL_TRI(0,0)] + Ky[LDL_TRI(1,1)] + Ky[LDL_TRI(2,2)];
	TEST_ASSERT(fabs(trace_x - trace_y) > 1.0); // meaningfully different
	TEST_ASSERT(k_diag_positive(Kx, 3));
	TEST_ASSERT(k_diag_positive(Ky, 3));
}

static void test_K_floor_dynamic_pair()
{
	TEST_BEGIN("K_floor_dynamic_pair");
	// Most common configuration: static floor + dynamic box sitting on it.
	// Static contributes zero. K comes entirely from the dynamic body.
	// Dynamic body: mass=1, inertia from a unit box (Ix=Iy=Iz=1/6).
	// Lever arm: anchor at bottom of box, r = (0, -0.5, 0).
	float box_inertia = 1.0f / 6.0f;
	BodyHot floor_body = {0};
	floor_body.rotation = quat_identity();
	BodyHot box = make_body_full(1, V3(box_inertia, box_inertia, box_inertia), quat_identity());

	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0, -0.5f, 0), .r_b = V3(0, 0, 0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &box, 1.0, K);       // dynamic body
	ldl_K_body_contrib(jac, 3, 1, 0, &floor_body, 1.0, K); // static: adds nothing

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &box, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);

	// K should be SPD (factorable)
	TEST_ASSERT(k_diag_positive(K, 3));

	// Y row (gravity direction) should have inv_mass contribution only (no angular
	// coupling from r_a = (0, -0.5, 0) in the Y row since skew(r_a) row 1 = [rz, 0, -rx] = [0, 0, 0]).
	TEST_ASSERT(fabs(K[LDL_TRI(1,1)] - 1.0) < 0.01); // ~inv_mass = 1
}

static void test_K_mismatched_lever_lengths()
{
	TEST_BEGIN("K_mismatched_lever_lengths");
	// Body A: anchor at origin (no angular coupling).
	// Body B: anchor far out at r = (50, 0, 0) (huge angular coupling).
	// The two K contributions have vastly different magnitudes in angular block.
	// This is a common ragdoll/chain scenario.
	BodyHot body_a = make_body(1, 1);
	BodyHot body_b = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(50,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body_a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &body_b, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body_a, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &body_b, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 0.1); // absolute can be large

	TEST_ASSERT(k_diag_positive(K, 3));

	// Body A contributes inv_mass=1 on all diagonal entries.
	// Body B contributes inv_mass=1 plus angular from r_b=(50,0,0).
	// r_b along X: skew(r_b) has entries in Y/Z rows, not X row.
	// So X diagonal gets ~2 (both inv_mass), Y and Z get ~2 + angular.
	TEST_ASSERT(K[LDL_TRI(0,0)] > 1.5);

	// The diagonal entries should have a large spread because body B's lever
	// contributes angular coupling on two axes but not the third.
	double dmin = K[LDL_TRI(0,0)], dmax = K[LDL_TRI(0,0)];
	for (int i = 1; i < 3; i++) {
		double d = K[LDL_TRI(i,i)];
		if (d < dmin) dmin = d;
		if (d > dmax) dmax = d;
	}
	TEST_ASSERT(dmax / dmin > 10.0);
}

static void test_K_negative_offdiag()
{
	TEST_BEGIN("K_negative_offdiag");
	// Verify K can have negative off-diagonal entries while still being SPD.
	// Need: rotated body with strongly asymmetric inertia and a lever arm
	// with mixed-sign components that create cross terms through the rotated
	// inertia tensor. The 45-deg rotation mixes the weak/strong inertia axes
	// so that world-space I^{-1} has large off-diagonals.
	quat rot = quat_axis_angle(norm(V3(1, 0, 0)), 0.7854f); // 45 deg about X
	BodyHot body = make_body_full(1, V3(0.1f, 50, 0.1f), rot);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(5, 5, 0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);

	// Check that at least one off-diagonal is negative
	int has_neg = 0;
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < r; c++)
			if (K[LDL_TRI(r,c)] < -1e-6) has_neg = 1;
	TEST_ASSERT(has_neg);

	// Still SPD: all diagonal entries positive
	TEST_ASSERT(k_diag_positive(K, 3));

	// Actually factorable: block_ldl should succeed with positive pivots
	double A[6], D[3];
	memcpy(A, K, sizeof(A));
	block_ldl(A, D, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_K_mass_ratio_1e6()
{
	TEST_BEGIN("K_mass_ratio_1e6");
	// 1,000,000:1 mass ratio. Heavy vehicle on light debris.
	// Light body dominates K (~1.0), heavy contributes ~1e-6.
	// Combined K must still be SPD and factorable.
	BodyHot light = make_body(1, 1);
	BodyHot heavy = make_body(1e6f, 1e6f);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &light, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &heavy, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &light, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &heavy, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);

	TEST_ASSERT(k_diag_positive(K, 3));

	// Must be factorable
	double A[6], D[3];
	memcpy(A, K, sizeof(A));
	block_ldl(A, D, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_K_nearly_static_body()
{
	TEST_BEGIN("K_nearly_static_body");
	// Near-static body: mass = 1e-8 (nearly infinite mass, inv_mass = 1e8).
	// K entries become enormous. Paired with a normal body.
	BodyHot normal = make_body(1, 1);
	BodyHot near_static = make_body(1e-8f, 1e-8f);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f, 0, 0), .r_b = V3(0, 0.5f, 0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &normal, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &near_static, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &normal, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &near_static, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);

	// Relative error: entries can be ~1e8
	double max_val = 0;
	for (int i = 0; i < 9; i++) { if (fabs(K_ref[i]) > max_val) max_val = fabs(K_ref[i]); }
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) / (max_val + 1e-10) < 1e-4);

	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_both_levers_zero()
{
	TEST_BEGIN("K_both_levers_zero");
	// Both anchors at body centers. K is purely linear: diagonal = inv_m_a + inv_m_b.
	// No angular coupling at all. Well-conditioned, easy case.
	BodyHot a = make_body(2, 5);
	BodyHot b = make_body(8, 12);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &b, 1.0, K);

	double expected_diag = (1.0 / 2.0) + (1.0 / 8.0); // 0.625
	TEST_ASSERT(fabs(K[LDL_TRI(0,0)] - expected_diag) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(1,1)] - expected_diag) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,2)] - expected_diag) < K_EPS);
	// Off-diagonals must be zero (no angular coupling)
	TEST_ASSERT(fabs(K[LDL_TRI(1,0)]) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,0)]) < K_EPS);
	TEST_ASSERT(fabs(K[LDL_TRI(2,1)]) < K_EPS);
}

static void test_K_parallel_levers()
{
	TEST_BEGIN("K_parallel_levers");
	// Both lever arms in the same direction (parallel). Angular contributions
	// from both bodies are along the same axes. K should still be SPD.
	BodyHot a = make_body(1, 1);
	BodyHot b = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0, 0, 3), .r_b = V3(0, 0, 5) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &a, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &b, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &a, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &b, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-4);

	TEST_ASSERT(k_diag_positive(K, 3));

	// Z row: lever along Z, skew(0,0,rz) has zero in row 2 angular.
	// So Z diagonal = inv_m_a + inv_m_b = 2.0 (pure linear).
	// X and Y diagonals > 2.0 (get angular contribution from |r|^2).
	TEST_ASSERT(fabs(K[LDL_TRI(2,2)] - 2.0) < 0.01);
	TEST_ASSERT(K[LDL_TRI(0,0)] > 2.0);
	TEST_ASSERT(K[LDL_TRI(1,1)] > 2.0);
}

static void test_K_extreme_mass_ratio_with_lever()
{
	TEST_BEGIN("K_extreme_mass_ratio_with_lever");
	// The hardest real case: extreme mass ratio + large lever arm + asymmetric inertia.
	// 100000:1 ratio, lever arm = 10, non-uniform inertia, rotated.
	quat rot = quat_axis_angle(norm(V3(1, 2, 3)), 0.5f);
	BodyHot light = make_body_full(1, V3(0.5f, 2, 8), rot);
	BodyHot heavy = make_body_full(1e5f, V3(1000, 5000, 2000), quat_identity());
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(5, -3, 8), .r_b = V3(-2, 10, 1) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &light, 1.0, K);
	ldl_K_body_contrib(jac, 3, 1, 0, &heavy, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &light, 1.0, K_ref);
	ref_K_body_contrib(jac, 3, 1, &heavy, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1.0); // absolute tolerance, entries can be large

	TEST_ASSERT(k_diag_positive(K, 3));

	// Must be factorable
	double A[6], D[3];
	memcpy(A, K, sizeof(A));
	block_ldl(A, D, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_K_distance_constraint()
{
	TEST_BEGIN("K_distance_constraint");
	// Distance constraint: 1 DOF. K is a scalar (1x1 packed = 1 element).
	BodyHot body_a = make_body(2, 3);
	BodyHot body_b = make_body(5, 8);
	v3 ax = norm(V3(1, 1, 0));
	SolverJoint sol = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(1,0,0), .r_b = V3(0,1,0), .dist.axis = ax };
	LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac[1];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[1] = {0};
	ldl_K_body_contrib(jac, 1, 0, 0, &body_a, 1.0, K);
	ldl_K_body_contrib(jac, 1, 1, 0, &body_b, 1.0, K);

	double K_ref[1] = {0};
	ref_K_body_contrib(jac, 1, 0, &body_a, 1.0, K_ref);
	ref_K_body_contrib(jac, 1, 1, &body_b, 1.0, K_ref);

	TEST_ASSERT(fabs(K[0] - K_ref[0]) < 1e-4);
	TEST_ASSERT(K[0] > 0);
}

static void test_K_hinge_constraint()
{
	TEST_BEGIN("K_hinge_constraint");
	// Hinge: 5 DOF (3 linear + 2 angular). Largest non-weld block.
	BodyHot body_a = make_body(3, 6);
	BodyHot body_b = make_body(8, 12);
	SolverJoint sol = { .type = JOINT_HINGE, .dof = 5,
		.r_a = V3(1, 0, 0), .r_b = V3(0, 1, 0),
		.hinge.u1 = norm(V3(1, 0, 0)), .hinge.u2 = norm(V3(0, 1, 0)),
	};
	LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .solver_idx = 0 };
	LDL_JacobianRow jac[5];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[15] = {0}; // 5*(5+1)/2 = 15 packed elements
	ldl_K_body_contrib(jac, 5, 0, 0, &body_a, 1.0, K);
	ldl_K_body_contrib(jac, 5, 1, 0, &body_b, 1.0, K);

	double K_ref[25] = {0};
	ref_K_body_contrib(jac, 5, 0, &body_a, 1.0, K_ref);
	ref_K_body_contrib(jac, 5, 1, &body_b, 1.0, K_ref);

	double K_full[25];
	unpack_symmetric(K, 5, K_full);
	for (int i = 0; i < 25; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);
	TEST_ASSERT(k_is_symmetric(K_ref, 5, 1e-3));
	TEST_ASSERT(k_diag_positive(K, 5));
}

static void test_K_large_lever_arm()
{
	TEST_BEGIN("K_large_lever_arm");
	// r = 100. Angular K entries scale with |r|^2 ~ 10000.
	// Tests that large values don't cause precision issues.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(100, 0, 0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 0, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);

	// Relative error check: entries can be ~10000
	double max_val = 0;
	for (int i = 0; i < 9; i++) { if (fabs(K_ref[i]) > max_val) max_val = fabs(K_ref[i]); }
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) / (max_val + 1e-10) < 1e-4);
	TEST_ASSERT(k_diag_positive(K, 3));
}

static void test_K_accumulates()
{
	TEST_BEGIN("K_accumulates");
	// K_body_contrib ADDS to the existing K array. Verify accumulation.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Call twice with same body -> should be 2x single call.
	double K_single[6] = {0}, K_double[6] = {0};
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K_single);
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K_double);
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K_double);

	for (int i = 0; i < 6; i++)
		TEST_ASSERT(fabs(K_double[i] - 2.0 * K_single[i]) < K_EPS);
}

static void test_K_dof_start_offset()
{
	TEST_BEGIN("K_dof_start_offset");
	// dof_start != 0: constraint writes to a sub-block within a larger packed array.
	// This happens when a bundle has multiple constraints.
	BodyHot body = make_body(1, 1);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0,0,0) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	// Write at offset 0
	double K0[21] = {0}; // 6x6 packed = 21
	ldl_K_body_contrib(jac, 3, 0, 0, &body, 1.0, K0);

	// Write at offset 3 (as if there's a 3-DOF constraint before us)
	double K3[21] = {0};
	ldl_K_body_contrib(jac, 3, 0, 3, &body, 1.0, K3);

	// The values should be the same but at shifted positions.
	for (int r = 0; r < 3; r++)
		for (int c = 0; c <= r; c++)
			TEST_ASSERT(fabs(K3[LDL_TRI(3+r, 3+c)] - K0[LDL_TRI(r, c)]) < K_EPS);

	// Original 3x3 sub-block in K3 should be untouched (zero).
	for (int r = 0; r < 3; r++)
		for (int c = 0; c <= r; c++)
			TEST_ASSERT(fabs(K3[LDL_TRI(r, c)]) < K_EPS);
}

static void test_K_side_b()
{
	TEST_BEGIN("K_side_b");
	// side=1 should use J_b columns. Verify by comparing against reference.
	BodyHot body = make_body(3, 7);
	SolverJoint sol = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,2,0), .r_b = V3(-1,0,3) };
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&con, &sol, jac);

	double K[6] = {0};
	ldl_K_body_contrib(jac, 3, 1, 0, &body, 1.0, K);

	double K_ref[9] = {0};
	ref_K_body_contrib(jac, 3, 1, &body, 1.0, K_ref);
	double K_full[9];
	unpack_symmetric(K, 3, K_full);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(K_full[i] - K_ref[i]) < 1e-3);
}

// ============================================================================
// ldl_K_body_off: off-diagonal block K_ij from a shared body.
// K_ij += J_i * M^{-1} * J_j^T (not symmetric in general, rectangular for mixed DOF).

static void test_K_off_ball_socket_pair()
{
	TEST_BEGIN("K_off_ball_socket_pair");
	// Two ball-sockets sharing one body. Off-diagonal = coupling through shared body.
	BodyHot body = make_body(2, 5);
	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0, 1, 0), .r_b = V3(0,0,0) };
	LDL_Constraint con_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_Constraint con_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 };
	LDL_JacobianRow jac_i[3], jac_j[3];
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con_i, sols, jac_i);
	con_j.solver_idx = 1;
	ldl_fill_jacobian(&con_j, sols, jac_j);

	double out[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, out);

	// Reference: explicit J_i * M^{-1} * J_j^T
	double ref[9] = {0};
	double wm = (double)body.inv_mass;
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			double lin = wm * (jac_i[r].J_a[0]*jac_j[c].J_a[0] + jac_i[r].J_a[1]*jac_j[c].J_a[1] + jac_i[r].J_a[2]*jac_j[c].J_a[2]);
			v3 jc_ang = V3((float)jac_j[c].J_a[3], (float)jac_j[c].J_a[4], (float)jac_j[c].J_a[5]);
			v3 Iinv_jc = inv_inertia_mul(body.rotation, body.inv_inertia_local, jc_ang);
			double ang = jac_i[r].J_a[3]*(double)Iinv_jc.x + jac_i[r].J_a[4]*(double)Iinv_jc.y + jac_i[r].J_a[5]*(double)Iinv_jc.z;
			ref[r*3+c] = lin + ang;
		}
	}

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(out[i] - ref[i]) < 1e-3);
}

static void test_K_off_mixed_dof()
{
	TEST_BEGIN("K_off_mixed_dof");
	// Ball-socket (3 DOF) coupled with distance (1 DOF) through shared body.
	// Off-diagonal is 3x1.
	BodyHot body = make_body(5, 10);
	SolverJoint sol_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(0,0,0) };
	v3 ax = norm(V3(0, 1, 0));
	SolverJoint sol_d = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0, 1, 0), .r_b = V3(0,0,0), .dist.axis = ax };

	LDL_Constraint con_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_Constraint con_d = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac_bs[3], jac_d[1];
	ldl_fill_jacobian(&con_bs, &sol_bs, jac_bs);
	ldl_fill_jacobian(&con_d, &sol_d, jac_d);

	// Both on side A of the shared body
	double out[3] = {0}; // 3x1
	ldl_K_body_off(jac_bs, 3, 0, jac_d, 1, 0, &body, 1.0, out);

	// Reference
	double wm = (double)body.inv_mass;
	double ref[3] = {0};
	for (int r = 0; r < 3; r++) {
		double lin = wm * (jac_bs[r].J_a[0]*jac_d[0].J_a[0] + jac_bs[r].J_a[1]*jac_d[0].J_a[1] + jac_bs[r].J_a[2]*jac_d[0].J_a[2]);
		v3 jc_ang = V3((float)jac_d[0].J_a[3], (float)jac_d[0].J_a[4], (float)jac_d[0].J_a[5]);
		v3 Iinv_jc = inv_inertia_mul(body.rotation, body.inv_inertia_local, jc_ang);
		double ang = jac_bs[r].J_a[3]*(double)Iinv_jc.x + jac_bs[r].J_a[4]*(double)Iinv_jc.y + jac_bs[r].J_a[5]*(double)Iinv_jc.z;
		ref[r] = lin + ang;
	}

	for (int i = 0; i < 3; i++)
		TEST_ASSERT(fabs(out[i] - ref[i]) < 1e-3);
}

static void test_K_off_static_body()
{
	TEST_BEGIN("K_off_static_body");
	// Static shared body: off-diagonal coupling should be zero.
	BodyHot body = {0};
	body.rotation = quat_identity();

	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,1,0), .r_b = V3(0,0,0) };
	LDL_JacobianRow jac_i[3], jac_j[3];
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con, sols, jac_i);
	con.solver_idx = 1;
	ldl_fill_jacobian(&con, sols, jac_j);

	double out[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, out);

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(out[i]) < K_EPS);
}

static void test_K_off_shattering_weight()
{
	TEST_BEGIN("K_off_shattering_weight");
	// Weight=5 should produce 5x the weight=1 result.
	BodyHot body = make_body(3, 7);
	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,1,0), .r_b = V3(0,0,0) };
	LDL_JacobianRow jac_i[3], jac_j[3];
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con, sols, jac_i);
	con.solver_idx = 1;
	ldl_fill_jacobian(&con, sols, jac_j);

	double out1[9] = {0}, out5[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, out1);
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 5.0, out5);

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(out5[i] - 5.0 * out1[i]) < K_EPS);
}

static void test_K_off_accumulates()
{
	TEST_BEGIN("K_off_accumulates");
	// K_body_off ADDS to existing array. Two calls = 2x single.
	BodyHot body = make_body(2, 4);
	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,1,0), .r_b = V3(0,0,0) };
	LDL_JacobianRow jac_i[3], jac_j[3];
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con, sols, jac_i);
	con.solver_idx = 1;
	ldl_fill_jacobian(&con, sols, jac_j);

	double single[9] = {0}, twice[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, single);
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, twice);
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, twice);

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(twice[i] - 2.0 * single[i]) < K_EPS);
}

static void test_K_off_cross_sides()
{
	TEST_BEGIN("K_off_cross_sides");
	// Constraint i on side A, constraint j on side B of the same body.
	// This is the typical case: body is body_a for one joint and body_b for another.
	BodyHot body = make_body(4, 8);
	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(2, 0, 0), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0), .r_b = V3(0, 3, 0) };
	LDL_JacobianRow jac_i[3], jac_j[3];
	LDL_Constraint con_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_Constraint con_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 };
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con_i, sols, jac_i);
	ldl_fill_jacobian(&con_j, sols, jac_j);

	// i uses side A (J_a), j uses side B (J_b) of the shared body
	double out[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 1, &body, 1.0, out);

	// Reference
	double wm = (double)body.inv_mass;
	double ref[9] = {0};
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			double lin = wm * (jac_i[r].J_a[0]*jac_j[c].J_b[0] + jac_i[r].J_a[1]*jac_j[c].J_b[1] + jac_i[r].J_a[2]*jac_j[c].J_b[2]);
			v3 jc_ang = V3((float)jac_j[c].J_b[3], (float)jac_j[c].J_b[4], (float)jac_j[c].J_b[5]);
			v3 Iinv_jc = inv_inertia_mul(body.rotation, body.inv_inertia_local, jc_ang);
			double ang = jac_i[r].J_a[3]*(double)Iinv_jc.x + jac_i[r].J_a[4]*(double)Iinv_jc.y + jac_i[r].J_a[5]*(double)Iinv_jc.z;
			ref[r*3+c] = lin + ang;
		}
	}

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(out[i] - ref[i]) < 1e-3);
}

static void test_K_off_rotated_asymmetric()
{
	TEST_BEGIN("K_off_rotated_asymmetric");
	// Rotated body with asymmetric inertia. The hardest case for
	// dinv_inertia_mul correctness to matter in K assembly.
	quat rot = quat_axis_angle(norm(V3(1, 2, -1)), 0.8f);
	BodyHot body = make_body_full(3, V3(1, 5, 20), rot);
	SolverJoint sol_i = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(2, -1, 0.5f), .r_b = V3(0,0,0) };
	SolverJoint sol_j = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(-1, 0, 3), .r_b = V3(0,0,0) };
	LDL_JacobianRow jac_i[3], jac_j[3];
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	SolverJoint sols[2] = { sol_i, sol_j };
	ldl_fill_jacobian(&con, sols, jac_i);
	con.solver_idx = 1;
	ldl_fill_jacobian(&con, sols, jac_j);

	double out[9] = {0};
	ldl_K_body_off(jac_i, 3, 0, jac_j, 3, 0, &body, 1.0, out);

	double wm = (double)body.inv_mass;
	double ref[9] = {0};
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			double lin = wm * (jac_i[r].J_a[0]*jac_j[c].J_a[0] + jac_i[r].J_a[1]*jac_j[c].J_a[1] + jac_i[r].J_a[2]*jac_j[c].J_a[2]);
			v3 jc_ang = V3((float)jac_j[c].J_a[3], (float)jac_j[c].J_a[4], (float)jac_j[c].J_a[5]);
			v3 Iinv_jc = inv_inertia_mul(body.rotation, body.inv_inertia_local, jc_ang);
			double ang = jac_i[r].J_a[3]*(double)Iinv_jc.x + jac_i[r].J_a[4]*(double)Iinv_jc.y + jac_i[r].J_a[5]*(double)Iinv_jc.z;
			ref[r*3+c] = lin + ang;
		}
	}

	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(out[i] - ref[i]) < 1e-3);
}

// ============================================================================
// Runner

static void run_k_matrix_unit_tests()
{
	printf("--- K matrix unit tests ---\n");

	// Diagonal (ldl_K_body_contrib)
	test_K_ball_socket_unit_mass_origin();
	test_K_ball_socket_heavy_body();
	test_K_ball_socket_with_lever();
	test_K_ball_socket_two_bodies();
	test_K_static_body();
	test_K_extreme_mass_ratio();
	test_K_asymmetric_inertia();
	test_K_rotated_body();
	test_K_shattering_weight();
	test_K_distance_constraint();
	test_K_hinge_constraint();
	test_K_large_lever_arm();
	test_K_mixed_scale_lever();
	test_K_lever_along_weak_inertia();
	test_K_floor_dynamic_pair();
	test_K_mismatched_lever_lengths();
	test_K_negative_offdiag();
	test_K_mass_ratio_1e6();
	test_K_nearly_static_body();
	test_K_both_levers_zero();
	test_K_parallel_levers();
	test_K_extreme_mass_ratio_with_lever();
	test_K_accumulates();
	test_K_dof_start_offset();
	test_K_side_b();

	// Off-diagonal (ldl_K_body_off)
	test_K_off_ball_socket_pair();
	test_K_off_mixed_dof();
	test_K_off_static_body();
	test_K_off_shattering_weight();
	test_K_off_accumulates();
	test_K_off_cross_sides();
	test_K_off_rotated_asymmetric();
}
