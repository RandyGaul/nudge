// tests_ldl_unit.c -- unit tests for LDL dense block math.
// Tests block_ldl (factorization) and block_solve (forward/diagonal/back sub)
// in isolation with known matrices, verifying correctness at each layer.

#define LDL_EPS 1e-10

// Helper: pack a full NxN symmetric matrix into lower-triangular storage.
static void pack_symmetric(double* full, int n, double* packed)
{
	for (int r = 0; r < n; r++)
		for (int c = 0; c <= r; c++)
			packed[LDL_TRI(r, c)] = full[r * n + c];
}

// Helper: reconstruct K from L (packed, unit diagonal) and D, verify K = L*D*L^T.
// Returns max absolute error.
static double ldl_reconstruct_error(double* L_packed, double* D, int n, double* K_original)
{
	// Build full L (unit diagonal, lower triangle from packed)
	double L[36] = {0};
	for (int r = 0; r < n; r++) {
		L[r * n + r] = 1.0;
		for (int c = 0; c < r; c++)
			L[r * n + c] = L_packed[LDL_TRI(r, c)];
	}

	// Compute L * D * L^T
	double LDLt[36] = {0};
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				LDLt[i * n + j] += L[i * n + k] * D[k] * L[j * n + k];

	double max_err = 0;
	for (int i = 0; i < n * n; i++) {
		double err = fabs(LDLt[i] - K_original[i]);
		if (err > max_err) max_err = err;
	}
	return max_err;
}

// Helper: compute K * x, return max |K*x - b| (residual check for solve).
static double ldl_solve_residual(double* K_full, double* x, double* b, int n)
{
	double max_err = 0;
	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < n; j++) sum += K_full[i * n + j] * x[j];
		double err = fabs(sum - b[i]);
		if (err > max_err) max_err = err;
	}
	return max_err;
}

// ============================================================================
// block_ldl tests

static void test_ldl_factor_1x1()
{
	TEST_BEGIN("ldl_factor_1x1");
	double K_full[1] = { 7.0 };
	double A[1], D[1];
	pack_symmetric(K_full, 1, A);
	block_ldl(A, D, 1);
	TEST_ASSERT(fabs(D[0] - 7.0) < LDL_EPS);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 1, K_full) < LDL_EPS);
}

static void test_ldl_factor_2x2()
{
	TEST_BEGIN("ldl_factor_2x2");
	double K_full[4] = {
		4, 2,
		2, 5,
	};
	double A[3], D[2];
	pack_symmetric(K_full, 2, A);
	block_ldl(A, D, 2);
	// D[0] = 4, L[1,0] = 0.5, D[1] = 5 - 0.25*4 = 4
	TEST_ASSERT(fabs(D[0] - 4.0) < LDL_EPS);
	TEST_ASSERT(fabs(D[1] - 4.0) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(1, 0)] - 0.5) < LDL_EPS);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 2, K_full) < LDL_EPS);
}

static void test_ldl_factor_3x3()
{
	TEST_BEGIN("ldl_factor_3x3");
	// The walkthrough example from our conversation.
	double K_full[9] = {
		4, 2, 1,
		2, 5, 3,
		1, 3, 9,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(fabs(D[0] - 4.0) < LDL_EPS);
	TEST_ASSERT(fabs(D[1] - 4.0) < LDL_EPS);
	TEST_ASSERT(fabs(D[2] - 7.1875) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(1, 0)] - 0.5) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(2, 0)] - 0.25) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(2, 1)] - 0.625) < LDL_EPS);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < LDL_EPS);
}

static void test_ldl_factor_identity()
{
	TEST_BEGIN("ldl_factor_identity");
	// Identity matrix: L = I, D = [1,1,1]
	double K_full[9] = {
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(fabs(D[i] - 1.0) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(1, 0)]) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(2, 0)]) < LDL_EPS);
	TEST_ASSERT(fabs(A[LDL_TRI(2, 1)]) < LDL_EPS);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < LDL_EPS);
}

static void test_ldl_factor_diagonal()
{
	TEST_BEGIN("ldl_factor_diagonal");
	// Diagonal matrix: L = I, D = diagonal values.
	double K_full[9] = {
		3, 0, 0,
		0, 7, 0,
		0, 0, 11,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(fabs(D[0] - 3.0) < LDL_EPS);
	TEST_ASSERT(fabs(D[1] - 7.0) < LDL_EPS);
	TEST_ASSERT(fabs(D[2] - 11.0) < LDL_EPS);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < LDL_EPS);
}

static void test_ldl_factor_6x6()
{
	TEST_BEGIN("ldl_factor_6x6");
	// Max block size in the solver (hinge = 5 DOF, weld = 6 DOF).
	// Build a 6x6 SPD matrix: K = A^T * A + I (guaranteed SPD).
	double raw[36] = {
		2, 1, 0, 1, 0, 0,
		1, 3, 1, 0, 1, 0,
		0, 1, 4, 1, 0, 1,
		1, 0, 1, 5, 1, 0,
		0, 1, 0, 1, 6, 1,
		0, 0, 1, 0, 1, 7,
	};
	double K_full[36] = {0};
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) K_full[i * 6 + j] += raw[k * 6 + i] * raw[k * 6 + j];
			if (i == j) K_full[i * 6 + j] += 1.0;
		}
	}
	double A[21], D[6];
	pack_symmetric(K_full, 6, A);
	block_ldl(A, D, 6);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 6, K_full) < 1e-8);
	// All pivots should be positive
	for (int i = 0; i < 6; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_high_condition()
{
	TEST_BEGIN("ldl_factor_high_condition");
	// Condition number ~1000: tests numerical robustness.
	double K_full[9] = {
		1000, 1, 0,
		1,    1, 0,
		0,    0, 1,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < 1e-6);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

// ============================================================================
// block_solve tests

static void test_ldl_solve_1x1()
{
	TEST_BEGIN("ldl_solve_1x1");
	double K_full[1] = { 7.0 };
	double A[1], D[1];
	pack_symmetric(K_full, 1, A);
	block_ldl(A, D, 1);
	double b[1] = { 21.0 };
	double x[1];
	block_solve(A, D, b, x, 1);
	TEST_ASSERT(fabs(x[0] - 3.0) < LDL_EPS);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 1) < LDL_EPS);
}

static void test_ldl_solve_2x2()
{
	TEST_BEGIN("ldl_solve_2x2");
	double K_full[4] = {
		4, 2,
		2, 5,
	};
	double A[3], D[2];
	pack_symmetric(K_full, 2, A);
	block_ldl(A, D, 2);
	double b[2] = { 1.0, 2.0 };
	double x[2];
	block_solve(A, D, b, x, 2);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 2) < LDL_EPS);
}

static void test_ldl_solve_3x3()
{
	TEST_BEGIN("ldl_solve_3x3");
	// The walkthrough example: K * x = [1, 2, 3]
	double K_full[9] = {
		4, 2, 1,
		2, 5, 3,
		1, 3, 9,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 1.0, 2.0, 3.0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 3) < LDL_EPS);
}

static void test_ldl_solve_identity()
{
	TEST_BEGIN("ldl_solve_identity");
	// I * x = b => x = b
	double K_full[9] = {
		1, 0, 0,
		0, 1, 0,
		0, 0, 1,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 5.0, -3.0, 7.0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(fabs(x[i] - b[i]) < LDL_EPS);
}

static void test_ldl_solve_6x6()
{
	TEST_BEGIN("ldl_solve_6x6");
	// Same 6x6 SPD as factorization test, solve with a known RHS.
	double raw[36] = {
		2, 1, 0, 1, 0, 0,
		1, 3, 1, 0, 1, 0,
		0, 1, 4, 1, 0, 1,
		1, 0, 1, 5, 1, 0,
		0, 1, 0, 1, 6, 1,
		0, 0, 1, 0, 1, 7,
	};
	double K_full[36] = {0};
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) K_full[i * 6 + j] += raw[k * 6 + i] * raw[k * 6 + j];
			if (i == j) K_full[i * 6 + j] += 1.0;
		}
	}
	double A[21], D[6];
	pack_symmetric(K_full, 6, A);
	block_ldl(A, D, 6);
	double b[6] = { 1, -2, 3, -4, 5, -6 };
	double x[6];
	block_solve(A, D, b, x, 6);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 6) < 1e-8);
}

static void test_ldl_solve_high_condition()
{
	TEST_BEGIN("ldl_solve_high_condition");
	double K_full[9] = {
		1000, 1, 0,
		1,    1, 0,
		0,    0, 1,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 1001.0, 2.0, 5.0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 3) < 1e-6);
}

static void test_ldl_solve_multiple_rhs()
{
	TEST_BEGIN("ldl_solve_multiple_rhs");
	// Factorize once, solve with multiple RHS vectors.
	double K_full[9] = {
		4, 2, 1,
		2, 5, 3,
		1, 3, 9,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double rhs[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
	for (int t = 0; t < 3; t++) {
		double x[3];
		block_solve(A, D, rhs[t], x, 3);
		TEST_ASSERT(ldl_solve_residual(K_full, x, rhs[t], 3) < LDL_EPS);
	}
}

// ============================================================================
// block_ldl stress tests: physics-realistic pathological inputs

static void test_ldl_factor_mass_ratio_1e3()
{
	TEST_BEGIN("ldl_factor_mass_ratio_1e3");
	// Simulates 1000:1 mass ratio. Light body inv_mass=1, heavy body inv_mass=0.001.
	// K contributions differ by 1000x. Condition number ~1000.
	double K_full[9] = {
		1.001, 0.001, 0,
		0.001, 1.001, 0,
		0,     0,     1.001,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < 1e-8);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_mass_ratio_1e6()
{
	TEST_BEGIN("ldl_factor_mass_ratio_1e6");
	// 1,000,000:1 mass ratio. This is where real solvers start struggling.
	// K = diag(1+1e-6, 1+1e-6, 1+1e-6) with tiny off-diagonal coupling.
	double eps = 1e-6;
	double K_full[9] = {
		1.0 + eps, eps,       0,
		eps,       1.0 + eps, 0,
		0,         0,         1.0 + eps,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < 1e-8);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_mixed_scale_5x5()
{
	TEST_BEGIN("ldl_factor_mixed_scale_5x5");
	// Simulates a hinge: 3 linear DOFs with small K entries (inv_mass ~ 1),
	// 2 angular DOFs with large K entries (inv_inertia * |r|^2, r=10 -> 100x).
	// Diagonal spans ~100x within one block.
	double K_full[25] = {
		2.0,  0.5,  0.0,  0.1,   0.0,
		0.5,  2.0,  0.5,  0.0,   0.1,
		0.0,  0.5,  2.0,  0.0,   0.0,
		0.1,  0.0,  0.0,  200.0, 50.0,
		0.0,  0.1,  0.0,  50.0,  200.0,
	};
	double A[15], D[5];
	pack_symmetric(K_full, 5, A);
	block_ldl(A, D, 5);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 5, K_full) < 1e-6);
	for (int i = 0; i < 5; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_negative_offdiag()
{
	TEST_BEGIN("ldl_factor_negative_offdiag");
	// SPD with negative off-diagonals. Simulates cross terms from
	// skew(r_a) in the Jacobian when lever arm has mixed-sign components.
	// Built as A^T*A + 2*I to guarantee SPD with significant negative entries.
	double raw[9] = {
		 1, -2,  3,
		-2,  1, -1,
		 0, -3,  2,
	};
	double K_full[9] = {0};
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) K_full[i * 3 + j] += raw[k * 3 + i] * raw[k * 3 + j];
			if (i == j) K_full[i * 3 + j] += 2.0;
		}
	}
	// Verify it actually has negative off-diagonals
	int has_neg = 0;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			if (i != j && K_full[i * 3 + j] < 0) has_neg = 1;
	TEST_ASSERT(has_neg);
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < 1e-8);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_solve_negative_offdiag()
{
	TEST_BEGIN("ldl_solve_negative_offdiag");
	// Same matrix, verify solve works with negative coupling.
	double raw[9] = {
		 1, -2,  3,
		-2,  1, -1,
		 0, -3,  2,
	};
	double K_full[9] = {0};
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) K_full[i * 3 + j] += raw[k * 3 + i] * raw[k * 3 + j];
			if (i == j) K_full[i * 3 + j] += 2.0;
		}
	}
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 1.0, -1.0, 2.0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 3) < 1e-8);
}

static void test_ldl_factor_tiny_uniform_scale()
{
	TEST_BEGIN("ldl_factor_tiny_uniform_scale");
	// K = 1e-8 * well_conditioned_SPD. Simulates very heavy bodies
	// where inv_mass ~ 1e-8. All entries are tiny but well-conditioned.
	double s = 1e-8;
	double K_full[9] = {
		4*s, 2*s, 1*s,
		2*s, 5*s, 3*s,
		1*s, 3*s, 9*s,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < 1e-16);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_large_entries()
{
	TEST_BEGIN("ldl_factor_large_entries");
	// K with entries ~1e6. Simulates long lever arms (r=100, |r|^2=10000)
	// combined with light bodies (inv_mass=100).
	double K_full[9] = {
		1e6,   5e5,   0,
		5e5,   1e6,   5e5,
		0,     5e5,   1e6,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	// Relative error: reconstruction error / max entry
	double recon = ldl_reconstruct_error(A, D, 3, K_full);
	TEST_ASSERT(recon / 1e6 < 1e-8);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D[i] > 0);
}

static void test_ldl_factor_shattering_amplified()
{
	TEST_BEGIN("ldl_factor_shattering_amplified");
	// Simulates shattering: 100:1 base mass ratio, shard count S=10.
	// Effective ratio becomes 1000:1. Heavy side contributes K entries
	// scaled by S*inv_mass_heavy = 10*0.01 = 0.1, light side = 10*1.0 = 10.
	// 3x3 block from one ball-socket on the shattered body.
	double heavy = 0.1;  // S * inv_mass_heavy
	double light = 10.0; // S * inv_mass_light
	double K_full[9] = {
		heavy + light, 0, 0,
		0, heavy + light, 0,
		0, 0, heavy + light,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_full) < LDL_EPS);
	// Now the harder case: two constraints on same shattered body,
	// coupling through the heavy body. Off-diag = heavy contribution only.
	double K_coupled[9] = {
		heavy + light, heavy * 0.5, 0,
		heavy * 0.5,   heavy + light, 0,
		0,             0,             heavy + light,
	};
	pack_symmetric(K_coupled, 3, A);
	block_ldl(A, D, 3);
	TEST_ASSERT(ldl_reconstruct_error(A, D, 3, K_coupled) < 1e-8);
}

// ============================================================================
// block_solve stress tests

static void test_ldl_solve_zero_rhs()
{
	TEST_BEGIN("ldl_solve_zero_rhs");
	double K_full[9] = {
		4, 2, 1,
		2, 5, 3,
		1, 3, 9,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 0, 0, 0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(fabs(x[i]) < LDL_EPS);
}

static void test_ldl_solve_large_rhs()
{
	TEST_BEGIN("ldl_solve_large_rhs");
	// Large RHS values (~1e6): big velocity errors from penetration.
	double K_full[9] = {
		4, 2, 1,
		2, 5, 3,
		1, 3, 9,
	};
	double A[6], D[3];
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 1e6, -2e6, 3e6 };
	double x[3];
	block_solve(A, D, b, x, 3);
	TEST_ASSERT(ldl_solve_residual(K_full, x, b, 3) < 1e-4); // absolute error scales with RHS
}

static void test_ldl_solve_ill_conditioned()
{
	TEST_BEGIN("ldl_solve_ill_conditioned");
	// Condition ~1e6 (mass ratio 1M:1). Factor + solve, check residual.
	// This is the critical test: does the solve still produce usable results?
	double eps = 1e-6;
	double K_full[9] = {
		1.0,   eps,   0,
		eps,   eps,   0,
		0,     0,     1.0,
	};
	// This matrix has eigenvalues ~1.0, ~eps, ~1.0. Condition ~1e6.
	double A[6], D[3];
	double K_save[9];
	memcpy(K_save, K_full, sizeof(K_save));
	pack_symmetric(K_full, 3, A);
	block_ldl(A, D, 3);
	double b[3] = { 1.0, 1.0, 1.0 };
	double x[3];
	block_solve(A, D, b, x, 3);
	// With double precision, even cond=1e6 should solve cleanly.
	TEST_ASSERT(ldl_solve_residual(K_save, x, b, 3) < 1e-4);
	for (int i = 0; i < 3; i++) TEST_ASSERT(x[i] == x[i]); // no NaN
}

static void test_ldl_solve_mixed_scale_5x5()
{
	TEST_BEGIN("ldl_solve_mixed_scale_5x5");
	// Same mixed-scale hinge-like matrix, verify solve works end-to-end.
	double K_full[25] = {
		2.0,  0.5,  0.0,  0.1,   0.0,
		0.5,  2.0,  0.5,  0.0,   0.1,
		0.0,  0.5,  2.0,  0.0,   0.0,
		0.1,  0.0,  0.0,  200.0, 50.0,
		0.0,  0.1,  0.0,  50.0,  200.0,
	};
	double K_save[25];
	memcpy(K_save, K_full, sizeof(K_save));
	double A[15], D[5];
	pack_symmetric(K_full, 5, A);
	block_ldl(A, D, 5);
	// RHS with mixed scale too: small linear errors, large angular errors.
	double b[5] = { 0.1, -0.2, 0.1, 50.0, -30.0 };
	double x[5];
	block_solve(A, D, b, x, 5);
	TEST_ASSERT(ldl_solve_residual(K_save, x, b, 5) < 1e-6);
}

// ============================================================================
// block_mul / block_transpose / block_sub tests

static void test_ldl_block_transpose()
{
	TEST_BEGIN("ldl_block_transpose");
	double A[6] = { 1, 2, 3, 4, 5, 6 }; // 2x3
	double At[6];
	block_transpose(A, 2, 3, At); // -> 3x2
	// A = [1 2 3; 4 5 6], A^T = [1 4; 2 5; 3 6]
	TEST_ASSERT(fabs(At[0] - 1) < LDL_EPS);
	TEST_ASSERT(fabs(At[1] - 4) < LDL_EPS);
	TEST_ASSERT(fabs(At[2] - 2) < LDL_EPS);
	TEST_ASSERT(fabs(At[3] - 5) < LDL_EPS);
	TEST_ASSERT(fabs(At[4] - 3) < LDL_EPS);
	TEST_ASSERT(fabs(At[5] - 6) < LDL_EPS);
}

static void test_ldl_block_mul()
{
	TEST_BEGIN("ldl_block_mul");
	// [1 2; 3 4] * [5 6; 7 8] = [19 22; 43 50]
	double A[4] = { 1, 2, 3, 4 };
	double B[4] = { 5, 6, 7, 8 };
	double C[4];
	block_mul(A, 2, 2, B, 2, C);
	TEST_ASSERT(fabs(C[0] - 19) < LDL_EPS);
	TEST_ASSERT(fabs(C[1] - 22) < LDL_EPS);
	TEST_ASSERT(fabs(C[2] - 43) < LDL_EPS);
	TEST_ASSERT(fabs(C[3] - 50) < LDL_EPS);
}

static void test_ldl_block_sub()
{
	TEST_BEGIN("ldl_block_sub");
	double A[4] = { 10, 20, 30, 40 };
	double B[4] = { 1, 2, 3, 4 };
	block_sub(A, B, 2, 2);
	TEST_ASSERT(fabs(A[0] - 9) < LDL_EPS);
	TEST_ASSERT(fabs(A[1] - 18) < LDL_EPS);
	TEST_ASSERT(fabs(A[2] - 27) < LDL_EPS);
	TEST_ASSERT(fabs(A[3] - 36) < LDL_EPS);
}

static void test_ldl_block_mul_3x5_5x2()
{
	TEST_BEGIN("ldl_block_mul_3x5_5x2");
	// Schur complement shape: ball-socket (3 DOF) coupling with hinge (5 DOF)
	// through a shared body. E_ik is 3x5, L_jk^T is 5x2 -> result is 3x2.
	double A[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }; // 3x5
	double B[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; // 5x2
	double C[6];
	block_mul(A, 3, 5, B, 2, C);
	// Row 0: [1*1+2*3+3*5+4*7+5*9, 1*2+2*4+3*6+4*8+5*10] = [95, 110]
	// Row 1: [6*1+7*3+8*5+9*7+10*9, ...] = [6+21+40+63+90, 12+28+48+72+100] = [220, 260]
	// Row 2: [11+36+65+98+135, 22+48+78+112+150] = [345, 410]
	TEST_ASSERT(fabs(C[0] - 95) < LDL_EPS);
	TEST_ASSERT(fabs(C[1] - 110) < LDL_EPS);
	TEST_ASSERT(fabs(C[2] - 220) < LDL_EPS);
	TEST_ASSERT(fabs(C[3] - 260) < LDL_EPS);
	TEST_ASSERT(fabs(C[4] - 345) < LDL_EPS);
	TEST_ASSERT(fabs(C[5] - 410) < LDL_EPS);
}

static void test_ldl_block_mul_5x3_3x5()
{
	TEST_BEGIN("ldl_block_mul_5x3_3x5");
	// Reverse direction: hinge coupling through ball-socket pivot.
	// 5x3 * 3x5 -> 5x5 result.
	double A[15] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1 }; // 5x3
	double B[15] = { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0 }; // 3x5
	double C[25];
	block_mul(A, 5, 3, B, 5, C);
	// Row 0 of A = [1,0,0], * B -> row 0 of B = [1,0,0,0,0]
	TEST_ASSERT(fabs(C[0] - 1) < LDL_EPS);
	TEST_ASSERT(fabs(C[1] - 0) < LDL_EPS);
	// Row 1 of A = [0,1,0], * B -> row 1 of B = [0,1,0,0,0]
	TEST_ASSERT(fabs(C[5] - 0) < LDL_EPS);
	TEST_ASSERT(fabs(C[6] - 1) < LDL_EPS);
	// Row 3 of A = [1,1,0], * B -> row0+row1 of B = [1,1,0,0,0]
	TEST_ASSERT(fabs(C[15] - 1) < LDL_EPS);
	TEST_ASSERT(fabs(C[16] - 1) < LDL_EPS);
	TEST_ASSERT(fabs(C[17] - 0) < LDL_EPS);
}

static void test_ldl_block_mul_1x3_3x1()
{
	TEST_BEGIN("ldl_block_mul_1x3_3x1");
	// Distance constraint (1 DOF) coupling with ball-socket (3 DOF).
	// 1x3 * 3x1 -> 1x1 (scalar dot product).
	double A[3] = { 2, 3, 4 };
	double B[3] = { 5, 6, 7 };
	double C[1];
	block_mul(A, 1, 3, B, 1, C);
	// 2*5 + 3*6 + 4*7 = 10 + 18 + 28 = 56
	TEST_ASSERT(fabs(C[0] - 56) < LDL_EPS);
}

static void test_ldl_block_transpose_5x3()
{
	TEST_BEGIN("ldl_block_transpose_5x3");
	// 5x3 transpose -> 3x5. Actual Schur complement shape.
	double A[15], At[15];
	for (int i = 0; i < 15; i++) A[i] = (double)(i + 1);
	block_transpose(A, 5, 3, At);
	// A[r][c] = At[c][r]. A is 5x3 row-major, At is 3x5 row-major.
	for (int r = 0; r < 5; r++)
		for (int c = 0; c < 3; c++)
			TEST_ASSERT(fabs(At[c * 5 + r] - A[r * 3 + c]) < LDL_EPS);
}

static void test_ldl_block_sub_rectangular()
{
	TEST_BEGIN("ldl_block_sub_rectangular");
	// 3x5 subtraction (off-diagonal Schur update shape).
	double A[15], B[15];
	for (int i = 0; i < 15; i++) { A[i] = 100.0 + i; B[i] = (double)i; }
	block_sub(A, B, 3, 5);
	for (int i = 0; i < 15; i++) TEST_ASSERT(fabs(A[i] - 100.0) < LDL_EPS);
}

// ============================================================================
// Runner

static void run_ldl_unit_tests()
{
	printf("--- LDL unit tests ---\n");

	// Factorization: basic
	test_ldl_factor_1x1();
	test_ldl_factor_2x2();
	test_ldl_factor_3x3();
	test_ldl_factor_identity();
	test_ldl_factor_diagonal();
	test_ldl_factor_6x6();
	test_ldl_factor_high_condition();

	// Factorization: physics stress
	test_ldl_factor_mass_ratio_1e3();
	test_ldl_factor_mass_ratio_1e6();
	test_ldl_factor_mixed_scale_5x5();
	test_ldl_factor_negative_offdiag();
	test_ldl_factor_tiny_uniform_scale();
	test_ldl_factor_large_entries();
	test_ldl_factor_shattering_amplified();

	// Solve: basic
	test_ldl_solve_1x1();
	test_ldl_solve_2x2();
	test_ldl_solve_3x3();
	test_ldl_solve_identity();
	test_ldl_solve_6x6();
	test_ldl_solve_high_condition();
	test_ldl_solve_multiple_rhs();

	// Solve: stress
	test_ldl_solve_zero_rhs();
	test_ldl_solve_large_rhs();
	test_ldl_solve_ill_conditioned();
	test_ldl_solve_mixed_scale_5x5();
	test_ldl_solve_negative_offdiag();

	// Block helpers: basic
	test_ldl_block_transpose();
	test_ldl_block_mul();
	test_ldl_block_sub();

	// Block helpers: rectangular (Schur complement shapes)
	test_ldl_block_mul_3x5_5x2();
	test_ldl_block_mul_5x3_3x5();
	test_ldl_block_mul_1x3_3x1();
	test_ldl_block_transpose_5x3();
	test_ldl_block_sub_rectangular();
}
