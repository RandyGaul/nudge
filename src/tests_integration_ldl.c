// tests_integration_ldl.c -- integration tests for the full sparse LDL pipeline.
// Compares sparse factorize + solve against a dense reference oracle.
// Every test: build physics -> sparse pipeline -> dense reference -> compare solutions.

// ============================================================================
// Dense reference oracle: scalar LDL factorize + solve on full n x n matrix.
// No blocks, no sparsity, no equilibration. Pure double precision.

static void dense_ldl_solve(double* K, double* rhs, double* x, int n)
{
	// Copy K for in-place factorization
	double* A = CK_ALLOC(n * n * sizeof(double));
	double* D = CK_ALLOC(n * sizeof(double));
	memcpy(A, K, n * n * sizeof(double));

	// LDL^T factorization (scalar, dense)
	for (int j = 0; j < n; j++) {
		double dj = A[j * n + j];
		for (int k = 0; k < j; k++) dj -= A[j * n + k] * A[j * n + k] * D[k];
		assert(dj > 0 && "dense_ldl_solve: non-positive pivot");
		D[j] = dj;
		double inv_dj = 1.0 / dj;
		for (int i = j + 1; i < n; i++) {
			double lij = A[i * n + j];
			for (int k = 0; k < j; k++) lij -= A[i * n + k] * A[j * n + k] * D[k];
			A[i * n + j] = lij * inv_dj;
		}
	}

	// Forward substitution
	for (int i = 0; i < n; i++) { x[i] = rhs[i]; for (int k = 0; k < i; k++) x[i] -= A[i * n + k] * x[k]; }
	// Diagonal solve
	for (int i = 0; i < n; i++) x[i] /= D[i];
	// Back substitution
	for (int i = n - 1; i >= 0; i--) for (int k = i + 1; k < n; k++) x[i] -= A[k * n + i] * x[k];

	CK_FREE(A);
	CK_FREE(D);
}

// ============================================================================
// Helper: run the full sparse pipeline on an LDL_Cache + WorldInternal.
// Returns the solution in x[]. Caller provides rhs[].

typedef struct IntegrationWorld
{
	BodyHot bodies[16];
	int body_count;
	SolverJoint sol_joints[16];
	int joint_count;
} IntegrationWorld;

static void integration_sparse_solve(LDL_Cache* c, IntegrationWorld* iw, double* rhs, double* x)
{
	// Build a minimal WorldInternal with just body_hot
	WorldInternal w;
	memset(&w, 0, sizeof(w));
	// Use ckit dynamic array for body_hot
	afit(w.body_hot, iw->body_count);
	asetlen(w.body_hot, iw->body_count);
	for (int i = 0; i < iw->body_count; i++) w.body_hot[i] = iw->bodies[i];

	ldl_build_bundles(c);
	ldl_build_topology(c, &w);
	ldl_numeric_factor(c, &w, iw->sol_joints);
	ldl_solve_topo(c->topo, c->diag_data, c->diag_D, c->L_factors, rhs, x);

	afree(w.body_hot);
}

// Helper: build dense K independently using Jacobians from the cache (after numeric_factor filled them).
static void integration_dense_solve(LDL_Cache* c, IntegrationWorld* iw, double* rhs, double* x)
{
	int n = c->n;
	int jc = c->joint_count;
	double* K = CK_ALLOC(n * n * sizeof(double));
	memset(K, 0, n * n * sizeof(double));

	// Use the Jacobians that numeric_factor computed (stored in c->jacobians)
	// but build K from scratch using tested functions.
	LDL_Topology* t = c->topo;

	for (int i = 0; i < jc; i++) {
		LDL_Constraint* ci = &c->constraints[i];
		if (ci->is_synthetic) continue;
		int oi = t->row_offset[ci->bundle_idx] + ci->bundle_offset;
		int di = ci->dof;
		LDL_JacobianRow* jac = &c->jacobians[ci->jacobian_start];
		BodyHot* a = &iw->bodies[ci->real_body_a];
		BodyHot* b = &iw->bodies[ci->real_body_b];

		// Diagonal: both bodies contribute
		double K_small[21] = {0};
		ldl_K_body_contrib(jac, di, 0, 0, a, ci->weight_a, K_small);
		ldl_K_body_contrib(jac, di, 1, 0, b, ci->weight_b, K_small);
		for (int r = 0; r < di; r++)
			for (int c2 = 0; c2 < di; c2++)
				K[(oi + r) * n + (oi + c2)] += K_small[LDL_TRI(r, c2)];

		// Regularization
		double trace = 0;
		for (int d = 0; d < di; d++) trace += K[(oi + d) * n + (oi + d)];
		double soft = 0;
		soft = iw->sol_joints[ci->solver_idx].softness;
		double reg = (soft > 0) ? soft : LDL_MIN_COMPLIANCE * trace / di;
		for (int d = 0; d < di; d++) K[(oi + d) * n + (oi + d)] += reg;
	}

	// Intra-bundle coupling
	for (int bi = 0; bi < c->bundle_count; bi++) {
		LDL_Bundle* bun = &c->bundles[bi];
		for (int ci = 0; ci < bun->count; ci++) {
			LDL_Constraint* con_i = &c->constraints[bun->start + ci];
			for (int cj = ci + 1; cj < bun->count; cj++) {
				LDL_Constraint* con_j = &c->constraints[bun->start + cj];
				int oi = t->row_offset[con_i->bundle_idx] + con_i->bundle_offset;
				int oj = t->row_offset[con_j->bundle_idx] + con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				LDL_JacobianRow* jac_i = &c->jacobians[con_i->jacobian_start];
				LDL_JacobianRow* jac_j = &c->jacobians[con_j->jacobian_start];
				for (int bs = 0; bs < 2; bs++) {
					int body = bs == 0 ? bun->body_a : bun->body_b;
					int real_b = bs == 0 ? con_i->real_body_a : con_i->real_body_b;
					double wt = (body == con_i->body_a) ? con_i->weight_a : con_i->weight_b;
					int side_i = (body == con_i->body_b) ? 1 : 0;
					int side_j = (body == con_j->body_b) ? 1 : 0;
					double buf[36] = {0};
					ldl_K_body_off(jac_i, di, side_i, jac_j, dj, side_j, &iw->bodies[real_b], wt, buf);
					for (int r = 0; r < di; r++)
						for (int c2 = 0; c2 < dj; c2++) {
							K[(oi + r) * n + (oj + c2)] += buf[r * dj + c2];
							K[(oj + c2) * n + (oi + r)] += buf[r * dj + c2];
						}
				}
			}
		}
	}

	// Inter-bundle coupling (shared body between different bundles)
	int fc = asize(t->couplings);
	for (int fi = 0; fi < fc; fi++) {
		LDL_Coupling* f = &t->couplings[fi];
		LDL_Bundle* bi_b = &c->bundles[f->block_i];
		LDL_Bundle* bj_b = &c->bundles[f->block_j];
		for (int ci = 0; ci < bi_b->count; ci++) {
			LDL_Constraint* con_i = &c->constraints[bi_b->start + ci];
			if (con_i->body_a != f->body && con_i->body_b != f->body) continue;
			int side_i = (f->body == con_i->body_b) ? 1 : 0;
			int real_b = side_i ? con_i->real_body_b : con_i->real_body_a;
			double wt = side_i ? con_i->weight_b : con_i->weight_a;
			LDL_JacobianRow* jac_i = &c->jacobians[con_i->jacobian_start];
			int oi = t->row_offset[con_i->bundle_idx] + con_i->bundle_offset;
			for (int cj = 0; cj < bj_b->count; cj++) {
				LDL_Constraint* con_j = &c->constraints[bj_b->start + cj];
				if (con_j->body_a != f->body && con_j->body_b != f->body) continue;
				int side_j = (f->body == con_j->body_b) ? 1 : 0;
				int oj = t->row_offset[con_j->bundle_idx] + con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				LDL_JacobianRow* jac_j = &c->jacobians[con_j->jacobian_start];
				double buf[36] = {0};
				ldl_K_body_off(jac_i, di, side_i, jac_j, dj, side_j, &iw->bodies[real_b], wt, buf);
				for (int r = 0; r < di; r++)
					for (int c2 = 0; c2 < dj; c2++) {
						K[(oi + r) * n + (oj + c2)] += buf[r * dj + c2];
						K[(oj + c2) * n + (oi + r)] += buf[r * dj + c2];
					}
			}
		}
	}

	dense_ldl_solve(K, rhs, x, n);
	CK_FREE(K);
}

// Helper: free cache resources
static void integration_cache_free(LDL_Cache* c)
{
	afree(c->constraints);
	afree(c->bundles);
	if (c->topo) { ldl_topology_free(c->topo); CK_FREE(c->topo); }
	afree(c->L_factors);
	afree(c->jacobians);
	afree(c->body_remap);
	afree(c->shard_counts);
}

// Helper: compare sparse vs dense solutions
static double max_solution_error(double* x_sparse, double* x_dense, int n)
{
	double max_err = 0;
	for (int i = 0; i < n; i++) {
		double err = fabs(x_sparse[i] - x_dense[i]);
		if (err > max_err) max_err = err;
	}
	return max_err;
}

// ============================================================================
// Tier 1: Structural correctness

static void test_integ_single_ball_socket()
{
	TEST_BEGIN("integ_single_ball_socket");
	IntegrationWorld iw = {0};
	iw.body_count = 2;
	iw.bodies[0] = make_body(2, 4);
	iw.bodies[1] = make_body(5, 10);
	// joint_count set below (was bs: 1
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(0,1,0), .body_a = 0, .body_b = 1 };

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	double rhs[3] = { 1, -0.5, 0.2 };
	double x_sparse[3], x_dense[3];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 3);
	TEST_ASSERT(err < 0.01);
	for (int i = 0; i < 3; i++) TEST_ASSERT(x_sparse[i] == x_sparse[i]); // no NaN

	integration_cache_free(&c);
}

static void test_integ_chain_2()
{
	TEST_BEGIN("integ_chain_2");
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	iw.bodies[0] = make_body(2, 4);
	iw.bodies[1] = make_body(5, 10);
	iw.bodies[2] = make_body(3, 6);
	// joint_count set below (was bs: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f,0,0), .r_b = V3(-0.5f,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f,0,0), .r_b = V3(-0.5f,0,0), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[6] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1 };
	double x_sparse[6], x_dense[6];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 6);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_chain_3()
{
	TEST_BEGIN("integ_chain_3");
	IntegrationWorld iw = {0};
	iw.body_count = 4;
	for (int i = 0; i < 4; i++) iw.bodies[i] = make_body((float)(2 + i), (float)(4 + i));
	// joint_count set below (was bs: 3
	for (int i = 0; i < 3; i++)
		iw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f,0,0), .r_b = V3(-0.5f,0,0), .body_a = i, .body_b = i + 1 };

	LDL_Cache c = {0};
	for (int i = 0; i < 3; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = i, .body_b = i + 1, .real_body_a = i, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_chain_5()
{
	TEST_BEGIN("integ_chain_5");
	IntegrationWorld iw = {0};
	iw.body_count = 6;
	for (int i = 0; i < 6; i++) iw.bodies[i] = make_body((float)(1 + i), (float)(2 + i));
	// joint_count set below (was bs: 5
	for (int i = 0; i < 5; i++)
		iw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f,0,0), .r_b = V3(-0.5f,0,0), .body_a = i, .body_b = i + 1 };

	LDL_Cache c = {0};
	for (int i = 0; i < 5; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = i, .body_b = i + 1, .real_body_a = i, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 5;

	double rhs[15];
	for (int i = 0; i < 15; i++) rhs[i] = (i % 3 == 0) ? 1.0 : (i % 3 == 1) ? -0.5 : 0.2;
	double x_sparse[15], x_dense[15];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 15);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_triangle()
{
	TEST_BEGIN("integ_triangle");
	// 3 bodies, 3 joints: 0-1, 1-2, 0-2. Every pair coupled.
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	iw.bodies[0] = make_body(2, 4);
	iw.bodies[1] = make_body(3, 6);
	iw.bodies[2] = make_body(4, 8);
	// joint_count set below (was bs: 3
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0.5f,0,0), .r_b = V3(-0.5f,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0.5f,0), .r_b = V3(0,-0.5f,0), .body_a = 1, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,0.5f), .r_b = V3(0,0,-0.5f), .body_a = 0, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 2, .real_body_a = 0, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 2 },
	};
	for (int i = 0; i < 3; i++) apush(c.constraints, cons[i]);
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_star_4()
{
	TEST_BEGIN("integ_star_4");
	// Hub body 0, leaves 1,2,3. Joints: 0-1, 0-2, 0-3.
	IntegrationWorld iw = {0};
	iw.body_count = 4;
	iw.bodies[0] = make_body(5, 10); // hub
	iw.bodies[1] = make_body(2, 4);
	iw.bodies[2] = make_body(3, 6);
	iw.bodies[3] = make_body(4, 8);
	// joint_count set below (was bs: 3
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,1,0), .r_b = V3(0,-1,0), .body_a = 0, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,1), .r_b = V3(0,0,-1), .body_a = 0, .body_b = 3 };

	LDL_Cache c = {0};
	for (int i = 0; i < 3; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_star_6()
{
	TEST_BEGIN("integ_star_6");
	// Hub body 0, 6 leaves. 6 joints = 18 DOF. Eliminating leaves creates
	// no fill, but the hub node has 6 neighbors -- all pair up in Schur updates
	// when the hub is eliminated (or when leaves reduce hub degree to tie).
	IntegrationWorld iw = {0};
	iw.body_count = 7;
	iw.bodies[0] = make_body(10, 20); // hub
	for (int i = 1; i <= 6; i++) iw.bodies[i] = make_body((float)(1 + i), (float)(2 + i));
	// joint_count set below (was bs: 6
	v3 dirs[6] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
	for (int i = 0; i < 6; i++)
		iw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = dirs[i], .r_b = scale(dirs[i], -1), .body_a = 0, .body_b = i + 1 };

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	double rhs[18];
	for (int i = 0; i < 18; i++) rhs[i] = (double)(i % 5) * 0.3 - 0.6;
	double x_sparse[18], x_dense[18];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 18);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_star_extreme_mass()
{
	TEST_BEGIN("integ_star_extreme_mass");
	// Heavy hub (mass=10000) + light leaves (mass=1). Hub contributes almost
	// nothing to K, so off-diagonal coupling through hub is tiny.
	IntegrationWorld iw = {0};
	iw.body_count = 5;
	iw.bodies[0] = make_body(10000, 10000); // heavy hub
	for (int i = 1; i <= 4; i++) iw.bodies[i] = make_body(1, 1);
	// joint_count set below (was bs: 4
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(2,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(-2,0,0), .r_b = V3(1,0,0), .body_a = 0, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,2,0), .r_b = V3(0,-1,0), .body_a = 0, .body_b = 3 };
	iw.sol_joints[3] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,-2,0), .r_b = V3(0,1,0), .body_a = 0, .body_b = 4 };

	LDL_Cache c = {0};
	for (int i = 0; i < 4; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 4;

	double rhs[12] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7, -0.1, 0.4, -0.6 };
	double x_sparse[12], x_dense[12];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 12);
	TEST_ASSERT(err < 0.1);
	for (int i = 0; i < 12; i++) TEST_ASSERT(x_sparse[i] == x_sparse[i]);

	integration_cache_free(&c);
}

static void test_integ_star_mixed_types()
{
	TEST_BEGIN("integ_star_mixed_types");
	// Hub with mixed spokes: ball-socket (3), distance (1), hinge (5).
	// 9 total DOF, 3 different block sizes coupling through hub.
	IntegrationWorld iw = {0};
	iw.body_count = 4;
	iw.bodies[0] = make_body(5, 10); // hub
	iw.bodies[1] = make_body(2, 4);
	iw.bodies[2] = make_body(3, 6);
	iw.bodies[3] = make_body(4, 8);
	// joint_count set below (was bs: 1
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,1,0), .r_b = V3(0,-1,0), .dist.axis = norm(V3(0,1,0)), .body_a = 0, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_HINGE, .dof = 5, .r_a = V3(0,0,1), .r_b = V3(0,0,-1), .hinge.u1 = V3(1,0,0), .hinge.u2 = V3(0,1,0), .body_a = 0, .body_b = 3 };
	iw.joint_count = 3;

	LDL_Cache c = {0};
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 2, .real_body_a = 0, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 0, .body_b = 3, .real_body_a = 0, .real_body_b = 3, .weight_a = 1, .weight_b = 1, .solver_idx = 2 },
	};
	for (int i = 0; i < 3; i++) apush(c.constraints, cons[i]);
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, 0.3, -0.3, 0.8, 0.1, 0.5, -0.2 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

// ============================================================================
// Tier 2: Physics stress

static void test_integ_chain_2_extreme_mass()
{
	TEST_BEGIN("integ_chain_2_extreme_mass");
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	iw.bodies[0] = make_body(1, 1);        // light
	iw.bodies[1] = make_body(10000, 10000); // heavy
	iw.bodies[2] = make_body(1, 1);         // light
	// joint_count set below (was bs: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[6] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1 };
	double x_sparse[6], x_dense[6];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 6);
	TEST_ASSERT(err < 0.1); // looser tolerance for extreme mass ratio
	for (int i = 0; i < 6; i++) TEST_ASSERT(x_sparse[i] == x_sparse[i]);

	integration_cache_free(&c);
}

static void test_integ_chain_2_large_levers()
{
	TEST_BEGIN("integ_chain_2_large_levers");
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	for (int i = 0; i < 3; i++) iw.bodies[i] = make_body(1, 1);
	// joint_count set below (was bs: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(50,0,0), .r_b = V3(-50,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,50,0), .r_b = V3(0,-50,0), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[6] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1 };
	double x_sparse[6], x_dense[6];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 6);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_chain_2_rotated_asymmetric()
{
	TEST_BEGIN("integ_chain_2_rotated_asymmetric");
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	quat r1 = quat_axis_angle(norm(V3(1, 2, -1)), 0.8f);
	quat r2 = quat_axis_angle(norm(V3(-1, 0, 3)), 1.2f);
	iw.bodies[0] = make_body_full(2, V3(0.5f, 5, 20), r1);
	iw.bodies[1] = make_body_full(8, V3(3, 10, 1), r2);
	iw.bodies[2] = make_body_full(1, V3(1, 1, 1), quat_identity());
	// joint_count set below (was bs: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(3,-1,2), .r_b = V3(-2,1,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,2,-1), .r_b = V3(0,-1,3), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[6] = { 2, -1, 0.5, -0.3, 1.5, -0.7 };
	double x_sparse[6], x_dense[6];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 6);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_chain_3_mixed_mass()
{
	TEST_BEGIN("integ_chain_3_mixed_mass");
	// 1:100:1 pattern: light-heavy-light-heavy
	IntegrationWorld iw = {0};
	iw.body_count = 4;
	iw.bodies[0] = make_body(1, 1);
	iw.bodies[1] = make_body(100, 100);
	iw.bodies[2] = make_body(1, 1);
	iw.bodies[3] = make_body(100, 100);
	// joint_count set below (was bs: 3
	for (int i = 0; i < 3; i++)
		iw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = i, .body_b = i + 1 };

	LDL_Cache c = {0};
	for (int i = 0; i < 3; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = i, .body_b = i + 1, .real_body_a = i, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.1);

	integration_cache_free(&c);
}

static void test_integ_triangle_extreme_mass()
{
	TEST_BEGIN("integ_triangle_extreme_mass");
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	iw.bodies[0] = make_body(1, 1);
	iw.bodies[1] = make_body(10000, 10000);
	iw.bodies[2] = make_body(1, 1);
	// joint_count set below (was bs: 3
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,1,0), .r_b = V3(0,-1,0), .body_a = 1, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(0,0,1), .r_b = V3(0,0,-1), .body_a = 0, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 2, .real_body_a = 0, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 2 },
	};
	for (int i = 0; i < 3; i++) apush(c.constraints, cons[i]);
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, -0.3, 0.8, 0.1, 0.5, -0.2, 0.7 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.1);

	integration_cache_free(&c);
}

// ============================================================================
// Tier 3: Mixed constraint types

static void test_integ_chain_distance()
{
	TEST_BEGIN("integ_chain_distance");
	// Two distance constraints (1 DOF each). Total DOF = 2.
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	for (int i = 0; i < 3; i++) iw.bodies[i] = make_body(2, 4);
	// joint_count set below (was dist: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .dist.axis = norm(V3(1,0,0)), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .dist.axis = norm(V3(1,0,0)), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[2] = { 0.5, -0.3 };
	double x_sparse[2], x_dense[2];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 2);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_chain_hinge()
{
	TEST_BEGIN("integ_chain_hinge");
	// Two hinge constraints (5 DOF each). Total DOF = 10.
	IntegrationWorld iw = {0};
	iw.body_count = 3;
	for (int i = 0; i < 3; i++) iw.bodies[i] = make_body(3, 6);
	// joint_count set below (was hinge: 2
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_HINGE, .dof = 5, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .hinge.u1 = V3(1,0,0), .hinge.u2 = V3(0,1,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_HINGE, .dof = 5, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .hinge.u1 = V3(1,0,0), .hinge.u2 = V3(0,1,0), .body_a = 1, .body_b = 2 };

	LDL_Cache c = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(c.constraints, cons[i]);
	c.joint_count = 2;

	double rhs[10] = { 1, -0.5, 0.2, 0.3, -0.1, -0.3, 0.8, 0.1, 0.5, -0.2 };
	double x_sparse[10], x_dense[10];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 10);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

static void test_integ_mixed_all()
{
	TEST_BEGIN("integ_mixed_all");
	// BS(3) + Distance(1) + Hinge(5) = 9 DOF, 3 constraints on 4 bodies.
	IntegrationWorld iw = {0};
	iw.body_count = 4;
	for (int i = 0; i < 4; i++) iw.bodies[i] = make_body((float)(2 + i), (float)(4 + i * 2));
	// joint_count set below (was bs: 1
	iw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1,0,0), .r_b = V3(-1,0,0), .body_a = 0, .body_b = 1 };
	iw.sol_joints[1] = (SolverJoint){ .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(0,1,0), .r_b = V3(0,-1,0), .dist.axis = norm(V3(0,1,0)), .body_a = 1, .body_b = 2 };
	iw.sol_joints[2] = (SolverJoint){ .type = JOINT_HINGE, .dof = 5, .r_a = V3(0,0,1), .r_b = V3(0,0,-1), .hinge.u1 = V3(1,0,0), .hinge.u2 = V3(0,1,0), .body_a = 2, .body_b = 3 };
	iw.joint_count = 3;

	LDL_Cache c = {0};
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 2, .body_b = 3, .real_body_a = 2, .real_body_b = 3, .weight_a = 1, .weight_b = 1, .solver_idx = 2 },
	};
	for (int i = 0; i < 3; i++) apush(c.constraints, cons[i]);
	c.joint_count = 3;

	double rhs[9] = { 1, -0.5, 0.2, 0.3, -0.3, 0.8, 0.1, 0.5, -0.2 };
	double x_sparse[9], x_dense[9];
	integration_sparse_solve(&c, &iw, rhs, x_sparse);
	integration_dense_solve(&c, &iw, rhs, x_dense);

	double err = max_solution_error(x_sparse, x_dense, 9);
	TEST_ASSERT(err < 0.01);

	integration_cache_free(&c);
}

// ============================================================================
// Runner

static void run_integration_ldl_tests()
{
	printf("--- LDL integration tests ---\n");

	// Tier 1: structural
	test_integ_single_ball_socket();
	test_integ_chain_2();
	test_integ_chain_3();
	test_integ_chain_5();
	test_integ_triangle();
	test_integ_star_4();

	test_integ_star_6();
	test_integ_star_extreme_mass();
	test_integ_star_mixed_types();

	// Tier 2: physics stress
	test_integ_chain_2_extreme_mass();
	test_integ_chain_2_large_levers();
	test_integ_chain_2_rotated_asymmetric();
	test_integ_chain_3_mixed_mass();
	test_integ_triangle_extreme_mass();

	// Tier 3: mixed types
	test_integ_chain_distance();
	test_integ_chain_hinge();
	test_integ_mixed_all();
}
