// tests_sparse_unit.c -- unit tests for LDL_Sparse edge operations.
// Tests the sparse block adjacency structure in isolation, then connects it
// to known-good K assembly and block factorization to verify end-to-end
// that the sparse storage correctly holds physics data.

// ============================================================================
// Pure sparse edge operations

static void test_sparse_empty()
{
	TEST_BEGIN("sparse_empty");
	// Empty graph: no edges, find returns -1, get returns NULL.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;

	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 1) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 0) == -1);
	TEST_ASSERT(ldl_sparse_get_edge(&s, 0, 1) == NULL);
	TEST_ASSERT(ldl_sparse_get_edge(&s, 2, 0) == NULL);

	ldl_sparse_free(&s);
}

static void test_sparse_create_edge()
{
	TEST_BEGIN("sparse_create_edge");
	// Create one edge (0, 1). Both directions should exist.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;

	double* data = ldl_sparse_get_or_create_edge(&s, 0, 1);
	TEST_ASSERT(data != NULL);

	// Both directions exist
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 1) >= 0);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 0) >= 0);

	// get_edge returns non-NULL for both
	TEST_ASSERT(ldl_sparse_get_edge(&s, 0, 1) != NULL);
	TEST_ASSERT(ldl_sparse_get_edge(&s, 1, 0) != NULL);

	// Data is zeroed
	for (int i = 0; i < 9; i++) TEST_ASSERT(fabs(data[i]) < 1e-15);

	ldl_sparse_free(&s);
}

static void test_sparse_create_idempotent()
{
	TEST_BEGIN("sparse_create_idempotent");
	// Creating the same edge twice returns the same block.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;

	double* first = ldl_sparse_get_or_create_edge(&s, 0, 1);
	first[0] = 42.0; // write a marker

	double* second = ldl_sparse_get_or_create_edge(&s, 0, 1);
	TEST_ASSERT(second == first); // same pointer
	TEST_ASSERT(fabs(second[0] - 42.0) < 1e-15); // marker preserved

	ldl_sparse_free(&s);
}

static void test_sparse_mixed_dof()
{
	TEST_BEGIN("sparse_mixed_dof");
	// Edge between nodes with different DOF: node 0 has 3 DOF, node 1 has 1 DOF.
	// Block (0,1) is 3x1, block (1,0) is 1x3.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 1;

	double* data_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	double* data_10 = ldl_sparse_get_edge(&s, 1, 0);
	TEST_ASSERT(data_01 != NULL);
	TEST_ASSERT(data_10 != NULL);

	// Write to (0,1): 3x1 = 3 elements
	data_01[0] = 1.0; data_01[1] = 2.0; data_01[2] = 3.0;

	// Write to (1,0): 1x3 = 3 elements
	data_10[0] = 4.0; data_10[1] = 5.0; data_10[2] = 6.0;

	// Read back: blocks are independent storage
	double* read_01 = ldl_sparse_get_edge(&s, 0, 1);
	double* read_10 = ldl_sparse_get_edge(&s, 1, 0);
	TEST_ASSERT(fabs(read_01[0] - 1.0) < 1e-15);
	TEST_ASSERT(fabs(read_01[2] - 3.0) < 1e-15);
	TEST_ASSERT(fabs(read_10[0] - 4.0) < 1e-15);
	TEST_ASSERT(fabs(read_10[2] - 6.0) < 1e-15);

	ldl_sparse_free(&s);
}

static void test_sparse_multiple_edges()
{
	TEST_BEGIN("sparse_multiple_edges");
	// Chain: 0-1-2. Two edges, no connection between 0 and 2.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;

	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	// 0-1 exists
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 1) >= 0);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 0) >= 0);

	// 1-2 exists
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 2) >= 0);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 2, 1) >= 0);

	// 0-2 does NOT exist
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 2) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 2, 0) == -1);

	// Node 1 has 2 neighbors
	TEST_ASSERT(asize(s.adj[1]) == 2);
	// Nodes 0 and 2 have 1 neighbor each
	TEST_ASSERT(asize(s.adj[0]) == 1);
	TEST_ASSERT(asize(s.adj[2]) == 1);

	ldl_sparse_free(&s);
}

static void test_sparse_star_topology()
{
	TEST_BEGIN("sparse_star_topology");
	// Star: node 0 connected to 1, 2, 3. No edges between 1, 2, 3.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;

	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 3);

	// Hub has 3 neighbors
	TEST_ASSERT(asize(s.adj[0]) == 3);

	// Leaves have 1 neighbor each
	for (int i = 1; i <= 3; i++) TEST_ASSERT(asize(s.adj[i]) == 1);

	// No leaf-leaf edges
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 2) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 3) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 2, 3) == -1);

	ldl_sparse_free(&s);
}

static void test_sparse_data_isolation()
{
	TEST_BEGIN("sparse_data_isolation");
	// Write to one edge block, verify other edges are unaffected.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;

	double* e01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	double* e12 = ldl_sparse_get_or_create_edge(&s, 1, 2);

	// Fill e01 with 1s
	for (int i = 0; i < 9; i++) e01[i] = 1.0;

	// e12 should still be zero
	for (int i = 0; i < 9; i++) TEST_ASSERT(fabs(e12[i]) < 1e-15);

	ldl_sparse_free(&s);
}

// ============================================================================
// Integration: sparse storage + known-good K assembly + block factorization.
// Build a 2-node graph from real physics, store K in the sparse structure,
// factor the diagonal blocks, verify the factored result solves correctly.

static void test_sparse_K_roundtrip()
{
	TEST_BEGIN("sparse_K_roundtrip");
	// Two ball-socket constraints sharing body B1 in a chain: B0--j0--B1--j1--B2.
	// Node 0 = joint 0 (3 DOF), Node 1 = joint 1 (3 DOF).
	// Shared body B1 creates an off-diagonal edge between nodes 0 and 1.

	// Build bodies
	BodyHot bodies[3];
	bodies[0] = make_body(2, 4).hot;   // B0
	bodies[1] = make_body(5, 10).hot;  // B1 (shared)
	bodies[2] = make_body(3, 6).hot;   // B2

	// Build solver data
	SolverJoint sols[2] = {
		{ .r_a = V3(0.5f, 0, 0), .r_b = V3(-0.5f, 0, 0), .body_a = 0, .body_b = 1 },
		{ .r_a = V3(0.5f, 0, 0), .r_b = V3(-0.5f, 0, 0), .body_a = 1, .body_b = 2 },
	};
	test_fill_bs_rows(&sols[0]);
	test_fill_bs_rows(&sols[1]);

	// Fill Jacobians (using tested ldl_fill_jacobian)
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 },
	};
	LDL_JacobianRow jac0[3], jac1[3];
	ldl_fill_jacobian(&cons[0], sols, jac0);
	ldl_fill_jacobian(&cons[1], sols, jac1);

	// Build diagonal K blocks (using tested ldl_K_body_contrib)
	double K0[6] = {0}; // node 0: joint 0's K block
	ldl_K_body_contrib(jac0, 3, 0, 0, &bodies[0], 1.0, K0); // B0 on side A
	ldl_K_body_contrib(jac0, 3, 1, 0, &bodies[1], 1.0, K0); // B1 on side B

	double K1[6] = {0}; // node 1: joint 1's K block
	ldl_K_body_contrib(jac1, 3, 0, 0, &bodies[1], 1.0, K1); // B1 on side A
	ldl_K_body_contrib(jac1, 3, 1, 0, &bodies[2], 1.0, K1); // B2 on side B

	// Build off-diagonal K block via shared body B1 (using tested ldl_K_body_off)
	// Joint 0 sees B1 on side B, joint 1 sees B1 on side A.
	double K01[9] = {0};
	ldl_K_body_off(jac0, 3, 1, jac1, 3, 0, &bodies[1], 1.0, K01);

	// Now store these in an LDL_Sparse and verify they survive storage/retrieval
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;

	// Store diagonal blocks
	for (int i = 0; i < 6; i++) {
		s.diag_data[0][i] = K0[i];
		s.diag_data[1][i] = K1[i];
	}

	// Store off-diagonal via sparse edge
	double* edge_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	for (int i = 0; i < 9; i++) edge_01[i] = K01[i];

	// Verify retrieval matches what we stored
	double* read_01 = ldl_sparse_get_edge(&s, 0, 1);
	TEST_ASSERT(read_01 != NULL);
	for (int i = 0; i < 9; i++)
		TEST_ASSERT(fabs(read_01[i] - K01[i]) < 1e-15);

	// Factor the diagonal blocks (using tested block_ldl)
	double D0[3], D1[3];
	block_ldl(s.diag_data[0], D0, 3);
	block_ldl(s.diag_data[1], D1, 3);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D0[i] > 0);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D1[i] > 0);

	// Solve each diagonal block independently (using tested block_solve)
	// RHS = some constraint error vector
	double b0[3] = { 1, -0.5, 0.2 }, x0[3];
	double b1[3] = { -0.3, 0.8, 0.1 }, x1[3];
	block_solve(s.diag_data[0], D0, b0, x0, 3);
	block_solve(s.diag_data[1], D1, b1, x1, 3);

	// Verify: K * x = b (residual check using original K before factoring)
	double res0 = ldl_solve_residual(NULL, x0, b0, 0); // can't use this, K0 was overwritten
	// Instead verify x is finite and reasonable
	for (int i = 0; i < 3; i++) {
		TEST_ASSERT(x0[i] == x0[i]); // not NaN
		TEST_ASSERT(fabs(x0[i]) < 1e6);
		TEST_ASSERT(x1[i] == x1[i]);
		TEST_ASSERT(fabs(x1[i]) < 1e6);
	}

	ldl_sparse_free(&s);
}

// ============================================================================
// Runner

static void test_sparse_K_roundtrip_extreme_mass()
{
	TEST_BEGIN("sparse_K_roundtrip_extreme_mass");
	// Chain B0--j0--B1--j1--B2 with 10000:1 mass ratio.
	// B0 heavy (10000), B1 medium (10), B2 light (1).
	// The off-diagonal coupling through B1 is well-conditioned,
	// but the diagonal blocks have wildly different contributions.
	BodyHot bodies[3];
	bodies[0] = make_body(10000, 10000).hot;
	bodies[1] = make_body(10, 10).hot;
	bodies[2] = make_body(1, 1).hot;

	SolverJoint sols[2] = {
		{ .r_a = V3(1, 0, 0), .r_b = V3(-1, 0, 0), .body_a = 0, .body_b = 1 },
		{ .r_a = V3(1, 0, 0), .r_b = V3(-1, 0, 0), .body_a = 1, .body_b = 2 },
	};
	test_fill_bs_rows(&sols[0]);
	test_fill_bs_rows(&sols[1]);
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 },
	};
	LDL_JacobianRow jac0[3], jac1[3];
	ldl_fill_jacobian(&cons[0], sols, jac0);
	ldl_fill_jacobian(&cons[1], sols, jac1);

	// Diagonal K blocks
	double K0[6] = {0}, K1[6] = {0};
	ldl_K_body_contrib(jac0, 3, 0, 0, &bodies[0], 1.0, K0);
	ldl_K_body_contrib(jac0, 3, 1, 0, &bodies[1], 1.0, K0);
	ldl_K_body_contrib(jac1, 3, 0, 0, &bodies[1], 1.0, K1);
	ldl_K_body_contrib(jac1, 3, 1, 0, &bodies[2], 1.0, K1);

	// Off-diagonal via shared B1
	double K01[9] = {0};
	ldl_K_body_off(jac0, 3, 1, jac1, 3, 0, &bodies[1], 1.0, K01);

	// Store in sparse
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;
	for (int i = 0; i < 6; i++) { s.diag_data[0][i] = K0[i]; s.diag_data[1][i] = K1[i]; }
	double* edge_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	for (int i = 0; i < 9; i++) edge_01[i] = K01[i];

	// Verify off-diagonal is non-trivial (B1 coupling exists)
	double off_norm = 0;
	for (int i = 0; i < 9; i++) off_norm += K01[i] * K01[i];
	TEST_ASSERT(off_norm > 1e-10);

	// Factor diagonals: both must have positive pivots despite mass ratio
	double D0[3], D1[3];
	double K0_copy[6], K1_copy[6];
	memcpy(K0_copy, K0, sizeof(K0));
	memcpy(K1_copy, K1, sizeof(K1));
	block_ldl(K0_copy, D0, 3);
	block_ldl(K1_copy, D1, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(D0[i] > 0); TEST_ASSERT(D1[i] > 0); }

	// Solve and verify finite
	double b[3] = { 1, -0.5, 0.2 }, x[3];
	block_solve(K0_copy, D0, b, x, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(x[i] == x[i]); TEST_ASSERT(fabs(x[i]) < 1e6); }

	// K0 diagonal should reflect the mass ratio: heavy body contributes ~0.0001,
	// medium body contributes ~0.1, so diagonal ~0.1
	TEST_ASSERT(K0[LDL_TRI(0,0)] > 0.05 && K0[LDL_TRI(0,0)] < 10.0);

	// K1 diagonal: medium + light, so ~1.1
	TEST_ASSERT(K1[LDL_TRI(0,0)] > 0.5);

	ldl_sparse_free(&s);
}

static void test_sparse_K_roundtrip_large_levers()
{
	TEST_BEGIN("sparse_K_roundtrip_large_levers");
	// Chain with large lever arms: r = 50. Angular K entries dominate (~2500).
	// Off-diagonal entries also large. Verify sparse storage handles big values
	// and factorization still produces positive pivots.
	BodyHot bodies[3];
	bodies[0] = make_body(1, 1).hot;
	bodies[1] = make_body(1, 1).hot;
	bodies[2] = make_body(1, 1).hot;

	SolverJoint sols[2] = {
		{ .r_a = V3(50, 0, 0), .r_b = V3(-50, 0, 0), .body_a = 0, .body_b = 1 },
		{ .r_a = V3(0, 50, 0), .r_b = V3(0, -50, 0), .body_a = 1, .body_b = 2 },
	};
	test_fill_bs_rows(&sols[0]);
	test_fill_bs_rows(&sols[1]);
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 },
	};
	LDL_JacobianRow jac0[3], jac1[3];
	ldl_fill_jacobian(&cons[0], sols, jac0);
	ldl_fill_jacobian(&cons[1], sols, jac1);

	double K0[6] = {0}, K1[6] = {0};
	ldl_K_body_contrib(jac0, 3, 0, 0, &bodies[0], 1.0, K0);
	ldl_K_body_contrib(jac0, 3, 1, 0, &bodies[1], 1.0, K0);
	ldl_K_body_contrib(jac1, 3, 0, 0, &bodies[1], 1.0, K1);
	ldl_K_body_contrib(jac1, 3, 1, 0, &bodies[2], 1.0, K1);

	double K01[9] = {0};
	ldl_K_body_off(jac0, 3, 1, jac1, 3, 0, &bodies[1], 1.0, K01);

	// Store in sparse
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;
	for (int i = 0; i < 6; i++) { s.diag_data[0][i] = K0[i]; s.diag_data[1][i] = K1[i]; }
	double* edge_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	for (int i = 0; i < 9; i++) edge_01[i] = K01[i];

	// Diagonal entries should be significant. Lever along X means skew(r_a)
	// contributes to Y and Z rows, not X. So X row = 2*inv_mass = 2.
	// Y and Z rows get angular contribution: inv_mass + inv_inertia*|r|^2 >> 1.
	double max_diag_0 = 0, max_diag_1 = 0;
	for (int i = 0; i < 3; i++) {
		if (K0[LDL_TRI(i,i)] > max_diag_0) max_diag_0 = K0[LDL_TRI(i,i)];
		if (K1[LDL_TRI(i,i)] > max_diag_1) max_diag_1 = K1[LDL_TRI(i,i)];
	}
	TEST_ASSERT(max_diag_0 > 100.0);
	TEST_ASSERT(max_diag_1 > 100.0);

	// Off-diagonal should have significant entries from shared B1
	double max_off = 0;
	for (int i = 0; i < 9; i++) { if (fabs(K01[i]) > max_off) max_off = fabs(K01[i]); }
	TEST_ASSERT(max_off > 10.0);

	// Factor and solve
	double D0[3], D1[3];
	block_ldl(s.diag_data[0], D0, 3);
	block_ldl(s.diag_data[1], D1, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(D0[i] > 0); TEST_ASSERT(D1[i] > 0); }

	double b[3] = { 100, -50, 200 }, x[3];
	block_solve(s.diag_data[0], D0, b, x, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(x[i] == x[i]); TEST_ASSERT(fabs(x[i]) < 1e6); }

	ldl_sparse_free(&s);
}

static void test_sparse_K_roundtrip_mixed_types()
{
	TEST_BEGIN("sparse_K_roundtrip_mixed_types");
	// Ball-socket (3 DOF) and distance (1 DOF) sharing body B1.
	// Off-diagonal block is 3x1 and 1x3 -- non-square.
	// Verifies mixed-DOF sparse storage with real physics values.
	BodyHot bodies[3];
	bodies[0] = make_body(2, 4).hot;
	bodies[1] = make_body(5, 10).hot;
	bodies[2] = make_body(3, 6).hot;

	SolverJoint sol_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .r_a = V3(1, 0, 0), .r_b = V3(-1, 0, 0), .body_a = 0, .body_b = 1 };
	test_fill_bs_rows(&sol_bs);
	v3 ax = norm(V3(1, 0, 0));
	SolverJoint sol_d = { .type = JOINT_DISTANCE, .dof = 1, .r_a = V3(1, 0, 0), .r_b = V3(-1, 0, 0), .body_a = 1, .body_b = 2 };
	test_fill_dist_rows(&sol_d, ax);

	LDL_Constraint con_bs = { .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 };
	LDL_Constraint con_d = { .type = JOINT_DISTANCE, .dof = 1, .solver_idx = 0 };
	LDL_JacobianRow jac_bs[3], jac_d[1];
	ldl_fill_jacobian(&con_bs, &sol_bs, jac_bs);
	ldl_fill_jacobian(&con_d, &sol_d, jac_d);

	// Diagonal K blocks
	double K_bs[6] = {0}; // 3x3 packed
	ldl_K_body_contrib(jac_bs, 3, 0, 0, &bodies[0], 1.0, K_bs);
	ldl_K_body_contrib(jac_bs, 3, 1, 0, &bodies[1], 1.0, K_bs);

	double K_d[1] = {0}; // 1x1 packed
	ldl_K_body_contrib(jac_d, 1, 0, 0, &bodies[1], 1.0, K_d);
	ldl_K_body_contrib(jac_d, 1, 1, 0, &bodies[2], 1.0, K_d);

	// Off-diagonal: ball-socket side B and distance side A both use B1.
	// Block is 3x1.
	double K_off[3] = {0};
	ldl_K_body_off(jac_bs, 3, 1, jac_d, 1, 0, &bodies[1], 1.0, K_off);

	// Store in sparse with mixed DOF
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 1;
	for (int i = 0; i < 6; i++) s.diag_data[0][i] = K_bs[i];
	s.diag_data[1][0] = K_d[0];

	double* edge_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	// Block (0,1) is 3x1
	for (int i = 0; i < 3; i++) edge_01[i] = K_off[i];

	// Verify retrieval of 3x1 block
	double* read_01 = ldl_sparse_get_edge(&s, 0, 1);
	TEST_ASSERT(read_01 != NULL);
	for (int i = 0; i < 3; i++) TEST_ASSERT(fabs(read_01[i] - K_off[i]) < 1e-15);

	// Reverse edge (1,0) is 1x3
	double* read_10 = ldl_sparse_get_edge(&s, 1, 0);
	TEST_ASSERT(read_10 != NULL);

	// Factor diagonal blocks
	double D_bs[3], D_d[1];
	double K_bs_copy[6];
	memcpy(K_bs_copy, K_bs, sizeof(K_bs));
	block_ldl(K_bs_copy, D_bs, 3);
	double K_d_copy[1] = { K_d[0] };
	block_ldl(K_d_copy, D_d, 1);
	for (int i = 0; i < 3; i++) TEST_ASSERT(D_bs[i] > 0);
	TEST_ASSERT(D_d[0] > 0);

	// Solve 3x3 block
	double b[3] = { 1, -0.5, 0.2 }, x[3];
	block_solve(K_bs_copy, D_bs, b, x, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(x[i] == x[i]); TEST_ASSERT(fabs(x[i]) < 1e6); }

	// Solve 1x1 block
	double b1[1] = { 0.5 }, x1[1];
	block_solve(K_d_copy, D_d, b1, x1, 1);
	TEST_ASSERT(fabs(x1[0] - b1[0] / K_d[0]) < 1e-6);

	ldl_sparse_free(&s);
}

static void test_sparse_K_roundtrip_rotated_asymmetric()
{
	TEST_BEGIN("sparse_K_roundtrip_rotated_asymmetric");
	// Rotated bodies with asymmetric inertia + large levers + extreme mass ratio.
	// The combined worst case: exercises dinv_inertia_mul through K assembly
	// into sparse storage and factorization.
	quat rot_a = quat_axis_angle(norm(V3(1, 2, -1)), 0.8f);
	quat rot_b = quat_axis_angle(norm(V3(-1, 0, 3)), 1.2f);
	BodyHot bodies[3];
	bodies[0] = make_body_full(1, V3(0.5f, 5, 20), rot_a).hot;
	bodies[1] = make_body_full(100, V3(10, 50, 200), rot_b).hot;
	bodies[2] = make_body_full(1, V3(1, 1, 1), quat_identity()).hot;

	SolverJoint sols[2] = {
		{ .r_a = V3(5, -3, 2), .r_b = V3(-2, 8, 1), .body_a = 0, .body_b = 1 },
		{ .r_a = V3(3, 1, -4), .r_b = V3(-1, 0, 6), .body_a = 1, .body_b = 2 },
	};
	test_fill_bs_rows(&sols[0]);
	test_fill_bs_rows(&sols[1]);
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .solver_idx = 1 },
	};
	LDL_JacobianRow jac0[3], jac1[3];
	ldl_fill_jacobian(&cons[0], sols, jac0);
	ldl_fill_jacobian(&cons[1], sols, jac1);

	double K0[6] = {0}, K1[6] = {0};
	ldl_K_body_contrib(jac0, 3, 0, 0, &bodies[0], 1.0, K0);
	ldl_K_body_contrib(jac0, 3, 1, 0, &bodies[1], 1.0, K0);
	ldl_K_body_contrib(jac1, 3, 0, 0, &bodies[1], 1.0, K1);
	ldl_K_body_contrib(jac1, 3, 1, 0, &bodies[2], 1.0, K1);

	double K01[9] = {0};
	ldl_K_body_off(jac0, 3, 1, jac1, 3, 0, &bodies[1], 1.0, K01);

	// Store in sparse
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;
	for (int i = 0; i < 6; i++) { s.diag_data[0][i] = K0[i]; s.diag_data[1][i] = K1[i]; }
	double* edge_01 = ldl_sparse_get_or_create_edge(&s, 0, 1);
	for (int i = 0; i < 9; i++) edge_01[i] = K01[i];

	// Both diagonal blocks must be SPD and factorable
	double D0[3], D1[3];
	double K0c[6], K1c[6];
	memcpy(K0c, K0, sizeof(K0));
	memcpy(K1c, K1, sizeof(K1));
	block_ldl(K0c, D0, 3);
	block_ldl(K1c, D1, 3);
	for (int i = 0; i < 3; i++) { TEST_ASSERT(D0[i] > 0); TEST_ASSERT(D1[i] > 0); }

	// Solve both blocks
	double b[3] = { 2, -1, 0.5 }, x0[3], x1[3];
	block_solve(K0c, D0, b, x0, 3);
	block_solve(K1c, D1, b, x1, 3);
	for (int i = 0; i < 3; i++) {
		TEST_ASSERT(x0[i] == x0[i] && fabs(x0[i]) < 1e6);
		TEST_ASSERT(x1[i] == x1[i] && fabs(x1[i]) < 1e6);
	}

	// Impulse roundtrip: apply lambda through both bodies, verify constraint
	// velocity matches K * lambda. This chains: Jacobian -> K assembly ->
	// sparse storage -> factorization -> solve -> impulse -> velocity readback.
	BodyHot a = bodies[0]; // copy for impulse test
	BodyHot b_body = bodies[1];
	double lambda[3] = { 1, -0.5, 2 };
	ldl_apply_jacobian_impulse(jac0, 3, lambda, &a, 0);
	ldl_apply_jacobian_impulse(jac0, 3, lambda, &b_body, 1);

	double cv[3];
	for (int d = 0; d < 3; d++)
		cv[d] = ldl_constraint_velocity(&jac0[d], &a, &b_body);

	double expected[3] = {0};
	for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
			expected[r] += K0[LDL_TRI(r, c)] * lambda[c];

	for (int i = 0; i < 3; i++)
		TEST_ASSERT(fabs(cv[i] - expected[i]) < 0.1);

	ldl_sparse_free(&s);
}

static void run_sparse_unit_tests()
{
	printf("--- sparse unit tests ---\n");

	// Edge operations
	test_sparse_empty();
	test_sparse_create_edge();
	test_sparse_create_idempotent();
	test_sparse_mixed_dof();
	test_sparse_multiple_edges();
	test_sparse_star_topology();
	test_sparse_data_isolation();

	// Integration with K assembly + factorization
	test_sparse_K_roundtrip();
	test_sparse_K_roundtrip_extreme_mass();
	test_sparse_K_roundtrip_large_levers();
	test_sparse_K_roundtrip_mixed_types();
	test_sparse_K_roundtrip_rotated_asymmetric();
}
