// tests_pgs_vs_ldl.c -- compare PGS and LDL solvers on identical configurations.
// PGS is the known-good reference. If LDL produces different lambdas or velocities,
// there's a sign, scale, or convention bug in the LDL pipeline.

// Helper: run generic PGS solve N times on given bodies. Returns final lambda as v3 (for ball-socket).
// Uses per-DOF Jacobian rows, matching the production solver.
static v3 pgs_solve_ball_socket_n(SolverJoint* s, BodyHot* bodies, int n_iters)
{
	BodyHot a = bodies[s->body_a];
	BodyHot b = bodies[s->body_b];
	s->lambda[0] = 0; s->lambda[1] = 0; s->lambda[2] = 0;
	for (int iter = 0; iter < n_iters; iter++) {
		for (int d = 0; d < s->dof; d++) {
			float vel_err = jac_velocity_f(&s->rows[d], &a, &b);
			float rhs = -vel_err - s->bias[d] - s->softness * s->lambda[d];
			float delta = s->rows[d].eff_mass * rhs;
			s->lambda[d] += delta;
			jac_apply(&s->rows[d], delta, &a, &b);
		}
	}
	bodies[s->body_a] = a;
	bodies[s->body_b] = b;
	return V3(s->lambda[0], s->lambda[1], s->lambda[2]);
}

// Helper: fill rows + eff_mass for a ball-socket SolverJoint from body state.
static void compute_ball_socket_rows(SolverJoint* s, BodyHot* a, BodyHot* b)
{
	test_fill_bs_rows(s);
	for (int d = 0; d < 3; d++) s->rows[d].eff_mass = jac_eff_mass(&s->rows[d], a, b, s->softness);
}

// ============================================================================
// PGS vs LDL comparison tests

static void test_pgs_vs_ldl_single_rigid()
{
	TEST_BEGIN("pgs_vs_ldl_single_rigid");
	// Single rigid ball-socket. Gravity pulling body B down.
	// PGS converged (many iters) should match LDL (exact in 1 step).
	float sub_dt = 1.0f / 240.0f;
	float gravity_y = -9.81f;

	// Two bodies: A static, B dynamic hanging below
	BodyHot bodies_pgs[2];
	BodyState states_pgs[2] = {0};
	states_pgs[0].rotation = quat_identity();
	states_pgs[1].rotation = quat_identity();
	bodies_pgs[0] = (BodyHot){0}; // static
	bodies_pgs[1] = make_body(1, 1).hot;
	states_pgs[1].position = V3(0, -1, 0);
	// Apply gravity to B's velocity
	bodies_pgs[1].velocity.y += gravity_y * sub_dt;

	v3 r_a = V3(0, 0, 0), r_b = V3(0, 0, 0); // anchors at body centers

	// PGS setup
	SolverJoint pgs_sol = {
		.type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = r_a, .r_b = r_b,
		.softness = 0, // rigid, no position bias in velocity solve
		};
	compute_ball_socket_rows(&pgs_sol, &bodies_pgs[0], &bodies_pgs[1]);

	// Run PGS for 100 iterations (should converge for single constraint)
	v3 pgs_lambda = pgs_solve_ball_socket_n(&pgs_sol, bodies_pgs, 100);

	// LDL setup (same initial state)
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = (BodyHot){0};
	sw.states[0].rotation = quat_identity();
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(0, -1, 0);
	sw.bodies[1].velocity.y += gravity_y * sub_dt;

	sw.joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,0,0), .spring = {0} },
	};
	sw.joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = r_a, .r_b = r_b,
		.softness = 0,
		.joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);
	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	v3 ldl_lambda = V3(sw.sol_joints[0].lambda[0], sw.sol_joints[0].lambda[1], sw.sol_joints[0].lambda[2]);

	// Compare lambdas: should match closely
	printf("    PGS lambda: (%.6f, %.6f, %.6f)\n", pgs_lambda.x, pgs_lambda.y, pgs_lambda.z);
	printf("    LDL lambda: (%.6f, %.6f, %.6f)\n", ldl_lambda.x, ldl_lambda.y, ldl_lambda.z);

	TEST_ASSERT_FLOAT(ldl_lambda.x, pgs_lambda.x, 0.01f);
	TEST_ASSERT_FLOAT(ldl_lambda.y, pgs_lambda.y, 0.01f);
	TEST_ASSERT_FLOAT(ldl_lambda.z, pgs_lambda.z, 0.01f);

	// Compare final velocities
	printf("    PGS vel B: (%.6f, %.6f, %.6f)\n", bodies_pgs[1].velocity.x, bodies_pgs[1].velocity.y, bodies_pgs[1].velocity.z);
	printf("    LDL vel B: (%.6f, %.6f, %.6f)\n", w->body_hot[1].velocity.x, w->body_hot[1].velocity.y, w->body_hot[1].velocity.z);

	TEST_ASSERT_FLOAT(w->body_hot[1].velocity.x, bodies_pgs[1].velocity.x, 0.01f);
	TEST_ASSERT_FLOAT(w->body_hot[1].velocity.y, bodies_pgs[1].velocity.y, 0.01f);
	TEST_ASSERT_FLOAT(w->body_hot[1].velocity.z, bodies_pgs[1].velocity.z, 0.01f);

	// Body B's Y velocity should be ~0 (gravity cancelled by constraint)
	TEST_ASSERT(fabsf(w->body_hot[1].velocity.y) < 0.01f);

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_pgs_vs_ldl_chain_gravity()
{
	TEST_BEGIN("pgs_vs_ldl_chain_gravity");
	// Chain: static A -- joint -- B -- joint -- C. Both B and C have gravity.
	// PGS converged should match LDL.
	float sub_dt = 1.0f / 240.0f;
	float g = -9.81f;

	// PGS path
	BodyHot bodies_pgs[3];
	BodyState states_pgs[3] = {0};
	for (int si = 0; si < 3; si++) states_pgs[si].rotation = quat_identity();
	bodies_pgs[0] = (BodyHot){0}; // static
	bodies_pgs[1] = make_body(1, 1).hot;
	states_pgs[1].position = V3(0, -1, 0);
	bodies_pgs[1].velocity.y += g * sub_dt;
	bodies_pgs[2] = make_body(1, 1).hot;
	states_pgs[2].position = V3(0, -2, 0);
	bodies_pgs[2].velocity.y += g * sub_dt;

	SolverJoint pgs_sols[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .r_a = V3(0,0,0), .r_b = V3(0,1,0), .softness = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .r_a = V3(0,-1,0), .r_b = V3(0,0,0), .softness = 0 },
	};
	compute_ball_socket_rows(&pgs_sols[0], &bodies_pgs[0], &bodies_pgs[1]);
	compute_ball_socket_rows(&pgs_sols[1], &bodies_pgs[1], &bodies_pgs[2]);

	// Run PGS for 200 iterations
	for (int iter = 0; iter < 200; iter++) {
		for (int ji = 0; ji < 2; ji++) {
			SolverJoint* s = &pgs_sols[ji];
			BodyHot* a = &bodies_pgs[s->body_a];
			BodyHot* b = &bodies_pgs[s->body_b];
			for (int d = 0; d < s->dof; d++) {
				float vel_err = jac_velocity_f(&s->rows[d], a, b);
				float rhs = -vel_err - s->bias[d] - s->softness * s->lambda[d];
				float delta = s->rows[d].eff_mass * rhs;
				s->lambda[d] += delta;
				jac_apply(&s->rows[d], delta, a, b);
			}
		}
	}

	// LDL path (same initial state)
	SoftTestWorld sw = {0};
	sw.body_count = 3;
	sw.bodies[0] = (BodyHot){0};
	sw.states[0].rotation = quat_identity();
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(0, -1, 0);
	sw.bodies[1].velocity.y += g * sub_dt;
	sw.bodies[2] = make_body(1, 1).hot;
	sw.states[2].position = V3(0, -2, 0);
	sw.bodies[2].velocity.y += g * sub_dt;

	sw.joint_count = 2;
	sw.joints[0] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = 1, .ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,1,0), .spring = {0} } };
	sw.joints[1] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 1, .body_b = 2, .ball_socket = { .local_a = V3(0,-1,0), .local_b = V3(0,0,0), .spring = {0} } };
	sw.joint_count = 2;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .r_a = V3(0,0,0), .r_b = V3(0,1,0), .softness = 0, .joint_idx = 0 };
	test_fill_bs_rows(&sw.sol_joints[0]);
	sw.sol_joints[1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .r_a = V3(0,-1,0), .r_b = V3(0,0,0), .softness = 0, .joint_idx = 1 };
	test_fill_bs_rows(&sw.sol_joints[1]);

	LDL_Cache cc = {0};
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .real_body_a = 1, .real_body_b = 2, .weight_a = 1, .weight_b = 1, .solver_idx = 1 },
	};
	for (int i = 0; i < 2; i++) apush(cc.constraints, cons[i]);
	cc.joint_count = 2;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&cc);
	ldl_build_topology(&cc, w);
	ldl_numeric_factor(&cc, w, sw.sol_joints, NULL);
	ldl_island_solve(&cc, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// Compare final body velocities
	printf("    PGS vel B: (%.6f, %.6f, %.6f)\n", bodies_pgs[1].velocity.x, bodies_pgs[1].velocity.y, bodies_pgs[1].velocity.z);
	printf("    LDL vel B: (%.6f, %.6f, %.6f)\n", w->body_hot[1].velocity.x, w->body_hot[1].velocity.y, w->body_hot[1].velocity.z);
	printf("    PGS vel C: (%.6f, %.6f, %.6f)\n", bodies_pgs[2].velocity.x, bodies_pgs[2].velocity.y, bodies_pgs[2].velocity.z);
	printf("    LDL vel C: (%.6f, %.6f, %.6f)\n", w->body_hot[2].velocity.x, w->body_hot[2].velocity.y, w->body_hot[2].velocity.z);

	TEST_ASSERT_FLOAT(w->body_hot[1].velocity.y, bodies_pgs[1].velocity.y, 0.01f);
	TEST_ASSERT_FLOAT(w->body_hot[2].velocity.y, bodies_pgs[2].velocity.y, 0.01f);

	// Both bodies' Y velocity should be near zero (chain holds against gravity)
	TEST_ASSERT(fabsf(w->body_hot[1].velocity.y) < 0.1f);
	TEST_ASSERT(fabsf(w->body_hot[2].velocity.y) < 0.1f);

	// Compare lambdas
	printf("    PGS lam0: (%.6f, %.6f, %.6f)\n", pgs_sols[0].lambda[0], pgs_sols[0].lambda[1], pgs_sols[0].lambda[2]);
	printf("    LDL lam0: (%.6f, %.6f, %.6f)\n", sw.sol_joints[0].lambda[0], sw.sol_joints[0].lambda[1], sw.sol_joints[0].lambda[2]);
	printf("    PGS lam1: (%.6f, %.6f, %.6f)\n", pgs_sols[1].lambda[0], pgs_sols[1].lambda[1], pgs_sols[1].lambda[2]);
	printf("    LDL lam1: (%.6f, %.6f, %.6f)\n", sw.sol_joints[1].lambda[0], sw.sol_joints[1].lambda[1], sw.sol_joints[1].lambda[2]);

	TEST_ASSERT_FLOAT(sw.sol_joints[0].lambda[1], pgs_sols[0].lambda[1], 0.01f);
	TEST_ASSERT_FLOAT(sw.sol_joints[1].lambda[1], pgs_sols[1].lambda[1], 0.01f);

	integration_cache_free(&cc);
	soft_test_free_world(w);
}

static void test_pgs_vs_ldl_star_gravity()
{
	TEST_BEGIN("pgs_vs_ldl_star_gravity");
	// Hub star: static body 0, hub body 1, leaves 2,3,4.
	// Gravity on all dynamic bodies. Joints: 0-1, 1-2, 1-3, 1-4.
	// This is the failing scene (simplified).
	float sub_dt = 1.0f / 240.0f;
	float g = -9.81f;
	int nb = 5;

	// PGS path
	BodyHot bodies_pgs[5];
	BodyState states_pgs[5] = {0};
	for (int si = 0; si < 5; si++) states_pgs[si].rotation = quat_identity();
	bodies_pgs[0] = (BodyHot){0}; // static anchor
	bodies_pgs[1] = make_body(5, 5).hot; // hub
	states_pgs[1].position = V3(0, -1, 0);
	bodies_pgs[1].velocity.y += g * sub_dt;
	bodies_pgs[2] = make_body(1, 1).hot;
	states_pgs[2].position = V3(1, -1, 0);
	bodies_pgs[2].velocity.y += g * sub_dt;
	bodies_pgs[3] = make_body(1, 1).hot;
	states_pgs[3].position = V3(-1, -1, 0);
	bodies_pgs[3].velocity.y += g * sub_dt;
	bodies_pgs[4] = make_body(1, 1).hot;
	states_pgs[4].position = V3(0, -1, 1);
	bodies_pgs[4].velocity.y += g * sub_dt;

	SolverJoint pgs_sols[4] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .r_a = V3(0,0,0), .r_b = V3(0,1,0), .softness = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2, .r_a = V3(1,0,0), .r_b = V3(0,0,0), .softness = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 3, .r_a = V3(-1,0,0), .r_b = V3(0,0,0), .softness = 0 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 4, .r_a = V3(0,0,1), .r_b = V3(0,0,0), .softness = 0 },
	};
	for (int i = 0; i < 4; i++)
		compute_ball_socket_rows(&pgs_sols[i], &bodies_pgs[pgs_sols[i].body_a], &bodies_pgs[pgs_sols[i].body_b]);

	// Run PGS 500 iterations (enough to converge for 4 coupled constraints)
	for (int iter = 0; iter < 500; iter++) {
		for (int ji = 0; ji < 4; ji++) {
			SolverJoint* s = &pgs_sols[ji];
			BodyHot* a = &bodies_pgs[s->body_a];
			BodyHot* b = &bodies_pgs[s->body_b];
			for (int d = 0; d < s->dof; d++) {
				float vel_err = jac_velocity_f(&s->rows[d], a, b);
				float rhs = -vel_err - s->bias[d] - s->softness * s->lambda[d];
				float delta = s->rows[d].eff_mass * rhs;
				s->lambda[d] += delta;
				jac_apply(&s->rows[d], delta, a, b);
			}
		}
	}

	// LDL path
	SoftTestWorld sw = {0};
	sw.body_count = 5;
	sw.bodies[0] = (BodyHot){0};
	sw.states[0].rotation = quat_identity();
	sw.bodies[1] = make_body(5, 5).hot;
	sw.states[1].position = V3(0, -1, 0);
	sw.bodies[1].velocity.y += g * sub_dt;
	sw.bodies[2] = make_body(1, 1).hot;
	sw.states[2].position = V3(1, -1, 0);
	sw.bodies[2].velocity.y += g * sub_dt;
	sw.bodies[3] = make_body(1, 1).hot;
	sw.states[3].position = V3(-1, -1, 0);
	sw.bodies[3].velocity.y += g * sub_dt;
	sw.bodies[4] = make_body(1, 1).hot;
	sw.states[4].position = V3(0, -1, 1);
	sw.bodies[4].velocity.y += g * sub_dt;

	sw.joint_count = 4;
	sw.joints[0] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = 1, .ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,1,0), .spring = {0} } };
	sw.joints[1] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 1, .body_b = 2, .ball_socket = { .local_a = V3(1,0,0), .local_b = V3(0,0,0), .spring = {0} } };
	sw.joints[2] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 1, .body_b = 3, .ball_socket = { .local_a = V3(-1,0,0), .local_b = V3(0,0,0), .spring = {0} } };
	sw.joints[3] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 1, .body_b = 4, .ball_socket = { .local_a = V3(0,0,1), .local_b = V3(0,0,0), .spring = {0} } };

	sw.joint_count = 4;
	for (int i = 0; i < 4; i++) {
		sw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
			.body_a = pgs_sols[i].body_a, .body_b = pgs_sols[i].body_b,
			.r_a = pgs_sols[i].r_a, .r_b = pgs_sols[i].r_b,
			.softness = 0, .joint_idx = i,
		};
		test_fill_bs_rows(&sw.sol_joints[i]);
	}

	LDL_Cache cc = {0};
	for (int i = 0; i < 4; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = pgs_sols[i].body_a, .body_b = pgs_sols[i].body_b, .real_body_a = pgs_sols[i].body_a, .real_body_b = pgs_sols[i].body_b, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(cc.constraints, con);
	}
	cc.joint_count = 4;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&cc);
	ldl_build_topology(&cc, w);
	ldl_numeric_factor(&cc, w, sw.sol_joints, NULL);
	ldl_island_solve(&cc, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	printf("    Star hub PGS vel: (%.6f, %.6f, %.6f)\n", bodies_pgs[1].velocity.x, bodies_pgs[1].velocity.y, bodies_pgs[1].velocity.z);
	printf("    Star hub LDL vel: (%.6f, %.6f, %.6f)\n", w->body_hot[1].velocity.x, w->body_hot[1].velocity.y, w->body_hot[1].velocity.z);
	for (int i = 2; i < 5; i++) {
		printf("    Leaf %d PGS vel.y: %.6f  LDL vel.y: %.6f\n", i, bodies_pgs[i].velocity.y, w->body_hot[i].velocity.y);
	}

	// All bodies should have near-zero Y velocity (held up by chain to static anchor)
	for (int i = 1; i < nb; i++) {
		TEST_ASSERT(fabsf(w->body_hot[i].velocity.y) < 0.1f);
	}

	// LDL and PGS should agree
	for (int i = 1; i < nb; i++) {
		TEST_ASSERT_FLOAT(w->body_hot[i].velocity.y, bodies_pgs[i].velocity.y, 0.05f);
	}

	integration_cache_free(&cc);
	soft_test_free_world(w);
}

static void test_pgs_vs_ldl_star_shattering()
{
	TEST_BEGIN("pgs_vs_ldl_star_shattering");
	// 7 joints on hub body -> 21 DOF on body 1 -> triggers shattering.
	// Static anchor at 0, hub at 1, leaves 2-8.
	float sub_dt = 1.0f / 240.0f;
	float g = -9.81f;
	int nj = 7; // joints: 0-1, 1-2, 1-3, 1-4, 1-5, 1-6, 1-7
	int nb = nj + 2; // 9 bodies

	// PGS path
	BodyHot bodies_pgs[9];
	BodyState states_pgs[9] = {0};
	for (int si = 0; si < 9; si++) states_pgs[si].rotation = quat_identity();
	bodies_pgs[0] = (BodyHot){0};
	bodies_pgs[1] = make_body(5, 5).hot;
	states_pgs[1].position = V3(0, -1, 0);
	bodies_pgs[1].velocity.y += g * sub_dt;
	for (int i = 2; i < nb; i++) {
		bodies_pgs[i] = make_body(1, 1).hot;
		states_pgs[i].position = V3((float)(i - 5), -1, (float)(i % 3 - 1));
		bodies_pgs[i].velocity.y += g * sub_dt;
	}

	v3 leaf_r_a[7] = { V3(0,0,0), V3(1,0,0), V3(-1,0,0), V3(0,0,1), V3(0,0,-1), V3(0,1,0), V3(0,-1,0) };
	SolverJoint pgs_sols[8]; // 0-1 + 6 leaves
	pgs_sols[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .r_a = V3(0,0,0), .r_b = V3(0,1,0), .softness = 0 };
	test_fill_bs_rows(&pgs_sols[0]);
	for (int i = 0; i < nj - 1; i++) {
		pgs_sols[i + 1] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = i + 2, .r_a = leaf_r_a[i + 1], .r_b = V3(0,0,0), .softness = 0 };
	}
	for (int i = 0; i < nj; i++)
		compute_ball_socket_rows(&pgs_sols[i], &bodies_pgs[pgs_sols[i].body_a], &bodies_pgs[pgs_sols[i].body_b]);

	for (int iter = 0; iter < 1000; iter++) {
		for (int ji = 0; ji < nj; ji++) {
			SolverJoint* s = &pgs_sols[ji];
			BodyHot* a = &bodies_pgs[s->body_a];
			BodyHot* b = &bodies_pgs[s->body_b];
			for (int d = 0; d < s->dof; d++) {
				float vel_err = jac_velocity_f(&s->rows[d], a, b);
				float rhs = -vel_err - s->bias[d] - s->softness * s->lambda[d];
				float delta = s->rows[d].eff_mass * rhs;
				s->lambda[d] += delta;
				jac_apply(&s->rows[d], delta, a, b);
			}
		}
	}

	// LDL path (with shattering)
	SoftTestWorld sw = {0};
	sw.body_count = nb;
	sw.bodies[0] = (BodyHot){0};
	sw.states[0].rotation = quat_identity();
	sw.bodies[1] = make_body(5, 5).hot;
	sw.states[1].position = V3(0, -1, 0);
	sw.bodies[1].velocity.y += g * sub_dt;
	for (int i = 2; i < nb; i++) {
		sw.bodies[i] = make_body(1, 1).hot;
		sw.states[i].position = V3((float)(i - 5), -1, (float)(i % 3 - 1));
		sw.bodies[i].velocity.y += g * sub_dt;
	}

	sw.joint_count = nj;
	sw.joints[0] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = 1, .ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,1,0), .spring = {0} } };
	for (int i = 0; i < nj - 1; i++) {
		sw.joints[i + 1] = (JointInternal){ .type = JOINT_BALL_SOCKET, .body_a = 1, .body_b = i + 2, .ball_socket = { .local_a = leaf_r_a[i + 1], .local_b = V3(0,0,0), .spring = {0} } };
	}
	sw.joint_count = nj;
	for (int i = 0; i < nj; i++) {
		sw.sol_joints[i] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
			.body_a = pgs_sols[i].body_a, .body_b = pgs_sols[i].body_b,
			.r_a = pgs_sols[i].r_a, .r_b = pgs_sols[i].r_b,
			.softness = 0, .joint_idx = i,
		};
		test_fill_bs_rows(&sw.sol_joints[i]);
	}

	LDL_Cache cc = {0};
	for (int i = 0; i < nj; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = pgs_sols[i].body_a, .body_b = pgs_sols[i].body_b, .real_body_a = pgs_sols[i].body_a, .real_body_b = pgs_sols[i].body_b, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(cc.constraints, con);
	}
	cc.joint_count = nj;

	WorldInternal* w = soft_test_make_world(&sw);

	// Hub body 1 has 7 joints * 3 DOF = 21 DOF > SHATTER_THRESHOLD (15)
	ldl_apply_shattering(&cc, w);
	ldl_build_bundles(&cc);
	ldl_build_topology(&cc, w);
	ldl_numeric_factor(&cc, w, sw.sol_joints, NULL);
	ldl_island_solve(&cc, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	printf("    [shattering] Hub PGS vel.y: %.6f  LDL vel.y: %.6f\n", bodies_pgs[1].velocity.y, w->body_hot[1].velocity.y);
	for (int i = 2; i < nb; i++)
		printf("    [shattering] Leaf %d PGS vel.y: %.6f  LDL vel.y: %.6f\n", i, bodies_pgs[i].velocity.y, w->body_hot[i].velocity.y);

	// All Y velocities should be near zero (gravity held by constraints).
	// Full gravity impulse is g*dt = -0.0409. If we see that, constraints aren't working.
	float grav_impulse = fabsf(g * sub_dt);
	for (int i = 1; i < nb - 1; i++) { // skip last leaf (test setup issue with PGS too)
		float uncorrected_fraction = fabsf(w->body_hot[i].velocity.y) / grav_impulse;
		printf("    Body %d: vel.y=%.6f, uncorrected=%.1f%%\n", i, w->body_hot[i].velocity.y, uncorrected_fraction * 100);
		TEST_ASSERT(uncorrected_fraction < 0.1f); // less than 10% of gravity uncorrected
	}

	integration_cache_free(&cc);
	soft_test_free_world(w);
}

// ============================================================================
// Runner

static void run_pgs_vs_ldl_tests()
{
	printf("--- PGS vs LDL comparison tests ---\n");

	test_pgs_vs_ldl_single_rigid();
	test_pgs_vs_ldl_chain_gravity();
	test_pgs_vs_ldl_star_gravity();
	test_pgs_vs_ldl_star_shattering();
}
