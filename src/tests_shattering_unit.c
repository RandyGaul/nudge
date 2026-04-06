// tests_shattering_unit.c -- unit tests for ldl_apply_shattering and the full
// sparse pipeline with shattering enabled. Reproduces the hub star crash.

// ============================================================================
// ldl_apply_shattering tests

static void test_shatter_below_threshold()
{
	TEST_BEGIN("shatter_below_threshold");
	// 2 ball-sockets on hub = 6 DOF <= SHATTER_THRESHOLD (6). No shattering.
	BodyHot bodies[3];
	for (int i = 0; i < 3; i++) bodies[i] = make_body((float)(1 + i), (float)(2 + i));
	WorldInternal w = {0};
	afit(w.body_hot, 3); asetlen(w.body_hot, 3);
	for (int i = 0; i < 3; i++) w.body_hot[i] = bodies[i];

	LDL_Cache c = {0};
	for (int i = 0; i < 2; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 2;

	ldl_apply_shattering(&c, &w);

	// No shattering: all weights should still be 1, no virtual bodies
	TEST_ASSERT(c.virtual_body_count == 0);
	for (int i = 0; i < c.joint_count; i++) {
		TEST_ASSERT(c.constraints[i].weight_a == 1.0);
		TEST_ASSERT(c.constraints[i].weight_b == 1.0);
		TEST_ASSERT(!c.constraints[i].is_synthetic);
	}
	TEST_ASSERT(c.joint_count == 2); // no synthetic welds added

	afree(c.constraints); afree(c.body_remap); afree(c.shard_counts);
	afree(w.body_hot);
}

static void test_shatter_above_threshold()
{
	TEST_BEGIN("shatter_above_threshold");
	// 6 ball-sockets on hub body 0 = 18 DOF > 15. Triggers shattering.
	// S = ceil(18 / 6) = 3 shards.
	BodyHot bodies[7];
	for (int i = 0; i < 7; i++) bodies[i] = make_body((float)(1 + i), (float)(2 + i));
	WorldInternal w = {0};
	afit(w.body_hot, 7); asetlen(w.body_hot, 7);
	for (int i = 0; i < 7; i++) w.body_hot[i] = bodies[i];

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	ldl_apply_shattering(&c, &w);

	// Shattering should have occurred
	TEST_ASSERT(c.virtual_body_count > 0);
	int S = c.shard_counts[0];
	TEST_ASSERT(S >= 2);

	// Original constraints touching hub should have weight_a = S (hub is body_a)
	for (int i = 0; i < 6; i++) {
		TEST_ASSERT(c.constraints[i].weight_a == (double)S);
		TEST_ASSERT(c.constraints[i].real_body_a == 0); // real body preserved
		TEST_ASSERT(c.constraints[i].body_a >= 7); // remapped to virtual shard
	}

	// Synthetic welds should have been added
	int synth_count = 0;
	for (int i = 0; i < c.joint_count; i++)
		if (c.constraints[i].is_synthetic) synth_count++;
	TEST_ASSERT(synth_count == S); // wrap-around chain: S welds for S shards

	afree(c.constraints); afree(c.body_remap); afree(c.shard_counts);
	afree(w.body_hot);
}

static void test_shatter_static_hub_excluded()
{
	TEST_BEGIN("shatter_static_hub_excluded");
	// Static hub (inv_mass = 0): should NOT be shattered even if DOF > threshold.
	BodyHot bodies[7];
	bodies[0] = (BodyHot){0}; // static
	bodies[0].rotation = quat_identity();
	for (int i = 1; i < 7; i++) bodies[i] = make_body(1, 1);
	WorldInternal w = {0};
	afit(w.body_hot, 7); asetlen(w.body_hot, 7);
	for (int i = 0; i < 7; i++) w.body_hot[i] = bodies[i];

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	ldl_apply_shattering(&c, &w);

	TEST_ASSERT(c.virtual_body_count == 0);
	TEST_ASSERT(c.joint_count == 6); // no synthetic welds

	afree(c.constraints); afree(c.body_remap); afree(c.shard_counts);
	afree(w.body_hot);
}

// ============================================================================
// Full pipeline with shattering: reproduce hub star crash

static void test_shatter_hub_star_pipeline()
{
	TEST_BEGIN("shatter_hub_star_pipeline");
	// 6 ball-sockets on dynamic hub. Full pipeline: shattering -> bundles ->
	// topology -> numeric_factor -> solve. This is the crash scenario.
	SoftTestWorld sw = {0};
	sw.body_count = 7;
	sw.bodies[0] = make_body(5, 10); // hub
	for (int i = 1; i <= 6; i++) sw.bodies[i] = make_body(1, 1); // leaves
	sw.bs_count = 6;
	v3 dirs[6] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
	for (int i = 0; i < 6; i++) {
		sw.sol_bs[i] = (SolverBallSocket){ .r_a = dirs[i], .r_b = scale(dirs[i], -1), .body_a = 0, .body_b = i + 1, .joint_idx = i };
	}

	sw.joint_count = 6;
	for (int i = 0; i < 6; i++) {
		sw.joints[i] = (JointInternal){
			.type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = i + 1,
			.ball_socket = { .local_a = dirs[i], .local_b = scale(dirs[i], -1), .spring = { 0 } },
		};
	}

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	WorldInternal* w = soft_test_make_world(&sw);

	// This is where the crash happens: shattering changes weights and adds
	// synthetic welds, then numeric_factor builds K with those weights.
	ldl_apply_shattering(&c, w);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_bs, sw.sol_dist, sw.sol_hinge);

	// Verify all pivots are positive (the assert that was crashing)
	LDL_Topology* t = c.topo;
	for (int i = 0; i < t->node_count; i++)
		for (int d = 0; d < t->dof[i]; d++)
			TEST_ASSERT(c.diag_D[i][d] > 0);

	// If we get here without assert, try solving
	float sub_dt = 1.0f / 240.0f;
	ldl_island_solve(&c, w, sw.sol_bs, sw.bs_count, sw.sol_dist, sw.dist_count, sw.sol_hinge, sw.hinge_count, sub_dt);

	// Verify no NaN in body velocities
	for (int i = 0; i < sw.body_count; i++) {
		TEST_ASSERT(w->body_hot[i].velocity.x == w->body_hot[i].velocity.x);
		TEST_ASSERT(w->body_hot[i].velocity.y == w->body_hot[i].velocity.y);
		TEST_ASSERT(w->body_hot[i].velocity.z == w->body_hot[i].velocity.z);
	}

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_shatter_hub_star_extreme_mass()
{
	TEST_BEGIN("shatter_hub_star_extreme_mass");
	// Heavy hub (1000) + light leaves (1). Shattering weights amplify the ratio.
	SoftTestWorld sw = {0};
	sw.body_count = 7;
	sw.bodies[0] = make_body(1000, 1000); // heavy hub
	for (int i = 1; i <= 6; i++) sw.bodies[i] = make_body(1, 1);
	sw.bs_count = 6;
	v3 dirs[6] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
	for (int i = 0; i < 6; i++) {
		sw.sol_bs[i] = (SolverBallSocket){ .r_a = dirs[i], .r_b = scale(dirs[i], -1), .body_a = 0, .body_b = i + 1, .joint_idx = i };
	}
	sw.joint_count = 6;
	for (int i = 0; i < 6; i++) {
		sw.joints[i] = (JointInternal){
			.type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = i + 1,
			.ball_socket = { .local_a = dirs[i], .local_b = scale(dirs[i], -1), .spring = { 0 } },
		};
	}

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_apply_shattering(&c, w);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_bs, sw.sol_dist, sw.sol_hinge);

	// Check pivots
	LDL_Topology* t = c.topo;
	for (int i = 0; i < t->node_count; i++)
		for (int d = 0; d < t->dof[i]; d++)
			TEST_ASSERT(c.diag_D[i][d] > 0);

	float sub_dt = 1.0f / 240.0f;
	ldl_island_solve(&c, w, sw.sol_bs, sw.bs_count, sw.sol_dist, sw.dist_count, sw.sol_hinge, sw.hinge_count, sub_dt);

	for (int i = 0; i < sw.body_count; i++) {
		TEST_ASSERT(w->body_hot[i].velocity.x == w->body_hot[i].velocity.x);
	}

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_shatter_hub_with_soft_springs()
{
	TEST_BEGIN("shatter_hub_with_soft_springs");
	// Hub star with soft springs on all joints. Shattering + softness.
	SoftTestWorld sw = {0};
	sw.body_count = 7;
	sw.bodies[0] = make_body(5, 10);
	sw.bodies[0].position = V3(0, 0, 0);
	for (int i = 1; i <= 6; i++) {
		sw.bodies[i] = make_body(1, 1);
		sw.bodies[i].position = V3((float)(i % 3 - 1) * 2, (float)(i / 3) * 2, 0);
	}
	sw.bs_count = 6;
	float sub_dt = 1.0f / 240.0f;
	v3 dirs[6] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
	for (int i = 0; i < 6; i++) {
		float ptv, soft;
		SpringParams sp = { .frequency = 30, .damping_ratio = 1 };
		spring_compute(sp, sub_dt, &ptv, &soft);
		v3 err = sub(sw.bodies[i + 1].position, sw.bodies[0].position);
		sw.sol_bs[i] = (SolverBallSocket){
			.r_a = dirs[i], .r_b = scale(dirs[i], -1),
			.body_a = 0, .body_b = i + 1,
			.softness = soft, .bias = scale(err, ptv),
			.lambda = V3(0,0,0), .joint_idx = i,
		};
	}
	sw.joint_count = 6;
	for (int i = 0; i < 6; i++) {
		sw.joints[i] = (JointInternal){
			.type = JOINT_BALL_SOCKET, .body_a = 0, .body_b = i + 1,
			.ball_socket = { .local_a = dirs[i], .local_b = scale(dirs[i], -1), .spring = { .frequency = 30, .damping_ratio = 1 } },
		};
	}

	LDL_Cache c = {0};
	for (int i = 0; i < 6; i++) {
		LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = i + 1, .real_body_a = 0, .real_body_b = i + 1, .weight_a = 1, .weight_b = 1, .solver_idx = i };
		apush(c.constraints, con);
	}
	c.joint_count = 6;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_apply_shattering(&c, w);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_bs, sw.sol_dist, sw.sol_hinge);

	LDL_Topology* t = c.topo;
	for (int i = 0; i < t->node_count; i++)
		for (int d = 0; d < t->dof[i]; d++)
			TEST_ASSERT(c.diag_D[i][d] > 0);

	ldl_island_solve(&c, w, sw.sol_bs, sw.bs_count, sw.sol_dist, sw.dist_count, sw.sol_hinge, sw.hinge_count, sub_dt);

	// Velocities should be finite and bounded
	for (int i = 0; i < sw.body_count; i++) {
		float speed = len(w->body_hot[i].velocity);
		TEST_ASSERT(speed == speed); // not NaN
		TEST_ASSERT(speed < 1000.0f); // not shooting off
	}

	integration_cache_free(&c);
	soft_test_free_world(w);
}

// ============================================================================
// Runner

static void run_shattering_unit_tests()
{
	printf("--- shattering unit tests ---\n");

	test_shatter_below_threshold();
	test_shatter_above_threshold();
	test_shatter_static_hub_excluded();
	test_shatter_hub_star_pipeline();
	test_shatter_hub_star_extreme_mass();
	test_shatter_hub_with_soft_springs();
}
