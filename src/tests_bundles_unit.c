// tests_bundles_unit.c -- unit tests for ldl_build_bundles.
// Verifies: sorting by body pair, grouping into bundles, DOF accumulation,
// bundle splitting at 6 DOF, and bundle_idx/bundle_offset assignment.

// Helper: make a minimal LDL_Cache with a given constraint list.
// Caller must free c->constraints and c->bundles after test.
static void bundles_setup(LDL_Cache* c, LDL_Constraint* cons, int count)
{
	memset(c, 0, sizeof(*c));
	c->joint_count = count;
	c->constraints = NULL;
	c->bundles = NULL;
	afit(c->constraints, count);
	asetlen(c->constraints, count);
	for (int i = 0; i < count; i++) c->constraints[i] = cons[i];
}

static void bundles_teardown(LDL_Cache* c)
{
	afree(c->constraints);
	afree(c->bundles);
}

// ============================================================================
// Tests

static void test_bundles_single_constraint()
{
	TEST_BEGIN("bundles_single_constraint");
	// One ball-socket: bodies (0, 1), 3 DOF. Should produce 1 bundle.
	LDL_Constraint cons[1] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 1);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 1);
	TEST_ASSERT(c.bundles[0].body_a == 0);
	TEST_ASSERT(c.bundles[0].body_b == 1);
	TEST_ASSERT(c.bundles[0].dof == 3);
	TEST_ASSERT(c.bundles[0].start == 0);
	TEST_ASSERT(c.bundles[0].count == 1);
	TEST_ASSERT(c.constraints[0].bundle_idx == 0);
	TEST_ASSERT(c.constraints[0].bundle_offset == 0);

	bundles_teardown(&c);
}

static void test_bundles_two_different_pairs()
{
	TEST_BEGIN("bundles_two_different_pairs");
	// Two constraints on different body pairs. Should produce 2 bundles.
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 1, .body_b = 2 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 2);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 2);
	TEST_ASSERT(c.bundles[0].count == 1);
	TEST_ASSERT(c.bundles[1].count == 1);
	TEST_ASSERT(c.constraints[0].bundle_idx == 0);
	TEST_ASSERT(c.constraints[1].bundle_idx == 1);

	bundles_teardown(&c);
}

static void test_bundles_same_pair_merged()
{
	TEST_BEGIN("bundles_same_pair_merged");
	// Two constraints on the same body pair (ball-socket + distance = 3+1 = 4 DOF).
	// Should merge into one bundle.
	LDL_Constraint cons[2] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 2);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 1);
	TEST_ASSERT(c.bundles[0].dof == 4);
	TEST_ASSERT(c.bundles[0].count == 2);
	TEST_ASSERT(c.constraints[0].bundle_idx == 0);
	TEST_ASSERT(c.constraints[0].bundle_offset == 0);
	TEST_ASSERT(c.constraints[1].bundle_idx == 0);
	TEST_ASSERT(c.constraints[1].bundle_offset == 3); // after 3-DOF ball-socket

	bundles_teardown(&c);
}

static void test_bundles_split_at_6_dof()
{
	TEST_BEGIN("bundles_split_at_6_dof");
	// Same pair: ball-socket (3) + ball-socket (3) + distance (1) = 7 DOF.
	// First bundle takes 3+3=6, second bundle takes 1.
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 3);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 2);
	TEST_ASSERT(c.bundles[0].dof == 6);
	TEST_ASSERT(c.bundles[0].count == 2);
	TEST_ASSERT(c.bundles[1].dof == 1);
	TEST_ASSERT(c.bundles[1].count == 1);

	// Offsets within first bundle
	TEST_ASSERT(c.constraints[0].bundle_offset == 0);
	TEST_ASSERT(c.constraints[1].bundle_offset == 3);

	// Second bundle
	TEST_ASSERT(c.constraints[2].bundle_idx == 1);
	TEST_ASSERT(c.constraints[2].bundle_offset == 0);

	bundles_teardown(&c);
}

static void test_bundles_hinge_alone()
{
	TEST_BEGIN("bundles_hinge_alone");
	// Single hinge: 5 DOF. Fits in one bundle (< 6).
	LDL_Constraint cons[1] = {
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 2, .body_b = 3 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 1);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 1);
	TEST_ASSERT(c.bundles[0].dof == 5);

	bundles_teardown(&c);
}

static void test_bundles_hinge_plus_distance_split()
{
	TEST_BEGIN("bundles_hinge_plus_distance_split");
	// Same pair: hinge (5) + distance (1) = 6 DOF. Should fit in one bundle.
	LDL_Constraint cons[2] = {
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 2);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 1);
	TEST_ASSERT(c.bundles[0].dof == 6);
	TEST_ASSERT(c.bundles[0].count == 2);
	TEST_ASSERT(c.constraints[1].bundle_offset == 5);

	bundles_teardown(&c);
}

static void test_bundles_hinge_plus_ball_socket_split()
{
	TEST_BEGIN("bundles_hinge_plus_ball_socket_split");
	// Same pair: hinge (5) + ball-socket (3) = 8 DOF > 6. Must split.
	// First bundle takes hinge (5), second takes ball-socket (3).
	LDL_Constraint cons[2] = {
		{ .type = JOINT_HINGE, .dof = 5, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 2);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 2);
	TEST_ASSERT(c.bundles[0].dof == 5);
	TEST_ASSERT(c.bundles[0].count == 1);
	TEST_ASSERT(c.bundles[1].dof == 3);
	TEST_ASSERT(c.bundles[1].count == 1);
	TEST_ASSERT(c.constraints[0].bundle_idx == 0);
	TEST_ASSERT(c.constraints[1].bundle_idx == 1);

	bundles_teardown(&c);
}

static void test_bundles_sorting()
{
	TEST_BEGIN("bundles_sorting");
	// Constraints given in reverse body-pair order. Build_bundles should sort
	// them by (body_a, body_b) before grouping.
	LDL_Constraint cons[3] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 2, .body_b = 3 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 3);
	ldl_build_bundles(&c);

	// After sorting: (0,1) constraints first, then (2,3).
	TEST_ASSERT(c.bundle_count == 2);
	TEST_ASSERT(c.bundles[0].body_a == 0 && c.bundles[0].body_b == 1);
	TEST_ASSERT(c.bundles[0].dof == 4); // 3 + 1 merged
	TEST_ASSERT(c.bundles[0].count == 2);
	TEST_ASSERT(c.bundles[1].body_a == 2 && c.bundles[1].body_b == 3);
	TEST_ASSERT(c.bundles[1].dof == 3);

	bundles_teardown(&c);
}

static void test_bundles_many_pairs()
{
	TEST_BEGIN("bundles_many_pairs");
	// 5 constraints on 4 different body pairs, one pair has 2 constraints.
	LDL_Constraint cons[5] = {
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 3, .body_b = 4 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 1, .body_b = 2 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 1, .body_b = 2 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 2, .body_b = 3 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 5);
	ldl_build_bundles(&c);

	// Sorted pairs: (0,1), (1,2), (1,2), (2,3), (3,4)
	// Bundles: (0,1) 3 DOF, (1,2) 2 DOF (1+1 merged), (2,3) 3 DOF, (3,4) 3 DOF
	TEST_ASSERT(c.bundle_count == 4);
	TEST_ASSERT(c.bundles[0].body_a == 0 && c.bundles[0].body_b == 1);
	TEST_ASSERT(c.bundles[1].body_a == 1 && c.bundles[1].body_b == 2);
	TEST_ASSERT(c.bundles[1].dof == 2);
	TEST_ASSERT(c.bundles[1].count == 2);
	TEST_ASSERT(c.bundles[2].body_a == 2 && c.bundles[2].body_b == 3);
	TEST_ASSERT(c.bundles[3].body_a == 3 && c.bundles[3].body_b == 4);

	bundles_teardown(&c);
}

static void test_bundles_empty()
{
	TEST_BEGIN("bundles_empty");
	// Zero constraints. Should produce zero bundles.
	LDL_Cache c;
	memset(&c, 0, sizeof(c));
	c.joint_count = 0;
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 0);
	TEST_ASSERT(c.bundles == NULL || asize(c.bundles) == 0);

	afree(c.bundles);
}

static void test_bundles_offset_accumulation()
{
	TEST_BEGIN("bundles_offset_accumulation");
	// Three constraints on same pair: dist(1) + dist(1) + bs(3) = 5 DOF.
	// All fit in one bundle. Offsets: 0, 1, 2.
	LDL_Constraint cons[3] = {
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_DISTANCE, .dof = 1, .body_a = 0, .body_b = 1 },
		{ .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1 },
	};
	LDL_Cache c;
	bundles_setup(&c, cons, 3);
	ldl_build_bundles(&c);

	TEST_ASSERT(c.bundle_count == 1);
	TEST_ASSERT(c.bundles[0].dof == 5);
	TEST_ASSERT(c.constraints[0].bundle_offset == 0);
	TEST_ASSERT(c.constraints[1].bundle_offset == 1);
	TEST_ASSERT(c.constraints[2].bundle_offset == 2);

	bundles_teardown(&c);
}

// ============================================================================
// Runner

static void run_bundles_unit_tests()
{
	printf("--- bundles unit tests ---\n");

	test_bundles_single_constraint();
	test_bundles_two_different_pairs();
	test_bundles_same_pair_merged();
	test_bundles_split_at_6_dof();
	test_bundles_hinge_alone();
	test_bundles_hinge_plus_distance_split();
	test_bundles_hinge_plus_ball_socket_split();
	test_bundles_sorting();
	test_bundles_many_pairs();
	test_bundles_empty();
	test_bundles_offset_accumulation();
}
