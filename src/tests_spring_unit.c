// tests_spring_unit.c -- unit tests for spring_compute + ldl_island_solve with soft constraints.
// Tests the full soft constraint path: spring parameters -> bias/softness -> RHS -> solve -> impulse.
// Physical correctness: verify the solver pulls stretched bodies together, not apart.

// ============================================================================
// spring_compute tests

static void test_spring_rigid()
{
	TEST_BEGIN("spring_rigid");
	// frequency = 0: rigid constraint. ptv = 0, softness = 0.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = 0, .damping_ratio = 1 }, 1.0f / 60.0f, &ptv, &soft);
	TEST_ASSERT_FLOAT(ptv, 0.0f, 1e-10f);
	TEST_ASSERT_FLOAT(soft, 0.0f, 1e-10f);
}

static void test_spring_negative_freq()
{
	TEST_BEGIN("spring_negative_freq");
	// Negative frequency: treated as rigid.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = -5, .damping_ratio = 1 }, 1.0f / 60.0f, &ptv, &soft);
	TEST_ASSERT_FLOAT(ptv, 0.0f, 1e-10f);
	TEST_ASSERT_FLOAT(soft, 0.0f, 1e-10f);
}

static void test_spring_basic()
{
	TEST_BEGIN("spring_basic");
	// 30 Hz, critically damped, dt = 1/60.
	// omega = 2*pi*30 ~= 188.5
	// d = 2 * 1.0 * 188.5 = 377
	// k = 188.5^2 = 35531
	// hd = dt * d = 6.283, hk = dt * k = 592.2, hhk = dt * hk = 9.87
	// denom = hd + hhk = 16.15
	// softness = 1/denom ~= 0.0619
	// ptv = hk / denom = 592.2/16.15 ~= 36.67
	float ptv, soft;
	float dt = 1.0f / 60.0f;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 1 }, dt, &ptv, &soft);

	TEST_ASSERT(soft > 0.01f && soft < 1.0f);
	TEST_ASSERT(ptv > 10.0f && ptv < 100.0f);
	// ptv / soft = hk, which is omega^2 * dt = large
	TEST_ASSERT(ptv / soft > 100.0f);
}

static void test_spring_high_frequency()
{
	TEST_BEGIN("spring_high_frequency");
	// 240 Hz, critically damped, dt = 1/240.
	// Very stiff spring: softness should be small, ptv should be large.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = 240, .damping_ratio = 1 }, 1.0f / 240.0f, &ptv, &soft);

	TEST_ASSERT(soft > 0.0f);
	TEST_ASSERT(ptv > 0.0f);
	// Stiffer spring = larger ptv (more aggressive correction)
	float ptv_30, soft_30;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 1 }, 1.0f / 240.0f, &ptv_30, &soft_30);
	TEST_ASSERT(ptv > ptv_30);
}

static void test_spring_overdamped()
{
	TEST_BEGIN("spring_overdamped");
	// damping_ratio = 10: heavily overdamped. More damping = more softness.
	float ptv_crit, soft_crit, ptv_over, soft_over;
	float dt = 1.0f / 60.0f;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 1 }, dt, &ptv_crit, &soft_crit);
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 10 }, dt, &ptv_over, &soft_over);

	// More damping -> larger denom -> smaller softness? No -- larger d -> larger hd -> larger denom -> smaller softness.
	// But also larger denom means ptv = hk/denom is smaller.
	TEST_ASSERT(soft_over < soft_crit); // more damping = stiffer (less compliance)
	TEST_ASSERT(ptv_over < ptv_crit);   // but also less position correction per step
}

static void test_spring_underdamped()
{
	TEST_BEGIN("spring_underdamped");
	// damping_ratio = 0.1: bouncy spring.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 0.1f }, 1.0f / 60.0f, &ptv, &soft);
	TEST_ASSERT(soft > 0.0f);
	TEST_ASSERT(ptv > 0.0f);
}

static void test_spring_tiny_dt()
{
	TEST_BEGIN("spring_tiny_dt");
	// Very small timestep (substep of substep). Should not produce NaN/inf.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 1 }, 1e-6f, &ptv, &soft);
	TEST_ASSERT(ptv == ptv); // not NaN
	TEST_ASSERT(soft == soft);
	TEST_ASSERT(soft >= 0.0f);
}

static void test_spring_large_dt()
{
	TEST_BEGIN("spring_large_dt");
	// Large timestep (low framerate). Should still produce valid values.
	float ptv, soft;
	spring_compute((SpringParams){ .frequency = 30, .damping_ratio = 1 }, 1.0f / 10.0f, &ptv, &soft);
	TEST_ASSERT(ptv == ptv);
	TEST_ASSERT(soft == soft);
	TEST_ASSERT(soft > 0.0f);
	TEST_ASSERT(ptv > 0.0f);
}

// ============================================================================
// Physical correctness: ldl_island_solve with soft ball-socket.
// Key test: a stretched spring should produce a lambda that pulls bodies together.

// Helper: set up a minimal WorldInternal with joints array for soft constraint solve.
// This is needed because ldl_island_solve reads w->joints for soft bias recomputation.
typedef struct SoftTestWorld
{
	BodyHot bodies[16];
	BodyState states[16];
	int body_count;
	JointInternal joints[16];
	int joint_count;
	SolverJoint sol_joints[16];
	int sol_joint_count;
} SoftTestWorld;

static WorldInternal* soft_test_make_world(SoftTestWorld* sw)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	afit(w->body_hot, sw->body_count);
	asetlen(w->body_hot, sw->body_count);
	for (int i = 0; i < sw->body_count; i++) w->body_hot[i] = sw->bodies[i];
	afit(w->body_state, sw->body_count);
	asetlen(w->body_state, sw->body_count);
	for (int i = 0; i < sw->body_count; i++) w->body_state[i] = sw->states[i];
	afit(w->joints, sw->joint_count);
	asetlen(w->joints, sw->joint_count);
	for (int i = 0; i < sw->joint_count; i++) w->joints[i] = sw->joints[i];
	return w;
}

static void soft_test_free_world(WorldInternal* w)
{
	afree(w->body_hot);
	afree(w->body_state);
	afree(w->joints);
	CK_FREE(w);
}

static void test_soft_spring_pulls_together()
{
	TEST_BEGIN("soft_spring_pulls_together");
	// Two bodies separated along X. Soft ball-socket spring connecting them.
	// Body A at origin, body B at (3, 0, 0). Anchors at body centers.
	// The spring should produce a lambda that pushes A toward +X and B toward -X.
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = make_body(1, 1).hot;
	sw.states[0].position = V3(0, 0, 0);
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(3, 0, 0);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,0,0), .spring = { .frequency = 30, .damping_ratio = 1 } },
	};

	float sub_dt = 1.0f / 240.0f; // 60Hz * 4 substeps
	float ptv, soft;
	spring_compute(sw.joints[0].ball_socket.spring, sub_dt, &ptv, &soft);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = V3(0,0,0), .r_b = V3(0,0,0),
		.softness = soft,
		 // error * ptv
		.joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);
	{ v3 _b = scale(sub(sw.states[1].position, sw.states[0].position), ptv); sw.sol_joints[0].bias[0] = _b.x; sw.sol_joints[0].bias[1] = _b.y; sw.sol_joints[0].bias[2] = _b.z; }

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);

	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);
	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// After solve: body A should have gained +X velocity, body B should have gained -X velocity.
	// Both moving toward each other = spring pulling.
	TEST_ASSERT(w->body_hot[0].velocity.x > 0.01f); // A moves right
	TEST_ASSERT(w->body_hot[1].velocity.x < -0.01f); // B moves left

	// Y and Z should be near zero (separation is along X only)
	TEST_ASSERT(fabsf(w->body_hot[0].velocity.y) < 0.01f);
	TEST_ASSERT(fabsf(w->body_hot[0].velocity.z) < 0.01f);

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_soft_spring_no_overshoot()
{
	TEST_BEGIN("soft_spring_no_overshoot");
	// Same setup but verify the velocity magnitude is reasonable.
	// For a critically damped spring at 30Hz, the correction should be gentle.
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = make_body(1, 1).hot;
	sw.states[0].position = V3(0, 0, 0);
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(1, 0, 0); // 1 unit separation

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,0,0), .spring = { .frequency = 30, .damping_ratio = 1 } },
	};

	float sub_dt = 1.0f / 240.0f;
	float ptv, soft;
	spring_compute(sw.joints[0].ball_socket.spring, sub_dt, &ptv, &soft);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = V3(0,0,0), .r_b = V3(0,0,0),
		.softness = soft,

		.joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);
	{ v3 _b = scale(sub(sw.states[1].position, sw.states[0].position), ptv); sw.sol_joints[0].bias[0] = _b.x; sw.sol_joints[0].bias[1] = _b.y; sw.sol_joints[0].bias[2] = _b.z; }

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);
	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// Velocity should be moderate -- not shooting off. At sub_dt = 1/240,
	// even aggressive correction shouldn't exceed a few m/s for 1m separation.
	float speed_a = len(w->body_hot[0].velocity);
	float speed_b = len(w->body_hot[1].velocity);
	TEST_ASSERT(speed_a < 50.0f);
	TEST_ASSERT(speed_b < 50.0f);
	// Both should be nonzero (spring is doing work)
	TEST_ASSERT(speed_a > 0.001f);
	TEST_ASSERT(speed_b > 0.001f);

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_soft_spring_heavy_light()
{
	TEST_BEGIN("soft_spring_heavy_light");
	// Heavy (1000) + light (1) with soft spring. Light body should move much more.
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = make_body(1000, 1000).hot;
	sw.states[0].position = V3(0, 0, 0);
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(2, 0, 0);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,0,0), .spring = { .frequency = 30, .damping_ratio = 1 } },
	};

	float sub_dt = 1.0f / 240.0f;
	float ptv, soft;
	spring_compute(sw.joints[0].ball_socket.spring, sub_dt, &ptv, &soft);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = V3(0,0,0), .r_b = V3(0,0,0),
		.softness = soft,

		.joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);
	{ v3 _b = scale(sub(sw.states[1].position, sw.states[0].position), ptv); sw.sol_joints[0].bias[0] = _b.x; sw.sol_joints[0].bias[1] = _b.y; sw.sol_joints[0].bias[2] = _b.z; }

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);
	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// Both should move toward each other
	TEST_ASSERT(w->body_hot[0].velocity.x > 0.0f);
	TEST_ASSERT(w->body_hot[1].velocity.x < 0.0f);

	// Light body should move ~1000x faster than heavy body
	float ratio = fabsf(w->body_hot[1].velocity.x) / (fabsf(w->body_hot[0].velocity.x) + 1e-10f);
	TEST_ASSERT(ratio > 100.0f);

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_rigid_constraint_zeroes_velocity()
{
	TEST_BEGIN("rigid_constraint_zeroes_velocity");
	// Rigid ball-socket (softness = 0). Both bodies separating along X.
	// The solver should produce a lambda that exactly cancels the separation velocity.
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = make_body(1, 1).hot;
	sw.states[0].position = V3(0, 0, 0);
	sw.bodies[0].velocity = V3(-1, 0, 0); // moving apart
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(1, 0, 0);
	sw.bodies[1].velocity = V3(1, 0, 0); // moving apart

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(0,0,0), .local_b = V3(0,0,0), .spring = { .frequency = 0, .damping_ratio = 0 } },
	};

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = V3(0,0,0), .r_b = V3(0,0,0),
		.softness = 0, // rigid: no position bias in velocity solve
		.joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);

	float sub_dt = 1.0f / 240.0f;

	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);
	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// After solve: constraint velocity should be near zero.
	// Bodies were separating at 2 m/s relative; the rigid constraint should cancel that.
	// With equal masses and zero lever arms, each body gets half the correction.
	// Final: both should have ~0 relative velocity.
	float rel_vx = w->body_hot[1].velocity.x - w->body_hot[0].velocity.x;
	TEST_ASSERT(fabsf(rel_vx) < 0.01f);

	integration_cache_free(&c);
	soft_test_free_world(w);
}

static void test_rigid_constraint_with_lever()
{
	TEST_BEGIN("rigid_constraint_with_lever");
	// Rigid ball-socket with lever arm. Body A rotating, creating constraint velocity.
	// Solver should cancel the constraint velocity.
	SoftTestWorld sw = {0};
	sw.body_count = 2;
	sw.bodies[0] = make_body(1, 1).hot;
	sw.states[0].position = V3(0, 0, 0);
	sw.bodies[0].angular_velocity = V3(0, 0, 2); // spinning about Z
	sw.bodies[1] = make_body(1, 1).hot;
	sw.states[1].position = V3(2, 0, 0);

	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.joints[0] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = 0, .body_b = 1,
		.ball_socket = { .local_a = V3(1,0,0), .local_b = V3(0,0,0), .spring = { 0 } },
	};

	v3 r_a = V3(1, 0, 0); // local_a rotated by identity
	sw.joint_count = 1;
	sw.sol_joint_count = 1;
	sw.sol_joints[0] = (SolverJoint){ .type = JOINT_BALL_SOCKET, .dof = 3,
		.body_a = 0, .body_b = 1,
		.r_a = r_a, .r_b = V3(0,0,0),
		.softness = 0, .joint_idx = 0,
	};
	test_fill_bs_rows(&sw.sol_joints[0]);

	float sub_dt = 1.0f / 240.0f;
	LDL_Cache c = {0};
	LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = 0, .body_b = 1, .real_body_a = 0, .real_body_b = 1, .weight_a = 1, .weight_b = 1, .solver_idx = 0 };
	apush(c.constraints, con);
	c.joint_count = 1;

	WorldInternal* w = soft_test_make_world(&sw);
	ldl_build_bundles(&c);
	ldl_build_topology(&c, w);
	ldl_numeric_factor(&c, w, sw.sol_joints, NULL);

	// Compute constraint velocity before solve
	LDL_JacobianRow jac[3];
	ldl_fill_jacobian(&c.constraints[0], sw.sol_joints, jac);
	double cv_before[3];
	for (int d = 0; d < 3; d++) cv_before[d] = ldl_constraint_velocity(&jac[d], &w->body_hot[0], &w->body_hot[1]);

	ldl_island_solve(&c, w, sw.sol_joints, sw.sol_joint_count, sub_dt);

	// Recompute constraint velocity after solve
	double cv_after[3];
	for (int d = 0; d < 3; d++) cv_after[d] = ldl_constraint_velocity(&jac[d], &w->body_hot[0], &w->body_hot[1]);

	// Constraint velocity should be much closer to zero after solve
	double cv_mag_before = fabs(cv_before[0]) + fabs(cv_before[1]) + fabs(cv_before[2]);
	double cv_mag_after = fabs(cv_after[0]) + fabs(cv_after[1]) + fabs(cv_after[2]);
	TEST_ASSERT(cv_mag_before > 0.1); // there was constraint velocity before
	TEST_ASSERT(cv_mag_after < cv_mag_before * 0.1); // reduced by at least 10x

	integration_cache_free(&c);
	soft_test_free_world(w);
}

// ============================================================================
// Runner

static void run_spring_unit_tests()
{
	printf("--- spring + soft constraint tests ---\n");

	// spring_compute
	test_spring_rigid();
	test_spring_negative_freq();
	test_spring_basic();
	test_spring_high_frequency();
	test_spring_overdamped();
	test_spring_underdamped();
	test_spring_tiny_dt();
	test_spring_large_dt();

	// Physical correctness with ldl_island_solve
	test_soft_spring_pulls_together();
	test_soft_spring_no_overshoot();
	test_soft_spring_heavy_light();
	test_rigid_constraint_zeroes_velocity();
	test_rigid_constraint_with_lever();
}
