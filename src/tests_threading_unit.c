// tests_threading_unit.c -- cross-thread-count determinism gate.
//
// Per-island parallel dispatch (phase 2) + big-island per-color dispatch
// (phase 3) must produce bit-identical body_hot / body_state across thread
// counts. This is a hard requirement: islands are disjoint by construction
// so dispatch order cannot affect per-body update sequences, and within a
// color graph coloring guarantees disjoint body writes. If this test fails,
// some path has introduced cross-island or cross-color body writes.

static void threading_build_scene(World w)
{
	WorldInternal* wi = (WorldInternal*)w.id;
	// Keep sleep off so every body stays active every step -- maximizes the
	// parallel dispatch path we want to exercise.
	wi->sleep_enabled = 0;

	// Floor.
	Body floor_body = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_body, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(20, 0.5f, 20) });

	// 5x5x3 pile -> one big contact island. Stressed enough that phase-3
	// per-color dispatch fires at n_workers >= 2.
	float half = 0.5f;
	for (int row = 0; row < 5; row++) {
		for (int col = 0; col < 5; col++) {
			float x = (float)(col - 2) * (2.0f * half + 0.01f);
			float z = (float)(row - 2) * (2.0f * half + 0.01f);
			for (int h = 0; h < 3; h++) {
				float y = 1.0f + half + (float)h * (2.0f * half + 0.01f);
				Body b = create_body(w, (BodyParams){ .position = V3(x, y, z), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(half, half, half) });
			}
		}
	}

	// Four pendulum chains far from the pile -- four small joint-linked
	// islands, exercise the whole-island parallel claim path.
	for (int chain = 0; chain < 4; chain++) {
		float cx = 8.0f + (float)chain * 2.0f;
		Body anchor = create_body(w, (BodyParams){ .position = V3(cx, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body prev = anchor;
		for (int i = 0; i < 4; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(cx, 10 - (i + 1) * 0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0, -0.4f, 0), .local_offset_b = V3(0, 0.4f, 0) });
			prev = b;
		}
	}
}

static void threading_run_and_snapshot(int thread_count, int frames, BodyHot** out_hot, BodyState** out_state, int* out_count)
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });
	WorldInternal* wi = (WorldInternal*)w.id;
	threading_build_scene(w);
	wi->thread_count = thread_count;

	for (int f = 0; f < frames; f++) world_step(w, 1.0f / 60.0f);

	int n = asize(wi->body_hot);
	BodyHot* hot = (BodyHot*)CK_ALLOC(n * sizeof(BodyHot));
	BodyState* state = (BodyState*)CK_ALLOC(n * sizeof(BodyState));
	memcpy(hot, wi->body_hot, n * sizeof(BodyHot));
	memcpy(state, wi->body_state, n * sizeof(BodyState));

	destroy_world(w);

	*out_hot = hot;
	*out_state = state;
	*out_count = n;
}

static int threading_compare_snapshots(BodyHot* a_hot, BodyState* a_state, BodyHot* b_hot, BodyState* b_state, int n, int threads_a, int threads_b)
{
	int mismatches = 0;
	for (int i = 0; i < n; i++) {
		if (memcmp(&a_hot[i], &b_hot[i], sizeof(BodyHot)) != 0) {
			mismatches++;
			if (mismatches <= 3) {
				printf("    body[%d] hot differs between threads=%d and threads=%d:\n", i, threads_a, threads_b);
				printf("      vel a=(%.9f,%.9f,%.9f) b=(%.9f,%.9f,%.9f)\n",
					a_hot[i].velocity.x, a_hot[i].velocity.y, a_hot[i].velocity.z,
					b_hot[i].velocity.x, b_hot[i].velocity.y, b_hot[i].velocity.z);
				printf("      ang a=(%.9f,%.9f,%.9f) b=(%.9f,%.9f,%.9f)\n",
					a_hot[i].angular_velocity.x, a_hot[i].angular_velocity.y, a_hot[i].angular_velocity.z,
					b_hot[i].angular_velocity.x, b_hot[i].angular_velocity.y, b_hot[i].angular_velocity.z);
			}
		}
		if (memcmp(&a_state[i], &b_state[i], sizeof(BodyState)) != 0) {
			mismatches++;
			if (mismatches <= 3) {
				printf("    body[%d] state differs between threads=%d and threads=%d:\n", i, threads_a, threads_b);
				printf("      pos a=(%.9f,%.9f,%.9f) b=(%.9f,%.9f,%.9f)\n",
					a_state[i].position.x, a_state[i].position.y, a_state[i].position.z,
					b_state[i].position.x, b_state[i].position.y, b_state[i].position.z);
			}
		}
	}
	return mismatches;
}

static void test_threading_determinism_fixed_count(int threads_a, int threads_b, int frames)
{
	char name[64];
	snprintf(name, sizeof(name), "threading_determinism_%d_vs_%d", threads_a, threads_b);
	TEST_BEGIN(name);

	BodyHot *hot_a = NULL, *hot_b = NULL;
	BodyState *state_a = NULL, *state_b = NULL;
	int n_a = 0, n_b = 0;

	threading_run_and_snapshot(threads_a, frames, &hot_a, &state_a, &n_a);
	threading_run_and_snapshot(threads_b, frames, &hot_b, &state_b, &n_b);

	TEST_ASSERT(n_a == n_b);
	if (n_a == n_b) {
		int mismatches = threading_compare_snapshots(hot_a, state_a, hot_b, state_b, n_a, threads_a, threads_b);
		TEST_ASSERT(mismatches == 0);
		if (mismatches == 0) {
			printf("  [threading-det] threads=%d vs threads=%d: %d bodies bit-identical after %d frames\n", threads_a, threads_b, n_a, frames);
		} else {
			printf("  [threading-det] threads=%d vs threads=%d: %d / %d bodies DIVERGED after %d frames\n", threads_a, threads_b, mismatches, n_a * 2, frames);
		}
	}

	CK_FREE(hot_a); CK_FREE(hot_b);
	CK_FREE(state_a); CK_FREE(state_b);
}

static void run_threading_unit_tests()
{
	printf("--- threading determinism tests ---\n");
	// Cover: frames before contacts, frames during first contact, steady state.
	// Cover: single-threaded vs multi-threaded, plus multi-vs-multi.
	test_threading_determinism_fixed_count(1, 2, 60);
	test_threading_determinism_fixed_count(1, 4, 60);
	test_threading_determinism_fixed_count(1, 8, 60);
	test_threading_determinism_fixed_count(2, 4, 60);
	test_threading_determinism_fixed_count(4, 8, 60);
}
