// tests_heightfield_unit.c -- static heightfield shape: creation, collision,
// raycast, per-cell material palette, snapshot roundtrip.

// Build an N*N flat heightfield at height h0. All cells share the same height;
// the resulting surface is a flat plane useful for basic rest / drop tests.
static Heightfield* hf_build_flat(int N, float cell_size, float h0)
{
	float* heights = (float*)malloc(sizeof(float) * (size_t)N * (size_t)N);
	for (int k = 0; k < N * N; k++) heights[k] = h0;
	Heightfield* hf = heightfield_create(heights, N, cell_size);
	free(heights);
	return hf;
}

// Build an N*N ramp rising along +X. heights[j*N + i] = slope * i * cell_size.
static Heightfield* hf_build_ramp(int N, float cell_size, float slope)
{
	float* heights = (float*)malloc(sizeof(float) * (size_t)N * (size_t)N);
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) {
			heights[j * N + i] = slope * (float)i * cell_size;
		}
	}
	Heightfield* hf = heightfield_create(heights, N, cell_size);
	free(heights);
	return hf;
}

static void test_heightfield_create_and_free()
{
	TEST_BEGIN("heightfield_create_and_free");
	Heightfield* hf = hf_build_flat(8, 1.0f, 0.0f);
	TEST_ASSERT(hf != NULL);
	TEST_ASSERT(heightfield_tri_count(hf) == 7 * 7 * 2);
	heightfield_set_name(hf, "flat8");
	TEST_ASSERT(heightfield_get_name(hf) != NULL);
	heightfield_free(hf);
}

static void test_heightfield_sphere_rests_on_flat()
{
	TEST_BEGIN("heightfield_sphere_rests_on_flat");
	// 8x8 flat field at y=0, centered-ish on origin.
	Heightfield* hf = hf_build_flat(8, 1.0f, 0.0f);
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	// Place the heightfield so the grid center lies on the origin.
	Body floor_b = create_body(w, (BodyParams){ .position = V3(-3.5f, 0, -3.5f), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
	Body ball = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, ball);
	v3 vel = body_get_velocity(w, ball);
	printf("  heightfield sphere final y=%.4f vel=%.4f\n", p.y, v3_len(vel));
	TEST_ASSERT(p.y > 0.4f);
	TEST_ASSERT(p.y < 0.6f);
	TEST_ASSERT(v3_len(vel) < 0.1f);
	destroy_world(w);
	heightfield_free(hf);
}

static void test_heightfield_box_rests_on_flat()
{
	TEST_BEGIN("heightfield_box_rests_on_flat");
	Heightfield* hf = hf_build_flat(8, 1.0f, 0.0f);
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(-3.5f, 0, -3.5f), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
	Body box = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, box);
	v3 vel = body_get_velocity(w, box);
	printf("  heightfield box final y=%.4f vel=%.4f\n", p.y, v3_len(vel));
	TEST_ASSERT(p.y > 0.3f);
	TEST_ASSERT(p.y < 0.5f);
	TEST_ASSERT(v3_len(vel) < 0.1f);
	destroy_world(w);
	heightfield_free(hf);
}

static void test_heightfield_raycast()
{
	TEST_BEGIN("heightfield_raycast");
	// Flat field at y = 2. Ray cast from above should hit around y = 2.
	Heightfield* hf = hf_build_flat(4, 1.0f, 2.0f);
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(-1.5f, 0, -1.5f), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
	RayHit hit = {0};
	int ok = world_raycast(w, V3(0, 10, 0), V3(0, -1, 0), 50.0f, &hit);
	TEST_ASSERT(ok);
	TEST_ASSERT_FLOAT(hit.point.y, 2.0f, 0.01f);
	// Normal should point up on the flat surface.
	TEST_ASSERT(hit.normal.y > 0.95f);
	destroy_world(w);
	heightfield_free(hf);
}

static void test_heightfield_ramp_slides_ball()
{
	TEST_BEGIN("heightfield_ramp_slides_ball");
	// Ball on a ramp should end up with x > start_x (rolling downhill only
	// if slope falls toward -x). Here slope rises with +x, so ball released
	// near the top should move toward -x.
	Heightfield* hf = hf_build_ramp(8, 1.0f, 0.4f);
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(-3.5f, 0, -3.5f), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
	// Drop ball slightly up the slope.
	Body ball = create_body(w, (BodyParams){ .position = V3(1.0f, 3.0f, 0), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.1f });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.4f });
	for (int i = 0; i < 120; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, ball);
	printf("  heightfield ramp ball final x=%.4f (started at 1.0)\n", p.x);
	TEST_ASSERT(p.x < 0.9f); // has slid back toward the lower end
	destroy_world(w);
	heightfield_free(hf);
}

static void test_heightfield_per_cell_material()
{
	TEST_BEGIN("heightfield_per_cell_material");
	Heightfield* hf = hf_build_flat(4, 1.0f, 0.0f);
	// Set every cell's material to 42.
	uint8_t mats[3 * 3];
	for (int c = 0; c < 3 * 3; c++) mats[c] = 42;
	heightfield_set_material_ids(hf, mats);
	TEST_ASSERT(heightfield_get_material_id(hf, 0) == 42);
	TEST_ASSERT(heightfield_get_material_id(hf, 8) == 42);

	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(-1.5f, 0, -1.5f), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
	Body ball = create_body(w, (BodyParams){ .position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	body_set_material_id(w, ball, 99);
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n >= 1);
	int found = 0;
	for (int i = 0; i < n; i++) {
		uint32_t sub = s[i].sub_a ? s[i].sub_a : s[i].sub_b;
		if (sub == 0) continue;
		uint8_t surface_mat = s[i].sub_a ? s[i].material_a : s[i].material_b;
		uint8_t other_mat   = s[i].sub_a ? s[i].material_b : s[i].material_a;
		TEST_ASSERT(surface_mat == 42);
		TEST_ASSERT(other_mat == 99);
		found = 1;
	}
	TEST_ASSERT(found == 1);
	destroy_world(w);
	heightfield_free(hf);
}

static void test_heightfield_snapshot_roundtrip()
{
	TEST_BEGIN("heightfield_snapshot_roundtrip");
	const char* path = "test_hf_snapshot.dat";
	Heightfield* hf = hf_build_ramp(6, 1.0f, 0.2f);
	heightfield_set_name(hf, "ramp6");
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
		world_register_heightfield(w, hf);
		Body floor_b = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_HEIGHTFIELD, .heightfield.hf = hf });
		TEST_ASSERT(world_save_snapshot(w, path));
		destroy_world(w);
	}
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
		world_register_heightfield(w, hf);
		TEST_ASSERT(world_load_snapshot_into(w, path));
		Body bodies[4]; int nb = world_get_bodies(w, bodies, 4);
		TEST_ASSERT(nb == 1);
		destroy_world(w);
	}
	heightfield_free(hf);
	remove(path);
}

static void run_heightfield_unit_tests()
{
	printf("--- heightfield unit tests ---\n");
	test_heightfield_create_and_free();
	test_heightfield_sphere_rests_on_flat();
	test_heightfield_box_rests_on_flat();
	test_heightfield_raycast();
	test_heightfield_ramp_slides_ball();
	test_heightfield_per_cell_material();
	test_heightfield_snapshot_roundtrip();
}
