// tests_multishape_unit.c -- multi-shape body (local_rot + composite inertia).
// Covers: local_rot defaulting, rotated-child inertia composition, multi-shape
// collision dispatch emitting per-child manifolds.

// Compute body inertia via internal path. Returns v3 of body_inv_inertia_local
// after body_add_shape calls. 1 / component = actual diagonal inertia.
static v3 multishape_get_inv_inertia_local(World w_handle, Body body)
{
	WorldInternal* w = (WorldInternal*)w_handle.id;
	int idx = handle_index(body);
	return body_inv_inertia_local(w, idx);
}

static void test_multishape_identity_rot_matches_baseline()
{
	TEST_BEGIN("multishape_identity_rot_single_box");
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	Body b = create_body(w, (BodyParams){ .mass = 6.0f, .rotation = quat_identity() });
	// Local_rot omitted from designated init -> zero quat -> normalised to identity.
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(1, 2, 3) });
	v3 inv_i = multishape_get_inv_inertia_local(w, b);
	// Expected diagonal inertia for h=(1,2,3), m=6: (26, 20, 10).
	TEST_ASSERT_FLOAT(1.0f / inv_i.x, 26.0f, 1e-4f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.y, 20.0f, 1e-4f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.z, 10.0f, 1e-4f);
	destroy_world(w);
}

static void test_multishape_90deg_about_x_swaps_yz()
{
	TEST_BEGIN("multishape_90deg_about_x_swaps_yz");
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	Body b = create_body(w, (BodyParams){ .mass = 6.0f, .rotation = quat_identity() });
	// 90 deg about X: y-axis <-> z-axis. Box local inertia (26, 20, 10) becomes
	// body-frame diagonal (26, 10, 20).
	quat r = quat_axis_angle(V3(1, 0, 0), 3.14159265f * 0.5f);
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .local_rot = r, .box.half_extents = V3(1, 2, 3) });
	v3 inv_i = multishape_get_inv_inertia_local(w, b);
	TEST_ASSERT_FLOAT(1.0f / inv_i.x, 26.0f, 1e-3f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.y, 10.0f, 1e-3f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.z, 20.0f, 1e-3f);
	destroy_world(w);
}

static void test_multishape_two_spheres_parallel_axis()
{
	TEST_BEGIN("multishape_two_spheres_parallel_axis");
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	Body b = create_body(w, (BodyParams){ .mass = 2.0f, .rotation = quat_identity() });
	// Two unit spheres at (-1,0,0) and (+1,0,0). Equal volume -> m_per = 1 each.
	// Each: local I diag = 0.4. After parallel axis with d=(+/-1, 0, 0):
	//   I_xx = 0.4 + 0 = 0.4
	//   I_yy = 0.4 + 1 = 1.4
	//   I_zz = 0.4 + 1 = 1.4
	// Two spheres -> (0.8, 2.8, 2.8).
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .local_pos = V3(-1, 0, 0), .sphere.radius = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .local_pos = V3( 1, 0, 0), .sphere.radius = 1.0f });
	v3 inv_i = multishape_get_inv_inertia_local(w, b);
	TEST_ASSERT_FLOAT(1.0f / inv_i.x, 0.8f, 1e-3f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.y, 2.8f, 1e-3f);
	TEST_ASSERT_FLOAT(1.0f / inv_i.z, 2.8f, 1e-3f);
	destroy_world(w);
}

static void test_multishape_collides_per_child()
{
	TEST_BEGIN("multishape_collides_per_child");
	// Body A: two small boxes at (-2, 0, 0) and (+2, 0, 0) -- local shape.
	// Body B: static floor. Start penetrating so broadphase definitely matches.
	// After step, confirm debug contact list has contacts from both child pairs
	// (at least 2 distinct contact points with x<0 and x>0).
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body a = create_body(w, (BodyParams){ .mass = 2.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_BOX, .local_pos = V3(-2, 0, 0), .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_BOX, .local_pos = V3( 2, 0, 0), .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	// world_get_contacts must be called once before the step to allocate the
	// debug contact array (nudge.c populates it only when already non-NULL).
	const Contact* contacts = NULL;
	(void)world_get_contacts(w, &contacts);
	world_step(w, 1.0f / 60.0f);
	int count = world_get_contacts(w, &contacts);
	int has_neg = 0, has_pos = 0;
	for (int i = 0; i < count; i++) {
		if (contacts[i].point.x < -1.0f) has_neg = 1;
		if (contacts[i].point.x >  1.0f) has_pos = 1;
	}
	TEST_ASSERT(count >= 2);
	TEST_ASSERT(has_neg);
	TEST_ASSERT(has_pos);
	destroy_world(w);
}

static void run_multishape_unit_tests()
{
	printf("--- multi-shape unit tests ---\n");
	test_multishape_identity_rot_matches_baseline();
	test_multishape_90deg_about_x_swaps_yz();
	test_multishape_two_spheres_parallel_axis();
	test_multishape_collides_per_child();
}
