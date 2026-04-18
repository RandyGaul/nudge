// tests_contacts_unit.c -- ContactSummary build + sort + dedupe + per-body listener dispatch.

typedef struct ContactsTestCapture
{
	int fires;
	int total_pairs;
	Body last_self;
	ContactSummary last_pair;
} ContactsTestCapture;

static void contacts_test_cb(Body self, const ContactSummary* pairs, int count, void* ud)
{
	ContactsTestCapture* cap = (ContactsTestCapture*)ud;
	cap->fires++;
	cap->total_pairs += count;
	cap->last_self = self;
	if (count > 0) cap->last_pair = pairs[0];
}

static void test_contacts_single_pair_basic()
{
	TEST_BEGIN("contacts_single_pair_basic");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n == 1);
	TEST_ASSERT(s[0].a.id < s[0].b.id);               // canonical order
	TEST_ASSERT(s[0].radius >= 0.0f);
	TEST_ASSERT(s[0].depth >= 0.0f);
	// Contact normal for a ball resting on floor points from 'a' to 'b'; whichever
	// body got the smaller handle determines the sign. Just check it's unit-ish.
	float nlen = v3_len(s[0].normal);
	TEST_ASSERT(nlen > 0.9f && nlen < 1.1f);
	destroy_world(w);
}

static void test_contacts_sorted_by_body_id()
{
	TEST_BEGIN("contacts_sorted_by_body_id");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	for (int i = 0; i < 5; i++) {
		Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3((float)i * 1.5f - 3, 0.4f, 0), .rotation = quat_identity() });
		body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	}
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n == 5);
	for (int i = 0; i < n; i++) TEST_ASSERT(s[i].a.id < s[i].b.id);
	for (int i = 1; i < n; i++) {
		TEST_ASSERT(s[i-1].a.id < s[i].a.id || (s[i-1].a.id == s[i].a.id && s[i-1].b.id < s[i].b.id));
	}
	destroy_world(w);
}

static void test_contacts_listener_fires_normalized()
{
	TEST_BEGIN("contacts_listener_fires_normalized");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	ContactsTestCapture cap = {0};
	body_set_contact_listener(w, ball, contacts_test_cb, &cap);
	world_step(w, 1.0f / 60.0f);
	TEST_ASSERT(cap.fires == 1);
	TEST_ASSERT(cap.total_pairs == 1);
	TEST_ASSERT(cap.last_self.id == ball.id);
	TEST_ASSERT(cap.last_pair.a.id == ball.id);       // normalized into `a` slot
	TEST_ASSERT(cap.last_pair.b.id == floor.id);
	// Normal should point from ball toward floor = roughly (0, -1, 0).
	TEST_ASSERT(cap.last_pair.normal.y < -0.8f);
	destroy_world(w);
}

static void test_contacts_listener_clear()
{
	TEST_BEGIN("contacts_listener_clear");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	ContactsTestCapture cap = {0};
	body_set_contact_listener(w, ball, contacts_test_cb, &cap);
	world_step(w, 1.0f / 60.0f);
	TEST_ASSERT(cap.fires == 1);
	body_set_contact_listener(w, ball, NULL, NULL);
	world_step(w, 1.0f / 60.0f);
	TEST_ASSERT(cap.fires == 1); // still 1, listener cleared
	destroy_world(w);
}

static void test_contacts_no_collision_zero_summaries()
{
	TEST_BEGIN("contacts_no_collision_zero_summaries");
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	Body a = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(-5, 0, 0), .rotation = quat_identity() });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	Body b = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3( 5, 0, 0), .rotation = quat_identity() });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	world_step(w, 1.0f / 60.0f);
	int n = -1;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n == 0);
	(void)s;
	destroy_world(w);
}

static void test_contacts_destroy_body_clears_listener()
{
	TEST_BEGIN("contacts_destroy_body_clears_listener");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	ContactsTestCapture cap = {0};
	body_set_contact_listener(w, ball, contacts_test_cb, &cap);
	destroy_body(w, ball);
	// Next step must not fire the (freed) callback.
	world_step(w, 1.0f / 60.0f);
	TEST_ASSERT(cap.fires == 0);
	destroy_world(w);
}

static void test_contacts_body_material_id_flows_through()
{
	TEST_BEGIN("contacts_body_material_id_flows_through");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	world_set_material(w, 7, (Material){ .friction = 0.9f, .restitution = 0.1f, .user_data = 0xABCD });
	world_set_material(w, 3, (Material){ .friction = 0.2f, .restitution = 0.5f, .user_data = 0x1234 });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	body_set_material_id(w, floor, 7);
	Body ball = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	body_set_material_id(w, ball, 3);
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n == 1);
	// Summary is canonical (a.id < b.id). Whichever body got the lower handle
	// wears its own material on the `a` side.
	uint8_t exp_a = s[0].a.id == ball.id ? 3 : 7;
	uint8_t exp_b = s[0].b.id == ball.id ? 3 : 7;
	TEST_ASSERT(s[0].material_a == exp_a);
	TEST_ASSERT(s[0].material_b == exp_b);
	// Palette readback.
	Material m = world_get_material(w, 7);
	TEST_ASSERT_FLOAT(m.friction, 0.9f, 1e-6f);
	TEST_ASSERT(m.user_data == 0xABCD);
	TEST_ASSERT(body_get_material_id(w, floor) == 7);
	destroy_world(w);
}

static void test_contacts_trimesh_per_tri_material()
{
	TEST_BEGIN("contacts_trimesh_per_tri_material");
	// Quad floor = 2 triangles; assign distinct material ids per triangle.
	// A sphere dropped on the diagonal should report a mesh-side material
	// matching one of those ids (whichever triangle its contact landed on).
	v3 verts[4] = { V3(-2, 0, -2), V3(2, 0, -2), V3(2, 0, 2), V3(-2, 0, 2) };
	uint32_t indices[6] = { 0, 2, 1, 0, 3, 2 };
	TriMesh* mesh = trimesh_create(verts, 4, indices, 2);
	uint8_t tri_mats[2] = { 11, 22 };
	trimesh_set_material_ids(mesh, tri_mats);
	TEST_ASSERT(trimesh_get_material_id(mesh, 0) == 11);
	TEST_ASSERT(trimesh_get_material_id(mesh, 1) == 22);
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body sphere_b = create_body(w, (BodyParams){ .position = V3(0.5f, 0.4f, 0.5f), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, sphere_b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	body_set_material_id(w, sphere_b, 99);
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n >= 1);
	// Find the summary whose sub index resolves to the mesh side. Canonical
	// order: `a.id < b.id`, so the mesh side is whichever has the nonzero
	// sub index.
	int found = 0;
	for (int i = 0; i < n; i++) {
		uint32_t sub = s[i].sub_a ? s[i].sub_a : s[i].sub_b;
		if (sub == 0) continue;
		uint8_t mesh_mat = s[i].sub_a ? s[i].material_a : s[i].material_b;
		uint8_t other_mat = s[i].sub_a ? s[i].material_b : s[i].material_a;
		TEST_ASSERT(mesh_mat == 11 || mesh_mat == 22);
		TEST_ASSERT(other_mat == 99);
		found = 1;
	}
	TEST_ASSERT(found == 1);
	destroy_world(w);
	trimesh_free(mesh);
}

static void test_contacts_multishape_dedupes()
{
	TEST_BEGIN("contacts_multishape_dedupes");
	// Two small boxes on one body vs floor. Narrowphase emits two manifolds
	// (one per child pair); ContactSummary must coalesce to exactly one row.
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body a = create_body(w, (BodyParams){ .mass = 2.0f, .position = V3(0, 0.4f, 0), .rotation = quat_identity() });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_BOX, .local_pos = V3(-2, 0, 0), .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_BOX, .local_pos = V3( 2, 0, 0), .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	Body floor = create_body(w, (BodyParams){ .mass = 0.0f, .position = V3(0, -0.5f, 0), .rotation = quat_identity() });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	world_step(w, 1.0f / 60.0f);
	int n = 0;
	const ContactSummary* s = world_contact_summaries(w, &n);
	TEST_ASSERT(n == 1);
	(void)s;
	destroy_world(w);
}

static void test_contacts_materials_snapshot_roundtrip()
{
	TEST_BEGIN("contacts_materials_snapshot_roundtrip");
	const char* path = "test_materials_snapshot.dat";
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		world_set_material(w, 5, (Material){ .friction = 0.77f, .restitution = 0.22f, .user_data = 0xDEADBEEF });
		world_set_material(w, 100, (Material){ .friction = 0.10f, .restitution = 0.90f, .user_data = 0xCAFEF00D });
		Body b = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 1, 0), .rotation = quat_identity() });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
		body_set_material_id(w, b, 5);
		TEST_ASSERT(world_save_snapshot(w, path));
		destroy_world(w);
	}
	{
		World w = world_load_snapshot(path);
		TEST_ASSERT(w.id != 0);
		Material m5 = world_get_material(w, 5);
		TEST_ASSERT_FLOAT(m5.friction, 0.77f, 1e-6f);
		TEST_ASSERT_FLOAT(m5.restitution, 0.22f, 1e-6f);
		TEST_ASSERT(m5.user_data == 0xDEADBEEF);
		Material m100 = world_get_material(w, 100);
		TEST_ASSERT(m100.user_data == 0xCAFEF00D);
		Body bodies[4]; int nb = world_get_bodies(w, bodies, 4);
		TEST_ASSERT(nb == 1);
		TEST_ASSERT(body_get_material_id(w, bodies[0]) == 5);
		destroy_world(w);
	}
	remove(path);
}

static void test_contacts_materials_rewind_roundtrip()
{
	TEST_BEGIN("contacts_materials_rewind_roundtrip");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	world_rewind_init(w, (RewindParams){ .max_frames = 8, .auto_capture = 1 });
	Body b = create_body(w, (BodyParams){ .mass = 1.0f, .position = V3(0, 1, 0), .rotation = quat_identity() });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	world_set_material(w, 9, (Material){ .friction = 0.5f, .restitution = 0.1f, .user_data = 111 });
	body_set_material_id(w, b, 9);
	world_step(w, 1.0f / 60.0f);                 // captures the "before" state
	uint64_t anchor_id = world_rewind_capture(w); // manual capture of known state

	world_set_material(w, 9, (Material){ .friction = 0.1f, .restitution = 0.9f, .user_data = 999 });
	body_set_material_id(w, b, 200);
	world_step(w, 1.0f / 60.0f);

	Material cur = world_get_material(w, 9);
	TEST_ASSERT(cur.user_data == 999);
	TEST_ASSERT(body_get_material_id(w, b) == 200);

	TEST_ASSERT(world_rewind_to_frame(w, anchor_id));
	Material restored = world_get_material(w, 9);
	TEST_ASSERT_FLOAT(restored.friction, 0.5f, 1e-6f);
	TEST_ASSERT(restored.user_data == 111);
	TEST_ASSERT(body_get_material_id(w, b) == 9);
	destroy_world(w);
}

static void run_contacts_unit_tests()
{
	printf("--- contacts unit tests ---\n");
	test_contacts_single_pair_basic();
	test_contacts_sorted_by_body_id();
	test_contacts_listener_fires_normalized();
	test_contacts_listener_clear();
	test_contacts_no_collision_zero_summaries();
	test_contacts_destroy_body_clears_listener();
	test_contacts_body_material_id_flows_through();
	test_contacts_trimesh_per_tri_material();
	test_contacts_multishape_dedupes();
	test_contacts_materials_snapshot_roundtrip();
	test_contacts_materials_rewind_roundtrip();
}
