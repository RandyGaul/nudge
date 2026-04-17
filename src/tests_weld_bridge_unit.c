// tests_weld_bridge_unit.c -- regression test for the Weld Bridge scene.
//
// Replicates the core topology of scenes.c::scene_weld_bridge_setup: a chain
// of boxes connected by fixed joints with static bodies at both endpoints.
// No spheres dropped; gravity only. Measures how far the middle bodies drift
// from their initial positions over N frames.
//
// This is a regression gate: a rigid bridge with both ends fixed should
// hold position indefinitely under gravity (the chain is over-constrained
// but structurally sound). Drift > ~0.2 m in the first second means the
// LDL / PGS solver isn't maintaining the joint constraint.

static void test_weld_bridge_drift(int use_ldl)
{
	char name[64];
	snprintf(name, sizeof(name), "weld_bridge_drift_%s", use_ldl ? "LDL" : "PGS");
	TEST_BEGIN(name);

	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = use_ldl;
	wi->sleep_enabled = 0; // keep bridge active so we can measure drift continuously

	// Bridge topology matching scene_weld_bridge_setup: n=8, link_len=1.0, half=0.4,
	// endpoints static (mass=0), interior mass=2.
	int n = 8;
	float link_len = 1.0f;
	float half = 0.4f;
	float start_x = -(n - 1) * link_len * 0.5f;

	CK_DYNA Body* bodies = NULL;
	Body prev = (Body){0};
	for (int i = 0; i < n; i++) {
		float x = start_x + i * link_len;
		float mass = (i == 0 || i == n - 1) ? 0.0f : 2.0f;
		Body b = create_body(w, (BodyParams){
			.position = V3(x, 4.0f, 0),
			.rotation = quat_identity(),
			.mass = mass,
		});
		body_add_shape(w, b, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(half, 0.15f, half),
		});
		if (i > 0) {
			create_fixed(w, (FixedParams){
				.body_a = prev, .body_b = b,
				.local_offset_a = V3(link_len * 0.5f, 0, 0),
				.local_offset_b = V3(-link_len * 0.5f, 0, 0),
			});
		}
		apush(bodies, b);
		prev = b;
	}

	// Record initial positions for interior dynamic bodies (index 1..n-2).
	v3 initial[16] = {0};
	for (int i = 1; i < n - 1; i++) initial[i] = body_get_position(w, bodies[i]);

	// Step 60 frames (1 second at 60 Hz) -- before any sphere impact in the
	// real scene. The chain is purely hanging under gravity with static anchors.
	for (int f = 0; f < 60; f++) world_step(w, 1.0f / 60.0f);

	// Measure drift. Interior bodies should stay very close to their initial
	// positions (the constraint is rigid). Report max drift across all axes.
	float max_x_drift = 0.0f, max_y_drift = 0.0f, max_z_drift = 0.0f;
	int worst_body = -1;
	float worst_x = 0.0f;
	for (int i = 1; i < n - 1; i++) {
		v3 p = body_get_position(w, bodies[i]);
		float dx = fabsf(p.x - initial[i].x);
		float dy = fabsf(p.y - initial[i].y);
		float dz = fabsf(p.z - initial[i].z);
		if (dx > max_x_drift) { max_x_drift = dx; worst_body = i; worst_x = p.x - initial[i].x; }
		if (dy > max_y_drift) max_y_drift = dy;
		if (dz > max_z_drift) max_z_drift = dz;
	}

	printf("  [weld-bridge %s] max_drift x=%.3f y=%.3f z=%.3f  worst=body%d dx=%+.3f (initial x=%.3f)\n",
		use_ldl ? "LDL" : "PGS", (double)max_x_drift, (double)max_y_drift, (double)max_z_drift,
		worst_body, (double)worst_x, worst_body >= 0 ? (double)initial[worst_body].x : 0.0);

	// Bridge is over-constrained but structurally rigid. Expect < 0.2 m drift
	// after 1 second. Anything above that is solver drift, not physics.
	TEST_ASSERT(max_x_drift < 0.2f);
	TEST_ASSERT(max_y_drift < 0.2f);
	TEST_ASSERT(max_z_drift < 0.2f);

	afree(bodies);
	destroy_world(w);
}

static void run_weld_bridge_unit_tests()
{
	printf("--- weld bridge drift tests ---\n");
	test_weld_bridge_drift(0); // PGS baseline
	test_weld_bridge_drift(1); // LDL under test
}
