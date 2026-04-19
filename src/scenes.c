// scenes.c -- scene setup functions for nudge demo app
// Included from main.c (unity build). Accesses globals defined there:
//   g_world, g_draw_list, DrawEntry, MESH_SPHERE, MESH_BOX,
//   g_test_hull, g_mesh_capsule, g_mesh_hull, CAP_HALF_H, CAP_RADIUS,
//   render_create_capsule_mesh, render_create_hull_mesh

// Scene-specific globals (shape showcase)
#define CHAIN_LEN 5
static Body g_chain[CHAIN_LEN];
static Body g_chain_anchor;
static Body g_spring_a, g_spring_b;

// Helper: create floor and add to draw list
static Body add_floor()
{
	Body floor = create_body(g_world, (BodyParams){
		.position = V3(0, -1, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, floor, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(10, 1, 10),
	});
	apush(g_draw_list, ((DrawEntry){ floor, MESH_BOX, V3(10, 1, 10), V3(0.4f, 0.4f, 0.45f) }));
	return floor;
}

// ---------------------------------------------------------------------------
// Scene: Shape Showcase (all shape types, chain, spring)
// ---------------------------------------------------------------------------
static void scene_showcase_setup()
{
	add_floor();

	// Dynamic sphere (bouncy)
	Body sphere = create_body(g_world, (BodyParams){
		.position = V3(-3, 5, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.restitution = 0.5f,
	});
	body_add_shape(g_world, sphere, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.5f,
	});
	apush(g_draw_list, ((DrawEntry){ sphere, MESH_SPHERE, V3(0.5f, 0.5f, 0.5f), V3(0.9f, 0.3f, 0.2f) }));

	// Dynamic capsule
	Body capsule = create_body(g_world, (BodyParams){
		.position = V3(-1, 6, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.restitution = 0.5f,
	});
	body_add_shape(g_world, capsule, (ShapeParams){
		.type = SHAPE_CAPSULE,
		.capsule = { .half_height = CAP_HALF_H, .radius = CAP_RADIUS },
	});
	apush(g_draw_list, ((DrawEntry){ capsule, g_mesh_capsule, V3(1, 1, 1), V3(0.2f, 0.8f, 0.3f) }));

	// Dynamic box
	Body box = create_body(g_world, (BodyParams){
		.position = V3(1, 7, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.restitution = 0.5f,
	});
	body_add_shape(g_world, box, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.4f, 0.4f, 0.4f),
	});
	apush(g_draw_list, ((DrawEntry){ box, MESH_BOX, V3(0.4f, 0.4f, 0.4f), V3(0.3f, 0.5f, 0.9f) }));

	// Dynamic hull
	Body hull_body = create_body(g_world, (BodyParams){
		.position = V3(3, 8, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.restitution = 0.5f,
	});
	body_add_shape(g_world, hull_body, (ShapeParams){
		.type = SHAPE_HULL,
		.hull = { .hull = g_test_hull, .scale = V3(1, 1, 1) },
	});
	apush(g_draw_list, ((DrawEntry){ hull_body, g_mesh_hull, V3(1, 1, 1), V3(0.9f, 0.7f, 0.2f) }));

	// Dynamic cylinder
	Body cyl = create_body(g_world, (BodyParams){
		.position = V3(5, 9, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.restitution = 0.5f,
	});
	body_add_shape(g_world, cyl, (ShapeParams){
		.type = SHAPE_CYLINDER,
		.cylinder = { .half_height = CYL_HALF_H, .radius = CYL_RADIUS },
	});
	apush(g_draw_list, ((DrawEntry){ cyl, g_mesh_cylinder, V3(1, 1, 1), V3(0.8f, 0.4f, 0.6f) }));

	// --- Pendulum chain (ball sockets) ---
	g_chain_anchor = create_body(g_world, (BodyParams){
		.position = V3(0, 8, -4),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, g_chain_anchor, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.15f,
	});

	float link_len = 1.0f;
	Body prev = g_chain_anchor;
	for (int i = 0; i < CHAIN_LEN; i++) {
		// Spawn horizontally from anchor -- zero initial joint error.
		// Gravity swings the chain down.
		g_chain[i] = create_body(g_world, (BodyParams){
			.position = V3((i + 1) * link_len, 8, -4),
			.rotation = quat_identity(),
			.mass = 0.5f,
		});
		body_add_shape(g_world, g_chain[i], (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = 0.2f,
		});
		create_distance(g_world, (DistanceParams){
			.body_a = prev,
			.body_b = g_chain[i],
			.rest_length = link_len,
		});
		apush(g_draw_list, ((DrawEntry){ g_chain[i], MESH_SPHERE, V3(0.2f, 0.2f, 0.2f), V3(0.8f, 0.4f, 0.9f) }));
		prev = g_chain[i];
	}

	// --- Spring-connected pair (distance joint) ---
	g_spring_a = create_body(g_world, (BodyParams){
		.position = V3(7, 4, -4),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_spring_a, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.3f, 0.3f, 0.3f),
	});
	g_spring_b = create_body(g_world, (BodyParams){
		.position = V3(5, 6, -4),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, g_spring_b, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.15f,
	});
	create_distance(g_world, (DistanceParams){
		.body_a = g_spring_a,
		.body_b = g_spring_b,
		.rest_length = 0,
		.spring = { .frequency = 3.0f, .damping_ratio = 0.3f },
	});
	apush(g_draw_list, ((DrawEntry){ g_spring_a, MESH_BOX, V3(0.3f, 0.3f, 0.3f), V3(0.2f, 0.7f, 0.9f) }));
	apush(g_draw_list, ((DrawEntry){ g_spring_b, MESH_SPHERE, V3(0.15f, 0.15f, 0.15f), V3(0.7f, 0.7f, 0.7f) }));
}

// ---------------------------------------------------------------------------
// Scene: Box Pyramid
// ---------------------------------------------------------------------------
static void scene_pyramid_setup()
{
	add_floor();

	float box_size = 0.5f; // half-extent
	float spacing = 1.05f;
	int base = 5;
	v3 colors[] = { V3(0.9f, 0.3f, 0.2f), V3(0.9f, 0.6f, 0.2f), V3(0.9f, 0.9f, 0.3f), V3(0.3f, 0.8f, 0.3f), V3(0.3f, 0.5f, 0.9f) };

	for (int layer = 0; layer < base; layer++) {
		int count = base - layer;
		float offset = -(count - 1) * 0.5f * spacing;
		float y = box_size + layer * spacing;
		v3 col = colors[layer % 5];
		for (int r = 0; r < count; r++) {
			for (int c = 0; c < count; c++) {
				Body b = create_body(g_world, (BodyParams){
					.position = V3(offset + c * spacing, y, offset + r * spacing),
					.rotation = quat_identity(),
					.mass = 1.0f,
				});
				body_add_shape(g_world, b, (ShapeParams){
					.type = SHAPE_BOX,
					.box.half_extents = V3(box_size, box_size, box_size),
				});
				apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(box_size, box_size, box_size), col }));
			}
		}
	}
}

// ---------------------------------------------------------------------------
// Scene: 2D Pyramid (single row deep, easier to debug contact normals)
// ---------------------------------------------------------------------------
static void scene_pyramid_2d_setup()
{
	add_floor();

	float half = 0.5f;
	int base = 15;
	v3 colors[] = { V3(0.9f, 0.3f, 0.2f), V3(0.9f, 0.6f, 0.2f), V3(0.9f, 0.9f, 0.3f), V3(0.3f, 0.8f, 0.3f), V3(0.3f, 0.5f, 0.9f) };

	for (int row = 0; row < base; row++) {
		int count = base - row;
		float startX = -(count - 1) * 0.5f;
		v3 col = colors[row % 5];
		for (int i = 0; i < count; i++) {
			Body b = create_body(g_world, (BodyParams){
				.position = V3(startX + i, half + row, 0),
				.rotation = quat_identity(),
				.mass = 1.0f,
				.friction = 0.6f,
			});
			body_add_shape(g_world, b, (ShapeParams){
				.type = SHAPE_BOX,
				.box.half_extents = V3(half, half, half),
			});
			apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(half, half, half), col }));
		}
	}
}

// ---------------------------------------------------------------------------
// Scene: Varied Stacks (columns of boxes with different heights)
// ---------------------------------------------------------------------------
static void scene_stacks_setup()
{
	add_floor();

	v3 col_colors[] = { V3(0.9f, 0.3f, 0.2f), V3(0.2f, 0.8f, 0.3f), V3(0.3f, 0.5f, 0.9f), V3(0.9f, 0.7f, 0.2f), V3(0.8f, 0.3f, 0.8f) };
	float col_x[] = { -4.0f, -2.0f, 0.0f, 2.0f, 4.0f };
	float heights[][6] = {
		{ 0.3f, 0.5f, 0.2f, 0.8f, 0.4f, 0.0f },
		{ 0.6f, 0.3f, 0.6f, 0.3f, 0.0f, 0.0f },
		{ 0.2f, 0.2f, 0.7f, 0.2f, 0.5f, 0.3f },
		{ 0.8f, 0.4f, 0.4f, 0.0f, 0.0f, 0.0f },
		{ 0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f },
	};
	float half_w = 0.45f;

	for (int col = 0; col < 5; col++) {
		float y = 0.0f;
		for (int row = 0; row < 6; row++) {
			float h = heights[col][row];
			if (h == 0.0f) break;
			y += h; // y at center of this box
			Body b = create_body(g_world, (BodyParams){
				.position = V3(col_x[col], y, 0),
				.rotation = quat_identity(),
				.mass = 1.0f,
			});
			body_add_shape(g_world, b, (ShapeParams){
				.type = SHAPE_BOX,
				.box.half_extents = V3(half_w, h, half_w),
			});
			apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(half_w, h, half_w), col_colors[col] }));
			y += h; // y at top of this box
		}
	}
}

// ---------------------------------------------------------------------------
// Scene: Friction Test -- rows of boxes with varying initial velocities.
// Row 0: linear slide (increasing speed along +X)
// Row 1: spin-only (angular velocity around Y, no linear)
// Row 2: combined linear + spin
// ---------------------------------------------------------------------------
static void scene_friction_setup()
{
	add_floor();

	float half = 0.4f;
	v3 colors[] = {
		V3(0.9f, 0.4f, 0.2f), V3(0.8f, 0.6f, 0.2f), V3(0.6f, 0.8f, 0.2f),
		V3(0.2f, 0.8f, 0.4f), V3(0.2f, 0.6f, 0.8f), V3(0.4f, 0.3f, 0.9f),
	};
	int ncols = 6;

	for (int row = 0; row < 3; row++) {
		float z = (row - 1) * 3.0f;
		for (int col = 0; col < ncols; col++) {
			float x = (col - (ncols - 1) * 0.5f) * 2.0f;
			Body b = create_body(g_world, (BodyParams){
				.position = V3(x, half, z),
				.rotation = quat_identity(),
				.mass = 1.0f,
			});
			body_add_shape(g_world, b, (ShapeParams){
				.type = SHAPE_BOX,
				.box.half_extents = V3(half, half, half),
			});
			apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(half, half, half), colors[col] }));

			float speed = (col + 1) * 1.5f;
			if (row == 0)
				body_set_velocity(g_world, b, V3(speed, 0, 0));
			else if (row == 1)
				body_set_angular_velocity(g_world, b, V3(0, speed * 3.0f, 0));
			else {
				body_set_velocity(g_world, b, V3(speed, 0, 0));
				body_set_angular_velocity(g_world, b, V3(0, speed * 2.0f, 0));
			}
		}
	}
}

// ---------------------------------------------------------------------------
// --- Heavy Chain: long chain with massive weight on the end ---
// PGS struggles with high mass ratio across many joints.
// LDL should solve the joint system exactly.
#define HEAVY_CHAIN_LEN 10
static Body g_hchain[HEAVY_CHAIN_LEN];
static Body g_hchain_anchor;

static void scene_heavy_chain_setup()
{

	// Static anchor at top
	g_hchain_anchor = create_body(g_world, (BodyParams){
		.position = V3(0, 10, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, g_hchain_anchor, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.15f,
	});

	float link_len = 0.8f;
	Body prev = g_hchain_anchor;
	for (int i = 0; i < HEAVY_CHAIN_LEN; i++) {
		int last = (i == HEAVY_CHAIN_LEN - 1);
		float mass = last ? 500.0f : 1.0f;
		float radius = last ? 2.0f : 0.15f;
		v3 color = last ? V3(0.9f, 0.2f, 0.2f) : V3(0.6f, 0.6f, 0.9f);

		g_hchain[i] = create_body(g_world, (BodyParams){
			.position = V3((i + 1) * link_len, 10, 0),
			.rotation = quat_identity(),
			.mass = mass,
		});
		body_add_shape(g_world, g_hchain[i], (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = radius,
		});
		create_distance(g_world, (DistanceParams){
			.body_a = prev,
			.body_b = g_hchain[i],
			.rest_length = link_len,
		});
		apush(g_draw_list, ((DrawEntry){ g_hchain[i], MESH_SPHERE, V3(radius, radius, radius), color }));
		prev = g_hchain[i];
	}
}

// --- Joint Demo: suspension bridge with hinged planks and hanging cables ---
static void scene_joint_demo_setup()
{
	add_floor();

	int plank_count = 12;
	float plank_w = 0.6f;   // half-width (Z)
	float plank_h = 0.06f;  // half-height (Y)
	float plank_d = 0.25f;  // half-depth (X, along bridge)
	float spacing = plank_d * 2.0f + 0.05f;
	float bridge_len = plank_count * spacing;
	float bridge_y = 5.0f;
	float tower_h = 4.0f;

	// Two static towers
	float tower_x_l = -bridge_len * 0.5f - 1.0f;
	float tower_x_r =  bridge_len * 0.5f + 1.0f;
	Body tower_l = create_body(g_world, (BodyParams){ .position = V3(tower_x_l, bridge_y + tower_h * 0.5f, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, tower_l, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, tower_h * 0.5f, 0.3f) });
	apush(g_draw_list, ((DrawEntry){ tower_l, MESH_BOX, V3(0.3f, tower_h * 0.5f, 0.3f), V3(0.5f, 0.45f, 0.4f) }));

	Body tower_r = create_body(g_world, (BodyParams){ .position = V3(tower_x_r, bridge_y + tower_h * 0.5f, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, tower_r, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, tower_h * 0.5f, 0.3f) });
	apush(g_draw_list, ((DrawEntry){ tower_r, MESH_BOX, V3(0.3f, tower_h * 0.5f, 0.3f), V3(0.5f, 0.45f, 0.4f) }));

	// Bridge planks connected by two distance joints each (near Z corners)
	Body planks[12];
	for (int i = 0; i < plank_count; i++) {
		float x = -bridge_len * 0.5f + spacing * 0.5f + i * spacing;
		planks[i] = create_body(g_world, (BodyParams){ .position = V3(x, bridge_y, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(g_world, planks[i], (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(plank_d, plank_h, plank_w) });
		apush(g_draw_list, ((DrawEntry){ planks[i], MESH_BOX, V3(plank_d, plank_h, plank_w), V3(0.7f, 0.55f, 0.3f) }));
	}
	// Connect consecutive planks at both Z edges (soft rope).
	// 8 Hz + damping 1.0 is loose enough for PGS/CR convergence while looking solid.
	// Higher frequency (15+) causes oscillation with iterative solvers on long chains.
	SpringParams rope = { .frequency = 8.0f, .damping_ratio = 1.0f };
	for (int i = 0; i < plank_count - 1; i++) {
		for (int side = -1; side <= 1; side += 2) {
			create_distance(g_world, (DistanceParams){ .body_a = planks[i], .body_b = planks[i+1],
				.local_offset_a = V3(plank_d, 0, side * plank_w), .local_offset_b = V3(-plank_d, 0, side * plank_w),
				.rest_length = spacing - plank_d * 2, .spring = rope });
		}
	}
	// Anchor end planks to towers at both Z edges
	for (int side = -1; side <= 1; side += 2) {
		create_distance(g_world, (DistanceParams){ .body_a = tower_l, .body_b = planks[0],
			.local_offset_a = V3(0.3f, -tower_h * 0.5f, side * plank_w),
			.local_offset_b = V3(-plank_d, 0, side * plank_w),
			.rest_length = spacing - plank_d * 2, .spring = rope });
		create_distance(g_world, (DistanceParams){ .body_a = planks[plank_count-1], .body_b = tower_r,
			.local_offset_a = V3(plank_d, 0, side * plank_w),
			.local_offset_b = V3(-0.3f, -tower_h * 0.5f, side * plank_w),
			.rest_length = spacing - plank_d * 2, .spring = rope });
	}

	// A heavy ball sitting on the bridge to stress it
	Body ball = create_body(g_world, (BodyParams){ .position = V3(0, bridge_y + 1.5f, 0), .rotation = quat_identity(), .mass = 20.0f });
	body_add_shape(g_world, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.6f });
	apush(g_draw_list, ((DrawEntry){ ball, MESH_SPHERE, V3(0.6f, 0.6f, 0.6f), V3(0.9f, 0.2f, 0.2f) }));
}

// Scene: Mass Ratio -- tiny box at bottom, each box above is larger and heavier.
// Stress test for solver stability under extreme mass ratios.
// ---------------------------------------------------------------------------
static void scene_mass_ratio_setup()
{
	add_floor();

	float sizes[] = { 0.15f, 0.25f, 0.4f, 0.55f, 0.75f };
	float masses[] = { 0.5f, 2.0f, 8.0f, 30.0f, 100.0f };
	v3 colors[] = {
		V3(0.9f, 0.2f, 0.2f), V3(0.9f, 0.5f, 0.2f), V3(0.9f, 0.8f, 0.2f),
		V3(0.4f, 0.8f, 0.3f), V3(0.3f, 0.5f, 0.9f), V3(0.5f, 0.3f, 0.9f),
	};
	float gap = 0.6f;
	float y = 0.0f;
	for (int i = 0; i < 5; i++) {
		float h = sizes[i];
		y += h + gap;
		Body b = create_body(g_world, (BodyParams){
			.position = V3(0, y, 0),
			.rotation = quat_identity(),
			.mass = masses[i],
		});
		body_add_shape(g_world, b, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(h, h, h),
		});
		apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(h, h, h), colors[i] }));
		y += h;
	}
}

// ---------------------------------------------------------------------------
// Scene: Hull Pile -- various convex hull shapes dropped in a heap.
// Stress test for hull-hull collision and feature ID stability.
// ---------------------------------------------------------------------------
static Hull* g_hull_tet;
static Hull* g_hull_wedge;
static Hull* g_hull_rock;
static int g_mesh_tet, g_mesh_wedge, g_mesh_rock;

static void hull_pile_init_shapes()
{
	static int done = 0;
	if (done) return;
	done = 1;

	// Tetrahedron
	v3 tet[] = { {0,0.6f,0}, {0.5f,-0.3f,0.3f}, {-0.5f,-0.3f,0.3f}, {0,-0.3f,-0.5f} };
	g_hull_tet = quickhull(tet, 4);
	g_mesh_tet = render_create_hull_mesh(g_hull_tet, V3(1,1,1));

	// Wedge / ramp shape
	v3 wedge[] = {
		{-0.5f,-0.3f,-0.4f}, {0.5f,-0.3f,-0.4f}, {-0.5f,-0.3f,0.4f}, {0.5f,-0.3f,0.4f},
		{-0.5f, 0.3f,-0.4f}, {0.5f, 0.3f,-0.4f},
	};
	g_hull_wedge = quickhull(wedge, 6);
	g_mesh_wedge = render_create_hull_mesh(g_hull_wedge, V3(1,1,1));

	// Irregular rock (perturbed icosahedron-ish)
	v3 rock[] = {
		{0, 0.55f, 0}, {0, -0.5f, 0},
		{0.45f, 0.1f, 0.3f}, {-0.4f, 0.15f, 0.35f},
		{0.35f, 0.1f, -0.4f}, {-0.3f, 0.1f, -0.45f},
		{0.5f, -0.15f, -0.1f}, {-0.5f, -0.1f, 0.05f},
		{0.1f, -0.2f, 0.55f}, {-0.15f, -0.2f, -0.5f},
	};
	g_hull_rock = quickhull(rock, 10);
	g_mesh_rock = render_create_hull_mesh(g_hull_rock, V3(1,1,1));
}

static void scene_hull_pile_setup()
{
	hull_pile_init_shapes();
	add_floor();

	// Hull type table: shape, mesh, scale, color
	struct { Hull* hull; int mesh; v3 scale; v3 color; } types[] = {
		{ g_hull_tet,   g_mesh_tet,   V3(0.8f,0.8f,0.8f), V3(0.9f,0.3f,0.2f) },
		{ g_hull_wedge, g_mesh_wedge, V3(0.7f,0.7f,0.7f), V3(0.2f,0.8f,0.3f) },
		{ g_hull_rock,  g_mesh_rock,  V3(0.6f,0.6f,0.6f), V3(0.3f,0.5f,0.9f) },
		{ g_test_hull,  g_mesh_hull,  V3(0.7f,0.7f,0.7f), V3(0.9f,0.7f,0.2f) },
	};
	int ntypes = sizeof(types) / sizeof(types[0]);

	// Drop a grid of hulls from various heights
	int nx = 3, nz = 3, ny = 3;
	float spacing = 1.2f;
	int idx = 0;
	for (int layer = 0; layer < ny; layer++) {
		for (int ix = 0; ix < nx; ix++) {
			for (int iz = 0; iz < nz; iz++) {
				float x = (ix - nx/2) * spacing + ((layer % 2) ? 0.3f : 0.0f);
				float z = (iz - nz/2) * spacing + ((layer % 2) ? 0.3f : 0.0f);
				float y = 1.0f + layer * 2.0f;
				int t = idx % ntypes;

				// Slightly different rotation per body
				float a = (float)idx * 1.1f;
				v3 ax = norm(V3(sinf(a), cosf(a), sinf(a*0.7f)));
				float ha = a * 0.4f;
				quat rot = (quat){ ax.x*sinf(ha), ax.y*sinf(ha), ax.z*sinf(ha), cosf(ha) };

				Body b = create_body(g_world, (BodyParams){
					.position = V3(x, y, z),
					.rotation = rot,
					.mass = 1.0f,
				});
				body_add_shape(g_world, b, (ShapeParams){
					.type = SHAPE_HULL,
					.hull = { .hull = types[t].hull, .scale = types[t].scale },
				});
				apush(g_draw_list, ((DrawEntry){ b, types[t].mesh, types[t].scale, types[t].color }));
				idx++;
			}
		}
	}

	// Toss a few boxes into the mix
	for (int i = 0; i < 4; i++) {
		float x = (i - 2) * 0.9f;
		float y = 8.0f + i * 0.5f;
		Body b = create_body(g_world, (BodyParams){
			.position = V3(x, y, 0),
			.rotation = quat_identity(),
			.mass = 1.0f,
		});
		float h = 0.3f + i * 0.05f;
		body_add_shape(g_world, b, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(h, h, h),
		});
		v3 col = V3(0.7f + i*0.05f, 0.4f, 0.8f - i*0.1f);
		apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(h, h, h), col }));
	}
}

// ---------------------------------------------------------------------------
// Scene: Weld Bridge -- fixed joints create a rigid bridge that breaks under load
// ---------------------------------------------------------------------------
static void scene_weld_bridge_setup()
{
	add_floor();

	// Bridge: chain of boxes connected by fixed joints, anchored at both ends
	int n = 8;
	float link_len = 1.0f;
	float half = 0.4f;
	float start_x = -(n - 1) * link_len * 0.5f;
	Body prev = (Body){0};
	for (int i = 0; i < n; i++) {
		float x = start_x + i * link_len;
		float mass = (i == 0 || i == n - 1) ? 0.0f : 2.0f; // endpoints static
		Body b = create_body(g_world, (BodyParams){
			.position = V3(x, 4, 0),
			.rotation = quat_identity(),
			.mass = mass,
		});
		body_add_shape(g_world, b, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(half, 0.15f, half),
		});
		v3 col = mass > 0 ? V3(0.8f, 0.6f, 0.3f) : V3(0.5f, 0.5f, 0.55f);
		apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(half, 0.15f, half), col }));
		if (i > 0) {
			create_fixed(g_world, (FixedParams){
				.body_a = prev, .body_b = b,
				.local_offset_a = V3(link_len * 0.5f, 0, 0),
				.local_offset_b = V3(-link_len * 0.5f, 0, 0),
			});
		}
		prev = b;
	}

	// Drop some heavy spheres onto the bridge
	for (int i = 0; i < 3; i++) {
		float x = start_x + (2 + i * 2) * link_len;
		Body s = create_body(g_world, (BodyParams){
			.position = V3(x, 7 + i, 0),
			.rotation = quat_identity(),
			.mass = 5.0f,
		});
		body_add_shape(g_world, s, (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = 0.4f,
		});
		apush(g_draw_list, ((DrawEntry){ s, MESH_SPHERE, V3(0.4f, 0.4f, 0.4f), V3(0.9f, 0.2f, 0.2f) }));
	}
}

// ---------------------------------------------------------------------------
// Scene: Slider Crane -- prismatic joints create rails for sliding bodies
// ---------------------------------------------------------------------------
// Persistent handles so the tick function can drive the motors over time.
static Joint g_crane_trolley_joint;
static Joint g_crane_elevator_joint;
static float g_crane_time;

static void scene_slider_crane_tick(float dt)
{
	g_crane_time += dt;
	// Trolley: sweeps -2..+2 m/s roughly every 4 seconds so the crane paces
	// back and forth along the rail. max_impulse must be large enough to
	// overcome chain drag at peak speeds.
	float trolley_speed = 2.5f * sinf(g_crane_time * 1.4f);
	joint_set_prismatic_motor(g_world, g_crane_trolley_joint, trolley_speed, 200.0f);
	// Elevator: slower oscillation to raise/lower the platform.
	float elev_speed = 1.5f * sinf(g_crane_time * 0.9f);
	joint_set_prismatic_motor(g_world, g_crane_elevator_joint, elev_speed, 200.0f);
}

static void scene_slider_crane_setup()
{
	g_crane_time = 0.0f;
	add_floor();

	// --- Crane assembly (compound_id 1) ---
	// Rail, trolley, hook arm, and hook all share compound_id so they don't
	// generate fighting contacts when they overlap (e.g. trolley sliding under
	// rail, hook swinging through the rail). Cargo is NOT in the compound,
	// so it can be picked up.
	const uint32_t crane_cid = 1;

	// Horizontal rail along X.
	Body rail = create_body(g_world, (BodyParams){
		.position = V3(0, 6, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(g_world, rail, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(5, 0.10f, 0.10f) });
	apush(g_draw_list, ((DrawEntry){ rail, MESH_BOX, V3(5, 0.10f, 0.10f), V3(0.50f, 0.50f, 0.60f) }));
	body_set_compound_id(g_world, rail, crane_cid);

	// Trolley slides along the rail; motor lets you drive it back and forth.
	Body trolley = create_body(g_world, (BodyParams){
		.position = V3(-3, 6, 0), .rotation = quat_identity(), .mass = 3.0f,
	});
	body_add_shape(g_world, trolley, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.40f, 0.30f, 0.30f) });
	apush(g_draw_list, ((DrawEntry){ trolley, MESH_BOX, V3(0.40f, 0.30f, 0.30f), V3(0.30f, 0.70f, 0.30f) }));
	body_set_compound_id(g_world, trolley, crane_cid);
	// Motor-driven bodies must stay awake so the oscillating target_speed
	// actually drives them through the zero-crossing each cycle.
	body_set_sleep_allowed(g_world, trolley, 0);
	g_crane_trolley_joint = create_prismatic(g_world, (PrismaticParams){
		.body_a = rail, .body_b = trolley,
		.local_offset_a = V3(-3, 0, 0), .local_offset_b = V3(0, 0, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(1, 0, 0),
	});
	joint_set_prismatic_motor(g_world, g_crane_trolley_joint, 0.0f, 40.0f);

	// Hook arm: a chain of 3 segments hanging from the trolley.
	const int chain_n = 3;
	const float seg_hh = 0.35f;
	Body chain[3];
	Body parent = trolley;
	float anchor_y = -0.30f; // bottom of trolley
	for (int i = 0; i < chain_n; i++) {
		float seg_y = 6.0f - 0.30f - (2 * i + 1) * seg_hh - i * 0.05f;
		chain[i] = create_body(g_world, (BodyParams){ .position = V3(-3, seg_y, 0), .rotation = quat_identity(), .mass = 0.5f });
		body_add_shape(g_world, chain[i], (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.06f, seg_hh, 0.06f) });
		apush(g_draw_list, ((DrawEntry){ chain[i], MESH_BOX, V3(0.06f, seg_hh, 0.06f), V3(0.85f, 0.65f, 0.30f) }));
		body_set_compound_id(g_world, chain[i], crane_cid);
		create_swing_twist(g_world, (SwingTwistParams){
			.body_a = parent, .body_b = chain[i],
			.local_offset_a = V3(0, anchor_y, 0),
			.local_offset_b = V3(0, seg_hh, 0),
			.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
			.cone_half_angle = 0.6f, .twist_min = -0.3f, .twist_max = 0.3f,
		});
		parent = chain[i];
		anchor_y = -seg_hh;
	}

	// Hook (hung off the bottom of the chain) -- still in the crane compound.
	Body hook = create_body(g_world, (BodyParams){ .position = V3(-3, 6.0f - 0.30f - 2 * chain_n * seg_hh - chain_n * 0.05f - 0.20f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(g_world, hook, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.30f, 0.10f, 0.30f) });
	apush(g_draw_list, ((DrawEntry){ hook, MESH_BOX, V3(0.30f, 0.10f, 0.30f), V3(0.85f, 0.40f, 0.30f) }));
	body_set_compound_id(g_world, hook, crane_cid);
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = chain[chain_n - 1], .body_b = hook,
		.local_offset_a = V3(0, -seg_hh, 0), .local_offset_b = V3(0, 0.10f, 0),
		.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
		.cone_half_angle = 0.4f, .twist_min = -0.3f, .twist_max = 0.3f,
	});

	// --- Cargo (NO compound) ---
	// Three crates of varying mass for the crane to pick up. They do collide with
	// the hook (different compound) and with each other and the floor.
	for (int i = 0; i < 3; i++) {
		float cargo_x = 1.5f + i * 0.85f;
		float cargo_size = 0.30f;
		float mass = 1.0f + i * 1.5f;
		Body crate = create_body(g_world, (BodyParams){ .position = V3(cargo_x, cargo_size, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(g_world, crate, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(cargo_size, cargo_size, cargo_size) });
		v3 col = V3(0.30f + 0.20f * i, 0.50f - 0.10f * i, 0.20f);
		apush(g_draw_list, ((DrawEntry){ crate, MESH_BOX, V3(cargo_size, cargo_size, cargo_size), col }));
	}

	// --- Elevator (compound_id 2) ---
	// Pole + platform share their own compound; cargo on top is independent.
	const uint32_t elev_cid = 2;

	Body pole = create_body(g_world, (BodyParams){ .position = V3(4, 3, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, pole, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.08f, 3, 0.08f) });
	apush(g_draw_list, ((DrawEntry){ pole, MESH_BOX, V3(0.08f, 3, 0.08f), V3(0.40f, 0.40f, 0.45f) }));
	body_set_compound_id(g_world, pole, elev_cid);

	Body platform = create_body(g_world, (BodyParams){ .position = V3(4, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(g_world, platform, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.80f, 0.10f, 0.80f) });
	apush(g_draw_list, ((DrawEntry){ platform, MESH_BOX, V3(0.80f, 0.10f, 0.80f), V3(0.30f, 0.50f, 0.80f) }));
	body_set_compound_id(g_world, platform, elev_cid);
	body_set_sleep_allowed(g_world, platform, 0);
	g_crane_elevator_joint = create_prismatic(g_world, (PrismaticParams){
		.body_a = pole, .body_b = platform,
		.local_offset_a = V3(0, 2, 0), .local_offset_b = V3(0, 0, 0),
		.local_axis_a = V3(0, 1, 0), .local_axis_b = V3(0, 1, 0),
	});
	joint_set_prismatic_motor(g_world, g_crane_elevator_joint, 0.0f, 40.0f);

	// Box sitting on the platform off to one side (pole runs through the
	// platform center; load sits clear of it). No compound -- can fall off.
	Body load = create_body(g_world, (BodyParams){ .position = V3(4.5f, 5.6f, 0.4f), .rotation = quat_identity(), .mass = 0.5f });
	body_add_shape(g_world, load, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.25f, 0.25f, 0.25f) });
	apush(g_draw_list, ((DrawEntry){ load, MESH_BOX, V3(0.25f, 0.25f, 0.25f), V3(0.90f, 0.80f, 0.20f) }));
}

// ---------------------------------------------------------------------------
// Scene: Hinge Limits -- hinged doors and pendulums with angle limits
// ---------------------------------------------------------------------------
static void scene_hinge_limits_setup()
{
	add_floor();

	// Row of pendulums with increasing angle limits
	for (int i = 0; i < 5; i++) {
		float x = -4.0f + i * 2.0f;
		float limit = 0.3f + i * 0.3f; // 0.3 to 1.5 rad
		Body anchor = create_body(g_world, (BodyParams){
			.position = V3(x, 5, 0),
			.rotation = quat_identity(),
			.mass = 0,
		});
		body_add_shape(g_world, anchor, (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = 0.15f,
		});
		apush(g_draw_list, ((DrawEntry){ anchor, MESH_SPHERE, V3(0.15f, 0.15f, 0.15f), V3(0.5f, 0.5f, 0.55f) }));

		Body arm = create_body(g_world, (BodyParams){
			.position = V3(x, 4.3f, 0),  // anchor at y=5, arm top at y=5-0.7=4.3
			.rotation = quat_identity(),
			.mass = 1.0f,
		});
		body_add_shape(g_world, arm, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(0.1f, 0.7f, 0.1f),
		});
		float t = (float)i / 4.0f;
		v3 col = V3(0.3f + t * 0.6f, 0.7f - t * 0.4f, 0.3f);
		apush(g_draw_list, ((DrawEntry){ arm, MESH_BOX, V3(0.1f, 0.7f, 0.1f), col }));

		Joint h = create_hinge(g_world, (HingeParams){
			.body_a = anchor, .body_b = arm,
			.local_offset_a = V3(0, 0, 0),
			.local_offset_b = V3(0, 0.7f, 0),
			.local_axis_a = V3(0, 0, 1),
			.local_axis_b = V3(0, 0, 1),
		});
		joint_set_hinge_limits(g_world, h, -limit, limit);
	}

	// A door: hinge with limits on a heavy box
	Body door_frame = create_body(g_world, (BodyParams){
		.position = V3(0, 2, 3),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, door_frame, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.1f, 1.5f, 0.1f),
	});
	apush(g_draw_list, ((DrawEntry){ door_frame, MESH_BOX, V3(0.1f, 1.5f, 0.1f), V3(0.5f, 0.5f, 0.55f) }));

	Body door = create_body(g_world, (BodyParams){
		.position = V3(1, 2, 3),
		.rotation = quat_identity(),
		.mass = 3.0f,
	});
	body_add_shape(g_world, door, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(1.0f, 1.5f, 0.08f),
	});
	apush(g_draw_list, ((DrawEntry){ door, MESH_BOX, V3(1.0f, 1.5f, 0.08f), V3(0.6f, 0.35f, 0.15f) }));

	Joint door_hinge = create_hinge(g_world, (HingeParams){
		.body_a = door_frame, .body_b = door,
		.local_offset_a = V3(0, 0, 0),
		.local_offset_b = V3(-1.0f, 0, 0),
		.local_axis_a = V3(0, 1, 0),
		.local_axis_b = V3(0, 1, 0),
	});
	joint_set_hinge_limits(g_world, door_hinge, -1.57f, 1.57f); // 90 degrees each way
}

// ---------------------------------------------------------------------------
// Scene: Cylinder Playground
// Mixed cylinders, boxes, spheres, and capsules for manual interaction.
// Designed to exercise all native cyl-* narrowphase pairs visually.
// ---------------------------------------------------------------------------
static void scene_cylinder_playground_setup()
{
	add_floor();

	v3 cyl_color = V3(0.8f, 0.4f, 0.6f);
	v3 box_color = V3(0.3f, 0.5f, 0.9f);
	v3 sph_color = V3(0.9f, 0.3f, 0.2f);
	v3 cap_color = V3(0.2f, 0.8f, 0.3f);

	// --- Upright cylinders (test cap-on-floor, cyl-cyl stacking) ---
	for (int i = 0; i < 3; i++) {
		Body c = create_body(g_world, (BodyParams){
			.position = V3(-3.0f + i * 1.2f, 0.5f + i * 1.2f, 0),
			.rotation = quat_identity(),
			.mass = 1.0f,
			.friction = 0.6f,
		});
		body_add_shape(g_world, c, (ShapeParams){
			.type = SHAPE_CYLINDER,
			.cylinder = { .half_height = 0.5f, .radius = 0.4f },
		});
		apush(g_draw_list, ((DrawEntry){ c, g_mesh_cylinder, V3(1, 1, 1), cyl_color }));
	}

	// --- Sideways cylinder (test side-on-floor rolling) ---
	{
		float ang = 3.14159265f * 0.5f;
		quat rot_z90 = { 0, 0, sinf(ang * 0.5f), cosf(ang * 0.5f) };
		Body c = create_body(g_world, (BodyParams){
			.position = V3(2, 0.5f, 0),
			.rotation = rot_z90,
			.mass = 1.0f,
			.friction = 0.5f,
		});
		body_add_shape(g_world, c, (ShapeParams){
			.type = SHAPE_CYLINDER,
			.cylinder = { .half_height = 0.8f, .radius = 0.4f },
		});
		int side_mesh = render_create_cylinder_mesh(0.4f, 0.8f);
		apush(g_draw_list, ((DrawEntry){ c, side_mesh, V3(1, 1, 1), V3(0.9f, 0.6f, 0.3f) }));
	}

	// --- Boxes (test cyl-box pairs) ---
	for (int i = 0; i < 3; i++) {
		Body b = create_body(g_world, (BodyParams){
			.position = V3(-2.0f + i * 2.0f, 3.0f, -2.0f),
			.rotation = quat_identity(),
			.mass = 1.0f,
			.friction = 0.5f,
		});
		body_add_shape(g_world, b, (ShapeParams){
			.type = SHAPE_BOX,
			.box.half_extents = V3(0.4f, 0.4f, 0.4f),
		});
		apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(0.4f, 0.4f, 0.4f), box_color }));
	}

	// --- Spheres (test cyl-sphere pairs) ---
	for (int i = 0; i < 3; i++) {
		Body s = create_body(g_world, (BodyParams){
			.position = V3(-1.5f + i * 1.5f, 4.0f, 2.0f),
			.rotation = quat_identity(),
			.mass = 0.5f,
			.restitution = 0.4f,
		});
		body_add_shape(g_world, s, (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = 0.35f,
		});
		apush(g_draw_list, ((DrawEntry){ s, MESH_SPHERE, V3(0.35f, 0.35f, 0.35f), sph_color }));
	}

	// --- Capsules (test cyl-capsule pairs) ---
	for (int i = 0; i < 2; i++) {
		Body c = create_body(g_world, (BodyParams){
			.position = V3(3.0f, 2.0f + i * 1.5f, -1.0f + i * 2.0f),
			.rotation = quat_identity(),
			.mass = 0.8f,
		});
		body_add_shape(g_world, c, (ShapeParams){
			.type = SHAPE_CAPSULE,
			.capsule = { .half_height = CAP_HALF_H, .radius = CAP_RADIUS },
		});
		apush(g_draw_list, ((DrawEntry){ c, g_mesh_capsule, V3(1, 1, 1), cap_color }));
	}

	// --- A big cylinder as a "barrel" to push things into ---
	{
		Body barrel = create_body(g_world, (BodyParams){
			.position = V3(0, 1.0f, 4.0f),
			.rotation = quat_identity(),
			.mass = 5.0f,
			.friction = 0.4f,
		});
		body_add_shape(g_world, barrel, (ShapeParams){
			.type = SHAPE_CYLINDER,
			.cylinder = { .half_height = 1.0f, .radius = 0.7f },
		});
		int barrel_mesh = render_create_cylinder_mesh(0.7f, 1.0f);
		apush(g_draw_list, ((DrawEntry){ barrel, barrel_mesh, V3(1, 1, 1), V3(0.5f, 0.3f, 0.15f) }));
	}

	// --- Hull body (test cyl-hull pair) ---
	{
		Body h = create_body(g_world, (BodyParams){
			.position = V3(-4, 3, 1),
			.rotation = quat_identity(),
			.mass = 1.0f,
			.restitution = 0.3f,
		});
		body_add_shape(g_world, h, (ShapeParams){
			.type = SHAPE_HULL,
			.hull = { .hull = g_test_hull, .scale = V3(0.8f, 0.8f, 0.8f) },
		});
		int hm = render_create_hull_mesh(g_test_hull, V3(0.8f, 0.8f, 0.8f));
		apush(g_draw_list, ((DrawEntry){ h, hm, V3(1, 1, 1), V3(0.9f, 0.7f, 0.2f) }));
	}
}

// ---------------------------------------------------------------------------
// Scene: Capsule Test (isolated capsule behavior on flat ground)

static Body add_big_floor()
{
	Body floor = create_body(g_world, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
	apush(g_draw_list, ((DrawEntry){ floor, MESH_BOX, V3(10, 1, 10), V3(0.4f, 0.4f, 0.45f) }));
	return floor;
}

static void scene_capsule_test_setup()
{
	add_big_floor();

	// Single cylinder for isolation
	Body c = create_body(g_world, (BodyParams){
		.position = V3(0, CYL_HALF_H + CYL_RADIUS, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.friction = 0.5f,
		.restitution = 0.0f,
	});
	body_add_shape(g_world, c, (ShapeParams){
		.type = SHAPE_CYLINDER,
		.cylinder = { .half_height = CYL_HALF_H, .radius = CYL_RADIUS },
	});
	apush(g_draw_list, ((DrawEntry){ c, g_mesh_cylinder, V3(1, 1, 1), V3(0.8f, 0.4f, 0.6f) }));
}

// ---------------------------------------------------------------------------
// Scene: Joint Gallery -- one clear example of every joint type, in a row.
// Each station is anchored to a static post on top so the joint behavior is
// obvious as the dynamic part swings/spins/extends under gravity.
// ---------------------------------------------------------------------------
static void scene_joint_gallery_setup()
{
	add_big_floor();

	// Ceiling rail of static posts. 9 stations, each ~1.6m apart.
	const int n_stations = 9;
	const float spacing = 1.8f;
	const float row_x0 = -(spacing * (n_stations - 1)) * 0.5f;
	const float ceiling_y = 5.5f;

	v3 col_post = V3(0.40f, 0.40f, 0.45f);

	#define MAKE_POST(_var, _xpos, _ypos) \
		Body _var = create_body(g_world, (BodyParams){ .position = V3(_xpos, _ypos, 0), .rotation = quat_identity(), .mass = 0 }); \
		body_add_shape(g_world, _var, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.10f, 0.10f, 0.10f) }); \
		apush(g_draw_list, ((DrawEntry){ _var, MESH_BOX, V3(0.10f, 0.10f, 0.10f), col_post }))

	#define DROP_BOX(_var, _xpos, _ypos, _hx, _hy, _hz, _mass, _color) \
		Body _var = create_body(g_world, (BodyParams){ .position = V3(_xpos, _ypos, 0), .rotation = quat_identity(), .mass = _mass }); \
		body_add_shape(g_world, _var, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(_hx, _hy, _hz) }); \
		apush(g_draw_list, ((DrawEntry){ _var, MESH_BOX, V3(_hx, _hy, _hz), _color }))

	int i = 0;

	// Station 0: BALL SOCKET -- free 3-DOF swing pendulum.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(bob, x, ceiling_y - 1.2f, 0.20f, 0.20f, 0.20f, 1.0f, V3(0.85f, 0.55f, 0.30f));
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = post, .body_b = bob,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 1.10f, 0),
		});
		i++;
	}

	// Station 1: DISTANCE -- rope (rest length, soft spring).
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		Body bob = create_body(g_world, (BodyParams){ .position = V3(x, ceiling_y - 1.4f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(g_world, bob, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.20f });
		apush(g_draw_list, ((DrawEntry){ bob, MESH_SPHERE, V3(0.20f, 0.20f, 0.20f), V3(0.30f, 0.85f, 0.55f) }));
		create_distance(g_world, (DistanceParams){
			.body_a = post, .body_b = bob,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0, 0),
			.rest_length = 1.3f, .spring = { .frequency = 4.0f, .damping_ratio = 0.7f },
		});
		i++;
	}

	// Station 2: HINGE -- arm limited to +-30 deg, swings like a clock pendulum.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(arm, x, ceiling_y - 0.8f, 0.08f, 0.7f, 0.08f, 1.0f, V3(0.55f, 0.55f, 0.85f));
		Joint h = create_hinge(g_world, (HingeParams){
			.body_a = post, .body_b = arm,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0.7f, 0),
			.local_axis_a = V3(0, 0, 1), .local_axis_b = V3(0, 0, 1),
		});
		joint_set_hinge_limits(g_world, h, -0.5f, 0.5f);
		i++;
	}

	// Station 3: FIXED -- rigid weld; bob hangs straight down with no rotation freedom.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(arm, x, ceiling_y - 0.8f, 0.08f, 0.7f, 0.08f, 1.0f, V3(0.85f, 0.85f, 0.30f));
		create_fixed(g_world, (FixedParams){
			.body_a = post, .body_b = arm,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0.7f, 0),
		});
		i++;
	}

	// Station 4: PRISMATIC + MOTOR -- elevator with motor that pushes upward.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(plat, x, ceiling_y - 1.2f, 0.40f, 0.05f, 0.40f, 1.5f, V3(0.55f, 0.85f, 0.85f));
		Joint p = create_prismatic(g_world, (PrismaticParams){
			.body_a = post, .body_b = plat,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0, 0),
			.local_axis_a = V3(0, 1, 0), .local_axis_b = V3(0, 1, 0),
		});
		joint_set_prismatic_motor(g_world, p, 0.0f, 30.0f);
		i++;
	}

	// Station 5: ANGULAR MOTOR -- spinning rotor about Y axis.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		// Anchor sphere holds rotor in place via ball socket; angular motor adds spin.
		Body rotor = create_body(g_world, (BodyParams){ .position = V3(x, ceiling_y - 0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(g_world, rotor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.6f, 0.05f, 0.10f) });
		apush(g_draw_list, ((DrawEntry){ rotor, MESH_BOX, V3(0.6f, 0.05f, 0.10f), V3(0.85f, 0.30f, 0.55f) }));
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = post, .body_b = rotor,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0, 0),
		});
		create_angular_motor(g_world, (AngularMotorParams){
			.body_a = post, .body_b = rotor,
			.local_axis_a = V3(0, 1, 0), .local_axis_b = V3(0, 1, 0),
			.target_speed = 4.0f, .max_impulse = 50.0f,
		});
		i++;
	}

	// Station 6: CONE LIMIT -- pendulum that can swing within a cone but not past.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(arm, x, ceiling_y - 1.0f, 0.06f, 0.9f, 0.06f, 1.0f, V3(0.85f, 0.55f, 0.55f));
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = post, .body_b = arm,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0.9f, 0),
		});
		// Arm's local +Y points up toward post when at rest; limit deviation to ~25 deg.
		create_cone_limit(g_world, (ConeLimitParams){
			.body_a = post, .body_b = arm,
			.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
			.half_angle = 0.45f,
		});
		i++;
	}

	// Station 7: TWIST LIMIT -- bob hangs and can spin around the rope, but only +-45 deg.
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(bob, x, ceiling_y - 1.2f, 0.50f, 0.10f, 0.10f, 1.0f, V3(0.55f, 0.30f, 0.85f));
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = post, .body_b = bob,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 1.10f, 0),
		});
		create_twist_limit(g_world, (TwistLimitParams){
			.body_a = post, .body_b = bob,
			.local_axis_a = V3(0, 1, 0), .local_axis_b = V3(0, 1, 0),
			.limit_min = -0.8f, .limit_max = 0.8f,
		});
		i++;
	}

	// Station 8: SWING-TWIST -- single ragdoll-style joint (cone + twist + ball socket).
	{
		float x = row_x0 + i * spacing;
		MAKE_POST(post, x, ceiling_y);
		DROP_BOX(bone, x, ceiling_y - 1.0f, 0.10f, 0.7f, 0.10f, 1.0f, V3(0.30f, 0.55f, 0.85f));
		create_swing_twist(g_world, (SwingTwistParams){
			.body_a = post, .body_b = bone,
			.local_offset_a = V3(0, -0.10f, 0), .local_offset_b = V3(0, 0.7f, 0),
			.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
			.cone_half_angle = 0.6f, .twist_min = -0.5f, .twist_max = 0.5f,
		});
	}

	#undef MAKE_POST
	#undef DROP_BOX
}

// ---------------------------------------------------------------------------
// Scene: Ragdoll
// Multiple humanoid figures built from swing-twist (ball + cone + twist) joints
// at the spine/neck/shoulders/hips and hinge joints at the elbows/knees.
// Each ragdoll uses a unique compound_id so its parts don't self-collide,
// while different ragdolls (and the floor) collide normally.
// ---------------------------------------------------------------------------

// Bone dimensions. Pelvis is a sphere (just a ball at the hips); everything
// else is capsules so limbs don't snag on each other.
static const float RAG_PELVIS_R  = 0.18f;              // sphere radius
static const float RAG_CHEST_HH  = 0.24f, RAG_CHEST_R  = 0.16f;
static const float RAG_HEAD_R    = 0.13f;
static const float RAG_UARM_HH   = 0.18f, RAG_LIMB_R   = 0.06f;
static const float RAG_FARM_HH   = 0.16f;
static const float RAG_THIGH_HH  = 0.22f, RAG_THIGH_R  = 0.08f;
static const float RAG_SHIN_HH   = 0.20f;

static int g_rag_mesh_chest = -1;
static int g_rag_mesh_uarm = -1, g_rag_mesh_farm  = -1;
static int g_rag_mesh_thigh = -1, g_rag_mesh_shin = -1;

static void ragdoll_init_meshes()
{
	if (g_rag_mesh_chest >= 0) return;
	g_rag_mesh_chest = render_create_capsule_mesh(RAG_CHEST_R, RAG_CHEST_HH);
	g_rag_mesh_uarm  = render_create_capsule_mesh(RAG_LIMB_R,  RAG_UARM_HH);
	g_rag_mesh_farm  = render_create_capsule_mesh(RAG_LIMB_R,  RAG_FARM_HH);
	g_rag_mesh_thigh = render_create_capsule_mesh(RAG_THIGH_R, RAG_THIGH_HH);
	g_rag_mesh_shin  = render_create_capsule_mesh(RAG_THIGH_R, RAG_SHIN_HH);
}

// Spawn one ragdoll at world position `origin`. compound_id must be nonzero
// and unique per ragdoll instance. Color tints the torso/limbs slightly so
// multiple ragdolls in the same scene are easy to tell apart.
//
// Bones are capsules (smooth cylindrical caps) so limbs don't snag on each
// other or on the floor; head is a sphere.
static void spawn_ragdoll(v3 origin, uint32_t compound_id, v3 tint)
{
	ragdoll_init_meshes();

	// Local short names for the module-static dims.
	const float PR  = RAG_PELVIS_R; // pelvis sphere radius
	const float CHH = RAG_CHEST_HH,  CR = RAG_CHEST_R;
	const float HR  = RAG_HEAD_R;
	const float UHH = RAG_UARM_HH,   LR = RAG_LIMB_R;
	const float FHH = RAG_FARM_HH;
	const float THH = RAG_THIGH_HH,  TR = RAG_THIGH_R;
	const float SHH = RAG_SHIN_HH;

	float pelvis_y = origin.y;
	float chest_y  = pelvis_y + PR + CHH + CR - 0.10f;
	float head_y   = chest_y  + CHH + CR + HR + 0.03f;
	float thigh_y  = pelvis_y - PR - THH;
	float shin_y   = thigh_y  - THH - TR - SHH - TR + 0.02f;
	float ua_y     = chest_y  + CHH - UHH;
	float fa_y     = ua_y     - UHH - LR - FHH;
	float shoulder_dx = 0.28f;
	float hip_dx      = 0.10f;

	v3 col_torso = V3(0.55f * tint.x, 0.40f * tint.y, 0.35f * tint.z);
	v3 col_head  = V3(0.85f * tint.x, 0.70f * tint.y, 0.55f * tint.z);
	v3 col_arm   = V3(0.45f * tint.x, 0.55f * tint.y, 0.75f * tint.z);
	v3 col_leg   = V3(0.40f * tint.x, 0.35f * tint.y, 0.55f * tint.z);

	// Capsule bone: physics SHAPE_CAPSULE along local +Y, drawn with a
	// pre-sized capsule mesh (scale=1 since meshes were baked to these exact
	// dimensions in ragdoll_init_meshes).
	#define MAKE_CAPSULE(_var, _pos_rel, _hh, _r, _mass, _color, _mesh) \
		Body _var = create_body(g_world, (BodyParams){ .position = add(origin, _pos_rel), .rotation = quat_identity(), .mass = (_mass) }); \
		body_add_shape(g_world, _var, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = (_hh), .radius = (_r) } }); \
		apush(g_draw_list, ((DrawEntry){ _var, _mesh, V3(1, 1, 1), _color }))

	Body pelvis = create_body(g_world, (BodyParams){ .position = V3(origin.x, pelvis_y, origin.z), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(g_world, pelvis, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = PR });
	apush(g_draw_list, ((DrawEntry){ pelvis, MESH_SPHERE, V3(PR, PR, PR), col_torso }));
	MAKE_CAPSULE(chest,  V3(0, chest_y  - origin.y, 0), CHH,  CR,  6.0f, col_torso, g_rag_mesh_chest);

	Body head = create_body(g_world, (BodyParams){ .position = V3(origin.x, head_y, origin.z), .rotation = quat_identity(), .mass = 2.0f });
	body_add_shape(g_world, head, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = HR });
	apush(g_draw_list, ((DrawEntry){ head, MESH_SPHERE, V3(HR, HR, HR), col_head }));

	MAKE_CAPSULE(ua_l, V3(-shoulder_dx, ua_y - origin.y, 0), UHH, LR, 1.8f, col_arm, g_rag_mesh_uarm);
	MAKE_CAPSULE(ua_r, V3( shoulder_dx, ua_y - origin.y, 0), UHH, LR, 1.8f, col_arm, g_rag_mesh_uarm);
	MAKE_CAPSULE(fa_l, V3(-shoulder_dx, fa_y - origin.y, 0), FHH, LR, 1.3f, col_arm, g_rag_mesh_farm);
	MAKE_CAPSULE(fa_r, V3( shoulder_dx, fa_y - origin.y, 0), FHH, LR, 1.3f, col_arm, g_rag_mesh_farm);
	MAKE_CAPSULE(thigh_l, V3(-hip_dx, thigh_y - origin.y, 0), THH, TR, 4.0f, col_leg, g_rag_mesh_thigh);
	MAKE_CAPSULE(thigh_R, V3( hip_dx, thigh_y - origin.y, 0), THH, TR, 4.0f, col_leg, g_rag_mesh_thigh);
	MAKE_CAPSULE(shin_l,  V3(-hip_dx, shin_y  - origin.y, 0), SHH,  TR, 2.5f, col_leg, g_rag_mesh_shin);
	MAKE_CAPSULE(shin_r,  V3( hip_dx, shin_y  - origin.y, 0), SHH,  TR, 2.5f, col_leg, g_rag_mesh_shin);
	#undef MAKE_CAPSULE

	// One compound_id per ragdoll: parts of the SAME ragdoll skip pairwise
	// collision; parts of DIFFERENT ragdolls (and world geometry) collide
	// normally. Scales to unlimited ragdolls -- no bit budget.
	Body parts[] = { pelvis, chest, head, ua_l, ua_r, fa_l, fa_r, thigh_l, thigh_R, shin_l, shin_r };
	for (int i = 0; i < (int)(sizeof(parts) / sizeof(parts[0])); i++)
		body_set_compound_id(g_world, parts[i], compound_id);

	const v3 axis_up = V3(0, 1, 0);

	// Spine: pelvis <-> chest. Anchors at the capsule tips (along +/- local Y
	// at half_height, since the capsule cap sphere extends past that by
	// radius -- we anchor just past the sphere to avoid joint-driven penetration).
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = pelvis, .body_b = chest,
		.local_offset_a = V3(0,  PR - 0.04f, 0),
		.local_offset_b = V3(0, -CHH  - CR  + 0.04f, 0),
		.local_axis_a = axis_up, .local_axis_b = axis_up,
		.cone_half_angle = 0.6f, .twist_min = -0.45f, .twist_max = 0.45f,
	});
	// Neck.
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = chest, .body_b = head,
		.local_offset_a = V3(0,  CHH + CR - 0.02f, 0),
		.local_offset_b = V3(0, -HR, 0),
		.local_axis_a = axis_up, .local_axis_b = axis_up,
		.cone_half_angle = 0.7f, .twist_min = -0.6f, .twist_max = 0.6f,
	});
	// Shoulders.
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = chest, .body_b = ua_l,
		.local_offset_a = V3(-shoulder_dx,  CHH - 0.05f, 0),
		.local_offset_b = V3(0, UHH, 0),
		.local_axis_a = V3(-1, 0, 0), .local_axis_b = V3(0, 1, 0),
		.cone_half_angle = 1.4f, .twist_min = -0.6f, .twist_max = 0.6f,
	});
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = chest, .body_b = ua_r,
		.local_offset_a = V3( shoulder_dx,  CHH - 0.05f, 0),
		.local_offset_b = V3(0, UHH, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(0, 1, 0),
		.cone_half_angle = 1.4f, .twist_min = -0.6f, .twist_max = 0.6f,
	});
	// Elbows.
	Joint elbow_l = create_hinge(g_world, (HingeParams){
		.body_a = ua_l, .body_b = fa_l,
		.local_offset_a = V3(0, -UHH, 0), .local_offset_b = V3(0,  FHH, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(1, 0, 0),
	});
	joint_set_hinge_limits(g_world, elbow_l, -2.2f, 0.05f);
	Joint elbow_r = create_hinge(g_world, (HingeParams){
		.body_a = ua_r, .body_b = fa_r,
		.local_offset_a = V3(0, -UHH, 0), .local_offset_b = V3(0,  FHH, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(1, 0, 0),
	});
	joint_set_hinge_limits(g_world, elbow_r, -2.2f, 0.05f);
	// Hips.
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = pelvis, .body_b = thigh_l,
		.local_offset_a = V3(-hip_dx, -PR + 0.04f, 0),
		.local_offset_b = V3(0,  THH, 0),
		.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
		.cone_half_angle = 1.0f, .twist_min = -0.4f, .twist_max = 0.4f,
	});
	create_swing_twist(g_world, (SwingTwistParams){
		.body_a = pelvis, .body_b = thigh_R,
		.local_offset_a = V3( hip_dx, -PR + 0.04f, 0),
		.local_offset_b = V3(0,  THH, 0),
		.local_axis_a = V3(0, -1, 0), .local_axis_b = V3(0, -1, 0),
		.cone_half_angle = 1.0f, .twist_min = -0.4f, .twist_max = 0.4f,
	});
	// Knees.
	Joint knee_l = create_hinge(g_world, (HingeParams){
		.body_a = thigh_l, .body_b = shin_l,
		.local_offset_a = V3(0, -THH, 0), .local_offset_b = V3(0,  SHH, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(1, 0, 0),
	});
	joint_set_hinge_limits(g_world, knee_l, -0.05f, 2.2f);
	Joint knee_r = create_hinge(g_world, (HingeParams){
		.body_a = thigh_R, .body_b = shin_r,
		.local_offset_a = V3(0, -THH, 0), .local_offset_b = V3(0,  SHH, 0),
		.local_axis_a = V3(1, 0, 0), .local_axis_b = V3(1, 0, 0),
	});
	joint_set_hinge_limits(g_world, knee_r, -0.05f, 2.2f);
}

static void scene_ragdoll_setup()
{
	add_big_floor();

	// 3x3 grid of ragdolls: each gets a unique compound_id so its own parts
	// don't collide, but different ragdolls do collide with each other.
	v3 tints[9] = {
		V3(1.0f, 1.0f, 1.0f),  // default
		V3(1.4f, 0.6f, 0.6f),  // red-ish
		V3(0.6f, 1.4f, 0.6f),  // green-ish
		V3(0.6f, 0.6f, 1.4f),  // blue-ish
		V3(1.4f, 1.4f, 0.6f),  // yellow-ish
		V3(1.4f, 0.6f, 1.4f),  // magenta-ish
		V3(0.6f, 1.4f, 1.4f),  // cyan-ish
		V3(1.2f, 1.0f, 0.7f),  // warm
		V3(0.7f, 1.0f, 1.2f),  // cool
	};
	int idx = 0;
	for (int row = -1; row <= 1; row++) {
		for (int col = -1; col <= 1; col++) {
			v3 origin = V3((float)col * 1.6f, 4.5f + (float)row * 0.05f, (float)row * 1.6f);
			spawn_ragdoll(origin, (uint32_t)(idx + 1), tints[idx]);
			idx++;
		}
	}
}

// ---------------------------------------------------------------------------
// Scene: Trimesh Terrain -- a procedural low-poly hill mesh with a V-groove.
// Drop boxes, spheres, capsules, hulls, cylinders onto it. Exercises:
//   - flat edges between coplanar triangles (no ghost collisions)
//   - convex ridges (Gauss map Voronoi classification)
//   - concave valleys (snap normal to owning face)
//   - per-triangle manifold simultaneity (shapes wedging in V-groove)
// ---------------------------------------------------------------------------
static TriMesh* g_scene_trimesh = NULL;
static int g_mesh_trimesh = -1;

static void scene_trimesh_terrain_setup()
{
	// Scene setup runs after world destroy; bodies have been freed but the
	// trimesh heap block is ours to clean up. Rebuild fresh each time.
	if (g_scene_trimesh) { trimesh_free(g_scene_trimesh); g_scene_trimesh = NULL; }

	// Grid parameters. 13x13 vertices -> 12x12 quads -> 288 triangles. Apply a
	// procedural height field with a ridge, a concave dip, and a flat plateau.
	#define TM_N 13
	#define TM_EXTENT 6.0f
	v3 verts[TM_N * TM_N];
	for (int z = 0; z < TM_N; z++) {
		for (int x = 0; x < TM_N; x++) {
			float fx = ((float)x / (TM_N - 1)) * 2.0f - 1.0f;  // [-1, 1]
			float fz = ((float)z / (TM_N - 1)) * 2.0f - 1.0f;
			float px = fx * TM_EXTENT;
			float pz = fz * TM_EXTENT;
			// Two hills (convex ridge) + a V-shaped trench through the middle.
			float hill = 1.2f * expf(-8.0f * ((fx + 0.45f) * (fx + 0.45f) + (fz - 0.45f) * (fz - 0.45f)))
			           + 0.8f * expf(-10.0f * ((fx - 0.5f) * (fx - 0.5f) + (fz + 0.4f) * (fz + 0.4f)));
			float trench = 0.6f * expf(-40.0f * fz * fz) * (1.0f - fabsf(fx));
			float py = hill - trench;
			verts[z * TM_N + x] = V3(px, py, pz);
		}
	}

	uint32_t indices[(TM_N - 1) * (TM_N - 1) * 6];
	int ti = 0;
	for (int z = 0; z < TM_N - 1; z++) {
		for (int x = 0; x < TM_N - 1; x++) {
			int i00 = z * TM_N + x;
			int i10 = i00 + 1;
			int i01 = i00 + TM_N;
			int i11 = i01 + 1;
			// Two triangles per quad, CCW when viewed from +Y (normal up).
			indices[ti++] = i00; indices[ti++] = i01; indices[ti++] = i10;
			indices[ti++] = i10; indices[ti++] = i01; indices[ti++] = i11;
		}
	}

	g_scene_trimesh = trimesh_create(verts, TM_N * TM_N, indices, (TM_N - 1) * (TM_N - 1) * 2);
	g_mesh_trimesh = render_create_trimesh_mesh(verts, TM_N * TM_N, indices, (TM_N - 1) * (TM_N - 1) * 2);

	Body floor_b = create_body(g_world, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, floor_b, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = g_scene_trimesh });
	apush(g_draw_list, ((DrawEntry){ floor_b, g_mesh_trimesh, V3(1, 1, 1), V3(0.35f, 0.5f, 0.35f) }));

	// A pile of mixed shapes dropped from varied starting positions.
	v3 drop_positions[] = {
		V3(-2.0f, 6.0f,  1.5f), V3( 1.5f, 5.5f, -1.5f),
		V3( 0.0f, 7.0f,  0.0f), V3(-3.5f, 4.5f, -0.5f),
		V3( 2.5f, 5.0f,  2.5f), V3(-1.0f, 6.5f, -2.5f),
		V3( 3.0f, 5.5f, -2.0f), V3(-2.5f, 4.0f,  2.0f),
	};
	v3 colors[] = {
		V3(0.9f, 0.3f, 0.2f), V3(0.2f, 0.8f, 0.3f),
		V3(0.3f, 0.5f, 0.9f), V3(0.9f, 0.7f, 0.2f),
		V3(0.8f, 0.4f, 0.6f), V3(0.5f, 0.9f, 0.9f),
		V3(0.95f, 0.55f, 0.35f), V3(0.7f, 0.5f, 0.8f),
	};
	for (int i = 0; i < 8; i++) {
		Body b = create_body(g_world, (BodyParams){ .position = drop_positions[i], .rotation = quat_identity(), .mass = 1.0f, .friction = 0.6f });
		switch (i % 5) {
			case 0:
				body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.35f });
				apush(g_draw_list, ((DrawEntry){ b, MESH_SPHERE, V3(0.35f, 0.35f, 0.35f), colors[i] }));
				break;
			case 1:
				body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, 0.3f, 0.3f) });
				apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(0.3f, 0.3f, 0.3f), colors[i] }));
				break;
			case 2:
				body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = CAP_HALF_H, .radius = CAP_RADIUS } });
				apush(g_draw_list, ((DrawEntry){ b, g_mesh_capsule, V3(1, 1, 1), colors[i] }));
				break;
			case 3:
				body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_HULL, .hull = { .hull = g_test_hull, .scale = V3(1, 1, 1) } });
				apush(g_draw_list, ((DrawEntry){ b, g_mesh_hull, V3(1, 1, 1), colors[i] }));
				break;
			case 4:
				body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_CYLINDER, .cylinder = { .half_height = 0.4f, .radius = 0.3f } });
				apush(g_draw_list, ((DrawEntry){ b, g_mesh_cylinder, V3(1, 1, 1), colors[i] }));
				break;
		}
	}
	#undef TM_N
	#undef TM_EXTENT
}

// ---------------------------------------------------------------------------
// Scene: Trimesh Stress -- many bodies dropped onto a large fine-tessellated
// terrain mesh. Exercises the SIMD-batched narrowphase + pair sort under load.
// ---------------------------------------------------------------------------
static void scene_trimesh_stress_setup()
{
	if (g_scene_trimesh) { trimesh_free(g_scene_trimesh); g_scene_trimesh = NULL; }

	#define TMS_N 25
	#define TMS_EXTENT 12.0f
	v3 verts[TMS_N * TMS_N];
	for (int z = 0; z < TMS_N; z++) {
		for (int x = 0; x < TMS_N; x++) {
			float fx = ((float)x / (TMS_N - 1)) * 2.0f - 1.0f;
			float fz = ((float)z / (TMS_N - 1)) * 2.0f - 1.0f;
			float px = fx * TMS_EXTENT;
			float pz = fz * TMS_EXTENT;
			// Gentle undulating terrain with a few bumps.
			float py = 0.4f * sinf(px * 0.6f) * cosf(pz * 0.6f)
			         + 0.2f * sinf(px * 1.7f + pz * 1.3f);
			verts[z * TMS_N + x] = V3(px, py, pz);
		}
	}
	uint32_t indices[(TMS_N - 1) * (TMS_N - 1) * 6];
	int ti = 0;
	for (int z = 0; z < TMS_N - 1; z++) {
		for (int x = 0; x < TMS_N - 1; x++) {
			int i00 = z * TMS_N + x;
			int i10 = i00 + 1;
			int i01 = i00 + TMS_N;
			int i11 = i01 + 1;
			indices[ti++] = i00; indices[ti++] = i01; indices[ti++] = i10;
			indices[ti++] = i10; indices[ti++] = i01; indices[ti++] = i11;
		}
	}

	g_scene_trimesh = trimesh_create(verts, TMS_N * TMS_N, indices, (TMS_N - 1) * (TMS_N - 1) * 2);
	int mesh_render = render_create_trimesh_mesh(verts, TMS_N * TMS_N, indices, (TMS_N - 1) * (TMS_N - 1) * 2);

	Body floor_b = create_body(g_world, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, floor_b, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = g_scene_trimesh });
	apush(g_draw_list, ((DrawEntry){ floor_b, mesh_render, V3(1, 1, 1), V3(0.3f, 0.45f, 0.32f) }));

	// Drop 7x7 = 49 mixed bodies from varying heights.
	int grid = 7;
	float spacing = 1.3f;
	float x0 = -(spacing * (grid - 1)) * 0.5f;
	float z0 = x0;
	int idx = 0;
	for (int gz = 0; gz < grid; gz++) {
		for (int gx = 0; gx < grid; gx++) {
			float px = x0 + gx * spacing + 0.1f * sinf((float)(gz * 7 + gx));
			float pz = z0 + gz * spacing + 0.1f * cosf((float)(gz * 5 + gx * 3));
			float py = 4.0f + 0.3f * (float)((gx + gz) % 3);
			v3 pos = V3(px, py, pz);
			v3 color = V3(0.4f + 0.5f * ((idx * 37) % 7) / 7.0f,
			              0.4f + 0.5f * ((idx * 19) % 11) / 11.0f,
			              0.4f + 0.5f * ((idx * 53) % 13) / 13.0f);
			Body b = create_body(g_world, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f, .friction = 0.6f });
			switch (idx % 5) {
				case 0:
					body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
					apush(g_draw_list, ((DrawEntry){ b, MESH_SPHERE, V3(0.3f, 0.3f, 0.3f), color }));
					break;
				case 1:
					body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, 0.3f, 0.3f) });
					apush(g_draw_list, ((DrawEntry){ b, MESH_BOX, V3(0.3f, 0.3f, 0.3f), color }));
					break;
				case 2:
					body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = CAP_HALF_H, .radius = CAP_RADIUS } });
					apush(g_draw_list, ((DrawEntry){ b, g_mesh_capsule, V3(1, 1, 1), color }));
					break;
				case 3:
					body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_HULL, .hull = { .hull = g_test_hull, .scale = V3(1, 1, 1) } });
					apush(g_draw_list, ((DrawEntry){ b, g_mesh_hull, V3(1, 1, 1), color }));
					break;
				case 4:
					body_add_shape(g_world, b, (ShapeParams){ .type = SHAPE_CYLINDER, .cylinder = { .half_height = 0.35f, .radius = 0.25f } });
					apush(g_draw_list, ((DrawEntry){ b, g_mesh_cylinder, V3(1, 1, 1), color }));
					break;
			}
			idx++;
		}
	}
	#undef TMS_N
	#undef TMS_EXTENT
}
