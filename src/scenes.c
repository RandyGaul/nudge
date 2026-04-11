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

// --- Minimal Chain: static anchor + 2 dynamic bodies for debugging ---
static Body g_mini_anchor;
static Body g_mini_a, g_mini_b;

static Body g_mini_c;

static void scene_mini_chain_setup()
{
	float link_len = 0.8f;

	g_mini_anchor = create_body(g_world, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(g_world, g_mini_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });

	g_mini_a = create_body(g_world, (BodyParams){ .position = V3(link_len, 8, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(g_world, g_mini_a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });

	g_mini_b = create_body(g_world, (BodyParams){ .position = V3(link_len * 2, 8, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(g_world, g_mini_b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });

	g_mini_c = create_body(g_world, (BodyParams){ .position = V3(link_len * 3, 8, 0), .rotation = quat_identity(), .mass = 100.0f });
	body_add_shape(g_world, g_mini_c, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	create_distance(g_world, (DistanceParams){ .body_a = g_mini_anchor, .body_b = g_mini_a, .rest_length = link_len });
	create_distance(g_world, (DistanceParams){ .body_a = g_mini_a, .body_b = g_mini_b, .rest_length = link_len });
	create_distance(g_world, (DistanceParams){ .body_a = g_mini_b, .body_b = g_mini_c, .rest_length = link_len });

	apush(g_draw_list, ((DrawEntry){ g_mini_a, MESH_SPHERE, V3(0.2f, 0.2f, 0.2f), V3(0.6f, 0.6f, 0.9f) }));
	apush(g_draw_list, ((DrawEntry){ g_mini_b, MESH_SPHERE, V3(0.2f, 0.2f, 0.2f), V3(0.6f, 0.6f, 0.9f) }));
	apush(g_draw_list, ((DrawEntry){ g_mini_c, MESH_SPHERE, V3(0.5f, 0.5f, 0.5f), V3(0.9f, 0.2f, 0.2f) }));
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

// --- Hub Star: central body with 8 radial joints (exercises body shattering) ---
#define HUB_STAR_ARMS 8
static Body g_hub_center;
static Body g_hub_arms[HUB_STAR_ARMS];

static void scene_hub_star_setup()
{
	// Central hub body (dynamic, moderate mass)
	g_hub_center = create_body(g_world, (BodyParams){
		.position = V3(0, 8, 0),
		.rotation = quat_identity(),
		.mass = 5.0f,
	});
	body_add_shape(g_world, g_hub_center, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.4f,
	});
	apush(g_draw_list, ((DrawEntry){ g_hub_center, MESH_SPHERE, V3(0.4f, 0.4f, 0.4f), V3(1, 0.8f, 0.2f) }));

	// Static anchor above hub
	Body anchor = create_body(g_world, (BodyParams){
		.position = V3(0, 10, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	create_ball_socket(g_world, (BallSocketParams){
		.body_a = anchor, .body_b = g_hub_center,
		.local_offset_a = V3(0, -1, 0), .local_offset_b = V3(0, 1, 0),
	});

	// 8 radial arms
	float arm_len = 1.5f;
	for (int i = 0; i < HUB_STAR_ARMS; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / HUB_STAR_ARMS;
		float cx = cosf(angle), cz = sinf(angle);
		v3 dir = V3(cx, 0, cz);
		v3 arm_pos = add(V3(0, 8, 0), scale(dir, arm_len));

		float mass = (i == 0) ? 20.0f : 1.0f; // one heavy arm
		float radius = (i == 0) ? 0.35f : 0.2f;
		v3 color = (i == 0) ? V3(0.9f, 0.2f, 0.2f) : V3(0.4f, 0.7f, 0.9f);

		g_hub_arms[i] = create_body(g_world, (BodyParams){
			.position = arm_pos,
			.rotation = quat_identity(),
			.mass = mass,
		});
		body_add_shape(g_world, g_hub_arms[i], (ShapeParams){
			.type = SHAPE_SPHERE,
			.sphere.radius = radius,
		});
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = g_hub_center, .body_b = g_hub_arms[i],
			.local_offset_a = scale(dir, 0.5f),
			.local_offset_b = scale(dir, -arm_len + 0.5f),
		});
		apush(g_draw_list, ((DrawEntry){ g_hub_arms[i], MESH_SPHERE, V3(radius, radius, radius), color }));
	}
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
static void scene_slider_crane_setup()
{
	add_floor();

	// Horizontal rail (static box)
	Body rail = create_body(g_world, (BodyParams){
		.position = V3(0, 6, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, rail, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(5, 0.1f, 0.1f),
	});
	apush(g_draw_list, ((DrawEntry){ rail, MESH_BOX, V3(5, 0.1f, 0.1f), V3(0.5f, 0.5f, 0.6f) }));

	// Trolley slides along the rail
	Body trolley = create_body(g_world, (BodyParams){
		.position = V3(-3, 6, 0),
		.rotation = quat_identity(),
		.mass = 3.0f,
	});
	body_add_shape(g_world, trolley, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.4f, 0.3f, 0.3f),
	});
	apush(g_draw_list, ((DrawEntry){ trolley, MESH_BOX, V3(0.4f, 0.3f, 0.3f), V3(0.3f, 0.7f, 0.3f) }));
	create_prismatic(g_world, (PrismaticParams){
		.body_a = rail, .body_b = trolley,
		.local_offset_a = V3(-3, 0, 0),
		.local_offset_b = V3(0, 0, 0),
		.local_axis_a = V3(1, 0, 0),
		.local_axis_b = V3(1, 0, 0),
	});

	// Hang a payload from the trolley via distance joint (rope)
	Body payload = create_body(g_world, (BodyParams){
		.position = V3(-3, 3, 0),
		.rotation = quat_identity(),
		.mass = 2.0f,
	});
	body_add_shape(g_world, payload, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});
	apush(g_draw_list, ((DrawEntry){ payload, MESH_BOX, V3(0.5f, 0.5f, 0.5f), V3(0.8f, 0.5f, 0.2f) }));
	create_distance(g_world, (DistanceParams){
		.body_a = trolley, .body_b = payload,
		.local_offset_a = V3(0, -0.3f, 0),
		.local_offset_b = V3(0, 0.5f, 0),
		.rest_length = 2.5f,
	});

	// Second vertical slider: an elevator platform
	Body pole = create_body(g_world, (BodyParams){
		.position = V3(4, 3, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, pole, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.08f, 3, 0.08f),
	});
	apush(g_draw_list, ((DrawEntry){ pole, MESH_BOX, V3(0.08f, 3, 0.08f), V3(0.4f, 0.4f, 0.45f) }));

	Body platform = create_body(g_world, (BodyParams){
		.position = V3(4, 5, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, platform, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.8f, 0.1f, 0.8f),
	});
	apush(g_draw_list, ((DrawEntry){ platform, MESH_BOX, V3(0.8f, 0.1f, 0.8f), V3(0.3f, 0.5f, 0.8f) }));
	Joint elev = create_prismatic(g_world, (PrismaticParams){
		.body_a = pole, .body_b = platform,
		.local_offset_a = V3(0, 2, 0),
		.local_offset_b = V3(0, 0, 0),
		.local_axis_a = V3(0, 1, 0),
		.local_axis_b = V3(0, 1, 0),
	});
	// Motor holds the platform at the anchor with gentle upward push
	joint_set_prismatic_motor(g_world, elev, 0.0f, 20.0f);

	// Box sitting on the platform
	Body cargo = create_body(g_world, (BodyParams){
		.position = V3(4, 5.6f, 0),
		.rotation = quat_identity(),
		.mass = 0.5f,
	});
	body_add_shape(g_world, cargo, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.3f, 0.3f, 0.3f),
	});
	apush(g_draw_list, ((DrawEntry){ cargo, MESH_BOX, V3(0.3f, 0.3f, 0.3f), V3(0.9f, 0.8f, 0.2f) }));
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
