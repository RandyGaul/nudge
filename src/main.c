// See LICENSE for licensing info.
#define CKIT_IMPLEMENTATION
#include "ckit.h"

#include <SDL3/SDL.h>
#include <SDL3/SDL_opengl.h>

#include "dcimgui.h"
#include "dcimgui_impl_sdl3.h"
#include "dcimgui_impl_opengl3.h"

#include "nudge.h"
#include "nudge.c"
#include "render.c"

// -----------------------------------------------------------------------------
// App settings & platform layer.

typedef struct AppSettings
{
	const char* title;
	int w;
	int h;
} AppSettings;

static SDL_Window* g_window;
static SDL_GLContext g_glctx;
static int g_running;
static int g_width;
static int g_height;

void start_app(AppSettings settings)
{
	SDL_Init(SDL_INIT_VIDEO);

	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);

	g_window = SDL_CreateWindow(settings.title, settings.w, settings.h, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	g_glctx = SDL_GL_CreateContext(g_window);
	SDL_GL_SetSwapInterval(1);

	g_width = settings.w;
	g_height = settings.h;
	g_running = 1;
}

// App interface -- user implements these.
void init();
void update();
void draw();

static void platform_shutdown()
{
	SDL_GL_DestroyContext(g_glctx);
	SDL_DestroyWindow(g_window);
	SDL_Quit();
}

int main(int argc, char* argv[])
{
	(void)argc; (void)argv;

	init();

	while (g_running) {
		SDL_Event ev;
		while (SDL_PollEvent(&ev)) {
			cImGui_ImplSDL3_ProcessEvent(&ev);
			if (ev.type == SDL_EVENT_QUIT) g_running = 0;
			if (ev.type == SDL_EVENT_WINDOW_RESIZED) {
				g_width = ev.window.data1;
				g_height = ev.window.data2;
				glViewport(0, 0, g_width, g_height);
			}
		}
		cImGui_ImplOpenGL3_NewFrame();
		cImGui_ImplSDL3_NewFrame();
		ImGui_NewFrame();
		update();
		draw();
		ImGui_Render();
		cImGui_ImplOpenGL3_RenderDrawData(ImGui_GetDrawData());
		SDL_GL_SwapWindow(g_window);
	}

	cImGui_ImplOpenGL3_Shutdown();
	cImGui_ImplSDL3_Shutdown();
	ImGui_DestroyContext(NULL);
	platform_shutdown();
	return 0;
}

// -----------------------------------------------------------------------------
// App implementation

static World g_world;
static Hull* g_test_hull;
static int g_mesh_capsule;
static int g_mesh_hull;
static bool g_show_contacts = true;
static bool g_show_joints = true;
static bool g_show_bvh = true;
static bool g_show_sleep = true;
static bool g_show_shadows = true;
static bool g_sleep_enabled = true;
static int g_friction_model = FRICTION_PATCH;
static int g_solver_type = SOLVER_SOFT_STEP;
static bool g_ldl_enabled = false;
static bool g_ldl_debug = false;
static bool g_paused = false;
static bool g_step_once = false;

// Capsule rendering params (baked into mesh)
static const float CAP_RADIUS = 0.3f;
static const float CAP_HALF_H = 0.5f;

// Scene system: each scene has a name, setup, and optional extra draw (joints etc.)
typedef struct DrawEntry { Body body; int mesh; v3 scale; v3 color; } DrawEntry;
static DrawEntry* g_draw_list; // ckit dynamic array

typedef struct Scene {
	const char* name;
	void (*setup)();
	void (*draw_extras)();
} Scene;

static int g_scene_index = 0;

// Scene-specific globals (shape showcase)
#define CHAIN_LEN 5
static Body g_chain[CHAIN_LEN];
static Body g_chain_anchor;
static Body g_spring_a, g_spring_b;

// Forward declarations
static void setup_scene();
static void scene_showcase_setup();
static void scene_showcase_draw_extras();
static void scene_pyramid_setup();
static void scene_stacks_setup();
static void scene_friction_setup();

static void scene_mass_ratio_setup();
static void scene_heavy_chain_setup();
static void scene_heavy_chain_draw_extras();
static void scene_hub_star_setup();
static void scene_hub_star_draw_extras();
static void scene_hull_pile_setup();

static Scene g_scenes[] = {
	{ "Shape Showcase",  scene_showcase_setup,  scene_showcase_draw_extras },
	{ "Box Pyramid",     scene_pyramid_setup,   NULL },
	{ "Varied Stacks",   scene_stacks_setup,    NULL },
	{ "Friction Test",   scene_friction_setup,  NULL },
	{ "Mass Ratio",      scene_mass_ratio_setup, NULL },
	{ "Heavy Chain",     scene_heavy_chain_setup, scene_heavy_chain_draw_extras },
	{ "Hub Star",        scene_hub_star_setup,  scene_hub_star_draw_extras },
	{ "Hull Pile",       scene_hull_pile_setup,  NULL },
};
#define SCENE_COUNT (sizeof(g_scenes) / sizeof(g_scenes[0]))

// Maya-style orbit camera: yaw/pitch angles, quaternion rebuilt each frame.
// Y-locked: up is always world (0,1,0), no roll.
static v3    g_cam_focus = { 0, 1, 0 };
static float g_cam_yaw   = 0.0f;       // radians around world Y
static float g_cam_pitch = -0.25f;     // radians around local X (clamped)
static float g_cam_dist  = 15.0f;

static quat cam_orientation()
{
	// Yaw around Y, then pitch around X. No roll.
	float sy = sinf(g_cam_yaw * 0.5f),   cy = cosf(g_cam_yaw * 0.5f);
	float sp = sinf(g_cam_pitch * 0.5f), cp = cosf(g_cam_pitch * 0.5f);
	quat qy = { 0, sy, 0, cy };
	quat qp = { sp, 0, 0, cp };
	return quat_mul(qy, qp);
}

static void cam_init() {}

static mat4 cam_view_matrix()
{
	quat rot = cam_orientation();
	v3 offset = quat_rotate(rot, V3(0, 0, g_cam_dist));
	v3 eye = v3_add(g_cam_focus, offset);
	return mat4_look_at(eye, g_cam_focus, V3(0, 1, 0));
}

static void cam_orbit(float dx, float dy)
{
	float sensitivity = 0.005f;
	g_cam_yaw   += -dx * sensitivity;
	g_cam_pitch += -dy * sensitivity;
	// Clamp pitch to avoid flipping at poles
	float limit = 1.5f; // ~86 degrees
	if (g_cam_pitch < -limit) g_cam_pitch = -limit;
	if (g_cam_pitch >  limit) g_cam_pitch =  limit;
}

static void cam_pan(float dx, float dy)
{
	float sensitivity = 0.005f * g_cam_dist;
	quat rot = cam_orientation();
	v3 right = quat_rotate(rot, V3(1, 0, 0));
	v3 up    = V3(0, 1, 0); // pan in world Y, not camera Y
	g_cam_focus = v3_add(g_cam_focus, v3_scale(right, -dx * sensitivity));
	g_cam_focus = v3_add(g_cam_focus, v3_scale(up,     dy * sensitivity));
}

static void cam_zoom(float delta)
{
	g_cam_dist *= 1.0f - delta * 0.1f;
	if (g_cam_dist < 0.5f) g_cam_dist = 0.5f;
	if (g_cam_dist > 200.0f) g_cam_dist = 200.0f;
}

void init()
{
	AppSettings settings = {
		.title = "nudge",
		.w = 1280,
		.h = 720,
	};
	start_app(settings);
	render_init();
	cam_init();

	ImGui_CreateContext(NULL);
	ImGui_StyleColorsDark(NULL);
	cImGui_ImplSDL3_InitForOpenGL(g_window, g_glctx);
	cImGui_ImplOpenGL3_InitEx("#version 330");

	// Build test hull (bipyramid)
	v3 hull_pts[] = {
		{  0.0f,  0.7f,  0.0f },
		{  0.5f,  0.0f,  0.5f },
		{ -0.5f,  0.0f,  0.5f },
		{  0.5f,  0.0f, -0.5f },
		{ -0.5f,  0.0f, -0.5f },
		{  0.0f, -0.7f,  0.0f },
	};
	g_test_hull = quickhull(hull_pts, 6);

	// Register custom meshes
	g_mesh_capsule = render_create_capsule_mesh(CAP_RADIUS, CAP_HALF_H);
	g_mesh_hull = render_create_hull_mesh(g_test_hull, V3(1, 1, 1));

	setup_scene();
}

static void setup_scene()
{
	if (g_world.id) destroy_world(g_world);
	aclear(g_draw_list);
	g_ldl_debug_info.valid = 0;

	g_world = create_world((WorldParams){
		.gravity = V3(0, -9.81f, 0),
		.broadphase = BROADPHASE_BVH,
	});

	((WorldInternal*)g_world.id)->sleep_enabled = g_sleep_enabled;
	((WorldInternal*)g_world.id)->ldl_enabled = g_ldl_enabled;
	world_set_friction_model(g_world, (FrictionModel)g_friction_model);
	world_set_solver_type(g_world, (SolverType)g_solver_type);
	g_scenes[g_scene_index].setup();
}

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
	});
	body_add_shape(g_world, hull_body, (ShapeParams){
		.type = SHAPE_HULL,
		.hull = { .hull = g_test_hull, .scale = V3(1, 1, 1) },
	});
	apush(g_draw_list, ((DrawEntry){ hull_body, g_mesh_hull, V3(1, 1, 1), V3(0.9f, 0.7f, 0.2f) }));

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
		// Spawn horizontally from anchor — zero initial joint error.
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
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = prev,
			.body_b = g_chain[i],
			.local_offset_a = V3(link_len * 0.5f, 0, 0),
			.local_offset_b = V3(-link_len * 0.5f, 0, 0),
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

static void scene_showcase_draw_extras()
{
	if (!g_show_joints) return;
	v3 jcol = V3(1.0f, 0.4f, 0.1f);
	for (int i = 0; i < CHAIN_LEN; i++) {
		v3 p = body_get_position(g_world, g_chain[i]);
		v3 above = i == 0 ? body_get_position(g_world, g_chain_anchor) : body_get_position(g_world, g_chain[i-1]);
		render_debug_line(above, p, jcol);
	}
	v3 sa = body_get_position(g_world, g_spring_a);
	v3 sb = body_get_position(g_world, g_spring_b);
	render_debug_line(sa, sb, V3(0.1f, 0.9f, 1.0f));
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
		float mass = last ? 100.0f : 1.0f;
		float radius = last ? 0.5f : 0.15f;
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
		create_ball_socket(g_world, (BallSocketParams){
			.body_a = prev,
			.body_b = g_hchain[i],
			.local_offset_a = V3(link_len * 0.5f, 0, 0),
			.local_offset_b = V3(-link_len * 0.5f, 0, 0),
		});
		apush(g_draw_list, ((DrawEntry){ g_hchain[i], MESH_SPHERE, V3(radius, radius, radius), color }));
		prev = g_hchain[i];
	}
}

static void scene_heavy_chain_draw_extras()
{
	if (!g_show_joints) return;
	v3 prev_pos = body_get_position(g_world, g_hchain_anchor);
	for (int i = 0; i < HEAVY_CHAIN_LEN; i++) {
		v3 pos = body_get_position(g_world, g_hchain[i]);
		render_debug_line(prev_pos, pos, V3(1, 1, 0));
		prev_pos = pos;
	}
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

static void scene_hub_star_draw_extras()
{
	if (!g_show_joints) return;
	v3 hub_pos = body_get_position(g_world, g_hub_center);
	for (int i = 0; i < HUB_STAR_ARMS; i++) {
		v3 arm_pos = body_get_position(g_world, g_hub_arms[i]);
		render_debug_line(hub_pos, arm_pos, V3(1, 1, 0));
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

// --- LDL Debug Visualizer ---
extern LDL_DebugInfo g_ldl_debug_info;
extern int g_ldl_debug_enabled;

static ImU32 ldl_heat_color(float val, float max_val)
{
	if (max_val < 1e-12f) return ImGui_GetColorU32ImVec4((ImVec4){0.1f, 0.1f, 0.1f, 1.0f});
	float t = fabsf(val) / max_val;
	if (t > 1.0f) t = 1.0f;
	float r, g, b;
	if (t < 0.33f) { float s = t / 0.33f; r = 0; g = 0; b = s; }
	else if (t < 0.66f) { float s = (t - 0.33f) / 0.33f; r = s; g = 0; b = 1.0f - s; }
	else { float s = (t - 0.66f) / 0.34f; r = 1.0f; g = s; b = 0; }
	return ImGui_GetColorU32ImVec4((ImVec4){r, g, b, 1.0f});
}

// Build a short label for joint j (e.g. "BS0", "BS1", "D0").
static void ldl_joint_label(LDL_DebugInfo* info, int j, char* buf, int bufsize)
{
	if (info->block_types[j] == JOINT_BALL_SOCKET)
		snprintf(buf, bufsize, "BS%d", j);
	else
		snprintf(buf, bufsize, "D%d", j - info->bs_count);
}

// Find which joint owns a given DOF row index.
static int ldl_joint_for_row(LDL_DebugInfo* info, int row)
{
	for (int j = info->joint_count - 1; j >= 0; j--)
		if (row >= info->block_rows[j]) return j;
	return 0;
}

// Help marker: gray "(?) " that shows tooltip on hover.
static void ldl_help(const char* text)
{
	ImGui_SameLine();
	ImGui_TextDisabled("(?)");
	if (ImGui_IsItemHovered(0)) ImGui_SetTooltipUnformatted(text);
}

static void draw_ldl_debug()
{
	LDL_DebugInfo* info = &g_ldl_debug_info;
	if (!info->valid) { ImGui_Begin("LDL Debug", NULL, 0); ImGui_Text("No data yet"); ImGui_End(); return; }

	ImGui_Begin("LDL Debug", NULL, 0);
	int n = info->n;

	// --- Stats ---
	ImGui_SeparatorText("Stats");
	ImGui_Text("DOFs: %d  Joints: %d (BS:%d D:%d)", n, info->joint_count, info->bs_count, info->dist_count);
	ldl_help("Total scalar degrees of freedom in the joint system.\nBS = ball socket (3 DOF each), D = distance (1 DOF each).");

	float d_min = 1e18f, d_max = -1e18f;
	for (int i = 0; i < n; i++) { if (info->D[i] < d_min) d_min = info->D[i]; if (info->D[i] > d_max) d_max = info->D[i]; }
	float d_ratio = d_max / (d_min > 1e-12f ? d_min : 1e-12f);
	ImGui_Text("Pivot D: min=%.4g  max=%.4g  ratio=%.1f", d_min, d_max, d_ratio);
	ldl_help("Diagonal pivots from LDL^T factorization.\nHigh ratio = ill-conditioned system (large mass ratios\nor long chains). PGS struggles most when ratio is high.\nLDL solves exactly regardless of conditioning.");

	float max_delta = 0, total_corr = 0;
	for (int i = 0; i < n; i++) {
		float d = fabsf(info->lambda_ldl[i] - info->lambda_pgs[i]);
		if (d > max_delta) max_delta = d;
		total_corr += d * d;
	}
	ImGui_Text("Max delta: %.4g  RMS correction: %.4g", max_delta, sqrtf(total_corr / (n > 0 ? n : 1)));
	ldl_help("How much LDL changed the PGS result.\nSmall values = PGS was already accurate.\nLarge values = PGS was struggling, LDL is helping.");

	// --- Matrix Heatmap ---
	ImGui_SeparatorText("Constraint Matrix A");
	ldl_help("A = J * M^-1 * J^T -- the joint coupling matrix.\nDiagonal blocks = each joint's effective mass.\nOff-diagonal = coupling between joints sharing a body.\nBrighter = larger magnitude. Hover cells for values.");

	if (n > 0 && n <= LDL_MAX_DOF) {
		float cell = 10.0f;
		float margin = 30.0f; // space for labels
		ImDrawList* dl = ImGui_GetWindowDrawList();
		ImVec2 cursor = ImGui_GetCursorScreenPos();
		float ox = cursor.x + margin; // matrix origin x (after labels)
		float oy = cursor.y + margin; // matrix origin y (after labels)
		ImU32 text_col = ImGui_GetColorU32ImVec4((ImVec4){0.7f, 0.7f, 0.7f, 1.0f});

		// Find max magnitude for color scaling
		float a_max = 0;
		for (int i = 0; i < n * n; i++) { float v = fabsf(info->A[i]); if (v > a_max) a_max = v; }

		// Top labels (one per joint block, centered over its columns)
		for (int j = 0; j < info->joint_count; j++) {
			char lbl[8]; ldl_joint_label(info, j, lbl, sizeof(lbl));
			float bx = ox + info->block_rows[j] * cell + info->block_dofs[j] * cell * 0.5f - 8;
			ImDrawList_AddText(dl, (ImVec2){bx, cursor.y}, text_col, lbl);
		}
		// Left labels (one per joint block, centered vertically)
		for (int j = 0; j < info->joint_count; j++) {
			char lbl[8]; ldl_joint_label(info, j, lbl, sizeof(lbl));
			float by = oy + info->block_rows[j] * cell + info->block_dofs[j] * cell * 0.5f - 6;
			ImDrawList_AddText(dl, (ImVec2){cursor.x, by}, text_col, lbl);
		}

		// Draw cells
		int hover_row = -1, hover_col = -1;
		ImVec2 mouse = ImGui_GetMousePos();
		for (int row = 0; row < n; row++) {
			for (int col = 0; col < n; col++) {
				float x0 = ox + col * cell;
				float y0 = oy + row * cell;
				ImU32 c = ldl_heat_color(info->A[row * n + col], a_max);
				ImDrawList_AddRectFilled(dl, (ImVec2){x0, y0}, (ImVec2){x0 + cell - 1, y0 + cell - 1}, c);
				if (mouse.x >= x0 && mouse.x < x0 + cell && mouse.y >= y0 && mouse.y < y0 + cell) {
					hover_row = row; hover_col = col;
				}
			}
		}

		// Block boundary lines
		ImU32 line_col = ImGui_GetColorU32ImVec4((ImVec4){1, 1, 1, 0.3f});
		for (int j = 0; j < info->joint_count; j++) {
			float px = ox + info->block_rows[j] * cell;
			float py = oy + info->block_rows[j] * cell;
			ImDrawList_AddLine(dl, (ImVec2){px, oy}, (ImVec2){px, oy + n * cell}, line_col);
			ImDrawList_AddLine(dl, (ImVec2){ox, py}, (ImVec2){ox + n * cell, py}, line_col);
		}
		// Outer border
		ImDrawList_AddRect(dl, (ImVec2){ox, oy}, (ImVec2){ox + n * cell, oy + n * cell}, line_col);

		ImGui_Dummy((ImVec2){margin + n * cell + 4, margin + n * cell + 4});

		// Hover tooltip for matrix cell
		if (hover_row >= 0 && ImGui_IsMouseHoveringRect((ImVec2){ox, oy}, (ImVec2){ox + n * cell, oy + n * cell})) {
			int jr = ldl_joint_for_row(info, hover_row);
			int jc = ldl_joint_for_row(info, hover_col);
			char lr[8], lc[8];
			ldl_joint_label(info, jr, lr, sizeof(lr));
			ldl_joint_label(info, jc, lc, sizeof(lc));
			int dr = hover_row - info->block_rows[jr];
			int dc = hover_col - info->block_rows[jc];
			const char* axes = "xyz";
			if (ImGui_BeginTooltip()) {
				ImGui_Text("A[%d,%d] = %.6g", hover_row, hover_col, (double)info->A[hover_row * n + hover_col]);
				if (jr == jc)
					ImGui_Text("Diagonal: %s (DOF %c,%c)", lr, axes[dr % 3], axes[dc % 3]);
				else
					ImGui_Text("Coupling: %s.%c <-> %s.%c (shared body)", lr, axes[dr % 3], lc, axes[dc % 3]);
				ImGui_EndTooltip();
			}
		}

		// Color legend
		ImGui_TextDisabled("Color: black=0  blue=small  red=medium  yellow=large");
	}

	// --- Lambda Comparison ---
	ImGui_SeparatorText("Lambda Impulses");
	ldl_help("Accumulated constraint impulses per DOF.\nGreen = warm-start (previous frame).\nYellow = LDL exact result.\nLarge differences = system is changing rapidly.");

	if (n > 0) {
		float bar_h = 10.0f, bar_max_w = 150.0f, row_h = bar_h * 2 + 6, label_w = 50.0f, value_x = label_w + bar_max_w + 8;
		ImDrawList* dl = ImGui_GetWindowDrawList();
		ImVec2 cursor = ImGui_GetCursorScreenPos();

		float l_max = 1e-6f;
		for (int i = 0; i < n; i++) { float v = fabsf(info->lambda_pgs[i]); if (v > l_max) l_max = v; v = fabsf(info->lambda_ldl[i]); if (v > l_max) l_max = v; }

		ImU32 pgs_col = ImGui_GetColorU32ImVec4((ImVec4){0.2f, 0.8f, 0.2f, 0.9f});
		ImU32 ldl_col = ImGui_GetColorU32ImVec4((ImVec4){1.0f, 0.9f, 0.1f, 0.9f});
		ImU32 text_col = ImGui_GetColorU32ImVec4((ImVec4){1, 1, 1, 1});
		ImU32 dim_col = ImGui_GetColorU32ImVec4((ImVec4){0.5f, 0.5f, 0.5f, 1.0f});
		int dof_idx = 0;
		for (int j = 0; j < info->joint_count; j++) {
			for (int d = 0; d < info->block_dofs[j]; d++) {
				float y = cursor.y + dof_idx * row_h;
				int ri = info->block_rows[j] + d;
				char label[16];
				if (info->block_types[j] == JOINT_BALL_SOCKET)
					snprintf(label, sizeof(label), "BS%d.%c", j, "xyz"[d]);
				else
					snprintf(label, sizeof(label), "D%d", j - info->bs_count);
				ImDrawList_AddText(dl, (ImVec2){cursor.x, y}, text_col, label);

				// PGS bar (green)
				float pgs_val = info->lambda_pgs[ri];
				float pgs_w = (fabsf(pgs_val) / l_max) * bar_max_w;
				ImDrawList_AddRectFilled(dl, (ImVec2){cursor.x + label_w, y}, (ImVec2){cursor.x + label_w + pgs_w, y + bar_h}, pgs_col);

				// LDL bar (yellow)
				float ldl_val = info->lambda_ldl[ri];
				float ldl_w = (fabsf(ldl_val) / l_max) * bar_max_w;
				ImDrawList_AddRectFilled(dl, (ImVec2){cursor.x + label_w, y + bar_h + 2}, (ImVec2){cursor.x + label_w + ldl_w, y + bar_h + 2 + bar_h}, ldl_col);

				// Numeric values
				char vals[32];
				snprintf(vals, sizeof(vals), "%.3f", (double)ldl_val);
				ImDrawList_AddText(dl, (ImVec2){cursor.x + value_x, y + bar_h * 0.5f - 4}, dim_col, vals);

				dof_idx++;
			}
		}

		ImGui_Dummy((ImVec2){value_x + 60, dof_idx * row_h + 4});
	}

	ImGui_End();
}

void update()
{
	// Camera input (skip when imgui wants the mouse)
	ImGuiIO* io = ImGui_GetIO();
	if (!io->WantCaptureMouse) {
		if (io->MouseDown[0])
			cam_orbit(io->MouseDelta.x, io->MouseDelta.y);
		if (io->MouseDown[2])
			cam_pan(io->MouseDelta.x, io->MouseDelta.y);
		if (io->MouseWheel != 0.0f)
			cam_zoom(io->MouseWheel);
	}

	if (!g_paused || g_step_once) { world_step(g_world, 1.0f / 60.0f); g_step_once = false; }

	// Debug panel
	ImGui_Begin("Debug", NULL, 0);

	// Scene selector -- fixed-width name so > button doesn't shift
	ImGui_SeparatorText("Scene");
	if (ImGui_Button("<<")) { g_scene_index = (g_scene_index + SCENE_COUNT - 1) % SCENE_COUNT; setup_scene(); }
	ImGui_SameLine();
	ImGui_Text("%-18s", g_scenes[g_scene_index].name);
	ImGui_SameLine();
	if (ImGui_Button(">>")) { g_scene_index = (g_scene_index + 1) % SCENE_COUNT; setup_scene(); }
	if (ImGui_Button("Restart")) setup_scene();
	ImGui_SameLine();
	ImGui_Checkbox("Pause", &g_paused);
	if (g_paused) { ImGui_SameLine(); if (ImGui_Button("Step")) g_step_once = true; }

	// Resolve after scene buttons -- setup_scene() may destroy/recreate the world
	WorldInternal* dbg_w = (WorldInternal*)g_world.id;

	// Systems
	ImGui_SeparatorText("Systems");
	if (ImGui_Checkbox("Sleep", &g_sleep_enabled)) {
		dbg_w->sleep_enabled = g_sleep_enabled;
		if (!g_sleep_enabled) {
			for (int i = 0; i < asize(dbg_w->islands); i++) {
				if ((dbg_w->island_gen[i] & 1) && !dbg_w->islands[i].awake)
					island_wake(dbg_w, i);
			}
		}
	}
	if (ImGui_Combo("Friction", &g_friction_model, "Coulomb\0Patch\0"))
		world_set_friction_model(g_world, (FrictionModel)g_friction_model);
	if (ImGui_Combo("Solver", &g_solver_type, "Soft Step\0SI Soft\0SI\0Block\0AVBD\0"))
		world_set_solver_type(g_world, (SolverType)g_solver_type);
	if (g_solver_type != SOLVER_AVBD) {
		if (ImGui_Checkbox("LDL Joints", &g_ldl_enabled))
			dbg_w->ldl_enabled = g_ldl_enabled;
		if (g_ldl_enabled) {
			ImGui_SameLine();
			ImGui_Checkbox("Debug##ldl", &g_ldl_debug);
			g_ldl_debug_enabled = g_ldl_debug;
		} else {
			g_ldl_debug = false;
			g_ldl_debug_enabled = 0;
		}
	} else {
		g_ldl_enabled = false;
		g_ldl_debug = false;
		dbg_w->ldl_enabled = 0;
		g_ldl_debug_enabled = 0;
	}

	// Visualization
	ImGui_SeparatorText("Visualization");
	ImGui_Checkbox("Contacts", &g_show_contacts);
	ImGui_Checkbox("Joints", &g_show_joints);
	ImGui_Checkbox("BVH", &g_show_bvh);
	ImGui_Checkbox("Sleeping bodies", &g_show_sleep);
	ImGui_Checkbox("Shadows", &g_show_shadows);

	// Stats
	ImGui_SeparatorText("Stats");
	const Contact* contacts;
	int ncontacts = world_get_contacts(g_world, &contacts);
	ImGui_Text("Broadphase: %s", dbg_w->broadphase_type == BROADPHASE_BVH ? "BVH" : "N^2");
	ImGui_Text("Contacts: %d", ncontacts);
	{
		int n_islands = 0, n_sleeping = 0;
		for (int i = 0; i < asize(dbg_w->islands); i++) {
			if (!(dbg_w->island_gen[i] & 1)) continue;
			n_islands++;
			if (!dbg_w->islands[i].awake) n_sleeping++;
		}
		ImGui_Text("Islands: %d (%d sleeping)", n_islands, n_sleeping);
	}
	ImGui_Text("Bodies: %d", asize(g_draw_list));
	ImGui_End();

	if (g_ldl_debug && g_ldl_enabled)
		draw_ldl_debug();
}

static void draw_body_mesh(int mesh, Body body, v3 sc, v3 color)
{
	v3 pos = body_get_position(g_world, body);
	quat rot = body_get_rotation(g_world, body);
	float opacity = 1.0f;
	if (g_show_sleep && body_is_asleep(g_world, body)) {
		color = V3(0.3f, 0.35f, 0.5f); // desaturated blue tint
		opacity = 0.6f;
	}
	render_push(mesh, mat4_trs(pos, rot, sc), color, opacity);
}

static void draw_aabb_wireframe(v3 lo, v3 hi, v3 color)
{
	v3 c[8] = {
		V3(lo.x, lo.y, lo.z), V3(hi.x, lo.y, lo.z), V3(hi.x, hi.y, lo.z), V3(lo.x, hi.y, lo.z),
		V3(lo.x, lo.y, hi.z), V3(hi.x, lo.y, hi.z), V3(hi.x, hi.y, hi.z), V3(lo.x, hi.y, hi.z),
	};
	// Bottom face
	render_debug_line(c[0], c[1], color); render_debug_line(c[1], c[2], color);
	render_debug_line(c[2], c[3], color); render_debug_line(c[3], c[0], color);
	// Top face
	render_debug_line(c[4], c[5], color); render_debug_line(c[5], c[6], color);
	render_debug_line(c[6], c[7], color); render_debug_line(c[7], c[4], color);
	// Verticals
	render_debug_line(c[0], c[4], color); render_debug_line(c[1], c[5], color);
	render_debug_line(c[2], c[6], color); render_debug_line(c[3], c[7], color);
}

static void bvh_debug_draw_cb(v3 mn, v3 mx, int depth, int is_leaf, void* user)
{
	(void)user;
	// Color by depth: cycle through distinct hues.
	float hues[][3] = { {1,0.3f,0.3f}, {0.3f,1,0.3f}, {0.3f,0.3f,1}, {1,1,0.3f}, {1,0.3f,1}, {0.3f,1,1} };
	int hi = depth % 6;
	float a = is_leaf ? 1.0f : 0.5f;
	v3 col = V3(hues[hi][0] * a, hues[hi][1] * a, hues[hi][2] * a);
	draw_aabb_wireframe(mn, mx, col);
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	render_draw_bg(V3(0.35f, 0.38f, 0.42f), V3(0.14f, 0.14f, 0.16f));

	float aspect = (float)g_width / (float)g_height;
	mat4 proj = mat4_perspective(1.0f, aspect, 0.1f, 500.0f);
	mat4 view = cam_view_matrix();
	mat4 vp = mul(proj, view);

	render_set_shadows(g_show_shadows);
	render_begin(vp);

	// Ground grid (XZ plane, y=0)
	{
		v3 grid_color = V3(0.3f, 0.3f, 0.3f);
		float extent = 10.0f;
		float step = 1.0f;
		for (float x = -extent; x <= extent; x += step) {
			render_debug_line(V3(x, 0, -extent), V3(x, 0, extent), grid_color);
			render_debug_line(V3(-extent, 0, x), V3(extent, 0, x), grid_color);
		}
	}

	// Draw all bodies from scene draw list
	for (int i = 0; i < asize(g_draw_list); i++) {
		DrawEntry* e = &g_draw_list[i];
		draw_body_mesh(e->mesh, e->body, e->scale, e->color);
	}

	// Scene-specific extras (joint lines, etc.)
	if (g_scenes[g_scene_index].draw_extras) g_scenes[g_scene_index].draw_extras();

	// Debug: contact points and normals
	if (g_show_contacts) {
		const Contact* contacts;
		int ncontacts = world_get_contacts(g_world, &contacts);
		for (int i = 0; i < ncontacts; i++) {
			v3 p = contacts[i].point;
			v3 n = contacts[i].normal;
			float nlen = 0.3f; // visual length of normal arrow

			// Contact point as small cross
			float s = 0.04f;
			v3 yellow = V3(1.0f, 1.0f, 0.0f);
			render_debug_line(V3(p.x-s, p.y, p.z), V3(p.x+s, p.y, p.z), yellow);
			render_debug_line(V3(p.x, p.y-s, p.z), V3(p.x, p.y+s, p.z), yellow);
			render_debug_line(V3(p.x, p.y, p.z-s), V3(p.x, p.y, p.z+s), yellow);

			// Normal arrow
			v3 tip = add(p, scale(n, nlen));
			render_debug_line(p, tip, V3(0.0f, 1.0f, 0.5f));
		}
	}

	if (g_show_bvh) world_debug_bvh(g_world, bvh_debug_draw_cb, NULL);

	render_end();
}
