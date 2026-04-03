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
static bool g_sleep_enabled = true;
static int g_friction_model = FRICTION_PATCH;
static int g_solver_type = SOLVER_SOFT_STEP;
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

static Scene g_scenes[] = {
	{ "Shape Showcase",  scene_showcase_setup,  scene_showcase_draw_extras },
	{ "Box Pyramid",     scene_pyramid_setup,   NULL },
	{ "Varied Stacks",   scene_stacks_setup,    NULL },
	{ "Friction Test",   scene_friction_setup,  NULL },
	{ "Mass Ratio",      scene_mass_ratio_setup, NULL },
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

	g_world = create_world((WorldParams){
		.gravity = V3(0, -9.81f, 0),
		.broadphase = BROADPHASE_BVH,
	});

	((WorldInternal*)g_world.id)->sleep_enabled = g_sleep_enabled;
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

	// Visualization
	ImGui_SeparatorText("Visualization");
	ImGui_Checkbox("Contacts", &g_show_contacts);
	ImGui_Checkbox("Joints", &g_show_joints);
	ImGui_Checkbox("BVH", &g_show_bvh);
	ImGui_Checkbox("Sleeping bodies", &g_show_sleep);

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
