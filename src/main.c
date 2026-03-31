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
// App implementation -- test scene with all shape types.

static World g_world;
static Body g_floor;
static Body g_sphere;
static Body g_capsule;
static Body g_box;
static Body g_hull_body;
static Hull* g_test_hull;
static int g_mesh_capsule;
static int g_mesh_hull;
static bool g_show_contacts = true;

// Capsule rendering params (baked into mesh)
static const float CAP_RADIUS = 0.3f;
static const float CAP_HALF_H = 0.5f;

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

	g_world = create_world((WorldParams){
		.gravity = V3(0, -9.81f, 0),
	});

	// Static floor
	g_floor = create_body(g_world, (BodyParams){
		.position = V3(0, -1, 0),
		.rotation = quat_identity(),
		.mass = 0,
	});
	body_add_shape(g_world, g_floor, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(10, 1, 10),
	});

	// Dynamic sphere
	g_sphere = create_body(g_world, (BodyParams){
		.position = V3(-3, 5, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_sphere, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.5f,
	});

	// Dynamic capsule
	g_capsule = create_body(g_world, (BodyParams){
		.position = V3(-1, 6, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_capsule, (ShapeParams){
		.type = SHAPE_CAPSULE,
		.capsule = { .half_height = CAP_HALF_H, .radius = CAP_RADIUS },
	});

	// Dynamic box
	g_box = create_body(g_world, (BodyParams){
		.position = V3(1, 7, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_box, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(0.4f, 0.4f, 0.4f),
	});

	// Dynamic hull
	g_hull_body = create_body(g_world, (BodyParams){
		.position = V3(3, 8, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_hull_body, (ShapeParams){
		.type = SHAPE_HULL,
		.hull = { .hull = g_test_hull, .scale = V3(1, 1, 1) },
	});
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

	world_step(g_world, 1.0f / 60.0f);

	// Debug panel
	ImGui_Begin("Debug", NULL, 0);
	ImGui_Checkbox("Show contacts", &g_show_contacts);
	const Contact* contacts;
	int ncontacts = world_get_contacts(g_world, &contacts);
	ImGui_Text("Contacts: %d", ncontacts);
	ImGui_End();
}

static void draw_body_mesh(int mesh, Body body, v3 sc, v3 color)
{
	v3 pos = body_get_position(g_world, body);
	quat rot = body_get_rotation(g_world, body);
	render_push(mesh, mat4_trs(pos, rot, sc), color, 1.0f);
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

	draw_body_mesh(MESH_BOX,      g_floor,     V3(10, 1, 10),           V3(0.4f, 0.4f, 0.45f));
	draw_body_mesh(MESH_SPHERE,   g_sphere,    V3(0.5f, 0.5f, 0.5f),   V3(0.9f, 0.3f, 0.2f));
	draw_body_mesh(g_mesh_capsule,g_capsule,   V3(1, 1, 1),            V3(0.2f, 0.8f, 0.3f));
	draw_body_mesh(MESH_BOX,      g_box,       V3(0.4f, 0.4f, 0.4f),  V3(0.3f, 0.5f, 0.9f));
	draw_body_mesh(g_mesh_hull,   g_hull_body, V3(1, 1, 1),            V3(0.9f, 0.7f, 0.2f));

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

	render_end();
}
