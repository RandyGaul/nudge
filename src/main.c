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

// Capsule rendering params (baked into mesh)
static const float CAP_RADIUS = 0.3f;
static const float CAP_HALF_H = 0.5f;

void init()
{
	AppSettings settings = {
		.title = "nudge",
		.w = 1280,
		.h = 720,
	};
	start_app(settings);
	render_init();

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
	world_step(g_world, 1.0f / 60.0f);
}

static void draw_body_mesh(int mesh, Body body, v3 sc, v3 color)
{
	v3 pos = body_get_position(g_world, body);
	quat rot = body_get_rotation(g_world, body);
	render_push(mesh, mat4_trs(pos, rot, sc), color, 1.0f);
}

void draw()
{
	glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	float aspect = (float)g_width / (float)g_height;
	mat4 proj = mat4_perspective(1.0f, aspect, 0.1f, 100.0f);
	mat4 view = mat4_look_at(V3(0, 5, 15), V3(0, 1, 0), V3(0, 1, 0));
	mat4 vp = mul(proj, view);

	render_begin(vp);

	draw_body_mesh(MESH_BOX,      g_floor,     V3(10, 1, 10),           V3(0.4f, 0.4f, 0.45f));
	draw_body_mesh(MESH_SPHERE,   g_sphere,    V3(0.5f, 0.5f, 0.5f),   V3(0.9f, 0.3f, 0.2f));
	draw_body_mesh(g_mesh_capsule,g_capsule,   V3(1, 1, 1),            V3(0.2f, 0.8f, 0.3f));
	draw_body_mesh(MESH_BOX,      g_box,       V3(0.4f, 0.4f, 0.4f),  V3(0.3f, 0.5f, 0.9f));
	draw_body_mesh(g_mesh_hull,   g_hull_body, V3(1, 1, 1),            V3(0.9f, 0.7f, 0.2f));

	render_end();
}
