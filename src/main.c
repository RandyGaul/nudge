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
// App implementation (placeholder).

static World g_world;
static Body g_floor;
static Body g_ball;

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

	g_world = create_world((WorldParams){
		.gravity = V3(0, -9.81f, 0),
});

	// static floor
	g_floor = create_body(g_world, (BodyParams){
		.position = V3(0, -1, 0),
		.rotation = quat_identity(),
		.mass = 0, // static
	});
	body_add_shape(g_world, g_floor, (ShapeParams){
		.type = SHAPE_BOX,
		.box.half_extents = V3(10, 1, 10),
	});

	// dynamic sphere
	g_ball = create_body(g_world, (BodyParams){
		.position = V3(0, 5, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
	});
	body_add_shape(g_world, g_ball, (ShapeParams){
		.type = SHAPE_SPHERE,
		.sphere.radius = 0.5f,
	});
}

void update()
{
	world_step(g_world, 1.0f / 60.0f);
}

void draw()
{
	glClearColor(0.1f, 0.1f, 0.12f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	float aspect = (float)g_width / (float)g_height;
	mat4 proj = mat4_perspective(1.0f, aspect, 0.1f, 100.0f);
	mat4 view = mat4_look_at(V3(0, 5, 15), V3(0, 0, 0), V3(0, 1, 0));
	mat4 vp = mul(proj, view);

	render_begin(vp);

	// Floor box
	v3 floor_pos = body_get_position(g_world, g_floor);
	quat floor_rot = body_get_rotation(g_world, g_floor);
	mat4 floor_model = mat4_trs(floor_pos, floor_rot, V3(10, 1, 10));
	render_push(MESH_BOX, floor_model, V3(0.4f, 0.4f, 0.45f), 1.0f);

	// Sphere
	v3 ball_pos = body_get_position(g_world, g_ball);
	quat ball_rot = body_get_rotation(g_world, g_ball);
	mat4 ball_model = mat4_trs(ball_pos, ball_rot, V3(0.5f, 0.5f, 0.5f));
	render_push(MESH_SPHERE, ball_model, V3(0.9f, 0.3f, 0.2f), 1.0f);

	render_end();
}
