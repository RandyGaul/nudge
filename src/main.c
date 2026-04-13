// See LICENSE for licensing info.
#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#undef small
#undef near
#undef far

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

// Forward declarations for debug server (defined in debug_server.c, included later).
static void debug_server_init();
static void debug_server_poll();
static void debug_server_shutdown();

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

	debug_server_shutdown();
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
static int g_mesh_cylinder;
static bool g_show_contacts = true;
static bool g_show_joints = true;
static bool g_show_bvh = true;
static bool g_show_proxies = false;
static bool g_show_sleep = true;
static bool g_show_shadows = true;
static bool g_translucent_shapes = false;
static bool g_sleep_enabled = true;
static bool g_sat_hint = true;
static bool g_sat_hillclimb = true;
static bool g_box_use_hull = false;
static bool g_warm_start = true;
static bool g_incremental_np = true;
// Coulomb friction removed -- patch friction is the only mode.
static int g_solver_type = SOLVER_SOFT_STEP;
static bool g_ldl_enabled = true;
static int g_ldl_inspect_island = -1;   // selected island for LDL inspector (-1 = none)
static int g_ldl_inspect_step = 0;      // factorization step slider
static int g_ldl_hover_body = -1;       // body highlighted by matrix hover (-1 = none)
static bool g_paused = false;
static bool g_step_once = false;
static float g_time_scale = 1.0f; // slow-mo: 0.1 = 10x slow, 1 = normal
static int g_run_frames = 0;      // >0: run this many frames then auto-pause

// Debug replay: highlights, label overlay
#define MAX_HIGHLIGHTS 64
static struct { int body_idx; v3 color; int active; } g_highlights[MAX_HIGHLIGHTS];
static int g_highlight_count;
static char g_label_text[256];
static float g_label_timer; // seconds remaining to show label

// Capsule rendering params (baked into mesh)
static const float CAP_RADIUS = 0.3f;
static const float CAP_HALF_H = 0.5f;
static const float CYL_RADIUS = 0.4f;
static const float CYL_HALF_H = 0.5f;

// Scene system: each scene has a name, setup, and optional extra draw (joints etc.)
typedef struct DrawEntry { Body body; int mesh; v3 scale; v3 color; } DrawEntry;
static DrawEntry* g_draw_list; // ckit dynamic array

typedef struct Scene {
	const char* name;
	void (*setup)();
} Scene;

static int g_scene_index = 0; // Shape Showcase

// Mouse constraint state (right-click drag to interact with bodies)
static Body g_mouse_body;         // body being dragged ({0} if none)
static Body g_mouse_anchor;      // hidden static body at mouse target
static Joint g_mouse_joint;      // soft ball-socket connecting anchor to body
static v3 g_mouse_local_hit;     // hit point in body's local space
static float g_mouse_ray_dist;   // camera distance at pick time
static int g_mouse_dragging;     // 1 while right-click is held (even if no body hit)

// --- Mouse input recorder ---
// Records mouse drag events per frame for deterministic test replay.
// Toggle with F5. Saves to mouse_recording.c on stop.
typedef struct
{
	int type;  // 0=step (no drag), 1=begin, 2=update, 3=end
	int body_draw_idx;   // index into g_draw_list at pick time (-1 if none)
	v3 local_hit;        // hit point in body-local space
	v3 anchor_pos;       // world-space anchor position this frame
	float ray_dist;      // camera distance at pick time
} RecordedFrame;

static CK_DYNA RecordedFrame* g_recorded_frames;
static int g_recording;
static int g_rec_body_draw_idx;  // draw list index of picked body
static int g_rec_start_frame;    // world frame counter when recording started

static void recording_save();

// Forward declarations
static void setup_scene();

#include "scenes.c"

static Scene g_scenes[] = {
	{ "Shape Showcase",  scene_showcase_setup },
	{ "Box Pyramid",     scene_pyramid_setup },
	{ "Pyramid 2D",      scene_pyramid_2d_setup },
	{ "Varied Stacks",   scene_stacks_setup },
	{ "Friction Test",   scene_friction_setup },
	{ "Mass Ratio",      scene_mass_ratio_setup },
	{ "Heavy Chain",     scene_heavy_chain_setup },
	{ "Mini Chain",      scene_mini_chain_setup },
	{ "Joint Demo",      scene_joint_demo_setup },
	{ "Hub Star",        scene_hub_star_setup },
	{ "Hull Pile",       scene_hull_pile_setup },
	{ "Weld Bridge",     scene_weld_bridge_setup },
	{ "Slider Crane",    scene_slider_crane_setup },
	{ "Hinge Limits",    scene_hinge_limits_setup },
	{ "Cylinder Playground", scene_cylinder_playground_setup },
	{ "Capsule Test", scene_capsule_test_setup },
};
#define SCENE_COUNT (sizeof(g_scenes) / sizeof(g_scenes[0]))

static void recording_save()
{
	int n = asize(g_recorded_frames);
	if (n == 0) return;
	FILE* f = fopen("mouse_recording.bin", "wb");
	if (!f) { printf("Failed to open mouse_recording.bin for writing\n"); return; }
	int scene = g_scene_index;
	fwrite(&scene, sizeof(int), 1, f);
	fwrite(&g_rec_start_frame, sizeof(int), 1, f);
	fwrite(&n, sizeof(int), 1, f);
	fwrite(g_recorded_frames, sizeof(RecordedFrame), n, f);
	fclose(f);
	printf("Saved %d frames to mouse_recording.bin (scene=%d start_frame=%d nframes=%d)\n", n, scene, g_rec_start_frame, n);
}

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

// -----------------------------------------------------------------------------
// Mouse constraint: right-click to pick and drag bodies.

static void screen_to_ray(float mx, float my, v3* origin, v3* dir)
{
	float aspect = (float)g_width / (float)g_height;
	quat rot = cam_orientation();
	*origin = add(g_cam_focus, rotate(rot, V3(0, 0, g_cam_dist)));

	v3 fwd = neg(rotate(rot, V3(0, 0, 1)));
	v3 right = rotate(rot, V3(1, 0, 0));
	v3 up = rotate(rot, V3(0, 1, 0));

	float fov = 1.0f; // must match mat4_perspective call in draw()
	float half_h = tanf(fov * 0.5f);
	float half_w = half_h * aspect;

	float ndc_x = 2.0f * mx / g_width - 1.0f;
	float ndc_y = 1.0f - 2.0f * my / g_height;

	*dir = norm(add(fwd, add(scale(right, ndc_x * half_w), scale(up, ndc_y * half_h))));
}

// Ray-sphere intersection, returns t >= 0 or -1 on miss.
static float ray_sphere_simple(v3 origin, v3 dir, v3 center, float radius)
{
	v3 oc = sub(origin, center);
	float b = dot(oc, dir);
	float c = dot(oc, oc) - radius * radius;
	float disc = b * b - c;
	if (disc < 0) return -1.0f;
	float t = -b - sqrtf(disc);
	if (t < 0) t = -b + sqrtf(disc);
	return t >= 0 ? t : -1.0f;
}

static void mouse_begin_drag(float mx, float my)
{
	v3 origin, dir;
	screen_to_ray(mx, my, &origin, &dir);

	WorldInternal* w = (WorldInternal*)g_world.id;
	Body best = {0};
	float best_t = 1e30f;

	for (int i = 0; i < asize(g_draw_list); i++) {
		DrawEntry* e = &g_draw_list[i];
		int idx = handle_index(e->body);
		if (body_inv_mass(w, idx) == 0.0f) continue; // skip static

		v3 pos = body_get_position(g_world, e->body);
		float radius = fmaxf(fmaxf(e->scale.x, e->scale.y), e->scale.z);
		float t = ray_sphere_simple(origin, dir, pos, radius);
		if (t >= 0 && t < best_t) {
			best_t = t;
			best = e->body;
		}
	}

	if (!best.id) return;

	g_mouse_body = best;
	g_mouse_ray_dist = best_t;

	v3 hit_world = add(origin, scale(dir, best_t));

	// Compute hit point in body's local space
	v3 body_pos = body_get_position(g_world, best);
	quat body_rot = body_get_rotation(g_world, best);
	g_mouse_local_hit = rotate(inv(body_rot), sub(hit_world, body_pos));

	// Create hidden static anchor at hit point
	g_mouse_anchor = create_body(g_world, (BodyParams){
		.position = hit_world,
		.rotation = quat_identity(),
		.mass = 0,
	});

	// Soft ball-socket: pulls body toward anchor
	g_mouse_joint = create_ball_socket(g_world, (BallSocketParams){
		.body_a = g_mouse_anchor,
		.body_b = g_mouse_body,
		.local_offset_a = V3(0, 0, 0),
		.local_offset_b = g_mouse_local_hit,
		.spring = { .frequency = 5.0f, .damping_ratio = 0.7f },
	});
}

static void mouse_update_drag(float mx, float my)
{
	if (!g_mouse_body.id) return;

	v3 origin, dir;
	screen_to_ray(mx, my, &origin, &dir);

	v3 target = add(origin, scale(dir, g_mouse_ray_dist));

	// Move anchor body to new target
	WorldInternal* w = (WorldInternal*)g_world.id;
	int anchor_idx = handle_index(g_mouse_anchor);
	body_pos(w, anchor_idx) = target;

	// Keep dragged body's island awake
	int joint_idx = handle_index(g_mouse_joint);
	int island_id = w->joints[joint_idx].island_id;
	if (island_id >= 0 && !w->islands[island_id].awake)
		island_wake(w, island_id);
}

static void mouse_end_drag()
{
	if (!g_mouse_body.id) return;
	destroy_joint(g_world, g_mouse_joint);
	destroy_body(g_world, g_mouse_anchor);
	g_mouse_body = (Body){0};
	g_mouse_anchor = (Body){0};
	g_mouse_joint = (Joint){0};
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
	g_mesh_cylinder = render_create_cylinder_mesh(CYL_RADIUS, CYL_HALF_H);

	setup_scene();

	debug_server_init();
}

static void setup_scene()
{
	// Clear mouse drag state before destroying world
	g_mouse_body = (Body){0};
	g_mouse_anchor = (Body){0};
	g_mouse_joint = (Joint){0};
	g_mouse_dragging = 0;

	if (g_world.id) destroy_world(g_world);
	aclear(g_draw_list);
	g_ldl_debug_info.valid = 0;

	g_world = create_world((WorldParams){
		.gravity = V3(0, -9.81f, 0),
		.broadphase = BROADPHASE_BVH,
	});

	((WorldInternal*)g_world.id)->sleep_enabled = g_sleep_enabled;
	((WorldInternal*)g_world.id)->ldl_enabled = g_ldl_enabled;
	((WorldInternal*)g_world.id)->sat_hint_enabled = g_sat_hint;
	((WorldInternal*)g_world.id)->sat_hillclimb_enabled = g_sat_hillclimb;
	((WorldInternal*)g_world.id)->box_use_hull = g_box_use_hull;
	((WorldInternal*)g_world.id)->warm_start_enabled = g_warm_start;
	((WorldInternal*)g_world.id)->incremental_np_enabled = g_incremental_np;
	world_set_solver_type(g_world, (SolverType)g_solver_type);
	g_scenes[g_scene_index].setup();
}

// Debug UI panels (LDL Inspector, joint/BVH debug callbacks).
// Must be before update()/draw() which call these functions.
// Requires draw_aabb_wireframe, so forward-declare it.
static void draw_aabb_wireframe(v3 lo, v3 hi, v3 color);

#include "debug_ui.c"

static bool g_npv_mode = false;
#include "np_viz.c"
#include "debug_server.c"

void update()
{
	if (g_npv_mode) { npv_update(); return; }

	// Camera input (skip when imgui wants the mouse)
	ImGuiIO* io = ImGui_GetIO();
	if (!io->WantCaptureMouse) {
		if (io->MouseDown[0] && !g_mouse_dragging && !io->KeyCtrl)
			cam_orbit(io->MouseDelta.x, io->MouseDelta.y);
		if (io->MouseDown[2])
			cam_pan(io->MouseDelta.x, io->MouseDelta.y);
		if (io->MouseWheel != 0.0f)
			cam_zoom(io->MouseWheel);

		// Mouse constraint: right-click drag
		if (io->MouseDown[1]) {
			if (!g_mouse_dragging) {
				mouse_begin_drag(io->MousePos.x, io->MousePos.y);
				g_mouse_dragging = 1;
				if (g_recording && g_mouse_body.id) {
					g_rec_body_draw_idx = -1;
					for (int di = 0; di < asize(g_draw_list); di++) {
						if (g_draw_list[di].body.id == g_mouse_body.id) { g_rec_body_draw_idx = di; break; }
					}
					WorldInternal* rw = (WorldInternal*)g_world.id;
					int ai = handle_index(g_mouse_anchor);
					RecordedFrame rf = { .type = 1, .body_draw_idx = g_rec_body_draw_idx, .local_hit = g_mouse_local_hit, .anchor_pos = body_pos(rw, ai), .ray_dist = g_mouse_ray_dist };
					apush(g_recorded_frames, rf);
				} else if (g_recording) {
					apush(g_recorded_frames, ((RecordedFrame){ .type = 0 }));
				}
			} else {
				mouse_update_drag(io->MousePos.x, io->MousePos.y);
				if (g_recording && g_mouse_body.id) {
					WorldInternal* rw = (WorldInternal*)g_world.id;
					int ai = handle_index(g_mouse_anchor);
					RecordedFrame rf = { .type = 2, .anchor_pos = body_pos(rw, ai) };
					apush(g_recorded_frames, rf);
				} else if (g_recording) {
					apush(g_recorded_frames, ((RecordedFrame){ .type = 0 }));
				}
			}
		} else if (g_mouse_dragging) {
			mouse_end_drag();
			g_mouse_dragging = 0;
			if (g_recording) apush(g_recorded_frames, ((RecordedFrame){ .type = 3 }));
		} else if (g_recording) {
			apush(g_recorded_frames, ((RecordedFrame){ .type = 0 }));
		}
	} else if (g_mouse_dragging) {
		// Release if mouse moves over UI while dragging
		mouse_end_drag();
		g_mouse_dragging = 0;
		if (g_recording) apush(g_recorded_frames, ((RecordedFrame){ .type = 3 }));
	} else if (g_recording) {
		apush(g_recorded_frames, ((RecordedFrame){ .type = 0 }));
	}

	// LDL Inspector: Ctrl+left-click picks island
	if (!io->WantCaptureMouse && io->KeyCtrl && ImGui_IsMouseClicked(0)) {
		v3 origin, dir;
		screen_to_ray(io->MousePos.x, io->MousePos.y, &origin, &dir);
		WorldInternal* w = (WorldInternal*)g_world.id;
		int best_body = -1;
		float best_t = 1e30f;
		for (int i = 0; i < asize(g_draw_list); i++) {
			DrawEntry* e = &g_draw_list[i];
			int idx = handle_index(e->body);
			v3 pos = body_get_position(g_world, e->body);
			float radius = fmaxf(fmaxf(e->scale.x, e->scale.y), e->scale.z);
			if (radius < 0.3f) radius = 0.3f; // minimum pick radius for small bodies
			float t = ray_sphere_simple(origin, dir, pos, radius);
			if (t >= 0 && t < best_t) {
				best_t = t;
				best_body = idx;
			}
		}
		if (best_body >= 0) {
			int isl = w->body_cold[best_body].island_id;
			if (isl >= 0 && (w->island_gen[isl] & 1)) {
				g_ldl_inspect_island = isl;
				g_ldl_inspect_step = 0;
			}
		} else {
			g_ldl_inspect_island = -1;
		}
	}

	// F5: toggle recording
	if (ImGui_IsKeyPressedEx(ImGuiKey_F5, false)) {
		if (!g_recording) {
			afree(g_recorded_frames);
			g_recorded_frames = NULL;
			g_recording = 1;
			g_rec_start_frame = ((WorldInternal*)g_world.id)->frame;
			printf("[REC] Recording started (scene=%d '%s' frame=%d)\n", g_scene_index, g_scenes[g_scene_index].name, g_rec_start_frame);
		} else {
			g_recording = 0;
			recording_save();
		}
	}

	// Playback mode: drive mouse events from recording
	static int playback_frame = 0;
	if (g_recording == 2) {
		int n = asize(g_recorded_frames);
		if (playback_frame < n) {
			RecordedFrame* rf = &g_recorded_frames[playback_frame];
			WorldInternal* pw = (WorldInternal*)g_world.id;
			if (rf->type == 1 && rf->body_draw_idx >= 0 && rf->body_draw_idx < asize(g_draw_list)) {
				Body target = g_draw_list[rf->body_draw_idx].body;
				g_mouse_anchor = create_body(g_world, (BodyParams){ .position = rf->anchor_pos, .rotation = quat_identity(), .mass = 0 });
				g_mouse_joint = create_ball_socket(g_world, (BallSocketParams){ .body_a = g_mouse_anchor, .body_b = target, .local_offset_a = V3(0,0,0), .local_offset_b = rf->local_hit, .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
				g_mouse_body = target;
			} else if (rf->type == 2 && g_mouse_anchor.id) {
				body_pos(pw, handle_index(g_mouse_anchor)) = rf->anchor_pos;
			} else if (rf->type == 3 && g_mouse_anchor.id) {
				mouse_end_drag();
			}
			playback_frame++;
		} else {
			g_recording = 0;
			playback_frame = 0;
			printf("[PLAY] done\n");
		}
	}

	if (!g_paused || g_step_once) {
		world_step(g_world, (1.0f / 60.0f) * g_time_scale);
		g_step_once = false;
		if (g_run_frames > 0 && --g_run_frames == 0) g_paused = true;
	}

	debug_server_poll();

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
	if (ImGui_Button("NP Viz")) g_npv_mode = true;
	ImGui_Checkbox("Pause", &g_paused);
	if (g_paused) { ImGui_SameLine(); if (ImGui_Button("Step")) g_step_once = true; }
	if (g_recording) { ImGui_SameLine(); ImGui_TextColored((ImVec4){1,0.2f,0.2f,1}, "REC %d", asize(g_recorded_frames)); }
	if (ImGui_Button("Play Recording")) {
		FILE* rf = fopen("mouse_recording.bin", "rb");
		if (rf) {
			int rec_scene, rec_start, rec_n;
			fread(&rec_scene, sizeof(int), 1, rf);
			fread(&rec_start, sizeof(int), 1, rf);
			fread(&rec_n, sizeof(int), 1, rf);
			afree(g_recorded_frames);
			g_recorded_frames = NULL;
			afit(g_recorded_frames, rec_n);
			asetlen(g_recorded_frames, rec_n);
			fread(g_recorded_frames, sizeof(RecordedFrame), rec_n, rf);
			fclose(rf);
			g_rec_start_frame = rec_start;
			g_scene_index = rec_scene;
			setup_scene();
			// Run the exact settle frames
			WorldInternal* pw = (WorldInternal*)g_world.id;
			for (int i = pw->frame; i < rec_start; i++) world_step(g_world, 1.0f / 60.0f);
			g_recording = 2; // 2 = playback mode
			g_rec_body_draw_idx = -1;
			playback_frame = 0;
			printf("[PLAY] loaded %d frames, scene=%d start=%d\n", rec_n, rec_scene, rec_start);
		} else {
			printf("[PLAY] no mouse_recording.bin found\n");
		}
	}

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
	// Friction model: patch friction only (Coulomb removed).
	if (ImGui_Combo("Solver", &g_solver_type, "Soft Step\0SI Soft\0SI\0"))
		world_set_solver_type(g_world, (SolverType)g_solver_type);
	if (ImGui_Checkbox("LDL Joints", &g_ldl_enabled)) {
		dbg_w->ldl_enabled = g_ldl_enabled;
	}
	if (ImGui_Checkbox("SAT Axis Hints", &g_sat_hint)) dbg_w->sat_hint_enabled = g_sat_hint;
	if (ImGui_Checkbox("SAT Hill-Climb", &g_sat_hillclimb)) dbg_w->sat_hillclimb_enabled = g_sat_hillclimb;
	if (ImGui_Checkbox("Box via Hull", &g_box_use_hull)) dbg_w->box_use_hull = g_box_use_hull;
	if (ImGui_Checkbox("Warm Start", &g_warm_start)) dbg_w->warm_start_enabled = g_warm_start;
	if (ImGui_Checkbox("Incremental NP", &g_incremental_np)) dbg_w->incremental_np_enabled = g_incremental_np;
	if (!g_ldl_enabled) {
		g_ldl_inspect_island = -1;
	}

	// Visualization
	ImGui_SeparatorText("Visualization");
	ImGui_Checkbox("Contacts", &g_show_contacts);
	ImGui_Checkbox("Joints", &g_show_joints);
	ImGui_Checkbox("BVH", &g_show_bvh);
	ImGui_SameLine(); ImGui_Checkbox("Proxies", &g_show_proxies);
	ImGui_Checkbox("Sleeping bodies", &g_show_sleep);
	ImGui_Checkbox("Shadows", &g_show_shadows);
	ImGui_Checkbox("Translucent shapes", &g_translucent_shapes);

	// Stats
	ImGui_SeparatorText("Stats");
	const Contact* contacts;
	int ncontacts = world_get_contacts(g_world, &contacts);
	{ int bp = dbg_w->broadphase_type;
	  if (ImGui_Combo("Broadphase", &bp, "N^2\0BVH\0")) dbg_w->broadphase_type = bp; }
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
	if (g_ldl_inspect_island >= 0) {
		ImGui_TextDisabled("Inspecting island %d", g_ldl_inspect_island);
	} else {
		ImGui_TextDisabled("Ctrl+click body to inspect island");
	}
	ImGui_End();

	{
		extern int g_ldl_debug_island;
		g_ldl_debug_island = (g_ldl_enabled && g_ldl_inspect_island >= 0) ? g_ldl_inspect_island : -1;
		draw_ldl_inspector();
	}

	// Debug replay label overlay
	if (g_label_timer > 0) {
		g_label_timer -= 1.0f / 60.0f;
		ImGui_SetNextWindowPosEx((ImVec2){g_width * 0.5f, 40}, ImGuiCond_Always, (ImVec2){0.5f, 0});
		ImGui_SetNextWindowBgAlpha(0.7f);
		ImGui_Begin("##label", NULL, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoNav);
		ImGui_TextColored((ImVec4){1, 1, 0.3f, 1}, "%s", g_label_text);
		ImGui_End();
	}
}

static void draw_body_mesh(int mesh, Body body, v3 sc, v3 color)
{
	v3 pos = body_get_position(g_world, body);
	quat rot = body_get_rotation(g_world, body);
	float opacity = g_translucent_shapes ? 0.3f : 1.0f;
	if (g_show_sleep && body_is_asleep(g_world, body)) {
		color = V3(0.3f, 0.35f, 0.5f);
		opacity = g_translucent_shapes ? 0.15f : 0.6f;
	}
	// Debug highlight override
	int idx = handle_index(body);
	for (int h = 0; h < g_highlight_count; h++) {
		if (g_highlights[h].active && g_highlights[h].body_idx == idx) {
			color = g_highlights[h].color;
			opacity = 1.0f;
			break;
		}
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

void draw()
{
	if (g_npv_mode) { npv_draw(); return; }

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	render_draw_bg(V3(0.35f, 0.38f, 0.42f), V3(0.14f, 0.14f, 0.16f));

	float aspect = (float)g_width / (float)g_height;
	mat4 proj = mat4_perspective(1.0f, aspect, 0.1f, 500.0f);
	mat4 view = cam_view_matrix();
	mat4 vp = mul(proj, view);

	render_set_shadows(g_show_shadows);
	render_set_no_depth_write(g_translucent_shapes);
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

	// Draw all bodies from scene draw list (with island highlighting)
	WorldInternal* draw_w = (WorldInternal*)g_world.id;
	for (int i = 0; i < asize(g_draw_list); i++) {
		DrawEntry* e = &g_draw_list[i];
		v3 color = e->color;
		if (g_ldl_inspect_island >= 0) {
			int idx = handle_index(e->body);
			int isl = draw_w->body_cold[idx].island_id;
			if (isl == g_ldl_inspect_island) {
				// Warm orange tint for selected island
				color.x = color.x * 0.6f + 1.0f * 0.4f;
				color.y = color.y * 0.6f + 0.6f * 0.4f;
				color.z = color.z * 0.6f + 0.1f * 0.4f;
				// Cyan override for hovered body (from matrix tooltip)
				if (idx == g_ldl_hover_body) {
					color = V3(0.0f, 0.9f, 1.0f);
				}
			}
		}
		draw_body_mesh(e->mesh, e->body, e->scale, color);
	}

	// Draw joint lines for selected island
	if (g_ldl_inspect_island >= 0 && (draw_w->island_gen[g_ldl_inspect_island] & 1)) {
		Island* isl = &draw_w->islands[g_ldl_inspect_island];
		int ji = isl->head_joint;
		while (ji >= 0) {
			JointInternal* j = &draw_w->joints[ji];
			int ba = j->body_a, bb = j->body_b;
			v3 pa = body_pos(draw_w, ba);
			v3 pb = body_pos(draw_w, bb);
			render_debug_line(pa, pb, V3(1.0f, 0.7f, 0.2f));
			ji = j->island_next;
		}
	}

	// Generic joint debug rendering
	if (g_show_joints) {
		world_debug_joints(g_world, draw_joint_debug, NULL);
	}

	// Mouse constraint visual
	if (g_mouse_body.id) {
		WorldInternal* w = (WorldInternal*)g_world.id;
		int anchor_idx = handle_index(g_mouse_anchor);
		v3 target = body_pos(w, anchor_idx);
		v3 body_pos = body_get_position(g_world, g_mouse_body);
		quat body_rot = body_get_rotation(g_world, g_mouse_body);
		v3 attach = add(body_pos, rotate(body_rot, g_mouse_local_hit));
		render_debug_line(target, attach, V3(1.0f, 0.5f, 0.0f));
	}

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

	if (g_show_bvh) {
		world_debug_bvh(g_world, bvh_debug_draw_cb, NULL);
		// Validate BVH leaves every frame when debug draw is on
		WorldInternal* vw = (WorldInternal*)g_world.id;
		int stale = bvh_validate_leaves(vw->bvh_dynamic, vw);
		if (stale > 0) printf("[BVH] %d stale dynamic leaves at draw time (frame %d)\n", stale, vw->frame);
	}

	if (g_show_proxies) {
		// Draw per-body fat AABB (green) and tight AABB (yellow) for dynamic bodies.
		// Stale proxies (tight outside fat) show as red tight AABB.
		WorldInternal* pw = (WorldInternal*)g_world.id;
		for (int i = 0; i < asize(pw->body_hot); i++) {
			if (!split_alive(pw->body_gen, i)) continue;
			if (asize(pw->body_cold[i].shapes) == 0) continue;
			if (body_inv_mass(pw, i) == 0.0f) continue;
			int leaf = pw->body_cold[i].bvh_leaf;
			if (leaf < 0) continue;
			AABB tight_bb = body_aabb(&pw->body_state[i], &pw->body_cold[i]);
			v3 fmin = pw->bvh_dynamic->leaves[leaf].fat_min;
			v3 fmax = pw->bvh_dynamic->leaves[leaf].fat_max;
			int stale = tight_bb.min.x < fmin.x || tight_bb.min.y < fmin.y || tight_bb.min.z < fmin.z || tight_bb.max.x > fmax.x || tight_bb.max.y > fmax.y || tight_bb.max.z > fmax.z;
			draw_aabb_wireframe(fmin, fmax, V3(0.2f, 0.8f, 0.2f)); // fat = green
			draw_aabb_wireframe(tight_bb.min, tight_bb.max, stale ? V3(1, 0.2f, 0.2f) : V3(0.9f, 0.9f, 0.3f)); // tight = yellow (red if stale)
		}
	}

	render_end();
}
