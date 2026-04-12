// nudge_dll.c -- unity build for nudge as a shared library.
// Thin C wrapper with flat parameter lists for easy P/Invoke / FFI.

#define CKIT_IMPLEMENTATION
#include "ckit.h"
#include "nudge.h"
#include "nudge.c"

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

// --- World ---

EXPORT uint64_t nudge_create_world(float gx, float gy, float gz, int solver_type, int sub_steps, int velocity_iters) {
	WorldParams p = {0};
	p.gravity = V3(gx, gy, gz);
	p.broadphase = BROADPHASE_BVH;
	p.solver_type = solver_type;
	p.sub_steps = sub_steps > 0 ? sub_steps : 2;
	p.velocity_iters = velocity_iters > 0 ? velocity_iters : 8;
	return create_world(p).id;
}

EXPORT void nudge_destroy_world(uint64_t world) {
	destroy_world((World){world});
}

EXPORT void nudge_step(uint64_t world, float dt) {
	world_step((World){world}, dt);
}

EXPORT double nudge_get_step_time(uint64_t world) {
	return world_get_perf((World){world}).total;
}

// --- Body ---

EXPORT uint64_t nudge_create_body(uint64_t world, float px, float py, float pz, float mass, float friction, float restitution) {
	BodyParams p = {0};
	p.position = V3(px, py, pz);
	p.rotation = (quat){0, 0, 0, 1};
	p.mass = mass;
	p.friction = friction > 0 ? friction : 0.5f;
	p.restitution = restitution;
	return create_body((World){world}, p).id;
}

EXPORT void nudge_body_add_box(uint64_t world, uint64_t body, float hx, float hy, float hz) {
	ShapeParams s = {0};
	s.type = SHAPE_BOX;
	s.box.half_extents = V3(hx, hy, hz);
	body_add_shape((World){world}, (Body){body}, s);
}

EXPORT void nudge_body_add_sphere(uint64_t world, uint64_t body, float radius) {
	ShapeParams s = {0};
	s.type = SHAPE_SPHERE;
	s.sphere.radius = radius;
	body_add_shape((World){world}, (Body){body}, s);
}

EXPORT void nudge_body_add_capsule(uint64_t world, uint64_t body, float half_height, float radius) {
	ShapeParams s = {0};
	s.type = SHAPE_CAPSULE;
	s.capsule.half_height = half_height;
	s.capsule.radius = radius;
	body_add_shape((World){world}, (Body){body}, s);
}

EXPORT void nudge_get_position(uint64_t world, uint64_t body, float* out) {
	v3 p = body_get_position((World){world}, (Body){body});
	out[0] = p.x; out[1] = p.y; out[2] = p.z;
}

EXPORT void nudge_get_rotation(uint64_t world, uint64_t body, float* out) {
	quat q = body_get_rotation((World){world}, (Body){body});
	out[0] = q.x; out[1] = q.y; out[2] = q.z; out[3] = q.w;
}

EXPORT int nudge_body_is_asleep(uint64_t world, uint64_t body) {
	return body_is_asleep((World){world}, (Body){body});
}

EXPORT void nudge_set_sleep_enabled(uint64_t world, int enabled) {
	world_set_sleep_enabled((World){world}, enabled);
}

EXPORT int nudge_get_sleep_enabled(uint64_t world) {
	return world_get_sleep_enabled((World){world});
}

EXPORT void nudge_body_set_sleep_allowed(uint64_t world, uint64_t body, int allowed) {
	body_set_sleep_allowed((World){world}, (Body){body}, allowed);
}

EXPORT void nudge_get_perf(uint64_t world, double* out) {
	PerfTimers t = world_get_perf((World){world});
	// Top-level phases [0..6]
	out[0] = t.broadphase;
	out[1] = t.pre_solve;
	out[2] = t.pgs_solve;
	out[3] = t.position_correct;
	out[4] = t.integrate;
	out[5] = t.islands;
	out[6] = t.total;
	// PGS sub-timers [7..16]
	out[7]  = t.pgs.pre_solve;
	out[8]  = t.pgs.warm_start;
	out[9]  = t.pgs.graph_color;
	out[10] = t.pgs.iterations;
	out[11] = t.pgs.joint_limits;
	out[12] = t.pgs.ldl;
	out[13] = t.pgs.relax;
	out[14] = t.pgs.pos_contacts;
	out[15] = t.pgs.pos_joints;
	out[16] = t.pgs.post_solve;
}

EXPORT void nudge_debug_sleep(uint64_t world) {
	WorldInternal* w = (WorldInternal*)world;
	int count = asize(w->body_hot);
	int awake = 0, jittery = 0;
	float max_v2 = 0;
	int max_idx = -1;
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (w->body_hot[i].inv_mass == 0) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && (w->island_gen[isl] & 1) && w->islands[isl].awake) awake++;
		float v2 = len2(w->body_hot[i].velocity) + len2(w->body_hot[i].angular_velocity);
		if (v2 > SLEEP_VEL_THRESHOLD) jittery++;
		if (v2 > max_v2) { max_v2 = v2; max_idx = i; }
	}
	if (max_idx >= 0) {
		v3 v = w->body_hot[max_idx].velocity;
		v3 av = w->body_hot[max_idx].angular_velocity;
		fprintf(stderr, "sleep: awake=%d jittery=%d max_v2=%.6f (body %d: v=(%.4f,%.4f,%.4f) w=(%.4f,%.4f,%.4f) sleep_t=%.2f)\n", awake, jittery, max_v2, max_idx, v.x, v.y, v.z, av.x, av.y, av.z, w->body_hot[max_idx].sleep_time);
	}
}

EXPORT void nudge_body_set_velocity(uint64_t world, uint64_t body, float vx, float vy, float vz) {
	body_set_velocity((World){world}, (Body){body}, V3(vx, vy, vz));
}
