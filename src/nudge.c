// See LICENSE for licensing info.
// nudge.c -- physics world implementation

#include "nudge_internal.h"
#include "gjk.c"
#include "quickhull.c"
#include "bvh.c"
#include "collision.c"
#include "inertia.c"
#include "solver.c"
#include "joints.c"
#include "islands.c"
#include "solver_avbd.c"

// -----------------------------------------------------------------------------
// World.

World create_world(WorldParams params)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	w->gravity = params.gravity;
	w->broadphase_type = params.broadphase;
	w->friction_model = params.friction_model;
	w->solver_type = params.solver_type;
	w->sleep_enabled = 1;
	w->velocity_iters = params.velocity_iters > 0 ? params.velocity_iters : SOLVER_VELOCITY_ITERS;
	w->position_iters = params.position_iters > 0 ? params.position_iters : SOLVER_POSITION_ITERS;
	w->contact_hertz = params.contact_hertz > 0.0f ? params.contact_hertz : 60.0f;
	w->contact_damping_ratio = params.contact_damping_ratio > 0.0f ? params.contact_damping_ratio : 3.0f;
	w->max_push_velocity = params.max_push_velocity > 0.0f ? params.max_push_velocity : 3.0f;
	w->sub_steps = params.sub_steps > 0 ? params.sub_steps : 4;
	w->avbd_alpha = 0.99f;
	w->avbd_beta_lin = 10000.0f;
	w->avbd_beta_ang = 100.0f;
	w->avbd_gamma = 0.999f;
	w->avbd_iterations = 20;
	w->bvh_static = CK_ALLOC(sizeof(BVHTree));
	w->bvh_dynamic = CK_ALLOC(sizeof(BVHTree));
	bvh_init(w->bvh_static);
	bvh_init(w->bvh_dynamic);
	return (World){ (uint64_t)w };
}

void destroy_world(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	for (int i = 0; i < asize(w->body_cold); i++) {
		afree(w->body_cold[i].shapes);
	}
	afree(w->debug_contacts);
	map_free(w->warm_cache);
	map_free(w->avbd_warm_cache);
	afree(w->avbd_prev_velocity);
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	CK_FREE(w);
}

// Integrate velocities for a sub-step (gravity + damping).
static void integrate_velocities(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
		if (h->linear_damping > 0.0f)
			h->velocity = scale(h->velocity, 1.0f / (1.0f + h->linear_damping * dt));
		if (h->angular_damping > 0.0f)
			h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + h->angular_damping * dt));
	}
}

// Integrate positions and rotations for a sub-step.
static void integrate_positions(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;

		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL)
			h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL)
			h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));

		h->position = add(h->position, scale(h->velocity, dt));

		h->angular_velocity = solve_gyroscopic(h->rotation, h->inv_inertia_local, h->angular_velocity, dt);

		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, h->rotation);
		h->rotation.x += 0.5f * dt * dq.x;
		h->rotation.y += 0.5f * dt * dq.y;
		h->rotation.z += 0.5f * dt * dq.z;
		h->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y
			+ h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql;
		h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
	}
}

void world_step(World world, float dt)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->frame++;
	int n_sub = w->sub_steps;
	float sub_dt = dt / (float)n_sub;

	// Age warm cache once per frame (AVBD does its own in avbd_solve)
	if (w->solver_type != SOLVER_AVBD)
		warm_cache_age_and_evict(w);

	// --- Collision detection (once per frame) ---
	// Dual solvers need velocity integration before collision for first substep.
	// AVBD handles gravity internally via inertial position — skip.
	if (w->solver_type != SOLVER_AVBD)
		integrate_velocities(w, sub_dt);

	CK_DYNA InternalManifold* manifolds = NULL;
	broadphase_and_collide(w, &manifolds);
	islands_update_contacts(w, manifolds, asize(manifolds));

	aclear(w->debug_contacts);
	for (int i = 0; i < asize(manifolds); i++)
		for (int c = 0; c < manifolds[i].m.count; c++)
			apush(w->debug_contacts, manifolds[i].m.contacts[c]);

	int manifold_count = asize(manifolds);

	// AVBD takes a completely different path (primal-dual position solver)
	if (w->solver_type == SOLVER_AVBD) {
		avbd_solve(w, manifolds, manifold_count, dt);
		if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
		afree(manifolds);
		return;
	}

	// --- Pre-solve (once per frame, using sub_dt for softness/bias) ---
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, sub_dt);

	SolverBallSocket* sol_bs = NULL;
	SolverDistance*    sol_dist = NULL;
	joints_pre_solve(w, sub_dt, &sol_bs, &sol_dist);
	joints_warm_start(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	// --- Graph color (once per frame) ---
	int count = asize(w->body_hot);
	CK_DYNA ConstraintRef* crefs = NULL;
	int sm_count = asize(sm);
	for (int i = 0; i < sm_count; i++) {
		ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
			.body_a = sm[i].body_a, .body_b = sm[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_bs); i++) {
		ConstraintRef r = { .type = CTYPE_BALL_SOCKET, .index = i,
			.body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_dist); i++) {
		ConstraintRef r = { .type = CTYPE_DISTANCE, .index = i,
			.body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b };
		apush(crefs, r);
	}

	int cref_count = asize(crefs);
	int batch_starts[65] = {0};
	int color_count = 0;
	if (cref_count > 0)
		color_constraints(crefs, cref_count, count, batch_starts, &color_count);

	// --- Sub-step loop: velocity solve + position integrate ---
	// First sub-step: velocities already integrated above.
	// Joint bias applied on first sub-step only (error is stale after position integration).
	for (int sub = 0; sub < n_sub; sub++) {
		if (sub > 0)
			integrate_velocities(w, sub_dt);

		for (int iter = 0; iter < w->velocity_iters; iter++)
			for (int c = 0; c < color_count; c++)
				for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
					solve_constraint(w, &crefs[i], sm, sc, sol_bs, sol_dist);

		integrate_positions(w, sub_dt);

		// Relax contacts: refresh separation/bias from updated positions
		if (w->solver_type == SOLVER_SOFT_STEP || w->solver_type == SOLVER_BLOCK)
			solver_relax_contacts(w, sm, asize(sm), sc, sub_dt);

		// Zero rigid joint bias after first sub-step (position error is stale)
		if (sub == 0) {
			for (int i = 0; i < asize(sol_bs); i++)
				if (sol_bs[i].softness == 0.0f) sol_bs[i].bias = V3(0, 0, 0);
			for (int i = 0; i < asize(sol_dist); i++)
				if (sol_dist[i].softness == 0.0f) sol_dist[i].bias = 0.0f;
		}
	}

	afree(crefs);

	// Position correction: NGS for hard SI; also for soft modes when contact_hertz is off
	if (w->solver_type == SOLVER_SI)
		solver_position_correct(w, sm, asize(sm), sc);
	else if (w->contact_hertz <= 0.0f)
		solver_position_correct(w, sm, asize(sm), sc);

	// Post-solve (once per frame)
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	if (w->sleep_enabled) islands_evaluate_sleep(w, dt);

	afree(manifolds);
}

void world_set_friction_model(World world, FrictionModel model)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->friction_model = model;
}

void world_set_solver_type(World world, SolverType type)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->solver_type = type;
}

// -----------------------------------------------------------------------------
// Body.

Body create_body(World world, BodyParams params)
{
	assert(is_valid(params.position) && "create_body: position is NaN/inf");
	assert(is_valid(params.rotation) && "create_body: rotation is NaN/inf");
	assert(is_valid(params.mass) && params.mass >= 0.0f && "create_body: mass must be >= 0 and finite");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	split_add(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);

	w->body_cold[idx] = (BodyCold){
		.mass = params.mass,
		.shapes = NULL,
		.bvh_leaf = -1,
		.island_id = -1,
		.island_prev = -1,
		.island_next = -1,
	};
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	float ang_damp = params.angular_damping;
	if (ang_damp == 0.0f) ang_damp = 0.03f; // default: 3%/s (BEPU-style)
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
		.friction = fric,
		.restitution = params.restitution,
		.linear_damping = params.linear_damping,
		.angular_damping = ang_damp,
	};

	return split_handle(Body, w->body_gen, idx);
}

void destroy_body(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	// Remove from island
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl)) {
		// Remove all joints connected to this body
		int ji = w->islands[isl].head_joint;
		while (ji >= 0) {
			int next = w->joints[ji].island_next;
			if (w->joints[ji].body_a == idx || w->joints[ji].body_b == idx) {
				island_remove_joint(w, isl, ji);
				w->islands[isl].constraint_remove_count++;
			}
			ji = next;
		}
		island_remove_body(w, isl, idx);
		w->islands[isl].constraint_remove_count++;
	}
	if (w->body_cold[idx].bvh_leaf >= 0) {
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		int moved_body = bvh_remove(tree, w->body_cold[idx].bvh_leaf);
		if (moved_body >= 0) w->body_cold[moved_body].bvh_leaf = w->body_cold[idx].bvh_leaf;
	}
	afree(w->body_cold[idx].shapes);
	split_del(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);
}

void body_add_shape(World world, Body body, ShapeParams params)
{
	assert(is_valid(params.local_pos) && "body_add_shape: local_pos is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));

	ShapeInternal s = {0};
	s.type = params.type;
	s.local_pos = params.local_pos;
	switch (params.type) {
	case SHAPE_SPHERE:  s.sphere.radius = params.sphere.radius; break;
	case SHAPE_CAPSULE: s.capsule.half_height = params.capsule.half_height;
	                    s.capsule.radius = params.capsule.radius; break;
	case SHAPE_BOX:     s.box.half_extents = params.box.half_extents; break;
	case SHAPE_HULL:    s.hull.hull = params.hull.hull;
	                    s.hull.scale = params.hull.scale; break;
	}
	apush(w->body_cold[idx].shapes, s);
	recompute_body_inertia(w, idx);

	// Insert into BVH on first shape add.
	if (w->broadphase_type == BROADPHASE_BVH && asize(w->body_cold[idx].shapes) == 1) {
		AABB box = aabb_expand(body_aabb(&w->body_hot[idx], &w->body_cold[idx]), BVH_AABB_MARGIN);
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		w->body_cold[idx].bvh_leaf = bvh_insert(tree, idx, box);
	}
}

v3 body_get_position(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].position;
}

quat body_get_rotation(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].rotation;
}

void body_wake(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_velocity(World world, Body body, v3 vel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].velocity = vel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_angular_velocity(World world, Body body, v3 avel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].angular_velocity = avel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

int body_is_asleep(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl < 0 || !island_alive(w, isl)) return 0;
	return !w->islands[isl].awake;
}

// -----------------------------------------------------------------------------
// Joints.

Joint create_ball_socket(World world, BallSocketParams params)
{
	assert(is_valid(params.local_offset_a) && "create_ball_socket: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_ball_socket: local_offset_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	// Grow joint arrays manually (no split_add -- joints don't need hot/cold split)
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1); // odd = alive
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = ba, .body_b = bb,
		.ball_socket = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_distance(World world, DistanceParams params)
{
	assert(is_valid(params.local_offset_a) && "create_distance: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_distance: local_offset_b is NaN/inf");
	assert(is_valid(params.rest_length) && "create_distance: rest_length is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	// Auto-compute rest length if not specified
	float rest = params.rest_length;
	if (rest <= 0.0f) {
		BodyHot* a = &w->body_hot[ba];
		BodyHot* b = &w->body_hot[bb];
		v3 wa = add(a->position, rotate(a->rotation, params.local_offset_a));
		v3 wb = add(b->position, rotate(b->rotation, params.local_offset_b));
		rest = len(sub(wb, wa));
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_DISTANCE,
		.body_a = ba, .body_b = bb,
		.distance = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.rest_length = rest,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	unlink_joint_from_island(w, idx);
	memset(&w->joints[idx], 0, sizeof(JointInternal));
	w->joint_gen[idx]++; // even = dead
	apush(w->joint_free, idx);
}

static void bvh_debug_walk(BVHTree* t, int ni, int depth, BVHDebugFn fn, void* user)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		fn(c->min, c->max, depth, bvh_child_is_leaf(c), user);
		if (bvh_child_is_internal(c)) bvh_debug_walk(t, c->index, depth + 1, fn, user);
	}
}

void world_debug_bvh(World world, BVHDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (w->bvh_dynamic->root >= 0) bvh_debug_walk(w->bvh_dynamic, w->bvh_dynamic->root, 0, fn, user);
	if (w->bvh_static->root >= 0) bvh_debug_walk(w->bvh_static, w->bvh_static->root, 0, fn, user);
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}
