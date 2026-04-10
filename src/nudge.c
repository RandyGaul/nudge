// See LICENSE for licensing info.
// nudge.c -- physics world implementation

#include "perf.h"
#include "nudge_internal.h"
#include "gjk.c"
#include "gjk_batch.c"
#include "quickhull.c"
#include "bvh.c"
#include "collision.c"
#include "inertia.c"
#include "solver_pgs.c"
#include "joints.c"
#include "solver_ldl.c"
#include "islands.c"

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
	w->ldl_correction_iter = -2; // -2 = auto: velocity_iters/2 (mid-loop, PGS can recover after LDL)
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
	afree(w->body_vel);
	map_free(w->warm_cache);
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	for (int i = 0; i < asize(w->islands); i++) ldl_cache_free(&w->islands[i].ldl);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	map_free(w->joint_pairs);
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

static void precompute_world_inertia(WorldInternal* w)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		body_compute_inv_inertia_world(&w->body_hot[i]);
	}
}

static int perf_initialized;

void world_step(World world, float dt)
{
	if (!perf_initialized) { perf_init(); perf_initialized = 1; }
	double t_total = perf_now();

	WorldInternal* w = (WorldInternal*)world.id;
	w->frame++;
	int n_sub = w->sub_steps;
	float sub_dt = dt / (float)n_sub;

	double t0 = perf_now();
	warm_cache_age_and_evict(w);
	integrate_velocities(w, sub_dt);
	precompute_world_inertia(w);
	w->perf.integrate = perf_now() - t0;

	double t1 = perf_now();
	CK_DYNA InternalManifold* manifolds = NULL;
	broadphase_and_collide(w, &manifolds);
	islands_update_contacts(w, manifolds, asize(manifolds));
	w->perf.broadphase = perf_now() - t1;

	aclear(w->debug_contacts);
	for (int i = 0; i < asize(manifolds); i++)
		for (int c = 0; c < manifolds[i].m.count; c++)
			apush(w->debug_contacts, manifolds[i].m.contacts[c]);

	int manifold_count = asize(manifolds);

	// --- Pre-solve (once per frame, using sub_dt for softness/bias) ---
	double t2 = perf_now();
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, sub_dt);

	SolverJoint* sol_joints = NULL;
	joints_pre_solve(w, sub_dt, &sol_joints);

	// LDL is a direct solver -- rigid joints don't need warm-start. Stale
	// warm-start impulses inject energy when lever arms rotate between frames.
	// Only zero bilateral DOFs that LDL handles; preserve limit/motor DOF warm-start.
	if (w->ldl_enabled) {
		for (int i = 0; i < asize(sol_joints); i++) {
			if (sol_joints[i].softness != 0.0f) continue;
			JointInternal* j = &w->joints[sol_joints[i].joint_idx];
			int ldl_dof = sol_joints[i].dof;
			if (j->type == JOINT_HINGE && (j->hinge.limit_min != 0 || j->hinge.limit_max != 0 || j->hinge.motor_max_impulse > 0)) ldl_dof = 5;
			if (j->type == JOINT_PRISMATIC && j->prismatic.motor_max_impulse > 0) ldl_dof = 5;
			for (int d = 0; d < ldl_dof; d++) sol_joints[i].lambda[d] = 0;
		}
	}
	w->perf.pgs.pre_solve = perf_now() - t2;

	double t_ws = perf_now();
	joints_warm_start(w, sol_joints, asize(sol_joints));
	w->perf.pgs.warm_start = perf_now() - t_ws;

	// --- Graph color (once per frame) ---
	// When LDL enabled: rigid joints excluded from PGS (they get diagonal GS + LDL K^-1).
	// Soft spring joints stay in PGS (LDL handles them via softness term).
	double t_gc = perf_now();
	int count = asize(w->body_hot);
	CK_DYNA ConstraintRef* crefs = NULL;
	int sm_count = asize(sm);
	for (int i = 0; i < sm_count; i++) {
		ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
			.body_a = sm[i].body_a, .body_b = sm[i].body_b };
		apush(crefs, r);
	}
	if (!w->ldl_enabled) {
		// No LDL: all joints go into PGS
		for (int i = 0; i < asize(sol_joints); i++) {
			ConstraintRef r = { .type = CTYPE_JOINT, .index = i,
				.body_a = sol_joints[i].body_a, .body_b = sol_joints[i].body_b };
			apush(crefs, r);
		}
	} else {
		// LDL enabled: LDL handles all joints (bilateral DOFs).
		// Limit DOFs are solved separately after PGS via joints_solve_limits().
	}

	int cref_count = asize(crefs);
	int batch_starts[65] = {0};
	int color_count = 0;
	if (cref_count > 0)
		color_constraints(crefs, cref_count, count, batch_starts, &color_count);
	w->perf.pgs.graph_color = perf_now() - t_gc;

	w->perf.pre_solve = perf_now() - t2;

	// --- Sub-step loop ---
	// Unified path: PGS iterates all constraints (contacts + joints).
	// When LDL enabled, K is factored once at substep start, and a mid-loop
	// K^-1 residual correction is applied at the configured iteration.
	double t_pgs = 0, t_pos = 0, t_int_sub = 0;
	double t_jlim = 0, t_ldl = 0, t_relax = 0, t_posJ = 0;
	for (int sub = 0; sub < n_sub; sub++) {
		if (sub > 0) {
			double ti = perf_now();
			integrate_velocities(w, sub_dt);
			precompute_world_inertia(w);
			t_int_sub += perf_now() - ti;
			// Refresh joint Jacobians/limits from current body state (positions
			// changed by integrate_positions last substep, velocities just updated).
			joints_refresh_substep(w, sol_joints, asize(sol_joints), sub_dt);
		}

		int has_ldl = w->ldl_enabled && asize(sol_joints) > 0;

		// Resolve LDL correction iteration: -2 = auto (velocity_iters/2), -1 = after loop
		int ldl_iter = w->ldl_correction_iter;
		if (ldl_iter == -2) ldl_iter = w->velocity_iters / 2;

		// LDL: factor K once at start of substep (topology + numeric)
		double tl0 = perf_now();
		if (has_ldl)
			ldl_factor(w, sol_joints, asize(sol_joints), sub, sub_dt);

		// LDL: direct solve for joints (velocity correction).
		// Runs before PGS so joint impulses are already applied.
		if (has_ldl) {
			if (sub > 0) {
				for (int i = 0; i < asize(sol_joints); i++) if (sol_joints[i].softness > 0.0f) for (int d = 0; d < sol_joints[i].dof; d++) sol_joints[i].lambda[d] = 0;
			}
			ldl_velocity_correct(w, sol_joints, asize(sol_joints), sub_dt);
		}
		t_ldl += perf_now() - tl0;

		// PGS: iterate all constraints (contacts, and joints when LDL is off).
		// Fast path: patch friction contacts with no joints use compact SolverBodyVel
		// (32 bytes/body instead of 120 bytes — fits more bodies in cache).
		int use_body_vel = (w->friction_model == FRICTION_PATCH && asize(sol_joints) == 0);
		if (use_body_vel) solver_sync_vel_in(w);

		double tp = perf_now();
		if (use_body_vel) {
			for (int iter = 0; iter < w->velocity_iters; iter++) {
				for (int c = 0; c < color_count; c++)
					for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
						if (crefs[i].type == CTYPE_CONTACT)
							solve_contact_patch_sv(w->body_vel, &sm[crefs[i].index], sc);
				double tjl = perf_now();
				joints_solve_limits(w, sol_joints, asize(sol_joints));
				t_jlim += perf_now() - tjl;
			}
		} else {
			for (int iter = 0; iter < w->velocity_iters; iter++) {
				for (int c = 0; c < color_count; c++)
					for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
						solve_constraint(w, &crefs[i], sm, sc, sol_joints);
				double tjl = perf_now();
				joints_solve_limits(w, sol_joints, asize(sol_joints));
				t_jlim += perf_now() - tjl;
			}
		}
		t_pgs += perf_now() - tp;

		if (use_body_vel) solver_sync_vel_out(w);

		double ti2 = perf_now();
		integrate_positions(w, sub_dt);
		t_int_sub += perf_now() - ti2;

		// Relax contacts: refresh separation/bias from updated positions
		double tr = perf_now();
		if (w->solver_type == SOLVER_SOFT_STEP)
			solver_relax_contacts(w, sm, asize(sm), sc, sub_dt);
		t_relax += perf_now() - tr;

		// Position correction after integration.
		// LDL: direct solve projects out bulk joint error, NGS cleans up residual.
		double tl1 = perf_now();
		if (has_ldl)
			ldl_position_correct(w, sol_joints, asize(sol_joints), sub_dt);
		t_ldl += perf_now() - tl1;

		double tpj = perf_now();
		joints_position_correct(w, sol_joints, asize(sol_joints), w->position_iters);
		t_posJ += perf_now() - tpj;

		t_pos += (perf_now() - tr);
	}

	w->perf.pgs_solve = t_pgs;
	w->perf.integrate += t_int_sub;
	w->perf.position_correct = t_pos;
	w->perf.pgs.iterations = t_pgs - t_jlim;
	w->perf.pgs.joint_limits = t_jlim;
	w->perf.pgs.ldl = t_ldl;
	w->perf.pgs.relax = t_relax;
	w->perf.pgs.pos_joints = t_posJ;

	afree(crefs);

	// Position correction: NGS for hard SI; also for soft modes when contact_hertz is off
	double t3 = perf_now();
	if (w->solver_type == SOLVER_SI)
		solver_position_correct(w, sm, asize(sm), sc);
	else if (w->contact_hertz <= 0.0f)
		solver_position_correct(w, sm, asize(sm), sc);
	w->perf.pgs.pos_contacts = perf_now() - t3;
	w->perf.position_correct += perf_now() - t3;

	// Post-solve (once per frame)
	double t_ps = perf_now();
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_joints, asize(sol_joints));
	w->perf.pgs.post_solve = perf_now() - t_ps;

	double t4 = perf_now();
	if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
	w->perf.islands = perf_now() - t4;

	afree(manifolds);
	w->perf.total = perf_now() - t_total;
}

PerfTimers world_get_perf(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	return w->perf;
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
	w->ldl_topo_version++;
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
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_hinge(World world, HingeParams params)
{
	assert(is_valid(params.local_offset_a) && "create_hinge: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_hinge: local_offset_b is NaN/inf");
	assert(is_valid(params.local_axis_a) && "create_hinge: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_hinge: local_axis_b is NaN/inf");

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

	v3 axis_a_local = norm(params.local_axis_a);
	v3 axis_b_local = norm(params.local_axis_b);

	// Compute reference directions perpendicular to the hinge axis for angle measurement.
	// local_ref_a: arbitrary unit vector perpendicular to axis_a in body A's local space.
	// local_ref_b: chosen so measured angle = 0 at the initial configuration.
	v3 ref_a, ref_a_t2;
	hinge_tangent_basis(axis_a_local, &ref_a, &ref_a_t2);
	// Transform ref_a into world, then into body B's local space
	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
	v3 ref_a_world = rotate(q_a, ref_a);
	v3 ref_b = rotate(inv(q_b), ref_a_world);

	w->joints[idx] = (JointInternal){
		.type = JOINT_HINGE,
		.body_a = ba, .body_b = bb,
		.hinge = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_axis_a = axis_a_local,
			.local_axis_b = axis_b_local,
			.local_ref_a = ref_a,
			.local_ref_b = ref_b,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_fixed(World world, FixedParams params)
{
	assert(is_valid(params.local_offset_a) && "create_fixed: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_fixed: local_offset_b is NaN/inf");

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

	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
	w->joints[idx] = (JointInternal){
		.type = JOINT_FIXED,
		.body_a = ba, .body_b = bb,
		.fixed = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_rel_quat = mul(inv(q_a), q_b),
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_prismatic(World world, PrismaticParams params)
{
	assert(is_valid(params.local_offset_a) && "create_prismatic: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_prismatic: local_offset_b is NaN/inf");
	assert(is_valid(params.local_axis_a) && "create_prismatic: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_prismatic: local_axis_b is NaN/inf");

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

	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
	w->joints[idx] = (JointInternal){
		.type = JOINT_PRISMATIC,
		.body_a = ba, .body_b = bb,
		.prismatic = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.local_rel_quat = mul(inv(q_a), q_b),
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
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
	w->ldl_topo_version++;
	apush(w->joint_free, idx);
}

void joint_set_hinge_limits(World world, Joint joint, float min_angle, float max_angle)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_HINGE && "joint_set_hinge_limits: not a hinge joint");
	w->joints[idx].hinge.limit_min = min_angle;
	w->joints[idx].hinge.limit_max = max_angle;
}

void joint_set_distance_limits(World world, Joint joint, float min_distance, float max_distance)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_DISTANCE && "joint_set_distance_limits: not a distance joint");
	w->joints[idx].distance.limit_min = min_distance;
	w->joints[idx].distance.limit_max = max_distance;
}

void joint_clear_limits(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	if (w->joints[idx].type == JOINT_HINGE) {
		w->joints[idx].hinge.limit_min = 0;
		w->joints[idx].hinge.limit_max = 0;
	} else if (w->joints[idx].type == JOINT_DISTANCE) {
		w->joints[idx].distance.limit_min = 0;
		w->joints[idx].distance.limit_max = 0;
	}
}

void joint_set_hinge_motor(World world, Joint joint, float speed, float max_impulse)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_HINGE);
	w->joints[idx].hinge.motor_speed = speed;
	w->joints[idx].hinge.motor_max_impulse = max_impulse;
}

void joint_set_prismatic_motor(World world, Joint joint, float speed, float max_impulse)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_PRISMATIC);
	w->joints[idx].prismatic.motor_speed = speed;
	w->joints[idx].prismatic.motor_max_impulse = max_impulse;
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

void world_debug_joints(World world, JointDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int jcount = asize(w->joints);
	for (int i = 0; i < jcount; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];
		JointDebugInfo info = {0};
		info.type = j->type;
		if (j->type == JOINT_BALL_SOCKET) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->ball_socket.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->ball_socket.local_b));
			info.is_soft = j->ball_socket.spring.frequency > 0;
		} else if (j->type == JOINT_DISTANCE) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->distance.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->distance.local_b));
			info.is_soft = j->distance.spring.frequency > 0;
		} else if (j->type == JOINT_HINGE) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->hinge.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->hinge.local_b));
			info.axis_a = norm(rotate(a->rotation, j->hinge.local_axis_a));
			info.is_soft = j->hinge.spring.frequency > 0;
			info.motor_speed = j->hinge.motor_speed;
			info.motor_max_impulse = j->hinge.motor_max_impulse;
			info.limit_min = j->hinge.limit_min;
			info.limit_max = j->hinge.limit_max;
			info.ref_a = rotate(a->rotation, j->hinge.local_ref_a);
			info.ref_b = rotate(b->rotation, j->hinge.local_ref_b);
			float angle = atan2f(dot(cross(info.ref_a, info.ref_b), info.axis_a), dot(info.ref_a, info.ref_b));
			info.current_angle = angle;
			info.limit_active = (j->hinge.limit_min != 0 && angle <= j->hinge.limit_min) || (j->hinge.limit_max != 0 && angle >= j->hinge.limit_max);
		} else if (j->type == JOINT_FIXED) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->fixed.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->fixed.local_b));
			info.is_soft = j->fixed.spring.frequency > 0;
		} else if (j->type == JOINT_PRISMATIC) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->prismatic.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->prismatic.local_b));
			info.axis_a = norm(rotate(a->rotation, j->prismatic.local_axis_a));
			info.is_soft = j->prismatic.spring.frequency > 0;
			info.motor_speed = j->prismatic.motor_speed;
			info.motor_max_impulse = j->prismatic.motor_max_impulse;
		}
		fn(info, user);
	}
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}

// -----------------------------------------------------------------------------
// World queries.

int world_query_aabb(World world, v3 lo, v3 hi, Body* results, int max_results)
{
	WorldInternal* w = (WorldInternal*)world.id;
	AABB query = { lo, hi };
	CK_DYNA int* candidates = NULL;
	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_query_aabb(w->bvh_dynamic, query, &candidates);
		bvh_query_aabb(w->bvh_static, query, &candidates);
	} else {
		int count = asize(w->body_hot);
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->body_gen, i)) continue;
			if (asize(w->body_cold[i].shapes) == 0) continue;
			apush(candidates, i);
		}
	}
	int total = 0;
	for (int i = 0; i < asize(candidates); i++) {
		int idx = candidates[i];
		AABB b = body_aabb(&w->body_hot[idx], &w->body_cold[idx]);
		if (!aabb_overlaps(query, b)) continue;
		if (total < max_results)
			results[total] = split_handle(Body, w->body_gen, idx);
		total++;
	}
	afree(candidates);
	return total;
}

int world_raycast(World world, v3 origin, v3 direction, float max_distance, RayHit* hit)
{
	WorldInternal* w = (WorldInternal*)world.id;
	float dl = len(direction);
	if (dl < 1e-12f) return 0;
	v3 dir = scale(direction, 1.0f / dl);
	v3 inv_dir = V3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);

	CK_DYNA int* candidates = NULL;
	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_query_ray(w->bvh_dynamic, origin, inv_dir, max_distance, &candidates);
		bvh_query_ray(w->bvh_static, origin, inv_dir, max_distance, &candidates);
	} else {
		int count = asize(w->body_hot);
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->body_gen, i)) continue;
			if (asize(w->body_cold[i].shapes) == 0) continue;
			apush(candidates, i);
		}
	}

	float best_t = max_distance;
	v3 best_n = {0};
	int best_idx = -1;
	for (int i = 0; i < asize(candidates); i++) {
		int idx = candidates[i];
		float t; v3 n;
		if (ray_body(w, idx, origin, dir, best_t, &t, &n)) {
			best_t = t; best_n = n; best_idx = idx;
		}
	}
	afree(candidates);

	if (best_idx < 0) return 0;
	if (hit) {
		hit->body = split_handle(Body, w->body_gen, best_idx);
		hit->point = add(origin, scale(dir, best_t));
		hit->normal = best_n;
		hit->distance = best_t;
	}
	return 1;
}
