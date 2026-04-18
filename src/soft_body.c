// soft_body.c -- Native soft-body solver (XPBD-style PGS).
//
// Each SoftBody owns an independent particle network (3-DOF nodes + 1-DOF
// scalar distance links). The solver is projected Gauss-Seidel: every
// substep, for each PGS iteration, every link computes a delta-lambda
// impulse from the current velocity/position error and applies it to its
// two endpoints. Warm-start lambdas persist across iterations, substeps,
// and frames.
//
// Why PGS and not LDL: soft bodies have user-authored topology that is
// routinely over-constrained (e.g. surface mesh + interior bracing). A
// direct LDL solve on rank-deficient K amplifies null-space lambdas by
// 1 / compliance and blows up on contact. PGS doesn't care -- every
// constraint is solved independently, and redundant constraints just
// converge slightly slower.
//
// For a distance link with axis a = (p_j - p_i) / |p_j - p_i| and nodes
// (i, j), the velocity-level constraint is dC/dt = a . (v_j - v_i) = 0.
// PGS iteration for one link:
//     rhs = -jv - bias - softness * lambda_acc
//     delta_lambda = rhs / (inv_m_i + inv_m_j + softness)
//     lambda_acc += delta_lambda
//     v_i -= axis * inv_m_i * delta_lambda
//     v_j += axis * inv_m_j * delta_lambda
//
// Pins set effective inv_mass to 0. Pinned nodes don't feel impulses and
// don't integrate; their position is snapped to the pin target each
// substep so the user can drag a pin by updating its world_pos.
//
// v1 scope: internal dynamics + static pins + node-as-sphere collision
// against static rigid bodies. No dynamic-rigid coupling, no snapshot.

// Diagnostics: updated each soft_body_step_world, read by bench/debug UI.
double g_soft_body_max_lambda;
double g_soft_body_max_rhs;

#define SB_DEFAULT_ITERATIONS 12
#define SB_BAUMGARTE_GAIN     0.2f  // beta for rigid links' position drift correction
#define SB_MAX_PUSH_VEL       3.0f  // per-substep clamp on Baumgarte bias magnitude (m/s)

// -----------------------------------------------------------------------------
// Handle <-> pointer resolution.

static SoftBodyInternal* sb_lookup(WorldInternal* w, SoftBody sb)
{
	int idx = handle_index(sb);
	assert(idx >= 0 && idx < asize(w->soft_body_gen));
	assert(w->soft_body_gen[idx] == handle_gen(sb));
	return &w->soft_bodies[idx];
}

static SoftBodyInternal* sb_lookup_safe(WorldInternal* w, SoftBody sb)
{
	int idx = handle_index(sb);
	if (idx < 0 || idx >= asize(w->soft_body_gen)) return NULL;
	if (w->soft_body_gen[idx] != handle_gen(sb)) return NULL;
	return &w->soft_bodies[idx];
}

static int sb_find_pin(SoftBodyInternal* sb, int node)
{
	int pc = asize(sb->pins);
	for (int i = 0; i < pc; i++) if (sb->pins[i].node == node) return i;
	return -1;
}

static void sb_free_storage(SoftBodyInternal* sb)
{
	afree(sb->node_pos);
	afree(sb->node_vel);
	afree(sb->node_ext_force);
	afree(sb->node_inv_mass);
	afree(sb->node_user_inv_mass);
	afree(sb->links);
	afree(sb->pins);
	memset(sb, 0, sizeof(*sb));
}

// -----------------------------------------------------------------------------
// Lifecycle.

SoftBody create_soft_body(World world, SoftBodyParams params)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	if (asize(w->soft_body_free) > 0) {
		idx = apop(w->soft_body_free);
		w->soft_body_gen[idx]++;
		memset(&w->soft_bodies[idx], 0, sizeof(SoftBodyInternal));
	} else {
		idx = asize(w->soft_bodies);
		apush(w->soft_bodies, (SoftBodyInternal){0});
		apush(w->soft_body_gen, 1);
	}
	SoftBodyInternal* sb = &w->soft_bodies[idx];
	sb->params = params;
	if (sb->params.linear_damping <= 0.0f) sb->params.linear_damping = 0.02f;
	if (sb->params.iterations <= 0) sb->params.iterations = SB_DEFAULT_ITERATIONS;
	return (SoftBody){ handle_make(idx, w->soft_body_gen[idx]) };
}

void destroy_soft_body(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup_safe(w, sb_h);
	if (!sb) return;
	int idx = handle_index(sb_h);
	sb_free_storage(sb);
	w->soft_body_gen[idx]++;
	apush(w->soft_body_free, idx);
}

int soft_body_is_valid(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(sb_h);
	if (idx < 0 || idx >= asize(w->soft_body_gen)) return 0;
	return w->soft_body_gen[idx] == handle_gen(sb_h);
}

// -----------------------------------------------------------------------------
// Assembly.

int soft_body_add_node(World world, SoftBody sb_h, v3 pos, float mass)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(!sb->built);
	int idx = asize(sb->node_pos);
	apush(sb->node_pos, pos);
	apush(sb->node_vel, V3(0, 0, 0));
	apush(sb->node_ext_force, V3(0, 0, 0));
	float inv_m = (mass > 0.0f) ? (1.0f / mass) : 0.0f;
	apush(sb->node_inv_mass, inv_m);
	apush(sb->node_user_inv_mass, inv_m);
	return idx;
}

void soft_body_add_link(World world, SoftBody sb_h, int node_i, int node_j, float rest_length, SpringParams spring)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(!sb->built);
	assert(node_i != node_j);
	int nc = asize(sb->node_pos);
	assert(node_i >= 0 && node_i < nc && node_j >= 0 && node_j < nc);
	if (rest_length < 0.0f) {
		v3 d = sub(sb->node_pos[node_j], sb->node_pos[node_i]);
		rest_length = len(d);
	}
	// v1: per-link SpringParams reserved -- compliance is taken from
	// SoftBodyParams.default_spring for the whole body.
	(void)spring;
	SoftLink lk = (SoftLink){
		.node_i = node_i, .node_j = node_j,
		.rest_length = rest_length,
		.axis = V3(0, 1, 0),
		.lambda = 0.0f,
	};
	apush(sb->links, lk);
}

void soft_body_build(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(!sb->built);
	sb->built = 1;
}

void soft_body_pin_static(World world, SoftBody sb_h, int node, v3 world_pos)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(node >= 0 && node < asize(sb->node_pos));
	int pi = sb_find_pin(sb, node);
	if (pi < 0) {
		SoftPin p = (SoftPin){ .node = node, .world_pos = world_pos };
		apush(sb->pins, p);
	} else {
		sb->pins[pi].world_pos = world_pos;
	}
	sb->node_pos[node] = world_pos;
	sb->node_vel[node] = V3(0, 0, 0);
	sb->node_inv_mass[node] = 0.0f;
}

void soft_body_unpin(World world, SoftBody sb_h, int node)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	int pi = sb_find_pin(sb, node);
	if (pi < 0) return;
	int last = asize(sb->pins) - 1;
	if (pi != last) sb->pins[pi] = sb->pins[last];
	apop(sb->pins);
	sb->node_inv_mass[node] = sb->node_user_inv_mass[node];
}

// -----------------------------------------------------------------------------
// Gameplay pokes.

void soft_body_apply_force(World world, SoftBody sb_h, int node, v3 force)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(node >= 0 && node < asize(sb->node_ext_force));
	sb->node_ext_force[node] = add(sb->node_ext_force[node], force);
}

void soft_body_apply_impulse(World world, SoftBody sb_h, int node, v3 impulse)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(node >= 0 && node < asize(sb->node_vel));
	if (sb->node_inv_mass[node] == 0.0f) return;
	sb->node_vel[node] = add(sb->node_vel[node], scale(impulse, sb->node_inv_mass[node]));
}

// -----------------------------------------------------------------------------
// Readback.

int soft_body_node_count(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	return asize(sb->node_pos);
}

const v3* soft_body_node_positions(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	return sb->node_pos;
}

const v3* soft_body_node_velocities(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	return sb->node_vel;
}

int soft_body_link_count(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	return asize(sb->links);
}

int soft_body_get_links(World world, SoftBody sb_h, int* out_pairs, int max)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	int L = asize(sb->links);
	int n = L < max ? L : max;
	for (int k = 0; k < n; k++) {
		out_pairs[2 * k + 0] = sb->links[k].node_i;
		out_pairs[2 * k + 1] = sb->links[k].node_j;
	}
	return L;
}

int world_get_soft_body_count(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int n = 0;
	for (int i = 0; i < asize(w->soft_body_gen); i++) if (split_alive(w->soft_body_gen, i)) n++;
	return n;
}

int world_get_soft_bodies(World world, SoftBody* out, int max)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int n = 0, total = 0;
	for (int i = 0; i < asize(w->soft_body_gen); i++) {
		if (!split_alive(w->soft_body_gen, i)) continue;
		if (n < max) out[n++] = (SoftBody){ handle_make(i, w->soft_body_gen[i]) };
		total++;
	}
	return total;
}

// -----------------------------------------------------------------------------
// PBD-style node-vs-static-rigid collision. v1: static bodies only, SHAPE_BOX
// and SHAPE_SPHERE. Each soft-body node is treated as a sphere of
// sb->params.node_radius. Corrections applied inline -- no PGS contact
// constraint, no manifold cache, no ContactSummary. Dynamic rigid bodies are
// ignored since proper coupling needs a K-matrix border block (phase 2+).
static void sb_collide_node_with_world(WorldInternal* w, SoftBodyInternal* sb, int node_idx)
{
	if (sb->node_inv_mass[node_idx] == 0.0f) return; // pinned don't collide
	float r = sb->params.node_radius;
	if (r <= 0.0f) return;

	Sphere node_sphere = { .center = sb->node_pos[node_idx], .radius = r };
	int bc = asize(w->body_hot);
	for (int bi = 0; bi < bc; bi++) {
		if (!split_alive(w->body_gen, bi)) continue;
		if (w->body_hot[bi].inv_mass > 0.0f) continue; // v1: static only
		BodyState* bs = &w->body_state[bi];
		BodyCold* bcold = &w->body_cold[bi];
		int sc = asize(bcold->shapes);
		for (int si = 0; si < sc; si++) {
			ShapeInternal* sh = &bcold->shapes[si];
			v3 w_pos = add(bs->position, rotate(bs->rotation, sh->local_pos));
			quat w_rot = quat_mul(bs->rotation, sh->local_rot);
			Manifold m = (Manifold){ 0 };
			int hit = 0;
			switch (sh->type) {
				case SHAPE_BOX: {
					Box b = { .center = w_pos, .rotation = w_rot, .half_extents = sh->box.half_extents };
					hit = collide_sphere_box(node_sphere, b, &m);
					break;
				}
				case SHAPE_SPHERE: {
					Sphere s = { .center = w_pos, .radius = sh->sphere.radius };
					hit = collide_sphere_sphere(node_sphere, s, &m);
					break;
				}
				default: break;
			}
			if (!hit) continue;
			// Normal points from A (node) toward B (world body). If the node
			// penetrates, push it along -normal by penetration; cancel the
			// inward (into-body) velocity component; apply Coulomb-clamped
			// tangential damping.
			for (int ci = 0; ci < m.count; ci++) {
				Contact* c = &m.contacts[ci];
				sb->node_pos[node_idx] = sub(sb->node_pos[node_idx], scale(c->normal, c->penetration));
				node_sphere.center = sb->node_pos[node_idx];
				v3 v = sb->node_vel[node_idx];
				float vn = dot(v, c->normal);
				if (vn <= 0.0f) continue; // already separating
				v3 vt = sub(v, scale(c->normal, vn));
				float mu = 0.4f;
				float max_friction = mu * vn;
				float vt_len = len(vt);
				v3 new_v;
				if (vt_len > max_friction && vt_len > 1e-6f) {
					new_v = scale(vt, (vt_len - max_friction) / vt_len);
				} else {
					new_v = V3(0, 0, 0);
				}
				sb->node_vel[node_idx] = new_v;
			}
		}
	}
}

// -----------------------------------------------------------------------------
// Solver substep.

static void soft_body_substep(WorldInternal* w, SoftBodyInternal* sb, v3 gravity, float dt)
{
	if (!sb->built) return;
	int N = asize(sb->node_pos);
	int L = asize(sb->links);

	// 1) Integrate velocity: gravity + ext_force + damping.
	float damp_mul = 1.0f / (1.0f + sb->params.linear_damping * dt);
	for (int n = 0; n < N; n++) {
		if (sb->node_inv_mass[n] == 0.0f) { sb->node_ext_force[n] = V3(0, 0, 0); continue; }
		v3 a = add(gravity, scale(sb->node_ext_force[n], sb->node_inv_mass[n]));
		sb->node_vel[n] = add(sb->node_vel[n], scale(a, dt));
		sb->node_vel[n] = scale(sb->node_vel[n], damp_mul);
		sb->node_ext_force[n] = V3(0, 0, 0);
	}

	if (L > 0) {
		// 2) Refresh per-link state: axis, softness/bias, inv_eff_mass.
		//    Identical for every link in this substep (uses the one shared
		//    SoftBodyParams.default_spring). Bias is clamped up front so a
		//    sudden large position error (mouse drag, deep collision) can't
		//    produce an explosive impulse.
		SpringParams sp = sb->params.default_spring;
		float softness = 0.0f, bias_gain = 0.0f;
		if (sp.frequency > 0.0f) {
			float omega = 6.2831853f * sp.frequency;
			float zeta = sp.damping_ratio > 0.0f ? sp.damping_ratio : 1.0f;
			float hk = dt * omega * omega;
			float d = 2.0f * zeta * omega;
			float denom = dt * (d + hk);
			if (denom > 1e-12f) {
				softness = 1.0f / denom;
				bias_gain = hk * softness;
			}
		} else {
			// Rigid: Baumgarte position bias in velocity space, no compliance.
			bias_gain = SB_BAUMGARTE_GAIN / dt;
		}
		for (int k = 0; k < L; k++) {
			SoftLink* lk = &sb->links[k];
			v3 d = sub(sb->node_pos[lk->node_j], sb->node_pos[lk->node_i]);
			float l2 = len2(d);
			float length = l2 > 1e-16f ? sqrtf(l2) : 1e-8f;
			lk->axis = scale(d, 1.0f / length);
			float err = length - lk->rest_length;
			float b = bias_gain * err;
			if (b >  SB_MAX_PUSH_VEL) b =  SB_MAX_PUSH_VEL;
			if (b < -SB_MAX_PUSH_VEL) b = -SB_MAX_PUSH_VEL;
			lk->bias = b;
			lk->softness = softness;
			float mi = sb->node_inv_mass[lk->node_i];
			float mj = sb->node_inv_mass[lk->node_j];
			float denom = mi + mj + softness;
			lk->inv_eff_mass = denom > 1e-10f ? (1.0f / denom) : 0.0f;
		}

		// 3) Warm-start: apply the previously-converged lambda as an initial
		//    impulse. Subsequent iterations compute delta-lambdas that adjust
		//    from this warm state.
		for (int k = 0; k < L; k++) {
			SoftLink* lk = &sb->links[k];
			if (lk->inv_eff_mass == 0.0f) continue;
			v3 imp = scale(lk->axis, lk->lambda);
			float mi = sb->node_inv_mass[lk->node_i];
			float mj = sb->node_inv_mass[lk->node_j];
			sb->node_vel[lk->node_i] = sub(sb->node_vel[lk->node_i], scale(imp, mi));
			sb->node_vel[lk->node_j] = add(sb->node_vel[lk->node_j], scale(imp, mj));
		}

		// 4) PGS iteration. Each iteration: for each link, compute the impulse
		//    that drives its constraint toward satisfaction, apply it, and
		//    accumulate the total into lk->lambda.
		int iters = sb->params.iterations;
		double local_max_lambda = 0, local_max_rhs = 0;
		for (int it = 0; it < iters; it++) {
			for (int k = 0; k < L; k++) {
				SoftLink* lk = &sb->links[k];
				if (lk->inv_eff_mass == 0.0f) continue; // both endpoints pinned
				v3 va = sb->node_vel[lk->node_i];
				v3 vb = sb->node_vel[lk->node_j];
				float jv = dot(sub(vb, va), lk->axis);
				float rhs = -jv - lk->bias - lk->softness * lk->lambda;
				float delta_lambda = rhs * lk->inv_eff_mass;
				lk->lambda += delta_lambda;
				v3 imp = scale(lk->axis, delta_lambda);
				float mi = sb->node_inv_mass[lk->node_i];
				float mj = sb->node_inv_mass[lk->node_j];
				sb->node_vel[lk->node_i] = sub(sb->node_vel[lk->node_i], scale(imp, mi));
				sb->node_vel[lk->node_j] = add(sb->node_vel[lk->node_j], scale(imp, mj));
				float la = lk->lambda < 0 ? -lk->lambda : lk->lambda;
				float rm = rhs < 0 ? -rhs : rhs;
				if (la > local_max_lambda) local_max_lambda = la;
				if (rm > local_max_rhs) local_max_rhs = rm;
			}
		}
		if (local_max_lambda > g_soft_body_max_lambda) g_soft_body_max_lambda = local_max_lambda;
		if (local_max_rhs > g_soft_body_max_rhs) g_soft_body_max_rhs = local_max_rhs;
	}

	// 5) Integrate positions.
	for (int n = 0; n < N; n++) {
		if (sb->node_inv_mass[n] == 0.0f) continue;
		sb->node_pos[n] = add(sb->node_pos[n], scale(sb->node_vel[n], dt));
	}

	// 6) Collide each node against world static rigid bodies.
	for (int n = 0; n < N; n++) {
		sb_collide_node_with_world(w, sb, n);
	}

	// 7) Snap pinned nodes to pin target. Updating the pin target via
	//    soft_body_pin_static is the drag mechanism.
	for (int p = 0; p < asize(sb->pins); p++) {
		sb->node_pos[sb->pins[p].node] = sb->pins[p].world_pos;
		sb->node_vel[sb->pins[p].node] = V3(0, 0, 0);
	}
}

// -----------------------------------------------------------------------------
// World-level entry points.

static void soft_body_step_world(WorldInternal* w, float dt)
{
	g_soft_body_max_lambda = 0.0;
	g_soft_body_max_rhs = 0.0;
	int count = asize(w->soft_bodies);
	if (count == 0) return;
	int n_sub = w->sub_steps > 0 ? w->sub_steps : 1;
	float sub_dt = dt / (float)n_sub;
	for (int s = 0; s < n_sub; s++) {
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->soft_body_gen, i)) continue;
			soft_body_substep(w, &w->soft_bodies[i], w->gravity, sub_dt);
		}
	}
}

static void soft_body_free_all(WorldInternal* w)
{
	int count = asize(w->soft_bodies);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->soft_body_gen, i)) continue;
		sb_free_storage(&w->soft_bodies[i]);
	}
	afree(w->soft_bodies);
	afree(w->soft_body_gen);
	afree(w->soft_body_free);
}
