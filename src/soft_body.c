// soft_body.c -- Native soft-body solver.
//
// Each SoftBody owns an independent particle network (3-DOF nodes + 1-DOF
// distance links) and its own dense LDL factorization. Topology is frozen at
// soft_body_build(); the K matrix is refactored every substep because link
// axes change with positions. K is the constraint-space effective mass:
//
//     K[k,l] = J_k * M^-1 * J_l^T  (+ compliance * delta_kl)
//
// For a distance link with axis a and endpoints (i, j):
//   J_k = [ ..., -a, ..., +a, ...]  (row vector, -a at column i, +a at column j)
//
// Off-diagonal K[k,l] when links k and l share a node s:
//   sign_k(s) = -1 if s is the "first" node of k, +1 if s is the "second"
//   K[k,l] += sign_k(s) * sign_l(s) * inv_m[s] * (a_k . a_l)
//
// Pinned nodes use effective inv_mass = 0 -- they don't feel impulses and
// don't move during integrate. We still snap pinned positions each substep
// so the user can translate a pin target by calling soft_body_pin_static
// again with a new world_pos.
//
// v1 scope: internal dynamics + static pins. No collision, no coupling.

#define SB_COMPLIANCE_FLOOR 1e-8   // tiny regularization for rigid links to keep K SPD
#define SB_BAUMGARTE_GAIN   0.2f   // beta for rigid links' position drift correction
#define SB_POS_CLAMP        100.0f // sanity clamp on lambda solution

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
	afree(sb->K);
	afree(sb->D);
	afree(sb->rhs);
	afree(sb->lambda_sol);
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
		w->soft_body_gen[idx]++; // bump to odd (alive)
		memset(&w->soft_bodies[idx], 0, sizeof(SoftBodyInternal));
	} else {
		idx = asize(w->soft_bodies);
		apush(w->soft_bodies, (SoftBodyInternal){0});
		apush(w->soft_body_gen, 1); // start alive
	}
	SoftBodyInternal* sb = &w->soft_bodies[idx];
	sb->params = params;
	if (sb->params.linear_damping <= 0.0f) sb->params.linear_damping = 0.02f;
	return (SoftBody){ handle_make(idx, w->soft_body_gen[idx]) };
}

void destroy_soft_body(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup_safe(w, sb_h);
	if (!sb) return;
	int idx = handle_index(sb_h);
	sb_free_storage(sb);
	w->soft_body_gen[idx]++; // bump to even (dead)
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
	SoftLink lk = (SoftLink){
		.node_i = node_i, .node_j = node_j,
		.rest_length = rest_length,
		.compliance = 0.0f, .pos_to_vel = 0.0f, .lambda = 0.0f,
		.axis = V3(0, 1, 0),
	};
	// Stash SpringParams into the link via the compliance/pos_to_vel pair; we
	// recompute per-substep from the user's SpringParams. But we don't store
	// SpringParams directly on the hot struct -- instead we derive compliance
	// once at build. For v1, let users override via default_spring on params
	// and a later call; simplest path: stash frequency + damping into the link
	// as extra fields OR just treat spring.frequency == 0 as rigid.
	//
	// Implementation choice: keep SpringParams out of the hot struct, record
	// compliance only. Compliance for soft links is `softness` from
	// spring_compute; for rigid links it stays at SB_COMPLIANCE_FLOOR (set at
	// build). Per-link spring is stored as an auxiliary array, not here.
	// For v1 simplicity we go with a single global spring from params.default_spring.
	(void)spring;
	apush(sb->links, lk);
}

void soft_body_build(World world, SoftBody sb_h)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(!sb->built);
	int L = asize(sb->links);
	sb->link_count = L;
	if (L > 0) {
		afit(sb->K, L * L);          asetlen(sb->K, L * L);
		afit(sb->D, L);              asetlen(sb->D, L);
		afit(sb->rhs, L);            asetlen(sb->rhs, L);
		afit(sb->lambda_sol, L);     asetlen(sb->lambda_sol, L);
	}
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
	// Snap effective inv_mass to 0 (and remember user's original for unpin).
	sb->node_inv_mass[node] = 0.0f;
	sb->node_pos[node] = world_pos;
	sb->node_vel[node] = V3(0, 0, 0);
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
// Per-substep solver.
//
// Dense LDL: factor K = L D L^T in place inside sb->K (lower triangle holds L,
// unit diagonal implicit), D separate. Then forward/diagonal/back substitute
// against sb->rhs into sb->lambda_sol.

// Returns 0 on success, -1 if a non-positive pivot is encountered (singular).
static int sb_ldl_factor(double* K, double* D, int n)
{
	for (int k = 0; k < n; k++) {
		double d = K[k * n + k];
		for (int j = 0; j < k; j++) d -= K[k * n + j] * K[k * n + j] * D[j];
		if (d <= 1e-14) return -1;
		D[k] = d;
		for (int i = k + 1; i < n; i++) {
			double s = K[i * n + k];
			for (int j = 0; j < k; j++) s -= K[i * n + j] * K[k * n + j] * D[j];
			K[i * n + k] = s / d;
		}
	}
	return 0;
}

static void sb_ldl_solve(const double* L, const double* D, const double* rhs, double* x, int n)
{
	// L y = rhs (L unit lower)
	for (int k = 0; k < n; k++) {
		double v = rhs[k];
		for (int j = 0; j < k; j++) v -= L[k * n + j] * x[j];
		x[k] = v;
	}
	// D z = y
	for (int k = 0; k < n; k++) x[k] /= D[k];
	// L^T xf = z (back substitute)
	for (int k = n - 1; k >= 0; k--) {
		double v = x[k];
		for (int i = k + 1; i < n; i++) v -= L[i * n + k] * x[i];
		x[k] = v;
	}
}

// PBD-style node-vs-static-rigid collision. v1: static bodies only, SHAPE_BOX
// and SHAPE_SPHERE. Each soft-body node is treated as a sphere of
// sb->params.node_radius. Corrections are applied inline -- no PGS, no
// manifold cache, no ContactSummary. Dynamic rigid bodies are ignored since
// proper coupling needs a K-matrix border block (phase 2+).
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
			// Manifold normal points from A (node sphere) toward B (world).
			// Push node along -normal by penetration; cancel inward velocity;
			// apply simple Coulomb-clamped tangential damping.
			for (int ci = 0; ci < m.count; ci++) {
				Contact* c = &m.contacts[ci];
				sb->node_pos[node_idx] = sub(sb->node_pos[node_idx], scale(c->normal, c->penetration));
				node_sphere.center = sb->node_pos[node_idx];
				v3 v = sb->node_vel[node_idx];
				float vn = dot(v, c->normal);
				if (vn >= 0.0f) continue; // already separating
				v3 vn_vec = scale(c->normal, vn);
				v3 vt = sub(v, vn_vec);
				v3 new_v = vt; // inelastic: zero out normal
				// Friction: clamp |vt| by mu * |vn|
				float mu = 0.4f;
				float vt_len = len(vt);
				float max_friction = mu * (-vn);
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

// Single soft-body substep. Called once per rigid substep from world_step.
static void soft_body_substep(WorldInternal* w, SoftBodyInternal* sb, v3 gravity, float dt)
{
	if (!sb->built) return;
	int N = asize(sb->node_pos);
	int L = sb->link_count;

	// Shared spring params from soft-body defaults.
	float ptv = 0.0f, soft_compliance = 0.0f;
	{
		SpringParams sp = sb->params.default_spring;
		if (sp.frequency > 0.0f) {
			float omega = 6.2831853f * sp.frequency;
			float zeta = sp.damping_ratio > 0.0f ? sp.damping_ratio : 1.0f;
			float hk = dt * omega * omega;
			float d = 2.0f * zeta * omega;
			float denom = dt * (d + hk);
			if (denom > 1e-12f) {
				soft_compliance = 1.0f / denom;
				ptv = hk * soft_compliance;
			}
		}
	}
	int rigid = (soft_compliance == 0.0f); // user asked for rigid default
	float rigid_compliance = SB_COMPLIANCE_FLOOR;
	float rigid_ptv = SB_BAUMGARTE_GAIN / dt;

	// 1) Integrate velocity: gravity + ext_force + damping.
	float damp_mul = 1.0f / (1.0f + sb->params.linear_damping * dt);
	for (int n = 0; n < N; n++) {
		if (sb->node_inv_mass[n] == 0.0f) { sb->node_ext_force[n] = V3(0, 0, 0); continue; }
		v3 a = add(gravity, scale(sb->node_ext_force[n], sb->node_inv_mass[n]));
		sb->node_vel[n] = add(sb->node_vel[n], scale(a, dt));
		sb->node_vel[n] = scale(sb->node_vel[n], damp_mul);
		sb->node_ext_force[n] = V3(0, 0, 0);
	}

	if (L == 0) {
		// No constraints -- just integrate positions and snap pins.
		for (int n = 0; n < N; n++) {
			if (sb->node_inv_mass[n] == 0.0f) continue;
			sb->node_pos[n] = add(sb->node_pos[n], scale(sb->node_vel[n], dt));
		}
		for (int p = 0; p < asize(sb->pins); p++) {
			sb->node_pos[sb->pins[p].node] = sb->pins[p].world_pos;
			sb->node_vel[sb->pins[p].node] = V3(0, 0, 0);
		}
		return;
	}

	// 2) Refresh link axes and rhs.
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		v3 d = sub(sb->node_pos[lk->node_j], sb->node_pos[lk->node_i]);
		float l2 = len2(d);
		float length = l2 > 1e-16f ? sqrtf(l2) : 1e-8f;
		lk->axis = scale(d, 1.0f / length);
		float err = length - lk->rest_length;

		v3 va = sb->node_vel[lk->node_i];
		v3 vb = sb->node_vel[lk->node_j];
		float jv = dot(sub(vb, va), lk->axis);

		float compliance = rigid ? rigid_compliance : soft_compliance;
		float bias_gain  = rigid ? rigid_ptv : ptv;
		lk->compliance = compliance;
		lk->pos_to_vel = bias_gain;

		double bias = (double)bias_gain * (double)err;
		sb->rhs[k] = -(double)jv - bias - (double)compliance * (double)lk->lambda;
	}

	// 3) Build K (dense). K is symmetric, but we fill the full matrix because
	//    sb_ldl_factor reads both halves during elimination.
	double* K = sb->K;
	memset(K, 0, (size_t)L * (size_t)L * sizeof(double));
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		double mi = (double)sb->node_inv_mass[lk->node_i];
		double mj = (double)sb->node_inv_mass[lk->node_j];
		K[k * L + k] = mi + mj + (double)lk->compliance;
	}
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		for (int l = 0; l < k; l++) {
			SoftLink* ll = &sb->links[l];
			double dot_ax = (double)dot(lk->axis, ll->axis);
			double sum = 0.0;
			// sign of a node in a link: -1 for node_i, +1 for node_j.
			// contribution: sign_k(n) * sign_l(n) * inv_m[n] * (a_k . a_l)
			if (lk->node_i == ll->node_i) sum += (double)sb->node_inv_mass[lk->node_i] * dot_ax; // (-)(-)=+
			if (lk->node_i == ll->node_j) sum -= (double)sb->node_inv_mass[lk->node_i] * dot_ax; // (-)(+)=-
			if (lk->node_j == ll->node_i) sum -= (double)sb->node_inv_mass[lk->node_j] * dot_ax; // (+)(-)=-
			if (lk->node_j == ll->node_j) sum += (double)sb->node_inv_mass[lk->node_j] * dot_ax; // (+)(+)=+
			K[k * L + l] = sum;
			K[l * L + k] = sum;
		}
	}

	// 4) Factor + solve.
	if (sb_ldl_factor(K, sb->D, L) != 0) {
		// Singular: bail on this substep. Impulses and positions still advance
		// from step 1's integration-up-to-this-point, so the worst case is a
		// frame of drift, not a divergence.
		goto integrate_and_snap;
	}
	sb_ldl_solve(K, sb->D, sb->rhs, sb->lambda_sol, L);

	// 5) Apply impulses: v += M^-1 J^T lambda.
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		double lam = sb->lambda_sol[k];
		// Clamp paranoia: if factorization gave something preposterous, ignore.
		if (!(lam > -1e12 && lam < 1e12)) continue;
		v3 imp = scale(lk->axis, (float)lam);
		float mi = sb->node_inv_mass[lk->node_i];
		float mj = sb->node_inv_mass[lk->node_j];
		sb->node_vel[lk->node_i] = sub(sb->node_vel[lk->node_i], scale(imp, mi));
		sb->node_vel[lk->node_j] = add(sb->node_vel[lk->node_j], scale(imp, mj));
		// Warm cache: soft links accumulate for spring inertia; rigid SET (since
		// rigid_compliance is tiny, the contribution to next-frame RHS is negligible).
		if (rigid) lk->lambda = (float)lam;
		else       lk->lambda += (float)lam;
	}

integrate_and_snap:
	// 6) Integrate positions.
	for (int n = 0; n < N; n++) {
		if (sb->node_inv_mass[n] == 0.0f) continue;
		sb->node_pos[n] = add(sb->node_pos[n], scale(sb->node_vel[n], dt));
	}

	// 7) Collide each node against world static rigid bodies.
	for (int n = 0; n < N; n++) {
		sb_collide_node_with_world(w, sb, n);
	}

	// 8) Snap pinned nodes to their target world position. Lets the user move
	//    a pin by re-calling soft_body_pin_static -- the node teleports to the
	//    new target each substep.
	for (int p = 0; p < asize(sb->pins); p++) {
		sb->node_pos[sb->pins[p].node] = sb->pins[p].world_pos;
		sb->node_vel[sb->pins[p].node] = V3(0, 0, 0);
	}
}

// World-level entry: advance every live soft body by dt using w->sub_steps
// internal substeps. Decoupled from the rigid substep loop -- soft bodies
// don't touch rigid state in v1, and the world-step entry keeps them stepping
// even on rigid fast-paths (all-asleep shortcut).
static void soft_body_step_world(WorldInternal* w, float dt)
{
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

// Free storage for every live soft body. Called from destroy_world.
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
