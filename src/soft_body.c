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

// Diagnostics: updated each soft_body_step_world, read by test harness/debug UI.
double g_soft_body_max_lambda;
double g_soft_body_max_rhs;
double g_soft_body_min_K_diag;
int    g_soft_body_trace;

#define SB_COMPLIANCE_FLOOR      5e-5  // regularization floor for rigid links, scaled by trace(K)/n
#define SB_OVER_CONSTRAINT_SCALE 1.0   // compliance *= (1 + scale * redundant_dofs) for over-constrained islands
#define SB_BAUMGARTE_GAIN        0.2f  // beta for rigid links' position drift correction
#define SB_MAX_PUSH_VEL          3.0f  // per-substep clamp on Baumgarte bias magnitude (m/s)

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
	sb->k_dirty = 1; // force K build + factor on first step
}

void soft_body_pin_static(World world, SoftBody sb_h, int node, v3 world_pos)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SoftBodyInternal* sb = sb_lookup(w, sb_h);
	assert(node >= 0 && node < asize(sb->node_pos));
	int pi = sb_find_pin(sb, node);
	int was_pinned = pi >= 0;
	if (pi < 0) {
		SoftPin p = (SoftPin){ .node = node, .world_pos = world_pos };
		apush(sb->pins, p);
	} else {
		sb->pins[pi].world_pos = world_pos;
	}
	sb->node_pos[node] = world_pos;
	sb->node_vel[node] = V3(0, 0, 0);
	// inv_mass transitions to 0 only on first pin for this node; K diagonal
	// changes then, so mark dirty. Re-pinning the same node (moving the pin
	// target) doesn't change K -- pos/vel snap each substep is enough.
	if (!was_pinned) {
		sb->node_inv_mass[node] = 0.0f;
		sb->k_dirty = 1;
	}
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
	sb->k_dirty = 1; // inv_mass changed -> K diagonal changes
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

// Modified-Cholesky factorization: LDL with pivot floor. Pivots below
// `min_pivot` are raised to that floor; this prevents 1/D amplification of
// null-space RHS components in over-constrained / rank-deficient systems.
// Returns 0 on success, -1 only if a pivot is non-finite.
static int sb_ldl_factor(double* K, double* D, int n, double min_pivot)
{
	for (int k = 0; k < n; k++) {
		double d = K[k * n + k];
		for (int j = 0; j < k; j++) d -= K[k * n + j] * K[k * n + j] * D[j];
		if (!(d == d)) return -1; // NaN guard
		if (d < min_pivot) d = min_pivot;
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
			// Manifold normal points from A (node sphere) toward B (world body).
			// For a node penetrating a static body, `normal` points into the body,
			// so:
			//   - dot(v, normal) > 0 => node moving INTO the surface (cancel)
			//   - dot(v, normal) < 0 => node moving away (already separating)
			// Position correction: push the node back along -normal by penetration.
			for (int ci = 0; ci < m.count; ci++) {
				Contact* c = &m.contacts[ci];
				sb->node_pos[node_idx] = sub(sb->node_pos[node_idx], scale(c->normal, c->penetration));
				node_sphere.center = sb->node_pos[node_idx];
				v3 v = sb->node_vel[node_idx];
				float vn = dot(v, c->normal);
				if (vn <= 0.0f) continue; // already separating
				// Remove the inward (into-body) component.
				v3 vt = sub(v, scale(c->normal, vn));
				// Coulomb-clamped tangential damping against max_friction = mu * |vn|.
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

// Once-per-frame: refresh link axes, compute per-link compliance/bias, build K
// from current state, factor K. The factored K is used for every substep this
// frame -- factorization is still the dominant cost, so a future optimization
// pass should avoid rebuilding/refactoring when K is approximately unchanged
// (stiff body undergoing only rigid motion; see project_soft_bodies.md).
//
// Returns 0 on success, -1 if K is singular (caller skips solve but still
// integrates + collides).
static int sb_prepare_frame(SoftBodyInternal* sb, float sub_dt)
{
	int L = sb->link_count;
	if (L == 0) return 0;

	// Shared spring params for rigid vs soft default.
	float soft_compliance = 0.0f, soft_ptv = 0.0f;
	SpringParams sp = sb->params.default_spring;
	if (sp.frequency > 0.0f) {
		float omega = 6.2831853f * sp.frequency;
		float zeta = sp.damping_ratio > 0.0f ? sp.damping_ratio : 1.0f;
		float hk = sub_dt * omega * omega;
		float d = 2.0f * zeta * omega;
		float denom = sub_dt * (d + hk);
		if (denom > 1e-12f) {
			soft_compliance = 1.0f / denom;
			soft_ptv = hk * soft_compliance;
		}
	}
	int rigid_default = (soft_compliance == 0.0f);
	float rigid_ptv = SB_BAUMGARTE_GAIN / sub_dt;

	// Refresh axes.
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		v3 d = sub(sb->node_pos[lk->node_j], sb->node_pos[lk->node_i]);
		float l2 = len2(d);
		float length = l2 > 1e-16f ? sqrtf(l2) : 1e-8f;
		lk->axis = scale(d, 1.0f / length);
	}

	// Build K (symmetric, full storage). Diagonal first.
	double* K = sb->K;
	memset(K, 0, (size_t)L * (size_t)L * sizeof(double));
	double trace = 0.0;
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		double mi = (double)sb->node_inv_mass[lk->node_i];
		double mj = (double)sb->node_inv_mass[lk->node_j];
		K[k * L + k] = mi + mj;
		trace += mi + mj;
	}
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		for (int l = 0; l < k; l++) {
			SoftLink* ll = &sb->links[l];
			double dot_ax = (double)dot(lk->axis, ll->axis);
			double sum = 0.0;
			if (lk->node_i == ll->node_i) sum += (double)sb->node_inv_mass[lk->node_i] * dot_ax;
			if (lk->node_i == ll->node_j) sum -= (double)sb->node_inv_mass[lk->node_i] * dot_ax;
			if (lk->node_j == ll->node_i) sum -= (double)sb->node_inv_mass[lk->node_j] * dot_ax;
			if (lk->node_j == ll->node_j) sum += (double)sb->node_inv_mass[lk->node_j] * dot_ax;
			K[k * L + l] = sum;
			K[l * L + k] = sum;
		}
	}

	// Regularization. Count dynamic DOFs vs constraint DOFs; if over-constrained,
	// scale compliance by (1 + scale * redundant_dofs) -- K is rank-deficient in
	// that case and the null space needs meaningful regularization to be SPD.
	// Without this, tiny compliance leaves a near-singular K that produces
	// garbage lambdas (ball explodes, tunnels through floor, etc.).
	int dyn_dof = 0;
	for (int n = 0; n < asize(sb->node_pos); n++) {
		if (sb->node_inv_mass[n] > 0.0f) dyn_dof += 3;
	}
	int redundant = L - dyn_dof; // positive if over-constrained
	double over_scale = redundant > 0 ? (1.0 + SB_OVER_CONSTRAINT_SCALE * (double)redundant) : 1.0;
	double diag_avg = trace > 0.0 ? trace / (double)L : 1.0;
	double rigid_compliance = SB_COMPLIANCE_FLOOR * diag_avg * over_scale;

	// Apply compliance to K diagonal + store per-link params for RHS.
	// K diagonal and RHS warm-start term use the SAME compliance value
	// (standard XPBD: `K*lambda + alpha*lambda_warm = -Jv - alpha/dt * C`).
	for (int k = 0; k < L; k++) {
		SoftLink* lk = &sb->links[k];
		double c = rigid_default ? rigid_compliance : (double)soft_compliance;
		K[k * L + k] += c;
		lk->compliance = (float)c;
		lk->pos_to_vel = rigid_default ? rigid_ptv : soft_ptv;
	}

	// Pivot floor for modified-Cholesky: keeps null-space pivots from
	// collapsing to ~compliance, which would amplify lambdas by 1/compliance.
	// Scaling by diag_avg keeps the floor mass-adaptive.
	double min_pivot = 0.1 * diag_avg;
	return sb_ldl_factor(K, sb->D, L, min_pivot);
}

// Per-substep: integrate velocity, refresh axes from current positions,
// compute RHS + back-substitute against the factored K, apply impulses,
// integrate positions, collide, snap pins. k_factored = result of
// sb_prepare_frame (0 = usable, -1 = skip the solve).
//
// Axes are refreshed from current positions every substep (cheap: one sqrt +
// scale per link) because the RHS projection jv = axis . (vb - va) is not
// rotation-invariant -- using stale axes from frame start produces spurious
// "velocity error" during rigid-body rotation of the soft body and fights the
// rotation. K's off-diagonals ARE invariant under rigid rotation (dot
// products preserved), so factoring once per frame stays valid even though
// axes advance within the frame.
static void soft_body_substep(WorldInternal* w, SoftBodyInternal* sb, v3 gravity, float dt, int k_factored)
{
	if (!sb->built) return;
	int N = asize(sb->node_pos);
	int L = sb->link_count;

	// 1) Integrate velocity: gravity + ext_force + damping.
	float damp_mul = 1.0f / (1.0f + sb->params.linear_damping * dt);
	for (int n = 0; n < N; n++) {
		if (sb->node_inv_mass[n] == 0.0f) { sb->node_ext_force[n] = V3(0, 0, 0); continue; }
		v3 a = add(gravity, scale(sb->node_ext_force[n], sb->node_inv_mass[n]));
		sb->node_vel[n] = add(sb->node_vel[n], scale(a, dt));
		sb->node_vel[n] = scale(sb->node_vel[n], damp_mul);
		sb->node_ext_force[n] = V3(0, 0, 0);
	}

	if (L > 0 && k_factored == 0) {
		// 2) Refresh axes from CURRENT positions (not the frame-start cache).
		//    Then compute RHS with those fresh axes. The Baumgarte position
		//    bias `pos_to_vel * err` is clamped to +/-SB_MAX_PUSH_VEL per
		//    substep so a large sudden position error (mouse drag teleport,
		//    deep collision penetration) doesn't produce an impulse that
		//    explodes the body.
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
			float bias = lk->pos_to_vel * err;
			if (bias >  SB_MAX_PUSH_VEL) bias =  SB_MAX_PUSH_VEL;
			if (bias < -SB_MAX_PUSH_VEL) bias = -SB_MAX_PUSH_VEL;
			sb->rhs[k] = -(double)jv - (double)bias - (double)lk->compliance * (double)lk->lambda;
		}

		// 3) Back-substitute against the factored K.
		sb_ldl_solve(sb->K, sb->D, sb->rhs, sb->lambda_sol, L);
		// Diagnostics.
		double lam_peak = 0, rhs_peak = 0;
		int lam_peak_k = 0;
		for (int k = 0; k < L; k++) {
			double mag = sb->lambda_sol[k]; if (mag < 0) mag = -mag;
			if (mag > lam_peak) { lam_peak = mag; lam_peak_k = k; }
			double rm = sb->rhs[k]; if (rm < 0) rm = -rm;
			if (rm > rhs_peak) rhs_peak = rm;
			if (sb->D[k] < g_soft_body_min_K_diag) g_soft_body_min_K_diag = sb->D[k];
		}
		if (lam_peak > g_soft_body_max_lambda) g_soft_body_max_lambda = lam_peak;
		if (rhs_peak > g_soft_body_max_rhs) g_soft_body_max_rhs = rhs_peak;

		(void)lam_peak_k;

		// 4) Apply impulses along CURRENT axes.
		int rigid_default = (sb->params.default_spring.frequency <= 0.0f);
		for (int k = 0; k < L; k++) {
			SoftLink* lk = &sb->links[k];
			double lam = sb->lambda_sol[k];
			if (!(lam > -1e12 && lam < 1e12)) continue;
			v3 imp = scale(lk->axis, (float)lam);
			float mi = sb->node_inv_mass[lk->node_i];
			float mj = sb->node_inv_mass[lk->node_j];
			sb->node_vel[lk->node_i] = sub(sb->node_vel[lk->node_i], scale(imp, mi));
			sb->node_vel[lk->node_j] = add(sb->node_vel[lk->node_j], scale(imp, mj));
			if (rigid_default) lk->lambda = (float)lam;
			else               lk->lambda += (float)lam;
		}
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

	// 7) Snap pinned nodes to pin target (idempotent; lets the user move a
	//    pin by calling soft_body_pin_static again with a new target).
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
	g_soft_body_max_lambda = 0.0;
	g_soft_body_max_rhs = 0.0;
	g_soft_body_min_K_diag = 1e30;
	int count = asize(w->soft_bodies);
	if (count == 0) return;
	int n_sub = w->sub_steps > 0 ? w->sub_steps : 1;
	float sub_dt = dt / (float)n_sub;

	// TODO(perf): today we rebuild + factor K once per frame. That cost
	// dominates for ~100+ link bodies. Optimization passes worth trying:
	//   1. Skip factor when K is approximately unchanged (dirty flag for
	//      topology/pin change; works well for stiff bodies with small
	//      deformation, breaks for squishy jelly).
	//   2. Sparse LDL reusing the helpers in solver_ldl.c -- K is highly
	//      sparse for soft bodies (each link only adjacent to a few others).
	//   3. Deferred refactor: reuse the factorization for N frames when the
	//      soft body is not deforming much (detectable from lambda magnitudes).
	// Keep correctness simple for now: factor every frame.
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->soft_body_gen, i)) continue;
		SoftBodyInternal* sb = &w->soft_bodies[i];
		if (!sb->built) continue;
		if (sb->link_count == 0) { sb->k_factor_ok = 0; continue; }
		int rc = sb_prepare_frame(sb, sub_dt);
		sb->k_factor_ok = (rc == 0);
		sb->k_sub_dt = sub_dt;
		sb->k_dirty = 0;
	}

	// Substeps re-use the factored K; they refresh axes, compute RHS, back-sub,
	// apply, integrate, collide, snap pins.
	for (int s = 0; s < n_sub; s++) {
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->soft_body_gen, i)) continue;
			SoftBodyInternal* sb = &w->soft_bodies[i];
			soft_body_substep(w, sb, w->gravity, sub_dt, sb->k_factor_ok ? 0 : -1);
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
