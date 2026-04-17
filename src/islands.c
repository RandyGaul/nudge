// islands.c -- island connectivity, sleep/wake management

static int island_create(WorldInternal* w)
{
	int idx;
	if (asize(w->island_free) > 0) {
		idx = apop(w->island_free);
		w->island_gen[idx]++;
	} else {
		idx = asize(w->islands);
		Island zero = {0};
		apush(w->islands, zero);
		apush(w->island_gen, 0);
	}
	w->islands[idx] = (Island){
		.head_body = -1, .tail_body = -1, .body_count = 0,
		.head_joint = -1, .tail_joint = -1, .joint_count = 0,
		.constraint_remove_count = 0, .awake = 1,
	};
	w->island_gen[idx] |= 1; // odd = alive
	return idx;
}

static void island_destroy(WorldInternal* w, int id)
{
	w->island_gen[id]++;
	if (w->island_gen[id] & 1) w->island_gen[id]++; // ensure even = dead
	apush(w->island_free, id);
}

static int island_alive(WorldInternal* w, int id)
{
	return id >= 0 && id < asize(w->islands) && (w->island_gen[id] & 1);
}

static void island_add_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	c->island_id = island_id;
	c->island_prev = isl->tail_body;
	c->island_next = -1;
	if (isl->tail_body >= 0)
		w->body_cold[isl->tail_body].island_next = body_idx;
	else
		isl->head_body = body_idx;
	isl->tail_body = body_idx;
	isl->body_count++;
}

static void island_remove_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	if (c->island_prev >= 0)
		w->body_cold[c->island_prev].island_next = c->island_next;
	else
		isl->head_body = c->island_next;
	if (c->island_next >= 0)
		w->body_cold[c->island_next].island_prev = c->island_prev;
	else
		isl->tail_body = c->island_prev;
	c->island_id = -1;
	c->island_prev = -1;
	c->island_next = -1;
	isl->body_count--;
}

static void island_add_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	j->island_id = island_id;
	j->island_prev = isl->tail_joint;
	j->island_next = -1;
	if (isl->tail_joint >= 0)
		w->joints[isl->tail_joint].island_next = joint_idx;
	else
		isl->head_joint = joint_idx;
	isl->tail_joint = joint_idx;
	isl->joint_count++;
}

static void island_remove_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_prev >= 0)
		w->joints[j->island_prev].island_next = j->island_next;
	else
		isl->head_joint = j->island_next;
	if (j->island_next >= 0)
		w->joints[j->island_next].island_prev = j->island_prev;
	else
		isl->tail_joint = j->island_prev;
	j->island_id = -1;
	j->island_prev = -1;
	j->island_next = -1;
	isl->joint_count--;
}

// Merge two islands. Keep the bigger one, walk the smaller one remapping bodies/joints.
static int island_merge(WorldInternal* w, int id_a, int id_b)
{
	if (id_a == id_b) return id_a;
	Island* a = &w->islands[id_a];
	Island* b = &w->islands[id_b];
	// Keep the bigger island
	if (a->body_count + a->joint_count < b->body_count + b->joint_count) {
		int tmp = id_a; id_a = id_b; id_b = tmp;
		a = &w->islands[id_a]; b = &w->islands[id_b];
	}
	// Remap smaller's bodies to bigger
	int bi = b->head_body;
	while (bi >= 0) {
		int next = w->body_cold[bi].island_next;
		w->body_cold[bi].island_id = id_a;
		bi = next;
	}
	// Splice body lists: append b's list to a's tail
	if (b->head_body >= 0) {
		if (a->tail_body >= 0) {
			w->body_cold[a->tail_body].island_next = b->head_body;
			w->body_cold[b->head_body].island_prev = a->tail_body;
		} else {
			a->head_body = b->head_body;
		}
		a->tail_body = b->tail_body;
		a->body_count += b->body_count;
	}
	// Remap smaller's joints
	int ji = b->head_joint;
	while (ji >= 0) {
		int next = w->joints[ji].island_next;
		w->joints[ji].island_id = id_a;
		ji = next;
	}
	// Splice joint lists
	if (b->head_joint >= 0) {
		if (a->tail_joint >= 0) {
			w->joints[a->tail_joint].island_next = b->head_joint;
			w->joints[b->head_joint].island_prev = a->tail_joint;
		} else {
			a->head_joint = b->head_joint;
		}
		a->tail_joint = b->tail_joint;
		a->joint_count += b->joint_count;
	}
	a->constraint_remove_count += b->constraint_remove_count;
	// If either was awake, result is awake
	if (b->awake) a->awake = 1;
	island_destroy(w, id_b);
	return id_a;
}

// Transfer a body's BVH leaf from one tree to another.
static void bvh_transfer_body(WorldInternal* w, int body_idx, BVH_Tree* from, BVH_Tree* to)
{
	int leaf = w->body_cold[body_idx].bvh_leaf;
	if (leaf < 0) return;
	assert(leaf < asize(from->leaves) && "bvh_transfer: leaf index out of range");
	assert(from->leaves[leaf].body_idx == body_idx && "bvh_transfer: leaf doesn't belong to this body");
	int moved = bvh_remove(from, leaf);
	if (moved >= 0) w->body_cold[moved].bvh_leaf = leaf;
	AABB box = aabb_expand(body_aabb(&w->body_state[body_idx], &w->body_cold[body_idx]), BVH_AABB_MARGIN);
	w->body_cold[body_idx].bvh_leaf = bvh_insert(to, body_idx, box);
}

static void island_wake(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 1;
	int bi = isl->head_body;
	while (bi >= 0) {
		body_sleep_time(w, bi) = 0.0f;
		bi = w->body_cold[bi].island_next;
	}
}

static void island_sleep(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 0;
	ldl_cache_sleep(&isl->ldl);
	int bi = isl->head_body;
	while (bi >= 0) {
		body_vel(w, bi) = V3(0, 0, 0);
		body_angvel(w, bi) = V3(0, 0, 0);
		bi = w->body_cold[bi].island_next;
	}
}

// Link a newly created joint to islands: merge/create as needed, wake if sleeping.
static void link_joint_to_islands(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	int ba = j->body_a, bb = j->body_b;
	int is_static_a = body_inv_mass(w, ba) == 0.0f;
	int is_static_b = body_inv_mass(w, bb) == 0.0f;
	// Static bodies never get islands
	if (is_static_a && is_static_b) return;
	int isl_a = is_static_a ? -1 : w->body_cold[ba].island_id;
	int isl_b = is_static_b ? -1 : w->body_cold[bb].island_id;
	int target;
	if (isl_a >= 0 && isl_b >= 0) {
		target = island_merge(w, isl_a, isl_b);
	} else if (isl_a >= 0) {
		target = isl_a;
	} else if (isl_b >= 0) {
		target = isl_b;
	} else {
		target = island_create(w);
	}
	// Add dynamic bodies that aren't in the target island yet
	if (!is_static_a && w->body_cold[ba].island_id != target)
		island_add_body(w, target, ba);
	if (!is_static_b && w->body_cold[bb].island_id != target)
		island_add_body(w, target, bb);
	island_add_joint(w, target, joint_idx);
	// Wake if sleeping
	if (!w->islands[target].awake)
		island_wake(w, target);
}

// Unlink a joint from its island before destruction.
static void unlink_joint_from_island(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_id < 0) return;
	int isl = j->island_id;
	island_remove_joint(w, isl, joint_idx);
	w->islands[isl].constraint_remove_count++;
}

// Update islands from this frame's contact manifolds.
// Merge/create islands for touching dynamic body pairs.
// Detect lost contacts for island splitting.
static void islands_update_contacts(WorldInternal* w, InternalManifold* manifolds, int manifold_count)
{
	CK_MAP(uint8_t) curr_touching = NULL;

	for (int i = 0; i < manifold_count; i++) {
		int a = manifolds[i].body_a, b = manifolds[i].body_b;
		int is_static_a = body_inv_mass(w, a) == 0.0f;
		int is_static_b = body_inv_mass(w, b) == 0.0f;
		// Skip static-static (shouldn't happen, but guard)
		if (is_static_a && is_static_b) continue;

		uint64_t key = body_pair_key(a, b);
		map_set(curr_touching, key, (uint8_t)1);

		// Static bodies never get island_id
		int isl_a = is_static_a ? -1 : w->body_cold[a].island_id;
		int isl_b = is_static_b ? -1 : w->body_cold[b].island_id;
		// Capture awake state before merge may destroy the source island
		int was_awake_a = isl_a >= 0 && w->islands[isl_a].awake;
		int was_awake_b = isl_b >= 0 && w->islands[isl_b].awake;

		int target;
		if (isl_a >= 0 && isl_b >= 0) {
			target = island_merge(w, isl_a, isl_b);
		} else if (isl_a >= 0) {
			target = isl_a;
		} else if (isl_b >= 0) {
			target = isl_b;
		} else {
			target = island_create(w);
		}
		int added_a = 0, added_b = 0;
		if (!is_static_a && w->body_cold[a].island_id != target) {
			island_add_body(w, target, a);
			added_a = 1;
		}
		if (!is_static_b && w->body_cold[b].island_id != target) {
			island_add_body(w, target, b);
			added_b = 1;
		}
		// Only wake a sleeping island if an AWAKE body or NEW body joins.
		// Sleeping bodies re-merging after a split should not wake the island.
		if (!w->islands[target].awake) {
			int should_wake = 0;
			if (added_a && !is_static_a && (isl_a < 0 || was_awake_a)) should_wake = 1;
			if (added_b && !is_static_b && (isl_b < 0 || was_awake_b)) should_wake = 1;
			if (!added_a && !added_b) {
				if (was_awake_a && isl_a != target) should_wake = 1;
				if (was_awake_b && isl_b != target) should_wake = 1;
			}
			if (should_wake) island_wake(w, target);
		}
	}

	// Detect lost contacts: pairs in prev_touching but not in curr_touching.
	// Increment constraint_remove_count on their island so island_try_split fires.
	map_each(w->prev_touching, pi) {
		uint64_t pk = map_key(w->prev_touching, pi);
		if (map_get_ptr(curr_touching, pk)) continue; // still touching
		// Contact lost. Find which island this pair belongs to.
		int pa = (int)(pk >> 32), pb = (int)(pk & 0xFFFFFFFF);
		if (pa >= asize(w->body_hot) || pb >= asize(w->body_hot)) continue;
		if (!split_alive(w->body_gen, pa) || !split_alive(w->body_gen, pb)) continue;
		int isl = w->body_cold[pa].island_id;
		if (isl < 0) isl = w->body_cold[pb].island_id;
		if (isl >= 0 && (w->island_gen[isl] & 1)) w->islands[isl].constraint_remove_count++;
	}

	// Swap: prev_touching = curr_touching
	map_free(w->prev_touching);
	w->prev_touching = curr_touching;
}

// Return the island that owns a body pair. Under normal operation (sleep
// enabled), both dynamic bodies share an island after islands_update_contacts /
// link_joint_to_island; one body may be static (island_id == -1), in which
// case the other body's island owns the pair. When sleep is disabled, contact
// pairs don't link islands, so two dynamic bodies can briefly be in different
// islands -- pick the lower id so bucket placement stays deterministic.
// Both-static pairs return -1.
static int pair_owning_island(WorldInternal* w, int a, int b)
{
	int isl_a = body_inv_mass(w, a) == 0.0f ? -1 : w->body_cold[a].island_id;
	int isl_b = body_inv_mass(w, b) == 0.0f ? -1 : w->body_cold[b].island_id;
	if (isl_a >= 0 && isl_b >= 0) return isl_a < isl_b ? isl_a : isl_b;
	if (isl_a >= 0) return isl_a;
	if (isl_b >= 0) return isl_b;
	return -1;
}

static int manifold_owning_island(WorldInternal* w, InternalManifold* m)
{
	return pair_owning_island(w, m->body_a, m->body_b);
}

// Partition manifolds by owning island. Arrays are arena-allocated; callers
// must not afree them.
//
//   *out_offsets -- int[n_islands + 1]. Exclusive prefix sum: island i's
//                   manifolds live in perm[offsets[i] .. offsets[i+1]).
//   *out_perm    -- int[offsets[n_islands]]. Permutation of manifold indices;
//                   manifolds with two static bodies are omitted.
//
// Determinism: within each island, manifold indices appear in ascending
// original order (stable partition over BVH pair-emission order).
//
// Returns n_islands (== asize(w->islands)).
static int islands_bucket_manifolds(WorldInternal* w, InternalManifold* manifolds, int manifold_count, Arena* arena, int** out_offsets, int** out_perm)
{
	int n_islands = asize(w->islands);

	// No islands to bucket into. Return NULLs; callers check island_manifold_perm
	// before using it. Without this guard the zero-length stretchy-buffer path
	// below hits asetlen(NULL, 0), which dereferences a NULL header.
	if (n_islands == 0) {
		*out_offsets = NULL;
		*out_perm = NULL;
		return 0;
	}

	int* offsets = NULL;
	arena_fit(arena, offsets, n_islands + 1);
	asetlen(offsets, n_islands + 1);
	for (int i = 0; i <= n_islands; i++) offsets[i] = 0;

	// Count pass. Store count in offsets[i+1] so the prefix sum below produces
	// exclusive-prefix offsets[] in one forward scan.
	for (int m = 0; m < manifold_count; m++) {
		int isl = manifold_owning_island(w, &manifolds[m]);
		if (isl < 0) continue;
		offsets[isl + 1]++;
	}
	for (int i = 1; i <= n_islands; i++) offsets[i] += offsets[i - 1];
	int total = offsets[n_islands];

	int* perm = NULL;
	if (total > 0) {
		arena_fit(arena, perm, total);
		asetlen(perm, total);
	}

	// Running cursor per island, initialized from offsets[]. A copy avoids
	// clobbering the exclusive-prefix offsets we return to the caller.
	int* cursor = NULL;
	arena_fit(arena, cursor, n_islands);
	asetlen(cursor, n_islands);
	for (int i = 0; i < n_islands; i++) cursor[i] = offsets[i];

	for (int m = 0; m < manifold_count; m++) {
		int isl = manifold_owning_island(w, &manifolds[m]);
		if (isl < 0) continue;
		perm[cursor[isl]++] = m;
	}

	*out_offsets = offsets;
	*out_perm = perm;
	return n_islands;
}

// Same shape as islands_bucket_manifolds, but over sol_joints (the solver's
// per-frame joint array populated by joints_pre_solve). Active joints always
// have at least one dynamic body, which owns the joint's island.
static int islands_bucket_sol_joints(WorldInternal* w, SolverJoint* sol_joints, int joint_count, Arena* arena, int** out_offsets, int** out_perm)
{
	int n_islands = asize(w->islands);

	if (n_islands == 0) {
		*out_offsets = NULL;
		*out_perm = NULL;
		return 0;
	}

	int* offsets = NULL;
	arena_fit(arena, offsets, n_islands + 1);
	asetlen(offsets, n_islands + 1);
	for (int i = 0; i <= n_islands; i++) offsets[i] = 0;

	for (int j = 0; j < joint_count; j++) {
		int isl = pair_owning_island(w, sol_joints[j].body_a, sol_joints[j].body_b);
		if (isl < 0) continue;
		offsets[isl + 1]++;
	}
	for (int i = 1; i <= n_islands; i++) offsets[i] += offsets[i - 1];
	int total = offsets[n_islands];

	int* perm = NULL;
	if (total > 0) {
		arena_fit(arena, perm, total);
		asetlen(perm, total);
	}

	int* cursor = NULL;
	arena_fit(arena, cursor, n_islands);
	asetlen(cursor, n_islands);
	for (int i = 0; i < n_islands; i++) cursor[i] = offsets[i];

	for (int j = 0; j < joint_count; j++) {
		int isl = pair_owning_island(w, sol_joints[j].body_a, sol_joints[j].body_b);
		if (isl < 0) continue;
		perm[cursor[isl]++] = j;
	}

	*out_offsets = offsets;
	*out_perm = perm;
	return n_islands;
}

static void island_try_split(WorldInternal* w, int island_id); // forward decl

// Split islands that lost contacts. Must run every frame regardless of sleep.
static void islands_try_splits(WorldInternal* w)
{
	int island_count = asize(w->islands);
	for (int i = 0; i < island_count; i++) {
		if (!(w->island_gen[i] & 1)) continue;
		if (!w->islands[i].awake) continue;
		if (w->islands[i].constraint_remove_count > 0)
			island_try_split(w, i);
		// Re-index after split (island_create may realloc w->islands)
		w->islands[i].constraint_remove_count = 0;
	}
}

// Evaluate sleep for all awake islands after post-solve.
static void islands_evaluate_sleep(WorldInternal* w, float dt)
{
	int island_count = asize(w->islands);
	for (int i = 0; i < island_count; i++) {
		if (!(w->island_gen[i] & 1)) continue; // dead slot
		Island* isl = &w->islands[i];
		if (!isl->awake) continue;

		int all_sleepy = 1;
		int bi = isl->head_body;
		while (bi >= 0) {
			BodyHot* h = &w->body_hot[bi];
			BodyState* s = &w->body_state[bi];
			if (!s->sleep_allowed) { all_sleepy = 0; bi = w->body_cold[bi].island_next; continue; }
			float v2 = len2(h->velocity) + len2(h->angular_velocity);
			if (v2 > SLEEP_VEL_THRESHOLD) {
				s->sleep_time = fmaxf(s->sleep_time - dt * 2.0f, 0.0f);
				all_sleepy = 0;
			} else {
				s->sleep_time += dt;
				if (s->sleep_time < SLEEP_TIME_THRESHOLD)
					all_sleepy = 0;
			}
			bi = w->body_cold[bi].island_next;
		}

		if (all_sleepy)
			island_sleep(w, i);
	}
}

// DFS island split when contacts have been lost.
static void island_try_split(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];

	// Collect all body indices from island
	CK_DYNA int* bodies = NULL;
	int bi = isl->head_body;
	while (bi >= 0) {
		apush(bodies, bi);
		bi = w->body_cold[bi].island_next;
	}
	int n = asize(bodies);
	if (n <= 1) { afree(bodies); return; } // single body, nothing to split

	// Visited array and component assignment
	CK_DYNA int* component = NULL;
	afit(component, n);
	asetlen(component, n);
	for (int i = 0; i < n; i++) component[i] = -1;

	// Map body_idx -> local index for fast lookup
	CK_MAP(int) body_to_local = NULL;
	for (int i = 0; i < n; i++) {
		map_set(body_to_local, (uint64_t)bodies[i], i);
	}

	int num_components = 0;
	CK_DYNA int* stack = NULL;

	for (int start = 0; start < n; start++) {
		if (component[start] >= 0) continue;
		int comp_id = num_components++;

		// DFS from start
		aclear(stack);
		apush(stack, start);
		component[start] = comp_id;

		while (asize(stack) > 0) {
			int cur_local = apop(stack);
			int cur_body = bodies[cur_local];

			// Neighbors from joints in this island
			int ji = isl->head_joint;
			while (ji >= 0) {
				JointInternal* j = &w->joints[ji];
				int other = -1;
				if (j->body_a == cur_body) other = j->body_b;
				else if (j->body_b == cur_body) other = j->body_a;
				if (other >= 0) {
					int* loc = map_get_ptr(body_to_local, (uint64_t)other);
					if (loc && component[*loc] < 0) {
						component[*loc] = comp_id;
						apush(stack, *loc);
					}
				}
				ji = w->joints[ji].island_next;
			}

			// Neighbors from current contacts (prev_touching map)
			for (int j = 0; j < n; j++) {
				if (component[j] >= 0) continue;
				uint64_t key = body_pair_key(cur_body, bodies[j]);
				uint8_t* val = map_get_ptr(w->prev_touching, key);
				if (val) {
					component[j] = comp_id;
					apush(stack, j);
				}
			}
		}
	}

	if (num_components <= 1) {
		// Still one connected component, no split needed
	} else {
		// First component keeps the original island. Remove all bodies/joints,
		// then re-add per component.
		// Detach all bodies from the original island
		for (int i = 0; i < n; i++) {
			w->body_cold[bodies[i]].island_id = -1;
			w->body_cold[bodies[i]].island_prev = -1;
			w->body_cold[bodies[i]].island_next = -1;
		}
		// Detach all joints
		CK_DYNA int* joint_list = NULL;
		int ji = isl->head_joint;
		while (ji >= 0) {
			apush(joint_list, ji);
			int next = w->joints[ji].island_next;
			w->joints[ji].island_id = -1;
			w->joints[ji].island_prev = -1;
			w->joints[ji].island_next = -1;
			ji = next;
		}
		// Reset the original island
		*isl = (Island){ .head_body = -1, .tail_body = -1, .head_joint = -1, .tail_joint = -1, .awake = 1 };

		// Create islands per component (reuse original for component 0)
		CK_DYNA int* comp_islands = NULL;
		afit(comp_islands, num_components);
		asetlen(comp_islands, num_components);
		comp_islands[0] = island_id;
		for (int c = 1; c < num_components; c++) {
			comp_islands[c] = island_create(w);
		}

		// Assign bodies to their component's island
		for (int i = 0; i < n; i++) {
			island_add_body(w, comp_islands[component[i]], bodies[i]);
		}

		// Assign joints to the island of body_a (both endpoints should be in same component)
		for (int i = 0; i < asize(joint_list); i++) {
			int jidx = joint_list[i];
			int ba = w->joints[jidx].body_a;
			int* loc = map_get_ptr(body_to_local, (uint64_t)ba);
			int comp = loc ? component[*loc] : 0;
			island_add_joint(w, comp_islands[comp], jidx);
		}

		afree(comp_islands);
		afree(joint_list);
	}

	afree(stack);
	map_free(body_to_local);
	afree(component);
	afree(bodies);
}
