// See LICENSE for licensing info.
// broadphase.c -- broad-phase collision detection (SAP, BVH, joint filtering)

static uint64_t body_pair_key(int a, int b)
{
	uint32_t lo = a < b ? a : b;
	uint32_t hi = a < b ? b : a;
	return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static int jointed_pair_skip(CK_MAP(uint8_t) joint_pairs, int a, int b)
{
	if (!joint_pairs) return 0;
	uint64_t key = body_pair_key(a, b);
	return map_get_ptr(joint_pairs, key) != NULL;
}

// Collision filter: two bodies collide iff (a.group & b.mask) && (b.group & a.mask).
// Cheap bitwise check; skips the narrowphase entirely for filtered pairs.
static int filter_pair_skip(WorldInternal* w, int a, int b)
{
	uint32_t ga = w->body_cold[a].collision_group, ma = w->body_cold[a].collision_mask;
	uint32_t gb = w->body_cold[b].collision_group, mb = w->body_cold[b].collision_mask;
	return ((ga & mb) == 0) || ((gb & ma) == 0);
}

static void broadphase_n2(WorldInternal* w, InternalManifold** manifolds)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (asize(w->body_cold[i].shapes) == 0) continue;
		for (int j = i + 1; j < count; j++) {
			if (!split_alive(w->body_gen, j)) continue;
			if (asize(w->body_cold[j].shapes) == 0) continue;
			// Skip pair if both bodies are inactive (static or sleeping)
			int isl_i = w->body_cold[i].island_id, isl_j = w->body_cold[j].island_id;
			int inactive_i = body_inv_mass(w, i) == 0.0f || (isl_i >= 0 && (w->island_gen[isl_i] & 1) && !w->islands[isl_i].awake);
			int inactive_j = body_inv_mass(w, j) == 0.0f || (isl_j >= 0 && (w->island_gen[isl_j] & 1) && !w->islands[isl_j].awake);
			if (inactive_i && inactive_j) continue;
			if (jointed_pair_skip(w->joint_pairs, i, j)) continue;
			if (filter_pair_skip(w, i, j)) continue;
			narrowphase_pair(w, i, j, manifolds);
		}
	}
}

// Sweep-and-prune entry for axis-sorted broadphase.
typedef struct SAP_Entry
{
	int body_idx;
	float min_val, max_val; // along chosen axis
} SAP_Entry;

static int sap_cmp(const void* a, const void* b) { float d = ((SAP_Entry*)a)->min_val - ((SAP_Entry*)b)->min_val; return (d > 0) - (d < 0); }

// Insertion sort: O(n) for nearly-sorted data (settled physics), O(n^2) worst case.
static void sap_insertion_sort(SAP_Entry* arr, int n)
{
	for (int i = 1; i < n; i++) {
		SAP_Entry key = arr[i];
		int j = i - 1;
		while (j >= 0 && arr[j].min_val > key.min_val) { arr[j + 1] = arr[j]; j--; }
		arr[j + 1] = key;
	}
}

// Broadphase sub-phase timing accumulators (seconds, summed across frames).
double bp_refit_acc, bp_precomp_acc, bp_sweep_acc, bp_cross_acc;
int bp_frame_count;

static void broadphase_bvh(WorldInternal* w, InternalManifold** manifolds)
{
	double t0 = perf_now();
	bvh_refit(w->bvh_dynamic, w);
	bp_refit_acc += perf_now() - t0;

	// Build tight AABBs + sweep-and-prune entries for awake dynamic bodies only.
	// Sleeping bodies live in bvh_sleeping and are excluded from SAP.
	double t1 = perf_now();
	int body_count = asize(w->body_hot);
	AABB tight_stack[256];
	AABB* tight = (body_count <= 256) ? tight_stack : CK_ALLOC(sizeof(AABB) * body_count);
	AABB scene_bounds = aabb_empty();
	CK_DYNA SAP_Entry* sap = NULL;
	// Compute tight AABBs. Skip sleeping dynamic bodies (they don't move).
	CK_DYNA int* sleeping_bodies = NULL;
	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i) || asize(w->body_cold[i].shapes) == 0) { tight[i] = aabb_empty(); continue; }
		if (body_inv_mass(w, i) > 0.0f) {
			int isl = w->body_cold[i].island_id;
			if (isl >= 0 && (w->island_gen[isl] & 1) && !w->islands[isl].awake) {
				tight[i] = body_aabb(&w->body_state[i], &w->body_cold[i]); // needed for wake detection
				apush(sleeping_bodies, i);
				continue;
			}
		}
		tight[i] = body_aabb(&w->body_state[i], &w->body_cold[i]);
		if (body_inv_mass(w, i) > 0.0f) { scene_bounds = aabb_merge(scene_bounds, tight[i]); }
	}
	v3 extent = sub(scene_bounds.max, scene_bounds.min);
	int axis = (extent.y > extent.x && extent.y > extent.z) ? 1 : (extent.z > extent.x) ? 2 : 0;
	for (int i = 0; i < body_count; i++) {
		if (!split_alive(w->body_gen, i) || asize(w->body_cold[i].shapes) == 0) continue;
		if (body_inv_mass(w, i) == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && (w->island_gen[isl] & 1) && !w->islands[isl].awake) continue;
		apush(sap, ((SAP_Entry){ i, ((float*)&tight[i].min)[axis], ((float*)&tight[i].max)[axis] }));
	}

	// Sort awake dynamic bodies by chosen axis AABB min.
	int sap_count = asize(sap);
	if (sap_count > 1) sap_insertion_sort(sap, sap_count);
	bp_precomp_acc += perf_now() - t1;

	// Sweep: test overlapping awake-awake pairs.
	double t2 = perf_now();
	CK_DYNA BroadPair* dd_pairs = NULL;
	for (int i = 0; i < sap_count; i++) {
		float max_val = sap[i].max_val;
		int a = sap[i].body_idx;
		AABB ta = tight[a];
		for (int j = i + 1; j < sap_count && sap[j].min_val <= max_val; j++) {
			int b = sap[j].body_idx;
			if (!aabb_overlaps(ta, tight[b])) continue;
			if (jointed_pair_skip(w->joint_pairs, a, b)) continue;
			if (filter_pair_skip(w, a, b)) continue;
			apush(dd_pairs, ((BroadPair){ a, b }));
		}
	}
	// Wake detection: test each awake body against sleeping bodies.
	// O(awake * sleeping) AABB tests, but fast with SIMD branchless overlap.
	int n_sleeping = asize(sleeping_bodies);
	for (int i = 0; i < sap_count; i++) {
		int a = sap[i].body_idx;
		AABB ta = tight[a];
		for (int si = 0; si < n_sleeping; si++) {
			int b = sleeping_bodies[si];
			if (!aabb_overlaps(ta, tight[b])) continue;
			if (jointed_pair_skip(w->joint_pairs, a, b)) continue;
			if (filter_pair_skip(w, a, b)) continue;
			apush(dd_pairs, ((BroadPair){ a, b }));
		}
	}
	afree(sleeping_bodies);

	// Collect awake-dynamic vs static pairs.
	bp_sweep_acc += perf_now() - t2;
	double t3 = perf_now();
	CK_DYNA BroadPair* pairs = NULL;
	bvh_cross_test(w->bvh_dynamic, w->bvh_static, &pairs);
	for (int i = 0; i < asize(pairs); i++) {
		int a = pairs[i].a, b = pairs[i].b;
		if (!aabb_overlaps(tight[a], tight[b])) continue;
		// Skip sleeping dynamic vs static
		int isl_a = w->body_cold[a].island_id;
		if (isl_a >= 0 && (w->island_gen[isl_a] & 1) && !w->islands[isl_a].awake) continue;
		apush(dd_pairs, ((BroadPair){ a, b }));
	}
	afree(pairs);
	bp_cross_acc += perf_now() - t3;

	// Narrowphase on all collected pairs.
	// If parallel dispatch is available (thread_count > 1), output pairs for external dispatch.
	// Otherwise run narrowphase inline.
	if (w->thread_count > 1 && w->np_pairs_out) {
		// Output pairs for parallel narrowphase in nudge.c
		CK_DYNA BroadPair** out = (CK_DYNA BroadPair**)w->np_pairs_out;
		for (int i = 0; i < asize(dd_pairs); i++) apush(*out, dd_pairs[i]);
	} else {
		int dd_count = asize(dd_pairs);
		for (int i = 0; i < dd_count; i++) {
			narrowphase_pair(w, dd_pairs[i].a, dd_pairs[i].b, manifolds);
		}
	}
	afree(dd_pairs);
	bp_frame_count++;

	if (tight != tight_stack) CK_FREE(tight);
	afree(sap);
}

// Rebuild joint_pairs map when joint topology changes (create/destroy joint).
// O(joint_count) build, O(1) lookup per broadphase pair.
static void joint_pairs_refresh(WorldInternal* w)
{
	if (w->joint_pairs_version == w->ldl_topo_version) return;
	map_free(w->joint_pairs);
	w->joint_pairs = NULL;
	int jcount = asize(w->joints);
	for (int i = 0; i < jcount; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		uint64_t key = body_pair_key(j->body_a, j->body_b);
		map_set(w->joint_pairs, key, (uint8_t)1);
	}
	w->joint_pairs_version = w->ldl_topo_version;
}

static void broadphase_and_collide(WorldInternal* w, InternalManifold** manifolds)
{
	joint_pairs_refresh(w);
	if (w->broadphase_type == BROADPHASE_BVH) broadphase_bvh(w, manifolds);
	else broadphase_n2(w, manifolds);
}
