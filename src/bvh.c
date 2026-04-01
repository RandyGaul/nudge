// See LICENSE for licensing info.
// bvh.c -- binary BVH broadphase (Bepu-inspired, 64-byte cache-line nodes).

// -----------------------------------------------------------------------------
// AABB type and helpers.

typedef struct AABB { v3 min, max; } AABB;

static AABB aabb_empty() { return (AABB){ V3(1e18f, 1e18f, 1e18f), V3(-1e18f, -1e18f, -1e18f) }; }

static AABB aabb_merge(AABB a, AABB b) { return (AABB){ v3_min(a.min, b.min), v3_max(a.max, b.max) }; }

static AABB aabb_from_point(v3 p) { return (AABB){ p, p }; }

static AABB aabb_expand(AABB a, float margin) { v3 m = V3(margin, margin, margin); return (AABB){ sub(a.min, m), add(a.max, m) }; }

static float aabb_surface_area(AABB a) { v3 d = sub(a.max, a.min); return d.x*d.y + d.y*d.z + d.z*d.x; }

static int aabb_overlaps(AABB a, AABB b) { return a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y && a.min.z <= b.max.z && a.max.z >= b.min.z; }

// Compute world-space AABB for a single shape on a body.
static AABB shape_aabb(BodyHot* h, ShapeInternal* s)
{
	v3 world_pos = add(h->position, rotate(h->rotation, s->local_pos));
	switch (s->type) {
	case SHAPE_SPHERE: {
		v3 r = V3(s->sphere.radius, s->sphere.radius, s->sphere.radius);
		return (AABB){ sub(world_pos, r), add(world_pos, r) };
	}
	case SHAPE_CAPSULE: {
		v3 up = rotate(h->rotation, V3(0, s->capsule.half_height, 0));
		v3 p = sub(world_pos, up), q = add(world_pos, up);
		v3 r = V3(s->capsule.radius, s->capsule.radius, s->capsule.radius);
		return (AABB){ sub(v3_min(p, q), r), add(v3_max(p, q), r) };
	}
	case SHAPE_BOX: {
		// OBB -> AABB: project rotated half-extents onto each world axis.
		v3 e = s->box.half_extents;
		v3 ax = rotate(h->rotation, V3(e.x, 0, 0));
		v3 ay = rotate(h->rotation, V3(0, e.y, 0));
		v3 az = rotate(h->rotation, V3(0, 0, e.z));
		v3 half = V3(fabsf(ax.x) + fabsf(ay.x) + fabsf(az.x), fabsf(ax.y) + fabsf(ay.y) + fabsf(az.y), fabsf(ax.z) + fabsf(ay.z) + fabsf(az.z));
		return (AABB){ sub(world_pos, half), add(world_pos, half) };
	}
	case SHAPE_HULL: {
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		AABB box = aabb_from_point(add(world_pos, rotate(h->rotation, V3(hull->verts[0].x*sc.x, hull->verts[0].y*sc.y, hull->verts[0].z*sc.z))));
		for (int i = 1; i < hull->vert_count; i++) {
			v3 v = add(world_pos, rotate(h->rotation, V3(hull->verts[i].x*sc.x, hull->verts[i].y*sc.y, hull->verts[i].z*sc.z)));
			box.min = v3_min(box.min, v);
			box.max = v3_max(box.max, v);
		}
		return box;
	}
	}
	return aabb_empty();
}

// Compute world-space AABB for an entire body (union of all shapes).
static AABB body_aabb(BodyHot* h, BodyCold* c)
{
	AABB box = shape_aabb(h, &c->shapes[0]);
	for (int i = 1; i < asize(c->shapes); i++)
		box = aabb_merge(box, shape_aabb(h, &c->shapes[i]));
	return box;
}

// -----------------------------------------------------------------------------
// BVH node layout: 64 bytes = one cache line.
// Two children inline. Metadata (parent pointers) stored separately.

typedef struct BVHChild
{
	v3 min;          // 12
	int32_t index;   // 4 -- negative: encoded leaf (~idx), non-negative: child node
	v3 max;          // 12
	int32_t leaf_count; // 4 -- 0=empty, 1=leaf, >1=internal subtree
} BVHChild;          // 32 bytes

typedef struct BVHNode
{
	BVHChild a; // 32
	BVHChild b; // 32
} BVHNode;      // 64 bytes

typedef struct BVHMeta
{
	int parent;     // parent node index, -1 for root
	int child_slot; // 0=A, 1=B
} BVHMeta;

typedef struct BVHLeaf
{
	int body_idx;
	int node_idx;
	int child_slot; // 0=A, 1=B
	v3 fat_min, fat_max; // expanded AABB for motion threshold
} BVHLeaf;

typedef struct BVHTree
{
	CK_DYNA BVHNode* nodes;
	CK_DYNA BVHMeta* meta;
	CK_DYNA BVHLeaf* leaves;
	CK_DYNA int* node_free;
	CK_DYNA int* leaf_free;
	int root; // -1 = empty
} BVHTree;

typedef struct BroadPair { int a, b; } BroadPair;

#define BVH_AABB_MARGIN 0.05f
#define BVH_FAT_MARGIN  0.2f

// -----------------------------------------------------------------------------
// Tree lifecycle.

static void bvh_init(BVHTree* t) { memset(t, 0, sizeof(*t)); t->root = -1; }

static void bvh_free(BVHTree* t) { afree(t->nodes); afree(t->meta); afree(t->leaves); afree(t->node_free); afree(t->leaf_free); }

// Alloc/free with freelist.
static int bvh_alloc_node(BVHTree* t)
{
	int idx;
	if (asize(t->node_free) > 0) { idx = apop(t->node_free); }
	else { idx = asize(t->nodes); BVHNode z = {0}; apush(t->nodes, z); BVHMeta zm = {-1, 0}; apush(t->meta, zm); }
	t->nodes[idx] = (BVHNode){0};
	t->meta[idx] = (BVHMeta){-1, 0};
	return idx;
}

static void bvh_free_node(BVHTree* t, int idx) { apush(t->node_free, idx); }

static int bvh_alloc_leaf(BVHTree* t)
{
	int idx;
	if (asize(t->leaf_free) > 0) { idx = apop(t->leaf_free); }
	else { idx = asize(t->leaves); BVHLeaf z = {0}; apush(t->leaves, z); }
	return idx;
}

static void bvh_free_leaf(BVHTree* t, int idx) { apush(t->leaf_free, idx); }

// -----------------------------------------------------------------------------
// Child helpers.

static BVHChild* bvh_child(BVHNode* n, int slot) { return slot == 0 ? &n->a : &n->b; }

static AABB bvh_child_aabb(BVHChild* c) { return (AABB){ c->min, c->max }; }

static void bvh_child_set_aabb(BVHChild* c, AABB box) { c->min = box.min; c->max = box.max; }

static void bvh_child_set_leaf(BVHChild* c, AABB box, int leaf_idx) { bvh_child_set_aabb(c, box); c->index = ~leaf_idx; c->leaf_count = 1; }

// Place a leaf with fat AABB. The node child stores the fat AABB for conservative overlap.
static void bvh_place_leaf(BVHTree* t, int ni, int slot, AABB tight, int leaf_idx)
{
	AABB fat = aabb_expand(tight, BVH_FAT_MARGIN);
	bvh_child_set_leaf(bvh_child(&t->nodes[ni], slot), fat, leaf_idx);
	t->leaves[leaf_idx].node_idx = ni;
	t->leaves[leaf_idx].child_slot = slot;
	t->leaves[leaf_idx].fat_min = fat.min;
	t->leaves[leaf_idx].fat_max = fat.max;
}

static void bvh_child_set_node(BVHChild* c, AABB box, int node_idx, int lcount) { bvh_child_set_aabb(c, box); c->index = node_idx; c->leaf_count = lcount; }

static void bvh_child_set_empty(BVHChild* c) { *c = (BVHChild){0}; }

static int bvh_child_is_leaf(BVHChild* c) { return c->leaf_count == 1; }
static int bvh_child_is_internal(BVHChild* c) { return c->leaf_count > 1; }
static int bvh_child_is_empty(BVHChild* c) { return c->leaf_count == 0; }
static int bvh_child_leaf_idx(BVHChild* c) { return ~c->index; }

static AABB bvh_node_aabb(BVHNode* n)
{
	if (n->a.leaf_count == 0) return bvh_child_aabb(&n->b);
	if (n->b.leaf_count == 0) return bvh_child_aabb(&n->a);
	return aabb_merge(bvh_child_aabb(&n->a), bvh_child_aabb(&n->b));
}

// Count leaves reachable from a node (for validation).
static int bvh_count_leaves(BVHTree* t, int ni)
{
	BVHNode* n = &t->nodes[ni];
	int count = 0;
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		if (bvh_child_is_leaf(c)) count++;
		else count += bvh_count_leaves(t, c->index);
	}
	return count;
}

// -----------------------------------------------------------------------------
// Tree rotation: swap a grandchild with its uncle to reduce SAH cost.
// Evaluates 4 candidates at a node and applies the best improvement.

static void bvh_try_rotate(BVHTree* t, int ni)
{
	BVHNode* n = &t->nodes[ni];
	if (bvh_child_is_empty(&n->a) || bvh_child_is_empty(&n->b)) return;

	float base_cost = aabb_surface_area(bvh_child_aabb(&n->a)) + aabb_surface_area(bvh_child_aabb(&n->b));
	float best_cost = base_cost;
	int best_rot = -1; // 0-3: which rotation to apply

	// For each internal child, try swapping each of its grandchildren with the other child.
	for (int side = 0; side < 2; side++) {
		BVHChild* parent_child = bvh_child(n, side);
		BVHChild* uncle = bvh_child(n, 1 - side);
		if (!bvh_child_is_internal(parent_child)) continue;

		BVHNode* pn = &t->nodes[parent_child->index];
		for (int gc = 0; gc < 2; gc++) {
			BVHChild* grandchild = bvh_child(pn, gc);
			BVHChild* sibling = bvh_child(pn, 1 - gc);
			if (bvh_child_is_empty(grandchild) || bvh_child_is_empty(sibling)) continue;

			// After rotation: parent_child would contain (sibling, uncle_old), uncle slot gets grandchild.
			AABB new_parent_aabb = aabb_merge(bvh_child_aabb(sibling), bvh_child_aabb(uncle));
			float new_parent_sa = aabb_surface_area(new_parent_aabb);
			float new_uncle_sa = aabb_surface_area(bvh_child_aabb(grandchild));
			float cost = new_parent_sa + new_uncle_sa;

			if (cost < best_cost - 1e-6f) {
				best_cost = cost;
				best_rot = side * 2 + gc;
			}
		}
	}

	if (best_rot < 0) return;

	// Apply the rotation.
	int side = best_rot / 2;
	int gc = best_rot % 2;
	BVHChild* parent_child = bvh_child(n, side);
	BVHChild* uncle_slot = bvh_child(n, 1 - side);
	int pn_idx = parent_child->index;
	BVHNode* pn = &t->nodes[pn_idx];

	// Save grandchild and uncle.
	BVHChild gc_saved = *bvh_child(pn, gc);
	BVHChild uncle_saved = *uncle_slot;

	// Grandchild moves to uncle's slot in N.
	*uncle_slot = gc_saved;
	// Uncle moves into grandchild's slot in parent_child's node.
	*bvh_child(pn, gc) = uncle_saved;

	// Update parent_child's AABB and leaf_count in N.
	bvh_child_set_aabb(parent_child, bvh_node_aabb(pn));
	parent_child->leaf_count = pn->a.leaf_count + pn->b.leaf_count;

	// Fix back-pointers for the swapped children.
	// The grandchild (now in uncle_slot at node ni):
	if (bvh_child_is_leaf(uncle_slot)) {
		t->leaves[bvh_child_leaf_idx(uncle_slot)].node_idx = ni;
		t->leaves[bvh_child_leaf_idx(uncle_slot)].child_slot = 1 - side;
	} else if (bvh_child_is_internal(uncle_slot)) {
		t->meta[uncle_slot->index].parent = ni;
		t->meta[uncle_slot->index].child_slot = 1 - side;
	}
	// The uncle (now in pn at slot gc):
	BVHChild* new_gc = bvh_child(pn, gc);
	if (bvh_child_is_leaf(new_gc)) {
		t->leaves[bvh_child_leaf_idx(new_gc)].node_idx = pn_idx;
		t->leaves[bvh_child_leaf_idx(new_gc)].child_slot = gc;
	} else if (bvh_child_is_internal(new_gc)) {
		t->meta[new_gc->index].parent = pn_idx;
		t->meta[new_gc->index].child_slot = gc;
	}
}

// -----------------------------------------------------------------------------
// Insertion: SAH-guided top-down.

static int bvh_insert(BVHTree* t, int body_idx, AABB bounds)
{
	int li = bvh_alloc_leaf(t);
	t->leaves[li] = (BVHLeaf){ .body_idx = body_idx };

	// Empty tree: create root with leaf in child A.
	if (t->root == -1) {
		int ni = bvh_alloc_node(t);
		t->meta[ni] = (BVHMeta){-1, 0};
		bvh_place_leaf(t, ni, 0, bounds, li);
		t->root = ni;
		return li;
	}

	// Single leaf in tree: root has child A occupied and B empty (or vice versa).
	BVHNode* root = &t->nodes[t->root];
	if (bvh_child_is_empty(&root->b) && bvh_child_is_leaf(&root->a)) {
		bvh_place_leaf(t, t->root, 1, bounds, li);
		return li;
	}
	if (bvh_child_is_empty(&root->a) && bvh_child_is_leaf(&root->b)) {
		bvh_place_leaf(t, t->root, 0, bounds, li);
		return li;
	}

	// Walk tree from root to find best insertion point.
	int cur = t->root;
	for (;;) {
		BVHNode* node = &t->nodes[cur];
		// Evaluate cost of inserting into child A vs B.
		AABB ma = aabb_merge(bvh_child_aabb(&node->a), bounds);
		AABB mb = aabb_merge(bvh_child_aabb(&node->b), bounds);
		float cost_a = aabb_surface_area(ma) - aabb_surface_area(bvh_child_aabb(&node->a));
		float cost_b = aabb_surface_area(mb) - aabb_surface_area(bvh_child_aabb(&node->b));

		int slot = cost_a <= cost_b ? 0 : 1;
		BVHChild* child = bvh_child(node, slot);

		if (bvh_child_is_leaf(child)) {
			// Split: create new node containing existing leaf + new leaf.
			int existing_li = bvh_child_leaf_idx(child);
			AABB existing_aabb = bvh_child_aabb(child);
			AABB merged = slot == 0 ? ma : mb;

			int new_ni = bvh_alloc_node(t);
			// NOTE: alloc may realloc t->nodes, so re-fetch pointers.

			// Existing leaf keeps its current (fat) AABB.
			bvh_child_set_leaf(&t->nodes[new_ni].a, existing_aabb, existing_li);
			t->leaves[existing_li].node_idx = new_ni;
			t->leaves[existing_li].child_slot = 0;

			// New leaf gets fat AABB.
			bvh_place_leaf(t, new_ni, 1, bounds, li);

			// Replace child in parent with internal node pointer (re-fetch after alloc).
			bvh_child_set_node(bvh_child(&t->nodes[cur], slot), merged, new_ni, 2);
			t->meta[new_ni] = (BVHMeta){ cur, slot };

			// Walk back up refitting bounds and applying rotations.
			int ri = cur;
			while (ri != -1) {
				BVHNode* rn = &t->nodes[ri];
				if (bvh_child_is_internal(&rn->a)) { BVHNode* cn = &t->nodes[rn->a.index]; rn->a.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->a, bvh_node_aabb(cn)); }
				if (bvh_child_is_internal(&rn->b)) { BVHNode* cn = &t->nodes[rn->b.index]; rn->b.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->b, bvh_node_aabb(cn)); }
				bvh_try_rotate(t, ri);
				ri = t->meta[ri].parent;
			}
			return li;
		}

		if (bvh_child_is_internal(child)) {
			// Descend into the child subtree.
			cur = child->index;
			continue;
		}

		// Empty slot: place leaf directly.
		bvh_place_leaf(t, cur, slot, bounds, li);
		return li;
	}
}

// -----------------------------------------------------------------------------
// Removal.

static void bvh_remove(BVHTree* t, int leaf_idx)
{
	BVHLeaf* leaf = &t->leaves[leaf_idx];
	int ni = leaf->node_idx;
	int slot = leaf->child_slot;

	// Clear the child slot.
	bvh_child_set_empty(bvh_child(&t->nodes[ni], slot));

	// Get the sibling.
	int sib_slot = 1 - slot;
	BVHChild sib = *bvh_child(&t->nodes[ni], sib_slot);

	int parent = t->meta[ni].parent;
	if (parent == -1) {
		// ni is root. If sibling is empty, tree becomes empty. Otherwise keep root with one child.
		if (bvh_child_is_empty(&sib)) { bvh_free_node(t, ni); t->root = -1; }
		// else: root has one child remaining, that's fine.
	} else {
		// Replace ni in parent with sibling.
		int parent_slot = t->meta[ni].child_slot;
		BVHChild* pc = bvh_child(&t->nodes[parent], parent_slot);
		*pc = sib;

		// Update sibling's back-pointers.
		if (bvh_child_is_leaf(pc)) {
			int sib_li = bvh_child_leaf_idx(pc);
			t->leaves[sib_li].node_idx = parent;
			t->leaves[sib_li].child_slot = parent_slot;
		} else if (bvh_child_is_internal(pc)) {
			t->meta[pc->index].parent = parent;
			t->meta[pc->index].child_slot = parent_slot;
		}

		bvh_free_node(t, ni);

		// Refit upward.
		int ri = parent;
		while (ri != -1) {
			BVHNode* rn = &t->nodes[ri];
			if (bvh_child_is_internal(&rn->a)) { BVHNode* cn = &t->nodes[rn->a.index]; rn->a.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->a, bvh_node_aabb(cn)); }
			if (bvh_child_is_internal(&rn->b)) { BVHNode* cn = &t->nodes[rn->b.index]; rn->b.leaf_count = cn->a.leaf_count + cn->b.leaf_count; bvh_child_set_aabb(&rn->b, bvh_node_aabb(cn)); }
			ri = t->meta[ri].parent;
		}
	}

	bvh_free_leaf(t, leaf_idx);
}

// -----------------------------------------------------------------------------
// Binned SAH builder: top-down build from a set of leaf indices.
// lut maps leaf_index -> AABB (indexed by leaf index, not array position).

#define BVH_SAH_BINS 8

// Build a subtree from lis[0..count-1]. Each element is a leaf index.
// lut[leaf_idx] gives the AABB for that leaf.
// If count==1, returns ~leaf_idx (encoded as a "leaf result"). Caller checks sign.
// If count>=2, returns a node index.
static int bvh_binned_build(BVHTree* t, int* lis, AABB* lut, int count)
{
	assert(count >= 2);

	if (count == 2) {
		int ni = bvh_alloc_node(t);
		bvh_child_set_leaf(&t->nodes[ni].a, lut[lis[0]], lis[0]);
		bvh_child_set_leaf(&t->nodes[ni].b, lut[lis[1]], lis[1]);
		t->leaves[lis[0]].node_idx = ni; t->leaves[lis[0]].child_slot = 0;
		t->leaves[lis[1]].node_idx = ni; t->leaves[lis[1]].child_slot = 1;
		return ni;
	}

	// Compute centroid bounds.
	v3 cmin = V3(1e18f, 1e18f, 1e18f), cmax = V3(-1e18f, -1e18f, -1e18f);
	for (int i = 0; i < count; i++) {
		v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
		cmin = v3_min(cmin, c); cmax = v3_max(cmax, c);
	}

	// Find best split across all 3 axes.
	float best_cost = 1e18f;
	int best_axis = 0, best_split = 0;

	for (int axis = 0; axis < 3; axis++) {
		float amin = (&cmin.x)[axis], amax = (&cmax.x)[axis];
		float extent = amax - amin;
		if (extent < 1e-7f) continue;

		int bin_count[BVH_SAH_BINS] = {0};
		AABB bin_aabb[BVH_SAH_BINS];
		for (int b = 0; b < BVH_SAH_BINS; b++) bin_aabb[b] = aabb_empty();

		float inv_ext = (float)BVH_SAH_BINS / extent;
		for (int i = 0; i < count; i++) {
			v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
			int b = (int)(((&c.x)[axis] - amin) * inv_ext);
			if (b >= BVH_SAH_BINS) b = BVH_SAH_BINS - 1;
			bin_count[b]++; bin_aabb[b] = aabb_merge(bin_aabb[b], lut[lis[i]]);
		}

		// Prefix scan left and right.
		AABB la[BVH_SAH_BINS - 1], ra[BVH_SAH_BINS - 1];
		int lc[BVH_SAH_BINS - 1], rc[BVH_SAH_BINS - 1];
		AABB run = aabb_empty(); int rn = 0;
		for (int i = 0; i < BVH_SAH_BINS - 1; i++) { run = aabb_merge(run, bin_aabb[i]); rn += bin_count[i]; la[i] = run; lc[i] = rn; }
		run = aabb_empty(); rn = 0;
		for (int i = BVH_SAH_BINS - 1; i > 0; i--) { run = aabb_merge(run, bin_aabb[i]); rn += bin_count[i]; ra[i-1] = run; rc[i-1] = rn; }

		for (int i = 0; i < BVH_SAH_BINS - 1; i++) {
			if (lc[i] == 0 || rc[i] == 0) continue;
			float cost = aabb_surface_area(la[i]) * lc[i] + aabb_surface_area(ra[i]) * rc[i];
			if (cost < best_cost) { best_cost = cost; best_axis = axis; best_split = i; }
		}
	}

	// Partition by best split.
	float amin = (&cmin.x)[best_axis], amax = (&cmax.x)[best_axis];
	float extent = amax - amin;
	float inv_ext = extent > 1e-7f ? (float)BVH_SAH_BINS / extent : 0.0f;

	int mid = 0;
	for (int i = 0; i < count; i++) {
		v3 c = scale(add(lut[lis[i]].min, lut[lis[i]].max), 0.5f);
		int b = inv_ext > 0 ? (int)(((&c.x)[best_axis] - amin) * inv_ext) : 0;
		if (b >= BVH_SAH_BINS) b = BVH_SAH_BINS - 1;
		if (b <= best_split) { int tmp = lis[mid]; lis[mid] = lis[i]; lis[i] = tmp; mid++; }
	}
	if (mid == 0 || mid == count) mid = count / 2;

	// Build children. Each half has count >= 1.
	int ni = bvh_alloc_node(t);

	// Left child.
	if (mid == 1) {
		bvh_child_set_leaf(&t->nodes[ni].a, lut[lis[0]], lis[0]);
		t->leaves[lis[0]].node_idx = ni; t->leaves[lis[0]].child_slot = 0;
	} else {
		int left = bvh_binned_build(t, lis, lut, mid);
		BVHNode* ln = &t->nodes[left];
		bvh_child_set_node(&t->nodes[ni].a, bvh_node_aabb(ln), left, ln->a.leaf_count + ln->b.leaf_count);
		t->meta[left].parent = ni; t->meta[left].child_slot = 0;
	}

	// Right child.
	int rcount = count - mid;
	if (rcount == 1) {
		bvh_child_set_leaf(&t->nodes[ni].b, lut[lis[mid]], lis[mid]);
		t->leaves[lis[mid]].node_idx = ni; t->leaves[lis[mid]].child_slot = 1;
	} else {
		int right = bvh_binned_build(t, lis + mid, lut, rcount);
		BVHNode* rnode = &t->nodes[right];
		bvh_child_set_node(&t->nodes[ni].b, bvh_node_aabb(rnode), right, rnode->a.leaf_count + rnode->b.leaf_count);
		t->meta[right].parent = ni; t->meta[right].child_slot = 1;
	}

	return ni;
}

// Collect all leaf indices reachable from node ni into lis array. Frees internal nodes.
static void bvh_collect_and_free(BVHTree* t, int ni, CK_DYNA int** lis)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		if (bvh_child_is_leaf(c)) apush(*lis, bvh_child_leaf_idx(c));
		else bvh_collect_and_free(t, c->index, lis);
	}
	bvh_free_node(t, ni);
}

// Rebuild the subtree at node ni using binned SAH.
// parent_ni and parent_slot describe where ni is attached (-1 for root).
static void bvh_refine_subtree(BVHTree* t, int ni, AABB* lut)
{
	int parent_ni = t->meta[ni].parent;
	int parent_slot = t->meta[ni].child_slot;

	// Collect leaves and free old internal nodes.
	CK_DYNA int* lis = NULL;
	bvh_collect_and_free(t, ni, &lis);
	int count = asize(lis);

	if (count < 2) {
		// Edge case: 0 or 1 leaves. Just put the leaf back.
		if (count == 1 && parent_ni >= 0) {
			bvh_child_set_leaf(bvh_child(&t->nodes[parent_ni], parent_slot), lut[lis[0]], lis[0]);
			t->leaves[lis[0]].node_idx = parent_ni;
			t->leaves[lis[0]].child_slot = parent_slot;
		}
		afree(lis);
		return;
	}

	// Rebuild.
	int new_root = bvh_binned_build(t, lis, lut, count);

	// Splice into parent.
	if (parent_ni == -1) {
		t->root = new_root;
		t->meta[new_root].parent = -1;
	} else {
		BVHNode* cn = &t->nodes[new_root];
		AABB box = bvh_node_aabb(cn);
		bvh_child_set_node(bvh_child(&t->nodes[parent_ni], parent_slot), box, new_root, cn->a.leaf_count + cn->b.leaf_count);
		t->meta[new_root].parent = parent_ni;
		t->meta[new_root].child_slot = parent_slot;
	}

	afree(lis);
}

// Per-frame incremental refinement: rebuild ~fraction of the tree using binned SAH.
static void bvh_incremental_refine(BVHTree* t, AABB* lut, float fraction)
{
	if (t->root == -1) return;
	int total_leaves = bvh_count_leaves(t, t->root);
	if (total_leaves < 4) return; // too small to bother

	int target_size = (int)sqrtf((float)total_leaves);
	if (target_size < 2) target_size = 2;

	// Walk nodes to find a subtree near the target size.
	int best_ni = -1, best_diff = 0x7fffffff;
	for (int i = 0; i < asize(t->nodes); i++) {
		BVHNode* n = &t->nodes[i];
		int lc = n->a.leaf_count + n->b.leaf_count;
		if (lc < 2) continue;
		int diff = abs(lc - target_size);
		if (diff < best_diff) { best_diff = diff; best_ni = i; }
	}

	if (best_ni >= 0) bvh_refine_subtree(t, best_ni, lut);
	(void)fraction;
}

// -----------------------------------------------------------------------------
// DFS cache reorder: compact live nodes into DFS traversal order for cache locality.

static void bvh_dfs_collect(BVHTree* t, int ni, int* remap, BVHNode* dst_nodes, BVHMeta* dst_meta, int* cursor)
{
	int new_idx = (*cursor)++;
	remap[ni] = new_idx;
	dst_nodes[new_idx] = t->nodes[ni];
	dst_meta[new_idx] = t->meta[ni];
	BVHNode* n = &t->nodes[ni];
	if (bvh_child_is_internal(&n->a)) bvh_dfs_collect(t, n->a.index, remap, dst_nodes, dst_meta, cursor);
	if (bvh_child_is_internal(&n->b)) bvh_dfs_collect(t, n->b.index, remap, dst_nodes, dst_meta, cursor);
}

static void bvh_cache_reorder(BVHTree* t)
{
	if (t->root == -1) return;

	int cap = asize(t->nodes);
	int* remap = CK_ALLOC(sizeof(int) * cap);
	for (int i = 0; i < cap; i++) remap[i] = -1;

	BVHNode* new_nodes = CK_ALLOC(sizeof(BVHNode) * cap);
	BVHMeta* new_meta = CK_ALLOC(sizeof(BVHMeta) * cap);
	int cursor = 0;

	bvh_dfs_collect(t, t->root, remap, new_nodes, new_meta, &cursor);
	int live_count = cursor;

	// Fix child indices in new_nodes.
	for (int i = 0; i < live_count; i++) {
		BVHNode* n = &new_nodes[i];
		if (bvh_child_is_internal(&n->a)) n->a.index = remap[n->a.index];
		if (bvh_child_is_internal(&n->b)) n->b.index = remap[n->b.index];
	}

	// Fix parent indices in new_meta.
	for (int i = 0; i < live_count; i++) {
		if (new_meta[i].parent >= 0) new_meta[i].parent = remap[new_meta[i].parent];
	}

	// Fix leaf back-pointers.
	for (int i = 0; i < asize(t->leaves); i++) {
		if (remap[t->leaves[i].node_idx] >= 0) t->leaves[i].node_idx = remap[t->leaves[i].node_idx];
	}

	// Swap into the ckit dynamic arrays. Copy live nodes into t->nodes, resize to live_count.
	// We can't just swap pointers because ckit tracks capacity in the array header.
	aclear(t->nodes); aclear(t->meta);
	for (int i = 0; i < live_count; i++) { apush(t->nodes, new_nodes[i]); apush(t->meta, new_meta[i]); }

	// Rebuild freelist: empty since all slots are live and packed.
	aclear(t->node_free);

	t->root = remap[t->root];

	CK_FREE(remap);
	CK_FREE(new_nodes);
	CK_FREE(new_meta);
}

// -----------------------------------------------------------------------------
// Refit: update all leaf AABBs from current body poses, propagate up.

// Returns 1 if anything changed in this subtree, 0 if all leaves stayed within fat AABBs.
static int bvh_refit_recurse(BVHTree* t, int ni, WorldInternal* w)
{
	BVHNode* n = &t->nodes[ni];
	int changed = 0;
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		if (bvh_child_is_leaf(c)) {
			int li = bvh_child_leaf_idx(c);
			int bi = t->leaves[li].body_idx;
			AABB tight = aabb_expand(body_aabb(&w->body_hot[bi], &w->body_cold[bi]), BVH_AABB_MARGIN);
			// Check if tight AABB is still within fat AABB.
			v3 fmin = t->leaves[li].fat_min, fmax = t->leaves[li].fat_max;
			if (tight.min.x >= fmin.x && tight.min.y >= fmin.y && tight.min.z >= fmin.z && tight.max.x <= fmax.x && tight.max.y <= fmax.y && tight.max.z <= fmax.z)
				continue; // still within fat -- skip
			// Re-expand fat AABB.
			AABB fat = aabb_expand(tight, BVH_FAT_MARGIN);
			t->leaves[li].fat_min = fat.min;
			t->leaves[li].fat_max = fat.max;
			bvh_child_set_aabb(c, fat);
			changed = 1;
		} else {
			int sub = bvh_refit_recurse(t, c->index, w);
			if (sub) {
				BVHNode* cn = &t->nodes[c->index];
				bvh_child_set_aabb(c, bvh_node_aabb(cn));
				c->leaf_count = cn->a.leaf_count + cn->b.leaf_count;
				changed = 1;
			}
		}
	}
	return changed;
}

static void bvh_refit(BVHTree* t, WorldInternal* w)
{
	if (t->root == -1) return;
	bvh_refit_recurse(t, t->root, w);
	bvh_cache_reorder(t);
}

// -----------------------------------------------------------------------------
// Self-test: find all overlapping leaf pairs within a single tree.

static void bvh_test_children(BVHTree* t, BVHChild* a, BVHChild* b, CK_DYNA BroadPair** pairs)
{
	if (bvh_child_is_empty(a) || bvh_child_is_empty(b)) return;
	if (!aabb_overlaps(bvh_child_aabb(a), bvh_child_aabb(b))) return;

	if (bvh_child_is_leaf(a) && bvh_child_is_leaf(b)) {
		BroadPair p = { t->leaves[bvh_child_leaf_idx(a)].body_idx, t->leaves[bvh_child_leaf_idx(b)].body_idx };
		apush(*pairs, p);
		return;
	}

	// Expand the larger child (by leaf count).
	if (bvh_child_is_leaf(a) || (bvh_child_is_internal(b) && b->leaf_count > a->leaf_count)) {
		BVHNode* bn = &t->nodes[b->index];
		bvh_test_children(t, a, &bn->a, pairs);
		bvh_test_children(t, a, &bn->b, pairs);
	} else {
		BVHNode* an = &t->nodes[a->index];
		bvh_test_children(t, &an->a, b, pairs);
		bvh_test_children(t, &an->b, b, pairs);
	}
}

static void bvh_self_test_node(BVHTree* t, int ni, CK_DYNA BroadPair** pairs)
{
	BVHNode* n = &t->nodes[ni];
	// Test child A vs child B.
	bvh_test_children(t, &n->a, &n->b, pairs);
	// Recurse into internal children.
	if (bvh_child_is_internal(&n->a)) bvh_self_test_node(t, n->a.index, pairs);
	if (bvh_child_is_internal(&n->b)) bvh_self_test_node(t, n->b.index, pairs);
}

static void bvh_self_test(BVHTree* t, CK_DYNA BroadPair** pairs)
{
	if (t->root == -1) return;
	bvh_self_test_node(t, t->root, pairs);
}

// -----------------------------------------------------------------------------
// Cross-test: find all overlapping leaf pairs between two trees.

static void bvh_cross_children(BVHTree* ta, BVHChild* a, BVHTree* tb, BVHChild* b, CK_DYNA BroadPair** pairs)
{
	if (bvh_child_is_empty(a) || bvh_child_is_empty(b)) return;
	if (!aabb_overlaps(bvh_child_aabb(a), bvh_child_aabb(b))) return;

	if (bvh_child_is_leaf(a) && bvh_child_is_leaf(b)) {
		BroadPair p = { ta->leaves[bvh_child_leaf_idx(a)].body_idx, tb->leaves[bvh_child_leaf_idx(b)].body_idx };
		apush(*pairs, p);
		return;
	}

	if (bvh_child_is_leaf(a) || (bvh_child_is_internal(b) && b->leaf_count > a->leaf_count)) {
		BVHNode* bn = &tb->nodes[b->index];
		bvh_cross_children(ta, a, tb, &bn->a, pairs);
		bvh_cross_children(ta, a, tb, &bn->b, pairs);
	} else {
		BVHNode* an = &ta->nodes[a->index];
		bvh_cross_children(ta, &an->a, tb, b, pairs);
		bvh_cross_children(ta, &an->b, tb, b, pairs);
	}
}

static void bvh_cross_test(BVHTree* ta, BVHTree* tb, CK_DYNA BroadPair** pairs)
{
	if (ta->root == -1 || tb->root == -1) return;
	BVHNode* ra = &ta->nodes[ta->root];
	BVHNode* rb = &tb->nodes[tb->root];
	// Test all 4 root child combinations.
	bvh_cross_children(ta, &ra->a, tb, &rb->a, pairs);
	bvh_cross_children(ta, &ra->a, tb, &rb->b, pairs);
	bvh_cross_children(ta, &ra->b, tb, &rb->a, pairs);
	bvh_cross_children(ta, &ra->b, tb, &rb->b, pairs);
}
