// See LICENSE for licensing info.
// epa.c -- Expanding Polytope Algorithm for penetration depth / contact manifold.
//
// Alternative narrowphase backend: run GJK to termination, seed EPA from the
// resulting simplex, expand outward until the closest origin-face converges on
// the penetration axis, and read out a single world-space contact via
// barycentric projection.
//
// Compared to SAT, EPA yields one contact per frame. Incremental manifolds
// accumulate contacts over frames in a per-pair cache (see nudge_internal.h
// EpaManifold).
//
// Data structures are local to this file. Polytope is hand-rolled (not
// quickhull) but mirrors quickhull.c's horizon walk / weld / validation
// patterns (see CLAUDE.md plan for the reference table).
//
// Half-edge convention (same as quickhull / nudge Hull):
//   Edges stored in twin pairs at 2k / 2k+1. Each edge stores its origin
//   (tail) vertex, and the head is origin(next(e)). Face cycle: three edges
//   per triangle connected via next[].

// -----------------------------------------------------------------------------
// Constants.

#define EPA_MAX_VERTS       96
#define EPA_MAX_FACES       192
#define EPA_MAX_EDGES       (EPA_MAX_FACES * 3)
#define EPA_MAX_ITERATIONS  48
#define EPA_TOLERANCE       1e-4f
// EPA_DEEP_THRESHOLD: gjk_distance() returns positive values for separated
// shapes (distance between witnesses). Below this, we treat the pair as
// overlapping and the EPA expansion path is eligible. Combined with
// EPA_SHALLOW_MIN_DIST below, this defines the 3-way dispatch:
//   r.distance > EPA_DEEP_THRESHOLD          → no hit (separated)
//   r.distance in (EPA_SHALLOW_MIN_DIST, threshold] → shallow path (GJK witness)
//   r.distance <= EPA_SHALLOW_MIN_DIST       → EPA expansion (true depth)
#define EPA_DEEP_THRESHOLD  1e-4f
#define EPA_DEGENERATE_EPS  1e-8f    // face-normal length^2 floor

// -----------------------------------------------------------------------------
// Polytope types.

typedef struct EpaVert
{
	v3 mink;   // Minkowski difference: pb - pa
	v3 pa;     // support on A (world space)
	v3 pb;     // support on B (world space)
	int fa;    // feature ID on A
	int fb;    // feature ID on B
} EpaVert;

typedef struct EpaFace
{
	int v[3];     // vertex indices
	int e[3];     // half-edge indices (one per triangle side)
	v3 normal;    // outward plane normal
	float dist;   // dot(normal, v[0]) -- distance from origin along outward normal
	int alive;    // 1 = alive, 0 = removed
} EpaFace;

typedef struct EpaEdge
{
	int v0;       // origin (tail) vertex
	int face;     // owning face
	int next;     // next half-edge in face cycle
	// Twin is implicit: twin(i) = i ^ 1.
} EpaEdge;

typedef struct EpaPoly
{
	EpaVert verts[EPA_MAX_VERTS];
	int vert_count;
	EpaFace faces[EPA_MAX_FACES];
	int face_count;
	EpaEdge edges[EPA_MAX_EDGES];
	int edge_count;
	float epsilon;   // weld / orientation epsilon
} EpaPoly;

// EPA output (single contact, world-space).
typedef struct EpaHit
{
	v3 point_a;     // world-space contact on A
	v3 point_b;     // world-space contact on B
	v3 normal;      // from A toward B (i.e. direction along which A is pushed)
	float depth;
	uint32_t feature_id;
} EpaHit;

// -----------------------------------------------------------------------------
// Polytope helpers.

static inline int epa_twin(int e) { return e ^ 1; }

// Allocate 2 adjacent half-edge slots for a new twin pair. Returns first slot, or -1.
static int epa_alloc_edge_pair(EpaPoly* p)
{
	if (p->edge_count + 2 > EPA_MAX_EDGES) return -1;
	int e = p->edge_count;
	p->edge_count += 2;
	return e;
}

// -----------------------------------------------------------------------------
// Polytope seed: build an initial tetrahedron from the GJK termination simplex.
//
// Inputs: 1-4 Minkowski-difference points in `simplex`. We promote to 4.
// Support functions are provided via callbacks so this works across
// shape-pair backends.

typedef v3 (*EpaSupportFn)(void* ctx, v3 dir, int* out_feat);

typedef struct EpaSupportPair
{
	EpaSupportFn fa;
	void* ctxa;
	EpaSupportFn fb;
	void* ctxb;
} EpaSupportPair;

static void epa_support_pair(EpaSupportPair* sp, v3 dir, EpaVert* out)
{
	// Minkowski difference M = B - A. Support in dir = support(B, dir) - support(A, -dir).
	int fa, fb;
	v3 pa = sp->fa(sp->ctxa, neg(dir), &fa);
	v3 pb = sp->fb(sp->ctxb, dir, &fb);
	out->pa = pa;
	out->pb = pb;
	out->mink = sub(pb, pa);
	out->fa = fa;
	out->fb = fb;
}

// Build a tetrahedron containing the origin. If `warm_dirs` is non-NULL and
// warm_valid is set, seed from those four Minkowski-space directions (they
// approximate the previous terminating face + interior) so the initial tetra
// is already close to convergence. Falls back to regular-tetrahedron cardinal
// directions on degeneracy.
static int epa_build_tetra(EpaPoly* p, EpaSupportPair* sp, const v3* warm_dirs, int warm_valid)
{
	p->vert_count = 0;

	// Four regular-tetrahedron directions (vertices of an inscribed tetra in a cube).
	static const v3 tetra_dirs[4] = {
		{ 1,  1,  1},
		{ 1, -1, -1},
		{-1,  1, -1},
		{-1, -1,  1},
	};

	// Select seed directions: warm or regular.
	v3 seed_dirs[4];
	int used_warm = 0;
	if (warm_valid && warm_dirs) {
		int ok = 1;
		for (int i = 0; i < 4; i++) {
			if (len2(warm_dirs[i]) < 1e-12f) { ok = 0; break; }
			seed_dirs[i] = warm_dirs[i];
		}
		if (ok) { used_warm = 1; }
	}
	if (!used_warm) {
		for (int i = 0; i < 4; i++) seed_dirs[i] = tetra_dirs[i];
	}

	for (int i = 0; i < 4; i++) {
		EpaVert v;
		epa_support_pair(sp, seed_dirs[i], &v);
		// Weld check: if duplicates one existing vertex, try an alternative direction.
		int dup = 0;
		for (int j = 0; j < p->vert_count; j++) {
			if (len2(sub(v.mink, p->verts[j].mink)) <= p->epsilon * p->epsilon) { dup = 1; break; }
		}
		if (dup) {
			// Try a randomized direction by combining with cardinal axes.
			v3 alt_dirs[6] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
			for (int k = 0; k < 6 && dup; k++) {
				epa_support_pair(sp, alt_dirs[k], &v);
				dup = 0;
				for (int j = 0; j < p->vert_count; j++) {
					if (len2(sub(v.mink, p->verts[j].mink)) <= p->epsilon * p->epsilon) { dup = 1; break; }
				}
			}
			if (dup) return 0;
		}
		p->verts[p->vert_count++] = v;
	}

	// Ensure the 4 vertices form a non-degenerate tetrahedron. If warm-start
	// produced a degenerate shape, caller retries with warm_valid=0.
	{
		v3 a = p->verts[0].mink, b = p->verts[1].mink, c = p->verts[2].mink, d = p->verts[3].mink;
		v3 ab = sub(b, a), ac = sub(c, a), ad = sub(d, a);
		float vol = dot(ab, cross(ac, ad));
		if (fabsf(vol) < 1e-8f) return 0;
	}

	// Build tetrahedron faces. The tetra interior must contain the origin for
	// EPA to work (otherwise GJK would not have reported penetration).
	//
	// Use consistent winding across all 4 faces. The natural winding of
	// {0,1,2},{0,2,3},{0,3,1},{1,3,2} is outward-consistent iff vertex 3 lies
	// on the negative side of face {0,1,2}'s geometric normal. If not, flip
	// the whole tetrahedron's winding by swapping verts 2 and 3 (which flips
	// every face's winding).
	v3 va0 = p->verts[0].mink, vb0 = p->verts[1].mink, vc0 = p->verts[2].mink, vd0 = p->verts[3].mink;
	v3 n012 = cross(sub(vb0, va0), sub(vc0, va0));
	if (dot(n012, sub(vd0, va0)) > 0.0f) {
		// Swap verts 2 and 3 so the standard face table has consistent outward winding.
		EpaVert tmp = p->verts[2]; p->verts[2] = p->verts[3]; p->verts[3] = tmp;
	}

	int face_verts[4][3] = {
		{ 0, 1, 2 },
		{ 0, 2, 3 },
		{ 0, 3, 1 },
		{ 1, 3, 2 },
	};

	p->face_count = 0;
	p->edge_count = 0;

	// Phase 1: set up face vertex indices and compute planes.
	for (int fi = 0; fi < 4; fi++) {
		EpaFace* f = &p->faces[fi];
		f->v[0] = face_verts[fi][0];
		f->v[1] = face_verts[fi][1];
		f->v[2] = face_verts[fi][2];
		f->alive = 1;
		v3 va = p->verts[f->v[0]].mink, vb = p->verts[f->v[1]].mink, vc = p->verts[f->v[2]].mink;
		v3 n = cross(sub(vb, va), sub(vc, va));
		float l2 = len2(n);
		if (l2 < EPA_DEGENERATE_EPS) return 0;
		n = scale(n, 1.0f / sqrtf(l2));
		float d = dot(n, va);
		// Ensure outward (dist >= 0). EPA expects all face normals to point away from origin.
		if (d < 0.0f) { n = neg(n); d = -d; }
		f->normal = n;
		f->dist = d;
	}
	p->face_count = 4;

	// Phase 2: allocate edges (3 per face = 12 total). Since each undirected
	// edge is shared by 2 faces, we get 6 twin pairs. We allocate 12 edge
	// slots (6 pairs), then for each face find the pair that contains its
	// directed edge.
	//
	// Strategy: first pass, collect all 12 directed edges as (u,v,face,slot)
	// without committing twin indices. Second pass, group (u,v) with (v,u)
	// and assign them to adjacent pair-slots (2k,2k+1).

	struct { int u, v, face, fedge; } dirs[12];
	int dn = 0;
	for (int fi = 0; fi < 4; fi++) {
		EpaFace* f = &p->faces[fi];
		for (int k = 0; k < 3; k++) {
			dirs[dn].u = f->v[k];
			dirs[dn].v = f->v[(k + 1) % 3];
			dirs[dn].face = fi;
			dirs[dn].fedge = k;
			dn++;
		}
	}

	// Assign each directed edge to a half-edge slot. For each pair (u,v)/(v,u),
	// place them at 2p and 2p+1.
	int assigned[12] = { 0 };
	int next_pair = 0;
	for (int i = 0; i < 12; i++) {
		if (assigned[i]) continue;
		int partner = -1;
		for (int j = i + 1; j < 12; j++) {
			if (assigned[j]) continue;
			if (dirs[j].u == dirs[i].v && dirs[j].v == dirs[i].u) { partner = j; break; }
		}
		if (partner < 0) return 0;
		int slot_i = next_pair * 2;
		int slot_j = next_pair * 2 + 1;
		next_pair++;
		assigned[i] = 1; assigned[partner] = 1;
		p->edges[slot_i] = (EpaEdge){ .v0 = dirs[i].u, .face = dirs[i].face, .next = -1 };
		p->edges[slot_j] = (EpaEdge){ .v0 = dirs[partner].u, .face = dirs[partner].face, .next = -1 };
		p->faces[dirs[i].face].e[dirs[i].fedge] = slot_i;
		p->faces[dirs[partner].face].e[dirs[partner].fedge] = slot_j;
	}
	p->edge_count = next_pair * 2;

	// Phase 3: hook up `next` pointers for each face's 3 edges.
	for (int fi = 0; fi < 4; fi++) {
		EpaFace* f = &p->faces[fi];
		p->edges[f->e[0]].next = f->e[1];
		p->edges[f->e[1]].next = f->e[2];
		p->edges[f->e[2]].next = f->e[0];
	}

	return 1;
}

// -----------------------------------------------------------------------------
// EPA expansion loop.
//
// Each iteration: find alive face closest to origin, query support in its
// outward normal, if not further than current best by tolerance we are done,
// else mark visible faces dead and stitch new triangle fan from horizon edges
// to the new vertex. Mirrors qh_horizon_and_cone (see CLAUDE.md).

// Pick the alive face with smallest distance from origin. Returns -1 if none.
static int epa_pick_closest_face(EpaPoly* p)
{
	int best = -1;
	float best_d = 1e30f;
	for (int i = 0; i < p->face_count; i++) {
		if (!p->faces[i].alive) continue;
		if (p->faces[i].dist < best_d) { best_d = p->faces[i].dist; best = i; }
	}
	return best;
}

// Horizon walk + cone stitch: mark all alive faces visible to `new_vert` dead;
// collect the boundary (horizon) half-edges; create a new triangle per
// horizon edge connecting to the new vertex; set up twin links along the
// cone's interior edges. Returns 1 on success.
//
// Implementation: reuses the horizon-edge half-edge slots for the new-face
// boundary edges. Since the horizon's twin is on a live face (unchanged),
// implicit twin(i) = i ^ 1 is preserved automatically for the boundary.
//
// For cone-internal edges (v->new_vert and new_vert->u), we allocate one
// twin pair per horizon segment: the "down" edge for horizon h twins with
// the "up" edge for horizon h+1 (mod H).
static int epa_expand(EpaPoly* p, int new_vert)
{
	v3 pv = p->verts[new_vert].mink;

	// Step 1: flag visible faces (dot(n, pv - v[0]) > 0 means new_vert is outside).
	for (int i = 0; i < p->face_count; i++) {
		if (!p->faces[i].alive) continue;
		EpaFace* f = &p->faces[i];
		if (dot(f->normal, sub(pv, p->verts[f->v[0]].mink)) > 0.0f) {
			f->alive = 0;
		}
	}

	// Step 2: collect horizon half-edges (dead-face side; their twin is on an alive face).
	// The horizon forms a closed loop around the deleted region.
	int horizon[EPA_MAX_EDGES];
	int horizon_count = 0;
	for (int i = 0; i < p->edge_count; i++) {
		int fi = p->edges[i].face;
		if (fi < 0 || fi >= p->face_count) continue;
		if (p->faces[fi].alive) continue;
		int twin = epa_twin(i);
		int tfi = p->edges[twin].face;
		if (tfi < 0 || tfi >= p->face_count) continue;
		if (!p->faces[tfi].alive) continue;
		if (horizon_count >= EPA_MAX_EDGES) return 0;
		horizon[horizon_count++] = i;
	}
	if (horizon_count < 3) return 0;

	// Step 3: sort horizon edges into a cycle (v[i]->v[i+1] boundary loop).
	// Start from horizon[0], then repeatedly find the horizon edge whose v0
	// matches the previous edge's head.
	int ordered[EPA_MAX_EDGES];
	int ordered_count = 0;
	int used[EPA_MAX_EDGES] = { 0 };
	ordered[ordered_count++] = horizon[0];
	used[0] = 1;
	while (ordered_count < horizon_count) {
		int prev = ordered[ordered_count - 1];
		// Head of prev edge = origin of edges[prev.next]. But prev.next points
		// to the next edge on the SAME (dead) face, which may not be another
		// horizon edge. So we compute head directly: the twin's v0 (since twin
		// goes head->tail on live face).
		int head = p->edges[epa_twin(prev)].v0;
		int found = -1;
		for (int h = 0; h < horizon_count; h++) {
			if (used[h]) continue;
			if (p->edges[horizon[h]].v0 == head) { found = h; break; }
		}
		if (found < 0) return 0; // horizon not a simple cycle — bail
		used[found] = 1;
		ordered[ordered_count++] = horizon[found];
	}

	// Step 4: allocate new faces and cone-internal edge pairs.
	int face_idx[EPA_MAX_EDGES];
	int pair_base[EPA_MAX_EDGES];
	for (int h = 0; h < ordered_count; h++) {
		if (p->face_count >= EPA_MAX_FACES) return 0;
		face_idx[h] = p->face_count++;
	}
	for (int h = 0; h < ordered_count; h++) {
		int slot = epa_alloc_edge_pair(p);
		if (slot < 0) return 0;
		pair_base[h] = slot;
	}

	// Step 5: fill faces.
	// For horizon h (u -> v on dead side; twin v -> u on live side):
	//   boundary edge (reuses horizon[h] slot): u -> v, twin = horizon[h]^1 unchanged.
	//   down edge (new alloc, pair_base[h]+0): v -> new_vert.
	//   up edge (new alloc, pair_base[(h-1+H)%H]+1): new_vert -> u, twins with the
	//     down edge of the PREVIOUS horizon (whose v is our u).
	//
	// Twin pairing: down of h is pair_base[h]+0; its sibling up belongs to
	// horizon (h+1)%H, whose u = head_of_h = v. So we place that horizon's
	// up edge at pair_base[h]+1. Thus up[h+1] = pair_base[h]+1.
	for (int h = 0; h < ordered_count; h++) {
		int he = ordered[h];
		int u = p->edges[he].v0;
		// Head: twin's v0 (twin goes head->tail on live side).
		int v = p->edges[epa_twin(he)].v0;

		int e_boundary = he;
		int e_down = pair_base[h] + 0;
		int e_up = pair_base[(h - 1 + ordered_count) % ordered_count] + 1;
		int fi = face_idx[h];

		p->edges[e_boundary].v0 = u;
		p->edges[e_boundary].face = fi;
		p->edges[e_boundary].next = e_down;

		p->edges[e_down].v0 = v;
		p->edges[e_down].face = fi;
		p->edges[e_down].next = e_up;

		p->edges[e_up].v0 = new_vert;
		p->edges[e_up].face = fi;
		p->edges[e_up].next = e_boundary;

		EpaFace* f = &p->faces[fi];
		f->v[0] = u;
		f->v[1] = v;
		f->v[2] = new_vert;
		f->e[0] = e_boundary;
		f->e[1] = e_down;
		f->e[2] = e_up;
		f->alive = 1;

		// Compute plane; outward normal points AWAY from origin (dist >= 0).
		v3 a = p->verts[u].mink, b = p->verts[v].mink, c = p->verts[new_vert].mink;
		v3 n = cross(sub(b, a), sub(c, a));
		float l2 = len2(n);
		if (l2 < EPA_DEGENERATE_EPS) {
			f->alive = 0;
			f->normal = V3(1, 0, 0);
			f->dist = 1e30f;
			continue;
		}
		n = scale(n, 1.0f / sqrtf(l2));
		float d = dot(n, a);
		if (d < 0.0f) { n = neg(n); d = -d; }
		f->normal = n;
		f->dist = d;
	}

	return 1;
}

// -----------------------------------------------------------------------------
// Run EPA. If `warm_dirs`/`warm_valid` are provided, seed the initial tetra
// from those Minkowski-space directions (prior-frame terminating face + a
// centroid direction) rather than the cardinal regular-tetra directions. On
// success and when `out_warm_dirs` is non-NULL, writes the terminating face's
// Minkowski directions back for next-frame reuse.
//
// `used_warm_out` (optional) is set to 1 if the warm seed was actually used.
// `iters_out` (optional) is set to the iteration count actually executed.
// `cap_hit_out` (optional) is set to 1 if the loop ran to EPA_MAX_ITERATIONS.

static int epa_run(EpaSupportPair* sp, const v3* warm_dirs, int warm_valid, EpaHit* out, v3* out_warm_dirs, int* used_warm_out, int* iters_out, int* cap_hit_out)
{
	EpaPoly poly;
	memset(&poly, 0, sizeof(poly));

	// Epsilon from an initial cardinal-direction sample (we don't have a seed
	// simplex any more, but a single support pair gives an adequate extent).
	{
		EpaVert probe;
		epa_support_pair(sp, V3(1, 0, 0), &probe);
		float ax = fabsf(probe.pa.x) + fabsf(probe.pa.y) + fabsf(probe.pa.z);
		float bx = fabsf(probe.pb.x) + fabsf(probe.pb.y) + fabsf(probe.pb.z);
		float mx = ax > bx ? ax : bx;
		poly.epsilon = 3.0f * (mx + 1.0f) * FLT_EPSILON;
		if (poly.epsilon < 1e-7f) poly.epsilon = 1e-7f;
	}

	// First attempt: warm-seeded tetra if available. Fall back to cardinal.
	int used_warm = 0;
	if (!epa_build_tetra(&poly, sp, warm_dirs, warm_valid)) {
		memset(&poly, 0, sizeof(poly));
		poly.epsilon = 1e-7f;
		if (!epa_build_tetra(&poly, sp, NULL, 0)) return 0;
	} else {
		used_warm = warm_valid ? 1 : 0;
	}
	if (used_warm_out) *used_warm_out = used_warm;

	// Main loop.
	int closest = -1;
	float closest_dist = 0.0f;
	v3 closest_normal = V3(0, 1, 0);
	int iter = 0;
	int cap_hit = 0;

	for (iter = 0; iter < EPA_MAX_ITERATIONS; iter++) {
		int fi = epa_pick_closest_face(&poly);
		if (fi < 0) {
			if (iters_out) *iters_out = iter;
			if (cap_hit_out) *cap_hit_out = 0;
			return 0;
		}
		EpaFace* f = &poly.faces[fi];

		// Query support in face normal.
		EpaVert nv;
		epa_support_pair(sp, f->normal, &nv);

		float new_dist = dot(f->normal, nv.mink);

		closest = fi;
		closest_dist = f->dist;
		closest_normal = f->normal;

		if (new_dist - f->dist < EPA_TOLERANCE) break;

		// Weld: don't add a vertex that duplicates an existing one.
		int dup = 0;
		for (int i = 0; i < poly.vert_count; i++) {
			if (len2(sub(nv.mink, poly.verts[i].mink)) <= poly.epsilon * poly.epsilon) { dup = 1; break; }
		}
		if (dup) break;

		if (poly.vert_count >= EPA_MAX_VERTS) break;
		int vi = poly.vert_count++;
		poly.verts[vi] = nv;

		if (!epa_expand(&poly, vi)) break;
	}
	if (iter >= EPA_MAX_ITERATIONS) cap_hit = 1;
	if (iters_out) *iters_out = iter;
	if (cap_hit_out) *cap_hit_out = cap_hit;

	if (closest < 0) return 0;
	EpaFace* f = &poly.faces[closest];

	// Barycentric projection of origin onto face plane to get witness points.
	v3 a = poly.verts[f->v[0]].mink, b = poly.verts[f->v[1]].mink, c = poly.verts[f->v[2]].mink;
	v3 p_on_face = scale(closest_normal, closest_dist); // projection of origin onto face plane

	v3 v0 = sub(b, a), v1 = sub(c, a), v2 = sub(p_on_face, a);
	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);
	float denom = d00 * d11 - d01 * d01;
	float wa, wb, wc;
	if (fabsf(denom) > 1e-20f) {
		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0f - v - w;
		wa = u; wb = v; wc = w;
	} else {
		wa = 1.0f; wb = 0.0f; wc = 0.0f;
	}

	v3 pa = add(add(scale(poly.verts[f->v[0]].pa, wa), scale(poly.verts[f->v[1]].pa, wb)), scale(poly.verts[f->v[2]].pa, wc));
	v3 pb = add(add(scale(poly.verts[f->v[0]].pb, wa), scale(poly.verts[f->v[1]].pb, wb)), scale(poly.verts[f->v[2]].pb, wc));

	out->point_a = pa;
	out->point_b = pb;
	// Normal convention: from A toward B. The MD = B - A contains origin; EPA's
	// outward normal n is the direction to move the MD-origin to the boundary.
	// To separate, we translate B by -n*d (equivalently, the contact pushes A
	// in -n and B in +n). So the contact normal (A->B) is -n.
	out->normal = neg(closest_normal);
	out->depth = closest_dist;

	// Feature id: pack the three vertex feature pairs into a single uint32.
	// Simple scheme: xor of fa/fb values (adequate for matching incremental contacts).
	uint32_t fid = 0;
	for (int i = 0; i < 3; i++) {
		fid ^= (uint32_t)poly.verts[f->v[i]].fa * 2654435761u;
		fid ^= (uint32_t)poly.verts[f->v[i]].fb * 2246822519u;
	}
	out->feature_id = fid;

	// Warm-dir export: three terminating-face Minkowski points + one "interior"
	// direction (centroid as a vector). These seed next frame's support queries.
	if (out_warm_dirs) {
		out_warm_dirs[0] = poly.verts[f->v[0]].mink;
		out_warm_dirs[1] = poly.verts[f->v[1]].mink;
		out_warm_dirs[2] = poly.verts[f->v[2]].mink;
		v3 centroid = scale(add(add(out_warm_dirs[0], out_warm_dirs[1]), out_warm_dirs[2]), 1.0f / 3.0f);
		// Prefer a direction that points INTO the polytope (opposite the terminating face normal)
		// so the four support queries span the MD shape well.
		out_warm_dirs[3] = neg(centroid);
		// Guard against zero-length centroid.
		if (len2(out_warm_dirs[3]) < 1e-12f) out_warm_dirs[3] = neg(closest_normal);
	}

	return 1;
}

// -----------------------------------------------------------------------------
// Shape-pair adapters: wrap GJK_Shape into a support callback for EPA.

typedef struct EpaGjkShapeCtx
{
	GJK_Shape* shape;
} EpaGjkShapeCtx;

static v3 epa_support_gjk_shape(void* ctx_, v3 dir, int* out_feat)
{
	EpaGjkShapeCtx* ctx = (EpaGjkShapeCtx*)ctx_;
	GJK_Shape* shape = ctx->shape;
	int feat = 0;
	v3 p = V3(0, 0, 0);
	gjk_support(shape, dir, &feat, p);
	// Apply radius for sphere/capsule.
	if (shape->radius != 0.0f) {
		float dl = len(dir);
		if (dl > 1e-12f) p = add(p, scale(dir, shape->radius / dl));
	}
	*out_feat = feat;
	return p;
}

// -----------------------------------------------------------------------------
// Public entry points.
//
// Each runs GJK first; if distance > threshold, reports no hit. Otherwise runs
// EPA and returns a single contact. Optional `em` carries per-pair warm-start
// state (warm_dirs/warm_valid) and is updated on success. Optional `stats`
// accumulates telemetry; pass NULL for callers without a world context.

// Shallow-path floor: only skip EPA for GJK distances genuinely above the
// containment noise. Below this we treat the pair as touching and run EPA,
// which returns a real penetration depth even when GJK saw ~1e-7 jitter.
#define EPA_SHALLOW_MIN_DIST 1e-5f

static int epa_hit_from_gjk_shapes_ex(GJK_Shape* ga, GJK_Shape* gb, EpaHit* out, EpaManifold* em, EpaStats* stats)
{
	GJK_Result r = gjk_distance(ga, gb, NULL);
	if (r.distance > EPA_DEEP_THRESHOLD) return 0;
	// Shallow / touching: derive contact directly from GJK witness. Require a
	// meaningful separation so GJK jitter inside a contained configuration
	// falls through to EPA for a real depth.
	if (r.distance > EPA_SHALLOW_MIN_DIST) {
		v3 d = sub(r.point2, r.point1);
		float dist = r.distance;
		if (dist < 1e-8f) dist = 1e-8f;
		out->point_a = r.point1;
		out->point_b = r.point2;
		out->normal = scale(d, 1.0f / dist);
		out->depth = 0.0f;
		out->feature_id = (uint32_t)r.feat1 * 2654435761u ^ (uint32_t)r.feat2 * 2246822519u;
		return 1;
	}

	// Deep: run EPA.
	EpaGjkShapeCtx ca = { .shape = ga };
	EpaGjkShapeCtx cb = { .shape = gb };
	EpaSupportPair sp = { .fa = epa_support_gjk_shape, .ctxa = &ca, .fb = epa_support_gjk_shape, .ctxb = &cb };

	const v3* warm_dirs = (em && em->warm_valid) ? em->warm_dirs : NULL;
	int warm_valid = (em && em->warm_valid) ? 1 : 0;
	v3 new_warm[4];
	int used_warm = 0;
	int iters = 0;
	int cap_hit = 0;
	int ok = epa_run(&sp, warm_dirs, warm_valid, out, em ? new_warm : NULL, &used_warm, &iters, &cap_hit);

	if (stats) {
		stats->queries++;
		stats->total_iters += iters;
		if (cap_hit) stats->iter_cap_hits++;
		if (used_warm) stats->warm_reseeds++;
	}
	if (ok && em) {
		for (int i = 0; i < 4; i++) em->warm_dirs[i] = new_warm[i];
		em->warm_valid = 1;
	}
	return ok;
}

// Cacheless variant used by the direct shape-pair entry points (bench / tests).
static int epa_hit_from_gjk_shapes(GJK_Shape* ga, GJK_Shape* gb, EpaHit* out)
{
	return epa_hit_from_gjk_shapes_ex(ga, gb, out, NULL, NULL);
}

// Hull-hull EPA entry point.
static int epa_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold)
{
	v3 sa[MAX_HULL_VERTS], sb[MAX_HULL_VERTS]; float soa_a[MAX_HULL_VERTS*3], soa_b[MAX_HULL_VERTS*3];
	GJK_Shape ga = gjk_hull_scaled(a.hull, a.center, a.rotation, a.scale, sa, soa_a);
	GJK_Shape gb = gjk_hull_scaled(b.hull, b.center, b.rotation, b.scale, sb, soa_b);
	EpaHit hit;
	if (!epa_hit_from_gjk_shapes(&ga, &gb, &hit)) return 0;
	if (!manifold) return 1;
	v3 contact_point = scale(add(hit.point_a, hit.point_b), 0.5f);
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = contact_point,
		.normal = hit.normal,
		.penetration = hit.depth,
		.feature_id = hit.feature_id,
	};
	return 1;
}

// Sphere-hull EPA entry point.
static int epa_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold)
{
	v3 sb[MAX_HULL_VERTS]; float soa_b[MAX_HULL_VERTS*3];
	GJK_Shape ga = gjk_sphere(a.center, a.radius);
	GJK_Shape gb = gjk_hull_scaled(b.hull, b.center, b.rotation, b.scale, sb, soa_b);
	EpaHit hit;
	if (!epa_hit_from_gjk_shapes(&ga, &gb, &hit)) return 0;
	if (!manifold) return 1;
	v3 contact_point = scale(add(hit.point_a, hit.point_b), 0.5f);
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = contact_point,
		.normal = hit.normal,
		.penetration = hit.depth,
		.feature_id = hit.feature_id,
	};
	return 1;
}

// Capsule-hull EPA entry point.
static int epa_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold)
{
	v3 sb[MAX_HULL_VERTS]; float soa_b[MAX_HULL_VERTS*3];
	GJK_Shape ga = gjk_capsule(a.p, a.q, a.radius);
	GJK_Shape gb = gjk_hull_scaled(b.hull, b.center, b.rotation, b.scale, sb, soa_b);
	EpaHit hit;
	if (!epa_hit_from_gjk_shapes(&ga, &gb, &hit)) return 0;
	if (!manifold) return 1;
	v3 contact_point = scale(add(hit.point_a, hit.point_b), 0.5f);
	manifold->count = 1;
	manifold->contacts[0] = (Contact){
		.point = contact_point,
		.normal = hit.normal,
		.penetration = hit.depth,
		.feature_id = hit.feature_id,
	};
	return 1;
}

// -----------------------------------------------------------------------------
// Incremental manifold merge: accumulate 1-per-frame EPA contacts over time.
//
// Flow per frame:
//   1. Age existing contacts.
//   2. Re-validate: transform local -> world, drop if slipped / separated too far.
//   3. Run EPA once; if hit, merge new contact into cache (replace-by-feature
//      or append, evicting oldest if full).
//   4. Emit accumulated contacts to Manifold* in world space.

// Transform a local-space point to world space using position + rotation.
static inline v3 epa_local_to_world(v3 local, v3 pos, quat rot)
{
	return add(pos, rotate(rot, local));
}

static inline v3 epa_world_to_local(v3 world, v3 pos, quat rot)
{
	return rotate(inv(rot), sub(world, pos));
}

// Sum of squared pairwise distances — spread metric for eviction policy. We
// prefer configurations that maximize this across the candidate 4-point set so
// the manifold stays dispersed across the contact patch.
static float epa_spread_sum_sq(const v3* pts, int n)
{
	float s = 0.0f;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			v3 d = sub(pts[i], pts[j]);
			s += dot(d, d);
		}
	}
	return s;
}

// Merge a new EPA contact into the manifold.
// Rule: if any existing contact has the same feature_id, REPLACE in place and
// reset age. Else APPEND; if full, evict the contact whose removal leaves the
// largest spatial spread among the remaining + new contact (Box2D/Bullet
// pattern). Age-FIFO falls out only on ties.
static void epa_merge_contact(EpaManifold* em, const EpaContact* c_new)
{
	for (int i = 0; i < em->count; i++) {
		if (em->contacts[i].feature_id == c_new->feature_id) {
			em->contacts[i] = *c_new;
			em->contacts[i].age = 0;
			return;
		}
	}
	if (em->count < MAX_CONTACTS) {
		em->contacts[em->count++] = *c_new;
		return;
	}

	// Full manifold: evaluate spread for each eviction candidate. Use
	// point_a_local for the metric — any single body-local frame is fine since
	// we're comparing relative distances.
	int best_evict = 0;
	float best_spread = -1.0f;
	int best_age = -1;
	for (int ei = 0; ei < MAX_CONTACTS; ei++) {
		v3 pts[MAX_CONTACTS];
		int n = 0;
		for (int i = 0; i < MAX_CONTACTS; i++) {
			if (i == ei) continue;
			pts[n++] = em->contacts[i].point_a_local;
		}
		pts[n++] = c_new->point_a_local;
		float s = epa_spread_sum_sq(pts, n);
		int age = em->contacts[ei].age;
		if (s > best_spread + 1e-8f) {
			best_spread = s;
			best_evict = ei;
			best_age = age;
		} else if (fabsf(s - best_spread) <= 1e-8f && age > best_age) {
			// Tie on spread: fall back to age-FIFO (evict oldest).
			best_evict = ei;
			best_age = age;
		}
	}
	em->contacts[best_evict] = *c_new;
}

// Re-validate existing contacts: drop those whose world-space geometry has
// drifted more than tolerance along the contact normal or tangentially across
// the original patch.
static void epa_validate_contacts(EpaManifold* em, v3 pa, quat ra, v3 pb, quat rb)
{
	int write = 0;
	for (int i = 0; i < em->count; i++) {
		EpaContact* c = &em->contacts[i];
		v3 wa = epa_local_to_world(c->point_a_local, pa, ra);
		v3 wb = epa_local_to_world(c->point_b_local, pb, rb);
		v3 wn = rotate(ra, c->normal_local_a);
		v3 d = sub(wb, wa);
		float along = dot(d, wn);
		// Normal-axis separation: drop if geometry moved apart along the contact normal.
		if (along > EPA_CONTACT_REFRESH_DIST) continue;
		// Tangential drift: the original witness points A and B should coincide
		// along the normal; any remaining offset is tangential. If that grew
		// large, the contact patch slid, so the cached pair is stale.
		v3 tang = sub(d, scale(wn, along));
		if (len2(tang) > EPA_CONTACT_TANGENTIAL_DRIFT * EPA_CONTACT_TANGENTIAL_DRIFT) continue;
		// Age-out.
		if (c->age > EPA_MAX_CONTACT_AGE) continue;
		// Update penetration from latest geometry (negated: if along<0, penetration=-along).
		c->penetration = -along;
		em->contacts[write++] = *c;
	}
	em->count = write;
}

// EPA narrowphase entry for a body pair. Handles caching, validation, merging,
// and manifold emission. Returns 1 if any contact emitted.
//
// Callers must ensure the shape types are hull-involved (SHAPE_BOX | SHAPE_HULL
// for A and B, or sphere/capsule + hull).
static int epa_narrowphase_pair(struct WorldInternal* w, int body_a_idx, int body_b_idx, ShapeInternal* s0, ShapeInternal* s1, BodyState* bs0, BodyState* bs1, Manifold* manifold)
{
	(void)body_a_idx; (void)body_b_idx;
	// Build GJK shapes for each body.
	v3 sa_buf[MAX_HULL_VERTS], sb_buf[MAX_HULL_VERTS]; float soa_a[MAX_HULL_VERTS*3], soa_b[MAX_HULL_VERTS*3];
	GJK_Shape ga, gb;
	if (s0->type == SHAPE_SPHERE) {
		v3 center = add(bs0->position, rotate(bs0->rotation, s0->local_pos));
		ga = gjk_sphere(center, s0->sphere.radius);
	} else if (s0->type == SHAPE_CAPSULE) {
		v3 lp = add(s0->local_pos, V3(0, -s0->capsule.half_height, 0));
		v3 lq = add(s0->local_pos, V3(0, s0->capsule.half_height, 0));
		v3 P = add(bs0->position, rotate(bs0->rotation, lp));
		v3 Q = add(bs0->position, rotate(bs0->rotation, lq));
		ga = gjk_capsule(P, Q, s0->capsule.radius);
	} else if (s0->type == SHAPE_BOX) {
		extern const Hull s_unit_box_hull;
		ga = gjk_hull_scaled(&s_unit_box_hull, bs0->position, bs0->rotation, s0->box.half_extents, sa_buf, soa_a);
	} else if (s0->type == SHAPE_HULL) {
		ga = gjk_hull_scaled(s0->hull.hull, bs0->position, bs0->rotation, s0->hull.scale, sa_buf, soa_a);
	} else {
		return 0; // unsupported
	}
	if (s1->type == SHAPE_SPHERE) {
		v3 center = add(bs1->position, rotate(bs1->rotation, s1->local_pos));
		gb = gjk_sphere(center, s1->sphere.radius);
	} else if (s1->type == SHAPE_CAPSULE) {
		v3 lp = add(s1->local_pos, V3(0, -s1->capsule.half_height, 0));
		v3 lq = add(s1->local_pos, V3(0, s1->capsule.half_height, 0));
		v3 P = add(bs1->position, rotate(bs1->rotation, lp));
		v3 Q = add(bs1->position, rotate(bs1->rotation, lq));
		gb = gjk_capsule(P, Q, s1->capsule.radius);
	} else if (s1->type == SHAPE_BOX) {
		extern const Hull s_unit_box_hull;
		gb = gjk_hull_scaled(&s_unit_box_hull, bs1->position, bs1->rotation, s1->box.half_extents, sb_buf, soa_b);
	} else if (s1->type == SHAPE_HULL) {
		gb = gjk_hull_scaled(s1->hull.hull, bs1->position, bs1->rotation, s1->hull.scale, sb_buf, soa_b);
	} else {
		return 0;
	}

	// Look up or create cache entry.
	uint64_t key = body_pair_key(body_a_idx, body_b_idx);
	EpaManifold* em = map_get_ptr(w->epa_cache, key);
	if (!em) {
		EpaManifold zero = { 0 };
		map_set(w->epa_cache, key, zero);
		em = map_get_ptr(w->epa_cache, key);
	}

	// Age existing contacts and re-validate.
	for (int i = 0; i < em->count; i++) em->contacts[i].age++;
	epa_validate_contacts(em, bs0->position, bs0->rotation, bs1->position, bs1->rotation);

	// Dormant-pair detection: if the manifold has gone empty (all contacts
	// aged / drifted out), the Minkowski geometry we had cached is no longer
	// trustworthy as a seed. Drop warm directions so the next EPA call starts
	// cold from cardinal directions.
	if (em->count == 0) em->warm_valid = 0;

	// Run EPA with this pair's warm state. On success, the cache writes back
	// fresh warm_dirs for next frame.
	EpaHit hit;
	int got = epa_hit_from_gjk_shapes_ex(&ga, &gb, &hit, em, &w->epa_stats);
	if (got) {
		EpaContact c_new;
		v3 world_point = scale(add(hit.point_a, hit.point_b), 0.5f);
		c_new.point_a_local = epa_world_to_local(hit.point_a, bs0->position, bs0->rotation);
		c_new.point_b_local = epa_world_to_local(hit.point_b, bs1->position, bs1->rotation);
		c_new.normal_local_a = rotate(inv(bs0->rotation), hit.normal);
		c_new.penetration = hit.depth;
		c_new.feature_id = hit.feature_id;
		c_new.age = 0;
		epa_merge_contact(em, &c_new);
		(void)world_point;
	}

	em->stale = 0;

	// Emit accumulated contacts to manifold (world space).
	if (em->count == 0) return 0;
	w->epa_stats.pair_count++;
	w->epa_stats.contacts_emitted += em->count;
	if (!manifold) return 1;
	manifold->count = em->count;
	for (int i = 0; i < em->count; i++) {
		EpaContact* c = &em->contacts[i];
		v3 wa = epa_local_to_world(c->point_a_local, bs0->position, bs0->rotation);
		v3 wb = epa_local_to_world(c->point_b_local, bs1->position, bs1->rotation);
		v3 wn = rotate(bs0->rotation, c->normal_local_a);
		manifold->contacts[i] = (Contact){
			.point = scale(add(wa, wb), 0.5f),
			.normal = wn,
			.penetration = c->penetration,
			.feature_id = c->feature_id,
		};
	}
	return 1;
}
