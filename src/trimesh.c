// See LICENSE for licensing info.
// trimesh.c -- static triangle-mesh collision geometry with smooth-normal
// (ghost-edge) handling, Box2D-chain-shape style.
//
// Architecture
// ------------
// Build-time:
//   - Copy verts/indices into a single owning block.
//   - Compute per-triangle face normals.
//   - Walk every undirected edge, pair it with at most one neighbor triangle.
//     Manifold input assumed (>= 3 triangles sharing an edge -> assert).
//   - Classify each half-edge as BOUNDARY / CONVEX (ridge) / FLAT / CONCAVE
//     from the dihedral between the two face normals. Pack two bits per edge.
//   - Build an immutable internal BVH over triangle AABBs using the shared
//     binned SAH builder.
//
// Collision time:
//   - Broadphase reports a (mesh body, convex body) pair.
//   - Transform convex into mesh-local space; query mesh BVH.
//   - For each candidate triangle, call collide_<convex>_triangle_ctx, passing
//     a TriNeighbor context with neighbor face normals and edge flags. The
//     per-pair routine applies the same Skip/Admit/Snap rule Box2D uses:
//       * BOUNDARY edge            -> admit contact as-is (no neighbor exists).
//       * FLAT edge                -> admit only if this tri is edge owner
//                                     (lower tri index wins) to avoid dupes.
//       * CONVEX (ridge) edge      -> Voronoi-compare this face normal vs the
//                                     neighbor's; admit if this one "owns" the
//                                     contact-normal direction, else skip.
//       * CONCAVE (valley) edge    -> snap contact normal to this face normal.
//     This preserves the invariant that sliding a convex across a flat or
//     concave seam cannot produce a spurious edge contact pointing outside
//     the triangle face plane.
//
// v1 scope:
//   - static mesh only (mass=0 body); movement is allowed but not dynamic mass.
//   - convex-vs-triangle pairs for sphere first; capsule/box/hull/cylinder
//     land in subsequent phases wired through the same TriNeighbor context.
//   - raycast, smooth-vertex handling beyond the edge case, welding, and
//     streaming are not implemented and not on the v1 list.

// -----------------------------------------------------------------------------
// Edge classification tags (2 bits per half-edge, packed into tri_edge_flag).

enum
{
	TRI_EDGE_BOUNDARY = 0, // no neighbor; emit contact as-is
	TRI_EDGE_CONVEX   = 1, // ridge; Voronoi-compare with neighbor normal
	TRI_EDGE_FLAT     = 2, // coplanar neighbor; emit only from the owner side
	TRI_EDGE_CONCAVE  = 3, // valley; snap contact normal to this face
};

// -----------------------------------------------------------------------------
// Internal TriMesh struct. Single heap block per mesh: the struct header
// followed by all arrays laid out end-to-end.

struct TriMesh
{
	v3 aabb_min;
	v3 aabb_max;

	// Geometry.
	const v3* verts;              // [vert_count]
	const uint32_t* indices;      // [3 * tri_count]
	const v3* tri_normal;         // [tri_count] unit face normals
	const int* tri_edge_adj;      // [3 * tri_count] neighbor tri idx per edge, -1 if boundary
	const uint8_t* tri_edge_flag; // [tri_count] two bits per edge, packed (edge i in bits (2*i)..(2*i+1))

	int vert_count;
	int tri_count;

	// Internal BVH over triangle AABBs.
	BVH_Tree bvh;
};

// -----------------------------------------------------------------------------
// TriNeighbor: context passed into every convex-vs-triangle pair routine so
// that smooth-normal classification happens inside the pair (matching Box2D's
// b2CollideChainSegmentAndPolygon taking b2ChainSegmentParams).

typedef struct TriNeighbor
{
	v3 v0, v1, v2;             // triangle vertices (in the space the pair expects, typically mesh-local)
	v3 normal;                 // face normal (unit)
	int tri_idx;               // own triangle index (for flat-edge ownership tiebreak)
	int adj_idx[3];            // neighbor tri per edge i (edge i = v_i -> v_{i+1})
	v3 adj_normal[3];          // neighbor face normal (valid iff flag != BOUNDARY)
	uint8_t edge_flag[3];      // TRI_EDGE_* per edge
} TriNeighbor;

// -----------------------------------------------------------------------------
// Smooth-normal classifier (3D analog of Box2D's b2ClassifyNormal).
//
// Given a candidate contact whose closest triangle feature is `edge_i` and
// whose world-space contact normal points FROM triangle TOWARD the convex,
// return one of:
//   TRI_CLASSIFY_ADMIT -- use the contact as-is
//   TRI_CLASSIFY_SKIP  -- discard; the neighbor triangle will generate it
//   TRI_CLASSIFY_SNAP  -- force the normal to the face normal before use
//
// Convention: contact_normal is the direction the convex should be pushed out
// along. For a sphere sitting on a triangle this is very close to the face
// normal.

enum
{
	TRI_CLASSIFY_ADMIT = 0,
	TRI_CLASSIFY_SKIP  = 1,
	TRI_CLASSIFY_SNAP  = 2,
};

// Closest-feature classification for a point on a triangle. Returns:
//   -1 = interior
//    0..2 = closest to edge i (vertices v_i and v_{i+1})
// Based on signed barycentrics; the closest edge is the one with the most
// negative barycentric coordinate after clamping into the triangle plane.
// Vertex-dominant contacts are treated as edge contacts on whichever of the
// two incident edges has the more negative bary (conservative: emits smooth
// classification against at least one of the two neighbors).
static int tri_closest_feature(v3 p, v3 v0, v3 v1, v3 v2)
{
	// Signed barycentric of p projected onto the triangle plane.
	v3 e01 = sub(v1, v0);
	v3 e12 = sub(v2, v1);
	v3 e20 = sub(v0, v2);
	v3 n = cross(e01, sub(v2, v0));
	float denom = dot(n, n);
	if (denom <= 1e-20f) return -1; // degenerate; shouldn't happen post-build
	// Bary of p w.r.t. (v0, v1, v2) is (w0, w1, w2) with w_i proportional to
	// signed area of sub-triangle opposite vertex i.
	float w0 = dot(cross(e12, sub(p, v1)), n) / denom;
	float w1 = dot(cross(e20, sub(p, v2)), n) / denom;
	float w2 = 1.0f - w0 - w1;
	// Edge i is opposite vertex (i+2)%3 -- so edge 0 (v0->v1) is opposite v2.
	// The bary most below zero identifies the edge farthest from p: that is
	// the edge opposite, which is the edge the point crossed to leave the
	// triangle. When all w_i > 0 the point is inside; return interior.
	if (w0 >= 0.0f && w1 >= 0.0f && w2 >= 0.0f) return -1;
	// Pick the edge opposite the most-negative bary.
	if (w0 < w1 && w0 < w2) return 1; // opposite v0 = edge v1->v2 = edge 1
	if (w1 < w2) return 2;            // opposite v1 = edge v2->v0 = edge 2
	return 0;                         // opposite v2 = edge v0->v1 = edge 0
}

static int tri_classify_edge(const TriNeighbor* ctx, int edge_i, v3 contact_normal)
{
	uint8_t flag = ctx->edge_flag[edge_i];
	if (flag == TRI_EDGE_BOUNDARY) return TRI_CLASSIFY_ADMIT;
	if (flag == TRI_EDGE_FLAT) {
		// Ownership tiebreak: the lower-index triangle emits.
		return (ctx->tri_idx < ctx->adj_idx[edge_i]) ? TRI_CLASSIFY_ADMIT : TRI_CLASSIFY_SKIP;
	}
	if (flag == TRI_EDGE_CONCAVE) return TRI_CLASSIFY_SNAP;
	// CONVEX: Voronoi slice on the Gauss map. This triangle owns normals that
	// lie closer to its face normal than to the neighbor's face normal.
	// Small slack lets borderline normals pass rather than skip/skip both sides.
	const float sin_tol = 0.01f;
	float d_self  = dot(contact_normal, ctx->normal);
	float d_other = dot(contact_normal, ctx->adj_normal[edge_i]);
	if (d_self >= d_other - sin_tol) return TRI_CLASSIFY_ADMIT;
	return TRI_CLASSIFY_SKIP;
}

// -----------------------------------------------------------------------------
// trimesh_create -- build pipeline.
//
// Layout of the single heap block:
//   [ TriMesh header | verts[vc] | indices[3*tc] | tri_normal[tc]
//     | tri_edge_adj[3*tc] | tri_edge_flag[tc] ]
// Followed by a separately-allocated BVH_Tree (its arrays are ckit dynamic
// arrays, owned but not pulled into the header block).

static inline uint8_t trimesh_pack_flag(uint8_t tag, int edge_i)
{
	return (uint8_t)((tag & 3u) << (edge_i * 2));
}

static inline uint8_t trimesh_unpack_flag(uint8_t packed, int edge_i)
{
	return (uint8_t)((packed >> (edge_i * 2)) & 3u);
}

// Undirected edge key from two vertex ids.
static inline uint64_t trimesh_edge_key(uint32_t a, uint32_t b)
{
	uint32_t lo = a < b ? a : b;
	uint32_t hi = a < b ? b : a;
	return ((uint64_t)hi << 32) | (uint64_t)lo;
}

// Slot recorded per undirected edge during adjacency building.
typedef struct TriMeshEdgeSlot
{
	int tri_a;    // first triangle sharing this edge
	int edge_a;   // which half-edge (0..2) on tri_a
	int tri_b;    // second tri, or -1 if boundary so far
	int edge_b;
} TriMeshEdgeSlot;

TriMesh* trimesh_create(const v3* verts, int vert_count, const uint32_t* indices, int tri_count)
{
	assert(vert_count > 0 && tri_count > 0);
	assert(verts && indices);

	size_t verts_bytes    = (size_t)vert_count * sizeof(v3);
	size_t indices_bytes  = (size_t)tri_count * 3 * sizeof(uint32_t);
	size_t normals_bytes  = (size_t)tri_count * sizeof(v3);
	size_t edge_adj_bytes = (size_t)tri_count * 3 * sizeof(int);
	size_t edge_flag_bytes = (size_t)tri_count * sizeof(uint8_t);

	// 16-byte alignment for v3 arrays.
	size_t header = (sizeof(TriMesh) + 15) & ~(size_t)15;
	size_t off_verts    = header;
	size_t off_indices  = off_verts + ((verts_bytes + 15) & ~(size_t)15);
	size_t off_normals  = off_indices + ((indices_bytes + 15) & ~(size_t)15);
	size_t off_adj      = off_normals + ((normals_bytes + 15) & ~(size_t)15);
	size_t off_flag     = off_adj + ((edge_adj_bytes + 15) & ~(size_t)15);
	size_t total        = off_flag + ((edge_flag_bytes + 15) & ~(size_t)15);

	char* block = (char*)calloc(1, total);
	TriMesh* m = (TriMesh*)block;

	v3*       mverts     = (v3*)(block + off_verts);
	uint32_t* mindices   = (uint32_t*)(block + off_indices);
	v3*       mnormals   = (v3*)(block + off_normals);
	int*      medge_adj  = (int*)(block + off_adj);
	uint8_t*  medge_flag = (uint8_t*)(block + off_flag);

	m->verts          = mverts;
	m->indices        = mindices;
	m->tri_normal     = mnormals;
	m->tri_edge_adj   = medge_adj;
	m->tri_edge_flag  = medge_flag;
	m->vert_count     = vert_count;
	m->tri_count      = tri_count;

	memcpy(mverts, verts, verts_bytes);
	memcpy(mindices, indices, indices_bytes);

	// Compute per-triangle face normals and global AABB.
	v3 amin = V3(1e18f, 1e18f, 1e18f), amax = V3(-1e18f, -1e18f, -1e18f);
	for (int t = 0; t < tri_count; t++) {
		uint32_t i0 = mindices[3*t + 0], i1 = mindices[3*t + 1], i2 = mindices[3*t + 2];
		assert((int)i0 < vert_count && (int)i1 < vert_count && (int)i2 < vert_count);
		v3 a = mverts[i0], b = mverts[i1], c = mverts[i2];
		v3 n = cross(sub(b, a), sub(c, a));
		float len2 = dot(n, n);
		assert(len2 > 1e-20f); // degenerate triangle
		float inv_len = 1.0f / sqrtf(len2);
		mnormals[t] = scale(n, inv_len);
		amin = v3_min(amin, v3_min(a, v3_min(b, c)));
		amax = v3_max(amax, v3_max(a, v3_max(b, c)));
	}
	m->aabb_min = amin;
	m->aabb_max = amax;

	// Init adjacency to boundary.
	for (int i = 0; i < 3 * tri_count; i++) medge_adj[i] = -1;
	for (int i = 0; i < tri_count; i++) medge_flag[i] = trimesh_pack_flag(TRI_EDGE_BOUNDARY, 0) | trimesh_pack_flag(TRI_EDGE_BOUNDARY, 1) | trimesh_pack_flag(TRI_EDGE_BOUNDARY, 2);

	// Build edge adjacency: for each half-edge, record (tri_a, edge_a). When
	// the second triangle for an undirected edge arrives, wire both directions.
	CK_MAP(TriMeshEdgeSlot) edge_map = NULL;
	for (int t = 0; t < tri_count; t++) {
		for (int e = 0; e < 3; e++) {
			uint32_t va = mindices[3*t + e];
			uint32_t vb = mindices[3*t + (e + 1) % 3];
			uint64_t key = trimesh_edge_key(va, vb);
			TriMeshEdgeSlot* slot = map_get_ptr(edge_map, key);
			if (!slot) {
				TriMeshEdgeSlot fresh = { .tri_a = t, .edge_a = e, .tri_b = -1, .edge_b = -1 };
				map_set(edge_map, key, fresh);
			} else {
				assert(slot->tri_b == -1 && "non-manifold mesh: edge shared by >2 triangles");
				slot->tri_b = t;
				slot->edge_b = e;
			}
		}
	}

	// Resolve flags from paired slots.
	const float flat_cos = 0.9998f;   // dihedral < ~1.1 degrees from coplanar
	const float convex_sin_tol = 0.01f;
	map_each(edge_map, mi) {
		TriMeshEdgeSlot* slot = &map_val(edge_map, mi);
		if (slot->tri_b == -1) continue; // boundary stays as default
		int ta = slot->tri_a, ea = slot->edge_a;
		int tb = slot->tri_b, eb = slot->edge_b;
		v3 na = mnormals[ta];
		v3 nb = mnormals[tb];

		medge_adj[3*ta + ea] = tb;
		medge_adj[3*tb + eb] = ta;

		uint8_t tag;
		float d = dot(na, nb);
		if (d >= flat_cos) {
			tag = TRI_EDGE_FLAT;
		} else {
			// Convex vs concave: sign of cross(na, nb) . edge_dir tells us
			// whether the neighbor folds away from or toward the face.
			uint32_t via = mindices[3*ta + ea];
			uint32_t vib = mindices[3*ta + (ea + 1) % 3];
			v3 edge_dir = sub(mverts[vib], mverts[via]);
			v3 fold = cross(na, nb);
			float sign = dot(fold, edge_dir);
			// Convention: CCW winding means a convex ridge has fold aligned with
			// edge_dir (right-hand rule); concave valley has opposite sign.
			if (sign > convex_sin_tol) tag = TRI_EDGE_CONVEX;
			else if (sign < -convex_sin_tol) tag = TRI_EDGE_CONCAVE;
			else tag = TRI_EDGE_FLAT; // near-coplanar fallback
		}

		// Pack tag into each tri's flag byte for its respective edge slot.
		uint8_t a_flag = medge_flag[ta];
		a_flag = (uint8_t)((a_flag & ~(3u << (ea * 2))) | trimesh_pack_flag(tag, ea));
		medge_flag[ta] = a_flag;

		uint8_t b_flag = medge_flag[tb];
		b_flag = (uint8_t)((b_flag & ~(3u << (eb * 2))) | trimesh_pack_flag(tag, eb));
		medge_flag[tb] = b_flag;
	}
	map_free(edge_map);

	// Build the internal BVH over triangle AABBs. bvh_binned_build expects
	// pre-allocated leaf slots; each leaf's body_idx holds the triangle index.
	bvh_init(&m->bvh);
	AABB* tri_aabb = (AABB*)malloc(sizeof(AABB) * tri_count);
	CK_DYNA int* lis = NULL;
	for (int t = 0; t < tri_count; t++) {
		v3 a = mverts[mindices[3*t + 0]];
		v3 b = mverts[mindices[3*t + 1]];
		v3 c = mverts[mindices[3*t + 2]];
		int li = bvh_alloc_leaf(&m->bvh);
		m->bvh.leaves[li].body_idx = t;
		tri_aabb[li] = (AABB){ v3_min(a, v3_min(b, c)), v3_max(a, v3_max(b, c)) };
		apush(lis, li);
	}
	if (tri_count == 1) {
		// Single-leaf degenerate tree: root is a node whose child A is the leaf.
		int ni = bvh_alloc_node(&m->bvh);
		bvh_child_set_leaf(&m->bvh.nodes[ni].a, tri_aabb[lis[0]], lis[0]);
		m->bvh.leaves[lis[0]].node_idx = ni;
		m->bvh.leaves[lis[0]].child_slot = 0;
		m->bvh.root = ni;
	} else {
		m->bvh.root = bvh_binned_build(&m->bvh, lis, tri_aabb, tri_count);
		m->bvh.meta[m->bvh.root].parent = -1;
	}
	afree(lis);
	free(tri_aabb);

	return m;
}

void trimesh_free(TriMesh* m)
{
	if (!m) return;
	bvh_free(&m->bvh);
	free(m);
}

int trimesh_tri_count(const TriMesh* m)
{
	return m ? m->tri_count : 0;
}

// Precomputed mesh-local AABB. Used by bvh.c's shape_aabb to transform the
// mesh body's bounds for broadphase.
static AABB trimesh_local_aabb(const TriMesh* m)
{
	return (AABB){ m->aabb_min, m->aabb_max };
}

// -----------------------------------------------------------------------------
// Raycast against mesh. Uses the internal BVH to prune candidate triangles,
// Moeller-Trumbore for ray-triangle intersection. Returns nearest hit.

static int ray_triangle(v3 ro, v3 rd, v3 v0, v3 v1, v3 v2, float max_t, float* t_out)
{
	v3 e1 = sub(v1, v0);
	v3 e2 = sub(v2, v0);
	v3 h = cross(rd, e2);
	float a = dot(e1, h);
	if (fabsf(a) < 1e-12f) return 0; // ray parallel to triangle
	float f = 1.0f / a;
	v3 s = sub(ro, v0);
	float u = f * dot(s, h);
	if (u < 0.0f || u > 1.0f) return 0;
	v3 qv = cross(s, e1);
	float v = f * dot(rd, qv);
	if (v < 0.0f || u + v > 1.0f) return 0;
	float t = f * dot(e2, qv);
	if (t < 0.0f || t > max_t) return 0;
	*t_out = t;
	return 1;
}

// BVH ray traversal: populate candidate triangle indices whose AABB the ray pierces.
static void trimesh_query_ray_node(const BVH_Tree* t, int ni, v3 ro, v3 inv_dir, float max_t, CK_DYNA int** out)
{
	const BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		const BVH_Child* c = s == 0 ? &n->a : &n->b;
		if (c->leaf_count == 0) continue;
		// Slab test.
		v3 tmin = V3((c->min.x - ro.x) * inv_dir.x, (c->min.y - ro.y) * inv_dir.y, (c->min.z - ro.z) * inv_dir.z);
		v3 tmax = V3((c->max.x - ro.x) * inv_dir.x, (c->max.y - ro.y) * inv_dir.y, (c->max.z - ro.z) * inv_dir.z);
		float t0x = tmin.x < tmax.x ? tmin.x : tmax.x, t1x = tmin.x < tmax.x ? tmax.x : tmin.x;
		float t0y = tmin.y < tmax.y ? tmin.y : tmax.y, t1y = tmin.y < tmax.y ? tmax.y : tmin.y;
		float t0z = tmin.z < tmax.z ? tmin.z : tmax.z, t1z = tmin.z < tmax.z ? tmax.z : tmin.z;
		float tmn = fmaxf(fmaxf(t0x, t0y), t0z);
		float tmx = fminf(fminf(t1x, t1y), t1z);
		if (tmx < 0.0f || tmn > tmx || tmn > max_t) continue;
		if (c->index < 0) apush(*out, t->leaves[~c->index].body_idx);
		else trimesh_query_ray_node(t, c->index, ro, inv_dir, max_t, out);
	}
}

// World-space ray vs mesh. Mesh body pose (pos/rot) applied; ray transformed
// into mesh-local frame for traversal. Returns 1 on hit with t and world-space
// normal (from mesh toward exterior).
static int ray_mesh(v3 ro, v3 rd, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, float max_t, float* t_out, v3* n_out)
{
	quat inv_rot = inv(mesh_rot);
	v3 lro = rotate(inv_rot, sub(ro, mesh_pos));
	v3 lrd = rotate(inv_rot, rd);
	v3 inv_dir = V3(fabsf(lrd.x) > 1e-20f ? 1.0f / lrd.x : 1e20f,
	                 fabsf(lrd.y) > 1e-20f ? 1.0f / lrd.y : 1e20f,
	                 fabsf(lrd.z) > 1e-20f ? 1.0f / lrd.z : 1e20f);

	CK_DYNA int* cands = NULL;
	if (mesh->bvh.root >= 0) trimesh_query_ray_node(&mesh->bvh, mesh->bvh.root, lro, inv_dir, max_t, &cands);

	float best_t = max_t;
	int best_tri = -1;
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		uint32_t i0 = mesh->indices[3*t + 0], i1 = mesh->indices[3*t + 1], i2 = mesh->indices[3*t + 2];
		float tt;
		if (ray_triangle(lro, lrd, mesh->verts[i0], mesh->verts[i1], mesh->verts[i2], best_t, &tt) && tt < best_t) {
			best_t = tt; best_tri = t;
		}
	}
	afree(cands);

	if (best_tri < 0) return 0;
	*t_out = best_t;
	*n_out = rotate(mesh_rot, mesh->tri_normal[best_tri]);
	return 1;
}

// -----------------------------------------------------------------------------
// Shape-local helpers: fetch triangle vertices and pack a TriNeighbor.

static inline void trimesh_load_ctx(const TriMesh* m, int t, TriNeighbor* ctx)
{
	uint32_t i0 = m->indices[3*t + 0], i1 = m->indices[3*t + 1], i2 = m->indices[3*t + 2];
	ctx->v0 = m->verts[i0];
	ctx->v1 = m->verts[i1];
	ctx->v2 = m->verts[i2];
	ctx->normal = m->tri_normal[t];
	ctx->tri_idx = t;
	for (int e = 0; e < 3; e++) {
		int adj = m->tri_edge_adj[3*t + e];
		ctx->adj_idx[e] = adj;
		ctx->adj_normal[e] = adj >= 0 ? m->tri_normal[adj] : V3(0,0,0);
		ctx->edge_flag[e] = trimesh_unpack_flag(m->tri_edge_flag[t], e);
	}
}

// -----------------------------------------------------------------------------
// Sphere vs triangle (mesh-local frame).
// Finds the closest point on the triangle to the sphere center, tests overlap,
// classifies the contact against smooth-normal rules. Returns 1 on accepted
// contact and writes {point, normal, penetration} into *out.

// Closest point on a triangle to point p. Returns the closest point and
// optionally the triangle feature (-1=interior, 0..2=edge i, 10..12=vertex (i-10)).
static v3 closest_point_triangle(v3 p, v3 a, v3 b, v3 c, int* out_feature)
{
	v3 ab = sub(b, a);
	v3 ac = sub(c, a);
	v3 ap = sub(p, a);
	float d1 = dot(ab, ap);
	float d2 = dot(ac, ap);
	if (d1 <= 0.0f && d2 <= 0.0f) { if (out_feature) *out_feature = 10; return a; }

	v3 bp = sub(p, b);
	float d3 = dot(ab, bp);
	float d4 = dot(ac, bp);
	if (d3 >= 0.0f && d4 <= d3) { if (out_feature) *out_feature = 11; return b; }

	float vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
		float t = d1 / (d1 - d3);
		if (out_feature) *out_feature = 0; // edge v0-v1
		return add(a, scale(ab, t));
	}

	v3 cp = sub(p, c);
	float d5 = dot(ab, cp);
	float d6 = dot(ac, cp);
	if (d6 >= 0.0f && d5 <= d6) { if (out_feature) *out_feature = 12; return c; }

	float vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
		float t = d2 / (d2 - d6);
		if (out_feature) *out_feature = 2; // edge v2-v0
		return add(a, scale(ac, t));
	}

	float va = d3 * d6 - d5 * d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
		float t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		if (out_feature) *out_feature = 1; // edge v1-v2
		return add(b, scale(sub(c, b), t));
	}

	float denom = 1.0f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	if (out_feature) *out_feature = -1;
	return add(a, add(scale(ab, v), scale(ac, w)));
}

// Returns 1 on emitted contact. `center` and out->point are in the frame the
// triangle lives in (mesh-local for the mesh pair entry point).
static int collide_sphere_triangle_ctx(v3 center, float radius, const TriNeighbor* ctx, Contact* out)
{
	int feature;
	v3 cp = closest_point_triangle(center, ctx->v0, ctx->v1, ctx->v2, &feature);
	v3 d = sub(center, cp);
	float d2 = dot(d, d);
	float r = radius;
	if (d2 > r * r) return 0;

	v3 normal;
	float pen;
	if (d2 > 1e-12f) {
		float len = sqrtf(d2);
		normal = scale(d, 1.0f / len);
		pen = r - len;
	} else {
		// Sphere center coincides with triangle surface: fall back to face normal.
		normal = ctx->normal;
		pen = r;
	}

	// Smooth-normal classification for edge-dominant contacts.
	if (feature >= 0 && feature <= 2) {
		int verdict = tri_classify_edge(ctx, feature, normal);
		if (verdict == TRI_CLASSIFY_SKIP) return 0;
		if (verdict == TRI_CLASSIFY_SNAP) {
			normal = ctx->normal;
			pen = r - dot(sub(center, cp), ctx->normal);
			if (pen <= 0.0f) return 0;
		}
	} else if (feature >= 10 && feature <= 12) {
		// Vertex-dominant: classify against both incident edges. Skip only if
		// BOTH neighbors skip (conservative admit). Snap if either snaps.
		int v = feature - 10;
		int e_prev = (v + 2) % 3; // edge leading into vertex v
		int e_next = v;           // edge leaving vertex v
		int a = tri_classify_edge(ctx, e_prev, normal);
		int b = tri_classify_edge(ctx, e_next, normal);
		if (a == TRI_CLASSIFY_SKIP && b == TRI_CLASSIFY_SKIP) return 0;
		if (a == TRI_CLASSIFY_SNAP || b == TRI_CLASSIFY_SNAP) {
			normal = ctx->normal;
			pen = r - dot(sub(center, cp), ctx->normal);
			if (pen <= 0.0f) return 0;
		}
	}

	out->point = cp;
	out->normal = normal;
	out->penetration = pen;
	// Feature id: pack triangle index + feature slot for warm matching.
	out->feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)(feature & 0xFF) << 24);
	return 1;
}

// -----------------------------------------------------------------------------
// Mesh-level BVH traversal: collect candidate triangle indices whose AABB
// overlaps the query AABB (expressed in mesh-local space).

static void trimesh_query_aabb_node(const BVH_Tree* t, int ni, AABB q, CK_DYNA int** out)
{
	const BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		const BVH_Child* c = s == 0 ? &n->a : &n->b;
		if (c->leaf_count == 0) continue;
		AABB cb = { c->min, c->max };
		if (!aabb_overlaps(cb, q)) continue;
		if (c->index < 0) apush(*out, t->leaves[~c->index].body_idx);
		else trimesh_query_aabb_node(t, c->index, q, out);
	}
}

static void trimesh_query_aabb(const TriMesh* m, AABB q, CK_DYNA int** out)
{
	if (m->bvh.root < 0) return;
	trimesh_query_aabb_node(&m->bvh, m->bvh.root, q, out);
}

// -----------------------------------------------------------------------------
// Capsule vs triangle (mesh-local frame).
//
// Strategy: compute closest point pair between the capsule segment and the
// triangle. Three candidate sources:
//   - each capsule endpoint vs triangle (point-vs-triangle closest)
//   - capsule segment vs each triangle edge (segment-vs-segment closest)
//   - segment pierces triangle plane inside the triangle (distance zero)
// Pick the minimum. Classify the closest triangle feature and apply the
// Skip/Admit/Snap smooth-normal rule. Emit a single contact per triangle;
// rotational stability across a flat area comes from the capsule touching
// multiple triangles (one manifold each).

// Segment-triangle closest pair. Returns squared distance. Writes:
//   *cs = closest point on segment,
//   *ct = closest point on triangle,
//   *feat = -1 interior, 0..2 edge i, 10..12 vertex (i-10).
static float segment_triangle_closest(v3 p, v3 q, v3 v0, v3 v1, v3 v2, v3 tri_n, v3* cs, v3* ct, int* feat)
{
	float best_d2 = 1e30f;
	int best_feat = -1;
	v3 best_cs = p, best_ct = v0;

	// (a) Segment pierces triangle plane inside it -> distance 0.
	v3 pq = sub(q, p);
	float denom = dot(tri_n, pq);
	if (fabsf(denom) > 1e-12f) {
		float t = (dot(tri_n, v0) - dot(tri_n, p)) / denom;
		if (t >= 0.0f && t <= 1.0f) {
			v3 pierce = add(p, scale(pq, t));
			// Inside-triangle test via signed barycentrics.
			v3 e01 = sub(v1, v0), e12 = sub(v2, v1), e20 = sub(v0, v2);
			int ok = dot(cross(e01, sub(pierce, v0)), tri_n) >= 0.0f
			      && dot(cross(e12, sub(pierce, v1)), tri_n) >= 0.0f
			      && dot(cross(e20, sub(pierce, v2)), tri_n) >= 0.0f;
			if (ok) {
				*cs = pierce; *ct = pierce; *feat = -1;
				return 0.0f;
			}
		}
	}

	// (b) Each segment endpoint vs triangle (point-vs-triangle).
	for (int e = 0; e < 2; e++) {
		v3 pt = e == 0 ? p : q;
		int f;
		v3 cp = closest_point_triangle(pt, v0, v1, v2, &f);
		v3 d = sub(cp, pt);
		float d2 = dot(d, d);
		if (d2 < best_d2) { best_d2 = d2; best_cs = pt; best_ct = cp; best_feat = f; }
	}

	// (c) Segment vs each triangle edge.
	v3 ev[3][2] = { { v0, v1 }, { v1, v2 }, { v2, v0 } };
	for (int i = 0; i < 3; i++) {
		v3 ca, cb;
		segments_closest_points(p, q, ev[i][0], ev[i][1], &ca, &cb);
		v3 d = sub(cb, ca);
		float d2 = dot(d, d);
		if (d2 < best_d2) { best_d2 = d2; best_cs = ca; best_ct = cb; best_feat = i; }
	}

	*cs = best_cs; *ct = best_ct; *feat = best_feat;
	return best_d2;
}

// Contact: normal from triangle toward capsule (caller flips to match A=capsule, B=mesh).
static int collide_capsule_triangle_ctx(v3 p, v3 q, float radius, const TriNeighbor* ctx, Contact* out)
{
	v3 cs, ct;
	int feat;
	float d2 = segment_triangle_closest(p, q, ctx->v0, ctx->v1, ctx->v2, ctx->normal, &cs, &ct, &feat);
	if (d2 > radius * radius) return 0;

	v3 d = sub(cs, ct);
	v3 normal;
	float pen;
	if (d2 > 1e-12f) {
		float len = sqrtf(d2);
		normal = scale(d, 1.0f / len);
		pen = radius - len;
	} else {
		// Segment touches or pierces the triangle: fall back to face normal.
		normal = ctx->normal;
		pen = radius;
	}

	if (feat >= 0 && feat <= 2) {
		int verdict = tri_classify_edge(ctx, feat, normal);
		if (verdict == TRI_CLASSIFY_SKIP) return 0;
		if (verdict == TRI_CLASSIFY_SNAP) {
			normal = ctx->normal;
			pen = radius - dot(sub(cs, ct), ctx->normal);
			if (pen <= 0.0f) return 0;
		}
	} else if (feat >= 10 && feat <= 12) {
		int v = feat - 10;
		int e_prev = (v + 2) % 3;
		int e_next = v;
		int a = tri_classify_edge(ctx, e_prev, normal);
		int b = tri_classify_edge(ctx, e_next, normal);
		if (a == TRI_CLASSIFY_SKIP && b == TRI_CLASSIFY_SKIP) return 0;
		if (a == TRI_CLASSIFY_SNAP || b == TRI_CLASSIFY_SNAP) {
			normal = ctx->normal;
			pen = radius - dot(sub(cs, ct), ctx->normal);
			if (pen <= 0.0f) return 0;
		}
	}

	out->point = ct;
	out->normal = normal;
	out->penetration = pen;
	out->feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)(feat & 0xFF) << 24);
	return 1;
}

// -----------------------------------------------------------------------------
// Multi-contact box vs triangle. Iterates the box's 8 corners; any corner at
// or below the triangle plane whose projection falls inside the triangle
// becomes a contact. Up to 4 contacts per triangle (MAX_CONTACTS).
//
// This gives rotational stability for a box resting flat on a triangle face
// -- single-contact-per-triangle fails because a single point can't resist
// tipping. For corners that fall outside the triangle (overhanging the
// edge), the neighboring triangle's manifold handles them via the smooth-
// normal rules. No edge/vertex smooth-normal classification needed here
// since we only emit contacts for interior projections.
static int collide_box_triangle_multi(Box box, const TriNeighbor* ctx, Manifold* out)
{
	int count = 0;
	v3 he = box.half_extents;
	for (int i = 0; i < 8 && count < MAX_CONTACTS; i++) {
		v3 local = V3((i & 1) ? he.x : -he.x, (i & 2) ? he.y : -he.y, (i & 4) ? he.z : -he.z);
		v3 world = add(box.center, rotate(box.rotation, local));
		float signed_dist = dot(ctx->normal, sub(world, ctx->v0));
		if (signed_dist > LINEAR_SLOP) continue;
		v3 projected = sub(world, scale(ctx->normal, signed_dist));
		int feat;
		closest_point_triangle(projected, ctx->v0, ctx->v1, ctx->v2, &feat);
		if (feat != -1) continue; // projection outside triangle; neighbor owns it
		out->contacts[count].point = projected;
		out->contacts[count].normal = ctx->normal;
		out->contacts[count].penetration = -signed_dist;
		out->contacts[count].feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)i << 24);
		count++;
	}
	if (count == 0) return 0;
	out->count = count;
	return 1;
}

// -----------------------------------------------------------------------------
// SIMD-batched sample-points-vs-4-triangles. The sample-point approach
// unifies box corners, hull vertices, cylinder cap rim points, and capsule
// endpoints: each shape pre-computes a list of surface-or-near-surface
// reference points, and this routine projects them onto 4 triangles in
// parallel.
//
// Semantics (radius = 0 for box/hull/cylinder; radius = r for capsule):
//   For each sample point P and triangle (v0, v1, v2, n):
//     d_p     = dot(n, P - v0)                  // P's signed distance to plane
//     skip if d_p > radius + LINEAR_SLOP        // out of range
//     proj    = P - n * d_p                     // P projected onto tri plane
//     inside  = all three bary-signs >= 0       // proj within triangle
//     emit contact at proj with pen = radius - d_p, normal = tri_normal
//
// Each shape's sample set captures the geometric "corners" that matter for
// rotational stability:
//   box      -- 8 vertices
//   hull     -- up to vert_count vertices
//   cylinder -- 4 rim points per cap (8 total)
//   capsule  -- 2 endpoint points (radius = capsule.radius)
//
// Returns a 4-bit mask: bit i set iff lane i emitted at least one contact.
static int collide_samples_tris_batch4(const v3* samples, int n_samples, float radius, const TriNeighbor ctx[4], Manifold out_mfs[4])
{
	v3w v0 = v3w_load4(ctx[0].v0, ctx[1].v0, ctx[2].v0, ctx[3].v0);
	v3w v1 = v3w_load4(ctx[0].v1, ctx[1].v1, ctx[2].v1, ctx[3].v1);
	v3w v2 = v3w_load4(ctx[0].v2, ctx[1].v2, ctx[2].v2, ctx[3].v2);
	v3w nml = v3w_load4(ctx[0].normal, ctx[1].normal, ctx[2].normal, ctx[3].normal);
	v3w e01 = v3w_sub(v1, v0);
	v3w e12 = v3w_sub(v2, v1);
	v3w e20 = v3w_sub(v0, v2);

	simd4f range = simd_set1(radius + LINEAR_SLOP);
	simd4f zero = simd_set1(0.0f);

	int any = 0;
	out_mfs[0].count = out_mfs[1].count = out_mfs[2].count = out_mfs[3].count = 0;

	for (int s = 0; s < n_samples; s++) {
		v3w cw = v3w_set1(samples[s]);
		v3w cv0 = v3w_sub(cw, v0);
		simd4f d_p = v3w_dot(nml, cv0);
		simd4f within = simd_cmple(d_p, range);
		v3w proj = v3w_sub(cw, v3w_scale(nml, d_p));
		simd4f s0 = v3w_dot(v3w_cross(e01, v3w_sub(proj, v0)), nml);
		simd4f s1 = v3w_dot(v3w_cross(e12, v3w_sub(proj, v1)), nml);
		simd4f s2 = v3w_dot(v3w_cross(e20, v3w_sub(proj, v2)), nml);
		simd4f inside = simd_and(simd_and(simd_cmpge(s0, zero), simd_cmpge(s1, zero)), simd_cmpge(s2, zero));
		simd4f emit = simd_and(within, inside);
		int mask = simd_movemask(emit);
		if (!mask) continue;

		float dp[4]; simd_store(dp, d_p);
		float px[4]; simd_store(px, proj.x);
		float py[4]; simd_store(py, proj.y);
		float pz[4]; simd_store(pz, proj.z);

		for (int lane = 0; lane < 4; lane++) {
			if (!(mask & (1 << lane))) continue;
			if (out_mfs[lane].count >= MAX_CONTACTS) continue;
			int ci = out_mfs[lane].count++;
			out_mfs[lane].contacts[ci].point = V3(px[lane], py[lane], pz[lane]);
			out_mfs[lane].contacts[ci].normal = ctx[lane].normal;
			out_mfs[lane].contacts[ci].penetration = radius - dp[lane];
			out_mfs[lane].contacts[ci].feature_id = (uint32_t)(ctx[lane].tri_idx & 0x00FFFFFFu) | ((uint32_t)s << 24);
			any |= (1 << lane);
		}
	}
	return any;
}

// Multi-contact capsule vs triangle. Samples both segment endpoints (the
// capsule's "deepest surface points" along -tri_normal are endpoint - n*radius).
// For a horizontal capsule flat on a triangle, both endpoints emit contacts
// so the capsule doesn't jitter about its long axis. If neither endpoint
// projects inside the triangle, falls back to the single-contact closest-
// point path (capsule poking mid-segment into an edge/corner).
static int collide_capsule_triangle_multi(v3 p, v3 q, float radius, const TriNeighbor* ctx, Manifold* out)
{
	int count = 0;
	v3 endpoints[2] = { p, q };
	for (int i = 0; i < 2 && count < MAX_CONTACTS; i++) {
		v3 pt = endpoints[i];
		float proj_dist = dot(ctx->normal, sub(pt, ctx->v0));
		float surface_dist = proj_dist - radius;
		if (surface_dist > LINEAR_SLOP) continue;
		// Surface point is endpoint offset by radius toward the plane.
		v3 surface_pt = sub(pt, scale(ctx->normal, radius));
		v3 projected = sub(surface_pt, scale(ctx->normal, surface_dist));
		int feat;
		closest_point_triangle(projected, ctx->v0, ctx->v1, ctx->v2, &feat);
		if (feat != -1) continue;
		out->contacts[count].point = projected;
		out->contacts[count].normal = ctx->normal;
		out->contacts[count].penetration = -surface_dist;
		out->contacts[count].feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)i << 24);
		count++;
	}
	if (count == 0) {
		// Fallback to single-contact (angled capsule where neither endpoint
		// projects inside triangle but mid-segment contacts an edge/vertex).
		Contact c;
		if (collide_capsule_triangle_ctx(p, q, radius, ctx, &c)) {
			out->contacts[0] = c;
			out->count = 1;
			return 1;
		}
		return 0;
	}
	out->count = count;
	return 1;
}

// Multi-contact cylinder vs triangle. Samples 4 rim points on each cap
// (perpendicular to cylinder axis in local X and Z), 8 samples total. Gives
// proper rotational stability for:
//   - cylinder on side: bottom-rim points at each end contact the plane
//   - cylinder on end: cap rim in the -tri_normal direction contacts the plane
// Falls back to GJK closest-point for edge-dominant contacts.
static int collide_cylinder_triangle_multi(Cylinder cyl, const TriNeighbor* ctx, Manifold* out)
{
	int count = 0;
	v3 axis = rotate(cyl.rotation, V3(0, 1, 0));
	v3 perp1 = rotate(cyl.rotation, V3(1, 0, 0));
	v3 perp2 = rotate(cyl.rotation, V3(0, 0, 1));
	v3 c0 = sub(cyl.center, scale(axis, cyl.half_height));
	v3 c1 = add(cyl.center, scale(axis, cyl.half_height));

	// Sample 4 rim points per cap: ±perp1, ±perp2 scaled by radius.
	v3 caps[2] = { c0, c1 };
	for (int cap = 0; cap < 2; cap++) {
		for (int k = 0; k < 4 && count < MAX_CONTACTS; k++) {
			float cx = (k == 0) ? 1.0f : ((k == 2) ? -1.0f : 0.0f);
			float cy = (k == 1) ? 1.0f : ((k == 3) ? -1.0f : 0.0f);
			v3 rim_pt = add(caps[cap], add(scale(perp1, cx * cyl.radius), scale(perp2, cy * cyl.radius)));
			float signed_dist = dot(ctx->normal, sub(rim_pt, ctx->v0));
			if (signed_dist > LINEAR_SLOP) continue;
			v3 projected = sub(rim_pt, scale(ctx->normal, signed_dist));
			int feat;
			closest_point_triangle(projected, ctx->v0, ctx->v1, ctx->v2, &feat);
			if (feat != -1) continue;
			out->contacts[count].point = projected;
			out->contacts[count].normal = ctx->normal;
			out->contacts[count].penetration = -signed_dist;
			out->contacts[count].feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)((cap * 4 + k)) << 24);
			count++;
		}
	}
	if (count == 0) return 0;
	out->count = count;
	return 1;
}

// Multi-contact hull vs triangle. Same idea as the box version but iterates
// hull vertices (in world space, after applying scale + rotation + center).
// For small hulls this is cheap; for large hulls the vertex count grows linearly.
// TODO: gather only support-direction vertices via hull adjacency for speed.
static int collide_hull_triangle_multi(ConvexHull hull, const TriNeighbor* ctx, Manifold* out)
{
	int count = 0;
	const Hull* h = hull.hull;
	v3 sc = hull.scale;
	for (int i = 0; i < h->vert_count && count < MAX_CONTACTS; i++) {
		v3 local = V3(h->verts[i].x * sc.x, h->verts[i].y * sc.y, h->verts[i].z * sc.z);
		v3 world = add(hull.center, rotate(hull.rotation, local));
		float signed_dist = dot(ctx->normal, sub(world, ctx->v0));
		if (signed_dist > LINEAR_SLOP) continue;
		v3 projected = sub(world, scale(ctx->normal, signed_dist));
		int feat;
		closest_point_triangle(projected, ctx->v0, ctx->v1, ctx->v2, &feat);
		if (feat != -1) continue;
		out->contacts[count].point = projected;
		out->contacts[count].normal = ctx->normal;
		out->contacts[count].penetration = -signed_dist;
		out->contacts[count].feature_id = (uint32_t)(ctx->tri_idx & 0x00FFFFFFu) | ((uint32_t)i << 24);
		count++;
	}
	if (count == 0) return 0;
	out->count = count;
	return 1;
}

// -----------------------------------------------------------------------------
// collide_sphere_mesh_emit -- one InternalManifold per contacted triangle.
//
// Dispatch order: mesh is always body B (SHAPE_MESH is the highest enum value,
// so canonical dispatch swaps it to s1). Normal convention therefore points
// from sphere (A) toward mesh (B). Since the sphere-triangle routine produces
// a normal pointing FROM the triangle TOWARD the sphere (d = center - cp), we
// flip it before emitting.
//
// Per-triangle warm cache: each (body_pair, tri_idx+1) gets its own WarmManifold
// so independent triangles can warm-start independently.

static void collide_sphere_mesh_emit(WorldInternal* w, int body_a, int body_b, Sphere sphere_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	// Transform sphere into mesh-local space (cheaper than transforming N tris
	// into world space). Since sphere is a point+radius and orientation doesn't
	// affect it, rotation can be dropped from the transform.
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(sphere_world.center, mesh_pos));

	v3 r = V3(sphere_world.radius, sphere_world.radius, sphere_world.radius);
	AABB q = { sub(local_center, r), add(local_center, r) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);

	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		TriNeighbor ctx;
		trimesh_load_ctx(mesh, t, &ctx);
		Contact c;
		if (!collide_sphere_triangle_ctx(local_center, sphere_world.radius, &ctx, &c)) continue;

		// Flip normal to match A=sphere, B=mesh dispatch convention.
		c.normal = neg(c.normal);

		// Transform contact back to world space.
		c.point = add(mesh_pos, rotate(mesh_rot, c.point));
		c.normal = rotate(mesh_rot, c.normal);

		uint32_t sub_id = (uint32_t)(t + 1);
		InternalManifold im = { .body_a = body_a, .body_b = body_b, .sub_id = sub_id };
		im.m.contacts[0] = c;
		im.m.count = 1;
		if (w->warm_start_enabled) {
			uint64_t key = warm_cache_key(body_a, body_b, sub_id);
			im.warm = map_get_ptr(w->warm_cache, key);
		}
		apush(*manifolds, im);
	}

	afree(cands);
}

// Helper: transform a box manifold's contacts from mesh-local to world space
// and push the InternalManifold. Normal is flipped to match dispatch (A=box,
// B=mesh) so it points from box toward mesh.
static inline void trimesh_push_box_manifold(WorldInternal* w, int body_a, int body_b, int tri_idx, const Manifold* local_m, v3 mesh_pos, quat mesh_rot, InternalManifold** manifolds)
{
	uint32_t sub_id = (uint32_t)(tri_idx + 1);
	InternalManifold im = { .body_a = body_a, .body_b = body_b, .sub_id = sub_id };
	im.m = *local_m;
	for (int c = 0; c < im.m.count; c++) {
		im.m.contacts[c].normal = neg(im.m.contacts[c].normal);
		im.m.contacts[c].point = add(mesh_pos, rotate(mesh_rot, im.m.contacts[c].point));
		im.m.contacts[c].normal = rotate(mesh_rot, im.m.contacts[c].normal);
	}
	if (w->warm_start_enabled) {
		uint64_t key = warm_cache_key(body_a, body_b, sub_id);
		im.warm = map_get_ptr(w->warm_cache, key);
	}
	apush(*manifolds, im);
}

// Box-mesh: multi-contact-per-triangle with SIMD batch over 4 triangles at a
// time (when candidate count >= 4). Scalar tail for the remainder. Box's 8
// corner positions are computed once in mesh-local space and broadcast.
static void collide_box_mesh_emit(WorldInternal* w, int body_a, int body_b, Box box_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(box_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, box_world.rotation);
	Box local_box = { local_center, local_rot, box_world.half_extents };

	v3 he = box_world.half_extents;
	v3 ax = rotate(local_rot, V3(he.x, 0, 0));
	v3 ay = rotate(local_rot, V3(0, he.y, 0));
	v3 az = rotate(local_rot, V3(0, 0, he.z));
	v3 half = V3(fabsf(ax.x) + fabsf(ay.x) + fabsf(az.x), fabsf(ax.y) + fabsf(ay.y) + fabsf(az.y), fabsf(ax.z) + fabsf(ay.z) + fabsf(az.z));
	v3 slop = V3(LINEAR_SLOP, LINEAR_SLOP, LINEAR_SLOP);
	AABB q = { sub(sub(local_center, half), slop), add(add(local_center, half), slop) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	int n = asize(cands);

	// Precompute 8 box corners in mesh-local space, shared across all triangles.
	v3 corners[8];
	for (int k = 0; k < 8; k++) {
		v3 local_corner = V3((k & 1) ? he.x : -he.x, (k & 2) ? he.y : -he.y, (k & 4) ? he.z : -he.z);
		corners[k] = add(local_center, rotate(local_rot, local_corner));
	}

	int i = 0;
	// SIMD batch: 4 triangles at a time (gated for perf comparison).
	if (w->trimesh_simd_enabled) {
		for (; i + 4 <= n; i += 4) {
			TriNeighbor ctx[4];
			trimesh_load_ctx(mesh, cands[i + 0], &ctx[0]);
			trimesh_load_ctx(mesh, cands[i + 1], &ctx[1]);
			trimesh_load_ctx(mesh, cands[i + 2], &ctx[2]);
			trimesh_load_ctx(mesh, cands[i + 3], &ctx[3]);
			Manifold mfs[4];
			int hit = collide_samples_tris_batch4(corners, 8, 0.0f, ctx, mfs);
			if (!hit) continue;
			for (int lane = 0; lane < 4; lane++) {
				if (!(hit & (1 << lane))) continue;
				trimesh_push_box_manifold(w, body_a, body_b, cands[i + lane], &mfs[lane], mesh_pos, mesh_rot, manifolds);
			}
		}
	}
	// Scalar tail (also handles everything when SIMD is disabled).
	for (; i < n; i++) {
		int t = cands[i];
		TriNeighbor ctx;
		trimesh_load_ctx(mesh, t, &ctx);
		Manifold m = {0};
		if (!collide_box_triangle_multi(local_box, &ctx, &m)) continue;
		trimesh_push_box_manifold(w, body_a, body_b, t, &m, mesh_pos, mesh_rot, manifolds);
	}

	afree(cands);
}

// Hull-mesh: multi-contact via hull-vertex iteration. Vertex positions are
// precomputed once in mesh-local space and shared across 4-triangle SIMD
// batches. Scalar tail handles the last (< 4) triangles.
// Stack cap at MAX_HULL_VERTS to avoid heap alloc for typical hulls.
static void collide_hull_mesh_emit(WorldInternal* w, int body_a, int body_b, ConvexHull hull_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(hull_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, hull_world.rotation);
	ConvexHull local_hull = { hull_world.hull, local_center, local_rot, hull_world.scale };

	const Hull* h = hull_world.hull;
	v3 sc = hull_world.scale;
	int vn = h->vert_count;
	if (vn > MAX_HULL_VERTS) vn = MAX_HULL_VERTS;
	v3 samples[MAX_HULL_VERTS];
	v3 v0 = add(local_center, rotate(local_rot, V3(h->verts[0].x * sc.x, h->verts[0].y * sc.y, h->verts[0].z * sc.z)));
	samples[0] = v0;
	AABB local_aabb = { v0, v0 };
	for (int i = 1; i < vn; i++) {
		v3 vi = add(local_center, rotate(local_rot, V3(h->verts[i].x * sc.x, h->verts[i].y * sc.y, h->verts[i].z * sc.z)));
		samples[i] = vi;
		local_aabb.min = v3_min(local_aabb.min, vi);
		local_aabb.max = v3_max(local_aabb.max, vi);
	}
	v3 slop = V3(LINEAR_SLOP, LINEAR_SLOP, LINEAR_SLOP);
	AABB q = { sub(local_aabb.min, slop), add(local_aabb.max, slop) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	int n = asize(cands);

	int i = 0;
	if (w->trimesh_simd_enabled) {
		for (; i + 4 <= n; i += 4) {
			TriNeighbor ctx[4];
			trimesh_load_ctx(mesh, cands[i + 0], &ctx[0]);
			trimesh_load_ctx(mesh, cands[i + 1], &ctx[1]);
			trimesh_load_ctx(mesh, cands[i + 2], &ctx[2]);
			trimesh_load_ctx(mesh, cands[i + 3], &ctx[3]);
			Manifold mfs[4];
			int hit = collide_samples_tris_batch4(samples, vn, 0.0f, ctx, mfs);
			if (!hit) continue;
			for (int lane = 0; lane < 4; lane++) {
				if (!(hit & (1 << lane))) continue;
				trimesh_push_box_manifold(w, body_a, body_b, cands[i + lane], &mfs[lane], mesh_pos, mesh_rot, manifolds);
			}
		}
	}
	for (; i < n; i++) {
		int t = cands[i];
		TriNeighbor ctx;
		trimesh_load_ctx(mesh, t, &ctx);
		Manifold m = {0};
		if (!collide_hull_triangle_multi(local_hull, &ctx, &m)) continue;
		trimesh_push_box_manifold(w, body_a, body_b, t, &m, mesh_pos, mesh_rot, manifolds);
	}

	afree(cands);
}

// Cylinder-mesh: multi-contact via 8-rim-point sampling. SIMD-batched.
static void collide_cylinder_mesh_emit(WorldInternal* w, int body_a, int body_b, Cylinder cyl_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(cyl_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, cyl_world.rotation);
	Cylinder local_cyl = { local_center, local_rot, cyl_world.half_height, cyl_world.radius };

	v3 axis = rotate(local_rot, V3(0, 1, 0));
	v3 perp1 = rotate(local_rot, V3(1, 0, 0));
	v3 perp2 = rotate(local_rot, V3(0, 0, 1));
	float hh = cyl_world.half_height, r = cyl_world.radius;
	v3 caps[2] = { sub(local_center, scale(axis, hh)), add(local_center, scale(axis, hh)) };

	// Build 8 rim samples: 4 per cap (±perp1, ±perp2 at radius).
	v3 samples[8];
	int s = 0;
	for (int cap = 0; cap < 2; cap++) {
		samples[s++] = add(caps[cap], scale(perp1,  r));
		samples[s++] = add(caps[cap], scale(perp2,  r));
		samples[s++] = add(caps[cap], scale(perp1, -r));
		samples[s++] = add(caps[cap], scale(perp2, -r));
	}

	v3 abs_axis = V3(fabsf(axis.x), fabsf(axis.y), fabsf(axis.z));
	v3 half = V3(abs_axis.x * hh + r, abs_axis.y * hh + r, abs_axis.z * hh + r);
	v3 slop = V3(LINEAR_SLOP, LINEAR_SLOP, LINEAR_SLOP);
	AABB q = { sub(sub(local_center, half), slop), add(add(local_center, half), slop) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	int n = asize(cands);

	int i = 0;
	if (w->trimesh_simd_enabled) {
		for (; i + 4 <= n; i += 4) {
			TriNeighbor ctx[4];
			trimesh_load_ctx(mesh, cands[i + 0], &ctx[0]);
			trimesh_load_ctx(mesh, cands[i + 1], &ctx[1]);
			trimesh_load_ctx(mesh, cands[i + 2], &ctx[2]);
			trimesh_load_ctx(mesh, cands[i + 3], &ctx[3]);
			Manifold mfs[4];
			int hit = collide_samples_tris_batch4(samples, 8, 0.0f, ctx, mfs);
			if (!hit) continue;
			for (int lane = 0; lane < 4; lane++) {
				if (!(hit & (1 << lane))) continue;
				trimesh_push_box_manifold(w, body_a, body_b, cands[i + lane], &mfs[lane], mesh_pos, mesh_rot, manifolds);
			}
		}
	}
	for (; i < n; i++) {
		int t = cands[i];
		TriNeighbor ctx;
		trimesh_load_ctx(mesh, t, &ctx);
		Manifold m = {0};
		if (!collide_cylinder_triangle_multi(local_cyl, &ctx, &m)) continue;
		trimesh_push_box_manifold(w, body_a, body_b, t, &m, mesh_pos, mesh_rot, manifolds);
	}

	afree(cands);
}

// Capsule-mesh: 2 endpoint samples with radius offset = capsule.radius.
// SIMD-batched over 4 triangles, with scalar tail fallback that also handles
// mid-segment edge contacts via segment-triangle closest point (rare case).
static void collide_capsule_mesh_emit(WorldInternal* w, int body_a, int body_b, Capsule capsule_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 lp = rotate(inv_rot, sub(capsule_world.p, mesh_pos));
	v3 lq = rotate(inv_rot, sub(capsule_world.q, mesh_pos));
	float rad = capsule_world.radius;
	v3 samples[2] = { lp, lq };

	v3 r = V3(rad, rad, rad);
	v3 slop = V3(LINEAR_SLOP, LINEAR_SLOP, LINEAR_SLOP);
	AABB q = { sub(sub(v3_min(lp, lq), r), slop), add(add(v3_max(lp, lq), r), slop) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	int n = asize(cands);

	int i = 0;
	if (w->trimesh_simd_enabled) {
		for (; i + 4 <= n; i += 4) {
			TriNeighbor ctx[4];
			trimesh_load_ctx(mesh, cands[i + 0], &ctx[0]);
			trimesh_load_ctx(mesh, cands[i + 1], &ctx[1]);
			trimesh_load_ctx(mesh, cands[i + 2], &ctx[2]);
			trimesh_load_ctx(mesh, cands[i + 3], &ctx[3]);
			Manifold mfs[4];
			int hit = collide_samples_tris_batch4(samples, 2, rad, ctx, mfs);
			// Missing-lane fallback: for any lane that didn't hit via endpoint
			// sampling, try the scalar single-contact path (mid-segment edge case).
			for (int lane = 0; lane < 4; lane++) {
				if (hit & (1 << lane)) {
					trimesh_push_box_manifold(w, body_a, body_b, cands[i + lane], &mfs[lane], mesh_pos, mesh_rot, manifolds);
					continue;
				}
				Manifold m = {0};
				Contact c;
				if (!collide_capsule_triangle_ctx(lp, lq, rad, &ctx[lane], &c)) continue;
				m.contacts[0] = c;
				m.count = 1;
				trimesh_push_box_manifold(w, body_a, body_b, cands[i + lane], &m, mesh_pos, mesh_rot, manifolds);
			}
		}
	}
	for (; i < n; i++) {
		int t = cands[i];
		TriNeighbor ctx;
		trimesh_load_ctx(mesh, t, &ctx);
		Manifold m = {0};
		if (!collide_capsule_triangle_multi(lp, lq, rad, &ctx, &m)) continue;

		for (int c = 0; c < m.count; c++) {
			m.contacts[c].normal = neg(m.contacts[c].normal);
			m.contacts[c].point = add(mesh_pos, rotate(mesh_rot, m.contacts[c].point));
			m.contacts[c].normal = rotate(mesh_rot, m.contacts[c].normal);
		}

		uint32_t sub_id = (uint32_t)(t + 1);
		InternalManifold im = { .body_a = body_a, .body_b = body_b, .sub_id = sub_id };
		im.m = m;
		if (w->warm_start_enabled) {
			uint64_t key = warm_cache_key(body_a, body_b, sub_id);
			im.warm = map_get_ptr(w->warm_cache, key);
		}
		apush(*manifolds, im);
	}

	afree(cands);
}
