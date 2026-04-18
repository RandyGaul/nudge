// See LICENSE for licensing info.
// trimesh.c -- static triangle-mesh collision geometry.
//
// Each triangle is a degenerate closed Hull (zero-volume, double-sided disc:
// 3 verts, 2 faces with opposing normals, 6 half-edges wired into proper twin
// pairs). Every convex-vs-triangle pair reuses the engine's existing
// convex-vs-hull routines -- no custom trimesh narrowphase.
//
// Layout:
//   * Half-edge topology is identical for every triangle; shared via static
//     module-level arrays (6 half-edges, 2 faces).
//   * Per-triangle: 3 verts + 2 planes + one Hull header pointing at those
//     plus the shared topology.
//   * TriMesh owns a single heap block with all per-tri arrays contiguous.
//
// Collision: dispatch iterates AABB-BVH candidates; each candidate triangle
// becomes a ConvexHull (pointing at its stored Hull) passed to the existing
// sphere-hull / capsule-hull / hull-hull / cylinder-hull routines. Output
// is one InternalManifold per (body_pair, triangle).

// -----------------------------------------------------------------------------
// Shared half-edge topology for all triangles.
//
// Twin pairs at adjacent indices (2k, 2k+1) -- engine-wide convention required
// by sat_query_edges' undirected-edge iteration.
//   pair 0 (v0-v1): 0 front v0->v1 / 1 back v1->v0
//   pair 1 (v1-v2): 2 front v1->v2 / 3 back v2->v1
//   pair 2 (v2-v0): 4 front v2->v0 / 5 back v0->v2
// Front face (face 0, +n) CCW: 0 -> 2 -> 4. Back face (face 1, -n) CCW: 5 -> 3 -> 1.

static const int s_tri_edge_twin[6]   = { 1, 0, 3, 2, 5, 4 };
static const int s_tri_edge_next[6]   = { 2, 5, 4, 1, 0, 3 };
static const int s_tri_edge_origin[6] = { 0, 1, 1, 2, 2, 0 };
static const int s_tri_edge_face[6]   = { 0, 1, 0, 1, 0, 1 };
static const HullFace s_tri_faces[2]  = { { .edge = 0 }, { .edge = 5 } };

// -----------------------------------------------------------------------------
// TriMesh layout. Single heap block with:
//   [TriMesh header][verts[V]][indices[3T]][tri_normal[T]]
//   [hull_verts[3T]][hull_planes[2T]][hulls[T]]

// Internal-edge / ghost-contact handling: Bepu-style MeshReduction.
//
// After narrowphase produces one manifold per contacted triangle for a
// body-vs-mesh pair, we run a reduction pass that tests each manifold's
// representative (deepest) contact against every OTHER triangle in the same
// pair. If a contact is physically "inside" a neighbor triangle's face region
// AND its normal is infringing on that neighbor's outward face (or outward
// edge-planes for edge-region contacts), the source manifold is a ghost and
// gets deleted. In the wedge case (both manifolds block each other) the
// source's normal is replaced with the neighbor's face normal.
//
// No mesh-bake precompute -- runs purely per-frame on collected manifolds.
// Reference: BepuPhysics/CollisionDetection/MeshReduction.cs.
struct TriMesh
{
	v3 aabb_min;
	v3 aabb_max;

	const v3* verts;
	const uint32_t* indices;
	const v3* tri_normal;

	const v3* hull_verts;         // [3 * tri_count]
	const HullPlane* hull_planes; // [2 * tri_count]
	const Hull* hulls;            // [tri_count]

	int vert_count;
	int tri_count;

	BVH_Tree bvh;

	const char* name;  // sinterned name for snapshot identification; NULL = unnamed

	// Per-triangle material ids; NULL when unset (callers get body default).
	// Heap-allocated outside the main block so it remains optional without
	// forcing a rebuild of the mesh. Length equals tri_count when present.
	uint8_t* material_ids;
};

void trimesh_set_name(TriMesh* mesh, const char* name)
{
	mesh->name = name ? sintern(name) : NULL;
}

const char* trimesh_get_name(const TriMesh* mesh) { return mesh->name; }

// -----------------------------------------------------------------------------
// trimesh_create -- build per-triangle Hull structs + BVH.

TriMesh* trimesh_create(const v3* verts, int vert_count, const uint32_t* indices, int tri_count)
{
	assert(vert_count > 0 && tri_count > 0);
	assert(verts && indices);

	size_t verts_bytes       = (size_t)vert_count * sizeof(v3);
	size_t indices_bytes     = (size_t)tri_count * 3 * sizeof(uint32_t);
	size_t tri_normal_bytes  = (size_t)tri_count * sizeof(v3);
	size_t hull_verts_bytes  = (size_t)tri_count * 3 * sizeof(v3);
	size_t hull_planes_bytes = (size_t)tri_count * 2 * sizeof(HullPlane);
	size_t hulls_bytes       = (size_t)tri_count * sizeof(Hull);

	#define ALN16(x) (((x) + 15) & ~(size_t)15)
	size_t off = ALN16(sizeof(TriMesh));
	size_t off_verts       = off;                                  off += ALN16(verts_bytes);
	size_t off_indices     = off;                                  off += ALN16(indices_bytes);
	size_t off_tri_normal  = off;                                  off += ALN16(tri_normal_bytes);
	size_t off_hull_verts  = off;                                  off += ALN16(hull_verts_bytes);
	size_t off_hull_planes = off;                                  off += ALN16(hull_planes_bytes);
	size_t off_hulls       = off;                                  off += ALN16(hulls_bytes);
	size_t total = off;
	#undef ALN16

	char* block = (char*)calloc(1, total);
	TriMesh* m = (TriMesh*)block;

	v3* mverts            = (v3*)(block + off_verts);
	uint32_t* mindices    = (uint32_t*)(block + off_indices);
	v3* mnormals          = (v3*)(block + off_tri_normal);
	v3* mhull_verts       = (v3*)(block + off_hull_verts);
	HullPlane* mhull_planes = (HullPlane*)(block + off_hull_planes);
	Hull* mhulls          = (Hull*)(block + off_hulls);
	m->verts       = mverts;
	m->indices     = mindices;
	m->tri_normal  = mnormals;
	m->hull_verts  = mhull_verts;
	m->hull_planes = mhull_planes;
	m->hulls       = mhulls;
	m->vert_count  = vert_count;
	m->tri_count   = tri_count;

	memcpy(mverts, verts, verts_bytes);
	memcpy(mindices, indices, indices_bytes);

	v3 amin = V3(1e18f, 1e18f, 1e18f), amax = V3(-1e18f, -1e18f, -1e18f);
	for (int t = 0; t < tri_count; t++) {
		uint32_t i0 = mindices[3*t + 0], i1 = mindices[3*t + 1], i2 = mindices[3*t + 2];
		assert((int)i0 < vert_count && (int)i1 < vert_count && (int)i2 < vert_count);
		v3 a = mverts[i0], b = mverts[i1], c = mverts[i2];
		v3 n = cross(sub(b, a), sub(c, a));
		float len2 = dot(n, n);
		assert(len2 > 1e-20f);
		v3 nu = scale(n, 1.0f / sqrtf(len2));
		mnormals[t] = nu;

		mhull_verts[3*t + 0] = a;
		mhull_verts[3*t + 1] = b;
		mhull_verts[3*t + 2] = c;

		mhull_planes[2*t + 0] = (HullPlane){ .normal = nu,      .offset =  dot(nu, a) };
		mhull_planes[2*t + 1] = (HullPlane){ .normal = neg(nu), .offset = -dot(nu, a) };

		v3 centroid = scale(add(add(a, b), c), 1.0f / 3.0f);
		float abs_sum = fmaxf(fabsf(a.x), fmaxf(fabsf(b.x), fabsf(c.x)))
		              + fmaxf(fabsf(a.y), fmaxf(fabsf(b.y), fabsf(c.y)))
		              + fmaxf(fabsf(a.z), fmaxf(fabsf(b.z), fabsf(c.z)));
		mhulls[t] = (Hull){
			.centroid    = centroid,
			.verts       = &mhull_verts[3*t],
			.soa_verts   = NULL,
			.edge_twin   = s_tri_edge_twin,
			.edge_next   = s_tri_edge_next,
			.edge_origin = s_tri_edge_origin,
			.edge_face   = s_tri_edge_face,
			.faces       = s_tri_faces,
			.planes      = &mhull_planes[2*t],
			.vert_count  = 3,
			.edge_count  = 6,
			.face_count  = 2,
			.epsilon     = 3.0f * abs_sum * FLT_EPSILON,
			.maxoutside  = 0.0f,
		};

		amin = v3_min(amin, v3_min(a, v3_min(b, c)));
		amax = v3_max(amax, v3_max(a, v3_max(b, c)));
	}
	m->aabb_min = amin;
	m->aabb_max = amax;

	// Internal BVH over triangle AABBs.
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
	free(m->material_ids);
	free(m);
}

int trimesh_tri_count(const TriMesh* m) { return m ? m->tri_count : 0; }

void trimesh_set_material_ids(TriMesh* mesh, const uint8_t* ids)
{
	free(mesh->material_ids);
	mesh->material_ids = NULL;
	if (!ids) return;
	size_t bytes = (size_t)mesh->tri_count * sizeof(uint8_t);
	mesh->material_ids = (uint8_t*)malloc(bytes);
	memcpy(mesh->material_ids, ids, bytes);
}

uint8_t trimesh_get_material_id(const TriMesh* mesh, int tri_index)
{
	if (!mesh->material_ids) return 0;
	assert(tri_index >= 0 && tri_index < mesh->tri_count);
	return mesh->material_ids[tri_index];
}

static AABB trimesh_local_aabb(const TriMesh* m) { return (AABB){ m->aabb_min, m->aabb_max }; }

// -----------------------------------------------------------------------------
// Raycast against mesh. BVH-pruned, Moeller-Trumbore per leaf.

static int ray_triangle(v3 ro, v3 rd, v3 v0, v3 v1, v3 v2, float max_t, float* t_out)
{
	v3 e1 = sub(v1, v0);
	v3 e2 = sub(v2, v0);
	v3 h = cross(rd, e2);
	float a = dot(e1, h);
	if (fabsf(a) < 1e-12f) return 0;
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

static void trimesh_query_ray_node(const BVH_Tree* t, int ni, v3 ro, v3 inv_dir, float max_t, CK_DYNA int** out)
{
	const BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		const BVH_Child* c = s == 0 ? &n->a : &n->b;
		if (c->leaf_count == 0) continue;
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
// BVH AABB traversal (candidate triangle gather).

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
// Emit helpers. Each mesh emit:
//   1. Transform convex into mesh-local space (cheaper than moving N hulls).
//   2. BVH-query the convex's mesh-local AABB.
//   3. For each candidate triangle, call the matching convex-vs-hull routine
//      using the pre-built triangle Hull, emit one InternalManifold per hit.

// --- Bepu-style MeshReduction -----------------------------------------------
//
// Per-triangle test data: face normal + 3 outward edge-plane normals. A
// contact at point P is "inside" this triangle's face region when P lies
// close to the face plane AND on the interior side of all three edge planes.
// If additionally the contact's normal infringes on this triangle (points
// into the face, or on an edge the contact is touching, points out of that
// edge plane), the contact is a ghost from a neighbor and should be killed.
typedef struct TriTest
{
	v3 anchor;    // verts[0], reference point for plane distances
	v3 face_n;    // unit face normal
	v3 edge_n[3]; // outward edge-plane normals (in-plane, perpendicular to each edge)
	float thr;    // distance threshold -- scale-relative (~1e-3 * triangle bbox)
} TriTest;

static void tri_test_build(TriTest* tt, v3 v0, v3 v1, v3 v2, v3 face_n)
{
	tt->anchor = v0;
	tt->face_n = face_n;
	v3 verts[3] = { v0, v1, v2 };
	for (int k = 0; k < 3; k++) {
		v3 a = verts[k], b = verts[(k + 1) % 3];
		v3 e = sub(b, a);
		// outward edge-plane normal: perpendicular to edge, in the face plane,
		// pointing AWAY from triangle interior. cross(face_n, edge) gives this
		// when face_n is the CCW outward normal and the edge is traversed in
		// CCW order; the triangle interior sits on the negative side.
		v3 en = cross(face_n, e);
		float l2 = dot(en, en);
		tt->edge_n[k] = l2 > 1e-20f ? scale(en, 1.0f / sqrtf(l2)) : V3(0, 0, 0);
	}
	// Scale-relative threshold: use the longest edge's length, 1e-3 factor.
	float max_e2 = 0.0f;
	for (int k = 0; k < 3; k++) {
		v3 e = sub(verts[(k + 1) % 3], verts[k]);
		float l2 = dot(e, e);
		if (l2 > max_e2) max_e2 = l2;
	}
	tt->thr = 1e-3f * sqrtf(max_e2);
	if (tt->thr < 1e-5f) tt->thr = 1e-5f;
}

// True if `cp` (in mesh-local space) is "inside" the face region of `target`
// AND contact_n (A->B direction, body->mesh) is infringing on target. The
// caller passes a candidate representative contact of the SOURCE manifold.
//
// Face infringement: dot(contact_n, target.face_n) < -thr (normal points
// into the face -- i.e. source says "the body is on the OTHER side of this
// face", which is a ghost signal).
//
// Edge classification: for each edge plane, distance to plane (positive =
// outside triangle, negative = inside). Contact is "on edge k" if
// |edgeDist[k]| <= thr. It is "inside the triangle face region" if all
// edgeDist <= thr.
//
// Edge infringement: for each edge the contact is on, dot(contact_n,
// edge_n[k]) must be positive (contact's normal points OUT of that edge
// plane -- a ghost telling you "the body is past this edge"). The primary
// (closest) edge gets a strict > eps test; secondary edges just need
// parallel (>= -1e-2) to avoid gaps at shared vertices.
//
// Returns 1 if target blocks source.
static int tri_test_blocks(const TriTest* target, v3 cp, v3 cn)
{
	// Face plane distance.
	float face_d = dot(sub(cp, target->anchor), target->face_n);
	if (face_d > target->thr || face_d < -target->thr) return 0; // not in face plane

	// Edge plane distances (signed, positive = outside).
	float ed[3];
	for (int k = 0; k < 3; k++) ed[k] = dot(sub(cp, target->anchor), target->edge_n[k]);
	// "Near the triangle face region": all edge distances <= thr.
	if (ed[0] > target->thr || ed[1] > target->thr || ed[2] > target->thr) return 0;

	// Infringement test.
	float edge_on_thr = -target->thr * 0.01f; // "on edge k" when ed[k] >= edge_on_thr
	int on_edge_count = 0;
	int primary_edge = -1;
	float primary_ed = -1e18f;
	for (int k = 0; k < 3; k++) {
		if (ed[k] >= edge_on_thr) {
			on_edge_count++;
			if (ed[k] > primary_ed) { primary_ed = ed[k]; primary_edge = k; }
		}
	}

	if (on_edge_count == 0) {
		// Strictly inside -- face infringement.
		return dot(cn, target->face_n) < -1e-3f;
	}

	// On one or more edges. Primary edge: strict positive normal dot with
	// outward edge plane. Secondary edges: just parallel (>= small neg tol).
	if (dot(cn, target->edge_n[primary_edge]) <= 1e-6f) return 0;
	for (int k = 0; k < 3; k++) {
		if (k == primary_edge || ed[k] < edge_on_thr) continue;
		if (dot(cn, target->edge_n[k]) < -1e-2f) return 0;
	}
	return 1;
}

// Source record for one triangle manifold during reduction.
typedef struct ReduceRec
{
	int tri_idx;
	Manifold m;               // contacts in mesh-local space
	TriTest tt;               // precomputed test data for THIS triangle
	int rep_contact;          // index of the deepest contact in m.contacts
	uint8_t killed;           // deleted by a neighbor
	int8_t corrected_by;      // >=0: index of neighbor whose face_n replaces this normal
} ReduceRec;

// Reduce an array of per-triangle manifolds for one body-mesh pair in-place.
// Deletes ghost manifolds (sets m.count = 0) and, in the wedge case, rewrites
// the manifold normal with the neighbor's face normal.
static void trimesh_reduce(ReduceRec* recs, int n)
{
	// One-sided-mesh enforcement (Bepu/Jolt default). A triangle-hull is built
	// with two opposing faces so SAT can work, but physically the mesh is a
	// surface, not a volume. When a body crosses through, SAT picks the back
	// face and returns a normal pointing from the body *toward* the body's
	// new side of the plane (A->B convention with mesh=B, body below plane ->
	// normal = +face_n). The solver then applies -normal to the body, which
	// is along -face_n -- pushing the body *further away from the surface*
	// in the tunneled direction, NOT back through to where it came from.
	// Result: catastrophic tunneling.
	//
	// Fix: for every contact whose normal has a positive component along the
	// source tri's face_n (i.e. back-face contact, body has tunneled), flip
	// the normal to -face_n so the solver pushes the body back toward the
	// front side. This is the "one-sided triangle" behaviour that Bepu and
	// Jolt use as their default. Penetration depth is kept as-is; fixing
	// the direction is what matters.
	for (int i = 0; i < n; i++) {
		ReduceRec* si = &recs[i];
		if (si->m.count == 0) continue;
		v3 front_facing = neg(si->tt.face_n); // A->B normal for a front-face contact
		for (int c = 0; c < si->m.count; c++) {
			v3 cn = si->m.contacts[c].normal;
			if (dot(cn, si->tt.face_n) > 0.05f) {
				// Back-face / tunneled contact. Override normal so solver
				// pushes body back toward the front side.
				si->m.contacts[c].normal = front_facing;
			}
		}
	}

	if (n < 2) return;
	// For each source i, find its deepest contact and test against every
	// other triangle j's TriTest.
	for (int i = 0; i < n; i++) {
		ReduceRec* si = &recs[i];
		if (si->m.count == 0) continue;
		// Pick deepest contact as representative.
		int best = 0;
		for (int c = 1; c < si->m.count; c++) {
			if (si->m.contacts[c].penetration > si->m.contacts[best].penetration) best = c;
		}
		si->rep_contact = best;
	}
	for (int i = 0; i < n; i++) {
		ReduceRec* si = &recs[i];
		if (si->m.count == 0) continue;
		v3 cp = si->m.contacts[si->rep_contact].point;
		v3 cn = si->m.contacts[si->rep_contact].normal;
		// FACE-CONTACT TRUST (Bepu has an explicit FaceCollisionFlag on contacts
		// set by narrowphase; we approximate by checking the source manifold's
		// own normal alignment with its own face normal). A contact whose normal
		// already points -face_n (A->B convention, body toward mesh surface) IS
		// a genuine face contact -- trust it, skip the block test.
		if (dot(neg(cn), si->tt.face_n) > 0.999f) continue;
		for (int j = 0; j < n; j++) {
			if (j == i) continue;
			ReduceRec* sj = &recs[j];
			if (sj->m.count == 0) continue;
			if (tri_test_blocks(&sj->tt, cp, cn)) {
				si->killed = 1;
				si->corrected_by = (int8_t)j;
				break;
			}
		}
	}
	// Resolve:
	//   killed and neighbor NOT killed by source -> delete source.
	//   both killed (wedge) -> replace source's normal with neighbor's face.
	for (int i = 0; i < n; i++) {
		ReduceRec* si = &recs[i];
		if (!si->killed) continue;
		int j = si->corrected_by;
		ReduceRec* sj = &recs[j];
		if (sj->killed && sj->corrected_by == i) {
			// Wedge: both block each other. Replace normal.
			v3 replacement_n = neg(sj->tt.face_n); // A->B convention (body toward mesh)
			for (int c = 0; c < si->m.count; c++) si->m.contacts[c].normal = replacement_n;
		} else {
			// One-sided block: kill source.
			si->m.count = 0;
		}
	}
}

// Push a single per-triangle manifold (contacts still in mesh-local space)
// into a local reduction buffer.
static inline void trimesh_collect(const TriMesh* mesh, int tri_idx, Manifold m_local, ReduceRec* recs, int* n)
{
	if (m_local.count == 0) return;
	ReduceRec* r = &recs[*n];
	r->tri_idx = tri_idx;
	r->m = m_local;
	r->killed = 0;
	r->corrected_by = -1;
	r->rep_contact = 0;
	uint32_t i0 = mesh->indices[3*tri_idx + 0];
	uint32_t i1 = mesh->indices[3*tri_idx + 1];
	uint32_t i2 = mesh->indices[3*tri_idx + 2];
	tri_test_build(&r->tt, mesh->verts[i0], mesh->verts[i1], mesh->verts[i2], mesh->tri_normal[tri_idx]);
	(*n)++;
}

// Flush survivors to the caller's manifold list, transforming to world space
// as we go.
static inline void trimesh_flush(WorldInternal* w, int body_a, int body_b, v3 mesh_pos, quat mesh_rot, uint32_t sub_id_base, ReduceRec* recs, int n, InternalManifold** manifolds)
{
	for (int i = 0; i < n; i++) {
		ReduceRec* r = &recs[i];
		if (r->m.count == 0) continue;
		uint32_t sub_id = sub_id_base | (uint32_t)(r->tri_idx + 1);
		InternalManifold im = { .body_a = body_a, .body_b = body_b, .sub_id = sub_id };
		im.m = r->m;
		for (int c = 0; c < im.m.count; c++) {
			im.m.contacts[c].point = add(mesh_pos, rotate(mesh_rot, im.m.contacts[c].point));
			im.m.contacts[c].normal = rotate(mesh_rot, im.m.contacts[c].normal);
		}
		apush(*manifolds, im);
	}
}

#define TRIMESH_MAX_LOCAL_MANIFOLDS 64

// Positioned wrapper: the convex is already in mesh-local space, so the
// triangle Hull sits at identity transform.
static inline ConvexHull trimesh_tri_as_hull_local(const TriMesh* mesh, int t)
{
	return (ConvexHull){ &mesh->hulls[t], V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
}

// -----------------------------------------------------------------------------
// Sphere vs mesh.
static void collide_sphere_mesh_emit(WorldInternal* w, int body_a, int body_b, Sphere sphere_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(sphere_world.center, mesh_pos));
	Sphere local = { local_center, sphere_world.radius };

	v3 rv = V3(sphere_world.radius, sphere_world.radius, sphere_world.radius);
	AABB q = { sub(local_center, rv), add(local_center, rv) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int i = 0; i < asize(cands) && n_recs < TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_sphere_hull(local, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_collect(mesh, t, m, recs, &n_recs);
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, mesh_pos, mesh_rot, sub_id_base, recs, n_recs, manifolds);
	afree(cands);
}

// -----------------------------------------------------------------------------
// Capsule vs mesh.
static void collide_capsule_mesh_emit(WorldInternal* w, int body_a, int body_b, Capsule capsule_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 lp = rotate(inv_rot, sub(capsule_world.p, mesh_pos));
	v3 lq = rotate(inv_rot, sub(capsule_world.q, mesh_pos));
	Capsule local = { lp, lq, capsule_world.radius };

	v3 rv = V3(capsule_world.radius, capsule_world.radius, capsule_world.radius);
	AABB q = { sub(v3_min(lp, lq), rv), add(v3_max(lp, lq), rv) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int i = 0; i < asize(cands) && n_recs < TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_capsule_hull(local, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_collect(mesh, t, m, recs, &n_recs);
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, mesh_pos, mesh_rot, sub_id_base, recs, n_recs, manifolds);
	afree(cands);
}

// -----------------------------------------------------------------------------
// Box vs mesh. Box is wrapped as a unit-box hull scaled by half_extents.
static void collide_box_mesh_emit(WorldInternal* w, int body_a, int body_b, Box box_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(box_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, box_world.rotation);

	v3 he = box_world.half_extents;
	v3 ax = rotate(local_rot, V3(he.x, 0, 0));
	v3 ay = rotate(local_rot, V3(0, he.y, 0));
	v3 az = rotate(local_rot, V3(0, 0, he.z));
	v3 half = V3(fabsf(ax.x) + fabsf(ay.x) + fabsf(az.x), fabsf(ax.y) + fabsf(ay.y) + fabsf(az.y), fabsf(ax.z) + fabsf(ay.z) + fabsf(az.z));
	AABB q = { sub(local_center, half), add(local_center, half) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	ConvexHull box_hull = { &s_unit_box_hull, local_center, local_rot, he };
	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int i = 0; i < asize(cands) && n_recs < TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_hull_hull(box_hull, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_collect(mesh, t, m, recs, &n_recs);
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, mesh_pos, mesh_rot, sub_id_base, recs, n_recs, manifolds);
	afree(cands);
}

// -----------------------------------------------------------------------------
// Hull vs mesh.
static void collide_hull_mesh_emit(WorldInternal* w, int body_a, int body_b, ConvexHull hull_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(hull_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, hull_world.rotation);
	ConvexHull local_hull = { hull_world.hull, local_center, local_rot, hull_world.scale };

	const Hull* h = hull_world.hull;
	v3 sc = hull_world.scale;
	v3 v0w = add(local_center, rotate(local_rot, V3(h->verts[0].x * sc.x, h->verts[0].y * sc.y, h->verts[0].z * sc.z)));
	AABB q = { v0w, v0w };
	for (int i = 1; i < h->vert_count; i++) {
		v3 vi = add(local_center, rotate(local_rot, V3(h->verts[i].x * sc.x, h->verts[i].y * sc.y, h->verts[i].z * sc.z)));
		q.min = v3_min(q.min, vi);
		q.max = v3_max(q.max, vi);
	}
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int i = 0; i < asize(cands) && n_recs < TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_hull_hull(local_hull, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_collect(mesh, t, m, recs, &n_recs);
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, mesh_pos, mesh_rot, sub_id_base, recs, n_recs, manifolds);
	afree(cands);
}

// -----------------------------------------------------------------------------
// Cylinder vs flat triangle (single-triangle collision, mesh-local space).
// Replaces collide_cylinder_hull for triangle-mesh because the generic hull
// routine treats a zero-volume triangle as a degenerate polyhedron and its
// edge-edge SAT family produces false-positive lateral contacts (edge axes
// treat the edge as an infinite 3D prism, not a 2D segment).
//
// Strategy:
//   1. Signed distance from cyl center to tri plane; reject if separated.
//   2. Contact normal = tri plane normal, oriented from cyl toward tri.
//   3. CAP-like contacts (cyl axis ~parallel to tri normal): emit up to 4
//      rim-point projections AND up to 3 tri-vertex projections that lie
//      inside the cylinder disk, clipping rim points to the triangle.
//   4. SIDE-like contacts (cyl axis ~perpendicular to tri normal): emit two
//      contacts at the cylinder axis endpoints offset by radius toward tri.
//   5. Depth = (cyl_extent_along_normal) - |signed_distance|.
//
// All contacts share the same face normal per triangle -- no false lateral
// normals from edge-prism SAT.
static int point_in_triangle_2d(v3 p, v3 v0, v3 v1, v3 v2, v3 n)
{
	// Barycentric via normal-projected cross products.
	v3 e0 = sub(v1, v0), e1 = sub(v2, v1), e2 = sub(v0, v2);
	v3 c0 = cross(e0, sub(p, v0));
	v3 c1 = cross(e1, sub(p, v1));
	v3 c2 = cross(e2, sub(p, v2));
	float d0 = dot(c0, n), d1 = dot(c1, n), d2 = dot(c2, n);
	return (d0 >= 0.0f && d1 >= 0.0f && d2 >= 0.0f) || (d0 <= 0.0f && d1 <= 0.0f && d2 <= 0.0f);
}

// Sample cylinder surface at representative points (both cap rims + axial side
// band), project each to triangle plane, keep those below the plane AND inside
// the triangle. Reduce to MAX_CONTACTS deepest. This unifies CAP/SIDE handling:
// vertical cyls dominated by cap-rim contacts, horizontal by axial, tilted get
// both -- no discontinuity at any tilt.
static int collide_cylinder_triangle_local(Cylinder cyl, v3 v0, v3 v1, v3 v2, v3 tri_n, Manifold* m)
{
	v3 cyl_axis = rotate(cyl.rotation, V3(0, 1, 0));
	float hh = cyl.half_height, r = cyl.radius;

	// Signed distance of cyl center from tri plane (positive = same side as tri_n).
	float plane_off = dot(tri_n, v0);
	float c_signed_d = dot(cyl.center, tri_n) - plane_off;

	// Cylinder extent along +/- tri_n.
	float an = dot(cyl_axis, tri_n);
	float axial = fabsf(an) * hh;
	float sin2 = 1.0f - an * an;
	float radial = sin2 > 0.0f ? sqrtf(sin2) * r : 0.0f;
	float cyl_extent = axial + radial;

	float c_abs_d = fabsf(c_signed_d);
	float gap = c_abs_d - cyl_extent;
	if (gap > LINEAR_SLOP) return 0;

	// Contact normal: cyl A toward tri B. Cyl on +tri_n side -> points -tri_n.
	v3 contact_n = c_signed_d >= 0.0f ? neg(tri_n) : tri_n;
	int on_plus_side = c_signed_d >= 0.0f;

	// Build orthonormal basis for cylinder: axis + two perpendicular directions.
	v3 probe = fabsf(cyl_axis.y) < 0.9f ? V3(0, 1, 0) : V3(1, 0, 0);
	v3 perp1 = norm(cross(cyl_axis, probe));
	v3 perp2 = cross(cyl_axis, perp1);

	// Candidate contact samples on the cylinder surface. Combine: cap rims at
	// both ends (8 directions each) + axial samples along the side (5 axial
	// positions at 2 opposing perpendicular directions = 10 samples).
	// Keep only those below tri plane AND inside triangle.
	#define CYL_TRI_MAX_SAMPLES 32
	v3 s_points[CYL_TRI_MAX_SAMPLES]; float s_depths[CYL_TRI_MAX_SAMPLES];
	int sc = 0;

	// Cap rims at +hh and -hh along cyl_axis (8 directions per cap).
	float cap_axials[2] = { -hh, +hh };
	int rim_dirs_n = 8;
	for (int cap = 0; cap < 2; cap++) {
		v3 cap_center = add(cyl.center, scale(cyl_axis, cap_axials[cap]));
		for (int i = 0; i < rim_dirs_n && sc < CYL_TRI_MAX_SAMPLES; i++) {
			float theta = (float)i / (float)rim_dirs_n * 6.28318530718f;
			v3 rim = add(cap_center, add(scale(perp1, cosf(theta) * r), scale(perp2, sinf(theta) * r)));
			float d_signed = dot(rim, tri_n) - plane_off;
			if (on_plus_side && d_signed > 0.0f) continue;
			if (!on_plus_side && d_signed < 0.0f) continue;
			v3 on_plane = sub(rim, scale(tri_n, d_signed));
			if (!point_in_triangle_2d(on_plane, v0, v1, v2, tri_n)) continue;
			s_points[sc] = on_plane;
			s_depths[sc] = fabsf(d_signed);
			sc++;
		}
	}

	// Axial side samples along cylinder's deepest contact line: at each of
	// 5 axial positions, the cyl surface point nearest the tri plane.
	int axial_samples = 5;
	for (int i = 0; i < axial_samples && sc < CYL_TRI_MAX_SAMPLES; i++) {
		float t = (float)i / (float)(axial_samples - 1) * 2.0f - 1.0f;
		v3 on_axis = add(cyl.center, scale(cyl_axis, t * hh));
		// Deepest perpendicular offset: project contact_n onto axis-perp plane, scale to r.
		float ap = dot(contact_n, cyl_axis);
		v3 perp = sub(contact_n, scale(cyl_axis, ap));
		float pl2 = len2(perp);
		v3 offset = pl2 > 1e-8f ? scale(perp, r / sqrtf(pl2)) : V3(0, 0, 0);
		v3 on_surface = add(on_axis, offset);
		float d_signed = dot(on_surface, tri_n) - plane_off;
		if (on_plus_side && d_signed > 0.0f) continue;
		if (!on_plus_side && d_signed < 0.0f) continue;
		v3 on_plane = sub(on_surface, scale(tri_n, d_signed));
		if (!point_in_triangle_2d(on_plane, v0, v1, v2, tri_n)) continue;
		s_points[sc] = on_plane;
		s_depths[sc] = fabsf(d_signed);
		sc++;
	}

	if (sc == 0) return 0;

	// Reduce to MAX_CONTACTS by picking the 4 deepest and well-separated points.
	// Simple greedy: sort by depth desc, dedupe by spatial distance.
	int cp = 0;
	v3 points[MAX_CONTACTS]; float depths[MAX_CONTACTS];
	for (int pass = 0; pass < MAX_CONTACTS; pass++) {
		int best = -1; float best_d = -1e18f;
		for (int i = 0; i < sc; i++) {
			if (s_depths[i] < 0.0f) continue; // already consumed
			int too_close = 0;
			for (int j = 0; j < cp; j++) {
				if (len2(sub(s_points[i], points[j])) < 1e-4f) { too_close = 1; break; }
			}
			if (too_close) { s_depths[i] = -1.0f; continue; }
			if (s_depths[i] > best_d) { best_d = s_depths[i]; best = i; }
		}
		if (best < 0) break;
		points[cp] = s_points[best]; depths[cp] = s_depths[best];
		s_depths[best] = -1.0f; cp++;
	}
	if (cp == 0) return 0;

	m->count = cp;
	for (int i = 0; i < cp; i++) m->contacts[i] = (Contact){ .point = points[i], .normal = contact_n, .penetration = depths[i] };
	return 1;
}

// Cylinder vs mesh.
static void collide_cylinder_mesh_emit(WorldInternal* w, int body_a, int body_b, Cylinder cyl_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(cyl_world.center, mesh_pos));
	quat local_rot = mul(inv_rot, cyl_world.rotation);
	Cylinder local = { local_center, local_rot, cyl_world.half_height, cyl_world.radius };

	v3 axis = rotate(local_rot, V3(0, 1, 0));
	float hh = cyl_world.half_height, r = cyl_world.radius;
	v3 abs_axis = V3(fabsf(axis.x), fabsf(axis.y), fabsf(axis.z));
	v3 half = V3(abs_axis.x * hh + r, abs_axis.y * hh + r, abs_axis.z * hh + r);
	AABB q = { sub(local_center, half), add(local_center, half) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int i = 0; i < asize(cands) && n_recs < TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
		int t = cands[i];
		uint32_t i0 = mesh->indices[3*t + 0], i1 = mesh->indices[3*t + 1], i2 = mesh->indices[3*t + 2];
		v3 v0 = mesh->verts[i0], v1 = mesh->verts[i1], v2 = mesh->verts[i2];
		Manifold m = {0};
		if (!collide_cylinder_triangle_local(local, v0, v1, v2, mesh->tri_normal[t], &m)) continue;
		trimesh_collect(mesh, t, m, recs, &n_recs);
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, mesh_pos, mesh_rot, sub_id_base, recs, n_recs, manifolds);
	afree(cands);
}
