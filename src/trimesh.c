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
// Front face edges (CCW around +n): 0: v0->v1, 1: v1->v2, 2: v2->v0.
// Back face edges (CCW around -n):  3: v0->v2, 4: v2->v1, 5: v1->v0.
// Twin pairs: (0, 5), (1, 4), (2, 3) -- opposite half-edges across the disc.

static const int s_tri_edge_twin[6]   = { 5, 4, 3, 2, 1, 0 };
static const int s_tri_edge_next[6]   = { 1, 2, 0, 4, 5, 3 };
static const int s_tri_edge_origin[6] = { 0, 1, 2, 0, 2, 1 };
static const int s_tri_edge_face[6]   = { 0, 0, 0, 1, 1, 1 };
static const HullFace s_tri_faces[2]  = { { .edge = 0 }, { .edge = 3 } };

// -----------------------------------------------------------------------------
// TriMesh layout. Single heap block with:
//   [TriMesh header][verts[V]][indices[3T]][tri_normal[T]]
//   [hull_verts[3T]][hull_planes[2T]][hulls[T]]

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
};

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
	free(m);
}

int trimesh_tri_count(const TriMesh* m) { return m ? m->tri_count : 0; }

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
	// DIAGNOSTIC: bypass BVH, linear scan every triangle.
	for (int t = 0; t < m->tri_count; t++) {
		uint32_t i0 = m->indices[3*t + 0], i1 = m->indices[3*t + 1], i2 = m->indices[3*t + 2];
		v3 a = m->verts[i0], b = m->verts[i1], c = m->verts[i2];
		AABB ta = { v3_min(a, v3_min(b, c)), v3_max(a, v3_max(b, c)) };
		if (aabb_overlaps(ta, q)) apush(*out, t);
	}
}

// -----------------------------------------------------------------------------
// Emit helpers. Each mesh emit:
//   1. Transform convex into mesh-local space (cheaper than moving N hulls).
//   2. BVH-query the convex's mesh-local AABB.
//   3. For each candidate triangle, call the matching convex-vs-hull routine
//      using the pre-built triangle Hull, emit one InternalManifold per hit.

// Push a per-triangle manifold with local-space contacts transformed to world.
// collide_{sphere,capsule,hull,cylinder}_hull already returns Contact.normal
// in A-toward-B convention (convex toward triangle-hull), which matches the
// dispatch's expected A=convex, B=mesh orientation -- no normal flip needed.
static inline void trimesh_push_mesh_manifold(WorldInternal* w, int body_a, int body_b, int tri_idx, Manifold m_local, v3 mesh_pos, quat mesh_rot, InternalManifold** manifolds)
{
	if (m_local.count == 0) return;
	uint32_t sub_id = (uint32_t)(tri_idx + 1);
	InternalManifold im = { .body_a = body_a, .body_b = body_b, .sub_id = sub_id };
	im.m = m_local;
	for (int c = 0; c < im.m.count; c++) {
		im.m.contacts[c].point = add(mesh_pos, rotate(mesh_rot, im.m.contacts[c].point));
		im.m.contacts[c].normal = rotate(mesh_rot, im.m.contacts[c].normal);
	}
	// DIAGNOSTIC: mesh warm-start off to isolate warm-cache contribution.
	(void)sub_id;
	apush(*manifolds, im);
}

// Positioned wrapper: the convex is already in mesh-local space, so the
// triangle Hull sits at identity transform.
static inline ConvexHull trimesh_tri_as_hull_local(const TriMesh* mesh, int t)
{
	return (ConvexHull){ &mesh->hulls[t], V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
}

// -----------------------------------------------------------------------------
// Sphere vs mesh.
static void collide_sphere_mesh_emit(WorldInternal* w, int body_a, int body_b, Sphere sphere_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
{
	quat inv_rot = inv(mesh_rot);
	v3 local_center = rotate(inv_rot, sub(sphere_world.center, mesh_pos));
	Sphere local = { local_center, sphere_world.radius };

	v3 rv = V3(sphere_world.radius, sphere_world.radius, sphere_world.radius);
	AABB q = { sub(local_center, rv), add(local_center, rv) };
	if (!aabb_overlaps(q, (AABB){ mesh->aabb_min, mesh->aabb_max })) return;

	CK_DYNA int* cands = NULL;
	trimesh_query_aabb(mesh, q, &cands);
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_sphere_hull(local, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_push_mesh_manifold(w, body_a, body_b, t, m, mesh_pos, mesh_rot, manifolds);
	}
	afree(cands);
}

// -----------------------------------------------------------------------------
// Capsule vs mesh.
static void collide_capsule_mesh_emit(WorldInternal* w, int body_a, int body_b, Capsule capsule_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
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
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_capsule_hull(local, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_push_mesh_manifold(w, body_a, body_b, t, m, mesh_pos, mesh_rot, manifolds);
	}
	afree(cands);
}

// -----------------------------------------------------------------------------
// Box vs mesh. Box is wrapped as a unit-box hull scaled by half_extents.
static void collide_box_mesh_emit(WorldInternal* w, int body_a, int body_b, Box box_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
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
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_hull_hull(box_hull, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_push_mesh_manifold(w, body_a, body_b, t, m, mesh_pos, mesh_rot, manifolds);
	}
	afree(cands);
}

// -----------------------------------------------------------------------------
// Hull vs mesh.
static void collide_hull_mesh_emit(WorldInternal* w, int body_a, int body_b, ConvexHull hull_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
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
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_hull_hull(local_hull, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_push_mesh_manifold(w, body_a, body_b, t, m, mesh_pos, mesh_rot, manifolds);
	}
	afree(cands);
}

// -----------------------------------------------------------------------------
// Cylinder vs mesh.
static void collide_cylinder_mesh_emit(WorldInternal* w, int body_a, int body_b, Cylinder cyl_world, v3 mesh_pos, quat mesh_rot, const TriMesh* mesh, InternalManifold** manifolds)
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
	for (int i = 0; i < asize(cands); i++) {
		int t = cands[i];
		Manifold m = {0};
		if (!collide_cylinder_hull(local, trimesh_tri_as_hull_local(mesh, t), &m)) continue;
		trimesh_push_mesh_manifold(w, body_a, body_b, t, m, mesh_pos, mesh_rot, manifolds);
	}
	afree(cands);
}
