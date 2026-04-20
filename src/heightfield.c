// heightfield.c -- static regular-grid heightfield shape. See nudge.h for the
// user-facing contract. Pipeline mirrors trimesh.c: per-convex-shape emit
// functions transform the convex into heightfield-local space, iterate the
// cells overlapping the query AABB, build one Hull per candidate triangle
// on the stack, collide it with the convex (reusing the existing convex-vs-
// hull routines), then funnel the local manifolds through trimesh.c's
// manifold reducer to kill internal-edge ghost contacts at cell seams.
//
// Storage: float heights[N*N] in row-major (j-major) layout with cell_size
// spacing on X and Z. No quantization, no holes in v1.
//
// Triangle numbering: cell (i, j) owns tris (cell * 2) and (cell * 2 + 1)
// where cell = j * (N-1) + i. Tri 0 = diagonal-positive half:
// (i,j)-(i+1,j+1)-(i+1,j). Tri 1 = diagonal-negative half: (i,j)-(i,j+1)-(i+1,j+1).
// Both wind CCW from +Y so face normals point up. The sub-id packing used
// by InternalManifold reserves 16 bits for tri_idx + 1, so (N-1)*(N-1)*2
// must stay under 65535 -- enforced at heightfield_create (N <= 182).

struct Heightfield
{
	float*      heights;       // [N*N] row-major: heights[j*N + i]
	int         N;             // grid dimension (vertex count per axis)
	float       cell_size;     // world-space spacing per cell
	v3          aabb_min;      // local-space AABB (computed at create)
	v3          aabb_max;
	uint8_t*    material_ids;  // optional [(N-1)*(N-1)]; NULL = unset
	const char* name;          // sinterned snapshot id; NULL = unnamed
};

static inline int hf_cell_count(const Heightfield* hf) { return (hf->N - 1) * (hf->N - 1); }

Heightfield* heightfield_create(const float* heights, int N, float cell_size)
{
	assert(N >= 2 && "heightfield_create: N must be >= 2");
	assert(cell_size > 0.0f && "heightfield_create: cell_size must be > 0");
	assert((N - 1) * (N - 1) * 2 <= 65535 && "heightfield_create: N too large (sub-id limit); max N = 182");
	assert(heights && "heightfield_create: heights is NULL");

	Heightfield* hf = (Heightfield*)calloc(1, sizeof(Heightfield));
	hf->N = N;
	hf->cell_size = cell_size;

	size_t bytes = (size_t)N * (size_t)N * sizeof(float);
	hf->heights = (float*)malloc(bytes);
	memcpy(hf->heights, heights, bytes);

	float hmin =  1e30f, hmax = -1e30f;
	for (int k = 0; k < N * N; k++) {
		float h = hf->heights[k];
		if (h < hmin) hmin = h;
		if (h > hmax) hmax = h;
	}
	float extent = (float)(N - 1) * cell_size;
	hf->aabb_min = V3(0.0f,   hmin, 0.0f);
	hf->aabb_max = V3(extent, hmax, extent);
	return hf;
}

void heightfield_free(Heightfield* hf)
{
	if (!hf) return;
	free(hf->heights);
	free(hf->material_ids);
	free(hf);
}

int heightfield_tri_count(const Heightfield* hf) { return hf ? hf_cell_count(hf) * 2 : 0; }

static AABB heightfield_local_aabb(const Heightfield* hf) { return (AABB){ hf->aabb_min, hf->aabb_max }; }

void heightfield_set_material_ids(Heightfield* hf, const uint8_t* ids)
{
	free(hf->material_ids);
	hf->material_ids = NULL;
	if (!ids) return;
	size_t n = (size_t)hf_cell_count(hf);
	hf->material_ids = (uint8_t*)malloc(n);
	memcpy(hf->material_ids, ids, n);
}

uint8_t heightfield_get_material_id(const Heightfield* hf, int cell_index)
{
	if (!hf->material_ids) return 0;
	assert(cell_index >= 0 && cell_index < hf_cell_count(hf));
	return hf->material_ids[cell_index];
}

void heightfield_set_name(Heightfield* hf, const char* name) { hf->name = name ? sintern(name) : NULL; }
const char* heightfield_get_name(const Heightfield* hf) { return hf->name; }

// -----------------------------------------------------------------------------
// Local geometry helpers.

static inline v3 hf_vert(const Heightfield* hf, int i, int j)
{
	return V3((float)i * hf->cell_size, hf->heights[j * hf->N + i], (float)j * hf->cell_size);
}

// Fill v0,v1,v2 + face normal for triangle t_local in [0, 2) of cell (i, j).
// Winding: CCW from +Y (normal points up for flat ground).
static inline void hf_tri_verts(const Heightfield* hf, int i, int j, int t_local, v3* v0, v3* v1, v3* v2, v3* n_out)
{
	v3 v00 = hf_vert(hf, i, j);
	v3 v11 = hf_vert(hf, i + 1, j + 1);
	if (t_local == 0) {
		*v0 = v00;
		*v1 = v11;
		*v2 = hf_vert(hf, i + 1, j);
	} else {
		*v0 = v00;
		*v1 = hf_vert(hf, i, j + 1);
		*v2 = v11;
	}
	v3 e1 = sub(*v1, *v0), e2 = sub(*v2, *v0);
	v3 n = cross(e1, e2);
	*n_out = v3_norm(n);
}

// Build a stack-allocated Hull for one triangle. Reuses the static edge/face
// tables from trimesh.c (s_tri_edge_*, s_tri_faces) since they're pure POD.
// verts[3] and planes[2] must outlive the returned hull.
static inline Hull hf_tri_hull(v3* verts, HullPlane* planes, v3 a, v3 b, v3 c, v3 n)
{
	verts[0] = a;
	verts[1] = b;
	verts[2] = c;
	planes[0] = (HullPlane){ .normal = n,      .offset =  dot(n, a) };
	planes[1] = (HullPlane){ .normal = neg(n), .offset = -dot(n, a) };
	v3 centroid = scale(add(add(a, b), c), 1.0f / 3.0f);
	float abs_sum = fmaxf(fabsf(a.x), fmaxf(fabsf(b.x), fabsf(c.x)))
	              + fmaxf(fabsf(a.y), fmaxf(fabsf(b.y), fabsf(c.y)))
	              + fmaxf(fabsf(a.z), fmaxf(fabsf(b.z), fabsf(c.z)));
	return (Hull){
		.centroid    = centroid,
		.verts       = verts,
		.soa_verts   = NULL,
		.edge_twin   = s_tri_edge_twin,
		.edge_next   = s_tri_edge_next,
		.edge_origin = s_tri_edge_origin,
		.edge_face   = s_tri_edge_face,
		.faces       = s_tri_faces,
		.planes      = planes,
		.vert_count  = 3,
		.edge_count  = 6,
		.face_count  = 2,
		.epsilon     = 3.0f * abs_sum * FLT_EPSILON,
		.maxoutside  = 0.0f,
	};
}

// Clip a local-space AABB to the grid's cell-index range.
// Writes [i0, i1) x [j0, j1) cell bounds. Returns 0 if the AABB misses.
static int hf_cell_range(const Heightfield* hf, AABB q, int* i0, int* i1, int* j0, int* j1)
{
	float cs = hf->cell_size;
	int N = hf->N;
	if (q.max.x < 0.0f || q.max.z < 0.0f) return 0;
	if (q.min.x > (float)(N - 1) * cs) return 0;
	if (q.min.z > (float)(N - 1) * cs) return 0;
	int lo_i = (int)floorf(q.min.x / cs);
	int hi_i = (int)floorf(q.max.x / cs);
	int lo_j = (int)floorf(q.min.z / cs);
	int hi_j = (int)floorf(q.max.z / cs);
	if (lo_i < 0) lo_i = 0;
	if (lo_j < 0) lo_j = 0;
	if (hi_i > N - 2) hi_i = N - 2;
	if (hi_j > N - 2) hi_j = N - 2;
	if (lo_i > hi_i || lo_j > hi_j) return 0;
	*i0 = lo_i; *i1 = hi_i + 1;
	*j0 = lo_j; *j1 = hi_j + 1;
	return 1;
}

// Per-candidate collect: same job as trimesh_collect but using raw vertices
// since heightfield triangles aren't backed by an indices[] array. Uses the
// existing ReduceRec + tri_test_build + trimesh_reduce machinery from
// trimesh.c (we're part of the same unity TU, included after trimesh.c).
static inline void hf_collect(int tri_idx, v3 v0, v3 v1, v3 v2, v3 n, Manifold m_local, ReduceRec* recs, int* n_recs)
{
	if (m_local.count == 0) return;
	ReduceRec* r = &recs[*n_recs];
	r->tri_idx = tri_idx;
	r->m = m_local;
	r->killed = 0;
	r->corrected_by = -1;
	r->rep_contact = 0;
	tri_test_build(&r->tt, v0, v1, v2, n);
	(*n_recs)++;
}

// -----------------------------------------------------------------------------
// Narrowphase: one emit routine per convex shape, each mirroring the trimesh
// counterpart's 3-step flow (local-transform, cell-walk, reduce+flush). The
// convex argument is in world space; we rotate it into heightfield-local
// space, iterate overlapping cells, generate 2 triangles per cell, collide,
// and hand off to the shared manifold reducer.

static void collide_sphere_heightfield_emit(WorldInternal* w, int body_a, int body_b, Sphere sphere_world, v3 hf_pos, quat hf_rot, const Heightfield* hf, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(hf_rot);
	v3 local_center = rotate(inv_rot, sub(sphere_world.center, hf_pos));
	Sphere local = { local_center, sphere_world.radius };

	v3 rv = V3(sphere_world.radius, sphere_world.radius, sphere_world.radius);
	AABB q = { sub(local_center, rv), add(local_center, rv) };
	if (!aabb_overlaps(q, (AABB){ hf->aabb_min, hf->aabb_max })) return;
	int i0, i1, j0, j1;
	if (!hf_cell_range(hf, q, &i0, &i1, &j0, &j1)) return;

	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int j = j0; j < j1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; j++) {
		for (int i = i0; i < i1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
			int cell = j * (hf->N - 1) + i;
			for (int t = 0; t < 2; t++) {
				v3 a, b, c, n;
				hf_tri_verts(hf, i, j, t, &a, &b, &c, &n);
				v3 tv[3]; HullPlane tp[2];
				Hull th = hf_tri_hull(tv, tp, a, b, c, n);
				ConvexHull tch = { &th, V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
				Manifold m = {0};
				if (!collide_sphere_hull(local, tch, &m)) continue;
				hf_collect(cell * 2 + t, a, b, c, n, m, recs, &n_recs);
			}
		}
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, hf_pos, hf_rot, sub_id_base, recs, n_recs, manifolds);
}

static void collide_capsule_heightfield_emit(WorldInternal* w, int body_a, int body_b, Capsule capsule_world, v3 hf_pos, quat hf_rot, const Heightfield* hf, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(hf_rot);
	v3 lp = rotate(inv_rot, sub(capsule_world.p, hf_pos));
	v3 lq = rotate(inv_rot, sub(capsule_world.q, hf_pos));
	Capsule local = { lp, lq, capsule_world.radius };

	v3 rv = V3(capsule_world.radius, capsule_world.radius, capsule_world.radius);
	AABB q = { sub(v3_min(lp, lq), rv), add(v3_max(lp, lq), rv) };
	if (!aabb_overlaps(q, (AABB){ hf->aabb_min, hf->aabb_max })) return;
	int i0, i1, j0, j1;
	if (!hf_cell_range(hf, q, &i0, &i1, &j0, &j1)) return;

	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int j = j0; j < j1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; j++) {
		for (int i = i0; i < i1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
			int cell = j * (hf->N - 1) + i;
			for (int t = 0; t < 2; t++) {
				v3 a, b, c, n;
				hf_tri_verts(hf, i, j, t, &a, &b, &c, &n);
				v3 tv[3]; HullPlane tp[2];
				Hull th = hf_tri_hull(tv, tp, a, b, c, n);
				ConvexHull tch = { &th, V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
				Manifold m = {0};
				if (!collide_capsule_hull(local, tch, &m)) continue;
				hf_collect(cell * 2 + t, a, b, c, n, m, recs, &n_recs);
			}
		}
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, hf_pos, hf_rot, sub_id_base, recs, n_recs, manifolds);
}

static void collide_box_heightfield_emit(WorldInternal* w, int body_a, int body_b, Box box_world, v3 hf_pos, quat hf_rot, const Heightfield* hf, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(hf_rot);
	v3 local_center = rotate(inv_rot, sub(box_world.center, hf_pos));
	quat local_rot = mul(inv_rot, box_world.rotation);

	// OBB AABB bound: project each oriented half-extent axis onto world axes.
	v3 he = box_world.half_extents;
	v3 ax = rotate(local_rot, V3(he.x, 0, 0));
	v3 ay = rotate(local_rot, V3(0, he.y, 0));
	v3 az = rotate(local_rot, V3(0, 0, he.z));
	v3 half = V3(fabsf(ax.x) + fabsf(ay.x) + fabsf(az.x), fabsf(ax.y) + fabsf(ay.y) + fabsf(az.y), fabsf(ax.z) + fabsf(ay.z) + fabsf(az.z));
	AABB q = { sub(local_center, half), add(local_center, half) };
	if (!aabb_overlaps(q, (AABB){ hf->aabb_min, hf->aabb_max })) return;
	int i0, i1, j0, j1;
	if (!hf_cell_range(hf, q, &i0, &i1, &j0, &j1)) return;

	ConvexHull box_hull = { &s_unit_box_hull, local_center, local_rot, he };

	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int j = j0; j < j1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; j++) {
		for (int i = i0; i < i1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
			int cell = j * (hf->N - 1) + i;
			for (int t = 0; t < 2; t++) {
				v3 a, b, c, n;
				hf_tri_verts(hf, i, j, t, &a, &b, &c, &n);
				v3 tv[3]; HullPlane tp[2];
				Hull th = hf_tri_hull(tv, tp, a, b, c, n);
				ConvexHull tch = { &th, V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
				Manifold m = {0};
				if (!collide_hull_hull(box_hull, tch, &m)) continue;
				hf_collect(cell * 2 + t, a, b, c, n, m, recs, &n_recs);
			}
		}
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, hf_pos, hf_rot, sub_id_base, recs, n_recs, manifolds);
}

static void collide_hull_heightfield_emit(WorldInternal* w, int body_a, int body_b, ConvexHull hull_world, v3 hf_pos, quat hf_rot, const Heightfield* hf, uint32_t sub_id_base, InternalManifold** manifolds)
{
	quat inv_rot = inv(hf_rot);
	v3 local_center = rotate(inv_rot, sub(hull_world.center, hf_pos));
	quat local_rot = mul(inv_rot, hull_world.rotation);
	ConvexHull local = { hull_world.hull, local_center, local_rot, hull_world.scale };

	const Hull* h = hull_world.hull;
	v3 sc = hull_world.scale;
	v3 v0w = add(local_center, rotate(local_rot, V3(h->verts[0].x * sc.x, h->verts[0].y * sc.y, h->verts[0].z * sc.z)));
	AABB q = { v0w, v0w };
	for (int vi = 1; vi < h->vert_count; vi++) {
		v3 vw = add(local_center, rotate(local_rot, V3(h->verts[vi].x * sc.x, h->verts[vi].y * sc.y, h->verts[vi].z * sc.z)));
		q.min = v3_min(q.min, vw);
		q.max = v3_max(q.max, vw);
	}
	if (!aabb_overlaps(q, (AABB){ hf->aabb_min, hf->aabb_max })) return;
	int i0, i1, j0, j1;
	if (!hf_cell_range(hf, q, &i0, &i1, &j0, &j1)) return;

	ReduceRec recs[TRIMESH_MAX_LOCAL_MANIFOLDS]; int n_recs = 0;
	for (int j = j0; j < j1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; j++) {
		for (int i = i0; i < i1 && n_recs + 2 <= TRIMESH_MAX_LOCAL_MANIFOLDS; i++) {
			int cell = j * (hf->N - 1) + i;
			for (int t = 0; t < 2; t++) {
				v3 a, b, c, n;
				hf_tri_verts(hf, i, j, t, &a, &b, &c, &n);
				v3 tv[3]; HullPlane tp[2];
				Hull th = hf_tri_hull(tv, tp, a, b, c, n);
				ConvexHull tch = { &th, V3(0, 0, 0), quat_identity(), V3(1, 1, 1) };
				Manifold m = {0};
				if (!collide_hull_hull(local, tch, &m)) continue;
				hf_collect(cell * 2 + t, a, b, c, n, m, recs, &n_recs);
			}
		}
	}
	trimesh_reduce(recs, n_recs);
	trimesh_flush(w, body_a, body_b, hf_pos, hf_rot, sub_id_base, recs, n_recs, manifolds);
}

// -----------------------------------------------------------------------------
// Raycast: walk cells along the ray, Moeller-Trumbore on each triangle.

static int ray_heightfield(v3 ro, v3 rd, v3 hf_pos, quat hf_rot, const Heightfield* hf, float max_t, float* t_out, v3* n_out)
{
	quat inv_rot = inv(hf_rot);
	v3 lo = rotate(inv_rot, sub(ro, hf_pos));
	v3 ld = rotate(inv_rot, rd);

	// Coarse AABB clip against heightfield extent.
	v3 inv_d = V3(ld.x != 0.0f ? 1.0f / ld.x : 1e30f,
	              ld.y != 0.0f ? 1.0f / ld.y : 1e30f,
	              ld.z != 0.0f ? 1.0f / ld.z : 1e30f);
	float tmin = 0.0f, tmax = max_t;
	for (int ax = 0; ax < 3; ax++) {
		float o = ((float*)&lo)[ax];
		float d = ((float*)&inv_d)[ax];
		float a = (((float*)&hf->aabb_min)[ax] - o) * d;
		float b = (((float*)&hf->aabb_max)[ax] - o) * d;
		if (a > b) { float tmp = a; a = b; b = tmp; }
		if (a > tmin) tmin = a;
		if (b < tmax) tmax = b;
		if (tmin > tmax) return 0;
	}

	// Brute-force per-triangle. Acceptable for v1 -- a grid DDA walk is
	// straightforward once profiling calls for it.
	float best_t = max_t;
	v3 best_n = V3(0, 1, 0);
	int hit = 0;
	int cells = hf_cell_count(hf);
	for (int cell = 0; cell < cells; cell++) {
		int i = cell % (hf->N - 1);
		int j = cell / (hf->N - 1);
		for (int t = 0; t < 2; t++) {
			v3 a, b, c, n;
			hf_tri_verts(hf, i, j, t, &a, &b, &c, &n);
			v3 e1 = sub(b, a), e2 = sub(c, a);
			v3 h = cross(ld, e2);
			float det = dot(e1, h);
			if (fabsf(det) < 1e-12f) continue;
			float f = 1.0f / det;
			v3 s = sub(lo, a);
			float u = f * dot(s, h);
			if (u < 0.0f || u > 1.0f) continue;
			v3 qv = cross(s, e1);
			float vv = f * dot(ld, qv);
			if (vv < 0.0f || u + vv > 1.0f) continue;
			float tt = f * dot(e2, qv);
			if (tt < 0.0f || tt > best_t) continue;
			best_t = tt;
			best_n = n;
			hit = 1;
		}
	}
	if (!hit) return 0;
	*t_out = best_t;
	*n_out = rotate(hf_rot, best_n);
	return 1;
}
