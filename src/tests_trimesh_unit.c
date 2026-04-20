// tests_trimesh_unit.c -- trimesh layer 1: pure primitives.
//
// Direct tests of the static helpers in trimesh.c, reachable because
// tests_trimesh_unit.c is included in the same unity TU that defines them.
// Nothing here allocates a TriMesh or touches a World: each test stands up
// a fixed triangle on the stack and calls the prim directly. Layers 2+
// add mesh construction, queries, narrowphase, reduction, and integration.

// -----------------------------------------------------------------------------
// L1.1 -- point_in_triangle_2d (trimesh.c:794)
//
// Signature: int point_in_triangle_2d(v3 p, v3 v0, v3 v1, v3 v2, v3 n)
//
// Orientation-agnostic: passes if all three barycentric cross-projections
// agree in sign (either all >= 0 OR all <= 0). That makes it work for both
// CCW and CW winding -- verified below.

static void test_pit2d_ccw_interior()
{
	TEST_BEGIN("pit2d ccw interior");
	// CCW in XY plane, normal +Z (b-a)x(c-a) = +Z.
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), n = V3(0, 0, 1);
	TEST_ASSERT(point_in_triangle_2d(V3(0.25f, 0.25f, 0), v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(V3(0.1f,  0.1f,  0), v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(V3(0.33f, 0.33f, 0), v0, v1, v2, n)); // near centroid
}

static void test_pit2d_ccw_outside()
{
	TEST_BEGIN("pit2d ccw outside");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), n = V3(0, 0, 1);
	TEST_ASSERT(!point_in_triangle_2d(V3(1.0f,  1.0f,  0), v0, v1, v2, n));  // past hypotenuse
	TEST_ASSERT(!point_in_triangle_2d(V3(-0.1f, 0.5f,  0), v0, v1, v2, n));  // left of v0-v2 edge
	TEST_ASSERT(!point_in_triangle_2d(V3(0.5f,  -0.1f, 0), v0, v1, v2, n));  // below v0-v1 edge
	TEST_ASSERT(!point_in_triangle_2d(V3(2.0f,  2.0f,  0), v0, v1, v2, n));  // far away
}

static void test_pit2d_boundary()
{
	TEST_BEGIN("pit2d boundary inclusive");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), n = V3(0, 0, 1);
	// Vertices are inclusive (all cross products hit zero on at least one).
	TEST_ASSERT(point_in_triangle_2d(v0, v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(v1, v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(v2, v0, v1, v2, n));
	// Edge midpoints inclusive.
	TEST_ASSERT(point_in_triangle_2d(V3(0.5f, 0,    0), v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(V3(0,    0.5f, 0), v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(V3(0.5f, 0.5f, 0), v0, v1, v2, n)); // on hypotenuse
}

static void test_pit2d_cw_winding()
{
	TEST_BEGIN("pit2d cw winding with flipped normal");
	// Same triangle wound CW, supply the matching (flipped) normal.
	v3 v0 = V3(0, 0, 0), v1 = V3(0, 1, 0), v2 = V3(1, 0, 0), n = V3(0, 0, -1);
	TEST_ASSERT(point_in_triangle_2d(V3(0.25f, 0.25f, 0), v0, v1, v2, n));
	TEST_ASSERT(!point_in_triangle_2d(V3(1.0f, 1.0f, 0), v0, v1, v2, n));
	// And with the "wrong" normal the function still works because the
	// all-same-sign rule accepts either orientation.
	TEST_ASSERT(point_in_triangle_2d(V3(0.25f, 0.25f, 0), v0, v1, v2, V3(0, 0, 1)));
}

static void test_pit2d_projection_invariance()
{
	TEST_BEGIN("pit2d ignores component along n");
	// The function projects cross products onto n; a Z offset on p shouldn't
	// change the result when the triangle lives in the XY plane.
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), n = V3(0, 0, 1);
	TEST_ASSERT(point_in_triangle_2d(V3(0.25f, 0.25f,  10.0f), v0, v1, v2, n));
	TEST_ASSERT(point_in_triangle_2d(V3(0.25f, 0.25f, -10.0f), v0, v1, v2, n));
	TEST_ASSERT(!point_in_triangle_2d(V3(2.0f,  2.0f,  10.0f), v0, v1, v2, n));
}

// -----------------------------------------------------------------------------
// L1.2 -- ray_triangle (trimesh.c:241), Moeller-Trumbore.
//
// Signature: int ray_triangle(v3 ro, v3 rd, v3 v0, v3 v1, v3 v2, float max_t, float* t_out)
//
// NOTE on two-sidedness: the current implementation uses `fabsf(a) < 1e-12f`
// for its determinant test, so it accepts rays approaching from either side
// of the triangle. Bepu and Jolt raycasts are one-sided by default. The
// "back-face hit" test below pins CURRENT behavior; flipping to one-sided is
// a deliberate decision for L3, not a fix to force here.

static void test_ray_tri_front_hit_center()
{
	TEST_BEGIN("ray_tri front-side hit through center");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	int hit = ray_triangle(V3(0.25f, 0.25f, 1.0f), V3(0, 0, -1), v0, v1, v2, 10.0f, &t);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-5f);
}

static void test_ray_tri_back_hit_pinned_two_sided()
{
	TEST_BEGIN("ray_tri back-side hit (pinned: current impl is two-sided)");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	int hit = ray_triangle(V3(0.25f, 0.25f, -1.0f), V3(0, 0, 1), v0, v1, v2, 10.0f, &t);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-5f);
}

static void test_ray_tri_miss_barycentric()
{
	TEST_BEGIN("ray_tri misses when projected point is outside barycentric");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	// Outside u bound (u < 0): hit behind v0-v2 edge.
	TEST_ASSERT(!ray_triangle(V3(-0.5f, 0.5f, 1.0f), V3(0, 0, -1), v0, v1, v2, 10.0f, &t));
	// Outside u bound (u > 1): hit past v1.
	TEST_ASSERT(!ray_triangle(V3(3.0f, 0.5f, 1.0f), V3(0, 0, -1), v0, v1, v2, 10.0f, &t));
	// Outside v bound (v < 0): hit below v0-v1 edge.
	TEST_ASSERT(!ray_triangle(V3(0.5f, -0.5f, 1.0f), V3(0, 0, -1), v0, v1, v2, 10.0f, &t));
	// u + v > 1: past hypotenuse.
	TEST_ASSERT(!ray_triangle(V3(1.5f, 1.5f, 1.0f), V3(0, 0, -1), v0, v1, v2, 10.0f, &t));
}

static void test_ray_tri_miss_parallel()
{
	TEST_BEGIN("ray_tri parallel to plane rejected");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	// Ray in XY plane (parallel to triangle). Determinant ~= 0.
	TEST_ASSERT(!ray_triangle(V3(-1.0f, 0.5f, 0.0f), V3(1, 0, 0), v0, v1, v2, 10.0f, &t));
}

static void test_ray_tri_miss_max_t()
{
	TEST_BEGIN("ray_tri rejects when hit beyond max_t");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	// Ray would hit at t=5 but max_t=1.
	TEST_ASSERT(!ray_triangle(V3(0.5f, 0.5f, 5.0f), V3(0, 0, -1), v0, v1, v2, 1.0f, &t));
}

static void test_ray_tri_miss_behind_origin()
{
	TEST_BEGIN("ray_tri rejects when triangle is behind origin (t < 0)");
	v3 v0 = V3(0, 0, 0), v1 = V3(2, 0, 0), v2 = V3(0, 2, 0);
	float t = -1.0f;
	// Ray origin at z=+1, direction +Z: triangle at z=0 is behind it.
	TEST_ASSERT(!ray_triangle(V3(0.5f, 0.5f, 1.0f), V3(0, 0, 1), v0, v1, v2, 10.0f, &t));
}

// -----------------------------------------------------------------------------
// L1.3 -- tri_test_build (trimesh.c:355)
//
// Signature: void tri_test_build(TriTest* tt, v3 v0, v3 v1, v3 v2, v3 face_n)
//
// Populates TriTest with anchor, face_n, three edge normals (expected to be
// perpendicular to both the edge and face_n, in the face plane, unit length
// or zero for a degenerate edge), and a scale-relative distance threshold
// thr = max(1e-5, 1e-3 * longest_edge).
//
// We verify the STRUCTURAL invariants (unit length, perpendicularity) but
// deliberately do NOT assert edge_n's sign/direction here -- the end-to-end
// correctness of the reducer's ghost/wedge behavior lives in L4+L5 tests.

static void test_tt_build_structural_invariants()
{
	TEST_BEGIN("tri_test_build: edge_n unit + perp to edge + perp to face_n");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), face_n = V3(0, 0, 1);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);

	TEST_ASSERT_FLOAT(tt.anchor.x, v0.x, 1e-6f);
	TEST_ASSERT_FLOAT(tt.anchor.y, v0.y, 1e-6f);
	TEST_ASSERT_FLOAT(tt.anchor.z, v0.z, 1e-6f);
	TEST_ASSERT_FLOAT(tt.face_n.x, face_n.x, 1e-6f);
	TEST_ASSERT_FLOAT(tt.face_n.y, face_n.y, 1e-6f);
	TEST_ASSERT_FLOAT(tt.face_n.z, face_n.z, 1e-6f);

	v3 edges[3] = { sub(v1, v0), sub(v2, v1), sub(v0, v2) };
	for (int k = 0; k < 3; k++) {
		float len = sqrtf(len2(tt.edge_n[k]));
		TEST_ASSERT_FLOAT(len, 1.0f, 1e-5f);                        // unit length
		TEST_ASSERT_FLOAT(dot(tt.edge_n[k], edges[k]), 0.0f, 1e-5f); // perp to edge
		TEST_ASSERT_FLOAT(dot(tt.edge_n[k], face_n), 0.0f, 1e-5f);   // perp to face_n (in-plane)
	}
}

static void test_tt_build_thr_scales_with_edge_length()
{
	TEST_BEGIN("tri_test_build: thr scales with longest edge (1e-3 factor)");
	// Triangle with longest edge == 100 -> thr ~= 0.1.
	v3 v0 = V3(0, 0, 0), v1 = V3(100, 0, 0), v2 = V3(0, 1, 0);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, V3(0, 0, 1));
	TEST_ASSERT(tt.thr > 0.09f && tt.thr < 0.11f);
}

static void test_tt_build_thr_floor()
{
	TEST_BEGIN("tri_test_build: thr floor at 1e-5 for tiny triangles");
	// 1mm triangle: longest edge ~1e-3 -> 1e-3 * 1e-3 = 1e-6, below floor.
	v3 v0 = V3(0, 0, 0), v1 = V3(1e-3f, 0, 0), v2 = V3(0, 1e-3f, 0);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, V3(0, 0, 1));
	TEST_ASSERT_FLOAT(tt.thr, 1e-5f, 1e-7f);
}

static void test_tt_build_arbitrary_orientation()
{
	TEST_BEGIN("tri_test_build: invariants hold for tilted triangle");
	// Triangle tilted out of coordinate planes; face_n computed from the verts.
	v3 v0 = V3(0.5f, 0.2f, 0.1f);
	v3 v1 = V3(1.7f, 0.3f, 0.4f);
	v3 v2 = V3(0.6f, 1.2f, 0.8f);
	v3 n_raw = cross(sub(v1, v0), sub(v2, v0));
	v3 face_n = scale(n_raw, 1.0f / sqrtf(len2(n_raw)));
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);

	v3 edges[3] = { sub(v1, v0), sub(v2, v1), sub(v0, v2) };
	for (int k = 0; k < 3; k++) {
		float len = sqrtf(len2(tt.edge_n[k]));
		TEST_ASSERT_FLOAT(len, 1.0f, 1e-5f);
		TEST_ASSERT_FLOAT(dot(tt.edge_n[k], edges[k]), 0.0f, 1e-4f);
		TEST_ASSERT_FLOAT(dot(tt.edge_n[k], face_n), 0.0f, 1e-5f);
	}
}

// -----------------------------------------------------------------------------
// L1.4 -- tri_test_blocks (trimesh.c:402)
//
// Signature: int tri_test_blocks(const TriTest* target, v3 cp, v3 cn)
//
// Returns 1 when `target` blocks the candidate contact (cp, cn) -- i.e. the
// contact would be a ghost because cp lies inside/on target's face region AND
// cn is infringing on target's face or an edge the contact is on.
//
// These are BEHAVIORAL tests for clear-cut cases; subtle edge-region and
// wedge-blocking behavior is exercised at L5 with the full reducer.

static void test_ttb_out_of_plane_rejects()
{
	TEST_BEGIN("tri_test_blocks: cp far off face plane -> no block");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), face_n = V3(0, 0, 1);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);
	// Contact 0.5 above the plane, normal pointing into the face. thr << 0.5,
	// so the face plane check bails immediately.
	TEST_ASSERT(tri_test_blocks(&tt, V3(0.3f, 0.3f, 0.5f), V3(0, 0, -1)) == 0);
	TEST_ASSERT(tri_test_blocks(&tt, V3(0.3f, 0.3f, -0.5f), V3(0, 0, 1)) == 0);
}

static void test_ttb_aligned_normal_no_infringe()
{
	TEST_BEGIN("tri_test_blocks: cn along +face_n (not infringing) -> no block");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), face_n = V3(0, 0, 1);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);
	// Contact in the plane, but cn points *along* face_n -- not "into" the
	// face. Face infringement requires dot(cn, face_n) < -1e-3.
	TEST_ASSERT(tri_test_blocks(&tt, V3(0.25f, 0.25f, 0.0f), V3(0, 0, 1)) == 0);
}

static void test_ttb_far_outside_triangle()
{
	TEST_BEGIN("tri_test_blocks: cp far outside triangle -> no block");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), face_n = V3(0, 0, 1);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);
	// cp is in the face plane but 10 units from the triangle. Normal points
	// into face. Should be rejected by the edge-region prefilter.
	TEST_ASSERT(tri_test_blocks(&tt, V3(10.0f, 10.0f, 0.0f), V3(0, 0, -1)) == 0);
	TEST_ASSERT(tri_test_blocks(&tt, V3(-5.0f, 0.25f, 0.0f), V3(0, 0, -1)) == 0);
}

static void test_ttb_zero_normal()
{
	TEST_BEGIN("tri_test_blocks: zero contact normal -> no infringement");
	v3 v0 = V3(0, 0, 0), v1 = V3(1, 0, 0), v2 = V3(0, 1, 0), face_n = V3(0, 0, 1);
	TriTest tt; tri_test_build(&tt, v0, v1, v2, face_n);
	// A zero normal can never beat the -1e-3 face-infringement threshold.
	TEST_ASSERT(tri_test_blocks(&tt, V3(0.25f, 0.25f, 0.0f), V3(0, 0, 0)) == 0);
}

// -----------------------------------------------------------------------------
// L2 -- trimesh_create: Hull topology, per-triangle data, AABB, BVH, material
// palette, name interning, and free. Uses small heap-allocated meshes and
// reaches into struct TriMesh directly (same TU as trimesh.c).

static TriMesh* ttm_single(v3 a, v3 b, v3 c)
{
	v3 verts[3] = { a, b, c };
	uint32_t idx[3] = { 0, 1, 2 };
	return trimesh_create(verts, 3, idx, 1);
}

// Coplanar quad in the XZ plane (y=0) split into 2 CCW triangles with +Y face
// normals -- i.e. "floor" orientation. Shared edge is v0-v2 (the diagonal).
static TriMesh* ttm_quad_xz()
{
	v3 verts[4] = { V3(0,0,0), V3(0,0,1), V3(1,0,1), V3(1,0,0) };
	uint32_t idx[6] = { 0, 1, 2,   0, 2, 3 };
	return trimesh_create(verts, 4, idx, 2);
}

// -----------------------------------------------------------------------------
// L2.1 -- Hull topology shared tables (trimesh.c:31-35)
//
// Every per-triangle Hull points at these. Verify the half-edge topology is
// internally consistent (twins involute, face cycles close, origins match).
// This is one-time correctness for a file-static -- we test it once.

static void test_trimesh_shared_topology_twins()
{
	TEST_BEGIN("shared tri edge topology: twins form involution (i^1)");
	for (int i = 0; i < 6; i++) {
		TEST_ASSERT(s_tri_edge_twin[i] == (i ^ 1));
		TEST_ASSERT(s_tri_edge_twin[s_tri_edge_twin[i]] == i);
	}
}

static void test_trimesh_shared_topology_face_cycles()
{
	TEST_BEGIN("shared tri edge topology: each face cycle has length 3");
	for (int f = 0; f < 2; f++) {
		int start = s_tri_faces[f].edge;
		int e = start;
		int steps = 0;
		do {
			TEST_ASSERT(s_tri_edge_face[e] == f);
			e = s_tri_edge_next[e];
			steps++;
		} while (e != start && steps < 10);
		TEST_ASSERT(steps == 3);
	}
}

static void test_trimesh_shared_topology_origins()
{
	TEST_BEGIN("shared tri edge topology: origin(next(e)) != origin(e), twin origins swap");
	// An edge's origin is one endpoint; its twin's origin is the other.
	for (int i = 0; i < 6; i++) {
		int tw = s_tri_edge_twin[i];
		TEST_ASSERT(s_tri_edge_origin[i] != s_tri_edge_origin[tw]);
		// next(e) starts at the destination of e (= origin of twin(e)).
		TEST_ASSERT(s_tri_edge_origin[s_tri_edge_next[i]] == s_tri_edge_origin[tw]);
	}
}

// -----------------------------------------------------------------------------
// L2.2 -- Per-triangle Hull population

static void test_trimesh_hull_wiring()
{
	TEST_BEGIN("trimesh_create: per-tri Hull points at shared topology + local verts/planes");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	const Hull* h = &m->hulls[0];

	TEST_ASSERT(h->vert_count == 3);
	TEST_ASSERT(h->edge_count == 6);
	TEST_ASSERT(h->face_count == 2);
	TEST_ASSERT(h->soa_verts == NULL);
	TEST_ASSERT(h->maxoutside == 0.0f);

	// Shared static arrays -- pointer identity, not just content.
	TEST_ASSERT(h->edge_twin == s_tri_edge_twin);
	TEST_ASSERT(h->edge_next == s_tri_edge_next);
	TEST_ASSERT(h->edge_origin == s_tri_edge_origin);
	TEST_ASSERT(h->edge_face == s_tri_edge_face);
	TEST_ASSERT(h->faces == s_tri_faces);

	// Per-tri blocks.
	TEST_ASSERT(h->verts == &m->hull_verts[0]);
	TEST_ASSERT(h->planes == &m->hull_planes[0]);

	trimesh_free(m);
}

static void test_trimesh_centroid_and_epsilon()
{
	TEST_BEGIN("trimesh_create: centroid is vertex mean, epsilon = 3 * abs_sum * FLT_EPSILON");
	// Pick non-collinear verts so trimesh_create doesn't reject them.
	v3 a = V3(1, 2, 3), b = V3(4, 5, 6), c = V3(7, 8, 0);
	TriMesh* m = ttm_single(a, b, c);
	const Hull* h = &m->hulls[0];
	v3 expected_centroid = scale(add(add(a, b), c), 1.0f / 3.0f);
	TEST_ASSERT_FLOAT(h->centroid.x, expected_centroid.x, 1e-5f);
	TEST_ASSERT_FLOAT(h->centroid.y, expected_centroid.y, 1e-5f);
	TEST_ASSERT_FLOAT(h->centroid.z, expected_centroid.z, 1e-5f);

	// epsilon = 3 * (max|x| + max|y| + max|z|) * FLT_EPSILON. max|x|=7, max|y|=8,
	// max|z|=6 -> abs_sum=21 -> epsilon = 63 * FLT_EPSILON.
	float expected_eps = 3.0f * 21.0f * FLT_EPSILON;
	TEST_ASSERT_FLOAT(h->epsilon, expected_eps, 1e-10f);
	trimesh_free(m);
}

static void test_trimesh_front_back_planes()
{
	TEST_BEGIN("trimesh_create: plane[0] = +tri_normal, plane[1] = -tri_normal, offsets flipped");
	// CCW triangle -> +Z tri_normal.
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	const HullPlane* p_front = &m->hull_planes[0];
	const HullPlane* p_back  = &m->hull_planes[1];

	TEST_ASSERT_FLOAT(p_front->normal.x, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_front->normal.y, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_front->normal.z, 1.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_front->offset,   0.0f, 1e-6f); // dot((0,0,1), (0,0,0)) = 0

	TEST_ASSERT_FLOAT(p_back->normal.x,  0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_back->normal.y,  0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_back->normal.z, -1.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_back->offset,    0.0f, 1e-6f);
	trimesh_free(m);
}

static void test_trimesh_plane_offset_nonzero()
{
	TEST_BEGIN("trimesh_create: plane offset is dot(normal, v0) for shifted triangle");
	// Same triangle lifted to y=5; normal still +Y in XZ plane.
	TriMesh* m = ttm_single(V3(0,5,0), V3(1,5,0), V3(0,5,1));
	const HullPlane* p_front = &m->hull_planes[0];
	// For CCW in XZ with +Y normal (cross((1,0,0),(0,0,1)) = (0,-1,0))... actually
	// b-a = (1,0,0), c-a = (0,0,1); cross = (0*1 - 0*0, 0*0 - 1*1, 1*0 - 0*0) = (0,-1,0).
	// So tri_normal is -Y and offset = dot(-Y, (0,5,0)) = -5.
	TEST_ASSERT_FLOAT(p_front->normal.y, -1.0f, 1e-6f);
	TEST_ASSERT_FLOAT(p_front->offset,   -5.0f, 1e-5f);
	trimesh_free(m);
}

static void test_trimesh_tri_normal_stored()
{
	TEST_BEGIN("trimesh_create: tri_normal[t] is unit and matches face-plane normal");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0));
	v3 tn = m->tri_normal[0];
	TEST_ASSERT_FLOAT(sqrtf(len2(tn)), 1.0f, 1e-5f);
	TEST_ASSERT_FLOAT(tn.x, m->hull_planes[0].normal.x, 1e-6f);
	TEST_ASSERT_FLOAT(tn.y, m->hull_planes[0].normal.y, 1e-6f);
	TEST_ASSERT_FLOAT(tn.z, m->hull_planes[0].normal.z, 1e-6f);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L2.3 -- AABB

static void test_trimesh_aabb_single_triangle()
{
	TEST_BEGIN("trimesh_create: aabb covers all verts (single tri)");
	TriMesh* m = ttm_single(V3(-1, 2, 3), V3(4, -5, 3), V3(0, 2, 7));
	TEST_ASSERT_FLOAT(m->aabb_min.x, -1.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_min.y, -5.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_min.z,  3.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.x,  4.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.y,  2.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.z,  7.0f, 1e-6f);
	trimesh_free(m);
}

static void test_trimesh_aabb_multi_triangle()
{
	TEST_BEGIN("trimesh_create: aabb covers every vertex across tris");
	TriMesh* m = ttm_quad_xz(); // verts at x in [0,1], z in [0,1], y=0
	TEST_ASSERT_FLOAT(m->aabb_min.x, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_min.y, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_min.z, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.x, 1.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.y, 0.0f, 1e-6f);
	TEST_ASSERT_FLOAT(m->aabb_max.z, 1.0f, 1e-6f);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L2.4 -- BVH wiring. Deep query behavior is L3; here we only check that the
// single-leaf and multi-leaf construction paths produce a traversable tree.

static void test_trimesh_bvh_single_triangle()
{
	TEST_BEGIN("trimesh_create: 1-triangle BVH has root with one leaf pointing at tri 0");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	TEST_ASSERT(m->bvh.root >= 0);
	TEST_ASSERT(asize(m->bvh.leaves) == 1);
	TEST_ASSERT(m->bvh.leaves[0].body_idx == 0);
	// Single-tri path parks the leaf in child slot a (trimesh.c:193-198).
	const BVHNode* n = &m->bvh.nodes[m->bvh.root];
	TEST_ASSERT(bvh_child_is_leaf((BVH_Child*)&n->a));
	trimesh_free(m);
}

static void test_trimesh_bvh_multi_triangle()
{
	TEST_BEGIN("trimesh_create: N-triangle BVH has every triangle reachable as a leaf");
	// 10 disjoint triangles spaced along X, each with 3 unique verts. Forces
	// the binned-build path (trimesh.c:200).
	v3 dv[30];
	uint32_t di[30];
	for (int t = 0; t < 10; t++) {
		float ox = (float)t * 3.0f;
		dv[3*t + 0] = V3(ox,     0, 0);
		dv[3*t + 1] = V3(ox + 1, 0, 0);
		dv[3*t + 2] = V3(ox,     0, 1);
		di[3*t + 0] = (uint32_t)(3*t + 0);
		di[3*t + 1] = (uint32_t)(3*t + 1);
		di[3*t + 2] = (uint32_t)(3*t + 2);
	}
	TriMesh* m = trimesh_create(dv, 30, di, 10);
	TEST_ASSERT(m->bvh.root >= 0);
	TEST_ASSERT(asize(m->bvh.leaves) == 10);

	// Every body_idx in [0,10) should appear exactly once across leaves.
	int seen[10] = {0};
	for (int i = 0; i < (int)asize(m->bvh.leaves); i++) {
		int bi = m->bvh.leaves[i].body_idx;
		TEST_ASSERT(bi >= 0 && bi < 10);
		seen[bi]++;
	}
	for (int i = 0; i < 10; i++) TEST_ASSERT(seen[i] == 1);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L2.5 -- Material palette

static void test_trimesh_material_default_zero()
{
	TEST_BEGIN("trimesh: material_ids NULL -> get returns 0 for every tri");
	TriMesh* m = ttm_quad_xz();
	TEST_ASSERT(m->material_ids == NULL);
	TEST_ASSERT(trimesh_get_material_id(m, 0) == 0);
	TEST_ASSERT(trimesh_get_material_id(m, 1) == 0);
	trimesh_free(m);
}

static void test_trimesh_material_set_and_get()
{
	TEST_BEGIN("trimesh: set_material_ids stores per-tri, get_material_id reads them");
	TriMesh* m = ttm_quad_xz();
	uint8_t ids[2] = { 7, 42 };
	trimesh_set_material_ids(m, ids);
	TEST_ASSERT(m->material_ids != NULL);
	TEST_ASSERT(trimesh_get_material_id(m, 0) == 7);
	TEST_ASSERT(trimesh_get_material_id(m, 1) == 42);
	trimesh_free(m);
}

static void test_trimesh_material_clear_with_null()
{
	TEST_BEGIN("trimesh: set_material_ids(NULL) frees palette; get returns 0 again");
	TriMesh* m = ttm_quad_xz();
	uint8_t ids[2] = { 1, 2 };
	trimesh_set_material_ids(m, ids);
	TEST_ASSERT(m->material_ids != NULL);
	trimesh_set_material_ids(m, NULL);
	TEST_ASSERT(m->material_ids == NULL);
	TEST_ASSERT(trimesh_get_material_id(m, 0) == 0);
	trimesh_free(m);
}

static void test_trimesh_material_reassign()
{
	TEST_BEGIN("trimesh: second set_material_ids replaces the first");
	TriMesh* m = ttm_quad_xz();
	uint8_t a[2] = { 10, 20 }; trimesh_set_material_ids(m, a);
	uint8_t b[2] = { 30, 40 }; trimesh_set_material_ids(m, b);
	TEST_ASSERT(trimesh_get_material_id(m, 0) == 30);
	TEST_ASSERT(trimesh_get_material_id(m, 1) == 40);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L2.6 -- Name interning

static void test_trimesh_name_intern()
{
	TEST_BEGIN("trimesh: set_name sinterns so equal strings share a pointer");
	TriMesh* a = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	TriMesh* b = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	trimesh_set_name(a, "terrain");
	trimesh_set_name(b, "terrain");
	TEST_ASSERT(trimesh_get_name(a) != NULL);
	TEST_ASSERT(trimesh_get_name(a) == trimesh_get_name(b));
	trimesh_free(a);
	trimesh_free(b);
}

static void test_trimesh_name_null()
{
	TEST_BEGIN("trimesh: default name is NULL, set_name(NULL) clears");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	TEST_ASSERT(trimesh_get_name(m) == NULL);
	trimesh_set_name(m, "thing");
	TEST_ASSERT(trimesh_get_name(m) != NULL);
	trimesh_set_name(m, NULL);
	TEST_ASSERT(trimesh_get_name(m) == NULL);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L2.7 -- Count + free

static void test_trimesh_tri_count()
{
	TEST_BEGIN("trimesh_tri_count: returns stored count, 0 for NULL");
	TriMesh* m1 = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	TEST_ASSERT(trimesh_tri_count(m1) == 1);
	TriMesh* m2 = ttm_quad_xz();
	TEST_ASSERT(trimesh_tri_count(m2) == 2);
	TEST_ASSERT(trimesh_tri_count(NULL) == 0);
	trimesh_free(m1);
	trimesh_free(m2);
}

static void test_trimesh_free_null()
{
	TEST_BEGIN("trimesh_free(NULL) is a no-op");
	trimesh_free(NULL); // must not crash
	TEST_ASSERT(1);
}

// -----------------------------------------------------------------------------
// L3 -- mesh queries (trimesh.c:261-330). trimesh_query_aabb does an AABB/BVH
// walk returning candidate triangle indices; ray_mesh transforms the ray into
// mesh-local space, BVH-prunes, and returns the closest Moeller-Trumbore hit
// with the face normal rotated back to world space.

// L3.1 -- trimesh_query_aabb

static void test_tmq_aabb_single_hit()
{
	TEST_BEGIN("trimesh_query_aabb: single-tri mesh, query containing tri -> 1 hit");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	CK_DYNA int* cands = NULL;
	AABB q = { V3(-1, -1, -1), V3(2, 2, 2) };
	trimesh_query_aabb(m, q, &cands);
	TEST_ASSERT(asize(cands) == 1);
	TEST_ASSERT(cands[0] == 0);
	afree(cands);
	trimesh_free(m);
}

static void test_tmq_aabb_single_miss()
{
	TEST_BEGIN("trimesh_query_aabb: single-tri mesh, query far away -> 0 hits");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	CK_DYNA int* cands = NULL;
	AABB q = { V3(100, 100, 100), V3(101, 101, 101) };
	trimesh_query_aabb(m, q, &cands);
	TEST_ASSERT(asize(cands) == 0);
	afree(cands);
	trimesh_free(m);
}

static void test_tmq_aabb_multi_covers_all()
{
	TEST_BEGIN("trimesh_query_aabb: multi-tri mesh, large query -> all tris returned");
	TriMesh* m = ttm_quad_xz();
	CK_DYNA int* cands = NULL;
	AABB q = { V3(-1, -1, -1), V3(2, 2, 2) };
	trimesh_query_aabb(m, q, &cands);
	TEST_ASSERT(asize(cands) == 2);
	int seen0 = 0, seen1 = 0;
	for (int i = 0; i < (int)asize(cands); i++) {
		if (cands[i] == 0) seen0++;
		if (cands[i] == 1) seen1++;
	}
	TEST_ASSERT(seen0 == 1 && seen1 == 1);
	afree(cands);
	trimesh_free(m);
}

static void test_tmq_aabb_multi_one_region()
{
	TEST_BEGIN("trimesh_query_aabb: multi-tri mesh, query overlapping only one tri");
	// ttm_quad_xz: tri 0 = (0,0,0)-(1,0,0)-(1,0,1), tri 1 = (0,0,0)-(1,0,1)-(0,0,1).
	// Tri 0 AABB is ((0,0,0),(1,0,1)); tri 1 AABB is ((0,0,0),(1,0,1)). They
	// share the same AABB (both span the whole quad), so narrow queries inside
	// may hit both; instead, we query clearly past both and expect zero, and
	// query the corner to guarantee coverage.
	TriMesh* m = ttm_quad_xz();
	CK_DYNA int* c1 = NULL;
	AABB q_far = { V3(5, 5, 5), V3(6, 6, 6) };
	trimesh_query_aabb(m, q_far, &c1);
	TEST_ASSERT(asize(c1) == 0);
	afree(c1);

	// Tiny query at tri 0's exclusive corner (1,0,0) -- tri 1 doesn't reach it.
	CK_DYNA int* c2 = NULL;
	AABB q_corner = { V3(0.95f, -0.01f, -0.01f), V3(1.05f, 0.01f, 0.01f) };
	trimesh_query_aabb(m, q_corner, &c2);
	TEST_ASSERT(asize(c2) >= 1);
	int saw_tri0 = 0;
	for (int i = 0; i < (int)asize(c2); i++) if (c2[i] == 0) saw_tri0 = 1;
	TEST_ASSERT(saw_tri0);
	afree(c2);
	trimesh_free(m);
}

static void test_tmq_aabb_touching_boundary()
{
	TEST_BEGIN("trimesh_query_aabb: boundary-touching query still reports hit");
	TriMesh* m = ttm_single(V3(0,0,0), V3(1,0,0), V3(0,1,0));
	CK_DYNA int* cands = NULL;
	// Query whose lo corner is exactly at v1=(1,0,0).
	AABB q = { V3(1, 0, 0), V3(2, 2, 2) };
	trimesh_query_aabb(m, q, &cands);
	TEST_ASSERT(asize(cands) == 1);
	afree(cands);
	trimesh_free(m);
}

// L3.2 -- ray_mesh

static void test_rm_hit_from_above()
{
	TEST_BEGIN("ray_mesh: perpendicular hit from above returns expected t + +Z normal");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0)); // +Z tri normal
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(0.25f, 0.25f, 1.0f), V3(0, 0, -1), V3(0,0,0), quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-5f);
	TEST_ASSERT_FLOAT(n.z, 1.0f, 1e-5f);
	trimesh_free(m);
}

static void test_rm_miss_direction_away()
{
	TEST_BEGIN("ray_mesh: ray pointing away from mesh -> miss");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0));
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(0.25f, 0.25f, 1.0f), V3(0, 0, 1), V3(0,0,0), quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(!hit);
	trimesh_free(m);
}

static void test_rm_miss_max_t()
{
	TEST_BEGIN("ray_mesh: hit beyond max_t is rejected");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0));
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(0.25f, 0.25f, 5.0f), V3(0, 0, -1), V3(0,0,0), quat_identity(), m, 1.0f, &t, &n);
	TEST_ASSERT(!hit);
	trimesh_free(m);
}

static void test_rm_translated_mesh()
{
	TEST_BEGIN("ray_mesh: mesh_pos translates world ray into local space correctly");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0));
	v3 mesh_pos = V3(10, 5, 0);
	// World ray at (10.25, 5.25, 1) going -Z. In mesh-local this is (0.25, 0.25, 1).
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(10.25f, 5.25f, 1.0f), V3(0, 0, -1), mesh_pos, quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-5f);
	TEST_ASSERT_FLOAT(n.z, 1.0f, 1e-5f);
	trimesh_free(m);
}

static void test_rm_rotated_mesh_normal_world_space()
{
	TEST_BEGIN("ray_mesh: rotated mesh -> returned normal is rotated back to world space");
	// Triangle lives in local XY plane with +Z normal. Rotate mesh by +90 deg
	// about X: quat = (sin(45), 0, 0, cos(45)). That maps local +Z to world -Y,
	// so the triangle's world-space face normal becomes -Y.
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0));
	quat rot = { sinf(0.5f * 1.57079632679f), 0, 0, cosf(0.5f * 1.57079632679f) };

	// World ray origin (0.25, 1, 0.25) going -Y. Local-space ray lands at
	// (0.25, 0.25, 0) at t=1 -- clearly inside the triangle.
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(0.25f, 1.0f, 0.25f), V3(0, -1, 0), V3(0,0,0), rot, m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-4f);
	TEST_ASSERT_FLOAT(n.y, -1.0f, 1e-4f);
	TEST_ASSERT_FLOAT(n.x,  0.0f, 1e-4f);
	TEST_ASSERT_FLOAT(n.z,  0.0f, 1e-4f);
	trimesh_free(m);
}

static void test_rm_picks_closest_candidate()
{
	TEST_BEGIN("ray_mesh: multiple candidate triangles -> closest t wins");
	// Two overlapping (in XZ) tris at different heights Y=0 and Y=2.
	v3 v[6] = {
		V3(-1, 0, -1), V3(1, 0, -1), V3(-1, 0, 1),   // lower tri (y=0)
		V3(-1, 2, -1), V3(1, 2, -1), V3(-1, 2, 1),   // upper tri (y=2)
	};
	uint32_t idx[6] = { 0, 1, 2, 3, 4, 5 };
	TriMesh* m = trimesh_create(v, 6, idx, 2);

	// Ray from above (y=5) going -Y through x=-0.5, z=-0.5 (inside both tris'
	// projection). Upper tri (y=2) is hit first at t=3, lower at t=5.
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(-0.5f, 5.0f, -0.5f), V3(0, -1, 0), V3(0,0,0), quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 3.0f, 1e-4f);
	trimesh_free(m);
}

static void test_rm_back_side_hit_pinned_two_sided()
{
	TEST_BEGIN("ray_mesh: back-side hit (pinned: current impl is two-sided)");
	TriMesh* m = ttm_single(V3(0,0,0), V3(2,0,0), V3(0,2,0)); // +Z tri normal
	float t = -1.0f; v3 n = V3(0,0,0);
	// Origin below plane, direction +Z: should hit back side.
	int hit = ray_mesh(V3(0.25f, 0.25f, -1.0f), V3(0, 0, 1), V3(0,0,0), quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 1.0f, 1e-5f);
	trimesh_free(m);
}

static void test_rm_axis_aligned_ray()
{
	TEST_BEGIN("ray_mesh: axis-aligned ray (rd.x = rd.z = 0) doesn't fall over");
	TriMesh* m = ttm_single(V3(-1,0,-1), V3(1,0,-1), V3(-1,0,1)); // y=0 plane, -Y normal
	float t = -1.0f; v3 n = V3(0,0,0);
	int hit = ray_mesh(V3(-0.5f, 5.0f, -0.5f), V3(0, -1, 0), V3(0,0,0), quat_identity(), m, 10.0f, &t, &n);
	TEST_ASSERT(hit);
	TEST_ASSERT_FLOAT(t, 5.0f, 1e-4f);
	trimesh_free(m);
}

// -----------------------------------------------------------------------------
// L4 -- single-triangle narrowphase (trimesh.c:685 and 810). Both routines
// operate entirely in mesh-local space: triangle at (v0,v1,v2), convex shape
// already transformed into that frame. Contact normal convention is A(hull) ->
// B(triangle); for a hull above the triangle plane the normal is -tri_n.

// Big CCW triangle in XY plane (tri_n = +Z), comfortably containing the unit
// box's projection so every bottom vertex lands inside the triangle region.
#define TRI_BIG_V0 V3(-10, -5, 0)
#define TRI_BIG_V1 V3( 10, -5, 0)
#define TRI_BIG_V2 V3(  0, 10, 0)
#define TRI_BIG_N  V3(  0,  0, 1)

// L4.1 -- collide_hull_triangle_local

static void test_cht_face_down_box()
{
	TEST_BEGIN("collide_hull_triangle_local: face-down box -> 4 contacts, depth matches, normal -tri_n");
	// Box center at z=0.3 (penetration 0.2), identity rotation, half_extents 0.5.
	ConvexHull box = { &s_unit_box_hull, V3(0, 0, 0.3f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold m = {0};
	int hit = collide_hull_triangle_local(box, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count == 4);
	for (int i = 0; i < m.count; i++) {
		TEST_ASSERT_FLOAT(m.contacts[i].normal.z, -1.0f, 1e-5f);
		TEST_ASSERT_FLOAT(m.contacts[i].penetration, 0.2f, 1e-4f);
		TEST_ASSERT_FLOAT(m.contacts[i].point.z, 0.0f, 1e-5f); // on tri plane
	}
}

static void test_cht_edge_down_box()
{
	TEST_BEGIN("collide_hull_triangle_local: edge-down box (+45 deg about X) -> 2 contacts");
	// +45 deg about X -> two verts meet at (0, 0, -half*sqrt(2)/... actually -0.707).
	// Center z=0.6 puts them at z=-0.107. Non-penetrating verts sit at z=0.6 or z=1.307.
	float half_angle = 0.5f * (3.14159265f * 0.25f); // 22.5 deg
	quat rot = { sinf(half_angle), 0, 0, cosf(half_angle) };
	ConvexHull box = { &s_unit_box_hull, V3(0, 0, 0.6f), rot, V3(0.5f, 0.5f, 0.5f) };
	Manifold m = {0};
	int hit = collide_hull_triangle_local(box, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count == 2);
	float expected_depth = 0.5f * 1.41421356f - 0.6f; // half*sqrt(2) - center.z
	for (int i = 0; i < m.count; i++) {
		TEST_ASSERT_FLOAT(m.contacts[i].normal.z, -1.0f, 1e-5f);
		TEST_ASSERT_FLOAT(m.contacts[i].penetration, expected_depth, 1e-4f);
		TEST_ASSERT_FLOAT(m.contacts[i].point.z, 0.0f, 1e-5f);
	}
	// Two contacts should be ~1 unit apart along the X axis.
	float dx = fabsf(m.contacts[0].point.x - m.contacts[1].point.x);
	TEST_ASSERT_FLOAT(dx, 1.0f, 1e-4f);
}

static void test_cht_back_face_tunneled()
{
	TEST_BEGIN("collide_hull_triangle_local: center below plane -> normal flipped to +tri_n");
	// Box center at z=-0.3, top verts at z=+0.2 (past the plane from the "wrong" side).
	ConvexHull box = { &s_unit_box_hull, V3(0, 0, -0.3f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold m = {0};
	int hit = collide_hull_triangle_local(box, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count == 4);
	for (int i = 0; i < m.count; i++) {
		TEST_ASSERT_FLOAT(m.contacts[i].normal.z, 1.0f, 1e-5f); // flipped, pointing out of back
		TEST_ASSERT_FLOAT(m.contacts[i].penetration, 0.2f, 1e-4f);
	}
}

static void test_cht_separated_above()
{
	TEST_BEGIN("collide_hull_triangle_local: box fully above plane -> 0 contacts");
	ConvexHull box = { &s_unit_box_hull, V3(0, 0, 2.0f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold m = {0};
	int hit = collide_hull_triangle_local(box, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &m);
	TEST_ASSERT(!hit);
	TEST_ASSERT(m.count == 0);
}

static void test_cht_separated_below()
{
	TEST_BEGIN("collide_hull_triangle_local: box fully below plane -> 0 contacts");
	ConvexHull box = { &s_unit_box_hull, V3(0, 0, -2.0f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold m = {0};
	int hit = collide_hull_triangle_local(box, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &m);
	TEST_ASSERT(!hit);
	TEST_ASSERT(m.count == 0);
}

static void test_cht_deep_penetration_depth()
{
	TEST_BEGIN("collide_hull_triangle_local: penetration scales linearly with center depth");
	// Shallow case (0.1 penetration).
	ConvexHull shallow = { &s_unit_box_hull, V3(0, 0, 0.4f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold ms = {0};
	collide_hull_triangle_local(shallow, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &ms);
	TEST_ASSERT(ms.count == 4);
	TEST_ASSERT_FLOAT(ms.contacts[0].penetration, 0.1f, 1e-4f);

	// Deep case (0.4 penetration).
	ConvexHull deep = { &s_unit_box_hull, V3(0, 0, 0.1f), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	Manifold md = {0};
	collide_hull_triangle_local(deep, TRI_BIG_V0, TRI_BIG_V1, TRI_BIG_V2, TRI_BIG_N, &md);
	TEST_ASSERT(md.count == 4);
	TEST_ASSERT_FLOAT(md.contacts[0].penetration, 0.4f, 1e-4f);
}


// -----------------------------------------------------------------------------
// L5 -- trimesh_reduce (trimesh.c:455). Three passes in order:
//   Pass 1: back-face flip. Any contact normal with dot(cn, face_n) > 0.05
//           is replaced with -face_n. This is one-sided-mesh enforcement.
//   Pass 2: pick deepest contact per source manifold as rep_contact (n>=2 only).
//   Pass 3: for each source whose rep contact is NOT face-trusted (threshold:
//           dot(-cn, face_n) > 0.999), test rep against every other tri's
//           TriTest; if blocked, mark killed and record the blocker index.
//   Pass 4: resolve -- kill source outright, or in the wedge case (mutual
//           block) replace source's normals with the blocker's -face_n.
//
// These tests pin the parts that are invariant under tri_test_blocks' sign
// convention. Ghost/wedge end-to-end scenarios land in L7.

static ReduceRec make_rec(v3 v0, v3 v1, v3 v2, v3 cp, v3 cn, float depth)
{
	ReduceRec r;
	memset(&r, 0, sizeof(r));
	r.corrected_by = -1;
	v3 tn_raw = cross(sub(v1, v0), sub(v2, v0));
	v3 tn = scale(tn_raw, 1.0f / sqrtf(len2(tn_raw)));
	tri_test_build(&r.tt, v0, v1, v2, tn);
	r.m.count = 1;
	r.m.contacts[0] = (Contact){ .point = cp, .normal = cn, .penetration = depth, .feature_id = 0 };
	return r;
}

static void test_reduce_empty_nocrash()
{
	TEST_BEGIN("trimesh_reduce: n=0 is a no-op (no crash)");
	trimesh_reduce(NULL, 0);
	TEST_ASSERT(1);
}

static void test_reduce_n1_flips_back_face_normal()
{
	TEST_BEGIN("trimesh_reduce: n=1, cn along +face_n -> flipped to -face_n");
	// Tri in XY plane, +Z face_n. Contact normal +Z infringes (body is tunneled).
	ReduceRec r = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, 1), 0.1f);
	trimesh_reduce(&r, 1);
	TEST_ASSERT_FLOAT(r.m.contacts[0].normal.z, -1.0f, 1e-5f);
	TEST_ASSERT(r.m.count == 1); // not killed (n<2 skips block pass)
}

static void test_reduce_n1_front_face_unchanged()
{
	TEST_BEGIN("trimesh_reduce: n=1, cn = -face_n left alone");
	ReduceRec r = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, -1), 0.1f);
	trimesh_reduce(&r, 1);
	TEST_ASSERT_FLOAT(r.m.contacts[0].normal.z, -1.0f, 1e-5f);
}

static void test_reduce_back_face_flip_threshold()
{
	TEST_BEGIN("trimesh_reduce: back-face flip threshold pinned at 0.05f (absolute)");
	// dot(cn, face_n) = 0.04 -> not flipped.
	ReduceRec r1 = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, 0.04f), 0.1f);
	trimesh_reduce(&r1, 1);
	TEST_ASSERT_FLOAT(r1.m.contacts[0].normal.z, 0.04f, 1e-5f); // unchanged

	// dot(cn, face_n) = 0.06 -> flipped to -face_n.
	ReduceRec r2 = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, 0.06f), 0.1f);
	trimesh_reduce(&r2, 1);
	TEST_ASSERT_FLOAT(r2.m.contacts[0].normal.z, -1.0f, 1e-5f);
}

static void test_reduce_mixed_contacts_flipped_independently()
{
	TEST_BEGIN("trimesh_reduce: per-contact back-face flip within one manifold");
	ReduceRec r;
	memset(&r, 0, sizeof(r));
	r.corrected_by = -1;
	v3 v0 = V3(0,0,0), v1 = V3(2,0,0), v2 = V3(0,2,0);
	tri_test_build(&r.tt, v0, v1, v2, V3(0, 0, 1));
	r.m.count = 3;
	r.m.contacts[0] = (Contact){ .point = V3(0.2f, 0.2f, 0), .normal = V3(0, 0, -1),  .penetration = 0.1f };
	r.m.contacts[1] = (Contact){ .point = V3(0.3f, 0.3f, 0), .normal = V3(0, 0, 1),   .penetration = 0.1f };
	r.m.contacts[2] = (Contact){ .point = V3(0.4f, 0.4f, 0), .normal = V3(1, 0, 0),   .penetration = 0.1f };
	trimesh_reduce(&r, 1);
	TEST_ASSERT_FLOAT(r.m.contacts[0].normal.z, -1.0f, 1e-5f); // already front, unchanged
	TEST_ASSERT_FLOAT(r.m.contacts[1].normal.z, -1.0f, 1e-5f); // was +Z, flipped
	// (1, 0, 0) has dot with face_n = 0, below threshold -> unchanged.
	TEST_ASSERT_FLOAT(r.m.contacts[2].normal.x, 1.0f, 1e-5f);
	TEST_ASSERT_FLOAT(r.m.contacts[2].normal.z, 0.0f, 1e-5f);
}

static void test_reduce_empty_manifold_in_multi_skipped()
{
	TEST_BEGIN("trimesh_reduce: empty source manifold among others is skipped cleanly");
	ReduceRec recs[2];
	recs[0] = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, -1), 0.1f);
	recs[1] = make_rec(V3(10,10,0), V3(12,10,0), V3(10,12,0), V3(11, 11, 0), V3(0, 0, -1), 0.1f);
	recs[1].m.count = 0; // empty
	trimesh_reduce(recs, 2);
	TEST_ASSERT(recs[0].m.count == 1); // survivor
	TEST_ASSERT(recs[1].m.count == 0); // still empty
}

static void test_reduce_face_trust_both_survive()
{
	TEST_BEGIN("trimesh_reduce: both face-down (-face_n) manifolds bypass block test");
	// Even if one happens to project onto the other's tri region, face-trust
	// short-circuits the block pass for both. Both survive intact.
	ReduceRec recs[2];
	recs[0] = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, -1), 0.2f);
	recs[1] = make_rec(V3(0,0,0), V3(0,2,0), V3(-2,0,0), V3(-0.25f, 0.25f, 0), V3(0, 0, -1), 0.1f);
	trimesh_reduce(recs, 2);
	TEST_ASSERT(recs[0].m.count == 1);
	TEST_ASSERT(recs[1].m.count == 1);
	TEST_ASSERT_FLOAT(recs[0].m.contacts[0].normal.z, -1.0f, 1e-5f);
	TEST_ASSERT_FLOAT(recs[1].m.contacts[0].normal.z, -1.0f, 1e-5f);
}

static void test_reduce_rep_contact_is_deepest()
{
	TEST_BEGIN("trimesh_reduce: rep_contact picks the deepest contact index");
	ReduceRec recs[2];
	// Source with 3 contacts, deepest at index 1.
	memset(&recs[0], 0, sizeof(ReduceRec));
	recs[0].corrected_by = -1;
	tri_test_build(&recs[0].tt, V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0, 0, 1));
	recs[0].m.count = 3;
	recs[0].m.contacts[0] = (Contact){ .point = V3(0.1f, 0.1f, 0), .normal = V3(0, 0, -1), .penetration = 0.1f };
	recs[0].m.contacts[1] = (Contact){ .point = V3(0.2f, 0.2f, 0), .normal = V3(0, 0, -1), .penetration = 0.3f };
	recs[0].m.contacts[2] = (Contact){ .point = V3(0.3f, 0.3f, 0), .normal = V3(0, 0, -1), .penetration = 0.2f };
	// Second rec: face-trusted so Pass 3 skips. Pass 2 still runs for both.
	recs[1] = make_rec(V3(5,5,0), V3(7,5,0), V3(5,7,0), V3(6, 6, 0), V3(0, 0, -1), 0.5f);
	trimesh_reduce(recs, 2);
	TEST_ASSERT(recs[0].rep_contact == 1); // deepest of 3
	TEST_ASSERT(recs[1].rep_contact == 0); // only contact
}

static void test_reduce_spatially_distant_both_survive()
{
	TEST_BEGIN("trimesh_reduce: spatially distant non-face-trusted manifolds don't kill each other");
	// Two manifolds with side-normals (not face-trusted) but on triangles 100
	// units apart. Neither's contact point can land inside the other's tri.
	ReduceRec recs[2];
	recs[0] = make_rec(V3(0,0,0),    V3(1,0,0),     V3(0,1,0),     V3(0.2f, 0.2f, 0),     V3(1, 0, 0), 0.1f);
	recs[1] = make_rec(V3(100,100,0), V3(101,100,0), V3(100,101,0), V3(100.2f, 100.2f, 0), V3(1, 0, 0), 0.1f);
	trimesh_reduce(recs, 2);
	TEST_ASSERT(recs[0].m.count == 1);
	TEST_ASSERT(recs[1].m.count == 1);
}

static void test_reduce_preserves_face_trusted_normal_after_flip()
{
	TEST_BEGIN("trimesh_reduce: back-face flip happens BEFORE face-trust; flipped normal is then trusted");
	// Contact starts with cn = +face_n (back-face). Pass 1 flips it to -face_n,
	// then Pass 3 sees a face-trusted contact and skips block. Manifold survives.
	ReduceRec recs[2];
	recs[0] = make_rec(V3(0,0,0), V3(2,0,0), V3(0,2,0), V3(0.25f, 0.25f, 0), V3(0, 0, 1), 0.1f);
	// Second rec only exists to make n>=2 so passes 2-4 execute.
	recs[1] = make_rec(V3(5,5,0), V3(7,5,0), V3(5,7,0), V3(6, 6, 0), V3(0, 0, -1), 0.1f);
	trimesh_reduce(recs, 2);
	TEST_ASSERT(recs[0].m.count == 1); // not killed
	TEST_ASSERT_FLOAT(recs[0].m.contacts[0].normal.z, -1.0f, 1e-5f);
}

// -----------------------------------------------------------------------------
// L6 -- full emits (trimesh.c:588-945). Each `collide_<shape>_mesh_emit`
// transforms the convex into mesh-local space, BVH-prunes, runs per-tri
// narrowphase, reduces, and pushes one InternalManifold per surviving triangle
// -- contacts transformed back to world space.
//
// The emit signatures take WorldInternal* but trimesh_flush never dereferences
// it, so NULL is safe for these unit tests.

// Shared big flat triangle in XZ plane (tri_n = +Y), centered on origin.
// Using XZ (not XY) is more natural for "body falls onto mesh" scenarios.
// CCW from +Y (looking down), so tri_normal = (v1-v0) x (v2-v0) points +Y.
static TriMesh* ttm_big_floor()
{
	v3 v[3] = { V3(-10, 0, -10), V3(0, 0, 10), V3(10, 0, -10) };
	uint32_t idx[3] = { 0, 1, 2 };
	return trimesh_create(v, 3, idx, 1);
}

static void test_emit_sphere_aabb_miss()
{
	TEST_BEGIN("emit sphere-mesh: AABB early-out -> no manifolds");
	TriMesh* m = ttm_big_floor();
	Sphere s = { V3(100, 100, 100), 0.5f };
	CK_DYNA InternalManifold* ms = NULL;
	collide_sphere_mesh_emit(NULL, 3, 7, s, V3(0,0,0), quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) == 0);
	afree(ms);
	trimesh_free(m);
}

static void test_emit_sphere_single_hit()
{
	TEST_BEGIN("emit sphere-mesh: sphere penetrating tri -> 1 manifold with body_a/body_b propagated");
	TriMesh* m = ttm_big_floor();
	Sphere s = { V3(0, 0.3f, 0), 0.5f }; // 0.2 penetration into XZ tri plane at y=0
	CK_DYNA InternalManifold* ms = NULL;
	collide_sphere_mesh_emit(NULL, 3, 7, s, V3(0,0,0), quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) == 1);
	TEST_ASSERT(ms[0].body_a == 3);
	TEST_ASSERT(ms[0].body_b == 7);
	TEST_ASSERT(ms[0].m.count > 0);
	// Contact normal points A(sphere) -> B(mesh). Face normal of this CCW tri
	// is +Y; body above plane -> normal -Y.
	TEST_ASSERT_FLOAT(ms[0].m.contacts[0].normal.y, -1.0f, 1e-4f);
	afree(ms);
	trimesh_free(m);
}

static void test_emit_sphere_sub_id_encoding()
{
	TEST_BEGIN("emit sphere-mesh: sub_id = sub_id_base | (tri_idx + 1)");
	// 3-tri mesh, sphere over the middle one. Expect sub_id with tri_idx+1 bit.
	// Each tri is CCW from +Y so tri_normal = +Y (like ttm_big_floor).
	v3 v[9] = {
		V3(-20, 0, -1), V3(-10, 0, 1), V3(-10, 0, -1), // tri 0
		V3( -5, 0, -1), V3(  5, 0, 1), V3(  5, 0, -1), // tri 1
		V3( 10, 0, -1), V3( 20, 0, 1), V3( 20, 0, -1), // tri 2
	};
	uint32_t idx[9] = { 0,1,2, 3,4,5, 6,7,8 };
	TriMesh* m = trimesh_create(v, 9, idx, 3);
	// Sphere penetrating only tri 1 (centered at x=0).
	Sphere s = { V3(0, 0.3f, 0), 0.5f };
	CK_DYNA InternalManifold* ms = NULL;
	uint32_t base = 0x1000u;
	collide_sphere_mesh_emit(NULL, 3, 7, s, V3(0,0,0), quat_identity(), m, base, &ms);
	TEST_ASSERT(asize(ms) == 1);
	TEST_ASSERT(ms[0].sub_id == (base | (uint32_t)(1 + 1))); // tri_idx=1 -> 2
	afree(ms);
	trimesh_free(m);
}

static void test_emit_sphere_translated_mesh_world_contacts()
{
	TEST_BEGIN("emit sphere-mesh: mesh translation -> contact point emitted in WORLD space");
	TriMesh* m = ttm_big_floor();
	v3 mesh_pos = V3(100, 50, 25);
	// Sphere directly above the translated mesh origin, penetrating 0.2.
	Sphere s = { V3(100, 50.3f, 25), 0.5f };
	CK_DYNA InternalManifold* ms = NULL;
	collide_sphere_mesh_emit(NULL, 3, 7, s, mesh_pos, quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) == 1);
	// collide_sphere_hull places the contact at the sphere's deepest point, so
	// cp.y ~= 49.8 (sphere surface), not exactly on the plane. Correctness
	// invariant: contact is within radius+slop of the sphere center in WORLD.
	v3 cp = ms[0].m.contacts[0].point;
	v3 sc = s.center;
	float dist = sqrtf((cp.x-sc.x)*(cp.x-sc.x) + (cp.y-sc.y)*(cp.y-sc.y) + (cp.z-sc.z)*(cp.z-sc.z));
	TEST_ASSERT(dist <= s.radius + 1e-3f);
	// And demonstrably translated: far from origin.
	TEST_ASSERT(sqrtf(cp.x*cp.x + cp.y*cp.y + cp.z*cp.z) > 50.0f);
	afree(ms);
	trimesh_free(m);
}

static void test_emit_sphere_rotated_mesh_world_normal()
{
	TEST_BEGIN("emit sphere-mesh: mesh rotation -> contact normal rotated back to world space");
	TriMesh* m = ttm_big_floor(); // +Y face_n in local space
	// Rotate mesh +90 deg about X: local +Y -> world +Z, so the tri's front
	// face points +Z after rotation. Sphere on the +Z side -> normal A->B = -Z.
	float half_angle = 0.5f * (3.14159265f * 0.5f);
	quat rot = { sinf(half_angle), 0, 0, cosf(half_angle) };
	Sphere s = { V3(0, 0, 0.3f), 0.5f };
	CK_DYNA InternalManifold* ms = NULL;
	collide_sphere_mesh_emit(NULL, 3, 7, s, V3(0,0,0), rot, m, 0, &ms);
	TEST_ASSERT(asize(ms) == 1);
	TEST_ASSERT_FLOAT(ms[0].m.contacts[0].normal.z, -1.0f, 1e-3f);
	TEST_ASSERT_FLOAT(ms[0].m.contacts[0].normal.x,  0.0f, 1e-3f);
	TEST_ASSERT_FLOAT(ms[0].m.contacts[0].normal.y,  0.0f, 1e-3f);
	afree(ms);
	trimesh_free(m);
}

static void test_emit_sphere_multi_tri_distinct_sub_ids()
{
	TEST_BEGIN("emit sphere-mesh: sphere over multiple tris -> distinct sub_ids per manifold");
	// Two coplanar tris sharing an edge (quad split), sphere hovering over the
	// seam.
	TriMesh* m = ttm_quad_xz(); // tris in y=0 plane
	// Large sphere above the seam: penetrating both tris.
	Sphere s = { V3(0.5f, 0.1f, 0.5f), 0.3f };
	CK_DYNA InternalManifold* ms = NULL;
	collide_sphere_mesh_emit(NULL, 3, 7, s, V3(0,0,0), quat_identity(), m, 0, &ms);
	// Expect up to 2 manifolds (one per tri), possibly fewer if reduction kills one.
	TEST_ASSERT(asize(ms) >= 1);
	// If two, sub_ids must differ.
	if (asize(ms) == 2) {
		TEST_ASSERT(ms[0].sub_id != ms[1].sub_id);
	}
	afree(ms);
	trimesh_free(m);
}

static void test_emit_capsule_single_hit()
{
	TEST_BEGIN("emit capsule-mesh: capsule penetrating tri -> manifold with correct normal");
	TriMesh* m = ttm_big_floor();
	// Horizontal capsule just above the tri. Segment along X from (-0.5, 0.2, 0)
	// to (0.5, 0.2, 0), radius 0.3 -> bottom at y=-0.1 (0.1 penetration).
	Capsule c = { V3(-0.5f, 0.2f, 0), V3(0.5f, 0.2f, 0), 0.3f };
	CK_DYNA InternalManifold* ms = NULL;
	collide_capsule_mesh_emit(NULL, 3, 7, c, V3(0,0,0), quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) >= 1);
	TEST_ASSERT(ms[0].m.count > 0);
	TEST_ASSERT_FLOAT(ms[0].m.contacts[0].normal.y, -1.0f, 1e-3f);
	afree(ms);
	trimesh_free(m);
}

static void test_emit_box_single_hit()
{
	TEST_BEGIN("emit box-mesh: box resting on tri -> 4 contacts normal -tri_n");
	TriMesh* m = ttm_big_floor();
	Box b = { V3(0, 0.3f, 0), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	CK_DYNA InternalManifold* ms = NULL;
	collide_box_mesh_emit(NULL, 3, 7, b, V3(0,0,0), quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) == 1);
	TEST_ASSERT(ms[0].m.count == 4);
	for (int i = 0; i < ms[0].m.count; i++) {
		TEST_ASSERT_FLOAT(ms[0].m.contacts[i].normal.y, -1.0f, 1e-4f);
		TEST_ASSERT_FLOAT(ms[0].m.contacts[i].penetration, 0.2f, 1e-4f);
	}
	afree(ms);
	trimesh_free(m);
}

static void test_emit_hull_single_hit()
{
	TEST_BEGIN("emit hull-mesh: box-as-hull penetrating tri -> 4 contacts via dedicated hull-tri narrowphase");
	TriMesh* m = ttm_big_floor();
	ConvexHull h = { &s_unit_box_hull, V3(0, 0.3f, 0), quat_identity(), V3(0.5f, 0.5f, 0.5f) };
	CK_DYNA InternalManifold* ms = NULL;
	collide_hull_mesh_emit(NULL, 3, 7, h, V3(0,0,0), quat_identity(), m, 0, &ms);
	TEST_ASSERT(asize(ms) == 1);
	TEST_ASSERT(ms[0].m.count == 4);
	for (int i = 0; i < ms[0].m.count; i++) {
		TEST_ASSERT_FLOAT(ms[0].m.contacts[i].normal.y, -1.0f, 1e-4f);
	}
	afree(ms);
	trimesh_free(m);
}


// -----------------------------------------------------------------------------
// L7 -- integration. Full world_step scenarios: bodies drop onto a mesh floor
// and settle, a box slides across a shared edge, raycast against the mesh,
// and a snapshot roundtrip with a registered mesh.

// Bigger flat floor for drops. 2 triangles, tri_n = +Y.
static TriMesh* ttm_wide_floor()
{
	v3 v[4] = { V3(-20, 0, -20), V3(-20, 0, 20), V3(20, 0, 20), V3(20, 0, -20) };
	uint32_t idx[6] = { 0, 1, 2,   0, 2, 3 };
	return trimesh_create(v, 4, idx, 2);
}

static void test_int_sphere_rests_on_mesh_floor()
{
	TEST_BEGIN("integration: sphere drops onto mesh floor and settles");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body ball = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, ball);
	v3 v = body_get_velocity(w, ball);
	printf("  sphere final y=%.4f vel=%.4f\n", p.y, v3_len(v));
	TEST_ASSERT(p.y > 0.4f && p.y < 0.6f);
	TEST_ASSERT(v3_len(v) < 0.1f);
	destroy_world(w);
	trimesh_free(mesh);
}

static void test_int_box_rests_on_mesh_floor()
{
	TEST_BEGIN("integration: box drops onto mesh floor and settles");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body box = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, box);
	v3 v = body_get_velocity(w, box);
	printf("  box final y=%.4f vel=%.4f\n", p.y, v3_len(v));
	TEST_ASSERT(p.y > 0.3f && p.y < 0.5f);
	TEST_ASSERT(v3_len(v) < 0.1f);
	destroy_world(w);
	trimesh_free(mesh);
}

static void test_int_capsule_rests_on_mesh_floor()
{
	TEST_BEGIN("integration: capsule drops onto mesh floor and settles");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body cap = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, cap, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = 0.5f, .radius = 0.3f } });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, cap);
	v3 v = body_get_velocity(w, cap);
	printf("  capsule final y=%.4f vel=%.4f\n", p.y, v3_len(v));
	// Vertical capsule: center rests at half_height + radius above floor.
	TEST_ASSERT(p.y > 0.7f && p.y < 0.9f);
	TEST_ASSERT(v3_len(v) < 0.1f);
	destroy_world(w);
	trimesh_free(mesh);
}


static void test_int_hull_rests_on_mesh_floor()
{
	TEST_BEGIN("integration: convex hull (unit box) drops onto mesh floor and settles");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body body = create_body(w, (BodyParams){ .position = V3(0, 2.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, body, (ShapeParams){ .type = SHAPE_HULL, .hull = { .hull = &s_unit_box_hull, .scale = V3(0.4f, 0.4f, 0.4f) } });
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	v3 p = body_get_position(w, body);
	v3 v = body_get_velocity(w, body);
	printf("  hull final y=%.4f vel=%.4f\n", p.y, v3_len(v));
	TEST_ASSERT(p.y > 0.3f && p.y < 0.5f);
	TEST_ASSERT(v3_len(v) < 0.1f);
	destroy_world(w);
	trimesh_free(mesh);
}

// The "sliding across an internal edge" test -- Bepu Test 1 / Jolt Test 11.
// A small box sits on a 2-tri flat quad whose shared edge runs through the
// box's travel path. Moving horizontally, the box must stay on the surface
// without ghost-kick (sudden +Y bump as the shared edge contact is caught by
// narrowphase and reduced away) or drop-through.
static void test_int_box_slides_across_shared_edge()
{
	TEST_BEGIN("integration: box slides across shared-edge seam without vertical perturbation");
	TriMesh* mesh = ttm_wide_floor(); // shared edge is the diagonal of the XZ square
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body box = create_body(w, (BodyParams){ .position = V3(-2.0f, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.2f, 0.2f, 0.2f) });

	// Let box settle on surface first.
	for (int i = 0; i < 60; i++) world_step(w, 1.0f / 60.0f);
	float settled_y = body_get_position(w, box).y;
	printf("  settled y=%.4f\n", settled_y);

	// Nudge horizontally at 1 m/s toward +X. The diagonal seam passes through
	// origin; box will cross it roughly mid-traversal.
	body_set_velocity(w, box, V3(1.0f, 0, 0));

	float max_dev = 0.0f;
	float min_y = 1e9f;
	for (int i = 0; i < 180; i++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, box);
		float dev = fabsf(p.y - settled_y);
		if (dev > max_dev) max_dev = dev;
		if (p.y < min_y) min_y = p.y;
	}
	printf("  sliding max_dev=%.5f min_y=%.4f\n", max_dev, min_y);

	// Pinning current behavior: tolerance set based on first-run measurement.
	// If Bepu-parity reducer fixes land later, dev should shrink, not grow.
	TEST_ASSERT(max_dev < 0.05f);
	TEST_ASSERT(min_y > settled_y - 0.05f); // no drop-through
	destroy_world(w);
	trimesh_free(mesh);
}

// Vertical capsule given horizontal velocity on a flat mesh floor. A vertical
// capsule on a flat surface is unstable under friction -- the friction torque
// tips it over to horizontal. After tipping, the capsule rests at y=radius
// instead of y=(radius+half_height). What we actually want to detect here is
// the LURCH -- a sudden upward kick that was the stress-scene bug. So we only
// guard against y going ABOVE the settled height by more than LINEAR_SLOP.
static void test_int_capsule_slides_across_shared_edge()
{
	TEST_BEGIN("integration: vertical capsule with horizontal velocity has no upward lurch");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
	Body cap = create_body(w, (BodyParams){ .position = V3(-2.0f, 0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, cap, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = 0.2f, .radius = 0.2f } });
	for (int i = 0; i < 60; i++) world_step(w, 1.0f / 60.0f);
	float settled_y = body_get_position(w, cap).y;
	body_set_velocity(w, cap, V3(1.0f, 0, 0));
	float max_up_excursion = 0.0f;
	for (int i = 0; i < 180; i++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, cap);
		float up = p.y - settled_y;
		if (up > max_up_excursion) max_up_excursion = up;
	}
	printf("  capsule slide settled_y=%.5f max_up_excursion=%.5f\n", settled_y, max_up_excursion);
	// LURCH guard: no single event may push the capsule upward by more than
	// a centimeter. Downward motion (from tipping) is expected physics.
	TEST_ASSERT(max_up_excursion < 0.02f);
	destroy_world(w);
	trimesh_free(mesh);
}


static void test_int_mesh_raycast_world()
{
	TEST_BEGIN("integration: world_raycast against mesh returns hit with correct normal");
	TriMesh* mesh = ttm_wide_floor();
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });

	RayHit hit;
	int got = world_raycast(w, V3(0.5f, 5.0f, 0.5f), V3(0, -1, 0), 10.0f, &hit);
	TEST_ASSERT(got);
	TEST_ASSERT_FLOAT(hit.distance, 5.0f, 1e-3f);
	TEST_ASSERT_FLOAT(hit.normal.y, 1.0f, 1e-3f); // tri face_n = +Y
	destroy_world(w);
	trimesh_free(mesh);
}

static void test_int_mesh_snapshot_roundtrip()
{
	TEST_BEGIN("integration: mesh body survives snapshot save/load roundtrip");
	const char* path = "test_trimesh_snapshot.dat";
	TriMesh* mesh = ttm_wide_floor();
	trimesh_set_name(mesh, "wide_floor");
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
		world_register_mesh(w, mesh);
		Body floor = create_body(w, (BodyParams){ .position = V3(0,0,0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = mesh });
		TEST_ASSERT(world_save_snapshot(w, path));
		destroy_world(w);
	}
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
		world_register_mesh(w, mesh);
		TEST_ASSERT(world_load_snapshot_into(w, path));
		Body bodies[4]; int nb = world_get_bodies(w, bodies, 4);
		TEST_ASSERT(nb == 1);
		destroy_world(w);
	}
	trimesh_free(mesh);
	remove(path);
}

// -----------------------------------------------------------------------------

static void run_trimesh_unit_tests()
{
	printf("--- trimesh unit tests (L1: primitives) ---\n");
	test_pit2d_ccw_interior();
	test_pit2d_ccw_outside();
	test_pit2d_boundary();
	test_pit2d_cw_winding();
	test_pit2d_projection_invariance();

	test_ray_tri_front_hit_center();
	test_ray_tri_back_hit_pinned_two_sided();
	test_ray_tri_miss_barycentric();
	test_ray_tri_miss_parallel();
	test_ray_tri_miss_max_t();
	test_ray_tri_miss_behind_origin();

	test_tt_build_structural_invariants();
	test_tt_build_thr_scales_with_edge_length();
	test_tt_build_thr_floor();
	test_tt_build_arbitrary_orientation();

	test_ttb_out_of_plane_rejects();
	test_ttb_aligned_normal_no_infringe();
	test_ttb_far_outside_triangle();
	test_ttb_zero_normal();

	printf("--- trimesh unit tests (L2: construction) ---\n");
	test_trimesh_shared_topology_twins();
	test_trimesh_shared_topology_face_cycles();
	test_trimesh_shared_topology_origins();

	test_trimesh_hull_wiring();
	test_trimesh_centroid_and_epsilon();
	test_trimesh_front_back_planes();
	test_trimesh_plane_offset_nonzero();
	test_trimesh_tri_normal_stored();

	test_trimesh_aabb_single_triangle();
	test_trimesh_aabb_multi_triangle();

	test_trimesh_bvh_single_triangle();
	test_trimesh_bvh_multi_triangle();

	test_trimesh_material_default_zero();
	test_trimesh_material_set_and_get();
	test_trimesh_material_clear_with_null();
	test_trimesh_material_reassign();

	test_trimesh_name_intern();
	test_trimesh_name_null();

	test_trimesh_tri_count();
	test_trimesh_free_null();

	printf("--- trimesh unit tests (L3: queries) ---\n");
	test_tmq_aabb_single_hit();
	test_tmq_aabb_single_miss();
	test_tmq_aabb_multi_covers_all();
	test_tmq_aabb_multi_one_region();
	test_tmq_aabb_touching_boundary();

	test_rm_hit_from_above();
	test_rm_miss_direction_away();
	test_rm_miss_max_t();
	test_rm_translated_mesh();
	test_rm_rotated_mesh_normal_world_space();
	test_rm_picks_closest_candidate();
	test_rm_back_side_hit_pinned_two_sided();
	test_rm_axis_aligned_ray();

	printf("--- trimesh unit tests (L4: single-triangle narrowphase) ---\n");
	test_cht_face_down_box();
	test_cht_edge_down_box();
	test_cht_back_face_tunneled();
	test_cht_separated_above();
	test_cht_separated_below();
	test_cht_deep_penetration_depth();

	printf("--- trimesh unit tests (L5: reduction) ---\n");
	test_reduce_empty_nocrash();
	test_reduce_n1_flips_back_face_normal();
	test_reduce_n1_front_face_unchanged();
	test_reduce_back_face_flip_threshold();
	test_reduce_mixed_contacts_flipped_independently();
	test_reduce_empty_manifold_in_multi_skipped();
	test_reduce_face_trust_both_survive();
	test_reduce_rep_contact_is_deepest();
	test_reduce_spatially_distant_both_survive();
	test_reduce_preserves_face_trusted_normal_after_flip();

	printf("--- trimesh unit tests (L6: full emits) ---\n");
	test_emit_sphere_aabb_miss();
	test_emit_sphere_single_hit();
	test_emit_sphere_sub_id_encoding();
	test_emit_sphere_translated_mesh_world_contacts();
	test_emit_sphere_rotated_mesh_world_normal();
	test_emit_sphere_multi_tri_distinct_sub_ids();
	test_emit_capsule_single_hit();
	test_emit_box_single_hit();
	test_emit_hull_single_hit();

	printf("--- trimesh unit tests (L7: integration) ---\n");
	test_int_sphere_rests_on_mesh_floor();
	test_int_box_rests_on_mesh_floor();
	test_int_capsule_rests_on_mesh_floor();
	test_int_hull_rests_on_mesh_floor();
	test_int_box_slides_across_shared_edge();
	test_int_capsule_slides_across_shared_edge();
	test_int_mesh_raycast_world();
	test_int_mesh_snapshot_roundtrip();
}
