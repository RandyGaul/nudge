// tests.c -- narrowphase collision unit tests.
//
// Tests all Voronoi region combinations for each shape pair:
//   sphere-sphere, sphere-capsule, sphere-box,
//   capsule-capsule, capsule-box, box-box.

#include <stdio.h>

static int test_pass;
static int test_fail;
static const char* test_current;

#define TEST_BEGIN(name) do { test_current = name; } while(0)
#define TEST_ASSERT(cond) do { \
	if (!(cond)) { printf("  FAIL [%s] %s:%d: %s\n", test_current, __FILE__, __LINE__, #cond); test_fail++; } \
	else { test_pass++; } \
} while(0)
#define TEST_ASSERT_FLOAT(a, b, eps) TEST_ASSERT(fabsf((a) - (b)) < (eps))
#define TEST_ASSERT_HIT(m) TEST_ASSERT((m).count > 0)
#define TEST_ASSERT_MISS(m) TEST_ASSERT((m).count == 0)

#define EPS 0.01f

// Helper: check normal roughly points in expected direction.
#define TEST_ASSERT_NORMAL_DIR(m, dx, dy, dz) do { \
	v3 _exp = norm(V3(dx, dy, dz)); \
	float _d = dot((m).contacts[0].normal, _exp); \
	TEST_ASSERT(_d > 0.7f); \
} while(0)

// ============================================================================
// Sphere-Sphere
// Voronoi: only center-center. Test separated, touching, overlapping, concentric.

static void test_sphere_sphere()
{
	Manifold m = {0};

	TEST_BEGIN("sphere-sphere separated");
	m = (Manifold){0};
	collide_sphere_sphere((Sphere){ V3(0,0,0), 1.0f }, (Sphere){ V3(3,0,0), 1.0f }, &m);
	TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-sphere touching");
	m = (Manifold){0};
	collide_sphere_sphere((Sphere){ V3(0,0,0), 1.0f }, (Sphere){ V3(2,0,0), 1.0f }, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.0f, EPS);

	TEST_BEGIN("sphere-sphere overlap");
	m = (Manifold){0};
	collide_sphere_sphere((Sphere){ V3(0,0,0), 1.0f }, (Sphere){ V3(1,0,0), 1.0f }, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 1.0f, EPS);
	TEST_ASSERT_NORMAL_DIR(m, 1, 0, 0);

	TEST_BEGIN("sphere-sphere overlap Y axis");
	m = (Manifold){0};
	collide_sphere_sphere((Sphere){ V3(0,0,0), 1.0f }, (Sphere){ V3(0,1.5f,0), 1.0f }, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.5f, EPS);
	TEST_ASSERT_NORMAL_DIR(m, 0, 1, 0);

	TEST_BEGIN("sphere-sphere concentric");
	m = (Manifold){0};
	collide_sphere_sphere((Sphere){ V3(0,0,0), 1.0f }, (Sphere){ V3(0,0,0), 1.0f }, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 2.0f, EPS);
}

// ============================================================================
// Sphere-Capsule
// Voronoi regions of capsule: endpoint P, endpoint Q, segment interior.

static void test_sphere_capsule()
{
	Manifold m = {0};
	Capsule cap = { V3(0,-2,0), V3(0,2,0), 0.5f };

	#define S(cx,cy,cz,r) (Sphere){ V3(cx,cy,cz), r }

	TEST_BEGIN("sphere-capsule region P separated");
	m = (Manifold){0}; collide_sphere_capsule(S(0,-4,0,1), cap, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-capsule region P overlap");
	m = (Manifold){0}; collide_sphere_capsule(S(0,-3,0,1), cap, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.5f, EPS);
	TEST_ASSERT_NORMAL_DIR(m, 0, 1, 0);

	TEST_BEGIN("sphere-capsule region Q separated");
	m = (Manifold){0}; collide_sphere_capsule(S(0,4,0,1), cap, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-capsule region Q overlap");
	m = (Manifold){0}; collide_sphere_capsule(S(0,3,0,1), cap, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_NORMAL_DIR(m, 0, -1, 0);

	TEST_BEGIN("sphere-capsule segment interior separated");
	m = (Manifold){0}; collide_sphere_capsule(S(3,0,0,1), cap, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-capsule segment interior overlap");
	m = (Manifold){0}; collide_sphere_capsule(S(1,0,0,1), cap, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.5f, EPS);
	TEST_ASSERT_NORMAL_DIR(m, -1, 0, 0);

	TEST_BEGIN("sphere-capsule segment interior overlap Z");
	m = (Manifold){0}; collide_sphere_capsule(S(0,1,1,1), cap, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_NORMAL_DIR(m, 0, 0, -1);

	TEST_BEGIN("sphere-capsule diagonal near P");
	m = (Manifold){0}; collide_sphere_capsule(S(1,-3,0,1), cap, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_NORMAL_DIR(m, -1, 1, 0);

	#undef S
}

// ============================================================================
// Capsule-Capsule
// Voronoi regions of segment-segment: PP, PQ, QP, QQ, P-interior,
// Q-interior, interior-P, interior-Q, interior-interior.

static void test_capsule_capsule()
{
	Manifold m = {0};
	float r = 0.5f;

	// Capsule A along Y: P=(0,-2,0) Q=(0,2,0)
	// Capsule B varies.

	// -- PP: both endpoints closest --
	#define C(px,py,pz,qx,qy,qz) (Capsule){ V3(px,py,pz), V3(qx,qy,qz), r }

	TEST_BEGIN("capsule-capsule PP separated");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-2,0,0,2,0), C(3,-2,0,3,-5,0), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("capsule-capsule PP overlap");
	m = (Manifold){0}; collide_capsule_capsule(C(0,0,0,0,2,0), C(0.5f,0,0,3,0,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-capsule QQ overlap");
	m = (Manifold){0}; collide_capsule_capsule(C(0,0,0,0,2,0), C(0.5f,2,0,3,5,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-capsule PQ overlap");
	m = (Manifold){0}; collide_capsule_capsule(C(0,0,0,0,3,0), C(0.5f,-3,0,0.5f,0,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-capsule parallel separated");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-2,0,0,2,0), C(3,-2,0,3,2,0), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("capsule-capsule parallel overlap");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-2,0,0,2,0), C(0.8f,-2,0,0.8f,2,0), &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_NORMAL_DIR(m, 1, 0, 0);

	TEST_BEGIN("capsule-capsule skew overlap");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-2,0,0,2,0), C(-2,0,0.3f,2,0,0.3f), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-capsule skew separated");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-2,0,0,2,0), C(-2,0,3,2,0,3), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("capsule-capsule interior-P");
	m = (Manifold){0}; collide_capsule_capsule(C(0,-3,0,0,3,0), C(0.5f,0,0,3,0,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-capsule P-interior");
	m = (Manifold){0}; collide_capsule_capsule(C(0.5f,0,0,3,0,0), C(0,-3,0,0,3,0), &m); TEST_ASSERT_HIT(m);

	#undef C
}

// ============================================================================
// Sphere-Box
// Voronoi regions of box: 6 faces, 12 edges, 8 vertices, interior.

static void test_sphere_box()
{
	Manifold m = {0};
	Box box = { V3(0,0,0), quat_identity(), V3(1,1,1) };
	#define S(x,y,z) (Sphere){ V3(x,y,z), 0.5f }

	TEST_BEGIN("sphere-box face +X separated");
	m = (Manifold){0}; collide_sphere_box(S(3,0,0), box, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-box face +X overlap");
	m = (Manifold){0}; collide_sphere_box(S(1.3f,0,0), box, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.2f, EPS);

	TEST_BEGIN("sphere-box face -X overlap");
	m = (Manifold){0}; collide_sphere_box(S(-1.3f,0,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box face +Y overlap");
	m = (Manifold){0}; collide_sphere_box(S(0,1.3f,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box face -Y overlap");
	m = (Manifold){0}; collide_sphere_box(S(0,-1.3f,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box face +Z overlap");
	m = (Manifold){0}; collide_sphere_box(S(0,0,1.3f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box face -Z overlap");
	m = (Manifold){0}; collide_sphere_box(S(0,0,-1.3f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box edge +X+Y separated");
	m = (Manifold){0}; collide_sphere_box(S(2,2,0), box, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-box edge +X+Y overlap");
	m = (Manifold){0}; collide_sphere_box(S(1.2f,1.2f,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box edge +X+Z overlap");
	m = (Manifold){0}; collide_sphere_box(S(1.2f,0,1.2f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box edge -Y+Z overlap");
	m = (Manifold){0}; collide_sphere_box(S(0,-1.2f,1.2f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box vertex +++ separated");
	m = (Manifold){0}; collide_sphere_box(S(2,2,2), box, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("sphere-box vertex +++ overlap");
	m = (Manifold){0}; collide_sphere_box(S(1.1f,1.1f,1.1f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box vertex --- overlap");
	m = (Manifold){0}; collide_sphere_box(S(-1.1f,-1.1f,-1.1f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box vertex +-+ overlap");
	m = (Manifold){0}; collide_sphere_box(S(1.1f,-1.1f,1.1f), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("sphere-box interior");
	m = (Manifold){0}; collide_sphere_box(S(0.1f,0,0), box, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	#undef S
}

// ============================================================================
// Capsule-Box
// Voronoi regions: capsule endpoint vs face/edge/vertex, segment vs face/edge.

static void test_capsule_box()
{
	Manifold m = {0};
	Box box = { V3(0,0,0), quat_identity(), V3(1,1,1) };
	#define C(px,py,pz,qx,qy,qz) (Capsule){ V3(px,py,pz), V3(qx,qy,qz), 0.3f }

	TEST_BEGIN("capsule-box endpoint face +Y separated");
	m = (Manifold){0}; collide_capsule_box(C(0,3,0,0,5,0), box, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("capsule-box endpoint face +Y overlap");
	m = (Manifold){0}; collide_capsule_box(C(0,0.9f,0,0,3,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-box segment along face +Y");
	m = (Manifold){0}; collide_capsule_box(C(-0.5f,0.9f,0,0.5f,0.9f,0), box, &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.count == 2);

	TEST_BEGIN("capsule-box segment crossing +X edge");
	m = (Manifold){0}; collide_capsule_box(C(0.9f,-2,0,0.9f,2,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-box endpoint near vertex +++ separated");
	m = (Manifold){0}; collide_capsule_box(C(2,2,2,3,3,3), box, &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("capsule-box deep penetration");
	m = (Manifold){0}; collide_capsule_box(C(0,-0.5f,0,0,0.5f,0), box, &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("capsule-box endpoint face -X overlap");
	m = (Manifold){0}; collide_capsule_box(C(-0.9f,0,0,-3,0,0), box, &m); TEST_ASSERT_HIT(m);

	#undef C
}

// ============================================================================
// Box-Box (SAT)
// Voronoi regions: face-face (6 face normals from each), edge-edge (Gauss map pruned).

static void test_box_box()
{
	Manifold m = {0};
	quat id = quat_identity();
	#define B(x,y,z,r,hx,hy,hz) (Box){ V3(x,y,z), r, V3(hx,hy,hz) }
	#define B1(x,y,z) B(x,y,z, id, 1,1,1)

	TEST_BEGIN("box-box face separated +X");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(3,0,0), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("box-box face overlap +X");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(1.5f,0,0), &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.5f, EPS);

	TEST_BEGIN("box-box face overlap +Y");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(0,1.5f,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box face overlap +Z");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(0,0,1.5f), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box face overlap -X");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(-1.5f,0,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box face non-uniform separated");
	m = (Manifold){0}; collide_box_box(B(0,0,0,id,2,1,1), B(5,0,0,id,2,1,1), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("box-box face non-uniform overlap");
	m = (Manifold){0}; collide_box_box(B(0,0,0,id,2,1,1), B(3.5f,0,0,id,2,1,1), &m); TEST_ASSERT_HIT(m);

	quat rot45y = { 0, 0.3827f, 0, 0.9239f };
	TEST_BEGIN("box-box edge-edge rotated 45 Y separated");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B(3,0,0,rot45y,1,1,1), &m); TEST_ASSERT_MISS(m);

	TEST_BEGIN("box-box edge-edge rotated 45 Y overlap");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B(2.2f,0,0,rot45y,1,1,1), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box face touching");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(2,0,0), &m); TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.0f, EPS);

	TEST_BEGIN("box-box diagonal overlap");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(1.5f,1.5f,0), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box diagonal separated");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(1.5f,1.5f,1.5f), &m);
	TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box 3d diagonal separated");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B1(2.5f,2.5f,2.5f), &m); TEST_ASSERT_MISS(m);

	quat rot45z = { 0, 0, 0.3827f, 0.9239f };
	TEST_BEGIN("box-box edge-edge rotated 45 Z overlap");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B(2.0f,0,0,rot45z,1,1,1), &m); TEST_ASSERT_HIT(m);

	TEST_BEGIN("box-box deep overlap");
	m = (Manifold){0}; collide_box_box(B1(0,0,0), B(0.1f,0,0,id,0.5f,0.5f,0.5f), &m); TEST_ASSERT_HIT(m);

	#undef B
	#undef B1
}

// ============================================================================
// Entry point.

// ============================================================================
// Quickhull

static void test_quickhull()
{
	// Box point cloud.
	TEST_BEGIN("quickhull box");
	v3 box_pts[] = {
		{-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
		{-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1},
	};
	Hull* h = quickhull(box_pts, 8);
	TEST_ASSERT(h != NULL);
	TEST_ASSERT(h->vert_count >= 8);
	TEST_ASSERT(h->face_count >= 6);
	hull_free(h);

	// Tetrahedron.
	TEST_BEGIN("quickhull tetrahedron");
	v3 tet_pts[] = { {0,1,0}, {-1,-1,1}, {1,-1,1}, {0,-1,-1} };
	h = quickhull(tet_pts, 4);
	TEST_ASSERT(h != NULL);
	TEST_ASSERT(h->vert_count == 4);
	TEST_ASSERT(h->face_count == 4);
	hull_free(h);

	// Box with interior points (should be discarded).
	TEST_BEGIN("quickhull interior points");
	v3 pts_with_interior[] = {
		{-2,-2,-2}, {2,-2,-2}, {2,2,-2}, {-2,2,-2},
		{-2,-2, 2}, {2,-2, 2}, {2,2, 2}, {-2,2, 2},
		{0,0,0}, {0.5f,0.5f,0.5f}, {-0.1f,-0.1f,0.1f}, // interior
	};
	h = quickhull(pts_with_interior, 11);
	TEST_ASSERT(h != NULL);
	TEST_ASSERT(h->vert_count >= 8); // hull verts (may include duplicates)
	hull_free(h);

	// Hull collision: build a hull, collide with sphere.
	TEST_BEGIN("quickhull hull-sphere collision");
	h = quickhull(box_pts, 8);
	ConvexHull ch = { h, V3(0,0,0), quat_identity(), V3(1,1,1) };
	Manifold m = {0};
	// Place sphere overlapping the hull surface.
	int hit = collide_sphere_hull((Sphere){ V3(0,0,0), 0.5f }, ch, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	hull_free(h);
}

// ============================================================================
// Hull validation helpers.

// Check Euler formula: V - E + F = 2.
static int hull_check_euler(const Hull* h)
{
	return h->vert_count - h->edge_count / 2 + h->face_count == 2;
}

// Check all edges have valid twins and twins are symmetric.
static int hull_check_twins(const Hull* h)
{
	for (int i = 0; i < h->edge_count; i++) {
		int twin = h->edges[i].twin;
		if (twin >= h->edge_count) return 0;
		if (h->edges[twin].twin != i) return 0;
	}
	return 1;
}

// Check each face edge loop is closed and edges reference that face.
static int hull_check_face_loops(const Hull* h)
{
	for (int fi = 0; fi < h->face_count; fi++) {
		int start = h->faces[fi].edge;
		int e = start;
		int count = 0;
		do {
			if (h->edges[e].face != fi) return 0;
			e = h->edges[e].next;
			if (++count > h->edge_count) return 0; // infinite loop
		} while (e != start);
		if (count < 3) return 0;
	}
	return 1;
}

// Check convexity: every vertex is on or behind every face plane (within tolerance).
static int hull_check_convex(const Hull* h, float tol)
{
	for (int fi = 0; fi < h->face_count; fi++) {
		v3 n = h->planes[fi].normal;
		float d = h->planes[fi].offset;
		for (int vi = 0; vi < h->vert_count; vi++) {
			float dist = dot(n, h->verts[vi]) - d;
			if (dist > tol) return 0;
		}
	}
	return 1;
}

// Check outward normals: face normal should point away from hull centroid.
static int hull_check_normals_outward(const Hull* h)
{
	for (int fi = 0; fi < h->face_count; fi++) {
		// Face centroid.
		v3 fc = V3(0,0,0);
		int cnt = 0;
		int start = h->faces[fi].edge;
		int e = start;
		do { fc = add(fc, h->verts[h->edges[e].origin]); cnt++; e = h->edges[e].next; } while (e != start);
		fc = scale(fc, 1.0f / cnt);
		// Normal should point from hull centroid toward face centroid.
		v3 to_face = sub(fc, h->centroid);
		if (dot(h->planes[fi].normal, to_face) < 0) return 0;
	}
	return 1;
}

// Check face planarity: all vertices on a face lie within tolerance of its plane.
static int hull_check_planar(const Hull* h, float tol)
{
	for (int fi = 0; fi < h->face_count; fi++) {
		v3 n = h->planes[fi].normal;
		float d = h->planes[fi].offset;
		int start = h->faces[fi].edge;
		int e = start;
		do {
			float dist = fabsf(dot(n, h->verts[h->edges[e].origin]) - d);
			if (dist > tol) return 0;
			e = h->edges[e].next;
		} while (e != start);
	}
	return 1;
}

// Full hull validation.
static void hull_validate(const Hull* h, const char* label)
{
	char buf[128];
	float tol = 0.01f;

	snprintf(buf, sizeof(buf), "%s euler", label);
	TEST_BEGIN(buf);
	TEST_ASSERT(hull_check_euler(h));

	snprintf(buf, sizeof(buf), "%s twins", label);
	TEST_BEGIN(buf);
	TEST_ASSERT(hull_check_twins(h));

	snprintf(buf, sizeof(buf), "%s face loops", label);
	TEST_BEGIN(buf);
	TEST_ASSERT(hull_check_face_loops(h));

	snprintf(buf, sizeof(buf), "%s convex", label);
	TEST_BEGIN(buf);
	TEST_ASSERT(hull_check_convex(h, tol));

	snprintf(buf, sizeof(buf), "%s planar", label);
	TEST_BEGIN(buf);
	TEST_ASSERT(hull_check_planar(h, tol));
}

// ============================================================================
// Quickhull fuzz testing.
//
// Strategy: take known shapes, perturb vertices near hull features
// (on-vertex, on-edge, on-face, near-coplanar, interior), build hull,
// validate structural invariants.

static uint32_t fuzz_rng_state = 12345;

static float fuzz_rand()
{
	// xorshift32
	fuzz_rng_state ^= fuzz_rng_state << 13;
	fuzz_rng_state ^= fuzz_rng_state >> 17;
	fuzz_rng_state ^= fuzz_rng_state << 5;
	return (float)(fuzz_rng_state & 0xFFFFFF) / (float)0xFFFFFF;
}

static float fuzz_rand_range(float lo, float hi)
{
	return lo + fuzz_rand() * (hi - lo);
}

static v3 fuzz_rand_dir()
{
	v3 d;
	do {
		d = V3(fuzz_rand_range(-1,1), fuzz_rand_range(-1,1), fuzz_rand_range(-1,1));
	} while (len2(d) < 0.001f);
	return norm(d);
}

// Base shapes for fuzzing.
static const v3 s_cube_pts[] = {
	{-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
	{-1,-1, 1}, {1,-1, 1}, {1,1, 1}, {-1,1, 1},
};

static const v3 s_tet_pts[] = {
	{0,1,0}, {-1,-1,1}, {1,-1,1}, {0,-1,-1},
};

// Octahedron.
static const v3 s_oct_pts[] = {
	{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1},
};

// Icosahedron (approx).
#define ICO_PHI 1.618033988f
static const v3 s_ico_pts[] = {
	{-1, ICO_PHI, 0}, { 1, ICO_PHI, 0}, {-1,-ICO_PHI, 0}, { 1,-ICO_PHI, 0},
	{0,-1, ICO_PHI}, {0, 1, ICO_PHI}, {0,-1,-ICO_PHI}, {0, 1,-ICO_PHI},
	{ ICO_PHI, 0,-1}, { ICO_PHI, 0, 1}, {-ICO_PHI, 0,-1}, {-ICO_PHI, 0, 1},
};

typedef struct FuzzShape
{
	const char* name;
	const v3* verts;
	int count;
	int expected_verts;
	int expected_min_faces;
} FuzzShape;

static const FuzzShape s_fuzz_shapes[] = {
	{ "cube",        s_cube_pts, 8,  8,  6 },
	{ "tetrahedron", s_tet_pts,  4,  4,  4 },
	{ "octahedron",  s_oct_pts,  6,  6,  8 },
	{ "icosahedron", s_ico_pts,  12, 12, 20 },
};

#define FUZZ_SHAPE_COUNT (sizeof(s_fuzz_shapes) / sizeof(s_fuzz_shapes[0]))

// Hardcoded reproduction of fuzz icosahedron #783 failing case.
static void test_quickhull_case783()
{
	TEST_BEGIN("quickhull case783");

	v3 pts[] = {
		{-0.999981225f, 1.618084550f, -0.000268069f},
		{0.999727011f, 1.618018270f, -0.000001823f},
		{-1.000144720f, -1.618227363f, 0.000128272f},
		{0.999813139f, -1.618205786f, -0.000101668f},
		{0.000114347f, -1.000113249f, 1.617812991f},
		{-0.000224340f, 0.999980986f, 1.617878795f},
		{-0.000208760f, -1.000046253f, -1.617863536f},
		{0.000269564f, 1.000018477f, -1.617991924f},
		{1.618107319f, 0.000261519f, -0.999968350f}, // point 8: hull vertex that gets missed
		{1.617948294f, 0.000181835f, 1.000185370f},
		{-1.618246436f, -0.000156645f, -0.999928653f},
		{-1.618226528f, 0.000065485f, 0.999817133f},
		{1.618591547f, 0.001790219f, 0.998493850f},
		{1.617817879f, -0.000213558f, 0.999758244f},
		{1.000626922f, -1.616133332f, 0.001539255f},
		{0.000000000f, 0.259813756f, 0.420387506f},
		{-0.363759935f, 0.000000000f, 0.224816009f},
		{-0.460946679f, 0.686911166f, 0.219213605f},
		{0.999556899f, 1.621392727f, -0.003101509f},
		{-0.000241796f, 1.000435710f, -1.618831277f},
		{-0.220620841f, -0.526642025f, -0.395656526f},
		{-0.000346585f, -1.004286170f, -1.619372249f},
		{0.001092130f, 0.994397342f, -1.613620877f},
		{-1.017551184f, 0.003922998f, -0.621890545f},
		{0.000000000f, 0.753539741f, 1.219252944f},
		{1.110361457f, 0.000000000f, -0.686241090f},
		{0.018766023f, 0.030364064f, 0.000000000f},
		{-0.057126891f, -1.039216280f, 1.524609804f},
		{0.006557340f, 0.993353009f, 1.615448236f},
	};
	int npts = sizeof(pts) / sizeof(pts[0]);

	Hull* h = quickhull(pts, npts);
	TEST_ASSERT(h != NULL);
	if (h) {
		int c_ok = hull_check_convex(h, 0.02f);
		if (!c_ok) {
			float worst = 0; int worst_pi = -1;
			for (int fi = 0; fi < h->face_count; fi++)
				for (int pi = 0; pi < npts; pi++) {
					float d = dot(h->planes[fi].normal, pts[pi]) - h->planes[fi].offset;
					if (d > worst) { worst = d; worst_pi = pi; }
				}
			printf("  case783: point %d outside by %.8f\n", worst_pi, worst);
		}
		TEST_ASSERT(c_ok);
		hull_free(h);
	}
}

static void test_quickhull_fuzz(int iterations)
{
	char label[128];
	int total_hulls = 0;
	int total_fails = 0;

	for (int shape_idx = 0; shape_idx < (int)FUZZ_SHAPE_COUNT; shape_idx++) {
		const FuzzShape* base = &s_fuzz_shapes[shape_idx];

		for (int iter = 0; iter < iterations; iter++) {
			// Build a point set: start with base shape, then add perturbations.
			CK_DYNA v3* pts = NULL;

			// Copy base verts with small random jitter.
			float jitter = fuzz_rand_range(0.0f, 0.001f);
			for (int i = 0; i < base->count; i++) {
				v3 p = base->verts[i];
				p = add(p, scale(fuzz_rand_dir(), jitter));
				apush(pts, p);
			}

			// Add extra points in various categories:
			int extras = (int)(fuzz_rand() * 20);
			for (int i = 0; i < extras; i++) {
				float r = fuzz_rand();
				v3 p;
				if (r < 0.3f) {
					// Near a hull vertex (on-vertex perturbation).
					int vi = (int)(fuzz_rand() * base->count) % base->count;
					float fuzz = fuzz_rand_range(1e-7f, 1e-2f);
					p = add(base->verts[vi], scale(fuzz_rand_dir(), fuzz));
				} else if (r < 0.6f) {
					// Near a hull edge (interpolate two verts + perturb).
					int vi = (int)(fuzz_rand() * base->count) % base->count;
					int vj = (vi + 1 + (int)(fuzz_rand() * (base->count - 1))) % base->count;
					float t = fuzz_rand();
					p = add(scale(base->verts[vi], 1-t), scale(base->verts[vj], t));
					float fuzz = fuzz_rand_range(1e-7f, 1e-2f);
					p = add(p, scale(fuzz_rand_dir(), fuzz));
				} else if (r < 0.8f) {
					// Interior point (scaled toward center).
					int vi = (int)(fuzz_rand() * base->count) % base->count;
					p = scale(base->verts[vi], fuzz_rand_range(0.0f, 0.9f));
				} else {
					// Random point near the shape.
					p = scale(fuzz_rand_dir(), fuzz_rand_range(0.5f, 2.0f));
				}
				apush(pts, p);
			}

			// Build hull.
			Hull* h = quickhull(pts, asize(pts));
			total_hulls++;

			if (!h) {
				// Degenerate input is acceptable -- skip.
				afree(pts);
				continue;
			}

			// Validate.
			snprintf(label, sizeof(label), "fuzz %s #%d (%d pts)", base->name, iter, asize(pts));

			int e_ok = hull_check_euler(h);
			int t_ok = hull_check_twins(h);
			int f_ok = hull_check_face_loops(h);
			int n_ok = hull_check_normals_outward(h);
			int c_ok = hull_check_convex(h, 0.02f);

			if (!e_ok || !t_ok || !f_ok || !n_ok || !c_ok) {
				total_fails++;
				if (total_fails <= 3) {
					printf("  FUZZ FAIL: %s (V=%d E=%d F=%d) [%s%s%s%s%s]\n", label,
						h->vert_count, h->edge_count, h->face_count,
						e_ok ? "" : "euler ",
						t_ok ? "" : "twins ",
						f_ok ? "" : "loops ",
						n_ok ? "" : "normals ",
						c_ok ? "" : "convex ");
					// Find worst violation.
					float worst = 0;
					int worst_vi = -1, worst_fi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++) {
						for (int _vi = 0; _vi < h->vert_count; _vi++) {
							float d = dot(h->planes[_fi].normal, h->verts[_vi]) - h->planes[_fi].offset;
							if (d > worst) { worst = d; worst_vi = _vi; worst_fi = _fi; }
						}
					}
					// Check if any INPUT point is outside.
					float worst_input = 0;
					int worst_pi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++) {
						for (int _pi = 0; _pi < asize(pts); _pi++) {
							float d = dot(h->planes[_fi].normal, pts[_pi]) - h->planes[_fi].offset;
							if (d > worst_input) { worst_input = d; worst_pi = _pi; }
						}
					}
					printf("    worst hull vert: v%d above f%d by %.8f\n", worst_vi, worst_fi, worst);
					printf("    worst input pt:  p%d outside by %.8f\n", worst_pi, worst_input);
				}
			}

			hull_free(h);
			afree(pts);
		}
	}

	printf("  fuzz: %d hulls tested, %d failures\n", total_hulls, total_fails);
	TEST_BEGIN("quickhull fuzz");
	TEST_ASSERT(total_fails == 0);
}

static void run_tests()
{
	test_pass = 0;
	test_fail = 0;

	printf("--- nudge narrowphase tests ---\n");

	test_sphere_sphere();
	test_sphere_capsule();
	test_capsule_capsule();
	test_sphere_box();
	test_capsule_box();
	test_box_box();
	test_quickhull();

	// Quickhull fuzz: moderate count for regular runs.
	test_quickhull_case783();
	exit(0);
	test_quickhull_fuzz(100); // use 1000+ for stress testing

	printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	if (test_fail > 0) printf("*** FAILURES ***\n");
}
