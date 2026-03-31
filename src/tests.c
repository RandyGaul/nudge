// See LICENSE for licensing info.
// tests.c -- narrowphase collision unit tests.
//
// Tests all Voronoi region combinations for each shape pair:
//   sphere-sphere, sphere-capsule, sphere-box,
//   capsule-capsule, capsule-box, box-box.

#include <stdio.h>
#include <time.h>

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
	TEST_ASSERT(hull_check_convex(h, h->epsilon));
}

// ============================================================================
// Quickhull fuzz testing.
//
// Strategy: take known shapes, perturb vertices near hull features
// (on-vertex, on-edge, on-face, near-coplanar, interior), build hull,
// validate structural invariants.

static uint32_t fuzz_rng_state = 12345;
static uint32_t fuzz_rng_seed  = 12345;

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

// Edge pairs (vertex index pairs) for each base shape.
static const int s_cube_edges[][2] = {
	{0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4}, {0,4},{1,5},{2,6},{3,7},
};
static const int s_tet_edges[][2] = {
	{0,1},{0,2},{0,3},{1,2},{1,3},{2,3},
};
static const int s_oct_edges[][2] = {
	{0,2},{0,3},{0,4},{0,5},{1,2},{1,3},{1,4},{1,5},{2,4},{4,3},{3,5},{5,2},
};
static const int s_ico_edges[][2] = {
	{0,1},{0,5},{0,7},{0,10},{0,11},{1,5},{1,7},{1,8},{1,9},
	{2,3},{2,4},{2,10},{2,11},{3,4},{3,8},{3,9},
	{4,5},{4,11},{5,9},{6,7},{6,8},{6,10},{7,10},
	{8,3},{8,6},{9,3},{9,4},{10,6},{10,11},{11,4},
};

// Face triangles for each base shape (triangulated).
static const int s_cube_tris[][3] = {
	{0,2,1},{0,3,2}, {4,5,6},{4,6,7}, {0,1,5},{0,5,4},
	{2,3,7},{2,7,6}, {0,4,7},{0,7,3}, {1,2,6},{1,6,5},
};
static const int s_tet_tris[][3] = {
	{0,1,2},{0,3,1},{0,2,3},{1,3,2},
};
static const int s_oct_tris[][3] = {
	{0,4,2},{0,2,5},{0,5,3},{0,3,4},{1,2,4},{1,5,2},{1,3,5},{1,4,3},
};
static const int s_ico_tris[][3] = {
	{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},
	{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
	{3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},
	{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1},
};

typedef struct FuzzShape
{
	const char* name;
	const v3* verts;
	int vert_count;
	const int (*edges)[2];
	int edge_count;
	const int (*tris)[3];
	int tri_count;
} FuzzShape;

#define COUNTOF(a) (int)(sizeof(a)/sizeof(a[0]))

static const FuzzShape s_fuzz_shapes[] = {
	{ "cube",        s_cube_pts, 8,  s_cube_edges, COUNTOF(s_cube_edges), s_cube_tris, COUNTOF(s_cube_tris) },
	{ "tetrahedron", s_tet_pts,  4,  s_tet_edges,  COUNTOF(s_tet_edges),  s_tet_tris,  COUNTOF(s_tet_tris) },
	{ "octahedron",  s_oct_pts,  6,  s_oct_edges,  COUNTOF(s_oct_edges),  s_oct_tris,  COUNTOF(s_oct_tris) },
	{ "icosahedron", s_ico_pts,  12, s_ico_edges,  COUNTOF(s_ico_edges),  s_ico_tris,  COUNTOF(s_ico_tris) },
};

#define FUZZ_SHAPE_COUNT COUNTOF(s_fuzz_shapes)

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
		int c_ok = hull_check_convex(h, h->epsilon);
		if (!c_ok) {
			float worst = 0; int worst_pi = -1;
			for (int fi = 0; fi < h->face_count; fi++)
				for (int pi = 0; pi < npts; pi++) {
					float d = dot(h->planes[fi].normal, pts[pi]) - h->planes[fi].offset;
					if (d > worst) { worst = d; worst_pi = pi; }
				}
			printf("  case783: point %d outside by %.8f (epsilon=%.2e)\n", worst_pi, worst, h->epsilon);
		}
		TEST_ASSERT(c_ok);
		hull_free(h);
	}
}

// Fuzz icosahedron #7037: euler/twin topology corruption.
static void test_quickhull_ico7037()
{
	TEST_BEGIN("quickhull ico7037");
	v3 pts[] = {
		{-9.998392463e-01f, 1.618750095e+00f, 5.350695574e-04f},
		{1.000362515e+00f, 1.618721724e+00f, 4.695687385e-04f},
		{-9.991761446e-01f, -1.617793798e+00f, 2.973712108e-04f},
		{1.000283957e+00f, -1.617188692e+00f, 1.721137232e-04f},
		{5.902615376e-04f, -1.000672579e+00f, 1.618189454e+00f},
		{-2.665717329e-05f, 9.995228648e-01f, 1.617261648e+00f},
		{-8.071046905e-04f, -9.999146461e-01f, -1.617626429e+00f},
		{-7.048695697e-04f, 9.997880459e-01f, -1.618566036e+00f},
		{1.617915750e+00f, -7.171047037e-04f, -1.000544667e+00f},
		{1.618297100e+00f, 7.590552559e-04f, 1.000423670e+00f},
		{-1.618417263e+00f, -6.951298565e-04f, -1.000441313e+00f},
		{-1.617776036e+00f, 8.026060095e-05f, 9.991329312e-01f},
		{-6.373115350e-03f, -9.993319511e-01f, -1.619585514e+00f},
		{1.213525295e+00f, -2.499996871e-01f, -1.154508233e+00f},
		{8.090168238e-01f, -4.999997318e-01f, -1.309017301e+00f},
		{4.045084119e-01f, -7.499998212e-01f, -1.463525772e+00f},
		{-1.399685621e+00f, 5.716435909e-01f, 4.926615953e-02f},
		{-9.301985502e-01f, 3.140034676e-01f, -8.295849562e-01f},
		{-1.348361850e+00f, -1.666672230e-01f, 1.103005290e+00f},
		{-1.078689456e+00f, -3.333331943e-01f, 1.206010938e+00f},
		{-8.090164661e-01f, -4.999998212e-01f, 1.309017062e+00f},
		{-5.393444896e-01f, -6.666666865e-01f, 1.412021995e+00f},
		{-2.696720958e-01f, -8.333333731e-01f, 1.515028000e+00f},
		{-7.183008790e-01f, -5.533111095e-01f, -1.048815489e+00f},
		{6.177458545e-07f, -3.333335817e-01f, 1.618033528e+00f},
		{6.292121952e-07f, 3.333334029e-01f, 1.618034244e+00f},
		{8.570642471e-01f, -1.676853746e-01f, -1.162361279e-01f},
		{-7.177171111e-01f, -5.562684536e-01f, 1.343872428e+00f},
		{-1.278867960e+00f, 0.000000000e+00f, 7.903838754e-01f},
		{-2.423702329e-01f, 8.502069712e-01f, -1.525456786e+00f},
		{-1.602371335e-01f, 5.667126775e-01f, 5.072154403e-01f},
		{-1.454859734e+00f, -4.271956086e-01f, -5.881867409e-01f},
		{9.983093739e-01f, 1.616120100e+00f, 9.734939667e-04f},
	};
	int npts = sizeof(pts) / sizeof(pts[0]);
	Hull* h = quickhull(pts, npts);
	TEST_ASSERT(h != NULL);
	if (h) {
		hull_validate(h, "ico7037");
		hull_free(h);
	}
}

// Fuzz tetrahedron #8098: topology corruption from collinear edge clusters.
static void test_quickhull_tet8098()
{
	TEST_BEGIN("quickhull tet8098");
	v3 pts[] = {
		{-5.337806651e-04f, 1.000092030e+00f, -6.938962615e-04f},
		{-9.994471669e-01f, -9.997874498e-01f, 1.000651240e+00f},
		{9.997070432e-01f, -9.999066591e-01f, 9.991751313e-01f},
		{-1.580310782e-04f, -9.997072220e-01f, -9.991850257e-01f},
		{0.000000000e+00f, 8.300971389e-01f, 0.000000000e+00f},
		{1.044850349e-01f, -8.912307620e-01f, 9.455661178e-01f},
		{-2.500004470e-01f, 5.000004768e-01f, 2.500002980e-01f},
		{-4.999998212e-01f, -2.450137799e-07f, 5.000001788e-01f},
		{-7.499995828e-01f, -5.000000000e-01f, 7.499994636e-01f},
		{1.666666269e-01f, 6.666663289e-01f, 1.666665822e-01f},
		{3.333336711e-01f, 3.333330750e-01f, 3.333334625e-01f},
		{4.999996424e-01f, 1.129004161e-07f, 5.000006557e-01f},
		{6.666668057e-01f, -3.333334029e-01f, 6.666667461e-01f},
		{8.333331347e-01f, -6.666662693e-01f, 8.333335519e-01f},
		{6.678278446e-01f, -2.809335291e-01f, 6.951886415e-01f},
		{1.490461230e-01f, -1.664296389e-01f, -2.851225734e-01f},
		{6.641381979e-01f, -6.641381979e-01f, 6.641381979e-01f},
		{-8.298630118e-01f, -6.601527333e-01f, 8.301952481e-01f},
		{-7.142857909e-01f, -1.000000000e+00f, 1.000000119e+00f},
		{-4.285713434e-01f, -9.999998212e-01f, 9.999994040e-01f},
		{-1.428574622e-01f, -9.999997020e-01f, 9.999996424e-01f},
		{1.428571790e-01f, -1.000000358e+00f, 9.999997616e-01f},
		{4.285711646e-01f, -1.000000238e+00f, 9.999998808e-01f},
		{7.142857313e-01f, -1.000000000e+00f, 1.000000000e+00f},
		{0.000000000e+00f, -6.698234081e-01f, -6.698234081e-01f},
		{-2.959918082e-01f, -1.000000000e+00f, 3.595360518e-01f},
		{-6.061693430e-01f, -2.126072943e-01f, 6.061887145e-01f},
		{-1.145148743e-02f, -1.000101686e+00f, 1.000135541e+00f},
		{0.000000000e+00f, 8.311111927e-01f, -8.444441110e-02f},
		{-6.775360703e-01f, -2.215145826e-01f, 7.540838122e-01f},
		{1.416199088e+00f, -6.233173609e-01f, -1.041991234e+00f},
		{0.000000000e+00f, 1.000000000e+00f, 0.000000000e+00f},
		{8.559031487e-01f, -9.999999404e-01f, 8.832371235e-01f},
		{1.002887964e+00f, -9.991275668e-01f, 9.971886873e-01f},
		{0.000000000e+00f, 1.000000000e+00f, 0.000000000e+00f},
		{5.826290324e-02f, 7.816259265e-01f, 7.234503981e-03f},
		{1.392189413e-01f, -1.392189413e-01f, 1.392189413e-01f},
		{1.777941585e-01f, 3.107546568e-01f, 3.445743918e-01f},
		{3.966283798e-02f, 4.953246713e-01f, -3.316636682e-01f},
		{-1.010569483e-01f, -1.010569483e-01f, 1.010569483e-01f},
		{-4.551982880e-01f, -3.737179339e-01f, 2.235377431e-01f},
		{-8.439664841e-01f, -8.512523174e-01f, 9.256261587e-01f},
		{-3.355232179e-01f, -3.355232179e-01f, 3.355232179e-01f},
		{1.730938554e-01f, -9.999105930e-01f, -1.821184158e-02f},
		{1.000000000e+00f, -1.000000000e+00f, 1.000000000e+00f},
		{2.957146168e-01f, 2.861931026e-01f, 3.986243904e-01f},
		{-4.653079212e-01f, 6.629047990e-01f, 7.649911046e-01f},
		{9.992291033e-02f, -9.992291033e-02f, 9.992291033e-02f},
		{0.000000000e+00f, 4.690935612e-01f, 0.000000000e+00f},
		{-1.554469913e-01f, -9.999763966e-01f, -1.179022640e-01f},
		{2.778943777e-01f, 3.717072010e-01f, -2.266341448e-01f},
		{6.300083478e-04f, 5.262714624e-01f, -2.376070619e-01f},
		{-8.241868019e-02f, -1.000000000e+00f, 9.856159687e-01f},
		{1.666664630e-01f, 6.666671634e-01f, 1.666671485e-01f},
		{3.333337903e-01f, 3.333327174e-01f, 3.333327472e-01f},
		{5.000000596e-01f, -1.742947973e-07f, 5.000001192e-01f},
		{6.666667461e-01f, -3.333328962e-01f, 6.666672826e-01f},
		{8.333339095e-01f, -6.666662097e-01f, 8.333337903e-01f},
		{3.736703396e-01f, -1.000000000e+00f, -2.526593208e-01f},
		{9.833589196e-01f, -1.000387073e+00f, 9.995107651e-01f},
		{-3.333334923e-01f, 3.333329260e-01f, 3.333332837e-01f},
		{-6.666665673e-01f, -3.333339095e-01f, 6.666662097e-01f},
		{-7.641195059e-01f, 6.257113218e-01f, -2.739081383e-01f},
		{0.000000000e+00f, -3.903293312e-01f, -3.903293312e-01f},
		{2.721109092e-01f, 4.672608078e-01f, 1.395917773e+00f},
		{5.039944053e-01f, -1.000681520e+00f, 8.123803884e-03f},
		{0.000000000e+00f, 6.808243990e-01f, 0.000000000e+00f},
		{8.260955215e-01f, -8.592320681e-01f, 7.225213051e-01f},
		{5.945055648e-08f, 5.000000596e-01f, -2.499999851e-01f},
		{4.738698749e-07f, -1.294760352e-07f, -4.999996722e-01f},
		{-9.222519566e-07f, -4.999997020e-01f, -7.500002384e-01f},
		{7.923394442e-02f, -1.000063539e+00f, -5.345927477e-01f},
		{2.671395838e-01f, 4.657208622e-01f, 2.671395838e-01f},
		{7.841181159e-01f, -1.000000000e+00f, 8.455902338e-01f},
		{-6.013367772e-01f, -1.000067711e+00f, 9.996876717e-01f},
		{9.940116405e-01f, -9.970571399e-01f, 1.003938913e+00f},
		{-2.634056611e-03f, -1.003241777e+00f, -9.956597090e-01f},
		{1.965457201e-01f, -4.393271804e-01f, 7.196121216e-01f},
		{-7.210817337e-01f, 3.009600639e-01f, -3.347507119e-01f},
		{6.105791032e-02f, -1.000000000e+00f, -1.498298407e+00f},
		{4.062898755e-01f, -9.999670386e-01f, 2.205693722e-02f},
		{4.144960046e-01f, -2.809199095e-01f, 6.404599547e-01f},
		{3.797146678e-01f, -3.797146678e-01f, 3.797146678e-01f},
		{5.042100674e-04f, -7.175825238e-01f, -8.593770266e-01f},
		{8.309749688e-08f, 6.666666865e-01f, -1.666668504e-01f},
		{-6.299016491e-07f, 3.333329558e-01f, -3.333338797e-01f},
		{5.598515713e-07f, -2.059473019e-07f, -4.999996126e-01f},
		{-4.236425752e-07f, -3.333335221e-01f, -6.666664481e-01f},
		{2.755549247e-07f, -6.666669250e-01f, -8.333336115e-01f},
		{2.957612742e-03f, 9.946585298e-01f, -1.287558000e-03f},
		{0.000000000e+00f, 1.000000000e+00f, 0.000000000e+00f},
		{0.000000000e+00f, 2.191416919e-01f, -3.904291689e-01f},
		{-7.459217310e-02f, 1.496296763e+00f, 9.896406531e-02f},
		{7.533715963e-01f, -7.533715963e-01f, 7.533715963e-01f},
		{0.000000000e+00f, -1.000000000e+00f, -1.000000000e+00f},
		{-5.177285075e-01f, -9.584468007e-01f, 9.792233706e-01f},
		{-2.989606261e-01f, 1.172264934e+00f, 6.840536594e-01f},
		{4.209921360e-01f, -1.056362629e+00f, 1.028181314e+00f},
		{1.000000000e+00f, -1.000000000e+00f, 1.000000000e+00f},
		{-9.972440600e-01f, -1.007318377e+00f, 9.942823648e-01f},
		{7.500000596e-01f, -9.999997020e-01f, 4.999999404e-01f},
		{5.000000596e-01f, -9.999999404e-01f, -3.031544438e-08f},
		{2.499999255e-01f, -1.000000119e+00f, -4.999999106e-01f},
		{-7.235236117e-04f, 8.018327355e-01f, -9.863330424e-02f},
		{-4.490792155e-01f, -6.626325250e-01f, 6.670039892e-02f},
		{1.970293466e-03f, 1.003574491e+00f, -1.577983960e-03f},
		{7.273958325e-01f, -7.273958325e-01f, 7.273958325e-01f},
		{5.083889961e-01f, -1.050412416e+00f, 1.025206327e+00f},
		{-7.142858505e-01f, -1.000000000e+00f, 1.000000238e+00f},
		{-4.285714328e-01f, -9.999997616e-01f, 1.000000238e+00f},
		{-1.428570747e-01f, -1.000000000e+00f, 1.000000000e+00f},
		{1.428574026e-01f, -9.999998212e-01f, 9.999998808e-01f},
		{4.285716116e-01f, -9.999999404e-01f, 9.999998808e-01f},
		{7.142857909e-01f, -9.999998212e-01f, 1.000000715e+00f},
		{0.000000000e+00f, 1.000000000e+00f, 0.000000000e+00f},
		{-3.928686976e-01f, 2.162244916e-01f, 3.917980492e-01f},
		{3.523567226e-03f, 9.966949821e-01f, -1.270852983e-03f},
		{1.435524784e-03f, 9.999296665e-01f, -9.782379493e-04f},
		{3.327154517e-01f, -1.000495911e+00f, -3.347327411e-01f},
		{-4.169500619e-02f, 1.300986111e-02f, 1.240348935e+00f},
		{-1.905424595e-01f, -4.132390022e-01f, -3.255345821e-01f},
		{9.973774552e-01f, -1.000233173e+00f, 9.978656173e-01f},
		{-9.240421057e-01f, -9.262731075e-01f, 8.849738240e-01f},
	};
	int npts = sizeof(pts) / sizeof(pts[0]);
	Hull* h = quickhull(pts, npts);
	TEST_ASSERT(h != NULL);
	if (h) {
		TEST_ASSERT(hull_check_euler(h));
		TEST_ASSERT(hull_check_twins(h));
		hull_free(h);
	}
}

// Check that ALL input points lie inside (or on) every hull face plane.
static int hull_check_contains_inputs(const Hull* h, const v3* pts, int npts, float tol)
{
	for (int fi = 0; fi < h->face_count; fi++) {
		v3 n = h->planes[fi].normal;
		float d = h->planes[fi].offset;
		for (int pi = 0; pi < npts; pi++) {
			if (dot(n, pts[pi]) - d > tol) return 0;
		}
	}
	return 1;
}

static int fuzz_rand_int(int n) { return (int)(fuzz_rand() * n) % n; }

static void test_quickhull_fuzz(int iterations)
{
	char label[128];
	int total_hulls = 0;
	int total_fails = 0;

	for (int shape_idx = 0; shape_idx < (int)FUZZ_SHAPE_COUNT; shape_idx++) {
		const FuzzShape* base = &s_fuzz_shapes[shape_idx];

		for (int iter = 0; iter < iterations; iter++) {
			CK_DYNA v3* pts = NULL;

			// Random uniform scale to stress epsilon at different magnitudes.
			float hull_scale = 1.0f;
			float sr = fuzz_rand();
			if      (sr < 0.1f) hull_scale = fuzz_rand_range(0.001f, 0.01f);
			else if (sr < 0.2f) hull_scale = fuzz_rand_range(10.0f, 1000.0f);

			// Copy base verts with small random jitter.
			float jitter = fuzz_rand_range(0.0f, 0.001f) * hull_scale;
			for (int i = 0; i < base->vert_count; i++) {
				v3 p = scale(base->verts[i], hull_scale);
				p = add(p, scale(fuzz_rand_dir(), jitter));
				apush(pts, p);
			}

			// Add extra points in various categories.
			// Vary count: usually 0-20, sometimes 50-100.
			int extras = (int)(fuzz_rand() * 20);
			if (fuzz_rand() < 0.1f) extras += 50 + (int)(fuzz_rand() * 50);

			for (int i = 0; i < extras; i++) {
				float r = fuzz_rand();
				v3 p;
				if (r < 0.15f) {
					// Near a hull vertex.
					int vi = fuzz_rand_int(base->vert_count);
					float fuzz = fuzz_rand_range(1e-7f, 1e-2f) * hull_scale;
					p = add(scale(base->verts[vi], hull_scale), scale(fuzz_rand_dir(), fuzz));
				} else if (r < 0.20f) {
					// Exact duplicate of a base vertex.
					p = scale(base->verts[fuzz_rand_int(base->vert_count)], hull_scale);
				} else if (r < 0.35f) {
					// On an actual edge (interpolate two adjacent verts).
					int ei = fuzz_rand_int(base->edge_count);
					v3 a = scale(base->verts[base->edges[ei][0]], hull_scale);
					v3 b = scale(base->verts[base->edges[ei][1]], hull_scale);
					float t = fuzz_rand();
					p = add(scale(a, 1-t), scale(b, t));
					// Sometimes add small perturbation.
					if (fuzz_rand() < 0.5f) {
						float fuzz = fuzz_rand_range(1e-7f, 1e-3f) * hull_scale;
						p = add(p, scale(fuzz_rand_dir(), fuzz));
					}
				} else if (r < 0.45f) {
					// Collinear cluster: multiple points along one edge.
					int ei = fuzz_rand_int(base->edge_count);
					v3 a = scale(base->verts[base->edges[ei][0]], hull_scale);
					v3 b = scale(base->verts[base->edges[ei][1]], hull_scale);
					int cluster = 2 + (int)(fuzz_rand() * 5);
					for (int c = 0; c < cluster; c++) {
						float t = (float)(c + 1) / (float)(cluster + 1);
						float fuzz = fuzz_rand_range(0, 1e-6f) * hull_scale;
						apush(pts, add(add(scale(a, 1-t), scale(b, t)), scale(fuzz_rand_dir(), fuzz)));
					}
					continue; // already pushed
				} else if (r < 0.60f) {
					// On a face (barycentric interpolation of a triangle).
					int ti = fuzz_rand_int(base->tri_count);
					v3 a = scale(base->verts[base->tris[ti][0]], hull_scale);
					v3 b = scale(base->verts[base->tris[ti][1]], hull_scale);
					v3 c = scale(base->verts[base->tris[ti][2]], hull_scale);
					float u = fuzz_rand(), v = fuzz_rand();
					if (u + v > 1.0f) { u = 1.0f - u; v = 1.0f - v; }
					p = add(add(scale(a, 1-u-v), scale(b, u)), scale(c, v));
					// Sometimes perturb along face normal.
					if (fuzz_rand() < 0.3f) {
						v3 fn = norm(cross(sub(b, a), sub(c, a)));
						float fuzz = fuzz_rand_range(-1e-4f, 1e-4f) * hull_scale;
						p = add(p, scale(fn, fuzz));
					}
				} else if (r < 0.70f) {
					// Coplanar: project a random point onto a face plane.
					int ti = fuzz_rand_int(base->tri_count);
					v3 a = scale(base->verts[base->tris[ti][0]], hull_scale);
					v3 b = scale(base->verts[base->tris[ti][1]], hull_scale);
					v3 c = scale(base->verts[base->tris[ti][2]], hull_scale);
					v3 fn = norm(cross(sub(b, a), sub(c, a)));
					float plane_d = dot(fn, a);
					p = scale(fuzz_rand_dir(), hull_scale * 1.5f);
					p = sub(p, scale(fn, dot(fn, p) - plane_d));
				} else if (r < 0.85f) {
					// Interior point.
					int vi = fuzz_rand_int(base->vert_count);
					p = scale(scale(base->verts[vi], hull_scale), fuzz_rand_range(0.0f, 0.9f));
				} else {
					// Random point near the shape.
					p = scale(fuzz_rand_dir(), fuzz_rand_range(0.5f, 2.0f) * hull_scale);
				}
				apush(pts, p);
			}

			// Build hull.
			Hull* h = quickhull(pts, asize(pts));
			total_hulls++;

			if (!h) {
				afree(pts);
				continue;
			}

			// Validate.
			snprintf(label, sizeof(label), "fuzz %s #%d (%d pts, scale=%.3g)",
				base->name, iter, asize(pts), hull_scale);

			int e_ok = hull_check_euler(h);
			int t_ok = hull_check_twins(h);
			int f_ok = hull_check_face_loops(h);
			int n_ok = hull_check_normals_outward(h);
			int c_ok = hull_check_convex(h, h->epsilon);
			int i_ok = hull_check_contains_inputs(h, pts, asize(pts), h->epsilon + h->maxoutside);

			if (!e_ok || !t_ok || !f_ok || !n_ok || !c_ok || !i_ok) {
				total_fails++;
				if (total_fails <= 3) {
					printf("  FUZZ FAIL: %s (V=%d E=%d F=%d) [%s%s%s%s%s%s]\n", label,
						h->vert_count, h->edge_count, h->face_count,
						e_ok ? "" : "euler ",
						t_ok ? "" : "twins ",
						f_ok ? "" : "loops ",
						n_ok ? "" : "normals ",
						c_ok ? "" : "convex ",
						i_ok ? "" : "inputs ");
					// Find worst violations.
					float worst = 0;
					int worst_vi = -1, worst_fi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++)
						for (int _vi = 0; _vi < h->vert_count; _vi++) {
							float d = dot(h->planes[_fi].normal, h->verts[_vi]) - h->planes[_fi].offset;
							if (d > worst) { worst = d; worst_vi = _vi; worst_fi = _fi; }
						}
					float worst_input = 0;
					int worst_pi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++)
						for (int _pi = 0; _pi < asize(pts); _pi++) {
							float d = dot(h->planes[_fi].normal, pts[_pi]) - h->planes[_fi].offset;
							if (d > worst_input) { worst_input = d; worst_pi = _pi; }
						}
					printf("    worst hull vert: v%d above f%d by %.8f\n", worst_vi, worst_fi, worst);
					printf("    worst input pt:  p%d outside by %.8f\n", worst_pi, worst_input);
					printf("    epsilon=%.2e maxoutside=%.2e\n", h->epsilon, h->maxoutside);
					// Dump input array as C code for reproduction.
					printf("    // --- copy-paste for unit test ---\n");
					printf("    v3 pts[] = {\n");
					for (int _pi = 0; _pi < asize(pts); _pi++)
						printf("        {%.9ef, %.9ef, %.9ef},\n",
							pts[_pi].x, pts[_pi].y, pts[_pi].z);
					printf("    };\n");
					printf("    // --- end ---\n");
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

// ============================================================================
// Cylinder stress test for face merging.
//
// A discretized cylinder has N coplanar triangles on each cap and N coplanar
// quads (2 tris each) on the barrel. This hammers the merge code because:
//  - Caps must merge N triangles into a single N-gon
//  - Barrel quads are nearly coplanar with their neighbors at high N
//  - Arbitrary rotation prevents axis-aligned shortcuts
//
// After quickhull + merging, the ideal hull has N+2 faces (N barrel quads +
// 2 caps), but we accept any topologically valid convex result.

static quat quat_axis_angle(v3 axis, float angle)
{
	float s = sinf(angle * 0.5f);
	float c = cosf(angle * 0.5f);
	return (quat){ axis.x * s, axis.y * s, axis.z * s, c };
}

static void test_quickhull_cylinder()
{
	int segs[] = { 8, 12, 16, 24, 32, 64 };
	int nseg_counts = sizeof(segs) / sizeof(segs[0]);

	// Rotation/scale combos: {axis, angle, scaleXYZ}
	struct { v3 axis; float angle; v3 scl; const char* tag; } xforms[] = {
		{ {0,1,0},  0.0f,        {1,1,1},           "identity"   },
		{ {1,0,0},  1.5707963f,  {1,1,1},           "rot90X"     },
		{ {0,0,1},  0.7853982f,  {1,1,1},           "rot45Z"     },
		{ {0,1,0},  2.3561945f,  {1,1,1},           "rot135Y"    },
		{ {1,1,1},  1.0471976f,  {1,1,1},           "rot60diag"  },
		{ {0,1,0},  0.0f,        {0.01f,0.01f,0.01f}, "tiny"     },
		{ {0,1,0},  0.0f,        {100,100,100},     "large"      },
		{ {1,0,0},  0.3926991f,  {2,0.5f,1},        "aniso+rot"  },
		{ {0,0,1},  2.0943951f,  {0.1f,10,0.1f},    "pancake"    },
		{ {1,1,0},  0.5235988f,  {5,5,0.01f},       "flat disc"  },
	};
	int nxforms = sizeof(xforms) / sizeof(xforms[0]);

	float radius = 1.0f;
	float half_h = 1.0f;
	int total = 0, fails = 0;

	for (int si = 0; si < nseg_counts; si++) {
		int n = segs[si];
		for (int xi = 0; xi < nxforms; xi++) {
			char label[128];
			snprintf(label, sizeof(label), "cylinder n=%d %s", n, xforms[xi].tag);

			quat q = quat_axis_angle(norm(xforms[xi].axis), xforms[xi].angle);
			v3 s = xforms[xi].scl;

			// Generate cylinder points: ring on top cap, ring on bottom cap.
			CK_DYNA v3* pts = NULL;
			for (int i = 0; i < n; i++) {
				float theta = 2.0f * 3.14159265f * (float)i / (float)n;
				float cx = radius * cosf(theta);
				float cz = radius * sinf(theta);

				// Top cap vertex.
				v3 pt = V3(cx * s.x, half_h * s.y, cz * s.z);
				apush(pts, rotate(q, pt));

				// Bottom cap vertex.
				v3 pb = V3(cx * s.x, -half_h * s.y, cz * s.z);
				apush(pts, rotate(q, pb));
			}

			Hull* h = quickhull(pts, asize(pts));
			total++;

			if (!h) {
				snprintf(label, sizeof(label), "cylinder n=%d %s build", n, xforms[xi].tag);
				TEST_BEGIN(label);
				TEST_ASSERT(h != NULL);
				afree(pts);
				continue;
			}

			int e_ok = hull_check_euler(h);
			int t_ok = hull_check_twins(h);
			int f_ok = hull_check_face_loops(h);
			int n_ok = hull_check_normals_outward(h);
			int c_ok = hull_check_convex(h, h->epsilon);
			int i_ok = hull_check_contains_inputs(h, pts, asize(pts), h->epsilon + h->maxoutside);

			char buf[128];
			snprintf(buf, sizeof(buf), "%s euler", label);
			TEST_BEGIN(buf); TEST_ASSERT(e_ok);
			snprintf(buf, sizeof(buf), "%s twins", label);
			TEST_BEGIN(buf); TEST_ASSERT(t_ok);
			snprintf(buf, sizeof(buf), "%s loops", label);
			TEST_BEGIN(buf); TEST_ASSERT(f_ok);
			snprintf(buf, sizeof(buf), "%s normals", label);
			TEST_BEGIN(buf); TEST_ASSERT(n_ok);
			snprintf(buf, sizeof(buf), "%s convex", label);
			TEST_BEGIN(buf); TEST_ASSERT(c_ok);
			snprintf(buf, sizeof(buf), "%s inputs", label);
			TEST_BEGIN(buf); TEST_ASSERT(i_ok);

			if (!e_ok || !t_ok || !f_ok || !n_ok || !c_ok || !i_ok) {
				fails++;
				if (fails <= 5) {
					printf("  CYLINDER FAIL: %s V=%d E=%d F=%d [%s%s%s%s%s%s]\n", label,
						h->vert_count, h->edge_count, h->face_count,
						e_ok ? "" : "euler ",
						t_ok ? "" : "twins ",
						f_ok ? "" : "loops ",
						n_ok ? "" : "normals ",
						c_ok ? "" : "convex ",
						i_ok ? "" : "inputs ");
					printf("    v3 pts[] = {\n");
					for (int pi = 0; pi < asize(pts); pi++)
						printf("        {%.9ef, %.9ef, %.9ef},\n",
							pts[pi].x, pts[pi].y, pts[pi].z);
					printf("    };\n");
				}
			}

			hull_free(h);
			afree(pts);
		}
	}

	printf("  cylinder: %d hulls tested, %d failures\n", total, fails);
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

	test_quickhull_case783();
	test_quickhull_ico7037();
	test_quickhull_tet8098();
	test_quickhull_cylinder();

	// Quickhull fuzz: moderate count for regular runs.
	test_quickhull_fuzz(100); // use 1000+ for stress testing

	printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	if (test_fail > 0) printf("*** FAILURES ***\n");
}

// Soak test: infinite-loop fuzz with crash logging to file.
static void test_quickhull_soak()
{
	const int ITERS_PER_SHAPE = 200;
	uint32_t base_seed = (uint32_t)time(NULL);
	int round = 0;
	int total_hulls = 0;
	FILE* logf = NULL;

	printf("=== quickhull soak test (seed base %u) ===\n", base_seed);
	printf("    logging failures to crash_log.txt\n");
	printf("    press Ctrl+C to stop\n\n");
	fflush(stdout);

	for (;;) {
		uint32_t round_seed = base_seed ^ (uint32_t)(round * 2654435761u);
		fuzz_rng_state = round_seed;
		fuzz_rng_seed  = round_seed;

		for (int shape_idx = 0; shape_idx < (int)FUZZ_SHAPE_COUNT; shape_idx++) {
			const FuzzShape* base = &s_fuzz_shapes[shape_idx];

			for (int iter = 0; iter < ITERS_PER_SHAPE; iter++) {
				CK_DYNA v3* pts = NULL;

				float hull_scale = 1.0f;
				float sr = fuzz_rand();
				if      (sr < 0.1f) hull_scale = fuzz_rand_range(0.001f, 0.01f);
				else if (sr < 0.2f) hull_scale = fuzz_rand_range(10.0f, 1000.0f);

				float jitter = fuzz_rand_range(0.0f, 0.001f) * hull_scale;
				for (int i = 0; i < base->vert_count; i++) {
					v3 p = scale(base->verts[i], hull_scale);
					p = add(p, scale(fuzz_rand_dir(), jitter));
					apush(pts, p);
				}

				int extras = (int)(fuzz_rand() * 20);
				if (fuzz_rand() < 0.1f) extras += 50 + (int)(fuzz_rand() * 50);

				for (int i = 0; i < extras; i++) {
					float r = fuzz_rand();
					v3 p;
					if (r < 0.15f) {
						int vi = fuzz_rand_int(base->vert_count);
						float fuzz = fuzz_rand_range(1e-7f, 1e-2f) * hull_scale;
						p = add(scale(base->verts[vi], hull_scale), scale(fuzz_rand_dir(), fuzz));
					} else if (r < 0.20f) {
						p = scale(base->verts[fuzz_rand_int(base->vert_count)], hull_scale);
					} else if (r < 0.35f) {
						int ei = fuzz_rand_int(base->edge_count);
						v3 a = scale(base->verts[base->edges[ei][0]], hull_scale);
						v3 b = scale(base->verts[base->edges[ei][1]], hull_scale);
						float t = fuzz_rand();
						p = add(scale(a, 1-t), scale(b, t));
						if (fuzz_rand() < 0.5f) {
							float fuzz = fuzz_rand_range(1e-7f, 1e-3f) * hull_scale;
							p = add(p, scale(fuzz_rand_dir(), fuzz));
						}
					} else if (r < 0.45f) {
						int ei = fuzz_rand_int(base->edge_count);
						v3 a = scale(base->verts[base->edges[ei][0]], hull_scale);
						v3 b = scale(base->verts[base->edges[ei][1]], hull_scale);
						int cluster = 2 + (int)(fuzz_rand() * 5);
						for (int c = 0; c < cluster; c++) {
							float t = (float)(c + 1) / (float)(cluster + 1);
							float fuzz = fuzz_rand_range(0, 1e-6f) * hull_scale;
							apush(pts, add(add(scale(a, 1-t), scale(b, t)), scale(fuzz_rand_dir(), fuzz)));
						}
						continue;
					} else if (r < 0.60f) {
						int ti = fuzz_rand_int(base->tri_count);
						v3 a = scale(base->verts[base->tris[ti][0]], hull_scale);
						v3 b = scale(base->verts[base->tris[ti][1]], hull_scale);
						v3 c = scale(base->verts[base->tris[ti][2]], hull_scale);
						float u = fuzz_rand(), v = fuzz_rand();
						if (u + v > 1.0f) { u = 1.0f - u; v = 1.0f - v; }
						p = add(add(scale(a, 1-u-v), scale(b, u)), scale(c, v));
						if (fuzz_rand() < 0.3f) {
							v3 fn = norm(cross(sub(b, a), sub(c, a)));
							float fuzz = fuzz_rand_range(-1e-4f, 1e-4f) * hull_scale;
							p = add(p, scale(fn, fuzz));
						}
					} else if (r < 0.70f) {
						int ti = fuzz_rand_int(base->tri_count);
						v3 a = scale(base->verts[base->tris[ti][0]], hull_scale);
						v3 b = scale(base->verts[base->tris[ti][1]], hull_scale);
						v3 c = scale(base->verts[base->tris[ti][2]], hull_scale);
						v3 fn = norm(cross(sub(b, a), sub(c, a)));
						float plane_d = dot(fn, a);
						p = scale(fuzz_rand_dir(), hull_scale * 1.5f);
						p = sub(p, scale(fn, dot(fn, p) - plane_d));
					} else if (r < 0.85f) {
						int vi = fuzz_rand_int(base->vert_count);
						p = scale(scale(base->verts[vi], hull_scale), fuzz_rand_range(0.0f, 0.9f));
					} else {
						p = scale(fuzz_rand_dir(), fuzz_rand_range(0.5f, 2.0f) * hull_scale);
					}
					apush(pts, p);
				}

				Hull* h = quickhull(pts, asize(pts));
				total_hulls++;

				if (!h) {
					afree(pts);
					continue;
				}

				int e_ok = hull_check_euler(h);
				int t_ok = hull_check_twins(h);
				int f_ok = hull_check_face_loops(h);
				int n_ok = hull_check_normals_outward(h);
				int c_ok = hull_check_convex(h, h->epsilon);
				int i_ok = hull_check_contains_inputs(h, pts, asize(pts), h->epsilon + h->maxoutside);

				if (!e_ok || !t_ok || !f_ok || !n_ok || !c_ok || !i_ok) {
					printf("\n!!! SOAK FAILURE at round %d, %s #%d (seed %u, %d pts, scale=%.3g)\n",
						round, base->name, iter, round_seed, asize(pts), hull_scale);
					printf("    V=%d E=%d F=%d [%s%s%s%s%s%s]\n",
						h->vert_count, h->edge_count, h->face_count,
						e_ok ? "" : "euler ",
						t_ok ? "" : "twins ",
						f_ok ? "" : "loops ",
						n_ok ? "" : "normals ",
						c_ok ? "" : "convex ",
						i_ok ? "" : "inputs ");

					float worst = 0;
					int worst_vi = -1, worst_fi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++)
						for (int _vi = 0; _vi < h->vert_count; _vi++) {
							float d = dot(h->planes[_fi].normal, h->verts[_vi]) - h->planes[_fi].offset;
							if (d > worst) { worst = d; worst_vi = _vi; worst_fi = _fi; }
						}
					float worst_input = 0;
					int worst_pi = -1;
					for (int _fi = 0; _fi < h->face_count; _fi++)
						for (int _pi = 0; _pi < asize(pts); _pi++) {
							float d = dot(h->planes[_fi].normal, pts[_pi]) - h->planes[_fi].offset;
							if (d > worst_input) { worst_input = d; worst_pi = _pi; }
						}
					printf("    worst hull vert: v%d above f%d by %.8f\n", worst_vi, worst_fi, worst);
					printf("    worst input pt:  p%d outside by %.8f\n", worst_pi, worst_input);
					printf("    epsilon=%.2e maxoutside=%.2e\n", h->epsilon, h->maxoutside);

					// Log to file.
					logf = fopen("crash_log.txt", "a");
					if (logf) {
						fprintf(logf, "// SOAK FAILURE: round %d, %s #%d, seed %u\n",
							round, base->name, iter, round_seed);
						fprintf(logf, "// V=%d E=%d F=%d [%s%s%s%s%s%s]\n",
							h->vert_count, h->edge_count, h->face_count,
							e_ok ? "" : "euler ",
							t_ok ? "" : "twins ",
							f_ok ? "" : "loops ",
							n_ok ? "" : "normals ",
							c_ok ? "" : "convex ",
							i_ok ? "" : "inputs ");
						fprintf(logf, "// worst hull vert: v%d above f%d by %.8f\n", worst_vi, worst_fi, worst);
						fprintf(logf, "// worst input pt:  p%d outside by %.8f\n", worst_pi, worst_input);
						fprintf(logf, "// epsilon=%.2e maxoutside=%.2e\n", h->epsilon, h->maxoutside);
						fprintf(logf, "// %d points, scale=%.3g\n", asize(pts), hull_scale);
						fprintf(logf, "v3 pts[] = {\n");
						for (int _pi = 0; _pi < asize(pts); _pi++)
							fprintf(logf, "\t{%.9ef, %.9ef, %.9ef},\n",
								pts[_pi].x, pts[_pi].y, pts[_pi].z);
						fprintf(logf, "};\n\n");
						fclose(logf);
						logf = NULL;
					}

					printf("    -> logged to crash_log.txt\n");
					printf("    -> stopping soak test.\n");
					fflush(stdout);

					hull_free(h);
					afree(pts);
					return;
				}

				hull_free(h);
				afree(pts);
			}
		}

		round++;
		if (round % 10 == 0) {
			printf("  soak: round %d done (%d hulls so far, no failures)\n", round, total_hulls);
			fflush(stdout);
		}
	}
}
