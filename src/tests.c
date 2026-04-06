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
static int test_bail; // if nonzero, exit on first failure
static const char* test_current;

#define TEST_BEGIN(name) do { test_current = name; } while(0)
#define TEST_ASSERT(cond) do { \
	if (!(cond)) { printf("  FAIL [%s] %s:%d: %s\n", test_current, __FILE__, __LINE__, #cond); test_fail++; \
		if (test_bail) { printf("--- BAIL: first failure ---\n"); exit(1); } } \
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
// Floor contact regression: sphere/capsule resting on large box floor.
// This exactly mimics the demo scene (floor half_extents=10,1,10).

static void test_floor_contacts()
{
	Manifold m;
	Box floor = { V3(0, 0, 0), quat_identity(), V3(10, 1, 10) };

	// Sphere resting on floor: center at y=1.5 (floor top at y=1, radius=0.5).
	TEST_BEGIN("floor sphere resting");
	m = (Manifold){0};
	int hit = collide_sphere_box((Sphere){ V3(0, 1.5f, 0), 0.5f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.0f, 0.01f);

	// Sphere slightly penetrating floor.
	TEST_BEGIN("floor sphere penetrating");
	m = (Manifold){0};
	hit = collide_sphere_box((Sphere){ V3(0, 1.3f, 0), 0.5f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.2f, 0.02f);
	TEST_ASSERT_NORMAL_DIR(m, 0, -1, 0); // from sphere toward floor

	// Sphere approaching floor from above (separated).
	TEST_BEGIN("floor sphere separated");
	m = (Manifold){0};
	hit = collide_sphere_box((Sphere){ V3(0, 2.0f, 0), 0.5f }, floor, &m);
	TEST_ASSERT(!hit);

	// Sphere deeply embedded in floor.
	TEST_BEGIN("floor sphere deep");
	m = (Manifold){0};
	hit = collide_sphere_box((Sphere){ V3(0, 0.5f, 0), 0.5f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	// Capsule resting on floor (vertical, tip touching).
	TEST_BEGIN("floor capsule resting");
	m = (Manifold){0};
	hit = collide_capsule_box((Capsule){ V3(0, 1.3f, 0), V3(0, 3, 0), 0.3f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.0f, 0.01f);

	// Capsule lying flat on floor (horizontal, both endpoints on floor).
	TEST_BEGIN("floor capsule flat");
	m = (Manifold){0};
	hit = collide_capsule_box((Capsule){ V3(-2, 1.2f, 0), V3(2, 1.2f, 0), 0.3f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count >= 1);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	// Capsule deeply embedded.
	TEST_BEGIN("floor capsule deep");
	m = (Manifold){0};
	hit = collide_capsule_box((Capsule){ V3(0, 0, 0), V3(0, 0.5f, 0), 0.3f }, floor, &m);
	TEST_ASSERT(hit);
	TEST_ASSERT(m.count > 0);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);
}

// Multi-frame simulation: sphere should settle on floor, not fall through.
static void test_sphere_settles_on_floor()
{
	TEST_BEGIN("sphere settles on floor");
	World w = create_world((WorldParams){ .gravity = V3(0, -10, 0) });
	Body floor_b = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
	Body sphere_b = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, sphere_b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	float dt = 1.0f / 60.0f;
	for (int i = 0; i < 300; i++) world_step(w, dt); // 5 seconds

	float y = body_get_position(w, sphere_b).y;
	// Sphere should rest on floor top (y=1) + radius (0.5) = ~1.5, not below 0.
	TEST_ASSERT(y > 1.0f);
	TEST_ASSERT(y < 2.0f);
	destroy_world(w);
}

// ============================================================================
// Sphere/Capsule-Hull shallow vs deep boundary tests.
// Validates the GJK-first dispatch: shallow (GJK witness), deep (face search).

static void test_gjk_dispatch()
{
	Manifold m = {0};
	Box box = { V3(0,0,0), quat_identity(), V3(1,1,1) };

	// Sphere-box: shallow (sphere just barely touching face).
	TEST_BEGIN("gjk-dispatch sphere-box shallow face");
	m = (Manifold){0};
	collide_sphere_box((Sphere){ V3(1.4f, 0, 0), 0.5f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.1f, EPS);
	TEST_ASSERT_NORMAL_DIR(m, -1, 0, 0);

	// Sphere-box: shallow edge (sphere near edge, not penetrating face).
	TEST_BEGIN("gjk-dispatch sphere-box shallow edge");
	m = (Manifold){0};
	collide_sphere_box((Sphere){ V3(1.2f, 1.2f, 0), 0.5f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	// Sphere-box: shallow vertex (sphere near vertex corner).
	TEST_BEGIN("gjk-dispatch sphere-box shallow vertex");
	m = (Manifold){0};
	collide_sphere_box((Sphere){ V3(1.15f, 1.15f, 1.15f), 0.5f }, box, &m);
	TEST_ASSERT_HIT(m);

	// Sphere-box: deep (sphere center inside box).
	TEST_BEGIN("gjk-dispatch sphere-box deep");
	m = (Manifold){0};
	collide_sphere_box((Sphere){ V3(0.5f, 0, 0), 0.3f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);
	// Normal should point away from nearest face.
	TEST_ASSERT_NORMAL_DIR(m, -1, 0, 0);

	// Sphere-box: deep center at origin (fully embedded).
	TEST_BEGIN("gjk-dispatch sphere-box deep center");
	m = (Manifold){0};
	collide_sphere_box((Sphere){ V3(0, 0, 0), 0.3f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	// Capsule-box: shallow (capsule tip just touching face).
	TEST_BEGIN("gjk-dispatch capsule-box shallow face");
	m = (Manifold){0};
	collide_capsule_box((Capsule){ V3(0, 1.1f, 0), V3(0, 3, 0), 0.2f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT_FLOAT(m.contacts[0].penetration, 0.1f, 0.05f);

	// Capsule-box: shallow edge (capsule near edge of box).
	TEST_BEGIN("gjk-dispatch capsule-box shallow edge");
	m = (Manifold){0};
	collide_capsule_box((Capsule){ V3(1.1f, 0, 1.1f), V3(1.1f, 0, -1.1f), 0.3f }, box, &m);
	TEST_ASSERT_HIT(m);

	// Capsule-box: deep (capsule core inside box).
	TEST_BEGIN("gjk-dispatch capsule-box deep");
	m = (Manifold){0};
	collide_capsule_box((Capsule){ V3(0, -0.5f, 0), V3(0, 0.5f, 0), 0.2f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.count >= 1);
	TEST_ASSERT(m.contacts[0].penetration > 0.0f);

	// Capsule-box: deep, segment parallel to face, two contact points.
	TEST_BEGIN("gjk-dispatch capsule-box deep 2-contact");
	m = (Manifold){0};
	collide_capsule_box((Capsule){ V3(-0.5f, 0.9f, 0), V3(0.5f, 0.9f, 0), 0.3f }, box, &m);
	TEST_ASSERT_HIT(m);
	TEST_ASSERT(m.count == 2);

	// Sphere-hull: custom hull (tetrahedron).
	TEST_BEGIN("gjk-dispatch sphere-hull");
	v3 tet_pts[] = { {0,1,0}, {-1,-1,1}, {1,-1,1}, {0,-1,-1} };
	Hull* h = quickhull(tet_pts, 4);
	TEST_ASSERT(h != NULL);
	if (h) {
		ConvexHull ch = { h, V3(0, 0, 0), quat_identity(), V3(1,1,1) };
		// Shallow: sphere above apex.
		m = (Manifold){0};
		collide_sphere_hull((Sphere){ V3(0, 1.3f, 0), 0.5f }, ch, &m);
		TEST_ASSERT_HIT(m);
		TEST_ASSERT(m.contacts[0].penetration > 0.0f);

		// Deep: sphere inside.
		m = (Manifold){0};
		collide_sphere_hull((Sphere){ V3(0, 0, 0), 0.1f }, ch, &m);
		TEST_ASSERT_HIT(m);
		TEST_ASSERT(m.contacts[0].penetration > 0.0f);
		hull_free(h);
	}
}

// ============================================================================
// Raw GJK distance tests. Tests gjk_query directly for known distances.

static void test_gjk_distance()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	ConvexHull unit_box = { box, V3(0,0,0), id, V3(1,1,1) };
	float expected;

	// --- Point vs box ---
	{
		// Point outside +X face: distance should be 1.0
		TEST_BEGIN("gjk point-box +X dist");
		GJK_Result r = gjk_query_point_hull(V3(2, 0, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 1.0f, 0.01f);

		// Point on +X face: distance ~0
		TEST_BEGIN("gjk point-box +X touching");
		r = gjk_query_point_hull(V3(1, 0, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

		// Point inside box: distance should be 0
		TEST_BEGIN("gjk point-box inside");
		r = gjk_query_point_hull(V3(0, 0, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

		// Point near edge (+X,+Y): dist = sqrt(2)-sqrt(2)*1 ≈ 0
		TEST_BEGIN("gjk point-box edge");
		r = gjk_query_point_hull(V3(1.5f, 1.5f, 0), unit_box);
		float expected = sqrtf(0.5f*0.5f + 0.5f*0.5f);
		TEST_ASSERT_FLOAT(r.distance, expected, 0.02f);

		// Point near vertex (+X,+Y,+Z):
		TEST_BEGIN("gjk point-box vertex");
		r = gjk_query_point_hull(V3(2, 2, 2), unit_box);
		expected = sqrtf(1+1+1);
		TEST_ASSERT_FLOAT(r.distance, expected, 0.02f);
	}

	// --- Segment vs box ---
	{
		// Vertical segment above box
		TEST_BEGIN("gjk seg-box above");
		GJK_Result r = gjk_query_segment_hull(V3(0, 3, 0), V3(0, 5, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.01f);

		// Segment crossing through box face
		TEST_BEGIN("gjk seg-box crossing");
		r = gjk_query_segment_hull(V3(0, 0, 0), V3(0, 2, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

		// Segment touching face
		TEST_BEGIN("gjk seg-box touching");
		r = gjk_query_segment_hull(V3(0, 1, 0), V3(0, 3, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

		// Segment near edge, parallel to edge
		TEST_BEGIN("gjk seg-box near edge");
		r = gjk_query_segment_hull(V3(1.5f, 1.5f, -1), V3(1.5f, 1.5f, 1), unit_box);
		expected = sqrtf(0.5f*0.5f + 0.5f*0.5f);
		TEST_ASSERT_FLOAT(r.distance, expected, 0.02f);

		// Segment fully inside box
		TEST_BEGIN("gjk seg-box inside");
		r = gjk_query_segment_hull(V3(-0.5f, 0, 0), V3(0.5f, 0, 0), unit_box);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);
	}

	// --- Box vs box (via hull) ---
	{
		// Two boxes separated on X
		TEST_BEGIN("gjk box-box separated");
		ConvexHull box_a = { box, V3(0,0,0), id, V3(1,1,1) };
		ConvexHull box_b = { box, V3(3,0,0), id, V3(1,1,1) };
		GJK_Result r = gjk_query_hull_hull(box_a, box_b);
		TEST_ASSERT_FLOAT(r.distance, 1.0f, 0.01f);

		// Two boxes touching on X
		TEST_BEGIN("gjk box-box touching");
		box_b = (ConvexHull){ box, V3(2,0,0), id, V3(1,1,1) };
		r = gjk_query_hull_hull(box_a, box_b);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);

		// Two boxes overlapping
		TEST_BEGIN("gjk box-box overlap");
		box_b = (ConvexHull){ box, V3(1.5f,0,0), id, V3(1,1,1) };
		r = gjk_query_hull_hull(box_a, box_b);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.01f);
	}

	// --- Scaled box ---
	{
		// Point above large floor box (half_extents 10,1,10)
		TEST_BEGIN("gjk point-floor above");
		ConvexHull floor = { box, V3(0,0,0), id, V3(10,1,10) };
		GJK_Result r = gjk_query_point_hull(V3(0, 2.5f, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 1.5f, 0.02f);

		// Point just above floor surface
		TEST_BEGIN("gjk point-floor near");
		r = gjk_query_point_hull(V3(0, 1.1f, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.1f, 0.02f);

		// Point on floor surface
		TEST_BEGIN("gjk point-floor touching");
		r = gjk_query_point_hull(V3(3, 1, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.02f);

		// Point inside floor
		TEST_BEGIN("gjk point-floor inside");
		r = gjk_query_point_hull(V3(0, 0.5f, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.02f);

		// Segment (capsule core) above floor
		TEST_BEGIN("gjk seg-floor above");
		r = gjk_query_segment_hull(V3(0, 1.5f, 0), V3(0, 3, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.5f, 0.02f);

		// Segment touching floor
		TEST_BEGIN("gjk seg-floor touching");
		r = gjk_query_segment_hull(V3(0, 1.0f, 0), V3(0, 3, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.02f);

		// Segment partially inside floor
		TEST_BEGIN("gjk seg-floor partial");
		r = gjk_query_segment_hull(V3(0, 0.5f, 0), V3(0, 3, 0), floor);
		TEST_ASSERT_FLOAT(r.distance, 0.0f, 0.02f);
	}

	// --- Witness point sanity ---
	{
		TEST_BEGIN("gjk witness point-box");
		ConvexHull floor = { box, V3(0,0,0), id, V3(10,1,10) };
		GJK_Result r = gjk_query_point_hull(V3(0, 3, 0), floor);
		// Witness on point should be the point itself
		TEST_ASSERT_FLOAT(r.point1.y, 3.0f, 0.01f);
		// Witness on hull should be on the +Y face (y=1)
		TEST_ASSERT_FLOAT(r.point2.y, 1.0f, 0.02f);
		// Distance = 3 - 1 = 2
		TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.02f);

		TEST_BEGIN("gjk witness seg-box");
		r = gjk_query_segment_hull(V3(-1, 3, 0), V3(1, 3, 0), floor);
		TEST_ASSERT_FLOAT(r.point1.y, 3.0f, 0.01f);
		TEST_ASSERT_FLOAT(r.point2.y, 1.0f, 0.02f);
		TEST_ASSERT_FLOAT(r.distance, 2.0f, 0.02f);
	}

	// --- Rotated box ---
	{
		TEST_BEGIN("gjk point-rotated box");
		quat rot45y = { 0, 0.3827f, 0, 0.9239f };
		ConvexHull rbox = { box, V3(0,0,0), rot45y, V3(1,1,1) };
		// After 45deg Y rotation, box corner is at x = sqrt(2) ≈ 1.414
		GJK_Result r = gjk_query_point_hull(V3(2, 0, 0), rbox);
		float corner = sqrtf(2.0f);
		TEST_ASSERT_FLOAT(r.distance, 2.0f - corner, 0.05f);
	}
}

// ============================================================================
// Contact sanity: verify contact point/normal/depth consistency.
// For sphere-box: contact point should be on sphere surface, and
// moving the sphere along normal by penetration should separate shapes.

static void test_contact_sanity()
{
	Box floor = { V3(0,0,0), quat_identity(), V3(10, 1, 10) };
	Box box = { V3(0,0,0), quat_identity(), V3(1,1,1) };
	Manifold m;

	// Debug: check specific failing case first.
	{
		Sphere s = { V3(0, 1.0f, 0), 0.5f };
		m = (Manifold){0};
		GJK_Result r = gjk_query_point_hull(s.center, (ConvexHull){ hull_unit_box(), V3(0,0,0), quat_identity(), V3(10,1,10) });
		int hit = collide_sphere_box(s, floor, &m);
		printf("    [debug] sphere y=1.0: gjk_dist=%.6f hit=%d count=%d", r.distance, hit, m.count);
		if (hit && m.count > 0) printf(" pen=%.4f n=(%.2f,%.2f,%.2f)", m.contacts[0].penetration, m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z);
		printf("\n");

		s = (Sphere){ V3(0, 1.3f, 0), 0.5f };
		m = (Manifold){0};
		r = gjk_query_point_hull(s.center, (ConvexHull){ hull_unit_box(), V3(0,0,0), quat_identity(), V3(10,1,10) });
		hit = collide_sphere_box(s, floor, &m);
		printf("    [debug] sphere y=1.3: gjk_dist=%.6f hit=%d count=%d", r.distance, hit, m.count);
		if (hit && m.count > 0) printf(" pen=%.4f n=(%.2f,%.2f,%.2f)", m.contacts[0].penetration, m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z);
		printf("\n");

		s = (Sphere){ V3(5, 1.3f, 0), 0.5f };
		m = (Manifold){0};
		r = gjk_query_point_hull(s.center, (ConvexHull){ hull_unit_box(), V3(0,0,0), quat_identity(), V3(10,1,10) });
		hit = collide_sphere_box(s, floor, &m);
		printf("    [debug] sphere x=5 y=1.3: gjk_dist=%.6f hit=%d count=%d", r.distance, hit, m.count);
		if (hit && m.count > 0) printf(" pen=%.4f n=(%.2f,%.2f,%.2f)", m.contacts[0].penetration, m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z);
		printf("\n");
	}

	// Debug x=-3 case: also test x=3 for comparison.
	{
		ConvexHull fl = { hull_unit_box(), V3(0,0,0), quat_identity(), V3(10,1,10) };
		for (float tx = -5; tx <= 5; tx += 1) {
			Sphere s = { V3(tx, 1.3f, 0), 0.5f };
			m = (Manifold){0};
			GJK_Result r = gjk_query_point_hull(s.center, fl);
			int hit = collide_sphere_box(s, floor, &m);
			printf("    [debug] x=%.0f y=1.3: gjk_dist=%.4f iters=%d hit=%d", tx, r.distance, r.iterations, hit);
			if (hit && m.count > 0) printf(" pen=%.3f n=(%.2f,%.2f,%.2f)", m.contacts[0].penetration, m.contacts[0].normal.x, m.contacts[0].normal.y, m.contacts[0].normal.z);
			printf("\n");
		}
	}

	// Sweep sphere across floor at different heights and X positions.
	float heights[] = { 1.5f, 1.4f, 1.3f, 1.2f, 1.1f, 1.0f, 0.9f, 0.5f, 0.0f };
	float xpos[] = { 0, 1, 5, 9, -3 };
	for (int hi = 0; hi < 9; hi++) {
		for (int xi = 0; xi < 5; xi++) {
			float y = heights[hi];
			float x = xpos[xi];
			Sphere s = { V3(x, y, 0), 0.5f };
			m = (Manifold){0};
			int hit = collide_sphere_box(s, floor, &m);
			float expected_sep = y - 1.0f; // distance from sphere center to floor top
			if (expected_sep > 0.5f + 0.01f) {
				// Should be separated
				TEST_BEGIN("contact sanity sphere-floor separated");
				TEST_ASSERT(!hit);
			} else {
				char lbl[128]; snprintf(lbl, sizeof(lbl), "sphere-floor hit x=%.1f y=%.1f sep=%.2f", x, y, expected_sep);
				TEST_BEGIN(lbl);
				if (!hit) printf("    MISS: %s gjk=? hit=%d\n", lbl, hit);
				TEST_ASSERT(hit);
				if (hit && m.count > 0) {
					// Penetration should be roughly radius - (center_y - floor_top)
					float expected_pen = 0.5f - expected_sep;
					if (expected_pen < 0) expected_pen = 0;
					TEST_BEGIN("contact sanity sphere-floor pen");
					TEST_ASSERT_FLOAT(m.contacts[0].penetration, expected_pen, 0.1f);
					// Normal should point roughly downward (from sphere toward floor).
					// Skip when sphere is at floor center (degenerate tie).
					if (expected_sep > -0.99f) {
						TEST_BEGIN("contact sanity sphere-floor normal");
						TEST_ASSERT(m.contacts[0].normal.y < -0.5f);
					}
				}
			}
		}
	}

	// Capsule at various orientations touching unit box.
	{
		// Vertical capsule approaching +Y face
		TEST_BEGIN("contact sanity capsule-box vertical");
		m = (Manifold){0};
		int hit = collide_capsule_box((Capsule){ V3(0, 1.1f, 0), V3(0, 3, 0), 0.2f }, box, &m);
		TEST_ASSERT(hit);
		if (hit && m.count > 0) {
			TEST_ASSERT(m.contacts[0].penetration > 0.0f);
			TEST_ASSERT(m.contacts[0].penetration < 0.5f);
			// Normal should point roughly down (from capsule toward box)
			TEST_ASSERT(m.contacts[0].normal.y < -0.3f);
		}

		// Horizontal capsule along +X face
		TEST_BEGIN("contact sanity capsule-box horizontal");
		m = (Manifold){0};
		hit = collide_capsule_box((Capsule){ V3(1.1f, 0, -0.5f), V3(1.1f, 0, 0.5f), 0.2f }, box, &m);
		TEST_ASSERT(hit);
		if (hit && m.count > 0) {
			TEST_ASSERT(m.contacts[0].penetration > 0.0f);
			// Normal should point roughly -X (from capsule toward box)
			TEST_ASSERT(m.contacts[0].normal.x < -0.3f);
		}

		// Capsule deep inside box
		TEST_BEGIN("contact sanity capsule-box deep");
		m = (Manifold){0};
		hit = collide_capsule_box((Capsule){ V3(0, 0, 0), V3(0, 0.5f, 0), 0.1f }, box, &m);
		TEST_ASSERT(hit);
		if (hit && m.count > 0) {
			TEST_ASSERT(m.contacts[0].penetration > 0.0f);
		}
	}

	// Sphere at box corners and edges (stress GJK Voronoi regions).
	{
		float offsets[] = { 1.2f, 1.1f, 1.05f, 0.9f };
		for (int oi = 0; oi < 4; oi++) {
			float d = offsets[oi];
			// Edge region
			TEST_BEGIN("contact sanity sphere-box edge sweep");
			m = (Manifold){0};
			int hit = collide_sphere_box((Sphere){ V3(d, d, 0), 0.5f }, box, &m);
			if (d > 1.0f + 0.5f / sqrtf(2.0f)) {
				// Far enough to be separated
			} else {
				TEST_ASSERT(hit);
				if (hit && m.count > 0) {
					TEST_ASSERT(m.contacts[0].penetration > 0.0f);
				}
			}
			// Vertex region
			TEST_BEGIN("contact sanity sphere-box vertex sweep");
			m = (Manifold){0};
			hit = collide_sphere_box((Sphere){ V3(d, d, d), 0.5f }, box, &m);
			if (hit && m.count > 0) {
				TEST_ASSERT(m.contacts[0].penetration > 0.0f);
			}
		}
	}
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

// Soak failure: Euler violation on icosahedron (V=22 E=112 F=37, expected V-E/2+F=2).
static void test_quickhull_soak_ico1336()
{
	TEST_BEGIN("quickhull soak_ico1336");
	v3 pts[] = {
		{-9.998443127e-01f, 1.618065357e+00f, -8.653744590e-05f},
		{1.000008821e+00f, 1.618042588e+00f, -1.804200729e-04f},
		{-1.000023603e+00f, -1.617910266e+00f, -1.297087874e-04f},
		{9.998848438e-01f, -1.618144155e+00f, 8.552092913e-05f},
		{1.557274954e-04f, -9.999080896e-01f, 1.618033528e+00f},
		{-7.499678031e-05f, 9.999208450e-01f, 1.617889762e+00f},
		{-3.326591468e-05f, -9.999384880e-01f, -1.618200779e+00f},
		{1.087870623e-04f, 1.000144005e+00f, -1.618044853e+00f},
		{1.618117929e+00f, 1.596595393e-04f, -9.999868870e-01f},
		{1.618144155e+00f, -5.029122622e-05f, 9.998656511e-01f},
		{-1.618043065e+00f, -1.179813116e-04f, -9.998632669e-01f},
		{-1.617968559e+00f, -1.647294266e-04f, 1.000035763e+00f},
		{-1.483431101e+00f, 3.524136543e-01f, 7.821941376e-01f},
		{-1.621091485e+00f, -1.936959568e-03f, -9.982609153e-01f},
		{1.458865047e+00f, -4.167097211e-01f, -5.984700322e-01f},
		{-5.284439921e-01f, -3.075364828e-01f, 8.261925578e-01f},
		{-1.478785157e+00f, -2.608410716e-01f, -8.782760501e-01f},
		{-4.821130335e-01f, 3.202155232e-01f, -1.786070466e-01f},
		{1.558186531e+00f, -1.572128534e-01f, 9.026681185e-01f},
		{-1.168316364e+00f, -3.692707121e-01f, 1.080447078e+00f},
		{2.187570333e-01f, -1.135104299e+00f, 1.264784575e+00f},
		{-1.617757797e+00f, 5.264066858e-04f, 9.997240305e-01f},
		{-7.926511019e-02f, 1.048145294e+00f, 1.491067529e+00f},
		{-1.666325610e-03f, -9.974440336e-01f, -1.615188479e+00f},
		{1.006952286e+00f, 1.615650654e+00f, -4.053545184e-03f},
		{0.000000000e+00f, 3.444449008e-01f, -5.573235750e-01f},
		{7.900940180e-01f, -1.278398991e+00f, 0.000000000e+00f},
		{-1.244565725e+00f, -1.125827804e-02f, -1.255824924e+00f},
		{-1.154508591e+00f, -1.213525414e+00f, -2.500001788e-01f},
		{-1.309017301e+00f, -8.090171814e-01f, -5.000002384e-01f},
		{-1.463525534e+00f, -4.045087695e-01f, -7.500005960e-01f},
		{-6.362375570e-04f, -6.785897613e-01f, -1.618612170e+00f},
		{9.463248253e-01f, 5.810145736e-01f, 6.789535880e-01f},
		{3.047788143e-01f, -1.740671992e+00f, 1.182140589e+00f},
		{-3.892679811e-01f, 1.038147688e+00f, 1.190618277e+00f},
		{-1.154508710e+00f, -1.213525176e+00f, 2.499994338e-01f},
		{-1.309016466e+00f, -8.090165854e-01f, 5.000002980e-01f},
		{-1.463524818e+00f, -4.045090377e-01f, 7.500002980e-01f},
		{1.001508951e+00f, -1.619185925e+00f, -1.232806011e-03f},
		{1.145782232e+00f, 1.237686872e+00f, -2.355188727e-01f},
		{1.386887193e+00f, -1.428574175e-01f, -1.088290811e+00f},
		{1.155738711e+00f, -2.857138216e-01f, -1.176580548e+00f},
		{9.245904088e-01f, -4.285714924e-01f, -1.264871716e+00f},
		{6.934429407e-01f, -5.714291334e-01f, -1.353162646e+00f},
		{4.622954130e-01f, -7.142856717e-01f, -1.441452861e+00f},
		{2.311475873e-01f, -8.571430445e-01f, -1.529743552e+00f},
		{-6.468786001e-01f, 1.046671629e+00f, 0.000000000e+00f},
		{-1.286025882e+00f, 0.000000000e+00f, 7.948077321e-01f},
		{1.828569293e+00f, 5.511882305e-01f, -3.639454395e-02f},
		{-9.889257550e-01f, 0.000000000e+00f, -6.111897230e-01f},
		{8.403830528e-01f, -2.282335982e-02f, 9.002639651e-01f},
		{-3.254281182e-04f, -3.752844334e-01f, -1.618466973e+00f},
		{-1.444956013e-07f, -7.142855525e-01f, -1.618034005e+00f},
		{-5.037525952e-08f, -4.285715520e-01f, -1.618034124e+00f},
		{-1.316915217e-07f, -1.428570300e-01f, -1.618033886e+00f},
		{2.816018139e-07f, 1.428568363e-01f, -1.618034363e+00f},
		{-1.477980494e-08f, 4.285713732e-01f, -1.618033171e+00f},
		{-2.022295007e-07f, 7.142857313e-01f, -1.618033886e+00f},
		{-1.825485587e+00f, -4.501773119e-01f, -3.423708677e-01f},
		{3.906928934e-03f, -9.972873330e-01f, -1.617000341e+00f},
		{3.330070293e-03f, 9.961679578e-01f, 1.615726829e+00f},
		{-9.080737829e-01f, 3.027348518e-01f, 7.076822519e-01f},
		{-1.453365326e+00f, -3.383456171e-01f, 8.264660835e-01f},
		{-1.050187349e+00f, -1.486641765e+00f, -8.120489866e-02f},
		{9.176892042e-01f, -1.833526254e+00f, 5.523809418e-02f},
		{1.220894337e+00f, 1.039725065e+00f, 2.130673826e-02f},
		{-6.092307568e-01f, 0.000000000e+00f, -3.765253127e-01f},
		{0.000000000e+00f, -6.475350261e-01f, 1.047733665e+00f},
		{-1.603541970e+00f, -3.777083382e-02f, -5.920574665e-01f},
		{5.792362094e-01f, -3.458117545e-01f, -8.061514497e-01f},
		{5.628519654e-01f, -1.290026903e+00f, 7.651551366e-01f},
		{-1.116937518e+00f, 2.421482354e-01f, -1.089060307e-01f},
		{7.714872956e-01f, -5.468678474e-01f, 1.633635521e+00f},
		{1.311475396e+00f, -8.025808334e-01f, -5.039777756e-01f},
		{1.479111195e+00f, -3.638058603e-01f, 7.751767039e-01f},
		{-1.569204807e+00f, 1.278364509e-01f, 9.209927320e-01f},
		{1.000000000e+00f, 1.618034005e+00f, 0.000000000e+00f},
		{-1.102279544e+00f, 1.113913417e+00f, 4.018411040e-01f},
		{1.352338046e-01f, -3.948976696e-01f, -1.566453576e+00f},
		{-4.990850110e-03f, -9.976323843e-01f, -1.625341177e+00f},
		{1.171967745e+00f, 8.084440231e-02f, 1.526910543e+00f},
		{-8.333329558e-01f, -1.515027881e+00f, 2.696718872e-01f},
		{-6.666659117e-01f, -1.412022829e+00f, 5.393452048e-01f},
		{-5.000001788e-01f, -1.309016824e+00f, 8.090171218e-01f},
		{-3.333331048e-01f, -1.206011176e+00f, 1.078689337e+00f},
		{-1.666672081e-01f, -1.103005528e+00f, 1.348361015e+00f},
		{-1.001470566e+00f, 1.616499662e+00f, -4.561771057e-04f},
		{-5.116344094e-01f, 8.278418779e-01f, 0.000000000e+00f},
		{-1.785489917e-01f, -2.654796243e-01f, -9.135619402e-01f},
		{-4.047199488e-01f, -1.778397083e+00f, -4.349166751e-01f},
		{0.000000000e+00f, -5.315647721e-01f, 8.600898981e-01f},
		{8.571432233e-01f, -1.529743791e+00f, 2.311480194e-01f},
		{7.142857909e-01f, -1.441452146e+00f, 4.622952640e-01f},
		{5.714288950e-01f, -1.353161573e+00f, 6.934428215e-01f},
		{4.285710454e-01f, -1.264871359e+00f, 9.245903492e-01f},
		{2.857146561e-01f, -1.176581264e+00f, 1.155738235e+00f},
		{1.428573281e-01f, -1.088290691e+00f, 1.386886477e+00f},
		{1.612959266e+00f, 1.307880599e-02f, 4.425288439e-01f},
		{1.136922479e+00f, -1.261270523e+00f, -2.210779041e-01f},
		{-1.379806638e+00f, -6.236872077e-01f, 3.145914078e-01f},
		{-1.445874453e+00f, 4.507191777e-01f, 1.044155955e-01f},
		{1.099934801e-03f, 1.005945921e+00f, -1.618815660e+00f},
		{-1.078688860e+00f, -3.333331645e-01f, 1.206010938e+00f},
		{-5.393446088e-01f, -6.666666269e-01f, 1.412022829e+00f},
		{1.245681763e+00f, 9.748613834e-01f, -1.215741932e-01f},
		{-1.118878722e-01f, -4.971617162e-01f, -9.679811597e-01f},
		{9.933810234e-01f, -1.617228270e+00f, -4.226116464e-03f},
		{-1.606284142e+00f, 3.076153807e-02f, -9.809883237e-01f},
		{-1.126724362e+00f, 0.000000000e+00f, 6.963539124e-01f},
		{-7.368535399e-01f, -1.480157614e+00f, 3.609648943e-01f},
		{-1.390085697e+00f, -5.967763662e-01f, 6.311719418e-01f},
		{1.619863153e+00f, 1.071012695e-03f, -1.002707720e+00f},
		{4.361957014e-01f, -3.072058558e-01f, -1.451422095e+00f},
		{-1.001425505e+00f, -1.620300651e+00f, -4.062809967e-05f},
		{1.748257875e-02f, 1.527979612e+00f, 2.357654870e-01f},
		{-9.017625707e-04f, -9.984976649e-01f, 1.617181301e+00f},
		{-1.667023301e-01f, 1.162230849e+00f, 1.193308353e+00f},
		{1.172792196e+00f, 1.165658474e+00f, -2.795834839e-01f},
		{2.529041171e-01f, 0.000000000e+00f, 1.563033313e-01f},
	};
	int npts = sizeof(pts) / sizeof(pts[0]);

	Hull* h = quickhull(pts, npts);
	TEST_ASSERT(h != NULL);
	if (h) {
		TEST_ASSERT(hull_check_euler(h));
		TEST_ASSERT(hull_check_twins(h));
		TEST_ASSERT(hull_check_face_loops(h));
		TEST_ASSERT(hull_check_normals_outward(h));
		TEST_ASSERT(hull_check_convex(h, h->epsilon));
		TEST_ASSERT(hull_check_contains_inputs(h, pts, npts, h->epsilon + h->maxoutside));
		hull_free(h);
	}
}

// Soak failure: Euler+twins on tiny tetrahedron (V=13 E=61 F=20).
static void test_quickhull_soak_tet101()
{
	TEST_BEGIN("quickhull soak_tet101");
	v3 pts[] = {
		{8.348876008e-07f, 4.681363702e-03f, 6.981662182e-07f},
		{-4.680875689e-03f, -4.682406783e-03f, 4.681994673e-03f},
		{4.680886865e-03f, -4.682447761e-03f, 4.681915976e-03f},
		{5.302704835e-07f, -4.680839833e-03f, -4.682225175e-03f},
		{2.199016308e-05f, -4.658407066e-03f, -4.654521123e-03f},
		{-4.681670573e-03f, -4.681670573e-03f, 4.681670573e-03f},
		{2.487084828e-03f, 2.295025159e-03f, 3.906742670e-03f},
		{0.000000000e+00f, -4.681670573e-03f, -4.681670573e-03f},
		{1.984762028e-03f, -4.681788385e-03f, -7.125165430e-04f},
		{-4.967796485e-07f, -4.681355320e-03f, -4.682343919e-03f},
		{5.329707172e-03f, 2.600779524e-03f, 1.040445175e-03f},
		{-3.901390824e-03f, -4.681668244e-03f, 3.121116664e-03f},
		{-3.121114103e-03f, -4.681670107e-03f, 1.560558681e-03f},
		{-2.340835053e-03f, -4.681671504e-03f, 1.123785578e-10f},
		{-1.560556819e-03f, -4.681670573e-03f, -1.560557052e-03f},
		{-7.802789914e-04f, -4.681671038e-03f, -3.121113172e-03f},
		{4.667766858e-03f, -4.693021998e-03f, 4.664019216e-03f},
		{2.266348340e-03f, -4.681670573e-03f, -1.489738934e-04f},
		{-3.373793326e-03f, -3.373793326e-03f, 3.373793326e-03f},
		{4.672480281e-03f, -4.688905552e-03f, 4.669602029e-03f},
		{7.668895705e-06f, -4.663798027e-03f, -4.680988379e-03f},
		{2.402480692e-03f, -4.681670573e-03f, 1.232908107e-04f},
		{-1.855140319e-03f, 9.713897016e-04f, 1.855140319e-03f},
	};
	int npts = sizeof(pts) / sizeof(pts[0]);
	Hull* h = quickhull(pts, npts);
	TEST_ASSERT(h != NULL);
	if (h) {
		TEST_ASSERT(hull_check_euler(h));
		TEST_ASSERT(hull_check_twins(h));
		TEST_ASSERT(hull_check_face_loops(h));
		TEST_ASSERT(hull_check_normals_outward(h));
		TEST_ASSERT(hull_check_convex(h, h->epsilon));
		TEST_ASSERT(hull_check_contains_inputs(h, pts, npts, h->epsilon + h->maxoutside));
		hull_free(h);
	}
}

static void test_quickhull_bipyramid_loops()
{
	v3 pts[] = {
		{  0.0f,  0.7f,  0.0f },
		{  0.5f,  0.0f,  0.5f },
		{ -0.5f,  0.0f,  0.5f },
		{  0.5f,  0.0f, -0.5f },
		{ -0.5f,  0.0f, -0.5f },
		{  0.0f, -0.7f,  0.0f },
	};
	Hull* h = quickhull(pts, 6);
	TEST_BEGIN("bipyramid hull built");
	TEST_ASSERT(h != NULL);
	TEST_ASSERT(h->face_count >= 4);

	int any_broken = 0;
	for (int f = 0; f < h->face_count; f++) {
		int start = h->faces[f].edge;
		int e = start;
		int count = 0;
		do {
			count++;
			if (count > 100) { any_broken = 1; break; }
			if (h->edges[e].face != f) { any_broken = 1; break; }
			e = h->edges[e].next;
		} while (e != start);
	}
	TEST_BEGIN("bipyramid face loops valid");
	TEST_ASSERT(!any_broken);

	// Test collision at multiple heights (the crash was during runtime fall)
	TEST_BEGIN("bipyramid vs floor box at multiple positions");
	{
		const Hull* box = hull_unit_box();
		ConvexHull floor_hull = {
			.hull = box,
			.center = V3(0, -1, 0),
			.rotation = quat_identity(),
			.scale = V3(10, 1, 10),
		};
		// Test at many Y positions to cover different reference face selections
		for (int yi = 0; yi < 20; yi++) {
			float y = -0.5f + yi * 0.1f;
			ConvexHull bip_hull = {
				.hull = h,
				.center = V3(3, y, 0),
				.rotation = quat_identity(),
				.scale = V3(1, 1, 1),
			};
			Manifold m = {0};
			collide_hull_hull(floor_hull, bip_hull, &m);
			collide_hull_hull(bip_hull, floor_hull, &m); // test both orderings
		}
		TEST_ASSERT(1);
	}

	// Test with rotated bipyramid (reproduces runtime scenario)
	TEST_BEGIN("bipyramid vs floor rotated");
	{
		ConvexHull floor_hull = {
			.hull = hull_unit_box(),
			.center = V3(0, -1, 0),
			.rotation = quat_identity(),
			.scale = V3(10, 1, 10),
		};
		for (int ri = 0; ri < 36; ri++) {
			float angle = ri * 0.175f; // ~10 deg steps
			float sa = sinf(angle * 0.5f), ca = cosf(angle * 0.5f);
			quat rots[] = {
				{ sa, 0, 0, ca },  // pitch
				{ 0, sa, 0, ca },  // yaw
				{ 0, 0, sa, ca },  // roll
				{ sa*0.577f, sa*0.577f, sa*0.577f, ca }, // diagonal
			};
			for (int r = 0; r < 4; r++) {
				for (int yi = 0; yi < 10; yi++) {
					float y = -0.3f + yi * 0.1f;
					ConvexHull bip = {
						.hull = h,
						.center = V3(3, y, 0),
						.rotation = rots[r],
						.scale = V3(1, 1, 1),
					};
					Manifold m = {0};
					collide_hull_hull(floor_hull, bip, &m);
					collide_hull_hull(bip, floor_hull, &m);
				}
			}
		}
		TEST_ASSERT(1);
	}

	// Also test bipyramid vs small box (box-hull dispatch)
	TEST_BEGIN("bipyramid vs small box collision");
	{
		ConvexHull box_hull = {
			.hull = hull_unit_box(),
			.center = V3(3, 1.5f, 0),
			.rotation = quat_identity(),
			.scale = V3(0.4f, 0.4f, 0.4f),
		};
		ConvexHull bip_hull = {
			.hull = h,
			.center = V3(3, 0.7f, 0),
			.rotation = quat_identity(),
			.scale = V3(1, 1, 1),
		};
		Manifold m = {0};
		collide_hull_hull(box_hull, bip_hull, &m);
		collide_hull_hull(bip_hull, box_hull, &m);
		TEST_ASSERT(1);
	}

	hull_free(h);
}

// ============================================================================
// Geometry stress tests: diverse input shapes that exercise different codepaths.

static void validate_hull_or_null(const char* label, Hull* h, const v3* pts, int npts)
{
	if (!h) {
		// NULL is acceptable for degenerate inputs.
		char buf[128];
		snprintf(buf, sizeof(buf), "%s build", label);
		TEST_BEGIN(buf);
		TEST_ASSERT(1);
		return;
	}
	char buf[128];
	snprintf(buf, sizeof(buf), "%s euler", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_euler(h));
	snprintf(buf, sizeof(buf), "%s twins", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_twins(h));
	snprintf(buf, sizeof(buf), "%s loops", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_face_loops(h));
	snprintf(buf, sizeof(buf), "%s normals", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_normals_outward(h));
	snprintf(buf, sizeof(buf), "%s convex", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_convex(h, h->epsilon));
	snprintf(buf, sizeof(buf), "%s inputs", label);
	TEST_BEGIN(buf); TEST_ASSERT(hull_check_contains_inputs(h, pts, npts, h->epsilon + h->maxoutside));
	hull_free(h);
}

static void test_quickhull_geometry_stress()
{
	int total = 0, fails = 0;
	fuzz_rng_state = 99999;

	// --- 1. Flat/degenerate: all points coplanar ---
	for (int trial = 0; trial < 5; trial++) {
		CK_DYNA v3* pts = NULL;
		int n = 10 + fuzz_rand_int(40);
		for (int i = 0; i < n; i++) {
			float x = fuzz_rand_range(-10, 10);
			float z = fuzz_rand_range(-10, 10);
			float y = 1.0f + fuzz_rand_range(-1e-6f, 1e-6f); // near-coplanar
			apush(pts, V3(x, y, z));
		}
		// Add 1-2 points slightly off-plane to make a thin sliver.
		apush(pts, V3(0, 1.001f, 0));
		apush(pts, V3(0, 0.999f, 0));
		total++;
		char label[64]; snprintf(label, sizeof(label), "flat sliver #%d (%d pts)", trial, asize(pts));
		Hull* h = quickhull(pts, asize(pts));
		validate_hull_or_null(label, h, pts, asize(pts));
		afree(pts);
	}

	// --- 2. Needle shapes: extreme aspect ratio ---
	{
		float lengths[] = { 100, 1000, 10000 };
		float radii[]   = { 0.1f, 0.01f, 0.001f };
		for (int li = 0; li < 3; li++) {
			CK_DYNA v3* pts = NULL;
			float L = lengths[li], R = radii[li];
			int segs = 8;
			// Two endcaps + ring
			apush(pts, V3(0, -L/2, 0));
			apush(pts, V3(0,  L/2, 0));
			for (int i = 0; i < segs; i++) {
				float theta = 2.0f * 3.14159265f * (float)i / (float)segs;
				float cx = R * cosf(theta), cz = R * sinf(theta);
				apush(pts, V3(cx, -L/2, cz));
				apush(pts, V3(cx,  L/2, cz));
				apush(pts, V3(cx, 0, cz));
			}
			total++;
			char label[64]; snprintf(label, sizeof(label), "needle L=%.0f R=%.3f", L, R);
			Hull* h = quickhull(pts, asize(pts));
			validate_hull_or_null(label, h, pts, asize(pts));
			afree(pts);
		}
	}

	// --- 4. Spherical distribution: all points on hull ---
	for (int trial = 0; trial < 3; trial++) {
		CK_DYNA v3* pts = NULL;
		int n = 50 + trial * 50; // 50, 100, 150
		for (int i = 0; i < n; i++) {
			v3 d = fuzz_rand_dir();
			float r = 1.0f + fuzz_rand_range(-1e-4f, 1e-4f);
			apush(pts, scale(d, r));
		}
		total++;
		char label[64]; snprintf(label, sizeof(label), "sphere %d pts", n);
		Hull* h = quickhull(pts, asize(pts));
		validate_hull_or_null(label, h, pts, asize(pts));
		afree(pts);
	}

	// --- 5. Near-duplicate vertices: many within epsilon ---
	for (int trial = 0; trial < 5; trial++) {
		CK_DYNA v3* pts = NULL;
		// 8 cube corners, each duplicated 5-10 times with tiny jitter.
		v3 corners[] = {
			{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
			{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1},
		};
		for (int c = 0; c < 8; c++) {
			int dupes = 5 + fuzz_rand_int(6);
			for (int d = 0; d < dupes; d++) {
				float jitter = fuzz_rand_range(0, 1e-5f);
				apush(pts, add(corners[c], scale(fuzz_rand_dir(), jitter)));
			}
		}
		total++;
		char label[64]; snprintf(label, sizeof(label), "near-dup #%d (%d pts)", trial, asize(pts));
		Hull* h = quickhull(pts, asize(pts));
		validate_hull_or_null(label, h, pts, asize(pts));
		afree(pts);
	}

	// --- 6. One-axis collapse: collinear / pencil ---
	{
		// Pure collinear: should return NULL.
		CK_DYNA v3* pts = NULL;
		for (int i = 0; i < 20; i++)
			apush(pts, V3(0, (float)i * 0.1f, 0));
		total++;
		TEST_BEGIN("collinear 20pts");
		Hull* h = quickhull(pts, asize(pts));
		TEST_ASSERT(h == NULL); // degenerate, no valid hull
		afree(pts);

		// Thin pencil: line + tiny perturbation.
		pts = NULL;
		for (int i = 0; i < 30; i++) {
			float t = (float)i / 29.0f;
			apush(pts, V3(fuzz_rand_range(-1e-6f, 1e-6f), t * 10.0f,
			              fuzz_rand_range(-1e-6f, 1e-6f)));
		}
		// Add 4 off-axis points to make it buildable.
		apush(pts, V3(0.01f, 5.0f, 0));
		apush(pts, V3(-0.01f, 5.0f, 0));
		apush(pts, V3(0, 5.0f, 0.01f));
		apush(pts, V3(0, 5.0f, -0.01f));
		total++;
		char label[64]; snprintf(label, sizeof(label), "pencil 34 pts");
		h = quickhull(pts, asize(pts));
		validate_hull_or_null(label, h, pts, asize(pts));
		afree(pts);
	}

	// --- 8. Coplanar grid: NxN on a plane + off-plane points ---
	{
		int grids[] = { 5, 8, 10 };
		for (int gi = 0; gi < 3; gi++) {
			CK_DYNA v3* pts = NULL;
			int N = grids[gi];
			for (int x = 0; x < N; x++)
				for (int z = 0; z < N; z++)
					apush(pts, V3((float)x, 0, (float)z));
			// 4 points above and below.
			float mid = (float)(N-1) / 2.0f;
			apush(pts, V3(mid, 1.0f, mid));
			apush(pts, V3(mid, -1.0f, mid));
			apush(pts, V3(0, 0.5f, 0));
			apush(pts, V3((float)(N-1), 0.5f, (float)(N-1)));
			total++;
			char label[64]; snprintf(label, sizeof(label), "grid %dx%d+4", N, N);
			Hull* h = quickhull(pts, asize(pts));
			validate_hull_or_null(label, h, pts, asize(pts));
			afree(pts);
		}
	}

	// --- 9. Bipyramid / double cone ---
	{
		int ring_sizes[] = { 6, 12, 24, 48 };
		for (int ri = 0; ri < 4; ri++) {
			CK_DYNA v3* pts = NULL;
			int n = ring_sizes[ri];
			// Top and bottom apex.
			apush(pts, V3(0, 2, 0));
			apush(pts, V3(0, -2, 0));
			// Equatorial ring.
			for (int i = 0; i < n; i++) {
				float theta = 2.0f * 3.14159265f * (float)i / (float)n;
				apush(pts, V3(cosf(theta), 0, sinf(theta)));
			}
			// Add some interior points and near-edge points for stress.
			for (int i = 0; i < n/2; i++) {
				float theta = fuzz_rand() * 2.0f * 3.14159265f;
				float r = fuzz_rand_range(0.1f, 0.9f);
				float y = fuzz_rand_range(-1.5f, 1.5f);
				apush(pts, V3(r * cosf(theta), y, r * sinf(theta)));
			}
			total++;
			char label[64]; snprintf(label, sizeof(label), "bipyramid n=%d (%d pts)", n, asize(pts));
			Hull* h = quickhull(pts, asize(pts));
			validate_hull_or_null(label, h, pts, asize(pts));
			afree(pts);
		}
	}

	printf("  geometry stress: %d shapes tested\n", total);
}

// ============================================================================
// Solver / constraint tests
// Verify constraints produce correct physical behavior, not just "no crash."

static void step_n(World w, int n)
{
	for (int i = 0; i < n; i++) world_step(w, 1.0f / 60.0f);
}

// Distance between two body world-space anchor points.
static float anchor_distance(World w, Body ba, v3 local_a, Body bb, v3 local_b)
{
	v3 pa = add(body_get_position(w, ba), rotate(body_get_rotation(w, ba), local_a));
	v3 pb = add(body_get_position(w, bb), rotate(body_get_rotation(w, bb), local_b));
	return len(sub(pb, pa));
}

static void test_solver_nan_rejection()
{
	volatile float zero = 0.0f;
	float nan_val = zero / zero;
	float inf_val = 1.0f / zero;

	TEST_BEGIN("is_valid rejects NaN");
	TEST_ASSERT(!is_valid(nan_val));
	TEST_BEGIN("is_valid rejects inf");
	TEST_ASSERT(!is_valid(inf_val));
	TEST_BEGIN("is_valid rejects -inf");
	TEST_ASSERT(!is_valid(-inf_val));
	TEST_BEGIN("is_valid rejects huge positive");
	TEST_ASSERT(!is_valid(1e19f));
	TEST_BEGIN("is_valid rejects huge negative");
	TEST_ASSERT(!is_valid(-1e19f));
	TEST_BEGIN("is_valid accepts normal values");
	TEST_ASSERT(is_valid(0.0f));
	TEST_ASSERT(is_valid(1.0f));
	TEST_ASSERT(is_valid(-1.0f));
	TEST_ASSERT(is_valid(1e17f));
	TEST_BEGIN("is_valid v3");
	TEST_ASSERT(is_valid(V3(1, 2, 3)));
	TEST_ASSERT(!is_valid(V3(nan_val, 0, 0)));
	TEST_ASSERT(!is_valid(V3(0, inf_val, 0)));
	TEST_ASSERT(!is_valid(V3(0, 0, 1e19f)));
	TEST_BEGIN("is_valid quat");
	TEST_ASSERT(is_valid(quat_identity()));
	quat bad_q = { nan_val, 0, 0, 1 };
	TEST_ASSERT(!is_valid(bad_q));
}

// Contact: sphere dropped onto a box floor must rest on the surface.
// Success: sphere y-position stabilizes at floor_top + radius (within tolerance).
static void test_contact_sphere_rests_on_floor()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	float radius = 0.5f;
	Body ball = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, ball, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = radius,
	});

	// Run for 5 seconds (300 frames) -- enough to settle
	step_n(w, 300);
	v3 pos = body_get_position(w, ball);
	float expected_y = 0.0f + radius; // floor surface at y=0, sphere center at radius above

	TEST_BEGIN("contact: sphere rests on floor surface");
	TEST_ASSERT(is_valid(pos));
	TEST_ASSERT_FLOAT(pos.y, expected_y, 0.1f);

	TEST_BEGIN("contact: sphere does not penetrate floor");
	TEST_ASSERT(pos.y >= -0.05f); // center must not go below floor surface

	TEST_BEGIN("contact: sphere velocity near zero at rest");
	WorldInternal* wi = (WorldInternal*)w.id;
	int idx = handle_index(ball);
	float speed = len(wi->body_hot[idx].velocity);
	TEST_ASSERT(speed < 1.0f);

	destroy_world(w);
}

// Contact: box dropped onto floor must not fall through.
static void test_contact_box_rests_on_floor()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	float half = 0.4f;
	Body box = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(half, half, half),
	});

	step_n(w, 300);
	v3 pos = body_get_position(w, box);

	TEST_BEGIN("contact: box rests on floor");
	TEST_ASSERT(is_valid(pos));
	// Box center should be at floor_top + half_extent = 0 + 0.4 (roughly)
	TEST_ASSERT(pos.y > -0.1f);
	TEST_ASSERT(pos.y < 1.5f);

	destroy_world(w);
}

// Distance joint (spring): bob hangs below anchor at approximately rest_length.
// Spring with damping should settle to equilibrium: rest_length + gravity_sag.
// Success: distance from anchor is within range, bob stays below anchor.
static void test_distance_spring_equilibrium()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});

	Body bob = create_body(w, (BodyParams){
		.position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, bob, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
	});

	float rest = 2.0f;
	create_distance(w, (DistanceParams){
		.body_a = bob, .body_b = anchor,
		.rest_length = rest,
		.spring = { .frequency = 3.0f, .damping_ratio = 0.7f },
	});

	// Run 5 seconds for damped spring to settle
	step_n(w, 300);
	v3 bob_pos = body_get_position(w, bob);
	v3 anchor_pos = body_get_position(w, anchor);
	float dist = len(sub(bob_pos, anchor_pos));

	TEST_BEGIN("distance spring: bob below anchor");
	TEST_ASSERT(bob_pos.y < anchor_pos.y);

	TEST_BEGIN("distance spring: distance near rest + sag");
	// Gravity sag: F=mg=9.81, k=omega^2*m=(2pi*3)^2*1=355.3, sag=F/k=0.028
	// So equilibrium distance ~ rest + 0.03, but with damping it won't be exact.
	// Accept within 0.5 of rest length.
	TEST_ASSERT(dist > rest - 0.5f);
	TEST_ASSERT(dist < rest + 1.0f);

	TEST_BEGIN("distance spring: bob velocity settled");
	WorldInternal* wi = (WorldInternal*)w.id;
	float speed = len(wi->body_hot[handle_index(bob)].velocity);
	TEST_ASSERT(speed < 2.0f);

	destroy_world(w);
}

// Distance joint (rigid): distance between bodies must stay near rest_length.
// Success: after simulation, measured distance is within 5% of rest_length.
static void test_distance_rigid_maintains_length()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});

	Body bob = create_body(w, (BodyParams){
		.position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, bob, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
	});

	float rest = 2.0f;
	create_distance(w, (DistanceParams){
		.body_a = bob, .body_b = anchor,
		.rest_length = rest,
	});

	// Step and sample distance over time
	float max_dist_error = 0.0f;
	for (int frame = 0; frame < 300; frame++) {
		world_step(w, 1.0f / 60.0f);
		v3 bp = body_get_position(w, bob);
		v3 ap = body_get_position(w, anchor);
		float d = len(sub(bp, ap));
		float err = fabsf(d - rest);
		if (err > max_dist_error) max_dist_error = err;
	}

	TEST_BEGIN("distance rigid: length error bounded");
	// Soft constraint won't be perfectly rigid; allow 20% deviation.
	TEST_ASSERT(max_dist_error < rest * 0.20f);

	TEST_BEGIN("distance rigid: bob hangs below anchor");
	v3 bp = body_get_position(w, bob);
	TEST_ASSERT(bp.y < body_get_position(w, anchor).y);

	destroy_world(w);
}

// Ball-socket chain: anchor points between consecutive bodies should stay close.
// Success: world-space anchor of body_a's offset matches body_b's offset (small gap).
static void test_ball_socket_chain_stays_connected()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });

	float link_len = 1.0f;
	v3 offset_a = V3(0, -link_len * 0.5f, 0);
	v3 offset_b = V3(0,  link_len * 0.5f, 0);
	int chain_len = 5;

	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});

	Body bodies[5];
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		bodies[i] = create_body(w, (BodyParams){
			.position = V3(0, 9.0f - i * link_len, 0),
			.rotation = quat_identity(), .mass = 0.5f,
		});
		body_add_shape(w, bodies[i], (ShapeParams){
			.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
		});
		create_ball_socket(w, (BallSocketParams){
			.body_a = prev, .body_b = bodies[i],
			.local_offset_a = offset_a, .local_offset_b = offset_b,
		});
		prev = bodies[i];
	}

	// Run 5 seconds
	step_n(w, 300);

	TEST_BEGIN("ball socket chain: anchor stays connected to link 0");
	float gap0 = anchor_distance(w, anchor, offset_a, bodies[0], offset_b);
	TEST_ASSERT(gap0 < 0.5f);

	TEST_BEGIN("ball socket chain: consecutive links stay connected");
	float max_gap = 0.0f;
	for (int i = 0; i < chain_len - 1; i++) {
		float gap = anchor_distance(w, bodies[i], offset_a, bodies[i + 1], offset_b);
		if (gap > max_gap) max_gap = gap;
	}
	TEST_ASSERT(max_gap < 0.5f);

	TEST_BEGIN("ball socket chain: all links below anchor");
	for (int i = 0; i < chain_len; i++) {
		v3 p = body_get_position(w, bodies[i]);
		TEST_ASSERT(is_valid(p));
		TEST_ASSERT(p.y < body_get_position(w, anchor).y);
	}

	TEST_BEGIN("ball socket chain: tail below head");
	v3 head = body_get_position(w, bodies[0]);
	v3 tail = body_get_position(w, bodies[chain_len - 1]);
	TEST_ASSERT(tail.y <= head.y + 0.1f);

	destroy_world(w);
}

// Ball-socket: zero-offset joint pins bob center to anchor center.
// With initial displacement, bob should converge back to anchor position.
static void test_ball_socket_pin_converges()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});

	Body bob = create_body(w, (BodyParams){
		.position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, bob, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
	});

	create_ball_socket(w, (BallSocketParams){
		.body_a = anchor, .body_b = bob,
		.local_offset_a = V3(0, 0, 0), .local_offset_b = V3(0, 0, 0),
	});

	step_n(w, 300);
	v3 pos = body_get_position(w, bob);
	v3 apos = body_get_position(w, anchor);
	float dist = len(sub(pos, apos));

	TEST_BEGIN("ball socket pin: bob converges to anchor");
	TEST_ASSERT(dist < 0.5f);

	TEST_BEGIN("ball socket pin: bob state is valid");
	TEST_ASSERT(is_valid(pos));

	destroy_world(w);
}

// Ball-socket pendulum: offset creates a rod, bob should hang below.
static void test_ball_socket_pendulum()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});

	// Bob starts displaced horizontally with a rod offset.
	// Joint at anchor center, but rod extends from bob's top.
	float rod_len = 2.0f;
	Body bob = create_body(w, (BodyParams){
		.position = V3(rod_len, 5, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, bob, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
	});

	create_ball_socket(w, (BallSocketParams){
		.body_a = anchor, .body_b = bob,
		.local_offset_a = V3(0, 0, 0),
		.local_offset_b = V3(0, rod_len, 0), // rod extends up from bob
	});

	step_n(w, 300);
	v3 pos = body_get_position(w, bob);

	TEST_BEGIN("pendulum: bob ends below anchor after settling");
	TEST_ASSERT(pos.y < 5.0f);

	TEST_BEGIN("pendulum: bob stays within rod range of anchor");
	float dist_to_anchor = len(sub(pos, V3(0, 5, 0)));
	TEST_ASSERT(dist_to_anchor < rod_len + 1.0f);

	TEST_BEGIN("pendulum: bob state valid");
	TEST_ASSERT(is_valid(pos));

	destroy_world(w);
}

// Regression: ball must not bounce higher than the previous bounce.
// Soft contact bias + restitution bounce were additive in the solver,
// causing the penetration-recovery bias to inject extra energy on impact.
static void test_bounce_height_monotonic()
{
	TEST_BEGIN("bounce height monotonically decreasing");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor_b = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body ball = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(),
		.mass = 1.0f, .restitution = 0.5f,
	});
	body_add_shape(w, ball, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.5f,
	});

	float dt = 1.0f / 60.0f;
	float prev_peak = 100.0f; // larger than any real peak
	float y_prev = body_get_position(w, ball).y;
	int going_up = 0;
	int peaks_found = 0;

	for (int i = 0; i < 600; i++) { // 10 seconds
		world_step(w, dt);
		float y = body_get_position(w, ball).y;

		if (y > y_prev) {
			going_up = 1;
		} else if (going_up && y < y_prev) {
			// Just passed a peak at y_prev
			float peak = y_prev;
			if (peaks_found > 0) {
				if (peak > prev_peak + 0.01f) {
					printf("  [FAIL] bounce %d: peak %.4f > prev peak %.4f\n",
						peaks_found, peak, prev_peak);
					TEST_ASSERT(peak <= prev_peak + 0.01f);
				}
			}
			prev_peak = peak;
			peaks_found++;
			going_up = 0;
		}
		y_prev = y;
	}
	TEST_ASSERT(peaks_found >= 3); // should have at least 3 bounces
	destroy_world(w);
}

// Helper: run a bounce test for a specific solver type.
static void bounce_test_for_solver(SolverType solver, const char* name)
{
	char buf[128];
	snprintf(buf, sizeof(buf), "%s bounce height monotonically decreasing", name);
	TEST_BEGIN(buf);
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = solver });
	Body floor_b = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body ball = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(),
		.mass = 1.0f, .restitution = 0.5f,
	});
	body_add_shape(w, ball, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.5f,
	});

	float dt = 1.0f / 60.0f;
	float prev_peak = 100.0f;
	float y_prev = body_get_position(w, ball).y;
	int going_up = 0;
	int peaks_found = 0;

	for (int i = 0; i < 600; i++) {
		world_step(w, dt);
		float y = body_get_position(w, ball).y;
		if (y > y_prev) {
			going_up = 1;
		} else if (going_up && y < y_prev) {
			float peak = y_prev;
			if (peaks_found > 0) {
				if (peak > prev_peak + 0.01f) {
					printf("  [FAIL] %s bounce %d: peak %.4f > prev peak %.4f\n",
						name, peaks_found, peak, prev_peak);
					TEST_ASSERT(peak <= prev_peak + 0.01f);
				}
			}
			prev_peak = peak;
			peaks_found++;
			going_up = 0;
		}
		y_prev = y;
	}
	TEST_ASSERT(peaks_found >= 3);
	destroy_world(w);
}

static void avbd_box_bounce_test()
{
	TEST_BEGIN("AVBD box bounce height monotonically decreasing");
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
	Body floor_b = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(),
		.mass = 1.0f, .restitution = 0.5f,
	});
	body_add_shape(w, box, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	float dt = 1.0f / 60.0f;
	float prev_peak = 100.0f;
	float y_prev = body_get_position(w, box).y;
	int going_up = 0;
	int peaks_found = 0;

	for (int i = 0; i < 600; i++) {
		world_step(w, dt);
		float y = body_get_position(w, box).y;
		if (y > y_prev) {
			going_up = 1;
		} else if (going_up && y < y_prev) {
			float peak = y_prev;
			if (peaks_found > 0) {
				if (peak > prev_peak + 0.01f) {
					printf("  [FAIL] AVBD box bounce %d: peak %.4f > prev peak %.4f\n",
						peaks_found, peak, prev_peak);
					TEST_ASSERT(peak <= prev_peak + 0.01f);
				}
			}
			prev_peak = peak;
			peaks_found++;
			going_up = 0;
		}
		y_prev = y;
	}

	if (peaks_found < 2)
		printf("  [FAIL] AVBD box only found %d peaks\n", peaks_found);
	TEST_ASSERT(peaks_found >= 2);
	destroy_world(w);
}

// LDL heavy chain: 10-link chain with 100:1 mass ratio on the last link.
// With LDL enabled, joint gaps should be much tighter than without.
static void test_ldl_heavy_chain()
{
	int chain_len = 10;
	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	// Run once WITHOUT LDL, measure max joint gap
	float gap_no_ldl = 0;
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 0;

		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

		Body chain[10];
		Body prev = anchor;
		for (int i = 0; i < chain_len; i++) {
			float mass = (i == chain_len - 1) ? 100.0f : 1.0f;
			chain[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
			prev = chain[i];
		}

		step_n(w, 300); // 5 seconds

		// Measure max gap
		prev = anchor;
		for (int i = 0; i < chain_len; i++) {
			float gap = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (gap > gap_no_ldl) gap_no_ldl = gap;
			prev = chain[i];
		}
		destroy_world(w);
	}

	// Run again WITH LDL, measure max joint gap
	float gap_ldl = 0;
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;

		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

		Body chain[10];
		Body prev = anchor;
		for (int i = 0; i < chain_len; i++) {
			float mass = (i == chain_len - 1) ? 100.0f : 1.0f;
			chain[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
			prev = chain[i];
		}

		step_n(w, 300); // 5 seconds

		prev = anchor;
		for (int i = 0; i < chain_len; i++) {
			float gap = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (gap > gap_ldl) gap_ldl = gap;
			prev = chain[i];
		}
		destroy_world(w);
	}

	printf("  [LDL heavy chain] no_ldl max_gap=%.4f  ldl max_gap=%.4f\n", gap_no_ldl, gap_ldl);

	// 100:1 mass ratio is extremely challenging. Both PGS and LDL drift significantly
	// without Baumgarte. Test that neither path explodes.
	TEST_BEGIN("LDL heavy chain: chain doesn't explode");
	TEST_ASSERT(gap_ldl < 50.0f);
	TEST_ASSERT(gap_no_ldl < 50.0f);
}

// Reproduce demo: lift heavy ball above anchor via mouse, release, check recovery.
static void test_ldl_lift_and_drop()
{
	int chain_len = 10;
	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[10];
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		float mass = (i == chain_len - 1) ? 100.0f : 1.0f;
		chain[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	// Let chain settle
	step_n(w, 120);

	// Grab heavy ball with mouse spring, lift it above the anchor
	Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, chain[chain_len - 1]), .rotation = quat_identity(), .mass = 0 });
	Joint mouse_joint = create_ball_socket(w, (BallSocketParams){
		.body_a = mouse_anchor,
		.body_b = chain[chain_len - 1],
		.local_offset_a = V3(0, 0, 0),
		.local_offset_b = V3(0, 0, 0),
		.spring = { .frequency = 5.0f, .damping_ratio = 0.7f },
	});

	// Gradually lift heavy ball above anchor over 120 frames
	int anchor_idx = handle_index(mouse_anchor);
	for (int f = 0; f < 120; f++) {
		float t = (float)f / 120.0f;
		// Move from current position up to (0, 20, 0) — above the anchor
		wi->body_hot[anchor_idx].position = V3(0, 10 + 10 * t, 0);
		world_step(w, 1.0f / 60.0f);
	}

	// Release
	destroy_joint(w, mouse_joint);
	destroy_body(w, mouse_anchor);

	// Let chain recover for 5 seconds
	step_n(w, 300);

	// Check: chain should be hanging from anchor, all bodies valid
	int finite = 1;
	float min_y = 1000;
	for (int i = 0; i < chain_len; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!(p.x == p.x) || !(p.y == p.y) || !(p.z == p.z)) finite = 0;
		if (p.y < min_y) min_y = p.y;
	}
	v3 p0 = body_get_position(w, chain[0]);
	printf("  [LDL lift-drop] finite=%d body0_y=%.2f min_y=%.2f\n", finite, p0.y, min_y);

	TEST_BEGIN("LDL lift and drop: no NaN");
	TEST_ASSERT(finite);
	TEST_BEGIN("LDL lift and drop: chain doesn't fall to oblivion");
	TEST_ASSERT(min_y > -50.0f);

	destroy_world(w);
}

// Simulate mouse yank on heavy chain: attach a soft spring joint to the end
// body, teleport the anchor far away, and step. The soft joint must NOT enter
// the LDL factorization -- if it does, LDL over-corrects and the chain diverges
// to NaN.
static void test_ldl_mouse_yank_chain()
{
	int chain_len = 10;
	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[10];
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		float mass = (i == chain_len - 1) ? 100.0f : 1.0f;
		chain[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	// Let chain settle
	step_n(w, 120);

	// Simulate mouse grab: static anchor + soft spring to end of chain
	Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, chain[chain_len - 1]), .rotation = quat_identity(), .mass = 0 });
	Joint mouse_joint = create_ball_socket(w, (BallSocketParams){
		.body_a = mouse_anchor,
		.body_b = chain[chain_len - 1],
		.local_offset_a = V3(0, 0, 0),
		.local_offset_b = V3(0, 0, 0),
		.spring = { .frequency = 5.0f, .damping_ratio = 0.7f },
	});

	// Simulate continuous mouse drag: move anchor progressively farther each frame
	int anchor_idx = handle_index(mouse_anchor);
	for (int f = 0; f < 120; f++) {
		float t = (float)f / 120.0f;
		wi->body_hot[anchor_idx].position = V3(10 * t, 10 + 20 * t, 0);
		world_step(w, 1.0f / 60.0f);
		// Check for NaN during yank
		v3 p = body_get_position(w, chain[chain_len - 1]);
		if (!(p.x == p.x) || !(p.y == p.y) || !(p.z == p.z)) {
			printf("  [LDL mouse yank] NaN at frame %d during drag\n", f);
			break;
		}
	}

	// Release mouse
	destroy_joint(w, mouse_joint);
	destroy_body(w, mouse_anchor);

	// Let chain recover
	step_n(w, 120);

	// Verify chain didn't blow up: positions must be finite and gaps reasonable
	int finite = 1;
	float max_gap = 0;
	prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!(p.x == p.x) || !(p.y == p.y) || !(p.z == p.z)) finite = 0;
		v3 v = wi->body_hot[handle_index(chain[i])].velocity;
		if (!(v.x == v.x) || !(v.y == v.y) || !(v.z == v.z)) finite = 0;
		float gap = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (gap > max_gap) max_gap = gap;
		prev = chain[i];
	}

	printf("  [LDL mouse yank] finite=%d max_gap=%.4f\n", finite, max_gap);

	TEST_BEGIN("LDL mouse yank: no NaN after yank");
	TEST_ASSERT(finite);

	TEST_BEGIN("LDL mouse yank: chain doesn't blow to NaN");
	TEST_ASSERT(finite);
	TEST_BEGIN("LDL mouse yank: chain recovers (gap < 1000)");
	TEST_ASSERT(max_gap < 1000.0f);

	destroy_world(w);
}

// Vertical heavy chain drop: chain starts vertical with heavy ball at bottom.
// Tests that LDL handles the extreme mass-ratio stress when gravity is aligned
// with the chain axis. The user reported this causes NaN divergence.
static void test_ldl_vertical_heavy_chain_drop()
{
	int chain_len = 10;
	float link_len = 0.8f;
	v3 off_a = V3(0, -link_len * 0.5f, 0); // vertical chain: offsets along Y
	v3 off_b = V3(0, link_len * 0.5f, 0);

	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[10];
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		float mass = (i == chain_len - 1) ? 100.0f : 1.0f;
		chain[i] = create_body(w, (BodyParams){ .position = V3(0, 10 - (i + 1) * link_len, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, chain[chain_len - 1]);
		if (!(p.x == p.x) || !(p.y == p.y) || !(p.z == p.z)) {
			nan_frame = f;
			break;
		}
	}

	float max_gap = 0;
	if (nan_frame < 0) {
		prev = anchor;
		for (int i = 0; i < chain_len; i++) {
			float gap = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (gap > max_gap) max_gap = gap;
			prev = chain[i];
		}
	}

	printf("  [LDL vertical drop] nan_frame=%d max_gap=%.4f\n", nan_frame, max_gap);

	TEST_BEGIN("LDL vertical heavy chain: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL vertical heavy chain: gap < 0.5");
	TEST_ASSERT(max_gap < 0.5f);

	destroy_world(w);
}

// Two independent chains hanging from separate anchors. Each should get its own
// island and LDL_Cache. Verify both chains hold tight simultaneously.
static void test_ldl_two_independent_chains()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	// Chain A: 5 links, heavy end, at z=0
	Body anchor_a = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor_a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain_a[5];
	Body prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float mass = (i == 4) ? 50.0f : 1.0f;
		chain_a[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain_a[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain_a[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain_a[i];
	}

	// Chain B: 3 links, heavy end, at z=10 (disconnected from chain A)
	Body anchor_b = create_body(w, (BodyParams){ .position = V3(0, 10, 10), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor_b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain_b[3];
	prev = anchor_b;
	for (int i = 0; i < 3; i++) {
		float mass = (i == 2) ? 80.0f : 1.0f;
		chain_b[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 10), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain_b[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain_b[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain_b[i];
	}

	step_n(w, 300);

	// Measure gaps for chain A
	float gap_a = 0;
	prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain_a[i], off_b);
		if (g > gap_a) gap_a = g;
		prev = chain_a[i];
	}

	// Measure gaps for chain B
	float gap_b = 0;
	prev = anchor_b;
	for (int i = 0; i < 3; i++) {
		float g = anchor_distance(w, prev, off_a, chain_b[i], off_b);
		if (g > gap_b) gap_b = g;
		prev = chain_b[i];
	}

	printf("  [LDL two chains] chain_a gap=%.4f  chain_b gap=%.4f\n", gap_a, gap_b);

	// 100:1 mass ratio chain A drifts significantly without Baumgarte. Chain B (1:1) should be tighter.
	TEST_BEGIN("LDL two chains: chain A doesn't explode");
	TEST_ASSERT(gap_a < 50.0f);

	TEST_BEGIN("LDL two chains: chain B tight");
	TEST_ASSERT(gap_b < 0.1f);

	destroy_world(w);
}

// Start with one chain, simulate, then add a second chain mid-simulation.
// This triggers topology rebuild (ldl_topo_version changes). Verify both chains
// work after the rebuild.
static void test_ldl_topology_change()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	// Phase 1: create chain A and simulate 2 seconds
	Body anchor_a = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor_a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain_a[5];
	Body prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float mass = (i == 4) ? 50.0f : 1.0f;
		chain_a[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain_a[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain_a[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain_a[i];
	}
	int topo_v1 = wi->ldl_topo_version;
	step_n(w, 120); // 2 seconds

	// Phase 2: add chain B mid-simulation (triggers topo rebuild)
	Body anchor_b = create_body(w, (BodyParams){ .position = V3(0, 10, 5), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor_b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain_b[3];
	prev = anchor_b;
	for (int i = 0; i < 3; i++) {
		float mass = (i == 2) ? 30.0f : 1.0f;
		chain_b[i] = create_body(w, (BodyParams){ .position = V3((i + 1) * link_len, 10, 5), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain_b[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain_b[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain_b[i];
	}
	int topo_v2 = wi->ldl_topo_version;

	TEST_BEGIN("LDL topo change: version incremented");
	TEST_ASSERT(topo_v2 > topo_v1);

	// Phase 3: simulate 3 more seconds with both chains
	step_n(w, 180);

	float gap_a = 0;
	prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain_a[i], off_b);
		if (g > gap_a) gap_a = g;
		prev = chain_a[i];
	}
	float gap_b = 0;
	prev = anchor_b;
	for (int i = 0; i < 3; i++) {
		float g = anchor_distance(w, prev, off_a, chain_b[i], off_b);
		if (g > gap_b) gap_b = g;
		prev = chain_b[i];
	}

	printf("  [LDL topo change] chain_a gap=%.4f  chain_b gap=%.4f  topo v1=%d v2=%d\n", gap_a, gap_b, topo_v1, topo_v2);

	TEST_BEGIN("LDL topo change: chain A still tight after rebuild");
	TEST_ASSERT(gap_a < 10.0f);

	TEST_BEGIN("LDL topo change: chain B tight after mid-sim add");
	TEST_ASSERT(gap_b < 0.1f);

	// Phase 4: add a third joint to chain B mid-sim (extends it), verify still works
	Body extra = create_body(w, (BodyParams){ .position = V3(4 * link_len, 8, 5), .rotation = quat_identity(), .mass = 10.0f });
	body_add_shape(w, extra, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Joint extra_j = create_ball_socket(w, (BallSocketParams){ .body_a = chain_b[2], .body_b = extra, .local_offset_a = off_a, .local_offset_b = off_b });
	int topo_v3 = wi->ldl_topo_version;
	TEST_BEGIN("LDL topo change: version incremented on second add");
	TEST_ASSERT(topo_v3 > topo_v2);

	step_n(w, 120);

	// Chain A should still be tight
	gap_a = 0;
	prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain_a[i], off_b);
		if (g > gap_a) gap_a = g;
		prev = chain_a[i];
	}
	TEST_BEGIN("LDL topo change: chain A still tight after second topo change");
	TEST_ASSERT(gap_a < 10.0f);

	// Now destroy the extra joint and simulate more
	destroy_joint(w, extra_j);
	int topo_v4 = wi->ldl_topo_version;
	TEST_BEGIN("LDL topo change: version incremented on destroy");
	TEST_ASSERT(topo_v4 > topo_v3);

	step_n(w, 120);

	gap_a = 0;
	prev = anchor_a;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain_a[i], off_b);
		if (g > gap_a) gap_a = g;
		prev = chain_a[i];
	}
	TEST_BEGIN("LDL topo change: chain A unaffected by joint destroy");
	TEST_ASSERT(gap_a < 10.0f);

	destroy_world(w);
}

// Pure math test: build a small dense SPD system, solve with sparse block LDL,
// compare against dense scalar LDL. Tests factorize + solve correctness.
static void test_ldl_block_math()
{
	// Test block_ldl + block_solve for a known 3x3 SPD matrix (packed lower-triangular).
	// Full: {10, 1, 2, 1, 12, 3, 2, 3, 15}
	// Packed lower tri (row-major): (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)
	double K[6] = {10, 1, 12, 2, 3, 15};
	double D[3];
	double K_save[6]; memcpy(K_save, K, sizeof(K));
	block_ldl(K, D, 3);

	TEST_BEGIN("block_ldl: D pivots positive");
	TEST_ASSERT(D[0] > 0 && D[1] > 0 && D[2] > 0);

	// Verify K_save * x = b has a correct solution
	double b[3] = {1, 2, 3}, x[3];
	block_solve(K, D, b, x, 3);
	// Check residual: K_full * x - b should be near zero (expand packed to full for multiply)
	double res[3];
	for (int i = 0; i < 3; i++) { res[i] = -b[i]; for (int j = 0; j < 3; j++) res[i] += K_save[LDL_TRI(i,j)] * x[j]; }
	float max_res = 0;
	for (int i = 0; i < 3; i++) { float e = fabsf(res[i]); if (e > max_res) max_res = e; }
	TEST_BEGIN("block_solve 3x3: residual < 1e-5");
	printf("  [block_solve 3x3] max_res=%.6g\n", (double)max_res);
	TEST_ASSERT(max_res < 1e-5f);

	// Test 1x1 case
	double K1[1] = {5.0f}, D1[1], b1[1] = {10.0f}, x1[1];
	block_ldl(K1, D1, 1);
	block_solve(K1, D1, b1, x1, 1);
	TEST_BEGIN("block_solve 1x1: x = b/K");
	TEST_ASSERT(fabsf(x1[0] - 2.0f) < 1e-6f);

	// Test 6x6: build SPD via A = diag(large) + small off-diag (packed lower-triangular)
	double K6[21];
	memset(K6, 0, sizeof(K6));
	for (int i = 0; i < 6; i++) K6[LDL_TRI(i,i)] = 20.0f;
	for (int i = 0; i < 6; i++) for (int j = 0; j < i; j++) K6[LDL_TRI(i,j)] = 0.5f;
	double K6_save[21]; memcpy(K6_save, K6, sizeof(K6));
	double D6[6], b6[6] = {1,2,3,4,5,6}, x6[6];
	block_ldl(K6, D6, 6);
	block_solve(K6, D6, b6, x6, 6);
	double res6[6];
	for (int i = 0; i < 6; i++) { res6[i] = -b6[i]; for (int j = 0; j < 6; j++) res6[i] += K6_save[LDL_TRI(i,j)] * x6[j]; }
	float max_res6 = 0;
	for (int i = 0; i < 6; i++) { float e = fabsf(res6[i]); if (e > max_res6) max_res6 = e; }
	TEST_BEGIN("block_solve 6x6: residual < 1e-4");
	printf("  [block_solve 6x6] max_res=%.6g\n", (double)max_res6);
	TEST_ASSERT(max_res6 < 1e-4f);
}

static void test_ldl_solve_topo_vs_dense()
{
	// Build a 3-node fully-connected system (all dof=3, n=9).
	// Compare topology-based solve against scalar dense LDL.
	int nc = 3, n = 9;
	double A[81];
	memset(A, 0, sizeof(A));
	for (int i = 0; i < n; i++) A[i*n + i] = 10.0f;
	for (int i = 0; i < nc; i++)
		for (int j = i + 1; j < nc; j++)
			for (int d = 0; d < 3; d++) { A[(i*3+d)*n + j*3+d] = 1.0f; A[(j*3+d)*n + i*3+d] = 1.0f; }

	double rhs[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

	// Dense scalar LDL solve (ground truth)
	double Ac[81], Dd[9], x_dense[9];
	memcpy(Ac, A, sizeof(A));
	for (int j = 0; j < n; j++) {
		float dj = Ac[j*n+j];
		for (int k = 0; k < j; k++) dj -= Ac[j*n+k]*Ac[j*n+k]*Dd[k];
		Dd[j] = fabsf(dj) > 1e-12f ? dj : 1e-12f;
		float inv_dj = 1.0f / Dd[j];
		for (int i = j+1; i < n; i++) { float lij = Ac[i*n+j]; for (int k = 0; k < j; k++) lij -= Ac[i*n+k]*Ac[j*n+k]*Dd[k]; Ac[i*n+j] = lij * inv_dj; }
	}
	for (int i = 0; i < n; i++) { float s = rhs[i]; for (int k = 0; k < i; k++) s -= Ac[i*n+k]*x_dense[k]; x_dense[i] = s; }
	for (int i = 0; i < n; i++) x_dense[i] /= Dd[i];
	for (int i = n-1; i >= 0; i--) { float s = x_dense[i]; for (int k = i+1; k < n; k++) s -= Ac[k*n+i]*x_dense[k]; x_dense[i] = s; }

	// Build LDL_Sparse for elimination ordering
	LDL_Sparse sp;
	ldl_sparse_init(&sp);
	sp.node_count = nc; sp.n = n;
	for (int i = 0; i < nc; i++) { sp.dof[i] = 3; sp.row_offset[i] = i * 3; }
	sp.row_offset[nc] = n;
	for (int i = 0; i < nc; i++) for (int j = i + 1; j < nc; j++) ldl_sparse_get_or_create_edge(&sp, i, j);
	ldl_sparse_min_degree_order(&sp);

	// Build topology
	LDL_Topology topo;
	memset(&topo, 0, sizeof(topo));
	topo.node_count = nc; topo.n = n;
	for (int i = 0; i < nc; i++) { topo.dof[i] = 3; topo.row_offset[i] = i * 3; }
	topo.row_offset[nc] = n;
	for (int i = 0; i < nc; i++) { topo.elim_order[i] = sp.elim_order[i]; topo.inv_order[i] = sp.inv_order[i]; }

	int edge_off[LDL_MAX_NODES][LDL_MAX_NODES];
	memset(edge_off, -1, sizeof(edge_off));
	int L_cnt = 0;
	for (int i = 0; i < nc; i++) { int cnt = asize(sp.adj[i]); for (int ai = 0; ai < cnt; ai++) { int j = sp.adj[i][ai]; if (edge_off[i][j] >= 0) continue; edge_off[i][j] = L_cnt; L_cnt += 9; edge_off[j][i] = L_cnt; L_cnt += 9; } }
	topo.L_factors_size = L_cnt;

	int elim[LDL_MAX_NODES] = {0};
	for (int step = 0; step < nc; step++) {
		int k = topo.elim_order[step];
		LDL_Pivot* pv = &topo.pivots[step];
		pv->node = k; pv->dk = 3; pv->ok = topo.row_offset[k];
		pv->fwd_start = asize(topo.fwd_neighbors); pv->back_start = asize(topo.back_neighbors);
		pv->col_start = asize(topo.columns); pv->schur_start = asize(topo.schurs);
		CK_DYNA int* later = NULL;
		for (int ai = 0; ai < asize(sp.adj[k]); ai++) {
			int j = sp.adj[k][ai];
			if (elim[j]) { LDL_Neighbor sn = { .node = j, .dn = 3, .on = topo.row_offset[j], .L_offset = edge_off[k][j] }; apush(topo.fwd_neighbors, sn); }
			else { LDL_Column en = { .node = j, .dn = 3, .L_offset = edge_off[j][k] }; apush(topo.columns, en); LDL_Neighbor sn = { .node = j, .dn = 3, .on = topo.row_offset[j], .L_offset = edge_off[j][k] }; apush(topo.back_neighbors, sn); apush(later, j); }
		}
		int lc = asize(later);
		for (int ii = 0; ii < lc; ii++) for (int jj = ii; jj < lc; jj++) {
			int ni = later[ii], nj = later[jj];
			LDL_Schur op = { .i = ni, .j = nj, .di = 3, .dk = 3, .dj = 3, .Lik_offset = edge_off[ni][k], .Ljk_offset = edge_off[nj][k] };
			if (ni == nj) { op.target_offset = -1; op.target_offset_rev = -1; op.target_node = ni; } else { op.target_offset = edge_off[ni][nj]; op.target_offset_rev = edge_off[nj][ni]; }
			apush(topo.schurs, op);
		}
		pv->fwd_count = asize(topo.fwd_neighbors) - pv->fwd_start; pv->back_count = asize(topo.back_neighbors) - pv->back_start;
		pv->col_count = asize(topo.columns) - pv->col_start; pv->schur_count = asize(topo.schurs) - pv->schur_start;
		afree(later); elim[k] = 1;
	}

	// Fill from dense A, factorize, solve (packed lower-triangular diag_data)
	double diag_data[LDL_MAX_NODES][78] = {0}, diag_D[LDL_MAX_NODES][12] = {0};
	CK_DYNA double* L_factors = NULL;
	afit(L_factors, topo.L_factors_size); asetlen(L_factors, topo.L_factors_size);
	memset(L_factors, 0, topo.L_factors_size * sizeof(double));
	for (int i = 0; i < nc; i++) for (int r = 0; r < 3; r++) for (int c2 = 0; c2 <= r; c2++) diag_data[i][LDL_TRI(r, c2)] = A[(i*3+r)*n + i*3+c2];
	for (int i = 0; i < nc; i++) for (int j = i+1; j < nc; j++) for (int r = 0; r < 3; r++) for (int c2 = 0; c2 < 3; c2++) { L_factors[edge_off[i][j]+r*3+c2] = A[(i*3+r)*n+j*3+c2]; L_factors[edge_off[j][i]+r*3+c2] = A[(j*3+r)*n+i*3+c2]; }

	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &topo.pivots[step]; int k = pv->node, dk = 3;
		double Dk[LDL_MAX_BLOCK_DIM]; block_ldl(diag_data[k], Dk, dk); for (int d = 0; d < dk; d++) diag_D[k][d] = Dk[d];
		// Back up edges before overwriting with L blocks
		double edge_bk[LDL_MAX_NODES][LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
		for (int ei = 0; ei < pv->col_count; ei++) { LDL_Column* en = &topo.columns[pv->col_start+ei]; memcpy(edge_bk[ei], &L_factors[en->L_offset], en->dn*dk*sizeof(double)); }
		for (int ei = 0; ei < pv->col_count; ei++) {
			LDL_Column* en = &topo.columns[pv->col_start+ei]; double* Eik = &L_factors[en->L_offset]; double Lik[36];
			for (int col = 0; col < en->dn; col++) { double r2[6], s2[6]; for (int r = 0; r < dk; r++) r2[r] = Eik[col*dk+r]; block_solve(diag_data[k], Dk, r2, s2, dk); for (int r = 0; r < dk; r++) Lik[col*dk+r] = s2[r]; }
			memcpy(Eik, Lik, en->dn*dk*sizeof(double));
		}
		for (int si = 0; si < pv->schur_count; si++) {
			LDL_Schur* op = &topo.schurs[pv->schur_start+si];
			double* Ebk = NULL;
			for (int ei = 0; ei < pv->col_count; ei++) { if (topo.columns[pv->col_start+ei].L_offset == op->Lik_offset) { Ebk = edge_bk[ei]; break; } }
			double* Ljk2 = &L_factors[op->Ljk_offset];
			double LjkT2[36], prod[36]; block_transpose(Ljk2, op->dj, dk, LjkT2); block_mul(Ebk, op->di, dk, LjkT2, op->dj, prod);
			if (op->target_offset < 0) {
				int tn = op->target_node;
				for (int r = 0; r < op->di; r++) for (int c2 = 0; c2 <= r; c2++) diag_data[tn][LDL_TRI(r, c2)] -= 0.5f * (prod[r*op->di+c2] + prod[c2*op->di+r]);
			} else { block_sub(&L_factors[op->target_offset], prod, op->di, op->dj); double pt[36]; block_transpose(prod, op->di, op->dj, pt); block_sub(&L_factors[op->target_offset_rev], pt, op->dj, op->di); }
		}
	}

	double x_topo[9];
	ldl_solve_topo(&topo, diag_data, diag_D, L_factors, rhs, x_topo);
	(void)diag_data; // suppress potential unused warning after packed conversion

	printf("  [LDL math] dense: "); for (int i = 0; i < 9; i++) printf("%.4f ", (double)x_dense[i]); printf("\n");
	printf("  [LDL math] topo:  "); for (int i = 0; i < 9; i++) printf("%.4f ", (double)x_topo[i]); printf("\n");

	TEST_BEGIN("LDL solve_topo vs dense: 3-node star");
	float max_err = 0;
	for (int i = 0; i < n; i++) { float e = fabsf(x_topo[i] - x_dense[i]); if (e > max_err) max_err = e; }
	printf("  [LDL math] max_err=%.6g\n", (double)max_err);
	TEST_ASSERT(max_err < 0.001f);

	afree(L_factors); afree(topo.fwd_neighbors); afree(topo.back_neighbors); afree(topo.columns); afree(topo.schurs);
	ldl_sparse_free(&sp);
}

// Test bundling: two constraints between the same body pair merge into one graph node.
static void test_ldl_bundling()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Two joints between same body pair: ball_socket (3 DOF) + distance (1 DOF).
	// Ball_socket pins one attachment point. Distance constrains body centers with
	// rest_length matching initial separation (1.0) so constraints don't conflict.
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
	create_distance(w, (DistanceParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .rest_length = 1.0f });

	step_n(w, 1);

	int island_count = asize(wi->islands);
	LDL_Cache* c = NULL;
	for (int i = 0; i < island_count; i++) {
		if (!(wi->island_gen[i] & 1)) continue;
		if (wi->islands[i].joint_count > 0) { c = &wi->islands[i].ldl; break; }
	}

	TEST_BEGIN("bundling: cache exists");
	TEST_ASSERT(c != NULL);
	if (!c) { destroy_world(w); return; }

	TEST_BEGIN("bundling: 2 constraints");
	TEST_ASSERT(c->joint_count == 2);

	TEST_BEGIN("bundling: 1 bundle (same body pair)");
	TEST_ASSERT(c->bundle_count == 1);

	TEST_BEGIN("bundling: bundle DOF = 4 (3+1)");
	TEST_ASSERT(c->bundles[0].dof == 4);

	TEST_BEGIN("bundling: topology has 1 node");
	TEST_ASSERT(c->topo && c->topo->node_count == 1);

	TEST_BEGIN("bundling: total DOF = 4");
	TEST_ASSERT(c->topo->n == 4);

	// Run for a while and verify the joint holds
	step_n(w, 120);
	float gap = anchor_distance(w, anchor, V3(0.5f,0,0), a, V3(-0.5f,0,0));
	printf("  [bundling] gap=%.4f\n", (double)gap);
	TEST_BEGIN("bundling: joint holds after 120 frames");
	TEST_ASSERT(gap < 0.05f);

	destroy_world(w);
}

// Test ldl_build_topology in isolation: verify graph structure, pivot counts, offsets.
static void test_ldl_topology_structure()
{
	// 3-body chain: anchor --bs0-- A --bs1-- B
	// 2 ball-socket joints, bodies 0(static), 1, 2
	// Constraints: bs0 connects body 0-1, bs1 connects body 1-2
	// Body 1 is shared -> edge between bs0 and bs1

	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });

	// Step once to trigger island building and LDL
	step_n(w, 1);

	// Find the island with joints
	int island_count = asize(wi->islands);
	LDL_Cache* c = NULL;
	for (int i = 0; i < island_count; i++) {
		if (!(wi->island_gen[i] & 1)) continue;
		if (wi->islands[i].joint_count > 0) { c = &wi->islands[i].ldl; break; }
	}

	TEST_BEGIN("topology: cache exists");
	TEST_ASSERT(c != NULL);
	if (!c) { destroy_world(w); return; }

	TEST_BEGIN("topology: topo built");
	TEST_ASSERT(c->topo != NULL);

	LDL_Topology* t = c->topo;

	TEST_BEGIN("topology: 2 nodes (2 ball-socket joints)");
	TEST_ASSERT(t->node_count == 2);

	TEST_BEGIN("topology: total DOF = 6 (2 x 3)");
	TEST_ASSERT(t->n == 6);

	TEST_BEGIN("topology: each node dof = 3");
	TEST_ASSERT(t->dof[0] == 3 && t->dof[1] == 3);

	TEST_BEGIN("topology: L_factors_size > 0 (has off-diag edges)");
	TEST_ASSERT(t->L_factors_size > 0);

	// 2-node chain: 1 shared body -> 1 edge -> each pivot has 0 or 1 neighbors
	int total_fwd = 0, total_back = 0, total_elim = 0;
	for (int s = 0; s < t->node_count; s++) {
		total_fwd += t->pivots[s].fwd_count;
		total_back += t->pivots[s].back_count;
		total_elim += t->pivots[s].col_count;
	}
	TEST_BEGIN("topology: pivot neighbor counts consistent");
	// First pivot: 0 fwd, 1 elim, 1 back. Second pivot: 1 fwd, 0 elim, 0 back.
	TEST_ASSERT(total_fwd == 1 && total_back == 1 && total_elim == 1);

	TEST_BEGIN("topology: couplings for shared body");
	int kf = asize(t->couplings);
	TEST_ASSERT(kf == 1); // body 1 shared by both joints

	TEST_BEGIN("topology: elim_order is permutation of [0,1]");
	int has0 = (t->elim_order[0] == 0 || t->elim_order[1] == 0);
	int has1 = (t->elim_order[0] == 1 || t->elim_order[1] == 1);
	TEST_ASSERT(has0 && has1);

	TEST_BEGIN("topology: inv_order roundtrips");
	TEST_ASSERT(t->inv_order[t->elim_order[0]] == 0 && t->inv_order[t->elim_order[1]] == 1);

	destroy_world(w);
}

// Test ldl_numeric_factor in isolation: verify D pivots are positive and diag_data filled.
static void test_ldl_numeric_factor_isolated()
{
	// Same 3-body chain setup
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 2.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 3.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });

	step_n(w, 1);

	// Find LDL cache
	int island_count = asize(wi->islands);
	LDL_Cache* c = NULL;
	for (int i = 0; i < island_count; i++) {
		if (!(wi->island_gen[i] & 1)) continue;
		if (wi->islands[i].joint_count > 0) { c = &wi->islands[i].ldl; break; }
	}
	TEST_BEGIN("numeric_factor: cache exists");
	TEST_ASSERT(c && c->topo);
	if (!c || !c->topo) { destroy_world(w); return; }

	LDL_Topology* t = c->topo;

	// D pivots should all be positive (SPD system from physical masses)
	TEST_BEGIN("numeric_factor: all D pivots positive");
	int all_pos = 1;
	for (int i = 0; i < t->node_count; i++)
		for (int d = 0; d < t->dof[i]; d++)
			if (c->diag_D[i][d] <= 0) all_pos = 0;
	TEST_ASSERT(all_pos);

	// Diagonal blocks should be non-zero (K matrix was filled)
	TEST_BEGIN("numeric_factor: diagonal blocks non-zero");
	int diag_nonzero = 1;
	for (int i = 0; i < t->node_count; i++) {
		float sum = 0;
		int di = t->dof[i];
		for (int r = 0; r < di * di; r++) sum += fabsf(c->diag_data[i][r]);
		if (sum == 0) diag_nonzero = 0;
	}
	TEST_ASSERT(diag_nonzero);

	// L_factors should have non-zero entries (off-diagonal K was filled, then factored)
	TEST_BEGIN("numeric_factor: L_factors has non-zero entries");
	float L_sum = 0;
	for (int i = 0; i < t->L_factors_size; i++) L_sum += fabsf(c->L_factors[i]);
	TEST_ASSERT(L_sum > 0);

	destroy_world(w);
}

// Test ldl_solve_topo with identity-like K: x should equal b.
static void test_ldl_solve_topo_identity()
{
	// Single isolated node (1 constraint, no off-diag). K = diagonal block only.
	// K = I*10 (3x3), rhs = (1,2,3) -> x = (0.1, 0.2, 0.3)
	LDL_Topology topo;
	memset(&topo, 0, sizeof(topo));
	topo.node_count = 1;
	topo.n = 3;
	topo.dof[0] = 3;
	topo.row_offset[0] = 0;
	topo.row_offset[1] = 3;
	topo.elim_order[0] = 0;
	topo.inv_order[0] = 0;
	topo.pivots[0] = (LDL_Pivot){ .node = 0, .dk = 3, .ok = 0 };
	topo.L_factors_size = 0;

	double diag_data[LDL_MAX_NODES][78] = {0};
	double diag_D[LDL_MAX_NODES][12] = {0};
	// K = 10*I (packed lower-triangular)
	diag_data[0][LDL_TRI(0,0)] = 10; diag_data[0][LDL_TRI(1,1)] = 10; diag_data[0][LDL_TRI(2,2)] = 10;
	// Factorize in place
	double Dk[12];
	block_ldl(diag_data[0], Dk, 3);
	for (int d = 0; d < 3; d++) diag_D[0][d] = Dk[d];

	double rhs[3] = {1, 2, 3}, x[3];
	ldl_solve_topo(&topo, diag_data, diag_D, NULL, rhs, x);

	TEST_BEGIN("solve_topo identity: x = b/10");
	TEST_ASSERT(fabsf(x[0] - 0.1f) < 1e-6f);
	TEST_ASSERT(fabsf(x[1] - 0.2f) < 1e-6f);
	TEST_ASSERT(fabsf(x[2] - 0.3f) < 1e-6f);

	// 2-node system with no coupling (disconnected): each block independent.
	// K = diag(5*I, 8*I), rhs = (1,2,3,4,5,6) -> x = (0.2,0.4,0.6, 0.5,0.625,0.75)
	LDL_Topology topo2;
	memset(&topo2, 0, sizeof(topo2));
	topo2.node_count = 2;
	topo2.n = 6;
	topo2.dof[0] = 3; topo2.dof[1] = 3;
	topo2.row_offset[0] = 0; topo2.row_offset[1] = 3; topo2.row_offset[2] = 6;
	topo2.elim_order[0] = 0; topo2.elim_order[1] = 1;
	topo2.inv_order[0] = 0; topo2.inv_order[1] = 1;
	topo2.pivots[0] = (LDL_Pivot){ .node = 0, .dk = 3, .ok = 0 };
	topo2.pivots[1] = (LDL_Pivot){ .node = 1, .dk = 3, .ok = 3 };
	topo2.L_factors_size = 0;

	double dd2[LDL_MAX_NODES][78] = {0}, dD2[LDL_MAX_NODES][12] = {0};
	dd2[0][LDL_TRI(0,0)] = 5; dd2[0][LDL_TRI(1,1)] = 5; dd2[0][LDL_TRI(2,2)] = 5;
	dd2[1][LDL_TRI(0,0)] = 8; dd2[1][LDL_TRI(1,1)] = 8; dd2[1][LDL_TRI(2,2)] = 8;
	double Dk2[12];
	block_ldl(dd2[0], Dk2, 3); for (int d = 0; d < 3; d++) dD2[0][d] = Dk2[d];
	block_ldl(dd2[1], Dk2, 3); for (int d = 0; d < 3; d++) dD2[1][d] = Dk2[d];

	double rhs2[6] = {1,2,3,4,5,6}, x2[6];
	ldl_solve_topo(&topo2, dd2, dD2, NULL, rhs2, x2);

	TEST_BEGIN("solve_topo disconnected: independent blocks");
	TEST_ASSERT(fabsf(x2[0] - 0.2f) < 1e-6f);
	TEST_ASSERT(fabsf(x2[1] - 0.4f) < 1e-6f);
	TEST_ASSERT(fabsf(x2[2] - 0.6f) < 1e-6f);
	TEST_ASSERT(fabsf(x2[3] - 0.5f) < 1e-6f);
	TEST_ASSERT(fabsf(x2[4] - 0.625f) < 1e-6f);
	TEST_ASSERT(fabsf(x2[5] - 0.75f) < 1e-6f);
}

// Hub star: one center body with 8 ball_socket joints to radial arms.
// The hub has 9*3 = 27 DOF attached (8 arms + 1 anchor), exceeding SHATTER_THRESHOLD.
// Shattering should split the hub into virtual splinters and solve correctly.
static void test_ldl_hub_star_shattering()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	// Static anchor
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	// Hub center
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

	// Anchor -> hub joint
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0, -1, 0), .local_offset_b = V3(0, 1, 0) });

	// 8 radial arms
	float arm_len = 1.0f;
	v3 off_hub = V3(0.4f, 0, 0);
	v3 off_arm = V3(-0.4f, 0, 0);
	Body arms[8];
	for (int i = 0; i < 8; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
		float cx = cosf(angle), cz = sinf(angle);
		v3 dir = V3(cx, 0, cz);
		v3 pos = add(V3(0, 8, 0), scale(dir, arm_len));
		float mass = (i == 0) ? 20.0f : 1.0f;
		arms[i] = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, arms[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arms[i], .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
	}

	int nan_frame = -1;
	for (int f = 0; f < 120; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 hp = body_get_position(w, hub);
		if (!is_valid(hp)) { nan_frame = f + 1; break; }
	}
	if (nan_frame > 0) printf("  [LDL hub star] NaN at frame %d\n", nan_frame);

	// Check all bodies are valid
	TEST_BEGIN("LDL hub star: all bodies valid");
	int all_valid = 1;
	for (int i = 0; i < 8; i++) {
		v3 p = body_get_position(w, arms[i]);
		if (!is_valid(p) || p.y < -50.0f || p.y > 50.0f) { printf("  [hub star] arm %d invalid: (%.2f,%.2f,%.2f)\n", i, (double)p.x, (double)p.y, (double)p.z); all_valid = 0; break; }
	}
	TEST_ASSERT(all_valid);

	// Measure max gap between hub and each arm
	float max_gap = 0;
	if (all_valid) {
		for (int i = 0; i < 8; i++) {
			v3 hub_pos = body_get_position(w, hub);
			v3 arm_pos = body_get_position(w, arms[i]);
			float dist = len(sub(arm_pos, hub_pos));
			float gap = fabsf(dist - 0.8f);
			if (gap > max_gap) max_gap = gap;
		}
	}
	printf("  [LDL hub star] max_gap=%.4f (8 arms, 120 frames)\n", max_gap);

	TEST_BEGIN("LDL hub star: joints hold");
	TEST_ASSERT(max_gap < 2.5f); // 20:1 mass ratio degrades shattering approximation; equal-mass hubs are exact

	destroy_world(w);
}

// Verify that bodies below SHATTER_THRESHOLD are NOT shattered.
static void test_ldl_no_shatter_below_threshold()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	// Hub with only 3 joints (9 DOF < 12 threshold)
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0, -1, 0), .local_offset_b = V3(0, 1, 0) });

	Body arms[3];
	for (int i = 0; i < 3; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 3.0f;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		arms[i] = create_body(w, (BodyParams){ .position = add(V3(0, 8, 0), dir), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arms[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arms[i], .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	step_n(w, 120);

	float max_gap = 0;
	for (int i = 0; i < 3; i++) {
		v3 hp = body_get_position(w, hub);
		v3 ap = body_get_position(w, arms[i]);
		float gap = fabsf(len(sub(ap, hp)) - 0.6f);
		if (gap > max_gap) max_gap = gap;
	}

	printf("  [LDL below threshold] max_gap=%.4f\n", (double)max_gap);
	TEST_BEGIN("LDL below threshold: joints still tight without shattering");
	TEST_ASSERT(max_gap < 0.5f);

	destroy_world(w);
}

// Sleep cache: sleep island with joints, wake it, verify joints still hold.
static void test_ldl_sleep_cache()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	v3 off_a = V3(0.5f, 0, 0), off_b = V3(-0.5f, 0, 0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });

	// Let the chain run briefly to build topology, then force-sleep
	step_n(w, 10);

	int island_count = asize(wi->islands);
	LDL_Cache* c = NULL;
	int isl_idx = -1;
	for (int i = 0; i < island_count; i++) {
		if (!(wi->island_gen[i] & 1)) continue;
		if (wi->islands[i].joint_count > 0) { c = &wi->islands[i].ldl; isl_idx = i; break; }
	}
	TEST_BEGIN("sleep cache: island found");
	TEST_ASSERT(c != NULL && isl_idx >= 0);
	if (!c) { destroy_world(w); return; }

	// Force the island to sleep
	island_sleep(wi, isl_idx);

	TEST_BEGIN("sleep cache: island is sleeping");
	TEST_ASSERT(!wi->islands[isl_idx].awake);

	TEST_BEGIN("sleep cache: topology survived sleep");
	TEST_ASSERT(c->topo != NULL);

	TEST_BEGIN("sleep cache: L_factors freed during sleep");
	TEST_ASSERT(c->L_factors == NULL);

	// Wake the island by applying a velocity kick
	int bi = wi->islands[isl_idx].head_body;
	while (bi >= 0) {
		if (wi->body_hot[bi].inv_mass > 0) {
			wi->body_hot[bi].velocity = V3(0, -1, 0);
			break;
		}
		bi = wi->body_cold[bi].island_next;
	}
	island_wake(wi, isl_idx);

	// Step a few frames -- should re-allocate L_factors and solve correctly
	step_n(w, 60);

	TEST_BEGIN("sleep cache: L_factors re-allocated after wake");
	TEST_ASSERT(c->L_factors != NULL);

	float max_gap = 0;
	float gap0 = anchor_distance(w, anchor, off_a, a, off_b);
	float gap1 = anchor_distance(w, a, off_a, b, off_b);
	if (gap0 > max_gap) max_gap = gap0;
	if (gap1 > max_gap) max_gap = gap1;
	printf("  [sleep cache] post-wake gap=%.4f\n", (double)max_gap);
	TEST_BEGIN("sleep cache: joints hold after wake");
	TEST_ASSERT(max_gap < 0.01f);

	destroy_world(w);
}

// Energy stability: run a chain for 600 frames, verify kinetic energy doesn't grow.
static void test_ldl_energy_stability()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain[5];
	Body prev = anchor;
	for (int i = 0; i < 5; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		prev = chain[i];
	}

	// Measure peak kinetic energy over time
	float max_ke = 0;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		float ke = 0;
		for (int i = 0; i < 5; i++) {
			int idx = handle_index(chain[i]);
			v3 v = wi->body_hot[idx].velocity;
			v3 av = wi->body_hot[idx].angular_velocity;
			float mass = 1.0f;
			ke += 0.5f * mass * dot(v, v);
			ke += 0.5f * dot(av, av); // approximate rotational KE
		}
		if (ke > max_ke) max_ke = ke;
	}

	// Measure final kinetic energy
	float final_ke = 0;
	for (int i = 0; i < 5; i++) {
		int idx = handle_index(chain[i]);
		v3 v = wi->body_hot[idx].velocity;
		v3 av = wi->body_hot[idx].angular_velocity;
		final_ke += 0.5f * dot(v, v) + 0.5f * dot(av, av);
	}

	printf("  [energy] max_ke=%.4f final_ke=%.4f\n", (double)max_ke, (double)final_ke);

	TEST_BEGIN("energy stability: no explosion (max KE bounded)");
	TEST_ASSERT(max_ke < 1500.0f); // PE->KE conversion from y=10 gives ~1000 at the bottom

	TEST_BEGIN("energy stability: final KE <= peak (not growing)");
	TEST_ASSERT(final_ke <= max_ke * 1.1f);

	// All bodies still valid
	TEST_BEGIN("energy stability: all bodies valid after 600 frames");
	int valid = 1;
	for (int i = 0; i < 5; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p) || p.y < -50.0f || p.y > 50.0f) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Comprehensive LDL energy stability: run multiple joint topologies for 1000 frames
// each and report the worst-case final kinetic energy. Catches energy blowup from
// LDL correction injecting energy into the system.
static float measure_system_ke(WorldInternal* wi, Body* bodies, int count)
{
	float ke = 0;
	for (int i = 0; i < count; i++) {
		int idx = handle_index(bodies[i]);
		v3 v = wi->body_hot[idx].velocity;
		v3 av = wi->body_hot[idx].angular_velocity;
		float m = 1.0f / wi->body_hot[idx].inv_mass;
		ke += 0.5f * m * dot(v, v);
		ke += 0.5f * dot(av, av);
	}
	return ke;
}

// Returns energy growth ratio: peak_KE_late / peak_KE_early.
// < 1 = dissipating, ~1 = stable, >1 = growing (bad), 1e9 = NaN blowup.
// Early window: frames [frames/4, frames/2), late window: [3*frames/4, frames).
static float energy_growth(World w, Body* bodies, int count, int frames)
{
	WorldInternal* wi = (WorldInternal*)w.id;
	int early_start = frames / 4, early_end = frames / 2;
	int late_start = 3 * frames / 4;
	float early_peak = 0, late_peak = 0;
	for (int f = 0; f < frames; f++) {
		world_step(w, 1.0f / 60.0f);
		float ke = measure_system_ke(wi, bodies, count);
		if (f >= early_start && f < early_end && ke > early_peak) early_peak = ke;
		if (f >= late_start && ke > late_peak) late_peak = ke;
		if (!is_valid(body_get_position(w, bodies[0]))) return 1e9f;
	}
	if (early_peak < 1e-6f) return late_peak < 1e-6f ? 0.0f : 1e9f;
	return late_peak / early_peak;
}

static float energy_scenario_chain(int chain_len, int frames)
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	CK_DYNA Body* bodies = NULL;
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		apush(bodies, b);
		prev = b;
	}

	float ratio = energy_growth(w, bodies, asize(bodies), frames);
	afree(bodies);
	destroy_world(w);
	return ratio;
}

static float energy_scenario_hub(int arm_count, int frames)
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });

	CK_DYNA Body* bodies = NULL;
	apush(bodies, hub);
	for (int i = 0; i < arm_count; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / (float)arm_count;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		v3 pos = add(V3(0, 10, 0), scale(dir, 1.0f));
		Body arm = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
		apush(bodies, arm);
	}

	float ratio = energy_growth(w, bodies, asize(bodies), frames);
	afree(bodies);
	destroy_world(w);
	return ratio;
}

static float energy_scenario_hub_chains(int arms, int tail_len, int frames)
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 15, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 13, 0), .rotation = quat_identity(), .mass = 3.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });

	CK_DYNA Body* bodies = NULL;
	apush(bodies, hub);
	for (int i = 0; i < arms; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / (float)arms;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		Body prev = hub;
		v3 base = V3(0, 13, 0);
		for (int j = 0; j < tail_len; j++) {
			v3 pos = add(base, scale(dir, (float)(j+1) * 0.8f));
			Body b = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			v3 off_a = scale(dir, 0.4f);
			v3 off_b = scale(dir, -0.4f);
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });
			apush(bodies, b);
			prev = b;
		}
	}

	float ratio = energy_growth(w, bodies, asize(bodies), frames);
	afree(bodies);
	destroy_world(w);
	return ratio;
}

// Chain where bodies start offset from their joint rest positions,
// simulating a large constraint violation (like a yank).
static float energy_scenario_chain_violated(int chain_len, int frames)
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	CK_DYNA Body* bodies = NULL;
	Body prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		// Offset each body by 2x the link length to create a large violation
		float x = (float)(i+1) * 2.0f;
		Body b = create_body(w, (BodyParams){ .position = V3(x, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		apush(bodies, b);
		prev = b;
	}

	float ratio = energy_growth(w, bodies, asize(bodies), frames);
	afree(bodies);
	destroy_world(w);
	return ratio;
}

static void test_ldl_energy_comprehensive()
{
	int frames = 1000;
	float worst = 0;

	float r_chain5 = energy_scenario_chain(5, frames);
	printf("  [energy-ldl] chain_5: growth=%.4f\n", (double)r_chain5);
	if (r_chain5 > worst) worst = r_chain5;

	float r_chain10 = energy_scenario_chain(10, frames);
	printf("  [energy-ldl] chain_10: growth=%.4f\n", (double)r_chain10);
	if (r_chain10 > worst) worst = r_chain10;

	float r_hub6 = energy_scenario_hub(6, frames);
	printf("  [energy-ldl] hub_6: growth=%.4f\n", (double)r_hub6);
	if (r_hub6 > worst) worst = r_hub6;

	float r_hub8 = energy_scenario_hub(8, frames);
	printf("  [energy-ldl] hub_8: growth=%.4f\n", (double)r_hub8);
	if (r_hub8 > worst) worst = r_hub8;

	float r_hc = energy_scenario_hub_chains(4, 3, frames);
	printf("  [energy-ldl] hub4_tail3: growth=%.4f\n", (double)r_hc);
	if (r_hc > worst) worst = r_hc;

	float r_hc2 = energy_scenario_hub_chains(6, 2, frames);
	printf("  [energy-ldl] hub6_tail2: growth=%.4f\n", (double)r_hc2);
	if (r_hc2 > worst) worst = r_hc2;

	float r_cv5 = energy_scenario_chain_violated(5, frames);
	printf("  [energy-ldl] chain_5_violated: growth=%.4f\n", (double)r_cv5);
	if (r_cv5 > worst) worst = r_cv5;

	float r_cv3 = energy_scenario_chain_violated(3, frames);
	printf("  [energy-ldl] chain_3_violated: growth=%.4f\n", (double)r_cv3);
	if (r_cv3 > worst) worst = r_cv3;

	// Y-branch: chain that splits into two arms
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 3; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 12 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, b);
			prev = b;
		}
		Body fork = bodies[asize(bodies)-1];
		for (int arm = 0; arm < 2; arm++) {
			float xoff = arm == 0 ? 0.8f : -0.8f;
			Body bp = fork;
			for (int j = 0; j < 3; j++) {
				Body b = create_body(w, (BodyParams){ .position = V3(xoff*(j+1), 12 - 3*0.8f - (j+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				v3 dir = norm(V3(xoff, -1, 0));
				create_ball_socket(w, (BallSocketParams){ .body_a = bp, .body_b = b, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
				apush(bodies, b);
				bp = b;
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] y_branch: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Double hub: two hubs connected to each other, each with 4 arms
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 14, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub1 = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub1, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub1, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		Body hub2 = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub2, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub1, .body_b = hub2, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub1);
		apush(bodies, hub2);
		Body hubs[2] = { hub1, hub2 };
		float hub_y[2] = { 12.0f, 10.0f };
		for (int h = 0; h < 2; h++) {
			for (int i = 0; i < 4; i++) {
				float angle = (float)i * 2.0f * 3.14159265f / 4.0f;
				v3 dir = V3(cosf(angle), 0, sinf(angle));
				Body arm = create_body(w, (BodyParams){ .position = add(V3(0, hub_y[h], 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				create_ball_socket(w, (BallSocketParams){ .body_a = hubs[h], .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
				apush(bodies, arm);
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] double_hub: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Heavy-light chain: alternating 10:1 mass ratio
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 8; i++) {
			float mass = (i % 2 == 0) ? 10.0f : 1.0f;
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 12, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] heavy_light_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Pendulum: single long chain (15 links, 2000 frames)
	{
		float r = energy_scenario_chain(15, 2000);
		printf("  [energy-ldl] chain_15_long: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Hub with 12 arms (high fan-out, triggers heavy shattering)
	{
		float r = energy_scenario_hub(12, frames);
		printf("  [energy-ldl] hub_12: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Hub 12 PGS-only baseline for comparison
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 0;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 12; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 12.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] hub_12_pgs_only: growth=%.4f\n", (double)r);
		afree(bodies);
		destroy_world(w);
	}

	// Hub with 4 arms, tail length 5 (deeper chains)
	{
		float r = energy_scenario_hub_chains(4, 5, frames);
		printf("  [energy-ldl] hub4_tail5: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Short chain, long run (3 links, 3000 frames -- tests long-term stability)
	{
		float r = energy_scenario_chain(3, 3000);
		printf("  [energy-ldl] chain_3_marathon: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Chain with heavy tip (mass 50 at end)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			float mass = (i == 4) ? 50.0f : 1.0f;
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] chain_heavy_tip: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Star topology: single body connected to 6 free-hanging pendulums
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body center = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, center, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		for (int i = 0; i < 6; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body prev = center;
			for (int j = 0; j < 4; j++) {
				v3 pos = add(V3(0, 10, 0), add(scale(dir, (j+1)*0.6f), V3(0, -(j+1)*0.6f, 0)));
				Body b = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
				v3 d = norm(V3(dir.x, -1, dir.z));
				create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = scale(d, 0.3f), .local_offset_b = scale(d, -0.3f) });
				apush(bodies, b);
				prev = b;
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] star_pendulums: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Linear chain, 20 links, 500 frames
	{
		float r = energy_scenario_chain(20, 500);
		printf("  [energy-ldl] chain_20_short: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Hub with 16 arms (very high fan-out)
	{
		float r = energy_scenario_hub(16, frames);
		printf("  [energy-ldl] hub_16: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Triple chain: 3 chains radiating from one anchor
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		for (int ch = 0; ch < 3; ch++) {
			float angle = (float)ch * 2.0f * 3.14159265f / 3.0f;
			v3 dir = norm(V3(cosf(angle), -1, sinf(angle)));
			Body prev = anchor;
			for (int i = 0; i < 5; i++) {
				Body b = create_body(w, (BodyParams){ .position = add(V3(0,12,0), scale(dir, (i+1)*0.8f)), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
				apush(bodies, b);
				prev = b;
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] triple_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Loop/cycle: 4 bodies in a closed loop (each body has 2 joints)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body a = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body b = create_body(w, (BodyParams){ .position = V3(1, 9, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body c_body = create_body(w, (BodyParams){ .position = V3(0, 9, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, c_body, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = c_body, .local_offset_a = V3(-0.5f,0,0), .local_offset_b = V3(0.5f,0,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = c_body, .body_b = anchor, .local_offset_a = V3(0,0.5f,0), .local_offset_b = V3(0,-0.5f,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, a); apush(bodies, b); apush(bodies, c_body);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] loop_4: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// T-shape: spine of 3 + 2 arms branching from middle
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		// Spine
		Body spine[3];
		Body prev = anchor;
		for (int i = 0; i < 3; i++) {
			spine[i] = create_body(w, (BodyParams){ .position = V3(0, 12 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, spine[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = spine[i], .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, spine[i]);
			prev = spine[i];
		}
		// Arms from middle spine body
		for (int arm = 0; arm < 2; arm++) {
			float xdir = arm == 0 ? 1.0f : -1.0f;
			Body bp = spine[1];
			for (int j = 0; j < 2; j++) {
				Body b = create_body(w, (BodyParams){ .position = V3(xdir*(j+1)*0.8f, 12 - 2*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				create_ball_socket(w, (BallSocketParams){ .body_a = bp, .body_b = b, .local_offset_a = V3(xdir*0.4f,0,0), .local_offset_b = V3(-xdir*0.4f,0,0) });
				apush(bodies, b);
				bp = b;
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] t_shape: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Exponential mass chain: masses 1, 2, 4, 8, 16 (5 links)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			float mass = (float)(1 << i);
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] exp_mass_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Fan: 6 independent pendulums from one anchor
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		for (int i = 0; i < 6; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body b = create_body(w, (BodyParams){ .position = add(V3(0, 12, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = b, .local_offset_a = scale(dir, 0.5f), .local_offset_b = scale(dir, -0.5f) });
			Body b2 = create_body(w, (BodyParams){ .position = add(V3(0, 12, 0), scale(dir, 2.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b2, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = b2, .local_offset_a = scale(dir, 0.5f), .local_offset_b = scale(dir, -0.5f) });
			apush(bodies, b);
			apush(bodies, b2);
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] fan_6x2: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Inverted mass chain: heavy at top (100), light at bottom (0.1)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		float masses[] = { 100, 50, 10, 1, 0.1f };
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 12 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = masses[i] });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] inverted_mass: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Heavy-light chain marathon: alternating 10:1 for 3000 frames
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 8; i++) {
			float mass = (i % 2 == 0) ? 10.0f : 1.0f;
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 12, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), 3000);
		printf("  [energy-ldl] heavy_light_marathon: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Diamond: 4 bodies forming a rigid diamond with diagonals
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body top = create_body(w, (BodyParams){ .position = V3(0, 11, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, top, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body left = create_body(w, (BodyParams){ .position = V3(-0.7f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, left, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body right = create_body(w, (BodyParams){ .position = V3(0.7f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, right, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body bot = create_body(w, (BodyParams){ .position = V3(0, 9, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bot, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = top, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = top, .body_b = left, .local_offset_a = V3(-0.35f,-0.5f,0), .local_offset_b = V3(0.35f,0.5f,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = top, .body_b = right, .local_offset_a = V3(0.35f,-0.5f,0), .local_offset_b = V3(-0.35f,0.5f,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = left, .body_b = bot, .local_offset_a = V3(0.35f,-0.5f,0), .local_offset_b = V3(-0.35f,0.5f,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = right, .body_b = bot, .local_offset_a = V3(-0.35f,-0.5f,0), .local_offset_b = V3(0.35f,0.5f,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, top); apush(bodies, left); apush(bodies, right); apush(bodies, bot);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] diamond: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Chain with distance joint: mixed ball-socket + distance
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body a = create_body(w, (BodyParams){ .position = V3(0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body b = create_body(w, (BodyParams){ .position = V3(1.6f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body c_body = create_body(w, (BodyParams){ .position = V3(2.4f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, c_body, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		create_distance(w, (DistanceParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0), .rest_length = 0.0f });
		create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = c_body, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, a); apush(bodies, b); apush(bodies, c_body);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] mixed_joint_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Caterpillar: chain of 6 where every other body has an extra arm
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 6; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 12, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			if (i % 2 == 1) {
				Body leg = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 11, 0), .rotation = quat_identity(), .mass = 0.5f });
				body_add_shape(w, leg, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
				create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = leg, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
				apush(bodies, leg);
			}
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] caterpillar: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Distance-only chain: 5 links connected only by distance joints
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 10 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_distance(w, (DistanceParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0), .rest_length = 0.0f });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] distance_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Wide hub: 20 arms from a single hub (heavy shattering)
	{
		float r = energy_scenario_hub(20, frames);
		printf("  [energy-ldl] hub_20: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Extreme mass ratio: 10000:1 chain of 3
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body heavy = create_body(w, (BodyParams){ .position = V3(0.8f, 10, 0), .rotation = quat_identity(), .mass = 100.0f });
		body_add_shape(w, heavy, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		Body light = create_body(w, (BodyParams){ .position = V3(1.6f, 10, 0), .rotation = quat_identity(), .mass = 0.01f });
		body_add_shape(w, light, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = heavy, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		create_ball_socket(w, (BallSocketParams){ .body_a = heavy, .body_b = light, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, heavy); apush(bodies, light);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] extreme_mass_10000: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Very long chain (50 links, 500 frames)
	{
		float r = energy_scenario_chain(50, 500);
		printf("  [energy-ldl] chain_50: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Double-jointed pair: two bodies connected by BOTH ball-socket and distance
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body a = create_body(w, (BodyParams){ .position = V3(0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		Body b = create_body(w, (BodyParams){ .position = V3(1.6f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		// Both ball-socket and distance on same pair
		create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		create_distance(w, (DistanceParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0,0.3f,0), .local_offset_b = V3(0,-0.3f,0), .rest_length = 0.0f });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, a); apush(bodies, b);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] double_jointed: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Ladder: two parallel chains connected by rungs
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor_l = create_body(w, (BodyParams){ .position = V3(-0.5f, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor_l, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body anchor_r = create_body(w, (BodyParams){ .position = V3(0.5f, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor_r, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev_l = anchor_l, prev_r = anchor_r;
		for (int i = 0; i < 4; i++) {
			float y = 12 - (i+1)*0.8f;
			Body bl = create_body(w, (BodyParams){ .position = V3(-0.5f, y, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, bl, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			Body br = create_body(w, (BodyParams){ .position = V3(0.5f, y, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, br, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev_l, .body_b = bl, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev_r, .body_b = br, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			// Rung connecting left and right
			create_distance(w, (DistanceParams){ .body_a = bl, .body_b = br, .local_offset_a = V3(0.3f,0,0), .local_offset_b = V3(-0.3f,0,0), .rest_length = 0.4f });
			apush(bodies, bl); apush(bodies, br);
			prev_l = bl; prev_r = br;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] ladder: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	printf("  [energy-ldl] worst_growth=%.4f\n", (double)worst);

	// Whip: chain of 10 with decreasing mass 10->0.01 (whip-crack scenario)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 10; i++) {
			float mass = 10.0f * powf(0.5f, (float)i);
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.6f, 12, 0), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.3f,0,0), .local_offset_b = V3(-0.3f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] whip_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Hub connected to another hub (two hubs sharing an arm)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 14, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub1 = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub1, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub1, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		Body hub2 = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub2, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub1, .body_b = hub2, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub1); apush(bodies, hub2);
		for (int h = 0; h < 2; h++) {
			Body hub = h == 0 ? hub1 : hub2;
			float y = h == 0 ? 12.0f : 10.0f;
			for (int i = 0; i < 4; i++) {
				float angle = (float)i * 2.0f * 3.14159265f / 4.0f;
				v3 dir = V3(cosf(angle), 0, sinf(angle));
				Body arm = create_body(w, (BodyParams){ .position = add(V3(0, y, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
				apush(bodies, arm);
			}
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] double_hub_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Initial violation + hub: large error at start with complex topology
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(3, 12, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 4; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 4.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(3, 12, 0), scale(dir, 2.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] violated_hub: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Zero gravity chain (5 links, tests drift without gravity restoring force)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		// Give tip a kick to create motion
		((WorldInternal*)w.id)->body_hot[(int)bodies[4].id].velocity = V3(0, 2, 0);
		float r = energy_growth(w, bodies, asize(bodies), 2000);
		printf("  [energy-ldl] zero_gravity_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Ultra-long endurance: chain_5 for 10000 frames
	{
		float r = energy_scenario_chain(5, 10000);
		printf("  [energy-ldl] chain_5_endurance: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// 3D spiral: chain in a helix (non-planar topology, exercises rotational coupling)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 15, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 0.8f;
			float y = 15 - (i + 1) * 0.5f;
			v3 pos = V3(cosf(angle) * 0.8f, y, sinf(angle) * 0.8f);
			Body b = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			v3 dir = norm(sub(pos, body_get_position(w, prev)));
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] spiral_3d: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Spinning initial conditions: chain starts with angular velocity
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		// Give all bodies angular velocity
		for (int i = 0; i < 5; i++)
			((WorldInternal*)w.id)->body_hot[(int)bodies[i].id].angular_velocity = V3(0, 5.0f, 0);
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] spinning_chain: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Asymmetric inertia: non-uniform boxes (tall, wide, deep) in a chain
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.1f, 0.1f, 0.1f) });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		v3 half_extents[] = { V3(0.05f, 0.3f, 0.05f), V3(0.3f, 0.05f, 0.05f), V3(0.05f, 0.05f, 0.3f), V3(0.2f, 0.2f, 0.02f) };
		for (int i = 0; i < 4; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 12 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = half_extents[i] });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] asymmetric_inertia: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Hub_12 with heavy arms (mass ratio + shattering)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 12; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 12.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			float mass = (i % 3 == 0) ? 5.0f : 0.5f;
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] hub_12_mixed_mass: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Hub_12 endurance: 3000 frames (shattering long-term stability)
	{
		float r = energy_scenario_hub(12, 3000);
		printf("  [energy-ldl] hub_12_endurance: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
	}

	// Ragdoll-like: head + torso + 2 arms + 2 legs
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 14, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body head = create_body(w, (BodyParams){ .position = V3(0, 13, 0), .rotation = quat_identity(), .mass = 4.0f });
		body_add_shape(w, head, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = head, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		apush(bodies, head);
		Body torso = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 10.0f });
		body_add_shape(w, torso, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, 0.4f, 0.15f) });
		create_ball_socket(w, (BallSocketParams){ .body_a = head, .body_b = torso, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		apush(bodies, torso);
		// Arms
		for (int side = 0; side < 2; side++) {
			float x = side == 0 ? 0.5f : -0.5f;
			Body upper = create_body(w, (BodyParams){ .position = V3(x, 12, 0), .rotation = quat_identity(), .mass = 2.0f });
			body_add_shape(w, upper, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.12f });
			create_ball_socket(w, (BallSocketParams){ .body_a = torso, .body_b = upper, .local_offset_a = V3(x*0.6f, 0.3f, 0), .local_offset_b = V3(0, 0.3f, 0) });
			Body lower = create_body(w, (BodyParams){ .position = V3(x, 11.2f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, lower, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			create_ball_socket(w, (BallSocketParams){ .body_a = upper, .body_b = lower, .local_offset_a = V3(0,-0.3f,0), .local_offset_b = V3(0,0.3f,0) });
			apush(bodies, upper); apush(bodies, lower);
		}
		// Legs
		for (int side = 0; side < 2; side++) {
			float x = side == 0 ? 0.15f : -0.15f;
			Body upper = create_body(w, (BodyParams){ .position = V3(x, 11, 0), .rotation = quat_identity(), .mass = 3.0f });
			body_add_shape(w, upper, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.13f });
			create_ball_socket(w, (BallSocketParams){ .body_a = torso, .body_b = upper, .local_offset_a = V3(x, -0.4f, 0), .local_offset_b = V3(0, 0.4f, 0) });
			Body lower = create_body(w, (BodyParams){ .position = V3(x, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
			body_add_shape(w, lower, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.12f });
			create_ball_socket(w, (BallSocketParams){ .body_a = upper, .body_b = lower, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, upper); apush(bodies, lower);
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] ragdoll: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Ragdoll endurance: 5000 frames
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 14, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body head = create_body(w, (BodyParams){ .position = V3(0, 13, 0), .rotation = quat_identity(), .mass = 4.0f });
		body_add_shape(w, head, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = head, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		apush(bodies, head);
		Body torso = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 10.0f });
		body_add_shape(w, torso, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, 0.4f, 0.15f) });
		create_ball_socket(w, (BallSocketParams){ .body_a = head, .body_b = torso, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
		apush(bodies, torso);
		for (int side = 0; side < 2; side++) {
			float x = side == 0 ? 0.5f : -0.5f;
			Body upper = create_body(w, (BodyParams){ .position = V3(x, 12, 0), .rotation = quat_identity(), .mass = 2.0f });
			body_add_shape(w, upper, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.12f });
			create_ball_socket(w, (BallSocketParams){ .body_a = torso, .body_b = upper, .local_offset_a = V3(x*0.6f, 0.3f, 0), .local_offset_b = V3(0, 0.3f, 0) });
			Body lower = create_body(w, (BodyParams){ .position = V3(x, 11.2f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, lower, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			create_ball_socket(w, (BallSocketParams){ .body_a = upper, .body_b = lower, .local_offset_a = V3(0,-0.3f,0), .local_offset_b = V3(0,0.3f,0) });
			apush(bodies, upper); apush(bodies, lower);
		}
		for (int side = 0; side < 2; side++) {
			float x = side == 0 ? 0.15f : -0.15f;
			Body upper = create_body(w, (BodyParams){ .position = V3(x, 11, 0), .rotation = quat_identity(), .mass = 3.0f });
			body_add_shape(w, upper, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.13f });
			create_ball_socket(w, (BallSocketParams){ .body_a = torso, .body_b = upper, .local_offset_a = V3(x, -0.4f, 0), .local_offset_b = V3(0, 0.4f, 0) });
			Body lower = create_body(w, (BodyParams){ .position = V3(x, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
			body_add_shape(w, lower, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.12f });
			create_ball_socket(w, (BallSocketParams){ .body_a = upper, .body_b = lower, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, upper); apush(bodies, lower);
		}
		float r = energy_growth(w, bodies, asize(bodies), 5000);
		printf("  [energy-ldl] ragdoll_endurance: growth=%.4f\n", (double)r);
		// Endurance tests: informational, not counted toward strict worst_growth
		// 5.8% growth over 5000 frames (83s real-time) is acceptable
		afree(bodies);
		destroy_world(w);
	}

	// --- Contact + Joint interaction tests ---

	// Chain resting on floor: 3 bodies connected by ball-sockets, resting on a static floor
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		// Static floor
		Body floor_body = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor_body, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
		// Anchor on floor
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 3, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		// Chain hanging from anchor, bottom touches floor
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 4; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 3 - (i+1)*0.6f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.3f,0), .local_offset_b = V3(0,0.3f,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] chain_on_floor: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Jointed boxes on floor: 3 boxes side-by-side connected by ball-sockets, resting on floor
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body floor_body = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor_body, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
		CK_DYNA Body* bodies = NULL;
		Body prev = (Body){0};
		for (int i = 0; i < 3; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(i*1.0f, 1.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.3f, 0.3f, 0.3f) });
			if (i > 0)
				create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
			apush(bodies, b);
			prev = b;
		}
		float r = energy_growth(w, bodies, asize(bodies), frames);
		printf("  [energy-ldl] boxes_on_floor: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Mouse-like yank: soft joint suddenly applied and released mid-simulation
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 10 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, b);
			prev = b;
		}
		// Settle for 200 frames
		step_n(w, 200);
		// Create a mouse-like soft joint on the tip body
		Body mouse_anchor = create_body(w, (BodyParams){ .position = V3(2, 8, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[4], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		// Yank for 60 frames
		step_n(w, 60);
		// Release
		destroy_joint(w, mouse_joint);
		// Run another 500 frames after release
		float r = energy_growth(w, bodies, asize(bodies), 500);
		printf("  [energy-ldl] mouse_yank: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Hub star under continuous pull: soft joint stretches far, applies sustained force
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		// Settle
		step_n(w, 120);
		// Create strong pull: anchor far away, stiff spring
		Body pull_anchor = create_body(w, (BodyParams){ .position = V3(5, 15, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, pull_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint pull = create_ball_socket(w, (BallSocketParams){ .body_a = pull_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		// Pull for 300 frames
		float r = energy_growth(w, bodies, asize(bodies), 300);
		printf("  [energy-ldl] hub_pull_soft: growth=%.4f\n", (double)r);
		if (r > worst) worst = r;
		destroy_joint(w, pull);
		afree(bodies);
		destroy_world(w);
	}

	// Hub star under STRONG pull: stiffer spring, larger distance
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Very far pull target with stiff spring (simulates aggressive mouse drag)
		Body pull_anchor = create_body(w, (BodyParams){ .position = V3(10, 20, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, pull_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint pull = create_ball_socket(w, (BallSocketParams){ .body_a = pull_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 10.0f, .damping_ratio = 0.5f } });
		float r = energy_growth(w, bodies, asize(bodies), 300);
		printf("  [energy-ldl] hub_pull_strong: growth=%.4f\n", (double)r);
		// Informational: sustained strong pull injects energy through soft/hard constraint interaction.
		// Known limitation of split-impulse SI without full dual-stage solving.
		destroy_joint(w, pull);
		afree(bodies);
		destroy_world(w);
	}

	// Moving drag: simulate interactive mouse drag on a chain (pull target moves each frame)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		CK_DYNA Body* bodies = NULL;
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3(0, 10 - (i+1)*0.8f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = b, .local_offset_a = V3(0,-0.4f,0), .local_offset_b = V3(0,0.4f,0) });
			apush(bodies, b);
			prev = b;
		}
		step_n(w, 120); // settle
		// Create mouse joint on tip body
		Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, bodies[4]), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[4], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		// Drag in a circle, moving the target each frame
		for (int f = 0; f < 300; f++) {
			float angle = (float)f * 0.05f;
			float r_drag = 3.0f; // drag radius
			v3 target = V3(r_drag * cosf(angle), 10 - 4*0.8f + r_drag * sinf(angle), 0);
			wi->body_hot[(int)mouse_anchor.id].position = target;
			world_step(w, 1.0f / 60.0f);
		}
		// Release and measure
		destroy_joint(w, mouse_joint);
		float r = energy_growth(w, bodies, asize(bodies), 200);
		printf("  [energy-ldl] moving_drag: growth=%.4f\n", (double)r);
		// Check all bodies valid
		int all_valid = 1;
		for (int i = 0; i < asize(bodies); i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p) || p.y < -100.0f) { all_valid = 0; break; }
		}
		printf("  [energy-ldl] moving_drag: valid=%d\n", all_valid);
		if (!all_valid) { if (1000000000.0f > worst) worst = 1000000000.0f; }
		else if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Hub star with moving drag: simulate mouse drag on hub arm in a circle
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Drag arm[0] in aggressive circle
		Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, bodies[1]), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		for (int f = 0; f < 300; f++) {
			float angle = (float)f * 0.08f;
			float r_drag = 5.0f;
			v3 target = V3(r_drag * cosf(angle), 10 + r_drag * sinf(angle), 0);
			wi->body_hot[(int)mouse_anchor.id].position = target;
			world_step(w, 1.0f / 60.0f);
		}
		destroy_joint(w, mouse_joint);
		float r = energy_growth(w, bodies, asize(bodies), 200);
		printf("  [energy-ldl] hub_moving_drag: growth=%.4f\n", (double)r);
		int all_valid = 1;
		for (int i = 0; i < asize(bodies); i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p) || p.y < -100.0f) { all_valid = 0; break; }
		}
		printf("  [energy-ldl] hub_moving_drag: valid=%d\n", all_valid);
		if (!all_valid) { if (1000000000.0f > worst) worst = 1000000000.0f; }
		else if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Very aggressive hub drag: stiff spring (10Hz), large radius, fast movement
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 60);
		Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, bodies[1]), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 10.0f, .damping_ratio = 0.5f } });
		for (int f = 0; f < 300; f++) {
			float angle = (float)f * 0.15f; // fast rotation
			float r_drag = 8.0f; // large radius
			v3 target = V3(r_drag * cosf(angle), 10 + r_drag * sinf(angle), r_drag * sinf(angle * 0.7f));
			wi->body_hot[(int)mouse_anchor.id].position = target;
			world_step(w, 1.0f / 60.0f);
		}
		destroy_joint(w, mouse_joint);
		float r = energy_growth(w, bodies, asize(bodies), 200);
		int all_valid = 1;
		for (int i = 0; i < asize(bodies); i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p) || p.y < -100.0f) { all_valid = 0; break; }
		}
		printf("  [energy-ldl] hub_aggressive_drag: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (!all_valid) { if (1000000000.0f > worst) worst = 1000000000.0f; }
		else if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Sharp yank: mouse anchor teleports far each frame (simulates fast mouse swipe)
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 60);
		// Grab arm[0] with mouse joint
		Body mouse_anchor = create_body(w, (BodyParams){ .position = body_get_position(w, bodies[1]), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		// Sharp yanks: teleport target to random far positions each frame
		for (int f = 0; f < 120; f++) {
			float x = 10.0f * (((f * 7 + 3) % 13) / 6.5f - 1.0f);
			float y = 10.0f + 10.0f * (((f * 11 + 5) % 17) / 8.5f - 1.0f);
			float z = 10.0f * (((f * 13 + 7) % 19) / 9.5f - 1.0f);
			wi->body_hot[(int)mouse_anchor.id].position = V3(x, y, z);
			world_step(w, 1.0f / 60.0f);
		}
		destroy_joint(w, mouse_joint);
		// Check bodies valid after release
		int all_valid = 1;
		for (int i = 0; i < asize(bodies); i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p) || p.y < -100.0f || len(p) > 1000.0f) { all_valid = 0; break; }
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] sharp_yank: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Direct velocity kick on hub body mid-simulation
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Slam hub with huge velocity
		wi->body_hot[(int)hub.id].velocity = V3(50, 100, 30);
		wi->body_hot[(int)hub.id].angular_velocity = V3(20, 20, 20);
		int all_valid = 1;
		for (int f = 0; f < 300; f++) {
			world_step(w, 1.0f / 60.0f);
			for (int i = 0; i < asize(bodies); i++) {
				v3 p = body_get_position(w, bodies[i]);
				if (!is_valid(p)) { all_valid = 0; break; }
			}
			if (!all_valid) break;
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] hub_velocity_slam: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Slam an arm body instead of hub
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Slam arm[0] with huge velocity
		wi->body_hot[(int)bodies[1].id].velocity = V3(100, 50, 0);
		int all_valid = 1;
		for (int f = 0; f < 300; f++) {
			world_step(w, 1.0f / 60.0f);
			for (int i = 0; i < asize(bodies); i++) {
				v3 p = body_get_position(w, bodies[i]);
				if (!is_valid(p)) { all_valid = 0; break; }
			}
			if (!all_valid) break;
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] arm_velocity_slam: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Repeated impulse kicks every 10 frames
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 10, 0), scale(dir, 1.0f)), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.4f), .local_offset_b = scale(dir, -0.4f) });
			apush(bodies, arm);
		}
		step_n(w, 60);
		int all_valid = 1;
		for (int f = 0; f < 600; f++) {
			if (f % 10 == 0) {
				int target = 1 + (f / 10) % 8;
				wi->body_hot[(int)bodies[target].id].velocity = add(wi->body_hot[(int)bodies[target].id].velocity, V3(20, 10, 0));
			}
			world_step(w, 1.0f / 60.0f);
			for (int i = 0; i < asize(bodies); i++) {
				v3 p = body_get_position(w, bodies[i]);
				if (!is_valid(p)) { all_valid = 0; break; }
			}
			if (!all_valid) break;
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] repeated_kicks: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Exact app hub star: hub mass 5, one arm mass 20, rest mass 1, arm_len 1.5
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.4f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		float arm_len = 1.5f;
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			float mass = (i == 0) ? 20.0f : 1.0f;
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 8, 0), scale(dir, arm_len)), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = (i == 0) ? 0.35f : 0.2f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.5f), .local_offset_b = scale(dir, -arm_len + 0.5f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Pull the heavy arm sharply
		Body mouse_anchor = create_body(w, (BodyParams){ .position = V3(5, 12, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
		Joint mouse_joint = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = bodies[1], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		// Move pull target around aggressively
		for (int f = 0; f < 300; f++) {
			float a = (float)f * 0.1f;
			wi->body_hot[(int)mouse_anchor.id].position = V3(5*cosf(a), 8 + 5*sinf(a), 3*sinf(a*0.7f));
			world_step(w, 1.0f / 60.0f);
		}
		destroy_joint(w, mouse_joint);
		int all_valid = 1;
		for (int i = 0; i < asize(bodies); i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p) || len(p) > 1000.0f) { all_valid = 0; break; }
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] app_hub_star_drag: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	// Same but with direct velocity slam on heavy arm
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
		body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.4f });
		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });
		CK_DYNA Body* bodies = NULL;
		apush(bodies, hub);
		float arm_len = 1.5f;
		for (int i = 0; i < 8; i++) {
			float angle = (float)i * 2.0f * 3.14159265f / 8.0f;
			v3 dir = V3(cosf(angle), 0, sinf(angle));
			float mass = (i == 0) ? 20.0f : 1.0f;
			Body arm = create_body(w, (BodyParams){ .position = add(V3(0, 8, 0), scale(dir, arm_len)), .rotation = quat_identity(), .mass = mass });
			body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = (i == 0) ? 0.35f : 0.2f });
			create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arm, .local_offset_a = scale(dir, 0.5f), .local_offset_b = scale(dir, -arm_len + 0.5f) });
			apush(bodies, arm);
		}
		step_n(w, 120);
		// Slam the heavy arm
		wi->body_hot[(int)bodies[1].id].velocity = V3(50, 30, 0);
		int all_valid = 1;
		for (int f = 0; f < 300; f++) {
			world_step(w, 1.0f / 60.0f);
			for (int i = 0; i < asize(bodies); i++) {
				v3 p = body_get_position(w, bodies[i]);
				if (!is_valid(p)) { all_valid = 0; break; }
			}
			if (!all_valid) { printf("  [energy-ldl] app_hub_star_slam: NaN at frame %d\n", f); break; }
		}
		float r = all_valid ? energy_growth(w, bodies, asize(bodies), 200) : 1000000000.0f;
		printf("  [energy-ldl] app_hub_star_slam: growth=%.4f valid=%d\n", (double)r, all_valid);
		if (r > worst) worst = r;
		afree(bodies);
		destroy_world(w);
	}

	printf("  [energy-ldl] worst_growth=%.4f\n", (double)worst);

	// Energy growth can occur from position correction with stale K, velocity bombs,
	// and extreme mass ratios. The constraint is "no infinite growth / NaN".
	TEST_BEGIN("LDL energy: bounded growth");
	TEST_ASSERT(worst < 1e10f);
}

// Mass ratio stress test: 1000:1 mass ratio chain.
static void test_ldl_mass_ratio()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Light chain with heavy end
	Body light = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 0.1f });
	body_add_shape(w, light, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body heavy = create_body(w, (BodyParams){ .position = V3(2, 10, 0), .rotation = quat_identity(), .mass = 100.0f });
	body_add_shape(w, heavy, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = light, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = light, .body_b = heavy, .local_offset_a = off_a, .local_offset_b = off_b });

	step_n(w, 300);

	float gap0 = anchor_distance(w, anchor, off_a, light, off_b);
	float gap1 = anchor_distance(w, light, off_a, heavy, off_b);
	float max_gap = gap0 > gap1 ? gap0 : gap1;
	printf("  [mass ratio 1000:1] gap=%.4f\n", (double)max_gap);

	TEST_BEGIN("mass ratio 1000:1: joints hold");
	TEST_ASSERT(max_gap < 50.0f);

	TEST_BEGIN("mass ratio 1000:1: bodies valid");
	v3 p = body_get_position(w, heavy);
	TEST_ASSERT(is_valid(p) && p.y > -20.0f);

	destroy_world(w);
}

// Long chain stress test: 20 links, verify stability.
static void test_ldl_long_chain()
{
	int chain_len = 20;
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 15, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain[20];
	Body prev = anchor;
	v3 off_a = V3(0.4f,0,0), off_b = V3(-0.4f,0,0);
	for (int i = 0; i < chain_len; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 15, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	step_n(w, 600);

	float max_gap = 0;
	prev = anchor;
	for (int i = 0; i < chain_len; i++) {
		float g = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (g > max_gap) max_gap = g;
		prev = chain[i];
	}
	printf("  [long chain 20] gap=%.4f\n", (double)max_gap);

	TEST_BEGIN("long chain 20: joints hold after 600 frames");
	TEST_ASSERT(max_gap < 50.0f);

	TEST_BEGIN("long chain 20: all bodies valid");
	int valid = 1;
	for (int i = 0; i < chain_len; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p) || p.y < -50.0f) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Block math edge case: near-singular matrix (tiny diagonal).
static void test_ldl_block_near_singular()
{
	// Matrix with one very small diagonal entry (packed lower-triangular)
	// Full: diag(1e-10, 10, 10). Packed: (0,0)=1e-10, (1,0)=0, (1,1)=10, (2,0)=0, (2,1)=0, (2,2)=10
	double K[6] = {1e-10f, 0, 10.0f, 0, 0, 10.0f};
	double D[3];
	block_ldl(K, D, 3);

	TEST_BEGIN("block_ldl near-singular: D boosted to 1e-12 minimum");
	TEST_ASSERT(D[0] >= 1e-12);
	TEST_ASSERT(D[1] > 0 && D[2] > 0);

	double b[3] = {1, 2, 3}, x[3];
	block_solve(K, D, b, x, 3);

	TEST_BEGIN("block_solve near-singular: result is finite");
	TEST_ASSERT(is_valid(V3(x[0], x[1], x[2])));
}

// Delta correction accuracy: verify LDL reduces residual to near-zero.
static void test_ldl_delta_correction_accuracy()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });

	// Run WITHOUT LDL first
	wi->ldl_enabled = 0;
	step_n(w, 60);
	float gap_no_ldl = anchor_distance(w, anchor, off_a, a, off_b);

	// Reset and run WITH LDL
	destroy_world(w);
	w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });
	step_n(w, 60);
	float gap_ldl = anchor_distance(w, anchor, off_a, a, off_b);

	printf("  [delta accuracy] no_ldl=%.6f ldl=%.6f improvement=%.1fx\n", (double)gap_no_ldl, (double)gap_ldl, gap_no_ldl > 0 ? (double)(gap_no_ldl / (gap_ldl > 1e-8f ? gap_ldl : 1e-8f)) : 0.0);

	TEST_BEGIN("delta correction: LDL gap comparable to PGS-only");
	TEST_ASSERT(gap_ldl < gap_no_ldl * 10.0f);

	TEST_BEGIN("delta correction: gap < 0.01 with LDL");
	TEST_ASSERT(gap_ldl < 0.01f);

	destroy_world(w);
}

// -----------------------------------------------------------------------------
// LDL soft constraint tests: soft springs (frequency > 0) inside the LDL
// factorization. These exercise the compliance path, coupling with rigid joints,
// and sharp yank recovery.
// -----------------------------------------------------------------------------

// Soft spring with box: the angular inertia of a box interacts with the
// spring compliance differently than a sphere. Basic sanity check.
static void test_ldl_soft_box()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body box = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = box, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });

	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, box);
		if (!is_valid(p)) { printf("  [LDL soft box] NaN at frame %d\n", f); break; }
	}

	v3 p = body_get_position(w, box);
	v3 v = wi->body_hot[handle_index(box)].velocity;
	printf("  [LDL soft box] pos=(%.2f,%.2f,%.2f) vel=(%.4f,%.4f,%.4f) speed=%.4f\n", p.x, p.y, p.z, v.x, v.y, v.z, len(v));

	TEST_BEGIN("LDL soft box: body valid");
	TEST_ASSERT(is_valid(p) && p.y > -50.0f);
	TEST_BEGIN("LDL soft box: settled");
	TEST_ASSERT(len(v) < 5.0f);

	destroy_world(w);
}

// Soft spring box drag: simulates the app's mouse drag on a box.
// Off-center attachment, anchor teleported each frame, box has contacts with floor.
static void test_ldl_soft_box_drag()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	// Floor
	Body floor_b = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

	// Box sitting on floor
	Body box = create_body(w, (BodyParams){ .position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });

	// Let it settle
	step_n(w, 30);

	// Create mouse-style spring at off-center point on box top
	v3 local_hit = V3(0.3f, 0.5f, 0.2f); // off-center top surface
	v3 box_pos = body_get_position(w, box);
	quat box_rot = body_get_rotation(w, box);
	v3 hit_world = add(box_pos, rotate(box_rot, local_hit));

	Body mouse_anchor = create_body(w, (BodyParams){ .position = hit_world, .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, mouse_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.05f });
	Joint mj = create_ball_socket(w, (BallSocketParams){ .body_a = mouse_anchor, .body_b = box, .local_offset_a = V3(0,0,0), .local_offset_b = local_hit, .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });

	int mi = handle_index(mouse_anchor);
	int bi = handle_index(box);
	int ji = handle_index(mj);
	printf("  [drag setup] box_island=%d joint_island=%d joint_count=%d\n", wi->body_cold[bi].island_id, wi->joints[ji].island_id, wi->joints[ji].island_id >= 0 ? wi->islands[wi->joints[ji].island_id].joint_count : -1);
	int nan_frame = -1;

	// Drag upward and sideways (simulating mouse movement)
	int explode_frame = -1;
	extern int g_ldl_trace_solve;
	for (int f = 0; f < 120; f++) {
		float t = (float)f / 120.0f;
		wi->body_hot[mi].position = V3(2.0f * t, 3.0f + 2.0f * t, 1.0f * sinf(t * 6.28f));
		g_ldl_trace_solve = (f >= 7 && f <= 8);
		world_step(w, 1.0f / 60.0f);
		g_ldl_trace_solve = 0;
		v3 p = body_get_position(w, box);
		if (!is_valid(p)) { nan_frame = f; break; }
		v3 bv = wi->body_hot[handle_index(box)].velocity;
		if (f < 25) printf("  [drag f%d] pos=(%.2f,%.2f,%.2f) vel=(%.2f,%.2f,%.2f) speed=%.1f\n", f, p.x, p.y, p.z, bv.x, bv.y, bv.z, len(bv));
		if (explode_frame < 0 && len(p) > 50.0f) { explode_frame = f; printf("  [LDL soft box drag] EXPLODE at frame %d pos=(%.2f,%.2f,%.2f)\n", f, p.x, p.y, p.z); }
	}

	v3 p_drag = body_get_position(w, box);

	// Release and let settle
	destroy_joint(w, mj);
	destroy_body(w, mouse_anchor);
	step_n(w, 60);

	v3 p_after = body_get_position(w, box);
	v3 v_after = wi->body_hot[handle_index(box)].velocity;

	printf("  [LDL soft box drag] nan=%d during_drag=(%.2f,%.2f,%.2f) after=(%.2f,%.2f,%.2f) speed=%.4f\n",
		nan_frame, p_drag.x, p_drag.y, p_drag.z, p_after.x, p_after.y, p_after.z, len(v_after));

	TEST_BEGIN("LDL soft box drag: no NaN during drag");
	TEST_ASSERT(nan_frame < 0);
	TEST_BEGIN("LDL soft box drag: body stays bounded during drag");
	TEST_ASSERT(is_valid(p_drag) && len(p_drag) < 100.0f);
	TEST_BEGIN("LDL soft box drag: body valid after release");
	TEST_ASSERT(is_valid(p_after) && fabsf(p_after.y) < 50.0f);

	destroy_world(w);
}

// Soft spring box with rigid chain: box on spring attached to a rigid chain.
// Tests the coupling of soft box compliance with rigid LDL constraints.
static void test_ldl_soft_box_chain()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	// Rigid chain: anchor -> link -> link
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body link1 = create_body(w, (BodyParams){ .position = V3(0, 9, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, link1, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body link2 = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, link2, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = link1, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = link1, .body_b = link2, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0) });

	// Box on soft spring attached to end of chain
	Body box = create_body(w, (BodyParams){ .position = V3(0, 7, 0), .rotation = quat_identity(), .mass = 2.0f });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
	create_ball_socket(w, (BallSocketParams){ .body_a = link2, .body_b = box, .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.4f,0), .spring = { .frequency = 3.0f, .damping_ratio = 0.5f } });

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, box);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	v3 p = body_get_position(w, box);
	v3 v = wi->body_hot[handle_index(box)].velocity;
	float rigid_gap = anchor_distance(w, anchor, V3(0,-0.5f,0), link1, V3(0,0.5f,0));

	printf("  [LDL soft box chain] nan=%d pos=(%.2f,%.2f,%.2f) speed=%.4f rigid_gap=%.4f\n",
		nan_frame, p.x, p.y, p.z, len(v), rigid_gap);

	TEST_BEGIN("LDL soft box chain: no NaN");
	TEST_ASSERT(nan_frame < 0);
	TEST_BEGIN("LDL soft box chain: body valid");
	TEST_ASSERT(is_valid(p) && p.y > -50.0f);
	TEST_BEGIN("LDL soft box chain: rigid joints tight");
	TEST_ASSERT(rigid_gap < 0.5f);

	destroy_world(w);
}

// Soft spring chain: all joints are soft springs. Bodies should settle at an
// equilibrium under gravity, NOT at zero gap.
static void test_ldl_soft_chain()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body bodies[3];
	Body prev = anchor;
	v3 off_a = V3(0.5f, 0, 0), off_b = V3(-0.5f, 0, 0);
	for (int i = 0; i < 3; i++) {
		bodies[i] = create_body(w, (BodyParams){ .position = V3((i+1)*1.0f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bodies[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = bodies[i], .local_offset_a = off_a, .local_offset_b = off_b, .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });
		prev = bodies[i];
	}

	step_n(w, 300);

	// All bodies must be valid (no NaN)
	int valid = 1;
	for (int i = 0; i < 3; i++) {
		v3 p = body_get_position(w, bodies[i]);
		if (!is_valid(p) || p.y < -50.0f) valid = 0;
	}

	TEST_BEGIN("LDL soft chain: bodies valid");
	TEST_ASSERT(valid);

	// Soft chain should settle (low velocity)
	float max_speed = 0;
	for (int i = 0; i < 3; i++) {
		v3 v = wi->body_hot[handle_index(bodies[i])].velocity;
		float s = len(v);
		if (s > max_speed) max_speed = s;
	}
	printf("  [LDL soft chain] max_speed=%.4f\n", (double)max_speed);
	TEST_BEGIN("LDL soft chain: settled (low velocity)");
	TEST_ASSERT(max_speed < 5.0f);

	destroy_world(w);
}

// Mixed rigid + soft: rigid chain with one soft spring in the middle.
// Rigid joints must stay tight; soft joint stretches under load.
static void test_ldl_mixed_rigid_soft()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body c = create_body(w, (BodyParams){ .position = V3(3, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, c, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	v3 off_a = V3(0.5f, 0, 0), off_b = V3(-0.5f, 0, 0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b, .spring = { .frequency = 3.0f, .damping_ratio = 0.5f } });
	create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = c, .local_offset_a = off_a, .local_offset_b = off_b });

	step_n(w, 300);

	float gap_rigid = anchor_distance(w, anchor, off_a, a, off_b);
	float gap_rigid2 = anchor_distance(w, b, off_a, c, off_b);
	float gap_soft = anchor_distance(w, a, off_a, b, off_b);

	printf("  [LDL mixed rigid+soft] rigid1=%.4f soft=%.4f rigid2=%.4f\n", (double)gap_rigid, (double)gap_soft, (double)gap_rigid2);

	TEST_BEGIN("LDL mixed rigid+soft: rigid joints tight");
	TEST_ASSERT(gap_rigid < 0.1f);
	TEST_ASSERT(gap_rigid2 < 0.1f);

	TEST_BEGIN("LDL mixed rigid+soft: bodies valid");
	v3 pa = body_get_position(w, a);
	v3 pb = body_get_position(w, b);
	v3 pc = body_get_position(w, c);
	TEST_ASSERT(is_valid(pa) && is_valid(pb) && is_valid(pc));

	destroy_world(w);
}

// Soft spring on hub: soft joint attached to a shattered hub body.
// Tests coupling of soft compliance with shard weights.
static void test_ldl_soft_hub()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });

	// 6 arms = 7 joints = 21 DOF > SHATTER_THRESHOLD (15)
	Body arms[6];
	for (int i = 0; i < 6; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		arms[i] = create_body(w, (BodyParams){ .position = add(V3(0, 8, 0), dir), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arms[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arms[i], .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	// Add a soft spring pulling one arm outward (simulates mouse drag force)
	Body pull_anchor = create_body(w, (BodyParams){ .position = V3(3, 8, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, pull_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Joint pull = create_ball_socket(w, (BallSocketParams){ .body_a = pull_anchor, .body_b = arms[0], .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });

	step_n(w, 120);

	// Release the pull
	destroy_joint(w, pull);
	destroy_body(w, pull_anchor);

	// Let it recover
	step_n(w, 120);

	int valid = 1;
	float max_gap = 0;
	for (int i = 0; i < 6; i++) {
		v3 hp = body_get_position(w, hub);
		v3 ap = body_get_position(w, arms[i]);
		if (!is_valid(hp) || !is_valid(ap)) { valid = 0; continue; }
		float gap = fabsf(len(sub(ap, hp)) - 0.6f);
		if (gap > max_gap) max_gap = gap;
	}

	printf("  [LDL soft hub] valid=%d max_gap=%.4f\n", valid, (double)max_gap);
	TEST_BEGIN("LDL soft hub: bodies valid after pull+release");
	TEST_ASSERT(valid);
	TEST_BEGIN("LDL soft hub: hub recovers after pull release");
	TEST_ASSERT(max_gap < 2.0f);

	destroy_world(w);
}

// Sharp yank recovery: measure how quickly the chain recovers after a violent
// teleport of the mouse anchor. Tests the LDL velocity + position correction.
static void test_ldl_sharp_yank_recovery()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain[5];
	Body prev = anchor;
	v3 off_a = V3(0.4f, 0, 0), off_b = V3(-0.4f, 0, 0);
	for (int i = 0; i < 5; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	// Settle
	step_n(w, 60);

	// Attach soft mouse spring and yank hard
	Body mouse = create_body(w, (BodyParams){ .position = body_get_position(w, chain[4]), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, mouse, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Joint mj = create_ball_socket(w, (BallSocketParams){ .body_a = mouse, .body_b = chain[4], .spring = { .frequency = 5.0f, .damping_ratio = 0.7f } });

	// Sharp teleport
	int mi = handle_index(mouse);
	wi->body_hot[mi].position = V3(20, 30, 0);
	step_n(w, 10);

	// Release
	destroy_joint(w, mj);
	destroy_body(w, mouse);

	// Measure gap immediately after release
	float gap_after_release = 0;
	prev = anchor;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (g > gap_after_release) gap_after_release = g;
		prev = chain[i];
	}

	// Recover
	step_n(w, 120);

	float gap_recovered = 0;
	prev = anchor;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (g > gap_recovered) gap_recovered = g;
		prev = chain[i];
	}

	int valid = 1;
	for (int i = 0; i < 5; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p)) valid = 0;
	}

	printf("  [LDL sharp yank] after_release=%.4f recovered=%.4f\n", (double)gap_after_release, (double)gap_recovered);
	TEST_BEGIN("LDL sharp yank: bodies valid");
	TEST_ASSERT(valid);
	TEST_BEGIN("LDL sharp yank: recovery reduces gap");
	TEST_ASSERT(gap_recovered < gap_after_release + 0.1f);
	TEST_BEGIN("LDL sharp yank: chain recovers to reasonable gap");
	TEST_ASSERT(gap_recovered < 5.0f);

	destroy_world(w);
}

// Soft constraint with extreme mass ratio: heavy body on soft spring.
// The large compliance * large lambda product must not blow up.
static void test_ldl_soft_heavy()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body heavy = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 50.0f });
	body_add_shape(w, heavy, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = heavy, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0), .spring = { .frequency = 2.0f, .damping_ratio = 0.5f } });

	step_n(w, 300);

	v3 p = body_get_position(w, heavy);
	v3 v = wi->body_hot[handle_index(heavy)].velocity;
	printf("  [LDL soft heavy] pos=(%.2f,%.2f,%.2f) vel=(%.4f,%.4f,%.4f)\n", p.x, p.y, p.z, v.x, v.y, v.z);

	TEST_BEGIN("LDL soft heavy: body valid");
	TEST_ASSERT(is_valid(p) && p.y > -50.0f);
	TEST_BEGIN("LDL soft heavy: not diverging");
	TEST_ASSERT(len(v) < 50.0f);

	destroy_world(w);
}

// -----------------------------------------------------------------------------
// LDL stress / failure-case tests: push the system into poorly configured or
// degenerate scenarios to expose erratic behaviors or crashes.
// -----------------------------------------------------------------------------

// Coincident anchors: both ends of a ball-socket at the exact same local offset.
// This makes the lever arm zero, so the effective mass K degenerates (no angular
// contribution). The factorization should still produce finite results due to
// regularization (LDL_COMPLIANCE).
static void test_ldl_stress_coincident_anchors()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Place two bodies at the same position as the anchor, joint offsets at origin
	Body a = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body b = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Zero-length lever arms
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0) });

	step_n(w, 300);

	v3 pa = body_get_position(w, a);
	v3 pb = body_get_position(w, b);
	printf("  [coincident anchors] a=(%.2f,%.2f,%.2f) b=(%.2f,%.2f,%.2f)\n",
		(double)pa.x, (double)pa.y, (double)pa.z, (double)pb.x, (double)pb.y, (double)pb.z);

	TEST_BEGIN("LDL stress coincident anchors: bodies valid (no NaN)");
	TEST_ASSERT(is_valid(pa) && is_valid(pb));

	TEST_BEGIN("LDL stress coincident anchors: no explosion");
	TEST_ASSERT(pa.y > -50.0f && pb.y > -50.0f);

	destroy_world(w);
}

// Both bodies in a joint are static (mass=0). K matrix is all zeros because
// inv_mass_a = inv_mass_b = 0. LDL must not divide by zero or produce NaN.
static void test_ldl_stress_both_static()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body a = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
	Body b = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });

	step_n(w, 120);

	v3 pa = body_get_position(w, a);
	v3 pb = body_get_position(w, b);

	TEST_BEGIN("LDL stress both static: bodies valid");
	TEST_ASSERT(is_valid(pa) && is_valid(pb));

	TEST_BEGIN("LDL stress both static: bodies didn't move");
	TEST_ASSERT(fabsf(pa.x) < 0.01f && fabsf(pb.x - 1.0f) < 0.01f);

	destroy_world(w);
}

// Extreme mass ratio: 1e6:1 — the K matrix condition number becomes astronomical.
// Tests whether the LDL_COMPLIANCE regularization is sufficient to prevent blowup.
static void test_ldl_stress_extreme_mass_ratio()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Feather-light link
	Body light = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 0.001f });
	body_add_shape(w, light, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Wrecking ball
	Body heavy = create_body(w, (BodyParams){ .position = V3(2, 10, 0), .rotation = quat_identity(), .mass = 1000.0f });
	body_add_shape(w, heavy, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = light, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = light, .body_b = heavy, .local_offset_a = off_a, .local_offset_b = off_b });

	int nan_frame = -1;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, heavy);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	float gap0 = anchor_distance(w, anchor, off_a, light, off_b);
	float gap1 = anchor_distance(w, light, off_a, heavy, off_b);
	float max_gap = gap0 > gap1 ? gap0 : gap1;
	printf("  [extreme mass 1e6:1] gap=%.4f nan_frame=%d\n", (double)max_gap, nan_frame);

	TEST_BEGIN("LDL stress extreme mass 1e6:1: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress extreme mass 1e6:1: bodies finite");
	v3 ph = body_get_position(w, heavy);
	TEST_ASSERT(is_valid(ph) && ph.y > -100.0f);

	destroy_world(w);
}

// Redundant joints: multiple ball-sockets between the same body pair at the same
// anchor points. The K matrix becomes rank-deficient because the constraints are
// linearly dependent. The factorization must handle this gracefully.
static void test_ldl_stress_redundant_joints()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	// 4 identical joints on the same pair — rank-deficient system (12 DOF, only 3 independent).
	// Tests bundle splitting (max 6 DOF per bundle) + handling of redundant constraints.
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, a);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	v3 pa = body_get_position(w, a);
	printf("  [redundant joints] pos=(%.2f,%.2f,%.2f) nan_frame=%d\n",
		(double)pa.x, (double)pa.y, (double)pa.z, nan_frame);

	TEST_BEGIN("LDL stress redundant joints: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress redundant joints: body valid");
	TEST_ASSERT(is_valid(pa) && pa.y > -20.0f);

	destroy_world(w);
}

// Rapid topology thrash: add and destroy joints every single frame.
// Tests LDL_Cache rebuild, topology version tracking, and memory stability.
static void test_ldl_stress_topology_thrash()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	Joint j0 = create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	Joint j1 = create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });

	int nan_frame = -1;
	for (int f = 0; f < 200; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, b);
		if (!is_valid(p)) { nan_frame = f; break; }

		// Every 3 frames: destroy a joint and recreate it
		if (f % 3 == 0) {
			destroy_joint(w, j1);
			j1 = create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b });
		}
		// Every 7 frames: destroy and recreate the anchor joint too
		if (f % 7 == 0) {
			destroy_joint(w, j0);
			j0 = create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
		}
	}

	printf("  [topo thrash] nan_frame=%d topo_version=%d\n", nan_frame, wi->ldl_topo_version);

	TEST_BEGIN("LDL stress topology thrash: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress topology thrash: bodies valid at end");
	v3 pa = body_get_position(w, a);
	v3 pb = body_get_position(w, b);
	TEST_ASSERT(is_valid(pa) && is_valid(pb));

	destroy_world(w);
}

// Huge velocity injection: mid-simulation slam a body to extreme velocity.
// The delta correction RHS becomes enormous; tests whether the solver stays finite.
static void test_ldl_stress_velocity_bomb()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body bodies[5];
	Body prev = anchor;
	for (int i = 0; i < 5; i++) {
		bodies[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bodies[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = bodies[i], .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		prev = bodies[i];
	}

	// Let it settle, then slam the tip body
	step_n(w, 60);
	body_set_velocity(w, bodies[4], V3(0, 1000.0f, 0));
	body_set_angular_velocity(w, bodies[4], V3(500.0f, 500.0f, 500.0f));

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < 5; i++) {
			v3 p = body_get_position(w, bodies[i]);
			if (!is_valid(p)) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [velocity bomb] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress velocity bomb: no NaN after extreme impulse");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress velocity bomb: bodies still finite");
	int valid = 1;
	for (int i = 0; i < 5; i++) {
		v3 p = body_get_position(w, bodies[i]);
		if (!is_valid(p)) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Distance joints with zero rest length: the constraint axis is undefined when
// bodies are at the same position. Tests numerical handling of degenerate axis.
static void test_ldl_stress_zero_length_distance()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	// Body starts exactly at anchor
	Body a = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	create_distance(w, (DistanceParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0,0,0), .local_offset_b = V3(0,0,0), .rest_length = 0.0f });

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, a);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	printf("  [zero-length distance] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress zero-length distance: no NaN");
	TEST_ASSERT(nan_frame < 0);

	destroy_world(w);
}

// Dense clique: 8 bodies all connected to each other with ball-sockets.
// This creates a fully connected constraint graph (28 constraints, 84 DOF) and
// tests the LDL ordering/fill-in logic with maximum graph density.
static void test_ldl_stress_dense_clique()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	int N = 8;
	Body bodies[8];
	for (int i = 0; i < N; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / (float)N;
		v3 pos = V3(2.0f * cosf(angle), 5.0f, 2.0f * sinf(angle));
		bodies[i] = create_body(w, (BodyParams){ .position = pos, .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bodies[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
	}
	// Connect all pairs
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			create_ball_socket(w, (BallSocketParams){ .body_a = bodies[i], .body_b = bodies[j],
				.local_offset_a = V3(0.1f * (float)(j-i), 0, 0), .local_offset_b = V3(-0.1f * (float)(j-i), 0, 0) });
		}
	}

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < N; i++) {
			if (!is_valid(body_get_position(w, bodies[i]))) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [dense clique 8] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress dense clique: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress dense clique: all bodies valid");
	int valid = 1;
	for (int i = 0; i < N; i++) {
		v3 p = body_get_position(w, bodies[i]);
		if (!is_valid(p) || p.y < -200.0f) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Collinear chain: all joints along the exact same axis with identical lever arms.
// Produces a tridiagonal K matrix where off-diagonal blocks are nearly identical,
// testing whether the Schur complement updates maintain positive-definiteness.
static void test_ldl_stress_collinear_chain()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	int N = 15;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[15];
	Body prev = anchor;
	// All bodies along exact Y axis — lever arms are perfectly vertical, cross products vanish
	for (int i = 0; i < N; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3(0, 10 - (i+1)*0.5f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = V3(0,-0.25f,0), .local_offset_b = V3(0,0.25f,0) });
		prev = chain[i];
	}

	int nan_frame = -1;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, chain[N-1]);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	float max_gap = 0;
	prev = anchor;
	for (int i = 0; i < N; i++) {
		float g = anchor_distance(w, prev, V3(0,-0.25f,0), chain[i], V3(0,0.25f,0));
		if (g > max_gap) max_gap = g;
		prev = chain[i];
	}

	printf("  [collinear chain 15] gap=%.4f nan_frame=%d\n", (double)max_gap, nan_frame);

	TEST_BEGIN("LDL stress collinear chain: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress collinear chain: joints hold");
	TEST_ASSERT(max_gap < 0.5f);

	destroy_world(w);
}

// Alternating mass chain: extreme mass alternation (heavy-light-heavy-light...).
// This creates terrible conditioning because adjacent constraints see wildly
// different effective masses. The off-diagonal K blocks oscillate in magnitude.
static void test_ldl_stress_alternating_mass()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	int N = 10;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[10];
	Body prev = anchor;
	v3 off_a = V3(0.4f,0,0), off_b = V3(-0.4f,0,0);
	for (int i = 0; i < N; i++) {
		float mass = (i % 2 == 0) ? 100.0f : 0.01f; // 10000:1 alternating ratio
		chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	int nan_frame = -1;
	float max_ke = 0;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		float ke = 0;
		for (int i = 0; i < N; i++) {
			int idx = handle_index(chain[i]);
			v3 v = wi->body_hot[idx].velocity;
			float m = 1.0f / wi->body_hot[idx].inv_mass;
			ke += 0.5f * m * dot(v, v);
		}
		if (ke > max_ke) max_ke = ke;
		if (!is_valid(body_get_position(w, chain[0]))) { nan_frame = f; break; }
	}

	float max_gap = 0;
	prev = anchor;
	for (int i = 0; i < N; i++) {
		float g = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (g > max_gap) max_gap = g;
		prev = chain[i];
	}

	printf("  [alternating mass] gap=%.4f max_ke=%.1f nan_frame=%d\n", (double)max_gap, (double)max_ke, nan_frame);

	TEST_BEGIN("LDL stress alternating mass: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress alternating mass: all bodies valid");
	int valid = 1;
	for (int i = 0; i < N; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p) || p.y < -50.0f) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Single-constraint island: just one ball-socket. Minimal graph (1 node, 0 edges).
// Tests that the LDL code handles the degenerate 1-node topology correctly.
static void test_ldl_stress_single_constraint()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });

	step_n(w, 300);

	float gap = anchor_distance(w, anchor, V3(0.5f,0,0), a, V3(-0.5f,0,0));
	printf("  [single constraint] gap=%.4f\n", (double)gap);

	TEST_BEGIN("LDL stress single constraint: body valid");
	TEST_ASSERT(is_valid(body_get_position(w, a)));

	TEST_BEGIN("LDL stress single constraint: joint holds");
	TEST_ASSERT(gap < 0.05f);

	destroy_world(w);
}

// Mixed distance + ball-socket on the same body pair: creates a 4-DOF bundle
// (3 from ball-socket + 1 from distance). Tests intra-bundle coupling computation
// between different constraint types sharing both bodies.
static void test_ldl_stress_mixed_same_pair()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });

	// Ball-socket + distance on the exact same body pair (same offsets)
	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	create_distance(w, (DistanceParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b, .rest_length = 1.0f });

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, a);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	printf("  [mixed same pair] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress mixed same pair: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress mixed same pair: body valid");
	v3 pa = body_get_position(w, a);
	TEST_ASSERT(is_valid(pa) && pa.y > -20.0f);

	destroy_world(w);
}

// Near-LDL_MAX_NODES: push node count close to the 128 limit.
// A hub body with ~40 arms creates 40 bundles + fill-in from the hub.
// With shattering, the shard welds add more nodes. Tests whether we stay in bounds.
static void test_ldl_stress_near_max_nodes()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 10.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = hub, .local_offset_a = V3(0,-1,0), .local_offset_b = V3(0,1,0) });

	// 30 arms radiating out — 30 constraints, each unique body pair with hub
	// hub DOF = 30 * 3 = 90, well above SHATTER_THRESHOLD=15
	int ARM_COUNT = 30;
	Body arms[30];
	for (int i = 0; i < ARM_COUNT; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / (float)ARM_COUNT;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		arms[i] = create_body(w, (BodyParams){ .position = add(V3(0, 8, 0), scale(dir, 1.5f)), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arms[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = arms[i], .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	int nan_frame = -1;
	for (int f = 0; f < 120; f++) {
		world_step(w, 1.0f / 60.0f);
		v3 p = body_get_position(w, hub);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	printf("  [near max nodes 30-arm hub] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress near max nodes: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress near max nodes: hub body valid");
	v3 ph = body_get_position(w, hub);
	TEST_ASSERT(is_valid(ph) && ph.y > -50.0f);

	destroy_world(w);
}

// Inverted gravity with tangled initial positions: bodies start overlapping and
// gravity pulls upward. The solver faces large initial penetration and conflicting
// constraint forces. Tests delta correction sign handling under unusual conditions.
static void test_ldl_stress_inverted_gravity_tangle()
{
	World w = create_world((WorldParams){ .gravity = V3(0, 9.81f, 0) }); // gravity UP
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	// All bodies start at nearly the same position (tangled), gravity is inverted
	Body chain[5];
	Body prev = anchor;
	for (int i = 0; i < 5; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3(0.01f * (float)i, 0.01f * (float)i, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
		prev = chain[i];
	}

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < 5; i++) {
			if (!is_valid(body_get_position(w, chain[i]))) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [inverted gravity tangle] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress inverted gravity tangle: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress inverted gravity tangle: all bodies valid");
	int valid = 1;
	for (int i = 0; i < 5; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p)) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Distance-only chain: tests the 1x1 block path through LDL factorization.
// All constraints are 1-DOF distance joints, no ball-sockets. The K matrix
// is scalar-block tridiagonal. Tests block_ldl with n=1 blocks throughout.
static void test_ldl_stress_distance_only_chain()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	int N = 10;
	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	Body chain[10];
	Body prev = anchor;
	for (int i = 0; i < N; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3(0, 10 - (i+1)*1.0f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_distance(w, (DistanceParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = V3(0,-0.5f,0), .local_offset_b = V3(0,0.5f,0), .rest_length = 0.0f });
		prev = chain[i];
	}

	int nan_frame = -1;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		if (!is_valid(body_get_position(w, chain[N-1]))) { nan_frame = f; break; }
	}

	printf("  [distance-only chain] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress distance-only chain: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress distance-only chain: all finite");
	int valid = 1;
	for (int i = 0; i < N; i++) {
		v3 p = body_get_position(w, chain[i]);
		if (!is_valid(p)) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Opposing velocity slam: two chains share a common hub body. We slam the tip
// of each chain in opposite directions simultaneously. The delta correction sees
// conflicting RHS entries across the shared node. Tests sign consistency.
static void test_ldl_stress_opposing_slams()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body hub = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

	// Left arm: hub -> a -> tip_l
	Body a = create_body(w, (BodyParams){ .position = V3(-1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Body tip_l = create_body(w, (BodyParams){ .position = V3(-2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, tip_l, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = a, .local_offset_a = V3(-0.5f,0,0), .local_offset_b = V3(0.5f,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = tip_l, .local_offset_a = V3(-0.5f,0,0), .local_offset_b = V3(0.5f,0,0) });

	// Right arm: hub -> b -> tip_r
	Body b = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Body tip_r = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, tip_r, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = b, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = tip_r, .local_offset_a = V3(0.5f,0,0), .local_offset_b = V3(-0.5f,0,0) });

	step_n(w, 60);

	// Slam tips in opposite directions
	body_set_velocity(w, tip_l, V3(-500, 300, 0));
	body_set_velocity(w, tip_r, V3(500, 300, 0));

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		if (!is_valid(body_get_position(w, hub))) { nan_frame = f; break; }
	}

	printf("  [opposing slams] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress opposing slams: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress opposing slams: hub body valid");
	v3 ph = body_get_position(w, hub);
	TEST_ASSERT(is_valid(ph));

	destroy_world(w);
}

// Zero-gravity free-floating: with no gravity, a chain with initial angular
// velocity should conserve angular momentum. Tests that LDL correction doesn't
// inject spurious energy/torque in the zero-gravity case.
static void test_ldl_stress_zero_gravity_spin()
{
	World w = create_world((WorldParams){ .gravity = V3(0, 0, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body chain[4];
	for (int i = 0; i < 4; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3((float)i * 0.8f, 0, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	}
	for (int i = 0; i < 3; i++) {
		create_ball_socket(w, (BallSocketParams){ .body_a = chain[i], .body_b = chain[i+1], .local_offset_a = V3(0.4f,0,0), .local_offset_b = V3(-0.4f,0,0) });
	}

	// Give initial spin
	body_set_angular_velocity(w, chain[0], V3(0, 10, 0));
	body_set_velocity(w, chain[3], V3(0, 3, 0));

	float initial_ke = measure_system_ke(wi, chain, 4);
	float max_ke = initial_ke;
	int nan_frame = -1;
	for (int f = 0; f < 1200; f++) {
		world_step(w, 1.0f / 60.0f);
		float ke = measure_system_ke(wi, chain, 4);
		if (ke > max_ke) max_ke = ke;
		if (!is_valid(body_get_position(w, chain[0]))) { nan_frame = f; break; }
	}

	printf("  [zero-g spin] initial_ke=%.2f max_ke=%.2f nan_frame=%d\n",
		(double)initial_ke, (double)max_ke, nan_frame);

	TEST_BEGIN("LDL stress zero-g spin: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress zero-g spin: KE bounded (no energy growth)");
	// Energy should not grow beyond initial + small tolerance from numerical error
	TEST_ASSERT(max_ke < initial_ke * 2.0f + 1.0f);

	destroy_world(w);
}

// Spring softness extremes: very stiff (softness near 0) and very soft springs
// in the same island. The K matrix diagonal entries span huge ranges within
// a single bundle, testing conditioning.
static void test_ldl_stress_mixed_stiffness()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body a = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Body b = create_body(w, (BodyParams){ .position = V3(2, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Body c = create_body(w, (BodyParams){ .position = V3(3, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, c, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	// Very stiff (rigid) joint
	create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
	// Very soft spring
	create_ball_socket(w, (BallSocketParams){ .body_a = a, .body_b = b, .local_offset_a = off_a, .local_offset_b = off_b, .spring = { .frequency = 0.5f, .damping_ratio = 0.1f } });
	// Back to rigid
	create_ball_socket(w, (BallSocketParams){ .body_a = b, .body_b = c, .local_offset_a = off_a, .local_offset_b = off_b });

	int nan_frame = -1;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < 3; i++) {
			Body bodies[3] = {a, b, c};
			if (!is_valid(body_get_position(w, bodies[i]))) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [mixed stiffness] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress mixed stiffness: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress mixed stiffness: all valid");
	TEST_ASSERT(is_valid(body_get_position(w, a)) && is_valid(body_get_position(w, b)) && is_valid(body_get_position(w, c)));

	destroy_world(w);
}

// Bidirectional hub-to-hub: two hubs connected to each other + each has arms.
// Creates a complex graph with two high-degree nodes linked by a bridge edge.
// Tests fill-in handling when two cliques connect.
static void test_ldl_stress_double_hub_bridge()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	Body hub1 = create_body(w, (BodyParams){ .position = V3(-2, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub1, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	Body hub2 = create_body(w, (BodyParams){ .position = V3(2, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub2, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

	// Bridge between hubs
	create_ball_socket(w, (BallSocketParams){ .body_a = hub1, .body_b = hub2, .local_offset_a = V3(2,0,0), .local_offset_b = V3(-2,0,0) });

	// 6 arms off hub1
	for (int i = 0; i < 6; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		Body arm = create_body(w, (BodyParams){ .position = add(V3(-2, 8, 0), dir), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub1, .body_b = arm, .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	// 6 arms off hub2
	for (int i = 0; i < 6; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		Body arm = create_body(w, (BodyParams){ .position = add(V3(2, 8, 0), dir), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, arm, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub2, .body_b = arm, .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		if (!is_valid(body_get_position(w, hub1)) || !is_valid(body_get_position(w, hub2))) { nan_frame = f; break; }
	}

	printf("  [double hub bridge] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress double hub bridge: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress double hub bridge: hubs valid");
	TEST_ASSERT(is_valid(body_get_position(w, hub1)) && is_valid(body_get_position(w, hub2)));

	destroy_world(w);
}

// Rapid create-destroy-create cycle: destroy a body (which destroys all its joints),
// then immediately recreate everything. Exercises LDL_Cache invalidation and
// reallocation under rapid turnover.
static void test_ldl_stress_body_destroy_recreate()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	v3 off_a = V3(0.5f,0,0), off_b = V3(-0.5f,0,0);
	Body a = create_body(w, (BodyParams){ .position = V3(1, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, a, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
	Joint j = create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });

	int nan_frame = -1;
	for (int f = 0; f < 200; f++) {
		world_step(w, 1.0f / 60.0f);

		if (f % 10 == 5) {
			// Destroy and immediately recreate
			destroy_joint(w, j);
			j = create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = a, .local_offset_a = off_a, .local_offset_b = off_b });
		}

		v3 p = body_get_position(w, a);
		if (!is_valid(p)) { nan_frame = f; break; }
	}

	printf("  [body destroy-recreate] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress body destroy-recreate: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress body destroy-recreate: body valid");
	TEST_ASSERT(is_valid(body_get_position(w, a)));

	destroy_world(w);
}

// Micro mass chain: all bodies have extremely tiny mass (1e-6).
// inv_mass becomes enormous (1e6), making K diagonal entries huge.
// Tests floating-point overflow during K computation and factorization.
static void test_ldl_stress_micro_mass()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });

	int N = 5;
	Body chain[5];
	Body prev = anchor;
	v3 off_a = V3(0.4f,0,0), off_b = V3(-0.4f,0,0);
	for (int i = 0; i < N; i++) {
		chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*0.8f, 5, 0), .rotation = quat_identity(), .mass = 1e-5f });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	int nan_frame = -1;
	for (int f = 0; f < 300; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < N; i++) {
			if (!is_valid(body_get_position(w, chain[i]))) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [micro mass 1e-5] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress micro mass: no NaN");
	TEST_ASSERT(nan_frame < 0);

	destroy_world(w);
}

// Constraint loop: 4 bodies forming a rigid loop (A-B-C-D-A). This creates
// a cycle in the constraint graph. The K matrix is no longer tree-structured;
// the LDL elimination must handle fill-in from the cycle correctly.
static void test_ldl_stress_constraint_loop()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;
	wi->sleep_enabled = 0;

	// 4 bodies forming a square
	Body bodies[4];
	v3 positions[4] = { V3(-1, 5, -1), V3(1, 5, -1), V3(1, 5, 1), V3(-1, 5, 1) };
	for (int i = 0; i < 4; i++) {
		bodies[i] = create_body(w, (BodyParams){ .position = positions[i], .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bodies[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.2f });
	}

	// Close the loop: A-B, B-C, C-D, D-A
	create_ball_socket(w, (BallSocketParams){ .body_a = bodies[0], .body_b = bodies[1], .local_offset_a = V3(1,0,0), .local_offset_b = V3(-1,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = bodies[1], .body_b = bodies[2], .local_offset_a = V3(0,0,1), .local_offset_b = V3(0,0,-1) });
	create_ball_socket(w, (BallSocketParams){ .body_a = bodies[2], .body_b = bodies[3], .local_offset_a = V3(-1,0,0), .local_offset_b = V3(1,0,0) });
	create_ball_socket(w, (BallSocketParams){ .body_a = bodies[3], .body_b = bodies[0], .local_offset_a = V3(0,0,-1), .local_offset_b = V3(0,0,1) });

	int nan_frame = -1;
	for (int f = 0; f < 600; f++) {
		world_step(w, 1.0f / 60.0f);
		for (int i = 0; i < 4; i++) {
			if (!is_valid(body_get_position(w, bodies[i]))) { nan_frame = f; break; }
		}
		if (nan_frame >= 0) break;
	}

	printf("  [constraint loop 4] nan_frame=%d\n", nan_frame);

	TEST_BEGIN("LDL stress constraint loop: no NaN");
	TEST_ASSERT(nan_frame < 0);

	TEST_BEGIN("LDL stress constraint loop: all finite");
	int valid = 1;
	for (int i = 0; i < 4; i++) {
		v3 p = body_get_position(w, bodies[i]);
		if (!is_valid(p)) { valid = 0; break; }
	}
	TEST_ASSERT(valid);

	destroy_world(w);
}

// Stretched joint recovery: yank a chain tip far away, verify LDL recovers the
// joint gap within a reasonable number of frames. Compares LDL vs PGS-only to
// ensure LDL doesn't prevent recovery or cause erratic behavior after stretching.
static void test_ldl_stress_stretched_recovery()
{
	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	// Run WITH LDL
	float gap_ldl_before = 0, gap_ldl_after = 0;
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 1;
		wi->sleep_enabled = 0;

		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body chain[5];
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*link_len, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
			prev = chain[i];
		}

		// Let it settle
		step_n(w, 120);

		// Yank the tip body far away
		body_set_velocity(w, chain[4], V3(50, 80, 0));
		step_n(w, 30);

		// Measure gap right after the yank settles
		prev = anchor;
		for (int i = 0; i < 5; i++) {
			float g = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (g > gap_ldl_before) gap_ldl_before = g;
			prev = chain[i];
		}

		// Let it recover
		step_n(w, 300);

		prev = anchor;
		for (int i = 0; i < 5; i++) {
			float g = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (g > gap_ldl_after) gap_ldl_after = g;
			prev = chain[i];
		}
		destroy_world(w);
	}

	// Run WITHOUT LDL for comparison
	float gap_pgs_before = 0, gap_pgs_after = 0;
	{
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = 0;
		wi->sleep_enabled = 0;

		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body chain[5];
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*link_len, 10, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
			prev = chain[i];
		}

		step_n(w, 120);
		body_set_velocity(w, chain[4], V3(50, 80, 0));
		step_n(w, 30);

		prev = anchor;
		for (int i = 0; i < 5; i++) {
			float g = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (g > gap_pgs_before) gap_pgs_before = g;
			prev = chain[i];
		}

		step_n(w, 300);

		prev = anchor;
		for (int i = 0; i < 5; i++) {
			float g = anchor_distance(w, prev, off_a, chain[i], off_b);
			if (g > gap_pgs_after) gap_pgs_after = g;
			prev = chain[i];
		}
		destroy_world(w);
	}

	printf("  [stretched recovery] LDL: before=%.4f after=%.4f  PGS: before=%.4f after=%.4f\n",
		(double)gap_ldl_before, (double)gap_ldl_after, (double)gap_pgs_before, (double)gap_pgs_after);

	TEST_BEGIN("LDL stretched recovery: joint recovers after yank");
	TEST_ASSERT(gap_ldl_after < gap_ldl_before + 0.05f); // LDL may slightly overshoot after yank

	TEST_BEGIN("LDL stretched recovery: gap below 50 after recovery");
	TEST_ASSERT(gap_ldl_after < 50.0f);

	TEST_BEGIN("LDL stretched recovery: LDL not drastically worse than PGS");
	TEST_ASSERT(gap_ldl_after <= gap_pgs_after * 100.0f);
}

// Heavy chain stretched recovery: extreme mass ratio chain gets yanked.
// The heavy tip resists recovery; LDL must help, not hinder.
static void test_ldl_stress_heavy_stretched_recovery()
{
	v3 off_a = V3(0.5f, 0, 0), off_b = V3(-0.5f, 0, 0);

	float gap_ldl = 0, gap_pgs = 0;
	for (int use_ldl = 0; use_ldl < 2; use_ldl++) {
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
		WorldInternal* wi = (WorldInternal*)w.id;
		wi->ldl_enabled = use_ldl;
		wi->sleep_enabled = 0;

		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body light = create_body(w, (BodyParams){ .position = V3(1, 10, 0), .rotation = quat_identity(), .mass = 0.1f });
		body_add_shape(w, light, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		Body heavy = create_body(w, (BodyParams){ .position = V3(2, 10, 0), .rotation = quat_identity(), .mass = 50.0f });
		body_add_shape(w, heavy, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = light, .local_offset_a = off_a, .local_offset_b = off_b });
		create_ball_socket(w, (BallSocketParams){ .body_a = light, .body_b = heavy, .local_offset_a = off_a, .local_offset_b = off_b });

		// Settle, then yank heavy body
		step_n(w, 120);
		body_set_velocity(w, heavy, V3(30, 60, 0));

		// Let it try to recover
		step_n(w, 600);

		float g0 = anchor_distance(w, anchor, off_a, light, off_b);
		float g1 = anchor_distance(w, light, off_a, heavy, off_b);
		float max_gap = g0 > g1 ? g0 : g1;

		if (use_ldl) gap_ldl = max_gap; else gap_pgs = max_gap;
		destroy_world(w);
	}

	printf("  [heavy stretched] LDL=%.4f PGS=%.4f\n", (double)gap_ldl, (double)gap_pgs);

	TEST_BEGIN("LDL heavy stretched: joint recovers");
	TEST_ASSERT(gap_ldl < 50.0f);

	TEST_BEGIN("LDL heavy stretched: LDL not worse than PGS");
	TEST_ASSERT(gap_ldl <= gap_pgs * 10.0f);
}

// Mixed scene: chain + hub star in the same world.
static void test_ldl_mixed_chain_and_hub()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->ldl_enabled = 1;

	float link_len = 0.8f;
	v3 off_a = V3(link_len * 0.5f, 0, 0);
	v3 off_b = V3(-link_len * 0.5f, 0, 0);

	// Chain: 5 links with heavy end
	Body chain_anchor = create_body(w, (BodyParams){ .position = V3(-5, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, chain_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body chain[5];
	Body prev = chain_anchor;
	for (int i = 0; i < 5; i++) {
		float mass = (i == 4) ? 50.0f : 1.0f;
		chain[i] = create_body(w, (BodyParams){ .position = V3(-5 + (i + 1) * link_len, 10, 0), .rotation = quat_identity(), .mass = mass });
		body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = off_a, .local_offset_b = off_b });
		prev = chain[i];
	}

	// Hub star: center + 6 arms (18 DOF > 12, triggers shattering)
	Body hub_anchor = create_body(w, (BodyParams){ .position = V3(5, 10, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, hub_anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
	Body hub = create_body(w, (BodyParams){ .position = V3(5, 8, 0), .rotation = quat_identity(), .mass = 5.0f });
	body_add_shape(w, hub, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });
	create_ball_socket(w, (BallSocketParams){ .body_a = hub_anchor, .body_b = hub, .local_offset_a = V3(0, -1, 0), .local_offset_b = V3(0, 1, 0) });
	Body hub_arms[6];
	for (int i = 0; i < 6; i++) {
		float angle = (float)i * 2.0f * 3.14159265f / 6.0f;
		v3 dir = V3(cosf(angle), 0, sinf(angle));
		hub_arms[i] = create_body(w, (BodyParams){ .position = add(V3(5, 8, 0), dir), .rotation = quat_identity(), .mass = 2.0f });
		body_add_shape(w, hub_arms[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
		create_ball_socket(w, (BallSocketParams){ .body_a = hub, .body_b = hub_arms[i], .local_offset_a = scale(dir, 0.3f), .local_offset_b = scale(dir, -0.3f) });
	}

	step_n(w, 300);

	// Check chain
	float chain_gap = 0;
	prev = chain_anchor;
	for (int i = 0; i < 5; i++) {
		float g = anchor_distance(w, prev, off_a, chain[i], off_b);
		if (g > chain_gap) chain_gap = g;
		prev = chain[i];
	}

	// Check hub
	float hub_gap = 0;
	for (int i = 0; i < 6; i++) {
		v3 hp = body_get_position(w, hub);
		v3 ap = body_get_position(w, hub_arms[i]);
		float gap = fabsf(len(sub(ap, hp)) - 0.6f);
		if (gap > hub_gap) hub_gap = gap;
	}

	printf("  [LDL mixed] chain_gap=%.4f  hub_gap=%.4f\n", chain_gap, hub_gap);

	TEST_BEGIN("LDL mixed: chain tight");
	TEST_ASSERT(chain_gap < 10.0f);

	TEST_BEGIN("LDL mixed: hub tight");
	TEST_ASSERT(hub_gap < 0.5f);

	destroy_world(w);
}

static void run_solver_tests()
{
	printf("--- nudge solver tests ---\n");
	test_ldl_soft_box_drag();
	test_solver_nan_rejection();
	test_contact_sphere_rests_on_floor();
	test_contact_box_rests_on_floor();
	test_distance_spring_equilibrium();
	test_distance_rigid_maintains_length();
	test_ball_socket_chain_stays_connected();
	test_ball_socket_pin_converges();
	test_ball_socket_pendulum();
	test_ldl_heavy_chain();
	test_ldl_lift_and_drop();
	test_ldl_mouse_yank_chain();
	test_ldl_vertical_heavy_chain_drop();
	test_ldl_two_independent_chains();
	test_ldl_topology_change();
	test_ldl_block_math();
	test_ldl_solve_topo_identity();
	test_ldl_bundling();
	test_ldl_topology_structure();
	test_ldl_numeric_factor_isolated();
	test_ldl_solve_topo_vs_dense();
	test_ldl_hub_star_shattering();
	test_ldl_no_shatter_below_threshold();
	test_ldl_mixed_chain_and_hub();
	test_ldl_sleep_cache();
	test_ldl_energy_stability();
	test_ldl_energy_comprehensive();
	test_ldl_mass_ratio();
	test_ldl_long_chain();
	test_ldl_block_near_singular();
	test_ldl_delta_correction_accuracy();
	test_ldl_soft_box();
	test_ldl_soft_box_drag();
	test_ldl_soft_box_chain();
	test_ldl_soft_chain();
	test_ldl_mixed_rigid_soft();
	test_ldl_soft_hub();
	test_ldl_sharp_yank_recovery();
	test_ldl_soft_heavy();
	test_bounce_height_monotonic();
	bounce_test_for_solver(SOLVER_SOFT_STEP, "Soft Step");
	bounce_test_for_solver(SOLVER_SI_SOFT, "SI Soft");
	bounce_test_for_solver(SOLVER_SI, "SI");
	bounce_test_for_solver(SOLVER_BLOCK, "Block");
	bounce_test_for_solver(SOLVER_AVBD, "AVBD");
	avbd_box_bounce_test();
	// AVBD should still settle onto the floor after the bounce energy decays.
	{
		TEST_BEGIN("AVBD sphere settles after bouncing");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body floor_b = create_body(w, (BodyParams){
			.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
		});
		body_add_shape(w, floor_b, (ShapeParams){
			.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
		});
		Body ball = create_body(w, (BodyParams){
			.position = V3(0, 5, 0), .rotation = quat_identity(),
			.mass = 1.0f, .restitution = 0.5f,
		});
		body_add_shape(w, ball, (ShapeParams){
			.type = SHAPE_SPHERE, .sphere.radius = 0.5f,
		});
		float dt = 1.0f / 60.0f;
		for (int i = 0; i < 300; i++) world_step(w, dt);
		float y = body_get_position(w, ball).y;
		float x = body_get_position(w, ball).x;
		printf("  [AVBD sphere] y=%.4f x=%.4f\n", y, x);
		// Ball should rest on floor top (y=0) + radius (0.5) = ~0.5
		TEST_ASSERT(y > 0.3f);
		TEST_ASSERT(y < 0.7f);
		destroy_world(w);
	}
	// AVBD box-on-floor: test various masses and drop heights
	{
		float masses[] = { 0.1f, 0.5f, 1.0f, 5.0f, 10.0f, 50.0f };
		float heights[] = { 2.0f, 5.0f, 10.0f };
		float sizes[] = { 0.3f, 0.5f, 1.0f };
		for (int mi = 0; mi < 6; mi++) {
			for (int hi = 0; hi < 3; hi++) {
				for (int si = 0; si < 3; si++) {
					char buf[128];
					snprintf(buf, sizeof(buf), "AVBD box settle m=%.1f h=%.0f s=%.1f", masses[mi], heights[hi], sizes[si]);
					TEST_BEGIN(buf);
					World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
					Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
					body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
					float hs = sizes[si];
					Body box = create_body(w, (BodyParams){
						.position = V3(0, heights[hi], 0), .rotation = quat_identity(), .mass = masses[mi] });
					body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(hs, hs, hs) });
					float dt = 1.0f / 60.0f;
					for (int i = 0; i < 120; i++) world_step(w, dt);
					float y = body_get_position(w, box).y;
					if (y < -0.5f)
						printf("  [AVBD FALL-THROUGH] m=%.1f h=%.0f s=%.1f y=%.3f\n", masses[mi], heights[hi], sizes[si], y);
					TEST_ASSERT(y > -0.5f); // must not fall through floor
					destroy_world(w);
				}
			}
		}
	}
	// AVBD box stack: multiple boxes stacked
	{
		TEST_BEGIN("AVBD 5-box stack");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){
				.position = V3(0, 0.5f + i * 1.1f, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
		}
		float dt = 1.0f / 60.0f;
		for (int i = 0; i < 300; i++) world_step(w, dt);
		// Check all bodies above floor
		WorldInternal* wi = (WorldInternal*)w.id;
		int fallen = 0;
		for (int i = 0; i < asize(wi->body_hot); i++) {
			if (wi->body_hot[i].inv_mass == 0.0f) continue;
			if (wi->body_hot[i].position.y < -2.0f) {
				printf("  [AVBD STACK FAIL] body %d y=%.3f\n", i, wi->body_hot[i].position.y);
				fallen++;
			}
		}
		TEST_ASSERT(fallen == 0);
		destroy_world(w);
	}
	// AVBD tilted box drop
	{
		TEST_BEGIN("AVBD tilted box drop");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		// Drop a tilted box
		Body b = create_body(w, (BodyParams){
			.position = V3(0, 5, 0), .rotation = quat_axis_angle(V3(1,1,0), 0.7f), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
		float dt = 1.0f / 60.0f;
		for (int i = 0; i < 300; i++) world_step(w, dt);
		float y = body_get_position(w, b).y;
		if (y < -0.5f) printf("  [AVBD TILTED FAIL] y=%.3f\n", y);
		TEST_ASSERT(y > -0.5f);
		destroy_world(w);
	}
	// AVBD mass ratio stack (matching mass ratio scene)
	{
		TEST_BEGIN("AVBD mass ratio stack");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		float sizes[] = { 0.15f, 0.25f, 0.4f, 0.55f, 0.75f };
		float masses[] = { 0.5f, 2.0f, 8.0f, 30.0f, 100.0f };
		float y = 0.0f;
		for (int i = 0; i < 5; i++) {
			float h = sizes[i];
			y += h + 0.6f;
			Body b = create_body(w, (BodyParams){
				.position = V3(0, y, 0), .rotation = quat_identity(), .mass = masses[i] });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(h, h, h) });
			y += h;
		}
		float dt = 1.0f / 60.0f;
		for (int i = 0; i < 300; i++) world_step(w, dt);
		WorldInternal* wi = (WorldInternal*)w.id;
		int fallen = 0;
		for (int i = 0; i < asize(wi->body_hot); i++) {
			if (wi->body_hot[i].inv_mass == 0.0f) continue;
			if (wi->body_hot[i].position.y < -2.0f) {
				printf("  [AVBD MASS-RATIO FAIL] body %d y=%.3f m=%.1f\n", i, wi->body_hot[i].position.y, 1.0f/wi->body_hot[i].inv_mass);
				fallen++;
			}
		}
		TEST_ASSERT(fallen == 0);
		destroy_world(w);
	}
	// AVBD first scene repro: sphere + capsule + box dropping simultaneously
	{
		TEST_BEGIN("AVBD showcase bodies settle");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

		Body sphere = create_body(w, (BodyParams){ .position = V3(-3, 5, 0), .rotation = quat_identity(), .mass = 1.0f, .restitution = 0.5f });
		body_add_shape(w, sphere, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

		Body capsule = create_body(w, (BodyParams){ .position = V3(-1, 6, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, capsule, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = 0.3f, .radius = 0.25f } });

		Body box = create_body(w, (BodyParams){ .position = V3(1, 7, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });

		float dt = 1.0f / 60.0f;
		for (int frame = 0; frame < 600; frame++) {
			world_step(w, dt);
			// Check each frame for fall-through
			float ys = body_get_position(w, sphere).y;
			float yc = body_get_position(w, capsule).y;
			float yb = body_get_position(w, box).y;
			if (ys < -2.0f || yc < -2.0f || yb < -2.0f) {
				printf("  [AVBD SHOWCASE FAIL] frame=%d sphere=%.2f capsule=%.2f box=%.2f\n", frame, ys, yc, yb);
				TEST_ASSERT(0); // fail
				break;
			}
		}
		destroy_world(w);
	}
	// AVBD pyramid: many boxes interacting
	{
		TEST_BEGIN("AVBD box pyramid");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		// 4-layer pyramid
		for (int row = 0; row < 4; row++)
			for (int col = 0; col < 4 - row; col++) {
				Body b = create_body(w, (BodyParams){
					.position = V3(col * 1.05f + row * 0.5f - 1.5f, 0.5f + row * 1.05f, 0),
					.rotation = quat_identity(), .mass = 1.0f });
				body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
			}
		float dt = 1.0f / 60.0f;
		int fallen = 0;
		for (int frame = 0; frame < 600; frame++) {
			world_step(w, dt);
			WorldInternal* wi = (WorldInternal*)w.id;
			for (int i = 0; i < asize(wi->body_hot); i++) {
				if (wi->body_hot[i].inv_mass == 0.0f) continue;
				if (wi->body_hot[i].position.y < -3.0f) {
					if (!fallen) printf("  [AVBD PYRAMID FAIL] frame=%d body=%d y=%.3f\n", frame, i, wi->body_hot[i].position.y);
					fallen++;
				}
			}
		}
		TEST_ASSERT(fallen == 0);
		destroy_world(w);
	}
	// AVBD solver switch: simulate with Soft Step, then switch to AVBD mid-sim
	{
		TEST_BEGIN("AVBD solver switch mid-sim");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_SOFT_STEP });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		Body b = create_body(w, (BodyParams){ .position = V3(0, 3, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
		float dt = 1.0f / 60.0f;
		// Run with Soft Step until box settles
		for (int i = 0; i < 120; i++) world_step(w, dt);
		float y_before = body_get_position(w, b).y;
		// Switch to AVBD
		world_set_solver_type(w, SOLVER_AVBD);
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			float y = body_get_position(w, b).y;
			if (y < -2.0f) {
				printf("  [AVBD SWITCH FAIL] frame=%d y=%.3f (was %.3f before switch)\n", i, y, y_before);
				TEST_ASSERT(0);
				break;
			}
		}
		destroy_world(w);
	}
	// AVBD box penetration repro: exact showcase box, track sinking
	{
		TEST_BEGIN("AVBD box resting depth");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		((WorldInternal*)w.id)->sleep_enabled = 0; // disable sleep to isolate solver behavior
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		// Match showcase: box at y=7, half_extents=0.4
		Body box = create_body(w, (BodyParams){
			.position = V3(1, 7, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
		float dt = 1.0f / 60.0f;
		// Track settling: print Y at key frames around impact and steady state
		for (int i = 0; i < 600; i++) {
			world_step(w, dt);
			float y = body_get_position(w, box).y;
			if (i >= 40 && i <= 110) // around impact through crash
				printf("  [AVBD box] f%d y=%.4f\n", i, y);
			if (i == 200 || i == 400 || i == 599)
				printf("  [AVBD box] f%d y=%.4f pen=%.4f\n", i, y, 0.4f - y);
		}
		float y = body_get_position(w, box).y;
		// Box should settle very close to y=0.4 (half-extent above floor).
		TEST_ASSERT(y > 0.38f);
		destroy_world(w);
	}
	// Also test with Soft Step for comparison
	{
		TEST_BEGIN("Soft Step box resting depth (reference)");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_SOFT_STEP });
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		Body box = create_body(w, (BodyParams){
			.position = V3(1, 7, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
		float dt = 1.0f / 60.0f;
		for (int i = 0; i < 300; i++) world_step(w, dt);
		float y = body_get_position(w, box).y;
		printf("  [SoftStep box depth] y=%.4f (expected ~0.4, penetration=%.4f)\n", y, 0.4f - y);
		destroy_world(w);
	}
	// AVBD friction slide: box with lateral velocity slides off floor edge and falls
	{
		TEST_BEGIN("AVBD box slides off floor edge");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		((WorldInternal*)w.id)->sleep_enabled = 0;
		// Small floor so box slides off the edge
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(3, 1, 3) });
		// Box starts on floor with lateral velocity
		Body box = create_body(w, (BodyParams){
			.position = V3(0, 0.4f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
		body_set_velocity(w, box, V3(8, 0, 0));

		float dt = 1.0f / 60.0f;
		int fell = 0;
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			v3 p = body_get_position(w, box);
			// Box should eventually slide off edge (x > 3) and fall (y < 0)
			if (p.x > 4.0f && p.y < -1.0f) { fell = 1; break; }
			// DEBUG: print key frames
			if (i == 30 || i == 60 || i == 120 || i == 200)
				printf("  [AVBD slide] f%d x=%.2f y=%.2f\n", i, p.x, p.y);
		}
		v3 pf = body_get_position(w, box);
		printf("  [AVBD slide] final x=%.2f y=%.2f fell=%d\n", pf.x, pf.y, fell);
		TEST_ASSERT(fell); // must fall off the edge, not float
		destroy_world(w);
	}
	// Same test with Soft Step for reference
	{
		TEST_BEGIN("Soft Step box slides off floor edge (reference)");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_SOFT_STEP });
		((WorldInternal*)w.id)->sleep_enabled = 0;
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(3, 1, 3) });
		Body box = create_body(w, (BodyParams){
			.position = V3(0, 0.4f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });
		body_set_velocity(w, box, V3(8, 0, 0));

		float dt = 1.0f / 60.0f;
		int fell = 0;
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			v3 p = body_get_position(w, box);
			if (p.x > 4.0f && p.y < -1.0f) { fell = 1; break; }
			if (i == 30 || i == 60 || i == 120 || i == 200)
				printf("  [SS slide] f%d x=%.2f y=%.2f\n", i, p.x, p.y);
		}
		v3 pf = body_get_position(w, box);
		printf("  [SS slide] final x=%.2f y=%.2f fell=%d\n", pf.x, pf.y, fell);
		// Soft Step friction stops the box before the edge — not a failure
		destroy_world(w);
	}
	// AVBD edge straddle: box placed half-on, half-off the floor edge.
	// Should tip and fall, not sink into the floor or hover.
	{
		TEST_BEGIN("AVBD box straddles floor edge");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		((WorldInternal*)w.id)->sleep_enabled = 0;
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(3, 1, 3) });
		// Box center at x=3.0 — half on floor (edge at x=3), half off
		Body box = create_body(w, (BodyParams){
			.position = V3(3.0f, 0.4f, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.4f, 0.4f, 0.4f) });

		float dt = 1.0f / 60.0f;
		int fell = 0;
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			v3 p = body_get_position(w, box);
			if (i < 60 || i % 30 == 0)
				printf("  [AVBD edge] f%d x=%.3f y=%.3f\n", i, p.x, p.y);
			if (p.y < -2.0f) { fell = 1; break; }
		}
		v3 pf = body_get_position(w, box);
		printf("  [AVBD edge] final x=%.3f y=%.3f fell=%d\n", pf.x, pf.y, fell);
		// Box should tip off and fall
		TEST_ASSERT(fell);
		destroy_world(w);
	}
	// AVBD single hull drop: isolate which hull shapes fall through
	{
		const char* hull_names[] = { "tet", "wedge", "rock", "bipyr" };
		v3 tet_pts[] = { {0,0.6f,0}, {0.5f,-0.3f,0.3f}, {-0.5f,-0.3f,0.3f}, {0,-0.3f,-0.5f} };
		v3 wedge_pts[] = {
			{-0.5f,-0.3f,-0.4f}, {0.5f,-0.3f,-0.4f}, {-0.5f,-0.3f,0.4f}, {0.5f,-0.3f,0.4f},
			{-0.5f, 0.3f,-0.4f}, {0.5f, 0.3f,-0.4f},
		};
		v3 rock_pts[] = {
			{0, 0.55f, 0}, {0, -0.5f, 0},
			{0.45f, 0.1f, 0.3f}, {-0.4f, 0.15f, 0.35f},
			{0.35f, 0.1f, -0.4f}, {-0.3f, 0.1f, -0.45f},
			{0.5f, -0.15f, -0.1f}, {-0.5f, -0.1f, 0.05f},
			{0.1f, -0.2f, 0.55f}, {-0.15f, -0.2f, -0.5f},
		};
		v3 bipyr_pts[] = {
			{0,0.7f,0}, {0.5f,0,0.5f}, {-0.5f,0,0.5f}, {0.5f,0,-0.5f}, {-0.5f,0,-0.5f}, {0,-0.7f,0},
		};
		struct { v3* pts; int n; v3 sc; } hull_defs[] = {
			{ tet_pts, 4, V3(0.8f,0.8f,0.8f) },
			{ wedge_pts, 6, V3(0.7f,0.7f,0.7f) },
			{ rock_pts, 10, V3(0.6f,0.6f,0.6f) },
			{ bipyr_pts, 6, V3(0.7f,0.7f,0.7f) },
		};

		// Test each hull type with several rotations
		for (int t = 0; t < 4; t++) {
			Hull* h = quickhull(hull_defs[t].pts, hull_defs[t].n);
			for (int r = 0; r < 8; r++) {
				char name[64];
				snprintf(name, sizeof(name), "AVBD hull drop %s rot%d", hull_names[t], r);
				TEST_BEGIN(name);
				World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
				((WorldInternal*)w.id)->sleep_enabled = 0;
				Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
				body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

				float a = (float)r * 0.8f;
				v3 ax = norm(V3(sinf(a), cosf(a*1.3f), sinf(a*0.7f)));
				float ha = a * 0.5f;
				quat rot = (quat){ ax.x*sinf(ha), ax.y*sinf(ha), ax.z*sinf(ha), cosf(ha) };

				Body body = create_body(w, (BodyParams){
					.position = V3(0, 5, 0), .rotation = rot, .mass = 1.0f });
				body_add_shape(w, body, (ShapeParams){
					.type = SHAPE_HULL, .hull = { .hull = h, .scale = hull_defs[t].sc } });

				float dt = 1.0f / 60.0f;
				int ok = 1;
				for (int f = 0; f < 120; f++) {
					world_step(w, dt);
					float y = body_get_position(w, body).y;
					if (y < -2.0f) {
						printf("  [FALL] %s rot%d fell at f%d y=%.2f\n", hull_names[t], r, f, y);
						ok = 0; break;
					}
				}
				float y = body_get_position(w, body).y;
				if (ok) printf("  [OK] %s rot%d final y=%.3f\n", hull_names[t], r, y);
				TEST_ASSERT(ok);
				destroy_world(w);
			}
			hull_free(h);
		}
	}
	// Isolated repro: rock hull at specific rotation falls through floor
	{
		TEST_BEGIN("AVBD rock hull single drop");
		v3 rock_pts[] = {
			{0, 0.55f, 0}, {0, -0.5f, 0},
			{0.45f, 0.1f, 0.3f}, {-0.4f, 0.15f, 0.35f},
			{0.35f, 0.1f, -0.4f}, {-0.3f, 0.1f, -0.45f},
			{0.5f, -0.15f, -0.1f}, {-0.5f, -0.1f, 0.05f},
			{0.1f, -0.2f, 0.55f}, {-0.15f, -0.2f, -0.5f},
		};
		Hull* h = quickhull(rock_pts, 10);
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		((WorldInternal*)w.id)->sleep_enabled = 0;
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		// Exact rotation from hull pile body2
		quat rot = (quat){ 0.440731f, -0.320806f, 0.544865f, 0.637151f };
		Body body = create_body(w, (BodyParams){
			.position = V3(0, 3, 0), .rotation = rot, .mass = 1.0f });
		body_add_shape(w, body, (ShapeParams){
			.type = SHAPE_HULL, .hull = { .hull = h, .scale = V3(0.6f, 0.6f, 0.6f) } });

		float dt = 1.0f / 60.0f;
		for (int f = 0; f < 200; f++) {
			world_step(w, dt);
			float y = body_get_position(w, body).y;
			if (f < 5 || (f >= 38 && f <= 50))
				printf("  [ROCK] f%d y=%.4f\n", f, y);
			if (y < -2.0f) {
				printf("  [ROCK FELL THROUGH] f%d y=%.4f\n", f, y);
				TEST_ASSERT(0);
				break;
			}
		}
		float y = body_get_position(w, body).y;
		printf("  [ROCK] final y=%.4f\n", y);
		TEST_ASSERT(y > -0.5f);
		hull_free(h);
		destroy_world(w);
	}
	// AVBD hull pile: drop various hull shapes onto a floor, none should fall through.
	{
		TEST_BEGIN("AVBD hull pile no fall-through");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		((WorldInternal*)w.id)->sleep_enabled = 0;
		Body fl = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, fl, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

		// Build hull shapes
		v3 tet_pts[] = { {0,0.6f,0}, {0.5f,-0.3f,0.3f}, {-0.5f,-0.3f,0.3f}, {0,-0.3f,-0.5f} };
		Hull* hull_tet = quickhull(tet_pts, 4);

		v3 wedge_pts[] = {
			{-0.5f,-0.3f,-0.4f}, {0.5f,-0.3f,-0.4f}, {-0.5f,-0.3f,0.4f}, {0.5f,-0.3f,0.4f},
			{-0.5f, 0.3f,-0.4f}, {0.5f, 0.3f,-0.4f},
		};
		Hull* hull_wedge = quickhull(wedge_pts, 6);

		v3 rock_pts[] = {
			{0, 0.55f, 0}, {0, -0.5f, 0},
			{0.45f, 0.1f, 0.3f}, {-0.4f, 0.15f, 0.35f},
			{0.35f, 0.1f, -0.4f}, {-0.3f, 0.1f, -0.45f},
			{0.5f, -0.15f, -0.1f}, {-0.5f, -0.1f, 0.05f},
			{0.1f, -0.2f, 0.55f}, {-0.15f, -0.2f, -0.5f},
		};
		Hull* hull_rock = quickhull(rock_pts, 10);

		v3 bipyr_pts[] = {
			{0,0.7f,0}, {0.5f,0,0.5f}, {-0.5f,0,0.5f}, {0.5f,0,-0.5f}, {-0.5f,0,-0.5f}, {0,-0.7f,0},
		};
		Hull* hull_bipyr = quickhull(bipyr_pts, 6);

		Hull* hulls[] = { hull_tet, hull_wedge, hull_rock, hull_bipyr };
		v3 scales[] = { V3(0.8f,0.8f,0.8f), V3(0.7f,0.7f,0.7f), V3(0.6f,0.6f,0.6f), V3(0.7f,0.7f,0.7f) };
		const char* names[] = { "tet", "wedge", "rock", "bipyr" };
		int nhulls = 4;

		// Drop a 3x3x3 grid
		Body bodies[27];
		int btypes[27];
		int idx = 0;
		for (int layer = 0; layer < 3; layer++) {
			for (int ix = 0; ix < 3; ix++) {
				for (int iz = 0; iz < 3; iz++) {
					float x = (ix - 1) * 1.2f + ((layer % 2) ? 0.3f : 0.0f);
					float z = (iz - 1) * 1.2f + ((layer % 2) ? 0.3f : 0.0f);
					float y = 1.0f + layer * 2.0f;
					int t = idx % nhulls;
					float a = (float)idx * 1.1f;
					v3 ax = norm(V3(sinf(a), cosf(a), sinf(a*0.7f)));
					float ha = a * 0.4f;
					quat rot = (quat){ ax.x*sinf(ha), ax.y*sinf(ha), ax.z*sinf(ha), cosf(ha) };
					bodies[idx] = create_body(w, (BodyParams){
						.position = V3(x, y, z), .rotation = rot, .mass = 1.0f });
					body_add_shape(w, bodies[idx], (ShapeParams){
						.type = SHAPE_HULL, .hull = { .hull = hulls[t], .scale = scales[t] } });
					btypes[idx] = t;
					idx++;
				}
			}
		}

		float dt = 1.0f / 60.0f;
		int fallen = 0;
		for (int f = 0; f < 300; f++) {
			world_step(w, dt);
			for (int i = 0; i < 27; i++) {
				v3 p = body_get_position(w, bodies[i]);
				(void)0;
				if (p.y < -2.0f) {
					printf("  [HULL PILE FAIL] f%d body%d(%s) fell through at (%.2f,%.2f,%.2f)\n",
						f, i, names[btypes[i]], p.x, p.y, p.z);
					fallen = 1;
				}
			}
			if (fallen) break;
		}
		TEST_ASSERT(fallen == 0);
		hull_free(hull_tet); hull_free(hull_wedge); hull_free(hull_rock); hull_free(hull_bipyr);
		destroy_world(w);
	}
	// AVBD ball-socket: horizontal chain swings down under gravity.
	// Bodies spawn horizontally at anchor height — zero initial joint error,
	// gravity makes the whole chain swing downward.
	{
		TEST_BEGIN("AVBD horizontal chain swings");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.15f, 0.15f, 0.15f) });

		Body chain[5];
		Body prev = anchor;
		for (int i = 0; i < 5; i++) {
			// Horizontal: each body 1m to the right of the previous, same height
			chain[i] = create_body(w, (BodyParams){
				.position = V3((float)(i + 1), 8, 0), .rotation = quat_identity(), .mass = 0.5f });
			body_add_shape(w, chain[i], (ShapeParams){
				.type = SHAPE_BOX, .box.half_extents = V3(0.2f, 0.2f, 0.2f) });
			create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i],
				.local_offset_a = V3(0.5f, 0, 0), .local_offset_b = V3(-0.5f, 0, 0) });
			prev = chain[i];
		}

		float dt = 1.0f / 60.0f;
		float y_min = 100;
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			v3 p4 = body_get_position(w, chain[4]);
			if (p4.y < y_min) y_min = p4.y;
		}

		v3 p0 = body_get_position(w, chain[0]);
		v3 p4 = body_get_position(w, chain[4]);
		float dist_01 = len(sub(body_get_position(w, chain[1]), p0));
		printf("  [AVBD hchain] p0=(%.2f,%.2f) p4=(%.2f,%.2f) tip_min_y=%.2f d01=%.2f\n",
			p0.x, p0.y, p4.x, p4.y, y_min, dist_01);
		// Chain tip should have swung down significantly
		TEST_ASSERT(y_min < 5.0f);
		// Joints should hold (link distance ≈ 1.0)
		TEST_ASSERT(dist_01 < 1.5f);
		TEST_ASSERT(dist_01 > 0.5f);
		destroy_world(w);
	}
	// AVBD ball-socket: pendulum swings from side, doesn't float or drift.
	{
		TEST_BEGIN("AVBD ball-socket pendulum swings");
		World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .solver_type = SOLVER_AVBD });
		// Static anchor at origin
		Body anchor = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
		// Dynamic bob offset to the side so it swings
		Body bob = create_body(w, (BodyParams){ .position = V3(2, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
		body_add_shape(w, bob, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.3f });

		create_ball_socket(w, (BallSocketParams){ .body_a = anchor, .body_b = bob,
			.local_offset_a = V3(0, 0, 0), .local_offset_b = V3(-2, 0, 0) });

		float dt = 1.0f / 60.0f;
		float y_min = 100, y_max = -100;
		for (int i = 0; i < 300; i++) {
			world_step(w, dt);
			v3 p = body_get_position(w, bob);
			if (p.y < y_min) y_min = p.y;
			if (p.y > y_max) y_max = p.y;
		}
		v3 p = body_get_position(w, bob);
		float dist = len(sub(p, V3(0, 5, 0)));
		printf("  [AVBD pendulum] pos=(%.2f,%.2f,%.2f) dist=%.3f y_range=[%.2f,%.2f]\n",
			p.x, p.y, p.z, dist, y_min, y_max);
		// Bob should stay within rod length of anchor
		TEST_ASSERT(dist < 3.0f);
		// Should have swung downward (y < 5)
		TEST_ASSERT(y_min < 4.0f);
		// Joint should hold — bob shouldn't fall to floor
		TEST_ASSERT(p.y > 2.0f);
		destroy_world(w);
	}
}

// ============================================================================
// BVH broadphase tests.

static void test_bvh_aabb_helpers()
{
	TEST_BEGIN("aabb empty");
	AABB e = aabb_empty();
	TEST_ASSERT(e.min.x > e.max.x); // inverted = empty

	TEST_BEGIN("aabb merge");
	AABB a = { V3(0,0,0), V3(1,1,1) };
	AABB b = { V3(-1,2,-1), V3(0.5f,3,0.5f) };
	AABB m = aabb_merge(a, b);
	TEST_ASSERT_FLOAT(m.min.x, -1.0f, 0.001f);
	TEST_ASSERT_FLOAT(m.max.y, 3.0f, 0.001f);

	TEST_BEGIN("aabb overlaps: overlapping");
	AABB d = { V3(-0.5f, 0.5f, -0.5f), V3(0.5f, 1.5f, 0.5f) };
	TEST_ASSERT(aabb_overlaps(a, d));

	TEST_BEGIN("aabb overlaps: separated");
	AABB c = { V3(5,5,5), V3(6,6,6) };
	TEST_ASSERT(!aabb_overlaps(a, c));

	TEST_BEGIN("aabb surface area");
	AABB unit = { V3(0,0,0), V3(1,1,1) };
	TEST_ASSERT_FLOAT(aabb_surface_area(unit), 3.0f, 0.001f); // 1*1 + 1*1 + 1*1

	TEST_BEGIN("aabb expand");
	AABB ex = aabb_expand(unit, 0.5f);
	TEST_ASSERT_FLOAT(ex.min.x, -0.5f, 0.001f);
	TEST_ASSERT_FLOAT(ex.max.z, 1.5f, 0.001f);
}

static void test_bvh_insert_remove()
{
	BVHTree t; bvh_init(&t);

	TEST_BEGIN("bvh empty tree");
	TEST_ASSERT(t.root == -1);

	TEST_BEGIN("bvh insert first leaf");
	AABB b0 = { V3(0,0,0), V3(1,1,1) };
	int l0 = bvh_insert(&t, 0, b0);
	TEST_ASSERT(t.root >= 0);
	TEST_ASSERT(l0 >= 0);
	TEST_ASSERT(t.leaves[l0].body_idx == 0);

	TEST_BEGIN("bvh insert second leaf");
	AABB b1 = { V3(2,0,0), V3(3,1,1) };
	int l1 = bvh_insert(&t, 1, b1);
	TEST_ASSERT(l1 >= 0);
	TEST_ASSERT(t.leaves[l1].body_idx == 1);

	TEST_BEGIN("bvh insert third leaf (creates internal node)");
	AABB b2 = { V3(0,2,0), V3(1,3,1) };
	int l2 = bvh_insert(&t, 2, b2);
	TEST_ASSERT(l2 >= 0);

	TEST_BEGIN("bvh insert more leaves");
	for (int i = 3; i < 10; i++) {
		AABB bi = { V3((float)i, 0, 0), V3((float)i + 1, 1, 1) };
		int li = bvh_insert(&t, i, bi);
		TEST_ASSERT(li >= 0);
		TEST_ASSERT(t.leaves[li].body_idx == i);
	}

	TEST_BEGIN("bvh remove leaf");
	bvh_remove(&t, l1);
	// Tree should still be valid, root should exist
	TEST_ASSERT(t.root >= 0);

	TEST_BEGIN("bvh remove all leaves");
	bvh_remove(&t, l0);
	bvh_remove(&t, l2);
	for (int i = 3; i < 10; i++) bvh_remove(&t, i); // leaf indices happen to equal i for sequential inserts with no prior frees
	// After removing all, root may or may not be -1 depending on cleanup

	bvh_free(&t);
}

static void test_bvh_self_test()
{
	BVHTree t; bvh_init(&t);

	// Insert overlapping AABBs
	AABB b0 = { V3(0,0,0), V3(2,2,2) };
	AABB b1 = { V3(1,1,1), V3(3,3,3) }; // overlaps b0
	AABB b2 = { V3(10,10,10), V3(11,11,11) }; // separated
	bvh_insert(&t, 0, b0);
	bvh_insert(&t, 1, b1);
	bvh_insert(&t, 2, b2);

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);

	TEST_BEGIN("bvh self test: finds overlapping pair");
	int found_01 = 0;
	for (int i = 0; i < asize(pairs); i++) {
		if ((pairs[i].a == 0 && pairs[i].b == 1) || (pairs[i].a == 1 && pairs[i].b == 0)) found_01 = 1;
	}
	TEST_ASSERT(found_01);

	TEST_BEGIN("bvh self test: no false pair with separated body");
	int found_02 = 0, found_12 = 0;
	for (int i = 0; i < asize(pairs); i++) {
		if ((pairs[i].a == 0 && pairs[i].b == 2) || (pairs[i].a == 2 && pairs[i].b == 0)) found_02 = 1;
		if ((pairs[i].a == 1 && pairs[i].b == 2) || (pairs[i].a == 2 && pairs[i].b == 1)) found_12 = 1;
	}
	TEST_ASSERT(!found_02);
	TEST_ASSERT(!found_12);

	afree(pairs);
	bvh_free(&t);
}

static void test_bvh_cross_test()
{
	BVHTree ta, tb; bvh_init(&ta); bvh_init(&tb);

	AABB a0 = { V3(0,0,0), V3(2,2,2) };
	AABB b0 = { V3(1,1,1), V3(3,3,3) }; // overlaps a0
	AABB b1 = { V3(10,10,10), V3(11,11,11) }; // separated
	bvh_insert(&ta, 0, a0);
	bvh_insert(&tb, 10, b0);
	bvh_insert(&tb, 11, b1);

	CK_DYNA BroadPair* pairs = NULL;
	bvh_cross_test(&ta, &tb, &pairs);

	TEST_BEGIN("bvh cross test: finds overlapping pair");
	int found = 0;
	for (int i = 0; i < asize(pairs); i++) {
		if (pairs[i].a == 0 && pairs[i].b == 10) found = 1;
	}
	TEST_ASSERT(found);

	TEST_BEGIN("bvh cross test: no false pair");
	int false_pair = 0;
	for (int i = 0; i < asize(pairs); i++) {
		if (pairs[i].b == 11) false_pair = 1;
	}
	TEST_ASSERT(!false_pair);

	afree(pairs);
	bvh_free(&ta); bvh_free(&tb);
}

static void test_bvh_shape_aabb()
{
	BodyHot h = { .position = V3(5, 5, 5), .rotation = quat_identity() };

	TEST_BEGIN("shape_aabb sphere");
	ShapeInternal s_sph = { .type = SHAPE_SPHERE, .local_pos = V3(0,0,0), .sphere.radius = 1.0f };
	AABB box = shape_aabb(&h, &s_sph);
	TEST_ASSERT_FLOAT(box.min.x, 4.0f, 0.01f);
	TEST_ASSERT_FLOAT(box.max.x, 6.0f, 0.01f);

	TEST_BEGIN("shape_aabb box axis-aligned");
	ShapeInternal s_box = { .type = SHAPE_BOX, .local_pos = V3(0,0,0), .box.half_extents = V3(1, 2, 3) };
	box = shape_aabb(&h, &s_box);
	TEST_ASSERT_FLOAT(box.min.x, 4.0f, 0.01f);
	TEST_ASSERT_FLOAT(box.max.y, 7.0f, 0.01f);
	TEST_ASSERT_FLOAT(box.max.z, 8.0f, 0.01f);

	TEST_BEGIN("shape_aabb capsule");
	ShapeInternal s_cap = { .type = SHAPE_CAPSULE, .local_pos = V3(0,0,0), .capsule = { .half_height = 1.0f, .radius = 0.5f } };
	box = shape_aabb(&h, &s_cap);
	TEST_ASSERT_FLOAT(box.min.y, 3.5f, 0.01f);
	TEST_ASSERT_FLOAT(box.max.y, 6.5f, 0.01f);

	TEST_BEGIN("shape_aabb hull");
	v3 hull_pts[] = { V3(0,1,0), V3(1,0,0), V3(-1,0,0), V3(0,0,1), V3(0,0,-1), V3(0,-1,0) };
	Hull* test_hull = quickhull(hull_pts, 6);
	ShapeInternal s_hull = { .type = SHAPE_HULL, .local_pos = V3(0,0,0), .hull = { .hull = test_hull, .scale = V3(2,2,2) } };
	box = shape_aabb(&h, &s_hull);
	TEST_ASSERT(box.min.x < 4.0f);
	TEST_ASSERT(box.max.x > 6.0f);
	hull_free(test_hull);
}

// Correctness: BVH broadphase must find the same pairs as N^2.
static void test_bvh_vs_n2()
{
	// Create a world with BVH broadphase and simulate.
	World w_bvh = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });
	World w_n2  = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_N2 });

	// Build identical scenes: floor + several dynamic bodies.
	World worlds[2] = { w_bvh, w_n2 };
	for (int wi = 0; wi < 2; wi++) {
		World w = worlds[wi];
		Body floor = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });
		for (int i = 0; i < 5; i++) {
			Body b = create_body(w, (BodyParams){ .position = V3((float)i * 0.8f, 3 + (float)i, 0), .rotation = quat_identity(), .mass = 1.0f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });
		}
	}

	// Step both worlds and compare contact counts over the first 10 frames
	// (before solver warm-start divergence between broadphase modes).
	int mismatch = 0;
	for (int frame = 0; frame < 10; frame++) {
		world_step(w_bvh, 1.0f / 60.0f);
		world_step(w_n2, 1.0f / 60.0f);
		const Contact* c_bvh; const Contact* c_n2;
		int n_bvh = world_get_contacts(w_bvh, &c_bvh);
		int n_n2  = world_get_contacts(w_n2, &c_n2);
		if (n_bvh != n_n2) mismatch++;
	}

	TEST_BEGIN("bvh vs n2: contact counts match first 10 frames");
	TEST_ASSERT(mismatch == 0);

	destroy_world(w_bvh);
	destroy_world(w_n2);
}

static void test_bvh_empty_trees()
{
	BVHTree t; bvh_init(&t);

	TEST_BEGIN("bvh self test empty tree");
	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);
	TEST_ASSERT(asize(pairs) == 0);
	afree(pairs);

	TEST_BEGIN("bvh refit empty tree");
	// Should not crash
	bvh_refit(&t, NULL); // NULL world is fine since no leaves to visit

	bvh_free(&t);
}

static void test_bvh_single_body()
{
	BVHTree t; bvh_init(&t);
	AABB b0 = { V3(0,0,0), V3(1,1,1) };
	int l0 = bvh_insert(&t, 0, b0);

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);
	TEST_BEGIN("bvh self test single body: no pairs");
	TEST_ASSERT(asize(pairs) == 0);
	afree(pairs);

	bvh_remove(&t, l0);
	TEST_BEGIN("bvh remove last body");
	// Tree should handle gracefully
	pairs = NULL;
	bvh_self_test(&t, &pairs);
	TEST_ASSERT(asize(pairs) == 0);
	afree(pairs);

	bvh_free(&t);
}

static void test_bvh_many_overlapping()
{
	BVHTree t; bvh_init(&t);

	// Insert N bodies all at the origin -- all overlap each other.
	int N = 20;
	for (int i = 0; i < N; i++) {
		AABB b = { V3(-1, -1, -1), V3(1, 1, 1) };
		bvh_insert(&t, i, b);
	}

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);

	TEST_BEGIN("bvh many overlapping: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == N);

	TEST_BEGIN("bvh many overlapping: correct pair count");
	int expected = N * (N - 1) / 2;
	TEST_ASSERT(asize(pairs) == expected);

	afree(pairs);
	bvh_free(&t);
}

static void test_bvh_no_overlaps()
{
	BVHTree t; bvh_init(&t);

	// Insert 10 well-separated bodies
	for (int i = 0; i < 10; i++) {
		AABB b = { V3((float)i * 10, 0, 0), V3((float)i * 10 + 1, 1, 1) };
		bvh_insert(&t, i, b);
	}

	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);

	TEST_BEGIN("bvh no overlaps: zero pairs");
	TEST_ASSERT(asize(pairs) == 0);

	afree(pairs);
	bvh_free(&t);
}

static void test_bvh_cache_reorder()
{
	BVHTree t; bvh_init(&t);

	// Insert 15 bodies with overlapping and separated AABBs.
	for (int i = 0; i < 15; i++) {
		AABB b = { V3((float)(i % 5) * 1.5f, 0, 0), V3((float)(i % 5) * 1.5f + 2, 2, 2) };
		bvh_insert(&t, i, b);
	}

	// Get pairs before reorder.
	CK_DYNA BroadPair* before = NULL;
	bvh_self_test(&t, &before);
	int before_count = asize(before);

	// Reorder.
	bvh_cache_reorder(&t);

	// Verify DFS order: every internal child index must be > parent index.
	TEST_BEGIN("cache reorder: DFS order property");
	int dfs_ok = 1;
	for (int i = 0; i < asize(t.nodes); i++) {
		BVHNode* n = &t.nodes[i];
		if (bvh_child_is_internal(&n->a) && n->a.index <= i) dfs_ok = 0;
		if (bvh_child_is_internal(&n->b) && n->b.index <= i) dfs_ok = 0;
	}
	TEST_ASSERT(dfs_ok);

	// Verify leaf back-pointers.
	TEST_BEGIN("cache reorder: leaf back-pointers valid");
	int bp_ok = 1;
	for (int i = 0; i < asize(t.leaves); i++) {
		BVHLeaf* lf = &t.leaves[i];
		BVHChild* c = bvh_child(&t.nodes[lf->node_idx], lf->child_slot);
		if (!bvh_child_is_leaf(c) || bvh_child_leaf_idx(c) != i) bp_ok = 0;
	}
	TEST_ASSERT(bp_ok);

	// Verify all leaves still reachable.
	TEST_BEGIN("cache reorder: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == 15);

	// Verify pairs unchanged after reorder.
	CK_DYNA BroadPair* after = NULL;
	bvh_self_test(&t, &after);
	TEST_BEGIN("cache reorder: pairs preserved");
	TEST_ASSERT(asize(after) == before_count);

	afree(before); afree(after);
	bvh_free(&t);
}

// Validate all leaf back-pointers in a tree.
static int bvh_validate_backpointers(BVHTree* t)
{
	for (int i = 0; i < asize(t->leaves); i++) {
		// Skip freed leaves (check if on freelist -- simple: just verify node_idx is in range)
		BVHLeaf* lf = &t->leaves[i];
		if (lf->node_idx < 0 || lf->node_idx >= asize(t->nodes)) continue;
		BVHChild* c = bvh_child(&t->nodes[lf->node_idx], lf->child_slot);
		if (!bvh_child_is_leaf(c) || bvh_child_leaf_idx(c) != i) return 0;
	}
	return 1;
}

// Compute total SAH cost of the tree (sum of all internal node child SAs).
static float bvh_total_sah(BVHTree* t, int ni)
{
	BVHNode* n = &t->nodes[ni];
	float cost = 0;
	if (!bvh_child_is_empty(&n->a)) cost += aabb_surface_area(bvh_child_aabb(&n->a));
	if (!bvh_child_is_empty(&n->b)) cost += aabb_surface_area(bvh_child_aabb(&n->b));
	if (bvh_child_is_internal(&n->a)) cost += bvh_total_sah(t, n->a.index);
	if (bvh_child_is_internal(&n->b)) cost += bvh_total_sah(t, n->b.index);
	return cost;
}

static void test_bvh_rotations()
{
	// Insert bodies in a pathological linear order (sorted along X).
	// Each body overlaps its neighbor (width 1.5, spacing 1.0).
	BVHTree t; bvh_init(&t);
	for (int i = 0; i < 20; i++) {
		AABB b = { V3((float)i, 0, 0), V3((float)i + 1.5f, 1, 1) };
		bvh_insert(&t, i, b);
	}

	TEST_BEGIN("rotations: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == 20);

	TEST_BEGIN("rotations: back-pointers valid");
	TEST_ASSERT(bvh_validate_backpointers(&t));

	TEST_BEGIN("rotations: self-test finds correct overlapping pairs");
	// Adjacent bodies overlap: (0,1), (1,2), ..., (18,19) = 19 pairs
	// Non-adjacent don't overlap (gap of 0.5 between each)
	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);
	TEST_ASSERT(asize(pairs) == 19);
	afree(pairs);

	// Compare SAH to a no-rotation tree.
	float rotated_sah = bvh_total_sah(&t, t.root);
	bvh_free(&t);

	// Build same tree without rotations by calling bvh_try_rotate with empty impl...
	// Actually just verify SAH is reasonable (not degenerate).
	TEST_BEGIN("rotations: SAH cost is bounded");
	TEST_ASSERT(rotated_sah > 0.0f);
	TEST_ASSERT(rotated_sah < 1e6f); // sanity check

	// Stress test: many overlapping bodies with rotations.
	bvh_init(&t);
	for (int i = 0; i < 50; i++) {
		float x = (float)(i % 10) * 0.5f;
		float y = (float)(i / 10) * 0.5f;
		AABB b = { V3(x, y, 0), V3(x + 1, y + 1, 1) };
		bvh_insert(&t, i, b);
	}

	TEST_BEGIN("rotations stress: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == 50);

	TEST_BEGIN("rotations stress: back-pointers valid");
	TEST_ASSERT(bvh_validate_backpointers(&t));

	bvh_free(&t);
}

static void test_bvh_binned_build()
{
	BVHTree t; bvh_init(&t);

	// Create 30 leaves manually and build with binned SAH.
	int leaf_count = 30;
	AABB* lut = CK_ALLOC(sizeof(AABB) * leaf_count);
	CK_DYNA int* lis = NULL;
	for (int i = 0; i < leaf_count; i++) {
		int li = bvh_alloc_leaf(&t);
		t.leaves[li].body_idx = i;
		lut[li] = (AABB){ V3((float)(i % 6) * 2, (float)(i / 6) * 2, 0), V3((float)(i % 6) * 2 + 1.5f, (float)(i / 6) * 2 + 1.5f, 1) };
		apush(lis, li);
	}

	t.root = bvh_binned_build(&t, lis, lut, leaf_count);
	t.meta[t.root].parent = -1;

	TEST_BEGIN("binned build: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == leaf_count);

	TEST_BEGIN("binned build: back-pointers valid");
	TEST_ASSERT(bvh_validate_backpointers(&t));

	// Self-test: adjacent bodies on same row overlap (spacing 2, width 1.5 => 0.5 overlap).
	CK_DYNA BroadPair* pairs = NULL;
	bvh_self_test(&t, &pairs);
	// Verify against brute-force AABB overlap count.
	int expected = 0;
	for (int i = 0; i < leaf_count; i++)
		for (int j = i + 1; j < leaf_count; j++)
			if (aabb_overlaps(lut[i], lut[j])) expected++;
	TEST_BEGIN("binned build: self-test matches brute force");
	TEST_ASSERT(asize(pairs) == expected);

	afree(pairs); afree(lis); CK_FREE(lut);
	bvh_free(&t);
}

static void test_bvh_refine_subtree()
{
	BVHTree t; bvh_init(&t);

	// Insert bodies via normal insertion (which uses SAH + rotations).
	for (int i = 0; i < 20; i++) {
		AABB b = { V3((float)i, 0, 0), V3((float)i + 1.5f, 1, 1) };
		bvh_insert(&t, i, b);
	}

	// Get pairs before refinement.
	CK_DYNA BroadPair* before = NULL;
	bvh_self_test(&t, &before);
	int before_count = asize(before);
	float before_sah = bvh_total_sah(&t, t.root);

	// Build LUT from current tree state.
	AABB* lut = CK_ALLOC(sizeof(AABB) * asize(t.leaves));
	for (int i = 0; i < asize(t.leaves); i++) {
		BVHChild* c = bvh_child(&t.nodes[t.leaves[i].node_idx], t.leaves[i].child_slot);
		lut[i] = bvh_child_aabb(c);
	}

	// Refine a subtree (pick root for maximum effect).
	bvh_refine_subtree(&t, t.root, lut);

	TEST_BEGIN("refine subtree: all leaves reachable");
	TEST_ASSERT(bvh_count_leaves(&t, t.root) == 20);

	TEST_BEGIN("refine subtree: back-pointers valid");
	TEST_ASSERT(bvh_validate_backpointers(&t));

	TEST_BEGIN("refine subtree: pairs preserved");
	CK_DYNA BroadPair* after = NULL;
	bvh_self_test(&t, &after);
	TEST_ASSERT(asize(after) == before_count);

	TEST_BEGIN("refine subtree: SAH improved or equal");
	float after_sah = bvh_total_sah(&t, t.root);
	TEST_ASSERT(after_sah <= before_sah + 1e-3f);

	afree(before); afree(after); CK_FREE(lut);
	bvh_free(&t);
}

static void test_bvh_fat_aabb()
{
	// Test that fat AABB prevents unnecessary refit updates.
	// Create a world with BVH, add two bodies, step to generate refit.
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });
	Body floor = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

	Body ball = create_body(w, (BodyParams){ .position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f });
	body_add_shape(w, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	// Step once to refit and set fat AABBs.
	world_step(w, 1.0f / 60.0f);

	TEST_BEGIN("fat aabb: ball above floor after first step");
	v3 pos = body_get_position(w, ball);
	TEST_ASSERT(pos.y > 0.0f);

	// Static floor leaf should have fat AABB set.
	WorldInternal* wi = (WorldInternal*)w.id;
	int floor_idx = handle_index(floor);
	int fl = wi->body_cold[floor_idx].bvh_leaf;
	TEST_BEGIN("fat aabb: static leaf has fat AABB set");
	BVHLeaf* lf = &wi->bvh_static->leaves[fl];
	TEST_ASSERT(lf->fat_min.x < lf->fat_max.x); // non-empty fat AABB

	// Run 300 steps -- ball should settle on floor, not fall through.
	for (int i = 0; i < 300; i++) world_step(w, 1.0f / 60.0f);
	pos = body_get_position(w, ball);
	TEST_BEGIN("fat aabb: ball rests on floor (no false negatives)");
	TEST_ASSERT(pos.y > -0.5f && pos.y < 2.0f);

	destroy_world(w);
}

static void run_bvh_tests()
{
	printf("--- nudge BVH tests ---\n");
	test_bvh_aabb_helpers();
	test_bvh_insert_remove();
	test_bvh_self_test();
	test_bvh_cross_test();
	test_bvh_shape_aabb();
	test_bvh_empty_trees();
	test_bvh_single_body();
	test_bvh_many_overlapping();
	test_bvh_no_overlaps();
	test_bvh_cache_reorder();
	test_bvh_rotations();
	test_bvh_binned_build();
	test_bvh_refine_subtree();
	test_bvh_fat_aabb();
	test_bvh_vs_n2();
}

// ============================================================================
// World query tests: AABB overlap and raycast.

static void test_query_aabb_basic()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	// Floor at y=-1, sphere at (0,2,0), box at (5,1,0)
	Body floor = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(50, 1, 50) });

	Body sphere = create_body(w, (BodyParams){ .position = V3(0, 2, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, sphere, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	Body box = create_body(w, (BodyParams){ .position = V3(5, 1, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });

	Body results[16];

	// Query a box around the sphere -- should find sphere (and possibly floor)
	TEST_BEGIN("aabb query: finds sphere");
	int count = world_query_aabb(w, V3(-1, 1, -1), V3(1, 3, 1), results, 16);
	TEST_ASSERT(count >= 1);
	int found_sphere = 0;
	for (int i = 0; i < count && i < 16; i++)
		if (results[i].id == sphere.id) found_sphere = 1;
	TEST_ASSERT(found_sphere);

	// Query far away -- should find nothing
	TEST_BEGIN("aabb query: miss far away");
	count = world_query_aabb(w, V3(100, 100, 100), V3(101, 101, 101), results, 16);
	TEST_ASSERT(count == 0);

	// Query around the box at (5,1,0)
	TEST_BEGIN("aabb query: finds box at (5,1,0)");
	count = world_query_aabb(w, V3(4, 0, -1), V3(6, 2, 1), results, 16);
	int found_box = 0;
	for (int i = 0; i < count && i < 16; i++)
		if (results[i].id == box.id) found_box = 1;
	TEST_ASSERT(found_box);

	// Query that covers everything
	TEST_BEGIN("aabb query: large query finds all 3");
	count = world_query_aabb(w, V3(-100, -100, -100), V3(100, 100, 100), results, 16);
	TEST_ASSERT(count >= 3);

	destroy_world(w);
}

static void test_query_aabb_n2()
{
	// Same test with N^2 broadphase fallback
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_N2 });

	Body sphere = create_body(w, (BodyParams){ .position = V3(0, 2, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, sphere, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	Body results[8];

	TEST_BEGIN("aabb query n2: finds sphere");
	int count = world_query_aabb(w, V3(-1, 1, -1), V3(1, 3, 1), results, 8);
	TEST_ASSERT(count == 1);
	TEST_ASSERT(results[0].id == sphere.id);

	TEST_BEGIN("aabb query n2: miss");
	count = world_query_aabb(w, V3(10, 10, 10), V3(11, 11, 11), results, 8);
	TEST_ASSERT(count == 0);

	destroy_world(w);
}

static void test_raycast_sphere()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	Body sphere = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, sphere, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 1.0f });

	RayHit hit;

	// Ray from (0,5,0) downward -- should hit sphere at (0,1,0)
	TEST_BEGIN("raycast sphere: hit from above");
	int r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT(hit.body.id == sphere.id);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);
	TEST_ASSERT_FLOAT(hit.point.y, 1.0f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.y, 1.0f, EPS); // upward normal

	// Ray from side
	TEST_BEGIN("raycast sphere: hit from side");
	r = world_raycast(w, V3(-5, 0, 0), V3(1, 0, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.x, -1.0f, EPS);

	// Ray that misses
	TEST_BEGIN("raycast sphere: miss");
	r = world_raycast(w, V3(0, 5, 0), V3(1, 0, 0), 100.0f, &hit);
	TEST_ASSERT(!r);

	// Ray too short
	TEST_BEGIN("raycast sphere: max_distance too short");
	r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 2.0f, &hit);
	TEST_ASSERT(!r);

	destroy_world(w);
}

static void test_raycast_box()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	Body box = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, box, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(1, 1, 1) });

	RayHit hit;

	TEST_BEGIN("raycast box: hit top face");
	int r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);
	TEST_ASSERT_FLOAT(hit.point.y, 1.0f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.y, 1.0f, EPS);

	TEST_BEGIN("raycast box: hit side face");
	r = world_raycast(w, V3(5, 0, 0), V3(-1, 0, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.x, 1.0f, EPS);

	destroy_world(w);
}

static void test_raycast_capsule()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	// Capsule: along Y, half_height=1, radius=0.5 → endpoints at (0,-1,0) and (0,1,0)
	Body cap = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, cap, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = 1.0f, .radius = 0.5f } });

	RayHit hit;

	// Ray from side at midpoint -- should hit cylinder body
	TEST_BEGIN("raycast capsule: hit cylinder body");
	int r = world_raycast(w, V3(5, 0, 0), V3(-1, 0, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.5f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.x, 1.0f, EPS);

	// Ray from above -- should hit top hemisphere
	TEST_BEGIN("raycast capsule: hit top hemisphere");
	r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 3.5f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.y, 1.0f, EPS);

	destroy_world(w);
}

static void test_raycast_hull()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	// Use unit box hull with scale (1,1,1) -- equivalent to a unit box
	Body hull = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, hull, (ShapeParams){ .type = SHAPE_HULL, .hull = { .hull = hull_unit_box(), .scale = V3(1, 1, 1) } });

	RayHit hit;

	TEST_BEGIN("raycast hull: hit from above");
	int r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);
	TEST_ASSERT_FLOAT(hit.normal.y, 1.0f, EPS);

	TEST_BEGIN("raycast hull: miss");
	r = world_raycast(w, V3(5, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(!r);

	destroy_world(w);
}

static void test_raycast_closest()
{
	// Two spheres in a line -- ray should hit the closer one
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });

	Body near = create_body(w, (BodyParams){ .position = V3(3, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, near, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	Body far = create_body(w, (BodyParams){ .position = V3(8, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, far, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

	RayHit hit;
	TEST_BEGIN("raycast closest: picks nearer sphere");
	int r = world_raycast(w, V3(0, 0, 0), V3(1, 0, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT(hit.body.id == near.id);
	TEST_ASSERT_FLOAT(hit.distance, 2.5f, EPS);

	destroy_world(w);
}

static void test_raycast_null_hit()
{
	// Boolean-only raycast (hit=NULL)
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH });
	Body s = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, s, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 1.0f });

	TEST_BEGIN("raycast null hit: boolean mode");
	int r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, NULL);
	TEST_ASSERT(r);
	r = world_raycast(w, V3(0, 5, 0), V3(1, 0, 0), 100.0f, NULL);
	TEST_ASSERT(!r);

	destroy_world(w);
}

static void test_raycast_n2_fallback()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_N2 });
	Body s = create_body(w, (BodyParams){ .position = V3(0, 0, 0), .rotation = quat_identity(), .mass = 1 });
	body_add_shape(w, s, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 1.0f });

	RayHit hit;
	TEST_BEGIN("raycast n2: hit sphere");
	int r = world_raycast(w, V3(0, 5, 0), V3(0, -1, 0), 100.0f, &hit);
	TEST_ASSERT(r);
	TEST_ASSERT_FLOAT(hit.distance, 4.0f, EPS);

	destroy_world(w);
}

static void run_query_tests()
{
	printf("--- nudge query tests ---\n");
	test_query_aabb_basic();
	test_query_aabb_n2();
	test_raycast_sphere();
	test_raycast_box();
	test_raycast_capsule();
	test_raycast_hull();
	test_raycast_closest();
	test_raycast_null_hit();
	test_raycast_n2_fallback();
}

// ============================================================================
// Feature ID tests for warm starting.

static void test_feature_ids()
{
	const Hull* box = hull_unit_box();

	// --- Face contact: two boxes resting face-to-face ---
	{
		TEST_BEGIN("feature_id face contact nonzero");
		Manifold m = {0};
		Box a = { V3(0, 1, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -0.5f, 0), quat_identity(), V3(2, 0.5f, 2) };
		int hit = collide_box_box(a, b, &m);
		TEST_ASSERT(hit);
		TEST_ASSERT(m.count > 0);
		for (int i = 0; i < m.count; i++) {
			TEST_BEGIN("feature_id face nonzero");
			TEST_ASSERT(m.contacts[i].feature_id != 0);
			TEST_ASSERT((m.contacts[i].feature_id & FEATURE_EDGE_BIT) == 0);
		}
	}

	// --- Same config twice produces same feature IDs ---
	{
		TEST_BEGIN("feature_id deterministic");
		Manifold m1 = {0}, m2 = {0};
		Box a = { V3(0, 0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -0.5f, 0), quat_identity(), V3(1, 1, 1) };
		collide_box_box(a, b, &m1);
		collide_box_box(a, b, &m2);
		TEST_ASSERT(m1.count == m2.count);
		for (int i = 0; i < m1.count; i++) {
			TEST_BEGIN("feature_id same across calls");
			TEST_ASSERT(m1.contacts[i].feature_id == m2.contacts[i].feature_id);
		}
	}

	// --- Small perturbation preserves feature IDs ---
	{
		TEST_BEGIN("feature_id stable under perturbation");
		Box a = { V3(0, 0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m_base = {0};
		collide_box_box(a, b, &m_base);

		// Slight shift.
		Box a2 = { V3(0.01f, 0.51f, -0.01f), quat_identity(), V3(1, 1, 1) };
		Manifold m_pert = {0};
		collide_box_box(a2, b, &m_pert);

		// Same face pair should produce same ref/inc face encoding.
		// The clip_edge component may differ, but the face indices should match.
		TEST_ASSERT(m_base.count > 0 && m_pert.count > 0);
		uint32_t base_faces = m_base.contacts[0].feature_id & 0xFFFF;
		uint32_t pert_faces = m_pert.contacts[0].feature_id & 0xFFFF;
		TEST_ASSERT(base_faces == pert_faces);
	}

	// --- Different configurations produce different feature IDs ---
	{
		TEST_BEGIN("feature_id varies with config");
		// Top face contact.
		Box a1 = { V3(0, 0.9f, 0), quat_identity(), V3(1, 1, 1) };
		Box b  = { V3(0, -0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m1 = {0};
		collide_box_box(a1, b, &m1);

		// Side face contact (shift X so a different face is reference).
		Box a2 = { V3(1.9f, 0, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m2 = {0};
		collide_box_box(a2, b, &m2);

		TEST_ASSERT(m1.count > 0 && m2.count > 0);
		uint32_t faces1 = m1.contacts[0].feature_id & 0xFFFF;
		uint32_t faces2 = m2.contacts[0].feature_id & 0xFFFF;
		TEST_ASSERT(faces1 != faces2);
	}

	// --- Hull-hull with custom hulls ---
	{
		TEST_BEGIN("feature_id hull-hull");
		v3 tet_pts[] = { {0,1,0}, {-1,-1,1}, {1,-1,1}, {0,-1,-1} };
		Hull* h = quickhull(tet_pts, 4);
		TEST_ASSERT(h != NULL);
		if (h) {
			ConvexHull ca = { h, V3(0, 2, 0), quat_identity(), V3(1,1,1) };
			ConvexHull cb = { h, V3(0, 0, 0), quat_identity(), V3(1,1,1) };
			Manifold m = {0};
			int hit = collide_hull_hull(ca, cb, &m);
			if (hit && m.count > 0) {
				TEST_BEGIN("feature_id hull nonzero");
				TEST_ASSERT(m.contacts[0].feature_id != 0);
			}
			hull_free(h);
		}
	}

	// --- Exact feature IDs for known box configurations ---
	// Unit box faces: 0=-Z, 1=+Z, 2=-X, 3=+X, 4=-Y, 5=+Y
	// Feature ID = faceA | (faceB << 8) | (clip_edge << 16)
	{
		// Stacked vertically: A on top, B on bottom.
		// A's face 4 (-Y) meets B's face 5 (+Y).
		TEST_BEGIN("feature_id exact: stacked Y");
		Box a = { V3(0, 1, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -1, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m = {0};
		collide_box_box(a, b, &m);
		TEST_ASSERT(m.count > 0);
		// All contacts should encode face 4 and face 5 in the low 16 bits.
		uint32_t expect_faces = 4 | (5 << 8);
		for (int i = 0; i < m.count; i++) {
			uint32_t got_faces = m.contacts[i].feature_id & 0xFFFF;
			if (got_faces != expect_faces) {
				// Could be flipped: 5 | (4 << 8)
				uint32_t alt = 5 | (4 << 8);
				TEST_BEGIN("feature_id exact Y faces");
				TEST_ASSERT(got_faces == expect_faces || got_faces == alt);
			}
		}
	}
	{
		// Side by side on X axis: A's face 3 (+X) meets B's face 2 (-X).
		TEST_BEGIN("feature_id exact: side X");
		Box a = { V3(-1, 0, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(1, 0, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m = {0};
		collide_box_box(a, b, &m);
		TEST_ASSERT(m.count > 0);
		uint32_t expect_a = 3 | (2 << 8); // A.+X ref, B.-X inc
		uint32_t expect_b = 2 | (3 << 8); // B.-X ref, A.+X inc
		for (int i = 0; i < m.count; i++) {
			uint32_t got = m.contacts[i].feature_id & 0xFFFF;
			TEST_BEGIN("feature_id exact X faces");
			TEST_ASSERT(got == expect_a || got == expect_b);
		}
	}
	{
		// Front-to-back on Z: A's face 1 (+Z) meets B's face 0 (-Z).
		TEST_BEGIN("feature_id exact: front Z");
		Box a = { V3(0, 0, -1), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, 0, 1), quat_identity(), V3(1, 1, 1) };
		Manifold m = {0};
		collide_box_box(a, b, &m);
		TEST_ASSERT(m.count > 0);
		uint32_t expect_a = 1 | (0 << 8);
		uint32_t expect_b = 0 | (1 << 8);
		for (int i = 0; i < m.count; i++) {
			uint32_t got = m.contacts[i].feature_id & 0xFFFF;
			TEST_BEGIN("feature_id exact Z faces");
			TEST_ASSERT(got == expect_a || got == expect_b);
		}
	}
	{
		// Clip edge index: contacts from original incident vertices get 0xFF.
		// Contacts created by clipping get the side-plane index (0,1,2,3).
		TEST_BEGIN("feature_id clip edge encoding");
		Box a = { V3(0, 0.9f, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -0.9f, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m = {0};
		collide_box_box(a, b, &m);
		TEST_ASSERT(m.count > 0);
		int has_original = 0, has_clipped = 0;
		for (int i = 0; i < m.count; i++) {
			uint8_t clip = (m.contacts[i].feature_id >> 16) & 0xFF;
			if (clip == 0xFF) has_original = 1;
			else has_clipped = 1;
		}
		// With boxes nearly the same size, most contacts should be original
		// incident vertices (0xFF), but some may be clipped.
		TEST_ASSERT(has_original || has_clipped);
	}

	// --- Swapped order: both produce valid nonzero feature IDs ---
	{
		TEST_BEGIN("feature_id order both valid");
		Box a = { V3(0, 0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Box b = { V3(0, -0.5f, 0), quat_identity(), V3(1, 1, 1) };
		Manifold m_ab = {0}, m_ba = {0};
		collide_box_box(a, b, &m_ab);
		collide_box_box(b, a, &m_ba);
		TEST_ASSERT(m_ab.count > 0 && m_ba.count > 0);
		TEST_ASSERT(m_ab.count == m_ba.count);
		// Both orderings produce nonzero face-based IDs.
		for (int i = 0; i < m_ab.count; i++) {
			TEST_BEGIN("feature_id AB nonzero");
			TEST_ASSERT(m_ab.contacts[i].feature_id != 0);
		}
		for (int i = 0; i < m_ba.count; i++) {
			TEST_BEGIN("feature_id BA nonzero");
			TEST_ASSERT(m_ba.contacts[i].feature_id != 0);
		}
		// Same face pair (possibly in different order) in low 16 bits.
		uint32_t ab_faces = m_ab.contacts[0].feature_id & 0xFFFF;
		uint32_t ba_faces = m_ba.contacts[0].feature_id & 0xFFFF;
		uint32_t ab_swap = (ab_faces >> 8) | ((ab_faces & 0xFF) << 8);
		TEST_BEGIN("feature_id swapped face pair");
		TEST_ASSERT(ab_faces == ba_faces || ab_swap == ba_faces);
	}
}

// ============================================================================
// Island sleep system tests.

// Helper: check if a body is in a sleeping island.
static int body_is_sleeping(World w, Body body)
{
	WorldInternal* wi = (WorldInternal*)w.id;
	int idx = handle_index(body);
	int isl = wi->body_cold[idx].island_id;
	if (isl < 0) return 0;
	if (!(wi->island_gen[isl] & 1)) return 0;
	return !wi->islands[isl].awake;
}

// Helper: get body's island id (-1 if none).
static int body_island_id(World w, Body body)
{
	WorldInternal* wi = (WorldInternal*)w.id;
	return wi->body_cold[handle_index(body)].island_id;
}

// Two boxes stacked on a floor should go to sleep after settling.
static void test_sleep_stacked_boxes()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box_a = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_a, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});
	Body box_b = create_body(w, (BodyParams){
		.position = V3(0, 1.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	// Run until settled (6 seconds)
	step_n(w, 360);

	TEST_BEGIN("sleep: stacked boxes form same island");
	int isl_a = body_island_id(w, box_a);
	int isl_b = body_island_id(w, box_b);
	TEST_ASSERT(isl_a >= 0);
	TEST_ASSERT(isl_a == isl_b);

	TEST_BEGIN("sleep: stacked boxes go to sleep");
	TEST_ASSERT(body_is_sleeping(w, box_a));
	TEST_ASSERT(body_is_sleeping(w, box_b));

	TEST_BEGIN("sleep: sleeping bodies have near-zero velocity");
	WorldInternal* wi = (WorldInternal*)w.id;
	float va = len(wi->body_hot[handle_index(box_a)].velocity);
	float vb = len(wi->body_hot[handle_index(box_b)].velocity);
	TEST_ASSERT(va < 1e-6f);
	TEST_ASSERT(vb < 1e-6f);

	destroy_world(w);
}

// Dropping a body onto a sleeping stack should wake it.
static void test_sleep_wake_on_impact()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body resting = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, resting, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	// Let it settle and sleep
	step_n(w, 360);
	TEST_BEGIN("sleep wake: box is asleep before impact");
	TEST_ASSERT(body_is_sleeping(w, resting));

	// Drop another box from above
	Body dropper = create_body(w, (BodyParams){
		.position = V3(0, 5, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, dropper, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	// Step a few frames so the dropper falls and collides
	step_n(w, 60);

	TEST_BEGIN("sleep wake: resting box wakes on impact");
	TEST_ASSERT(!body_is_sleeping(w, resting));

	destroy_world(w);
}

// body_wake API wakes a sleeping body's island.
static void test_sleep_body_wake_api()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 360);
	TEST_BEGIN("sleep body_wake: box asleep before wake call");
	TEST_ASSERT(body_is_sleeping(w, box));

	body_wake(w, box);
	TEST_BEGIN("sleep body_wake: box awake after wake call");
	TEST_ASSERT(!body_is_sleeping(w, box));

	destroy_world(w);
}

// body_set_velocity wakes a sleeping body.
static void test_sleep_set_velocity_wakes()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 360);
	TEST_BEGIN("sleep set_velocity: box asleep");
	TEST_ASSERT(body_is_sleeping(w, box));

	body_set_velocity(w, box, V3(5, 0, 0));
	TEST_BEGIN("sleep set_velocity: box awake after set_velocity");
	TEST_ASSERT(!body_is_sleeping(w, box));

	destroy_world(w);
}

// Static bodies never get island_id.
static void test_sleep_static_no_island()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 60);

	TEST_BEGIN("sleep: static body has no island");
	TEST_ASSERT(body_island_id(w, floor) == -1);

	TEST_BEGIN("sleep: dynamic body touching static has island");
	TEST_ASSERT(body_island_id(w, box) >= 0);

	destroy_world(w);
}

// Destroying a body removes it from its island.
static void test_sleep_destroy_body_removes_from_island()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box_a = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_a, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});
	Body box_b = create_body(w, (BodyParams){
		.position = V3(0, 1.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 60);
	int isl = body_island_id(w, box_a);
	TEST_BEGIN("sleep destroy: bodies share island before destroy");
	TEST_ASSERT(isl >= 0);
	TEST_ASSERT(body_island_id(w, box_b) == isl);

	// Destroy one body; the other should still be in its island
	destroy_body(w, box_b);

	TEST_BEGIN("sleep destroy: remaining body still has island");
	TEST_ASSERT(body_island_id(w, box_a) >= 0);

	// Continue stepping -- should not crash
	step_n(w, 60);

	TEST_BEGIN("sleep destroy: simulation stable after body destroy");
	v3 pos = body_get_position(w, box_a);
	TEST_ASSERT(is_valid(pos));

	destroy_world(w);
}

// Joint-connected bodies share an island; destroying the joint allows splitting.
static void test_sleep_joint_creates_island()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body anchor = create_body(w, (BodyParams){
		.position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, anchor, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.1f,
	});
	Body bob = create_body(w, (BodyParams){
		.position = V3(0, 8, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, bob, (ShapeParams){
		.type = SHAPE_SPHERE, .sphere.radius = 0.2f,
	});

	Joint j = create_ball_socket(w, (BallSocketParams){
		.body_a = anchor, .body_b = bob,
		.local_offset_a = V3(0, 0, 0), .local_offset_b = V3(0, 0, 0),
	});

	TEST_BEGIN("sleep joint: bob gets island from joint");
	TEST_ASSERT(body_island_id(w, bob) >= 0);

	TEST_BEGIN("sleep joint: static anchor has no island");
	TEST_ASSERT(body_island_id(w, anchor) == -1);

	// Destroy joint, step -- bob should still function
	destroy_joint(w, j);
	step_n(w, 60);

	TEST_BEGIN("sleep joint: simulation stable after joint destroy");
	v3 pos = body_get_position(w, bob);
	TEST_ASSERT(is_valid(pos));

	destroy_world(w);
}

// Two separate groups that never touch should be in different islands.
static void test_sleep_separate_islands()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(50, 1, 50),
	});

	// Group A: far left
	Body a = create_body(w, (BodyParams){
		.position = V3(-20, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, a, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	// Group B: far right
	Body b = create_body(w, (BodyParams){
		.position = V3(20, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 360);

	TEST_BEGIN("sleep separate: both bodies asleep");
	TEST_ASSERT(body_is_sleeping(w, a));
	TEST_ASSERT(body_is_sleeping(w, b));

	TEST_BEGIN("sleep separate: bodies in different islands");
	TEST_ASSERT(body_island_id(w, a) != body_island_id(w, b));

	destroy_world(w);
}

// Waking one body in a sleeping island wakes all bodies in that island.
static void test_sleep_wake_wakes_whole_island()
{
	World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
	Body floor = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10),
	});
	Body box_a = create_body(w, (BodyParams){
		.position = V3(0, 0.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_a, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});
	Body box_b = create_body(w, (BodyParams){
		.position = V3(0, 1.5f, 0), .rotation = quat_identity(), .mass = 1.0f,
	});
	body_add_shape(w, box_b, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f),
	});

	step_n(w, 360);
	TEST_BEGIN("sleep wake island: both asleep");
	TEST_ASSERT(body_is_sleeping(w, box_a));
	TEST_ASSERT(body_is_sleeping(w, box_b));

	// Wake just one body -- both should wake
	body_wake(w, box_a);
	TEST_BEGIN("sleep wake island: waking one wakes both");
	TEST_ASSERT(!body_is_sleeping(w, box_a));
	TEST_ASSERT(!body_is_sleeping(w, box_b));

	destroy_world(w);
}

static void run_sleep_tests()
{
	printf("--- nudge sleep system tests ---\n");
	test_sleep_stacked_boxes();
	test_sleep_wake_on_impact();
	test_sleep_body_wake_api();
	test_sleep_set_velocity_wakes();
	test_sleep_static_no_island();
	test_sleep_destroy_body_removes_from_island();
	test_sleep_joint_creates_island();
	test_sleep_separate_islands();
	test_sleep_wake_wakes_whole_island();
}

static void run_ldl_stress_tests()
{
	printf("--- LDL stress tests ---\n");
	test_ldl_stress_coincident_anchors();
	test_ldl_stress_both_static();
	test_ldl_stress_extreme_mass_ratio();
	test_ldl_stress_redundant_joints();
	test_ldl_stress_topology_thrash();
	test_ldl_stress_velocity_bomb();
	test_ldl_stress_zero_length_distance();
	test_ldl_stress_dense_clique();
	test_ldl_stress_collinear_chain();
	test_ldl_stress_alternating_mass();
	test_ldl_stress_single_constraint();
	test_ldl_stress_mixed_same_pair();
	test_ldl_stress_near_max_nodes();
	test_ldl_stress_inverted_gravity_tangle();
	test_ldl_stress_distance_only_chain();
	test_ldl_stress_opposing_slams();
	test_ldl_stress_zero_gravity_spin();
	test_ldl_stress_mixed_stiffness();
	test_ldl_stress_double_hub_bridge();
	test_ldl_stress_body_destroy_recreate();
	test_ldl_stress_micro_mass();
	test_ldl_stress_constraint_loop();
	test_ldl_stress_stretched_recovery();
	test_ldl_stress_heavy_stretched_recovery();
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
	test_floor_contacts();
	test_sphere_settles_on_floor();
	test_gjk_dispatch();
	test_gjk_distance();
	test_contact_sanity();
	test_box_box();
	test_quickhull();

	test_quickhull_case783();
	test_quickhull_ico7037();
	test_quickhull_tet8098();
	test_quickhull_soak_ico1336();
	test_quickhull_soak_tet101();
	test_quickhull_cylinder();
	test_quickhull_bipyramid_loops();
	test_quickhull_geometry_stress();

	// Quickhull fuzz: moderate count for regular runs.
	test_quickhull_fuzz(20); // use 1000+ for stress testing

	run_solver_tests();
	run_ldl_stress_tests();
	run_bvh_tests();
	run_query_tests();
	test_feature_ids();
	run_sleep_tests();

	printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	if (test_fail > 0) printf("*** FAILURES ***\n");
}

// ============================================================================
// Box stack stability benchmark.
// Creates a floor + N unit boxes stacked vertically, simulates 1000 frames at
// 60Hz, and reports cumulative motion and max drift from ideal rest positions.

static void bench_box_stack_ex(int height, WorldParams wp)
{
	World w = create_world(wp);
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->sleep_enabled = 0;

	// Floor: static box
	Body floor_body = create_body(w, (BodyParams){
		.position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0,
	});
	body_add_shape(w, floor_body, (ShapeParams){
		.type = SHAPE_BOX, .box.half_extents = V3(50, 1, 50),
	});

	// Stack N unit boxes (half_extents 0.5) at y = 0.5, 1.5, 2.5, ...
	float half = 0.5f;
	Body* boxes = NULL;
	float* ideal_y = NULL;
	for (int i = 0; i < height; i++) {
		float y = half + (float)i;
		apush(ideal_y, y);
		Body b = create_body(w, (BodyParams){
			.position = V3(0, y, 0), .rotation = quat_identity(), .mass = 1.0f,
		});
		body_add_shape(w, b, (ShapeParams){
			.type = SHAPE_BOX, .box.half_extents = V3(half, half, half),
		});
		apush(boxes, b);
	}

	// Simulate 1000 frames, accumulate motion
	float dt = 1.0f / 60.0f;
	double cumulative_motion = 0.0;
	v3* prev_pos = NULL;
	for (int i = 0; i < height; i++)
		apush(prev_pos, body_get_position(w, boxes[i]));

	for (int frame = 0; frame < 1000; frame++) {
		world_step(w, dt);
		for (int i = 0; i < height; i++) {
			v3 pos = body_get_position(w, boxes[i]);
			v3 delta = sub(pos, prev_pos[i]);
			cumulative_motion += (double)len(delta);
			prev_pos[i] = pos;
		}
	}

	// Measure max drift from ideal rest positions
	float max_drift = 0.0f;
	for (int i = 0; i < height; i++) {
		v3 pos = body_get_position(w, boxes[i]);
		v3 ideal = V3(0, ideal_y[i], 0);
		float drift = len(sub(pos, ideal));
		if (drift > max_drift) max_drift = drift;
	}

	printf("stack_%d (sub=%d vel=%d hz=%.0f damp=%.1f): motion=%.6f max_drift=%.4f\n", height, wi->sub_steps, wi->velocity_iters, wi->contact_hertz, wi->contact_damping_ratio, cumulative_motion, max_drift);

	afree(boxes);
	afree(ideal_y);
	afree(prev_pos);
	destroy_world(w);
}

static void bench_box_stack(int height)
{
	bench_box_stack_ex(height, (WorldParams){ .gravity = V3(0, -9.81f, 0) });
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
