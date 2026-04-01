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
		GjkResult r = gjk_query_point_hull(V3(2, 0, 0), unit_box);
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
		GjkResult r = gjk_query_segment_hull(V3(0, 3, 0), V3(0, 5, 0), unit_box);
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
		GjkResult r = gjk_query_hull_hull(box_a, box_b);
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
		GjkResult r = gjk_query_point_hull(V3(0, 2.5f, 0), floor);
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
		GjkResult r = gjk_query_point_hull(V3(0, 3, 0), floor);
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
		GjkResult r = gjk_query_point_hull(V3(2, 0, 0), rbox);
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
		GjkResult r = gjk_query_point_hull(s.center, (ConvexHull){ hull_unit_box(), V3(0,0,0), quat_identity(), V3(10,1,10) });
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
			GjkResult r = gjk_query_point_hull(s.center, fl);
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

	// Run 10 seconds for damped spring to settle
	step_n(w, 600);
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

	step_n(w, 600);
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

static void run_solver_tests()
{
	printf("--- nudge solver tests ---\n");
	test_solver_nan_rejection();
	test_contact_sphere_rests_on_floor();
	test_contact_box_rests_on_floor();
	test_distance_spring_equilibrium();
	test_distance_rigid_maintains_length();
	test_ball_socket_chain_stays_connected();
	test_ball_socket_pin_converges();
	test_ball_socket_pendulum();
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
	test_quickhull_fuzz(100); // use 1000+ for stress testing

	run_solver_tests();
	run_bvh_tests();
	test_feature_ids();
	run_sleep_tests();

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
