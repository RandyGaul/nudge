// See LICENSE for licensing info.
//
// Credits and references
// ----------------------
// BEPUphysics v2 -- broadphase, soft constraints, memory layout ideas
// Box2D -- GJK, sequential impulse solver reference
// Dirk Gregorius -- SAT with Gauss map pruning, contact clipping, direct enumeration LCP block solver
// Chris Giles -- Augmented Vertex Block Descent solver reference

#ifndef NUDGE_H
#define NUDGE_H

#include "vmath.h"
#include "split_store.h"

// -----------------------------------------------------------------------------
// Opaque handles.

typedef struct World { uint64_t id; } World;
typedef struct Body { uint64_t id; } Body;

// -----------------------------------------------------------------------------
// Shape types -- usable standalone for direct collision queries.

typedef struct Sphere
{
	v3 center;
	float radius;
} Sphere;

typedef struct Capsule
{
	v3 p, q;       // segment endpoints (world space)
	float radius;
} Capsule;

typedef struct Box
{
	v3 center;
	quat rotation;
	v3 half_extents;
} Box;

// Half-edge mesh for convex polyhedra.
// Edges stored in twin pairs: edge 2k and twin 2k+1.
typedef struct HalfEdge
{
	uint16_t twin;
	uint16_t next;
	uint16_t origin;
	uint16_t face;
} HalfEdge;

typedef struct HullPlane
{
	v3 normal;
	float offset;    // dot(normal, point_on_plane)
} HullPlane;

typedef struct HullFace
{
	uint16_t edge;    // first half-edge on this face
} HullFace;

// Convex hull -- half-edge mesh with precomputed face planes.
// Can point to static data (e.g. unit box) or dynamically built via quickhull.
typedef struct Hull
{
	v3 centroid;
	const v3*        verts;
	const HalfEdge*  edges;
	const HullFace*  faces;
	const HullPlane* planes;
	int vert_count;
	int edge_count;  // total half-edges (2x undirected edges)
	int face_count;
	float epsilon;      // build tolerance: 3*(max|x|+max|y|+max|z|)*FLT_EPSILON
	float maxoutside;   // max distance any vertex was widened beyond Newell plane
} Hull;

// Positioned hull for collision queries.
typedef struct ConvexHull
{
	const Hull* hull;
	v3 center;
	quat rotation;
	v3 scale;        // per-axis scale applied to unit hull verts
} ConvexHull;

// -----------------------------------------------------------------------------
// Contact manifold.

// Contact feature ID for warm starting.
// Encodes which geometric features (face/edge/vertex) produced the contact,
// enabling frame-to-frame matching by integer comparison (Box2D/BEPU style).
//
// Face contacts: ref_face | (inc_face << 8) | (clip_edge << 16)
//   ref_face:  index of reference face (the clipping face)
//   inc_face:  index of incident face
//   clip_edge: which side plane produced this vertex (0xFF = original vertex)
//
// Edge contacts: edge_a | (edge_b << 16) | FEATURE_EDGE_BIT

typedef struct Contact
{
	v3 point;
	v3 normal;         // from shape A toward shape B
	float penetration; // positive = overlapping
	uint32_t feature_id;
} Contact;

#define MAX_CONTACTS 4

typedef struct Manifold
{
	Contact contacts[MAX_CONTACTS];
	int count;
} Manifold;

// -----------------------------------------------------------------------------
// Collision queries.
//
// Each returns nonzero on hit. Pass manifold=NULL for boolean-only test.
// Normal convention: points from A toward B.

int collide_sphere_sphere(Sphere a, Sphere b, Manifold* manifold);
int collide_sphere_capsule(Sphere a, Capsule b, Manifold* manifold);
int collide_sphere_box(Sphere a, Box b, Manifold* manifold);
int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold);
int collide_capsule_capsule(Capsule a, Capsule b, Manifold* manifold);
int collide_capsule_box(Capsule a, Box b, Manifold* manifold);
int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold);
int collide_box_box(Box a, Box b, Manifold* manifold);
int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold);

// Built-in unit box hull (half-extents 1,1,1). Use with ConvexHull + scale for boxes.
const Hull* hull_unit_box();

// -----------------------------------------------------------------------------
// Quickhull -- build a convex hull from a point cloud.
//
// Returns a heap-allocated Hull. Caller frees with hull_free().
// The resulting hull has proper half-edge topology, face planes, and centroid.

Hull* quickhull(const v3* points, int count);
void hull_free(Hull* hull);

// -----------------------------------------------------------------------------
// Body params for world API.

typedef enum ShapeType
{
	SHAPE_SPHERE,
	SHAPE_CAPSULE,
	SHAPE_BOX,
	SHAPE_HULL,
} ShapeType;

typedef struct ShapeParams
{
	ShapeType type;
	v3 local_pos;   // offset from body origin
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule; // segment along local Y
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
	};
} ShapeParams;

typedef struct BodyParams
{
	v3 position;
	quat rotation;
	float mass;            // 0 = static/kinematic
	float friction;        // Coulomb mu (default 0.5)
	float restitution;     // bounce coefficient (default 0.0)
	float linear_damping;  // velocity decay coefficient (default 0.0)
	float angular_damping; // angular velocity decay coefficient (default 0.03)
} BodyParams;

typedef enum BroadphaseType { BROADPHASE_N2, BROADPHASE_BVH } BroadphaseType;

typedef enum FrictionModel
{
	FRICTION_COULOMB,  // per-point Coulomb (2 tangent rows per contact)
	FRICTION_PATCH,    // manifold-level 2D friction using patch area estimate
} FrictionModel;

typedef enum SolverType
{
	SOLVER_SOFT_STEP,  // soft contacts, relax each substep (default)
	SOLVER_SI_SOFT,    // soft contacts, no relax between substeps
	SOLVER_SI,         // hard constraints, NGS position correction
	SOLVER_BLOCK,      // direct LCP enumeration for 2-4 contact normals
	SOLVER_AVBD,       // augmented vertex block descent (primal-dual position solver)
} SolverType;

typedef struct WorldParams
{
	v3 gravity;
	BroadphaseType broadphase;
	FrictionModel friction_model;
	SolverType solver_type;
	int velocity_iters;  // 0 = default (10)
	int position_iters;  // 0 = default (4)
	float contact_hertz;          // 0 = default (30.0), soft contact frequency
	float contact_damping_ratio;  // 0 = default (10.0), heavily overdamped
	float max_push_velocity;      // 0 = default (3.0 m/s)
	int sub_steps;                // 0 = default (1)
} WorldParams;

// -----------------------------------------------------------------------------
// World API.

World create_world(WorldParams params);
void destroy_world(World world);
void world_step(World world, float dt);
void world_set_friction_model(World world, FrictionModel model);
void world_set_solver_type(World world, SolverType type);

Body create_body(World world, BodyParams params);
void destroy_body(World world, Body body);
void body_add_shape(World world, Body body, ShapeParams params);

v3 body_get_position(World world, Body body);
quat body_get_rotation(World world, Body body);

void body_wake(World world, Body body);
void body_set_velocity(World world, Body body, v3 vel);
void body_set_angular_velocity(World world, Body body, v3 avel);
int body_is_asleep(World world, Body body);

// Debug: contact points from last step. Returns count, *out valid until next step.
int world_get_contacts(World world, const Contact** out);

// -----------------------------------------------------------------------------
// World queries.

typedef struct RayHit
{
	Body body;
	v3 point;       // world-space hit point
	v3 normal;      // surface normal at hit (outward)
	float distance; // distance along ray
} RayHit;

// Find all bodies whose AABB overlaps the query box.
// Writes up to max_results body handles into results[]. Returns total hit count
// (may exceed max_results -- use to size a retry).
int world_query_aabb(World world, v3 lo, v3 hi, Body* results, int max_results);

// Cast a ray and find the closest body hit. Direction is normalized internally.
// Returns nonzero on hit; fills *hit (may be NULL for boolean-only test).
int world_raycast(World world, v3 origin, v3 direction, float max_distance, RayHit* hit);

// -----------------------------------------------------------------------------
// Joints.

typedef struct Joint { uint64_t id; } Joint;

typedef struct SpringParams
{
	float frequency;      // Hz (0 = rigid constraint)
	float damping_ratio;  // 1.0 = critically damped
} SpringParams;

typedef struct BallSocketParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	SpringParams spring;  // {0,0} = rigid
} BallSocketParams;

typedef struct DistanceParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	float rest_length;    // 0 = compute from initial body positions
	SpringParams spring;
} DistanceParams;

typedef struct HingeParams
{
	Body body_a, body_b;
	v3 local_offset_a;    // anchor in body A local space
	v3 local_offset_b;    // anchor in body B local space
	v3 local_axis_a;      // hinge axis in body A local space
	v3 local_axis_b;      // hinge axis in body B local space
	SpringParams spring;  // {0,0} = rigid
} HingeParams;

Joint create_ball_socket(World world, BallSocketParams params);
Joint create_distance(World world, DistanceParams params);
Joint create_hinge(World world, HingeParams params);
void destroy_joint(World world, Joint joint);

// Debug: iterate BVH nodes. Calls fn(min, max, depth, is_leaf, user) for each node child.
typedef void (*BVHDebugFn)(v3 min, v3 max, int depth, int is_leaf, void* user);
void world_debug_bvh(World world, BVHDebugFn fn, void* user);

// Debug: iterate all joints. Calls fn with world-space anchor points.
typedef struct JointDebugInfo
{
	int type;            // 0=ball_socket, 1=distance, 2=hinge
	v3 anchor_a;     // world-space anchor on body A
	v3 anchor_b;     // world-space anchor on body B
	v3 axis_a;       // hinge axis in world space (hinge only)
	int is_soft;     // 1 if spring.frequency > 0
} JointDebugInfo;
typedef void (*JointDebugFn)(JointDebugInfo info, void* user);
void world_debug_joints(World world, JointDebugFn fn, void* user);

#endif
