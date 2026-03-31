// See LICENSE for licensing info.
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

typedef struct Contact
{
	v3 point;
	v3 normal;         // from shape A toward shape B
	float penetration; // positive = overlapping
} Contact;

#define MAX_CONTACTS 5

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
	float mass;         // 0 = static/kinematic
	float friction;     // Coulomb mu (default 0.5)
	float restitution;  // bounce coefficient (default 0.0)
} BodyParams;

typedef struct WorldParams
{
	v3 gravity;
} WorldParams;

// -----------------------------------------------------------------------------
// World API.

World create_world(WorldParams params);
void destroy_world(World world);
void world_step(World world, float dt);

Body create_body(World world, BodyParams params);
void destroy_body(World world, Body body);
void body_add_shape(World world, Body body, ShapeParams params);

v3 body_get_position(World world, Body body);
quat body_get_rotation(World world, Body body);

// Debug: contact points from last step. Returns count, *out valid until next step.
int world_get_contacts(World world, const Contact** out);

#endif
