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

#include <stdint.h>
#include <stddef.h>

#include "vmath.h"
#include "split_store.h"

// -----------------------------------------------------------------------------
// DLL export / import macro. Every public declaration below is prefixed with
// NUDGE_API so the same header works for:
//
//   - static library / unity build (default): NUDGE_API expands to nothing.
//   - Windows DLL build of nudge:    -DNUDGE_BUILD_DLL -> __declspec(dllexport)
//   - Windows consumer of the DLL:   -DNUDGE_USE_DLL   -> __declspec(dllimport)
//   - Unix shared library build:     -DNUDGE_BUILD_DLL -> default visibility
//
// Override by defining NUDGE_API before including this header. Example:
//   #define NUDGE_API __attribute__((visibility("default")))
//   #include "nudge.h"

#ifndef NUDGE_API
	#if defined(_WIN32)
		#if defined(NUDGE_BUILD_DLL)
			#define NUDGE_API __declspec(dllexport)
		#elif defined(NUDGE_USE_DLL)
			#define NUDGE_API __declspec(dllimport)
		#else
			#define NUDGE_API
		#endif
	#else
		#if defined(NUDGE_BUILD_DLL)
			#define NUDGE_API __attribute__((visibility("default")))
		#else
			#define NUDGE_API
		#endif
	#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

// -----------------------------------------------------------------------------
// Performance timers (seconds per phase from last world_step).

typedef struct PGSTimers
{
	double pre_solve;      // build solver manifolds/contacts/joints
	double warm_start;     // apply warm start impulses
	double graph_color;    // constraint graph coloring
	double iterations;     // PGS velocity iterations (contacts + joints)
	double joint_limits;   // joint limit solving within PGS loop
	double ldl;            // LDL factor + velocity/position correction
	double relax;          // contact relaxation (soft step)
	double pos_contacts;   // NGS position correction for contacts
	double pos_joints;     // NGS position correction for joints
	double post_solve;     // warm cache update
} PGSTimers;

typedef struct PerfTimers
{
	double broadphase;
	double pre_solve;
	double pgs_solve;
	double position_correct;
	double integrate;
	double islands;
	double soft_body;      // soft-body LDL step (velocity solve + integrate + pin snap)
	double total;
	PGSTimers pgs;
} PerfTimers;

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

typedef struct Cylinder
{
	v3 center;
	quat rotation;
	float half_height; // along local Y
	float radius;
} Cylinder;

// Half-edge mesh for convex polyhedra.
// Edges stored in twin pairs: edge 2k and twin 2k+1.
typedef struct HalfEdge
{
	int twin;
	int next;
	int origin;
	int face;
} HalfEdge;

typedef struct HullPlane
{
	v3 normal;
	float offset;    // dot(normal, point_on_plane)
} HullPlane;

typedef struct HullFace
{
	int edge;    // first half-edge on this face
} HullFace;

// Convex hull -- half-edge mesh with precomputed face planes.
// Can point to static data (e.g. unit box) or dynamically built via quickhull.
//
// The `name` field (NULL by default) identifies the hull in snapshots. Tag a
// hull you want to save with `hull_set_name`, then `world_register_hull` so
// the loading world can resolve the name back to a live Hull* you provide.
typedef struct Hull
{
	v3 centroid;
	const v3*        verts;
	const float*     soa_verts; // SoA: x[0..n-1], y[0..n-1], z[0..n-1], 16-byte aligned
	const int*       edge_twin;   // SoA half-edge arrays
	const int*       edge_next;
	const int*       edge_origin;
	const int*       edge_face;
	const HullFace*  faces;
	const HullPlane* planes;
	int vert_count;
	int edge_count;  // total half-edges (2x undirected edges)
	int face_count;
	float epsilon;      // build tolerance: 3*(max|x|+max|y|+max|z|)*FLT_EPSILON
	float maxoutside;   // max distance any vertex was widened beyond Newell plane
	const char* name;   // stable sinterned tag for snapshots; NULL = not named
} Hull;

// Tag a hull with a name so it can appear in snapshot files. Interns the
// string internally; pass the same name on load to resolve it back.
NUDGE_API void hull_set_name(Hull* hull, const char* name);
NUDGE_API const char* hull_get_name(const Hull* hull);

// Compact hull variant -- uint8_t indices, single heap block, no SoA or HullFace.
// Caps: vert_count, edge_count, face_count all <= 256.
// Twin is implicit: twin(i) = i ^ 1 (edges stored in pair-adjacent slots).
// Built via hull_to_hull8(); returns NULL if the source exceeds any cap.
typedef struct Hull8
{
	v3 centroid;
	const v3*        verts;       // [vert_count]
	const uint8_t*   edge_next;   // [edge_count] (twin implicit: i ^ 1)
	const uint8_t*   edge_origin; // [edge_count]
	const uint8_t*   edge_face;   // [edge_count]
	const uint8_t*   face_edge;   // [face_count] first half-edge per face
	const HullPlane* planes;      // [face_count]
	int vert_count;
	int edge_count;  // always even (twin pairs)
	int face_count;
} Hull8;

// Positioned hull for collision queries.
typedef struct ConvexHull
{
	const Hull* hull;
	v3 center;
	quat rotation;
	v3 scale;        // per-axis scale applied to unit hull verts
} ConvexHull;

// -----------------------------------------------------------------------------
// Triangle mesh -- static (indexed) collision geometry with smooth-normal
// (ghost-edge) handling. Box2D's chain shape is the 2D analog.
//
// Input contract (caller responsibility, not validated):
//   - Mesh is a manifold surface: every interior edge is shared by at most two
//     triangles. Non-manifold meshes (>2 triangles per edge) will abort during
//     build. Boundary edges (exactly one triangle) are allowed.
//   - Triangle winding is consistent (CCW when viewed from the intended
//     collision side). Face normal is taken as cross(v1-v0, v2-v0).
//   - No degenerate triangles (zero-area); build asserts.
//   - Vertex positions are in the mesh-local frame. No welding is performed;
//     shared vertices must share their index.
//   - indices[] has length 3 * tri_count.
typedef struct TriMesh TriMesh;

NUDGE_API TriMesh* trimesh_create(const v3* verts, int vert_count, const uint32_t* indices, int tri_count);
NUDGE_API void trimesh_free(TriMesh* mesh);

// Debug: number of triangles in the mesh.
NUDGE_API int trimesh_tri_count(const TriMesh* mesh);

// Tag the mesh with a name for snapshot identification. Same pattern as
// hull_set_name: caller pre-registers the named mesh with both the saving
// world and the loading world.
NUDGE_API void trimesh_set_name(TriMesh* mesh, const char* name);
NUDGE_API const char* trimesh_get_name(const TriMesh* mesh);

// -----------------------------------------------------------------------------
// Heightfield -- static regular-grid collision surface. Cheaper than trimesh
// for terrain because the topology is implicit: N*N vertices on an N-by-N
// grid with uniform cell_size spacing, two triangles per cell with a fixed
// diagonal.
//
// Local frame: grid sits in the local X-Z plane with Y as height. Vertex
// (i, j) lives at world_local (i*cell_size, heights[j*N + i], j*cell_size),
// with i ranging over [0, N) along X and j over [0, N) along Z. The grid
// covers [0, (N-1)*cell_size] on both X and Z. Place the body's position
// to translate the whole field.
//
// Each cell (i, j) for i in [0, N-1), j in [0, N-1) produces two triangles:
// a diagonal from (i, j) to (i+1, j+1) splits each cell in half. Both
// triangles wind CCW when viewed from +Y so face normals point up.
//
// Static only: attach to a mass=0 body. Per-cell material ids are optional
// (heightfield_set_material_ids) and feed ContactSummary.material_* through
// the same palette body shapes use.
//
// v1 limits: plain float storage (quantization is a future optimization),
// fixed diagonal, no holes / no-collision sentinel. N must satisfy
// (N-1)*(N-1)*2 <= 65535 (i.e. N <= 182); the manifold sub-id reserves 16
// bits for the triangle index.
typedef struct Heightfield Heightfield;

NUDGE_API Heightfield* heightfield_create(const float* heights, int N, float cell_size);
NUDGE_API void heightfield_free(Heightfield* hf);

// Debug: number of triangles (= 2 * (N-1) * (N-1)).
NUDGE_API int heightfield_tri_count(const Heightfield* hf);

// Optional per-cell material ids; ids[] length must equal (N-1)*(N-1).
// Pass ids=NULL to clear. Summaries on a heightfield side read from this
// table when present; otherwise fall back to the body default material_id.
NUDGE_API void heightfield_set_material_ids(Heightfield* hf, const uint8_t* ids);
NUDGE_API uint8_t heightfield_get_material_id(const Heightfield* hf, int cell_index);

// Name tagging for snapshot identification (same pattern as hulls/meshes).
NUDGE_API void heightfield_set_name(Heightfield* hf, const char* name);
NUDGE_API const char* heightfield_get_name(const Heightfield* hf);

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

NUDGE_API int collide_sphere_sphere(Sphere a, Sphere b, Manifold* manifold);
NUDGE_API int collide_sphere_capsule(Sphere a, Capsule b, Manifold* manifold);
NUDGE_API int collide_sphere_box(Sphere a, Box b, Manifold* manifold);
NUDGE_API int collide_sphere_hull(Sphere a, ConvexHull b, Manifold* manifold);
NUDGE_API int collide_capsule_capsule(Capsule a, Capsule b, Manifold* manifold);
NUDGE_API int collide_capsule_box(Capsule a, Box b, Manifold* manifold);
NUDGE_API int collide_capsule_hull(Capsule a, ConvexHull b, Manifold* manifold);
NUDGE_API int collide_box_box(Box a, Box b, Manifold* manifold);
NUDGE_API int collide_hull_hull(ConvexHull a, ConvexHull b, Manifold* manifold);

// Native cylinder collision. Cylinder is always shape A; for cyl-as-B callers
// swap arguments and flip manifold normals after calling.
NUDGE_API int collide_cylinder_sphere(Cylinder a, Sphere b, Manifold* manifold);
NUDGE_API int collide_cylinder_capsule(Cylinder a, Capsule b, Manifold* manifold);
NUDGE_API int collide_cylinder_box(Cylinder a, Box b, Manifold* manifold);
NUDGE_API int collide_cylinder_hull(Cylinder a, ConvexHull b, Manifold* manifold);
NUDGE_API int collide_cylinder_cylinder(Cylinder a, Cylinder b, Manifold* manifold);

// Built-in unit box hull (half-extents 1,1,1). Use with ConvexHull + scale for boxes.
NUDGE_API const Hull* hull_unit_box();

// -----------------------------------------------------------------------------
// Quickhull -- build a convex hull from a point cloud.
//
// Returns a heap-allocated Hull. Caller frees with hull_free().
// The resulting hull has proper half-edge topology, face planes, and centroid.

NUDGE_API Hull* quickhull(const v3* points, int count);
NUDGE_API void hull_free(Hull* hull);

// Compact Hull8 converter. Returns NULL if src exceeds any uint8 cap
// (V/E/F > 256) or if src has inconsistent twin topology. The resulting
// Hull8 is a single heap block; free with hull8_free().
NUDGE_API Hull8* hull_to_hull8(const Hull* src);
NUDGE_API void hull8_free(Hull8* h);

// -----------------------------------------------------------------------------
// Body params for world API.

typedef enum ShapeType
{
	SHAPE_SPHERE,
	SHAPE_CAPSULE,
	SHAPE_BOX,
	SHAPE_HULL,
	SHAPE_CYLINDER,
	SHAPE_MESH,        // static triangle mesh; allowed only on mass=0 bodies
	SHAPE_HEIGHTFIELD, // static regular-grid heightfield; mass=0 only
} ShapeType;

typedef struct ShapeParams
{
	ShapeType type;
	v3 local_pos;   // offset from body origin
	quat local_rot; // rotation in body local frame; zero-quat = identity (default)
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule; // segment along local Y
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
		struct { float half_height; float radius; } cylinder; // segment along local Y, flat caps
		struct { const TriMesh* mesh; } mesh;                 // static only
		struct { const Heightfield* hf; } heightfield;        // static only
	};
} ShapeParams;

typedef struct BodyParams
{
	v3 position;
	quat rotation;
	float mass;            // 0 = static/kinematic
	float friction;        // Coulomb mu (default 0.5)
	float rolling_friction;// rolling-resistance mu: 1 angular row per manifold,
	                       // axis = snapshot of tangent-plane ang-vel, capped
	                       // by (rolling_friction * lambda_n) (default 0.0 =
	                       // off; try 0.02-0.1 for balls/cylinders that should
	                       // decelerate while rolling)
	float restitution;     // bounce coefficient (default 0.0)
	float linear_damping;  // velocity decay coefficient (default 0.0)
	float angular_damping; // angular velocity decay coefficient (default 0.03)
} BodyParams;

typedef enum BroadphaseType { BROADPHASE_N2, BROADPHASE_BVH } BroadphaseType;

// Narrowphase backend choice. Default SAT produces one-shot 4-contact manifolds
// with feature-ID warm cache. GJK+EPA produces a single contact per frame and
// accumulates contacts via an incremental per-pair manifold (see EpaManifold).
typedef enum NarrowphaseBackend
{
	NARROWPHASE_SAT,
	NARROWPHASE_GJK_EPA,
} NarrowphaseBackend;

typedef enum SolverType
{
	SOLVER_SOFT_STEP,  // soft contacts, relax each substep (default)
	SOLVER_SI_SOFT,    // soft contacts, no relax between substeps
	SOLVER_SI,         // hard constraints, NGS position correction
} SolverType;

typedef struct WorldParams
{
	v3 gravity;
	BroadphaseType broadphase;
	NarrowphaseBackend narrowphase_backend; // default NARROWPHASE_SAT
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

NUDGE_API World create_world(WorldParams params);
NUDGE_API void destroy_world(World world);
NUDGE_API void world_step(World world, float dt);
NUDGE_API void world_set_solver_type(World world, SolverType type);

NUDGE_API Body create_body(World world, BodyParams params);
NUDGE_API void destroy_body(World world, Body body);
NUDGE_API void body_add_shape(World world, Body body, ShapeParams params);

NUDGE_API v3 body_get_position(World world, Body body);
NUDGE_API void body_set_position(World world, Body body, v3 pos);
NUDGE_API quat body_get_rotation(World world, Body body);
NUDGE_API v3 body_get_velocity(World world, Body body);
NUDGE_API v3 body_get_angular_velocity(World world, Body body);

NUDGE_API void body_wake(World world, Body body);
NUDGE_API void body_set_velocity(World world, Body body, v3 vel);
NUDGE_API void body_set_angular_velocity(World world, Body body, v3 avel);
NUDGE_API int body_is_asleep(World world, Body body);

// Sleep control.
NUDGE_API void world_set_sleep_enabled(World world, int enabled);
NUDGE_API int world_get_sleep_enabled(World world);
NUDGE_API void body_set_sleep_allowed(World world, Body body, int allowed);  // per-body override: 0 = never sleep

// Debug: contact points from last step. Returns count, *out valid until next step.
NUDGE_API int world_get_contacts(World world, const Contact** out);

// -----------------------------------------------------------------------------
// Contact summaries -- user-facing, one entry per colliding body pair.
//
// Populated after narrowphase each world_step. The underlying buffer is:
//   - sorted by (a.id, b.id) with canonical order a.id < b.id
//   - deduped -- multi-triangle mesh pairs or multi-shape bodies collapse
//     into a single summary per body pair
//   - stable across a step (valid until the next world_step)
//
// Each summary reduces the manifold to a single contact patch: one normal,
// the patch centroid, the radius of the patch (max dist from centroid to any
// contact), and the deepest penetration. Sub-indices (sub_a/sub_b) carry
// the tri index + 1 when the body is a trimesh (0 otherwise); the material
// fields carry the resolved per-shape material id (0 until wired).

typedef struct ContactSummary
{
	Body a, b;
	v3 normal;          // world-space, from a toward b
	v3 point;           // patch centroid (world space)
	float radius;       // max dist from centroid to any contact in the patch
	float depth;        // deepest penetration in the patch (>= 0)
	uint8_t material_a; // resolved material id on side a (palette index)
	uint8_t material_b; // resolved material id on side b
	uint16_t _pad;
	uint32_t sub_a;     // tri index + 1 if a is a trimesh, else 0
	uint32_t sub_b;     // tri index + 1 if b is a trimesh, else 0
} ContactSummary;

// Fetch the current frame's contact summaries. Returned pointer + count are
// valid until the next world_step. Writes count into *out_count.
NUDGE_API const ContactSummary* world_contact_summaries(World world, int* out_count);

// -----------------------------------------------------------------------------
// Material palette. 256 palette entries per world; bodies and trimesh triangles
// carry a uint8 material_id that indexes this palette. Engine does not read
// friction/restitution from the palette -- per-body BodyParams.friction/
// restitution still drives the solver. The palette is reporting-only: it lets
// users pass opaque metadata (user_data) or simple tuning floats through the
// ContactSummary API without maintaining a side-table keyed by body+tri.
//
// Defaults: palette[0] = {0.5, 0.0, 0}. Bodies default to material_id 0.
// TriMesh triangles default to body's material_id until trimesh_set_material_ids.

typedef struct Material
{
	float friction;
	float restitution;
	uint32_t user_data;
} Material;

NUDGE_API void world_set_material(World world, uint8_t id, Material m);
NUDGE_API Material world_get_material(World world, uint8_t id);

// Set the default material id carried by a body. Summary.material_a/b reports
// this id for non-mesh sides and for mesh sides without per-tri material ids.
NUDGE_API void body_set_material_id(World world, Body body, uint8_t id);
NUDGE_API uint8_t body_get_material_id(World world, Body body);

// Attach a per-triangle material-id array to a trimesh. Array length must
// match trimesh_tri_count(mesh). The trimesh copies the array internally;
// caller may free after return. Pass ids=NULL to clear per-tri materials
// (summary falls back to the body's default material_id).
NUDGE_API void trimesh_set_material_ids(TriMesh* mesh, const uint8_t* ids);
// Read one triangle's material id. Returns 0 if no per-tri array is set.
NUDGE_API uint8_t trimesh_get_material_id(const TriMesh* mesh, int tri_index);

// Per-body contact listener. Fired once per body per step (only when the
// body has >= 1 contact that step). The pairs[] view is contiguous and
// normalized: `self` is always in the `a` slot -- when the underlying
// summary has `self == b`, the engine flips the normal and swaps the
// material/sub fields before the call so the callback always reads "from
// self's perspective." Pointer validity is only guaranteed for the
// duration of the call.
//
// Set fn=NULL to clear. Replacing overwrites the previous listener.
typedef void (*BodyContactListener)(Body self, const ContactSummary* pairs, int count, void* ud);
NUDGE_API void body_set_contact_listener(World world, Body body, BodyContactListener fn, void* ud);

// Per-frame EPA narrowphase stats. Reset at the top of each world_step.
// Only populated when narrowphase_backend == NARROWPHASE_GJK_EPA.
typedef struct WorldEpaStats
{
	int queries;
	int iter_cap_hits;
	int total_iters;
	int warm_reseeds;
	int contacts_emitted;
	int pair_count;
} WorldEpaStats;
NUDGE_API WorldEpaStats world_get_epa_stats(World world);

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
NUDGE_API int world_query_aabb(World world, v3 lo, v3 hi, Body* results, int max_results);

// Cast a ray and find the closest body hit. Direction is normalized internally.
// Returns nonzero on hit; fills *hit (may be NULL for boolean-only test).
NUDGE_API int world_raycast(World world, v3 origin, v3 direction, float max_distance, RayHit* hit);

// -----------------------------------------------------------------------------
// Rewind -- ring buffer of deterministic world snapshots.
//
// Captures everything the simulation reads: body hot/state/cold, joints,
// islands, warm-cache impulses, prev-touching + joint-pair maps, ldl topo
// version. Restoring + stepping forward is bit-identical to the original
// run from that point.
//
// Topology-change invalidation: any API call that creates or destroys a
// body / joint / shape flushes the ring buffer. Mutation replay is a v2.
//
// Not supported by v1:
//   - EPA backend (epa_cache is not captured; use the SAT backend).
//   - Mutating the world between capture and restore (flushes the buffer).

typedef struct RewindParams
{
	int max_frames;    // ring buffer depth (0 = disabled)
	int auto_capture;  // 1 = world_step captures at start automatically
} RewindParams;

NUDGE_API void world_rewind_init(World world, RewindParams params);
NUDGE_API void world_rewind_shutdown(World world);

// Manual capture. Returns monotonic frame_id, or 0 if rewind is disabled.
// Oldest frame is evicted when the buffer is full.
NUDGE_API uint64_t world_rewind_capture(World world);

// Restore world to a frame_id previously returned from capture.
// Returns 1 on success, 0 if the frame is no longer in the buffer.
// All frames captured after frame_id remain available for re-step and
// re-rewind; frames captured before are preserved as well.
NUDGE_API int world_rewind_to_frame(World world, uint64_t frame_id);

// Convenience: restore to the snapshot `n` steps before the newest one.
// n=0 restores the most recent snapshot. Returns 1 on success.
NUDGE_API int world_rewind_by_steps(World world, int n);

NUDGE_API int    world_rewind_frames_available(World world);
NUDGE_API size_t world_rewind_memory_used(World world);

// -----------------------------------------------------------------------------
// Snapshot save/load -- binary versioned world persistence.
//
// Save writes the current world's public state (world params + bodies +
// joints) to a file. Load reads the file, creates a new world, and
// recreates every body + joint. Backwards-compatible across engine
// versions via per-field version tags (see serialize.c "SV" framework).
//
// NOT captured by v1 snapshots (reconstructed on next world_step):
//   warm cache, island state, BVH, LDL factorization, rewind ring,
//   incremental narrowphase hints. Simulation after load may differ by a
//   few frames of warm-cache rebuild from an equivalent rewind-restore.
//
// Restrictions: SHAPE_HULL and SHAPE_MESH bodies are not supported in v1
// (their caller-owned Hull*/TriMesh* have no stable on-disk identity).
// Save asserts if the world contains either shape type.
//
// Returns 1 on success, 0 on I/O failure.
NUDGE_API int   world_save_snapshot(World world, const char* path);
// Returns a new World (like create_world) or (World){0} on failure.
// After a successful load, live body/joint indices match saved order, so
// world_get_bodies() + world_get_joints() return handles in the same order
// the caller supplied at save time. Use those to rebind gameplay-object
// references to the new handles.
//
// This variant asserts if the snapshot contains SHAPE_HULL or SHAPE_MESH
// bodies (no asset registry to resolve names). Use world_load_snapshot_into
// when the file references named hulls or meshes.
NUDGE_API World world_load_snapshot(const char* path);

// Load a snapshot into a pre-configured world. The world's registered hulls
// and meshes (via world_register_hull / world_register_mesh) are used to
// resolve SHAPE_HULL / SHAPE_MESH shape names in the file. The world should
// be empty (no existing bodies/joints/sensors) -- asserts otherwise. The
// file's WorldParams are ignored; the existing world's config is used.
NUDGE_API int world_load_snapshot_into(World world, const char* path);

// Asset registration. Attach a named hull/mesh to a world so snapshots can
// reference it by name. The Hull* or TriMesh* must have a non-NULL name set
// via hull_set_name / trimesh_set_name before registering.
NUDGE_API void world_register_hull(World world, const Hull* hull);
NUDGE_API void world_register_mesh(World world, const TriMesh* mesh);
NUDGE_API void world_register_heightfield(World world, const Heightfield* hf);

// Lookup by name (returns NULL if not registered). Lets callers keep their
// own "asset id -> Hull*" maps out of sight.
NUDGE_API const Hull*         world_find_hull(World world, const char* name);
NUDGE_API const TriMesh*      world_find_mesh(World world, const char* name);
NUDGE_API const Heightfield*  world_find_heightfield(World world, const char* name);

// Iterate all live bodies / joints in the world. Writes up to max handles
// into out[], returns the total count (may exceed max -- use to size a retry).
// Order is live-index ascending, which matches creation order in a fresh
// world (e.g. immediately after world_load_snapshot).
NUDGE_API int world_get_body_count(World world);
NUDGE_API int world_get_bodies(World world, Body* out, int max);
// (world_get_joint_count / world_get_joints declared below, after Joint.)

// -----------------------------------------------------------------------------
// Sensors -- read-only world-query volumes, owned by the world.
//
// A sensor is a compound of convex shapes with a transform. It never touches
// the solver, never enters the broadphase, never generates contacts. The
// only operation on a sensor is sensor_query, which walks the world and
// reports which bodies overlap the sensor volume.
//
// Sensors are world-owned so snapshot save/load and rewind persist them
// alongside bodies and joints. Use world_get_sensors() to recover handles
// after a rewind or load.
//
// Thread safety: sensor_query is read-only against the world. Concurrent
// queries on different sensors are safe as long as no thread is mutating
// the world. Do not call during world_step.
//
// Not supported: body shapes of type SHAPE_MESH are skipped during sensor
// queries (triangle mesh vs. convex overlap is not yet wired up).

typedef struct Sensor { uint64_t id; } Sensor;

typedef struct SensorParams
{
	v3 position;
	quat rotation;             // zero-quat = identity
	uint32_t collision_group;  // 0 = 0xFFFFFFFF (overlap everything)
	uint32_t collision_mask;   // 0 = 0xFFFFFFFF
} SensorParams;

NUDGE_API Sensor create_sensor(World world, SensorParams params);
NUDGE_API void destroy_sensor(World world, Sensor sensor);

// Attach a shape to the sensor. SHAPE_MESH is not allowed; asserts.
NUDGE_API void sensor_add_shape(World world, Sensor sensor, ShapeParams params);

// Move the sensor. No broadphase update (sensor is not in the broadphase).
NUDGE_API void sensor_set_transform(World world, Sensor sensor, v3 position, quat rotation);

// Query the world for bodies overlapping the sensor. Writes up to max_results
// body handles into results[]; returns the total overlap count (may exceed
// max_results -- use that to size a retry). Read-only against world.
NUDGE_API int sensor_query(World world, Sensor sensor, Body* results, int max_results);

NUDGE_API int world_get_sensor_count(World world);
NUDGE_API int world_get_sensors(World world, Sensor* out, int max);

// -----------------------------------------------------------------------------
// Handle revalidation.
//
// After a rewind or snapshot-load, any handle the caller held becomes either:
//   - still valid  (the underlying object existed at the rewound-to frame and
//                   its generation counter was restored to match)
//   - stale        (the object was created after the snapshot, so the slot's
//                   generation counter is now different from what the caller has)
// Use these predicates before dereferencing a held handle across a rewind or
// load boundary. Use world_get_bodies / _joints / _sensors to enumerate what
// the world currently contains.
NUDGE_API int body_is_valid(World world, Body body);
NUDGE_API int sensor_is_valid(World world, Sensor sensor);
// joint_is_valid declared below, after Joint.

// -----------------------------------------------------------------------------
// Joints.

typedef struct Joint { uint64_t id; } Joint;

NUDGE_API int world_get_joint_count(World world);
NUDGE_API int world_get_joints(World world, Joint* out, int max);
NUDGE_API int joint_is_valid(World world, Joint joint);

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

typedef struct FixedParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	SpringParams spring;  // {0,0} = rigid weld
} FixedParams;

typedef struct PrismaticParams
{
	Body body_a, body_b;
	v3 local_offset_a;
	v3 local_offset_b;
	v3 local_axis_a;      // slide axis in body A local space
	v3 local_axis_b;      // slide axis in body B local space
	SpringParams spring;
} PrismaticParams;

// Angular motor: drives relative angular velocity along a shared axis toward
// target_speed, bounded by max_impulse. No linear constraint, no position
// correction -- purely velocity-level actuation. Useful for powered doors,
// turrets, etc.
typedef struct AngularMotorParams
{
	Body body_a, body_b;
	v3 local_axis_a;      // motor axis in body A local space
	v3 local_axis_b;      // motor axis in body B local space (usually parallel to A's)
	float target_speed;   // relative angular velocity to drive toward, rad/s
	float max_impulse;    // |lambda| clamp per step
} AngularMotorParams;

// Twist limit: clamps the relative rotation about a shared axis to [min, max]
// radians, measured via quaternion swing-twist decomposition. Unilateral at
// whichever boundary is violated. No linear constraint.
typedef struct TwistLimitParams
{
	Body body_a, body_b;
	v3 local_axis_a;      // twist axis in body A local space
	v3 local_axis_b;      // twist axis in body B local space (parallel at zero twist)
	float limit_min;      // radians; must be <= 0
	float limit_max;      // radians; must be >= 0
	SpringParams spring;  // softness for the limit constraint
} TwistLimitParams;

// Cone limit: clamps the angle between two body-local axes to <= half_angle.
// Unilateral (lambda >= 0, can only push). No linear constraint.
typedef struct ConeLimitParams
{
	Body body_a, body_b;
	v3 local_axis_a;      // cone reference axis in body A local space
	v3 local_axis_b;      // cone reference axis in body B local space
	float half_angle;     // radians; max deviation between axes
	SpringParams spring;  // softness for the limit constraint
} ConeLimitParams;

// Swing-twist: ball socket (3-DOF linear) + cone limit (swing) + twist limit.
// The composite ragdoll joint: anchor two bodies at a pair of local points,
// limit how far body B's twist axis can deviate from body A's (swing cone),
// and limit twist rotation about that axis to [twist_min, twist_max].
typedef struct SwingTwistParams
{
	Body body_a, body_b;
	v3 local_offset_a;    // anchor in body A local space
	v3 local_offset_b;    // anchor in body B local space
	v3 local_axis_a;      // twist axis in body A local space
	v3 local_axis_b;      // twist axis in body B local space (parallel at rest)
	float cone_half_angle; // radians; max swing of axis_b from axis_a
	float twist_min;      // radians; min twist about the axis
	float twist_max;      // radians; max twist about the axis
	SpringParams spring;  // softness for anchor + limits
} SwingTwistParams;

NUDGE_API Joint create_ball_socket(World world, BallSocketParams params);
NUDGE_API Joint create_distance(World world, DistanceParams params);
NUDGE_API Joint create_hinge(World world, HingeParams params);
NUDGE_API Joint create_fixed(World world, FixedParams params);
NUDGE_API Joint create_prismatic(World world, PrismaticParams params);
NUDGE_API void body_set_collision_filter(World world, Body body, uint32_t group, uint32_t mask);
NUDGE_API void body_set_compound_id(World world, Body body, uint32_t compound_id);

NUDGE_API Joint create_angular_motor(World world, AngularMotorParams params);
NUDGE_API Joint create_twist_limit(World world, TwistLimitParams params);
NUDGE_API Joint create_cone_limit(World world, ConeLimitParams params);
NUDGE_API Joint create_swing_twist(World world, SwingTwistParams params);
NUDGE_API void destroy_joint(World world, Joint joint);
NUDGE_API void joint_set_hinge_limits(World world, Joint joint, float min_angle, float max_angle);
NUDGE_API void joint_set_distance_limits(World world, Joint joint, float min_distance, float max_distance);
NUDGE_API void joint_clear_limits(World world, Joint joint);
NUDGE_API void joint_set_hinge_motor(World world, Joint joint, float speed, float max_impulse);
NUDGE_API void joint_set_prismatic_motor(World world, Joint joint, float speed, float max_impulse);

// Performance: phase timings from the last world_step call.
NUDGE_API PerfTimers world_get_perf(World world);

// Debug: iterate BVH nodes. Calls fn(min, max, depth, is_leaf, user) for each node child.
typedef void (*BVHDebugFn)(v3 min, v3 max, int depth, int is_leaf, void* user);
NUDGE_API void world_debug_bvh(World world, BVHDebugFn fn, void* user);

// Debug: iterate all joints. Calls fn with world-space anchor points.
typedef struct JointDebugInfo
{
	int type;                // JOINT_BALL_SOCKET..JOINT_PRISMATIC
	v3 anchor_a, anchor_b;  // world-space anchor on body A/B
	v3 axis_a;              // hinge/prismatic axis in world space
	int is_soft;            // 1 if spring.frequency > 0
	float motor_speed;      // target speed (0 = no motor)
	float motor_max_impulse; // 0 = no motor
	float limit_min, limit_max; // 0/0 = no limits
	float current_angle;    // current hinge angle
	int limit_active;       // 1 if angle is at a limit boundary
	v3 ref_a, ref_b;        // hinge reference directions (world space, for arc drawing)
} JointDebugInfo;
typedef void (*JointDebugFn)(JointDebugInfo info, void* user);
NUDGE_API void world_debug_joints(World world, JointDebugFn fn, void* user);

// -----------------------------------------------------------------------------
// Soft bodies -- native particle networks solved with direct LDL.
//
// A SoftBody is a bag of 3-DOF particles connected by scalar distance links
// (XPBD-style compliance). Unlike chaining rigid bodies together, particles
// are first-class solver entities with their own storage: no Body handles,
// no per-particle inertia tensors, no rotations.
//
// Topology is baked at build time. The solver caches a dense LDL factorization
// of K = J M^-1 J^T per soft body and refactors each substep (lever arms
// change with positions). For the scale this targets (cobwebs, ropes, jelly
// blobs of up to a few hundred nodes), dense LDL is fast and has zero
// symbolic-analysis cost.
//
// Sweet spot: small-to-medium gameplay objects where topology is constant
// and stiffness matters -- rope, cobweb, cloth patch, squishy creature.
//
// Not for: tearing, runtime topology changes, thousands of particles,
// volumetric FEM with hyperelastic materials.
//
// v1 scope: internal dynamics + pin-to-static. No collision, no rigid-body
// coupling. Those land in follow-up phases once the solver path is proven.

typedef struct SoftBody { uint64_t id; } SoftBody;

typedef struct SoftBodyParams
{
	SpringParams default_spring; // compliance for all links in this soft body (0,0 = rigid)
	float node_radius;           // per-node sphere radius for collision (0 = no collision)
	float linear_damping;        // per-substep exponential decay (default 0.02)
	int iterations;              // PGS velocity iterations per substep (default 12)
	uint32_t collision_group;    // reserved for future collision filters
	uint32_t collision_mask;     // reserved for future collision filters
	uint8_t material_id;         // reserved for future contact summaries
} SoftBodyParams;

// Lifecycle. Nodes and links are added between create and build.
NUDGE_API SoftBody create_soft_body(World world, SoftBodyParams params);
NUDGE_API void destroy_soft_body(World world, SoftBody sb);
NUDGE_API int soft_body_is_valid(World world, SoftBody sb);

// Assembly. All add_* must happen before soft_body_build(). Returns the node
// index (0..n-1, local to this soft body). mass<=0 is treated as 0 = pinned.
NUDGE_API int  soft_body_add_node(World world, SoftBody sb, v3 pos, float mass);
// rest_length<0 = compute from current node positions.
NUDGE_API void soft_body_add_link(World world, SoftBody sb, int node_i, int node_j, float rest_length, SpringParams spring);

// Freeze topology. After this, subsequent add_* calls assert. The dense K
// matrix and lambda warm cache are allocated here.
NUDGE_API void soft_body_build(World world, SoftBody sb);

// Pinning. Static pin anchors a node to a fixed world position (the node's
// position is snapped each step). Safe to call after build.
NUDGE_API void soft_body_pin_static(World world, SoftBody sb, int node, v3 world_pos);
NUDGE_API void soft_body_unpin(World world, SoftBody sb, int node);

// Readback for rendering / gameplay. Pointers are engine-owned; valid until
// the next world_step or a topology change on this soft body.
NUDGE_API int           soft_body_node_count(World world, SoftBody sb);
NUDGE_API const v3*     soft_body_node_positions(World world, SoftBody sb);
NUDGE_API const v3*     soft_body_node_velocities(World world, SoftBody sb);
NUDGE_API int           soft_body_link_count(World world, SoftBody sb);
// Fills out[0..count-1] with {i, j} index pairs for each link. Return value
// is total link count (may exceed max; use to size a retry).
NUDGE_API int           soft_body_get_links(World world, SoftBody sb, int* out_pairs, int max);

// Gameplay pokes.
NUDGE_API void soft_body_apply_force(World world, SoftBody sb, int node, v3 force);
NUDGE_API void soft_body_apply_impulse(World world, SoftBody sb, int node, v3 impulse);

// Enumeration.
NUDGE_API int world_get_soft_body_count(World world);
NUDGE_API int world_get_soft_bodies(World world, SoftBody* out, int max);

#ifdef __cplusplus
}
#endif

#endif
