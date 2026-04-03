// nudge_internal.h -- internal types for nudge physics engine (unity build)
#ifndef NUDGE_INTERNAL_H
#define NUDGE_INTERNAL_H

// Portable float validity: rejects NaN, inf, and extreme magnitudes.
// NaN fails the self-equality test; inf and huge values fail the bound test.
static inline int float_valid(float f) { return f == f && f > -1e18f && f < 1e18f; }
static inline int v3_is_valid(v3 v) { return float_valid(v.x) && float_valid(v.y) && float_valid(v.z); }
static inline int quat_is_valid(quat q) { return float_valid(q.x) && float_valid(q.y) && float_valid(q.z) && float_valid(q.w); }
#define is_valid(x) _Generic((x), float: float_valid, v3: v3_is_valid, quat: quat_is_valid)(x)

// -----------------------------------------------------------------------------
// Internal types: cold (metadata) and hot (solver iteration) splits.

typedef struct ShapeInternal
{
	ShapeType type;
	v3 local_pos;
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule;
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
	};
} ShapeInternal;

// Cold: persistent metadata, topology, rarely touched during solve.
typedef struct BodyCold
{
	float mass;
	CK_DYNA ShapeInternal* shapes;
	int bvh_leaf;     // -1 if not in BVH
	int island_id;    // -1 = no island (static or unconnected)
	int island_prev;  // prev body in island body list, -1 = head
	int island_next;  // next body in island body list, -1 = tail
} BodyCold;

// Hot: solver working set, iterated every step, packed for cache.
typedef struct BodyHot
{
	v3 position;
	quat rotation;
	v3 velocity;
	v3 angular_velocity;
	float inv_mass;
	v3 inv_inertia_local; // diagonal of local-space inverse inertia tensor
	float friction;
	float restitution;
	float linear_damping;
	float angular_damping;
	float sleep_time; // accumulated seconds below velocity threshold
} BodyHot;

typedef struct WarmManifold WarmManifold; // forward decl for warm cache

// Joint persistent storage (handle-based, parallel arrays like bodies).
typedef enum JointType { JOINT_BALL_SOCKET, JOINT_DISTANCE } JointType;

typedef struct JointInternal
{
	JointType type;
	int body_a, body_b; // array indices (resolved from handles at creation)
	union {
		struct { v3 local_a, local_b; SpringParams spring; } ball_socket;
		struct { v3 local_a, local_b; float rest_length; SpringParams spring; } distance;
	};
	// Warm starting: accumulated impulses persisted across frames.
	union {
		v3 warm_lambda3;   // ball socket (3 DOF)
		float warm_lambda1; // distance (1 DOF)
	};
	// Island linked list fields.
	int island_id;    // -1 = none
	int island_prev;  // -1 = head
	int island_next;  // -1 = tail
} JointInternal;

typedef struct BVHTree BVHTree; // forward decl, defined in bvh.c

// Island: group of connected bodies that can sleep/wake together.
typedef struct Island
{
	int head_body, tail_body, body_count;
	int head_joint, tail_joint, joint_count;
	int constraint_remove_count;
	int awake; // 1 = awake, 0 = sleeping
} Island;

#define SLEEP_VEL_THRESHOLD   0.01f  // squared velocity magnitude threshold
#define SLEEP_TIME_THRESHOLD  0.5f   // seconds of stillness before sleep
#define LINEAR_SLOP           0.01f  // contact margin: keeps contacts alive near zero-penetration

typedef struct WorldInternal
{
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts;
	CK_MAP(WarmManifold) warm_cache;
	// Joints
	CK_DYNA JointInternal* joints;
	CK_DYNA uint32_t*      joint_gen;
	CK_DYNA int*           joint_free;
	// Broadphase
	BroadphaseType broadphase_type;
	BVHTree* bvh_static;
	BVHTree* bvh_dynamic;
	// Islands
	CK_DYNA Island*     islands;
	CK_DYNA uint32_t*   island_gen;
	CK_DYNA int*        island_free;
	CK_MAP(uint8_t)     prev_touching; // body_pair_key -> 1 for pairs touching last frame
	int sleep_enabled; // 1 = island sleep active (default)
	FrictionModel friction_model;
	SolverType solver_type;
	int velocity_iters;
	int position_iters;
	float contact_hertz;
	float contact_damping_ratio;
	float max_push_velocity;
	int sub_steps;
} WorldInternal;

// -----------------------------------------------------------------------------
// Solver types and constants.

#define SOLVER_VELOCITY_ITERS      10
#define SOLVER_POSITION_ITERS      4
#define SOLVER_BAUMGARTE           0.2f   // used by joints only; contacts use NGS
#define SOLVER_SLOP                0.005f // NGS position correction dead zone
#define SOLVER_RESTITUTION_THRESH  1.0f
#define SOLVER_POS_BAUMGARTE       0.2f
#define SOLVER_POS_MAX_CORRECTION  0.2f   // max position correction per step
#define SOLVER_MAX_LINEAR_VEL      500.0f
#define SOLVER_MAX_ANGULAR_VEL     100.0f

typedef struct SolverContact
{
	v3 r_a;              // contact offset from body A
	v3 r_b;              // contact offset from body B
	v3 normal;           // cached contact normal
	v3 tangent1;         // friction direction 1
	v3 tangent2;         // friction direction 2
	float eff_mass_n;    // 1/K for normal row
	float eff_mass_t1;   // 1/K for tangent1 row
	float eff_mass_t2;   // 1/K for tangent2 row
	float bias;          // velocity bias for penetration recovery
	float bounce;        // restitution velocity bias
	float lambda_n;      // accumulated normal impulse (>= 0)
	float lambda_t1;     // accumulated tangent1 impulse
	float lambda_t2;     // accumulated tangent2 impulse
	float softness;      // soft constraint regularization (0 = rigid/NGS)
	float penetration;   // cached for position correction pass
	uint32_t feature_id; // geometric feature key for warm starting
} SolverContact;

typedef struct SolverManifold
{
	int body_a, body_b;
	int contact_start;
	int contact_count;
	float friction;
	// Manifold-level patch friction data (FRICTION_PATCH only)
	v3 centroid_r_a;
	v3 centroid_r_b;
	v3 tangent1, tangent2;
	v3 normal;
	float eff_mass_t1, eff_mass_t2;
	float eff_mass_twist;
	float lambda_t1, lambda_t2;
	float lambda_twist;
	float patch_area;
	float patch_radius;
	// Block solver: normal coupling matrix A[i][j] (symmetric, max 4x4 = 10 unique)
	float K_nn[16];      // full NxN stored row-major (max 4x4)
	float K_nn_inv[16];  // inverse of K_nn
} SolverManifold;

// Warm starting: cached impulses from previous frame, keyed by body pair.
typedef struct WarmContact
{
	uint32_t feature_id; // geometric feature key for matching
	v3 r_a;              // body-A-relative position for spatial fallback matching
	float lambda_n;
	float lambda_t1;
	float lambda_t2;
} WarmContact;

struct WarmManifold
{
	WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale; // 0 = updated this frame, incremented each frame not touched, evicted at >1
	// Manifold-level friction warm data (FRICTION_PATCH)
	float manifold_lambda_t1;
	float manifold_lambda_t2;
	float manifold_lambda_twist;
};

// -----------------------------------------------------------------------------
// Joint solver types.

typedef struct SolverBallSocket
{
	int body_a, body_b;
	v3 r_a, r_b;
	float eff_mass[6]; // symmetric 3x3 inverse (xx,xy,xz,yy,yz,zz)
	v3 bias;
	float softness;
	v3 lambda;
	int joint_idx;
} SolverBallSocket;

typedef struct SolverDistance
{
	int body_a, body_b;
	v3 r_a, r_b;
	v3 axis;
	float eff_mass;
	float bias;
	float softness;
	float lambda;
	int joint_idx;
} SolverDistance;

// Constraint ref for graph coloring dispatch.
enum { CTYPE_CONTACT, CTYPE_BALL_SOCKET, CTYPE_DISTANCE };

typedef struct ConstraintRef
{
	uint8_t type;
	uint8_t color;
	int index;
	int body_a, body_b;
} ConstraintRef;

#endif // NUDGE_INTERNAL_H
