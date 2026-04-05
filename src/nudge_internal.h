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
typedef struct AVBD_WarmManifold AVBD_WarmManifold;

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
	// AVBD warm state (penalty/lambda for augmented Lagrangian)
	v3 avbd_penalty_lin;
	v3 avbd_lambda_lin;
	v3 avbd_C0_lin;        // constraint error at x- for stabilization
	// Island linked list fields.
	int island_id;    // -1 = none
	int island_prev;  // -1 = head
	int island_next;  // -1 = tail
} JointInternal;

typedef struct BVHTree BVHTree; // forward decl, defined in bvh.c

// LDL direct solver types (used by solver_ldl.c).
// Per-DOF Jacobian row: J_a and J_b each have 6 components [lin_x, lin_y, lin_z, ang_x, ang_y, ang_z].
// LDL computes K = J * M^{-1} * J^T, RHS = -J * v, and impulse = M^{-1} * J^T * lambda generically.
typedef struct LDL_JacobianRow
{
	float J_a[6];
	float J_b[6];
} LDL_JacobianRow;

typedef struct LDL_Constraint
{
	int type;         // JOINT_BALL_SOCKET or JOINT_DISTANCE (or -1 for synthetic)
	int dof;          // 3 or 1
	int body_a, body_b;         // body indices (may be virtual shard indices for graph topology)
	int real_body_a, real_body_b; // always real body indices (for physics data lookup)
	float weight_a, weight_b;   // shard weight (1.0 = normal, S = shattered into S shards)
	int solver_idx;   // index into sol_bs[] or sol_dist[] (-1 for synthetic)
	int is_synthetic;
	int bundle_idx;   // which bundle this constraint belongs to
	int bundle_offset; // local DOF offset within its bundle's DOF space
	int jacobian_start; // index into LDL_Cache.jacobians[] for first row
} LDL_Constraint;

// Group of constraints between the same body pair, forming one graph node.
typedef struct LDL_Bundle
{
	int body_a, body_b;     // body pair
	int dof;                // total DOF (sum of constituent constraints)
	int start, count;       // range into LDL_Cache.constraints[]
} LDL_Bundle;

// Sparse block matrix: adjacency-list format per constraint node.
// Each node i has a diagonal block (dof[i] x dof[i]) and off-diagonal
// blocks to its neighbors. Stored as flat float arrays.
#define LDL_MAX_NODES 128

typedef struct LDL_Sparse
{
	int node_count;
	int dof[LDL_MAX_NODES];          // DOF per node (3 or 1)
	int row_offset[LDL_MAX_NODES+1]; // cumulative DOF offset
	int n;                            // total scalar DOFs

	// Diagonal blocks: packed lower-triangular, dof*(dof+1)/2 floats per node.
	float diag_data[LDL_MAX_NODES][6]; // max 3x3 packed (6 elements)
	float diag_D[LDL_MAX_NODES][3];   // D pivots from block LDL (max 3 per 3x3 block)

	// Off-diagonal: per-node neighbor list.
	// adj[i][k] = neighbor node index, adj_data[i][k] = block data (dof[i]*dof[j] floats).
	CK_DYNA int* adj[LDL_MAX_NODES];       // neighbor indices
	CK_DYNA float* adj_data[LDL_MAX_NODES]; // packed block floats per neighbor

	// Elimination ordering (symbolic phase output).
	int elim_order[LDL_MAX_NODES]; // pivot sequence
	int inv_order[LDL_MAX_NODES];  // inv_order[node] = elimination step
} LDL_Sparse;

// Precomputed neighbor for forward/back substitution in solve phase.
typedef struct LDL_Neighbor
{
	int node;       // neighbor node index
	int dn;         // neighbor DOF
	int on;         // neighbor row_offset
	int L_offset;   // offset into L_factors for the L block
} LDL_Neighbor;

// Schur complement update during numeric factorization.
typedef struct LDL_Schur
{
	int i, j;           // target nodes being updated
	int Lik_offset;     // L_factors offset for L_{i,k}
	int Ljk_offset;     // L_factors offset for L_{j,k}
	int target_offset;  // L_factors offset for target edge (i,j) (-1 = diagonal)
	int target_offset_rev; // L_factors offset for reverse edge (j,i) (-1 if diagonal or same)
	int target_node;    // diagonal node index when target_offset == -1
	int di, dk, dj;     // block dimensions
} LDL_Schur;

// Non-eliminated neighbor for L-block computation during factorization.
typedef struct LDL_Column
{
	int node;       // non-eliminated neighbor
	int dn;         // its DOF
	int L_offset;   // where to write L_{node,pivot} in L_factors
} LDL_Column;

// Per-pivot precomputed index ranges.
typedef struct LDL_Pivot
{
	int node, dk, ok;                      // node index, DOF, row_offset
	int fwd_start, fwd_count;             // forward-sub neighbors (eliminated before)
	int back_start, back_count;           // back-sub neighbors (eliminated after)
	int col_start, col_count;             // L-column entries for factorization
	int schur_start, schur_count;         // Schur complement updates for this pivot
} LDL_Pivot;

// K-edge fill instruction: maps shared body to off-diagonal L_factors offsets.
typedef struct LDL_Coupling
{
	int body;                              // shared body index (real or virtual)
	int block_i, block_j;                 // constraint node indices
	int L_offset_ij, L_offset_ji;         // L_factors offsets for both directions
} LDL_Coupling;

// Cached constraint system topology: computed once per topology change,
// reused across substeps. Stores elimination ordering, sparsity pattern,
// and precomputed memory offsets into L_factors.
typedef struct LDL_Topology
{
	int node_count, n;
	int dof[LDL_MAX_NODES];
	int row_offset[LDL_MAX_NODES + 1];
	int elim_order[LDL_MAX_NODES];
	int inv_order[LDL_MAX_NODES];
	LDL_Pivot pivots[LDL_MAX_NODES];
	CK_DYNA LDL_Neighbor* fwd_neighbors;
	CK_DYNA LDL_Neighbor* back_neighbors;
	CK_DYNA LDL_Column* columns;
	CK_DYNA LDL_Schur* schurs;
	CK_DYNA LDL_Coupling* couplings;
	int L_factors_size;                    // total floats in L_factors
} LDL_Topology;

typedef struct LDL_Cache
{
	CK_DYNA LDL_Constraint* constraints;
	int joint_count;
	CK_DYNA LDL_Bundle* bundles;
	int bundle_count;
	int n;              // total DOFs
	LDL_Topology* topo; // cached topology (NULL until built)
	CK_DYNA float* L_factors; // contiguous off-diagonal L-factor blocks
	CK_DYNA LDL_JacobianRow* jacobians; // per-DOF Jacobian rows (filled each substep)
	CK_DYNA float* scale;               // diagonal equilibration: S[i] = 1/sqrt(K_ii) per DOF
	float diag_data[LDL_MAX_NODES][21]; // diagonal blocks: packed lower-triangular (max 6x6 = 21)
	float diag_D[LDL_MAX_NODES][6];    // D pivots (max 6 per block)
	int topo_version;   // world topo version when blocks were built

	// Shattering state (weight-based: no virtual body copies)
	int virtual_body_count;          // count of virtual shard indices (for graph topology)
	CK_DYNA int* body_remap;         // real_body_idx -> first virtual shard index (-1 if not shattered)
	CK_DYNA int* shard_counts;       // real_body_idx -> number of shards (0 if not shattered)
} LDL_Cache;

// Island: group of connected bodies that can sleep/wake together.
typedef struct Island
{
	int head_body, tail_body, body_count;
	int head_joint, tail_joint, joint_count;
	int constraint_remove_count;
	int awake; // 1 = awake, 0 = sleeping
	LDL_Cache ldl;
} Island;

#define SLEEP_VEL_THRESHOLD   0.01f  // squared velocity magnitude threshold
#define SLEEP_TIME_THRESHOLD  0.5f   // seconds of stillness before sleep
#define LINEAR_SLOP           0.01f  // contact margin: keeps contacts alive near zero-penetration

typedef struct WorldInternal
{
	int frame;         // monotonically increasing frame counter
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts;
	CK_MAP(WarmManifold) warm_cache;
	CK_MAP(AVBD_WarmManifold) avbd_warm_cache;
	CK_DYNA v3* avbd_prev_velocity; // per-body prev velocity for adaptive warm-start
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
	int ldl_enabled;   // 1 = LDL direct correction for joints (dual solvers only)
	int ldl_topo_version; // incremented on joint create/destroy
	FrictionModel friction_model;
	SolverType solver_type;
	int velocity_iters;
	int position_iters;
	float contact_hertz;
	float contact_damping_ratio;
	float max_push_velocity;
	int sub_steps;
	// AVBD parameters
	float avbd_alpha;       // stabilization (0.95-0.99)
	float avbd_beta_lin;    // penalty ramp, linear constraints
	float avbd_beta_ang;    // penalty ramp, angular constraints
	float avbd_gamma;       // warm-start decay
	int avbd_iterations;    // solver iterations
} WorldInternal;

// -----------------------------------------------------------------------------
// Solver types and constants.

#define SOLVER_VELOCITY_ITERS      10
#define SOLVER_POSITION_ITERS      4
// SOLVER_BAUMGARTE removed: joints use NGS or LDL position correction, not velocity bias.
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

// -----------------------------------------------------------------------------
// AVBD (Augmented Vertex Block Descent) solver types.

#define AVBD_STABLE_THRESH 0.05f // adaptive alpha: below this error, use full stabilization
#define AVBD_PENALTY_MIN  1.0f
#define AVBD_PENALTY_MAX  1e10f
#define AVBD_MARGIN       0.0f   // must be 0: nudge penetration is always positive, margin would hide it
#define AVBD_STICK_THRESH 0.00001f

typedef struct AVBD_Contact
{
	v3 r_a, r_b;       // contact offsets in body-local space
	v3 C0;             // constraint error at x- (normal, tangent1, tangent2)
	float normal_bias; // one-frame bias added directly to the normal C (restitution)
	v3 penalty;        // per-axis penalty parameters (ramp over iterations)
	v3 lambda;         // accumulated dual variables (clamped force)
	int stick;         // static friction flag
	uint32_t feature_id;
} AVBD_Contact;

typedef struct AVBD_Manifold
{
	int body_a, body_b;
	int contact_count;
	float friction;
	m3x3 basis;        // [normal; tangent1; tangent2] row-major
	AVBD_Contact contacts[MAX_CONTACTS];
} AVBD_Manifold;

// CSR adjacency: separate arrays for contacts and joints.
// body i's contact refs at ct_adj[ct_start[i]..ct_start[i+1]]
// body i's joint refs at jt_adj[jt_start[i]..jt_start[i+1]]

typedef struct AVBD_ContactAdj
{
	int manifold_idx;
	int contact_idx;
	int is_body_a;
} AVBD_ContactAdj;

typedef struct AVBD_JointAdj
{
	int joint_idx;
	int is_body_a;
} AVBD_JointAdj;

// AVBD warm cache: persists penalty/lambda across frames (keyed by body pair).
typedef struct AVBD_WarmContact
{
	uint32_t feature_id;
	v3 r_a;             // body-local contact offset A (also used for spatial matching)
	v3 r_b;             // body-local contact offset B
	v3 penalty;
	v3 lambda;
	int stick;
} AVBD_WarmContact;

typedef struct AVBD_WarmManifold
{
	AVBD_WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale;
	int touching;       // previous frame had an active/touching contact in this pair
	m3x3 basis;         // cached [normal; tangent1; tangent2] for synthetic contacts
	float friction;
} AVBD_WarmManifold;

// Per-body temporary state for AVBD (not stored in BodyHot)
typedef struct AVBD_BodyState
{
	v3 inertial_lin;
	v3 initial_lin;    // x- for velocity recovery
	quat inertial_ang;
	quat initial_ang;  // q- for angular velocity recovery
} AVBD_BodyState;

// -----------------------------------------------------------------------------
// LDL debug visualization data (populated by solver_ldl.c, read by UI).

#define LDL_MAX_DOF 128

typedef struct LDL_DebugInfo
{
	int n;                                // total DOFs
	float A[LDL_MAX_DOF * LDL_MAX_DOF];  // constraint matrix snapshot
	float D[LDL_MAX_DOF];                // diagonal pivots after factorization
	float rhs[LDL_MAX_DOF];              // RHS (constraint violations) before solve
	float lambda_ldl[LDL_MAX_DOF];       // exact lambdas after LDL solve
	int valid;
} LDL_DebugInfo;

#endif // NUDGE_INTERNAL_H
