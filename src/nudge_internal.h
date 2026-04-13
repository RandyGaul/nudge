// nudge_internal.h -- internal types for nudge physics engine (unity build)
#ifndef NUDGE_INTERNAL_H
#define NUDGE_INTERNAL_H

// Internal constant: feature ID bit flag for edge-edge contacts.
// Edge contacts: edge_a | (edge_b << 16) | FEATURE_EDGE_BIT
#define FEATURE_EDGE_BIT 0x80000000u

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
		struct { float half_height; float radius; } cylinder;
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

// Hot: solver working set, iterated every PGS step, packed for cache.
// Layout matches old SolverBodyVel so the SIMD gather/scatter macros work directly.
typedef struct BodyHot
{
	v3 velocity;
	v3 angular_velocity;
	v3 iw_diag;           // precomputed world-space inverse inertia (xx,yy,zz)
	v3 iw_off;            // precomputed world-space inverse inertia (xy,xz,yz)
	float inv_mass;
	float _pad[3];
} BodyHot;

// State: per-body data needed outside the inner PGS loop (integration, collision,
// joint geometry, sleep bookkeeping). Parallel array indexed the same as body_hot.
typedef struct BodyState
{
	v3 position;
	quat rotation;
	v3 inv_inertia_local; // diagonal of local-space inverse inertia tensor
	float friction;
	float restitution;
	float linear_damping;
	float angular_damping;
	float sleep_time;     // accumulated seconds below velocity threshold
	int sleep_allowed;    // 1 = can sleep (default), 0 = never sleep
} BodyState;

// Body field accessors: decouple call sites from tier assignment.
// Moving a field between BodyHot/BodyState = change one macro, not 50 files.
#define body_pos(w, i)              (w)->body_state[i].position
#define body_rot(w, i)              (w)->body_state[i].rotation
#define body_vel(w, i)              (w)->body_hot[i].velocity
#define body_angvel(w, i)           (w)->body_hot[i].angular_velocity
#define body_inv_mass(w, i)         (w)->body_hot[i].inv_mass
#define body_iw_diag(w, i)          (w)->body_hot[i].iw_diag
#define body_iw_off(w, i)           (w)->body_hot[i].iw_off
#define body_inv_inertia_local(w, i) (w)->body_state[i].inv_inertia_local
#define body_friction(w, i)         (w)->body_state[i].friction
#define body_restitution(w, i)      (w)->body_state[i].restitution
#define body_lin_damping(w, i)      (w)->body_state[i].linear_damping
#define body_ang_damping(w, i)      (w)->body_state[i].angular_damping
#define body_sleep_time(w, i)       (w)->body_state[i].sleep_time
#define body_sleep_allowed(w, i)    (w)->body_state[i].sleep_allowed

// Scalar body reference: bundles pointers to both tiers from a single index lookup.
// Prevents mismatched indices (body_hot[a] with body_state[b]).
// Usage: BodyRef r = body_ref(w, idx); r.hot->velocity; r.state->position;
typedef struct BodyRef { BodyHot* hot; BodyState* state; } BodyRef;
#define body_ref(w, i) ((BodyRef){ &(w)->body_hot[i], &(w)->body_state[i] })

// Cached narrowphase feature pair for incremental manifold refresh.
// type=0: cold (no cache), type=1: face-face, type=2: edge-edge.
typedef struct CachedFeaturePair
{
	int16_t type;      // 0=cold, 1=face-face, 2=edge-edge
	int16_t ref_body;  // 0=A is reference, 1=B is reference
	int16_t face_a;    // face index on hull A (or box face 0-5)
	int16_t face_b;    // face index on hull B
	int16_t edge_a;    // half-edge index (edge-edge only)
	int16_t edge_b;    // half-edge index (edge-edge only)
} CachedFeaturePair;

typedef struct WarmManifold WarmManifold; // forward decl for warm cache


// Joint persistent storage (handle-based, parallel arrays like bodies).
typedef enum JointType { JOINT_BALL_SOCKET, JOINT_DISTANCE, JOINT_HINGE, JOINT_FIXED, JOINT_PRISMATIC } JointType;

typedef struct JointInternal
{
	JointType type;
	int body_a, body_b; // array indices (resolved from handles at creation)
	union {
		struct { v3 local_a, local_b; SpringParams spring; } ball_socket;
		struct { v3 local_a, local_b; float rest_length; SpringParams spring; float limit_min, limit_max; } distance;
		struct { v3 local_a, local_b; v3 local_axis_a, local_axis_b; v3 local_ref_a, local_ref_b; SpringParams spring; float limit_min, limit_max; float motor_speed, motor_max_impulse; } hinge;
		struct { v3 local_a, local_b; quat local_rel_quat; SpringParams spring; } fixed;
		struct { v3 local_a, local_b; v3 local_axis_a, local_axis_b; quat local_rel_quat; SpringParams spring; float motor_speed, motor_max_impulse; } prismatic;
	};
	// Warm starting: accumulated impulses persisted across frames (up to 6 DOF).
	float warm_lambda[6];
	// Island linked list fields.
	int island_id;    // -1 = none
	int island_prev;  // -1 = head
	int island_next;  // -1 = tail
} JointInternal;

typedef struct BVH_Tree BVH_Tree; // forward decl, defined in bvh.c

// LDL direct solver types (used by solver_ldl.c).
// Per-DOF Jacobian row: J_a and J_b each have 6 components [lin_x, lin_y, lin_z, ang_x, ang_y, ang_z].
// LDL computes K = J * M^{-1} * J^T, RHS = -J * v, and impulse = M^{-1} * J^T * lambda generically.
typedef struct LDL_JacobianRow
{
	double J_a[6];
	double J_b[6];
} LDL_JacobianRow;

typedef struct LDL_Constraint
{
	int type;         // JOINT_BALL_SOCKET, JOINT_DISTANCE, or JOINT_HINGE (or -1 for synthetic)
	int dof;          // 3 or 1
	int body_a, body_b;         // body indices (may be virtual shard indices for graph topology)
	int real_body_a, real_body_b; // always real body indices (for physics data lookup)
	double weight_a, weight_b;   // shard weight (1.0 = normal, S = shattered into S shards)
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
#define LDL_MAX_BLOCK_DIM 12  // max DOF per bundle (12*13/2 = 78 packed elements)

typedef struct LDL_Sparse
{
	int node_count;
	int dof[LDL_MAX_NODES];          // DOF per node
	int row_offset[LDL_MAX_NODES+1]; // cumulative DOF offset
	int n;                            // total scalar DOFs

	// Diagonal blocks: packed lower-triangular, dof*(dof+1)/2 doubles per node.
	double diag_data[LDL_MAX_NODES][78]; // max 12x12 packed (78 elements)
	double diag_D[LDL_MAX_NODES][12];   // D pivots from block LDL (max 12 per block)

	// Off-diagonal: per-node neighbor list.
	// adj[i][k] = neighbor node index, adj_data[i][k] = block data (dof[i]*dof[j] doubles).
	CK_DYNA int* adj[LDL_MAX_NODES];       // neighbor indices
	CK_DYNA double* adj_data[LDL_MAX_NODES]; // packed block doubles per neighbor

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
	int node, dk, ok;             // node index, DOF, row_offset
	int fwd_start, fwd_count;     // forward-sub neighbors (eliminated before)
	int back_start, back_count;   // back-sub neighbors (eliminated after)
	int col_start, col_count;     // L-column entries for factorization
	int schur_start, schur_count; // Schur complement updates for this pivot
} LDL_Pivot;

// K-edge fill instruction: maps shared body to off-diagonal L_factors offsets.
typedef struct LDL_Coupling
{
	int body;                     // shared body index (real or virtual)
	int block_i, block_j;         // constraint node indices
	int L_offset_ij, L_offset_ji; // L_factors offsets for both directions
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
	int L_factors_size; // total floats in L_factors
} LDL_Topology;

typedef struct LDL_Cache
{
	CK_DYNA LDL_Constraint* constraints;
	int joint_count;
	CK_DYNA LDL_Bundle* bundles;
	int bundle_count;
	int n;              // total DOFs
	LDL_Topology* topo; // cached topology (NULL until built)
	CK_DYNA double* L_factors; // contiguous off-diagonal L-factor blocks (double precision)
	CK_DYNA LDL_JacobianRow* jacobians;  // per-DOF Jacobian rows (filled each substep)
	double diag_data[LDL_MAX_NODES][78]; // diagonal blocks: packed lower-triangular (max 12x12 = 78)
	double diag_D[LDL_MAX_NODES][12];    // D pivots (max 12 per block)
	int topo_version;   // world topo version when blocks were built

	// Preallocated solve scratch buffers (sized to n DOFs, reused across substeps)
	CK_DYNA double* solve_rhs;
	CK_DYNA double* solve_lambda;

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
	int frame; // monotonically increasing frame counter
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA BodyState*   body_state;
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
	BVH_Tree* bvh_static;
	BVH_Tree* bvh_dynamic;
	BVH_Tree* bvh_sleeping; // sleeping dynamic bodies (no refit, no SAP)
	// Islands
	CK_DYNA Island*     islands;
	CK_DYNA uint32_t*   island_gen;
	CK_DYNA int*        island_free;
	CK_MAP(uint8_t)     prev_touching; // body_pair_key -> 1 for pairs touching last frame
	CK_MAP(uint8_t)     joint_pairs;   // body_pair_key -> 1 for bodies connected by joints (skip collisions)
	int joint_pairs_version;           // ldl_topo_version when joint_pairs was last built
	int sleep_enabled;       // 1 = island sleep active (default)
	int sat_hint_enabled;    // 1 = warm-cache SAT axis hints fed to narrowphase
	int sat_hillclimb_enabled; // 1 = hill-climb face search for hulls (vs brute force)
	int box_use_hull;          // 1 = route box-box through hull-hull path (debug)
	int incremental_np_enabled; // 1 = incremental narrowphase (cached feature pair refresh)
	int warm_start_enabled;    // 1 = warm-start contact impulses from cache
	// Native cylinder narrowphase toggles (0 = route through hull-backed fallback).
	// Flipped on per-pair as native routines land; default 0 until Phase 6 cleanup.
	int cyl_native_sphere;
	int cyl_native_capsule;
	int cyl_native_box;
	int cyl_native_hull;
	int cyl_native_cyl;
	int thread_count;        // 0 or 1 = single-threaded, >1 = parallel PGS solver
	void* np_pairs_out; // CK_DYNA BroadPair** — when set, broadphase outputs pairs here instead of running narrowphase
	int ldl_enabled;         // 1 = LDL direct correction for joints (dual solvers only)
	int ldl_topo_version;    // incremented on joint create/destroy
	int ldl_correction_iter; // which PGS iter triggers LDL correction (-1 = after all PGS iters)
	SolverType solver_type;
	int velocity_iters;
	int position_iters;
	float contact_hertz;
	float contact_damping_ratio;
	float max_push_velocity;
	int sub_steps;
	PerfTimers perf;

	// Debug: solver arrays from last world_step, kept alive for remote viewer.
	// Freed at the start of the next world_step. Typed void* because SolverManifold
	// and SolverContact are defined later in the unity build.
	void *dbg_solver_manifolds; // CK_DYNA SolverManifold*
	void *dbg_solver_contacts;  // CK_DYNA SolverContact*
	void *dbg_solver_joints;    // CK_DYNA SolverJoint*
} WorldInternal;

// -----------------------------------------------------------------------------
// Solver types and constants.

#define SOLVER_VELOCITY_ITERS      10
#define SOLVER_POSITION_ITERS      2
// SOLVER_BAUMGARTE removed: joints use NGS or LDL position correction, not velocity bias.
#define SOLVER_SLOP                0.005f // NGS position correction dead zone
#define SOLVER_RESTITUTION_THRESH  1.0f
#define SOLVER_POS_BAUMGARTE       0.2f
#define SOLVER_POS_MAX_CORRECTION  0.2f   // max position correction per step
#define SOLVER_MAX_LINEAR_VEL      500.0f
#define SOLVER_MAX_ANGULAR_VEL     100.0f

typedef struct SolverContact
{
	// --- Hot fields for patch-mode PGS inner loop (first 92 bytes = 2 cache lines) ---
	v3 r_a;              // contact offset from body A
	v3 r_b;              // contact offset from body B
	v3 rn_a, rn_b;       // precomputed cross(r, normal)
	v3 w_n_a, w_n_b;     // precomputed I_w * cross(r, normal) — angular impulse direction
	float eff_mass_n;    // 1/K for normal row
	float bias;          // velocity bias for penetration recovery
	float bounce;        // restitution velocity bias
	float softness;      // soft constraint regularization (0 = rigid/NGS)
	float lambda_n;      // accumulated normal impulse (>= 0)
	// --- Cold fields (warm starting, position correction) ---
	v3 normal;           // cached contact normal
	float penetration;   // cached for position correction pass
	uint32_t feature_id; // geometric feature key for warm starting
} SolverContact;


typedef struct SolverManifold
{
	int body_a, body_b;
	int contact_start;
	int contact_count;
	float friction;
	float inv_mass_a, inv_mass_b; // cached from body (avoids body array lookup during iteration)
	// Manifold-level patch friction data (patch friction)
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
} SolverManifold;

// Warm starting: cached impulses from previous frame, keyed by body pair.
typedef struct WarmContact
{
	v3 r_a;              // body-A-relative position for spatial fallback matching
	float lambda_n;
	uint32_t feature_id; // geometric feature key for matching
} WarmContact;

struct WarmManifold
{
	WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale; // 0 = updated this frame, incremented each frame not touched, evicted at >1
	// Manifold-level friction warm data (patch friction)
	float manifold_lambda_t1;
	float manifold_lambda_t2;
	float manifold_lambda_twist;
	// Cached SAT axis for hill-climb warm-start (-1 = no cache).
	int sat_axis;
	// Cached feature pair for incremental narrowphase (skip full SAT when valid).
	CachedFeaturePair cached_pair;
};

// -----------------------------------------------------------------------------
// Joint solver types.

// Max DOF per joint: 6 for future 6-DOF joints (currently hinge = 5).
#define JOINT_MAX_DOF 6

// Per-DOF Jacobian row: J_a and J_b each have 6 components [lin_x, lin_y, lin_z, ang_x, ang_y, ang_z].
// Scalar effective mass = 1 / (J * M^-1 * J^T + softness).
typedef struct JacobianRow
{
	float J_a[6];
	float J_b[6];
	float eff_mass;
} JacobianRow;

// Unified solver joint: all joint types share one struct.
// PGS and LDL both read/write the common fields. Per-DOF Jacobian rows
// encode the constraint geometry; type-specific knowledge lives only in
// joint_fill_rows (called from joints_pre_solve and ldl_refresh_lever_arms).
typedef struct SolverJoint
{
	int body_a, body_b;
	int joint_idx;
	JointType type;
	int dof;              // 1 (distance), 3 (ball socket), 5 (hinge)
	float softness;

	v3 r_a, r_b;         // world-space lever arms
	float lambda[JOINT_MAX_DOF];
	float bias[JOINT_MAX_DOF];
	float pos_error[JOINT_MAX_DOF]; // raw position error per DOF (bias = ptv * pos_error)

	// Clamp bounds (for future joint limits / motors).
	// Default: lo=-FLT_MAX, hi=+FLT_MAX (bilateral, unclamped).
	float lo[JOINT_MAX_DOF];
	float hi[JOINT_MAX_DOF];

	JacobianRow rows[JOINT_MAX_DOF];
} SolverJoint;

// Constraint ref for graph coloring dispatch.
enum { CTYPE_CONTACT, CTYPE_JOINT };

typedef struct ConstraintRef
{
	uint8_t type;
	uint8_t color;
	int index;
	int body_a, body_b;
} ConstraintRef;

// -----------------------------------------------------------------------------
// LDL debug visualization data (populated by solver_ldl.c, read by UI).

#define LDL_MAX_DOF 128

typedef struct LDL_DebugInfo
{
	int n;                                // total DOFs
	double A[LDL_MAX_DOF * LDL_MAX_DOF];  // constraint matrix snapshot
	double D[LDL_MAX_DOF];                // diagonal pivots after factorization
	double rhs[LDL_MAX_DOF];              // RHS (constraint violations) before solve
	double lambda_ldl[LDL_MAX_DOF];       // exact lambdas after LDL solve
	int valid;
} LDL_DebugInfo;

#endif // NUDGE_INTERNAL_H
