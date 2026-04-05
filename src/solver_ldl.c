// solver_ldl.c -- Direct LDL joint correction for equality constraints.
// Reference: Mizerski, "Improving an Iterative Physics Solver Using a Direct Method", GDC 2020.
//
// Solves K * lambda = b where K = J * M^-1 * J^T (constraint-space effective mass),
// lambda = constraint impulses, b = velocity error + bias. K is symmetric positive
// (semi-)definite with block structure from the constraint graph.
//
// Architecture: symbolic/numeric split with cached topology.
//
//   Structs (nudge_internal.h)                              Lifetime
//   -----------------------------------------------------------------
//   LDL_Cache        Per-island root: constraints, topology,  island
//                    L_factors, diag_data, diag_D
//   LDL_Constraint   One per joint. Type, DOF, body pair,     topo change
//                    solver index, bundle offset.
//   LDL_Bundle       Same-body-pair constraint group. One     topo change
//                    graph node per bundle (reduces graph).
//   LDL_Topology     Cached symbolic analysis. Ordering,      topo change
//                    sparsity, precomputed offsets.
//   LDL_Pivot        Per elimination step. Index ranges       topo change
//                    into neighbor/column/schur arrays.
//   LDL_Neighbor     Adjacency for solve. Node, DOF,          topo change
//                    row offset, L_factors offset.
//   LDL_Column       L-block target during factorization.     topo change
//                    Node, DOF, L_factors write offset.
//   LDL_Schur        Schur complement update. Source and      topo change
//                    target L_factors offsets, dimensions.
//   LDL_Coupling     Shared-body K contribution. Body,        topo change
//                    node pair, L_factors offsets (both).
//   LDL_Sparse       Temporary adjacency-list graph for       transient
//                    min-degree ordering + fill-in.
//
//                    topology change
//                          |
//                          v
//   joints ------> LDL_Constraint[] ------> ldl_build_bundles
//                                                |
//                                           LDL_Bundle[] (same-pair grouping)
//                                                |
//                                                v
//                                           ldl_build_topology
//                                                |
//                                 +---> LDL_Sparse (temp)
//                                 |     min-degree ordering
//                                 |     fill-in edges
//                                 |     freed
//                                 |
//                                 +---> LDL_Topology (cached)
//                                         LDL_Pivot[]
//                                         LDL_Neighbor[] (fwd + back)
//                                         LDL_Column[]
//                                         LDL_Schur[]
//                                         LDL_Coupling[]
//                                              |
//                         each substep         |
//                             |                |
//                             v                v
//                       ldl_numeric_factor ----------> L_factors[]
//                         LDL_Coupling[] -> K fill     diag_data[]
//                         LDL_Column[]   -> L blocks   diag_D[]
//                         LDL_Schur[]    -> updates
//                                               |
//                                               v
//                                       ldl_solve_topo ------> delta_lambda[]
//                                         LDL_Neighbor[] (fwd)
//                                         diag_data[], diag_D[]
//                                         LDL_Neighbor[] (back)
//
//   Delta-correction mode:
//     PGS runs first with its own warm-starting and clamping.
//     LDL measures the residual error and computes a correction delta.
//     delta_lambda is ADDED to bodies on top of PGS result.
//     s->lambda accumulates PGS + LDL for next frame's warm-start.
//
//   Body shattering (currently disabled):
//     Splits hub bodies into virtual shards connected by synthetic rigid joints
//     to reduce fill-in during factorization. Needs greedy bin-packing, wrap-around
//     chain topology, and velocity sync averaging before re-enabling.

static LDL_DebugInfo g_ldl_debug_info;
int g_ldl_debug_enabled;
int g_ldl_debug_island = -1; // which island to capture debug data for (-1 = none)

#define SHATTER_THRESHOLD 15 // DOF threshold: 6+ ball-sockets (18 DOF) triggers shattering
#define SHARD_TARGET      6 // target DOF per shard
#define LDL_COMPLIANCE    5e-5f // regularization: compliance * trace(K) / dim on K diagonal + reg*lambda in RHS

// -----------------------------------------------------------------------------
// Small block math helpers (3x3 and 1x1 blocks stored as flat floats).

// In-place NxN LDL^T factorize. L stored in lower triangle of A (unit diagonal), D separate.
static void block_ldl(float* A, float* D, int n)
{
	// Track maximum original diagonal to set a relative pivot floor.
	// This prevents near-zero pivots from producing huge L entries in
	// overconstrained systems where Schur updates drive diagonals negative.
	float max_diag = 0;
	for (int j = 0; j < n; j++) {
		float ad = fabsf(A[j * n + j]);
		if (ad > max_diag) max_diag = ad;
	}
	float pivot_floor = max_diag * 1e-6f;
	if (pivot_floor < 1e-12f) pivot_floor = 1e-12f;

	for (int j = 0; j < n; j++) {
		float dj = A[j * n + j];
		for (int k = 0; k < j; k++) dj -= A[j * n + k] * A[j * n + k] * D[k];
		D[j] = fabsf(dj) > pivot_floor ? dj : pivot_floor;
		float inv_dj = 1.0f / D[j];
		for (int i = j + 1; i < n; i++) {
			float lij = A[i * n + j];
			for (int k = 0; k < j; k++) lij -= A[i * n + k] * A[j * n + k] * D[k];
			A[i * n + j] = lij * inv_dj;
		}
	}
}

// Solve NxN LDL system: L stored in lower triangle of A (unit diagonal), D separate.
static void block_solve(float* L, float* D, float* b, float* x, int n)
{
	// Forward: Ly = b
	for (int i = 0; i < n; i++) { x[i] = b[i]; for (int k = 0; k < i; k++) x[i] -= L[i * n + k] * x[k]; }
	// Diagonal
	for (int i = 0; i < n; i++) x[i] /= D[i];
	// Back: L^T x = z
	for (int i = n - 1; i >= 0; i--) for (int k = i + 1; k < n; k++) x[i] -= L[k * n + i] * x[k];
}

// Compute C = A * B where A is (ra x ca), B is (ca x cb). Out is (ra x cb).
static void block_mul(float* A, int ra, int ca, float* B, int cb, float* out)
{
	for (int i = 0; i < ra; i++) {
		for (int j = 0; j < cb; j++) {
			float sum = 0;
			for (int k = 0; k < ca; k++) sum += A[i*ca + k] * B[k*cb + j];
			out[i*cb + j] = sum;
		}
	}
}

// A -= B, both (r x c).
static void block_sub(float* A, float* B, int r, int c)
{
	for (int i = 0; i < r * c; i++) A[i] -= B[i];
}

// Transpose: out = A^T where A is (r x c), out is (c x r).
static void block_transpose(float* A, int r, int c, float* out)
{
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			out[j*r + i] = A[i*c + j];
		}
	}
}

// -----------------------------------------------------------------------------
// Sparse block matrix operations.

static void ldl_sparse_init(LDL_Sparse* s)
{
	memset(s, 0, sizeof(*s));
}

static void ldl_sparse_free(LDL_Sparse* s)
{
	for (int i = 0; i < s->node_count; i++) {
		afree(s->adj[i]);
		afree(s->adj_data[i]);
	}
	memset(s, 0, sizeof(*s));
}

// Find the off-diagonal block index for edge (i, j). Returns -1 if not found.
static int ldl_sparse_find_edge(LDL_Sparse* s, int i, int j)
{
	int cnt = asize(s->adj[i]);
	for (int k = 0; k < cnt; k++) {
		if (s->adj[i][k] == j) return k;
	}
	return -1;
}

// Get or create off-diagonal block for edge (i, j). Returns pointer to block data.
// Creates the block (zeroed) in both adj[i] and adj[j] if it doesn't exist.
static float* ldl_sparse_get_or_create_edge(LDL_Sparse* s, int i, int j)
{
	int ki = ldl_sparse_find_edge(s, i, j);
	if (ki >= 0) return &s->adj_data[i][ki * s->dof[i] * s->dof[j]];

	int di = s->dof[i], dj = s->dof[j];

	// Add j to adj[i]
	int block_size_ij = di * dj;
	int old_size_i = asize(s->adj_data[i]);
	apush(s->adj[i], j);
	afit(s->adj_data[i], old_size_i + block_size_ij);
	asetlen(s->adj_data[i], old_size_i + block_size_ij);
	float* data_ij = &s->adj_data[i][old_size_i];
	memset(data_ij, 0, block_size_ij * sizeof(float));

	// Add i to adj[j]
	int block_size_ji = dj * di;
	int old_size_j = asize(s->adj_data[j]);
	apush(s->adj[j], i);
	afit(s->adj_data[j], old_size_j + block_size_ji);
	asetlen(s->adj_data[j], old_size_j + block_size_ji);
	float* data_ji = &s->adj_data[j][old_size_j];
	memset(data_ji, 0, block_size_ji * sizeof(float));

	return data_ij;
}

// Get existing off-diagonal block data for edge (i, j). Returns NULL if not found.
static float* ldl_sparse_get_edge(LDL_Sparse* s, int i, int j)
{
	int ki = ldl_sparse_find_edge(s, i, j);
	if (ki < 0) return NULL;
	return &s->adj_data[i][ki * s->dof[i] * s->dof[j]];
}

// -----------------------------------------------------------------------------
// Diagonal block builders (K matrix computation).

static void ldl_ball_socket_K(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, float softness, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m };

	v3 ia = a->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		v3 e0 = inv_inertia_mul(a->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(a->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(a->rotation, ia, V3(-r_a.y, r_a.x, 0));
		K[0] += -r_a.z*e0.y + r_a.y*e0.z;
		K[1] += -r_a.z*e1.y + r_a.y*e1.z;
		K[2] += -r_a.z*e2.y + r_a.y*e2.z;
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z;
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z;
		K[5] += -r_a.y*e2.x + r_a.x*e2.y;
	}

	v3 ib = b->inv_inertia_local;
	if (ib.x > 0 || ib.y > 0 || ib.z > 0) {
		v3 e0 = inv_inertia_mul(b->rotation, ib, V3(0, -r_b.z, r_b.y));
		v3 e1 = inv_inertia_mul(b->rotation, ib, V3(r_b.z, 0, -r_b.x));
		v3 e2 = inv_inertia_mul(b->rotation, ib, V3(-r_b.y, r_b.x, 0));
		K[0] += -r_b.z*e0.y + r_b.y*e0.z;
		K[1] += -r_b.z*e1.y + r_b.y*e1.z;
		K[2] += -r_b.z*e2.y + r_b.y*e2.z;
		K[3] +=  r_b.z*e1.x - r_b.x*e1.z;
		K[4] +=  r_b.z*e2.x - r_b.x*e2.z;
		K[5] += -r_b.y*e2.x + r_b.x*e2.y;
	}

	K[0] += softness; K[3] += softness; K[5] += softness;
	// Write as full 3x3 row-major
	out[0] = K[0]; out[1] = K[1]; out[2] = K[2];
	out[3] = K[1]; out[4] = K[3]; out[5] = K[4];
	out[6] = K[2]; out[7] = K[4]; out[8] = K[5];
}

// Write symmetric 3x3 stored as 6 floats into full 3x3.
static void ldl_write_sym3x3(float* out, const float* K)
{
	out[0] = K[0]; out[1] = K[1]; out[2] = K[2];
	out[3] = K[1]; out[4] = K[3]; out[5] = K[4];
	out[6] = K[2]; out[7] = K[4]; out[8] = K[5];
}

// 6x6 K for synthetic weld (zero lever arm). Constrains both linear and angular velocity.
// K = diag(inv_m_sum * I_3, world_inv_inertia_sum)
static void ldl_weld_K(BodyHot* a, BodyHot* b, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	memset(out, 0, 36 * sizeof(float));
	// Linear block (top-left 3x3): (inv_m_a + inv_m_b) * I
	out[0] = inv_m; out[7] = inv_m; out[14] = inv_m;
	// Angular block (bottom-right 3x3): world_inv_inertia_a + world_inv_inertia_b
	// For diagonal inertia rotated to world: R * diag(I) * R^T
	// inv_inertia_mul(q, inv_I, v) = R * diag(inv_I) * R^T * v
	// We need the 3x3 matrix itself. Build it column by column.
	v3 ia = a->inv_inertia_local, ib = b->inv_inertia_local;
	v3 col0a = inv_inertia_mul(a->rotation, ia, V3(1,0,0));
	v3 col1a = inv_inertia_mul(a->rotation, ia, V3(0,1,0));
	v3 col2a = inv_inertia_mul(a->rotation, ia, V3(0,0,1));
	v3 col0b = inv_inertia_mul(b->rotation, ib, V3(1,0,0));
	v3 col1b = inv_inertia_mul(b->rotation, ib, V3(0,1,0));
	v3 col2b = inv_inertia_mul(b->rotation, ib, V3(0,0,1));
	out[21] = col0a.x + col0b.x; out[22] = col1a.x + col1b.x; out[23] = col2a.x + col2b.x;
	out[27] = col0a.y + col0b.y; out[28] = col1a.y + col1b.y; out[29] = col2a.y + col2b.y;
	out[33] = col0a.z + col0b.z; out[34] = col1a.z + col1b.z; out[35] = col2a.z + col2b.z;
}

// Apply a 6-DOF weld impulse (linear + angular) at body center (zero lever arm).
static void apply_weld_impulse(BodyHot* a, BodyHot* b, float* lambda6)
{
	v3 lin = V3(lambda6[0], lambda6[1], lambda6[2]);
	v3 ang = V3(lambda6[3], lambda6[4], lambda6[5]);
	a->velocity = sub(a->velocity, scale(lin, a->inv_mass));
	b->velocity = add(b->velocity, scale(lin, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, ang));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, ang));
}

// -----------------------------------------------------------------------------
// Off-diagonal block computation (same math as before, writes to block buffer).

// Ball socket x ball socket (3x3) for shared body B.
static void ldl_compute_off_diag_bs_bs(BodyHot* B, v3 ri, v3 rj, float sign, float* out)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	v3 e0 = inv_inertia_mul(B->rotation, ib, V3(0, -rj.z, rj.y));
	v3 e1 = inv_inertia_mul(B->rotation, ib, V3(rj.z, 0, -rj.x));
	v3 e2 = inv_inertia_mul(B->rotation, ib, V3(-rj.y, rj.x, 0));

	out[0] = sign * (inv_m + (-ri.z*e0.y + ri.y*e0.z));
	out[1] = sign * (        (-ri.z*e1.y + ri.y*e1.z));
	out[2] = sign * (        (-ri.z*e2.y + ri.y*e2.z));
	out[3] = sign * (        ( ri.z*e0.x - ri.x*e0.z));
	out[4] = sign * (inv_m + ( ri.z*e1.x - ri.x*e1.z));
	out[5] = sign * (        ( ri.z*e2.x - ri.x*e2.z));
	out[6] = sign * (        (-ri.y*e0.x + ri.x*e0.y));
	out[7] = sign * (        (-ri.y*e1.x + ri.x*e1.y));
	out[8] = sign * (inv_m + (-ri.y*e2.x + ri.x*e2.y));
}

// Ball socket x weld (3x6). Weld Jacobian is +-I_6, so coupling = sign * J_bs * M_B^-1.
// out is 3x6 row-major: [linear_3x3 | angular_3x3].
static void ldl_compute_off_diag_bs_weld(BodyHot* B, v3 ri, float sign, float* out)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	// Linear part (3x3): sign * inv_m * I + angular cross terms from r_i
	// Full formula: sign * [I | -skew(r_i)] * M^-1 = sign * [inv_m*I | -skew(r_i)*inv_I]
	// But we also need to account for the sign of the weld side (already folded into sign).
	memset(out, 0, 18 * sizeof(float));
	// Linear columns (0-2): sign * inv_m * I_3
	out[0] = sign * inv_m; out[7] = sign * inv_m; out[14] = sign * inv_m;
	// Angular columns (3-5): sign * (-skew(r_i) * inv_I)
	// -skew(r_i) * inv_I * e_j for each column j
	v3 c0 = inv_inertia_mul(B->rotation, ib, V3(1,0,0));
	v3 c1 = inv_inertia_mul(B->rotation, ib, V3(0,1,0));
	v3 c2 = inv_inertia_mul(B->rotation, ib, V3(0,0,1));
	// -skew(r) * col = -(r × col) = col × r
	v3 a0 = cross(c0, ri), a1 = cross(c1, ri), a2 = cross(c2, ri);
	out[3]  = sign * a0.x; out[4]  = sign * a1.x; out[5]  = sign * a2.x;
	out[9]  = sign * a0.y; out[10] = sign * a1.y; out[11] = sign * a2.y;
	out[15] = sign * a0.z; out[16] = sign * a1.z; out[17] = sign * a2.z;
}

// Distance x weld (1x6). out is 1x6.
static void ldl_compute_off_diag_dist_weld(BodyHot* B, v3 ri, v3 axis, float sign, float* out)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	// Linear part: sign * inv_m * axis^T (1x3)
	out[0] = sign * inv_m * axis.x;
	out[1] = sign * inv_m * axis.y;
	out[2] = sign * inv_m * axis.z;
	// Angular part: sign * (r × axis)^T * inv_I (1x3)
	v3 rxa = cross(ri, axis);
	v3 Irxa = inv_inertia_mul(B->rotation, ib, rxa);
	out[3] = sign * Irxa.x;
	out[4] = sign * Irxa.y;
	out[5] = sign * Irxa.z;
}

// Weld x weld (6x6). Two welds sharing a shard body. out is 6x6 row-major.
static void ldl_compute_off_diag_weld_weld(BodyHot* B, float sign, float* out)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	memset(out, 0, 36 * sizeof(float));
	// Linear block: sign * inv_m * I_3
	out[0] = sign * inv_m; out[7] = sign * inv_m; out[14] = sign * inv_m;
	// Angular block: sign * inv_I
	v3 c0 = inv_inertia_mul(B->rotation, ib, V3(1,0,0));
	v3 c1 = inv_inertia_mul(B->rotation, ib, V3(0,1,0));
	v3 c2 = inv_inertia_mul(B->rotation, ib, V3(0,0,1));
	out[21] = sign * c0.x; out[22] = sign * c1.x; out[23] = sign * c2.x;
	out[27] = sign * c0.y; out[28] = sign * c1.y; out[29] = sign * c2.y;
	out[33] = sign * c0.z; out[34] = sign * c1.z; out[35] = sign * c2.z;
}

// Ball socket x distance (3x1).
static void ldl_compute_off_diag_bs_dist(BodyHot* B, v3 ri, v3 rj, v3 axis, float sign, float* out)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	v3 Irj = inv_inertia_mul(B->rotation, ib, cross(rj, axis));
	out[0] = sign * (inv_m * axis.x + (-ri.z * Irj.y + ri.y * Irj.z));
	out[1] = sign * (inv_m * axis.y + ( ri.z * Irj.x - ri.x * Irj.z));
	out[2] = sign * (inv_m * axis.z + (-ri.y * Irj.x + ri.x * Irj.y));
}

// Distance x distance (1x1).
static float ldl_compute_off_diag_dist_dist(BodyHot* B, v3 ri, v3 rj, v3 axis_i, v3 axis_j, float sign)
{
	float inv_m = B->inv_mass;
	v3 ib = B->inv_inertia_local;
	v3 Irj = inv_inertia_mul(B->rotation, ib, cross(rj, axis_j));
	return sign * (inv_m * dot(axis_i, axis_j) + dot(cross(ri, axis_i), Irj));
}

// -----------------------------------------------------------------------------
// Sparse LDL factorize and solve.

// Count fill-in edges that would be created by eliminating node i.
static int ldl_sparse_fill_count(LDL_Sparse* s, int i, int* eliminated)
{
	int cnt = asize(s->adj[i]);
	int fill = 0;
	for (int ii = 0; ii < cnt; ii++) {
		int ni = s->adj[i][ii];
		if (eliminated[ni]) continue;
		for (int jj = ii + 1; jj < cnt; jj++) {
			int nj = s->adj[i][jj];
			if (eliminated[nj]) continue;
			if (ldl_sparse_find_edge(s, ni, nj) < 0) {
				fill++;
			}
		}
	}
	return fill;
}

// Minimum local fill elimination ordering.
// Picks the node whose elimination creates the fewest new edges (fill-in).
// Better than minimum degree for chain-ending-in-clique topologies (ragdolls).
static void ldl_sparse_min_fill_order(LDL_Sparse* s)
{
	int nc = s->node_count;
	int eliminated[LDL_MAX_NODES] = {0};
	int fill_cost[LDL_MAX_NODES];

	// Initial fill cost for each node
	for (int i = 0; i < nc; i++) {
		fill_cost[i] = ldl_sparse_fill_count(s, i, eliminated);
	}

	for (int step = 0; step < nc; step++) {
		// Pick node with minimum fill cost
		int best = -1, best_fill = 0x7FFFFFFF;
		for (int i = 0; i < nc; i++) {
			if (eliminated[i]) continue;
			if (fill_cost[i] < best_fill) { best_fill = fill_cost[i]; best = i; }
		}
		s->elim_order[step] = best;
		eliminated[best] = 1;

		// Eliminate best: add fill-in edges between its remaining neighbors
		int cnt = asize(s->adj[best]);
		for (int ii = 0; ii < cnt; ii++) {
			int ni = s->adj[best][ii];
			if (eliminated[ni]) continue;
			for (int jj = ii + 1; jj < cnt; jj++) {
				int nj = s->adj[best][jj];
				if (eliminated[nj]) continue;
				if (ldl_sparse_find_edge(s, ni, nj) < 0) {
					ldl_sparse_get_or_create_edge(s, ni, nj);
				}
			}
		}

		// Recompute fill cost only for neighbors of eliminated node
		for (int ii = 0; ii < cnt; ii++) {
			int ni = s->adj[best][ii];
			if (eliminated[ni]) continue;
			fill_cost[ni] = ldl_sparse_fill_count(s, ni, eliminated);
		}
	}

	// Build inverse permutation for O(1) lookups in solve
	for (int step = 0; step < nc; step++) {
		s->inv_order[s->elim_order[step]] = step;
	}
}

// Reorder pivots into depth-first order on the elimination tree.
// The elimination tree parent of step s is the earliest step t > s such that
// nodes elim_order[s] and elim_order[t] share an edge. DFS traversal groups
// related pivots consecutively for better cache coherence during factorization.
static void ldl_sparse_dfs_reorder(LDL_Sparse* s)
{
	int nc = s->node_count;
	if (nc <= 2) return;

	// Build elimination tree: parent[step] = first later step sharing an edge
	int parent[LDL_MAX_NODES];
	for (int i = 0; i < nc; i++) parent[i] = -1;

	for (int step = 0; step < nc; step++) {
		int node = s->elim_order[step];
		int cnt = asize(s->adj[node]);
		int best_step = nc; // sentinel: no parent
		for (int ai = 0; ai < cnt; ai++) {
			int neighbor = s->adj[node][ai];
			int nstep = s->inv_order[neighbor];
			if (nstep > step && nstep < best_step) {
				best_step = nstep;
			}
		}
		parent[step] = (best_step < nc) ? best_step : -1;
	}

	// Build child lists for DFS
	CK_DYNA int* children[LDL_MAX_NODES];
	memset(children, 0, sizeof(children));
	int roots[LDL_MAX_NODES], root_count = 0;
	for (int i = 0; i < nc; i++) {
		if (parent[i] >= 0) {
			apush(children[parent[i]], i);
		} else {
			roots[root_count++] = i;
		}
	}

	// Iterative DFS to produce new ordering
	int new_order[LDL_MAX_NODES];
	int out_idx = 0;
	int stack[LDL_MAX_NODES], sp = 0;
	for (int ri = root_count - 1; ri >= 0; ri--) {
		stack[sp++] = roots[ri];
	}
	while (sp > 0) {
		int step = stack[--sp];
		new_order[out_idx++] = s->elim_order[step];
		int cc = asize(children[step]);
		for (int ci = cc - 1; ci >= 0; ci--) {
			stack[sp++] = children[step][ci];
		}
	}

	// Apply new ordering
	for (int i = 0; i < nc; i++) {
		s->elim_order[i] = new_order[i];
	}

	// Rebuild inverse permutation
	for (int step = 0; step < nc; step++) {
		s->inv_order[s->elim_order[step]] = step;
	}

	for (int i = 0; i < nc; i++) {
		afree(children[i]);
	}
}

// -----------------------------------------------------------------------------
// Body shattering: split hub bodies into virtual shards.

// Apply shattering to the block list. Creates virtual bodies and synthetic joints.
// Modifies c->constraints (adds synthetic blocks, redirects body refs to virtual indices).
// Virtual body indices start at virtual_base = body_count (real bodies).
static void ldl_apply_shattering(LDL_Cache* c, WorldInternal* w)
{
	int jc = c->joint_count;
	int body_count = asize(w->body_hot);

	// Count total DOF per body
	int body_dof[4096] = {0}; // plenty for any reasonable body count
	assert(body_count < 4096);
	for (int i = 0; i < jc; i++) {
		body_dof[c->constraints[i].body_a] += c->constraints[i].dof;
		body_dof[c->constraints[i].body_b] += c->constraints[i].dof;
	}

	// Find hub bodies exceeding threshold
	afree(c->virtual_bodies);
	c->virtual_bodies = NULL;
	c->virtual_body_count = 0;

	// Map: real body -> first virtual body index (-1 if not shattered)
	afree(c->body_remap);
	c->body_remap = NULL;
	afit(c->body_remap, body_count);
	asetlen(c->body_remap, body_count);
	for (int b = 0; b < body_count; b++) c->body_remap[b] = -1;

	int virtual_base = body_count;
	typedef struct ShatterInfo { int body; int shard_count; int first_virtual; } ShatterInfo;
	CK_DYNA ShatterInfo* shatter_list = NULL;

	for (int b = 0; b < body_count; b++) {
		if (body_dof[b] <= SHATTER_THRESHOLD) continue;
		if (w->body_hot[b].inv_mass == 0.0f) continue; // don't shatter static bodies

		int S = (body_dof[b] + SHARD_TARGET - 1) / SHARD_TARGET;
		if (S < 2) S = 2;

		int first_virt = virtual_base + c->virtual_body_count;
		c->body_remap[b] = first_virt;

		// Create S virtual body entries with fractional mass
		for (int s = 0; s < S; s++) {
			BodyHot vb = w->body_hot[b];
			vb.inv_mass *= S;
			vb.inv_inertia_local = scale(vb.inv_inertia_local, (float)S);
			apush(c->virtual_bodies, vb);
			c->virtual_body_count++;
		}

		ShatterInfo si = { .body = b, .shard_count = S, .first_virtual = first_virt };
		apush(shatter_list, si);
	}

	if (asize(shatter_list) == 0) {
		afree(shatter_list);
		return; // nothing to shatter
	}

	// Greedy bin-packing: assign each constraint to the shard with smallest current DOF load.
	// Process higher-DOF constraints first for better balance.
	for (int si = 0; si < asize(shatter_list); si++) {
		ShatterInfo* info = &shatter_list[si];
		int hub = info->body;
		int S = info->shard_count;

		// Collect constraints touching this hub, sorted by DOF descending
		CK_DYNA int* touching = NULL;
		for (int i = 0; i < jc; i++) {
			if (c->constraints[i].body_a == hub || c->constraints[i].body_b == hub) {
				apush(touching, i);
			}
		}
		// Insertion sort by DOF descending
		int tc = asize(touching);
		for (int i = 1; i < tc; i++) {
			int tmp = touching[i];
			int j = i - 1;
			while (j >= 0 && c->constraints[touching[j]].dof < c->constraints[tmp].dof) {
				touching[j + 1] = touching[j];
				j--;
			}
			touching[j + 1] = tmp;
		}

		// Track DOF load per shard
		int shard_load[128] = {0};
		assert(S <= 128);
		for (int ti = 0; ti < tc; ti++) {
			int ci = touching[ti];
			// Find shard with smallest load
			int best = 0;
			for (int s = 1; s < S; s++) {
				if (shard_load[s] < shard_load[best]) best = s;
			}
			shard_load[best] += c->constraints[ci].dof;
			// Redirect the hub body reference to the chosen shard
			if (c->constraints[ci].body_a == hub) {
				c->constraints[ci].body_a = info->first_virtual + best;
			}
			if (c->constraints[ci].body_b == hub) {
				c->constraints[ci].body_b = info->first_virtual + best;
			}
		}
		afree(touching);
	}

	// Add synthetic rigid joints in a wrap-around chain between shards.
	// Wrapping the last shard back to the first creates a ring topology
	// that couples all shards more rigidly than a linear chain.
	for (int si = 0; si < asize(shatter_list); si++) {
		ShatterInfo* info = &shatter_list[si];
		for (int s = 0; s < info->shard_count; s++) {
			int next = (s + 1) % info->shard_count;
			LDL_Constraint synth = {
				.type = -1,
				.dof = 6,
				.body_a = info->first_virtual + s,
				.body_b = info->first_virtual + next,
				.solver_idx = -1,
				.is_synthetic = 1,
			};
			apush(c->constraints, synth);
			c->joint_count++;
			c->n += 6;
		}
	}

	afree(shatter_list);
}

// Get BodyHot for a body index (real or virtual).
static BodyHot* ldl_get_body(WorldInternal* w, LDL_Cache* c, int body_idx)
{
	int body_count = asize(w->body_hot);
	if (body_idx < body_count) return &w->body_hot[body_idx];
	return &c->virtual_bodies[body_idx - body_count];
}

// -----------------------------------------------------------------------------
// Topology building: cached symbolic analysis of constraint graph.

static void ldl_topology_free(LDL_Topology* t)
{
	afree(t->fwd_neighbors);
	afree(t->back_neighbors);
	afree(t->columns);
	afree(t->schurs);
	afree(t->couplings);
	memset(t, 0, sizeof(*t));
}

// Build the topology from the current block list. Uses a temporary LDL_Sparse
// to compute elimination ordering and fill-in edges, then extracts the
// precomputed index ranges and L_factors offsets.
static void ldl_build_topology(LDL_Cache* c, WorldInternal* w)
{
	int bc = c->bundle_count;
	if (bc == 0) return;

	// Allocate topology
	if (c->topo) { ldl_topology_free(c->topo); CK_FREE(c->topo); }
	c->topo = CK_ALLOC(sizeof(LDL_Topology));
	memset(c->topo, 0, sizeof(LDL_Topology));
	LDL_Topology* t = c->topo;

	t->node_count = bc;
	t->n = 0;
	for (int i = 0; i < bc; i++) {
		t->dof[i] = c->bundles[i].dof;
		t->row_offset[i] = t->n;
		t->n += t->dof[i];
	}
	t->row_offset[bc] = t->n;
	c->n = t->n;

	// Build temporary LDL_Sparse for elimination ordering + fill-in
	LDL_Sparse tmp;
	ldl_sparse_init(&tmp);
	tmp.node_count = bc;
	tmp.n = t->n;
	for (int i = 0; i < bc; i++) {
		tmp.dof[i] = t->dof[i];
		tmp.row_offset[i] = t->row_offset[i];
	}
	tmp.row_offset[bc] = t->n;

	// Build body adjacency to create graph edges (topology only, no numeric values)
	int real_body_count = asize(w->body_hot);
	int total_body_count = real_body_count + c->virtual_body_count;
	CK_DYNA int** body_adj = CK_ALLOC(total_body_count * sizeof(int*));
	memset(body_adj, 0, total_body_count * sizeof(int*));
	for (int i = 0; i < bc; i++) {
		apush(body_adj[c->bundles[i].body_a], i);
		apush(body_adj[c->bundles[i].body_b], i);
	}

	// Create edges for all pairs of constraints sharing a body
	for (int b = 0; b < total_body_count; b++) {
		int cnt = asize(body_adj[b]);
		for (int ii = 0; ii < cnt; ii++) {
			for (int jj = ii + 1; jj < cnt; jj++) {
				ldl_sparse_get_or_create_edge(&tmp, body_adj[b][ii], body_adj[b][jj]);
			}
		}
	}

	// Run elimination ordering (adds fill-in edges)
	ldl_sparse_min_fill_order(&tmp);
	ldl_sparse_dfs_reorder(&tmp);

	// Copy ordering
	for (int i = 0; i < bc; i++) {
		t->elim_order[i] = tmp.elim_order[i];
		t->inv_order[i] = tmp.inv_order[i];
	}

	// Assign L_factors offsets for every directed edge in the graph (including fill-in).
	// Each undirected edge (a,b) gets two slots: one for (a,b) direction, one for (b,a).
	int edge_L_offset[LDL_MAX_NODES][LDL_MAX_NODES];
	memset(edge_L_offset, -1, sizeof(edge_L_offset));
	int L_offset_counter = 0;
	for (int i = 0; i < bc; i++) {
		int di = tmp.dof[i];
		int cnt = asize(tmp.adj[i]);
		for (int ai = 0; ai < cnt; ai++) {
			int j = tmp.adj[i][ai];
			if (edge_L_offset[i][j] >= 0) continue; // already assigned
			int dj = tmp.dof[j];
			edge_L_offset[i][j] = L_offset_counter;
			L_offset_counter += di * dj;
			edge_L_offset[j][i] = L_offset_counter;
			L_offset_counter += dj * di;
		}
	}
	t->L_factors_size = L_offset_counter;

	// Walk elimination order to build pivot plans
	int eliminated[LDL_MAX_NODES] = {0};

	for (int step = 0; step < bc; step++) {
		int k = tmp.elim_order[step];
		int dk = tmp.dof[k];
		LDL_Pivot* pv = &t->pivots[step];
		pv->node = k;
		pv->dk = dk;
		pv->ok = t->row_offset[k];

		// Classify neighbors
		int cnt = asize(tmp.adj[k]);
		pv->fwd_start = asize(t->fwd_neighbors);
		pv->back_start = asize(t->back_neighbors);
		pv->col_start = asize(t->columns);
		pv->schur_start = asize(t->schurs);

		// Collect non-eliminated neighbors
		CK_DYNA int* later_nodes = NULL;
		for (int ai = 0; ai < cnt; ai++) {
			int j = tmp.adj[k][ai];
			int dj = tmp.dof[j];
			if (eliminated[j]) {
				// j eliminated before k -> forward-sub reads L_{k,j}
				LDL_Neighbor sn = { .node = j, .dn = dj, .on = t->row_offset[j], .L_offset = edge_L_offset[k][j] };
				apush(t->fwd_neighbors, sn);
			} else {
				// j not yet eliminated -> factorize writes L_{j,k}, back-sub reads it
				LDL_Column en = { .node = j, .dn = dj, .L_offset = edge_L_offset[j][k] };
				apush(t->columns, en);
				LDL_Neighbor sn = { .node = j, .dn = dj, .on = t->row_offset[j], .L_offset = edge_L_offset[j][k] };
				apush(t->back_neighbors, sn);
				apush(later_nodes, j);
			}
		}

		// Schur ops: all pairs of non-eliminated neighbors
		int later_count = asize(later_nodes);
		for (int ii = 0; ii < later_count; ii++) {
			int ni = later_nodes[ii];
			int di = tmp.dof[ni];
			for (int jj = ii; jj < later_count; jj++) {
				int nj = later_nodes[jj];
				int djj = tmp.dof[nj];
				LDL_Schur op = {0};
				op.i = ni;
				op.j = nj;
				op.di = di;
				op.dk = dk;
				op.dj = djj;
				op.Lik_offset = edge_L_offset[ni][k];
				op.Ljk_offset = edge_L_offset[nj][k];
				if (ni == nj) {
					op.target_offset = -1;
					op.target_offset_rev = -1;
					op.target_node = ni;
				} else {
					op.target_offset = edge_L_offset[ni][nj];
					op.target_offset_rev = edge_L_offset[nj][ni];
				}
				apush(t->schurs, op);
			}
		}

		pv->fwd_count = asize(t->fwd_neighbors) - pv->fwd_start;
		pv->back_count = asize(t->back_neighbors) - pv->back_start;
		pv->col_count = asize(t->columns) - pv->col_start;
		pv->schur_count = asize(t->schurs) - pv->schur_start;

		afree(later_nodes);
		eliminated[k] = 1;
	}

	// Build K-edge fill instructions from body adjacency
	for (int b = 0; b < total_body_count; b++) {
		int cnt = asize(body_adj[b]);
		if (cnt < 2) continue;
		for (int ii = 0; ii < cnt; ii++) {
			for (int jj = ii + 1; jj < cnt; jj++) {
				int bi = body_adj[b][ii], bj = body_adj[b][jj];
				// Find existing L_offset for this edge pair
				int off_ij = edge_L_offset[bi][bj];
				int off_ji = edge_L_offset[bj][bi];
				// These should exist since we created edges above
				assert(off_ij >= 0 && off_ji >= 0);
				LDL_Coupling fill = { .body = b, .block_i = bi, .block_j = bj, .L_offset_ij = off_ij, .L_offset_ji = off_ji };
				apush(t->couplings, fill);
			}
		}
	}

	// Free temporary data
	for (int b = 0; b < total_body_count; b++) afree(body_adj[b]);
	CK_FREE(body_adj);
	ldl_sparse_free(&tmp);

	// Allocate L_factors buffer
	afree(c->L_factors);
	c->L_factors = NULL;
	if (t->L_factors_size > 0) {
		afit(c->L_factors, t->L_factors_size);
		asetlen(c->L_factors, t->L_factors_size);
	}
}

// -----------------------------------------------------------------------------
// Numeric factorization using cached topology.

// Fill K matrix and factorize using precomputed topology.
static void ldl_numeric_factor(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, SolverDistance* sol_dist)
{
	LDL_Topology* t = c->topo;
	int nc = t->node_count;
	if (nc == 0) return;

	// Zero buffers
	memset(c->diag_data, 0, sizeof(c->diag_data));
	memset(c->diag_D, 0, sizeof(c->diag_D));
	if (c->L_factors) memset(c->L_factors, 0, t->L_factors_size * sizeof(float));

	// Fill diagonal K blocks: each constraint writes at its bundle_offset within its bundle's block.
	// Intra-bundle coupling between constraints sharing both bodies is also written here.
	int bc = c->bundle_count;
	for (int bi = 0; bi < bc; bi++) {
		LDL_Bundle* bun = &c->bundles[bi];
		int bdk = bun->dof;
		for (int ci = 0; ci < bun->count; ci++) {
			LDL_Constraint* con = &c->constraints[bun->start + ci];
			int off = con->bundle_offset;
			int dk = con->dof;
			// Compute constraint's own K into a temp buffer, then scatter into bundle diagonal
			float sub_K[36];
			if (con->is_synthetic) {
				BodyHot* a = ldl_get_body(w, c, con->body_a);
				BodyHot* b = ldl_get_body(w, c, con->body_b);
				ldl_weld_K(a, b, sub_K);
			} else if (con->type == JOINT_BALL_SOCKET) {
				SolverBallSocket* s = &sol_bs[con->solver_idx];
				BodyHot* a = ldl_get_body(w, c, con->body_a);
				BodyHot* b = ldl_get_body(w, c, con->body_b);
				ldl_ball_socket_K(a, b, s->r_a, s->r_b, s->softness, sub_K);
			} else {
				SolverDistance* s = &sol_dist[con->solver_idx];
				BodyHot* a = ldl_get_body(w, c, con->body_a);
				BodyHot* b = ldl_get_body(w, c, con->body_b);
				float inv_mass_sum = a->inv_mass + b->inv_mass;
				float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s->r_a, s->axis)), s->r_a), s->axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s->r_b, s->axis)), s->r_b), s->axis);
				k += s->softness;
				sub_K[0] = k;
			}
			// Trace-based regularization for genuine hard constraints.
			// reg = compliance * trace(K) / dim. Virtual welds skip (already PD).
			if (!con->is_synthetic) {
				float soft = con->type == JOINT_BALL_SOCKET ? sol_bs[con->solver_idx].softness : sol_dist[con->solver_idx].softness;
				if (soft == 0.0f) {
					float trace = 0;
					for (int d = 0; d < dk; d++) trace += sub_K[d * dk + d];
					float reg = LDL_COMPLIANCE * trace / dk;
					for (int d = 0; d < dk; d++) sub_K[d * dk + d] += reg;
				}
			}
			for (int r = 0; r < dk; r++) {
				for (int col = 0; col < dk; col++) {
					c->diag_data[bi][(off + r) * bdk + (off + col)] += sub_K[r * dk + col];
				}
			}
		}
		// Intra-bundle coupling: pairs of constraints within same bundle share both bodies.
		for (int ci = 0; ci < bun->count; ci++) {
			if (c->constraints[bun->start + ci].is_synthetic) continue;
			for (int cj = ci + 1; cj < bun->count; cj++) {
				LDL_Constraint* con_i = &c->constraints[bun->start + ci];
				LDL_Constraint* con_j = &c->constraints[bun->start + cj];
				if (con_j->is_synthetic) continue;
				int oi = con_i->bundle_offset, oj = con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				// Both bodies contribute coupling
				for (int body_side = 0; body_side < 2; body_side++) {
					int body = body_side == 0 ? bun->body_a : bun->body_b;
					BodyHot* B = ldl_get_body(w, c, body);
					float s_i = (body == con_i->body_b) ? 1.0f : -1.0f;
					float s_j = (body == con_j->body_b) ? 1.0f : -1.0f;
					float sign = s_i * s_j;
					v3 ri = V3(0,0,0), rj = V3(0,0,0);
					if (!con_i->is_synthetic) {
						if (con_i->type == JOINT_BALL_SOCKET) ri = (body == con_i->body_b) ? sol_bs[con_i->solver_idx].r_b : sol_bs[con_i->solver_idx].r_a;
						else ri = (body == con_i->body_b) ? sol_dist[con_i->solver_idx].r_b : sol_dist[con_i->solver_idx].r_a;
					}
					if (!con_j->is_synthetic) {
						if (con_j->type == JOINT_BALL_SOCKET) rj = (body == con_j->body_b) ? sol_bs[con_j->solver_idx].r_b : sol_bs[con_j->solver_idx].r_a;
						else rj = (body == con_j->body_b) ? sol_dist[con_j->solver_idx].r_b : sol_dist[con_j->solver_idx].r_a;
					}
					float buf[36];
					if (con_i->type == JOINT_BALL_SOCKET && con_j->type == JOINT_BALL_SOCKET) {
						ldl_compute_off_diag_bs_bs(B, ri, rj, sign, buf);
						for (int r = 0; r < di; r++) for (int col = 0; col < dj; col++) c->diag_data[bi][(oi+r)*bdk + (oj+col)] += buf[r*dj + col];
						float buf_T[36]; block_transpose(buf, di, dj, buf_T);
						for (int r = 0; r < dj; r++) for (int col = 0; col < di; col++) c->diag_data[bi][(oj+r)*bdk + (oi+col)] += buf_T[r*di + col];
					} else if (con_i->type == JOINT_BALL_SOCKET && con_j->type == JOINT_DISTANCE) {
						ldl_compute_off_diag_bs_dist(B, ri, rj, sol_dist[con_j->solver_idx].axis, sign, buf);
						for (int r = 0; r < di; r++) c->diag_data[bi][(oi+r)*bdk + oj] += buf[r];
						for (int r = 0; r < di; r++) c->diag_data[bi][oj*bdk + (oi+r)] += buf[r];
					} else if (con_i->type == JOINT_DISTANCE && con_j->type == JOINT_BALL_SOCKET) {
						ldl_compute_off_diag_bs_dist(B, rj, ri, sol_dist[con_i->solver_idx].axis, sign, buf);
						for (int r = 0; r < dj; r++) c->diag_data[bi][(oj+r)*bdk + oi] += buf[r];
						for (int r = 0; r < dj; r++) c->diag_data[bi][oi*bdk + (oj+r)] += buf[r];
					} else {
						float val = ldl_compute_off_diag_dist_dist(B, ri, rj, sol_dist[con_i->solver_idx].axis, sol_dist[con_j->solver_idx].axis, sign);
						c->diag_data[bi][oi*bdk + oj] += val;
						c->diag_data[bi][oj*bdk + oi] += val;
					}
				}
			}
		}
	}

	// Fill inter-bundle off-diagonal K blocks via precomputed couplings.
	// Each coupling identifies a shared body between two bundles. We iterate all
	// constraint pairs across the two bundles that touch that body.
	int fill_count = asize(t->couplings);
	for (int fi = 0; fi < fill_count; fi++) {
		LDL_Coupling* f = &t->couplings[fi];
		int bun_i = f->block_i, bun_j = f->block_j;
		LDL_Bundle* bi_b = &c->bundles[bun_i], *bj_b = &c->bundles[bun_j];
		int bdk_i = bi_b->dof, bdk_j = bj_b->dof;
		BodyHot* B = ldl_get_body(w, c, f->body);
		float* edge_ij = &c->L_factors[f->L_offset_ij];
		float* edge_ji = &c->L_factors[f->L_offset_ji];

		for (int ci = 0; ci < bi_b->count; ci++) {
			LDL_Constraint* con_i = &c->constraints[bi_b->start + ci];
			if (con_i->body_a != f->body && con_i->body_b != f->body) continue;
			float s_i = (f->body == con_i->body_b) ? 1.0f : -1.0f;
			v3 ri = V3(0,0,0);
			if (!con_i->is_synthetic) {
				if (con_i->type == JOINT_BALL_SOCKET) ri = (f->body == con_i->body_b) ? sol_bs[con_i->solver_idx].r_b : sol_bs[con_i->solver_idx].r_a;
				else ri = (f->body == con_i->body_b) ? sol_dist[con_i->solver_idx].r_b : sol_dist[con_i->solver_idx].r_a;
			}
			int is_weld_i = con_i->is_synthetic;
			for (int cj = 0; cj < bj_b->count; cj++) {
				LDL_Constraint* con_j = &c->constraints[bj_b->start + cj];
				if (con_j->body_a != f->body && con_j->body_b != f->body) continue;
				float s_j = (f->body == con_j->body_b) ? 1.0f : -1.0f;
				float sign = s_i * s_j;
				v3 rj = V3(0,0,0);
				if (!con_j->is_synthetic) {
					if (con_j->type == JOINT_BALL_SOCKET) rj = (f->body == con_j->body_b) ? sol_bs[con_j->solver_idx].r_b : sol_bs[con_j->solver_idx].r_a;
					else rj = (f->body == con_j->body_b) ? sol_dist[con_j->solver_idx].r_b : sol_dist[con_j->solver_idx].r_a;
				}
				int is_weld_j = con_j->is_synthetic;
				int oi = con_i->bundle_offset, oj = con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				float buf[36];
				// Dispatch by constraint type pair (including weld)
				if (is_weld_i && is_weld_j) {
					ldl_compute_off_diag_weld_weld(B, sign, buf);
				} else if (is_weld_i && !is_weld_j) {
					// weld(i) x constraint(j): compute j x weld then transpose
					if (con_j->type == JOINT_BALL_SOCKET) {
						float tmp[36]; ldl_compute_off_diag_bs_weld(B, rj, sign, tmp);
						block_transpose(tmp, dj, di, buf); // dj x di -> di x dj
					} else {
						float tmp[36]; ldl_compute_off_diag_dist_weld(B, rj, sol_dist[con_j->solver_idx].axis, sign, tmp);
						block_transpose(tmp, dj, di, buf);
					}
				} else if (!is_weld_i && is_weld_j) {
					if (con_i->type == JOINT_BALL_SOCKET) {
						ldl_compute_off_diag_bs_weld(B, ri, sign, buf);
					} else {
						ldl_compute_off_diag_dist_weld(B, ri, sol_dist[con_i->solver_idx].axis, sign, buf);
					}
				} else if (con_i->type == JOINT_BALL_SOCKET && con_j->type == JOINT_BALL_SOCKET) {
					ldl_compute_off_diag_bs_bs(B, ri, rj, sign, buf);
				} else if (con_i->type == JOINT_BALL_SOCKET && con_j->type == JOINT_DISTANCE) {
					ldl_compute_off_diag_bs_dist(B, ri, rj, sol_dist[con_j->solver_idx].axis, sign, buf);
				} else if (con_i->type == JOINT_DISTANCE && con_j->type == JOINT_BALL_SOCKET) {
					ldl_compute_off_diag_bs_dist(B, rj, ri, sol_dist[con_i->solver_idx].axis, sign, buf);
					// swap: this was computed as dj x di, need di x dj
					float tmp[36]; memcpy(tmp, buf, di * dj * sizeof(float));
					block_transpose(tmp, dj, di, buf);
				} else {
					float val = ldl_compute_off_diag_dist_dist(B, ri, rj, sol_dist[con_i->solver_idx].axis, sol_dist[con_j->solver_idx].axis, sign);
					buf[0] = val;
				}
				// Write to both directions
				for (int r = 0; r < di; r++) {
					for (int col = 0; col < dj; col++) {
						edge_ij[(oi+r)*bdk_j + (oj+col)] += buf[r*dj + col];
					}
				}
				float buf_T[36];
				block_transpose(buf, di, dj, buf_T);
				for (int r = 0; r < dj; r++) {
					for (int col = 0; col < di; col++) {
						edge_ji[(oj+r)*bdk_i + (oi+col)] += buf_T[r*di + col];
					}
				}
			}
		}
	}

	// Debug: capture dense A matrix before factorization
	int dbg = g_ldl_debug_enabled && t->n <= LDL_MAX_DOF && nc <= 64;
	if (dbg) {
		int n = t->n;
		memset(g_ldl_debug_info.A, 0, n * n * sizeof(float));
		for (int i = 0; i < nc; i++) {
			int di = t->dof[i], oi = t->row_offset[i];
			for (int r = 0; r < di; r++) {
				for (int col = 0; col < di; col++) {
					g_ldl_debug_info.A[(oi+r)*n + oi+col] = c->diag_data[i][r*di + col];
				}
			}
		}
		int fc = asize(t->couplings);
		for (int fi = 0; fi < fc; fi++) {
			LDL_Coupling* f = &t->couplings[fi];
			int bi = f->block_i, bj = f->block_j;
			int di = t->dof[bi], dj = t->dof[bj];
			int oi = t->row_offset[bi], oj = t->row_offset[bj];
			float* eij = &c->L_factors[f->L_offset_ij];
			float* eji = &c->L_factors[f->L_offset_ji];
			for (int r = 0; r < di; r++) {
				for (int col = 0; col < dj; col++) {
					g_ldl_debug_info.A[(oi+r)*n + oj+col] = eij[r*dj + col];
				}
			}
			for (int r = 0; r < dj; r++) {
				for (int col = 0; col < di; col++) {
					g_ldl_debug_info.A[(oj+r)*n + oi+col] = eji[r*di + col];
				}
			}
		}
	}

	// Factorize using precomputed pivot sequence
	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &t->pivots[step];
		int k = pv->node;
		int dk = pv->dk;

		// Factor diagonal block
		float Dk[6];
		block_ldl(c->diag_data[k], Dk, dk);
		for (int d = 0; d < dk; d++) c->diag_D[k][d] = Dk[d];

		// Compute L blocks for non-eliminated neighbors
		for (int ei = 0; ei < pv->col_count; ei++) {
			LDL_Column* en = &t->columns[pv->col_start + ei];
			int di = en->dn;
			float* Eik = &c->L_factors[en->L_offset]; // currently holds K_{i,k}, will be overwritten with L_{i,k}

			// L_{i,k} = E_{i,k} * N_k^{-1}: solve N_k * X = E^T column by column
			float Lik[36];
			for (int col = 0; col < di; col++) {
				float rhs[6], sol[6];
				for (int r = 0; r < dk; r++) rhs[r] = Eik[col * dk + r];
				block_solve(c->diag_data[k], Dk, rhs, sol, dk);
				for (int r = 0; r < dk; r++) Lik[col * dk + r] = sol[r];
			}
			memcpy(Eik, Lik, di * dk * sizeof(float));
		}

		// Reconstruct N_k for Schur complement
		float Nk[36];
		{
			float Lkk[36], LDk[36], LkkT[36];
			memset(Lkk, 0, dk * dk * sizeof(float));
			for (int r = 0; r < dk; r++) {
				Lkk[r * dk + r] = 1.0f;
				for (int cc = 0; cc < r; cc++) Lkk[r * dk + cc] = c->diag_data[k][r * dk + cc];
			}
			for (int r = 0; r < dk; r++) for (int cc = 0; cc < dk; cc++) LDk[r * dk + cc] = Lkk[r * dk + cc] * Dk[cc];
			block_transpose(Lkk, dk, dk, LkkT);
			block_mul(LDk, dk, dk, LkkT, dk, Nk);
		}

		// Apply Schur complement updates
		for (int si = 0; si < pv->schur_count; si++) {
			LDL_Schur* op = &t->schurs[pv->schur_start + si];
			float* Lik = &c->L_factors[op->Lik_offset];
			float* Ljk = &c->L_factors[op->Ljk_offset];
			int di = op->di, dj = op->dj;

			// product = L_{i,k} * N_k * L_{j,k}^T
			float LikNk[36], LjkT[36], product[36];
			block_mul(Lik, di, dk, Nk, dk, LikNk);
			block_transpose(Ljk, dj, dk, LjkT);
			block_mul(LikNk, di, dk, LjkT, dj, product);

			if (op->target_offset < 0) {
				// Diagonal update
				block_sub(c->diag_data[op->target_node], product, di, di);
			} else {
				// Off-diagonal update: write to (i,j) and transpose to (j,i)
				block_sub(&c->L_factors[op->target_offset], product, di, dj);
				float product_T[36];
				block_transpose(product, di, dj, product_T);
				block_sub(&c->L_factors[op->target_offset_rev], product_T, dj, di);
			}
		}
	}

	// Debug: capture D pivots and block info after factorization
	if (dbg) {
		g_ldl_debug_info.n = t->n;
		g_ldl_debug_info.joint_count = nc; // node (bundle) count
		g_ldl_debug_info.bs_count = c->joint_count; // total constraint count
		g_ldl_debug_info.dist_count = 0;
		for (int i = 0; i < nc; i++) {
			g_ldl_debug_info.block_dofs[i] = t->dof[i];
			g_ldl_debug_info.block_rows[i] = t->row_offset[i];
			g_ldl_debug_info.block_types[i] = c->bundles[i].count == 1 ? c->constraints[c->bundles[i].start].type : -1;
		}
		for (int i = 0; i < nc; i++) {
			int di = t->dof[i], oi = t->row_offset[i];
			for (int d = 0; d < di; d++) {
				g_ldl_debug_info.D[oi + d] = c->diag_D[i][d];
			}
		}
	}
}

// Forward/diagonal/back substitution using precomputed topology.
static void ldl_solve_topo(LDL_Topology* t, float diag_data[][36], float diag_D[][6], float* L_factors, float* rhs, float* x)
{
	int nc = t->node_count;
	memcpy(x, rhs, t->n * sizeof(float));

	// Forward substitution
	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &t->pivots[step];
		int ok = pv->ok, dk = pv->dk;
		for (int fi = 0; fi < pv->fwd_count; fi++) {
			LDL_Neighbor* sn = &t->fwd_neighbors[pv->fwd_start + fi];
			float* Lkj = &L_factors[sn->L_offset]; // L_{k,j}: dk x dj block
			int dj = sn->dn, oj = sn->on;
			for (int r = 0; r < dk; r++) {
				for (int c = 0; c < dj; c++) {
					x[ok + r] -= Lkj[r * dj + c] * x[oj + c];
				}
			}
		}
	}

	// Diagonal solve
	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &t->pivots[step];
		int k = pv->node, ok = pv->ok, dk = pv->dk;
		float tmp[6];
		for (int d = 0; d < dk; d++) tmp[d] = x[ok + d];
		block_solve(diag_data[k], diag_D[k], tmp, &x[ok], dk);
	}

	// Back substitution
	for (int step = nc - 1; step >= 0; step--) {
		LDL_Pivot* pv = &t->pivots[step];
		int ok = pv->ok, dk = pv->dk;
		for (int bi = 0; bi < pv->back_count; bi++) {
			LDL_Neighbor* sn = &t->back_neighbors[pv->back_start + bi];
			float* Ljk = &L_factors[sn->L_offset]; // L_{j,k}: dj x dk block. Apply L^T.
			int dj = sn->dn, oj = sn->on;
			for (int r = 0; r < dk; r++) {
				for (int c = 0; c < dj; c++) {
					x[ok + r] -= Ljk[c * dk + r] * x[oj + c];
				}
			}
		}
	}
}

// -----------------------------------------------------------------------------
// Cache management.

static void ldl_cache_free(LDL_Cache* c)
{
	afree(c->constraints);
	afree(c->bundles);
	if (c->topo) { ldl_topology_free(c->topo); CK_FREE(c->topo); }
	afree(c->L_factors);
	afree(c->virtual_bodies);
	afree(c->body_remap);
	memset(c, 0, sizeof(*c));
}

// Free numeric data but keep topology and constraint list for sleep cache.
// When the island wakes, topology is reused if topo_version still matches.
static void ldl_cache_sleep(LDL_Cache* c)
{
	afree(c->L_factors);
	c->L_factors = NULL;
	afree(c->virtual_bodies);
	c->virtual_bodies = NULL;
	c->virtual_body_count = 0;
}

// Build blocks for an island's joints.
static int ldl_cache_rebuild_blocks(LDL_Cache* c, WorldInternal* w, int island_idx, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count)
{
	afree(c->constraints);
	afree(c->bundles);
	afree(c->virtual_bodies);
	afree(c->body_remap);
	c->constraints = NULL;
	c->bundles = NULL;
	c->virtual_bodies = NULL;
	c->body_remap = NULL;
	c->virtual_body_count = 0;
	c->bundle_count = 0;
	c->n = 0;
	c->joint_count = 0;

	Island* isl = &w->islands[island_idx];
	int ji = isl->head_joint;
	while (ji >= 0) {
		JointInternal* j = &w->joints[ji];
		if (j->type == JOINT_BALL_SOCKET) {
			for (int i = 0; i < bs_count; i++) {
				if (sol_bs[i].joint_idx == ji) {
					if (sol_bs[i].softness != 0.0f) break; // soft joints handled by PGS only
					LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b, .solver_idx = i };
					apush(c->constraints, con);
					c->n += 3;
					c->joint_count++;
					break;
				}
			}
		} else if (j->type == JOINT_DISTANCE) {
			for (int i = 0; i < dist_count; i++) {
				if (sol_dist[i].joint_idx == ji) {
					if (sol_dist[i].softness != 0.0f) break; // soft joints handled by PGS only
					LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b, .solver_idx = i };
					apush(c->constraints, con);
					c->n += 1;
					c->joint_count++;
					break;
				}
			}
		}
		ji = j->island_next;
	}

	return c->joint_count;
}

// Group constraints by body pair into bundles. Each bundle becomes one graph node.
static void ldl_build_bundles(LDL_Cache* c)
{
	afree(c->bundles);
	c->bundles = NULL;
	c->bundle_count = 0;

	int jc = c->joint_count;
	if (jc == 0) return;

	// Sort constraints by body pair using insertion sort (small N, stable)
	for (int i = 1; i < jc; i++) {
		LDL_Constraint tmp = c->constraints[i];
		int j = i - 1;
		while (j >= 0 && (c->constraints[j].body_a > tmp.body_a || (c->constraints[j].body_a == tmp.body_a && c->constraints[j].body_b > tmp.body_b))) {
			c->constraints[j + 1] = c->constraints[j];
			j--;
		}
		c->constraints[j + 1] = tmp;
	}

	// Group consecutive constraints with same body pair.
	// Split bundles that would exceed 6 DOF to stay within diag_data[i][36] (max 6x6).
	int start = 0;
	while (start < jc) {
		int ba = c->constraints[start].body_a, bb = c->constraints[start].body_b;
		int end = start + 1;
		while (end < jc && c->constraints[end].body_a == ba && c->constraints[end].body_b == bb) end++;
		// Emit one or more bundles, splitting if cumulative DOF > 6
		int bstart = start;
		while (bstart < end) {
			int dof = 0;
			int bend = bstart;
			while (bend < end) {
				int next_dof = dof + c->constraints[bend].dof;
				if (next_dof > 6 && dof > 0) break; // split here (but always take at least one)
				dof = next_dof;
				bend++;
			}
			for (int i = bstart; i < bend; i++) {
				c->constraints[i].bundle_idx = c->bundle_count;
				c->constraints[i].bundle_offset = (i == bstart) ? 0 : c->constraints[i-1].bundle_offset + c->constraints[i-1].dof;
			}
			LDL_Bundle bundle = { .body_a = ba, .body_b = bb, .dof = dof, .start = bstart, .count = bend - bstart };
			apush(c->bundles, bundle);
			c->bundle_count++;
			bstart = bend;
		}
		start = end;
	}
}

// Solve: build RHS, forward/back-sub via topology, apply lambdas.
static void ldl_island_solve(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count)
{
	int jc = c->joint_count;
	int n = c->n;
	LDL_Topology* t = c->topo;

	// Refresh virtual body state from real hub bodies (they may have changed)
	if (c->virtual_body_count > 0) {
		int real_count = asize(w->body_hot);
		for (int b = 0; b < real_count; b++) {
			if (c->body_remap && c->body_remap[b] >= 0) {
				int first_v = c->body_remap[b] - real_count;
				for (int v = first_v; v < c->virtual_body_count; v++) {
					BodyHot* vb = &c->virtual_bodies[v];
					float S = vb->inv_mass / (w->body_hot[b].inv_mass > 0 ? w->body_hot[b].inv_mass : 1.0f);
					if (S < 1.5f) break;
					vb->position = w->body_hot[b].position;
					vb->rotation = w->body_hot[b].rotation;
					vb->velocity = w->body_hot[b].velocity;
					vb->angular_velocity = w->body_hot[b].angular_velocity;
				}
			}
		}
	}

	// Delta-correction approach (replaces the earlier undo+recompute approach).
	//
	// Why: the old approach undid PGS warm-start impulses, then had LDL recompute
	// exact lambdas from scratch. This required undo and apply to operate on the
	// same bodies. With shattering, PGS applies impulses to real bodies but LDL
	// operates on virtual shards -- the body sets don't match, causing energy
	// injection that makes the system explode.
	//
	// Delta-correction avoids this entirely: PGS runs first with its own warm-
	// starting. LDL then measures whatever residual error remains and computes
	// a small correction delta. No undo needed, no real/virtual body mismatch.
	float* rhs = CK_ALLOC(n * sizeof(float));
	float* lambda = CK_ALLOC(n * sizeof(float));
	memset(rhs, 0, n * sizeof(float));
	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
		if (con->is_synthetic) {
			BodyHot* a = ldl_get_body(w, c, con->body_a);
			BodyHot* b = ldl_get_body(w, c, con->body_b);
			v3 dv = sub(b->velocity, a->velocity);
			v3 dw = sub(b->angular_velocity, a->angular_velocity);
			rhs[oi] = -dv.x; rhs[oi+1] = -dv.y; rhs[oi+2] = -dv.z;
			rhs[oi+3] = -dw.x; rhs[oi+4] = -dw.y; rhs[oi+5] = -dw.z;
		} else if (con->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[con->solver_idx];
			BodyHot* a = ldl_get_body(w, c, con->body_a);
			BodyHot* b = ldl_get_body(w, c, con->body_b);
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			// Clamp warm-started lambda used in regularization term. If the joint was
			// stretched far, lambda may have grown huge from repeated impulses without
			// position recovery. Using a huge lambda in the RHS produces garbage output.
			v3 lam = s->lambda;
			float lam_mag = len(lam);
			if (lam_mag > 1e4f) lam = scale(lam, 1e4f / lam_mag);
			rhs[oi]   = -dv.x - LDL_COMPLIANCE * lam.x;
			rhs[oi+1] = -dv.y - LDL_COMPLIANCE * lam.y;
			rhs[oi+2] = -dv.z - LDL_COMPLIANCE * lam.z;
		} else if (con->type == JOINT_DISTANCE) {
			SolverDistance* s = &sol_dist[con->solver_idx];
			BodyHot* a = ldl_get_body(w, c, con->body_a);
			BodyHot* b = ldl_get_body(w, c, con->body_b);
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float lam = s->lambda;
			if (lam > 1e4f) lam = 1e4f; else if (lam < -1e4f) lam = -1e4f;
			rhs[oi] = -dot(dv, s->axis) - LDL_COMPLIANCE * lam;
		}
	}

	ldl_solve_topo(t, c->diag_data, c->diag_D, c->L_factors, rhs, lambda);

	// Safety: if the solve produced NaN/inf (e.g. from a severely overconstrained
	// or rank-deficient system), discard the correction to prevent simulation blow-up.
	int solve_valid = 1;
	for (int i = 0; i < n; i++) {
		if (!(lambda[i] == lambda[i]) || lambda[i] > 1e15f || lambda[i] < -1e15f) {
			solve_valid = 0;
			break;
		}
	}
	if (!solve_valid) {
		CK_FREE(rhs);
		CK_FREE(lambda);
		return;
	}

	// Debug
	int dbg = g_ldl_debug_enabled && n <= LDL_MAX_DOF;
	if (dbg) { memcpy(g_ldl_debug_info.lambda_ldl, lambda, n * sizeof(float)); g_ldl_debug_info.valid = 1; }

	// Apply delta correction to REAL bodies directly, bypassing virtual shards.
	// Synthetic weld deltas are skipped -- they only exist to couple shards in the
	// factorization and have no physical meaning for the real bodies.
	// s->lambda accumulates both PGS + LDL for next frame's warm-start.
	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		if (con->is_synthetic) continue;
		int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
		if (con->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[con->solver_idx];
			v3 delta = V3(lambda[oi], lambda[oi+1], lambda[oi+2]);
			s->lambda = add(s->lambda, delta);
			apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, delta);
		} else if (con->type == JOINT_DISTANCE) {
			SolverDistance* s = &sol_dist[con->solver_idx];
			float delta = lambda[oi];
			s->lambda += delta;
			apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, scale(s->axis, delta));
		}
	}

	CK_FREE(rhs);
	CK_FREE(lambda);
}

// -----------------------------------------------------------------------------
// Top-level entry point.

static void ldl_correction(WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, int sub)
{
	if (bs_count == 0 && dist_count == 0) return;

	int island_count = asize(w->islands);
	for (int ii = 0; ii < island_count; ii++) {
		if (!(w->island_gen[ii] & 1)) continue;
		Island* isl = &w->islands[ii];
		if (!isl->awake) continue;
		if (isl->joint_count == 0) continue;

		LDL_Cache* c = &isl->ldl;

		// Rebuild blocks on first substep or topology change
		if (sub == 0 || c->topo_version != w->ldl_topo_version) {
			ldl_cache_rebuild_blocks(c, w, ii, sol_bs, bs_count, sol_dist, dist_count);
			ldl_apply_shattering(c, w);
			ldl_build_bundles(c);
			ldl_build_topology(c, w);
			c->topo_version = w->ldl_topo_version;
		}
		if (c->n == 0) continue;

		// Enable debug capture only for the inspected island
		g_ldl_debug_enabled = (ii == g_ldl_debug_island);

		// Re-allocate L_factors if freed by sleep cache
		if (!c->L_factors && c->topo && c->topo->L_factors_size > 0) {
			afit(c->L_factors, c->topo->L_factors_size);
			asetlen(c->L_factors, c->topo->L_factors_size);
		}

		// Numeric factorization every substep (lever arms change after integrate_positions)
		ldl_numeric_factor(c, w, sol_bs, sol_dist);

		// Solve every substep
		ldl_island_solve(c, w, sol_bs, bs_count, sol_dist, dist_count);
	}
}
