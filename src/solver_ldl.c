// solver_ldl.c -- Direct LDL joint correction for equality constraints.
// Reference: Mizerski, "Improving an Iterative Physics Solver Using a Direct Method", GDC 2020.
//
// Solves K * lambda = b where K = J * M^-1 * J^T (constraint-space effective mass),
// lambda = constraint impulses, b = velocity/position error. K is symmetric positive
// definite (regularized) with block structure from the constraint graph.

// Double-precision vector/quat types for LDL internals.
typedef struct dv3 { double x, y, z; } dv3;
typedef struct dquat { double x, y, z, w; } dquat;
static inline dv3 DV3(double x, double y, double z) { return (dv3){x, y, z}; }
static inline dv3 dv3_from_v3(v3 v) { return (dv3){v.x, v.y, v.z}; }
static inline v3 dv3_to_v3(dv3 v) { return V3((float)v.x, (float)v.y, (float)v.z); }
static inline dv3 dv3_add(dv3 a, dv3 b) { return DV3(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline dv3 dv3_sub(dv3 a, dv3 b) { return DV3(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline dv3 dv3_scale(dv3 v, double s) { return DV3(v.x*s, v.y*s, v.z*s); }
static inline dv3 dv3_cross(dv3 a, dv3 b) { return DV3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
static inline double dv3_dot(dv3 a, dv3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double dv3_len(dv3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }

// Double-precision inv_inertia_mul: R * diag(inv_I) * R^T * v.
static dv3 dinv_inertia_mul(quat rot, v3 inv_i, dv3 v)
{
	// rotate_inv(q, v): R^T * v in double
	double qx = rot.x, qy = rot.y, qz = rot.z, qw = rot.w;
	double vx = v.x, vy = v.y, vz = v.z;
	double tx = 2.0*(qy*vz - qz*vy), ty = 2.0*(qz*vx - qx*vz), tz = 2.0*(qx*vy - qy*vx);
	// inv rotation: v - w*T + u x T (differs from forward only in sign of w*T)
	double rx = vx - qw*tx + (qy*tz - qz*ty);
	double ry = vy - qw*ty + (qz*tx - qx*tz);
	double rz = vz - qw*tz + (qx*ty - qy*tx);
	// scale by diagonal inertia
	double sx = rx * (double)inv_i.x, sy = ry * (double)inv_i.y, sz = rz * (double)inv_i.z;
	// rotate(q, scaled): R * scaled in double
	double tx2 = 2.0*(qy*sz - qz*sy), ty2 = 2.0*(qz*sx - qx*sz), tz2 = 2.0*(qx*sy - qy*sx);
	return DV3(sx + qw*tx2 + (qy*tz2 - qz*ty2), sy + qw*ty2 + (qz*tx2 - qx*tz2), sz + qw*tz2 + (qx*ty2 - qy*tx2));
}

static LDL_DebugInfo g_ldl_debug_info;
int g_ldl_debug_enabled;
int g_ldl_debug_island = -1; // which island to capture debug data for (-1 = none)
int g_ldl_trace_solve; // set to 1 to print detailed solve info for next LDL call

#define SHATTER_THRESHOLD 6 // DOF threshold: bodies with more DOF than this get shattered into virtual shards
#define SHARD_TARGET      6 // target DOF per shard
#define LDL_MIN_COMPLIANCE 5e-5 // min regularization floor: compliance * trace(K) / dim on K diagonal + compliance*lambda in RHS
#define LDL_SYNTH_COMPLIANCE 1e-3 // synthetic weld compliance: absolute (not trace-scaled) to break circular chain singularity

// -----------------------------------------------------------------------------
// Small block math helpers.
// Diagonal blocks use packed lower-triangular storage: n*(n+1)/2 elements.
// Element (r,c) with r >= c is at index r*(r+1)/2 + c. Symmetric by construction.

// Packed lower-triangular index. Always returns the lower-triangle position.
#define LDL_TRI(r, c) ((r) >= (c) ? (r)*((r)+1)/2 + (c) : (c)*((c)+1)/2 + (r))

// In-place NxN LDL^T factorize on packed lower-triangular storage.
// L stored in-place (unit diagonal implicit), D separate.
// Packed storage guarantees symmetry — no need to symmetrize.
// Block math operates in double precision for numerical stability.
// Returns 0 on success, -1 if a non-positive pivot is encountered (caller should skip solve).
static int block_ldl(double* A, double* D, int n)
{
	for (int j = 0; j < n; j++) {
		double dj = A[LDL_TRI(j,j)];
		for (int k = 0; k < j; k++) dj -= A[LDL_TRI(j,k)] * A[LDL_TRI(j,k)] * D[k];
		if (dj <= 0.0) {
			// Non-positive pivot: regularization bug or Schur complement issue.
			// Log and bail — caller skips LDL for this island.
			printf("[LDL] non-positive pivot dj=%.6e at j=%d/%d\n", dj, j, n);
			return -1;
		}
		D[j] = dj;
		double inv_dj = 1.0 / D[j];
		for (int i = j + 1; i < n; i++) {
			double lij = A[LDL_TRI(i,j)];
			for (int k = 0; k < j; k++) lij -= A[LDL_TRI(i,k)] * A[LDL_TRI(j,k)] * D[k];
			A[LDL_TRI(i,j)] = lij * inv_dj;
		}
	}
	return 0;
}

static void block_solve(double* L, double* D, double* b, double* x, int n)
{
	for (int i = 0; i < n; i++) { x[i] = b[i]; for (int k = 0; k < i; k++) x[i] -= L[LDL_TRI(i,k)] * x[k]; }
	for (int i = 0; i < n; i++) x[i] /= D[i];
	for (int i = n - 1; i >= 0; i--) for (int k = i + 1; k < n; k++) x[i] -= L[LDL_TRI(k,i)] * x[k];
}

static void block_mul(double* A, int ra, int ca, double* B, int cb, double* out)
{
	for (int i = 0; i < ra; i++) {
		for (int j = 0; j < cb; j++) {
			double sum = 0;
			for (int k = 0; k < ca; k++) sum += A[i*ca + k] * B[k*cb + j];
			out[i*cb + j] = sum;
		}
	}
}

static void block_sub(double* A, double* B, int r, int c)
{
	for (int i = 0; i < r * c; i++) A[i] -= B[i];
}

static void block_transpose(double* A, int r, int c, double* out)
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
static double* ldl_sparse_get_or_create_edge(LDL_Sparse* s, int i, int j)
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
	double* data_ij = &s->adj_data[i][old_size_i];
	memset(data_ij, 0, block_size_ij * sizeof(double));

	// Add i to adj[j]
	int block_size_ji = dj * di;
	int old_size_j = asize(s->adj_data[j]);
	apush(s->adj[j], i);
	afit(s->adj_data[j], old_size_j + block_size_ji);
	asetlen(s->adj_data[j], old_size_j + block_size_ji);
	double* data_ji = &s->adj_data[j][old_size_j];
	memset(data_ji, 0, block_size_ji * sizeof(double));

	return data_ij;
}

// Get existing off-diagonal block data for edge (i, j). Returns NULL if not found.
static double* ldl_sparse_get_edge(LDL_Sparse* s, int i, int j)
{
	int ki = ldl_sparse_find_edge(s, i, j);
	if (ki < 0) return NULL;
	return &s->adj_data[i][ki * s->dof[i] * s->dof[j]];
}

// Apply angular position correction: integrate rotation by angular velocity w for one step.
// Takes dv3 (double) from LDL internals, writes to quat body state.
static void apply_rotation_delta(quat* q, dv3 w)
{
	if (w.x == 0 && w.y == 0 && w.z == 0) return;
	dquat dq_spin = { w.x, w.y, w.z, 0.0 };
	dquat dq_cur = { q->x, q->y, q->z, q->w };
	// dq = spin * q (quaternion multiply in double)
	dquat dq = {
		dq_spin.w*dq_cur.x + dq_spin.x*dq_cur.w + dq_spin.y*dq_cur.z - dq_spin.z*dq_cur.y,
		dq_spin.w*dq_cur.y - dq_spin.x*dq_cur.z + dq_spin.y*dq_cur.w + dq_spin.z*dq_cur.x,
		dq_spin.w*dq_cur.z + dq_spin.x*dq_cur.y - dq_spin.y*dq_cur.x + dq_spin.z*dq_cur.w,
		dq_spin.w*dq_cur.w - dq_spin.x*dq_cur.x - dq_spin.y*dq_cur.y - dq_spin.z*dq_cur.z,
	};
	dq_cur.x += 0.5 * dq.x; dq_cur.y += 0.5 * dq.y;
	dq_cur.z += 0.5 * dq.z; dq_cur.w += 0.5 * dq.w;
	double ql = sqrt(dq_cur.x*dq_cur.x + dq_cur.y*dq_cur.y + dq_cur.z*dq_cur.z + dq_cur.w*dq_cur.w);
	double inv_ql = 1.0 / (ql > 1e-15 ? ql : 1.0);
	q->x = (float)(dq_cur.x * inv_ql); q->y = (float)(dq_cur.y * inv_ql);
	q->z = (float)(dq_cur.z * inv_ql); q->w = (float)(dq_cur.w * inv_ql);
}

// -----------------------------------------------------------------------------
// Jacobian-based K matrix computation.
// Each constraint provides Jacobian rows J_a, J_b (dof x 6 each).
// K = J * M^{-1} * J^T is computed generically for all constraint types.

// Fill Jacobian rows for a constraint. Writes dof rows starting at jac[0].
static void ldl_fill_jacobian(LDL_Constraint* con, SolverBallSocket* sol_bs, SolverDistance* sol_dist, SolverHinge* sol_hinge, LDL_JacobianRow* jac)
{
	int dof = con->dof;
	for (int d = 0; d < dof; d++) { memset(jac[d].J_a, 0, 6 * sizeof(double)); memset(jac[d].J_b, 0, 6 * sizeof(double)); }

	if (con->is_synthetic) {
		// Weld: J_a = -I_6, J_b = +I_6
		for (int d = 0; d < 6; d++) { jac[d].J_a[d] = -1.0f; jac[d].J_b[d] = 1.0f; }
	} else if (con->type == JOINT_BALL_SOCKET) {
		SolverBallSocket* s = &sol_bs[con->solver_idx];
		v3 ra = s->r_a, rb = s->r_b;
		// J_a = [-I_3, skew(r_a)],  J_b = [I_3, -skew(r_b)]
		jac[0].J_a[0] = -1; jac[0].J_a[4] = -ra.z; jac[0].J_a[5] =  ra.y;
		jac[1].J_a[1] = -1; jac[1].J_a[3] =  ra.z; jac[1].J_a[5] = -ra.x;
		jac[2].J_a[2] = -1; jac[2].J_a[3] = -ra.y; jac[2].J_a[4] =  ra.x;
		jac[0].J_b[0] =  1; jac[0].J_b[4] =  rb.z; jac[0].J_b[5] = -rb.y;
		jac[1].J_b[1] =  1; jac[1].J_b[3] = -rb.z; jac[1].J_b[5] =  rb.x;
		jac[2].J_b[2] =  1; jac[2].J_b[3] =  rb.y; jac[2].J_b[4] = -rb.x;
	} else if (con->type == JOINT_DISTANCE) {
		SolverDistance* s = &sol_dist[con->solver_idx];
		v3 ax = s->axis;
		v3 rxa = cross(s->r_a, ax), rxb = cross(s->r_b, ax);
		// J_a = [-axis^T, -(r_a x axis)^T],  J_b = [axis^T, (r_b x axis)^T]
		jac[0].J_a[0] = -ax.x; jac[0].J_a[1] = -ax.y; jac[0].J_a[2] = -ax.z;
		jac[0].J_a[3] = -rxa.x; jac[0].J_a[4] = -rxa.y; jac[0].J_a[5] = -rxa.z;
		jac[0].J_b[0] = ax.x; jac[0].J_b[1] = ax.y; jac[0].J_b[2] = ax.z;
		jac[0].J_b[3] = rxb.x; jac[0].J_b[4] = rxb.y; jac[0].J_b[5] = rxb.z;
	} else { // JOINT_HINGE
		SolverHinge* s = &sol_hinge[con->solver_idx];
		v3 ra = s->r_a, rb = s->r_b;
		// Linear rows 0-2: same as ball socket
		jac[0].J_a[0] = -1; jac[0].J_a[4] = -ra.z; jac[0].J_a[5] =  ra.y;
		jac[1].J_a[1] = -1; jac[1].J_a[3] =  ra.z; jac[1].J_a[5] = -ra.x;
		jac[2].J_a[2] = -1; jac[2].J_a[3] = -ra.y; jac[2].J_a[4] =  ra.x;
		jac[0].J_b[0] =  1; jac[0].J_b[4] =  rb.z; jac[0].J_b[5] = -rb.y;
		jac[1].J_b[1] =  1; jac[1].J_b[3] = -rb.z; jac[1].J_b[5] =  rb.x;
		jac[2].J_b[2] =  1; jac[2].J_b[3] =  rb.y; jac[2].J_b[4] = -rb.x;
		// Angular rows 3-4: J_a = [0, u_d], J_b = [0, -u_d]
		jac[3].J_a[3] = s->u1.x; jac[3].J_a[4] = s->u1.y; jac[3].J_a[5] = s->u1.z;
		jac[3].J_b[3] = -s->u1.x; jac[3].J_b[4] = -s->u1.y; jac[3].J_b[5] = -s->u1.z;
		jac[4].J_a[3] = s->u2.x; jac[4].J_a[4] = s->u2.y; jac[4].J_a[5] = s->u2.z;
		jac[4].J_b[3] = -s->u2.x; jac[4].J_b[4] = -s->u2.y; jac[4].J_b[5] = -s->u2.z;
	}
}

// Accumulate K contribution from one body: K += J * (w * M^{-1}) * J^T into packed lower-tri.
// jac has dof rows, side: 0 = J_a, 1 = J_b. weight scales inv_mass/inv_inertia (shattering).
static void ldl_K_body_contrib(LDL_JacobianRow* jac, int dof, int side, int dof_start, BodyHot* body, double weight, double* K_packed)
{
	double wm = (double)body->inv_mass * (double)weight;
	double W[6 * LDL_MAX_BLOCK_DIM];
	for (int d = 0; d < dof; d++) {
		double* J = side ? jac[d].J_b : jac[d].J_a;
		W[0*dof+d] = wm * J[0];
		W[1*dof+d] = wm * J[1];
		W[2*dof+d] = wm * J[2];
		dv3 j_ang = DV3(J[3], J[4], J[5]);
		dv3 w_ang = dv3_scale(dinv_inertia_mul(body->rotation, body->inv_inertia_local, j_ang), (double)weight);
		W[3*dof+d] = w_ang.x;
		W[4*dof+d] = w_ang.y;
		W[5*dof+d] = w_ang.z;
	}
	for (int r = 0; r < dof; r++)
		for (int c = 0; c <= r; c++) {
			double sum = 0;
			double* Jr = side ? jac[r].J_b : jac[r].J_a;
			for (int k = 0; k < 6; k++) sum += Jr[k] * W[k*dof + c];
			K_packed[LDL_TRI(dof_start+r, dof_start+c)] += sum;
		}
}

static void ldl_K_body_off(LDL_JacobianRow* jac_i, int di, int side_i, LDL_JacobianRow* jac_j, int dj, int side_j, BodyHot* body, double weight, double* out)
{
	double wm = (double)body->inv_mass * (double)weight;
	double W[6 * LDL_MAX_BLOCK_DIM];
	for (int d = 0; d < dj; d++) {
		double* J = side_j ? jac_j[d].J_b : jac_j[d].J_a;
		W[0*dj+d] = wm * J[0];
		W[1*dj+d] = wm * J[1];
		W[2*dj+d] = wm * J[2];
		dv3 j_ang = DV3(J[3], J[4], J[5]);
		dv3 w_ang = dv3_scale(dinv_inertia_mul(body->rotation, body->inv_inertia_local, j_ang), (double)weight);
		W[3*dj+d] = w_ang.x;
		W[4*dj+d] = w_ang.y;
		W[5*dj+d] = w_ang.z;
	}
	for (int r = 0; r < di; r++) {
		double* Jr = side_i ? jac_i[r].J_b : jac_i[r].J_a;
		for (int c = 0; c < dj; c++) {
			double sum = 0;
			for (int k = 0; k < 6; k++) sum += Jr[k] * W[k*dj + c];
			out[r*dj + c] += sum;
		}
	}
}

// Apply impulse generically: v += M^{-1} * J^T * lambda for one body.
// jac has dof rows, lambda has dof scalars. side: 0 = J_a, 1 = J_b.
// sign: +1 or -1 (convention: body_a gets -, body_b gets +... but sign is baked into J).
static void ldl_apply_jacobian_impulse(LDL_JacobianRow* jac, int dof, double* lambda, BodyHot* body, int side)
{
	for (int d = 0; d < dof; d++) {
		double* J = side ? jac[d].J_b : jac[d].J_a;
		double lam = lambda[d];
		double inv_m = (double)body->inv_mass;
		dv3 dv = DV3(inv_m * J[0] * lam, inv_m * J[1] * lam, inv_m * J[2] * lam);
		body->velocity.x += (float)dv.x;
		body->velocity.y += (float)dv.y;
		body->velocity.z += (float)dv.z;
		dv3 j_ang_d = DV3(J[3] * lam, J[4] * lam, J[5] * lam);
		dv3 dw = dinv_inertia_mul(body->rotation, body->inv_inertia_local, j_ang_d);
		body->angular_velocity.x += (float)dw.x;
		body->angular_velocity.y += (float)dw.y;
		body->angular_velocity.z += (float)dw.z;
	}
}

static double ldl_constraint_velocity(LDL_JacobianRow* jac, BodyHot* a, BodyHot* b)
{
	double* Ja = jac->J_a;
	double* Jb = jac->J_b;
	return Ja[0]*a->velocity.x + Ja[1]*a->velocity.y + Ja[2]*a->velocity.z + Ja[3]*a->angular_velocity.x + Ja[4]*a->angular_velocity.y + Ja[5]*a->angular_velocity.z + Jb[0]*b->velocity.x + Jb[1]*b->velocity.y + Jb[2]*b->velocity.z + Jb[3]*b->angular_velocity.x + Jb[4]*b->angular_velocity.y + Jb[5]*b->angular_velocity.z;
}

// -----------------------------------------------------------------------------
// Sparse LDL factorize and solve.

// Count non-eliminated neighbors of node i (current degree in elimination graph).
static int ldl_sparse_degree(LDL_Sparse* s, int i, int* eliminated)
{
	int cnt = asize(s->adj[i]);
	int deg = 0;
	for (int k = 0; k < cnt; k++)
		if (!eliminated[s->adj[i][k]]) deg++;
	return deg;
}

// Minimum degree elimination ordering.
// At each step, eliminates the node with the fewest non-eliminated neighbors.
// O(n^2) total for sparse graphs, O(n * d) per step where d = max degree.
// Produces fill-in edges that connect all remaining neighbors of the eliminated node.
static void ldl_sparse_min_degree_order(LDL_Sparse* s)
{
	int nc = s->node_count;
	int eliminated[LDL_MAX_NODES] = {0};
	int degree[LDL_MAX_NODES];

	for (int i = 0; i < nc; i++)
		degree[i] = ldl_sparse_degree(s, i, eliminated);

	for (int step = 0; step < nc; step++) {
		// Pick node with minimum degree, breaking ties by lowest DOF then lowest index
		int best = -1, best_deg = 0x7FFFFFFF, best_dof = 0x7FFFFFFF;
		for (int i = 0; i < nc; i++) {
			if (eliminated[i]) continue;
			if (degree[i] < best_deg || (degree[i] == best_deg && s->dof[i] < best_dof)) { best_deg = degree[i]; best_dof = s->dof[i]; best = i; }
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

		// Update degree for neighbors of eliminated node
		for (int ii = 0; ii < cnt; ii++) {
			int ni = s->adj[best][ii];
			if (eliminated[ni]) continue;
			degree[ni] = ldl_sparse_degree(s, ni, eliminated);
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

	// Iterative post-order DFS: children before parent preserves min-degree intent.
	// Pre-order would reverse elimination (root first = maximum fill-in = non-positive pivots).
	int new_order[LDL_MAX_NODES];
	int out_idx = 0;
	int stack[LDL_MAX_NODES], sp = 0;
	int child_idx[LDL_MAX_NODES];
	for (int ri = root_count - 1; ri >= 0; ri--) {
		stack[sp] = roots[ri];
		child_idx[sp] = 0;
		sp++;
	}
	while (sp > 0) {
		int step = stack[sp - 1];
		int cc = asize(children[step]);
		if (child_idx[sp - 1] < cc) {
			int ci = child_idx[sp - 1]++;
			stack[sp] = children[step][ci];
			child_idx[sp] = 0;
			sp++;
		} else {
			sp--;
			new_order[out_idx++] = s->elim_order[step];
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

// Apply shattering to the block list. Sets shard weights on constraints and adds
// synthetic welds. No virtual body copies — real bodies are used with weight scaling.
// Virtual shard indices are used only for graph topology (reducing fill-in).
static void ldl_apply_shattering(LDL_Cache* c, WorldInternal* w)
{
	int jc = c->joint_count;
	int body_count = asize(w->body_hot);

	// Count total DOF per body
	int body_dof[4096] = {0};
	assert(body_count < 4096);
	for (int i = 0; i < jc; i++) {
		body_dof[c->constraints[i].body_a] += c->constraints[i].dof;
		body_dof[c->constraints[i].body_b] += c->constraints[i].dof;
	}

	// Map: real body -> first virtual shard index (-1 if not shattered)
	afree(c->body_remap);
	c->body_remap = NULL;
	afit(c->body_remap, body_count);
	asetlen(c->body_remap, body_count);
	for (int b = 0; b < body_count; b++) c->body_remap[b] = -1;

	afree(c->shard_counts);
	c->shard_counts = NULL;
	afit(c->shard_counts, body_count);
	asetlen(c->shard_counts, body_count);
	for (int b = 0; b < body_count; b++) c->shard_counts[b] = 0;

	c->virtual_body_count = 0;

	int virtual_base = body_count;
	typedef struct ShatterInfo { int body; int shard_count; int first_virtual; } ShatterInfo;
	CK_DYNA ShatterInfo* shatter_list = NULL;

	for (int b = 0; b < body_count; b++) {
		if (body_dof[b] <= SHATTER_THRESHOLD) continue;
		if (w->body_hot[b].inv_mass == 0.0f) continue;

		int S = (body_dof[b] + SHARD_TARGET - 1) / SHARD_TARGET;
		if (S < 2) S = 2;

		int first_virt = virtual_base + c->virtual_body_count;
		c->body_remap[b] = first_virt;
		c->shard_counts[b] = S;
		c->virtual_body_count += S; // count only, no body copies

		ShatterInfo si = { .body = b, .shard_count = S, .first_virtual = first_virt };
		apush(shatter_list, si);
	}

	if (asize(shatter_list) == 0) {
		afree(shatter_list);
		return;
	}

	// Greedy bin-packing: assign each constraint to the shard with smallest current DOF load.
	for (int si = 0; si < asize(shatter_list); si++) {
		ShatterInfo* info = &shatter_list[si];
		int hub = info->body;
		int S = info->shard_count;

		CK_DYNA int* touching = NULL;
		for (int i = 0; i < jc; i++) {
			if (c->constraints[i].body_a == hub || c->constraints[i].body_b == hub)
				apush(touching, i);
		}
		int tc = asize(touching);
		for (int i = 1; i < tc; i++) {
			int tmp = touching[i];
			int j = i - 1;
			while (j >= 0 && c->constraints[touching[j]].dof < c->constraints[tmp].dof) { touching[j+1] = touching[j]; j--; }
			touching[j+1] = tmp;
		}

		int shard_load[128] = {0};
		assert(S <= 128);
		for (int ti = 0; ti < tc; ti++) {
			int ci = touching[ti];
			int best = 0;
			for (int s = 1; s < S; s++) { if (shard_load[s] < shard_load[best]) best = s; }
			shard_load[best] += c->constraints[ci].dof;
			// Redirect graph body index to virtual shard, set weight on real body
			if (c->constraints[ci].body_a == hub) {
				c->constraints[ci].body_a = info->first_virtual + best;
				// real_body_a stays as hub (set during rebuild), weight = S
				c->constraints[ci].weight_a = (double)S;
			}
			if (c->constraints[ci].body_b == hub) {
				c->constraints[ci].body_b = info->first_virtual + best;
				c->constraints[ci].weight_b = (double)S;
			}
		}
		afree(touching);
	}

	// Add synthetic weld joints in a wrap-around chain between shards.
	for (int si = 0; si < asize(shatter_list); si++) {
		ShatterInfo* info = &shatter_list[si];
		for (int s = 0; s < info->shard_count; s++) {
			int next = (s + 1) % info->shard_count;
			LDL_Constraint synth = {
				.type = -1, .dof = 6,
				.body_a = info->first_virtual + s,
				.body_b = info->first_virtual + next,
				.real_body_a = info->body,
				.real_body_b = info->body,
				.weight_a = (double)info->shard_count,
				.weight_b = (double)info->shard_count,
				.solver_idx = -1, .is_synthetic = 1,
			};
			apush(c->constraints, synth);
			c->joint_count++;
			c->n += 6;
		}
	}

	afree(shatter_list);
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
	ldl_sparse_min_degree_order(&tmp);
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

// Post-factorization condition number estimate: max(D) / min(D).
static double ldl_condition_check(LDL_Cache* c, LDL_Topology* t)
{
	double minD = 1e30, maxD = 0.0;
	for (int i = 0; i < t->node_count; i++) {
		for (int d = 0; d < t->dof[i]; d++) {
			double val = c->diag_D[i][d];
			if (val < minD) minD = val;
			if (val > maxD) maxD = val;
		}
	}
	return maxD / minD;
}


// -----------------------------------------------------------------------------
// Numeric factorization using cached topology.

// Fill K matrix and factorize using precomputed topology.
static void ldl_numeric_factor(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, SolverDistance* sol_dist, SolverHinge* sol_hinge)
{
	LDL_Topology* t = c->topo;
	int nc = t->node_count;
	if (nc == 0) return;

	// Zero buffers
	memset(c->diag_data, 0, sizeof(c->diag_data));
	memset(c->diag_D, 0, sizeof(c->diag_D));
	if (c->L_factors) memset(c->L_factors, 0, t->L_factors_size * sizeof(double));

	// Compute and store Jacobians for all constraints.
	int jc = c->joint_count;
	afree(c->jacobians);
	c->jacobians = NULL;
	afit(c->jacobians, c->n);
	asetlen(c->jacobians, c->n);
	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		con->jacobian_start = t->row_offset[con->bundle_idx] + con->bundle_offset;
		ldl_fill_jacobian(con, sol_bs, sol_dist, sol_hinge, &c->jacobians[con->jacobian_start]);
	}

	// Fill diagonal K blocks via generic K = J * M^{-1} * J^T.
	int bc = c->bundle_count;
	for (int bi = 0; bi < bc; bi++) {
		LDL_Bundle* bun = &c->bundles[bi];
		for (int ci = 0; ci < bun->count; ci++) {
			LDL_Constraint* con = &c->constraints[bun->start + ci];
			int off = con->bundle_offset;
			int dk = con->dof;
			BodyHot* a = &w->body_hot[con->real_body_a];
			BodyHot* b = &w->body_hot[con->real_body_b];
			LDL_JacobianRow* jac = &c->jacobians[con->jacobian_start];

			// K += J_a * (w_a * M_a^{-1}) * J_a^T + J_b * (w_b * M_b^{-1}) * J_b^T
			ldl_K_body_contrib(jac, dk, 0, off, a, con->weight_a, c->diag_data[bi]);
			ldl_K_body_contrib(jac, dk, 1, off, b, con->weight_b, c->diag_data[bi]);

			// Regularization: soft constraints add softness directly to K diagonal
			// (matching PGS convention). Rigid/synthetic use trace-scaled minimum compliance.
			// Synthetic welds MUST be regularized: a circular chain of S welds on S shards
			// is structurally singular (rank-deficient). Without regularization, the null-space
			// (uniform impulse vector) produces a zero pivot, violating SPD.
			if (con->is_synthetic) {
				// Synthetic welds: absolute compliance to break circular chain singularity.
				for (int d = 0; d < dk; d++) c->diag_data[bi][LDL_TRI(off+d, off+d)] += LDL_SYNTH_COMPLIANCE;
			} else {
				double soft = con->type == JOINT_BALL_SOCKET ? sol_bs[con->solver_idx].softness :
					con->type == JOINT_DISTANCE ? sol_dist[con->solver_idx].softness :
					sol_hinge[con->solver_idx].softness;
				if (soft > 0.0) {
					for (int d = 0; d < dk; d++) c->diag_data[bi][LDL_TRI(off+d, off+d)] += soft;
				} else {
					double trace = 0;
					for (int d = 0; d < dk; d++) trace += c->diag_data[bi][LDL_TRI(off+d, off+d)];
					double reg = LDL_MIN_COMPLIANCE * trace / dk;
					for (int d = 0; d < dk; d++) c->diag_data[bi][LDL_TRI(off+d, off+d)] += reg;
				}
			}
		}
		// Intra-bundle coupling: pairs of constraints within same bundle share both bodies.
		for (int ci = 0; ci < bun->count; ci++) {
			LDL_Constraint* con_i = &c->constraints[bun->start + ci];
			for (int cj = ci + 1; cj < bun->count; cj++) {
				LDL_Constraint* con_j = &c->constraints[bun->start + cj];
				int oi = con_i->bundle_offset, oj = con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				LDL_JacobianRow* jac_i = &c->jacobians[con_i->jacobian_start];
				LDL_JacobianRow* jac_j = &c->jacobians[con_j->jacobian_start];
				// Both bodies contribute coupling
				for (int bs = 0; bs < 2; bs++) {
					int body = bs == 0 ? bun->body_a : bun->body_b;
					int real_b = bs == 0 ? con_i->real_body_a : con_i->real_body_b;
					double wt = (body == con_i->body_a) ? con_i->weight_a : con_i->weight_b;
					BodyHot* B = &w->body_hot[real_b];
					int side_i = (body == con_i->body_b) ? 1 : 0;
					int side_j = (body == con_j->body_b) ? 1 : 0;
					double buf[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM] = {0};
					ldl_K_body_off(jac_i, di, side_i, jac_j, dj, side_j, B, wt, buf);
					for (int r = 0; r < di; r++)
						for (int col = 0; col < dj; col++)
							c->diag_data[bi][LDL_TRI(oi+r, oj+col)] += buf[r*dj + col];
				}
			}
		}
	}

	// Fill inter-bundle off-diagonal K blocks via precomputed couplings.
	// Generic: uses stored Jacobians, no type dispatch needed.
	int fill_count = asize(t->couplings);
	for (int fi = 0; fi < fill_count; fi++) {
		LDL_Coupling* f = &t->couplings[fi];
		int bun_i = f->block_i, bun_j = f->block_j;
		LDL_Bundle* bi_b = &c->bundles[bun_i], *bj_b = &c->bundles[bun_j];
		int bdk_i = bi_b->dof, bdk_j = bj_b->dof;
		double* edge_ij = &c->L_factors[f->L_offset_ij];
		double* edge_ji = &c->L_factors[f->L_offset_ji];

		for (int ci = 0; ci < bi_b->count; ci++) {
			LDL_Constraint* con_i = &c->constraints[bi_b->start + ci];
			if (con_i->body_a != f->body && con_i->body_b != f->body) continue;
			int side_i = (f->body == con_i->body_b) ? 1 : 0;
			int real_b = side_i ? con_i->real_body_b : con_i->real_body_a;
			double wt = side_i ? con_i->weight_b : con_i->weight_a;
			BodyHot* B = &w->body_hot[real_b];
			LDL_JacobianRow* jac_i = &c->jacobians[con_i->jacobian_start];
			for (int cj = 0; cj < bj_b->count; cj++) {
				LDL_Constraint* con_j = &c->constraints[bj_b->start + cj];
				if (con_j->body_a != f->body && con_j->body_b != f->body) continue;
				int side_j = (f->body == con_j->body_b) ? 1 : 0;
				int oi = con_i->bundle_offset, oj = con_j->bundle_offset;
				int di = con_i->dof, dj = con_j->dof;
				LDL_JacobianRow* jac_j = &c->jacobians[con_j->jacobian_start];
				double buf[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM] = {0};
				ldl_K_body_off(jac_i, di, side_i, jac_j, dj, side_j, B, wt, buf);
				// Write to both directions
				for (int r = 0; r < di; r++)
					for (int col = 0; col < dj; col++)
						edge_ij[(oi+r)*bdk_j + (oj+col)] += buf[r*dj + col];
				double buf_T[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
				block_transpose(buf, di, dj, buf_T);
				for (int r = 0; r < dj; r++)
					for (int col = 0; col < di; col++)
						edge_ji[(oj+r)*bdk_i + (oi+col)] += buf_T[r*di + col];
			}
		}
	}

	// Debug: capture dense A matrix before factorization
	int dbg = g_ldl_debug_enabled && t->n <= LDL_MAX_DOF && nc <= 64;
	if (dbg) {
		int n = t->n;
		memset(g_ldl_debug_info.A, 0, n * n * sizeof(g_ldl_debug_info.A[0]));
		for (int i = 0; i < nc; i++) {
			int di = t->dof[i], oi = t->row_offset[i];
			for (int r = 0; r < di; r++) {
				for (int col = 0; col < di; col++) {
					g_ldl_debug_info.A[(oi+r)*n + oi+col] = c->diag_data[i][LDL_TRI(r, col)];
				}
			}
		}
		int fc = asize(t->couplings);
		for (int fi = 0; fi < fc; fi++) {
			LDL_Coupling* f = &t->couplings[fi];
			int bi = f->block_i, bj = f->block_j;
			int di = t->dof[bi], dj = t->dof[bj];
			int oi = t->row_offset[bi], oj = t->row_offset[bj];
			double* eij = &c->L_factors[f->L_offset_ij];
			double* eji = &c->L_factors[f->L_offset_ji];
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

	// Factorize using precomputed pivot sequence.
	// Diagonal blocks are packed lower-triangular — symmetric by construction.
	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &t->pivots[step];
		int k = pv->node;
		int dk = pv->dk;

		// Factor diagonal block (packed storage)
		double Dk[LDL_MAX_BLOCK_DIM];
		if (block_ldl(c->diag_data[k], Dk, dk) < 0) { printf("[LDL] factorization failed at pivot step %d/%d (node %d, dim %d)\n", step, nc, k, dk); exit(-1); }
		for (int d = 0; d < dk; d++) c->diag_D[k][d] = Dk[d];

		// Back up edge matrices before overwriting with L blocks.
		// Schur updates use E_backup * L^T directly instead of L * N * L^T.
		double edge_backups[LDL_MAX_NODES][LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
		for (int ei = 0; ei < pv->col_count; ei++) {
			LDL_Column* en = &t->columns[pv->col_start + ei];
			memcpy(edge_backups[ei], &c->L_factors[en->L_offset], en->dn * dk * sizeof(double));
		}

		// Compute L blocks for non-eliminated neighbors
		for (int ei = 0; ei < pv->col_count; ei++) {
			LDL_Column* en = &t->columns[pv->col_start + ei];
			int di = en->dn;
			double* Eik = &c->L_factors[en->L_offset]; // currently holds K_{i,k}, will be overwritten with L_{i,k}

			// L_{i,k} = E_{i,k} * N_k^{-1}: solve N_k * X = E^T column by column
			double Lik[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
			for (int col = 0; col < di; col++) {
				double rhs[LDL_MAX_BLOCK_DIM], sol[LDL_MAX_BLOCK_DIM];
				for (int r = 0; r < dk; r++) rhs[r] = Eik[col * dk + r];
				block_solve(c->diag_data[k], Dk, rhs, sol, dk);
				for (int r = 0; r < dk; r++) Lik[col * dk + r] = sol[r];
			}
			memcpy(Eik, Lik, di * dk * sizeof(double));
		}

		// Apply Schur complement updates: product = E_{i,k}_backup * L_{j,k}^T
		for (int si = 0; si < pv->schur_count; si++) {
			LDL_Schur* op = &t->schurs[pv->schur_start + si];
			int di = op->di, dj = op->dj;

			// Find backup for E_{i,k}
			double* Eik_backup = NULL;
			for (int ei = 0; ei < pv->col_count; ei++) {
				if (t->columns[pv->col_start + ei].L_offset == op->Lik_offset) { Eik_backup = edge_backups[ei]; break; }
			}

			double* Ljk = &c->L_factors[op->Ljk_offset];
			double LjkT[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM], product[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
			block_transpose(Ljk, dj, dk, LjkT);
			block_mul(Eik_backup, di, dk, LjkT, dj, product);

			if (op->target_offset < 0) {
				// Diagonal Schur update — average (r,c) and (c,r) into packed storage.
				// E * L^T is symmetric in exact arithmetic; averaging absorbs float asymmetry.
				int tn = op->target_node;
				for (int r = 0; r < di; r++)
					for (int col = 0; col <= r; col++)
						c->diag_data[tn][LDL_TRI(r, col)] -= 0.5 * (product[r*di+col] + product[col*di+r]);
			} else {
				// Off-diagonal update: write to (i,j) and transpose to (j,i)
				block_sub(&c->L_factors[op->target_offset], product, di, dj);
				double product_T[LDL_MAX_BLOCK_DIM * LDL_MAX_BLOCK_DIM];
				block_transpose(product, di, dj, product_T);
				block_sub(&c->L_factors[op->target_offset_rev], product_T, dj, di);
			}
		}
	}

	// Post-factorization: verify all pivots are healthy
	ldl_condition_check(c, t);

	// Debug: capture D pivots after factorization
	if (dbg) {
		g_ldl_debug_info.n = t->n;
		for (int i = 0; i < nc; i++) {
			int di = t->dof[i], oi = t->row_offset[i];
			for (int d = 0; d < di; d++) {
				g_ldl_debug_info.D[oi + d] = c->diag_D[i][d];
			}
		}
	}
}

// Forward/diagonal/back substitution using precomputed topology.
static void ldl_solve_topo(LDL_Topology* t, double diag_data[][78], double diag_D[][12], double* L_factors, double* rhs, double* x)
{
	int nc = t->node_count;
	memcpy(x, rhs, t->n * sizeof(double));

	// Forward substitution
	for (int step = 0; step < nc; step++) {
		LDL_Pivot* pv = &t->pivots[step];
		int ok = pv->ok, dk = pv->dk;
		for (int fi = 0; fi < pv->fwd_count; fi++) {
			LDL_Neighbor* sn = &t->fwd_neighbors[pv->fwd_start + fi];
			double* Lkj = &L_factors[sn->L_offset]; // L_{k,j}: dk x dj block
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
		double tmp[LDL_MAX_BLOCK_DIM];
		for (int d = 0; d < dk; d++) tmp[d] = x[ok + d];
		block_solve(diag_data[k], diag_D[k], tmp, &x[ok], dk);
	}

	// Back substitution
	for (int step = nc - 1; step >= 0; step--) {
		LDL_Pivot* pv = &t->pivots[step];
		int ok = pv->ok, dk = pv->dk;
		for (int bi = 0; bi < pv->back_count; bi++) {
			LDL_Neighbor* sn = &t->back_neighbors[pv->back_start + bi];
			double* Ljk = &L_factors[sn->L_offset]; // L_{j,k}: dj x dk block. Apply L^T.
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
	afree(c->jacobians);
	afree(c->body_remap);
	afree(c->shard_counts);
	memset(c, 0, sizeof(*c));
}

// Free numeric data but keep topology and constraint list for sleep cache.
// When the island wakes, topology is reused if topo_version still matches.
static void ldl_cache_sleep(LDL_Cache* c)
{
	afree(c->L_factors);
	c->L_factors = NULL;
	afree(c->jacobians);
	c->jacobians = NULL;
}

// Build blocks for an island's joints.
static int ldl_cache_rebuild_blocks(LDL_Cache* c, WorldInternal* w, int island_idx, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count)
{
	afree(c->constraints);
	afree(c->bundles);
	afree(c->body_remap);
	afree(c->shard_counts);
	c->constraints = NULL;
	c->bundles = NULL;
	c->body_remap = NULL;
	c->shard_counts = NULL;
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
					LDL_Constraint con = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b, .real_body_a = sol_bs[i].body_a, .real_body_b = sol_bs[i].body_b, .weight_a = 1.0f, .weight_b = 1.0f, .solver_idx = i };
					assert(con.body_a != con.body_b && "Joint between body and itself");
					assert((w->body_hot[con.body_a].inv_mass > 0.0f || w->body_hot[con.body_b].inv_mass > 0.0f) && "Both bodies are static");
					apush(c->constraints, con);
					c->n += 3;
					c->joint_count++;
					break;
				}
			}
		} else if (j->type == JOINT_DISTANCE) {
			for (int i = 0; i < dist_count; i++) {
				if (sol_dist[i].joint_idx == ji) {
					LDL_Constraint con = { .type = JOINT_DISTANCE, .dof = 1, .body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b, .real_body_a = sol_dist[i].body_a, .real_body_b = sol_dist[i].body_b, .weight_a = 1.0f, .weight_b = 1.0f, .solver_idx = i };
					assert(con.body_a != con.body_b && "Joint between body and itself");
					assert((w->body_hot[con.body_a].inv_mass > 0.0f || w->body_hot[con.body_b].inv_mass > 0.0f) && "Both bodies are static");
					apush(c->constraints, con);
					c->n += 1;
					c->joint_count++;
					break;
				}
			}
		} else if (j->type == JOINT_HINGE) {
			for (int i = 0; i < hinge_count; i++) {
				if (sol_hinge[i].joint_idx == ji) {
					LDL_Constraint con = { .type = JOINT_HINGE, .dof = 5, .body_a = sol_hinge[i].body_a, .body_b = sol_hinge[i].body_b, .real_body_a = sol_hinge[i].body_a, .real_body_b = sol_hinge[i].body_b, .weight_a = 1.0f, .weight_b = 1.0f, .solver_idx = i };
					assert(con.body_a != con.body_b && "Joint between body and itself");
					assert((w->body_hot[con.body_a].inv_mass > 0.0f || w->body_hot[con.body_b].inv_mass > 0.0f) && "Both bodies are static");
					apush(c->constraints, con);
					c->n += 5;
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
	// Split bundles that would exceed LDL_MAX_BLOCK_DIM DOF to stay within packed storage limits.
	int start = 0;
	while (start < jc) {
		int ba = c->constraints[start].body_a, bb = c->constraints[start].body_b;
		int end = start + 1;
		while (end < jc && c->constraints[end].body_a == ba && c->constraints[end].body_b == bb) end++;
		// Emit one or more bundles, splitting if cumulative DOF > LDL_MAX_BLOCK_DIM
		int bstart = start;
		while (bstart < end) {
			int dof = 0;
			int bend = bstart;
			while (bend < end) {
				int next_dof = dof + c->constraints[bend].dof;
				if (next_dof > LDL_MAX_BLOCK_DIM && dof > 0) break; // split here (but always take at least one)
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

// Validate solve output: returns 0 if any NaN/inf/huge values found.
// Threshold 1e6: for a 1000 kg body at 60Hz/4 substeps, max reasonable
// lambda is ~mass*g*dt ≈ 40. Even extreme scenarios stay well under 1e6.
static int ldl_validate_lambda(double* lambda, int n)
{
	for (int i = 0; i < n; i++) {
		if (!(lambda[i] == lambda[i]) || lambda[i] > 1e6 || lambda[i] < -1e6) return 0;
	}
	return 1;
}

// Solve: build RHS, forward/back-sub via topology, apply velocity + position lambdas.
static void ldl_island_solve(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count, float sub_dt)
{
	int jc = c->joint_count;
	int n = c->n;
	LDL_Topology* t = c->topo;

	// --- Velocity solve ---
	// Generic RHS: vel_rhs[d] = -J_d * v - compliance * lambda_d for each DOF row.
	double* vel_rhs = CK_ALLOC(n * sizeof(double));
	double* vel_lambda = CK_ALLOC(n * sizeof(double));
	memset(vel_rhs, 0, n * sizeof(double));
	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
		BodyHot* a = &w->body_hot[con->real_body_a];
		BodyHot* b = &w->body_hot[con->real_body_b];
		LDL_JacobianRow* jac = &c->jacobians[con->jacobian_start];
		// Primary mode RHS: full constraint velocity error + compliance spring-back.
		double compliance = 0.0;
		if (!con->is_synthetic) {
			double soft = con->type == JOINT_BALL_SOCKET ? sol_bs[con->solver_idx].softness : con->type == JOINT_DISTANCE ? sol_dist[con->solver_idx].softness : sol_hinge[con->solver_idx].softness;
			compliance = soft > 0.0f ? soft : LDL_MIN_COMPLIANCE;
		}
		for (int d = 0; d < con->dof; d++) {
			double vel_err = ldl_constraint_velocity(&jac[d], a, b);
			double lam = 0.0;
			double bias = 0.0;
			if (!con->is_synthetic) {
				if (con->type == JOINT_BALL_SOCKET) {
					SolverBallSocket* s = &sol_bs[con->solver_idx];
					v3 l = s->lambda; lam = d == 0 ? l.x : d == 1 ? l.y : l.z;
					// For soft constraints, recompute bias from CURRENT positions each substep.
					// The pre-computed bias is stale after position integration.
					if (s->softness > 0.0f) {
						v3 ra = rotate(a->rotation, w->joints[s->joint_idx].ball_socket.local_a);
						v3 rb = rotate(b->rotation, w->joints[s->joint_idx].ball_socket.local_b);
						v3 err = sub(add(b->position, rb), add(a->position, ra));
						float ptv_f, soft_f;
						spring_compute(w->joints[s->joint_idx].ball_socket.spring, sub_dt, &ptv_f, &soft_f);
						v3 bb = scale(err, ptv_f);
						bias = d == 0 ? bb.x : d == 1 ? bb.y : bb.z;
					} else {
						v3 bb = s->bias; bias = d == 0 ? bb.x : d == 1 ? bb.y : bb.z;
					}
				} else if (con->type == JOINT_DISTANCE) {
					lam = sol_dist[con->solver_idx].lambda;
					bias = sol_dist[con->solver_idx].bias;
				} else {
					SolverHinge* sh = &sol_hinge[con->solver_idx];
					if (d < 3) { lam = d == 0 ? sh->lin_lambda.x : d == 1 ? sh->lin_lambda.y : sh->lin_lambda.z; bias = d == 0 ? sh->lin_bias.x : d == 1 ? sh->lin_bias.y : sh->lin_bias.z; }
					else { lam = sh->ang_lambda[d - 3]; bias = sh->ang_bias[d - 3]; }
				}
			}
			vel_rhs[oi + d] = -vel_err - bias - compliance * lam;
		}
	}

	ldl_solve_topo(t, c->diag_data, c->diag_D, c->L_factors, vel_rhs, vel_lambda);


	if (g_ldl_trace_solve && n <= 6) {
		printf("  [TRACE] n=%d jc=%d compliance=%.4g\n", n, jc, jc > 0 && !c->constraints[0].is_synthetic ? (c->constraints[0].type == JOINT_BALL_SOCKET ? (double)sol_bs[c->constraints[0].solver_idx].softness : 0.0) : 0.0);
		printf("  [TRACE] rhs=(%.4g, %.4g, %.4g)\n", vel_rhs[0], n>1?vel_rhs[1]:0, n>2?vel_rhs[2]:0);
		printf("  [TRACE] lam=(%.4g, %.4g, %.4g)\n", vel_lambda[0], n>1?vel_lambda[1]:0, n>2?vel_lambda[2]:0);
		printf("  [TRACE] D=(%.4g, %.4g, %.4g)\n", c->diag_D[0][0], c->diag_D[0][1], c->diag_D[0][2]);
		printf("  [TRACE] K_diag=(%.4g, %.4g, %.4g)\n", c->diag_data[0][LDL_TRI(0,0)], c->diag_data[0][LDL_TRI(1,1)], c->diag_data[0][LDL_TRI(2,2)]);
		LDL_Constraint* c0 = &c->constraints[0];
		BodyHot* ba = &w->body_hot[c0->real_body_a];
		BodyHot* bb = &w->body_hot[c0->real_body_b];
		printf("  [TRACE] vel_a=(%.4g,%.4g,%.4g) vel_b=(%.4g,%.4g,%.4g)\n", ba->velocity.x, ba->velocity.y, ba->velocity.z, bb->velocity.x, bb->velocity.y, bb->velocity.z);
		printf("  [TRACE] pos_b=(%.4g,%.4g,%.4g) rot_b=(%.4g,%.4g,%.4g,%.4g)\n", bb->position.x, bb->position.y, bb->position.z, bb->rotation.x, bb->rotation.y, bb->rotation.z, bb->rotation.w);
	}

	if (!ldl_validate_lambda(vel_lambda, n)) {
		CK_FREE(vel_rhs);
		CK_FREE(vel_lambda);
		return;
	}

	// Debug
	int dbg = g_ldl_debug_enabled && n <= LDL_MAX_DOF;
	if (dbg) { for (int i2 = 0; i2 < n; i2++) { g_ldl_debug_info.rhs[i2] = vel_rhs[i2]; g_ldl_debug_info.lambda_ldl[i2] = vel_lambda[i2]; }; g_ldl_debug_info.valid = 1; }

	// Apply velocity delta via generic Jacobian: v += M^{-1} * J^T * lambda.
	// Synthetic weld deltas are skipped -- they couple shards only.
	// Apply to REAL bodies (con->body_a/b may be virtual from shattering).
	// s->lambda accumulates PGS + LDL for next frame's warm-start.
	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		if (con->is_synthetic) continue;
		int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
		LDL_JacobianRow* jac = &c->jacobians[con->jacobian_start];
		// Use solver struct body indices (always real) for impulse application
		int real_a, real_b;
		if (con->type == JOINT_BALL_SOCKET) { real_a = sol_bs[con->solver_idx].body_a; real_b = sol_bs[con->solver_idx].body_b; }
		else if (con->type == JOINT_DISTANCE) { real_a = sol_dist[con->solver_idx].body_a; real_b = sol_dist[con->solver_idx].body_b; }
		else { real_a = sol_hinge[con->solver_idx].body_a; real_b = sol_hinge[con->solver_idx].body_b; }
		ldl_apply_jacobian_impulse(jac, con->dof, &vel_lambda[oi], &w->body_hot[real_a], 0);
		ldl_apply_jacobian_impulse(jac, con->dof, &vel_lambda[oi], &w->body_hot[real_b], 1);
		// LDL primary mode: vel_lambda is the full impulse, not a delta. SET lambda.
		if (con->type == JOINT_BALL_SOCKET) {
			sol_bs[con->solver_idx].lambda = V3((float)vel_lambda[oi], (float)vel_lambda[oi+1], (float)vel_lambda[oi+2]);
		} else if (con->type == JOINT_DISTANCE) {
			sol_dist[con->solver_idx].lambda = (float)vel_lambda[oi];
		} else {
			sol_hinge[con->solver_idx].lin_lambda = V3((float)vel_lambda[oi], (float)vel_lambda[oi+1], (float)vel_lambda[oi+2]);
			sol_hinge[con->solver_idx].ang_lambda[0] = (float)vel_lambda[oi+3];
			sol_hinge[con->solver_idx].ang_lambda[1] = (float)vel_lambda[oi+4];
		}
	}


	CK_FREE(vel_rhs);
	CK_FREE(vel_lambda);
}

// Position correction using factored K. Called AFTER integrate_positions to correct
// drift from velocity/rotation coupling. Reuses the K factored in ldl_numeric_factor.
static void ldl_island_position_correct(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, SolverDistance* sol_dist, SolverHinge* sol_hinge, float sub_dt)
{
	int jc = c->joint_count;
	int n = c->n;
	LDL_Topology* t = c->topo;

	double* pos_rhs = CK_ALLOC(n * sizeof(double));
	double* pos_lambda = CK_ALLOC(n * sizeof(double));
	memset(pos_rhs, 0, n * sizeof(double));
	double beta = SOLVER_POS_BAUMGARTE;
	double ptv = beta / sub_dt;

	for (int i = 0; i < jc; i++) {
		LDL_Constraint* con = &c->constraints[i];
		if (con->is_synthetic) continue;
		int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
		// Recompute lever arms from current rotation (positions changed since pre_solve)
		if (con->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[con->solver_idx];
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 ra = rotate(a->rotation, w->joints[s->joint_idx].ball_socket.local_a);
			v3 rb = rotate(b->rotation, w->joints[s->joint_idx].ball_socket.local_b);
			v3 err = sub(add(b->position, rb), add(a->position, ra));
			pos_rhs[oi]   = -ptv * err.x;
			pos_rhs[oi+1] = -ptv * err.y;
			pos_rhs[oi+2] = -ptv * err.z;
		} else if (con->type == JOINT_DISTANCE) {
			SolverDistance* s = &sol_dist[con->solver_idx];
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 ra = rotate(a->rotation, w->joints[s->joint_idx].distance.local_a);
			v3 rb = rotate(b->rotation, w->joints[s->joint_idx].distance.local_b);
			v3 d = sub(add(b->position, rb), add(a->position, ra));
			double dist_val = len(d);
			double err = dist_val - w->joints[s->joint_idx].distance.rest_length;
			pos_rhs[oi] = -ptv * err;
		} else if (con->type == JOINT_HINGE) {
			SolverHinge* s = &sol_hinge[con->solver_idx];
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 ra = rotate(a->rotation, w->joints[s->joint_idx].hinge.local_a);
			v3 rb = rotate(b->rotation, w->joints[s->joint_idx].hinge.local_b);
			v3 err = sub(add(b->position, rb), add(a->position, ra));
			pos_rhs[oi]   = -ptv * err.x;
			pos_rhs[oi+1] = -ptv * err.y;
			pos_rhs[oi+2] = -ptv * err.z;
			// Angular error: dot(t_d, axis_b) should be 0 when aligned
			v3 axis_a = norm(rotate(a->rotation, w->joints[s->joint_idx].hinge.local_axis_a));
			v3 axis_b = norm(rotate(b->rotation, w->joints[s->joint_idx].hinge.local_axis_b));
			v3 t1, t2;
			hinge_tangent_basis(axis_a, &t1, &t2);
			pos_rhs[oi+3] = -ptv * dot(t1, axis_b);
			pos_rhs[oi+4] = -ptv * dot(t2, axis_b);
		}
	}

	ldl_solve_topo(t, c->diag_data, c->diag_D, c->L_factors, pos_rhs, pos_lambda);

	// Apply position corrections simultaneously: accumulate per-body deltas, then apply once.
	// Sequential application would move shared bodies between constraint corrections,
	// breaking the coupled solution the LDL solve computed.
	if (ldl_validate_lambda(pos_lambda, n)) {
		int body_count = asize(w->body_hot);
		dv3* pos_delta = CK_ALLOC(body_count * sizeof(dv3));
		dv3* ang_delta = CK_ALLOC(body_count * sizeof(dv3));
		memset(pos_delta, 0, body_count * sizeof(dv3));
		memset(ang_delta, 0, body_count * sizeof(dv3));

		for (int i = 0; i < jc; i++) {
			LDL_Constraint* con = &c->constraints[i];
			if (con->is_synthetic) continue;
			int oi = t->row_offset[con->bundle_idx] + con->bundle_offset;
			LDL_JacobianRow* jac = &c->jacobians[con->jacobian_start];
			int real_a, real_b;
			if (con->type == JOINT_BALL_SOCKET) { real_a = sol_bs[con->solver_idx].body_a; real_b = sol_bs[con->solver_idx].body_b; }
			else if (con->type == JOINT_DISTANCE) { real_a = sol_dist[con->solver_idx].body_a; real_b = sol_dist[con->solver_idx].body_b; }
			else { real_a = sol_hinge[con->solver_idx].body_a; real_b = sol_hinge[con->solver_idx].body_b; }
			double inv_ma = (double)w->body_hot[real_a].inv_mass;
			double inv_mb = (double)w->body_hot[real_b].inv_mass;
			for (int d = 0; d < con->dof; d++) {
				double lam = pos_lambda[oi + d] * (double)sub_dt;
				pos_delta[real_a] = dv3_add(pos_delta[real_a], DV3(inv_ma * jac[d].J_a[0] * lam, inv_ma * jac[d].J_a[1] * lam, inv_ma * jac[d].J_a[2] * lam));
				dv3 j_ang_a = DV3(jac[d].J_a[3] * lam, jac[d].J_a[4] * lam, jac[d].J_a[5] * lam);
				dv3 dwa = dinv_inertia_mul(w->body_hot[real_a].rotation, w->body_hot[real_a].inv_inertia_local, j_ang_a);
				ang_delta[real_a] = dv3_add(ang_delta[real_a], dwa);
				pos_delta[real_b] = dv3_add(pos_delta[real_b], DV3(inv_mb * jac[d].J_b[0] * lam, inv_mb * jac[d].J_b[1] * lam, inv_mb * jac[d].J_b[2] * lam));
				dv3 j_ang_b = DV3(jac[d].J_b[3] * lam, jac[d].J_b[4] * lam, jac[d].J_b[5] * lam);
				dv3 dwb = dinv_inertia_mul(w->body_hot[real_b].rotation, w->body_hot[real_b].inv_inertia_local, j_ang_b);
				ang_delta[real_b] = dv3_add(ang_delta[real_b], dwb);
			}
		}

		for (int bi = 0; bi < body_count; bi++) {
			if (w->body_hot[bi].inv_mass == 0.0f) continue;
			w->body_hot[bi].position = add(w->body_hot[bi].position, dv3_to_v3(pos_delta[bi]));
			apply_rotation_delta(&w->body_hot[bi].rotation, ang_delta[bi]);
		}

		CK_FREE(pos_delta);
		CK_FREE(ang_delta);
	}

	CK_FREE(pos_rhs);
	CK_FREE(pos_lambda);
}

// -----------------------------------------------------------------------------
// Top-level entry point.

// Refresh solver-struct lever arms from current body rotations.
// Must be called before LDL factorization on substeps > 0, since integrate_positions
// rotates bodies but joints_pre_solve only computes lever arms once per frame.
static void ldl_refresh_lever_arms(WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count)
{
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &sol_bs[i];
		s->r_a = rotate(w->body_hot[s->body_a].rotation, w->joints[s->joint_idx].ball_socket.local_a);
		s->r_b = rotate(w->body_hot[s->body_b].rotation, w->joints[s->joint_idx].ball_socket.local_b);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &sol_dist[i];
		s->r_a = rotate(w->body_hot[s->body_a].rotation, w->joints[s->joint_idx].distance.local_a);
		s->r_b = rotate(w->body_hot[s->body_b].rotation, w->joints[s->joint_idx].distance.local_b);
		v3 anchor_a = add(w->body_hot[s->body_a].position, s->r_a);
		v3 anchor_b = add(w->body_hot[s->body_b].position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		float d = len(delta);
		s->axis = d > 1e-6f ? scale(delta, 1.0f / d) : V3(1, 0, 0);
	}
	for (int i = 0; i < hinge_count; i++) {
		SolverHinge* s = &sol_hinge[i];
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		s->r_a = rotate(a->rotation, w->joints[s->joint_idx].hinge.local_a);
		s->r_b = rotate(b->rotation, w->joints[s->joint_idx].hinge.local_b);
		v3 axis_a = norm(rotate(a->rotation, w->joints[s->joint_idx].hinge.local_axis_a));
		s->axis_b = norm(rotate(b->rotation, w->joints[s->joint_idx].hinge.local_axis_b));
		hinge_tangent_basis(axis_a, &s->t1, &s->t2);
		s->u1 = cross(s->t1, s->axis_b);
		s->u2 = cross(s->t2, s->axis_b);
	}
}

// Factor K for all islands (topology rebuild + numeric factorization). No solve.
static void ldl_factor(WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count, int sub)
{
	if (bs_count == 0 && dist_count == 0 && hinge_count == 0) return;

	// Refresh lever arms on substeps after the first (positions/rotations changed)
	if (sub > 0)
		ldl_refresh_lever_arms(w, sol_bs, bs_count, sol_dist, dist_count, sol_hinge, hinge_count);

	int island_count = asize(w->islands);
	for (int ii = 0; ii < island_count; ii++) {
		if (!(w->island_gen[ii] & 1)) continue;
		Island* isl = &w->islands[ii];
		if (!isl->awake) continue;
		if (isl->joint_count == 0) continue;

		LDL_Cache* c = &isl->ldl;

		// Rebuild blocks on first substep or topology change
		if (sub == 0 || c->topo_version != w->ldl_topo_version) {
			ldl_cache_rebuild_blocks(c, w, ii, sol_bs, bs_count, sol_dist, dist_count, sol_hinge, hinge_count);
			ldl_apply_shattering(c, w);
			ldl_build_bundles(c);
			// Only rebuild topology on actual graph change (skip if only solver indices changed)
			if (c->topo_version != w->ldl_topo_version)
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
		ldl_numeric_factor(c, w, sol_bs, sol_dist, sol_hinge);
	}
}

// Velocity correction using already-factored K. Solves and applies impulses.
static void ldl_velocity_correct(WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count, float sub_dt)
{
	if (bs_count == 0 && dist_count == 0 && hinge_count == 0) return;

	int island_count = asize(w->islands);
	for (int ii = 0; ii < island_count; ii++) {
		if (!(w->island_gen[ii] & 1)) continue;
		Island* isl = &w->islands[ii];
		if (!isl->awake) continue;
		if (isl->joint_count == 0) continue;
		LDL_Cache* c = &isl->ldl;
		if (c->n == 0 || !c->topo) continue;

		ldl_island_solve(c, w, sol_bs, bs_count, sol_dist, dist_count, sol_hinge, hinge_count, sub_dt);
	}
}

// Mass ratio compression exponent for position correction.
// 0.5 = sqrt compression (100:1 → 10:1). 1.0 = no compression (real masses).
#define LDL_POS_MASS_EXP 0.5f

// Position correction pass. Refactorizes K with compressed masses for better
// conditioning at extreme mass ratios. Position correction is stabilization,
// not dynamics — using compressed masses is safe.
static void ldl_position_correct(WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count, SolverHinge* sol_hinge, int hinge_count, float sub_dt)
{
	// Compress inv_mass and inv_inertia for position stage
	int body_count = asize(w->body_hot);
	float* save_inv_mass = CK_ALLOC(body_count * sizeof(float));
	v3* save_inv_inertia = CK_ALLOC(body_count * sizeof(v3));
	for (int i = 0; i < body_count; i++) {
		save_inv_mass[i] = w->body_hot[i].inv_mass;
		save_inv_inertia[i] = w->body_hot[i].inv_inertia_local;
		if (w->body_hot[i].inv_mass > 0.0f) {
			w->body_hot[i].inv_mass = powf(w->body_hot[i].inv_mass, LDL_POS_MASS_EXP);
			w->body_hot[i].inv_inertia_local = scale(w->body_hot[i].inv_inertia_local, w->body_hot[i].inv_mass / save_inv_mass[i]);
		}
	}

	// Refresh lever arms from current rotations (positions changed since pre-solve).
	// Without this, K uses stale lever arms while the RHS uses current positions.
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &sol_bs[i];
		s->r_a = rotate(w->body_hot[s->body_a].rotation, w->joints[s->joint_idx].ball_socket.local_a);
		s->r_b = rotate(w->body_hot[s->body_b].rotation, w->joints[s->joint_idx].ball_socket.local_b);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &sol_dist[i];
		s->r_a = rotate(w->body_hot[s->body_a].rotation, w->joints[s->joint_idx].distance.local_a);
		s->r_b = rotate(w->body_hot[s->body_b].rotation, w->joints[s->joint_idx].distance.local_b);
		v3 anchor_a = add(w->body_hot[s->body_a].position, s->r_a);
		v3 anchor_b = add(w->body_hot[s->body_b].position, s->r_b);
		v3 delta = sub(anchor_b, anchor_a);
		float d = len(delta);
		s->axis = d > 1e-6f ? scale(delta, 1.0f / d) : V3(1, 0, 0);
	}
	for (int i = 0; i < hinge_count; i++) {
		SolverHinge* s = &sol_hinge[i];
		BodyHot* a = &w->body_hot[s->body_a];
		BodyHot* b = &w->body_hot[s->body_b];
		s->r_a = rotate(a->rotation, w->joints[s->joint_idx].hinge.local_a);
		s->r_b = rotate(b->rotation, w->joints[s->joint_idx].hinge.local_b);
		v3 axis_a = norm(rotate(a->rotation, w->joints[s->joint_idx].hinge.local_axis_a));
		s->axis_b = norm(rotate(b->rotation, w->joints[s->joint_idx].hinge.local_axis_b));
		hinge_tangent_basis(axis_a, &s->t1, &s->t2);
		s->u1 = cross(s->t1, s->axis_b);
		s->u2 = cross(s->t2, s->axis_b);
	}

	int island_count = asize(w->islands);
	for (int ii = 0; ii < island_count; ii++) {
		if (!(w->island_gen[ii] & 1)) continue;
		Island* isl = &w->islands[ii];
		if (!isl->awake) continue;
		if (isl->joint_count == 0) continue;
		LDL_Cache* c = &isl->ldl;
		if (c->n == 0 || !c->topo) continue;
		// Refactorize K with compressed masses and current lever arms
		ldl_numeric_factor(c, w, sol_bs, sol_dist, sol_hinge);
		ldl_island_position_correct(c, w, sol_bs, sol_dist, sol_hinge, sub_dt);
	}

	// Restore real masses
	for (int i = 0; i < body_count; i++) {
		w->body_hot[i].inv_mass = save_inv_mass[i];
		w->body_hot[i].inv_inertia_local = save_inv_inertia[i];
	}
	CK_FREE(save_inv_mass);
	CK_FREE(save_inv_inertia);
}
