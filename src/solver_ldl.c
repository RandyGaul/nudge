// solver_ldl.c -- Direct LDL joint correction for equality constraints.
// Reference: Mizerski, "Improving an Iterative Physics Solver Using a Direct Method", GDC 2020.
//
// Sparse block LDL with body shattering:
//   Layer 1: Sparse block storage (adjacency-list per constraint node)
//   Layer 2: Sparse block LDL factorize/solve (min-degree ordering, fill-in)
//   Layer 3: Body shattering (virtual splinters for hub bodies)

static LDL_DebugInfo g_ldl_debug_info;
int g_ldl_debug_enabled;

#define SHATTER_THRESHOLD 9999 // shattering disabled pending star topology NaN debug
#define SPLINTER_TARGET   6  // target DOF per splinter

// =============================================================================
// Small block math helpers (3x3 and 1x1 blocks stored as flat floats).
// =============================================================================

// In-place 3x3 LDL^T factorize. L in lower triangle, D on diagonal.
static void block_ldl_3x3(float* A, float D[3])
{
	D[0] = A[0]; if (fabsf(D[0]) < 1e-12f) D[0] = 1e-12f;
	A[3] = A[3] / D[0]; A[6] = A[6] / D[0]; // L[1,0], L[2,0]
	D[1] = A[4] - A[3]*A[3]*D[0]; if (fabsf(D[1]) < 1e-12f) D[1] = 1e-12f;
	A[7] = (A[7] - A[6]*A[3]*D[0]) / D[1]; // L[2,1]
	D[2] = A[8] - A[6]*A[6]*D[0] - A[7]*A[7]*D[1]; if (fabsf(D[2]) < 1e-12f) D[2] = 1e-12f;
}

// Solve 3x3 LDL system: L stored in lower triangle of A (unit diagonal), D separate.
static void block_solve_3x3(float* L, float D[3], float* b, float* x)
{
	// Forward: Ly = b (L has unit diagonal, L[1,0]=L[3], L[2,0]=L[6], L[2,1]=L[7])
	x[0] = b[0];
	x[1] = b[1] - L[3]*x[0];
	x[2] = b[2] - L[6]*x[0] - L[7]*x[1];
	// Diagonal
	x[0] /= D[0]; x[1] /= D[1]; x[2] /= D[2];
	// Back: L^T x = z
	x[1] -= L[7]*x[2];
	x[0] -= L[3]*x[1] + L[6]*x[2];
}

// Compute C = A * B where A is (ra x ca), B is (ca x cb). Out is (ra x cb).
static void block_mul(float* A, int ra, int ca, float* B, int cb, float* out)
{
	for (int i = 0; i < ra; i++)
		for (int j = 0; j < cb; j++) {
			float sum = 0;
			for (int k = 0; k < ca; k++) sum += A[i*ca + k] * B[k*cb + j];
			out[i*cb + j] = sum;
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
	for (int i = 0; i < r; i++)
		for (int j = 0; j < c; j++)
			out[j*r + i] = A[i*c + j];
}

// =============================================================================
// Sparse block matrix operations.
// =============================================================================

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
	for (int k = 0; k < cnt; k++)
		if (s->adj[i][k] == j) return k;
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

// =============================================================================
// Diagonal block builders (K matrix computation).
// =============================================================================

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

// =============================================================================
// Off-diagonal block computation (same math as before, writes to block buffer).
// =============================================================================

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

// =============================================================================
// Sparse LDL factorize and solve.
// =============================================================================

// Minimum degree elimination ordering.
static void ldl_sparse_min_degree_order(LDL_Sparse* s)
{
	int nc = s->node_count;
	int eliminated[LDL_MAX_NODES] = {0};
	int degree[LDL_MAX_NODES];

	// Initial degree = number of neighbors
	for (int i = 0; i < nc; i++) degree[i] = asize(s->adj[i]);

	for (int step = 0; step < nc; step++) {
		// Find minimum-degree node among non-eliminated
		int best = -1, best_deg = 0x7FFFFFFF;
		for (int i = 0; i < nc; i++) {
			if (eliminated[i]) continue;
			if (degree[i] < best_deg) { best_deg = degree[i]; best = i; }
		}
		s->elim_order[step] = best;
		eliminated[best] = 1;

		// "Eliminate" best: add fill-in edges between its remaining neighbors
		int cnt = asize(s->adj[best]);
		for (int ii = 0; ii < cnt; ii++) {
			int ni = s->adj[best][ii];
			if (eliminated[ni]) continue;
			for (int jj = ii + 1; jj < cnt; jj++) {
				int nj = s->adj[best][jj];
				if (eliminated[nj]) continue;
				// Add fill-in edge (ni, nj) if it doesn't exist
				if (ldl_sparse_find_edge(s, ni, nj) < 0) {
					ldl_sparse_get_or_create_edge(s, ni, nj);
					degree[ni]++;
					degree[nj]++;
				}
			}
			// Decrement degree for losing neighbor 'best'
			degree[ni]--;
		}
	}
}

static void ldl_sparse_factorize(LDL_Sparse* s)
{
	int nc = s->node_count;

	for (int step = 0; step < nc; step++) {
		int k = s->elim_order[step];
		int dk = s->dof[k];

		// Factor diagonal block D_k
		float Dk[3];
		if (dk == 3) {
			block_ldl_3x3(s->diag_data[k], Dk);
			s->diag_D[k][0] = Dk[0]; s->diag_D[k][1] = Dk[1]; s->diag_D[k][2] = Dk[2];
		} else {
			float d = s->diag_data[k][0];
			Dk[0] = fabsf(d) > 1e-12f ? d : 1e-12f;
			s->diag_data[k][0] = Dk[0];
			s->diag_D[k][0] = Dk[0];
		}

		// For each neighbor i of k: compute L_{i,k} = E_{i,k} * D_k^{-1}
		int cnt_k = asize(s->adj[k]);
		for (int ai = 0; ai < cnt_k; ai++) {
			int i = s->adj[k][ai];
			int di = s->dof[i];

			// Get E_{i,k} (stored in adj[i]'s entry for k)
			float* Eik = ldl_sparse_get_edge(s, i, k);
			if (!Eik) continue;

			// L_{i,k} = E_{i,k} * D_k^{-1}
			// D_k^{-1} for LDL: D_k is factored, so D_k^{-1} = L_kk^{-T} * diag(Dk)^{-1} * L_kk^{-1}
			// But for simplicity: solve D_k * X^T = E_{i,k}^T column by column.
			// Actually, E_{i,k} * (L_kk * diag(Dk) * L_kk^T)^{-1} is complex.
			// Simpler: L_{i,k} = E_{i,k} * inv(diag_block_k) where diag_block_k is the original unfactored diagonal.
			// But we already factored it in place...
			//
			// The standard block LDL approach: L_{i,k} = E_{i,k} * N_k^{-1} where N_k is the ORIGINAL diagonal block.
			// Since we factored N_k in place (now contains L_k and D_k), we solve N_k * L_{i,k}^T = E_{i,k}^T.
			//
			// For 3x3: solve the 3x3 system for each column of L_{i,k}^T (i.e., each row of L_{i,k}).
			// For 1x1: L_{i,k} = E_{i,k} / D_k.

			if (dk == 1) {
				float inv_dk = 1.0f / Dk[0];
				for (int j = 0; j < di; j++) Eik[j] *= inv_dk;
				// Eik now stores L_{i,k}
			} else {
				// Solve N_k * X = E_{i,k}^T column by column, then transpose back.
				// N_k is stored in diag_data[k] as factored LDL.
				float Lik[9]; // di x dk block
				for (int col = 0; col < di; col++) {
					// Extract row 'col' of E_{i,k} as a dk-vector
					float rhs[3];
					for (int r = 0; r < dk; r++) rhs[r] = Eik[col * dk + r];
					float sol[3];
					block_solve_3x3(s->diag_data[k], Dk, rhs, sol);
					for (int r = 0; r < dk; r++) Lik[col * dk + r] = sol[r];
				}
				memcpy(Eik, Lik, di * dk * sizeof(float));
			}
		}

		// Schur complement: for each pair (i, j) of remaining neighbors of k:
		// off(i,j) -= L_{i,k} * D_k * L_{j,k}^T
		for (int ai = 0; ai < cnt_k; ai++) {
			int i = s->adj[k][ai];
			int di = s->dof[i];
			float* Lik = ldl_sparse_get_edge(s, i, k);
			if (!Lik) continue;

			// L_{i,k} * D_k: multiply each column of L_{i,k} by corresponding D_k element
			float LikDk[9]; // di x dk
			memcpy(LikDk, Lik, di * dk * sizeof(float));
			if (dk == 1) {
				for (int r = 0; r < di; r++) LikDk[r] *= Dk[0];
			} else {
				for (int r = 0; r < di; r++)
					for (int c = 0; c < dk; c++)
						LikDk[r * dk + c] *= Dk[c];
			}

			for (int aj = ai; aj < cnt_k; aj++) {
				int j = s->adj[k][aj];
				int dj = s->dof[j];
				float* Ljk = ldl_sparse_get_edge(s, j, k);
				if (!Ljk) continue;

				// Compute LikDk * Ljk^T (di x dj)
				float LjkT[9]; // dk x dj
				block_transpose(Ljk, dj, dk, LjkT);
				float product[9]; // di x dj
				block_mul(LikDk, di, dk, LjkT, dj, product);

				if (i == j) {
					block_sub(s->diag_data[i], product, di, di);
				} else {
					float* Sij = ldl_sparse_get_or_create_edge(s, i, j);
					block_sub(Sij, product, di, dj);
				}
			}
		}
	}
}

// Forward substitution, diagonal solve, back substitution using sparse block factors.
static void ldl_sparse_solve(LDL_Sparse* s, float* rhs, float* x)
{
	int nc = s->node_count;
	memcpy(x, rhs, s->n * sizeof(float));

	for (int step = 0; step < nc; step++) {
		int k = s->elim_order[step];
		int dk = s->dof[k];
		int ok = s->row_offset[k];

		// x_k -= sum over earlier pivots j: L_{k,j} * x_j
		// L_{k,j} is stored in adj[k] for j that was eliminated before k.
		// (We stored L_{i,k} in adj[i]'s entry for k during factorize.)
		// So L_{k,j} is in adj[k]'s entry for j.
		int cnt = asize(s->adj[k]);
		for (int ai = 0; ai < cnt; ai++) {
			int j = s->adj[k][ai];
			// Check if j was eliminated before k
			int j_step = -1;
			for (int ss = 0; ss < step; ss++) if (s->elim_order[ss] == j) { j_step = ss; break; }
			if (j_step < 0) continue; // j not yet eliminated

			int dj = s->dof[j];
			int oj = s->row_offset[j];
			float* Lkj = ldl_sparse_get_edge(s, k, j);
			if (!Lkj) continue;

			// x_k -= L_{k,j} * x_j
			for (int r = 0; r < dk; r++)
				for (int c = 0; c < dj; c++)
					x[ok + r] -= Lkj[r * dj + c] * x[oj + c];
		}
	}

	// Diagonal solve: x_k = D_k^{-1} * x_k
	for (int step = 0; step < nc; step++) {
		int k = s->elim_order[step];
		int dk = s->dof[k];
		int ok = s->row_offset[k];

		if (dk == 1) {
			x[ok] /= s->diag_D[k][0];
		} else {
			float tmp_in[3] = { x[ok], x[ok+1], x[ok+2] };
			block_solve_3x3(s->diag_data[k], s->diag_D[k], tmp_in, &x[ok]);
		}
	}

	// Back substitution: reverse elimination order
	for (int step = nc - 1; step >= 0; step--) {
		int k = s->elim_order[step];
		int dk = s->dof[k];
		int ok = s->row_offset[k];

		int cnt = asize(s->adj[k]);
		for (int ai = 0; ai < cnt; ai++) {
			int j = s->adj[k][ai];
			// Check if j was eliminated after k
			int j_step = -1;
			for (int ss = step + 1; ss < nc; ss++) if (s->elim_order[ss] == j) { j_step = ss; break; }
			if (j_step < 0) continue;

			int dj = s->dof[j];
			int oj = s->row_offset[j];
			float* Ljk = ldl_sparse_get_edge(s, j, k);
			if (!Ljk) continue;

			// x_k -= L_{j,k}^T * x_j
			for (int r = 0; r < dk; r++)
				for (int c = 0; c < dj; c++)
					x[ok + r] -= Ljk[c * dk + r] * x[oj + c];
		}
	}
}

// =============================================================================
// Body shattering: split hub bodies into virtual splinters.
// =============================================================================

// Apply shattering to the block list. Creates virtual bodies and synthetic joints.
// Modifies c->blocks (adds synthetic blocks, redirects body refs to virtual indices).
// Virtual body indices start at virtual_base = body_count (real bodies).
static void ldl_apply_shattering(LDL_Cache* c, WorldInternal* w)
{
	int jc = c->joint_count;
	int body_count = asize(w->body_hot);

	// Count total DOF per body
	int body_dof[4096] = {0}; // plenty for any reasonable body count
	assert(body_count < 4096);
	for (int i = 0; i < jc; i++) {
		body_dof[c->blocks[i].body_a] += c->blocks[i].dof;
		body_dof[c->blocks[i].body_b] += c->blocks[i].dof;
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
	typedef struct ShatterInfo { int body; int splinter_count; int first_virtual; } ShatterInfo;
	CK_DYNA ShatterInfo* shatter_list = NULL;

	for (int b = 0; b < body_count; b++) {
		if (body_dof[b] <= SHATTER_THRESHOLD) continue;
		if (w->body_hot[b].inv_mass == 0.0f) continue; // don't shatter static bodies

		int S = (body_dof[b] + SPLINTER_TARGET - 1) / SPLINTER_TARGET;
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

		ShatterInfo si = { .body = b, .splinter_count = S, .first_virtual = first_virt };
		apush(shatter_list, si);
	}

	if (asize(shatter_list) == 0) {
		afree(shatter_list);
		return; // nothing to shatter
	}

	// Redirect real joint body references to splinters (round-robin)
	// Track how many joints assigned to each splinter per hub body
	int assign_count[4096] = {0}; // per real body: running assignment counter
	for (int i = 0; i < jc; i++) {
		LDL_Block* blk = &c->blocks[i];

		int remap_a = c->body_remap[blk->body_a];
		if (remap_a >= 0) {
			// Find the ShatterInfo for this body
			for (int si = 0; si < asize(shatter_list); si++) {
				if (shatter_list[si].body == blk->body_a) {
					int S = shatter_list[si].splinter_count;
					int splinter = assign_count[blk->body_a] % S;
					blk->body_a = shatter_list[si].first_virtual + splinter;
					assign_count[blk->body_a < virtual_base ? blk->body_a : shatter_list[si].body]++;
					break;
				}
			}
		}

		int remap_b = c->body_remap[blk->body_b];
		if (remap_b >= 0) {
			for (int si = 0; si < asize(shatter_list); si++) {
				if (shatter_list[si].body == blk->body_b) {
					int S = shatter_list[si].splinter_count;
					int splinter = assign_count[blk->body_b] % S;
					blk->body_b = shatter_list[si].first_virtual + splinter;
					assign_count[blk->body_b < virtual_base ? blk->body_b : shatter_list[si].body]++;
					break;
				}
			}
		}
	}

	// Add synthetic rigid ball_socket joints between consecutive splinters
	for (int si = 0; si < asize(shatter_list); si++) {
		ShatterInfo* info = &shatter_list[si];
		for (int s = 0; s < info->splinter_count - 1; s++) {
			LDL_Block synth = {
				.type = JOINT_BALL_SOCKET,
				.dof = 3,
				.body_a = info->first_virtual + s,
				.body_b = info->first_virtual + s + 1,
				.solver_idx = -1,
				.is_synthetic = 1,
			};
			apush(c->blocks, synth);
			c->joint_count++;
			c->n += 3;
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

// =============================================================================
// Cache management.
// =============================================================================

static void ldl_cache_free(LDL_Cache* c)
{
	afree(c->blocks);
	if (c->sparse) { ldl_sparse_free(c->sparse); CK_FREE(c->sparse); }
	afree(c->virtual_bodies);
	afree(c->body_remap);
	memset(c, 0, sizeof(*c));
}

// Build blocks for an island's joints.
static int ldl_cache_rebuild_blocks(LDL_Cache* c, WorldInternal* w, int island_idx, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count)
{
	afree(c->blocks);
	afree(c->virtual_bodies);
	afree(c->body_remap);
	c->blocks = NULL;
	c->virtual_bodies = NULL;
	c->body_remap = NULL;
	c->virtual_body_count = 0;
	c->n = 0;
	c->joint_count = 0;

	Island* isl = &w->islands[island_idx];
	int ji = isl->head_joint;
	while (ji >= 0) {
		JointInternal* j = &w->joints[ji];
		if (j->type == JOINT_BALL_SOCKET) {
			for (int i = 0; i < bs_count; i++) {
				if (sol_bs[i].joint_idx == ji) {
					LDL_Block blk = { .type = JOINT_BALL_SOCKET, .dof = 3, .body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b, .solver_idx = i };
					apush(c->blocks, blk);
					c->n += 3;
					c->joint_count++;
					break;
				}
			}
		} else if (j->type == JOINT_DISTANCE) {
			for (int i = 0; i < dist_count; i++) {
				if (sol_dist[i].joint_idx == ji) {
					LDL_Block blk = { .type = JOINT_DISTANCE, .dof = 1, .body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b, .solver_idx = i };
					apush(c->blocks, blk);
					c->n += 1;
					c->joint_count++;
					break;
				}
			}
		}
		ji = j->island_next;
	}

	// Reallocate sparse matrix
	if (c->sparse) { ldl_sparse_free(c->sparse); CK_FREE(c->sparse); }
	c->sparse = NULL;

	return c->joint_count;
}

// Build the sparse block matrix from current body state and factorize.
static void ldl_cache_factorize(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, SolverDistance* sol_dist)
{
	int jc = c->joint_count;
	if (jc == 0) return;

	// Allocate fresh sparse matrix
	if (c->sparse) ldl_sparse_free(c->sparse);
	else c->sparse = CK_ALLOC(sizeof(LDL_Sparse));
	ldl_sparse_init(c->sparse);

	LDL_Sparse* sp = c->sparse;
	sp->node_count = jc;
	sp->n = 0;
	for (int i = 0; i < jc; i++) {
		sp->dof[i] = c->blocks[i].dof;
		sp->row_offset[i] = sp->n;
		sp->n += sp->dof[i];
	}
	sp->row_offset[jc] = sp->n;
	c->n = sp->n;

	// Fill diagonal blocks
	for (int i = 0; i < jc; i++) {
		LDL_Block* blk = &c->blocks[i];
		if (blk->is_synthetic) {
			// Synthetic rigid ball_socket: K from virtual splinter bodies
			BodyHot* a = ldl_get_body(w, c, blk->body_a);
			BodyHot* b = ldl_get_body(w, c, blk->body_b);
			// Synthetic joints are at body center (zero lever arm), so r_a = r_b = 0
			float softness = 0.0f;
			ldl_ball_socket_K(a, b, V3(0,0,0), V3(0,0,0), softness, sp->diag_data[i]);
		} else if (blk->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[blk->solver_idx];
			BodyHot* a = ldl_get_body(w, c, blk->body_a);
			BodyHot* b = ldl_get_body(w, c, blk->body_b);
			ldl_ball_socket_K(a, b, s->r_a, s->r_b, s->softness, sp->diag_data[i]);
		} else {
			SolverDistance* s = &sol_dist[blk->solver_idx];
			// Distance eff_mass was computed with real body mass -- recompute with virtual if shattered
			BodyHot* a = ldl_get_body(w, c, blk->body_a);
			BodyHot* b = ldl_get_body(w, c, blk->body_b);
			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s->r_a, s->axis)), s->r_a), s->axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s->r_b, s->axis)), s->r_b), s->axis);
			k += s->softness;
			sp->diag_data[i][0] = k;
		}
	}

	// Off-diagonal blocks via body adjacency
	// Body indices may include virtual bodies, so extend the adjacency range
	int real_body_count = asize(w->body_hot);
	int total_body_count = real_body_count + c->virtual_body_count;
	CK_DYNA int** body_adj = CK_ALLOC(total_body_count * sizeof(int*));
	memset(body_adj, 0, total_body_count * sizeof(int*));
	for (int i = 0; i < jc; i++) {
		apush(body_adj[c->blocks[i].body_a], i);
		apush(body_adj[c->blocks[i].body_b], i);
	}

	for (int b = 0; b < total_body_count; b++) {
		int cnt = asize(body_adj[b]);
		if (cnt < 2) continue;
		BodyHot* B = ldl_get_body(w, c, b);
		for (int ii = 0; ii < cnt; ii++) {
			for (int jj = ii + 1; jj < cnt; jj++) {
				int bi = body_adj[b][ii], bj = body_adj[b][jj];
				LDL_Block* blk_i = &c->blocks[bi], *blk_j = &c->blocks[bj];
				float s_i = (b == blk_i->body_b) ? 1.0f : -1.0f;
				float s_j = (b == blk_j->body_b) ? 1.0f : -1.0f;
				float sign = s_i * s_j;

				// Get lever arms (synthetic joints have zero lever arm)
				v3 ri = V3(0,0,0), rj = V3(0,0,0);
				if (!blk_i->is_synthetic) {
					if (blk_i->type == JOINT_BALL_SOCKET) ri = (b == blk_i->body_b) ? sol_bs[blk_i->solver_idx].r_b : sol_bs[blk_i->solver_idx].r_a;
					else ri = (b == blk_i->body_b) ? sol_dist[blk_i->solver_idx].r_b : sol_dist[blk_i->solver_idx].r_a;
				}
				if (!blk_j->is_synthetic) {
					if (blk_j->type == JOINT_BALL_SOCKET) rj = (b == blk_j->body_b) ? sol_bs[blk_j->solver_idx].r_b : sol_bs[blk_j->solver_idx].r_a;
					else rj = (b == blk_j->body_b) ? sol_dist[blk_j->solver_idx].r_b : sol_dist[blk_j->solver_idx].r_a;
				}

				float* edge_ij = ldl_sparse_get_or_create_edge(sp, bi, bj);
				float block_buf[9];
				if (blk_i->type == JOINT_BALL_SOCKET && blk_j->type == JOINT_BALL_SOCKET) {
					ldl_compute_off_diag_bs_bs(B, ri, rj, sign, block_buf);
					for (int r = 0; r < 9; r++) edge_ij[r] += block_buf[r];
					// Symmetric: also update (bj, bi)
					float* edge_ji = ldl_sparse_get_edge(sp, bj, bi);
					float block_T[9];
					block_transpose(block_buf, 3, 3, block_T);
					for (int r = 0; r < 9; r++) edge_ji[r] += block_T[r];
				} else if (blk_i->type == JOINT_BALL_SOCKET && blk_j->type == JOINT_DISTANCE) {
					ldl_compute_off_diag_bs_dist(B, ri, rj, sol_dist[blk_j->solver_idx].axis, sign, block_buf);
					for (int r = 0; r < 3; r++) edge_ij[r] += block_buf[r];
					float* edge_ji = ldl_sparse_get_edge(sp, bj, bi);
					for (int r = 0; r < 3; r++) edge_ji[r] += block_buf[r]; // 1x3 = transpose of 3x1
				} else if (blk_i->type == JOINT_DISTANCE && blk_j->type == JOINT_BALL_SOCKET) {
					v3 axis_i = sol_dist[blk_i->solver_idx].axis;
					ldl_compute_off_diag_bs_dist(B, rj, ri, axis_i, sign, block_buf);
					// This gives 3x1, but edge_ij is 1x3
					float* edge_ji = ldl_sparse_get_edge(sp, bj, bi);
					for (int r = 0; r < 3; r++) { edge_ji[r] += block_buf[r]; edge_ij[r] += block_buf[r]; }
				} else {
					float val = ldl_compute_off_diag_dist_dist(B, ri, rj, sol_dist[blk_i->solver_idx].axis, sol_dist[blk_j->solver_idx].axis, sign);
					edge_ij[0] += val;
					float* edge_ji = ldl_sparse_get_edge(sp, bj, bi);
					edge_ji[0] += val;
				}
			}
		}
	}

	for (int b = 0; b < total_body_count; b++) afree(body_adj[b]);
	CK_FREE(body_adj);

	// Debug snapshot
	int dbg = g_ldl_debug_enabled && sp->n <= LDL_MAX_DOF && jc <= 64;
	if (dbg) {
		// Rebuild dense A for debug visualization
		int n = sp->n;
		memset(g_ldl_debug_info.A, 0, n * n * sizeof(float));
		for (int i = 0; i < jc; i++) {
			int di = sp->dof[i], oi = sp->row_offset[i];
			for (int r = 0; r < di; r++)
				for (int col = 0; col < di; col++)
					g_ldl_debug_info.A[(oi+r)*n + oi+col] = sp->diag_data[i][r*di + col];
			int cnt = asize(sp->adj[i]);
			for (int ai = 0; ai < cnt; ai++) {
				int j = sp->adj[i][ai];
				int dj = sp->dof[j], oj = sp->row_offset[j];
				float* edata = &sp->adj_data[i][ai * di * dj];
				for (int r = 0; r < di; r++)
					for (int col = 0; col < dj; col++)
						g_ldl_debug_info.A[(oi+r)*n + oj+col] = edata[r*dj + col];
			}
		}
	}

	// Elimination ordering + factorize
	ldl_sparse_min_degree_order(sp);
	ldl_sparse_factorize(sp);

	if (dbg) {
		g_ldl_debug_info.n = sp->n;
		g_ldl_debug_info.joint_count = jc;
		g_ldl_debug_info.bs_count = 0;
		g_ldl_debug_info.dist_count = 0;
		for (int i = 0; i < jc; i++) {
			g_ldl_debug_info.block_dofs[i] = sp->dof[i];
			g_ldl_debug_info.block_rows[i] = sp->row_offset[i];
			g_ldl_debug_info.block_types[i] = c->blocks[i].type;
			if (c->blocks[i].type == JOINT_BALL_SOCKET) g_ldl_debug_info.bs_count++;
			else g_ldl_debug_info.dist_count++;
		}
		// Extract D values for debug
		for (int i = 0; i < jc; i++) {
			int di = sp->dof[i], oi = sp->row_offset[i];
			for (int d = 0; d < di; d++)
				g_ldl_debug_info.D[oi + d] = sp->diag_D[i][d];
		}
	}
}

// Solve: build RHS, sparse forward/back-sub, apply lambdas.
static void ldl_island_solve(LDL_Cache* c, WorldInternal* w, SolverBallSocket* sol_bs, int bs_count, SolverDistance* sol_dist, int dist_count)
{
	int jc = c->joint_count;
	int n = c->n;
	LDL_Sparse* sp = c->sparse;

	// Refresh virtual body state from real hub bodies (they may have changed)
	if (c->virtual_body_count > 0) {
		int real_count = asize(w->body_hot);
		for (int b = 0; b < real_count; b++) {
			if (c->body_remap && c->body_remap[b] >= 0) {
				int first_v = c->body_remap[b] - real_count;
				// Count how many splinters this body has
				for (int v = first_v; v < c->virtual_body_count; v++) {
					// Check if this virtual body belongs to this real body
					// (They're contiguous starting from first_v, same position/rotation/velocity)
					BodyHot* vb = &c->virtual_bodies[v];
					float S = vb->inv_mass / (w->body_hot[b].inv_mass > 0 ? w->body_hot[b].inv_mass : 1.0f);
					if (S < 1.5f) break; // not a splinter of this body
					vb->position = w->body_hot[b].position;
					vb->rotation = w->body_hot[b].rotation;
					vb->velocity = w->body_hot[b].velocity;
					vb->angular_velocity = w->body_hot[b].angular_velocity;
				}
			}
		}
	}

	// Undo warm-start impulses (skip synthetic joints)
	for (int i = 0; i < jc; i++) {
		LDL_Block* blk = &c->blocks[i];
		if (blk->is_synthetic) continue;
		if (blk->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[blk->solver_idx];
			if (s->lambda.x != 0 || s->lambda.y != 0 || s->lambda.z != 0)
				apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, neg(s->lambda));
		} else if (blk->type == JOINT_DISTANCE) {
			SolverDistance* s = &sol_dist[blk->solver_idx];
			if (s->lambda != 0)
				apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, scale(s->axis, -s->lambda));
		}
	}

	// Build RHS
	float* rhs = CK_ALLOC(n * sizeof(float));
	float* lambda = CK_ALLOC(n * sizeof(float));
	memset(rhs, 0, n * sizeof(float));
	for (int i = 0; i < jc; i++) {
		LDL_Block* blk = &c->blocks[i];
		int oi = sp->row_offset[i];
		if (blk->is_synthetic) {
			// Synthetic rigid joint: RHS = -(dv) where dv = velocity difference between splinters
			// (should be zero if splinters are synced, but bias corrects if not)
			BodyHot* a = ldl_get_body(w, c, blk->body_a);
			BodyHot* b = ldl_get_body(w, c, blk->body_b);
			v3 dv = sub(b->velocity, a->velocity);
			rhs[oi] = -dv.x; rhs[oi+1] = -dv.y; rhs[oi+2] = -dv.z;
		} else if (blk->type == JOINT_BALL_SOCKET) {
			SolverBallSocket* s = &sol_bs[blk->solver_idx];
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			v3 r = neg(add(dv, s->bias));
			rhs[oi] = r.x; rhs[oi+1] = r.y; rhs[oi+2] = r.z;
		} else if (blk->type == JOINT_DISTANCE) {
			SolverDistance* s = &sol_dist[blk->solver_idx];
			BodyHot* a = &w->body_hot[s->body_a];
			BodyHot* b = &w->body_hot[s->body_b];
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
			rhs[oi] = -dot(dv, s->axis) + s->bias;
		}
	}

	// Solve: use dense fallback for now (sparse solve has numerical issues with star topologies)
	// TODO: debug sparse solve for complete constraint graphs
	{
		// Rebuild dense A from the pre-factorize sparse structure
		// (sp was already factorized, so we can't read A from it. Rebuild from scratch.)
		float* A_dense = CK_ALLOC(n * n * sizeof(float));
		memset(A_dense, 0, n * n * sizeof(float));

		// Diagonal
		for (int i = 0; i < jc; i++) {
			LDL_Block* blk = &c->blocks[i];
			int oi = sp->row_offset[i], di = sp->dof[i];
			float K[9];
			if (blk->is_synthetic) {
				BodyHot* ba = ldl_get_body(w, c, blk->body_a);
				BodyHot* bb = ldl_get_body(w, c, blk->body_b);
				ldl_ball_socket_K(ba, bb, V3(0,0,0), V3(0,0,0), 0.0f, K);
			} else if (blk->type == JOINT_BALL_SOCKET) {
				SolverBallSocket* s = &sol_bs[blk->solver_idx];
				BodyHot* ba = ldl_get_body(w, c, blk->body_a);
				BodyHot* bb = ldl_get_body(w, c, blk->body_b);
				ldl_ball_socket_K(ba, bb, s->r_a, s->r_b, s->softness, K);
			} else {
				SolverDistance* s = &sol_dist[blk->solver_idx];
				BodyHot* ba = ldl_get_body(w, c, blk->body_a);
				BodyHot* bb = ldl_get_body(w, c, blk->body_b);
				float inv_m = ba->inv_mass + bb->inv_mass;
				float k = inv_m + dot(cross(inv_inertia_mul(ba->rotation, ba->inv_inertia_local, cross(s->r_a, s->axis)), s->r_a), s->axis) + dot(cross(inv_inertia_mul(bb->rotation, bb->inv_inertia_local, cross(s->r_b, s->axis)), s->r_b), s->axis) + s->softness;
				K[0] = k;
			}
			for (int r = 0; r < di; r++)
				for (int cc = 0; cc < di; cc++)
					A_dense[(oi+r)*n + oi+cc] = K[r*di + cc];
		}

		// Off-diagonal via body adjacency
		int real_bc = asize(w->body_hot);
		int total_bc = real_bc + c->virtual_body_count;
		CK_DYNA int** ba2 = CK_ALLOC(total_bc * sizeof(int*));
		memset(ba2, 0, total_bc * sizeof(int*));
		for (int i = 0; i < jc; i++) { apush(ba2[c->blocks[i].body_a], i); apush(ba2[c->blocks[i].body_b], i); }
		for (int b = 0; b < total_bc; b++) {
			int cnt = asize(ba2[b]);
			if (cnt < 2) continue;
			BodyHot* B = ldl_get_body(w, c, b);
			for (int ii = 0; ii < cnt; ii++) for (int jj2 = ii + 1; jj2 < cnt; jj2++) {
				int bi = ba2[b][ii], bj = ba2[b][jj2];
				LDL_Block* bi2 = &c->blocks[bi], *bj2 = &c->blocks[bj];
				float si2 = (b == bi2->body_b) ? 1.0f : -1.0f;
				float sj2 = (b == bj2->body_b) ? 1.0f : -1.0f;
				float sign = si2 * sj2;
				v3 ri = V3(0,0,0), rj2 = V3(0,0,0);
				if (!bi2->is_synthetic) { if (bi2->type == JOINT_BALL_SOCKET) ri = (b == bi2->body_b) ? sol_bs[bi2->solver_idx].r_b : sol_bs[bi2->solver_idx].r_a; else ri = (b == bi2->body_b) ? sol_dist[bi2->solver_idx].r_b : sol_dist[bi2->solver_idx].r_a; }
				if (!bj2->is_synthetic) { if (bj2->type == JOINT_BALL_SOCKET) rj2 = (b == bj2->body_b) ? sol_bs[bj2->solver_idx].r_b : sol_bs[bj2->solver_idx].r_a; else rj2 = (b == bj2->body_b) ? sol_dist[bj2->solver_idx].r_b : sol_dist[bj2->solver_idx].r_a; }
				int di = sp->dof[bi], dj = sp->dof[bj];
				int oi = sp->row_offset[bi], oj = sp->row_offset[bj];
				if (di == 3 && dj == 3) {
					float blk[9];
					ldl_compute_off_diag_bs_bs(B, ri, rj2, sign, blk);
					for (int r = 0; r < 3; r++) for (int cc = 0; cc < 3; cc++) { A_dense[(oi+r)*n+oj+cc] += blk[r*3+cc]; A_dense[(oj+cc)*n+oi+r] += blk[r*3+cc]; }
				} else if (di == 3 && dj == 1) {
					float blk[3];
					ldl_compute_off_diag_bs_dist(B, ri, rj2, sol_dist[bj2->solver_idx].axis, sign, blk);
					for (int r = 0; r < 3; r++) { A_dense[(oi+r)*n+oj] += blk[r]; A_dense[oj*n+oi+r] += blk[r]; }
				} else if (di == 1 && dj == 3) {
					float blk[3];
					ldl_compute_off_diag_bs_dist(B, rj2, ri, sol_dist[bi2->solver_idx].axis, sign, blk);
					for (int r = 0; r < 3; r++) { A_dense[(oj+r)*n+oi] += blk[r]; A_dense[oi*n+oj+r] += blk[r]; }
				} else {
					float val = ldl_compute_off_diag_dist_dist(B, ri, rj2, sol_dist[bi2->solver_idx].axis, sol_dist[bj2->solver_idx].axis, sign);
					A_dense[oi*n+oj] += val; A_dense[oj*n+oi] += val;
				}
			}
		}
		for (int b = 0; b < total_bc; b++) afree(ba2[b]);
		CK_FREE(ba2);

		// Dense LDL solve
		float* D_dense = CK_ALLOC(n * sizeof(float));
		// In-place LDL on A_dense
		for (int j = 0; j < n; j++) {
			float dj = A_dense[j*n+j];
			for (int k = 0; k < j; k++) dj -= A_dense[j*n+k]*A_dense[j*n+k]*D_dense[k];
			D_dense[j] = fabsf(dj) > 1e-12f ? dj : 1e-12f;
			float inv_dj = 1.0f / D_dense[j];
			for (int i = j+1; i < n; i++) {
				float lij = A_dense[i*n+j];
				for (int k = 0; k < j; k++) lij -= A_dense[i*n+k]*A_dense[j*n+k]*D_dense[k];
				A_dense[i*n+j] = lij * inv_dj;
			}
		}
		// Solve
		for (int i = 0; i < n; i++) { float s = rhs[i]; for (int k = 0; k < i; k++) s -= A_dense[i*n+k]*lambda[k]; lambda[i] = s; }
		for (int i = 0; i < n; i++) lambda[i] /= D_dense[i];
		for (int i = n-1; i >= 0; i--) { float s = lambda[i]; for (int k = i+1; k < n; k++) s -= A_dense[k*n+i]*lambda[k]; lambda[i] = s; }

		CK_FREE(A_dense);
		CK_FREE(D_dense);
	}


	// Debug
	int dbg = g_ldl_debug_enabled && n <= LDL_MAX_DOF;
	if (dbg) { memcpy(g_ldl_debug_info.lambda_ldl, lambda, n * sizeof(float)); g_ldl_debug_info.valid = 1; }

	// Apply exact lambdas
	if (c->virtual_body_count > 0) {
		// With shattering: apply ALL impulses (real + synthetic) to virtual/real bodies,
		// using the block's body references (which point to splinters for shattered bodies).
		int real_count = asize(w->body_hot);
		for (int i = 0; i < jc; i++) {
			LDL_Block* blk = &c->blocks[i];
			int oi = sp->row_offset[i];
			BodyHot* a = ldl_get_body(w, c, blk->body_a);
			BodyHot* b = ldl_get_body(w, c, blk->body_b);
			if (blk->is_synthetic) {
				// Synthetic rigid joint: apply to virtual splinters
				v3 impulse = V3(lambda[oi], lambda[oi+1], lambda[oi+2]);
				apply_impulse(a, b, V3(0,0,0), V3(0,0,0), impulse);
			} else if (blk->type == JOINT_BALL_SOCKET) {
				SolverBallSocket* s = &sol_bs[blk->solver_idx];
				v3 new_lambda = V3(lambda[oi], lambda[oi+1], lambda[oi+2]);
				s->lambda = new_lambda;
				apply_impulse(a, b, s->r_a, s->r_b, new_lambda);
			} else if (blk->type == JOINT_DISTANCE) {
				SolverDistance* s = &sol_dist[blk->solver_idx];
				s->lambda = lambda[oi];
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->axis, s->lambda));
			}
		}
		// Copy splinter velocities back to real hub bodies
		for (int b = 0; b < real_count; b++) {
			if (!c->body_remap || c->body_remap[b] < 0) continue;
			int first_v = c->body_remap[b] - real_count;
			w->body_hot[b].velocity = c->virtual_bodies[first_v].velocity;
			w->body_hot[b].angular_velocity = c->virtual_bodies[first_v].angular_velocity;
		}
	} else {
		// No shattering: apply real joint lambdas directly
		for (int i = 0; i < jc; i++) {
			LDL_Block* blk = &c->blocks[i];
			int oi = sp->row_offset[i];
			if (blk->type == JOINT_BALL_SOCKET) {
				SolverBallSocket* s = &sol_bs[blk->solver_idx];
				v3 new_lambda = V3(lambda[oi], lambda[oi+1], lambda[oi+2]);
				s->lambda = new_lambda;
				apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, new_lambda);
			} else if (blk->type == JOINT_DISTANCE) {
				SolverDistance* s = &sol_dist[blk->solver_idx];
				s->lambda = lambda[oi];
				apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b], s->r_a, s->r_b, scale(s->axis, s->lambda));
			}
		}
	}

	CK_FREE(rhs);
	CK_FREE(lambda);
}

// =============================================================================
// Top-level entry point.
// =============================================================================

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
			c->topo_version = w->ldl_topo_version;
		}
		if (c->n == 0) continue;

		// Factorize every substep (lever arms change after integrate_positions)
		ldl_cache_factorize(c, w, sol_bs, sol_dist);

		// Solve every substep
		ldl_island_solve(c, w, sol_bs, bs_count, sol_dist, dist_count);
	}
}
