// solver_cr.c -- Conjugate Residual
//
// Analysis: Global Krylov solver for contact velocity constraints. Builds the full island
// Delassus system A = J M^-1 J^T + Sigma and solves via conjugate residual with inequality
// projection. When LDL is active, CR handles contacts only (LDL direct-solves joints).
// When LDL is off, CR handles both contacts and joints.
//
// Evolution:
//   v1: whole-island CR with warm-start and post-clamp
//   v2: PGS warmup (3 sweeps) to discover active set before CR
//   v3: active-set masking (skip inactive contacts from PGS warmup)
//   v4: mid-solve reclamp every 5 iters with CR direction restart
//
// Failed experiments:
//   PGS inside CR: tried as SSOR preconditioning (every iter) and at reclamp points.
//   Both worse than plain clamping. PGS every iteration destroys conjugate directions,
//   reducing CR to steepest descent. PGS at reclamp created larger residual bumps than
//   simple constraint-space clamping.
//
//   Block-diagonal preconditioning: better spectral approximation than Jacobi, but
//   P^-1 A lost positive spectrum on high mass ratios, breaking CR's inner product.
//   Would need PCG to fix. Pivoted to symmetric scaling instead.
//
// Symmetric diagonal scaling (current approach):
//   S = diag(1/sqrt(A_dd)) transforms system to unit diagonal: A_tilde = SAS. Algebraically
//   identical to Jacobi preconditioning but preserves SPD (no sign breakdown). Cost: two
//   O(n) vector scales per mat-vec. Result: 2.3-3.7x convergence improvement on pyramids
//   (55-204 boxes), neutral on small systems.

#define CR_MIN_COMPLIANCE 5e-5f

static int g_cr_trace = 0; // set to 1 for per-island diagnostic output

// Per-DOF clamp types for post-solve projection.
enum
{
	CR_CLAMP_NONE,       // bilateral (joints, no clamp)
	CR_CLAMP_LO_ZERO,    // contact normal: lambda >= 0
	CR_CLAMP_FRICTION,   // friction tangent: first of a consecutive pair (t1, t2)
	CR_CLAMP_FRICTION2,  // friction tangent: second of pair (skipped, handled with FRICTION)
	CR_CLAMP_TWIST,      // torsional friction (PATCH mode): |lambda| <= mu * sum(normals) * radius
};

// Flat per-DOF descriptor for the CR system.
typedef struct CR_System
{
	int n;                   // total scalar DOFs
	JacobianRow* rows;       // [n] Jacobian per DOF (owned, built for contacts, copied for joints)
	int* body_a;             // [n] body A index per DOF
	int* body_b;             // [n] body B index per DOF
	float* softness;         // [n] compliance per DOF
	float* bias;             // [n] RHS bias term (includes bounce for contacts)
	float* lambda;           // [n] working solution vector (warm-started)
	float* lambda_warm;      // [n] snapshot of warm values before CR iteration
	uint8_t* clamp_type;     // [n] CR_CLAMP_*
	int* friction_parent;    // [n] for friction/twist: index of parent normal DOF
	int* normal_count;       // [n] for PATCH friction: how many consecutive normals to sum
	float* friction_coeff;   // [n] mu for friction/twist DOFs
	float* friction_radius;  // [n] patch radius for twist DOFs
	float** lambda_ptrs;     // [n] pointer to original lambda storage for writeback
} CR_System;

// Build a JacobianRow from a direction vector and contact offsets.
// Sign convention matches jac_apply: v_a += inv_mass_a * J_a * lambda,
// so J_a linear = -dir (body A pushed away from contact).
static JacobianRow cr_build_contact_jacobian(v3 dir, v3 r_a, v3 r_b)
{
	v3 cr_a = cross(r_a, dir);
	v3 cr_b = cross(r_b, dir);
	JacobianRow row;
	row.J_a[0] = -dir.x; row.J_a[1] = -dir.y; row.J_a[2] = -dir.z;
	row.J_a[3] = -cr_a.x; row.J_a[4] = -cr_a.y; row.J_a[5] = -cr_a.z;
	row.J_b[0] =  dir.x; row.J_b[1] =  dir.y; row.J_b[2] =  dir.z;
	row.J_b[3] =  cr_b.x; row.J_b[4] =  cr_b.y; row.J_b[5] =  cr_b.z;
	row.eff_mass = 0; // not used by CR (CR computes its own)
	return row;
}

// Build a JacobianRow for torsional friction (pure angular along normal).
static JacobianRow cr_build_twist_jacobian(v3 n)
{
	JacobianRow row;
	row.J_a[0] = 0; row.J_a[1] = 0; row.J_a[2] = 0;
	row.J_a[3] = -n.x; row.J_a[4] = -n.y; row.J_a[5] = -n.z;
	row.J_b[0] = 0; row.J_b[1] = 0; row.J_b[2] = 0;
	row.J_b[3] =  n.x; row.J_b[4] =  n.y; row.J_b[5] =  n.z;
	row.eff_mass = 0;
	return row;
}

// Add a DOF to the CR system.
static void cr_add_dof(CR_System* sys, JacobianRow jac, int ba, int bb, float soft, float b, float lam, float* lam_ptr, uint8_t clamp, int fric_parent, int norm_count, float fric_coeff, float fric_radius)
{
	int i = sys->n++;
	sys->rows[i] = jac;
	sys->body_a[i] = ba;
	sys->body_b[i] = bb;
	sys->softness[i] = soft > 0.0f ? soft : CR_MIN_COMPLIANCE;
	sys->bias[i] = b;
	sys->lambda[i] = lam;
	sys->lambda_warm[i] = lam;
	sys->clamp_type[i] = clamp;
	sys->friction_parent[i] = fric_parent;
	sys->normal_count[i] = norm_count;
	sys->friction_coeff[i] = fric_coeff;
	sys->friction_radius[i] = fric_radius;
	sys->lambda_ptrs[i] = lam_ptr;
}

// Count total DOFs for an island to pre-allocate.
// contacts_only: skip joints (when LDL handles them separately).
static int cr_count_island_dofs(WorldInternal* w, int island_idx, SolverManifold* sm, int sm_count, SolverContact* sc, SolverJoint* sol_joints, int joint_count, int contacts_only)
{
	Island* isl = &w->islands[island_idx];
	int n = 0;

	// Joints
	if (!contacts_only) {
	int ji = isl->head_joint;
	while (ji >= 0) {
		for (int i = 0; i < joint_count; i++) {
			if (sol_joints[i].joint_idx == ji) {
				n += sol_joints[i].dof;
				break;
			}
		}
		ji = w->joints[ji].island_next;
	}
	}

	// Contacts
	int patch = (w->friction_model == FRICTION_PATCH);
	for (int i = 0; i < sm_count; i++) {
		int isl_a = w->body_cold[sm[i].body_a].island_id;
		int isl_b = w->body_cold[sm[i].body_b].island_id;
		if (isl_a != island_idx && isl_b != island_idx) continue;
		if (patch) {
			n += sm[i].contact_count; // 1 normal per contact
			n += 3; // 2 tangent + 1 twist per manifold
		} else {
			n += sm[i].contact_count * 3; // normal + t1 + t2 per contact
		}
	}
	return n;
}

// Build the CR system for an island: collect joints and contacts into flat DOF arrays.
// contacts_only: skip joints (when LDL handles them separately).
static void cr_build_island_system(CR_System* sys, WorldInternal* w, int island_idx, SolverManifold* sm, int sm_count, SolverContact* sc, SolverJoint* sol_joints, int joint_count, int max_dofs, int contacts_only)
{
	// Allocate all arrays at max size
	sys->n = 0;
	sys->rows = CK_ALLOC(max_dofs * sizeof(JacobianRow));
	sys->body_a = CK_ALLOC(max_dofs * sizeof(int));
	sys->body_b = CK_ALLOC(max_dofs * sizeof(int));
	sys->softness = CK_ALLOC(max_dofs * sizeof(float));
	sys->bias = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda_warm = CK_ALLOC(max_dofs * sizeof(float));
	sys->clamp_type = CK_ALLOC(max_dofs * sizeof(uint8_t));
	sys->friction_parent = CK_ALLOC(max_dofs * sizeof(int));
	sys->normal_count = CK_ALLOC(max_dofs * sizeof(int));
	sys->friction_coeff = CK_ALLOC(max_dofs * sizeof(float));
	sys->friction_radius = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda_ptrs = CK_ALLOC(max_dofs * sizeof(float*));

	Island* isl = &w->islands[island_idx];
	int patch = (w->friction_model == FRICTION_PATCH);

	// --- Joints (skipped when LDL handles them) ---
	if (!contacts_only) {
	int ji = isl->head_joint;
	while (ji >= 0) {
		for (int i = 0; i < joint_count; i++) {
			if (sol_joints[i].joint_idx != ji) continue;
			SolverJoint* sj = &sol_joints[i];
			for (int d = 0; d < sj->dof; d++) {
				cr_add_dof(sys, sj->rows[d], sj->body_a, sj->body_b, sj->softness, sj->bias[d], sj->lambda[d], &sj->lambda[d], CR_CLAMP_NONE, -1, 0, 0, 0);
			}
			break;
		}
		ji = w->joints[ji].island_next;
	}
	}

	// --- Contacts ---
	for (int mi = 0; mi < sm_count; mi++) {
		SolverManifold* m = &sm[mi];
		int isl_a = w->body_cold[m->body_a].island_id;
		int isl_b = w->body_cold[m->body_b].island_id;
		if (isl_a != island_idx && isl_b != island_idx) continue;

		// Track where normals start for this manifold (for PATCH friction parent)
		int first_normal_idx = sys->n;
		int active_normals = 0;
		int mask = w->cr_active_set_mask;

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			// Active-set masking: PGS warmup already ran. If it left lambda_n = 0,
			// this contact is inactive (separating). Skip it — CR doesn't need to
			// solve for a DOF that's already at its bound.
			if (mask && s->lambda_n == 0.0f) continue;
			active_normals++;

			// Normal DOF
			JacobianRow jn = cr_build_contact_jacobian(s->normal, s->r_a, s->r_b);
			float bias_n = s->bias + s->bounce;
			cr_add_dof(sys, jn, m->body_a, m->body_b, s->softness, bias_n, s->lambda_n, &s->lambda_n, CR_CLAMP_LO_ZERO, -1, 0, 0, 0);

			if (!patch) {
				// COULOMB: per-contact tangent DOFs
				// Use same softness as parent normal for conditioning.
				int normal_idx = sys->n - 1;
				float fric_soft = s->softness;

				JacobianRow jt1 = cr_build_contact_jacobian(s->tangent1, s->r_a, s->r_b);
				cr_add_dof(sys, jt1, m->body_a, m->body_b, fric_soft, 0, s->lambda_t1, &s->lambda_t1, CR_CLAMP_FRICTION, normal_idx, 1, m->friction, 0);

				JacobianRow jt2 = cr_build_contact_jacobian(s->tangent2, s->r_a, s->r_b);
				cr_add_dof(sys, jt2, m->body_a, m->body_b, fric_soft, 0, s->lambda_t2, &s->lambda_t2, CR_CLAMP_FRICTION2, normal_idx, 1, m->friction, 0);
			}
		}

		if (patch) {
			// Skip manifold friction entirely if no normals are active (masking on)
			if (mask && active_normals == 0) continue;
			// PATCH: manifold-level tangent friction at centroid.
			// Use first active contact's softness for friction conditioning.
			int normal_count = active_normals;
			float patch_soft = (m->contact_count > 0) ? sc[m->contact_start].softness : 0;

			JacobianRow jt1 = cr_build_contact_jacobian(m->tangent1, m->centroid_r_a, m->centroid_r_b);
			cr_add_dof(sys, jt1, m->body_a, m->body_b, patch_soft, 0, m->lambda_t1, &m->lambda_t1, CR_CLAMP_FRICTION, first_normal_idx, normal_count, m->friction, 0);

			JacobianRow jt2 = cr_build_contact_jacobian(m->tangent2, m->centroid_r_a, m->centroid_r_b);
			cr_add_dof(sys, jt2, m->body_a, m->body_b, patch_soft, 0, m->lambda_t2, &m->lambda_t2, CR_CLAMP_FRICTION2, first_normal_idx, normal_count, m->friction, 0);

			// Torsional friction (pure angular)
			JacobianRow jtw = cr_build_twist_jacobian(m->normal);
			cr_add_dof(sys, jtw, m->body_a, m->body_b, patch_soft, 0, m->lambda_twist, &m->lambda_twist, CR_CLAMP_TWIST, first_normal_idx, normal_count, m->friction, m->patch_radius);
		}
	}
}

static void cr_free_system(CR_System* sys)
{
	CK_FREE(sys->rows);
	CK_FREE(sys->body_a);
	CK_FREE(sys->body_b);
	CK_FREE(sys->softness);
	CK_FREE(sys->bias);
	CK_FREE(sys->lambda);
	CK_FREE(sys->lambda_warm);
	CK_FREE(sys->clamp_type);
	CK_FREE(sys->friction_parent);
	CK_FREE(sys->normal_count);
	CK_FREE(sys->friction_coeff);
	CK_FREE(sys->friction_radius);
	CK_FREE(sys->lambda_ptrs);
	memset(sys, 0, sizeof(*sys));
}

// --- Vector ops ---

static float cr_dot(float* a, float* b, int n)
{
	float s = 0.0f;
	for (int i = 0; i < n; i++) s += a[i] * b[i];
	return s;
}

// --- Matrix-free mat-vec: out = (J M^{-1} J^T + Sigma) * v ---
//
// Uses temporary body-space accumulators (caller-provided for reuse).
// Does NOT touch body velocities -- purely operates on the constraint-space vector.
static void cr_matvec(CR_System* sys, WorldInternal* w, float* v, float* out, v3* f_lin, v3* f_ang, int body_count)
{
	int n = sys->n;

	// Zero body accumulators
	memset(f_lin, 0, body_count * sizeof(v3));
	memset(f_ang, 0, body_count * sizeof(v3));

	// Step 1: Scatter -- accumulate J^T * v into body-space forces
	for (int d = 0; d < n; d++) {
		float lam = v[d];
		if (lam == 0.0f) continue;
		JacobianRow* row = &sys->rows[d];
		int ba = sys->body_a[d], bb = sys->body_b[d];
		f_lin[ba].x += row->J_a[0] * lam;
		f_lin[ba].y += row->J_a[1] * lam;
		f_lin[ba].z += row->J_a[2] * lam;
		f_ang[ba].x += row->J_a[3] * lam;
		f_ang[ba].y += row->J_a[4] * lam;
		f_ang[ba].z += row->J_a[5] * lam;
		f_lin[bb].x += row->J_b[0] * lam;
		f_lin[bb].y += row->J_b[1] * lam;
		f_lin[bb].z += row->J_b[2] * lam;
		f_ang[bb].x += row->J_b[3] * lam;
		f_ang[bb].y += row->J_b[4] * lam;
		f_ang[bb].z += row->J_b[5] * lam;
	}

	// Step 2: Scale -- M^{-1} * f (in-place)
	for (int i = 0; i < body_count; i++) {
		BodyHot* bh = &w->body_hot[i];
		if (bh->inv_mass == 0.0f) {
			f_lin[i] = V3(0, 0, 0);
			f_ang[i] = V3(0, 0, 0);
			continue;
		}
		// Only process bodies that had force accumulated
		if (f_lin[i].x == 0 && f_lin[i].y == 0 && f_lin[i].z == 0 &&
			f_ang[i].x == 0 && f_ang[i].y == 0 && f_ang[i].z == 0)
			continue;
		f_lin[i] = scale(f_lin[i], bh->inv_mass);
		f_ang[i] = inv_inertia_mul(bh->rotation, bh->inv_inertia_local, f_ang[i]);
	}

	// Step 3: Gather -- out = J * a + softness * v
	for (int d = 0; d < n; d++) {
		JacobianRow* row = &sys->rows[d];
		int ba = sys->body_a[d], bb = sys->body_b[d];
		v3 al = f_lin[ba], aa = f_ang[ba];
		v3 bl = f_lin[bb], ba2 = f_ang[bb];
		float r = row->J_a[0]*al.x + row->J_a[1]*al.y + row->J_a[2]*al.z
		        + row->J_a[3]*aa.x + row->J_a[4]*aa.y + row->J_a[5]*aa.z
		        + row->J_b[0]*bl.x + row->J_b[1]*bl.y + row->J_b[2]*bl.z
		        + row->J_b[3]*ba2.x + row->J_b[4]*ba2.y + row->J_b[5]*ba2.z;
		out[d] = r + sys->softness[d] * v[d];
	}
}

// --- Post-solve clamp: project lambda onto feasible set ---
static void cr_post_clamp(CR_System* sys)
{
	int n = sys->n;
	for (int d = 0; d < n; d++) {
		switch (sys->clamp_type[d]) {
		case CR_CLAMP_NONE:
			break;
		case CR_CLAMP_LO_ZERO:
			if (sys->lambda[d] < 0.0f)
				sys->lambda[d] = 0.0f;
			break;
		case CR_CLAMP_FRICTION: {
			// Friction pair: d = tangent1, d+1 = tangent2 (FRICTION2)
			// Sum parent normal(s) for max friction force
			float normal_sum = 0.0f;
			int pi = sys->friction_parent[d];
			int nc = sys->normal_count[d];
			for (int k = 0; k < nc; k++)
				normal_sum += sys->lambda[pi + k];
			if (normal_sum < 0.0f) normal_sum = 0.0f;
			float max_f = sys->friction_coeff[d] * normal_sum;
			float t1 = sys->lambda[d];
			float t2 = sys->lambda[d + 1];
			float mag = sqrtf(t1 * t1 + t2 * t2);
			if (mag > max_f && mag > 1e-12f) {
				float s = max_f / mag;
				sys->lambda[d] = t1 * s;
				sys->lambda[d + 1] = t2 * s;
			}
			break;
		}
		case CR_CLAMP_FRICTION2:
			// Handled by the FRICTION case above (d-1)
			break;
		case CR_CLAMP_TWIST: {
			float normal_sum = 0.0f;
			int pi = sys->friction_parent[d];
			int nc = sys->normal_count[d];
			for (int k = 0; k < nc; k++)
				normal_sum += sys->lambda[pi + k];
			if (normal_sum < 0.0f) normal_sum = 0.0f;
			float max_tw = sys->friction_coeff[d] * normal_sum * sys->friction_radius[d];
			if (sys->lambda[d] > max_tw) sys->lambda[d] = max_tw;
			if (sys->lambda[d] < -max_tw) sys->lambda[d] = -max_tw;
			break;
		}
		}
	}
}

// --- CR solve for one island. Returns iteration count. ---
static int cr_island_solve(CR_System* sys, WorldInternal* w)
{
	int n = sys->n;
	if (n == 0) return 0;

	int body_count = asize(w->body_hot);
	int max_iters = w->cr_max_iters;
	float tol = w->cr_tolerance;

	// Allocate CR workspace + body-space accumulators
	float* r   = CK_ALLOC(n * sizeof(float));
	float* Ar  = CK_ALLOC(n * sizeof(float));
	float* p   = CK_ALLOC(n * sizeof(float));
	float* Ap  = CK_ALLOC(n * sizeof(float));
	v3* f_lin  = CK_ALLOC(body_count * sizeof(v3));
	v3* f_ang  = CK_ALLOC(body_count * sizeof(v3));

	// Symmetric mass scaling: S = diag(1/√A_dd). Transforms the system to
	// Ã = SAS with unit diagonal, removing the trivial mass-ratio component
	// from the condition number. All CR vectors live in scaled space.
	// s[d] = scaling factor, sinv[d] = 1/s[d] for unscaling.
	float* s    = NULL;
	float* sinv = NULL;
	int use_scale = w->cr_mass_scale;
	if (use_scale) {
		s    = CK_ALLOC(n * sizeof(float));
		sinv = CK_ALLOC(n * sizeof(float));
		for (int d = 0; d < n; d++) {
			JacobianRow* row = &sys->rows[d];
			BodyHot* ba = &w->body_hot[sys->body_a[d]];
			BodyHot* bb = &w->body_hot[sys->body_b[d]];
			float k = 0;
			k += (row->J_a[0]*row->J_a[0] + row->J_a[1]*row->J_a[1] + row->J_a[2]*row->J_a[2]) * ba->inv_mass;
			k += (row->J_b[0]*row->J_b[0] + row->J_b[1]*row->J_b[1] + row->J_b[2]*row->J_b[2]) * bb->inv_mass;
			v3 Ia_J = inv_inertia_mul(ba->rotation, ba->inv_inertia_local, V3(row->J_a[3], row->J_a[4], row->J_a[5]));
			k += row->J_a[3]*Ia_J.x + row->J_a[4]*Ia_J.y + row->J_a[5]*Ia_J.z;
			v3 Ib_J = inv_inertia_mul(bb->rotation, bb->inv_inertia_local, V3(row->J_b[3], row->J_b[4], row->J_b[5]));
			k += row->J_b[3]*Ib_J.x + row->J_b[4]*Ib_J.y + row->J_b[5]*Ib_J.z;
			k += sys->softness[d];
			float sq = (k > 1e-12f) ? sqrtf(k) : 1.0f;
			s[d] = 1.0f / sq;
			sinv[d] = sq;
		}
	}

	// Scaled mat-vec: computes Ã*v = S*A*S*v when scaling is on, or A*v otherwise.
	#define CR_MATVEC(in, out) do { \
		if (use_scale) { \
			for (int _d = 0; _d < n; _d++) out[_d] = s[_d] * in[_d]; \
			cr_matvec(sys, w, out, out, f_lin, f_ang, body_count); \
			for (int _d = 0; _d < n; _d++) out[_d] *= s[_d]; \
		} else { \
			cr_matvec(sys, w, in, out, f_lin, f_ang, body_count); \
		} \
	} while(0)

	// Initial residual: r₀ = -J·v_current - Σ·λ_warm - bias.
	// In scaled space: r̃ = S * r.
	for (int d = 0; d < n; d++) {
		BodyHot* ba = &w->body_hot[sys->body_a[d]];
		BodyHot* bb = &w->body_hot[sys->body_b[d]];
		float vel_err = jac_velocity_f(&sys->rows[d], ba, bb);
		r[d] = -vel_err - sys->bias[d] - sys->softness[d] * sys->lambda_warm[d];
		if (use_scale) r[d] *= s[d];
	}

	// Lambda in scaled space: λ̃ = S⁻¹ * λ (so that Ã*λ̃ + b̃ = S*(A*λ + b))
	if (use_scale) {
		for (int d = 0; d < n; d++) {
			sys->lambda[d] *= sinv[d];
			sys->lambda_warm[d] *= sinv[d];
		}
	}

	CR_MATVEC(r, Ar);
	memcpy(p, r, n * sizeof(float));
	memcpy(Ap, Ar, n * sizeof(float));

	float rAr = cr_dot(r, Ar, n);
	float tol_sq = tol * tol;

	if (g_cr_trace) printf("[CR] n=%d r0=%.6e rAr0=%.6e scale=%d\n", n, sqrtf(cr_dot(r, r, n)), rAr, use_scale);

	int final_iter = 0;
	for (int iter = 0; iter < max_iters; iter++) {
		final_iter = iter;
		float ApAp = cr_dot(Ap, Ap, n);
		if (ApAp < 1e-30f) break;

		float alpha = rAr / ApAp;

		for (int i = 0; i < n; i++) {
			sys->lambda[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		float r_norm_sq = cr_dot(r, r, n);
		if (g_cr_trace) printf("  iter %2d: |r|=%.6e\n", iter, sqrtf(r_norm_sq));
		if (r_norm_sq < tol_sq) break;

		// Mid-solve clamp: project inequality DOFs and restart CR directions.
		// Clamping operates in original space: unscale, clamp, rescale.
		int reclamp_interval = w->cr_reclamp_interval;
		if (reclamp_interval > 0 && (iter + 1) % reclamp_interval == 0) {
			// Unscale lambda for clamping
			if (use_scale)
				for (int d = 0; d < n; d++) sys->lambda[d] *= s[d];
			int clamped = 0;
			for (int d = 0; d < n; d++) {
				switch (sys->clamp_type[d]) {
				case CR_CLAMP_LO_ZERO:
					if (sys->lambda[d] < 0.0f) { sys->lambda[d] = 0.0f; clamped++; }
					break;
				case CR_CLAMP_FRICTION: {
					float normal_sum = 0.0f;
					int pi = sys->friction_parent[d];
					int nc = sys->normal_count[d];
					for (int k = 0; k < nc; k++) normal_sum += sys->lambda[pi + k];
					if (normal_sum < 0.0f) normal_sum = 0.0f;
					float max_f = sys->friction_coeff[d] * normal_sum;
					float t1 = sys->lambda[d], t2 = sys->lambda[d + 1];
					float mag = sqrtf(t1 * t1 + t2 * t2);
					if (mag > max_f && mag > 1e-12f) {
						float sc = max_f / mag;
						sys->lambda[d] = t1 * sc;
						sys->lambda[d + 1] = t2 * sc;
						clamped++;
					}
					break;
				}
				case CR_CLAMP_TWIST: {
					float normal_sum = 0.0f;
					int pi = sys->friction_parent[d];
					int nc = sys->normal_count[d];
					for (int k = 0; k < nc; k++) normal_sum += sys->lambda[pi + k];
					if (normal_sum < 0.0f) normal_sum = 0.0f;
					float max_tw = sys->friction_coeff[d] * normal_sum * sys->friction_radius[d];
					if (sys->lambda[d] > max_tw) { sys->lambda[d] = max_tw; clamped++; }
					if (sys->lambda[d] < -max_tw) { sys->lambda[d] = -max_tw; clamped++; }
					break;
				}
				default: break;
				}
			}
			// Rescale lambda back to scaled space
			if (use_scale)
				for (int d = 0; d < n; d++) sys->lambda[d] *= sinv[d];
			if (clamped > 0) {
				// Recompute residual from scratch in original space, then scale.
				// lambda/lambda_warm are already in original space after unscale+clamp+rescale,
				// but we need original-space values for the rhs. Temporarily unscale.
				if (use_scale)
					for (int d = 0; d < n; d++) { sys->lambda[d] *= s[d]; sys->lambda_warm[d] *= s[d]; }
				float* delta_lam = p;
				for (int d = 0; d < n; d++)
					delta_lam[d] = sys->lambda[d] - sys->lambda_warm[d];
				cr_matvec(sys, w, delta_lam, Ar, f_lin, f_ang, body_count);
				for (int d = 0; d < n; d++) {
					BodyHot* ba = &w->body_hot[sys->body_a[d]];
					BodyHot* bb = &w->body_hot[sys->body_b[d]];
					float vel_err = jac_velocity_f(&sys->rows[d], ba, bb);
					r[d] = -vel_err - sys->bias[d] - sys->softness[d] * sys->lambda_warm[d] - Ar[d];
					if (use_scale) r[d] *= s[d];
				}
				// Rescale lambda back
				if (use_scale)
					for (int d = 0; d < n; d++) { sys->lambda[d] *= sinv[d]; sys->lambda_warm[d] *= sinv[d]; }
				// Restart CR directions
				CR_MATVEC(r, Ar);
				memcpy(p, r, n * sizeof(float));
				memcpy(Ap, Ar, n * sizeof(float));
				rAr = cr_dot(r, Ar, n);
				if (g_cr_trace) printf("  iter %d: CLAMP clamped=%d |r|=%.6e\n", iter, clamped, sqrtf(cr_dot(r, r, n)));
				continue;
			}
		}

		// Normal CR direction update
		CR_MATVEC(r, Ar);
		float rAr_new = cr_dot(r, Ar, n);
		if (fabsf(rAr) < 1e-30f) break;
		float beta = rAr_new / rAr;
		rAr = rAr_new;

		for (int i = 0; i < n; i++) {
			p[i] = r[i] + beta * p[i];
			Ap[i] = Ar[i] + beta * Ap[i];
		}
	}

	// Unscale lambda back to original space for writeback
	if (use_scale) {
		for (int d = 0; d < n; d++) {
			sys->lambda[d] *= s[d];
			sys->lambda_warm[d] *= s[d];
		}
	}

	if (g_cr_trace) {
		int n_normal = 0, n_neg_normal = 0;
		float max_neg = 0;
		for (int d = 0; d < n; d++) {
			if (sys->clamp_type[d] == CR_CLAMP_LO_ZERO) {
				n_normal++;
				if (sys->lambda[d] < 0) { n_neg_normal++; if (sys->lambda[d] < max_neg) max_neg = sys->lambda[d]; }
			}
		}
		printf("[CR] done: %d iters, %d normals, %d negative (worst=%.3f)\n",
			final_iter + 1, n_normal, n_neg_normal, max_neg);
	}

	#undef CR_MATVEC
	CK_FREE(r);
	CK_FREE(Ar);
	CK_FREE(p);
	CK_FREE(Ap);
	CK_FREE(s);
	CK_FREE(sinv);
	CK_FREE(f_lin);
	CK_FREE(f_ang);
	return final_iter + 1;
}

// --- Writeback: apply delta impulses to body velocities ---
static void cr_writeback(CR_System* sys, WorldInternal* w)
{
	int n = sys->n;
	for (int d = 0; d < n; d++) {
		float delta = sys->lambda[d] - sys->lambda_warm[d];
		*sys->lambda_ptrs[d] = sys->lambda[d];
		if (delta == 0.0f) continue;
		BodyHot* ba = &w->body_hot[sys->body_a[d]];
		BodyHot* bb = &w->body_hot[sys->body_b[d]];
		jac_apply(&sys->rows[d], delta, ba, bb);
	}
}

// --- Outer driver: solve all awake islands ---
// contacts_only: when LDL handles joints, CR only solves contact DOFs.
static void cr_velocity_solve(WorldInternal* w,
	SolverManifold* sm, int sm_count, SolverContact* sc,
	SolverJoint* sol_joints, int joint_count, float sub_dt, int contacts_only)
{
	(void)sub_dt;
	int max_iters_this_solve = 0;
	int island_count = asize(w->islands);
	for (int ii = 0; ii < island_count; ii++) {
		if (!(w->island_gen[ii] & 1)) continue;
		Island* isl = &w->islands[ii];
		if (!isl->awake) continue;

		int max_dofs = cr_count_island_dofs(w, ii, sm, sm_count, sc, sol_joints, joint_count, contacts_only);
		if (max_dofs == 0) continue;

		CR_System sys = {0};
		cr_build_island_system(&sys, w, ii, sm, sm_count, sc, sol_joints, joint_count, max_dofs, contacts_only);

		int iters = cr_island_solve(&sys, w);
		if (iters > max_iters_this_solve) max_iters_this_solve = iters;
		cr_post_clamp(&sys);
		cr_writeback(&sys, w);

		cr_free_system(&sys);
	}
	w->cr_last_iters = max_iters_this_solve;
	if (max_iters_this_solve > w->cr_peak_iters) w->cr_peak_iters = max_iters_this_solve;
}
