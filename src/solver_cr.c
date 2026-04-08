// solver_cr.c -- Conjugate Residual
//
// Pure linear Krylov accelerator for the constraint-space Delassus system.
// CR handles normals, friction, and joints (when LDL is off) as a single
// unconstrained linear system A = J M^-1 J^T + Sigma. No clamping or
// inequality projection — PGS handles that before and after CR.
//
// Pipeline: PGS warmup (constraints) -> CR (global propagation)

#define CR_MIN_COMPLIANCE 5e-5f

static int g_cr_trace = 0; // set to 1 for per-island diagnostic output

typedef struct CR_System
{
	int n;
	JacobianRow* rows;   // [n] Jacobian per DOF
	int* body_a;         // [n] body A index
	int* body_b;         // [n] body B index
	float* softness;     // [n] compliance
	float* bias;         // [n] RHS bias
	float* lambda;       // [n] working solution
	float* lambda_warm;  // [n] snapshot before CR
	float** lambda_ptrs; // [n] writeback pointers
} CR_System;

static JacobianRow cr_build_contact_jacobian(v3 dir, v3 r_a, v3 r_b)
{
	v3 cr_a = cross(r_a, dir);
	v3 cr_b = cross(r_b, dir);
	JacobianRow row;
	row.J_a[0] = -dir.x; row.J_a[1] = -dir.y; row.J_a[2] = -dir.z;
	row.J_a[3] = -cr_a.x; row.J_a[4] = -cr_a.y; row.J_a[5] = -cr_a.z;
	row.J_b[0] =  dir.x; row.J_b[1] =  dir.y; row.J_b[2] =  dir.z;
	row.J_b[3] =  cr_b.x; row.J_b[4] =  cr_b.y; row.J_b[5] =  cr_b.z;
	row.eff_mass = 0;
	return row;
}


static void cr_add_dof(CR_System* sys, JacobianRow jac, int ba, int bb, float soft, float b, float lam, float* lam_ptr)
{
	int i = sys->n++;
	sys->rows[i] = jac;
	sys->body_a[i] = ba;
	sys->body_b[i] = bb;
	sys->softness[i] = soft > 0.0f ? soft : CR_MIN_COMPLIANCE;
	sys->bias[i] = b;
	sys->lambda[i] = lam;
	sys->lambda_warm[i] = lam;
	sys->lambda_ptrs[i] = lam_ptr;
}

static int cr_count_island_dofs(WorldInternal* w, int island_idx, SolverManifold* sm, int sm_count, SolverContact* sc, SolverJoint* sol_joints, int joint_count, int contacts_only)
{
	Island* isl = &w->islands[island_idx];
	int n = 0;
	if (!contacts_only) {
		int ji = isl->head_joint;
		while (ji >= 0) {
			for (int i = 0; i < joint_count; i++) {
				if (sol_joints[i].joint_idx == ji) { n += sol_joints[i].dof; break; }
			}
			ji = w->joints[ji].island_next;
		}
	}

	for (int i = 0; i < sm_count; i++) {
		int isl_a = w->body_cold[sm[i].body_a].island_id;
		int isl_b = w->body_cold[sm[i].body_b].island_id;
		if (isl_a != island_idx && isl_b != island_idx) continue;
		n += sm[i].contact_count;
	}
	return n;
}

static void cr_build_island_system(CR_System* sys, WorldInternal* w, int island_idx, SolverManifold* sm, int sm_count, SolverContact* sc, SolverJoint* sol_joints, int joint_count, int max_dofs, int contacts_only)
{
	sys->n = 0;
	sys->rows = CK_ALLOC(max_dofs * sizeof(JacobianRow));
	sys->body_a = CK_ALLOC(max_dofs * sizeof(int));
	sys->body_b = CK_ALLOC(max_dofs * sizeof(int));
	sys->softness = CK_ALLOC(max_dofs * sizeof(float));
	sys->bias = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda_warm = CK_ALLOC(max_dofs * sizeof(float));
	sys->lambda_ptrs = CK_ALLOC(max_dofs * sizeof(float*));

	Island* isl = &w->islands[island_idx];

	// Joints (skipped when LDL handles them)
	if (!contacts_only) {
		int ji = isl->head_joint;
		while (ji >= 0) {
			for (int i = 0; i < joint_count; i++) {
				if (sol_joints[i].joint_idx != ji) continue;
				SolverJoint* sj = &sol_joints[i];
				for (int d = 0; d < sj->dof; d++)
					cr_add_dof(sys, sj->rows[d], sj->body_a, sj->body_b, sj->softness, sj->bias[d], sj->lambda[d], &sj->lambda[d]);
				break;
			}
			ji = w->joints[ji].island_next;
		}
	}

	// Contact normals only. Friction handled by PGS.
	int mask = w->cr_active_set_mask;
	for (int mi = 0; mi < sm_count; mi++) {
		SolverManifold* m = &sm[mi];
		int isl_a = w->body_cold[m->body_a].island_id;
		int isl_b = w->body_cold[m->body_b].island_id;
		if (isl_a != island_idx && isl_b != island_idx) continue;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (mask && s->lambda_n == 0.0f) continue;
			cr_add_dof(sys, cr_build_contact_jacobian(s->normal, s->r_a, s->r_b), m->body_a, m->body_b, s->softness, s->bias + s->bounce, s->lambda_n, &s->lambda_n);
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
	CK_FREE(sys->lambda_ptrs);
	memset(sys, 0, sizeof(*sys));
}

static float cr_dot(float* a, float* b, int n)
{
	float s = 0.0f;
	for (int i = 0; i < n; i++) s += a[i] * b[i];
	return s;
}

// Matrix-free mat-vec: out = (J M^{-1} J^T + Sigma) * v
static void cr_matvec(CR_System* sys, WorldInternal* w, float* v, float* out, v3* f_lin, v3* f_ang, int body_count)
{
	int n = sys->n;
	memset(f_lin, 0, body_count * sizeof(v3));
	memset(f_ang, 0, body_count * sizeof(v3));

	for (int d = 0; d < n; d++) {
		float lam = v[d];
		if (lam == 0.0f) continue;
		JacobianRow* row = &sys->rows[d];
		int ba = sys->body_a[d], bb = sys->body_b[d];
		f_lin[ba].x += row->J_a[0] * lam; f_lin[ba].y += row->J_a[1] * lam; f_lin[ba].z += row->J_a[2] * lam;
		f_ang[ba].x += row->J_a[3] * lam; f_ang[ba].y += row->J_a[4] * lam; f_ang[ba].z += row->J_a[5] * lam;
		f_lin[bb].x += row->J_b[0] * lam; f_lin[bb].y += row->J_b[1] * lam; f_lin[bb].z += row->J_b[2] * lam;
		f_ang[bb].x += row->J_b[3] * lam; f_ang[bb].y += row->J_b[4] * lam; f_ang[bb].z += row->J_b[5] * lam;
	}

	for (int i = 0; i < body_count; i++) {
		BodyHot* bh = &w->body_hot[i];
		if (bh->inv_mass == 0.0f) { f_lin[i] = V3(0,0,0); f_ang[i] = V3(0,0,0); continue; }
		if (f_lin[i].x == 0 && f_lin[i].y == 0 && f_lin[i].z == 0 && f_ang[i].x == 0 && f_ang[i].y == 0 && f_ang[i].z == 0) continue;
		f_lin[i] = scale(f_lin[i], bh->inv_mass);
		f_ang[i] = inv_inertia_mul(bh->rotation, bh->inv_inertia_local, f_ang[i]);
	}

	for (int d = 0; d < n; d++) {
		JacobianRow* row = &sys->rows[d];
		int ba = sys->body_a[d], bb = sys->body_b[d];
		v3 al = f_lin[ba], aa = f_ang[ba], bl = f_lin[bb], ba2 = f_ang[bb];
		float r = row->J_a[0]*al.x + row->J_a[1]*al.y + row->J_a[2]*al.z + row->J_a[3]*aa.x + row->J_a[4]*aa.y + row->J_a[5]*aa.z + row->J_b[0]*bl.x + row->J_b[1]*bl.y + row->J_b[2]*bl.z + row->J_b[3]*ba2.x + row->J_b[4]*ba2.y + row->J_b[5]*ba2.z;
		out[d] = r + sys->softness[d] * v[d];
	}
}

// Pure unconstrained CR solve. No clamping — PGS handles constraints before/after.
static int cr_island_solve(CR_System* sys, WorldInternal* w)
{
	int n = sys->n;
	if (n == 0) return 0;

	int body_count = asize(w->body_hot);
	int max_iters = w->cr_max_iters;
	float tol = w->cr_tolerance;

	float* r   = CK_ALLOC(n * sizeof(float));
	float* Ar  = CK_ALLOC(n * sizeof(float));
	float* p   = CK_ALLOC(n * sizeof(float));
	float* Ap  = CK_ALLOC(n * sizeof(float));
	v3* f_lin  = CK_ALLOC(body_count * sizeof(v3));
	v3* f_ang  = CK_ALLOC(body_count * sizeof(v3));

	// Initial residual: r = -J*v - bias - sigma*lambda_warm
	for (int d = 0; d < n; d++) {
		BodyHot* ba = &w->body_hot[sys->body_a[d]];
		BodyHot* bb = &w->body_hot[sys->body_b[d]];
		r[d] = -jac_velocity_f(&sys->rows[d], ba, bb) - sys->bias[d] - sys->softness[d] * sys->lambda_warm[d];
	}

	cr_matvec(sys, w, r, Ar, f_lin, f_ang, body_count);
	memcpy(p, r, n * sizeof(float));
	memcpy(Ap, Ar, n * sizeof(float));

	float rAr = cr_dot(r, Ar, n);
	float tol_sq = tol * tol;

	if (g_cr_trace) printf("[CR] n=%d r0=%.6e rAr0=%.6e\n", n, sqrtf(cr_dot(r, r, n)), rAr);

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

		cr_matvec(sys, w, r, Ar, f_lin, f_ang, body_count);
		float rAr_new = cr_dot(r, Ar, n);
		if (fabsf(rAr) < 1e-30f) break;
		float beta = rAr_new / rAr;
		rAr = rAr_new;

		for (int i = 0; i < n; i++) {
			p[i] = r[i] + beta * p[i];
			Ap[i] = Ar[i] + beta * Ap[i];
		}
	}

	if (g_cr_trace) printf("[CR] done: %d iters\n", final_iter + 1);

	CK_FREE(r);
	CK_FREE(Ar);
	CK_FREE(p);
	CK_FREE(Ap);
	CK_FREE(f_lin);
	CK_FREE(f_ang);
	return final_iter + 1;
}

static void cr_writeback(CR_System* sys, WorldInternal* w)
{
	for (int d = 0; d < sys->n; d++) {
		float delta = sys->lambda[d] - sys->lambda_warm[d];
		*sys->lambda_ptrs[d] = sys->lambda[d];
		if (delta == 0.0f) continue;
		BodyHot* ba = &w->body_hot[sys->body_a[d]];
		BodyHot* bb = &w->body_hot[sys->body_b[d]];
		jac_apply(&sys->rows[d], delta, ba, bb);
	}
}

static void cr_velocity_solve(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, SolverJoint* sol_joints, int joint_count, float sub_dt, int contacts_only)
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
		cr_writeback(&sys, w);

		cr_free_system(&sys);
	}
	w->cr_last_iters = max_iters_this_solve;
	if (max_iters_this_solve > w->cr_peak_iters) w->cr_peak_iters = max_iters_this_solve;
}
