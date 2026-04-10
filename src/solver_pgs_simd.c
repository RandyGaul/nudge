// solver_pgs_simd.c -- SIMD PGS solver with AOS->SOA transpose gather/scatter.
// Packs 4 manifolds per SIMD group. Uses simd_transpose4 for body state
// gather (no scalar indexed loads). BEPU-style: recomputes cross+inertia inline.
// Uses v3w type and helpers from gjk_batch.c (included before this file).

// Wide 3x3 symmetric matrix: I_w = diag(xx,yy,zz) + off(xy,xz,yz).
typedef struct InertiaSIMD { v3w diag, off; } InertiaSIMD;

// Wide inertia multiply: I_w * v (macro for guaranteed inlining).
#define iw_mul(iw, v) ((v3w){ \
	simd_add(simd_add(simd_mul((iw).diag.x, (v).x), simd_mul((iw).off.x, (v).y)), simd_mul((iw).off.y, (v).z)), \
	simd_add(simd_add(simd_mul((iw).off.x, (v).x), simd_mul((iw).diag.y, (v).y)), simd_mul((iw).off.z, (v).z)), \
	simd_add(simd_add(simd_mul((iw).off.y, (v).x), simd_mul((iw).off.z, (v).y)), simd_mul((iw).diag.z, (v).z))})

// Wide cross product.
#define v3w_cross_m(a, b) ((v3w){ \
	simd_sub(simd_mul((a).y, (b).z), simd_mul((a).z, (b).y)), \
	simd_sub(simd_mul((a).z, (b).x), simd_mul((a).x, (b).z)), \
	simd_sub(simd_mul((a).x, (b).y), simd_mul((a).y, (b).x))})

// Gather a v3 field from 4 bodies using transpose (no scalar loads).
// Loads 4 v3.m (__m128) values, transposes to SOA v3w.
#define GATHER_V3(out, bodies, i0, i1, i2, i3, field) do { \
	simd4f _r0 = (bodies)[(i0)].field.m; \
	simd4f _r1 = (bodies)[(i1)].field.m; \
	simd4f _r2 = (bodies)[(i2)].field.m; \
	simd4f _r3 = (bodies)[(i3)].field.m; \
	simd_transpose4(&_r0, &_r1, &_r2, &_r3); \
	(out).x = _r0; (out).y = _r1; (out).z = _r2; \
} while(0)

// Scatter a v3w back to 4 bodies using reverse transpose.
#define SCATTER_V3(bodies, i0, i1, i2, i3, field, src) do { \
	simd4f _r0 = (src).x, _r1 = (src).y, _r2 = (src).z, _r3 = simd_zero(); \
	simd_transpose4(&_r0, &_r1, &_r2, &_r3); \
	(bodies)[(i0)].field.m = _r0; \
	(bodies)[(i1)].field.m = _r1; \
	(bodies)[(i2)].field.m = _r2; \
	(bodies)[(i3)].field.m = _r3; \
} while(0)

// Pack v3 from 4 sources into v3w using simd_set.
#define V3W_PACK(dst, v0, v1, v2, v3) do { (dst).x = simd_set((v0).x,(v1).x,(v2).x,(v3).x); (dst).y = simd_set((v0).y,(v1).y,(v2).y,(v3).y); (dst).z = simd_set((v0).z,(v1).z,(v2).z,(v3).z); } while(0)

// Lean SIMD group: 4 manifolds packed for parallel solve.
// No precomputed angular data — everything recomputed inline (BEPU style).
typedef struct SolverGroupSIMD
{
	int body_a[4];
	int body_b[4];
	int manifold_idx[4];
	int contact_start[4];
	int contact_count[4];
	int max_contacts;

	simd4f inv_mass_a, inv_mass_b;

	// Per-contact prestep (only what the solver reads): r_a, r_b, eff_mass, bias, bounce, softness, lambda
	struct { v3w r_a, r_b; simd4f eff_mass, bias, bounce, softness, lambda; } normals[MAX_CONTACTS];

	// Manifold-level friction prestep
	v3w tangent1, tangent2, normal;
	v3w centroid_r_a, centroid_r_b;
	simd4f eff_mass_t1, eff_mass_t2, eff_mass_twist;
	simd4f lambda_t1, lambda_t2, lambda_twist;
	simd4f friction, patch_radius;
} SolverGroupSIMD;

// Build SIMD groups from scalar manifolds within a color batch.
static SolverGroupSIMD* simd_prepare_groups(WorldInternal* w, SolverManifold* sm, SolverContact* sc, ConstraintRef* crefs, int start, int end, int* out_count)
{
	CK_DYNA SolverGroupSIMD* groups = NULL;

	CK_DYNA int* contact_crefs = NULL;
	for (int i = start; i < end; i++)
		if (crefs[i].type == CTYPE_CONTACT) apush(contact_crefs, crefs[i].index);
	while (asize(contact_crefs) % 4 != 0) apush(contact_crefs, -1);

	v3 zv = V3(0,0,0);
	for (int base = 0; base < asize(contact_crefs); base += 4) {
		SolverGroupSIMD g = {0};
		SolverManifold* ms[4] = {0};
		float ima[4] = {0}, imb[4] = {0};
		int max_cc = 0;

		for (int lane = 0; lane < 4; lane++) {
			int mi = contact_crefs[base + lane];
			if (mi < 0) { g.body_a[lane] = 0; g.body_b[lane] = 0; g.manifold_idx[lane] = -1; g.contact_start[lane] = -1; continue; }
			ms[lane] = &sm[mi];
			g.body_a[lane] = ms[lane]->body_a;
			g.body_b[lane] = ms[lane]->body_b;
			g.manifold_idx[lane] = mi;
			g.contact_start[lane] = ms[lane]->contact_start;
			g.contact_count[lane] = ms[lane]->contact_count;
			ima[lane] = ms[lane]->inv_mass_a;
			imb[lane] = ms[lane]->inv_mass_b;
			if (ms[lane]->contact_count > max_cc) max_cc = ms[lane]->contact_count;
		}
		g.inv_mass_a = simd_load(ima); g.inv_mass_b = simd_load(imb); g.max_contacts = max_cc;

		// Pack normal contact prestep
		for (int slot = 0; slot < max_cc; slot++) {
			v3 ra[4], rb[4]; float em[4]={0}, bi[4]={0}, bo[4]={0}, sf[4]={0}, lm[4]={0};
			for (int lane = 0; lane < 4; lane++) {
				if (!ms[lane] || slot >= ms[lane]->contact_count) { ra[lane]=rb[lane]=zv; continue; }
				SolverContact* s = &sc[ms[lane]->contact_start + slot];
				ra[lane]=s->r_a; rb[lane]=s->r_b;
				em[lane]=s->eff_mass_n; bi[lane]=s->bias; bo[lane]=s->bounce; sf[lane]=s->softness; lm[lane]=s->lambda_n;
			}
			V3W_PACK(g.normals[slot].r_a, ra[0], ra[1], ra[2], ra[3]);
			V3W_PACK(g.normals[slot].r_b, rb[0], rb[1], rb[2], rb[3]);
			g.normals[slot].eff_mass = simd_load(em); g.normals[slot].bias = simd_load(bi);
			g.normals[slot].bounce = simd_load(bo); g.normals[slot].softness = simd_load(sf);
			g.normals[slot].lambda = simd_load(lm);
		}

		// Pack manifold-level friction
		v3 t1[4], t2[4], nm[4], cra[4], crb[4];
		float emt1[4]={0}, emt2[4]={0}, emtw[4]={0}, lt1[4]={0}, lt2[4]={0}, ltw[4]={0}, fr[4]={0}, pr[4]={0};
		for (int lane = 0; lane < 4; lane++) {
			if (!ms[lane]) { t1[lane]=t2[lane]=nm[lane]=cra[lane]=crb[lane]=zv; continue; }
			SolverManifold* m = ms[lane];
			t1[lane]=m->tangent1; t2[lane]=m->tangent2; nm[lane]=m->normal;
			cra[lane]=m->centroid_r_a; crb[lane]=m->centroid_r_b;
			emt1[lane]=m->eff_mass_t1; emt2[lane]=m->eff_mass_t2; emtw[lane]=m->eff_mass_twist;
			lt1[lane]=m->lambda_t1; lt2[lane]=m->lambda_t2; ltw[lane]=m->lambda_twist;
			fr[lane]=m->friction; pr[lane]=m->patch_radius;
		}
		V3W_PACK(g.tangent1, t1[0], t1[1], t1[2], t1[3]);
		V3W_PACK(g.tangent2, t2[0], t2[1], t2[2], t2[3]);
		V3W_PACK(g.normal, nm[0], nm[1], nm[2], nm[3]);
		V3W_PACK(g.centroid_r_a, cra[0], cra[1], cra[2], cra[3]);
		V3W_PACK(g.centroid_r_b, crb[0], crb[1], crb[2], crb[3]);
		g.eff_mass_t1=simd_load(emt1); g.eff_mass_t2=simd_load(emt2); g.eff_mass_twist=simd_load(emtw);
		g.lambda_t1=simd_load(lt1); g.lambda_t2=simd_load(lt2); g.lambda_twist=simd_load(ltw);
		g.friction=simd_load(fr); g.patch_radius=simd_load(pr);

		apush(groups, g);
	}
	afree(contact_crefs);
	*out_count = asize(groups);
	return groups;
}

// Solve one SIMD group: 4 manifolds in parallel.
// Gathers body state via transpose (not scalar indexing).
// Recomputes cross+inertia inline (BEPU style).
static void simd_solve_group(SolverBodyVel* bodies, SolverGroupSIMD* g)
{
	int i0a = g->body_a[0], i1a = g->body_a[1], i2a = g->body_a[2], i3a = g->body_a[3];
	int i0b = g->body_b[0], i1b = g->body_b[1], i2b = g->body_b[2], i3b = g->body_b[3];

	// Gather body state via transpose (4 __m128 loads + 1 transpose per v3 field)
	v3w va, wa, vb, wb;
	InertiaSIMD iw_a, iw_b;
	GATHER_V3(va, bodies, i0a, i1a, i2a, i3a, velocity);
	GATHER_V3(wa, bodies, i0a, i1a, i2a, i3a, angular_velocity);
	GATHER_V3(iw_a.diag, bodies, i0a, i1a, i2a, i3a, iw_diag);
	GATHER_V3(iw_a.off, bodies, i0a, i1a, i2a, i3a, iw_off);
	GATHER_V3(vb, bodies, i0b, i1b, i2b, i3b, velocity);
	GATHER_V3(wb, bodies, i0b, i1b, i2b, i3b, angular_velocity);
	GATHER_V3(iw_b.diag, bodies, i0b, i1b, i2b, i3b, iw_diag);
	GATHER_V3(iw_b.off, bodies, i0b, i1b, i2b, i3b, iw_off);

	simd4f zero = simd_zero();
	simd4f ima = g->inv_mass_a, imb = g->inv_mass_b;
	simd4f inv_mass_sum = simd_add(ima, imb);
	v3w normal = g->normal;

	// Normal contacts — recompute cross+inertia inline
	simd4f linear_vn = v3w_dot(v3w_sub(vb, va), normal);
	simd4f total_lambda_n = zero;
	for (int slot = 0; slot < g->max_contacts; slot++) {
		v3w rn_a = v3w_cross_m(g->normals[slot].r_a, normal);
		v3w rn_b = v3w_cross_m(g->normals[slot].r_b, normal);
		simd4f vn = simd_add(linear_vn, simd_sub(v3w_dot(wb, rn_b), v3w_dot(wa, rn_a)));
		simd4f rhs = simd_add(simd_add(vn, g->normals[slot].bias), g->normals[slot].bounce);
		simd4f lambda_delta = simd_mul(g->normals[slot].eff_mass, simd_sub(simd_neg(rhs), simd_mul(g->normals[slot].softness, g->normals[slot].lambda)));
		simd4f old = g->normals[slot].lambda;
		g->normals[slot].lambda = simd_max(simd_add(old, lambda_delta), zero);
		simd4f delta = simd_sub(g->normals[slot].lambda, old);

		// Apply impulse — recompute angular inline
		va = v3w_sub(va, v3w_scale(normal, simd_mul(delta, ima)));
		vb = v3w_add(vb, v3w_scale(normal, simd_mul(delta, imb)));
		wa = v3w_sub(wa, v3w_scale(iw_mul(iw_a, rn_a), delta));
		wb = v3w_add(wb, v3w_scale(iw_mul(iw_b, rn_b), delta));
		linear_vn = simd_add(linear_vn, simd_mul(delta, inv_mass_sum));
		total_lambda_n = simd_add(total_lambda_n, g->normals[slot].lambda);
	}

	// Manifold-level friction
	simd4f max_f = simd_mul(g->friction, total_lambda_n);
	simd4f neg_max_f = simd_neg(max_f);

	// Tangent 1 — recompute
	v3w rct1_a = v3w_cross_m(g->centroid_r_a, g->tangent1);
	v3w rct1_b = v3w_cross_m(g->centroid_r_b, g->tangent1);
	simd4f vt1 = simd_add(v3w_dot(v3w_sub(vb, va), g->tangent1), simd_sub(v3w_dot(wb, rct1_b), v3w_dot(wa, rct1_a)));
	simd4f old_t1 = g->lambda_t1;
	g->lambda_t1 = simd_max(neg_max_f, simd_min(simd_add(old_t1, simd_mul(g->eff_mass_t1, simd_neg(vt1))), max_f));
	simd4f dt1 = simd_sub(g->lambda_t1, old_t1);
	va = v3w_sub(va, v3w_scale(g->tangent1, simd_mul(dt1, ima)));
	vb = v3w_add(vb, v3w_scale(g->tangent1, simd_mul(dt1, imb)));
	wa = v3w_sub(wa, v3w_scale(iw_mul(iw_a, rct1_a), dt1));
	wb = v3w_add(wb, v3w_scale(iw_mul(iw_b, rct1_b), dt1));

	// Tangent 2 — recompute
	v3w rct2_a = v3w_cross_m(g->centroid_r_a, g->tangent2);
	v3w rct2_b = v3w_cross_m(g->centroid_r_b, g->tangent2);
	simd4f vt2 = simd_add(v3w_dot(v3w_sub(vb, va), g->tangent2), simd_sub(v3w_dot(wb, rct2_b), v3w_dot(wa, rct2_a)));
	simd4f old_t2 = g->lambda_t2;
	g->lambda_t2 = simd_max(neg_max_f, simd_min(simd_add(old_t2, simd_mul(g->eff_mass_t2, simd_neg(vt2))), max_f));
	simd4f dt2 = simd_sub(g->lambda_t2, old_t2);
	va = v3w_sub(va, v3w_scale(g->tangent2, simd_mul(dt2, ima)));
	vb = v3w_add(vb, v3w_scale(g->tangent2, simd_mul(dt2, imb)));
	wa = v3w_sub(wa, v3w_scale(iw_mul(iw_a, rct2_a), dt2));
	wb = v3w_add(wb, v3w_scale(iw_mul(iw_b, rct2_b), dt2));

	// Torsional friction — recompute
	simd4f max_twist = simd_mul(max_f, g->patch_radius);
	simd4f w_rel = v3w_dot(v3w_sub(wb, wa), normal);
	simd4f old_tw = g->lambda_twist;
	g->lambda_twist = simd_max(simd_neg(max_twist), simd_min(simd_add(old_tw, simd_mul(g->eff_mass_twist, simd_neg(w_rel))), max_twist));
	simd4f dtw = simd_sub(g->lambda_twist, old_tw);
	wa = v3w_sub(wa, v3w_scale(iw_mul(iw_a, normal), dtw));
	wb = v3w_add(wb, v3w_scale(iw_mul(iw_b, normal), dtw));

	// Scatter body state via reverse transpose
	SCATTER_V3(bodies, i0a, i1a, i2a, i3a, velocity, va);
	SCATTER_V3(bodies, i0a, i1a, i2a, i3a, angular_velocity, wa);
	SCATTER_V3(bodies, i0b, i1b, i2b, i3b, velocity, vb);
	SCATTER_V3(bodies, i0b, i1b, i2b, i3b, angular_velocity, wb);
}

// Copy SIMD lambda values back to scalar structs for post_solve / warm cache.
static void simd_writeback_lambdas(SolverGroupSIMD* groups, int group_count, SolverManifold* sm, SolverContact* sc)
{
	for (int gi = 0; gi < group_count; gi++) {
		SolverGroupSIMD* g = &groups[gi];
		float lt1[4], lt2[4], ltw[4];
		simd_store(lt1, g->lambda_t1);
		simd_store(lt2, g->lambda_t2);
		simd_store(ltw, g->lambda_twist);
		for (int lane = 0; lane < 4; lane++) {
			int mi = g->manifold_idx[lane];
			if (mi < 0) continue;
			SolverManifold* m = &sm[mi];
			m->lambda_t1 = lt1[lane]; m->lambda_t2 = lt2[lane]; m->lambda_twist = ltw[lane];
			for (int c = 0; c < m->contact_count; c++) {
				float lam[4]; simd_store(lam, g->normals[c].lambda);
				sc[m->contact_start + c].lambda_n = lam[lane];
			}
		}
	}
}

// Refresh SIMD bias values from scalar contacts (after solver_relax_contacts).
static void simd_refresh_bias(SolverGroupSIMD* groups, int group_count, SolverManifold* sm, SolverContact* sc)
{
	for (int gi = 0; gi < group_count; gi++) {
		SolverGroupSIMD* g = &groups[gi];
		for (int slot = 0; slot < g->max_contacts; slot++) {
			float bi[4] = {0}, bo[4] = {0};
			for (int lane = 0; lane < 4; lane++) {
				int mi = g->manifold_idx[lane];
				if (mi < 0 || slot >= g->contact_count[lane]) continue;
				SolverContact* s = &sc[sm[mi].contact_start + slot];
				bi[lane] = s->bias; bo[lane] = s->bounce;
			}
			g->normals[slot].bias = simd_load(bi);
			g->normals[slot].bounce = simd_load(bo);
		}
	}
}
