// solver_pgs_simd.c -- SIMD-accelerated PGS solver for patch friction contacts.
// Packs 4 manifolds per SIMD group, solves normal + friction rows in parallel.
// Follows Box2D v3 pattern: SOA layout, gather/scatter body state per group.
// Uses v3w type and helpers from gjk_batch.c (included before this file).

// Per-contact-slot normal data (4 manifolds in parallel, 1 contact index).
typedef struct NormalSlotSIMD
{
	v3w normal;
	v3w rn_a, rn_b;      // cross(r, normal)
	v3w w_n_a, w_n_b;    // I_w * cross(r, normal)
	simd4f eff_mass;
	simd4f bias;
	simd4f bounce;
	simd4f softness;
	simd4f lambda;
} NormalSlotSIMD;

// SIMD group: 4 manifolds packed for parallel solve.
typedef struct SolverGroupSIMD
{
	int body_a[4];
	int body_b[4];
	int manifold_idx[4]; // index into scalar SolverManifold array (for lambda writeback)
	int contact_start[4]; // index into scalar SolverContact array

	simd4f inv_mass_a, inv_mass_b;
	int max_contacts; // max contact_count across all 4 manifolds

	NormalSlotSIMD normals[MAX_CONTACTS];

	// Manifold-level patch friction
	v3w tangent1, tangent2, normal;
	v3w rct1_a, rct1_b, rct2_a, rct2_b;
	v3w w_t1_a, w_t1_b, w_t2_a, w_t2_b;
	v3w w_tw_a, w_tw_b;
	simd4f eff_mass_t1, eff_mass_t2, eff_mass_twist;
	simd4f lambda_t1, lambda_t2, lambda_twist;
	simd4f friction;
	simd4f patch_radius;
} SolverGroupSIMD;

// Gather body velocities from body_hot[] into wide vectors.
static void simd_gather_velocities(BodyHot* bodies, int idx[4], v3w* vel, v3w* ang)
{
	float vx[4], vy[4], vz[4], wx[4], wy[4], wz[4];
	for (int i = 0; i < 4; i++) {
		BodyHot* b = &bodies[idx[i]];
		vx[i] = b->velocity.x; vy[i] = b->velocity.y; vz[i] = b->velocity.z;
		wx[i] = b->angular_velocity.x; wy[i] = b->angular_velocity.y; wz[i] = b->angular_velocity.z;
	}
	vel->x = simd_load(vx); vel->y = simd_load(vy); vel->z = simd_load(vz);
	ang->x = simd_load(wx); ang->y = simd_load(wy); ang->z = simd_load(wz);
}

// Scatter body velocities back to body_hot[].
static void simd_scatter_velocities(BodyHot* bodies, int idx[4], v3w vel, v3w ang)
{
	float vx[4], vy[4], vz[4], wx[4], wy[4], wz[4];
	simd_store(vx, vel.x); simd_store(vy, vel.y); simd_store(vz, vel.z);
	simd_store(wx, ang.x); simd_store(wy, ang.y); simd_store(wz, ang.z);
	for (int i = 0; i < 4; i++) {
		BodyHot* b = &bodies[idx[i]];
		b->velocity = V3(vx[i], vy[i], vz[i]);
		b->angular_velocity = V3(wx[i], wy[i], wz[i]);
	}
}

// Pack a v3 scalar into lane i of a v3w.
static inline void v3w_set_lane(v3w* w, int lane, v3 v)
{
	float buf[4];
	simd_store(buf, w->x); buf[lane] = v.x; w->x = simd_load(buf);
	simd_store(buf, w->y); buf[lane] = v.y; w->y = simd_load(buf);
	simd_store(buf, w->z); buf[lane] = v.z; w->z = simd_load(buf);
}

// Set lane of a simd4f.
static inline simd4f simd_set_lane(simd4f v, int lane, float val)
{
	float buf[4];
	simd_store(buf, v);
	buf[lane] = val;
	return simd_load(buf);
}

// Build SIMD groups from scalar manifolds within a color batch.
// Returns dynamic array of SolverGroupSIMD (caller frees).
static SolverGroupSIMD* simd_prepare_groups(WorldInternal* w, SolverManifold* sm, SolverContact* sc, ConstraintRef* crefs, int start, int end, int* out_count)
{
	CK_DYNA SolverGroupSIMD* groups = NULL;
	SolverGroupSIMD cur = {0};
	int lane = 0;

	// Sentinel body index: use index 0 (always exists as static floor or first body).
	// Inactive lanes get body_a=body_b=0 with inv_mass=0, so impulses are no-ops.
	for (int i = start; i <= end; i++) {
		if (i < end && crefs[i].type == CTYPE_CONTACT) {
			int mi = crefs[i].index;
			SolverManifold* m = &sm[mi];
			BodyHot* a = &w->body_hot[m->body_a];
			BodyHot* b = &w->body_hot[m->body_b];

			cur.body_a[lane] = m->body_a;
			cur.body_b[lane] = m->body_b;
			cur.manifold_idx[lane] = mi;
			cur.contact_start[lane] = m->contact_start;
			cur.inv_mass_a = simd_set_lane(cur.inv_mass_a, lane, a->inv_mass);
			cur.inv_mass_b = simd_set_lane(cur.inv_mass_b, lane, b->inv_mass);

			// Normal contacts
			for (int c = 0; c < m->contact_count; c++) {
				SolverContact* s = &sc[m->contact_start + c];
				NormalSlotSIMD* n = &cur.normals[c];
				v3w_set_lane(&n->normal, lane, s->normal);
				v3w_set_lane(&n->rn_a, lane, s->rn_a);
				v3w_set_lane(&n->rn_b, lane, s->rn_b);
				v3w_set_lane(&n->w_n_a, lane, s->w_n_a);
				v3w_set_lane(&n->w_n_b, lane, s->w_n_b);
				n->eff_mass = simd_set_lane(n->eff_mass, lane, s->eff_mass_n);
				n->bias = simd_set_lane(n->bias, lane, s->bias);
				n->bounce = simd_set_lane(n->bounce, lane, s->bounce);
				n->softness = simd_set_lane(n->softness, lane, s->softness);
				n->lambda = simd_set_lane(n->lambda, lane, s->lambda_n);
			}
			if (m->contact_count > cur.max_contacts) cur.max_contacts = m->contact_count;

			// Manifold-level friction
			v3w_set_lane(&cur.tangent1, lane, m->tangent1);
			v3w_set_lane(&cur.tangent2, lane, m->tangent2);
			v3w_set_lane(&cur.normal, lane, m->normal);
			v3w_set_lane(&cur.rct1_a, lane, m->rct1_a);
			v3w_set_lane(&cur.rct1_b, lane, m->rct1_b);
			v3w_set_lane(&cur.rct2_a, lane, m->rct2_a);
			v3w_set_lane(&cur.rct2_b, lane, m->rct2_b);
			v3w_set_lane(&cur.w_t1_a, lane, m->w_t1_a);
			v3w_set_lane(&cur.w_t1_b, lane, m->w_t1_b);
			v3w_set_lane(&cur.w_t2_a, lane, m->w_t2_a);
			v3w_set_lane(&cur.w_t2_b, lane, m->w_t2_b);
			v3w_set_lane(&cur.w_tw_a, lane, m->w_tw_a);
			v3w_set_lane(&cur.w_tw_b, lane, m->w_tw_b);
			cur.eff_mass_t1 = simd_set_lane(cur.eff_mass_t1, lane, m->eff_mass_t1);
			cur.eff_mass_t2 = simd_set_lane(cur.eff_mass_t2, lane, m->eff_mass_t2);
			cur.eff_mass_twist = simd_set_lane(cur.eff_mass_twist, lane, m->eff_mass_twist);
			cur.lambda_t1 = simd_set_lane(cur.lambda_t1, lane, m->lambda_t1);
			cur.lambda_t2 = simd_set_lane(cur.lambda_t2, lane, m->lambda_t2);
			cur.lambda_twist = simd_set_lane(cur.lambda_twist, lane, m->lambda_twist);
			cur.friction = simd_set_lane(cur.friction, lane, m->friction);
			cur.patch_radius = simd_set_lane(cur.patch_radius, lane, m->patch_radius);

			lane++;
		}

		// Flush group when full or at end of batch
		if (lane == 4 || (i == end && lane > 0)) {
			// Pad unused lanes with sentinel body 0 (inv_mass=0 → no-op)
			for (int p = lane; p < 4; p++) {
				cur.body_a[p] = 0;
				cur.body_b[p] = 0;
				cur.manifold_idx[p] = -1;
				cur.contact_start[p] = -1;
			}
			apush(groups, cur);
			memset(&cur, 0, sizeof(cur));
			lane = 0;
		}
	}

	*out_count = asize(groups);
	return groups;
}

// Solve one SIMD group: 4 manifolds in parallel.
static void simd_solve_group(BodyHot* bodies, SolverGroupSIMD* g)
{
	v3w va, wa, vb, wb;
	simd_gather_velocities(bodies, g->body_a, &va, &wa);
	simd_gather_velocities(bodies, g->body_b, &vb, &wb);

	simd4f zero = simd_zero();

	// Normal contacts
	for (int slot = 0; slot < g->max_contacts; slot++) {
		NormalSlotSIMD* n = &g->normals[slot];

		simd4f vn = simd_add(v3w_dot(v3w_sub(vb, va), n->normal), simd_sub(v3w_dot(wb, n->rn_b), v3w_dot(wa, n->rn_a)));
		simd4f rhs = simd_add(simd_add(vn, n->bias), n->bounce);
		simd4f lambda_delta = simd_mul(n->eff_mass, simd_sub(simd_neg(rhs), simd_mul(n->softness, n->lambda)));
		simd4f old = n->lambda;
		n->lambda = simd_max(simd_add(old, lambda_delta), zero);
		simd4f delta = simd_sub(n->lambda, old);

		va = v3w_sub(va, v3w_scale(n->normal, simd_mul(delta, g->inv_mass_a)));
		vb = v3w_add(vb, v3w_scale(n->normal, simd_mul(delta, g->inv_mass_b)));
		wa = v3w_sub(wa, v3w_scale(n->w_n_a, delta));
		wb = v3w_add(wb, v3w_scale(n->w_n_b, delta));
	}

	// Aggregate normal lambda for friction limit
	simd4f total_lambda_n = zero;
	for (int slot = 0; slot < g->max_contacts; slot++)
		total_lambda_n = simd_add(total_lambda_n, g->normals[slot].lambda);

	simd4f max_f = simd_mul(g->friction, total_lambda_n);
	simd4f neg_max_f = simd_neg(max_f);

	// Tangent 1
	simd4f vt1 = simd_add(v3w_dot(v3w_sub(vb, va), g->tangent1), simd_sub(v3w_dot(wb, g->rct1_b), v3w_dot(wa, g->rct1_a)));
	simd4f old_t1 = g->lambda_t1;
	g->lambda_t1 = simd_max(neg_max_f, simd_min(simd_add(old_t1, simd_mul(g->eff_mass_t1, simd_neg(vt1))), max_f));
	simd4f dt1 = simd_sub(g->lambda_t1, old_t1);
	va = v3w_sub(va, v3w_scale(g->tangent1, simd_mul(dt1, g->inv_mass_a)));
	vb = v3w_add(vb, v3w_scale(g->tangent1, simd_mul(dt1, g->inv_mass_b)));
	wa = v3w_sub(wa, v3w_scale(g->w_t1_a, dt1));
	wb = v3w_add(wb, v3w_scale(g->w_t1_b, dt1));

	// Tangent 2
	simd4f vt2 = simd_add(v3w_dot(v3w_sub(vb, va), g->tangent2), simd_sub(v3w_dot(wb, g->rct2_b), v3w_dot(wa, g->rct2_a)));
	simd4f old_t2 = g->lambda_t2;
	g->lambda_t2 = simd_max(neg_max_f, simd_min(simd_add(old_t2, simd_mul(g->eff_mass_t2, simd_neg(vt2))), max_f));
	simd4f dt2 = simd_sub(g->lambda_t2, old_t2);
	va = v3w_sub(va, v3w_scale(g->tangent2, simd_mul(dt2, g->inv_mass_a)));
	vb = v3w_add(vb, v3w_scale(g->tangent2, simd_mul(dt2, g->inv_mass_b)));
	wa = v3w_sub(wa, v3w_scale(g->w_t2_a, dt2));
	wb = v3w_add(wb, v3w_scale(g->w_t2_b, dt2));

	// Torsional friction
	simd4f max_twist = simd_mul(max_f, g->patch_radius);
	simd4f w_rel = v3w_dot(v3w_sub(wb, wa), g->normal);
	simd4f old_tw = g->lambda_twist;
	g->lambda_twist = simd_max(simd_neg(max_twist), simd_min(simd_add(old_tw, simd_mul(g->eff_mass_twist, simd_neg(w_rel))), max_twist));
	simd4f dtw = simd_sub(g->lambda_twist, old_tw);
	wa = v3w_sub(wa, v3w_scale(g->w_tw_a, dtw));
	wb = v3w_add(wb, v3w_scale(g->w_tw_b, dtw));

	simd_scatter_velocities(bodies, g->body_a, va, wa);
	simd_scatter_velocities(bodies, g->body_b, vb, wb);
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
			m->lambda_t1 = lt1[lane];
			m->lambda_t2 = lt2[lane];
			m->lambda_twist = ltw[lane];

			for (int c = 0; c < m->contact_count; c++) {
				float lam[4];
				simd_store(lam, g->normals[c].lambda);
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
		for (int lane = 0; lane < 4; lane++) {
			int mi = g->manifold_idx[lane];
			if (mi < 0) continue;
			SolverManifold* m = &sm[mi];
			for (int c = 0; c < m->contact_count; c++) {
				SolverContact* s = &sc[m->contact_start + c];
				g->normals[c].bias = simd_set_lane(g->normals[c].bias, lane, s->bias);
				g->normals[c].bounce = simd_set_lane(g->normals[c].bounce, lane, s->bounce);
			}
		}
	}
}
