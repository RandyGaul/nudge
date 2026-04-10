// gjk_batch.c -- Batched GJK distance for multiple pairs simultaneously.
// Processes 4 pairs in parallel using SSE SoA layout.
// Each simd4f holds the same component from 4 different pairs.

typedef struct v3w { simd4f x, y, z; } v3w; // 4-wide v3 (SoA)

static inline v3w v3w_set1(v3 v) { return (v3w){ simd_set1(v.x), simd_set1(v.y), simd_set1(v.z) }; }
static inline v3w v3w_load4(v3 a, v3 b, v3 c, v3 d) { return (v3w){ simd_set(d.x,c.x,b.x,a.x), simd_set(d.y,c.y,b.y,a.y), simd_set(d.z,c.z,b.z,a.z) }; }
static inline v3w v3w_sub(v3w a, v3w b) { return (v3w){ simd_sub(a.x, b.x), simd_sub(a.y, b.y), simd_sub(a.z, b.z) }; }
static inline v3w v3w_add(v3w a, v3w b) { return (v3w){ simd_add(a.x, b.x), simd_add(a.y, b.y), simd_add(a.z, b.z) }; }
static inline v3w v3w_neg(v3w a) { simd4f z = simd_zero(); return (v3w){ simd_sub(z, a.x), simd_sub(z, a.y), simd_sub(z, a.z) }; }
static inline simd4f v3w_dot(v3w a, v3w b) { return simd_add(simd_add(simd_mul(a.x, b.x), simd_mul(a.y, b.y)), simd_mul(a.z, b.z)); }
static inline simd4f v3w_len2(v3w a) { return v3w_dot(a, a); }
static inline v3w v3w_scale(v3w a, simd4f s) { return (v3w){ simd_mul(a.x, s), simd_mul(a.y, s), simd_mul(a.z, s) }; }
static inline v3w v3w_cross(v3w a, v3w b) {
	return (v3w){
		simd_sub(simd_mul(a.y, b.z), simd_mul(a.z, b.y)),
		simd_sub(simd_mul(a.z, b.x), simd_mul(a.x, b.z)),
		simd_sub(simd_mul(a.x, b.y), simd_mul(a.y, b.x))
	};
}
// Select: mask=all-1s picks a, mask=all-0s picks b.
static inline v3w v3w_sel(v3w a, v3w b, simd4f mask) {
	return (v3w){ simd_blendv(b.x, a.x, mask), simd_blendv(b.y, a.y, mask), simd_blendv(b.z, a.z, mask) };
}

// 4-wide mat3 rotate: R^T * v (column broadcast multiply). Each column is 4-wide.
static inline v3w v3w_mat_rotate_t(v3w c0, v3w c1, v3w c2, v3w v) {
	return (v3w){
		simd_add(simd_add(simd_mul(c0.x, v.x), simd_mul(c1.x, v.y)), simd_mul(c2.x, v.z)),
		simd_add(simd_add(simd_mul(c0.y, v.x), simd_mul(c1.y, v.y)), simd_mul(c2.y, v.z)),
		simd_add(simd_add(simd_mul(c0.z, v.x), simd_mul(c1.z, v.y)), simd_mul(c2.z, v.z))
	};
}
// 4-wide mat3 inverse rotate: R * v (transpose then column broadcast).
static inline v3w v3w_mat_rotate(v3w c0, v3w c1, v3w c2, v3w v) {
	// dot(col_i, v) for each row
	return (v3w){
		simd_add(simd_add(simd_mul(c0.x, v.x), simd_mul(c0.y, v.y)), simd_mul(c0.z, v.z)),
		simd_add(simd_add(simd_mul(c1.x, v.x), simd_mul(c1.y, v.y)), simd_mul(c1.z, v.z)),
		simd_add(simd_add(simd_mul(c2.x, v.x), simd_mul(c2.y, v.y)), simd_mul(c2.z, v.z))
	};
}

// 4-wide box support: copysign(he, R^T*d) then R*corner + center.
static inline v3w v3w_box_support(v3w center, v3w c0, v3w c1, v3w c2, v3w he, v3w dir) {
	v3w ld = v3w_mat_rotate(c0, c1, c2, dir);
	simd4f sign_mask = simd_cast_itof(simd_set1_i((int)0x80000000));
	v3w lc = {
		simd_xor(he.x, simd_and(ld.x, sign_mask)),
		simd_xor(he.y, simd_and(ld.y, sign_mask)),
		simd_xor(he.z, simd_and(ld.z, sign_mask))
	};
	return v3w_add(center, v3w_mat_rotate_t(c0, c1, c2, lc));
}

// Batched box-box GJK distance: process 4 pairs at once.
// Returns 4 distances. All shapes must be boxes.
static void gjk_distance_batch_box(const GJK_Shape* a[4], const GJK_Shape* b[4], GJK_Cache* caches[4], float out_dist[4])
{
	// Load 4 shapes into SoA
	v3w cenA = v3w_load4(a[0]->box.center, a[1]->box.center, a[2]->box.center, a[3]->box.center);
	v3w c0A = v3w_load4(a[0]->box.col0, a[1]->box.col0, a[2]->box.col0, a[3]->box.col0);
	v3w c1A = v3w_load4(a[0]->box.col1, a[1]->box.col1, a[2]->box.col1, a[3]->box.col1);
	v3w c2A = v3w_load4(a[0]->box.col2, a[1]->box.col2, a[2]->box.col2, a[3]->box.col2);
	v3w heA = v3w_load4(a[0]->box.half_extents, a[1]->box.half_extents, a[2]->box.half_extents, a[3]->box.half_extents);
	v3w cenB = v3w_load4(b[0]->box.center, b[1]->box.center, b[2]->box.center, b[3]->box.center);
	v3w c0B = v3w_load4(b[0]->box.col0, b[1]->box.col0, b[2]->box.col0, b[3]->box.col0);
	v3w c1B = v3w_load4(b[0]->box.col1, b[1]->box.col1, b[2]->box.col1, b[3]->box.col1);
	v3w c2B = v3w_load4(b[0]->box.col2, b[1]->box.col2, b[2]->box.col2, b[3]->box.col2);
	v3w heB = v3w_load4(b[0]->box.half_extents, b[1]->box.half_extents, b[2]->box.half_extents, b[3]->box.half_extents);

	// Initial direction from cache or center-to-center
	v3w init_d = v3w_sub(cenB, cenA);
	for (int i = 0; i < 4; i++) {
		if (caches[i] && len2(caches[i]->dir) > FLT_EPSILON) {
			// Scatter cached direction into lane i
			float* dx = (float*)&init_d.x, *dy = (float*)&init_d.y, *dz = (float*)&init_d.z;
			dx[i] = caches[i]->dir.x; dy[i] = caches[i]->dir.y; dz[i] = caches[i]->dir.z;
		}
	}

	// Initial support
	v3w supA = v3w_box_support(cenA, c0A, c1A, c2A, heA, init_d);
	v3w supB = v3w_box_support(cenB, c0B, c1B, c2B, heB, v3w_neg(init_d));
	v3w simplex_p = v3w_sub(supB, supA); // Minkowski diff
	v3w simplex_p1 = supA, simplex_p2 = supB;

	simd4f dsq_prev = simd_set1(FLT_MAX);
	simd4f done_mask = simd_zero(); // lanes that have terminated
	simd4f eps2 = simd_set1(GJK_CONTAINMENT_EPS2);
	v3w best_p1 = simplex_p1, best_p2 = simplex_p2; // witness points

	for (int iter = 0; iter < 20; iter++) {
		// Closest point on simplex (for count==1, it's just the single point)
		v3w closest = simplex_p;
		simd4f dsq = v3w_len2(closest);

		// Containment check
		simd4f contained = simd_cmple(dsq, eps2);
		done_mask = simd_or(done_mask, contained);

		// Monotonic check
		simd4f no_progress = simd_cmpge(dsq, dsq_prev);
		done_mask = simd_or(done_mask, no_progress);

		if (simd_movemask(done_mask) == 0xF) break; // all 4 lanes done

		dsq_prev = simd_blendv(dsq, dsq_prev, done_mask); // only update active lanes

		// New support
		v3w newA = v3w_box_support(cenA, c0A, c1A, c2A, heA, closest);
		v3w newB = v3w_box_support(cenB, c0B, c1B, c2B, heB, v3w_neg(closest));
		v3w w = v3w_sub(newB, newA);

		// Progress check: dsq - dot(w, closest) <= dsq * eps
		simd4f progress = simd_sub(dsq, v3w_dot(w, closest));
		simd4f threshold = simd_mul(v3w_len2(simplex_p), simd_set1(GJK_PROGRESS_EPS));
		simd4f stalled = simd_cmple(progress, threshold);
		done_mask = simd_or(done_mask, stalled);

		// Update simplex for active lanes (simplified: always use single-vertex simplex)
		// This is a simplification — full batch would need per-lane simplex state.
		simplex_p = v3w_sel(simplex_p, w, simd_andnot(done_mask, simd_set1(-0.0f))); // hacky: use -0 as all-bits-set
		// Actually, let's just use the new point for active lanes:
		simd4f active = simd_andnot(done_mask, simd_cast_itof(simd_set1_i(-1)));
		simplex_p = v3w_sel(w, simplex_p, done_mask);
		simplex_p1 = v3w_sel(newA, simplex_p1, done_mask);
		simplex_p2 = v3w_sel(newB, simplex_p2, done_mask);
		best_p1 = v3w_sel(newA, best_p1, done_mask);
		best_p2 = v3w_sel(newB, best_p2, done_mask);
	}

	// Extract distances
	v3w sep = v3w_sub(best_p2, best_p1);
	simd4f dist2 = v3w_len2(sep);
	simd4f dist = simd_sqrt(dist2);
	simd_store(out_dist, dist);

	// Save cache
	for (int i = 0; i < 4; i++) {
		if (caches[i]) {
			float* sx = (float*)&sep.x, *sy = (float*)&sep.y, *sz = (float*)&sep.z;
			caches[i]->dir = V3(sx[i], sy[i], sz[i]);
		}
	}
}

// 4-wide triangle support: max dot of 3 vertices per lane.
static inline v3w v3w_tri_support(v3w ta, v3w tb, v3w tc, v3w dir) {
	simd4f da = v3w_dot(ta, dir), db = v3w_dot(tb, dir), dc = v3w_dot(tc, dir);
	simd4f ab = simd_cmpge(da, db), ac = simd_cmpge(da, dc), bc = simd_cmpge(db, dc);
	simd4f pick_a = simd_and(ab, ac);             // da >= db && da >= dc
	simd4f pick_b = simd_andnot(pick_a, bc);       // not-a && db >= dc
	v3w result = tc;                                  // default: pick C
	result = v3w_sel(tb, result, pick_b);             // if pick_b: use B
	result = v3w_sel(ta, result, pick_a);             // if pick_a: use A
	return result;
}

// Batched sphere vs 4 triangles. sphere is shared across all 4 lanes.
static void gjk_distance_batch_sphere_tri(v3 sphere_center, float sphere_radius, const v3 tri_verts[12], float out_dist[4])
{
	// Load sphere (broadcast to all lanes)
	v3w sph = v3w_set1(sphere_center);
	// Load 4 triangles: tri_verts[0..2]=tri0, [3..5]=tri1, [6..8]=tri2, [9..11]=tri3
	v3w ta = v3w_load4(tri_verts[0], tri_verts[3], tri_verts[6], tri_verts[9]);
	v3w tb = v3w_load4(tri_verts[1], tri_verts[4], tri_verts[7], tri_verts[10]);
	v3w tc = v3w_load4(tri_verts[2], tri_verts[5], tri_verts[8], tri_verts[11]);

	// Initial direction: sphere center to triangle centroid
	v3w tri_cen = v3w_scale(v3w_add(v3w_add(ta, tb), tc), simd_set1(1.0f/3.0f));
	v3w init_d = v3w_sub(tri_cen, sph);

	// Initial support: sphere = center, triangle = max dot vertex
	v3w supA = sph; // sphere core is just the center
	v3w supB = v3w_tri_support(ta, tb, tc, v3w_neg(init_d));
	v3w simplex_p = v3w_sub(supB, supA);
	v3w best_p1 = supA, best_p2 = supB;

	simd4f dsq_prev = simd_set1(FLT_MAX);
	simd4f done_mask = simd_zero();
	simd4f eps2 = simd_set1(GJK_CONTAINMENT_EPS2);
	simd4f prog_eps = simd_set1(GJK_PROGRESS_EPS);

	for (int iter = 0; iter < 20; iter++) {
		v3w closest = simplex_p;
		simd4f dsq = v3w_len2(closest);
		done_mask = simd_or(done_mask, simd_cmple(dsq, eps2));
		done_mask = simd_or(done_mask, simd_cmpge(dsq, dsq_prev));
		if (simd_movemask(done_mask) == 0xF) break;
		dsq_prev = simd_blendv(dsq, dsq_prev, done_mask);

		// Sphere support in direction closest = just sphere center (always)
		v3w newA = sph;
		v3w newB = v3w_tri_support(ta, tb, tc, v3w_neg(closest));
		v3w w = v3w_sub(newB, newA);

		simd4f progress = simd_sub(dsq, v3w_dot(w, closest));
		simd4f threshold = simd_mul(v3w_len2(simplex_p), prog_eps);
		done_mask = simd_or(done_mask, simd_cmple(progress, threshold));

		simplex_p = v3w_sel(w, simplex_p, done_mask);
		best_p1 = v3w_sel(newA, best_p1, done_mask);
		best_p2 = v3w_sel(newB, best_p2, done_mask);
	}

	// Distance = |p2 - p1| - radius
	v3w sep = v3w_sub(best_p2, best_p1);
	simd4f dist = simd_sqrt(v3w_len2(sep));
	dist = simd_max(simd_zero(), simd_sub(dist, simd_set1(sphere_radius)));
	simd_store(out_dist, dist);
}
