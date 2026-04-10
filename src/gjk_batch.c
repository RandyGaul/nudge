// gjk_batch.c -- Batched GJK distance for multiple pairs simultaneously.
// Processes 4 pairs in parallel using SSE SoA layout.
// Each __m128 holds the same component from 4 different pairs.

typedef struct v3w { __m128 x, y, z; } v3w; // 4-wide v3 (SoA)

static inline v3w v3w_set1(v3 v) { return (v3w){ _mm_set1_ps(v.x), _mm_set1_ps(v.y), _mm_set1_ps(v.z) }; }
static inline v3w v3w_load4(v3 a, v3 b, v3 c, v3 d) { return (v3w){ _mm_set_ps(d.x,c.x,b.x,a.x), _mm_set_ps(d.y,c.y,b.y,a.y), _mm_set_ps(d.z,c.z,b.z,a.z) }; }
static inline v3w v3w_sub(v3w a, v3w b) { return (v3w){ _mm_sub_ps(a.x, b.x), _mm_sub_ps(a.y, b.y), _mm_sub_ps(a.z, b.z) }; }
static inline v3w v3w_add(v3w a, v3w b) { return (v3w){ _mm_add_ps(a.x, b.x), _mm_add_ps(a.y, b.y), _mm_add_ps(a.z, b.z) }; }
static inline v3w v3w_neg(v3w a) { __m128 z = _mm_setzero_ps(); return (v3w){ _mm_sub_ps(z, a.x), _mm_sub_ps(z, a.y), _mm_sub_ps(z, a.z) }; }
static inline __m128 v3w_dot(v3w a, v3w b) { return _mm_add_ps(_mm_add_ps(_mm_mul_ps(a.x, b.x), _mm_mul_ps(a.y, b.y)), _mm_mul_ps(a.z, b.z)); }
static inline __m128 v3w_len2(v3w a) { return v3w_dot(a, a); }
static inline v3w v3w_scale(v3w a, __m128 s) { return (v3w){ _mm_mul_ps(a.x, s), _mm_mul_ps(a.y, s), _mm_mul_ps(a.z, s) }; }
static inline v3w v3w_cross(v3w a, v3w b) {
	return (v3w){
		_mm_sub_ps(_mm_mul_ps(a.y, b.z), _mm_mul_ps(a.z, b.y)),
		_mm_sub_ps(_mm_mul_ps(a.z, b.x), _mm_mul_ps(a.x, b.z)),
		_mm_sub_ps(_mm_mul_ps(a.x, b.y), _mm_mul_ps(a.y, b.x))
	};
}
// Select: mask=all-1s picks a, mask=all-0s picks b.
static inline v3w v3w_sel(v3w a, v3w b, __m128 mask) {
	return (v3w){ _mm_blendv_ps(b.x, a.x, mask), _mm_blendv_ps(b.y, a.y, mask), _mm_blendv_ps(b.z, a.z, mask) };
}

// 4-wide mat3 rotate: R^T * v (column broadcast multiply). Each column is 4-wide.
static inline v3w v3w_mat_rotate_t(v3w c0, v3w c1, v3w c2, v3w v) {
	return (v3w){
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c0.x, v.x), _mm_mul_ps(c1.x, v.y)), _mm_mul_ps(c2.x, v.z)),
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c0.y, v.x), _mm_mul_ps(c1.y, v.y)), _mm_mul_ps(c2.y, v.z)),
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c0.z, v.x), _mm_mul_ps(c1.z, v.y)), _mm_mul_ps(c2.z, v.z))
	};
}
// 4-wide mat3 inverse rotate: R * v (transpose then column broadcast).
static inline v3w v3w_mat_rotate(v3w c0, v3w c1, v3w c2, v3w v) {
	// dot(col_i, v) for each row
	return (v3w){
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c0.x, v.x), _mm_mul_ps(c0.y, v.y)), _mm_mul_ps(c0.z, v.z)),
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c1.x, v.x), _mm_mul_ps(c1.y, v.y)), _mm_mul_ps(c1.z, v.z)),
		_mm_add_ps(_mm_add_ps(_mm_mul_ps(c2.x, v.x), _mm_mul_ps(c2.y, v.y)), _mm_mul_ps(c2.z, v.z))
	};
}

// 4-wide box support: copysign(he, R^T*d) then R*corner + center.
static inline v3w v3w_box_support(v3w center, v3w c0, v3w c1, v3w c2, v3w he, v3w dir) {
	v3w ld = v3w_mat_rotate(c0, c1, c2, dir);
	__m128 sign_mask = _mm_castsi128_ps(_mm_set1_epi32((int)0x80000000));
	v3w lc = {
		_mm_xor_ps(he.x, _mm_and_ps(ld.x, sign_mask)),
		_mm_xor_ps(he.y, _mm_and_ps(ld.y, sign_mask)),
		_mm_xor_ps(he.z, _mm_and_ps(ld.z, sign_mask))
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

	__m128 dsq_prev = _mm_set1_ps(FLT_MAX);
	__m128 done_mask = _mm_setzero_ps(); // lanes that have terminated
	__m128 eps2 = _mm_set1_ps(GJK_CONTAINMENT_EPS2);
	v3w best_p1 = simplex_p1, best_p2 = simplex_p2; // witness points

	for (int iter = 0; iter < 20; iter++) {
		// Closest point on simplex (for count==1, it's just the single point)
		v3w closest = simplex_p;
		__m128 dsq = v3w_len2(closest);

		// Containment check
		__m128 contained = _mm_cmple_ps(dsq, eps2);
		done_mask = _mm_or_ps(done_mask, contained);

		// Monotonic check
		__m128 no_progress = _mm_cmpge_ps(dsq, dsq_prev);
		done_mask = _mm_or_ps(done_mask, no_progress);

		if (_mm_movemask_ps(done_mask) == 0xF) break; // all 4 lanes done

		dsq_prev = _mm_blendv_ps(dsq, dsq_prev, done_mask); // only update active lanes

		// New support
		v3w newA = v3w_box_support(cenA, c0A, c1A, c2A, heA, closest);
		v3w newB = v3w_box_support(cenB, c0B, c1B, c2B, heB, v3w_neg(closest));
		v3w w = v3w_sub(newB, newA);

		// Progress check: dsq - dot(w, closest) <= dsq * eps
		__m128 progress = _mm_sub_ps(dsq, v3w_dot(w, closest));
		__m128 threshold = _mm_mul_ps(v3w_len2(simplex_p), _mm_set1_ps(GJK_PROGRESS_EPS));
		__m128 stalled = _mm_cmple_ps(progress, threshold);
		done_mask = _mm_or_ps(done_mask, stalled);

		// Update simplex for active lanes (simplified: always use single-vertex simplex)
		// This is a simplification — full batch would need per-lane simplex state.
		simplex_p = v3w_sel(simplex_p, w, _mm_andnot_ps(done_mask, _mm_set1_ps(-0.0f))); // hacky: use -0 as all-bits-set
		// Actually, let's just use the new point for active lanes:
		__m128 active = _mm_andnot_ps(done_mask, _mm_castsi128_ps(_mm_set1_epi32(-1)));
		simplex_p = v3w_sel(w, simplex_p, done_mask);
		simplex_p1 = v3w_sel(newA, simplex_p1, done_mask);
		simplex_p2 = v3w_sel(newB, simplex_p2, done_mask);
		best_p1 = v3w_sel(newA, best_p1, done_mask);
		best_p2 = v3w_sel(newB, best_p2, done_mask);
	}

	// Extract distances
	v3w sep = v3w_sub(best_p2, best_p1);
	__m128 dist2 = v3w_len2(sep);
	__m128 dist = _mm_sqrt_ps(dist2);
	_mm_storeu_ps(out_dist, dist);

	// Save cache
	for (int i = 0; i < 4; i++) {
		if (caches[i]) {
			float* sx = (float*)&sep.x, *sy = (float*)&sep.y, *sz = (float*)&sep.z;
			caches[i]->dir = V3(sx[i], sy[i], sz[i]);
		}
	}
}
