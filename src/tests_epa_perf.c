// See LICENSE for licensing info.
// tests_epa_perf.c -- per-query microbench for EPA vs SAT on canonical
// deep-contact configurations. Reuses qpc_now from tests_gjk_perf.c.

#define EPA_PERF_N 20000

static double epa_bench_box_box_deep()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	ConvexHull ha = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	ConvexHull hb = { box, V3(1.0f, 0, 0), id, V3(1, 1, 1) };
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (epa_hull_hull(ha, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

static double sat_bench_box_box_deep()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	ConvexHull ha = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	ConvexHull hb = { box, V3(1.0f, 0, 0), id, V3(1, 1, 1) };
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (collide_hull_hull(ha, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

static double epa_bench_sphere_hull_deep()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	ConvexHull hb = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	Sphere sa = { V3(0.5f, 0, 0), 1.0f }; // deep overlap
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (epa_sphere_hull(sa, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

static double sat_bench_sphere_hull_deep()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	ConvexHull hb = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	Sphere sa = { V3(0.5f, 0, 0), 1.0f };
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (collide_sphere_hull(sa, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

static double epa_bench_hull_hull_rotated()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	quat rot45y = { 0, 0.3827f, 0, 0.9239f };
	ConvexHull ha = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	ConvexHull hb = { box, V3(1.5f, 0, 0), rot45y, V3(1, 1, 1) };
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (epa_hull_hull(ha, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

static double sat_bench_hull_hull_rotated()
{
	const Hull* box = hull_unit_box();
	quat id = quat_identity();
	quat rot45y = { 0, 0.3827f, 0, 0.9239f };
	ConvexHull ha = { box, V3(0, 0, 0), id, V3(1, 1, 1) };
	ConvexHull hb = { box, V3(1.5f, 0, 0), rot45y, V3(1, 1, 1) };
	double t0 = qpc_now();
	int hits = 0;
	for (int i = 0; i < EPA_PERF_N; i++) {
		Manifold m;
		if (collide_hull_hull(ha, hb, &m)) hits++;
	}
	double t1 = qpc_now();
	(void)hits;
	return (t1 - t0) * 1e9 / EPA_PERF_N;
}

// World-level step bench: run a small box stack with each backend, measure
// ms/frame. Captures the incremental-manifold contact-growth effect on
// world_step rather than per-query cost.
static double bench_world_step_backend(int backend, int frames)
{
	WorldParams wp = { .gravity = V3(0, -9.81f, 0), .narrowphase_backend = (NarrowphaseBackend)backend };
	World w = create_world(wp);
	Body floor = create_body(w, (BodyParams){ .position = V3(0, -0.5f, 0), .rotation = quat_identity(), .mass = 0, .friction = 0.5f });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	for (int i = 0; i < 5; i++) {
		Body b = create_body(w, (BodyParams){ .position = V3(0, 1.0f + (float)i * 1.05f, 0), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.5f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	}
	// Settle a few frames before measuring.
	for (int i = 0; i < 20; i++) world_step(w, 1.0f / 60.0f);
	double t0 = qpc_now();
	for (int i = 0; i < frames; i++) world_step(w, 1.0f / 60.0f);
	double t1 = qpc_now();
	destroy_world(w);
	return (t1 - t0) * 1000.0 / frames; // ms/frame
}

// Forward decl: defined below.
static void bench_epa_scenes();

static void bench_epa_vs_sat()
{
	qpc_init();
	printf("\n=== EPA vs SAT per-query bench (%d iters per config) ===\n", EPA_PERF_N);
	printf("  %-30s  %10s  %10s  %6s\n", "config", "SAT ns/op", "EPA ns/op", "ratio");
	printf("  %-30s  %10s  %10s  %6s\n", "------", "---------", "---------", "-----");

	double sat_bb = sat_bench_box_box_deep();
	double epa_bb = epa_bench_box_box_deep();
	printf("  %-30s  %10.1f  %10.1f  %6.2fx\n", "box-box deep (overlap 1)", sat_bb, epa_bb, epa_bb / sat_bb);

	double sat_sh = sat_bench_sphere_hull_deep();
	double epa_sh = epa_bench_sphere_hull_deep();
	printf("  %-30s  %10.1f  %10.1f  %6.2fx\n", "sphere-hull deep", sat_sh, epa_sh, epa_sh / sat_sh);

	double sat_hh = sat_bench_hull_hull_rotated();
	double epa_hh = epa_bench_hull_hull_rotated();
	printf("  %-30s  %10.1f  %10.1f  %6.2fx\n", "hull-hull rotated deep", sat_hh, epa_hh, epa_hh / sat_hh);

	printf("\n=== World-step bench: 5-box stack on floor, 120 frames ===\n");
	double sat_ws = bench_world_step_backend(NARROWPHASE_SAT, 120);
	double epa_ws = bench_world_step_backend(NARROWPHASE_GJK_EPA, 120);
	printf("  SAT backend:     %7.3f ms/frame\n", sat_ws);
	printf("  GJK+EPA backend: %7.3f ms/frame   (%.2fx SAT)\n", epa_ws, epa_ws / sat_ws);

	printf("\n");
}

// -----------------------------------------------------------------------------
// Quality-metric scene bench (Phase D3/D4/D5). For each scene, step both
// backends and record:
//   - settle time: frames until max linear velocity drops below 0.01 m/s.
//   - peak penetration: max contact.penetration across all frames.
//   - avg contact count per pair per frame (last 20% of run).
//   - per-phase timing: narrowphase, solve, total (ms/frame).
//   - iteration-count histogram (EPA backend only).
//   - contact-count timeline for the first EPA pair.

typedef struct EpaSceneBuilder
{
	const char* name;
	void (*build)(World w);
	int frames;
} EpaSceneBuilder;

static void scene_3_box_stack(World w)
{
	Body floor = create_body(w, (BodyParams){ .position = V3(0,-0.5f,0), .rotation = quat_identity(), .mass = 0, .friction = 0.5f });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	for (int i = 0; i < 3; i++) {
		Body b = create_body(w, (BodyParams){ .position = V3(0, 1.0f + (float)i * 1.05f, 0), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.5f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	}
}

static void scene_10_box_stack(World w)
{
	Body floor = create_body(w, (BodyParams){ .position = V3(0,-0.5f,0), .rotation = quat_identity(), .mass = 0, .friction = 0.5f });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	for (int i = 0; i < 10; i++) {
		Body b = create_body(w, (BodyParams){ .position = V3(0, 1.0f + (float)i * 1.05f, 0), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.5f });
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	}
}

static void scene_5x5_pile(World w)
{
	Body floor = create_body(w, (BodyParams){ .position = V3(0,-0.5f,0), .rotation = quat_identity(), .mass = 0, .friction = 0.5f });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 0.5f, 10) });
	for (int z = 0; z < 5; z++) {
		for (int x = 0; x < 5; x++) {
			Body b = create_body(w, (BodyParams){ .position = V3(-2.0f + (float)x * 1.05f, 1.0f, -2.0f + (float)z * 1.05f), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.5f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
		}
	}
}

static void scene_small_pyramid(World w)
{
	// 10 base: 4 rows (10, 6, 3, 1)? Use simpler 4-3-2-1 triangular pyramid.
	Body floor = create_body(w, (BodyParams){ .position = V3(0,-0.5f,0), .rotation = quat_identity(), .mass = 0, .friction = 0.5f });
	body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(20, 0.5f, 20) });
	int rows[] = { 10, 6, 3, 1 };
	float row_y = 0.5f;
	for (int r = 0; r < 4; r++) {
		int n = rows[r];
		for (int i = 0; i < n; i++) {
			float x = ((float)i - 0.5f * ((float)n - 1.0f)) * 1.05f;
			Body b = create_body(w, (BodyParams){ .position = V3(x, row_y, 0), .rotation = quat_identity(), .mass = 1.0f, .friction = 0.5f });
			body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
		}
		row_y += 1.05f;
	}
}

typedef struct EpaSceneResult
{
	int settle_frame;       // -1 if never settled
	float peak_pen;
	float avg_contacts_per_pair; // averaged over last 20% of run
	double ms_total;        // averaged ms/frame
	double ms_np;           // averaged ms/frame narrowphase
	double ms_solve;        // averaged ms/frame solver (pre+pgs+pos_correct)
	int iter_hist[6];       // buckets: [1-4, 5-8, 9-16, 17-32, 33-48, cap]
	int total_queries;
	int total_warm;
	int total_iter_cap;
	long long total_iters;  // cumulative iterations
	int first_pair_ct[8];   // contact counts at frames [0,1,2,3,5,10,settled,final]
	int first_pair_frames[8];
} EpaSceneResult;

static EpaSceneResult run_scene(const EpaSceneBuilder* s, int backend)
{
	WorldParams wp = { .gravity = V3(0, -9.81f, 0), .narrowphase_backend = (NarrowphaseBackend)backend };
	World w = create_world(wp);
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->sleep_enabled = 0;
	s->build(w);

	EpaSceneResult r = { .settle_frame = -1, .peak_pen = 0.0f, .avg_contacts_per_pair = 0.0f, .ms_total = 0, .ms_np = 0, .ms_solve = 0 };
	for (int i = 0; i < 6; i++) r.iter_hist[i] = 0;
	for (int i = 0; i < 8; i++) { r.first_pair_ct[i] = 0; r.first_pair_frames[i] = -1; }

	int sample_start = (int)(0.8f * s->frames);
	int sample_count = s->frames - sample_start;
	double sample_contact_sum = 0;
	double sample_pair_sum = 0;

	double t0_total = qpc_now();
	double acc_np = 0, acc_solve = 0;
	int snapshot_idx = 0;
	int snapshot_frames[] = { 0, 1, 2, 3, 5, 10, -1, -1 }; // last two filled at end
	for (int frame = 0; frame < s->frames; frame++) {
		world_step(w, 1.0f / 60.0f);

		// Peak penetration from debug contacts.
		const Contact* contacts; int nc = world_get_contacts(w, &contacts);
		for (int i = 0; i < nc; i++) {
			if (contacts[i].penetration > r.peak_pen) r.peak_pen = contacts[i].penetration;
		}

		// Max linear velocity.
		float max_v2 = 0.0f;
		int dyn_count = 0;
		for (int bi = 0; bi < asize(wi->body_hot); bi++) {
			if (!split_alive(wi->body_gen, bi)) continue;
			if (body_inv_mass(wi, bi) <= 0.0f) continue;
			v3 v = wi->body_hot[bi].velocity;
			float v2 = dot(v, v);
			if (v2 > max_v2) max_v2 = v2;
			dyn_count++;
		}
		if (r.settle_frame < 0 && max_v2 < 0.0001f && frame > 10) {
			r.settle_frame = frame;
		}

		// Timing accumulators.
		acc_np    += wi->perf.broadphase + wi->perf.pre_solve;
		acc_solve += wi->perf.pgs_solve + wi->perf.position_correct;

		// EPA stats (backend-specific). Sample every frame regardless of backend;
		// counters are zero under SAT so the arithmetic is safe.
		WorldEpaStats es = world_get_epa_stats(w);
		if (backend == NARROWPHASE_GJK_EPA && es.queries > 0) {
			r.total_queries += es.queries;
			r.total_warm    += es.warm_reseeds;
			r.total_iter_cap += es.iter_cap_hits;
			r.total_iters    += (long long)es.total_iters;
			// Histogram: bucket by avg iters for the frame.
			float avg_iter = (float)es.total_iters / (float)es.queries;
			int bucket;
			if      (avg_iter <=  4) bucket = 0;
			else if (avg_iter <=  8) bucket = 1;
			else if (avg_iter <= 16) bucket = 2;
			else if (avg_iter <= 32) bucket = 3;
			else if (avg_iter <= 48) bucket = 4;
			else                     bucket = 5;
			r.iter_hist[bucket] += es.queries;
		}

		// Last 20% sampling for manifold quality.
		if (frame >= sample_start) {
			int pair_count = (int)map_size(wi->epa_cache);
			int contact_count = 0;
			for (int i = 0; i < pair_count; i++) contact_count += wi->epa_cache[i].count;
			if (backend == NARROWPHASE_GJK_EPA && pair_count > 0) {
				sample_contact_sum += (double)contact_count;
				sample_pair_sum    += (double)pair_count;
			}
		}

		// Contact-count timeline for the first EPA pair.
		if (backend == NARROWPHASE_GJK_EPA && map_size(wi->epa_cache) > 0 && snapshot_idx < 6) {
			if (frame == snapshot_frames[snapshot_idx]) {
				r.first_pair_ct[snapshot_idx] = wi->epa_cache[0].count;
				r.first_pair_frames[snapshot_idx] = frame;
				snapshot_idx++;
			}
		}
	}
	double t1_total = qpc_now();

	// Final snapshots: settled (or -1) and final.
	if (backend == NARROWPHASE_GJK_EPA && map_size(wi->epa_cache) > 0) {
		if (r.settle_frame >= 0) {
			r.first_pair_ct[6] = wi->epa_cache[0].count;
			r.first_pair_frames[6] = r.settle_frame;
		}
		r.first_pair_ct[7] = wi->epa_cache[0].count;
		r.first_pair_frames[7] = s->frames - 1;
	}

	r.ms_total = (t1_total - t0_total) * 1000.0 / (double)s->frames;
	r.ms_np    = acc_np * 1000.0 / (double)s->frames;
	r.ms_solve = acc_solve * 1000.0 / (double)s->frames;

	if (sample_pair_sum > 0) {
		r.avg_contacts_per_pair = (float)(sample_contact_sum / sample_pair_sum);
	}

	destroy_world(w);
	return r;
}

static void bench_epa_scenes()
{
	qpc_init();
	EpaSceneBuilder scenes[] = {
		{ "3-box stack",   scene_3_box_stack,  300 },
		{ "10-box stack",  scene_10_box_stack, 600 },
		{ "5x5 pile",      scene_5x5_pile,     400 },
		{ "small pyramid", scene_small_pyramid, 500 },
	};
	int n = (int)(sizeof(scenes) / sizeof(scenes[0]));

	printf("\n=== EPA scene quality bench (SAT vs GJK+EPA) ===\n");
	printf("  %-14s  %-9s  %8s  %8s  %9s  %9s  %9s  %9s\n",
		"scene", "backend", "settle", "peakpen", "avgCt/pr", "ms/np", "ms/solve", "ms/total");
	printf("  %-14s  %-9s  %8s  %8s  %9s  %9s  %9s  %9s\n",
		"-----", "-------", "------", "-------", "--------", "-----", "--------", "--------");
	for (int i = 0; i < n; i++) {
		EpaSceneResult sat = run_scene(&scenes[i], NARROWPHASE_SAT);
		EpaSceneResult epa = run_scene(&scenes[i], NARROWPHASE_GJK_EPA);
		printf("  %-14s  %-9s  %8d  %8.4f  %9.2f  %9.3f  %9.3f  %9.3f\n",
			scenes[i].name, "SAT",     sat.settle_frame, sat.peak_pen, sat.avg_contacts_per_pair,
			sat.ms_np, sat.ms_solve, sat.ms_total);
		printf("  %-14s  %-9s  %8d  %8.4f  %9.2f  %9.3f  %9.3f  %9.3f\n",
			"", "GJK+EPA", epa.settle_frame, epa.peak_pen, epa.avg_contacts_per_pair,
			epa.ms_np, epa.ms_solve, epa.ms_total);

		// EPA-only detail rows.
		float avg_iters = epa.total_queries > 0 ? (float)epa.total_iters / (float)epa.total_queries : 0.0f;
		printf("    EPA queries=%d  warm=%d (%.1f%%)  iter-cap=%d  avg-iters=%.2f\n",
			epa.total_queries, epa.total_warm,
			epa.total_queries > 0 ? 100.0f * (float)epa.total_warm / (float)epa.total_queries : 0.0f,
			epa.total_iter_cap, avg_iters);
		printf("    EPA iter histogram (queries per avg-iter bucket):\n");
		const char* bname[] = { "[1-4]", "[5-8]", "[9-16]", "[17-32]", "[33-48]", "[cap]" };
		for (int b = 0; b < 6; b++) {
			printf("      %-8s %d\n", bname[b], epa.iter_hist[b]);
		}
		printf("    EPA first-pair contact count timeline:\n");
		const char* tag[] = { "f0", "f1", "f2", "f3", "f5", "f10", "settled", "final" };
		for (int k = 0; k < 8; k++) {
			if (epa.first_pair_frames[k] >= 0) {
				printf("      %-8s (frame %d): %d contacts\n", tag[k], epa.first_pair_frames[k], epa.first_pair_ct[k]);
			}
		}
	}
	printf("\n");
}
