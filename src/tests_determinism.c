// See LICENSE for licensing info.
// tests_determinism.c -- cross-platform FP determinism hash check.
//
// Runs a fixed scene for a fixed number of steps with a fixed worker count,
// then hashes every body's transform + velocity. The expected hash is pinned
// per-platform (once) and checked in CI so any drift across compilers or ISAs
// is caught as a test failure. Uses deterministic transcendentals (nudge_*)
// and hardware sqrt only; build disables FMA contraction via -ffp-contract=off.

// FNV-1a 64-bit. Portable, no libm, one pass.
static uint64_t det_fnv1a_update(uint64_t h, const void* data, size_t n)
{
	const uint8_t* b = (const uint8_t*)data;
	for (size_t i = 0; i < n; i++) {
		h ^= b[i];
		h *= 0x100000001b3ull;
	}
	return h;
}

static uint64_t det_hash_world(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	uint64_t h = 0xcbf29ce484222325ull;
	int nb = asize(w->body_hot);
	for (int i = 0; i < nb; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyState* s = &w->body_state[i];
		BodyHot*   b = &w->body_hot[i];
		h = det_fnv1a_update(h, &s->position, sizeof(s->position));
		h = det_fnv1a_update(h, &s->rotation, sizeof(s->rotation));
		h = det_fnv1a_update(h, &b->velocity, sizeof(b->velocity));
		h = det_fnv1a_update(h, &b->angular_velocity, sizeof(b->angular_velocity));
	}
	return h;
}

// Build a canonical scene. Mixed shapes so we exercise several narrowphase
// pairs; small stack to exercise contact manifold accumulation; a tilt on
// the top boxes to give the PGS something to chew on.
static World det_build_scene()
{
	WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH };
	World w = create_world(wp);

	// Floor
	Body floor_b = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
	body_add_shape(w, floor_b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(20, 1, 20) });

	// 6-box stack, each tilted slightly so the rest isn't trivial.
	for (int i = 0; i < 6; i++) {
		float angle = 0.02f * (float)i;
		float s = nudge_sinf(angle * 0.5f);
		float c = nudge_cosf(angle * 0.5f);
		quat q = { 0, 0, s, c };
		Body b = create_body(w, (BodyParams){
			.position = V3(0, 0.5f + (float)i * 1.02f, 0),
			.rotation = q,
			.mass = 1.0f,
			.friction = 0.5f,
		});
		body_add_shape(w, b, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(0.5f, 0.5f, 0.5f) });
	}

	// A cylinder rolling in with lateral velocity so the cylinder pair and
	// friction code paths all contribute to the hash.
	Body cyl = create_body(w, (BodyParams){
		.position = V3(-4, 2, 0),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.friction = 0.5f,
	});
	body_add_shape(w, cyl, (ShapeParams){ .type = SHAPE_CYLINDER, .cylinder = { .half_height = 0.5f, .radius = 0.4f } });
	body_set_velocity(w, cyl, V3(3, 0, 0));

	// A capsule lofted over the stack, gives us sphere-capsule + capsule-box hits.
	Body cap = create_body(w, (BodyParams){
		.position = V3(2, 5, 0.3f),
		.rotation = quat_identity(),
		.mass = 1.0f,
		.friction = 0.3f,
	});
	body_add_shape(w, cap, (ShapeParams){ .type = SHAPE_CAPSULE, .capsule = { .half_height = 0.4f, .radius = 0.3f } });
	body_set_angular_velocity(w, cap, V3(0.5f, 0, 0));

	return w;
}

// Run the scene for N steps and return the hash.
static uint64_t det_run(int steps, int threads)
{
	World w = det_build_scene();
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->thread_count = threads;
	float dt = 1.0f / 60.0f;
	for (int i = 0; i < steps; i++) {
		world_step(w, dt);
	}
	uint64_t h = det_hash_world(w);
	destroy_world(w);
	return h;
}

// Build the scene, then print a hash at step 0, 1, 2, 10, 60, 240 so CI logs
// show where cross-arch divergence starts to accumulate. Threads-free so the
// trace is deterministic independent of the pool.
// Print the raw bits of a known `a*b+c` pattern. volatile locals prevent
// the compiler from constant-folding -- whatever the backend emits for a
// runtime mul-add is what we see. If macOS/WASM show different bits than
// x86, the backend is still fusing to FMA despite our flags.
static void det_dump_fma_probe()
{
	volatile float a = 1.0000001f;
	volatile float b = 3.0000003f;
	volatile float c = 5.0000005f;
	float r = a * b + c;
	uint32_t bits;
	memcpy(&bits, &r, 4);
	printf("det fma probe: a*b+c bits = 0x%08x (%.9g)\n", bits, (double)r);
}

// Run the scene single-threaded, printing the world hash at a handful of
// checkpoints so CI logs show where divergence starts accumulating.
static void det_trace()
{
	World w = det_build_scene();
	WorldInternal* wi = (WorldInternal*)w.id;
	wi->thread_count = 1;
	float dt = 1.0f / 60.0f;
	int checkpoints[] = { 0, 1, 2, 10, 60, 240 };
	int ci = 0;
	int next = checkpoints[ci];
	printf("det trace:\n");
	for (int i = 0; i <= 240; i++) {
		if (i == next) {
			printf("  step %3d: 0x%016llx\n", i, (unsigned long long)det_hash_world(w));
			ci++;
			if (ci < (int)(sizeof(checkpoints) / sizeof(checkpoints[0]))) next = checkpoints[ci];
			else next = -1;
		}
		if (i < 240) world_step(w, dt);
	}
	destroy_world(w);
}

// Per-arch pinned hashes. Within an arch every compiler (MSVC / GCC / Clang)
// and every SIMD backend (SSE / NEON / WASM SIMD128 / scalar) produces the
// exact same hash -- that's the cross-compiler guarantee this test enforces.
//
// Scene setup IS cross-arch bit-identical now (`step 0` hash matches on
// x86, aarch64, and wasm32) thanks to #pragma STDC FP_CONTRACT OFF and the
// -fno-{associative,reciprocal,unsafe}-math flag stack. Simulation steps
// still diverge cross-arch because the narrowphase makes discrete contact
// decisions on sub-ULP float comparisons, and a 1-ULP difference in SAT
// depth flips whether a contact fires. Fixing that would need hysteresis
// in the depth tests (e.g. `depth > -LINEAR_SLOP` instead of `depth > 0`),
// not just more FP flags -- tracked as a follow-up.
#if defined(__wasm__)
#define DET_EXPECTED_HASH 0xebbb50883882fe05ULL
#elif defined(__aarch64__) || defined(_M_ARM64)
#define DET_EXPECTED_HASH 0xc293410039033e05ULL
#else
#define DET_EXPECTED_HASH 0x86c44c829ce09c07ULL
#endif

// Runs the canonical scene twice (single-threaded and N-threaded) and checks:
//   1. Both runs produce the same hash -- threading does not affect output.
//      (Relies on disjoint-body graph coloring in the PGS solver.)
//   2. The hash matches DET_EXPECTED_HASH if pinned -- cross-platform check.
// Any failure returns nonzero so CI flags it.
static int run_determinism_test(int threads_hi)
{
	const int steps = 240;
	det_dump_fma_probe();
	det_trace();
	uint64_t h1 = det_run(steps, 1);
	printf("det hash (threads=1): 0x%016llx\n", (unsigned long long)h1);
	int fail = 0;
	// Threading determinism check only when threads_hi > 1 -- otherwise both
	// calls would be the same code path and the comparison is meaningless.
	// (On WASM we pass --threads 1 because pthreads aren't enabled there.)
	if (threads_hi > 1) {
		uint64_t hN = det_run(steps, threads_hi);
		printf("det hash (threads=%d): 0x%016llx\n", threads_hi, (unsigned long long)hN);
		if (h1 != hN) {
			fprintf(stderr, "FAIL: threading determinism -- threads=1 and threads=%d produced different hashes.\n", threads_hi);
			fail = 1;
		}
	}
	if (DET_EXPECTED_HASH != 0 && h1 != DET_EXPECTED_HASH) {
		fprintf(stderr, "FAIL: cross-platform determinism -- expected 0x%016llx, got 0x%016llx.\n",
			(unsigned long long)DET_EXPECTED_HASH, (unsigned long long)h1);
		fail = 1;
	}
	return fail;
}
