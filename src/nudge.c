// See LICENSE for licensing info.
// nudge.c -- physics world implementation

#include "perf.h"
#include "nudge_internal.h"

// Forward declarations for the thread-pool dispatch primitives. The full
// definitions live later in this TU; ldl per-island dispatch paths in
// solver_ldl.c (included below) call pool_dispatch, so the typedef + decl
// need to be visible before that include.
typedef void (*WorkFn)(void* ctx, int start, int count);
static void pool_dispatch(WorkFn fn, void* ctx, int total_items, int block_size, int n_workers);

#include "gjk.c"
#include "gjk_batch.c"
#include "quickhull.c"
#include "bvh.c"
#include "collision.c"
#include "broadphase.c"
#include "epa.c"
#include "inertia.c"
#include "solver_pgs.c"
#include "joints.c"
#include "solver_ldl.c"
#include "islands.c"

// -----------------------------------------------------------------------------
// World.

World create_world(WorldParams params)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	w->gravity = params.gravity;
	w->broadphase_type = params.broadphase;
	w->narrowphase_backend = (int)params.narrowphase_backend;
	w->solver_type = params.solver_type;
	w->sleep_enabled = 1;
	w->sat_hint_enabled = 1;
	w->sat_hillclimb_enabled = 1;
	w->warm_start_enabled = 1;
	w->incremental_np_enabled = 1;
	w->velocity_iters = params.velocity_iters > 0 ? params.velocity_iters : SOLVER_VELOCITY_ITERS;
	w->position_iters = params.position_iters > 0 ? params.position_iters : SOLVER_POSITION_ITERS;
	w->contact_hertz = params.contact_hertz > 0.0f ? params.contact_hertz : 60.0f;
	w->contact_damping_ratio = params.contact_damping_ratio > 0.0f ? params.contact_damping_ratio : 3.0f;
	w->max_push_velocity = params.max_push_velocity > 0.0f ? params.max_push_velocity : 3.0f;
	w->sub_steps = params.sub_steps > 0 ? params.sub_steps : 4;
	w->ldl_correction_iter = -2; // -2 = auto: velocity_iters/2 (mid-loop, PGS can recover after LDL)
	w->bvh_static = CK_ALLOC(sizeof(BVH_Tree));
	w->bvh_dynamic = CK_ALLOC(sizeof(BVH_Tree));
	w->bvh_sleeping = CK_ALLOC(sizeof(BVH_Tree));
	bvh_init(w->bvh_static);
	bvh_init(w->bvh_dynamic);
	bvh_init(w->bvh_sleeping);
	return (World){ (uint64_t)w };
}

void destroy_world(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	for (int i = 0; i < asize(w->body_cold); i++) {
		afree(w->body_cold[i].shapes);
	}
	afree(w->debug_contacts);
	afree(w->body_state);
	map_free(w->warm_cache);
	map_free(w->epa_cache);
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	bvh_free(w->bvh_sleeping); CK_FREE(w->bvh_sleeping);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	for (int i = 0; i < asize(w->islands); i++) ldl_cache_free(&w->islands[i].ldl);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	map_free(w->joint_pairs);
	for (int i = 0; i < asize(w->worker_arenas); i++) arena_free(&w->worker_arenas[i]);
	afree(w->worker_arenas);
	CK_FREE(w);
}

// Integrate velocities + precompute world inertia in one fused pass.
// Halves body_hot cache traversals vs separate loops (both touch same data).
static void integrate_velocities_and_inertia(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		BodyState* s = &w->body_state[i];
		body_compute_inv_inertia_world(h, s);
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
		if (s->linear_damping > 0.0f)
			h->velocity = scale(h->velocity, 1.0f / (1.0f + s->linear_damping * dt));
		if (s->angular_damping > 0.0f)
			h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + s->angular_damping * dt));
	}
}

// Integrate positions and rotations for a sub-step.
static void integrate_positions(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		BodyState* s = &w->body_state[i];

		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL)
			h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL)
			h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));

		s->position = add(s->position, scale(h->velocity, dt));

		// Skip gyroscopic for uniform inertia (cubes/spheres: cross(w,I*w)=0) or negligible angular vel.
		if (av2 > 0.01f && !(s->inv_inertia_local.x == s->inv_inertia_local.y && s->inv_inertia_local.y == s->inv_inertia_local.z))
			h->angular_velocity = solve_gyroscopic(s->rotation, s->inv_inertia_local, h->angular_velocity, dt);

		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, s->rotation);
		s->rotation.x += 0.5f * dt * dq.x;
		s->rotation.y += 0.5f * dt * dq.y;
		s->rotation.z += 0.5f * dt * dq.z;
		s->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(s->rotation.x*s->rotation.x + s->rotation.y*s->rotation.y
			+ s->rotation.z*s->rotation.z + s->rotation.w*s->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		s->rotation.x *= inv_ql; s->rotation.y *= inv_ql;
		s->rotation.z *= inv_ql; s->rotation.w *= inv_ql;
		// Recompute world-space inertia from updated rotation so substep 2+ has fresh values.
		body_compute_inv_inertia_world(h, s);
	}
}

static int perf_initialized;

// -----------------------------------------------------------------------------
// Thread pool for parallel PGS solver (Box2D/BEPU-style spin-wait).

#ifdef _WIN32
#include <windows.h>
#include <intrin.h>

#define SOLVER_MAX_THREADS 16
#define SOLVER_BLOCK_SIZE  8

typedef struct WorkBlock { int start, count; volatile long sync_index; } WorkBlock;
typedef struct WorkStage { WorkBlock* blocks; int block_count; volatile long completion_count; } WorkStage;

// Generic work function: called per block with user context + block range.
typedef void (*WorkFn)(void* ctx, int start, int count);

typedef struct ThreadPoolCtx
{
	volatile long sync_bits;  // (stage_index << 16) | sync_counter. UINT_MAX = exit.
	WorkStage* stage;         // current stage
	WorkFn fn;                // work function for current dispatch
	void* fn_ctx;             // user context passed to fn
	int worker_count;
	char _pad[40];
} ThreadPoolCtx;

// Execute blocks: claim via CAS, call work function.
static void pool_execute(ThreadPoolCtx* ctx, int prev_sync, int cur_sync)
{
	WorkStage* stage = ctx->stage;
	int completed = 0;
	for (int bi = 0; bi < stage->block_count; bi++) {
		WorkBlock* blk = &stage->blocks[bi];
		if (_InterlockedCompareExchange(&blk->sync_index, cur_sync, prev_sync) != prev_sync) continue;
		ctx->fn(ctx->fn_ctx, blk->start, blk->count);
		completed++;
	}
	for (int bi = stage->block_count - 1; bi >= 0; bi--) {
		WorkBlock* blk = &stage->blocks[bi];
		if (_InterlockedCompareExchange(&blk->sync_index, cur_sync, prev_sync) != prev_sync) continue;
		ctx->fn(ctx->fn_ctx, blk->start, blk->count);
		completed++;
	}
	if (completed > 0) _InterlockedExchangeAdd(&stage->completion_count, completed);
}

static DWORD WINAPI pool_worker_thread(LPVOID param)
{
	ThreadPoolCtx* ctx = (ThreadPoolCtx*)param;
	long last_sync = 0;
	for (;;) {
		long sync = ctx->sync_bits;
		if (sync == (long)0xFFFFFFFF) break;
		if (sync == last_sync) { _mm_pause(); continue; }
		int cur_sync = (sync >> 16) & 0xFFFF;
		pool_execute(ctx, cur_sync - 1, cur_sync);
		last_sync = sync;
	}
	return 0;
}

static HANDLE pool_threads[SOLVER_MAX_THREADS];
static ThreadPoolCtx pool_ctx;
static int pool_thread_count;

static void pool_ensure(int n_workers)
{
	if (n_workers <= 1) return;
	if (n_workers > SOLVER_MAX_THREADS) n_workers = SOLVER_MAX_THREADS;
	int needed = n_workers - 1;
	if (pool_thread_count >= needed) return;
	pool_ctx.worker_count = n_workers;
	pool_ctx.sync_bits = 0;
	for (int i = pool_thread_count; i < needed; i++) {
		pool_threads[i] = CreateThread(NULL, 0, pool_worker_thread, &pool_ctx, 0, NULL);
	}
	pool_thread_count = needed;
}

// Persistent dispatch state — survives across pool_dispatch calls so workers
// never access freed stack memory.
static int pool_sync_counter;
static WorkBlock pool_blocks[512];
static WorkStage pool_stage;

static void pool_dispatch(WorkFn fn, void* ctx, int total_items, int block_size, int n_workers)
{
	if (total_items <= 0) return;
	int n_blocks = (total_items + block_size - 1) / block_size;
	if (n_blocks > 512) n_blocks = 512;
	for (int i = 0; i < n_blocks; i++) { pool_blocks[i].start = i * block_size; int rem = total_items - pool_blocks[i].start; pool_blocks[i].count = rem < block_size ? rem : block_size; pool_blocks[i].sync_index = pool_sync_counter; }
	pool_stage.blocks = pool_blocks;
	pool_stage.block_count = n_blocks;
	pool_stage.completion_count = 0;
	pool_ctx.stage = &pool_stage;
	pool_ctx.fn = fn;
	pool_ctx.fn_ctx = ctx;
	pool_sync_counter++;
	long sync_bits = (long)((pool_sync_counter << 16) | 0);
	_InterlockedExchange(&pool_ctx.sync_bits, sync_bits);
	pool_execute(&pool_ctx, pool_sync_counter - 1, pool_sync_counter);
	while (pool_stage.completion_count < n_blocks) _mm_pause();
}

// --- PGS solver work function (per-color) ---
typedef struct PGS_WorkCtx { BodyHot* bodies; PGS_Batch4* batches; SolverManifold* sm; SolverContact* sc; int scatter; } PGS_WorkCtx;
static void pgs_work_fn(void* ctx, int start, int count)
{
	PGS_WorkCtx* p = (PGS_WorkCtx*)ctx;
	for (int i = start; i < start + count; i++) {
		solve_contact_batch4_sv(p->bodies, &p->batches[i]);
		if (p->scatter) {
			PGS_Batch4* bt = &p->batches[i];
			for (int j = 0; j < bt->lane_count; j++) {
				int mi = bt->manifold_idx[j];
				p->sm[mi].lambda_t1 = SIMD_LANE(bt->lambda_t1, j);
				p->sm[mi].lambda_t2 = SIMD_LANE(bt->lambda_t2, j);
				p->sm[mi].lambda_twist = SIMD_LANE(bt->lambda_twist, j);
				for (int cp2 = 0; cp2 < bt->max_contacts && cp2 < p->sm[mi].contact_count; cp2++) {
					p->sc[p->sm[mi].contact_start + cp2].lambda_n = SIMD_LANE(bt->cp[cp2].lambda_n, j);
				}
			}
		}
	}
}

// --- Per-island PGS work function (parallel island dispatch) ---
// Each worker claims a contiguous range of SolveIslands and runs the full
// velocity_iters x color_count x batch pipeline for each, end-to-end. Islands
// have disjoint bodies, manifolds, and joints, so per-body and per-manifold
// writes are non-overlapping across workers; within a color, graph coloring
// guarantees disjoint body writes across batches.
typedef struct IslandSolveCtx
{
	WorldInternal* w;
	SolveIsland* solve_islands;
	ConstraintRef* crefs;
	SolverManifold* sm;
	SolverContact* sc;
	SolverJoint* sol_joints;
	PGS_Batch4* simd_batches;
	int velocity_iters;
} IslandSolveCtx;

static void island_solve_work_fn(void* ctx, int start, int count)
{
	IslandSolveCtx* ic = (IslandSolveCtx*)ctx;
	WorldInternal* w = ic->w;
	int last_iter = ic->velocity_iters - 1;
	for (int idx = start; idx < start + count; idx++) {
		SolveIsland* island = &ic->solve_islands[idx];
		for (int iter = 0; iter <= last_iter; iter++) {
			int do_scatter = (iter == last_iter);
			for (int c = 0; c < island->color_count; c++) {
				int bs = island->simd_batch_start + island->simd_color_batch_starts[c];
				int be = island->simd_batch_start + island->simd_color_batch_starts[c + 1];
				for (int bi = bs; bi < be; bi++) {
					solve_contact_batch4_sv(w->body_hot, &ic->simd_batches[bi]);
					if (do_scatter) {
						PGS_Batch4* bt = &ic->simd_batches[bi];
						for (int j = 0; j < bt->lane_count; j++) {
							int mi = bt->manifold_idx[j];
							ic->sm[mi].lambda_t1 = SIMD_LANE(bt->lambda_t1, j);
							ic->sm[mi].lambda_t2 = SIMD_LANE(bt->lambda_t2, j);
							ic->sm[mi].lambda_twist = SIMD_LANE(bt->lambda_twist, j);
							for (int cp2 = 0; cp2 < bt->max_contacts && cp2 < ic->sm[mi].contact_count; cp2++) {
								ic->sc[ic->sm[mi].contact_start + cp2].lambda_n = SIMD_LANE(bt->cp[cp2].lambda_n, j);
							}
						}
					}
				}
				// Scalar joint solves for this color in this island.
				for (int i = island->contact_end[c]; i < island->batch_starts[c + 1]; i++) {
					solve_joint(w, &ic->sol_joints[ic->crefs[i].index]);
				}
			}
		}
	}
}

// --- Batch refresh work function (parallel batch refresh for substep 2+) ---
typedef struct RefreshCtx { PGS_Batch4* batches; SolverManifold* sm; SolverContact* sc; } RefreshCtx;
static void refresh_work_fn(void* ctx, int start, int count)
{
	RefreshCtx* r = (RefreshCtx*)ctx;
	for (int i = start; i < start + count; i++) pgs_batch4_refresh(&r->batches[i], r->sm, r->sc);
}

// --- Integrate work function (parallel body integration) ---
typedef struct IntegrateCtx { WorldInternal* w; float dt; int* body_indices; int mode; } IntegrateCtx;
static void integrate_vel_work_fn(void* ctx, int start, int count)
{
	IntegrateCtx* ic = (IntegrateCtx*)ctx;
	for (int k = start; k < start + count; k++) {
		int i = ic->body_indices ? ic->body_indices[k] : k;
		BodyHot* h = &ic->w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		BodyState* s = &ic->w->body_state[i];
		body_compute_inv_inertia_world(h, s);
		h->velocity = add(h->velocity, scale(ic->w->gravity, ic->dt));
		if (s->linear_damping > 0.0f) h->velocity = scale(h->velocity, 1.0f / (1.0f + s->linear_damping * ic->dt));
		if (s->angular_damping > 0.0f) h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + s->angular_damping * ic->dt));
	}
}

static void integrate_pos_work_fn(void* ctx, int start, int count)
{
	IntegrateCtx* ic = (IntegrateCtx*)ctx;
	float dt = ic->dt;
	for (int k = start; k < start + count; k++) {
		int i = ic->body_indices ? ic->body_indices[k] : k;
		BodyHot* h = &ic->w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		BodyState* s = &ic->w->body_state[i];
		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL) h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL) h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));
		s->position = add(s->position, scale(h->velocity, dt));
		// Skip gyroscopic for uniform inertia (cubes/spheres: cross(w,I*w)=0) or negligible angular vel.
		if (av2 > 0.01f && !(s->inv_inertia_local.x == s->inv_inertia_local.y && s->inv_inertia_local.y == s->inv_inertia_local.z))
			h->angular_velocity = solve_gyroscopic(s->rotation, s->inv_inertia_local, h->angular_velocity, dt);
		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, s->rotation);
		s->rotation.x += 0.5f * dt * dq.x; s->rotation.y += 0.5f * dt * dq.y; s->rotation.z += 0.5f * dt * dq.z; s->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(s->rotation.x*s->rotation.x + s->rotation.y*s->rotation.y + s->rotation.z*s->rotation.z + s->rotation.w*s->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		s->rotation.x *= inv_ql; s->rotation.y *= inv_ql; s->rotation.z *= inv_ql; s->rotation.w *= inv_ql;
		body_compute_inv_inertia_world(h, s);
	}
}

// --- Pre-solve work function (parallel manifold setup) ---
typedef struct PreSolveCtx { WorldInternal* w; InternalManifold* manifolds; SolverManifold* sm; SolverContact* sc; float dt; float soft_dd, bias_dd, soft_ds, bias_ds; } PreSolveCtx;
static void pre_solve_work_fn(void* ctx, int start, int count)
{
	PreSolveCtx* ps = (PreSolveCtx*)ctx;
	for (int i = start; i < start + count; i++) {
		pre_solve_manifold(ps->w, &ps->manifolds[i], i, ps->sm, ps->sc, ps->dt, ps->soft_dd, ps->bias_dd, ps->soft_ds, ps->bias_ds);
	}
}

// Parallel pre_solve: alloc fixed-stride → dispatch → sequential warm start.
static void solver_pre_solve_dispatch(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, float dt, WorkFn work_fn)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	afit(sm, manifold_count); asetlen(sm, manifold_count);
	int total_contacts = manifold_count * MAX_CONTACTS;
	afit(sc, total_contacts); asetlen(sc, total_contacts);
	memset(sm, 0, manifold_count * sizeof(SolverManifold));
	memset(sc, 0, total_contacts * sizeof(SolverContact));
	float soft_dd = 0, bias_dd = 0, soft_ds = 0, bias_ds = 0;
	if (w->solver_type != SOLVER_SI) {
		float h1 = w->contact_hertz, h2 = h1 * 2.0f, dr = w->contact_damping_ratio;
		float o1 = 6.28318530718f * h1, o2 = 6.28318530718f * h2;
		float den1 = dt*2*dr*o1 + dt*dt*o1*o1, den2 = dt*2*dr*o2 + dt*dt*o2*o2;
		if (den1 > 1e-12f) { soft_dd = 1.0f / den1; bias_dd = dt * o1*o1 * soft_dd; }
		if (den2 > 1e-12f) { soft_ds = 1.0f / den2; bias_ds = dt * o2*o2 * soft_ds; }
	}
	PreSolveCtx ps_ctx = { .w = w, .manifolds = manifolds, .sm = sm, .sc = sc, .dt = dt, .soft_dd = soft_dd, .bias_dd = bias_dd, .soft_ds = soft_ds, .bias_ds = bias_ds };
	pool_dispatch(work_fn, &ps_ctx, manifold_count, 32, w->thread_count);
	// Warm start (sequential — modifies shared body velocities)
	for (int i = 0; i < manifold_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue;
		BodyHot* a = &w->body_hot[m->body_a]; BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) { SolverContact* s = &sc[m->contact_start + ci]; if (s->lambda_n == 0.0f) continue; apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n)); }
		if (m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f) apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2)));
		if (m->lambda_twist != 0.0f) { v3 tw = scale(m->normal, m->lambda_twist); BodyState* sa = &w->body_state[m->body_a]; BodyState* sb = &w->body_state[m->body_b]; a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(sa->rotation, sa->inv_inertia_local, tw)); b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(sb->rotation, sb->inv_inertia_local, tw)); }
	}
	*out_sm = sm; *out_sc = sc;
}

// --- Narrowphase work function ---
typedef struct NP_WorkCtx
{
	WorldInternal* w;
	BroadPair* pairs;
	InternalManifold* output;       // pre-allocated output array (n_pairs slots)
} NP_WorkCtx;

// Narrowphase work: each pair writes to its fixed slot (output[i]) so output
// order matches serial pair-emission order regardless of worker scheduling --
// required for bit-identical cross-thread-count results. Non-hitting pairs
// leave their slot zero-initialized (contact_count=0); caller compacts.
//
// Calls the full narrowphase_pair so warm-cache lookups, SAT hints,
// incremental NP refresh, and EPA backend selection all match the serial
// path. Thread-safe because narrowphase_pair only reads the warm_cache map
// structure (no map_set / map_del) and only writes to each pair's own
// WarmManifold entry (disjoint across workers).
static void np_work_fn(void* ctx, int start, int count)
{
	NP_WorkCtx* np = (NP_WorkCtx*)ctx;
	CK_DYNA InternalManifold* local = NULL;
	for (int i = start; i < start + count; i++) {
		BroadPair* p = &np->pairs[i];
		aclear(local);
		narrowphase_pair(np->w, p->a, p->b, &local);
		if (asize(local) > 0) np->output[i] = local[0];
	}
	afree(local);
}

#endif // _WIN32

void world_step(World world, float dt)
{
	if (!perf_initialized) { perf_init(); perf_initialized = 1; }
	double t_total = perf_now();

	WorldInternal* w = (WorldInternal*)world.id;
	w->frame++;
	int n_sub = w->sub_steps;
	float sub_dt = dt / (float)n_sub;

	afree(w->dbg_solver_manifolds); w->dbg_solver_manifolds = NULL;
	afree(w->dbg_solver_contacts);  w->dbg_solver_contacts = NULL;
	afree(w->dbg_solver_joints);    w->dbg_solver_joints = NULL;

	double t0 = perf_now();
	warm_cache_age_and_evict(w);
	int n_workers = w->thread_count > 0 ? w->thread_count : 1;
#ifdef _WIN32
	if (n_workers > 1) pool_ensure(n_workers);
#endif
	// Ensure one scratch arena per worker, lazily. 4 MB each; grows capacity
	// only when thread_count rises. Reset (not freed) every step so solver-temp
	// arrays bump-allocate with zero malloc traffic.
	for (int wi = asize(w->worker_arenas); wi < n_workers; wi++) {
		Arena a = {0};
		arena_init(&a, 4 * 1024 * 1024);
		apush(w->worker_arenas, a);
	}
	for (int wi = 0; wi < n_workers; wi++) arena_reset(&w->worker_arenas[wi]);
	int body_count = asize(w->body_hot);
#ifdef _WIN32
	if (n_workers > 1 && body_count >= 256) {
		IntegrateCtx ic = { .w = w, .dt = sub_dt, .body_indices = NULL };
		pool_dispatch(integrate_vel_work_fn, &ic, body_count, 64, n_workers);
	} else
#endif
		integrate_velocities_and_inertia(w, sub_dt);
	w->perf.integrate = perf_now() - t0;

	// Ultra-fast path: skip broadphase + solver when all islands sleeping and no free bodies.
	// Check BEFORE broadphase to avoid BVH refit and SAP overhead.
	{
		int any_awake = 0;
		for (int ii = 0; ii < asize(w->islands) && !any_awake; ii++) {
			if ((w->island_gen[ii] & 1) && w->islands[ii].awake) any_awake = 1;
		}
		for (int bi2 = 0; bi2 < body_count && !any_awake; bi2++) {
			if (split_alive(w->body_gen, bi2) && body_inv_mass(w, bi2) > 0.0f && w->body_cold[bi2].island_id < 0) any_awake = 1;
		}
		if (!any_awake) {
			w->perf.broadphase = 0;
			w->perf.pre_solve = 0;
			w->perf.pgs_solve = 0;
			w->perf.position_correct = 0;
			w->perf.pgs = (PGSTimers){0};
			w->perf.islands = 0;
			w->perf.total = perf_now() - t_total;
			return;
		}
	}

	double t1 = perf_now();
	CK_DYNA InternalManifold* manifolds = NULL;
	// Parallel narrowphase: broadphase outputs pair list, we dispatch narrowphase here.
	CK_DYNA BroadPair* np_pairs = NULL;
#ifdef _WIN32
	if (n_workers > 1) {
		w->np_pairs_out = &np_pairs;
	} else
#endif
		w->np_pairs_out = NULL;

	broadphase_and_collide(w, &manifolds);
	w->np_pairs_out = NULL;

#ifdef _WIN32
	// Parallel narrowphase on collected pairs. Each pair writes to its fixed
	// slot (pair index) so output order matches the serial broadphase pair
	// order regardless of worker scheduling -- required for bit-identical
	// cross-thread-count output. Non-hitting pairs leave zero-initialized
	// slots which are compacted out afterwards in pair order.
	if (n_workers > 1 && asize(np_pairs) >= 32) {
		int existing = asize(manifolds);
		int np_count = asize(np_pairs);
		int total_cap = existing + np_count;
		afit(manifolds, total_cap); asetlen(manifolds, total_cap);
		memset(manifolds + existing, 0, (size_t)np_count * sizeof(InternalManifold));
		NP_WorkCtx np_ctx = { .w = w, .pairs = np_pairs, .output = manifolds + existing };
		pool_dispatch(np_work_fn, &np_ctx, np_count, 32, n_workers);
		// Compact out non-hitting slots, preserving pair-index order.
		int write = existing;
		for (int i = existing; i < total_cap; i++) {
			if (manifolds[i].m.count > 0) {
				if (write != i) manifolds[write] = manifolds[i];
				write++;
			}
		}
		asetlen(manifolds, write);
	} else if (asize(np_pairs) > 0) {
		for (int i = 0; i < asize(np_pairs); i++) {
			narrowphase_pair(w, np_pairs[i].a, np_pairs[i].b, &manifolds);
		}
	}
	afree(np_pairs);
#endif
	double t_iuc = perf_now();
	if (w->sleep_enabled) islands_update_contacts(w, manifolds, asize(manifolds));
	w->perf.pgs.pos_joints = perf_now() - t_iuc;
	w->perf.broadphase = perf_now() - t1;

	// Phase 1: bucket manifolds by owning island. Consumed by per-island
	// coloring / solve in later phases; computed now for wiring validation.
	// Arena is reset each step; arrays become dangling at next step start.
	int* island_manifold_offsets = NULL;
	int* island_manifold_perm = NULL;
	if (asize(manifolds) > 0) {
		islands_bucket_manifolds(w, manifolds, asize(manifolds), &w->worker_arenas[0], &island_manifold_offsets, &island_manifold_perm);
	}
	(void)island_manifold_offsets;

	// Only populate debug_contacts when the array was previously queried (non-NULL).
	// Avoids 48KB of per-frame apush when debug visualization is unused.
	if (w->debug_contacts) {
		aclear(w->debug_contacts);
		for (int i = 0; i < asize(manifolds); i++) {
			for (int c = 0; c < manifolds[i].m.count; c++) {
				apush(w->debug_contacts, manifolds[i].m.contacts[c]);
			}
		}
	}

	int manifold_count = asize(manifolds);

	// Fast path: skip entire solver when nothing to solve (all sleeping, no contacts, no joints).
	// Only safe when no islands are awake — free-falling bodies with 0 contacts still need integration.
	// When sleep is disabled all islands are awake, so skip the scan.
	int any_awake_island = !w->sleep_enabled;
	if (!any_awake_island) {
		for (int ii = 0; ii < asize(w->islands) && !any_awake_island; ii++) {
			if ((w->island_gen[ii] & 1) && w->islands[ii].awake) any_awake_island = 1;
		}
	}
	// Skip joint/body scans when we already know there's work (contacts or awake islands).
	int skip_solver = 0;
	if (manifold_count == 0 && !any_awake_island) {
		int any_active_joints = 0;
		for (int ji = 0; ji < asize(w->joints) && !any_active_joints; ji++) {
			if (split_alive(w->joint_gen, ji)) any_active_joints = 1;
		}
		int any_unisland_dynamic = 0;
		for (int bi2 = 0; bi2 < body_count && !any_unisland_dynamic; bi2++) {
			if (split_alive(w->body_gen, bi2) && body_inv_mass(w, bi2) > 0.0f && w->body_cold[bi2].island_id < 0) any_unisland_dynamic = 1;
		}
		skip_solver = !any_active_joints && !any_unisland_dynamic;
	}
	if (skip_solver) {
		w->perf.pre_solve = 0;
		w->perf.pgs_solve = 0;
		w->perf.position_correct = 0;
		w->perf.pgs = (PGSTimers){0};
		// No post-step bvh_refit needed — bodies didn't move (solver skipped).
		// The broadphase already refitted at start of this frame.
		double t4 = perf_now();
		if (w->sleep_enabled) {
			islands_try_splits(w);
			islands_evaluate_sleep(w, dt);
		}
		w->perf.islands = perf_now() - t4;
		afree(manifolds);
		w->perf.total = perf_now() - t_total;
		return;
	}

	// --- Pre-solve (once per frame, using sub_dt for softness/bias) ---
	double t2 = perf_now();
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	// Pre-solve: parallel when threading enabled (each manifold writes to fixed-stride slots).
#ifdef _WIN32
	if (n_workers > 1 && manifold_count >= 64)
		solver_pre_solve_dispatch(w, manifolds, manifold_count, &sm, &sc, sub_dt, pre_solve_work_fn);
	else
#endif
		solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, sub_dt);

	SolverJoint* sol_joints = NULL;
	joints_pre_solve(w, sub_dt, &sol_joints);

	// Phase 1: bucket sol_joints by owning island. Paired with
	// island_manifold_perm to drive a fully island-grouped crefs layout.
	int* island_joint_offsets = NULL;
	int* island_joint_perm = NULL;
	if (asize(sol_joints) > 0) {
		islands_bucket_sol_joints(w, sol_joints, asize(sol_joints), &w->worker_arenas[0], &island_joint_offsets, &island_joint_perm);
	}

	// LDL is a direct solver -- rigid joints don't need warm-start. Stale
	// warm-start impulses inject energy when lever arms rotate between frames.
	// Only zero bilateral DOFs that LDL handles; preserve limit/motor DOF warm-start.
	if (w->ldl_enabled) {
		for (int i = 0; i < asize(sol_joints); i++) {
			if (sol_joints[i].softness != 0.0f) continue;
			JointInternal* j = &w->joints[sol_joints[i].joint_idx];
			int ldl_dof = sol_joints[i].dof;
			if (j->type == JOINT_HINGE && (j->hinge.limit_min != 0 || j->hinge.limit_max != 0 || j->hinge.motor_max_impulse > 0)) ldl_dof = 5;
			if (j->type == JOINT_PRISMATIC && j->prismatic.motor_max_impulse > 0) ldl_dof = 5;
			for (int d = 0; d < ldl_dof; d++) sol_joints[i].lambda[d] = 0;
		}
	}
	w->perf.pgs.pre_solve = perf_now() - t2;

	double t_ws = perf_now();
	joints_warm_start(w, sol_joints, asize(sol_joints));
	w->perf.pgs.warm_start = perf_now() - t_ws;

	// --- Graph color (once per frame) ---
	// When LDL enabled: rigid joints excluded from PGS (they get diagonal GS + LDL K^-1).
	// Soft spring joints stay in PGS (LDL handles them via softness term).
	double t_gc = perf_now();
	int count = asize(w->body_hot);
	CK_DYNA ConstraintRef* crefs = NULL;
	int sm_count = asize(sm);
	int sj_count = asize(sol_joints);
	int n_islands = asize(w->islands);

	// Per-island solver state for this step. One entry per non-empty island
	// slice, plus optionally a trailing entry for orphan refs (manifolds/joints
	// whose bodies are outside any island). Drives per-island coloring, SIMD
	// batch building, and the per-island PGS iter loop. Arena-owned.
	CK_DYNA SolveIsland* solve_islands = NULL;
	Arena* main_arena = &w->worker_arenas[0];

	// Walk islands, emitting manifold-refs then joint-refs per island. Each
	// SolveIsland records its absolute slice into the flat crefs array.
	// Greedy coloring output remains bit-identical because islands have
	// disjoint bodies -- processing order across islands cannot affect
	// per-ref colors. Within an island, manifolds appear in ascending index
	// (stable partition) and joints follow in sol_joints order, matching the
	// prior flat ordering restricted to that island.
	//
	// PGS eligibility for joints: LDL is a direct solver for rigid equality
	// constraints only, so joints with any bounded DOF (limits / motors) are
	// shunted to PGS via joint_internal_has_limits.
	for (int isl = 0; isl < n_islands; isl++) {
		int slice_start = asize(crefs);
		if (island_manifold_perm) {
			int ms = island_manifold_offsets[isl];
			int me = island_manifold_offsets[isl + 1];
			for (int k = ms; k < me; k++) {
				int i = island_manifold_perm[k];
				if (sm[i].contact_count == 0) continue;
				ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
					.body_a = sm[i].body_a, .body_b = sm[i].body_b };
				apush(crefs, r);
			}
		}
		if (island_joint_perm) {
			int js = island_joint_offsets[isl];
			int je = island_joint_offsets[isl + 1];
			for (int k = js; k < je; k++) {
				int i = island_joint_perm[k];
				if (w->ldl_enabled && !joint_internal_has_limits(&w->joints[sol_joints[i].joint_idx])) continue;
				ConstraintRef r = { .type = CTYPE_JOINT, .index = i,
					.body_a = sol_joints[i].body_a, .body_b = sol_joints[i].body_b };
				apush(crefs, r);
			}
		}
		int slice_count = asize(crefs) - slice_start;
		if (slice_count > 0) {
			SolveIsland si = {0};
			si.island_id = isl;
			si.crefs_start = slice_start;
			si.crefs_count = slice_count;
			apush(solve_islands, si);
		}
	}
	// Orphan manifolds + joints (bodies outside any island) emit into a single
	// trailing slice. Matches prior behavior: manifolds first (linear scan
	// fallback when bucketing was skipped), then orphan joints in sol_joints
	// order.
	{
		int slice_start = asize(crefs);
		if (!island_manifold_perm) {
			for (int i = 0; i < sm_count; i++) {
				if (sm[i].contact_count == 0) continue;
				ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
					.body_a = sm[i].body_a, .body_b = sm[i].body_b };
				apush(crefs, r);
			}
		}
		for (int i = 0; i < sj_count; i++) {
			int isl = pair_owning_island(w, sol_joints[i].body_a, sol_joints[i].body_b);
			if (isl >= 0) continue;
			if (w->ldl_enabled && !joint_internal_has_limits(&w->joints[sol_joints[i].joint_idx])) continue;
			ConstraintRef r = { .type = CTYPE_JOINT, .index = i,
				.body_a = sol_joints[i].body_a, .body_b = sol_joints[i].body_b };
			apush(crefs, r);
		}
		int slice_count = asize(crefs) - slice_start;
		if (slice_count > 0) {
			SolveIsland si = {0};
			si.island_id = -1;
			si.crefs_start = slice_start;
			si.crefs_count = slice_count;
			apush(solve_islands, si);
		}
	}

	int cref_count = asize(crefs);
	int solve_island_count = asize(solve_islands);

	// Per-island coloring with a SHARED body_colors bitmap. Islands have
	// disjoint bodies so bits set by one island's coloring are never read by
	// another; bit-identical to isolated per-island coloring or the prior
	// flat coloring over the same ref order.
	uint64_t* body_colors = (uint64_t*)arena_alloc(main_arena, (size_t)count * sizeof(uint64_t), 16);
	memset(body_colors, 0, count * sizeof(uint64_t));
	for (int sii = 0; sii < solve_island_count; sii++) {
		SolveIsland* island = &solve_islands[sii];
		int rel_batch_starts[SOLVE_ISLAND_MAX_COLORS + 1];
		int island_color_count = 0;
		color_constraints_island(crefs + island->crefs_start, island->crefs_count, body_colors, main_arena, rel_batch_starts, &island_color_count);
		island->color_count = island_color_count;
		for (int c = 0; c <= island_color_count; c++) island->batch_starts[c] = island->crefs_start + rel_batch_starts[c];

		// Partition each color's refs so contacts precede joints within this
		// island. Swap preserves color validity (coloring is per-body-pair,
		// not per-type).
		for (int c = 0; c < island_color_count; c++) {
			int i = island->batch_starts[c], j = island->batch_starts[c + 1] - 1;
			while (i <= j) {
				if (crefs[i].type == CTYPE_CONTACT) { i++; }
				else if (crefs[j].type == CTYPE_JOINT) { j--; }
				else { ConstraintRef tmp = crefs[i]; crefs[i] = crefs[j]; crefs[j] = tmp; i++; j--; }
			}
			island->contact_end[c] = i;
		}
	}
	w->perf.pgs.graph_color = perf_now() - t_gc;

	w->perf.pre_solve = perf_now() - t2;

	// --- Sub-step loop ---
	// Unified path: PGS iterates all constraints (contacts + joints).
	// When LDL enabled, K is factored once at substep start, and a mid-loop
	// K^-1 residual correction is applied at the configured iteration.
	double t_pgs = 0, t_pos = 0, t_int_sub = 0;
	double t_jlim = 0, t_ldl = 0, t_relax = 0, t_posJ = 0;
#if SIMD_SSE
	// SIMD path handles mixed contacts + joints via per-color partition:
	// contacts batch4, joints scalar within the same color. Layout is
	// (solve_island, color, batch-in-color); offsets live on each SolveIsland.
	CK_DYNA PGS_Batch4* simd_batches = NULL;
	int simd_batch_count = 0;
#endif
	for (int sub = 0; sub < n_sub; sub++) {
		if (sub > 0) {
			double ti = perf_now();
#ifdef _WIN32
			if (n_workers > 1 && body_count >= 256) {
				IntegrateCtx ic = { .w = w, .dt = sub_dt, .body_indices = NULL };
				pool_dispatch(integrate_vel_work_fn, &ic, body_count, 64, n_workers);
			} else
#endif
				integrate_velocities_and_inertia(w, sub_dt);
			t_int_sub += perf_now() - ti;
			// Refresh joint Jacobians/limits from current body state (positions
			// changed by integrate_positions last substep, velocities just updated).
			joints_refresh_substep(w, sol_joints, asize(sol_joints), sub_dt);
		}

		int has_ldl = w->ldl_enabled && asize(sol_joints) > 0;

		// Resolve LDL correction iteration: -2 = auto (velocity_iters/2), -1 = after loop
		int ldl_iter = w->ldl_correction_iter;
		if (ldl_iter == -2) ldl_iter = w->velocity_iters / 2;

		// LDL: factor K once at start of substep (topology + numeric)
		double tl0 = perf_now();
		if (has_ldl)
			ldl_factor(w, sol_joints, asize(sol_joints), sub, sub_dt);

		// LDL: direct solve for joints (velocity correction).
		// Runs before PGS so joint impulses are already applied.
		if (has_ldl) {
			if (sub > 0) {
				for (int i = 0; i < asize(sol_joints); i++) if (sol_joints[i].softness > 0.0f) for (int d = 0; d < sol_joints[i].dof; d++) sol_joints[i].lambda[d] = 0;
			}
			ldl_velocity_correct(w, sol_joints, asize(sol_joints), sub_dt);
		}
		t_ldl += perf_now() - tl0;

		// PGS: iterate all constraints (contacts, and joints when LDL is off).
		// PGS iteration: BodyHot is now the lean solver array (same layout
		// as old SolverBodyVel), so no sync needed.

		// Build SIMD batches on first substep; refresh only bias/lambda on substep 2+.
		// Layout: batches are grouped by (solve_island, color, batch-in-color).
		// Each SolveIsland records simd_batch_start (absolute offset into the
		// flat simd_batches array) and simd_color_batch_starts[c] (relative
		// offsets within its own batch range).
#if SIMD_SSE
		if (sub == 0 && cref_count > 0) {
			int total_batches = 0;
			for (int sii = 0; sii < solve_island_count; sii++) {
				SolveIsland* island = &solve_islands[sii];
				island->simd_batch_start = total_batches;
				island->simd_color_batch_starts[0] = 0;
				for (int c = 0; c < island->color_count; c++) {
					int contact_n = island->contact_end[c] - island->batch_starts[c];
					int batches_in_color = (contact_n + 3) / 4;
					total_batches += batches_in_color;
					island->simd_color_batch_starts[c + 1] = island->simd_color_batch_starts[c] + batches_in_color;
				}
				island->simd_batch_count = island->simd_color_batch_starts[island->color_count];
			}
			if (total_batches > 0) {
				afit(simd_batches, total_batches); asetlen(simd_batches, total_batches);
				for (int sii = 0; sii < solve_island_count; sii++) {
					SolveIsland* island = &solve_islands[sii];
					for (int c = 0; c < island->color_count; c++) {
						int start = island->batch_starts[c], end = island->contact_end[c];
						int bi_out = island->simd_batch_start + island->simd_color_batch_starts[c];
						for (int i = start; i + 3 < end; i += 4) {
							int idx[4] = { crefs[i].index, crefs[i+1].index, crefs[i+2].index, crefs[i+3].index };
							pgs_batch4_prepare(&simd_batches[bi_out++], sm, idx, 4, sc);
						}
						int rem_start = start + ((end - start) / 4) * 4;
						if (rem_start < end) {
							int idx[4] = {0};
							int rem_count = end - rem_start;
							for (int j = 0; j < rem_count; j++) idx[j] = crefs[rem_start + j].index;
							pgs_batch4_prepare(&simd_batches[bi_out++], sm, idx, rem_count, sc);
						}
					}
				}
				simd_batch_count = asize(simd_batches);
			}
		}
#endif

		double tp = perf_now();
#if SIMD_SSE
		// Refresh on substep 2+. Batches are all addressed into sm/sc by
		// their stored indices, so a flat refresh pass works regardless of
		// island layout.
		if (sub > 0 && simd_batch_count > 0) {
#ifdef _WIN32
			if (n_workers > 1 && simd_batch_count >= n_workers * 2) {
				RefreshCtx rc = { .batches = simd_batches, .sm = sm, .sc = sc };
				pool_dispatch(refresh_work_fn, &rc, simd_batch_count, SOLVER_BLOCK_SIZE, n_workers);
			} else
#endif
			for (int bi = 0; bi < simd_batch_count; bi++) {
				pgs_batch4_refresh(&simd_batches[bi], sm, sc);
			}
		}

		// PGS dispatch: two strategies depending on island size.
		//
		// 1) Small islands -- whole-island parallel claim. Each worker owns an
		//    island end-to-end through velocity_iters x color_count x batch.
		//    Good when there are many islands (ragdolls, mixed scenes).
		//
		// 2) Big islands -- serial outer, per-color parallel dispatch. One
		//    giant island (tall stack, dense pile) can't be parallelized
		//    whole, but its color-c batches are disjoint-body so they split
		//    cleanly across workers per color. Matches the pre-phase-2
		//    perf profile for single-giant-island scenes.
		//
		// Partition swaps big islands to the tail. Small islands preserve
		// their original order; big islands don't, but both phases produce
		// bit-identical output because islands are disjoint (order across
		// islands cannot affect per-body update sequence).
		int first_big = solve_island_count;
#ifdef _WIN32
		if (n_workers > 1) {
			int i = 0;
			while (i < first_big) {
				if (solve_islands[i].simd_batch_count >= n_workers * 2) {
					first_big--;
					SolveIsland tmp = solve_islands[i];
					solve_islands[i] = solve_islands[first_big];
					solve_islands[first_big] = tmp;
				} else {
					i++;
				}
			}
		}
#endif
		int small_count = first_big;

		// Small islands: dispatch across workers (whole-island claim).
		IslandSolveCtx island_ctx = {
			.w = w, .solve_islands = solve_islands, .crefs = crefs,
			.sm = sm, .sc = sc, .sol_joints = sol_joints,
			.simd_batches = simd_batches, .velocity_iters = w->velocity_iters,
		};
#ifdef _WIN32
		if (n_workers > 1 && small_count >= 2) {
			pool_dispatch(island_solve_work_fn, &island_ctx, small_count, 1, n_workers);
		} else
#endif
		if (small_count > 0) {
			island_solve_work_fn(&island_ctx, 0, small_count);
		}

		// Big islands: serial outer, per-color parallel dispatch.
		int last_iter = w->velocity_iters - 1;
		for (int bi = first_big; bi < solve_island_count; bi++) {
			SolveIsland* island = &solve_islands[bi];
			for (int iter = 0; iter <= last_iter; iter++) {
				int do_scatter = (iter == last_iter);
				for (int c = 0; c < island->color_count; c++) {
					int bs = island->simd_batch_start + island->simd_color_batch_starts[c];
					int be = island->simd_batch_start + island->simd_color_batch_starts[c + 1];
					int n_color_batches = be - bs;
#ifdef _WIN32
					if (n_workers > 1 && n_color_batches >= n_workers * 2) {
						PGS_WorkCtx pgs_ctx = { .bodies = w->body_hot, .batches = simd_batches + bs, .sm = sm, .sc = sc, .scatter = do_scatter };
						pool_dispatch(pgs_work_fn, &pgs_ctx, n_color_batches, SOLVER_BLOCK_SIZE, n_workers);
					} else
#endif
					{
						for (int bi2 = bs; bi2 < be; bi2++) {
							solve_contact_batch4_sv(w->body_hot, &simd_batches[bi2]);
							if (do_scatter) {
								PGS_Batch4* bt = &simd_batches[bi2];
								for (int j = 0; j < bt->lane_count; j++) {
									int mi = bt->manifold_idx[j];
									sm[mi].lambda_t1 = SIMD_LANE(bt->lambda_t1, j);
									sm[mi].lambda_t2 = SIMD_LANE(bt->lambda_t2, j);
									sm[mi].lambda_twist = SIMD_LANE(bt->lambda_twist, j);
									for (int cp2 = 0; cp2 < bt->max_contacts && cp2 < sm[mi].contact_count; cp2++) {
										sc[sm[mi].contact_start + cp2].lambda_n = SIMD_LANE(bt->cp[cp2].lambda_n, j);
									}
								}
							}
						}
					}
					for (int i = island->contact_end[c]; i < island->batch_starts[c + 1]; i++) {
						solve_joint(w, &sol_joints[crefs[i].index]);
					}
				}
			}
		}
#else
#error "SIMD_SSE required: non-SSE solver path removed; add a SIMD_NEON backend if needed"
#endif
		t_pgs += perf_now() - tp;

		double ti2 = perf_now();
#ifdef _WIN32
		if (n_workers > 1 && body_count >= 256) {
			IntegrateCtx ic = { .w = w, .dt = sub_dt, .body_indices = NULL };
			pool_dispatch(integrate_pos_work_fn, &ic, body_count, 64, n_workers);
		} else
#endif
			integrate_positions(w, sub_dt);
		t_int_sub += perf_now() - ti2;

		// Relax contacts: refresh separation/bias from updated positions
		double tr = perf_now();
		if (w->solver_type == SOLVER_SOFT_STEP)
			solver_relax_contacts(w, sm, asize(sm), sc, sub_dt);
		t_relax += perf_now() - tr;

		// Position correction after integration.
		// LDL: direct solve projects out bulk joint error, NGS cleans up residual.
		double tl1 = perf_now();
		if (has_ldl)
			ldl_position_solve(w, sol_joints, asize(sol_joints), sub_dt);
		t_ldl += perf_now() - tl1;

		double tpj = perf_now();
		joints_position_correct(w, sol_joints, asize(sol_joints), w->position_iters);
		t_posJ += perf_now() - tpj;

		t_pos += (perf_now() - tr);
	}

#if SIMD_SSE
	afree(simd_batches);
#endif
	afree(solve_islands);

	w->perf.pgs_solve = t_pgs;
	w->perf.integrate += t_int_sub;
	w->perf.position_correct = t_pos;
	w->perf.pgs.iterations = t_pgs - t_jlim;
	w->perf.pgs.joint_limits = t_jlim;
	w->perf.pgs.ldl = t_ldl;
	w->perf.pgs.relax = t_relax;
	w->perf.pgs.pos_joints = t_posJ;

	afree(crefs);

	// Position correction: NGS for hard SI; also for soft modes when contact_hertz is off
	double t3 = perf_now();
	if (w->solver_type == SOLVER_SI)
		solver_position_correct(w, sm, asize(sm), sc);
	else if (w->contact_hertz <= 0.0f)
		solver_position_correct(w, sm, asize(sm), sc);
	w->perf.pgs.pos_contacts = perf_now() - t3;
	w->perf.position_correct += perf_now() - t3;

	// Post-solve (once per frame)
	double t_ps = perf_now();
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_joints, asize(sol_joints));
	w->perf.pgs.post_solve = perf_now() - t_ps;

	// Post-step BVH refit: update leaves after all substeps so they're correct
	// for next frame's broadphase. Without this, bodies accelerated by the solver
	// can escape their fat AABBs between frames.
	bvh_refit(w->bvh_dynamic, w);

	double t4 = perf_now();
	if (w->sleep_enabled) {
		islands_try_splits(w);
		islands_evaluate_sleep(w, dt);
	}
	w->perf.islands = perf_now() - t4;

	afree(manifolds);

	w->dbg_solver_manifolds = (void *)sm;
	w->dbg_solver_contacts = (void *)sc;
	w->dbg_solver_joints = (void *)sol_joints;

	w->perf.total = perf_now() - t_total;
}

PerfTimers world_get_perf(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	return w->perf;
}



void world_set_solver_type(World world, SolverType type)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->solver_type = type;
}

// -----------------------------------------------------------------------------
// Body.

Body create_body(World world, BodyParams params)
{
	assert(is_valid(params.position) && "create_body: position is NaN/inf");
	assert(is_valid(params.rotation) && "create_body: rotation is NaN/inf");
	assert(is_valid(params.mass) && params.mass >= 0.0f && "create_body: mass must be >= 0 and finite");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	split_add(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);
	split_ensure(w->body_state, idx);

	w->body_cold[idx] = (BodyCold){
		.mass = params.mass,
		.shapes = NULL,
		.bvh_leaf = -1,
		.island_id = -1,
		.island_prev = -1,
		.island_next = -1,
		.collision_group = 0xFFFFFFFFu,
		.collision_mask  = 0xFFFFFFFFu,
		.compound_id     = 0,
	};
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	float ang_damp = params.angular_damping;
	if (ang_damp == 0.0f) ang_damp = 0.03f; // default: 3%/s (BEPU-style)
	w->body_hot[idx] = (BodyHot){
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
	};
	w->body_state[idx] = (BodyState){
		.position = params.position,
		.rotation = params.rotation,
		.friction = fric,
		.restitution = params.restitution,
		.linear_damping = params.linear_damping,
		.angular_damping = ang_damp,
		.sleep_allowed = 1,
	};
	return split_handle(Body, w->body_gen, idx);
}

void destroy_body(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	// Remove from island
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl)) {
		// Remove all joints connected to this body
		int ji = w->islands[isl].head_joint;
		while (ji >= 0) {
			int next = w->joints[ji].island_next;
			if (w->joints[ji].body_a == idx || w->joints[ji].body_b == idx) {
				island_remove_joint(w, isl, ji);
				w->islands[isl].constraint_remove_count++;
			}
			ji = next;
		}
		island_remove_body(w, isl, idx);
		w->islands[isl].constraint_remove_count++;
	}
	if (w->body_cold[idx].bvh_leaf >= 0) {
		BVH_Tree* tree;
		if (body_inv_mass(w, idx) == 0.0f) tree = w->bvh_static;
		else if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) tree = w->bvh_sleeping;
		else tree = w->bvh_dynamic;
		// body_state[idx] will be cleared by split_del below
		int moved_body = bvh_remove(tree, w->body_cold[idx].bvh_leaf);
		if (moved_body >= 0) w->body_cold[moved_body].bvh_leaf = w->body_cold[idx].bvh_leaf;
	}
	afree(w->body_cold[idx].shapes);
	split_del(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);
}

void body_add_shape(World world, Body body, ShapeParams params)
{
	assert(is_valid(params.local_pos) && "body_add_shape: local_pos is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));

	ShapeInternal s = {0};
	s.type = params.type;
	s.local_pos = params.local_pos;
	switch (params.type) {
	case SHAPE_SPHERE:  s.sphere.radius = params.sphere.radius; break;
	case SHAPE_CAPSULE: s.capsule.half_height = params.capsule.half_height;
	                    s.capsule.radius = params.capsule.radius; break;
	case SHAPE_BOX:     s.box.half_extents = params.box.half_extents; break;
	case SHAPE_HULL:    s.hull.hull = params.hull.hull;
	                    s.hull.scale = params.hull.scale; break;
	case SHAPE_CYLINDER: s.cylinder.half_height = params.cylinder.half_height;
	                     s.cylinder.radius = params.cylinder.radius; break;
	}
	apush(w->body_cold[idx].shapes, s);
	if (params.type == SHAPE_CYLINDER && w->body_state[idx].angular_damping < 0.15f)
		w->body_state[idx].angular_damping = 0.15f;
	recompute_body_inertia(w, idx);

	// Insert into BVH on first shape add.
	if (w->broadphase_type == BROADPHASE_BVH && asize(w->body_cold[idx].shapes) == 1) {
		AABB box = aabb_expand(body_aabb(&w->body_state[idx], &w->body_cold[idx]), BVH_AABB_MARGIN);
		BVH_Tree* tree = body_inv_mass(w, idx) == 0.0f ? w->bvh_static : w->bvh_dynamic;
		w->body_cold[idx].bvh_leaf = bvh_insert(tree, idx, box);
	}
}

v3 body_get_position(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return body_pos(w, idx);
}

void body_set_position(World world, Body body, v3 pos)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	body_pos(w, idx) = pos;
}

quat body_get_rotation(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return body_rot(w, idx);
}

v3 body_get_velocity(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return body_vel(w, idx);
}

v3 body_get_angular_velocity(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return body_angvel(w, idx);
}

void body_wake(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_velocity(World world, Body body, v3 vel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	body_vel(w, idx) = vel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_angular_velocity(World world, Body body, v3 avel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	body_angvel(w, idx) = avel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

int body_is_asleep(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	if (body_inv_mass(w, idx) == 0.0f) return 1; // static bodies are always "asleep"
	int isl = w->body_cold[idx].island_id;
	if (isl < 0 || !island_alive(w, isl)) return 0;
	return !w->islands[isl].awake;
}

void world_set_sleep_enabled(World world, int enabled)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->sleep_enabled = enabled;
}

int world_get_sleep_enabled(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	return w->sleep_enabled;
}

// Collision filter: two bodies collide iff (a.group & b.mask) && (b.group & a.mask).
// Default on create is group=mask=0xFFFFFFFF (collide with everything).
// Useful for category-level filtering (terrain vs character vs projectile vs ...).
void body_set_collision_filter(World world, Body body, uint32_t group, uint32_t mask)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_cold[idx].collision_group = group;
	w->body_cold[idx].collision_mask = mask;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

// Compound/instance id. Bodies sharing the same nonzero compound_id never
// collide with each other (e.g. all parts of one ragdoll). Independent of the
// group/mask filter -- both must allow collision for a pair to reach
// narrowphase. Default 0 = no compound.
//
// Scales to unlimited instances: assign each ragdoll/vehicle/chain a unique
// nonzero id; their parts phase through each other but instances of the same
// class still collide normally.
void body_set_compound_id(World world, Body body, uint32_t compound_id)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_cold[idx].compound_id = compound_id;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_sleep_allowed(World world, Body body, int allowed)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	body_sleep_allowed(w, idx) = allowed;
	if (!allowed) {
		int isl = w->body_cold[idx].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
			island_wake(w, isl);
	}
}

// -----------------------------------------------------------------------------
// Joints.

Joint create_ball_socket(World world, BallSocketParams params)
{
	assert(is_valid(params.local_offset_a) && "create_ball_socket: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_ball_socket: local_offset_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	int idx = joint_alloc_slot(w);

	w->joints[idx] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = ba, .body_b = bb,
		.ball_socket = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_distance(World world, DistanceParams params)
{
	assert(is_valid(params.local_offset_a) && "create_distance: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_distance: local_offset_b is NaN/inf");
	assert(is_valid(params.rest_length) && "create_distance: rest_length is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	int idx = joint_alloc_slot(w);

	// Auto-compute rest length if not specified
	float rest = params.rest_length;
	if (rest <= 0.0f) {
		BodyState* sa = &w->body_state[ba];
		BodyState* sb = &w->body_state[bb];
		v3 wa = add(sa->position, rotate(sa->rotation, params.local_offset_a));
		v3 wb = add(sb->position, rotate(sb->rotation, params.local_offset_b));
		rest = len(sub(wb, wa));
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_DISTANCE,
		.body_a = ba, .body_b = bb,
		.distance = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.rest_length = rest,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_hinge(World world, HingeParams params)
{
	assert(is_valid(params.local_offset_a) && "create_hinge: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_hinge: local_offset_b is NaN/inf");
	assert(is_valid(params.local_axis_a) && "create_hinge: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_hinge: local_axis_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	int idx = joint_alloc_slot(w);

	v3 axis_a_local = norm(params.local_axis_a);
	v3 axis_b_local = norm(params.local_axis_b);

	// Compute reference directions perpendicular to the hinge axis for angle measurement.
	// local_ref_a: arbitrary unit vector perpendicular to axis_a in body A's local space.
	// local_ref_b: chosen so measured angle = 0 at the initial configuration.
	v3 ref_a, ref_a_t2;
	hinge_tangent_basis(axis_a_local, &ref_a, &ref_a_t2);
	// Transform ref_a into world, then into body B's local space
	quat q_a = body_rot(w, ba);
	quat q_b = body_rot(w, bb);
	v3 ref_a_world = rotate(q_a, ref_a);
	v3 ref_b = rotate(inv(q_b), ref_a_world);

	w->joints[idx] = (JointInternal){
		.type = JOINT_HINGE,
		.body_a = ba, .body_b = bb,
		.hinge = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_axis_a = axis_a_local,
			.local_axis_b = axis_b_local,
			.local_ref_a = ref_a,
			.local_ref_b = ref_b,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_fixed(World world, FixedParams params)
{
	assert(is_valid(params.local_offset_a) && "create_fixed: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_fixed: local_offset_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	int idx = joint_alloc_slot(w);

	quat q_a = body_rot(w, ba);
	quat q_b = body_rot(w, bb);
	w->joints[idx] = (JointInternal){
		.type = JOINT_FIXED,
		.body_a = ba, .body_b = bb,
		.fixed = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_rel_quat = mul(inv(q_a), q_b),
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_prismatic(World world, PrismaticParams params)
{
	assert(is_valid(params.local_offset_a) && "create_prismatic: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_prismatic: local_offset_b is NaN/inf");
	assert(is_valid(params.local_axis_a) && "create_prismatic: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_prismatic: local_axis_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	int idx = joint_alloc_slot(w);

	quat q_a = body_rot(w, ba);
	quat q_b = body_rot(w, bb);
	w->joints[idx] = (JointInternal){
		.type = JOINT_PRISMATIC,
		.body_a = ba, .body_b = bb,
		.prismatic = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.local_rel_quat = mul(inv(q_a), q_b),
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

// Common slot allocation used by all create_* joints. Keeps the joints[],
// joint_hot[] and joint_gen[] arrays in lockstep. Hot slot is zeroed so
// warm_lambda starts fresh (no stale impulse from a destroyed joint).
static int joint_alloc_slot(WorldInternal* w)
{
	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
		w->joint_hot[idx] = (JointHot){0};
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
		JointHot hot_zero = {0};
		apush(w->joint_hot, hot_zero);
	}
	return idx;
}

Joint create_angular_motor(World world, AngularMotorParams params)
{
	assert(is_valid(params.local_axis_a) && "create_angular_motor: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_angular_motor: local_axis_b is NaN/inf");
	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	assert(ba != bb && "joint requires two distinct bodies");
	assert((body_inv_mass(w, ba) > 0.0f || body_inv_mass(w, bb) > 0.0f) && "at least one body must be dynamic");

	int idx = joint_alloc_slot(w);
	w->joints[idx] = (JointInternal){
		.type = JOINT_ANGULAR_MOTOR,
		.body_a = ba, .body_b = bb,
		.angular_motor = {
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.target_speed = params.target_speed,
			.max_impulse = params.max_impulse,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_twist_limit(World world, TwistLimitParams params)
{
	assert(is_valid(params.local_axis_a) && "create_twist_limit: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_twist_limit: local_axis_b is NaN/inf");
	assert(params.limit_min <= 0.0f && params.limit_max >= 0.0f && "twist limits must bracket zero");
	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	assert(ba != bb && "joint requires two distinct bodies");
	assert((body_inv_mass(w, ba) > 0.0f || body_inv_mass(w, bb) > 0.0f) && "at least one body must be dynamic");

	int idx = joint_alloc_slot(w);
	w->joints[idx] = (JointInternal){
		.type = JOINT_TWIST_LIMIT,
		.body_a = ba, .body_b = bb,
		.twist_limit = {
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.limit_min = params.limit_min,
			.limit_max = params.limit_max,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_cone_limit(World world, ConeLimitParams params)
{
	assert(is_valid(params.local_axis_a) && "create_cone_limit: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_cone_limit: local_axis_b is NaN/inf");
	assert(params.half_angle > 0.0f && "cone half_angle must be positive");
	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	assert(ba != bb && "joint requires two distinct bodies");
	assert((body_inv_mass(w, ba) > 0.0f || body_inv_mass(w, bb) > 0.0f) && "at least one body must be dynamic");

	int idx = joint_alloc_slot(w);
	w->joints[idx] = (JointInternal){
		.type = JOINT_CONE_LIMIT,
		.body_a = ba, .body_b = bb,
		.cone_limit = {
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.half_angle = params.half_angle,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_swing_twist(World world, SwingTwistParams params)
{
	assert(is_valid(params.local_offset_a) && "create_swing_twist: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_swing_twist: local_offset_b is NaN/inf");
	assert(is_valid(params.local_axis_a) && "create_swing_twist: local_axis_a is NaN/inf");
	assert(is_valid(params.local_axis_b) && "create_swing_twist: local_axis_b is NaN/inf");
	assert(params.cone_half_angle > 0.0f && "swing_twist cone_half_angle must be positive");
	assert(params.twist_min <= 0.0f && params.twist_max >= 0.0f && "swing_twist twist limits must bracket zero");
	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));
	assert(ba != bb && "joint requires two distinct bodies");
	assert((body_inv_mass(w, ba) > 0.0f || body_inv_mass(w, bb) > 0.0f) && "at least one body must be dynamic");

	int idx = joint_alloc_slot(w);
	w->joints[idx] = (JointInternal){
		.type = JOINT_SWING_TWIST,
		.body_a = ba, .body_b = bb,
		.swing_twist = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.local_axis_a = norm(params.local_axis_a),
			.local_axis_b = norm(params.local_axis_b),
			.cone_half_angle = params.cone_half_angle,
			.twist_min = params.twist_min,
			.twist_max = params.twist_max,
			.spring = params.spring,
		},
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	w->ldl_topo_version++;
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	unlink_joint_from_island(w, idx);
	memset(&w->joints[idx], 0, sizeof(JointInternal));
	w->joint_hot[idx] = (JointHot){0};
	w->joint_gen[idx]++; // even = dead
	w->ldl_topo_version++;
	apush(w->joint_free, idx);
}

void joint_set_hinge_limits(World world, Joint joint, float min_angle, float max_angle)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_HINGE && "joint_set_hinge_limits: not a hinge joint");
	w->joints[idx].hinge.limit_min = min_angle;
	w->joints[idx].hinge.limit_max = max_angle;
}

void joint_set_distance_limits(World world, Joint joint, float min_distance, float max_distance)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_DISTANCE && "joint_set_distance_limits: not a distance joint");
	w->joints[idx].distance.limit_min = min_distance;
	w->joints[idx].distance.limit_max = max_distance;
}

void joint_clear_limits(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	if (w->joints[idx].type == JOINT_HINGE) {
		w->joints[idx].hinge.limit_min = 0;
		w->joints[idx].hinge.limit_max = 0;
	} else if (w->joints[idx].type == JOINT_DISTANCE) {
		w->joints[idx].distance.limit_min = 0;
		w->joints[idx].distance.limit_max = 0;
	}
}

void joint_set_hinge_motor(World world, Joint joint, float speed, float max_impulse)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_HINGE);
	w->joints[idx].hinge.motor_speed = speed;
	w->joints[idx].hinge.motor_max_impulse = max_impulse;
}

void joint_set_prismatic_motor(World world, Joint joint, float speed, float max_impulse)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	assert(w->joints[idx].type == JOINT_PRISMATIC);
	w->joints[idx].prismatic.motor_speed = speed;
	w->joints[idx].prismatic.motor_max_impulse = max_impulse;
}

static void bvh_debug_walk(BVH_Tree* t, int ni, int depth, BVHDebugFn fn, void* user)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVH_Child* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		fn(c->min, c->max, depth, bvh_child_is_leaf(c), user);
		if (bvh_child_is_internal(c)) bvh_debug_walk(t, c->index, depth + 1, fn, user);
	}
}

void world_debug_bvh(World world, BVHDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (w->bvh_dynamic->root >= 0) bvh_debug_walk(w->bvh_dynamic, w->bvh_dynamic->root, 0, fn, user);
	if (w->bvh_static->root >= 0) bvh_debug_walk(w->bvh_static, w->bvh_static->root, 0, fn, user);
	if (w->bvh_sleeping->root >= 0) bvh_debug_walk(w->bvh_sleeping, w->bvh_sleeping->root, 0, fn, user);
}

void world_debug_joints(World world, JointDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int jcount = asize(w->joints);
	for (int i = 0; i < jcount; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyState* sa = &w->body_state[j->body_a];
		BodyState* sb = &w->body_state[j->body_b];
		JointDebugInfo info = {0};
		info.type = j->type;
		if (j->type == JOINT_BALL_SOCKET) {
			info.anchor_a = add(sa->position, rotate(sa->rotation, j->ball_socket.local_a));
			info.anchor_b = add(sb->position, rotate(sb->rotation, j->ball_socket.local_b));
			info.is_soft = j->ball_socket.spring.frequency > 0;
		} else if (j->type == JOINT_DISTANCE) {
			info.anchor_a = add(sa->position, rotate(sa->rotation, j->distance.local_a));
			info.anchor_b = add(sb->position, rotate(sb->rotation, j->distance.local_b));
			info.is_soft = j->distance.spring.frequency > 0;
		} else if (j->type == JOINT_HINGE) {
			info.anchor_a = add(sa->position, rotate(sa->rotation, j->hinge.local_a));
			info.anchor_b = add(sb->position, rotate(sb->rotation, j->hinge.local_b));
			info.axis_a = norm(rotate(sa->rotation, j->hinge.local_axis_a));
			info.is_soft = j->hinge.spring.frequency > 0;
			info.motor_speed = j->hinge.motor_speed;
			info.motor_max_impulse = j->hinge.motor_max_impulse;
			info.limit_min = j->hinge.limit_min;
			info.limit_max = j->hinge.limit_max;
			info.ref_a = rotate(sa->rotation, j->hinge.local_ref_a);
			info.ref_b = rotate(sb->rotation, j->hinge.local_ref_b);
			float angle = atan2f(dot(cross(info.ref_a, info.ref_b), info.axis_a), dot(info.ref_a, info.ref_b));
			info.current_angle = angle;
			info.limit_active = (j->hinge.limit_min != 0 && angle <= j->hinge.limit_min) || (j->hinge.limit_max != 0 && angle >= j->hinge.limit_max);
		} else if (j->type == JOINT_FIXED) {
			info.anchor_a = add(sa->position, rotate(sa->rotation, j->fixed.local_a));
			info.anchor_b = add(sb->position, rotate(sb->rotation, j->fixed.local_b));
			info.is_soft = j->fixed.spring.frequency > 0;
		} else if (j->type == JOINT_PRISMATIC) {
			info.anchor_a = add(sa->position, rotate(sa->rotation, j->prismatic.local_a));
			info.anchor_b = add(sb->position, rotate(sb->rotation, j->prismatic.local_b));
			info.axis_a = norm(rotate(sa->rotation, j->prismatic.local_axis_a));
			info.is_soft = j->prismatic.spring.frequency > 0;
			info.motor_speed = j->prismatic.motor_speed;
			info.motor_max_impulse = j->prismatic.motor_max_impulse;
		}
		fn(info, user);
	}
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	// Ensure the array exists so world_step knows to populate it next frame.
	if (!w->debug_contacts) { afit(w->debug_contacts, 1); asetlen(w->debug_contacts, 0); }
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}

// -----------------------------------------------------------------------------
// World queries.

int world_query_aabb(World world, v3 lo, v3 hi, Body* results, int max_results)
{
	WorldInternal* w = (WorldInternal*)world.id;
	AABB query = { lo, hi };
	CK_DYNA int* candidates = NULL;
	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_query_aabb(w->bvh_dynamic, query, &candidates);
		bvh_query_aabb(w->bvh_static, query, &candidates);
		bvh_query_aabb(w->bvh_sleeping, query, &candidates);
	} else {
		int count = asize(w->body_hot);
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->body_gen, i)) continue;
			if (asize(w->body_cold[i].shapes) == 0) continue;
			apush(candidates, i);
		}
	}
	int total = 0;
	for (int i = 0; i < asize(candidates); i++) {
		int idx = candidates[i];
		AABB b = body_aabb(&w->body_state[idx], &w->body_cold[idx]);
		if (!aabb_overlaps(query, b)) continue;
		if (total < max_results)
			results[total] = split_handle(Body, w->body_gen, idx);
		total++;
	}
	afree(candidates);
	return total;
}

int world_raycast(World world, v3 origin, v3 direction, float max_distance, RayHit* hit)
{
	WorldInternal* w = (WorldInternal*)world.id;
	float dl = len(direction);
	if (dl < 1e-12f) return 0;
	v3 dir = scale(direction, 1.0f / dl);
	v3 inv_dir = V3(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z);

	CK_DYNA int* candidates = NULL;
	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_query_ray(w->bvh_dynamic, origin, inv_dir, max_distance, &candidates);
		bvh_query_ray(w->bvh_static, origin, inv_dir, max_distance, &candidates);
		bvh_query_ray(w->bvh_sleeping, origin, inv_dir, max_distance, &candidates);
	} else {
		int count = asize(w->body_hot);
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->body_gen, i)) continue;
			if (asize(w->body_cold[i].shapes) == 0) continue;
			apush(candidates, i);
		}
	}

	float best_t = max_distance;
	v3 best_n = {0};
	int best_idx = -1;
	for (int i = 0; i < asize(candidates); i++) {
		int idx = candidates[i];
		float t; v3 n;
		if (ray_body(w, idx, origin, dir, best_t, &t, &n)) {
			best_t = t; best_n = n; best_idx = idx;
		}
	}
	afree(candidates);

	if (best_idx < 0) return 0;
	if (hit) {
		hit->body = split_handle(Body, w->body_gen, best_idx);
		hit->point = add(origin, scale(dir, best_t));
		hit->normal = best_n;
		hit->distance = best_t;
	}
	return 1;
}
