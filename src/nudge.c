// See LICENSE for licensing info.
// nudge.c -- physics world implementation

#include "perf.h"
#include "nudge_internal.h"
#include "gjk.c"
#include "gjk_batch.c"
#include "quickhull.c"
#include "bvh.c"
#include "collision.c"
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
	w->solver_type = params.solver_type;
	w->sleep_enabled = 1;
	w->sat_hint_enabled = 1;
	w->sat_hillclimb_enabled = 1;
	w->warm_start_enabled = 1;
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
	afree(w->body_vel);
	map_free(w->warm_cache);
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	bvh_free(w->bvh_sleeping); CK_FREE(w->bvh_sleeping);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	for (int i = 0; i < asize(w->islands); i++) ldl_cache_free(&w->islands[i].ldl);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	map_free(w->joint_pairs);
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
		body_compute_inv_inertia_world(h);
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
		if (h->linear_damping > 0.0f)
			h->velocity = scale(h->velocity, 1.0f / (1.0f + h->linear_damping * dt));
		if (h->angular_damping > 0.0f)
			h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + h->angular_damping * dt));
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

		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL)
			h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL)
			h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));

		h->position = add(h->position, scale(h->velocity, dt));

		// Skip gyroscopic for uniform inertia (cubes/spheres: cross(w,I*w)=0) or negligible angular vel.
		if (av2 > 0.01f && !(h->inv_inertia_local.x == h->inv_inertia_local.y && h->inv_inertia_local.y == h->inv_inertia_local.z))
			h->angular_velocity = solve_gyroscopic(h->rotation, h->inv_inertia_local, h->angular_velocity, dt);

		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, h->rotation);
		h->rotation.x += 0.5f * dt * dq.x;
		h->rotation.y += 0.5f * dt * dq.y;
		h->rotation.z += 0.5f * dt * dq.z;
		h->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y
			+ h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql;
		h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
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
	for (int i = pool_thread_count; i < needed; i++)
		pool_threads[i] = CreateThread(NULL, 0, pool_worker_thread, &pool_ctx, 0, NULL);
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
typedef struct PGS_WorkCtx { SolverBodyVel* bodies; PGS_Batch4* batches; } PGS_WorkCtx;
static void pgs_work_fn(void* ctx, int start, int count)
{
	PGS_WorkCtx* p = (PGS_WorkCtx*)ctx;
	for (int i = start; i < start + count; i++)
		solve_contact_batch4_sv(p->bodies, &p->batches[i]);
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
		body_compute_inv_inertia_world(h);
		h->velocity = add(h->velocity, scale(ic->w->gravity, ic->dt));
		if (h->linear_damping > 0.0f) h->velocity = scale(h->velocity, 1.0f / (1.0f + h->linear_damping * ic->dt));
		if (h->angular_damping > 0.0f) h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + h->angular_damping * ic->dt));
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
		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL) h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL) h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));
		h->position = add(h->position, scale(h->velocity, dt));
		// Skip gyroscopic for uniform inertia (cubes/spheres: cross(w,I*w)=0) or negligible angular vel.
		if (av2 > 0.01f && !(h->inv_inertia_local.x == h->inv_inertia_local.y && h->inv_inertia_local.y == h->inv_inertia_local.z))
			h->angular_velocity = solve_gyroscopic(h->rotation, h->inv_inertia_local, h->angular_velocity, dt);
		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
		quat dq = mul(spin, h->rotation);
		h->rotation.x += 0.5f * dt * dq.x; h->rotation.y += 0.5f * dt * dq.y; h->rotation.z += 0.5f * dt * dq.z; h->rotation.w += 0.5f * dt * dq.w;
		float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y + h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
		if (ql < 1e-15f) ql = 1.0f;
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql; h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
	}
}

// --- Pre-solve work function (parallel manifold setup) ---
typedef struct PreSolveCtx { WorldInternal* w; InternalManifold* manifolds; SolverManifold* sm; SolverContact* sc; PatchContact* pc; float dt; } PreSolveCtx;
static void pre_solve_work_fn(void* ctx, int start, int count)
{
	PreSolveCtx* ps = (PreSolveCtx*)ctx;
	for (int i = start; i < start + count; i++)
		pre_solve_manifold(ps->w, &ps->manifolds[i], i, ps->sm, ps->sc, ps->pc, ps->dt);
}

// Parallel pre_solve: alloc fixed-stride → dispatch → sequential warm start.
static void solver_pre_solve_dispatch(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, CK_DYNA PatchContact** out_pc, float dt, WorkFn work_fn)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	CK_DYNA PatchContact*   pc = NULL;
	afit(sm, manifold_count); asetlen(sm, manifold_count);
	int total_contacts = manifold_count * MAX_CONTACTS;
	afit(sc, total_contacts); asetlen(sc, total_contacts);
	afit(pc, total_contacts); asetlen(pc, total_contacts);
	memset(sm, 0, manifold_count * sizeof(SolverManifold));
	memset(sc, 0, total_contacts * sizeof(SolverContact));
	memset(pc, 0, total_contacts * sizeof(PatchContact));
	PreSolveCtx ps_ctx = { .w = w, .manifolds = manifolds, .sm = sm, .sc = sc, .pc = pc, .dt = dt };
	pool_dispatch(work_fn, &ps_ctx, manifold_count, 32, w->thread_count);
	// Warm start (sequential — modifies shared body velocities)
	int patch_warm = 1;
	for (int i = 0; i < manifold_count; i++) {
		SolverManifold* m = &sm[i];
		if (m->contact_count == 0) continue;
		BodyHot* a = &w->body_hot[m->body_a]; BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) { SolverContact* s = &sc[m->contact_start + ci]; if (patch_warm) { if (s->lambda_n == 0.0f) continue; apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n)); } else { if (s->lambda_n == 0.0f && s->lambda_t1 == 0.0f && s->lambda_t2 == 0.0f) continue; apply_impulse(a, b, s->r_a, s->r_b, add(add(scale(s->normal, s->lambda_n), scale(s->tangent1, s->lambda_t1)), scale(s->tangent2, s->lambda_t2))); } }
		if (patch_warm && (m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f)) apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2)));
		if (patch_warm && m->lambda_twist != 0.0f) { v3 tw = scale(m->normal, m->lambda_twist); a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, tw)); b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, tw)); }
	}
	*out_sm = sm; *out_sc = sc;
	if (out_pc) *out_pc = pc; else afree(pc);
}

// --- Narrowphase work function ---
typedef struct NP_WorkCtx
{
	WorldInternal* w;
	BroadPair* pairs;
	CK_DYNA InternalManifold** per_thread_manifolds; // array of per-thread dyn arrays
} NP_WorkCtx;

// Thread-local narrowphase: each block writes to its own manifold array.
// We use the block start as a pseudo-thread-id for output routing.
static void np_work_fn(void* ctx, int start, int count)
{
	NP_WorkCtx* np = (NP_WorkCtx*)ctx;
	// Use a simple thread-local manifold list keyed by block index.
	// Each block appends to a shared per-thread array using the current thread.
	// Since blocks are claimed atomically, we use thread-local storage.
	CK_DYNA InternalManifold* local = NULL;
	for (int i = start; i < start + count; i++) {
		BroadPair* p = &np->pairs[i];
		ShapeInternal* s0 = &np->w->body_cold[p->a].shapes[0];
		ShapeInternal* s1 = &np->w->body_cold[p->b].shapes[0];
		int ia = p->a, ib = p->b;
		if (s0->type > s1->type || (s0->type == s1->type && ia > ib)) { int tmp = ia; ia = ib; ib = tmp; s0 = &np->w->body_cold[ia].shapes[0]; s1 = &np->w->body_cold[ib].shapes[0]; }
		BodyHot* h0 = &np->w->body_hot[ia];
		BodyHot* h1 = &np->w->body_hot[ib];
		InternalManifold im = { .body_a = ia, .body_b = ib };
		int hit = 0;
		if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_SPHERE) hit = collide_sphere_sphere(make_sphere(h0, s0), make_sphere(h1, s1), &im.m);
		else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CAPSULE) hit = collide_sphere_capsule(make_sphere(h0, s0), make_capsule(h1, s1), &im.m);
		else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_BOX) hit = collide_sphere_box(make_sphere(h0, s0), make_box(h1, s1), &im.m);
		else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CAPSULE) hit = collide_capsule_capsule(make_capsule(h0, s0), make_capsule(h1, s1), &im.m);
		else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_BOX) hit = collide_capsule_box(make_capsule(h0, s0), make_box(h1, s1), &im.m);
		else if (s0->type == SHAPE_BOX && s1->type == SHAPE_BOX) hit = collide_box_box_ex(make_box(h0, s0), make_box(h1, s1), &im.m, &(int){-1});
		else if (s0->type == SHAPE_BOX && s1->type == SHAPE_HULL) hit = collide_hull_hull((ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents }, make_convex_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_HULL) hit = collide_sphere_hull(make_sphere(h0, s0), make_convex_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_HULL) hit = collide_capsule_hull(make_capsule(h0, s0), make_convex_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_HULL && s1->type == SHAPE_HULL) hit = collide_hull_hull(make_convex_hull(h0, s0), make_convex_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_SPHERE && s1->type == SHAPE_CYLINDER) hit = collide_sphere_hull(make_sphere(h0, s0), make_cylinder_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_CAPSULE && s1->type == SHAPE_CYLINDER) hit = collide_capsule_hull(make_capsule(h0, s0), make_cylinder_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_BOX && s1->type == SHAPE_CYLINDER) hit = collide_hull_hull((ConvexHull){ &s_unit_box_hull, h0->position, h0->rotation, s0->box.half_extents }, make_cylinder_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_HULL && s1->type == SHAPE_CYLINDER) hit = collide_hull_hull(make_convex_hull(h0, s0), make_cylinder_hull(h1, s1), &im.m);
		else if (s0->type == SHAPE_CYLINDER && s1->type == SHAPE_CYLINDER) hit = collide_hull_hull(make_cylinder_hull(h0, s0), make_cylinder_hull(h1, s1), &im.m);
		if (hit) apush(local, im);
	}
	// Merge local manifolds into per-thread output (we just push to the global array with a lock).
	// Simple: use a global lock for the merge since blocks are large and merges are infrequent.
	if (asize(local) > 0) {
		static volatile long np_merge_lock;
		while (_InterlockedCompareExchange(&np_merge_lock, 1, 0) != 0) _mm_pause();
		for (int k = 0; k < asize(local); k++) apush(*np->per_thread_manifolds, local[k]);
		_InterlockedExchange(&np_merge_lock, 0);
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

	double t0 = perf_now();
	warm_cache_age_and_evict(w);
	int n_workers = w->thread_count > 0 ? w->thread_count : 1;
#ifdef _WIN32
	if (n_workers > 1) pool_ensure(n_workers);
#endif
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
		for (int ii = 0; ii < asize(w->islands) && !any_awake; ii++)
			if ((w->island_gen[ii] & 1) && w->islands[ii].awake) any_awake = 1;
		for (int bi2 = 0; bi2 < body_count && !any_awake; bi2++)
			if (split_alive(w->body_gen, bi2) && w->body_hot[bi2].inv_mass > 0.0f && w->body_cold[bi2].island_id < 0) any_awake = 1;
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
	// Parallel narrowphase on collected pairs.
	if (n_workers > 1 && asize(np_pairs) >= 32) {
		CK_DYNA InternalManifold* merged = manifolds;
		NP_WorkCtx np_ctx = { .w = w, .pairs = np_pairs, .per_thread_manifolds = &merged };
		pool_dispatch(np_work_fn, &np_ctx, asize(np_pairs), 32, n_workers);
		manifolds = merged;
	} else if (asize(np_pairs) > 0) {
		for (int i = 0; i < asize(np_pairs); i++)
			narrowphase_pair(w, np_pairs[i].a, np_pairs[i].b, &manifolds);
	}
	afree(np_pairs);
#endif
	double t_iuc = perf_now();
	islands_update_contacts(w, manifolds, asize(manifolds));
	w->perf.pgs.pos_joints = perf_now() - t_iuc; // hijack: store islands_update_contacts time
	w->perf.broadphase = perf_now() - t1;

	aclear(w->debug_contacts);
	for (int i = 0; i < asize(manifolds); i++)
		for (int c = 0; c < manifolds[i].m.count; c++)
			apush(w->debug_contacts, manifolds[i].m.contacts[c]);

	int manifold_count = asize(manifolds);

	// Fast path: skip entire solver when nothing to solve (all sleeping, no contacts, no joints).
	// Only safe when no islands are awake — free-falling bodies with 0 contacts still need integration.
	int any_awake_island = 0;
	for (int ii = 0; ii < asize(w->islands) && !any_awake_island; ii++)
		if ((w->island_gen[ii] & 1) && w->islands[ii].awake) any_awake_island = 1;
	int any_active_joints = 0;
	for (int ji = 0; ji < asize(w->joints) && !any_active_joints; ji++)
		if (split_alive(w->joint_gen, ji)) any_active_joints = 1;
	// Also check for bodies with no island (newly created, free-falling)
	int any_unisland_dynamic = 0;
	for (int bi2 = 0; bi2 < body_count && !any_unisland_dynamic; bi2++)
		if (split_alive(w->body_gen, bi2) && w->body_hot[bi2].inv_mass > 0.0f && w->body_cold[bi2].island_id < 0) any_unisland_dynamic = 1;
	if (manifold_count == 0 && !any_awake_island && !any_active_joints && !any_unisland_dynamic) {
		w->perf.pre_solve = 0;
		w->perf.pgs_solve = 0;
		w->perf.position_correct = 0;
		w->perf.pgs = (PGSTimers){0};
		// No post-step bvh_refit needed — bodies didn't move (solver skipped).
		// The broadphase already refitted at start of this frame.
		double t4 = perf_now();
		islands_try_splits(w);
		if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
		w->perf.islands = perf_now() - t4;
		afree(manifolds);
		w->perf.total = perf_now() - t_total;
		return;
	}

	// --- Pre-solve (once per frame, using sub_dt for softness/bias) ---
	double t2 = perf_now();
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	CK_DYNA PatchContact* pc = NULL;
	// Pre-solve: parallel when threading enabled (each manifold writes to fixed-stride slots).
#ifdef _WIN32
	if (n_workers > 1 && manifold_count >= 64)
		solver_pre_solve_dispatch(w, manifolds, manifold_count, &sm, &sc, &pc, sub_dt, pre_solve_work_fn);
	else
#endif
		solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, &pc, sub_dt);

	SolverJoint* sol_joints = NULL;
	joints_pre_solve(w, sub_dt, &sol_joints);

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
	for (int i = 0; i < sm_count; i++) {
		if (sm[i].contact_count == 0) continue; // skip empty manifolds (static-static filtered out)
		ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
			.body_a = sm[i].body_a, .body_b = sm[i].body_b };
		apush(crefs, r);
	}
	if (!w->ldl_enabled) {
		// No LDL: all joints go into PGS
		for (int i = 0; i < asize(sol_joints); i++) {
			ConstraintRef r = { .type = CTYPE_JOINT, .index = i,
				.body_a = sol_joints[i].body_a, .body_b = sol_joints[i].body_b };
			apush(crefs, r);
		}
	} else {
		// LDL enabled: LDL handles all joints (bilateral DOFs).
		// Limit DOFs are solved separately after PGS via joints_solve_limits().
	}

	int cref_count = asize(crefs);
	int batch_starts[65] = {0};
	int color_count = 0;
	if (cref_count > 0)
		color_constraints(crefs, cref_count, count, batch_starts, &color_count);
	w->perf.pgs.graph_color = perf_now() - t_gc;

	w->perf.pre_solve = perf_now() - t2;

	// --- Sub-step loop ---
	// Unified path: PGS iterates all constraints (contacts + joints).
	// When LDL enabled, K is factored once at substep start, and a mid-loop
	// K^-1 residual correction is applied at the configured iteration.
	double t_pgs = 0, t_pos = 0, t_int_sub = 0;
	double t_jlim = 0, t_ldl = 0, t_relax = 0, t_posJ = 0;
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
		// Fast path: contacts with no joints use compact SolverBodyVel
		// (32 bytes/body instead of 120 bytes — fits more bodies in cache).
		int use_body_vel = (asize(sol_joints) == 0);
		if (use_body_vel) solver_sync_vel_in(w);

		double tp = perf_now();
		if (use_body_vel) {
#if SIMD_SSE
			// Build persistent SoA batch array ONCE — reused across all iterations.
			// Eliminates per-iteration manifold→SoA conversion (was 40x redundant).
			int total_batches = 0;
			for (int c = 0; c < color_count; c++) total_batches += (batch_starts[c+1] - batch_starts[c] + 3) / 4;
			CK_DYNA PGS_Batch4* batches = NULL;
			afit(batches, total_batches);
			CK_DYNA int* color_batch_starts = NULL;
			apush(color_batch_starts, 0);
			for (int c = 0; c < color_count; c++) {
				int start = batch_starts[c], end = batch_starts[c + 1];
				for (int i = start; i + 3 < end; i += 4) {
					int idx[4] = { crefs[i].index, crefs[i+1].index, crefs[i+2].index, crefs[i+3].index };
					PGS_Batch4 bt;
					pgs_batch4_prepare(&bt, sm, idx, 4, pc);
					apush(batches, bt);
				}
				// Remainder manifolds: build partial batch
				int rem_start = start + ((end - start) / 4) * 4;
				if (rem_start < end) {
					int idx[4] = {0};
					int rem_count = end - rem_start;
					for (int j = 0; j < rem_count; j++) idx[j] = crefs[rem_start + j].index;
					PGS_Batch4 bt;
					pgs_batch4_prepare(&bt, sm, idx, rem_count, pc);
					apush(batches, bt);
				}
				apush(color_batch_starts, asize(batches));
			}
			int batch_count = asize(batches);

			// Iteration loop: dispatch per color within each iteration.
			// Colors must be sequential (graph coloring guarantee), but within each color
			// all blocks are parallel.
			// n_workers already computed at top of world_step
			for (int iter = 0; iter < w->velocity_iters; iter++) {
				for (int c = 0; c < color_count; c++) {
					int bs = color_batch_starts[c], be = color_batch_starts[c + 1];
					int n_color_batches = be - bs;
#ifdef _WIN32
					if (n_workers > 1 && n_color_batches >= n_workers * 2) {
						PGS_WorkCtx pgs_ctx = { .bodies = w->body_vel, .batches = batches + bs };
						pool_dispatch(pgs_work_fn, &pgs_ctx, n_color_batches, SOLVER_BLOCK_SIZE, n_workers);
					} else
#endif
					{
						for (int bi = bs; bi < be; bi++)
							solve_contact_batch4_sv(w->body_vel, &batches[bi]);
					}
				}
				double tjl = perf_now();
				joints_solve_limits(w, sol_joints, asize(sol_joints));
				t_jlim += perf_now() - tjl;
			}

			// Scatter lambdas back from persistent batches to SolverManifold/PatchContact.
			for (int bi = 0; bi < batch_count; bi++) {
				PGS_Batch4* bt = &batches[bi];
				// Manifold lambdas
				for (int c2 = 0; c2 < color_count; c2++) {
					if (bi >= color_batch_starts[c2] && bi < color_batch_starts[c2+1]) {
						int base = batch_starts[c2] + (bi - color_batch_starts[c2]) * 4;
						for (int j = 0; j < 4 && base + j < batch_starts[c2+1]; j++) {
							int mi = crefs[base + j].index;
							sm[mi].lambda_t1 = ((float*)&bt->lambda_t1)[j];
							sm[mi].lambda_t2 = ((float*)&bt->lambda_t2)[j];
							sm[mi].lambda_twist = ((float*)&bt->lambda_twist)[j];
						}
						break;
					}
				}
				// Contact lambdas
				for (int cp2 = 0; cp2 < bt->max_contacts; cp2++) {
					float nl[4]; _mm_storeu_ps(nl, bt->cp[cp2].lambda_n);
					for (int j = 0; j < 4; j++) {
						if (bt->body_a[j] == 0 && bt->body_b[j] == 0 && j > 0) continue; // padding lane
						// Find manifold index for this lane
						for (int c2 = 0; c2 < color_count; c2++) {
							if (bi >= color_batch_starts[c2] && bi < color_batch_starts[c2+1]) {
								int base = batch_starts[c2] + (bi - color_batch_starts[c2]) * 4;
								if (base + j < batch_starts[c2+1] && cp2 < sm[crefs[base+j].index].contact_count)
									pc[sm[crefs[base+j].index].contact_start + cp2].lambda_n = nl[j];
								break;
							}
						}
					}
				}
			}
			afree(batches);
			afree(color_batch_starts);
#else
			for (int iter = 0; iter < w->velocity_iters; iter++) {
				for (int c = 0; c < color_count; c++)
					for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
						if (crefs[i].type == CTYPE_CONTACT)
							solve_contact_patch_sv(w->body_vel, &sm[crefs[i].index], pc);
				double tjl = perf_now();
				joints_solve_limits(w, sol_joints, asize(sol_joints));
				t_jlim += perf_now() - tjl;
			}
#endif
		} else {
			for (int iter = 0; iter < w->velocity_iters; iter++) {
				for (int c = 0; c < color_count; c++)
					for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
						solve_constraint(w, &crefs[i], sm, sc, sol_joints);
				double tjl = perf_now();
				joints_solve_limits(w, sol_joints, asize(sol_joints));
				t_jlim += perf_now() - tjl;
			}
		}
		t_pgs += perf_now() - tp;

		if (use_body_vel) {
			solver_sync_vel_out(w);
			// Sync lambda_n from compact PatchContact back to SolverContact
			for (int ci2 = 0; ci2 < asize(pc); ci2++) sc[ci2].lambda_n = pc[ci2].lambda_n;
		}

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
	islands_try_splits(w);
	if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
	w->perf.islands = perf_now() - t4;

	afree(manifolds);
	afree(pc);
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

	w->body_cold[idx] = (BodyCold){
		.mass = params.mass,
		.shapes = NULL,
		.bvh_leaf = -1,
		.island_id = -1,
		.island_prev = -1,
		.island_next = -1,
	};
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	float ang_damp = params.angular_damping;
	if (ang_damp == 0.0f) ang_damp = 0.03f; // default: 3%/s (BEPU-style)
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
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
		if (w->body_hot[idx].inv_mass == 0.0f) tree = w->bvh_static;
		else if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) tree = w->bvh_sleeping;
		else tree = w->bvh_dynamic;
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
	recompute_body_inertia(w, idx);

	// Insert into BVH on first shape add.
	if (w->broadphase_type == BROADPHASE_BVH && asize(w->body_cold[idx].shapes) == 1) {
		AABB box = aabb_expand(body_aabb(&w->body_hot[idx], &w->body_cold[idx]), BVH_AABB_MARGIN);
		BVH_Tree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		w->body_cold[idx].bvh_leaf = bvh_insert(tree, idx, box);
	}
}

v3 body_get_position(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].position;
}

quat body_get_rotation(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].rotation;
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
	w->body_hot[idx].velocity = vel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_angular_velocity(World world, Body body, v3 avel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].angular_velocity = avel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

int body_is_asleep(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	if (w->body_hot[idx].inv_mass == 0.0f) return 1; // static bodies are always "asleep"
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

void body_set_sleep_allowed(World world, Body body, int allowed)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].sleep_allowed = allowed;
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
	int idx;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	// Grow joint arrays manually (no split_add -- joints don't need hot/cold split)
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1); // odd = alive
	}

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

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	// Auto-compute rest length if not specified
	float rest = params.rest_length;
	if (rest <= 0.0f) {
		BodyHot* a = &w->body_hot[ba];
		BodyHot* b = &w->body_hot[bb];
		v3 wa = add(a->position, rotate(a->rotation, params.local_offset_a));
		v3 wb = add(b->position, rotate(b->rotation, params.local_offset_b));
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

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	v3 axis_a_local = norm(params.local_axis_a);
	v3 axis_b_local = norm(params.local_axis_b);

	// Compute reference directions perpendicular to the hinge axis for angle measurement.
	// local_ref_a: arbitrary unit vector perpendicular to axis_a in body A's local space.
	// local_ref_b: chosen so measured angle = 0 at the initial configuration.
	v3 ref_a, ref_a_t2;
	hinge_tangent_basis(axis_a_local, &ref_a, &ref_a_t2);
	// Transform ref_a into world, then into body B's local space
	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
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

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
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

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	quat q_a = w->body_hot[ba].rotation;
	quat q_b = w->body_hot[bb].rotation;
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

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	unlink_joint_from_island(w, idx);
	memset(&w->joints[idx], 0, sizeof(JointInternal));
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
		BVHChild* c = bvh_child(n, s);
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
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];
		JointDebugInfo info = {0};
		info.type = j->type;
		if (j->type == JOINT_BALL_SOCKET) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->ball_socket.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->ball_socket.local_b));
			info.is_soft = j->ball_socket.spring.frequency > 0;
		} else if (j->type == JOINT_DISTANCE) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->distance.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->distance.local_b));
			info.is_soft = j->distance.spring.frequency > 0;
		} else if (j->type == JOINT_HINGE) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->hinge.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->hinge.local_b));
			info.axis_a = norm(rotate(a->rotation, j->hinge.local_axis_a));
			info.is_soft = j->hinge.spring.frequency > 0;
			info.motor_speed = j->hinge.motor_speed;
			info.motor_max_impulse = j->hinge.motor_max_impulse;
			info.limit_min = j->hinge.limit_min;
			info.limit_max = j->hinge.limit_max;
			info.ref_a = rotate(a->rotation, j->hinge.local_ref_a);
			info.ref_b = rotate(b->rotation, j->hinge.local_ref_b);
			float angle = atan2f(dot(cross(info.ref_a, info.ref_b), info.axis_a), dot(info.ref_a, info.ref_b));
			info.current_angle = angle;
			info.limit_active = (j->hinge.limit_min != 0 && angle <= j->hinge.limit_min) || (j->hinge.limit_max != 0 && angle >= j->hinge.limit_max);
		} else if (j->type == JOINT_FIXED) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->fixed.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->fixed.local_b));
			info.is_soft = j->fixed.spring.frequency > 0;
		} else if (j->type == JOINT_PRISMATIC) {
			info.anchor_a = add(a->position, rotate(a->rotation, j->prismatic.local_a));
			info.anchor_b = add(b->position, rotate(b->rotation, j->prismatic.local_b));
			info.axis_a = norm(rotate(a->rotation, j->prismatic.local_axis_a));
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
		AABB b = body_aabb(&w->body_hot[idx], &w->body_cold[idx]);
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
