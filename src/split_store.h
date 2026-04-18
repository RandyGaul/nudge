// See LICENSE for licensing info.
#ifndef SPLIT_STORE_H
#define SPLIT_STORE_H

#include <stdint.h>
#include <string.h>

// split_store.h -- Paired hot/cold storage with stable generational handles.
//
// PURPOSE: Cache-friendly data layout for tight iteration loops (solvers,
// broadphase, rendering). Inspired by Box2D's hot/cold body split.
//
// The idea: a fat struct with everything (position, velocity, mass, shapes,
// user data, topology links) destroys cache when the solver only needs
// position + velocity. Instead, split into two parallel arrays:
//
//   hot[]  -- only the fields the inner loop touches every frame (position,
//             velocity, inv_mass). Small struct, contiguous, cache-line friendly.
//   cold[] -- everything else (mass, shapes, topology, user data). Same indices
//             as hot[], but only accessed on-demand (creation, queries, etc).
//
// Both arrays share indices, so hot[i] and cold[i] describe the same object.
// A separate gen[] array tracks slot liveness (odd = alive, even = dead) and
// enables stable generational handles without polluting either data struct.
// A freelist[] recycles slots so indices remain stable across add/remove.
//
// CACHE OPTIMIZATION GUIDE:
//   - Put frequently-iterated solver fields in Hot (velocity, position, forces).
//   - Put rarely-read metadata in Cold (mass, shapes, user data, topology).
//   - Keep Hot structs small -- ideally <= 64 bytes (one cache line).
//   - The solver loop iterates hot[] linearly. Prefetch is automatic.
//   - Cold is only dereferenced when you need it (handle lookup, shape query).
//   - Gaps from deleted slots are zeroed and skipped via split_alive().
//     The branch is predictable and avoids the indirection-fixup cost of
//     swap-remove. Compact later if occupancy drops.
//
// Usage:
//     CK_DYNA MyCold*    cold = NULL;
//     CK_DYNA MyHot*     hot  = NULL;
//     CK_DYNA uint32_t*  gen  = NULL;
//     CK_DYNA int*       free_list = NULL;
//
//     int idx;
//     split_add(cold, hot, gen, free_list, idx);
//     cold[idx] = (MyCold){ ... };
//     hot[idx]  = (MyHot){ ... };
//     MyHandle h = split_handle(MyHandle, gen, idx);
//
//     assert(split_valid(gen, h));
//
//     for (int i = 0; i < asize(hot); i++) {
//         if (!split_alive(gen, i)) continue;
//         // iterate hot[i] ...
//     }
//
//     split_del(cold, hot, gen, free_list, idx);
//     split_free(cold, hot, gen, free_list);

// Allocate a slot. Reuses from freelist or grows all arrays.
// Sets out_idx to the allocated index. gen[out_idx] becomes odd (alive).
#define split_add(cold, hot, gen, freelist, out_idx) do { \
	if (asize(freelist) > 0) { \
		(out_idx) = apop(freelist); \
	} else { \
		(out_idx) = asize(cold); \
		afit((cold), (out_idx) + 1); \
		afit((hot), (out_idx) + 1); \
		afit((gen), (out_idx) + 1); \
		asetlen((cold), (out_idx) + 1); \
		asetlen((hot), (out_idx) + 1); \
		asetlen((gen), (out_idx) + 1); \
		(gen)[(out_idx)] = 0; \
	} \
	memset(&(cold)[(out_idx)], 0, sizeof(*(cold))); \
	memset(&(hot)[(out_idx)], 0, sizeof(*(hot))); \
	(gen)[(out_idx)]++; \
} while(0)

// Free a slot. Zeros cold+hot, bumps gen to even (dead), pushes to freelist.
#define split_del(cold, hot, gen, freelist, i) do { \
	memset(&(cold)[(i)], 0, sizeof(*(cold))); \
	memset(&(hot)[(i)], 0, sizeof(*(hot))); \
	(gen)[(i)]++; \
	apush((freelist), (i)); \
} while(0)

// Check if slot i is alive (odd generation).
#define split_alive(gen, i) ((gen)[(i)] & 1)

// Free all arrays.
#define split_free(cold, hot, gen, freelist) do { \
	afree(cold); afree(hot); afree(gen); afree(freelist); \
} while(0)

// Ensure an extra parallel array has capacity for index idx. Grows and zeros new slots.
// Use after split_add to keep additional arrays (e.g. body_vel) in sync:
//   split_add(cold, hot, gen, free, idx);
//   split_ensure(vel, idx);
#define split_ensure(arr, idx) do { \
	int _need = (idx) + 1; \
	if (asize(arr) < _need) { \
		int _old = asize(arr); \
		afit((arr), _need); \
		asetlen((arr), _need); \
		memset(&(arr)[_old], 0, (_need - _old) * sizeof(*(arr))); \
	} \
} while(0)

// Zero a slot of an extra parallel array (call alongside split_del).
#define split_clear(arr, idx) memset(&(arr)[(idx)], 0, sizeof(*(arr)))

// Handle encoding: low 32 bits = index, high 32 bits = generation.
#define handle_make(idx, gen)  ((uint64_t)(gen) << 32 | (uint32_t)(idx))
#define handle_index(h)        ((int)((h).id & 0xFFFFFFFF))
#define handle_gen(h)          ((uint32_t)((h).id >> 32))

// Build a handle from a split store's gen array and an index.
#define split_handle(T, gen, idx) ((T){ handle_make((idx), (gen)[(idx)]) })

// Validate a handle against the current gen array.
#define split_valid(gen, h) ((gen)[handle_index(h)] == handle_gen(h))

// --- 3-tier split store (hot + state + cold) ---
// Same as split_add/del/free but manages a third parallel "state" array.
// Use when you need two hot-path tiers (e.g. solver-hot + integration-warm)
// plus a cold tier for lifecycle/topology data.
//
// Example (bodies):
//   split3_add(body_cold, body_state, body_hot, body_gen, body_free, idx);
//   split3_del(body_cold, body_state, body_hot, body_gen, body_free, idx);
//   split3_free(body_cold, body_state, body_hot, body_gen, body_free);

#define split3_add(cold, state, hot, gen, freelist, out_idx) do { \
	split_add(cold, hot, gen, freelist, out_idx); \
	split_ensure(state, out_idx); \
	memset(&(state)[(out_idx)], 0, sizeof(*(state))); \
} while(0)

#define split3_del(cold, state, hot, gen, freelist, i) do { \
	split_del(cold, hot, gen, freelist, i); \
	split_clear(state, i); \
} while(0)

#define split3_free(cold, state, hot, gen, freelist) do { \
	split_free(cold, hot, gen, freelist); afree(state); \
} while(0)

#endif
