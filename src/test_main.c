// See LICENSE for licensing info.
// test_main.c -- unity build root for test executable.
// Includes only the physics engine + tests, no rendering/SDL/imgui.

#include <stdbool.h>
#ifdef _WIN32
	#define WIN32_LEAN_AND_MEAN
	#include <winsock2.h>
	#include <ws2tcpip.h>
	#pragma comment(lib, "ws2_32.lib")
	#undef small
	#undef near
	#undef far
#endif

#define CKIT_IMPLEMENTATION
#include "ckit.h"

#include "nudge.h"
#include "nudge.c"
// debug_server.c provides a tiny TCP server + reflection + the DBG_BREAK macro.
// In test mode (NUDGE_HOST_APP not defined) only inspection + break control are
// available. Tests opt in via --debug; the DBG_BREAK macro is a no-op otherwise.
#include "debug_server.c"
#include "tests.c"
#include "tests_ldl_unit.c"
#include "tests_inertia_unit.c"
#include "tests_jacobian_unit.c"
#include "tests_k_matrix_unit.c"
#include "tests_impulse_unit.c"
#include "tests_misc_ldl_unit.c"
#include "tests_sparse_unit.c"
#include "tests_ordering_unit.c"
#include "tests_dfs_reorder_unit.c"
#include "tests_bundles_unit.c"
#include "tests_integration_ldl.c"
#include "tests_spring_unit.c"
#include "tests_shattering_unit.c"
#include "tests_pgs_vs_ldl.c"
#include "tests_gjk_perf.c"
#include "tests_arena_unit.c"
#include "tests_threading_unit.c"
#include "tests_weld_bridge_unit.c"
#include "tests_multishape_unit.c"
#include "tests_contacts_unit.c"
#include "tests_heightfield_unit.c"
#include "tests_epa_debug.c"
#include "tests_epa_perf.c"
#include "tests_determinism.c"

int main(int argc, char* argv[])
{
	setvbuf(stdout, NULL, _IONBF, 0); // unbuffered so crashes don't hide output

	// --debug starts the TCP debug server so a viewer can attach. --break=<pat>
	// arms DBG_BREAK(name, world) sites whose names glob-match the pattern.
	// Pattern syntax: "*" matches all, "ldl_*" prefix, "*cant*" contains, "foo" exact.
	// Without --debug, DBG_BREAK is a no-op (predicate check only, ~free).
	const char* break_pat = NULL;
	int debug_enabled = 0;
	int debug_port = 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--debug") == 0) debug_enabled = 1;
		else if (strncmp(argv[i], "--break=", 8) == 0) { break_pat = argv[i] + 8; debug_enabled = 1; }
		else if (strncmp(argv[i], "--debug-port=", 13) == 0) debug_port = atoi(argv[i] + 13);
	}
	if (debug_enabled) {
		if (debug_port > 0) debug_server_set_port(debug_port);
		debug_server_set_break_filter(break_pat ? break_pat : "*");
		g_dbg_break_enabled = 1;
		debug_server_init();
		fprintf(stderr, "[dbg] test-mode debug server up. break filter=\"%s\". Connect with nudge_viewer.exe localhost %d.\n",
			break_pat ? break_pat : "*", debug_port > 0 ? debug_port : 9999);
	}

	for (int i = 1; i < argc; i++) if (strcmp(argv[i], "--debug-epa") == 0) { debug_epa(); return 0; }
	for (int i = 1; i < argc; i++) if (strcmp(argv[i], "--bench-epa") == 0) { bench_epa_vs_sat(); return 0; }
	for (int i = 1; i < argc; i++) if (strcmp(argv[i], "--bench-epa-scenes") == 0) { bench_epa_scenes(); return 0; }
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--determinism") == 0) {
			int det_threads = 1;
			for (int j = 1; j < argc; j++) {
				if (strcmp(argv[j], "--threads") == 0 && j + 1 < argc) det_threads = atoi(argv[j + 1]);
			}
			return run_determinism_test(det_threads);
		}
	}
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--fuzz-epa") == 0 && i + 1 < argc) {
			int n = atoi(argv[i + 1]);
			test_pass = 0; test_fail = 0;
			printf("--- epa fuzz (%d iterations per pair) ---\n", n);
			test_epa_fuzz(n);
			printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
			return test_fail > 0 ? 1 : 0;
		}
	}
	int fuzz_iters = 0;
	int soak = 0;
	int bench_stack = 0;
	int bench_pile = 0;
	int bench_suite_flag = 0;
	int bench_chaos = 0;
	int bench_qh = 0;
	int bench_ldl = 0;
	int ldl_chains = 10;
	int ldl_chain_len = 20;
	int ldl_frames = 200;
	int chaos_bodies = 500;
	int chaos_frames = 30;
	int chaos_churn = 10;
	int bench_softbody = 0;
	int bench_pyramid_base = 0;
	int bench_pile_grid = 10;
	int bench_pile_height = 5;
	int bench_pile_frames = 300;
	int bench_mixed = 0;
	int bench_ragdoll_flag = 0;
	int bench_incr_np = 0;
	int bench_planes = 0;
	int pyramid_test = 0;
	int pyramid_base = 5;
	int pyramid_frames = 600;
	int sub_steps = 0;
	int vel_iters = 0;
	int threads = 0;
	float hertz = 0;
	float damping = 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--fuzz") == 0 && i + 1 < argc)
			fuzz_iters = atoi(argv[++i]);
		else if (strcmp(argv[i], "--soak") == 0)
			soak = 1;
		else if (strcmp(argv[i], "--bench-stack") == 0)
			bench_stack = (i + 1 < argc && argv[i+1][0] != '-') ? atoi(argv[++i]) : 10;
		else if (strcmp(argv[i], "--bench-pile") == 0)
			bench_pile = 1;
		else if (strcmp(argv[i], "--bench-pyramid") == 0)
			bench_pyramid_base = (i + 1 < argc && argv[i+1][0] != '-') ? atoi(argv[++i]) : 20;
		else if (strcmp(argv[i], "--bench-suite") == 0)
			bench_suite_flag = 1;
		else if (strcmp(argv[i], "--bench-chaos") == 0)
			bench_chaos = 1;
		else if (strcmp(argv[i], "--bench-incr-np") == 0)
			bench_incr_np = 1;
		else if (strcmp(argv[i], "--bench-planes") == 0)
			bench_planes = 1;
		else if (strcmp(argv[i], "--bench-softbody") == 0)
			bench_softbody = (i + 1 < argc && argv[i+1][0] != '-') ? atoi(argv[++i]) : 1;
		else if (strcmp(argv[i], "--chaos-bodies") == 0 && i + 1 < argc)
			chaos_bodies = atoi(argv[++i]);
		else if (strcmp(argv[i], "--chaos-frames") == 0 && i + 1 < argc)
			chaos_frames = atoi(argv[++i]);
		else if (strcmp(argv[i], "--chaos-churn") == 0 && i + 1 < argc)
			chaos_churn = atoi(argv[++i]);
		else if (strcmp(argv[i], "--pile-grid") == 0 && i + 1 < argc)
			bench_pile_grid = atoi(argv[++i]);
		else if (strcmp(argv[i], "--pile-height") == 0 && i + 1 < argc)
			bench_pile_height = atoi(argv[++i]);
		else if (strcmp(argv[i], "--pile-frames") == 0 && i + 1 < argc)
			bench_pile_frames = atoi(argv[++i]);
		else if (strcmp(argv[i], "--sub-steps") == 0 && i + 1 < argc)
			sub_steps = atoi(argv[++i]);
		else if (strcmp(argv[i], "--vel-iters") == 0 && i + 1 < argc)
			vel_iters = atoi(argv[++i]);
		else if (strcmp(argv[i], "--hertz") == 0 && i + 1 < argc)
			hertz = (float)atof(argv[++i]);
		else if (strcmp(argv[i], "--damping") == 0 && i + 1 < argc)
			damping = (float)atof(argv[++i]);
		else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc)
			threads = atoi(argv[++i]);
		else if (strcmp(argv[i], "--pyramid") == 0)
			pyramid_test = 1;
		else if (strcmp(argv[i], "--pyramid-base") == 0 && i + 1 < argc)
			pyramid_base = atoi(argv[++i]);
		else if (strcmp(argv[i], "--pyramid-frames") == 0 && i + 1 < argc)
			pyramid_frames = atoi(argv[++i]);
		else if (strcmp(argv[i], "--bench-qh") == 0)
			bench_qh = 1;
		else if (strcmp(argv[i], "--bench-mixed") == 0)
			bench_mixed = 1;
		else if (strcmp(argv[i], "--bench-ragdoll") == 0)
			bench_ragdoll_flag = (i + 1 < argc && argv[i+1][0] != '-') ? atoi(argv[++i]) : 20;
		else if (strcmp(argv[i], "--bench-ldl") == 0)
			bench_ldl = 1;
		else if (strcmp(argv[i], "--ldl-chains") == 0 && i + 1 < argc)
			ldl_chains = atoi(argv[++i]);
		else if (strcmp(argv[i], "--ldl-chain-len") == 0 && i + 1 < argc)
			ldl_chain_len = atoi(argv[++i]);
		else if (strcmp(argv[i], "--ldl-frames") == 0 && i + 1 < argc)
			ldl_frames = atoi(argv[++i]);
	}

	extern int bench_thread_count;
	bench_thread_count = threads;

	if (pyramid_test) {
		test_pyramid_2d_jiggle(pyramid_base, pyramid_frames);
		printf("\n");
		test_pyramid_yank(pyramid_base, pyramid_frames);
		return 0;
	}

	if (bench_qh) {
		bench_quickhull();
		bench_quickhull_10k();
		return 0;
	}

	if (bench_ldl) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_ldl_joints(ldl_chains, ldl_chain_len, ldl_frames, wp);
		return 0;
	}

	if (bench_mixed) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_mixed_contacts_joints(8, 8, 20, 6, 300, wp);
		return 0;
	}

	if (bench_ragdoll_flag) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_ragdoll(bench_ragdoll_flag, 300, wp);
		return 0;
	}

	if (bench_chaos) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_hull_chaos(chaos_bodies, chaos_frames, chaos_churn, wp);
		return 0;
	}

	if (bench_softbody) {
		// Deterministic soft-body diagnostic. Prints per-frame positions,
		// velocities, link lengths, and total kinetic energy so we can see
		// exactly when collision fires, correction magnitudes, and whether
		// constraints hold.
		// bench_softbody=1 -> 3-node rigid chain (simplest case)
		// bench_softbody=2 -> 13-node icosahedron ball (what the demo scene does)
		World w = create_world((WorldParams){
			.gravity = V3(0, -9.81f, 0),
			.broadphase = BROADPHASE_BVH,
			.sub_steps = 4,
		});
		Body floor = create_body(w, (BodyParams){ .position = V3(0, -1, 0), .rotation = quat_identity(), .mass = 0 });
		body_add_shape(w, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = V3(10, 1, 10) });

		SoftBody sb = create_soft_body(w, (SoftBodyParams){
			.default_spring = { 0.0f, 0.0f },
			.node_radius = 0.1f,
			.linear_damping = 0.02f,
		});

		if (bench_softbody == 1) {
			// 3-node vertical chain, spacing 0.3, top at y=3.
			int a = soft_body_add_node(w, sb, V3(0.0f, 3.0f, 0.0f), 1.0f);
			int b = soft_body_add_node(w, sb, V3(0.0f, 2.7f, 0.0f), 1.0f);
			int c = soft_body_add_node(w, sb, V3(0.0f, 2.4f, 0.0f), 1.0f);
			soft_body_add_link(w, sb, a, b, -1.0f, (SpringParams){0});
			soft_body_add_link(w, sb, b, c, -1.0f, (SpringParams){0});
			(void)a; (void)c;
		} else {
			// 13-node icosahedron + center ball.
			float phi_f = (1.0f + sqrtf(5.0f)) * 0.5f;
			float nrm = 1.0f / sqrtf(1.0f + phi_f * phi_f);
			float aa = phi_f * nrm, bb = nrm;
			v3 ico[12] = {
				V3(-bb,  aa, 0), V3( bb,  aa, 0), V3(-bb, -aa, 0), V3( bb, -aa, 0),
				V3( 0, -bb,  aa), V3( 0,  bb,  aa), V3( 0, -bb, -aa), V3( 0,  bb, -aa),
				V3( aa, 0, -bb), V3( aa, 0,  bb), V3(-aa, 0, -bb), V3(-aa, 0,  bb),
			};
			int tris[20][3] = {
				{ 0,11, 5}, { 0, 5, 1}, { 0, 1, 7}, { 0, 7,10}, { 0,10,11},
				{ 1, 5, 9}, { 5,11, 4}, {11,10, 2}, {10, 7, 6}, { 7, 1, 8},
				{ 3, 9, 4}, { 3, 4, 2}, { 3, 2, 6}, { 3, 6, 8}, { 3, 8, 9},
				{ 4, 9, 5}, { 2, 4,11}, { 6, 2,10}, { 8, 6, 7}, { 9, 8, 1},
			};
			v3 center = V3(0, 3, 0);
			float R = 0.55f;
			int surf[12];
			for (int i = 0; i < 12; i++) surf[i] = soft_body_add_node(w, sb, add(center, scale(ico[i], R)), 0.1f);
			int ctr = soft_body_add_node(w, sb, center, 0.4f);
			int have[12][12] = {{0}};
			for (int t = 0; t < 20; t++) {
				int ti[3] = { tris[t][0], tris[t][1], tris[t][2] };
				for (int k = 0; k < 3; k++) {
					int u = ti[k], vv = ti[(k+1)%3];
					int lo = u<vv?u:vv, hi = u<vv?vv:u;
					if (have[lo][hi]) continue;
					have[lo][hi] = 1;
					soft_body_add_link(w, sb, surf[lo], surf[hi], -1.0f, (SpringParams){0});
				}
			}
			for (int i = 0; i < 12; i++) soft_body_add_link(w, sb, ctr, surf[i], -1.0f, (SpringParams){0});
		}

		soft_body_build(w, sb);

		int frames = 180;
		int N = soft_body_node_count(w, sb);
		int L = soft_body_link_count(w, sb);
		extern double g_soft_body_max_lambda;
		extern double g_soft_body_max_rhs;
		extern double g_soft_body_min_K_diag;
		extern int g_soft_body_trace;
		g_soft_body_trace = 1;
		printf("[softbody] N=%d L=%d\n", N, L);
		printf("frame | min_y | max_y | max|v| | max|lam| | max|rhs| | minD\n");
		for (int f = 0; f < frames; f++) {
			world_step(w, 1.0f / 60.0f);
			const v3* pos = soft_body_node_positions(w, sb);
			const v3* vel = soft_body_node_velocities(w, sb);
			float min_y = 1e9f, max_y = -1e9f;
			float max_v = 0;
			for (int n = 0; n < N; n++) {
				if (pos[n].y < min_y) min_y = pos[n].y;
				if (pos[n].y > max_y) max_y = pos[n].y;
				float vm = sqrtf(vel[n].x*vel[n].x + vel[n].y*vel[n].y + vel[n].z*vel[n].z);
				if (vm > max_v) max_v = vm;
			}
			printf("%5d | %7.3f | %7.3f | %7.2f | %8.2f | %8.2f | %6.3f\n",
				f, min_y, max_y, max_v,
				g_soft_body_max_lambda, g_soft_body_max_rhs, g_soft_body_min_K_diag);
		}
		destroy_world(w);
		return 0;
	}

	if (bench_suite_flag) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_suite(wp);
		return 0;
	}

	if (bench_pile) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_box_pile(bench_pile_grid, bench_pile_height, bench_pile_frames, wp);
		return 0;
	}
	if (bench_pyramid_base > 0) {
		bench_pyramid(bench_pyramid_base, 600);
		return 0;
	}
	if (bench_planes) {
		bench_plane_compute();
		return 0;
	}
	if (bench_incr_np) {
		bench_incremental_np(20, 8, 600);
		return 0;
	}
	// --bench-trimesh: stress test on 25x25 terrain mesh with 49 mixed bodies,
	// runs twice (SIMD on + off) and prints narrowphase ms per frame.
	for (int ai = 1; ai < argc; ai++) {
		if (strcmp(argv[ai], "--bench-trimesh") != 0) continue;
		bench_trimesh_stress();
		return 0;
	}

	if (bench_stack > 0) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_box_stack_ex(bench_stack, wp);
		return 0;
	}

	if (soak) {
		test_pass = 0;
		test_fail = 0;
		test_quickhull_soak();
		return test_fail > 0 ? 1 : 0;
	}

	if (fuzz_iters > 0) {
		test_pass = 0;
		test_fail = 0;
		printf("--- quickhull fuzz (%d iterations per shape) ---\n", fuzz_iters);
		test_quickhull_fuzz(fuzz_iters);
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--unit") == 0) {
		test_pass = 0;
		test_fail = 0;
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		test_aalign();
		run_arena_unit_tests();
		run_threading_unit_tests();
		run_weld_bridge_unit_tests();
		run_multishape_unit_tests();
		run_contacts_unit_tests();
		run_heightfield_unit_tests();
		run_ldl_unit_tests();
		run_inertia_unit_tests();
		run_jacobian_unit_tests();
		run_k_matrix_unit_tests();
		run_impulse_unit_tests();
		run_misc_ldl_unit_tests();
		run_sparse_unit_tests();
		run_ordering_unit_tests();
		run_dfs_reorder_unit_tests();
		run_bundles_unit_tests();
		run_integration_ldl_tests();
		run_spring_unit_tests();
		run_shattering_unit_tests();
		run_pgs_vs_ldl_tests();
		test_aalign();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--ldl-unit") == 0) {
		test_pass = 0;
		test_fail = 0;
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		run_ldl_unit_tests();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--replay") == 0) {
		// Replay a mouse_recording.bin captured with F5 in the app.
		// Usage: nudge_tests --replay [file.bin]
		const char* path = argc > 2 ? argv[2] : "mouse_recording.bin";
		test_pass = 0; test_fail = 0;
		test_replay_recording(path);
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--run") == 0 && argc > 2) {
		test_pass = 0; test_fail = 0;
		for (int i = 2; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) { test_bail = 1; continue; }
			// Call test by name (add entries as needed)
			if (strcmp(argv[i], "box-wall") == 0) test_box_wall_explosion();
			else if (strcmp(argv[i], "pull-down") == 0) test_ldl_pull_down_heavy_chain();
			else if (strcmp(argv[i], "trace") == 0) test_stretched_joint_trace();
			else if (strcmp(argv[i], "showcase-stretch") == 0) test_ldl_showcase_chain_stretch();
			else if (strcmp(argv[i], "lift-drop") == 0) test_ldl_lift_and_drop();
			else if (strcmp(argv[i], "heavy-chain") == 0) test_ldl_heavy_chain();
			else if (strcmp(argv[i], "mouse-yank") == 0) test_ldl_mouse_yank_chain();
			else printf("unknown test: %s\n", argv[i]);
		}
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--quick") == 0) {
		test_ldl_stress_single_constraint();
		test_ldl_heavy_chain();
		test_ldl_two_independent_chains();
		test_ldl_hub_star_shattering();
		test_ldl_topology_change();
		test_ldl_stress_dense_clique();
		test_ldl_stress_alternating_mass();
		test_ldl_mixed_chain_and_hub();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--gjk-perf") == 0) {
		test_pass = 0; test_fail = 0;
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		run_gjk_perf_tests();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--slow") == 0) {
		test_pass = 0; test_fail = 0;
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		run_solver_tests_slow();
		run_ldl_stress_tests();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else {
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
			if (strcmp(argv[i], "--timing") == 0) test_timing = 1;
		}
		run_tests();
		run_arena_unit_tests();
		run_threading_unit_tests();
		run_weld_bridge_unit_tests();
		run_multishape_unit_tests();
		run_contacts_unit_tests();
		run_heightfield_unit_tests();
		run_ldl_unit_tests();
		run_inertia_unit_tests();
		run_jacobian_unit_tests();
		run_k_matrix_unit_tests();
		run_impulse_unit_tests();
		run_misc_ldl_unit_tests();
		run_sparse_unit_tests();
		run_ordering_unit_tests();
		run_dfs_reorder_unit_tests();
		run_bundles_unit_tests();
		run_integration_ldl_tests();
		run_spring_unit_tests();
		run_shattering_unit_tests();
		run_pgs_vs_ldl_tests();
	}
	if (debug_enabled) debug_server_shutdown();
	return test_fail > 0 ? 1 : 0;
}
