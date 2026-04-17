// See LICENSE for licensing info.
// test_main.c -- unity build root for test executable.
// Includes only the physics engine + tests, no rendering/SDL/imgui.

#define CKIT_IMPLEMENTATION
#include "ckit.h"

#include "nudge.h"
#include "nudge.c"
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

int main(int argc, char* argv[])
{
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
	int bench_pyramid_base = 0;
	int bench_pile_grid = 10;
	int bench_pile_height = 5;
	int bench_pile_frames = 300;
	int bench_mixed = 0;
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

	if (bench_chaos) {
		WorldParams wp = { .gravity = V3(0, -9.81f, 0), .broadphase = BROADPHASE_BVH, .sub_steps = sub_steps, .velocity_iters = vel_iters, .contact_hertz = hertz, .contact_damping_ratio = damping };
		bench_hull_chaos(chaos_bodies, chaos_frames, chaos_churn, wp);
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
	return test_fail > 0 ? 1 : 0;
}
