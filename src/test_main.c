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

int main(int argc, char* argv[])
{
	int fuzz_iters = 0;
	int soak = 0;
	int bench_stack = 0;
	int sub_steps = 0;
	int vel_iters = 0;
	float hertz = 0;
	float damping = 0;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--fuzz") == 0 && i + 1 < argc)
			fuzz_iters = atoi(argv[++i]);
		else if (strcmp(argv[i], "--soak") == 0)
			soak = 1;
		else if (strcmp(argv[i], "--bench-stack") == 0)
			bench_stack = (i + 1 < argc && argv[i+1][0] != '-') ? atoi(argv[++i]) : 10;
		else if (strcmp(argv[i], "--sub-steps") == 0 && i + 1 < argc)
			sub_steps = atoi(argv[++i]);
		else if (strcmp(argv[i], "--vel-iters") == 0 && i + 1 < argc)
			vel_iters = atoi(argv[++i]);
		else if (strcmp(argv[i], "--hertz") == 0 && i + 1 < argc)
			hertz = (float)atof(argv[++i]);
		else if (strcmp(argv[i], "--damping") == 0 && i + 1 < argc)
			damping = (float)atof(argv[++i]);
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
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--ldl-unit") == 0) {
		test_pass = 0;
		test_fail = 0;
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		run_ldl_unit_tests();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else if (argc > 1 && strcmp(argv[1], "--quick") == 0) {
		test_ldl_stress_single_constraint();
		test_ldl_heavy_chain();
		test_ldl_two_independent_chains();
		test_ldl_hub_star_shattering();
		test_ldl_topology_change();
		test_ldl_energy_comprehensive();
		test_ldl_long_chain();
		test_ldl_stress_dense_clique();
		test_ldl_stress_alternating_mass();
		test_ldl_mixed_chain_and_hub();
		test_ldl_stress_stretched_recovery();
		test_ldl_stress_heavy_stretched_recovery();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else {
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
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
	}
	return test_fail > 0 ? 1 : 0;
}
