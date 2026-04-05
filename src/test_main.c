// See LICENSE for licensing info.
// test_main.c -- unity build root for test executable.
// Includes only the physics engine + tests, no rendering/SDL/imgui.

#define CKIT_IMPLEMENTATION
#include "ckit.h"

#include "nudge.h"
#include "nudge.c"
#include "tests.c"

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
	} else if (argc > 1 && strcmp(argv[1], "--quick") == 0) {
		// Diagnostic: 3-link chain with 100:1 on last link
		for (int use_ldl = 0; use_ldl <= 1; use_ldl++) {
			World w = create_world((WorldParams){ .gravity = V3(0, -9.81f, 0) });
			WorldInternal* wi = (WorldInternal*)w.id;
			wi->ldl_enabled = use_ldl;
			float ll = 0.8f;
			v3 oa = V3(ll*0.5f,0,0), ob = V3(-ll*0.5f,0,0);
			Body anchor = create_body(w, (BodyParams){ .position = V3(0, 10, 0), .rotation = quat_identity(), .mass = 0 });
			body_add_shape(w, anchor, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.1f });
			Body prev = anchor;
			Body chain[3];
			for (int i = 0; i < 3; i++) {
				float mass = (i == 2) ? 100.0f : 1.0f;
				chain[i] = create_body(w, (BodyParams){ .position = V3((i+1)*ll, 10, 0), .rotation = quat_identity(), .mass = mass });
				body_add_shape(w, chain[i], (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.15f });
				create_ball_socket(w, (BallSocketParams){ .body_a = prev, .body_b = chain[i], .local_offset_a = oa, .local_offset_b = ob });
				prev = chain[i];
			}
			for (int f = 0; f < 120; f++) {
				world_step(w, 1.0f/60.0f);
				if (f < 3 || f % 20 == 0) {
					float max_g = 0;
					prev = anchor;
					for (int i = 0; i < 3; i++) {
						float g = anchor_distance(w, prev, oa, chain[i], ob);
						if (g > max_g) max_g = g;
						prev = chain[i];
					}
					printf("  [%s f%d] max_gap=%.6f\n", use_ldl ? "ldl" : "pgs", f, (double)max_g);
				}
			}
			destroy_world(w);
		}
		test_ldl_stress_single_constraint();
		printf("--- results: %d passed, %d failed ---\n", test_pass, test_fail);
	} else {
		for (int i = 1; i < argc; i++) {
			if (strcmp(argv[i], "--bail") == 0) test_bail = 1;
		}
		run_tests();
	}
	return test_fail > 0 ? 1 : 0;
}
