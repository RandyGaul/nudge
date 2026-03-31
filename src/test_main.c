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
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--fuzz") == 0 && i + 1 < argc)
			fuzz_iters = atoi(argv[++i]);
		else if (strcmp(argv[i], "--soak") == 0)
			soak = 1;
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
	} else {
		run_tests();
	}
	return test_fail > 0 ? 1 : 0;
}
