// tests_arena_unit.c -- unit tests for the ckit Arena allocator + arena_push/arena_fit.

static void test_arena_init_reset_free()
{
	TEST_BEGIN("arena_init_reset_free");
	Arena a = {0};
	arena_init(&a, 64 * 1024);
	TEST_ASSERT(a.base != NULL);
	TEST_ASSERT(a.cap == 64 * 1024);
	TEST_ASSERT(a.used == 0);

	arena_alloc(&a, 1000, 16);
	TEST_ASSERT(a.used >= 1000);

	arena_reset(&a);
	TEST_ASSERT(a.used == 0);
	TEST_ASSERT(a.base != NULL);

	arena_free(&a);
	TEST_ASSERT(a.base == NULL);
	TEST_ASSERT(a.cap == 0);
}

static void test_arena_alloc_alignment()
{
	TEST_BEGIN("arena_alloc_alignment");
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	// Force an odd used offset, then request alignment.
	arena_alloc(&a, 1, 1);
	void* p16 = arena_alloc(&a, 32, 16);
	TEST_ASSERT(((uintptr_t)p16 & 15) == 0);

	arena_alloc(&a, 1, 1);
	void* p64 = arena_alloc(&a, 32, 64);
	TEST_ASSERT(((uintptr_t)p64 & 63) == 0);

	arena_free(&a);
}

static void test_arena_push_basic()
{
	TEST_BEGIN("arena_push_basic");
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	int* xs = NULL;
	arena_push(&a, xs, 10);
	arena_push(&a, xs, 20);
	arena_push(&a, xs, 30);

	TEST_ASSERT(asize(xs) == 3);
	TEST_ASSERT(xs[0] == 10);
	TEST_ASSERT(xs[1] == 20);
	TEST_ASSERT(xs[2] == 30);
	TEST_ASSERT(alast(xs) == 30);

	// Header marked is_static so afree is a no-op; arena owns the memory.
	TEST_ASSERT(CK_AHDR(xs)->is_static == 1);

	arena_free(&a);
}

static void test_arena_push_grow()
{
	TEST_BEGIN("arena_push_grow");
	Arena a = {0};
	arena_init(&a, 1 * 1024 * 1024);

	int* xs = NULL;
	// Push enough to force multiple grows. Initial capacity is 16; doubles each grow.
	int n = 10000;
	for (int i = 0; i < n; i++) arena_push(&a, xs, i);

	TEST_ASSERT(asize(xs) == n);
	for (int i = 0; i < n; i++) TEST_ASSERT(xs[i] == i);
	TEST_ASSERT(acap(xs) >= n);

	arena_free(&a);
}

static void test_arena_fit_preallocates()
{
	TEST_BEGIN("arena_fit_preallocates");
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	int* xs = NULL;
	arena_fit(&a, xs, 500);
	TEST_ASSERT(acap(xs) >= 500);
	TEST_ASSERT(asize(xs) == 0);

	// After arena_fit, subsequent arena_push should not need to grow (same pointer).
	int* before = xs;
	for (int i = 0; i < 500; i++) arena_push(&a, xs, i);
	TEST_ASSERT(xs == before);
	TEST_ASSERT(asize(xs) == 500);

	arena_free(&a);
}

static void test_arena_push_read_side_compat()
{
	TEST_BEGIN("arena_push_read_side_compat");
	// apop, aclear, adel, aend must all work on arena-backed arrays.
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	int* xs = NULL;
	for (int i = 0; i < 8; i++) arena_push(&a, xs, i);

	TEST_ASSERT(aend(xs) == xs + 8);
	TEST_ASSERT(apop(xs) == 7);
	TEST_ASSERT(asize(xs) == 7);

	adel(xs, 0); // swap-remove: xs[0] becomes 6 (last element)
	TEST_ASSERT(asize(xs) == 6);
	TEST_ASSERT(xs[0] == 6);

	aclear(xs);
	TEST_ASSERT(asize(xs) == 0);

	// Re-use after clear.
	arena_push(&a, xs, 99);
	TEST_ASSERT(asize(xs) == 1 && xs[0] == 99);

	arena_free(&a);
}

static void test_arena_reset_allows_reuse()
{
	TEST_BEGIN("arena_reset_allows_reuse");
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	int* xs = NULL;
	for (int i = 0; i < 100; i++) arena_push(&a, xs, i);
	size_t used_first = a.used;

	// After reset, xs is dangling; build a fresh array and expect used to climb
	// from zero again (bounded by same pattern).
	arena_reset(&a);
	TEST_ASSERT(a.used == 0);

	int* ys = NULL;
	for (int i = 0; i < 100; i++) arena_push(&a, ys, i * 2);
	TEST_ASSERT(asize(ys) == 100);
	for (int i = 0; i < 100; i++) TEST_ASSERT(ys[i] == i * 2);
	TEST_ASSERT(a.used == used_first); // same deterministic layout after reset

	arena_free(&a);
}

static void test_arena_push_struct_type()
{
	TEST_BEGIN("arena_push_struct_type");
	typedef struct P { float x, y, z; int tag; } P;
	Arena a = {0};
	arena_init(&a, 64 * 1024);

	P* ps = NULL;
	for (int i = 0; i < 200; i++) arena_push(&a, ps, ((P){ (float)i, (float)(i+1), (float)(i+2), i * 7 }));

	TEST_ASSERT(asize(ps) == 200);
	for (int i = 0; i < 200; i++) {
		TEST_ASSERT(ps[i].x == (float)i);
		TEST_ASSERT(ps[i].tag == i * 7);
	}

	// Data portion must be 16-byte aligned.
	TEST_ASSERT(((uintptr_t)ps & 15) == 0);

	arena_free(&a);
}

static void test_arena_multi_independent()
{
	TEST_BEGIN("arena_multi_independent");
	// Two arenas, two arrays. Writes to one must not touch the other.
	Arena a = {0}, b = {0};
	arena_init(&a, 64 * 1024);
	arena_init(&b, 64 * 1024);

	int* xa = NULL;
	int* xb = NULL;
	for (int i = 0; i < 1000; i++) { arena_push(&a, xa, i); arena_push(&b, xb, -i); }

	for (int i = 0; i < 1000; i++) {
		TEST_ASSERT(xa[i] == i);
		TEST_ASSERT(xb[i] == -i);
	}

	arena_free(&a);
	arena_free(&b);
}

static void run_arena_unit_tests()
{
	printf("--- arena unit tests ---\n");
	test_arena_init_reset_free();
	test_arena_alloc_alignment();
	test_arena_push_basic();
	test_arena_push_grow();
	test_arena_fit_preallocates();
	test_arena_push_read_side_compat();
	test_arena_reset_allows_reuse();
	test_arena_push_struct_type();
	test_arena_multi_independent();
}
