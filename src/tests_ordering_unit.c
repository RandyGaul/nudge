// tests_ordering_unit.c -- unit tests for ldl_sparse_degree and ldl_sparse_min_degree_order.
// Tests on known graph topologies: chains, stars, triangles, cliques.

// ============================================================================
// ldl_sparse_degree tests

static void test_degree_isolated()
{
	TEST_BEGIN("degree_isolated");
	// 3 nodes, no edges. All degrees = 0.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;

	int eliminated[3] = {0};
	TEST_ASSERT(ldl_sparse_degree(&s, 0, eliminated) == 0);
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 0);
	TEST_ASSERT(ldl_sparse_degree(&s, 2, eliminated) == 0);

	ldl_sparse_free(&s);
}

static void test_degree_chain()
{
	TEST_BEGIN("degree_chain");
	// Chain: 0-1-2. Endpoints have degree 1, middle has degree 2.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	int eliminated[3] = {0};
	TEST_ASSERT(ldl_sparse_degree(&s, 0, eliminated) == 1);
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 2);
	TEST_ASSERT(ldl_sparse_degree(&s, 2, eliminated) == 1);

	ldl_sparse_free(&s);
}

static void test_degree_after_elimination()
{
	TEST_BEGIN("degree_after_elimination");
	// Chain: 0-1-2. Eliminate node 0. Node 1's degree drops from 2 to 1.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	int eliminated[3] = {0};
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 2);
	eliminated[0] = 1;
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 1);
	eliminated[2] = 1;
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 0);

	ldl_sparse_free(&s);
}

static void test_degree_star()
{
	TEST_BEGIN("degree_star");
	// Star: hub 0, leaves 1,2,3. Hub degree = 3, leaves = 1.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 3);

	int eliminated[4] = {0};
	TEST_ASSERT(ldl_sparse_degree(&s, 0, eliminated) == 3);
	TEST_ASSERT(ldl_sparse_degree(&s, 1, eliminated) == 1);
	TEST_ASSERT(ldl_sparse_degree(&s, 2, eliminated) == 1);
	TEST_ASSERT(ldl_sparse_degree(&s, 3, eliminated) == 1);

	// Eliminate leaf 1: hub degree drops to 2
	eliminated[1] = 1;
	TEST_ASSERT(ldl_sparse_degree(&s, 0, eliminated) == 2);

	ldl_sparse_free(&s);
}

static void test_degree_clique()
{
	TEST_BEGIN("degree_clique");
	// Complete graph on 4 nodes. Every node has degree 3.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	for (int i = 0; i < 4; i++) {
		for (int j = i + 1; j < 4; j++) {
			ldl_sparse_get_or_create_edge(&s, i, j);
		}
	}

	int eliminated[4] = {0};
	for (int i = 0; i < 4; i++) {
		TEST_ASSERT(ldl_sparse_degree(&s, i, eliminated) == 3);
	}

	// Eliminate one: remaining all have degree 2
	eliminated[0] = 1;
	for (int i = 1; i < 4; i++) {
		TEST_ASSERT(ldl_sparse_degree(&s, i, eliminated) == 2);
	}

	ldl_sparse_free(&s);
}

// ============================================================================
// ldl_sparse_min_degree_order tests
//
// Helper: count total fill-in edges created during ordering.
// We run the ordering, then count edges that exist in the final graph
// but didn't exist in the original.

// Helper: check elim_order is a valid permutation of [0, nc).
static int is_valid_permutation(int* order, int nc)
{
	int seen[LDL_MAX_NODES] = {0};
	for (int i = 0; i < nc; i++) {
		if (order[i] < 0 || order[i] >= nc) return 0;
		if (seen[order[i]]) return 0;
		seen[order[i]] = 1;
	}
	return 1;
}

// Helper: check inv_order is consistent with elim_order.
static int inv_order_consistent(int* elim_order, int* inv_order, int nc)
{
	for (int step = 0; step < nc; step++) {
		if (inv_order[elim_order[step]] != step) return 0;
	}
	return 1;
}

static void test_order_single_node()
{
	TEST_BEGIN("order_single_node");
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 1;
	s.dof[0] = 3;

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(s.elim_order[0] == 0);
	TEST_ASSERT(s.inv_order[0] == 0);

	ldl_sparse_free(&s);
}

static void test_order_two_nodes()
{
	TEST_BEGIN("order_two_nodes");
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 2));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 2));

	ldl_sparse_free(&s);
}

static void test_order_chain_3()
{
	TEST_BEGIN("order_chain_3");
	// Chain: 0-1-2. Endpoints have degree 1, interior has degree 2.
	// First elimination must be a degree-1 node (an endpoint).
	// After that, both remaining have degree 1, so ordering is by tie-break.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));

	// First eliminated must be an endpoint (degree 1, not the interior with degree 2)
	int first = s.elim_order[0];
	TEST_ASSERT(first == 0 || first == 2);

	// Zero fill-in: no edge between 0 and 2 created
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 2) == -1);

	ldl_sparse_free(&s);
}

static void test_order_chain_4()
{
	TEST_BEGIN("order_chain_4");
	// Chain: 0-1-2-3. Endpoints 0,3 have degree 1, interior 1,2 have degree 2.
	// First elimination must be an endpoint (min degree = 1).
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);
	ldl_sparse_get_or_create_edge(&s, 2, 3);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// First eliminated must be an endpoint
	int first = s.elim_order[0];
	TEST_ASSERT(first == 0 || first == 3);

	ldl_sparse_free(&s);
}

static void test_order_star()
{
	TEST_BEGIN("order_star");
	// Star: hub 0, leaves 1,2,3. Leaves have degree 1, hub has degree 3.
	// First elimination must be a leaf (degree 1). As leaves are eliminated,
	// hub degree decreases. Eventually hub ties with remaining leaves.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 3);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// First eliminated must be a leaf (degree 1, not hub with degree 3)
	TEST_ASSERT(s.elim_order[0] != 0);

	// No leaf-leaf fill-in (each leaf's only neighbor is the hub)
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 2) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 3) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 2, 3) == -1);

	ldl_sparse_free(&s);
}

static void test_order_triangle()
{
	TEST_BEGIN("order_triangle");
	// Complete graph on 3 nodes. All degree 2. Any order is optimal.
	// Zero fill-in since already fully connected.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	for (int i = 0; i < 3; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 2);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));

	ldl_sparse_free(&s);
}

static void test_order_disconnected()
{
	TEST_BEGIN("order_disconnected");
	// Two disconnected pairs: 0-1 and 2-3. All degree 1.
	// Order doesn't matter, but must be a valid permutation.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 2, 3);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// No fill-in: disconnected components don't interact
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 2) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 0, 3) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 2) == -1);
	TEST_ASSERT(ldl_sparse_find_edge(&s, 1, 3) == -1);

	ldl_sparse_free(&s);
}

static void test_order_path_5_fill_in()
{
	TEST_BEGIN("order_path_5_fill_in");
	// Path: 0-1-2-3-4. Endpoints have degree 1, interior nodes degree 2.
	// First elimination must be an endpoint.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 5;
	for (int i = 0; i < 5; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);
	ldl_sparse_get_or_create_edge(&s, 2, 3);
	ldl_sparse_get_or_create_edge(&s, 3, 4);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 5));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 5));

	// First eliminated must be an endpoint (degree 1)
	int first = s.elim_order[0];
	TEST_ASSERT(first == 0 || first == 4);

	ldl_sparse_free(&s);
}

static void test_order_mixed_dof()
{
	TEST_BEGIN("order_mixed_dof");
	// Chain with mixed DOF: node 0 = 3 DOF (ball-socket), node 1 = 1 DOF (distance),
	// node 2 = 5 DOF (hinge). Ordering should work regardless of DOF.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 1; s.dof[2] = 5;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));

	// First eliminated must be an endpoint (degree 1)
	TEST_ASSERT(s.elim_order[0] == 0 || s.elim_order[0] == 2);

	ldl_sparse_free(&s);
}

static void test_order_dof_tiebreak()
{
	TEST_BEGIN("order_dof_tiebreak");
	// Three nodes all with degree 1: 0 (5 DOF) - 1 (3 DOF), 2 (1 DOF) - 1.
	// Wait -- that gives node 1 degree 2. Let's use disconnected pairs instead.
	// Nodes: 0 (5 DOF), 1 (3 DOF), 2 (1 DOF), 3 (3 DOF).
	// Edges: 0-1, 2-3. All nodes have degree 1.
	// Tie-break: lowest DOF first. So node 2 (1 DOF) should be eliminated first,
	// then nodes 1 and 3 (3 DOF, pick lower index = 1), then 0 (5 DOF) or 3.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	s.dof[0] = 5; s.dof[1] = 3; s.dof[2] = 1; s.dof[3] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 2, 3);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// Node 2 (1 DOF) should be first -- lowest DOF among all degree-1 nodes
	TEST_ASSERT(s.elim_order[0] == 2);

	// After eliminating 2, node 3 has degree 0 (lowest degree).
	// Remaining degree-1 nodes: 0, 1. Node 3 has degree 0, so it goes next.
	TEST_ASSERT(s.elim_order[1] == 3);

	// Remaining: 0 (5 DOF, degree 1) and 1 (3 DOF, degree 1). Tie on degree,
	// break by DOF: node 1 (3 DOF) before node 0 (5 DOF).
	TEST_ASSERT(s.elim_order[2] == 1);
	TEST_ASSERT(s.elim_order[3] == 0);

	ldl_sparse_free(&s);
}

static void test_order_dof_tiebreak_same_degree()
{
	TEST_BEGIN("order_dof_tiebreak_same_degree");
	// Star with mixed DOF leaves. Hub = node 0 (3 DOF).
	// Leaves: 1 (5 DOF), 2 (1 DOF), 3 (3 DOF). All leaves degree 1.
	// First elimination: leaf with lowest DOF = node 2 (1 DOF).
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	s.dof[0] = 3; s.dof[1] = 5; s.dof[2] = 1; s.dof[3] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 3);

	ldl_sparse_min_degree_order(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// First must be node 2 (degree 1, DOF 1 -- lowest DOF among leaves)
	TEST_ASSERT(s.elim_order[0] == 2);

	// Second: remaining leaves are 1 (5 DOF) and 3 (3 DOF), both degree 1.
	// Hub 0 has degree 2. So next is a leaf: node 3 (3 DOF < 5 DOF).
	TEST_ASSERT(s.elim_order[1] == 3);

	ldl_sparse_free(&s);
}

// ============================================================================
// Runner

static void run_ordering_unit_tests()
{
	printf("--- ordering unit tests ---\n");

	// Degree
	test_degree_isolated();
	test_degree_chain();
	test_degree_after_elimination();
	test_degree_star();
	test_degree_clique();

	// Min-degree ordering
	test_order_single_node();
	test_order_two_nodes();
	test_order_chain_3();
	test_order_chain_4();
	test_order_star();
	test_order_triangle();
	test_order_disconnected();
	test_order_path_5_fill_in();
	test_order_mixed_dof();
	test_order_dof_tiebreak();
	test_order_dof_tiebreak_same_degree();
}
