// tests_dfs_reorder_unit.c -- unit tests for ldl_sparse_dfs_reorder.
// Verifies that the DFS reorder produces a valid permutation, maintains
// inv_order consistency, and groups elimination tree subtrees contiguously.

// Helper: check that every node in the original ordering appears exactly once
// in the reordered output (same set, different sequence).
static int same_node_set(int* order_a, int* order_b, int nc)
{
	int seen_a[LDL_MAX_NODES] = {0}, seen_b[LDL_MAX_NODES] = {0};
	for (int i = 0; i < nc; i++) { seen_a[order_a[i]] = 1; seen_b[order_b[i]] = 1; }
	for (int i = 0; i < nc; i++) if (seen_a[i] != seen_b[i]) return 0;
	return 1;
}

static void test_dfs_single_node()
{
	TEST_BEGIN("dfs_single_node");
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 1;
	s.dof[0] = 3;
	s.elim_order[0] = 0;
	s.inv_order[0] = 0;

	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(s.elim_order[0] == 0);
	TEST_ASSERT(s.inv_order[0] == 0);

	ldl_sparse_free(&s);
}

static void test_dfs_two_nodes()
{
	TEST_BEGIN("dfs_two_nodes");
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 2;
	s.dof[0] = 3; s.dof[1] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);

	// Order: [0, 1]. Step 0 eliminates node 0, step 1 eliminates node 1.
	// Elimination tree: step 0's parent is step 1 (first later step sharing edge).
	// DFS from root (step 1): visit step 0 first, then step 1. But DFS visits
	// children before parent in the output? No -- DFS outputs in visit order.
	// Actually the DFS outputs nodes as it visits them (preorder).
	s.elim_order[0] = 0; s.elim_order[1] = 1;
	s.inv_order[0] = 0; s.inv_order[1] = 1;

	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 2));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 2));

	ldl_sparse_free(&s);
}

static void test_dfs_chain_3()
{
	TEST_BEGIN("dfs_chain_3");
	// Chain: 0-1-2. Run min-degree ordering first, then DFS reorder.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 3; s.dof[2] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	int order_before[3];
	ldl_sparse_min_degree_order(&s);
	for (int i = 0; i < 3; i++) order_before[i] = s.elim_order[i];

	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));
	TEST_ASSERT(same_node_set(order_before, s.elim_order, 3));

	ldl_sparse_free(&s);
}

static void test_dfs_chain_5()
{
	TEST_BEGIN("dfs_chain_5");
	// Path: 0-1-2-3-4. After min-degree + DFS, the ordering should group
	// parent-child pairs in the elimination tree consecutively.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 5;
	for (int i = 0; i < 5; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);
	ldl_sparse_get_or_create_edge(&s, 2, 3);
	ldl_sparse_get_or_create_edge(&s, 3, 4);

	ldl_sparse_min_degree_order(&s);
	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 5));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 5));

	ldl_sparse_free(&s);
}

static void test_dfs_star()
{
	TEST_BEGIN("dfs_star");
	// Star: hub 0, leaves 1,2,3. Min-degree eliminates leaves first.
	// Elimination tree: each leaf step's parent is the hub step.
	// DFS should visit all leaves (children of hub), then hub (root).
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 3);

	ldl_sparse_min_degree_order(&s);
	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	ldl_sparse_free(&s);
}

static void test_dfs_disconnected()
{
	TEST_BEGIN("dfs_disconnected");
	// Two disconnected pairs: 0-1 and 2-3. Each pair forms its own subtree.
	// DFS should visit each subtree contiguously.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 4;
	for (int i = 0; i < 4; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 2, 3);

	ldl_sparse_min_degree_order(&s);
	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 4));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 4));

	// Nodes from each pair should be adjacent in the ordering.
	// Find positions of nodes 0 and 1.
	int pos_0 = s.inv_order[0], pos_1 = s.inv_order[1];
	int pos_2 = s.inv_order[2], pos_3 = s.inv_order[3];
	TEST_ASSERT(abs(pos_0 - pos_1) == 1); // pair {0,1} adjacent
	TEST_ASSERT(abs(pos_2 - pos_3) == 1); // pair {2,3} adjacent

	ldl_sparse_free(&s);
}

static void test_dfs_binary_tree()
{
	TEST_BEGIN("dfs_binary_tree");
	// Binary tree: root 0, children 1 and 2, grandchildren 3,4 under 1 and 5,6 under 2.
	// 7 nodes, complex topology. Verify structural properties after reorder.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 7;
	for (int i = 0; i < 7; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 0, 2);
	ldl_sparse_get_or_create_edge(&s, 1, 3);
	ldl_sparse_get_or_create_edge(&s, 1, 4);
	ldl_sparse_get_or_create_edge(&s, 2, 5);
	ldl_sparse_get_or_create_edge(&s, 2, 6);

	int order_before[7];
	ldl_sparse_min_degree_order(&s);
	for (int i = 0; i < 7; i++) order_before[i] = s.elim_order[i];

	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 7));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 7));
	TEST_ASSERT(same_node_set(order_before, s.elim_order, 7));

	ldl_sparse_free(&s);
}

static void test_dfs_mixed_dof()
{
	TEST_BEGIN("dfs_mixed_dof");
	// Chain with mixed DOF: 0(3)-1(1)-2(5). DFS reorder should work regardless.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	s.dof[0] = 3; s.dof[1] = 1; s.dof[2] = 5;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);

	ldl_sparse_min_degree_order(&s);
	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));

	ldl_sparse_free(&s);
}

static void test_dfs_triangle()
{
	TEST_BEGIN("dfs_triangle");
	// Complete graph on 3 nodes. Fully connected, elimination tree is a chain.
	LDL_Sparse s;
	ldl_sparse_init(&s);
	s.node_count = 3;
	for (int i = 0; i < 3; i++) s.dof[i] = 3;
	ldl_sparse_get_or_create_edge(&s, 0, 1);
	ldl_sparse_get_or_create_edge(&s, 1, 2);
	ldl_sparse_get_or_create_edge(&s, 0, 2);

	ldl_sparse_min_degree_order(&s);
	ldl_sparse_dfs_reorder(&s);

	TEST_ASSERT(is_valid_permutation(s.elim_order, 3));
	TEST_ASSERT(inv_order_consistent(s.elim_order, s.inv_order, 3));

	ldl_sparse_free(&s);
}

// ============================================================================
// Runner

static void run_dfs_reorder_unit_tests()
{
	printf("--- DFS reorder unit tests ---\n");

	test_dfs_single_node();
	test_dfs_two_nodes();
	test_dfs_chain_3();
	test_dfs_chain_5();
	test_dfs_star();
	test_dfs_disconnected();
	test_dfs_binary_tree();
	test_dfs_mixed_dof();
	test_dfs_triangle();
}
