// debug_ui.c -- debug visualization panels for nudge demo app
// Included from main.c (unity build). Accesses globals defined there:
//   g_world, g_show_joints, g_show_bvh, g_ldl_inspect_island,
//   g_ldl_inspect_step, g_ldl_hover_body, WorldInternal, render functions, ImGui

// --- LDL Inspector ---
extern LDL_DebugInfo g_ldl_debug_info;
extern int g_ldl_debug_enabled;

static const char* ldl_constraint_type_name(int type)
{
	if (type == JOINT_BALL_SOCKET) return "BS";
	if (type == JOINT_DISTANCE) return "Dist";
	if (type == JOINT_HINGE) return "Hinge";
	if (type == JOINT_FIXED) return "Fixed";
	if (type == JOINT_PRISMATIC) return "Prism";
	return "Weld";
}

static void draw_ldl_overview(WorldInternal* w, Island* isl, LDL_Cache* c)
{
	LDL_Topology* t = c->topo;

	if (ImGui_TreeNodeEx("System", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui_Text("DOFs: %d  Nodes: %d  Constraints: %d", t ? t->n : 0, c->bundle_count, c->joint_count);
		ImGui_Text("Topo version: %d  L_factors: %d floats", c->topo_version, t ? t->L_factors_size : 0);
		if (c->virtual_body_count > 0) {
			ImGui_Text("Shattered: %d virtual bodies", c->virtual_body_count);
		}
		ImGui_TreePop();
	}

	if (ImGui_TreeNodeEx("Bodies", 0)) {
		int bi = isl->head_body;
		while (bi >= 0) {
			float mass = w->body_hot[bi].inv_mass > 0 ? 1.0f / w->body_hot[bi].inv_mass : 0;
			if (w->body_hot[bi].inv_mass == 0) {
				ImGui_Text("  [%d] static", bi);
			} else {
				v3 inv_I = w->body_hot[bi].inv_inertia_local;
				ImGui_Text("  [%d] mass=%.1f  inv_I=(%.2f, %.2f, %.2f)", bi, (double)mass, (double)inv_I.x, (double)inv_I.y, (double)inv_I.z);
			}
			if (ImGui_IsItemHovered(0)) {
				g_ldl_hover_body = bi;
			}
			bi = w->body_cold[bi].island_next;
		}
		ImGui_TreePop();
	}

	if (c->virtual_body_count > 0 && c->body_remap && ImGui_TreeNodeEx("Shattering", ImGuiTreeNodeFlags_DefaultOpen)) {
		int real_count = asize(w->body_hot);
		int bi = isl->head_body;
		while (bi >= 0) {
			if (c->body_remap[bi] >= 0) {
				int shard_n = c->shard_counts[bi];
				float real_mass = w->body_hot[bi].inv_mass > 0 ? 1.0f / w->body_hot[bi].inv_mass : 0;
				ImGui_Text("  Body %d (mass=%.1f) -> %d shards", bi, (double)real_mass, shard_n);
				if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
					ImGui_Text("Hub body %d was shattered into %d virtual shards.", bi, shard_n);
					ImGui_Text("Total DOF on this body exceeded %d.", SHATTER_THRESHOLD);
					ImGui_Separator();
					ImGui_Text("Each shard has mass = %.1f * %d = %.1f (inv_mass scaled up)", (double)real_mass, shard_n, (double)(real_mass * shard_n));
					ImGui_Text("Shards are connected by synthetic 6-DOF weld joints");
					ImGui_Text("in a wrap-around ring: s0-s1-s2-...-sN-s0");
					ImGui_Separator();
					ImGui_Text("Virtual body indices:");
					for (int v = 0; v < shard_n; v++) {
						int vid = c->body_remap[bi] + v;
						ImGui_Text("  shard %d = virtual body %d", v, vid);
					}
					ImGui_Separator();
					ImGui_Text("Real constraints were distributed across shards");
					ImGui_Text("using greedy bin-packing (largest DOF first).");
					ImGui_EndTooltip();
				}
				if (ImGui_IsItemHovered(0)) {
					g_ldl_hover_body = bi;
				}
			}
			bi = w->body_cold[bi].island_next;
		}

		// Show which constraints got redirected to which shards
		int synth_count = 0, real_redirected = 0;
		for (int i = 0; i < c->joint_count; i++) {
			if (c->constraints[i].is_synthetic) synth_count++;
			else if (c->constraints[i].body_a >= real_count || c->constraints[i].body_b >= real_count) real_redirected++;
		}
		ImGui_Text("  %d real constraints redirected to shards", real_redirected);
		ImGui_Text("  %d synthetic weld joints between shards", synth_count);
		ImGui_TreePop();
	}

	if (ImGui_TreeNodeEx("Constraints", 0)) {
		for (int i = 0; i < c->joint_count; i++) {
			LDL_Constraint* con = &c->constraints[i];
			if (con->is_synthetic) {
				ImGui_Text("  [%d] Weld (synthetic)  body %d <-> %d  dof=%d  bundle %d", i, con->body_a, con->body_b, con->dof, con->bundle_idx);
			} else {
				ImGui_Text("  [%d] %s  body %d <-> %d  dof=%d  bundle %d", i, ldl_constraint_type_name(con->type), con->body_a, con->body_b, con->dof, con->bundle_idx);
			}
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Constraint %d", i);
				ImGui_Text("Type: %s  DOF: %d", con->is_synthetic ? "Synthetic weld" : ldl_constraint_type_name(con->type), con->dof);
				ImGui_Text("Bodies: %d <-> %d", con->body_a, con->body_b);
				ImGui_Text("Bundle: %d  offset: %d", con->bundle_idx, con->bundle_offset);
				if (con->is_synthetic) {
					ImGui_Separator();
					ImGui_Text("Synthetic 6-DOF weld between virtual shards.");
					ImGui_Text("Constrains: v_b - v_a = 0 (linear, 3 DOF)");
					ImGui_Text("            w_b - w_a = 0 (angular, 3 DOF)");
					ImGui_Text("Forces shards to move as one rigid body.");
				} else {
					ImGui_Text("Solver index: %d", con->solver_idx);
				}
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}

	if (ImGui_TreeNodeEx("Bundles", 0)) {
		for (int i = 0; i < c->bundle_count; i++) {
			LDL_Bundle* b = &c->bundles[i];
			ImGui_Text("  [%d] body(%d,%d)  dof=%d  constraints=%d", i, b->body_a, b->body_b, b->dof, b->count);
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Bundle %d (graph node)", i);
				ImGui_Text("Body pair: %d <-> %d", b->body_a, b->body_b);
				ImGui_Text("Total DOF: %d  from %d constraint(s)", b->dof, b->count);
				ImGui_Text("Constraint range: [%d..%d)", b->start, b->start + b->count);
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}

	if (t && ImGui_TreeNodeEx("Elimination Order", 0)) {
		for (int s = 0; s < t->node_count; s++) {
			LDL_Pivot* pv = &t->pivots[s];
			ImGui_Text("  step %d: node %d (dof %d)  fwd=%d back=%d col=%d schur=%d", s, pv->node, pv->dk, pv->fwd_count, pv->back_count, pv->col_count, pv->schur_count);
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Elimination step %d", s);
				ImGui_Text("Eliminates node %d (bundle %d)", pv->node, pv->node);
				ImGui_Text("DOF: %d  Row offset: %d", pv->dk, pv->ok);
				ImGui_Separator();
				ImGui_Text("Forward-sub neighbors: %d (eliminated before this)", pv->fwd_count);
				ImGui_Text("Back-sub neighbors: %d (eliminated after this)", pv->back_count);
				ImGui_Text("L-column entries: %d (L blocks computed)", pv->col_count);
				ImGui_Text("Schur updates: %d (S -= L*N*L^T ops)", pv->schur_count);
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}
}

static ImU32 ldl_heat(float val, float max_val)
{
	if (max_val < 1e-12f) return ImGui_GetColorU32ImVec4((ImVec4){0.1f, 0.1f, 0.1f, 1.0f});
	float t = fabsf(val) / max_val;
	if (t > 1.0f) t = 1.0f;
	float r, g, b;
	if (t < 0.33f) { float s = t / 0.33f; r = 0; g = 0; b = s; }
	else if (t < 0.66f) { float s = (t - 0.33f) / 0.33f; r = s; g = 0; b = 1.0f - s; }
	else { float s = (t - 0.66f) / 0.34f; r = 1.0f; g = s; b = 0; }
	return ImGui_GetColorU32ImVec4((ImVec4){r, g, b, 1.0f});
}

static void draw_ldl_matrix(WorldInternal* w, Island* isl, LDL_Cache* c)
{
	LDL_Topology* t = c->topo;
	LDL_DebugInfo* info = &g_ldl_debug_info;
	if (!t || !info->valid || info->n == 0) {
		ImGui_Text("No matrix data (island may be sleeping)");
		return;
	}
	int n = info->n;
	int nc = t->node_count;
	if (n > LDL_MAX_DOF) { ImGui_Text("System too large to display (%d DOFs)", n); return; }

	float cell = 10.0f;
	float margin = 40.0f;
	ImDrawList* dl = ImGui_GetWindowDrawList();
	ImVec2 cursor = ImGui_GetCursorScreenPos();
	float ox = cursor.x + margin, oy = cursor.y + margin;
	ImU32 text_col = ImGui_GetColorU32ImVec4((ImVec4){0.7f, 0.7f, 0.7f, 1.0f});

	float a_max = 0;
	for (int i = 0; i < n * n; i++) {
		float v = fabsf(info->A[i]);
		if (v > a_max) a_max = v;
	}

	// Node labels (top and left). Synthetic weld nodes shown in cyan.
	ImU32 synth_col = ImGui_GetColorU32ImVec4((ImVec4){0.0f, 0.8f, 1.0f, 1.0f});
	for (int j = 0; j < nc; j++) {
		char lbl[16];
		LDL_Bundle* bun = &c->bundles[j];
		int is_synth = (bun->count == 1 && c->constraints[bun->start].is_synthetic);
		if (is_synth) {
			snprintf(lbl, sizeof(lbl), "W%d", j);
		} else if (bun->count == 1) {
			snprintf(lbl, sizeof(lbl), "%s%d", ldl_constraint_type_name(c->constraints[bun->start].type), j);
		} else {
			snprintf(lbl, sizeof(lbl), "B%d", j);
		}
		ImU32 lbl_col = is_synth ? synth_col : text_col;
		float bx = ox + t->row_offset[j] * cell + t->dof[j] * cell * 0.5f - 8;
		ImDrawList_AddText(dl, (ImVec2){bx, cursor.y}, lbl_col, lbl);
		float by = oy + t->row_offset[j] * cell + t->dof[j] * cell * 0.5f - 6;
		ImDrawList_AddText(dl, (ImVec2){cursor.x, by}, lbl_col, lbl);
	}

	// Draw cells
	int hover_row = -1, hover_col = -1;
	ImVec2 mouse = ImGui_GetMousePos();
	for (int row = 0; row < n; row++) {
		for (int col = 0; col < n; col++) {
			float x0 = ox + col * cell, y0 = oy + row * cell;
			ImU32 clr = ldl_heat(info->A[row * n + col], a_max);
			ImDrawList_AddRectFilled(dl, (ImVec2){x0, y0}, (ImVec2){x0 + cell - 1, y0 + cell - 1}, clr);
			if (mouse.x >= x0 && mouse.x < x0 + cell && mouse.y >= y0 && mouse.y < y0 + cell) {
				hover_row = row; hover_col = col;
			}
		}
	}

	// Block boundary lines
	ImU32 line_col = ImGui_GetColorU32ImVec4((ImVec4){1, 1, 1, 0.3f});
	for (int j = 0; j < nc; j++) {
		float px = ox + t->row_offset[j] * cell;
		float py = oy + t->row_offset[j] * cell;
		ImDrawList_AddLine(dl, (ImVec2){px, oy}, (ImVec2){px, oy + n * cell}, line_col);
		ImDrawList_AddLine(dl, (ImVec2){ox, py}, (ImVec2){ox + n * cell, py}, line_col);
	}
	ImDrawList_AddRect(dl, (ImVec2){ox, oy}, (ImVec2){ox + n * cell, oy + n * cell}, line_col);
	ImGui_Dummy((ImVec2){margin + n * cell + 4, margin + n * cell + 4});

	// Matrix hover tooltip with physics explanation
	if (hover_row >= 0 && ImGui_IsMouseHoveringRect((ImVec2){ox, oy}, (ImVec2){ox + n * cell, oy + n * cell})) {
		// Find which bundle owns each row/col
		int br = -1, bc_idx = -1;
		for (int j = nc - 1; j >= 0; j--) {
			if (hover_row >= t->row_offset[j]) { br = j; break; }
		}
		for (int j = nc - 1; j >= 0; j--) {
			if (hover_col >= t->row_offset[j]) { bc_idx = j; break; }
		}
		if (br >= 0 && bc_idx >= 0 && ImGui_BeginTooltip()) {
			float val = info->A[hover_row * n + hover_col];
			ImGui_Text("A[%d,%d] = %.6g", hover_row, hover_col, (double)val);
			ImGui_Separator();
			LDL_Bundle* bun_r = &c->bundles[br];
			LDL_Bundle* bun_c = &c->bundles[bc_idx];
			int r_synth = (bun_r->count == 1 && c->constraints[bun_r->start].is_synthetic);
			int c_synth = (bun_c->count == 1 && c->constraints[bun_c->start].is_synthetic);
			if (br == bc_idx) {
				if (r_synth) {
					ImGui_Text("Diagonal block: node %d (synthetic weld)", br);
					ImGui_Text("K = diag(inv_m_sum * I_3, inv_I_sum)");
					ImGui_TextDisabled("6-DOF weld between virtual shards %d, %d", bun_r->body_a, bun_r->body_b);
					ImGui_TextDisabled("Top-left 3x3: linear coupling (mass)");
					ImGui_TextDisabled("Bottom-right 3x3: angular coupling (inertia)");
				} else {
					ImGui_Text("Diagonal block: node %d", br);
					ImGui_Text("K_ii = J_i * M^-1 * J_i^T");
					ImGui_TextDisabled("Effective mass of constraint through bodies %d, %d", bun_r->body_a, bun_r->body_b);
				}
			} else {
				// Find shared body between the two bundles
				int shared = -1;
				if (bun_r->body_a == bun_c->body_a || bun_r->body_a == bun_c->body_b) shared = bun_r->body_a;
				else if (bun_r->body_b == bun_c->body_a || bun_r->body_b == bun_c->body_b) shared = bun_r->body_b;
				if (shared >= 0) {
					ImGui_Text("Off-diagonal: node %d <-> node %d", br, bc_idx);
					ImGui_Text("K_ij = J_i * M_%d^-1 * J_j^T", shared);
					ImGui_TextDisabled("Coupling through shared body %d", shared);
					g_ldl_hover_body = shared;
				} else if (fabsf(val) > 1e-12f) {
					ImGui_Text("Fill-in block: node %d <-> node %d", br, bc_idx);
					ImGui_TextDisabled("No direct body sharing (created by elimination)");
				} else {
					ImGui_Text("Zero block: node %d, node %d", br, bc_idx);
					ImGui_TextDisabled("No coupling (no shared body)");
				}
			}
			ImGui_EndTooltip();
		}
	}

	// Color legend
	{
		ImVec2 lc = ImGui_GetCursorScreenPos();
		ImDrawList* ldl = ImGui_GetWindowDrawList();
		float lw = 200.0f, lh = 12.0f;
		for (int i = 0; i < (int)lw; i++) {
			float t = (float)i / lw;
			ImU32 clr = ldl_heat(t, 1.0f);
			ImDrawList_AddRectFilled(ldl, (ImVec2){lc.x + i, lc.y}, (ImVec2){lc.x + i + 1, lc.y + lh}, clr);
		}
		ImU32 tc = ImGui_GetColorU32ImVec4((ImVec4){0.7f, 0.7f, 0.7f, 1.0f});
		ImDrawList_AddText(ldl, (ImVec2){lc.x, lc.y + lh + 2}, tc, "0");
		char maxlbl[16]; snprintf(maxlbl, sizeof(maxlbl), "%.3g", (double)a_max);
		ImDrawList_AddText(ldl, (ImVec2){lc.x + lw - 30, lc.y + lh + 2}, tc, maxlbl);
		ImGui_Dummy((ImVec2){lw, lh + 16});
	}
	if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
		ImGui_Text("|A[r,c]| mapped to color");
		ImGui_Text("Black = zero (no coupling)");
		ImGui_Text("Blue = small magnitude");
		ImGui_Text("Red = medium magnitude");
		ImGui_Text("Yellow = large magnitude (max = %.3g)", (double)a_max);
		ImGui_EndTooltip();
	}

	// D pivots
	ImGui_SeparatorText("D Pivots (LDL^T diagonal)");
	if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
		ImGui_Text("K = L * D * L^T");
		ImGui_Text("D contains the diagonal pivots from factorization.");
		ImGui_Text("Large ratio = ill-conditioned system.");
		ImGui_Text("PGS struggles with high ratios; LDL solves exactly.");
		ImGui_EndTooltip();
	}
	float d_min = 1e18f, d_max = -1e18f;
	for (int i = 0; i < n; i++) {
		float d = info->D[i];
		if (d < d_min) d_min = d;
		if (d > d_max) d_max = d;
	}
	float ratio = d_max / (d_min > 1e-12f ? d_min : 1e-12f);
	ImGui_Text("min=%.4g  max=%.4g  ratio=%.1f", (double)d_min, (double)d_max, (double)ratio);
}

static void draw_ldl_factorization(WorldInternal* w, Island* isl, LDL_Cache* c)
{
	LDL_Topology* t = c->topo;
	if (!t || t->node_count == 0) { ImGui_Text("No topology"); return; }

	int nc = t->node_count;
	if (g_ldl_inspect_step >= nc) g_ldl_inspect_step = nc - 1;
	if (g_ldl_inspect_step < 0) g_ldl_inspect_step = 0;

	ImGui_Text("Step %d / %d", g_ldl_inspect_step + 1, nc);
	ImGui_SameLine();
	if (ImGui_Button("<##fwd")) { if (g_ldl_inspect_step > 0) g_ldl_inspect_step--; }
	ImGui_SameLine();
	if (ImGui_Button(">##back")) { if (g_ldl_inspect_step < nc - 1) g_ldl_inspect_step++; }

	int step = g_ldl_inspect_step;
	LDL_Pivot* pv = &t->pivots[step];
	int node = pv->node;

	ImGui_Separator();

	// Node info
	LDL_Bundle* bun = &c->bundles[node];
	ImGui_Text("Eliminating: node %d  (%s, dof=%d)", node, bun->count == 1 ? ldl_constraint_type_name(c->constraints[bun->start].type) : "bundle", pv->dk);
	if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
		ImGui_Text("Gaussian elimination on the constraint graph.");
		ImGui_Text("Node %d is removed. Its neighbors become", node);
		ImGui_Text("connected to each other (fill-in edges).");
		ImGui_Text("L blocks are computed: L_{j,%d} = E_{j,%d} * K_%d^-1", node, node, node);
		ImGui_EndTooltip();
	}

	// D pivots for this node
	ImGui_Text("Pivot D:");
	ImGui_SameLine();
	for (int d = 0; d < pv->dk; d++) {
		ImGui_SameLine();
		ImGui_Text("%.3g", (double)c->diag_D[node][d]);
	}
	if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
		ImGui_Text("Diagonal pivots from block LDL factorization");
		ImGui_Text("of this node's effective mass matrix.");
		ImGui_Text("Small pivots indicate near-singularity.");
		ImGui_EndTooltip();
	}

	// L columns
	if (pv->col_count > 0 && ImGui_TreeNodeEx("L columns", ImGuiTreeNodeFlags_DefaultOpen)) {
		for (int ci = 0; ci < pv->col_count; ci++) {
			LDL_Column* col = &t->columns[pv->col_start + ci];
			ImGui_Text("  L_{%d,%d}: %dx%d block at L_factors[%d]", col->node, node, col->dn, pv->dk, col->L_offset);
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("L_{%d,%d} = E_{%d,%d} * K_%d^-1", col->node, node, col->node, node, node);
				ImGui_Text("This block of the lower-triangular factor L");
				ImGui_Text("encodes how constraint %d depends on %d.", col->node, node);
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}

	// Schur updates
	if (pv->schur_count > 0 && ImGui_TreeNodeEx("Schur updates", ImGuiTreeNodeFlags_DefaultOpen)) {
		for (int si = 0; si < pv->schur_count; si++) {
			LDL_Schur* op = &t->schurs[pv->schur_start + si];
			if (op->target_offset < 0) {
				ImGui_Text("  S_{%d,%d} -= L_{%d,%d} * N_%d * L_{%d,%d}^T  (diagonal)", op->i, op->j, op->i, node, node, op->j, node);
			} else {
				ImGui_Text("  S_{%d,%d} -= L_{%d,%d} * N_%d * L_{%d,%d}^T", op->i, op->j, op->i, node, node, op->j, node);
			}
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Schur complement update");
				ImGui_Text("When node %d is eliminated, nodes %d and %d", node, op->i, op->j);
				if (op->i == op->j) {
					ImGui_Text("get their diagonal block reduced.");
				} else {
					ImGui_Text("become coupled (fill-in if they weren't already).");
				}
				ImGui_Text("N_%d = original diagonal of node %d before factoring.", node, node);
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}
}

static void draw_ldl_solve(WorldInternal* w, Island* isl, LDL_Cache* c)
{
	LDL_Topology* t = c->topo;
	LDL_DebugInfo* info = &g_ldl_debug_info;
	if (!t || !info->valid) { ImGui_Text("No solve data"); return; }
	int n = info->n;

	if (ImGui_TreeNodeEx("RHS (constraint violations)", ImGuiTreeNodeFlags_DefaultOpen)) {
		for (int i = 0; i < c->bundle_count; i++) {
			LDL_Bundle* bun = &c->bundles[i];
			int oi = t->row_offset[i];
			char lbl[16];
			if (bun->count == 1) {
				snprintf(lbl, sizeof(lbl), "%s%d", ldl_constraint_type_name(c->constraints[bun->start].type), i);
			} else {
				snprintf(lbl, sizeof(lbl), "B%d", i);
			}
			if (bun->dof == 3) {
				ImGui_Text("  %s: (%.4f, %.4f, %.4f)", lbl, (double)info->rhs[oi], (double)info->rhs[oi+1], (double)info->rhs[oi+2]);
			} else if (bun->dof == 1) {
				ImGui_Text("  %s: %.4f", lbl, (double)info->rhs[oi]);
			} else {
				ImGui_Text("  %s: [%d DOF]", lbl, bun->dof);
			}
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Right-hand side: b = -(J*v + bias)");
				ImGui_Text("Measures how much this constraint is violated.");
				ImGui_Text("Zero = constraint perfectly satisfied.");
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}

	if (ImGui_TreeNodeEx("Lambda (exact impulses)", ImGuiTreeNodeFlags_DefaultOpen)) {
		for (int i = 0; i < c->bundle_count; i++) {
			LDL_Bundle* bun = &c->bundles[i];
			int oi = t->row_offset[i];
			char lbl[16];
			if (bun->count == 1) {
				snprintf(lbl, sizeof(lbl), "%s%d", ldl_constraint_type_name(c->constraints[bun->start].type), i);
			} else {
				snprintf(lbl, sizeof(lbl), "B%d", i);
			}
			if (bun->dof == 3) {
				ImGui_Text("  %s: (%.4f, %.4f, %.4f)", lbl, (double)info->lambda_ldl[oi], (double)info->lambda_ldl[oi+1], (double)info->lambda_ldl[oi+2]);
			} else if (bun->dof == 1) {
				ImGui_Text("  %s: %.4f", lbl, (double)info->lambda_ldl[oi]);
			} else {
				ImGui_Text("  %s: [%d DOF]", lbl, bun->dof);
			}
			if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
				ImGui_Text("Lambda = K^-1 * b");
				ImGui_Text("Exact constraint impulses from LDL solve.");
				ImGui_Text("Applied to bodies to correct velocity errors.");
				ImGui_EndTooltip();
			}
		}
		ImGui_TreePop();
	}

	ImGui_SeparatorText("Solve Residual");
	if (ImGui_IsItemHovered(0) && ImGui_BeginTooltip()) {
		ImGui_Text("r = A*lambda - rhs");
		ImGui_Text("Measures factorization/solve accuracy.");
		ImGui_Text("Zero = exact solution (expected for LDL).");
		ImGui_EndTooltip();
	}
	float max_res = 0, rms = 0;
	for (int i = 0; i < n; i++) {
		float Ax_i = 0;
		for (int j = 0; j < n; j++) Ax_i += info->A[i * n + j] * info->lambda_ldl[j];
		float r = fabsf(Ax_i - info->rhs[i]);
		if (r > max_res) max_res = r;
		rms += r * r;
	}
	rms = sqrtf(rms / (n > 0 ? n : 1));
	ImGui_Text("Max residual: %.4g  RMS: %.4g", (double)max_res, (double)rms);
}

static void draw_ldl_inspector()
{
	g_ldl_hover_body = -1; // reset each frame

	if (g_ldl_inspect_island < 0) return;
	WorldInternal* w = (WorldInternal*)g_world.id;
	if (!(w->island_gen[g_ldl_inspect_island] & 1)) { g_ldl_inspect_island = -1; return; }

	Island* isl = &w->islands[g_ldl_inspect_island];
	LDL_Cache* c = &isl->ldl;

	ImGui_Begin("LDL Inspector", NULL, 0);

	ImGui_Text("Island %d: %d bodies, %d joints, %s", g_ldl_inspect_island, isl->body_count, isl->joint_count, isl->awake ? "awake" : "sleeping");
	ImGui_SameLine();
	if (ImGui_Button("Deselect")) { g_ldl_inspect_island = -1; ImGui_End(); return; }

	if (!c->topo) {
		ImGui_Text("No LDL data (island has no joints or is sleeping)");
		ImGui_End();
		return;
	}

	if (ImGui_BeginTabBar("ldl_tabs", 0)) {
		if (ImGui_BeginTabItem("Overview", NULL, 0)) {
			draw_ldl_overview(w, isl, c);
			ImGui_EndTabItem();
		}
		if (ImGui_BeginTabItem("Matrix", NULL, 0)) {
			draw_ldl_matrix(w, isl, c);
			ImGui_EndTabItem();
		}
		if (ImGui_BeginTabItem("Factorization", NULL, 0)) {
			draw_ldl_factorization(w, isl, c);
			ImGui_EndTabItem();
		}
		if (ImGui_BeginTabItem("Solve", NULL, 0)) {
			draw_ldl_solve(w, isl, c);
			ImGui_EndTabItem();
		}
		ImGui_EndTabBar();
	}

	ImGui_End();
}

static void bvh_debug_draw_cb(v3 mn, v3 mx, int depth, int is_leaf, void* user)
{
	(void)user;
	// Color by depth: cycle through distinct hues.
	float hues[][3] = { {1,0.3f,0.3f}, {0.3f,1,0.3f}, {0.3f,0.3f,1}, {1,1,0.3f}, {1,0.3f,1}, {0.3f,1,1} };
	int hi = depth % 6;
	float a = is_leaf ? 1.0f : 0.5f;
	v3 col = V3(hues[hi][0] * a, hues[hi][1] * a, hues[hi][2] * a);
	draw_aabb_wireframe(mn, mx, col);
}

// Draw an arc of debug lines around axis at center, from angle0 to angle1.
// ref is the 0-angle direction (must be perpendicular to axis).
static void draw_debug_arc(v3 center, v3 axis, v3 ref, float angle0, float angle1, float radius, v3 color, int segments)
{
	v3 perp = cross(axis, ref);
	float step = (angle1 - angle0) / segments;
	for (int i = 0; i < segments; i++) {
		float a0 = angle0 + i * step, a1 = angle0 + (i + 1) * step;
		v3 p0 = add(center, add(scale(ref, cosf(a0) * radius), scale(perp, sinf(a0) * radius)));
		v3 p1 = add(center, add(scale(ref, cosf(a1) * radius), scale(perp, sinf(a1) * radius)));
		render_debug_line(p0, p1, color);
	}
}

static void draw_joint_debug(JointDebugInfo info, void* user)
{
	(void)user;
	v3 a = info.anchor_a, b = info.anchor_b;
	v3 mid = scale(add(a, b), 0.5f);

	if (info.type == JOINT_BALL_SOCKET) {
		v3 col = info.is_soft ? V3(0.2f, 0.8f, 1.0f) : V3(1.0f, 0.6f, 0.1f);
		render_debug_line(a, b, col);
	} else if (info.type == JOINT_DISTANCE) {
		v3 col = info.is_soft ? V3(0.2f, 1.0f, 0.4f) : V3(1.0f, 1.0f, 0.2f);
		render_debug_line(a, b, col);
	} else if (info.type == JOINT_HINGE) {
		v3 col = V3(1.0f, 0.3f, 0.8f);
		render_debug_line(a, b, col);
		// Axis indicator
		render_debug_line(mid, add(mid, scale(info.axis_a, 0.3f)), col);

		// Limits: draw arc showing allowed range
		if (info.limit_min != 0 || info.limit_max != 0) {
			v3 ref = norm(info.ref_a);
			float radius = 0.4f;
			// Limit range arc (dim)
			v3 limit_col = V3(0.5f, 0.5f, 0.3f);
			draw_debug_arc(mid, info.axis_a, ref, info.limit_min, info.limit_max, radius, limit_col, 16);
			// Limit boundary lines (bright)
			v3 perp = cross(info.axis_a, ref);
			v3 lo_dir = add(scale(ref, cosf(info.limit_min) * radius), scale(perp, sinf(info.limit_min) * radius));
			v3 hi_dir = add(scale(ref, cosf(info.limit_max) * radius), scale(perp, sinf(info.limit_max) * radius));
			v3 bound_col = info.limit_active ? V3(1.0f, 0.2f, 0.2f) : V3(0.8f, 0.8f, 0.3f);
			render_debug_line(mid, add(mid, lo_dir), bound_col);
			render_debug_line(mid, add(mid, hi_dir), bound_col);
			// Current angle indicator (white)
			v3 cur_dir = add(scale(ref, cosf(info.current_angle) * radius * 0.8f), scale(perp, sinf(info.current_angle) * radius * 0.8f));
			render_debug_line(mid, add(mid, cur_dir), V3(1, 1, 1));
		}

		// Motor: direction arrow
		if (info.motor_max_impulse > 0) {
			v3 motor_col = V3(0.2f, 1.0f, 0.6f);
			// Small arc arrow showing motor direction
			float sign = info.motor_speed >= 0 ? 1.0f : -1.0f;
			v3 ref = norm(info.ref_a);
			float r = 0.25f;
			draw_debug_arc(mid, info.axis_a, ref, info.current_angle, info.current_angle + sign * 0.5f, r, motor_col, 6);
			// Arrowhead at end of arc
			v3 perp = cross(info.axis_a, ref);
			float end_angle = info.current_angle + sign * 0.5f;
			v3 tip = add(mid, add(scale(ref, cosf(end_angle) * r), scale(perp, sinf(end_angle) * r)));
			v3 tangent = add(scale(ref, -sinf(end_angle) * sign), scale(perp, cosf(end_angle) * sign));
			render_debug_line(tip, add(tip, add(scale(tangent, -0.08f), scale(info.axis_a, 0.04f))), motor_col);
			render_debug_line(tip, add(tip, add(scale(tangent, -0.08f), scale(info.axis_a, -0.04f))), motor_col);
		}
	} else if (info.type == JOINT_FIXED) {
		// 6DOF weld: draw cross at anchor + line between anchors
		v3 col = V3(0.6f, 0.6f, 0.9f);
		render_debug_line(a, b, col);
		float s = 0.15f;
		render_debug_line(add(mid, V3(s, 0, 0)), add(mid, V3(-s, 0, 0)), col);
		render_debug_line(add(mid, V3(0, s, 0)), add(mid, V3(0, -s, 0)), col);
		render_debug_line(add(mid, V3(0, 0, s)), add(mid, V3(0, 0, -s)), col);
	} else if (info.type == JOINT_PRISMATIC) {
		// Slide axis rail + anchor line
		v3 col = V3(0.3f, 0.8f, 1.0f);
		render_debug_line(a, b, col);
		// Rail indicator: line along axis through anchor
		v3 rail_half = scale(info.axis_a, 0.5f);
		v3 rail_col = V3(0.2f, 0.5f, 0.7f);
		render_debug_line(sub(mid, rail_half), add(mid, rail_half), rail_col);
		// Cross-hatches perpendicular to axis at ends
		v3 t1, t2;
		float ax = fabsf(info.axis_a.x);
		v3 up = ax < 0.9f ? V3(1, 0, 0) : V3(0, 1, 0);
		t1 = norm(cross(info.axis_a, up));
		t2 = cross(info.axis_a, t1);
		float hs = 0.08f;
		v3 end1 = add(mid, rail_half), end2 = sub(mid, rail_half);
		render_debug_line(add(end1, scale(t1, hs)), add(end1, scale(t1, -hs)), rail_col);
		render_debug_line(add(end2, scale(t1, hs)), add(end2, scale(t1, -hs)), rail_col);

		// Motor: arrow along axis
		if (info.motor_max_impulse > 0) {
			v3 motor_col = V3(0.2f, 1.0f, 0.6f);
			float sign = info.motor_speed >= 0 ? 1.0f : -1.0f;
			v3 arrow_end = add(mid, scale(info.axis_a, 0.35f * sign));
			render_debug_line(mid, arrow_end, motor_col);
			// Arrowhead
			v3 head_back = scale(info.axis_a, -0.08f * sign);
			render_debug_line(arrow_end, add(arrow_end, add(head_back, scale(t1, 0.04f))), motor_col);
			render_debug_line(arrow_end, add(arrow_end, add(head_back, scale(t1, -0.04f))), motor_col);
		}
	}
}
