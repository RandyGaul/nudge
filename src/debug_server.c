// debug_server.c -- Debug driver + type reflection for nudge physics engine.
// Tiny TCP server for remote control + reflection data export.
// The viewer reads memory via ReadProcessMemory; this file provides the layout.
//
// Included from main.c (unity build).
// Winsock headers included from main.c before SDL (ordering requirement).

#define DBG_PORT_DEFAULT 9999
#define DBG_MAX_CLIENTS 4
#define DBG_RECV_BUF 4096

// Host-provided state. In app mode (nudge.exe, NUDGE_HOST_APP defined) these
// track the main scene. In test mode, tests publish the current World via
// dbg_set_world() and the server exposes only inspection + break control.
static World g_dbg_world;     // current world (paused or live)
static int g_dbg_port = DBG_PORT_DEFAULT;

// Breakpoint state -- shared by app + test mode. Host code calls DBG_BREAK()
// which routes through dbg_break_wait_impl() when g_dbg_break_enabled is set.
// Published so clients can see where the engine is paused.
volatile int g_dbg_break_enabled;            // set to 1 by host when --debug active
static const char* g_dbg_break_filter = "*"; // glob-ish pattern, "*" = all
static const char* g_dbg_break_name;         // name of current (active) break, NULL if running
static const char* g_dbg_break_file;
static int g_dbg_break_line;
static volatile int g_dbg_break_resume;      // set to 1 by "continue" command

#ifdef NUDGE_HOST_APP
// Remote drag state (mirrors the mouse constraint but driven via TCP).
static Body g_drag_body;
static Body g_drag_anchor;
static Joint g_drag_joint;
static v3 g_drag_local_hit;
#endif

// ============================================================================
// Reflection: type descriptors built from real engine structs via offsetof.
// ============================================================================

typedef enum R_FieldKind
{
	RFK_INT, RFK_UINT, RFK_FLOAT, RFK_DOUBLE, RFK_BOOL,
	RFK_U8, RFK_U16, RFK_U64,
	RFK_V3, RFK_QUAT,
	RFK_ENUM,
	RFK_PTR, RFK_ARRAY, RFK_MAP,
	RFK_STRUCT,
	RFK_PTR_RAW,
} R_FieldKind;

static const char *rfk_names[] = {
	"INT","UINT","FLOAT","DOUBLE","BOOL","U8","U16","U64",
	"V3","QUAT","ENUM","PTR","ARRAY","MAP","STRUCT","PTR_RAW",
};

typedef struct R_EnumEntry { const char *name; int value; } R_EnumEntry;

typedef struct R_FieldDesc
{
	const char *name;
	R_FieldKind kind;
	int offset;
	int size;
	const char *nested_type;         // type name for PTR/ARRAY/MAP/STRUCT
	const char *enum_name;           // enum table name for ENUM
} R_FieldDesc;

typedef struct R_TypeDesc
{
	const char *name;
	int size;
	int field_count;
	const R_FieldDesc *fields;
} R_TypeDesc;

typedef struct R_EnumTable
{
	const char *name;
	int count;
	const R_EnumEntry *entries;
} R_EnumTable;

// --- Field macros (use real offsetof on real engine types) ---
#define RF_INT(T, m)        { #m, RFK_INT,     offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, NULL }
#define RF_UINT(T, m)       { #m, RFK_UINT,    offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, NULL }
#define RF_FLOAT(T, m)      { #m, RFK_FLOAT,   offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, NULL }
#define RF_DOUBLE(T, m)     { #m, RFK_DOUBLE,  offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, NULL }
#define RF_BOOL(T, m)       { #m, RFK_BOOL,    offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, NULL }
#define RF_U64(T, m)        { #m, RFK_U64,     offsetof(T, m), 8, NULL, NULL }
#define RF_V3(T, m)         { #m, RFK_V3,      offsetof(T, m), 16, NULL, NULL }
#define RF_QUAT(T, m)       { #m, RFK_QUAT,    offsetof(T, m), 16, NULL, NULL }
#define RF_ENUM(T, m, e)    { #m, RFK_ENUM,    offsetof(T, m), (int)sizeof(((T*)0)->m), NULL, #e }
#define RF_PTR(T, m, nt)    { #m, RFK_PTR,     offsetof(T, m), 8, #nt, NULL }
#define RF_ARRAY(T, m, nt)  { #m, RFK_ARRAY,   offsetof(T, m), 8, #nt, NULL }
#define RF_MAP(T, m, nt)    { #m, RFK_MAP,     offsetof(T, m), 8, #nt, NULL }
#define RF_STRUCT(T, m, nt) { #m, RFK_STRUCT,  offsetof(T, m), (int)sizeof(((T*)0)->m), #nt, NULL }
#define RF_PTR_RAW(T, m)    { #m, RFK_PTR_RAW, offsetof(T, m), 8, NULL, NULL }

#define REFLECT(T, ...) \
	static const R_FieldDesc rfields_##T[] = { __VA_ARGS__ }; \
	static const R_TypeDesc rtype_##T = { #T, (int)sizeof(T), (int)(sizeof(rfields_##T)/sizeof(rfields_##T[0])), rfields_##T }

// --- Enum tables ---

static const R_EnumEntry renum_ShapeType_entries[] = { {"sphere",0}, {"capsule",1}, {"box",2}, {"hull",3}, {"cylinder",4}, {"mesh",5} };
static const R_EnumEntry renum_SolverType_entries[] = { {"soft_step",0}, {"si_soft",1}, {"si",2} };
static const R_EnumEntry renum_BroadphaseType_entries[] = { {"n2",0}, {"bvh",1} };
static const R_EnumEntry renum_JointType_entries[] = {
	{"ball_socket",0}, {"distance",1}, {"hinge",2}, {"fixed",3}, {"prismatic",4},
	{"angular_motor",5}, {"twist_limit",6}, {"cone_limit",7}, {"swing_twist",8}
};
static const R_EnumEntry renum_NarrowphaseBackend_entries[] = { {"sat",0}, {"gjk_epa",1} };

static const R_EnumTable g_enum_tables[] = {
	{ "ShapeType", 6, renum_ShapeType_entries },
	{ "SolverType", 3, renum_SolverType_entries },
	{ "BroadphaseType", 2, renum_BroadphaseType_entries },
	{ "JointType", 9, renum_JointType_entries },
	{ "NarrowphaseBackend", 2, renum_NarrowphaseBackend_entries },
};
#define NUM_ENUM_TABLES (int)(sizeof(g_enum_tables)/sizeof(g_enum_tables[0]))

// --- Type registrations (using REAL engine types, REAL offsetof) ---

REFLECT(ShapeInternal,
	RF_ENUM(ShapeInternal, type, ShapeType),
	RF_V3(ShapeInternal, local_pos),
	// Union: expose common subfields at known offsets
	// sphere.radius and capsule.half_height are at the same offset (start of union)
	{ "radius_or_half_height", RFK_FLOAT, offsetof(ShapeInternal, sphere.radius), 4, NULL, NULL },
	{ "capsule_radius", RFK_FLOAT, offsetof(ShapeInternal, capsule.radius), 4, NULL, NULL },
	{ "box_half_extents", RFK_V3, offsetof(ShapeInternal, box.half_extents), 16, NULL, NULL },
);

REFLECT(BodyCold,
	RF_FLOAT(BodyCold, mass),
	RF_ARRAY(BodyCold, shapes, ShapeInternal),
	RF_INT(BodyCold, bvh_leaf),
	RF_INT(BodyCold, island_id),
	RF_INT(BodyCold, island_prev),
	RF_INT(BodyCold, island_next),
	RF_UINT(BodyCold, collision_group),
	RF_UINT(BodyCold, collision_mask),
	RF_UINT(BodyCold, compound_id),
);

REFLECT(BodyHot,
	RF_V3(BodyHot, velocity),
	RF_V3(BodyHot, angular_velocity),
	RF_V3(BodyHot, iw_diag),
	RF_V3(BodyHot, iw_off),
	RF_FLOAT(BodyHot, inv_mass),
);

REFLECT(BodyState,
	RF_V3(BodyState, position),
	RF_QUAT(BodyState, rotation),
	RF_V3(BodyState, inv_inertia_local),
	RF_FLOAT(BodyState, friction),
	RF_FLOAT(BodyState, restitution),
	RF_FLOAT(BodyState, linear_damping),
	RF_FLOAT(BodyState, angular_damping),
	RF_FLOAT(BodyState, sleep_time),
	RF_INT(BodyState, sleep_allowed),
);

REFLECT(Contact,
	RF_V3(Contact, point),
	RF_V3(Contact, normal),
	RF_FLOAT(Contact, penetration),
	RF_UINT(Contact, feature_id),
);

REFLECT(SolverContact,
	RF_V3(SolverContact, r_a),
	RF_V3(SolverContact, r_b),
	RF_V3(SolverContact, rn_a),
	RF_V3(SolverContact, rn_b),
	RF_FLOAT(SolverContact, eff_mass_n),
	RF_FLOAT(SolverContact, bias),
	RF_FLOAT(SolverContact, bounce),
	RF_FLOAT(SolverContact, softness),
	RF_FLOAT(SolverContact, lambda_n),
	RF_V3(SolverContact, normal),
	RF_FLOAT(SolverContact, penetration),
	RF_UINT(SolverContact, feature_id),
);

REFLECT(SolverManifold,
	RF_INT(SolverManifold, body_a),
	RF_INT(SolverManifold, body_b),
	RF_INT(SolverManifold, contact_start),
	RF_INT(SolverManifold, contact_count),
	RF_FLOAT(SolverManifold, friction),
	RF_FLOAT(SolverManifold, inv_mass_a),
	RF_FLOAT(SolverManifold, inv_mass_b),
	RF_V3(SolverManifold, centroid_r_a),
	RF_V3(SolverManifold, centroid_r_b),
	RF_V3(SolverManifold, tangent1),
	RF_V3(SolverManifold, tangent2),
	RF_V3(SolverManifold, normal),
	RF_FLOAT(SolverManifold, eff_mass_t1),
	RF_FLOAT(SolverManifold, eff_mass_t2),
	RF_FLOAT(SolverManifold, eff_mass_twist),
	RF_FLOAT(SolverManifold, lambda_t1),
	RF_FLOAT(SolverManifold, lambda_t2),
	RF_FLOAT(SolverManifold, lambda_twist),
	RF_FLOAT(SolverManifold, patch_area),
	RF_FLOAT(SolverManifold, patch_radius),
);

REFLECT(WarmContact,
	RF_V3(WarmContact, r_a),
	RF_FLOAT(WarmContact, lambda_n),
	RF_UINT(WarmContact, feature_id),
);

REFLECT(WarmManifold,
	RF_INT(WarmManifold, count),
	RF_INT(WarmManifold, stale),
	RF_INT(WarmManifold, body_a),
	RF_INT(WarmManifold, body_b),
	RF_FLOAT(WarmManifold, manifold_lambda_t1),
	RF_FLOAT(WarmManifold, manifold_lambda_t2),
	RF_FLOAT(WarmManifold, manifold_lambda_twist),
	RF_INT(WarmManifold, sat_axis),
	// contacts[] is a fixed-size array of WarmContact; reflection lacks
	// fixed-array support so we expose each slot at its flat offset.
	{ "contact0_r_a",      RFK_V3,    offsetof(WarmManifold, contacts[0].r_a),        16, NULL, NULL },
	{ "contact0_lambda_n", RFK_FLOAT, offsetof(WarmManifold, contacts[0].lambda_n),    4, NULL, NULL },
	{ "contact0_feature",  RFK_UINT,  offsetof(WarmManifold, contacts[0].feature_id),  4, NULL, NULL },
	{ "contact1_r_a",      RFK_V3,    offsetof(WarmManifold, contacts[1].r_a),        16, NULL, NULL },
	{ "contact1_lambda_n", RFK_FLOAT, offsetof(WarmManifold, contacts[1].lambda_n),    4, NULL, NULL },
	{ "contact1_feature",  RFK_UINT,  offsetof(WarmManifold, contacts[1].feature_id),  4, NULL, NULL },
	{ "contact2_r_a",      RFK_V3,    offsetof(WarmManifold, contacts[2].r_a),        16, NULL, NULL },
	{ "contact2_lambda_n", RFK_FLOAT, offsetof(WarmManifold, contacts[2].lambda_n),    4, NULL, NULL },
	{ "contact2_feature",  RFK_UINT,  offsetof(WarmManifold, contacts[2].feature_id),  4, NULL, NULL },
	{ "contact3_r_a",      RFK_V3,    offsetof(WarmManifold, contacts[3].r_a),        16, NULL, NULL },
	{ "contact3_lambda_n", RFK_FLOAT, offsetof(WarmManifold, contacts[3].lambda_n),    4, NULL, NULL },
	{ "contact3_feature",  RFK_UINT,  offsetof(WarmManifold, contacts[3].feature_id),  4, NULL, NULL },
);

REFLECT(JointHot,
	{ "warm_lambda0", RFK_FLOAT, offsetof(JointHot, warm_lambda[0]), 4, NULL, NULL },
	{ "warm_lambda1", RFK_FLOAT, offsetof(JointHot, warm_lambda[1]), 4, NULL, NULL },
	{ "warm_lambda2", RFK_FLOAT, offsetof(JointHot, warm_lambda[2]), 4, NULL, NULL },
	{ "warm_lambda3", RFK_FLOAT, offsetof(JointHot, warm_lambda[3]), 4, NULL, NULL },
	{ "warm_lambda4", RFK_FLOAT, offsetof(JointHot, warm_lambda[4]), 4, NULL, NULL },
	{ "warm_lambda5", RFK_FLOAT, offsetof(JointHot, warm_lambda[5]), 4, NULL, NULL },
);

REFLECT(EpaContact,
	RF_V3(EpaContact, point_a_local),
	RF_V3(EpaContact, point_b_local),
	RF_V3(EpaContact, normal_local_a),
	RF_FLOAT(EpaContact, penetration),
	RF_UINT(EpaContact, feature_id),
	RF_INT(EpaContact, age),
);

REFLECT(EpaManifold,
	RF_INT(EpaManifold, count),
	RF_INT(EpaManifold, stale),
	RF_INT(EpaManifold, warm_valid),
);

REFLECT(EpaStats,
	RF_INT(EpaStats, queries),
	RF_INT(EpaStats, iter_cap_hits),
	RF_INT(EpaStats, total_iters),
	RF_INT(EpaStats, warm_reseeds),
	RF_INT(EpaStats, contacts_emitted),
	RF_INT(EpaStats, pair_count),
);

REFLECT(JointInternal,
	RF_ENUM(JointInternal, type, JointType),
	RF_INT(JointInternal, body_a),
	RF_INT(JointInternal, body_b),
	RF_INT(JointInternal, island_id),
	RF_INT(JointInternal, island_prev),
	RF_INT(JointInternal, island_next),
);

REFLECT(Island,
	RF_INT(Island, head_body),
	RF_INT(Island, tail_body),
	RF_INT(Island, body_count),
	RF_INT(Island, head_joint),
	RF_INT(Island, tail_joint),
	RF_INT(Island, joint_count),
	RF_INT(Island, constraint_remove_count),
	RF_INT(Island, awake),
);

REFLECT(BVH_Child,
	RF_V3(BVH_Child, min),
	RF_INT(BVH_Child, index),
	RF_V3(BVH_Child, max),
	RF_INT(BVH_Child, leaf_count),
);

REFLECT(BVHNode,
	RF_STRUCT(BVHNode, a, BVH_Child),
	RF_STRUCT(BVHNode, b, BVH_Child),
);

REFLECT(BVHLeaf,
	RF_INT(BVHLeaf, body_idx),
	RF_INT(BVHLeaf, node_idx),
	RF_INT(BVHLeaf, child_slot),
	RF_V3(BVHLeaf, fat_min),
	RF_V3(BVHLeaf, fat_max),
);

REFLECT(BVH_Tree,
	RF_ARRAY(BVH_Tree, nodes, BVHNode),
	RF_ARRAY(BVH_Tree, leaves, BVHLeaf),
	RF_INT(BVH_Tree, root),
	RF_INT(BVH_Tree, refine_cursor),
);

REFLECT(PGSTimers,
	RF_DOUBLE(PGSTimers, pre_solve),
	RF_DOUBLE(PGSTimers, warm_start),
	RF_DOUBLE(PGSTimers, graph_color),
	RF_DOUBLE(PGSTimers, iterations),
	RF_DOUBLE(PGSTimers, joint_limits),
	RF_DOUBLE(PGSTimers, ldl),
	RF_DOUBLE(PGSTimers, relax),
	RF_DOUBLE(PGSTimers, pos_contacts),
	RF_DOUBLE(PGSTimers, pos_joints),
	RF_DOUBLE(PGSTimers, post_solve),
);

REFLECT(SolverJoint,
	RF_INT(SolverJoint, body_a),
	RF_INT(SolverJoint, body_b),
	RF_INT(SolverJoint, joint_idx),
	RF_ENUM(SolverJoint, type, JointType),
	RF_INT(SolverJoint, dof),
	RF_FLOAT(SolverJoint, softness),
	RF_V3(SolverJoint, r_a),
	RF_V3(SolverJoint, r_b),
);

REFLECT(PerfTimers,
	RF_DOUBLE(PerfTimers, broadphase),
	RF_DOUBLE(PerfTimers, pre_solve),
	RF_DOUBLE(PerfTimers, pgs_solve),
	RF_DOUBLE(PerfTimers, position_correct),
	RF_DOUBLE(PerfTimers, integrate),
	RF_DOUBLE(PerfTimers, islands),
	RF_DOUBLE(PerfTimers, total),
	RF_STRUCT(PerfTimers, pgs, PGSTimers),
);

REFLECT(NP_DebugSnapshot,
	RF_INT(NP_DebugSnapshot, valid),
	RF_INT(NP_DebugSnapshot, frame),
	RF_INT(NP_DebugSnapshot, body_a),
	RF_INT(NP_DebugSnapshot, body_b),
	RF_ENUM(NP_DebugSnapshot, shape_a_kind, ShapeType),
	RF_ENUM(NP_DebugSnapshot, shape_b_kind, ShapeType),
	RF_INT(NP_DebugSnapshot, face_a_index),
	RF_FLOAT(NP_DebugSnapshot, face_a_sep),
	RF_INT(NP_DebugSnapshot, face_b_index),
	RF_FLOAT(NP_DebugSnapshot, face_b_sep),
	RF_INT(NP_DebugSnapshot, edge_a_index),
	RF_INT(NP_DebugSnapshot, edge_b_index),
	RF_FLOAT(NP_DebugSnapshot, edge_sep),
	RF_INT(NP_DebugSnapshot, winning_axis),
	RF_V3(NP_DebugSnapshot, contact_normal),
	RF_INT(NP_DebugSnapshot, contact_count),
	// contact_points / contact_pens are fixed-size C arrays. The reflection
	// system doesn't currently describe fixed arrays, so fall back to raw
	// field offsets so the viewer can read them via `raw`.
	{ "contact_points", RFK_V3, offsetof(NP_DebugSnapshot, contact_points[0]), 16, NULL, NULL },
	{ "contact_points_1", RFK_V3, offsetof(NP_DebugSnapshot, contact_points[1]), 16, NULL, NULL },
	{ "contact_points_2", RFK_V3, offsetof(NP_DebugSnapshot, contact_points[2]), 16, NULL, NULL },
	{ "contact_points_3", RFK_V3, offsetof(NP_DebugSnapshot, contact_points[3]), 16, NULL, NULL },
	{ "contact_pens_0", RFK_FLOAT, offsetof(NP_DebugSnapshot, contact_pens[0]), 4, NULL, NULL },
	{ "contact_pens_1", RFK_FLOAT, offsetof(NP_DebugSnapshot, contact_pens[1]), 4, NULL, NULL },
	{ "contact_pens_2", RFK_FLOAT, offsetof(NP_DebugSnapshot, contact_pens[2]), 4, NULL, NULL },
	{ "contact_pens_3", RFK_FLOAT, offsetof(NP_DebugSnapshot, contact_pens[3]), 4, NULL, NULL },
);

REFLECT(RewindIsland,
	RF_INT(RewindIsland, head_body),
	RF_INT(RewindIsland, tail_body),
	RF_INT(RewindIsland, body_count),
	RF_INT(RewindIsland, head_joint),
	RF_INT(RewindIsland, tail_joint),
	RF_INT(RewindIsland, joint_count),
	RF_INT(RewindIsland, constraint_remove_count),
	RF_INT(RewindIsland, awake),
);

REFLECT(RewindFrame,
	RF_U64(RewindFrame, frame_id),
	RF_INT(RewindFrame, sim_frame),
	RF_INT(RewindFrame, n_bodies),
	RF_INT(RewindFrame, is_keyframe),
	RF_INT(RewindFrame, n_dirty),
	RF_PTR_RAW(RewindFrame, dirty_indices),
	RF_PTR_RAW(RewindFrame, dirty_payload),
	RF_INT(RewindFrame, dirty_payload_size),
	RF_PTR_RAW(RewindFrame, dirty_payload_offsets),
	RF_PTR_RAW(RewindFrame, dirty_shapes),
	RF_PTR_RAW(RewindFrame, body_gen),
	RF_PTR_RAW(RewindFrame, body_free),
	RF_INT(RewindFrame, n_joints),
	RF_PTR_RAW(RewindFrame, joints),
	RF_PTR_RAW(RewindFrame, joint_hot),
	RF_PTR_RAW(RewindFrame, joint_gen),
	RF_PTR_RAW(RewindFrame, joint_free),
	RF_INT(RewindFrame, n_islands),
	RF_PTR_RAW(RewindFrame, islands),
	RF_PTR_RAW(RewindFrame, island_gen),
	RF_PTR_RAW(RewindFrame, island_free),
	RF_INT(RewindFrame, n_warm),
	RF_PTR_RAW(RewindFrame, warm_keys),
	RF_PTR_RAW(RewindFrame, warm_vals),
	RF_INT(RewindFrame, n_prev_touching),
	RF_INT(RewindFrame, n_joint_pairs),
	RF_INT(RewindFrame, joint_pairs_version),
	RF_INT(RewindFrame, ldl_topo_version),
);

REFLECT(RewindBuffer,
	RF_INT(RewindBuffer, max_frames),
	RF_INT(RewindBuffer, auto_capture),
	RF_U64(RewindBuffer, next_frame_id),
	RF_ARRAY(RewindBuffer, frames, RewindFrame),
	RF_INT(RewindBuffer, head),
	RF_INT(RewindBuffer, count),
	RF_INT(RewindBuffer, baseline_n_bodies),
	RF_PTR(RewindBuffer, baseline_hot, BodyHot),
	RF_PTR(RewindBuffer, baseline_state, BodyState),
	RF_PTR(RewindBuffer, baseline_cold, BodyCold),
	RF_PTR_RAW(RewindBuffer, baseline_shapes),
);

REFLECT(WorldInternal,
	RF_INT(WorldInternal, frame),
	RF_V3(WorldInternal, gravity),
	RF_ARRAY(WorldInternal, body_cold, BodyCold),
	RF_ARRAY(WorldInternal, body_hot, BodyHot),
	RF_ARRAY(WorldInternal, body_state, BodyState),
	RF_ARRAY(WorldInternal, body_gen, uint32_t),
	RF_PTR_RAW(WorldInternal, body_free),
	RF_ARRAY(WorldInternal, debug_contacts, Contact),
	RF_MAP(WorldInternal, warm_cache, WarmManifold),
	RF_ARRAY(WorldInternal, joints, JointInternal),
	RF_ARRAY(WorldInternal, joint_hot, JointHot),
	RF_ARRAY(WorldInternal, joint_gen, uint32_t),
	RF_PTR_RAW(WorldInternal, joint_free),
	RF_ENUM(WorldInternal, broadphase_type, BroadphaseType),
	RF_PTR(WorldInternal, bvh_static, BVH_Tree),
	RF_PTR(WorldInternal, bvh_dynamic, BVH_Tree),
	RF_PTR(WorldInternal, bvh_sleeping, BVH_Tree),
	RF_ARRAY(WorldInternal, islands, Island),
	RF_ARRAY(WorldInternal, island_gen, uint32_t),
	RF_PTR_RAW(WorldInternal, island_free),
	RF_INT(WorldInternal, joint_pairs_version),
	RF_INT(WorldInternal, sleep_enabled),
	RF_INT(WorldInternal, sat_hint_enabled),
	RF_INT(WorldInternal, sat_hillclimb_enabled),
	RF_INT(WorldInternal, box_use_hull),
	RF_INT(WorldInternal, incremental_np_enabled),
	RF_INT(WorldInternal, warm_start_enabled),
	RF_INT(WorldInternal, trimesh_simd_enabled),
	RF_ENUM(WorldInternal, narrowphase_backend, NarrowphaseBackend),
	RF_MAP(WorldInternal, epa_cache, EpaManifold),
	RF_STRUCT(WorldInternal, epa_stats, EpaStats),
	RF_INT(WorldInternal, cyl_native_sphere),
	RF_INT(WorldInternal, cyl_native_capsule),
	RF_INT(WorldInternal, cyl_native_box),
	RF_INT(WorldInternal, cyl_native_hull),
	RF_INT(WorldInternal, cyl_native_cyl),
	RF_INT(WorldInternal, thread_count),
	RF_INT(WorldInternal, ldl_enabled),
	RF_INT(WorldInternal, ldl_topo_version),
	RF_INT(WorldInternal, ldl_correction_iter),
	RF_ENUM(WorldInternal, solver_type, SolverType),
	RF_INT(WorldInternal, velocity_iters),
	RF_INT(WorldInternal, position_iters),
	RF_FLOAT(WorldInternal, contact_hertz),
	RF_FLOAT(WorldInternal, contact_damping_ratio),
	RF_FLOAT(WorldInternal, max_push_velocity),
	RF_INT(WorldInternal, sub_steps),
	RF_STRUCT(WorldInternal, perf, PerfTimers),
	RF_PTR_RAW(WorldInternal, dbg_solver_manifolds),
	RF_PTR_RAW(WorldInternal, dbg_solver_contacts),
	RF_PTR_RAW(WorldInternal, dbg_solver_joints),
	RF_INT(WorldInternal, np_debug_enabled),
	RF_INT(WorldInternal, np_debug_filter_body_a),
	RF_INT(WorldInternal, np_debug_filter_body_b),
	RF_STRUCT(WorldInternal, np_debug, NP_DebugSnapshot),
	RF_PTR(WorldInternal, rewind, RewindBuffer),
);

#ifdef _WIN32
// Forward declarations for send helpers (defined below).
typedef struct DbgClient DbgClient;
static int dbg_send_str(DbgClient *c, const char *s);

// Global type registry (all types the viewer can navigate).
static const R_TypeDesc *g_type_registry[] = {
	&rtype_WorldInternal,
	&rtype_BodyCold, &rtype_BodyHot, &rtype_BodyState,
	&rtype_ShapeInternal, &rtype_Contact,
	&rtype_SolverContact, &rtype_SolverManifold, &rtype_SolverJoint,
	&rtype_WarmContact, &rtype_WarmManifold,
	&rtype_JointInternal, &rtype_JointHot, &rtype_Island,
	&rtype_BVH_Child, &rtype_BVHNode, &rtype_BVHLeaf, &rtype_BVH_Tree,
	&rtype_PerfTimers, &rtype_PGSTimers, &rtype_NP_DebugSnapshot,
	&rtype_EpaContact, &rtype_EpaManifold, &rtype_EpaStats,
	&rtype_RewindIsland, &rtype_RewindFrame, &rtype_RewindBuffer,
};
#define NUM_TYPES (int)(sizeof(g_type_registry)/sizeof(g_type_registry[0]))

// Serialize all types as text for the viewer.
static void dbg_send_types(DbgClient *c)
{
	CK_SDYNA char *s = NULL;
	sfmt(s, "TYPES %d %d\n", NUM_TYPES, NUM_ENUM_TABLES);

	for (int t = 0; t < NUM_TYPES; t++) {
		const R_TypeDesc *td = g_type_registry[t];
		sfmt_append(s, "TYPE %s %d\n", td->name, td->size);
		for (int f = 0; f < td->field_count; f++) {
			const R_FieldDesc *fd = &td->fields[f];
			sfmt_append(s, "FIELD %s %s %d %d", fd->name, rfk_names[fd->kind], fd->offset, fd->size);
			if (fd->nested_type) sfmt_append(s, " %s", fd->nested_type);
			if (fd->enum_name) sfmt_append(s, " %s", fd->enum_name);
			sappend(s, "\n");
		}
		sappend(s, "END\n");
	}

	for (int e = 0; e < NUM_ENUM_TABLES; e++) {
		const R_EnumTable *et = &g_enum_tables[e];
		sfmt_append(s, "ENUM %s %d\n", et->name, et->count);
		for (int v = 0; v < et->count; v++)
			sfmt_append(s, "VALUE %s %d\n", et->entries[v].name, et->entries[v].value);
		sappend(s, "END\n");
	}
	sappend(s, "ENDTYPES\n");

	dbg_send_str(c, s);
	sfree(s);
}

// ============================================================================
// TCP server + command dispatch
// ============================================================================

struct DbgClient
{
	SOCKET sock;
	char recv_buf[DBG_RECV_BUF];
	int recv_len;
};

typedef struct DbgServer
{
	SOCKET listen_sock;
	DbgClient clients[DBG_MAX_CLIENTS];
	int client_count;
	bool initialized;
} DbgServer;

static DbgServer g_dbg;

static int dbg_send_str(DbgClient *c, const char *s)
{
	int len = (int)strlen(s);
	int total = 0;
	while (total < len) {
		int n = send(c->sock, s + total, len - total, 0);
		if (n == SOCKET_ERROR) return -1;
		total += n;
	}
	return total;
}

static void dbg_reply(DbgClient *c, CK_SDYNA char *s)
{
	dbg_send_str(c, s);
	sfree(s);
}

static char *str_trim(char *s)
{
	while (*s == ' ' || *s == '\t') s++;
	int len = (int)strlen(s);
	while (len > 0 && (s[len-1] == ' ' || s[len-1] == '\t' || s[len-1] == '\r' || s[len-1] == '\n')) s[--len] = '\0';
	return s;
}

static const char *str_next_word(const char *s, char *word, int word_cap)
{
	while (*s == ' ') s++;
	int i = 0;
	while (*s && *s != ' ' && i < word_cap - 1) word[i++] = *s++;
	word[i] = '\0';
	while (*s == ' ') s++;
	return s;
}

static void dbg_dispatch(DbgClient *c, char *line)
{
	line = str_trim(line);
	if (!*line) return;

	char cmd[32];
	const char *rest = str_next_word(line, cmd, sizeof(cmd));

	if (strcmp(cmd, "help") == 0) {
		dbg_send_str(c,
			"OK engine driver commands (sent via TCP):\n"
			"\n"
			"  Common (both app + test mode):\n"
			"    help                              this help\n"
			"    info                              PID, world pointer, frame, break state\n"
			"    types                             send type layout data for RPM viewer\n"
			"    continue | c                      resume from a break (no-op if not paused)\n"
			"    where                             show current break name + file:line\n"
			"    break-filter [pattern]            get/set runtime break filter (glob-ish)\n"
			"    np-debug <0|1> [body_a] [body_b]  toggle narrowphase snapshot capture (body indices: -1 = any)\n"
#ifdef NUDGE_HOST_APP
			"\n"
			"  App mode (nudge.exe):\n"
			"    pause                             toggle pause\n"
			"    step [n]                          advance n frames instantly (default 1)\n"
			"    run [n]                           run n frames at render speed, auto-pause\n"
			"    scene <name>                      load scene by name (prefix match ok)\n"
			"    scenes                            list all scene names\n"
			"    restart                           restart current scene\n"
			"    slow <factor>                     time scale: 0.1=slow-mo, 1=normal, 10=fast\n"
			"    playrecording                     replay mouse_recording.bin\n"
			"    push <idx> <fx fy fz> [rx ry rz]  impulse at point (default top of body)\n"
			"    drag <idx> <lx ly lz> <tx ty tz>  begin mouse-like drag (soft spring)\n"
			"    dragto <tx ty tz>                 update drag target position\n"
			"    release                           end drag\n"
			"    highlight <idx> <color|r g b>     tint body\n"
			"    unhighlight <idx|all>             remove tint\n"
			"    label <text>                      show overlay text for 5 seconds\n"
			"    npviz                             enter NP Viz mode\n"
			"    npvset <A|B> <kind> ...           set shape in NP Viz\n"
#else
			"\n"
			"  Test mode: simulation control is driven by the test. Use `continue`\n"
			"  to release a pause point; the test schedules the next step.\n"
#endif
		);
	} else if (strcmp(cmd, "info") == 0) {
		WorldInternal *w = g_dbg_world.id ? (WorldInternal *)g_dbg_world.id : NULL;
		CK_SDYNA char *s = NULL;
#ifdef NUDGE_HOST_APP
		sfmt(s, "OK pid=%lu world=0x%llX scene=\"%s\" frame=%d bodies=%d paused=%s\n",
			(unsigned long)GetCurrentProcessId(),
			(unsigned long long)g_dbg_world.id,
			g_scenes[g_scene_index].name,
			w ? w->frame : 0,
			asize(g_draw_list),
			g_paused ? "true" : "false");
#else
		sfmt(s, "OK pid=%lu world=0x%llX host=test frame=%d break=%s%s%s\n",
			(unsigned long)GetCurrentProcessId(),
			(unsigned long long)g_dbg_world.id,
			w ? w->frame : 0,
			g_dbg_break_name ? g_dbg_break_name : "(none)",
			g_dbg_break_name ? " at " : "",
			g_dbg_break_name ? g_dbg_break_file : "");
#endif
		dbg_reply(c, s);
	} else if (strcmp(cmd, "types") == 0) {
		dbg_send_types(c);
	} else if (strcmp(cmd, "continue") == 0 || strcmp(cmd, "c") == 0) {
		if (!g_dbg_break_name) { dbg_send_str(c, "ERR not at a break\n"); return; }
		g_dbg_break_resume = 1;
		dbg_send_str(c, "OK resuming\n");
	} else if (strcmp(cmd, "where") == 0) {
		if (!g_dbg_break_name) { dbg_send_str(c, "OK running (no active break)\n"); return; }
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK break=\"%s\" at %s:%d world=0x%llX\n",
			g_dbg_break_name, g_dbg_break_file, g_dbg_break_line,
			(unsigned long long)g_dbg_world.id);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "break-filter") == 0) {
		// Update the break filter at runtime. Empty => keep current.
		if (*rest) {
			static char filter_buf[128];
			strncpy(filter_buf, rest, sizeof(filter_buf) - 1);
			filter_buf[sizeof(filter_buf) - 1] = '\0';
			g_dbg_break_filter = filter_buf;
		}
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK filter=\"%s\" enabled=%d\n", g_dbg_break_filter ? g_dbg_break_filter : "", g_dbg_break_enabled);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "np-debug") == 0) {
		// np-debug <0|1> [filter_body_a] [filter_body_b]
		// Toggle narrowphase debug-snapshot capture on the current world.
		// -1 for either filter body = match any (narrowphase pair order-insensitive).
		if (!g_dbg_world.id) { dbg_send_str(c, "ERR no world\n"); return; }
		int on = -1, fa = -1, fb = -1;
		int nargs = sscanf(rest, "%d %d %d", &on, &fa, &fb);
		WorldInternal* w = (WorldInternal*)g_dbg_world.id;
		if (nargs >= 1) w->np_debug_enabled = on ? 1 : 0;
		if (nargs >= 2) w->np_debug_filter_body_a = fa;
		if (nargs >= 3) w->np_debug_filter_body_b = fb;
		if (w->np_debug_enabled) w->np_debug.valid = 0; // clear so next step re-captures
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK np_debug_enabled=%d filter=(%d,%d) last_valid=%d\n",
			w->np_debug_enabled, w->np_debug_filter_body_a, w->np_debug_filter_body_b, w->np_debug.valid);
		dbg_reply(c, s);
#ifdef NUDGE_HOST_APP
	} else if (strcmp(cmd, "pause") == 0) {
		g_paused = !g_paused;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK paused=%s\n", g_paused ? "true" : "false");
		dbg_reply(c, s);
	} else if (strcmp(cmd, "step") == 0) {
		int n = (*rest) ? atoi(rest) : 1;
		if (n < 1) n = 1;
		for (int i = 0; i < n; i++) world_step(g_dbg_world, 1.0f / 60.0f);
		WorldInternal *w = (WorldInternal *)g_dbg_world.id;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK frame=%d\n", w->frame);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "scene") == 0) {
		if (!*rest) { dbg_send_str(c, "ERR expected scene name\n"); return; }
		int found = -1;
		for (int i = 0; i < (int)SCENE_COUNT; i++) {
			if (strcmp(g_scenes[i].name, rest) == 0) { found = i; break; }
		}
		if (found < 0) {
			for (int i = 0; i < (int)SCENE_COUNT; i++) {
				if (strncmp(g_scenes[i].name, rest, strlen(rest)) == 0) { found = i; break; }
			}
		}
		if (found >= 0) {
			g_scene_index = found;
			setup_scene();
			WorldInternal *w = (WorldInternal *)g_dbg_world.id;
			CK_SDYNA char *s = NULL;
			sfmt(s, "OK scene=\"%s\" world=0x%llX frame=%d\n", g_scenes[found].name, (unsigned long long)g_dbg_world.id, w->frame);
			dbg_reply(c, s);
		} else {
			dbg_send_str(c, "ERR unknown scene\n");
		}
	} else if (strcmp(cmd, "scenes") == 0) {
		CK_SDYNA char *s = NULL;
		sset(s, "OK");
		for (int i = 0; i < (int)SCENE_COUNT; i++)
			sfmt_append(s, " \"%s\"", g_scenes[i].name);
		sappend(s, "\n");
		dbg_reply(c, s);
	} else if (strcmp(cmd, "restart") == 0) {
		setup_scene();
		WorldInternal *w = (WorldInternal *)g_dbg_world.id;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK world=0x%llX frame=%d\n", (unsigned long long)g_dbg_world.id, w->frame);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "push") == 0) {
		// push <body_idx> <fx> <fy> <fz> [rx ry rz]
		// Applies impulse J at point r (body-local offset, default = top of shape).
		// Linear: dv = J * inv_mass. Angular: dw = I_w^{-1} * cross(r_world, J).
		int idx; float fx, fy, fz;
		float rx = 0, ry = 0, rz = 0;
		int nargs = sscanf(rest, "%d %f %f %f %f %f %f", &idx, &fx, &fy, &fz, &rx, &ry, &rz);
		if (nargs < 4) {
			dbg_send_str(c, "ERR usage: push <body_idx> <fx> <fy> <fz> [rx ry rz]\n");
			return;
		}
		WorldInternal *w = (WorldInternal *)g_dbg_world.id;
		if (idx < 0 || idx >= asize(w->body_hot) || !(w->body_gen[idx] & 1)) {
			dbg_send_str(c, "ERR invalid body index\n");
			return;
		}
		if (body_inv_mass(w, idx) == 0.0f) {
			dbg_send_str(c, "ERR body is static\n");
			return;
		}
		v3 J = V3(fx, fy, fz);
		// Default application point: top of body (local Y = 0.5)
		v3 r_local = (nargs >= 7) ? V3(rx, ry, rz) : V3(0, 0.5f, 0);
		v3 r_world = quat_rotate(body_rot(w, idx), r_local);
		// Apply linear impulse
		body_vel(w, idx) = v3_add(body_vel(w, idx), v3_scale(J, body_inv_mass(w, idx)));
		// Apply angular impulse: dw = I_w^{-1} * cross(r, J)
		v3 torque_impulse = v3_cross(r_world, J);
		BodyState *st = &w->body_state[idx];
		v3 dw = inv_inertia_mul(st->rotation, st->inv_inertia_local, torque_impulse);
		body_angvel(w, idx) = v3_add(body_angvel(w, idx), dw);
		// Wake the body's island
		int isl = w->body_cold[idx].island_id;
		if (isl >= 0 && !w->islands[isl].awake) island_wake(w, isl);
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK pushed body %d: J=(%.1f,%.1f,%.1f) at r=(%.2f,%.2f,%.2f)\n", idx, fx, fy, fz, r_world.x, r_world.y, r_world.z);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "drag") == 0) {
		// drag <body_idx> <local_x> <local_y> <local_z> <target_x> <target_y> <target_z>
		// Begin dragging: creates soft ball-socket from body to a hidden anchor.
		int idx; float lx, ly, lz, tx, ty, tz;
		if (sscanf(rest, "%d %f %f %f %f %f %f", &idx, &lx, &ly, &lz, &tx, &ty, &tz) != 7) {
			dbg_send_str(c, "ERR usage: drag <body_idx> <local_xyz> <target_xyz>\n");
			return;
		}
		WorldInternal *w = (WorldInternal *)g_dbg_world.id;
		if (idx < 0 || idx >= asize(w->body_hot) || !(w->body_gen[idx] & 1)) {
			dbg_send_str(c, "ERR invalid body index\n");
			return;
		}
		if (g_drag_body.id) {
			destroy_joint(g_dbg_world, g_drag_joint);
			destroy_body(g_dbg_world, g_drag_anchor);
			g_drag_body = (Body){0};
		}
		// Find the Body handle from index
		Body target_body = {0};
		for (int i = 0; i < asize(g_draw_list); i++) {
			if (handle_index(g_draw_list[i].body) == idx) { target_body = g_draw_list[i].body; break; }
		}
		if (!target_body.id) { dbg_send_str(c, "ERR body not in draw list\n"); return; }
		g_drag_body = target_body;
		g_drag_local_hit = V3(lx, ly, lz);
		g_drag_anchor = create_body(g_dbg_world, (BodyParams){ .position = V3(tx, ty, tz), .rotation = quat_identity(), .mass = 0 });
		// Softer spring than mouse constraint (2Hz vs 5Hz) because TCP
		// commands move the target in discrete jumps, not smooth mouse motion.
		g_drag_joint = create_ball_socket(g_dbg_world, (BallSocketParams){
			.body_a = g_drag_anchor, .body_b = g_drag_body,
			.local_offset_a = V3(0, 0, 0), .local_offset_b = g_drag_local_hit,
			.spring = { .frequency = 2.0f, .damping_ratio = 1.0f },
		});
		int isl = w->body_cold[idx].island_id;
		if (isl >= 0 && !w->islands[isl].awake) island_wake(w, isl);
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK dragging body %d from local (%.2f,%.2f,%.2f) to (%.2f,%.2f,%.2f)\n", idx, lx, ly, lz, tx, ty, tz);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "dragto") == 0) {
		// dragto <target_x> <target_y> <target_z>
		// Update drag target position.
		float tx, ty, tz;
		if (sscanf(rest, "%f %f %f", &tx, &ty, &tz) != 3) {
			dbg_send_str(c, "ERR usage: dragto <target_xyz>\n");
			return;
		}
		if (!g_drag_body.id) { dbg_send_str(c, "ERR not dragging\n"); return; }
		WorldInternal *w = (WorldInternal *)g_dbg_world.id;
		int anchor_idx = handle_index(g_drag_anchor);
		body_pos(w, anchor_idx) = V3(tx, ty, tz);
		int joint_idx = handle_index(g_drag_joint);
		int isl = w->joints[joint_idx].island_id;
		if (isl >= 0 && !w->islands[isl].awake) island_wake(w, isl);
		dbg_send_str(c, "OK\n");
	} else if (strcmp(cmd, "release") == 0) {
		if (!g_drag_body.id) { dbg_send_str(c, "ERR not dragging\n"); return; }
		destroy_joint(g_dbg_world, g_drag_joint);
		destroy_body(g_dbg_world, g_drag_anchor);
		g_drag_body = (Body){0};
		g_drag_anchor = (Body){0};
		g_drag_joint = (Joint){0};
		dbg_send_str(c, "OK released\n");
	} else if (strcmp(cmd, "npviz") == 0) {
		g_npv_mode = 1;
		if (!npv_initialized) npv_init();
		dbg_send_str(c, "OK npviz mode\n");
	} else if (strcmp(cmd, "npvset") == 0) {
		// npvset <A|B> <kind> <px> <py> <pz> <qx> <qy> <qz> <qw> <radius> <hh> <hex> <hey> <hez>
		if (!npv_initialized) { g_npv_mode = 1; npv_init(); }
		char which; int kind;
		float px, py, pz, qx, qy, qz, qw, radius, hh, hex, hey, hez;
		if (sscanf(rest, "%c %d %f %f %f %f %f %f %f %f %f %f %f %f",
				&which, &kind, &px, &py, &pz, &qx, &qy, &qz, &qw, &radius, &hh, &hex, &hey, &hez) != 14) {
			dbg_send_str(c, "ERR usage: npvset <A|B> <kind> <px py pz> <qx qy qz qw> <radius> <hh> <hex hey hez>\n");
			return;
		}
		int idx = (which == 'B' || which == 'b') ? 1 : 0;
		NPV_Shape *s = &npv_shapes[idx];
		npv_free_hull(s);
		s->kind = kind;
		s->pos = V3(px, py, pz);
		s->rot = (quat){ qx, qy, qz, qw };
		s->radius = radius;
		s->half_height = hh;
		s->half_extents = V3(hex, hey, hez);
		if (s->kind >= NPV_HULL_TETRA) npv_rebuild_hull(s);
		npv_rebuild_mesh(s);
		CK_SDYNA char *r = NULL;
		sfmt(r, "OK shape %c set to kind=%d pos=(%.3f,%.3f,%.3f)\n", which, kind, px, py, pz);
		dbg_reply(c, r);
	} else if (strcmp(cmd, "playrecording") == 0) {
		FILE *rf = fopen("mouse_recording.bin", "rb");
		if (!rf) { dbg_send_str(c, "ERR no mouse_recording.bin\n"); return; }
		int rec_scene, rec_start, rec_n;
		fread(&rec_scene, sizeof(int), 1, rf);
		fread(&rec_start, sizeof(int), 1, rf);
		fread(&rec_n, sizeof(int), 1, rf);
		afree(g_recorded_frames);
		g_recorded_frames = NULL;
		afit(g_recorded_frames, rec_n);
		asetlen(g_recorded_frames, rec_n);
		fread(g_recorded_frames, sizeof(RecordedFrame), rec_n, rf);
		fclose(rf);
		g_rec_start_frame = rec_start;
		g_scene_index = rec_scene;
		setup_scene();
		WorldInternal *pw = (WorldInternal *)g_dbg_world.id;
		for (int i = pw->frame; i < rec_start; i++) world_step(g_dbg_world, 1.0f / 60.0f);
		g_recording = 2;
		g_rec_body_draw_idx = -1;
		g_playback_frame = 0;
		g_paused = false;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK playing %d frames (scene=%d start=%d) world=0x%llX\n", rec_n, rec_scene, rec_start, (unsigned long long)g_dbg_world.id);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "run") == 0) {
		int n = (*rest) ? atoi(rest) : 60;
		if (n < 1) n = 1;
		g_paused = false;
		g_run_frames = n;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK running %d frames\n", n);
		dbg_reply(c, s);
	} else if (strcmp(cmd, "highlight") == 0) {
		int idx; float r, g, b;
		if (sscanf(rest, "%d %f %f %f", &idx, &r, &g, &b) != 4) {
			// Named colors
			char color_name[32]; idx = -1;
			if (sscanf(rest, "%d %31s", &idx, color_name) == 2) {
				if (strcmp(color_name, "red") == 0) { r = 1; g = 0.2f; b = 0.2f; }
				else if (strcmp(color_name, "green") == 0) { r = 0.2f; g = 1; b = 0.2f; }
				else if (strcmp(color_name, "yellow") == 0) { r = 1; g = 1; b = 0.2f; }
				else if (strcmp(color_name, "cyan") == 0) { r = 0.2f; g = 1; b = 1; }
				else if (strcmp(color_name, "orange") == 0) { r = 1; g = 0.6f; b = 0.1f; }
				else { dbg_send_str(c, "ERR unknown color\n"); return; }
			} else { dbg_send_str(c, "ERR usage: highlight <body> <color|r g b>\n"); return; }
		}
		if (g_highlight_count < MAX_HIGHLIGHTS) {
			g_highlights[g_highlight_count++] = (typeof(g_highlights[0])){ idx, V3(r, g, b), 1 };
			dbg_send_str(c, "OK\n");
		} else dbg_send_str(c, "ERR too many highlights\n");
	} else if (strcmp(cmd, "unhighlight") == 0) {
		if (strcmp(rest, "all") == 0) { g_highlight_count = 0; dbg_send_str(c, "OK\n"); }
		else {
			int idx = atoi(rest);
			for (int i = 0; i < g_highlight_count; i++) {
				if (g_highlights[i].body_idx == idx) { g_highlights[i] = g_highlights[--g_highlight_count]; break; }
			}
			dbg_send_str(c, "OK\n");
		}
	} else if (strcmp(cmd, "label") == 0) {
		if (*rest) {
			strncpy(g_label_text, rest, sizeof(g_label_text) - 1);
			g_label_text[sizeof(g_label_text) - 1] = '\0';
			g_label_timer = 5.0f; // show for 5 seconds
			dbg_send_str(c, "OK\n");
		} else { g_label_text[0] = '\0'; g_label_timer = 0; dbg_send_str(c, "OK\n"); }
	} else if (strcmp(cmd, "slow") == 0) {
		float s = (float)atof(rest);
		if (s > 0 && s <= 10.0f) { g_time_scale = s; CK_SDYNA char *r = NULL; sfmt(r, "OK time_scale=%.2f\n", s); dbg_reply(c, r); }
		else dbg_send_str(c, "ERR expected 0 < factor <= 10\n");
#endif // NUDGE_HOST_APP
	} else {
		dbg_send_str(c, "ERR unknown command, try 'help'\n");
	}
}

// ============================================================================
// Server lifecycle
// ============================================================================

static void debug_server_init()
{
	WSADATA wsa;
	if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) return;

	g_dbg.listen_sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (g_dbg.listen_sock == INVALID_SOCKET) return;

	int opt = 1;
	setsockopt(g_dbg.listen_sock, SOL_SOCKET, SO_REUSEADDR, (const char *)&opt, sizeof(opt));

	u_long nonblock = 1;
	ioctlsocket(g_dbg.listen_sock, FIONBIO, &nonblock);

	struct sockaddr_in addr = {0};
	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = INADDR_ANY;
	addr.sin_port = htons((u_short)g_dbg_port);

	if (bind(g_dbg.listen_sock, (struct sockaddr *)&addr, sizeof(addr)) == SOCKET_ERROR) {
		closesocket(g_dbg.listen_sock); g_dbg.listen_sock = INVALID_SOCKET; return;
	}
	if (listen(g_dbg.listen_sock, 4) == SOCKET_ERROR) {
		closesocket(g_dbg.listen_sock); g_dbg.listen_sock = INVALID_SOCKET; return;
	}

	g_dbg.initialized = true;
	printf("[dbg] listening on port %d (PID %lu, world=0x%llX)\n",
		g_dbg_port, (unsigned long)GetCurrentProcessId(), (unsigned long long)g_dbg_world.id);
}

static void debug_server_poll()
{
	if (!g_dbg.initialized) return;

	while (g_dbg.client_count < DBG_MAX_CLIENTS) {
		SOCKET s = accept(g_dbg.listen_sock, NULL, NULL);
		if (s == INVALID_SOCKET) break;
		u_long nonblock = 1;
		ioctlsocket(s, FIONBIO, &nonblock);
		DbgClient *c = &g_dbg.clients[g_dbg.client_count++];
		memset(c, 0, sizeof(*c));
		c->sock = s;
		// Send banner + types automatically on connect
		CK_SDYNA char *banner = NULL;
		sfmt(banner, "nudge pid=%lu world=0x%llX\n",
			(unsigned long)GetCurrentProcessId(), (unsigned long long)g_dbg_world.id);
		dbg_send_str(c, banner);
		sfree(banner);
		dbg_send_types(c);
	}

	for (int i = 0; i < g_dbg.client_count; i++) {
		DbgClient *c = &g_dbg.clients[i];
		int space = DBG_RECV_BUF - c->recv_len - 1;
		if (space > 0) {
			int n = recv(c->sock, c->recv_buf + c->recv_len, space, 0);
			if (n > 0) {
				c->recv_len += n;
				c->recv_buf[c->recv_len] = '\0';
			} else if (n == 0 || (n == SOCKET_ERROR && WSAGetLastError() != WSAEWOULDBLOCK)) {
				closesocket(c->sock);
				g_dbg.clients[i] = g_dbg.clients[--g_dbg.client_count];
				i--;
				continue;
			}
		}
		char *start = c->recv_buf;
		char *nl;
		while ((nl = strchr(start, '\n')) != NULL) {
			*nl = '\0';
			dbg_dispatch(c, start);
			start = nl + 1;
		}
		int remaining = c->recv_len - (int)(start - c->recv_buf);
		if (remaining > 0 && start != c->recv_buf)
			memmove(c->recv_buf, start, remaining);
		c->recv_len = remaining;
	}
}

static void debug_server_shutdown()
{
	if (!g_dbg.initialized) return;
	for (int i = 0; i < g_dbg.client_count; i++)
		closesocket(g_dbg.clients[i].sock);
	g_dbg.client_count = 0;
	closesocket(g_dbg.listen_sock);
	WSACleanup();
	g_dbg.initialized = false;
}

// ============================================================================
// Host-facing setters / breakpoint API.
// ============================================================================

// Host publishes the World to inspect. In app mode, nudge.exe calls this once
// with the active world. In test mode, the DBG_BREAK macro calls it right
// before pausing so each test's World is surfaced to the viewer.
static void debug_server_set_world(World w) { g_dbg_world = w; }
static void debug_server_set_port(int port) { g_dbg_port = port > 0 ? port : DBG_PORT_DEFAULT; }
static void debug_server_set_break_filter(const char* pattern) { g_dbg_break_filter = pattern ? pattern : "*"; }

// Glob match with '*' wildcard (matches any run of characters, including empty).
// Simple greedy backtrack -- patterns are short (a few chars) so recursion depth is fine.
static int dbg_glob_match(const char* pat, const char* s)
{
	if (!pat) return 1;
	if (*pat == '\0') return *s == '\0';
	if (*pat == '*') {
		// Skip consecutive stars.
		while (*pat == '*') pat++;
		if (*pat == '\0') return 1;
		for (const char* p = s; *p; p++) if (dbg_glob_match(pat, p)) return 1;
		return dbg_glob_match(pat, s + strlen(s));
	}
	if (*s == '\0') return 0;
	if (*pat != *s) return 0;
	return dbg_glob_match(pat + 1, s + 1);
}

// Returns 1 if `name` matches the current break filter. Host code can call this
// before building expensive diagnostic state so breaks are truly zero-cost when
// the filter excludes them.
static int dbg_break_match(const char* name)
{
	if (!g_dbg_break_enabled) return 0;
	return dbg_glob_match(g_dbg_break_filter ? g_dbg_break_filter : "*", name);
}

// Actual implementation of DBG_BREAK. Publishes the World and break metadata,
// then polls the server until a `continue` command arrives. While paused, a
// connected viewer can issue inspection commands (get/summary/table/etc.)
// against the published World -- exactly like hitting a debugger breakpoint.
static void dbg_break_wait_impl(const char* name, const char* file, int line, World w)
{
	if (!g_dbg.initialized) return; // not running server -> no-op

	g_dbg_world = w;
	g_dbg_break_name = name;
	g_dbg_break_file = file;
	g_dbg_break_line = line;
	g_dbg_break_resume = 0;

	// Log locally so a user running the test without an attached viewer can
	// see exactly where it stopped. Async TCP pushes mix badly with the viewer's
	// command/response model, so we leave break-state discovery to the `where`
	// and `info` commands instead.
	fprintf(stderr, "[dbg] break \"%s\" at %s:%d (world=0x%llX). Connect viewer and send `continue`.\n",
		name, file, line, (unsigned long long)w.id);
	fflush(stderr);

	// Pump the TCP event loop until a viewer sends `continue`. Sleep briefly so
	// we don't spin at 100% CPU while waiting.
	while (!g_dbg_break_resume) {
		debug_server_poll();
		Sleep(5);
	}

	g_dbg_break_name = NULL;
	g_dbg_break_file = NULL;
	g_dbg_break_line = 0;
}

#else  // !_WIN32 -- debug server + break system are Windows-only for now.
static void debug_server_set_world(World w) { (void)w; g_dbg_world = w; }
static void debug_server_set_port(int port) { g_dbg_port = port > 0 ? port : DBG_PORT_DEFAULT; }
static void debug_server_set_break_filter(const char* p) { g_dbg_break_filter = p ? p : "*"; }
static void debug_server_init() {}
static void debug_server_shutdown() {}
static int dbg_break_match(const char* name) { (void)name; return 0; }
static void dbg_break_wait_impl(const char* name, const char* file, int line, World w) { (void)name; (void)file; (void)line; (void)w; }
#endif

// Pause-point macro for use in tests or engine code. Compiles to a single
// predicate check + function call; no-op when --debug is absent.
#define DBG_BREAK(name, world) do { if (dbg_break_match(name)) dbg_break_wait_impl((name), __FILE__, __LINE__, (world)); } while (0)
