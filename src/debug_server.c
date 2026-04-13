// debug_server.c -- Debug driver + type reflection for nudge physics engine.
// Tiny TCP server for remote control + reflection data export.
// The viewer reads memory via ReadProcessMemory; this file provides the layout.
//
// Included from main.c (unity build).
// Winsock headers included from main.c before SDL (ordering requirement).

#define DBG_PORT 9999
#define DBG_MAX_CLIENTS 4
#define DBG_RECV_BUF 4096

// Remote drag state (mirrors the mouse constraint but driven via TCP).
static Body g_drag_body;
static Body g_drag_anchor;
static Joint g_drag_joint;
static v3 g_drag_local_hit;

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

static const R_EnumEntry renum_ShapeType_entries[] = { {"sphere",0}, {"capsule",1}, {"box",2}, {"hull",3}, {"cylinder",4} };
static const R_EnumEntry renum_SolverType_entries[] = { {"soft_step",0}, {"si_soft",1}, {"si",2} };
static const R_EnumEntry renum_BroadphaseType_entries[] = { {"n2",0}, {"bvh",1} };
static const R_EnumEntry renum_JointType_entries[] = { {"ball_socket",0}, {"distance",1}, {"hinge",2}, {"fixed",3}, {"prismatic",4} };

static const R_EnumTable g_enum_tables[] = {
	{ "ShapeType", 5, renum_ShapeType_entries },
	{ "SolverType", 3, renum_SolverType_entries },
	{ "BroadphaseType", 2, renum_BroadphaseType_entries },
	{ "JointType", 5, renum_JointType_entries },
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

REFLECT(WarmManifold,
	RF_INT(WarmManifold, count),
	RF_INT(WarmManifold, stale),
	RF_FLOAT(WarmManifold, manifold_lambda_t1),
	RF_FLOAT(WarmManifold, manifold_lambda_t2),
	RF_FLOAT(WarmManifold, manifold_lambda_twist),
	RF_INT(WarmManifold, sat_axis),
);

REFLECT(JointInternal,
	RF_ENUM(JointInternal, type, JointType),
	RF_INT(JointInternal, body_a),
	RF_INT(JointInternal, body_b),
	RF_INT(JointInternal, island_id),
);

REFLECT(Island,
	RF_INT(Island, head_body),
	RF_INT(Island, tail_body),
	RF_INT(Island, body_count),
	RF_INT(Island, head_joint),
	RF_INT(Island, tail_joint),
	RF_INT(Island, joint_count),
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
	RF_ARRAY(WorldInternal, joint_gen, uint32_t),
	RF_ENUM(WorldInternal, broadphase_type, BroadphaseType),
	RF_PTR(WorldInternal, bvh_static, BVH_Tree),
	RF_PTR(WorldInternal, bvh_dynamic, BVH_Tree),
	RF_PTR(WorldInternal, bvh_sleeping, BVH_Tree),
	RF_ARRAY(WorldInternal, islands, Island),
	RF_ARRAY(WorldInternal, island_gen, uint32_t),
	RF_INT(WorldInternal, sleep_enabled),
	RF_INT(WorldInternal, sat_hint_enabled),
	RF_INT(WorldInternal, sat_hillclimb_enabled),
	RF_INT(WorldInternal, incremental_np_enabled),
	RF_INT(WorldInternal, warm_start_enabled),
	RF_INT(WorldInternal, ldl_enabled),
	RF_INT(WorldInternal, ldl_topo_version),
	RF_ENUM(WorldInternal, solver_type, SolverType),
	RF_INT(WorldInternal, velocity_iters),
	RF_INT(WorldInternal, position_iters),
	RF_FLOAT(WorldInternal, contact_hertz),
	RF_FLOAT(WorldInternal, contact_damping_ratio),
	RF_FLOAT(WorldInternal, max_push_velocity),
	RF_INT(WorldInternal, sub_steps),
	RF_STRUCT(WorldInternal, perf, PerfTimers),
	RF_ARRAY(WorldInternal, dbg_solver_manifolds, SolverManifold),
	RF_ARRAY(WorldInternal, dbg_solver_contacts, SolverContact),
	RF_ARRAY(WorldInternal, dbg_solver_joints, SolverJoint),
);

// Forward declarations for send helpers (defined below).
typedef struct DbgClient DbgClient;
static int dbg_send_str(DbgClient *c, const char *s);

// Global type registry (all types the viewer can navigate).
static const R_TypeDesc *g_type_registry[] = {
	&rtype_WorldInternal,
	&rtype_BodyCold, &rtype_BodyHot, &rtype_BodyState,
	&rtype_ShapeInternal, &rtype_Contact,
	&rtype_SolverContact, &rtype_SolverManifold, &rtype_SolverJoint,
	&rtype_WarmManifold, &rtype_JointInternal, &rtype_Island,
	&rtype_BVH_Child, &rtype_BVHNode, &rtype_BVHLeaf, &rtype_BVH_Tree,
	&rtype_PerfTimers, &rtype_PGSTimers,
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
			"OK commands:\n"
			"  help                          show this help\n"
			"  info                          PID, world pointer, scene, frame\n"
			"  types                         send type reflection data\n"
			"  pause / step [n]              pause, step n frames\n"
			"  scene <name> / scenes         load scene, list scenes\n"
			"  restart                       restart current scene\n"
			"  push <idx> <f> [r]            impulse at point\n"
			"  highlight <idx> <color|r g b> tint a body\n"
			"  unhighlight <idx|all>         remove tint\n"
			"  label <text>                  show overlay text\n"
			"  slow <factor>                 time scale (0.1-10)\n"
		);
	} else if (strcmp(cmd, "info") == 0) {
		WorldInternal *w = (WorldInternal *)g_world.id;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK pid=%lu world=0x%llX scene=\"%s\" frame=%d bodies=%d paused=%s\n",
			(unsigned long)GetCurrentProcessId(),
			(unsigned long long)g_world.id,
			g_scenes[g_scene_index].name,
			w->frame,
			asize(g_draw_list),
			g_paused ? "true" : "false");
		dbg_reply(c, s);
	} else if (strcmp(cmd, "types") == 0) {
		dbg_send_types(c);
	} else if (strcmp(cmd, "pause") == 0) {
		g_paused = !g_paused;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK paused=%s\n", g_paused ? "true" : "false");
		dbg_reply(c, s);
	} else if (strcmp(cmd, "step") == 0) {
		int n = (*rest) ? atoi(rest) : 1;
		if (n < 1) n = 1;
		for (int i = 0; i < n; i++) world_step(g_world, 1.0f / 60.0f);
		WorldInternal *w = (WorldInternal *)g_world.id;
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
			WorldInternal *w = (WorldInternal *)g_world.id;
			CK_SDYNA char *s = NULL;
			sfmt(s, "OK scene=\"%s\" world=0x%llX frame=%d\n", g_scenes[found].name, (unsigned long long)g_world.id, w->frame);
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
		WorldInternal *w = (WorldInternal *)g_world.id;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK world=0x%llX frame=%d\n", (unsigned long long)g_world.id, w->frame);
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
		WorldInternal *w = (WorldInternal *)g_world.id;
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
		WorldInternal *w = (WorldInternal *)g_world.id;
		if (idx < 0 || idx >= asize(w->body_hot) || !(w->body_gen[idx] & 1)) {
			dbg_send_str(c, "ERR invalid body index\n");
			return;
		}
		if (g_drag_body.id) {
			destroy_joint(g_world, g_drag_joint);
			destroy_body(g_world, g_drag_anchor);
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
		g_drag_anchor = create_body(g_world, (BodyParams){ .position = V3(tx, ty, tz), .rotation = quat_identity(), .mass = 0 });
		// Softer spring than mouse constraint (2Hz vs 5Hz) because TCP
		// commands move the target in discrete jumps, not smooth mouse motion.
		g_drag_joint = create_ball_socket(g_world, (BallSocketParams){
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
		WorldInternal *w = (WorldInternal *)g_world.id;
		int anchor_idx = handle_index(g_drag_anchor);
		body_pos(w, anchor_idx) = V3(tx, ty, tz);
		int joint_idx = handle_index(g_drag_joint);
		int isl = w->joints[joint_idx].island_id;
		if (isl >= 0 && !w->islands[isl].awake) island_wake(w, isl);
		dbg_send_str(c, "OK\n");
	} else if (strcmp(cmd, "release") == 0) {
		if (!g_drag_body.id) { dbg_send_str(c, "ERR not dragging\n"); return; }
		destroy_joint(g_world, g_drag_joint);
		destroy_body(g_world, g_drag_anchor);
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
		WorldInternal *pw = (WorldInternal *)g_world.id;
		for (int i = pw->frame; i < rec_start; i++) world_step(g_world, 1.0f / 60.0f);
		g_recording = 2;
		g_rec_body_draw_idx = -1;
		g_playback_frame = 0;
		g_paused = false;
		CK_SDYNA char *s = NULL;
		sfmt(s, "OK playing %d frames (scene=%d start=%d) world=0x%llX\n", rec_n, rec_scene, rec_start, (unsigned long long)g_world.id);
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
	addr.sin_port = htons(DBG_PORT);

	if (bind(g_dbg.listen_sock, (struct sockaddr *)&addr, sizeof(addr)) == SOCKET_ERROR) {
		closesocket(g_dbg.listen_sock); g_dbg.listen_sock = INVALID_SOCKET; return;
	}
	if (listen(g_dbg.listen_sock, 4) == SOCKET_ERROR) {
		closesocket(g_dbg.listen_sock); g_dbg.listen_sock = INVALID_SOCKET; return;
	}

	g_dbg.initialized = true;
	printf("[dbg] listening on port %d (PID %lu, world=0x%llX)\n",
		DBG_PORT, (unsigned long)GetCurrentProcessId(), (unsigned long long)g_world.id);
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
			(unsigned long)GetCurrentProcessId(), (unsigned long long)g_world.id);
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
