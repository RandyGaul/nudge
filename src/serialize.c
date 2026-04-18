// See LICENSE for licensing info.
// serialize.c -- binary versioned serialization.
//
// Port of Randy's "SV" (Save Version) design, originally inspired by Media
// Molecule. Each save file begins with a uint32 version; each serialized
// field is tagged with the version it was introduced, so older files can
// still be loaded by newer code (backwards-compatible). A removed field
// (via SV_REM) becomes a local variable on load so migrations can rewrite
// the old value into new state.
//
// Simple usage for a new type:
//
//     SV_SERIALIZABLE(MyType)
//     {
//         SV_ADD(SV_INITIAL, field_a);
//         SV_ADD(SV_ADD_B,   field_b);   // newer field
//     }
//
// Top-level save/load:
//
//     void save(World w, const char* path)
//     {
//         SV_SAVE_BEGIN(path);
//         SV_ADD_LOCAL(SV_INITIAL, some_data);
//         SV_SAVE_END();
//     }
//
// Add a new serializable type: declare it in SV_TYPES, bump the version
// enum, and write an SV_SERIALIZABLE routine.

// ---------------------------------------------------------------------------
// Version enum. Every change to any serialization routine bumps SV_LATEST.
// Add new entries just above SV_LATEST_PLUS_ONE, never reorder or remove.

enum
{
	SV_INITIAL,
	// Added deterministic internal world state: sleep + island linkage on
	// bodies, island linkage + warm lambda on joints, plus top-level frame /
	// version / island-scalars / warm_cache / prev_touching / joint_pairs
	// blocks so load+step reproduces the saved tick bit-for-bit.
	SV_DETERMINISTIC_WORLD,
	// Added sensors to snapshot save/load.
	SV_SENSORS_IN_SAVE,
	// --- insert new entries here ---
	SV_LATEST_PLUS_ONE
};

#define SV_LATEST (SV_LATEST_PLUS_ONE - 1)

// ---------------------------------------------------------------------------
// Types that can be serialized. Each needs an SV_SERIALIZABLE(T) routine.
// Types that never change (plain POD with no pointers) can go in
// SV_MEMCPY_SAFE_TYPES instead -- those get a single memcpy and arrays of
// them serialize in one call.

#define SV_TYPES(X) \
	X(BodyParams) \
	X(ShapeParams) \
	X(SpringParams) \
	X(BallSocketParams) \
	X(DistanceParams) \
	X(HingeParams) \
	X(FixedParams) \
	X(PrismaticParams) \
	X(AngularMotorParams) \
	X(TwistLimitParams) \
	X(ConeLimitParams) \
	X(SwingTwistParams) \
	X(WorldParams) \
	X(SavedBody) \
	X(SavedJoint) \
	X(SavedSensor) \
	X(CachedFeaturePair) \
	X(WarmContact) \
	X(WarmManifold) \
	X(Island)

#define SV_MEMCPY_SAFE_TYPES(X) \
	X(v3) \
	X(quat)

// Forward typedefs for all serializable user types.
#define SV_FORWARD(T) typedef struct T T;
SV_TYPES(SV_FORWARD)
#undef SV_FORWARD

// ---------------------------------------------------------------------------
// Context + open/close.

typedef struct SV_Context
{
	int version;        // file's schema version (== SV_LATEST when saving)
	int saving;
	int loading;
	FILE* file;
	int sync_counter;
} SV_Context;

// These open/close macros create a local `S` pointer so all SV_ADD calls
// in scope resolve to the same context.
#define SV_SAVE_BEGIN(path) SV_Context _sv_ctx = sv_make(path, 1), *S = &_sv_ctx
#define SV_SAVE_END()       sv_destroy(S)
#define SV_LOAD_BEGIN(path) SV_Context _sv_ctx = sv_make(path, 0), *S = &_sv_ctx
#define SV_LOAD_END()       sv_destroy(S)

static int sv_raw_write(SV_Context* S, const void* data, size_t n)
{
	return fwrite(data, 1, n, S->file) == n;
}

static int sv_raw_read(SV_Context* S, void* data, size_t n)
{
	return fread(data, 1, n, S->file) == n;
}

// ---------------------------------------------------------------------------
// Per-type declarations for the user-defined serializable types.

#define SV_DECL(T) static void sv_##T(SV_Context* S, T* o);
SV_TYPES(SV_DECL)
#undef SV_DECL

// ---------------------------------------------------------------------------
// Generic _Generic dispatch for SV_WRITE / SV_READ. Pointer-based so
// incomplete types can participate (pointers are always complete types).

#define SV_TYPE_FN(T) T*: sv_##T,
#define SV_MTYPE_FN(T) T*: sv_##T, T**: sv_##T##_array,

// Forward declarations for fundamental-type handlers.
static void sv_uint8 (SV_Context* S, uint8_t*  v);
static void sv_uint16(SV_Context* S, uint16_t* v);
static void sv_uint32(SV_Context* S, uint32_t* v);
static void sv_uint64(SV_Context* S, uint64_t* v);
static void sv_int8  (SV_Context* S, int8_t*   v);
static void sv_int16 (SV_Context* S, int16_t*  v);
static void sv_int32 (SV_Context* S, int32_t*  v);
static void sv_int64 (SV_Context* S, int64_t*  v);
static void sv_bool  (SV_Context* S, int*      v);
static void sv_float (SV_Context* S, float*    v);
static void sv_double(SV_Context* S, double*   v);
static void sv_cstr  (SV_Context* S, const char** v);
static void sv_v3    (SV_Context* S, v3*    v);
static void sv_quat  (SV_Context* S, quat*  v);
static void sv_v3_array  (SV_Context* S, v3**   a);
static void sv_quat_array(SV_Context* S, quat** a);

// Primary dispatch. Uses address-of so that _Generic sees pointer-to-type,
// which works for incomplete types (pointers are always complete).
#define SV_CALL(V) \
	_Generic(&(V), \
		SV_TYPES(SV_TYPE_FN) \
		SV_MEMCPY_SAFE_TYPES(SV_MTYPE_FN) \
		uint8_t*:      sv_uint8,  \
		uint16_t*:     sv_uint16, \
		uint32_t*:     sv_uint32, \
		uint64_t*:     sv_uint64, \
		int8_t*:       sv_int8,   \
		int16_t*:      sv_int16, \
		int32_t*:      sv_int32, \
		int64_t*:      sv_int64, \
		float*:        sv_float, \
		double*:       sv_double, \
		const char**:  sv_cstr, \
		char**:        sv_cstr \
	)(S, &(V))

// Internal: same call used for both save and load (the handler decides
// direction via S->saving / S->loading).
#define SV_WRITE(V) SV_CALL(V)
#define SV_READ(V)  SV_CALL(V)

// ---------------------------------------------------------------------------
// Field and local-variable macros.

// Serialize a struct member. On load, skip if the file's version predates
// when this field was introduced.
#define SV_ADD(VERSION, MEMBER) \
	do { \
		if (S->saving) { SV_CALL(o->MEMBER); } \
		else if (S->loading && S->version >= (VERSION)) { SV_CALL(o->MEMBER); } \
	} while (0)

// Serialize a local variable.
#define SV_ADD_LOCAL(VERSION, VAR) \
	do { \
		if (S->saving) { SV_CALL(VAR); } \
		else if (S->loading && S->version >= (VERSION)) { SV_CALL(VAR); } \
	} while (0)

// Handle a CK_DYNA member array: writes count then elements. On load
// resizes the live array and memsets new entries to zero (so SV_ADD
// overwrites from a clean slate).
#define SV_ADD_ARRAY(VERSION, ARRAY) \
	do { \
		if (S->saving || (S->loading && S->version >= (VERSION))) { \
			int _n = S->saving ? asize(o->ARRAY) : 0; \
			SV_ADD_LOCAL(VERSION, _n); \
			if (S->loading) { \
				afit(o->ARRAY, _n); asetlen(o->ARRAY, _n); \
				if (_n > 0) memset(o->ARRAY, 0, (size_t)_n * sizeof(o->ARRAY[0])); \
			} \
			for (int _i = 0; _i < _n; ++_i) { \
				if (S->saving) { SV_CALL(o->ARRAY[_i]); } \
				else { SV_CALL(o->ARRAY[_i]); } \
			} \
		} \
	} while (0)

// Local CK_DYNA array.
#define SV_ADD_LOCAL_ARRAY(VERSION, ARRAY) \
	do { \
		if (S->saving || (S->loading && S->version >= (VERSION))) { \
			int _n = S->saving ? asize(ARRAY) : 0; \
			SV_ADD_LOCAL(VERSION, _n); \
			if (S->loading) { \
				afit((ARRAY), _n); asetlen((ARRAY), _n); \
				if (_n > 0) memset((ARRAY), 0, (size_t)_n * sizeof((ARRAY)[0])); \
			} \
			for (int _i = 0; _i < _n; ++_i) { SV_CALL((ARRAY)[_i]); } \
		} \
	} while (0)

// Remove a field. Old saves' data for the field is read into a local
// variable named NAME with type T; caller can use it to migrate into new
// state. DEFAULT is the value NAME holds when the file predates the field.
//
//   SV_REM(SV_INITIAL, SV_REMOVED_FOO, int, old_x, 0);
//   o->new_x = old_x * 2;
#define SV_REM(VERSION_ADDED, VERSION_REMOVED, T, NAME, DEFAULT) \
	T NAME = (DEFAULT); \
	if (S->loading && S->version >= (VERSION_ADDED) && S->version < (VERSION_REMOVED)) { \
		SV_CALL(NAME); \
	}

// Optional sync-check: asserts on load if the per-routine counter diverges
// (indicating a missing SV_ADD/SV_REM somewhere). Place at end of a routine.
#define SV_SYNC() \
	do { \
		int _expected = S->sync_counter; \
		if (S->saving) { SV_CALL(_expected); } \
		else { int _got; SV_CALL(_got); assert(_got == _expected && "SV_SYNC mismatch"); } \
		S->sync_counter++; \
	} while (0)

// Define a serialization routine for user type T.
#define SV_SERIALIZABLE(T) static void sv_##T(SV_Context* S, T* o)

// ---------------------------------------------------------------------------
// Fundamental-type handlers.

static void sv_uint8 (SV_Context* S, uint8_t*  v) { if (S->saving) sv_raw_write(S, v, 1); else sv_raw_read(S, v, 1); }
static void sv_uint16(SV_Context* S, uint16_t* v) { if (S->saving) sv_raw_write(S, v, 2); else sv_raw_read(S, v, 2); }
static void sv_uint32(SV_Context* S, uint32_t* v) { if (S->saving) sv_raw_write(S, v, 4); else sv_raw_read(S, v, 4); }
static void sv_uint64(SV_Context* S, uint64_t* v) { if (S->saving) sv_raw_write(S, v, 8); else sv_raw_read(S, v, 8); }
static void sv_int8  (SV_Context* S, int8_t*   v) { if (S->saving) sv_raw_write(S, v, 1); else sv_raw_read(S, v, 1); }
static void sv_int16 (SV_Context* S, int16_t*  v) { if (S->saving) sv_raw_write(S, v, 2); else sv_raw_read(S, v, 2); }
static void sv_int32 (SV_Context* S, int32_t*  v) { if (S->saving) sv_raw_write(S, v, 4); else sv_raw_read(S, v, 4); }
static void sv_int64 (SV_Context* S, int64_t*  v) { if (S->saving) sv_raw_write(S, v, 8); else sv_raw_read(S, v, 8); }
static void sv_bool  (SV_Context* S, int*      v) { uint8_t b = S->saving ? (*v ? 1 : 0) : 0; sv_uint8(S, &b); if (S->loading) *v = (int)b; }
static void sv_float (SV_Context* S, float*    v) { if (S->saving) sv_raw_write(S, v, 4); else sv_raw_read(S, v, 4); }
static void sv_double(SV_Context* S, double*   v) { if (S->saving) sv_raw_write(S, v, 8); else sv_raw_read(S, v, 8); }
static void sv_cstr  (SV_Context* S, const char** v)
{
	uint32_t len;
	if (S->saving) {
		len = *v ? (uint32_t)strlen(*v) : 0;
		sv_uint32(S, &len);
		if (len > 0) sv_raw_write(S, *v, len);
	} else {
		sv_uint32(S, &len);
		char* buf = (char*)CK_ALLOC(len + 1);
		if (len > 0) sv_raw_read(S, buf, len);
		buf[len] = '\0';
		*v = sintern(buf);
		CK_FREE(buf);
	}
}

static void sv_v3   (SV_Context* S, v3*   v) { if (S->saving) sv_raw_write(S, v, sizeof(v3));   else sv_raw_read(S, v, sizeof(v3)); }
static void sv_quat (SV_Context* S, quat* v) { if (S->saving) sv_raw_write(S, v, sizeof(quat)); else sv_raw_read(S, v, sizeof(quat)); }

// Array handlers for memcpy-safe types. Single fwrite/fread of the entire array.
static void sv_v3_array(SV_Context* S, v3** a)
{
	if (S->saving) {
		int n = asize(*a);
		sv_int32(S, &n);
		if (n > 0) sv_raw_write(S, *a, (size_t)n * sizeof(v3));
	} else {
		int n = 0;
		sv_int32(S, &n);
		v3* p = *a;
		afit(p, n); asetlen(p, n);
		if (n > 0) sv_raw_read(S, p, (size_t)n * sizeof(v3));
		*a = p;
	}
}
static void sv_quat_array(SV_Context* S, quat** a)
{
	if (S->saving) {
		int n = asize(*a);
		sv_int32(S, &n);
		if (n > 0) sv_raw_write(S, *a, (size_t)n * sizeof(quat));
	} else {
		int n = 0;
		sv_int32(S, &n);
		quat* p = *a;
		afit(p, n); asetlen(p, n);
		if (n > 0) sv_raw_read(S, p, (size_t)n * sizeof(quat));
		*a = p;
	}
}

// ---------------------------------------------------------------------------
// Context lifecycle.

#define SV_MAGIC 0x5633444Eu  // "ND3V"

static SV_Context sv_make(const char* path, int saving)
{
	SV_Context ctx = { 0 };
	ctx.saving  = saving;
	ctx.loading = !saving;
	if (saving) {
		ctx.file = fopen(path, "wb");
		if (!ctx.file) return ctx;
		uint32_t magic = SV_MAGIC;
		fwrite(&magic, 4, 1, ctx.file);
		ctx.version = SV_LATEST;
		fwrite(&ctx.version, 4, 1, ctx.file);
	} else {
		ctx.file = fopen(path, "rb");
		if (!ctx.file) return ctx;
		uint32_t magic = 0;
		if (fread(&magic, 4, 1, ctx.file) != 1 || magic != SV_MAGIC) {
			fclose(ctx.file);
			ctx.file = NULL;
			return ctx;
		}
		if (fread(&ctx.version, 4, 1, ctx.file) != 1) {
			fclose(ctx.file);
			ctx.file = NULL;
			return ctx;
		}
		// Newer files can't be loaded by older code (forward incompat).
		if (ctx.version > SV_LATEST) {
			fclose(ctx.file);
			ctx.file = NULL;
			return ctx;
		}
	}
	return ctx;
}

static void sv_destroy(SV_Context* S)
{
	if (S->file) { fclose(S->file); S->file = NULL; }
}

// ---------------------------------------------------------------------------
// User-type serializers for public API structs. These describe the canonical
// scene-save format. Internal world state (warm cache, islands, BVH, LDL,
// rewind ring) is NOT saved -- it's reconstructed on load via fresh
// create_body / create_joint calls. Simulation after a load may differ
// from a pure-rewind-restore by a few frames of warm-cache rebuild; the
// scene topology and body poses are preserved exactly.

SV_SERIALIZABLE(SpringParams)
{
	SV_ADD(SV_INITIAL, frequency);
	SV_ADD(SV_INITIAL, damping_ratio);
}

SV_SERIALIZABLE(BodyParams)
{
	SV_ADD(SV_INITIAL, position);
	SV_ADD(SV_INITIAL, rotation);
	SV_ADD(SV_INITIAL, mass);
	SV_ADD(SV_INITIAL, friction);
	SV_ADD(SV_INITIAL, restitution);
	SV_ADD(SV_INITIAL, linear_damping);
	SV_ADD(SV_INITIAL, angular_damping);
}

SV_SERIALIZABLE(ShapeParams)
{
	int type = (int)o->type;
	SV_ADD_LOCAL(SV_INITIAL, type);
	if (S->loading) o->type = (ShapeType)type;
	SV_ADD(SV_INITIAL, local_pos);
	SV_ADD(SV_INITIAL, local_rot);
	switch (o->type) {
	case SHAPE_SPHERE:
		SV_ADD(SV_INITIAL, sphere.radius);
		break;
	case SHAPE_CAPSULE:
		SV_ADD(SV_INITIAL, capsule.half_height);
		SV_ADD(SV_INITIAL, capsule.radius);
		break;
	case SHAPE_BOX:
		SV_ADD(SV_INITIAL, box.half_extents);
		break;
	case SHAPE_CYLINDER:
		SV_ADD(SV_INITIAL, cylinder.half_height);
		SV_ADD(SV_INITIAL, cylinder.radius);
		break;
	case SHAPE_HULL:
	case SHAPE_MESH:
		// Hull and mesh shapes depend on caller-owned Hull* / TriMesh*
		// pointers; snapshots don't capture that geometry. Save/load of a
		// scene that contains them asserts.
		assert(0 && "SV: SHAPE_HULL and SHAPE_MESH are not supported in snapshots yet");
		break;
	}
}

SV_SERIALIZABLE(WorldParams)
{
	SV_ADD(SV_INITIAL, gravity);
	int bp = (int)o->broadphase;
	int nb = (int)o->narrowphase_backend;
	int st = (int)o->solver_type;
	SV_ADD_LOCAL(SV_INITIAL, bp);
	SV_ADD_LOCAL(SV_INITIAL, nb);
	SV_ADD_LOCAL(SV_INITIAL, st);
	if (S->loading) {
		o->broadphase = (BroadphaseType)bp;
		o->narrowphase_backend = (NarrowphaseBackend)nb;
		o->solver_type = (SolverType)st;
	}
	SV_ADD(SV_INITIAL, velocity_iters);
	SV_ADD(SV_INITIAL, position_iters);
	SV_ADD(SV_INITIAL, contact_hertz);
	SV_ADD(SV_INITIAL, contact_damping_ratio);
	SV_ADD(SV_INITIAL, max_push_velocity);
	SV_ADD(SV_INITIAL, sub_steps);
}

// SavedBody: one body's full serializable state. Velocity is captured so
// a load right after save restores dynamics as well as pose.
struct SavedBody
{
	BodyParams params;
	v3 velocity;
	v3 angular_velocity;
	uint32_t collision_group;
	uint32_t collision_mask;
	uint32_t compound_id;
	CK_DYNA ShapeParams* shapes;
	// Added in SV_DETERMINISTIC_WORLD so load+step matches save+step.
	float sleep_time;
	int sleep_allowed;
	int island_id;
	int island_prev;
	int island_next;
};

SV_SERIALIZABLE(SavedBody)
{
	SV_ADD(SV_INITIAL, params);
	SV_ADD(SV_INITIAL, velocity);
	SV_ADD(SV_INITIAL, angular_velocity);
	SV_ADD(SV_INITIAL, collision_group);
	SV_ADD(SV_INITIAL, collision_mask);
	SV_ADD(SV_INITIAL, compound_id);
	SV_ADD_ARRAY(SV_INITIAL, shapes);
	SV_ADD(SV_DETERMINISTIC_WORLD, sleep_time);
	SV_ADD(SV_DETERMINISTIC_WORLD, sleep_allowed);
	SV_ADD(SV_DETERMINISTIC_WORLD, island_id);
	SV_ADD(SV_DETERMINISTIC_WORLD, island_prev);
	SV_ADD(SV_DETERMINISTIC_WORLD, island_next);
}

// Joint params for each joint type. Each references bodies by saved-index
// (0..N-1) rather than by handle; the caller remaps on save and load.
SV_SERIALIZABLE(BallSocketParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(DistanceParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, rest_length);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(HingeParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(FixedParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(PrismaticParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(AngularMotorParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, target_speed);
	SV_ADD(SV_INITIAL, max_impulse);
}

SV_SERIALIZABLE(TwistLimitParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, limit_min);
	SV_ADD(SV_INITIAL, limit_max);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(ConeLimitParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, half_angle);
	SV_ADD(SV_INITIAL, spring);
}

SV_SERIALIZABLE(SwingTwistParams)
{
	SV_ADD(SV_INITIAL, body_a.id);
	SV_ADD(SV_INITIAL, body_b.id);
	SV_ADD(SV_INITIAL, local_offset_a);
	SV_ADD(SV_INITIAL, local_offset_b);
	SV_ADD(SV_INITIAL, local_axis_a);
	SV_ADD(SV_INITIAL, local_axis_b);
	SV_ADD(SV_INITIAL, cone_half_angle);
	SV_ADD(SV_INITIAL, twist_min);
	SV_ADD(SV_INITIAL, twist_max);
	SV_ADD(SV_INITIAL, spring);
}

// SavedJoint: typed union of every joint-params variant. Stored separately
// per joint rather than as a single flat array because the params are
// heterogeneous.
struct SavedJoint
{
	int type;  // JointType enum value
	union {
		BallSocketParams   ball_socket;
		DistanceParams     distance;
		HingeParams        hinge;
		FixedParams        fixed;
		PrismaticParams    prismatic;
		AngularMotorParams angular_motor;
		TwistLimitParams   twist_limit;
		ConeLimitParams    cone_limit;
		SwingTwistParams   swing_twist;
	};
	// Added in SV_DETERMINISTIC_WORLD.
	int island_id;
	int island_prev;
	int island_next;
	float warm_lambda[JOINT_MAX_DOF];  // JointHot[j].warm_lambda
};

SV_SERIALIZABLE(SavedJoint)
{
	SV_ADD(SV_INITIAL, type);
	switch (o->type) {
	case JOINT_BALL_SOCKET:   SV_ADD(SV_INITIAL, ball_socket);   break;
	case JOINT_DISTANCE:      SV_ADD(SV_INITIAL, distance);      break;
	case JOINT_HINGE:         SV_ADD(SV_INITIAL, hinge);         break;
	case JOINT_FIXED:         SV_ADD(SV_INITIAL, fixed);         break;
	case JOINT_PRISMATIC:     SV_ADD(SV_INITIAL, prismatic);     break;
	case JOINT_ANGULAR_MOTOR: SV_ADD(SV_INITIAL, angular_motor); break;
	case JOINT_TWIST_LIMIT:   SV_ADD(SV_INITIAL, twist_limit);   break;
	case JOINT_CONE_LIMIT:    SV_ADD(SV_INITIAL, cone_limit);    break;
	case JOINT_SWING_TWIST:   SV_ADD(SV_INITIAL, swing_twist);   break;
	}
	SV_ADD(SV_DETERMINISTIC_WORLD, island_id);
	SV_ADD(SV_DETERMINISTIC_WORLD, island_prev);
	SV_ADD(SV_DETERMINISTIC_WORLD, island_next);
	for (int i = 0; i < JOINT_MAX_DOF; i++)
		SV_ADD(SV_DETERMINISTIC_WORLD, warm_lambda[i]);
}

// Internal-state serializers -- used by snapshot.c to restore impulse caches,
// island linkage, and narrowphase hints for bit-identical determinism.

SV_SERIALIZABLE(CachedFeaturePair)
{
	SV_ADD(SV_DETERMINISTIC_WORLD, type);
	SV_ADD(SV_DETERMINISTIC_WORLD, ref_body);
	SV_ADD(SV_DETERMINISTIC_WORLD, face_a);
	SV_ADD(SV_DETERMINISTIC_WORLD, face_b);
	SV_ADD(SV_DETERMINISTIC_WORLD, edge_a);
	SV_ADD(SV_DETERMINISTIC_WORLD, edge_b);
}

SV_SERIALIZABLE(WarmContact)
{
	SV_ADD(SV_DETERMINISTIC_WORLD, r_a);
	SV_ADD(SV_DETERMINISTIC_WORLD, lambda_n);
	SV_ADD(SV_DETERMINISTIC_WORLD, feature_id);
}

SV_SERIALIZABLE(WarmManifold)
{
	SV_ADD(SV_DETERMINISTIC_WORLD, count);
	SV_ADD(SV_DETERMINISTIC_WORLD, stale);
	SV_ADD(SV_DETERMINISTIC_WORLD, body_a);
	SV_ADD(SV_DETERMINISTIC_WORLD, body_b);
	SV_ADD(SV_DETERMINISTIC_WORLD, manifold_lambda_t1);
	SV_ADD(SV_DETERMINISTIC_WORLD, manifold_lambda_t2);
	SV_ADD(SV_DETERMINISTIC_WORLD, manifold_lambda_twist);
	SV_ADD(SV_DETERMINISTIC_WORLD, sat_axis);
	for (int i = 0; i < MAX_CONTACTS; i++)
		SV_ADD(SV_DETERMINISTIC_WORLD, contacts[i]);
	SV_ADD(SV_DETERMINISTIC_WORLD, cached_pair);
}

SV_SERIALIZABLE(Island)
{
	SV_ADD(SV_DETERMINISTIC_WORLD, head_body);
	SV_ADD(SV_DETERMINISTIC_WORLD, tail_body);
	SV_ADD(SV_DETERMINISTIC_WORLD, body_count);
	SV_ADD(SV_DETERMINISTIC_WORLD, head_joint);
	SV_ADD(SV_DETERMINISTIC_WORLD, tail_joint);
	SV_ADD(SV_DETERMINISTIC_WORLD, joint_count);
	SV_ADD(SV_DETERMINISTIC_WORLD, constraint_remove_count);
	SV_ADD(SV_DETERMINISTIC_WORLD, awake);
	// LDL_Cache ldl -- recomputed next ldl_factor; don't save.
}

// SavedSensor: one world-owned sensor's full state.
struct SavedSensor
{
	v3 position;
	quat rotation;
	uint32_t collision_group;
	uint32_t collision_mask;
	CK_DYNA ShapeParams* shapes;
};

SV_SERIALIZABLE(SavedSensor)
{
	SV_ADD(SV_SENSORS_IN_SAVE, position);
	SV_ADD(SV_SENSORS_IN_SAVE, rotation);
	SV_ADD(SV_SENSORS_IN_SAVE, collision_group);
	SV_ADD(SV_SENSORS_IN_SAVE, collision_mask);
	SV_ADD_ARRAY(SV_SENSORS_IN_SAVE, shapes);
}
