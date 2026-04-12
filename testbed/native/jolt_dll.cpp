// jolt_dll.cpp -- thin C wrapper over JoltPhysics C++ API for P/Invoke.

#include <Jolt/Jolt.h>
#include <Jolt/RegisterTypes.h>
#include <Jolt/Core/Factory.h>
#include <Jolt/Core/TempAllocator.h>
#include <Jolt/Core/JobSystemThreadPool.h>
#include <Jolt/Physics/PhysicsSystem.h>
#include <Jolt/Physics/PhysicsSettings.h>
#include <Jolt/Physics/Collision/Shape/BoxShape.h>
#include <Jolt/Physics/Collision/Shape/SphereShape.h>
#include <Jolt/Physics/Collision/Shape/CapsuleShape.h>
#include <Jolt/Physics/Body/BodyCreationSettings.h>
#include <Jolt/Physics/Constraints/DistanceConstraint.h>
#include <chrono>
#include <vector>
#include <cstdio>
#include <cstdarg>

using namespace JPH;

// ---------------------------------------------------------------------------
// Layers (minimal two-layer setup).

namespace Layers {
	static constexpr ObjectLayer NON_MOVING = 0;
	static constexpr ObjectLayer MOVING = 1;
	static constexpr uint32_t NUM_LAYERS = 2;
}

namespace BroadPhaseLayers {
	static constexpr BroadPhaseLayer NON_MOVING(0);
	static constexpr BroadPhaseLayer MOVING(1);
	static constexpr uint32_t NUM_LAYERS = 2;
}

class BPLayerInterfaceImpl final : public BroadPhaseLayerInterface {
public:
	BPLayerInterfaceImpl() {
		mObjectToBroadPhase[Layers::NON_MOVING] = BroadPhaseLayers::NON_MOVING;
		mObjectToBroadPhase[Layers::MOVING] = BroadPhaseLayers::MOVING;
	}
	uint GetNumBroadPhaseLayers() const override { return BroadPhaseLayers::NUM_LAYERS; }
	BroadPhaseLayer GetBroadPhaseLayer(ObjectLayer inLayer) const override { return mObjectToBroadPhase[inLayer]; }
	const char* GetBroadPhaseLayerName(BroadPhaseLayer inLayer) const override { return inLayer == BroadPhaseLayers::NON_MOVING ? "NON_MOVING" : "MOVING"; }
private:
	BroadPhaseLayer mObjectToBroadPhase[Layers::NUM_LAYERS];
};

class ObjectVsBroadPhaseLayerFilterImpl : public ObjectVsBroadPhaseLayerFilter {
public:
	bool ShouldCollide(ObjectLayer inLayer1, BroadPhaseLayer inLayer2) const override {
		if (inLayer1 == Layers::NON_MOVING) return inLayer2 == BroadPhaseLayers::MOVING;
		return true;
	}
};

class ObjectLayerPairFilterImpl : public ObjectLayerPairFilter {
public:
	bool ShouldCollide(ObjectLayer inLayer1, ObjectLayer inLayer2) const override {
		if (inLayer1 == Layers::NON_MOVING && inLayer2 == Layers::NON_MOVING) return false;
		return true;
	}
};

// ---------------------------------------------------------------------------
// World handle.

struct JoltWorld {
	PhysicsSystem system;
	BPLayerInterfaceImpl bp_layer_interface;
	ObjectVsBroadPhaseLayerFilterImpl obj_vs_bp_filter;
	ObjectLayerPairFilterImpl obj_pair_filter;
	TempAllocatorImpl* temp_allocator;
	JobSystemThreadPool* job_system;
	std::vector<BodyID> bodies;
	std::vector<Body*> body_ptrs;
	double last_step_time;
};

static bool s_initialized = false;

static void trace_impl(const char* fmt, ...) {
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "\n");
	fflush(stderr);
	va_end(args);
}

#ifdef JPH_ENABLE_ASSERTS
static bool assert_impl(const char* expr, const char* msg, const char* file, uint line) {
	fprintf(stderr, "JOLT ASSERT: %s:%u: %s (%s)\n", file, line, expr, msg ? msg : "");
	fflush(stderr);
	return true;
}
#endif

static void ensure_init() {
	if (!s_initialized) {
		RegisterDefaultAllocator();
		Trace = trace_impl;
#ifdef JPH_ENABLE_ASSERTS
		AssertFailed = assert_impl;
#endif
		Factory::sInstance = new Factory();
		RegisterTypes();
		s_initialized = true;
	}
}

// ---------------------------------------------------------------------------
// Exported C API.

#ifdef _WIN32
#define EXPORT extern "C" __declspec(dllexport)
#else
#define EXPORT extern "C" __attribute__((visibility("default")))
#endif

// Shape type constants (match C# ShapeType enum)
enum { SHAPE_BOX = 0, SHAPE_SPHERE = 1, SHAPE_CAPSULE = 2 };

EXPORT void* jolt_create_world(float gx, float gy, float gz, int max_bodies) {
	ensure_init();
	auto* w = new JoltWorld();
	if (max_bodies < 1024) max_bodies = 1024;
	w->system.Init((uint)max_bodies, 0, (uint)max_bodies, (uint)max_bodies,
		w->bp_layer_interface, w->obj_vs_bp_filter, w->obj_pair_filter);
	w->system.SetGravity(Vec3(gx, gy, gz));
	w->temp_allocator = new TempAllocatorImpl(128 * 1024 * 1024);
	w->job_system = new JobSystemThreadPool(cMaxPhysicsJobs, cMaxPhysicsBarriers, 1);
	w->last_step_time = 0;
	return w;
}

EXPORT void jolt_destroy_world(void* world) {
	auto* w = (JoltWorld*)world;
	auto& bi = w->system.GetBodyInterface();
	for (auto& id : w->bodies) {
		bi.RemoveBody(id);
		bi.DestroyBody(id);
	}
	delete w->job_system;
	delete w->temp_allocator;
	delete w;
}

static int jolt_create_body_impl(JoltWorld* w, int shape_type, float s0, float s1, float s2, float px, float py, float pz, Quat rot, float mass, float friction, float restitution) {
	auto& bi = w->system.GetBodyInterface();

	const Shape* shape;
	switch (shape_type) {
		case SHAPE_BOX:     shape = new BoxShape(Vec3(s0, s1, s2)); break;
		case SHAPE_SPHERE:  shape = new SphereShape(s0); break;
		case SHAPE_CAPSULE: shape = new CapsuleShape(s0, s1); break;
		default:            shape = new BoxShape(Vec3(0.5f, 0.5f, 0.5f)); break;
	}

	EMotionType motion = (mass <= 0) ? EMotionType::Static : EMotionType::Dynamic;
	ObjectLayer layer = (mass <= 0) ? Layers::NON_MOVING : Layers::MOVING;
	EActivation activation = (mass <= 0) ? EActivation::DontActivate : EActivation::Activate;

	BodyCreationSettings settings(shape, RVec3(px, py, pz), rot, motion, layer);
	if (mass > 0) {
		settings.mOverrideMassProperties = EOverrideMassProperties::CalculateInertia;
		settings.mMassPropertiesOverride.mMass = mass;
	}
	settings.mFriction = (friction > 0) ? friction : 0.5f;
	settings.mRestitution = restitution;

	Body* body = bi.CreateBody(settings);
	bi.AddBody(body->GetID(), activation);
	int index = (int)w->bodies.size();
	w->bodies.push_back(body->GetID());
	w->body_ptrs.push_back(body);
	return index;
}

// Combined create body + shape. Returns body index.
EXPORT int jolt_create_body(void* world, int shape_type, float s0, float s1, float s2, float px, float py, float pz, float mass, float friction, float restitution) {
	return jolt_create_body_impl((JoltWorld*)world, shape_type, s0, s1, s2, px, py, pz, Quat::sIdentity(), mass, friction, restitution);
}

EXPORT int jolt_create_body_rotated(void* world, int shape_type, float s0, float s1, float s2, float px, float py, float pz, float qx, float qy, float qz, float qw, float mass, float friction, float restitution) {
	return jolt_create_body_impl((JoltWorld*)world, shape_type, s0, s1, s2, px, py, pz, Quat(qx, qy, qz, qw), mass, friction, restitution);
}

EXPORT void jolt_step(void* world, float dt) {
	auto* w = (JoltWorld*)world;
	auto start = std::chrono::high_resolution_clock::now();
	w->system.Update(dt, 1, w->temp_allocator, w->job_system);
	auto end = std::chrono::high_resolution_clock::now();
	w->last_step_time = std::chrono::duration<double>(end - start).count();
}

EXPORT void jolt_get_position(void* world, int body_index, float* out) {
	auto* w = (JoltWorld*)world;
	auto& bi = w->system.GetBodyInterface();
	RVec3 pos = bi.GetCenterOfMassPosition(w->bodies[body_index]);
	out[0] = (float)pos.GetX();
	out[1] = (float)pos.GetY();
	out[2] = (float)pos.GetZ();
}

EXPORT void jolt_get_rotation(void* world, int body_index, float* out) {
	auto* w = (JoltWorld*)world;
	auto& bi = w->system.GetBodyInterface();
	Quat rot = bi.GetRotation(w->bodies[body_index]);
	out[0] = rot.GetX();
	out[1] = rot.GetY();
	out[2] = rot.GetZ();
	out[3] = rot.GetW();
}

EXPORT int jolt_is_body_active(void* world, int body_index) {
	auto* w = (JoltWorld*)world;
	return w->system.GetBodyInterface().IsActive(w->bodies[body_index]) ? 1 : 0;
}

EXPORT int jolt_get_active_count(void* world) {
	auto* w = (JoltWorld*)world;
	auto& bi = w->system.GetBodyInterface();
	int count = 0;
	for (auto& id : w->bodies) {
		if (bi.IsActive(id)) count++;
	}
	return count;
}

EXPORT double jolt_get_step_time(void* world) {
	return ((JoltWorld*)world)->last_step_time;
}

EXPORT void jolt_set_velocity(void* world, int body_index, float vx, float vy, float vz) {
	auto* w = (JoltWorld*)world;
	w->system.GetBodyInterface().SetLinearVelocity(w->bodies[body_index], Vec3(vx, vy, vz));
}

EXPORT void jolt_optimize_broadphase(void* world) {
	((JoltWorld*)world)->system.OptimizeBroadPhase();
}

EXPORT void jolt_create_distance_joint(void* world, int body_a, int body_b, float ax, float ay, float az, float bx, float by, float bz, float rest_length) {
	auto* w = (JoltWorld*)world;

	DistanceConstraintSettings settings;
	settings.mSpace = EConstraintSpace::LocalToBodyCOM;
	settings.mPoint1 = Vec3(ax, ay, az);
	settings.mPoint2 = Vec3(bx, by, bz);
	settings.mMinDistance = rest_length;
	settings.mMaxDistance = rest_length;

	w->system.AddConstraint(settings.Create(*w->body_ptrs[body_a], *w->body_ptrs[body_b]));
}
