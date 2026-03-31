// nudge.c -- physics world implementation

// -----------------------------------------------------------------------------
// Internal types: cold (metadata) and hot (solver iteration) splits.

typedef struct ShapeInternal
{
	ShapeType type;
	v3 local_pos;
	union {
		struct { float radius; } sphere;
		struct { float half_height; float radius; } capsule;
		struct { v3 half_extents; } box;
		struct { const Hull* hull; v3 scale; } hull;
	};
} ShapeInternal;

// Cold: persistent metadata, topology, rarely touched during solve.
typedef struct BodyCold
{
	float mass;
	CK_DYNA ShapeInternal* shapes;
} BodyCold;

// Hot: solver working set, iterated every step, packed for cache.
typedef struct BodyHot
{
	v3 position;
	quat rotation;
	v3 velocity;
	v3 angular_velocity;
	float inv_mass;
} BodyHot;

typedef struct WorldInternal
{
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
} WorldInternal;

#include "gjk.c"
#include "quickhull.c"
#include "collision.c"

// -----------------------------------------------------------------------------
// World.

World create_world(WorldParams params)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	w->gravity = params.gravity;
	return (World){ (uint64_t)w };
}

void destroy_world(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	for (int i = 0; i < asize(w->body_cold); i++) {
		afree(w->body_cold[i].shapes);
	}
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	CK_FREE(w);
}

void world_step(World world, float dt)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int count = asize(w->body_hot);

	// Integrate velocities
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
	}

	// Collision detection
	CK_DYNA InternalManifold* manifolds = NULL;
	broadphase_and_collide(w, &manifolds);

	// Resolve contacts (placeholder: simple impulse)
	for (int i = 0; i < asize(manifolds); i++) {
		InternalManifold* im = &manifolds[i];
		BodyHot* a = &w->body_hot[im->body_a];
		BodyHot* b = &w->body_hot[im->body_b];
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		if (inv_mass_sum == 0.0f) continue;

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			v3 rel_vel = sub(b->velocity, a->velocity);
			float vn = dot(rel_vel, ct->normal);
			if (vn > 0.0f) continue; // separating

			float j = -vn / inv_mass_sum;
			v3 impulse = scale(ct->normal, j);
			a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
			b->velocity = add(b->velocity, scale(impulse, b->inv_mass));

			// Positional correction (Baumgarte)
			float slop = 0.01f;
			float percent = 0.4f;
			float corr = fmaxf(ct->penetration - slop, 0.0f) * percent / inv_mass_sum;
			v3 correction = scale(ct->normal, corr);
			a->position = sub(a->position, scale(correction, a->inv_mass));
			b->position = add(b->position, scale(correction, b->inv_mass));
		}
	}
	afree(manifolds);

	// Integrate positions
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		h->position = add(h->position, scale(h->velocity, dt));
	}
}

// -----------------------------------------------------------------------------
// Body.

Body create_body(World world, BodyParams params)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	split_add(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);

	w->body_cold[idx] = (BodyCold){
		.mass = params.mass,
		.shapes = NULL,
	};
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
	};

	return split_handle(Body, w->body_gen, idx);
}

void destroy_body(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	afree(w->body_cold[idx].shapes);
	split_del(w->body_cold, w->body_hot, w->body_gen, w->body_free, idx);
}

void body_add_shape(World world, Body body, ShapeParams params)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));

	ShapeInternal s = {0};
	s.type = params.type;
	s.local_pos = params.local_pos;
	switch (params.type) {
	case SHAPE_SPHERE:  s.sphere.radius = params.sphere.radius; break;
	case SHAPE_CAPSULE: s.capsule.half_height = params.capsule.half_height;
	                    s.capsule.radius = params.capsule.radius; break;
	case SHAPE_BOX:     s.box.half_extents = params.box.half_extents; break;
	case SHAPE_HULL:    s.hull.hull = params.hull.hull;
	                    s.hull.scale = params.hull.scale; break;
	}
	apush(w->body_cold[idx].shapes, s);
}

v3 body_get_position(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].position;
}

quat body_get_rotation(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	return w->body_hot[idx].rotation;
}
