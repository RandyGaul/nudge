// See LICENSE for licensing info.
// sensor.c -- read-only world-query volumes.
//
// Sensors are compounds of convex shapes with a transform. They never touch
// the solver, never enter the broadphase, and never generate contacts.
// sensor_query walks the world's BVHs, filters by collision group/mask, and
// runs shape-vs-shape narrowphase boolean tests against each candidate body.
//
// All state lives in a caller-owned heap allocation; the Sensor handle is a
// raw pointer. Lifetime errors (query against a destroyed sensor) crash, per
// the engine's non-defensive policy.

typedef struct SensorInternal
{
	v3 position;
	quat rotation;
	uint32_t collision_group;
	uint32_t collision_mask;
	CK_DYNA ShapeInternal* shapes;
} SensorInternal;

Sensor create_sensor(SensorParams params)
{
	assert(is_valid(params.position) && "create_sensor: position is NaN/inf");

	SensorInternal* s = (SensorInternal*)calloc(1, sizeof(SensorInternal));
	s->position = params.position;
	float r_m2 = params.rotation.x*params.rotation.x + params.rotation.y*params.rotation.y + params.rotation.z*params.rotation.z + params.rotation.w*params.rotation.w;
	s->rotation = (r_m2 < 0.5f) ? quat_identity() : params.rotation;
	s->collision_group = params.collision_group ? params.collision_group : 0xFFFFFFFFu;
	s->collision_mask  = params.collision_mask  ? params.collision_mask  : 0xFFFFFFFFu;
	return (Sensor){ (uint64_t)(uintptr_t)s };
}

void destroy_sensor(Sensor sensor)
{
	SensorInternal* s = (SensorInternal*)(uintptr_t)sensor.id;
	afree(s->shapes);
	free(s);
}

void sensor_add_shape(Sensor sensor, ShapeParams params)
{
	assert(is_valid(params.local_pos) && "sensor_add_shape: local_pos is NaN/inf");
	assert(params.type != SHAPE_MESH && "sensor_add_shape: SHAPE_MESH not supported on sensors");

	SensorInternal* s = (SensorInternal*)(uintptr_t)sensor.id;
	ShapeInternal sh = {0};
	sh.type = params.type;
	sh.local_pos = params.local_pos;
	float lr_m2 = params.local_rot.x*params.local_rot.x + params.local_rot.y*params.local_rot.y + params.local_rot.z*params.local_rot.z + params.local_rot.w*params.local_rot.w;
	sh.local_rot = (lr_m2 < 0.5f) ? quat_identity() : params.local_rot;
	switch (params.type) {
	case SHAPE_SPHERE:   sh.sphere.radius = params.sphere.radius; break;
	case SHAPE_CAPSULE:  sh.capsule.half_height = params.capsule.half_height;
	                     sh.capsule.radius = params.capsule.radius; break;
	case SHAPE_BOX:      sh.box.half_extents = params.box.half_extents; break;
	case SHAPE_HULL:     sh.hull.hull = params.hull.hull;
	                     sh.hull.scale = params.hull.scale; break;
	case SHAPE_CYLINDER: sh.cylinder.half_height = params.cylinder.half_height;
	                     sh.cylinder.radius = params.cylinder.radius; break;
	case SHAPE_MESH:     break; // unreachable, asserted above
	}
	apush(s->shapes, sh);
}

void sensor_set_transform(Sensor sensor, v3 position, quat rotation)
{
	assert(is_valid(position) && "sensor_set_transform: position is NaN/inf");
	SensorInternal* s = (SensorInternal*)(uintptr_t)sensor.id;
	s->position = position;
	float r_m2 = rotation.x*rotation.x + rotation.y*rotation.y + rotation.z*rotation.z + rotation.w*rotation.w;
	s->rotation = (r_m2 < 0.5f) ? quat_identity() : rotation;
}

// Boolean overlap test between two positioned shapes. Routes through the
// public collide_* API so no engine-internal state (warm cache, sat hint,
// world pointer) is touched.
static int sensor_shape_overlap(BodyState* bs_a, ShapeInternal* sh_a, BodyState* bs_b, ShapeInternal* sh_b)
{
	// Canonicalise so sh_a->type <= sh_b->type (matches collide_* API).
	if (sh_a->type > sh_b->type) {
		ShapeInternal* ts = sh_a; sh_a = sh_b; sh_b = ts;
		BodyState* tb = bs_a; bs_a = bs_b; bs_b = tb;
	}
	Manifold m = {0};
	switch (sh_a->type) {
	case SHAPE_SPHERE: {
		Sphere a = make_sphere(bs_a, sh_a);
		switch (sh_b->type) {
		case SHAPE_SPHERE:   return collide_sphere_sphere (a, make_sphere(bs_b, sh_b),      &m);
		case SHAPE_CAPSULE:  return collide_sphere_capsule(a, make_capsule(bs_b, sh_b),     &m);
		case SHAPE_BOX:      return collide_sphere_box    (a, make_box(bs_b, sh_b),         &m);
		case SHAPE_HULL:     return collide_sphere_hull   (a, make_convex_hull(bs_b, sh_b), &m);
		case SHAPE_CYLINDER: return collide_cylinder_sphere(make_cylinder(bs_b, sh_b), a,   &m);
		default: return 0;
		}
	}
	case SHAPE_CAPSULE: {
		Capsule a = make_capsule(bs_a, sh_a);
		switch (sh_b->type) {
		case SHAPE_CAPSULE:  return collide_capsule_capsule(a, make_capsule(bs_b, sh_b),     &m);
		case SHAPE_BOX:      return collide_capsule_box    (a, make_box(bs_b, sh_b),         &m);
		case SHAPE_HULL:     return collide_capsule_hull   (a, make_convex_hull(bs_b, sh_b), &m);
		case SHAPE_CYLINDER: return collide_cylinder_capsule(make_cylinder(bs_b, sh_b), a,   &m);
		default: return 0;
		}
	}
	case SHAPE_BOX: {
		Box ab = make_box(bs_a, sh_a);
		switch (sh_b->type) {
		case SHAPE_BOX:      return collide_box_box(ab, make_box(bs_b, sh_b), &m);
		case SHAPE_HULL:     return collide_hull_hull(
		                         (ConvexHull){ hull_unit_box(), ab.center, ab.rotation, ab.half_extents },
		                         make_convex_hull(bs_b, sh_b), &m);
		case SHAPE_CYLINDER: return collide_cylinder_box(make_cylinder(bs_b, sh_b), ab, &m);
		default: return 0;
		}
	}
	case SHAPE_HULL: {
		ConvexHull a = make_convex_hull(bs_a, sh_a);
		switch (sh_b->type) {
		case SHAPE_HULL:     return collide_hull_hull(a, make_convex_hull(bs_b, sh_b), &m);
		case SHAPE_CYLINDER: return collide_cylinder_hull(make_cylinder(bs_b, sh_b), a, &m);
		default: return 0;
		}
	}
	case SHAPE_CYLINDER: {
		if (sh_b->type == SHAPE_CYLINDER)
			return collide_cylinder_cylinder(make_cylinder(bs_a, sh_a), make_cylinder(bs_b, sh_b), &m);
		return 0;
	}
	default: return 0;
	}
}

int sensor_query(World world, Sensor sensor, Body* results, int max_results)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SensorInternal* s = (SensorInternal*)(uintptr_t)sensor.id;

	if (asize(s->shapes) == 0) return 0;

	// Build a transient BodyState so shape_aabb / make_* helpers work directly.
	BodyState sen_state = {0};
	sen_state.position = s->position;
	sen_state.rotation = s->rotation;

	// Union AABB of all sensor shapes.
	AABB sen_box = shape_aabb(&sen_state, &s->shapes[0]);
	for (int i = 1; i < asize(s->shapes); i++)
		sen_box = aabb_merge(sen_box, shape_aabb(&sen_state, &s->shapes[i]));

	CK_DYNA int* candidates = NULL;
	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_query_aabb(w->bvh_dynamic,  sen_box, &candidates);
		bvh_query_aabb(w->bvh_static,   sen_box, &candidates);
		bvh_query_aabb(w->bvh_sleeping, sen_box, &candidates);
	} else {
		int count = asize(w->body_hot);
		for (int i = 0; i < count; i++) {
			if (!split_alive(w->body_gen, i)) continue;
			if (asize(w->body_cold[i].shapes) == 0) continue;
			apush(candidates, i);
		}
	}

	int total = 0;
	int nshapes = asize(s->shapes);
	for (int ci = 0; ci < asize(candidates); ci++) {
		int idx = candidates[ci];
		BodyCold* bc = &w->body_cold[idx];
		BodyState* bs = &w->body_state[idx];

		// Collision filter (both sides must opt in).
		if (!((s->collision_group & bc->collision_mask) && (bc->collision_group & s->collision_mask)))
			continue;

		// Quick AABB reject (broadphase may be conservative, N^2 path never culls).
		AABB body_box = body_aabb(bs, bc);
		if (!aabb_overlaps(sen_box, body_box)) continue;

		int hit = 0;
		int body_nshapes = asize(bc->shapes);
		for (int sa = 0; sa < nshapes && !hit; sa++) {
			for (int sb = 0; sb < body_nshapes && !hit; sb++) {
				if (bc->shapes[sb].type == SHAPE_MESH) continue;
				if (sensor_shape_overlap(&sen_state, &s->shapes[sa], bs, &bc->shapes[sb]))
					hit = 1;
			}
		}
		if (hit) {
			if (total < max_results) results[total] = split_handle(Body, w->body_gen, idx);
			total++;
		}
	}
	afree(candidates);
	return total;
}
