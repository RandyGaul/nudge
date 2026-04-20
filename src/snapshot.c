// See LICENSE for licensing info.
// snapshot.c -- world save/load using SV (serialize.c).
//
// File format (all versioned via SV_Context):
//   [SV header: magic + version]
//   WorldParams
//   int body_count            -- alive bodies only
//   SavedBody[body_count]     -- indexed 0..body_count-1 = "saved index"
//   int joint_count
//   SavedJoint[joint_count]   -- joint bodies reference saved indices
//   SV_DETERMINISTIC_WORLD block:
//     int frame, ldl_topo_version, joint_pairs_version
//     int island_count
//     Island[island_count]        (scalar fields; LDL cache rebuilt on step)
//     uint32_t island_gen[]
//     CK_DYNA int island_free[]
//     int warm_count
//     uint64_t warm_keys[]
//     WarmManifold warm_vals[]
//     int prev_touching_count
//     uint64_t pt_keys[]; uint8_t pt_vals[]
//     int joint_pairs_count
//     uint64_t jp_keys[]; uint8_t jp_vals[]
//
// On load the world is created fresh, bodies + joints are recreated through
// the public API (so shape arrays, mass props, filters, and BVH leaves get
// set up correctly), then the SV_DETERMINISTIC_WORLD block is patched onto
// WorldInternal directly. Subsequent world_step calls produce bit-identical
// state to the saved run.
//
// Not captured (regenerated on next step):
//   - BVH trees (rebuilt from body AABBs during body_add_shape on load)
//   - LDL_Cache per island (rebuilt on next ldl_factor)
//   - EPA manifold cache (use SAT backend if you need determinism across save/load)
//   - Rewind ring (caller re-enables after load if desired)
//   - Sensors (caller-owned, outside the world)
//
// Restrictions:
//   - SHAPE_HULL and SHAPE_MESH are not supported yet (their Hull*/TriMesh*
//     pointers reference caller-owned data with no stable identity on disk).

static ShapeParams snapshot_shape_from_internal(const ShapeInternal* sh)
{
	ShapeParams sp = {0};
	sp.type = sh->type;
	sp.local_pos = sh->local_pos;
	sp.local_rot = sh->local_rot;
	switch (sh->type) {
	case SHAPE_SPHERE:
		sp.sphere.radius = sh->sphere.radius;
		break;
	case SHAPE_CAPSULE:
		sp.capsule.half_height = sh->capsule.half_height;
		sp.capsule.radius = sh->capsule.radius;
		break;
	case SHAPE_BOX:
		sp.box.half_extents = sh->box.half_extents;
		break;
	case SHAPE_HULL:
		sp.hull.hull = sh->hull.hull;   // SV_SERIALIZABLE reads the name field
		sp.hull.scale = sh->hull.scale;
		break;
	case SHAPE_MESH:
		sp.mesh.mesh = sh->mesh.mesh;
		break;
	case SHAPE_HEIGHTFIELD:
		sp.heightfield.hf = sh->heightfield.hf;
		break;
	}
	return sp;
}

// Deserialized shapes (from SV_SERIALIZABLE(ShapeParams) for SHAPE_HULL /
// SHAPE_MESH) carry a sinterned name in the pointer field. Resolve against
// the world's asset registry before the shape is applied (create_body +
// body_add_shape read the hull to compute AABBs / inertia).
static void snapshot_resolve_shape_name(WorldInternal* w, ShapeParams* sp)
{
	if (sp->type == SHAPE_HULL) {
		const char* name = (const char*)(uintptr_t)sp->hull.hull;
		if (!name) return;
		const Hull** slot = map_get_ptr(w->hull_registry, (uint64_t)(uintptr_t)name);
		assert(slot && "world_load_snapshot: SHAPE_HULL name not registered (use world_register_hull)");
		sp->hull.hull = *slot;
	} else if (sp->type == SHAPE_MESH) {
		const char* name = (const char*)(uintptr_t)sp->mesh.mesh;
		if (!name) return;
		const TriMesh** slot = map_get_ptr(w->mesh_registry, (uint64_t)(uintptr_t)name);
		assert(slot && "world_load_snapshot: SHAPE_MESH name not registered (use world_register_mesh)");
		sp->mesh.mesh = *slot;
	} else if (sp->type == SHAPE_HEIGHTFIELD) {
		const char* name = (const char*)(uintptr_t)sp->heightfield.hf;
		if (!name) return;
		const Heightfield** slot = map_get_ptr(w->heightfield_registry, (uint64_t)(uintptr_t)name);
		assert(slot && "world_load_snapshot: SHAPE_HEIGHTFIELD name not registered (use world_register_heightfield)");
		sp->heightfield.hf = *slot;
	}
}

int world_save_snapshot(World world, const char* path)
{
	WorldInternal* w = (WorldInternal*)world.id;
	SV_SAVE_BEGIN(path);
	if (!sv_ok(S)) return 0;

	// World parameters
	WorldParams wp = {
		.gravity = w->gravity,
		.broadphase = w->broadphase_type,
		.narrowphase_backend = (NarrowphaseBackend)w->narrowphase_backend,
		.solver_type = w->solver_type,
		.velocity_iters = w->velocity_iters,
		.position_iters = w->position_iters,
		.contact_hertz = w->contact_hertz,
		.contact_damping_ratio = w->contact_damping_ratio,
		.max_push_velocity = w->max_push_velocity,
		.sub_steps = w->sub_steps,
	};
	SV_ADD_LOCAL(SV_INITIAL, wp);

	// Map live body index -> saved index (0..N-1 over alive bodies only).
	int live_count = asize(w->body_hot);
	int* body_saved_idx = (int*)CK_ALLOC((size_t)(live_count > 0 ? live_count : 1) * sizeof(int));
	int saved_count = 0;
	for (int i = 0; i < live_count; i++) {
		if (split_alive(w->body_gen, i)) body_saved_idx[i] = saved_count++;
		else body_saved_idx[i] = -1;
	}
	SV_ADD_LOCAL(SV_INITIAL, saved_count);

	for (int i = 0; i < live_count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		SavedBody sb = {0};
		sb.params.position      = body_pos(w, i);
		sb.params.rotation      = body_rot(w, i);
		sb.params.mass          = w->body_cold[i].mass;
		sb.params.friction      = body_friction(w, i);
		sb.params.restitution   = body_restitution(w, i);
		sb.params.linear_damping  = body_lin_damping(w, i);
		sb.params.angular_damping = body_ang_damping(w, i);
		sb.velocity          = body_vel(w, i);
		sb.angular_velocity  = body_angvel(w, i);
		sb.collision_group   = w->body_cold[i].collision_group;
		sb.collision_mask    = w->body_cold[i].collision_mask;
		sb.compound_id       = w->body_cold[i].compound_id;
		sb.sleep_time        = body_sleep_time(w, i);
		sb.sleep_allowed     = body_sleep_allowed(w, i);
		sb.island_id         = w->body_cold[i].island_id;
		sb.island_prev       = w->body_cold[i].island_prev;
		sb.island_next       = w->body_cold[i].island_next;
		sb.material_id       = w->body_cold[i].material_id;
		int ns = asize(w->body_cold[i].shapes);
		for (int s = 0; s < ns; s++)
			apush(sb.shapes, snapshot_shape_from_internal(&w->body_cold[i].shapes[s]));
		SV_ADD_LOCAL(SV_INITIAL, sb);
		afree(sb.shapes);
	}

	// Map live joint index -> saved index for internal references.
	int live_joints = asize(w->joints);
	int* joint_saved_idx = (int*)CK_ALLOC((size_t)(live_joints > 0 ? live_joints : 1) * sizeof(int));
	int active_joints = 0;
	for (int j = 0; j < live_joints; j++) {
		if (w->joint_gen[j] & 1) joint_saved_idx[j] = active_joints++;
		else joint_saved_idx[j] = -1;
	}
	SV_ADD_LOCAL(SV_INITIAL, active_joints);

	for (int j = 0; j < live_joints; j++) {
		if (!(w->joint_gen[j] & 1)) continue;
		JointInternal* ji = &w->joints[j];
		SavedJoint sj = {0};
		sj.type = (int)ji->type;
		uint64_t sa = (body_saved_idx[ji->body_a] < 0) ? 0 : (uint64_t)body_saved_idx[ji->body_a];
		uint64_t sb_ = (body_saved_idx[ji->body_b] < 0) ? 0 : (uint64_t)body_saved_idx[ji->body_b];
		switch (ji->type) {
		case JOINT_BALL_SOCKET:
			sj.ball_socket.body_a.id = sa; sj.ball_socket.body_b.id = sb_;
			sj.ball_socket.local_offset_a = ji->ball_socket.local_a;
			sj.ball_socket.local_offset_b = ji->ball_socket.local_b;
			sj.ball_socket.spring = ji->ball_socket.spring;
			break;
		case JOINT_DISTANCE:
			sj.distance.body_a.id = sa; sj.distance.body_b.id = sb_;
			sj.distance.local_offset_a = ji->distance.local_a;
			sj.distance.local_offset_b = ji->distance.local_b;
			sj.distance.rest_length = ji->distance.rest_length;
			sj.distance.spring = ji->distance.spring;
			break;
		case JOINT_HINGE:
			sj.hinge.body_a.id = sa; sj.hinge.body_b.id = sb_;
			sj.hinge.local_offset_a = ji->hinge.local_a;
			sj.hinge.local_offset_b = ji->hinge.local_b;
			sj.hinge.local_axis_a = ji->hinge.local_axis_a;
			sj.hinge.local_axis_b = ji->hinge.local_axis_b;
			sj.hinge.spring = ji->hinge.spring;
			break;
		case JOINT_FIXED:
			sj.fixed.body_a.id = sa; sj.fixed.body_b.id = sb_;
			sj.fixed.local_offset_a = ji->fixed.local_a;
			sj.fixed.local_offset_b = ji->fixed.local_b;
			sj.fixed.spring = ji->fixed.spring;
			break;
		case JOINT_PRISMATIC:
			sj.prismatic.body_a.id = sa; sj.prismatic.body_b.id = sb_;
			sj.prismatic.local_offset_a = ji->prismatic.local_a;
			sj.prismatic.local_offset_b = ji->prismatic.local_b;
			sj.prismatic.local_axis_a = ji->prismatic.local_axis_a;
			sj.prismatic.local_axis_b = ji->prismatic.local_axis_b;
			sj.prismatic.spring = ji->prismatic.spring;
			break;
		case JOINT_ANGULAR_MOTOR:
			sj.angular_motor.body_a.id = sa; sj.angular_motor.body_b.id = sb_;
			sj.angular_motor.local_axis_a = ji->angular_motor.local_axis_a;
			sj.angular_motor.local_axis_b = ji->angular_motor.local_axis_b;
			sj.angular_motor.target_speed = ji->angular_motor.target_speed;
			sj.angular_motor.max_impulse = ji->angular_motor.max_impulse;
			break;
		case JOINT_TWIST_LIMIT:
			sj.twist_limit.body_a.id = sa; sj.twist_limit.body_b.id = sb_;
			sj.twist_limit.local_axis_a = ji->twist_limit.local_axis_a;
			sj.twist_limit.local_axis_b = ji->twist_limit.local_axis_b;
			sj.twist_limit.limit_min = ji->twist_limit.limit_min;
			sj.twist_limit.limit_max = ji->twist_limit.limit_max;
			sj.twist_limit.spring = ji->twist_limit.spring;
			break;
		case JOINT_CONE_LIMIT:
			sj.cone_limit.body_a.id = sa; sj.cone_limit.body_b.id = sb_;
			sj.cone_limit.local_axis_a = ji->cone_limit.local_axis_a;
			sj.cone_limit.local_axis_b = ji->cone_limit.local_axis_b;
			sj.cone_limit.half_angle = ji->cone_limit.half_angle;
			sj.cone_limit.spring = ji->cone_limit.spring;
			break;
		case JOINT_SWING_TWIST:
			sj.swing_twist.body_a.id = sa; sj.swing_twist.body_b.id = sb_;
			sj.swing_twist.local_offset_a = ji->swing_twist.local_a;
			sj.swing_twist.local_offset_b = ji->swing_twist.local_b;
			sj.swing_twist.local_axis_a = ji->swing_twist.local_axis_a;
			sj.swing_twist.local_axis_b = ji->swing_twist.local_axis_b;
			sj.swing_twist.cone_half_angle = ji->swing_twist.cone_half_angle;
			sj.swing_twist.twist_min = ji->swing_twist.twist_min;
			sj.swing_twist.twist_max = ji->swing_twist.twist_max;
			sj.swing_twist.spring = ji->swing_twist.spring;
			break;
		}
		sj.island_id   = ji->island_id;
		sj.island_prev = ji->island_prev;
		sj.island_next = ji->island_next;
		for (int d = 0; d < JOINT_MAX_DOF; d++) sj.warm_lambda[d] = w->joint_hot[j].warm_lambda[d];
		SV_ADD_LOCAL(SV_INITIAL, sj);
	}

	// SV_DETERMINISTIC_WORLD: internal state needed for bit-identical replay.
	// Guarded per-field by its VERSION in SV_ADD_LOCAL so older files simply
	// skip this block.
	int frame              = w->frame;
	int ldl_topo_version   = w->ldl_topo_version;
	int joint_pairs_version = w->joint_pairs_version;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, frame);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, ldl_topo_version);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, joint_pairs_version);

	// Islands: scalar fields + gen array + free list.
	int island_count = asize(w->islands);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, island_count);
	for (int i = 0; i < island_count; i++) {
		Island isl = w->islands[i];
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, isl);
	}
	for (int i = 0; i < island_count; i++) {
		uint32_t g = w->island_gen[i];
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, g);
	}
	int island_free_count = asize(w->island_free);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, island_free_count);
	for (int i = 0; i < island_free_count; i++) {
		int x = w->island_free[i];
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, x);
	}

	// Warm cache: dense (key, value) parallel arrays.
	int warm_count = map_size(w->warm_cache);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, warm_count);
	{
		uint64_t* keys = map_keys(w->warm_cache);
		for (int i = 0; i < warm_count; i++) {
			uint64_t k = keys[i];
			WarmManifold v = w->warm_cache[i];
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		}
	}

	int pt_count = map_size(w->prev_touching);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, pt_count);
	{
		uint64_t* keys = map_keys(w->prev_touching);
		for (int i = 0; i < pt_count; i++) {
			uint64_t k = keys[i];
			uint8_t v = w->prev_touching[i];
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		}
	}

	int jp_count = map_size(w->joint_pairs);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, jp_count);
	{
		uint64_t* keys = map_keys(w->joint_pairs);
		for (int i = 0; i < jp_count; i++) {
			uint64_t k = keys[i];
			uint8_t v = w->joint_pairs[i];
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
			SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		}
	}

	// Sensors: simple list. No index remapping needed -- sensor indices don't
	// appear in body/joint data.
	int sensor_live = asize(w->sensors);
	int sensor_count = 0;
	for (int i = 0; i < sensor_live; i++) if (split_alive(w->sensor_gen, i)) sensor_count++;
	SV_ADD_LOCAL(SV_SENSORS_IN_SAVE, sensor_count);
	for (int i = 0; i < sensor_live; i++) {
		if (!split_alive(w->sensor_gen, i)) continue;
		SensorInternal* ss = &w->sensors[i];
		SavedSensor ssv = {0};
		ssv.position = ss->position;
		ssv.rotation = ss->rotation;
		ssv.collision_group = ss->collision_group;
		ssv.collision_mask  = ss->collision_mask;
		int ns = asize(ss->shapes);
		for (int k = 0; k < ns; k++)
			apush(ssv.shapes, snapshot_shape_from_internal(&ss->shapes[k]));
		SV_ADD_LOCAL(SV_SENSORS_IN_SAVE, ssv);
		afree(ssv.shapes);
	}

	// Material palette: 256 fixed entries, written verbatim. Older files
	// simply don't have this block (SV_ADD_LOCAL skips when version < tag).
	for (int i = 0; i < 256; i++) {
		Material m = w->materials[i];
		SV_ADD_LOCAL(SV_MATERIALS, m);
	}

	CK_FREE(body_saved_idx);
	CK_FREE(joint_saved_idx);
	SV_SAVE_END();
	return 1;
}

// Shared body of world_load_snapshot + world_load_snapshot_into. If
// target.id == 0 creates a new world from the file's WorldParams; otherwise
// uses the existing target world and discards file's params.
static World snapshot_load_impl(SV_Context* S, World target)
{
	WorldParams wp = {0};
	SV_ADD_LOCAL(SV_INITIAL, wp);

	World world = target.id ? target : create_world(wp);
	WorldInternal* w = (WorldInternal*)world.id;

	int saved_count = 0;
	SV_ADD_LOCAL(SV_INITIAL, saved_count);

	// Bodies: create in saved order so live indices match saved indices.
	Body* body_table = (Body*)CK_ALLOC((size_t)(saved_count > 0 ? saved_count : 1) * sizeof(Body));
	SavedBody* saved_bodies = (SavedBody*)CK_ALLOC((size_t)(saved_count > 0 ? saved_count : 1) * sizeof(SavedBody));
	for (int i = 0; i < saved_count; i++) {
		SavedBody sb = {0};
		SV_ADD_LOCAL(SV_INITIAL, sb);
		Body b = create_body(world, sb.params);
		for (int s = 0; s < asize(sb.shapes); s++) {
			snapshot_resolve_shape_name(w, &sb.shapes[s]);
			body_add_shape(world, b, sb.shapes[s]);
		}
		body_set_velocity(world, b, sb.velocity);
		body_set_angular_velocity(world, b, sb.angular_velocity);
		body_set_collision_filter(world, b, sb.collision_group, sb.collision_mask);
		body_set_compound_id(world, b, sb.compound_id);
		body_set_sleep_allowed(world, b, sb.sleep_allowed);
		body_table[i] = b;
		saved_bodies[i] = sb;   // retain for post-create patching
		saved_bodies[i].shapes = NULL;
		afree(sb.shapes);
	}

	int active_joints = 0;
	SV_ADD_LOCAL(SV_INITIAL, active_joints);

	SavedJoint* saved_joints = (SavedJoint*)CK_ALLOC((size_t)(active_joints > 0 ? active_joints : 1) * sizeof(SavedJoint));
	Joint* joint_table = (Joint*)CK_ALLOC((size_t)(active_joints > 0 ? active_joints : 1) * sizeof(Joint));
	for (int j = 0; j < active_joints; j++) {
		SavedJoint sj = {0};
		SV_ADD_LOCAL(SV_INITIAL, sj);
		#define REMAP_BODY(P) do { \
			int sa = (int)sj.P.body_a.id; int sb2 = (int)sj.P.body_b.id; \
			sj.P.body_a = (sa >= 0 && sa < saved_count) ? body_table[sa] : (Body){0}; \
			sj.P.body_b = (sb2 >= 0 && sb2 < saved_count) ? body_table[sb2] : (Body){0}; \
		} while (0)
		Joint handle = {0};
		switch (sj.type) {
		case JOINT_BALL_SOCKET:   REMAP_BODY(ball_socket);   handle = create_ball_socket(world, sj.ball_socket); break;
		case JOINT_DISTANCE:      REMAP_BODY(distance);      handle = create_distance(world, sj.distance); break;
		case JOINT_HINGE:         REMAP_BODY(hinge);         handle = create_hinge(world, sj.hinge); break;
		case JOINT_FIXED:         REMAP_BODY(fixed);         handle = create_fixed(world, sj.fixed); break;
		case JOINT_PRISMATIC:     REMAP_BODY(prismatic);     handle = create_prismatic(world, sj.prismatic); break;
		case JOINT_ANGULAR_MOTOR: REMAP_BODY(angular_motor); handle = create_angular_motor(world, sj.angular_motor); break;
		case JOINT_TWIST_LIMIT:   REMAP_BODY(twist_limit);   handle = create_twist_limit(world, sj.twist_limit); break;
		case JOINT_CONE_LIMIT:    REMAP_BODY(cone_limit);    handle = create_cone_limit(world, sj.cone_limit); break;
		case JOINT_SWING_TWIST:   REMAP_BODY(swing_twist);   handle = create_swing_twist(world, sj.swing_twist); break;
		}
		#undef REMAP_BODY
		joint_table[j] = handle;
		saved_joints[j] = sj;
	}

	// SV_DETERMINISTIC_WORLD: only present in files >= that version. On older
	// files these fields remain at their defaults (zero) and the bodies will
	// settle via fresh simulation -- non-bit-identical but still correct.
	int frame               = 0;
	int ldl_topo_version    = 0;
	int joint_pairs_version = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, frame);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, ldl_topo_version);
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, joint_pairs_version);

	int island_count = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, island_count);

	// Free any islands create_body / create_* side-effects produced and
	// resize to the saved count. At this point no world_step has run so
	// islands should be empty, but LDL caches and the free list need cleanup.
	for (int i = 0; i < asize(w->islands); i++) ldl_cache_free(&w->islands[i].ldl);
	if (w->islands) asetlen(w->islands, 0);
	if (w->island_gen) asetlen(w->island_gen, 0);
	if (w->island_free) aclear(w->island_free);

	if (S->loading && S->version >= SV_DETERMINISTIC_WORLD && island_count > 0) {
		afit(w->islands, island_count); asetlen(w->islands, island_count);
		afit(w->island_gen, island_count); asetlen(w->island_gen, island_count);
		memset(w->islands, 0, (size_t)island_count * sizeof(Island));
	}

	for (int i = 0; i < island_count; i++) {
		Island isl = {0};
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, isl);
		// Overwrite scalar fields on the (already-zeroed) live island.
		w->islands[i].head_body             = isl.head_body;
		w->islands[i].tail_body             = isl.tail_body;
		w->islands[i].body_count            = isl.body_count;
		w->islands[i].head_joint            = isl.head_joint;
		w->islands[i].tail_joint            = isl.tail_joint;
		w->islands[i].joint_count           = isl.joint_count;
		w->islands[i].constraint_remove_count = isl.constraint_remove_count;
		w->islands[i].awake                 = isl.awake;
	}
	for (int i = 0; i < island_count; i++) {
		uint32_t g = 0;
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, g);
		w->island_gen[i] = g;
	}
	int island_free_count = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, island_free_count);
	for (int i = 0; i < island_free_count; i++) {
		int x = 0;
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, x);
		apush(w->island_free, x);
	}

	int warm_count = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, warm_count);
	map_clear(w->warm_cache);
	for (int i = 0; i < warm_count; i++) {
		uint64_t k = 0; WarmManifold v = {0};
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		map_set(w->warm_cache, k, v);
	}

	int pt_count = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, pt_count);
	map_clear(w->prev_touching);
	for (int i = 0; i < pt_count; i++) {
		uint64_t k = 0; uint8_t v = 0;
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		map_set(w->prev_touching, k, v);
	}

	int jp_count = 0;
	SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, jp_count);
	map_clear(w->joint_pairs);
	for (int i = 0; i < jp_count; i++) {
		uint64_t k = 0; uint8_t v = 0;
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, k);
		SV_ADD_LOCAL(SV_DETERMINISTIC_WORLD, v);
		map_set(w->joint_pairs, k, v);
	}

	// Apply body + joint patches (sleep state, island linkage, warm lambda).
	// Gated by version: SV_ADD(SV_DETERMINISTIC_WORLD, ...) read zero when
	// the file predates the feature.
	if (S->loading && S->version >= SV_DETERMINISTIC_WORLD) {
		for (int i = 0; i < saved_count; i++) {
			int bi = handle_index(body_table[i]);
			body_sleep_time(w, bi)    = saved_bodies[i].sleep_time;
			body_sleep_allowed(w, bi) = saved_bodies[i].sleep_allowed;
			w->body_cold[bi].island_id   = saved_bodies[i].island_id;
			w->body_cold[bi].island_prev = saved_bodies[i].island_prev;
			w->body_cold[bi].island_next = saved_bodies[i].island_next;
		}
		for (int j = 0; j < active_joints; j++) {
			int ji = handle_index(joint_table[j]);
			w->joints[ji].island_id   = saved_joints[j].island_id;
			w->joints[ji].island_prev = saved_joints[j].island_prev;
			w->joints[ji].island_next = saved_joints[j].island_next;
			for (int d = 0; d < JOINT_MAX_DOF; d++)
				w->joint_hot[ji].warm_lambda[d] = saved_joints[j].warm_lambda[d];
		}
		w->frame               = frame;
		w->ldl_topo_version    = ldl_topo_version;
		w->joint_pairs_version = joint_pairs_version;
	}

	// Sensors block (SV_SENSORS_IN_SAVE).
	int sensor_count = 0;
	SV_ADD_LOCAL(SV_SENSORS_IN_SAVE, sensor_count);
	for (int i = 0; i < sensor_count; i++) {
		SavedSensor ssv = {0};
		SV_ADD_LOCAL(SV_SENSORS_IN_SAVE, ssv);
		Sensor s = create_sensor(world, (SensorParams){
			.position = ssv.position, .rotation = ssv.rotation,
			.collision_group = ssv.collision_group, .collision_mask = ssv.collision_mask,
		});
		for (int k = 0; k < asize(ssv.shapes); k++) {
			snapshot_resolve_shape_name(w, &ssv.shapes[k]);
			sensor_add_shape(world, s, ssv.shapes[k]);
		}
		afree(ssv.shapes);
	}

	// Material palette + per-body material_id. Older files skip this block
	// entirely (SV_ADD_LOCAL is a no-op when version < tag); palette entries
	// keep their create_world defaults and body material_ids stay 0.
	if (S->version >= SV_MATERIALS) {
		for (int i = 0; i < saved_count; i++) {
			int bi = handle_index(body_table[i]);
			w->body_cold[bi].material_id = saved_bodies[i].material_id;
		}
	}
	for (int i = 0; i < 256; i++) {
		Material m = {0};
		SV_ADD_LOCAL(SV_MATERIALS, m);
		if (S->version >= SV_MATERIALS) w->materials[i] = m;
	}

	CK_FREE(body_table);
	CK_FREE(joint_table);
	CK_FREE(saved_bodies);
	CK_FREE(saved_joints);
	return world;
}

World world_load_snapshot(const char* path)
{
	SV_LOAD_BEGIN(path);
	if (!sv_ok(S)) { SV_LOAD_END(); return (World){0}; }
	World w = snapshot_load_impl(S, (World){0});
	SV_LOAD_END();
	return w;
}

int world_load_snapshot_into(World target, const char* path)
{
	assert(target.id && "world_load_snapshot_into: target world is null");
	SV_LOAD_BEGIN(path);
	if (!sv_ok(S)) { SV_LOAD_END(); return 0; }
	snapshot_load_impl(S, target);
	SV_LOAD_END();
	return 1;
}
