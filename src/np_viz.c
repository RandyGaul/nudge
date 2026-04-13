// np_viz.c -- Narrowphase Visualization Suite.
// Two shapes, run narrowphase, render full manifold, 3D manipulators.
// All behavioral code in one file (locality rule).

// ----------------------------------------------------------------------------
// Types.

typedef enum
{
	NPV_SPHERE,
	NPV_CAPSULE,
	NPV_BOX,
	NPV_CYLINDER,
	NPV_HULL_TETRA,
	NPV_HULL_OCTA,
	NPV_HULL_ICOSA,
	NPV_HULL_DODECA,
	NPV_HULL_RANDOM,
	NPV_SHAPE_COUNT,
} NPV_ShapeKind;

static const char* s_npv_shape_names = "Sphere\0Capsule\0Box\0Cylinder\0Tetrahedron\0Octahedron\0Icosahedron\0Dodecahedron\0Random Hull\0";

typedef struct NPV_Shape
{
	int kind; // NPV_ShapeKind as int for ImGui combo
	v3 pos;
	quat rot;
	float radius;
	float half_height;
	v3 half_extents;
	Hull* hull;       // owned, freed on change (NULL for non-hull types)
	int rand_verts;   // vertex count for random hull
	uint32_t rand_seed;
	int mesh_slot;    // allocated custom mesh slot index
	int mesh_valid;   // 1 if custom mesh has been generated
	v3 color;
} NPV_Shape;

// ----------------------------------------------------------------------------
// State.

static NPV_Shape npv_shapes[2];
static int npv_selected;           // 0 or 1
static int npv_gizmo_mode;        // 0=translate, 1=rotate
static int npv_drag_axis;         // -1=none, 0=X, 1=Y, 2=Z
static int npv_dragging;
static v3 npv_drag_plane_normal;
static v3 npv_drag_plane_point;
static float npv_drag_start_proj;
static v3 npv_drag_start_pos;
static float npv_drag_start_angle;
static quat npv_drag_start_rot;
static Manifold npv_manifold;
static int npv_colliding;
static int npv_initialized;
static float npv_gizmo_size = 1.5f;

static const v3 npv_axis_dirs[3] = { {1,0,0,0}, {0,1,0,0}, {0,0,1,0} };
static const v3 npv_axis_colors[3] = { {0.9f,0.2f,0.2f,0}, {0.2f,0.9f,0.2f,0}, {0.3f,0.3f,1.0f,0} };
static const v3 npv_axis_colors_bright[3] = { {1,0.5f,0.5f,0}, {0.5f,1,0.5f,0}, {0.5f,0.5f,1,0} };

// ----------------------------------------------------------------------------
// Hull preset builders.

static Hull* npv_build_tetrahedron()
{
	float s = 1.0f / sqrtf(3.0f);
	v3 pts[] = { V3(s,s,s), V3(s,-s,-s), V3(-s,s,-s), V3(-s,-s,s) };
	return quickhull(pts, 4);
}

static Hull* npv_build_octahedron()
{
	v3 pts[] = { V3(1,0,0), V3(-1,0,0), V3(0,1,0), V3(0,-1,0), V3(0,0,1), V3(0,0,-1) };
	return quickhull(pts, 6);
}

static Hull* npv_build_icosahedron()
{
	float p = 1.618033988f;
	float s = 1.0f / sqrtf(1.0f + p * p);
	float ps = p * s;
	v3 pts[] = {
		V3(-s, ps, 0), V3(s, ps, 0), V3(-s, -ps, 0), V3(s, -ps, 0),
		V3(0, -s, ps), V3(0, s, ps), V3(0, -s, -ps), V3(0, s, -ps),
		V3(ps, 0, -s), V3(ps, 0, s), V3(-ps, 0, -s), V3(-ps, 0, s),
	};
	return quickhull(pts, 12);
}

static Hull* npv_build_dodecahedron()
{
	float p = 1.618033988f;
	float ip = 1.0f / p;
	float s = 1.0f / sqrtf(3.0f);
	v3 pts[] = {
		V3(s,s,s), V3(s,s,-s), V3(s,-s,s), V3(s,-s,-s),
		V3(-s,s,s), V3(-s,s,-s), V3(-s,-s,s), V3(-s,-s,-s),
		V3(0, ip*s, p*s), V3(0, ip*s, -p*s), V3(0, -ip*s, p*s), V3(0, -ip*s, -p*s),
		V3(ip*s, p*s, 0), V3(ip*s, -p*s, 0), V3(-ip*s, p*s, 0), V3(-ip*s, -p*s, 0),
		V3(p*s, 0, ip*s), V3(p*s, 0, -ip*s), V3(-p*s, 0, ip*s), V3(-p*s, 0, -ip*s),
	};
	return quickhull(pts, 20);
}

static Hull* npv_build_random_hull(int n, uint32_t seed)
{
	if (n < 4) n = 4;
	if (n > 128) n = 128;
	v3* pts = NULL;
	afit(pts, n);
	uint32_t st = seed ? seed : 12345;
	for (int i = 0; i < n; i++) {
		st ^= st << 13; st ^= st >> 17; st ^= st << 5;
		float u = (float)(st & 0xFFFF) / 65535.0f;
		st ^= st << 13; st ^= st >> 17; st ^= st << 5;
		float v = (float)(st & 0xFFFF) / 65535.0f;
		float theta = 2.0f * 3.14159265f * u;
		float phi = acosf(1.0f - 2.0f * v);
		apush(pts, V3(sinf(phi)*cosf(theta), sinf(phi)*sinf(theta), cosf(phi)));
	}
	Hull* h = quickhull(pts, n);
	afree(pts);
	return h;
}

// ----------------------------------------------------------------------------
// Shape mesh rebuild.

static void npv_free_hull(NPV_Shape* s)
{
	if (s->hull) { hull_free(s->hull); s->hull = NULL; }
}

static void npv_rebuild_hull(NPV_Shape* s)
{
	npv_free_hull(s);
	switch (s->kind) {
	case NPV_HULL_TETRA: s->hull = npv_build_tetrahedron(); break;
	case NPV_HULL_OCTA:  s->hull = npv_build_octahedron();  break;
	case NPV_HULL_ICOSA: s->hull = npv_build_icosahedron(); break;
	case NPV_HULL_DODECA: s->hull = npv_build_dodecahedron(); break;
	case NPV_HULL_RANDOM: s->hull = npv_build_random_hull(s->rand_verts, s->rand_seed); break;
	default: break;
	}
}

static void npv_rebuild_mesh(NPV_Shape* s)
{
	// Sphere and box use built-in mesh types, no custom mesh needed.
	if (s->kind == NPV_SPHERE || s->kind == NPV_BOX) {
		if (s->mesh_valid) { mesh_destroy(&r_meshes[s->mesh_slot]); s->mesh_valid = 0; }
		return;
	}
	// Destroy old custom mesh.
	if (s->mesh_valid) { mesh_destroy(&r_meshes[s->mesh_slot]); s->mesh_valid = 0; }
	Mesh* m = &r_meshes[s->mesh_slot];
	switch (s->kind) {
	case NPV_CAPSULE:
		mesh_generate_capsule(m, s->radius, s->half_height);
		break;
	case NPV_CYLINDER:
		mesh_generate_cylinder(m, s->radius, s->half_height);
		break;
	default: // hull types
		if (s->hull) mesh_generate_hull(m, s->hull, V3(1, 1, 1));
		break;
	}
	r_mesh_ready[s->mesh_slot] = 1;
	s->mesh_valid = 1;
}

static int npv_mesh_type(NPV_Shape* s)
{
	if (s->kind == NPV_SPHERE) return MESH_SPHERE;
	if (s->kind == NPV_BOX) return MESH_BOX;
	return s->mesh_slot;
}

static v3 npv_mesh_scale(NPV_Shape* s)
{
	if (s->kind == NPV_SPHERE) return V3(s->radius, s->radius, s->radius);
	if (s->kind == NPV_BOX) return s->half_extents;
	return V3(1, 1, 1);
}

// ----------------------------------------------------------------------------
// Shape conversion -- build engine shape structs for collision.

static ShapeType npv_engine_type(int kind)
{
	switch (kind) {
	case NPV_SPHERE: return SHAPE_SPHERE;
	case NPV_CAPSULE: return SHAPE_CAPSULE;
	case NPV_BOX: return SHAPE_BOX;
	case NPV_CYLINDER: return SHAPE_CYLINDER;
	default: return SHAPE_HULL;
	}
}

static Sphere npv_make_sphere(NPV_Shape* s)
{
	return (Sphere){ .center = s->pos, .radius = s->radius };
}

static Capsule npv_make_capsule(NPV_Shape* s)
{
	v3 axis = quat_rotate(s->rot, V3(0, s->half_height, 0));
	return (Capsule){ .p = v3_sub(s->pos, axis), .q = v3_add(s->pos, axis), .radius = s->radius };
}

static Box npv_make_box(NPV_Shape* s)
{
	return (Box){ .center = s->pos, .rotation = s->rot, .half_extents = s->half_extents };
}

static Cylinder npv_make_cylinder(NPV_Shape* s)
{
	return (Cylinder){ .center = s->pos, .rotation = s->rot, .half_height = s->half_height, .radius = s->radius };
}

static ConvexHull npv_make_hull(NPV_Shape* s)
{
	return (ConvexHull){ .hull = s->hull, .center = s->pos, .rotation = s->rot, .scale = V3(1, 1, 1) };
}

static ConvexHull npv_box_as_hull(NPV_Shape* s)
{
	return (ConvexHull){ .hull = hull_unit_box(), .center = s->pos, .rotation = s->rot, .scale = s->half_extents };
}

// ----------------------------------------------------------------------------
// Collision dispatch.

static void npv_flip_manifold(Manifold* m)
{
	for (int i = 0; i < m->count; i++)
		m->contacts[i].normal = v3_neg(m->contacts[i].normal);
}

static void npv_run_collide()
{
	npv_manifold = (Manifold){0};
	NPV_Shape* sa = &npv_shapes[0];
	NPV_Shape* sb = &npv_shapes[1];
	ShapeType ta = npv_engine_type(sa->kind);
	ShapeType tb = npv_engine_type(sb->kind);

	// Hull types need a valid hull pointer.
	if (ta == SHAPE_HULL && !sa->hull) { npv_colliding = 0; return; }
	if (tb == SHAPE_HULL && !sb->hull) { npv_colliding = 0; return; }

	int hit = 0;

	if (ta == SHAPE_SPHERE && tb == SHAPE_SPHERE) {
		hit = collide_sphere_sphere(npv_make_sphere(sa), npv_make_sphere(sb), &npv_manifold);
	} else if (ta == SHAPE_SPHERE && tb == SHAPE_CAPSULE) {
		hit = collide_sphere_capsule(npv_make_sphere(sa), npv_make_capsule(sb), &npv_manifold);
	} else if (ta == SHAPE_SPHERE && tb == SHAPE_BOX) {
		hit = collide_sphere_box(npv_make_sphere(sa), npv_make_box(sb), &npv_manifold);
	} else if (ta == SHAPE_SPHERE && tb == SHAPE_HULL) {
		hit = collide_sphere_hull(npv_make_sphere(sa), npv_make_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_SPHERE && tb == SHAPE_CYLINDER) {
		hit = collide_cylinder_sphere(npv_make_cylinder(sb), npv_make_sphere(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_CAPSULE && tb == SHAPE_SPHERE) {
		hit = collide_sphere_capsule(npv_make_sphere(sb), npv_make_capsule(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_CAPSULE && tb == SHAPE_CAPSULE) {
		hit = collide_capsule_capsule(npv_make_capsule(sa), npv_make_capsule(sb), &npv_manifold);
	} else if (ta == SHAPE_CAPSULE && tb == SHAPE_BOX) {
		hit = collide_capsule_box(npv_make_capsule(sa), npv_make_box(sb), &npv_manifold);
	} else if (ta == SHAPE_CAPSULE && tb == SHAPE_HULL) {
		hit = collide_capsule_hull(npv_make_capsule(sa), npv_make_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_CAPSULE && tb == SHAPE_CYLINDER) {
		hit = collide_cylinder_capsule(npv_make_cylinder(sb), npv_make_capsule(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_BOX && tb == SHAPE_SPHERE) {
		hit = collide_sphere_box(npv_make_sphere(sb), npv_make_box(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_BOX && tb == SHAPE_CAPSULE) {
		hit = collide_capsule_box(npv_make_capsule(sb), npv_make_box(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_BOX && tb == SHAPE_BOX) {
		hit = collide_box_box(npv_make_box(sa), npv_make_box(sb), &npv_manifold);
	} else if (ta == SHAPE_BOX && tb == SHAPE_HULL) {
		hit = collide_hull_hull(npv_box_as_hull(sa), npv_make_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_BOX && tb == SHAPE_CYLINDER) {
		hit = collide_cylinder_box(npv_make_cylinder(sb), npv_make_box(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_HULL && tb == SHAPE_SPHERE) {
		hit = collide_sphere_hull(npv_make_sphere(sb), npv_make_hull(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_HULL && tb == SHAPE_CAPSULE) {
		hit = collide_capsule_hull(npv_make_capsule(sb), npv_make_hull(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_HULL && tb == SHAPE_BOX) {
		hit = collide_hull_hull(npv_make_hull(sa), npv_box_as_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_HULL && tb == SHAPE_HULL) {
		hit = collide_hull_hull(npv_make_hull(sa), npv_make_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_HULL && tb == SHAPE_CYLINDER) {
		hit = collide_cylinder_hull(npv_make_cylinder(sb), npv_make_hull(sa), &npv_manifold);
		npv_flip_manifold(&npv_manifold);
	} else if (ta == SHAPE_CYLINDER && tb == SHAPE_SPHERE) {
		hit = collide_cylinder_sphere(npv_make_cylinder(sa), npv_make_sphere(sb), &npv_manifold);
	} else if (ta == SHAPE_CYLINDER && tb == SHAPE_CAPSULE) {
		hit = collide_cylinder_capsule(npv_make_cylinder(sa), npv_make_capsule(sb), &npv_manifold);
	} else if (ta == SHAPE_CYLINDER && tb == SHAPE_BOX) {
		hit = collide_cylinder_box(npv_make_cylinder(sa), npv_make_box(sb), &npv_manifold);
	} else if (ta == SHAPE_CYLINDER && tb == SHAPE_HULL) {
		hit = collide_cylinder_hull(npv_make_cylinder(sa), npv_make_hull(sb), &npv_manifold);
	} else if (ta == SHAPE_CYLINDER && tb == SHAPE_CYLINDER) {
		hit = collide_cylinder_cylinder(npv_make_cylinder(sa), npv_make_cylinder(sb), &npv_manifold);
	}

	npv_colliding = hit;
}

// ----------------------------------------------------------------------------
// Gizmo math helpers.

// Closest approach between a ray and an infinite line.
// Returns the parameter along the line, and writes the closest distance.
static float npv_ray_line_closest(v3 ray_o, v3 ray_d, v3 line_o, v3 line_d, float* out_dist)
{
	v3 w = v3_sub(ray_o, line_o);
	float a = v3_dot(ray_d, ray_d);
	float b = v3_dot(ray_d, line_d);
	float c = v3_dot(line_d, line_d);
	float d = v3_dot(ray_d, w);
	float e = v3_dot(line_d, w);
	float denom = a * c - b * b;
	if (fabsf(denom) < 1e-8f) { *out_dist = 1e10f; return 0; }
	float t_ray = (b * e - c * d) / denom;
	float t_line = (a * e - b * d) / denom;
	if (t_ray < 0) t_ray = 0;
	v3 p_ray = v3_add(ray_o, v3_scale(ray_d, t_ray));
	v3 p_line = v3_add(line_o, v3_scale(line_d, t_line));
	*out_dist = v3_len(v3_sub(p_ray, p_line));
	return t_line;
}

// Ray-plane intersection. Returns t along ray, or -1 if parallel.
static float npv_ray_plane(v3 ray_o, v3 ray_d, v3 plane_pt, v3 plane_n)
{
	float denom = v3_dot(ray_d, plane_n);
	if (fabsf(denom) < 1e-8f) return -1.0f;
	return v3_dot(v3_sub(plane_pt, ray_o), plane_n) / denom;
}

// Build perpendicular vectors to an axis.
static void npv_perp_axes(v3 axis, v3* out_a, v3* out_b)
{
	v3 ref = (fabsf(axis.y) < 0.9f) ? V3(0, 1, 0) : V3(1, 0, 0);
	*out_a = v3_norm(v3_cross(axis, ref));
	*out_b = v3_cross(axis, *out_a);
}

// Quaternion from axis-angle.
static quat npv_quat_axis_angle(v3 axis, float angle)
{
	float s = sinf(angle * 0.5f);
	float c = cosf(angle * 0.5f);
	return (quat){ axis.x * s, axis.y * s, axis.z * s, c };
}

// Euler (radians, XYZ intrinsic) to quaternion.
static quat npv_euler_to_quat(float rx, float ry, float rz)
{
	float cx = cosf(rx * 0.5f), sx = sinf(rx * 0.5f);
	float cy = cosf(ry * 0.5f), sy = sinf(ry * 0.5f);
	float cz = cosf(rz * 0.5f), sz = sinf(rz * 0.5f);
	return (quat){
		sx*cy*cz - cx*sy*sz,
		cx*sy*cz + sx*cy*sz,
		cx*cy*sz - sx*sy*cz,
		cx*cy*cz + sx*sy*sz,
	};
}

// Quaternion to euler (radians, XYZ intrinsic).
static v3 npv_quat_to_euler(quat q)
{
	float sinr = 2.0f * (q.w * q.x + q.y * q.z);
	float cosr = 1.0f - 2.0f * (q.x * q.x + q.y * q.y);
	float rx = atan2f(sinr, cosr);
	float sinp = 2.0f * (q.w * q.y - q.z * q.x);
	float ry = (fabsf(sinp) >= 1.0f) ? copysignf(3.14159265f / 2.0f, sinp) : asinf(sinp);
	float siny = 2.0f * (q.w * q.z + q.x * q.y);
	float cosy = 1.0f - 2.0f * (q.y * q.y + q.z * q.z);
	float rz = atan2f(siny, cosy);
	return V3(rx, ry, rz);
}

// ----------------------------------------------------------------------------
// Gizmo hit testing and dragging.

// Test which translation axis the mouse ray is closest to. Returns 0-2 or -1.
static int npv_pick_translate_axis(v3 ray_o, v3 ray_d, v3 center, float threshold)
{
	int best = -1;
	float best_dist = threshold;
	for (int i = 0; i < 3; i++) {
		float dist;
		float t = npv_ray_line_closest(ray_o, ray_d, center, npv_axis_dirs[i], &dist);
		// Only consider the positive half of the axis (the arrow part).
		if (t < 0 || t > npv_gizmo_size * 1.2f) continue;
		if (dist < best_dist) { best_dist = dist; best = i; }
	}
	return best;
}

// Test which rotation circle the mouse ray is closest to. Returns 0-2 or -1.
static int npv_pick_rotate_axis(v3 ray_o, v3 ray_d, v3 center, float threshold)
{
	int best = -1;
	float best_err = threshold;
	for (int i = 0; i < 3; i++) {
		float t = npv_ray_plane(ray_o, ray_d, center, npv_axis_dirs[i]);
		if (t < 0) continue;
		v3 hit = v3_add(ray_o, v3_scale(ray_d, t));
		float dist = v3_len(v3_sub(hit, center));
		float err = fabsf(dist - npv_gizmo_size);
		if (err < best_err) { best_err = err; best = i; }
	}
	return best;
}

static void npv_begin_translate(v3 ray_o, v3 ray_d, int axis, NPV_Shape* s)
{
	npv_dragging = 1;
	npv_drag_axis = axis;
	npv_drag_start_pos = s->pos;
	// Build drag plane: contains the axis, perpendicular to camera-to-shape direction.
	v3 cam_dir = v3_norm(v3_sub(ray_o, s->pos));
	v3 axis_d = npv_axis_dirs[axis];
	v3 plane_n = v3_norm(v3_cross(axis_d, v3_cross(cam_dir, axis_d)));
	npv_drag_plane_normal = plane_n;
	npv_drag_plane_point = s->pos;
	float t = npv_ray_plane(ray_o, ray_d, npv_drag_plane_point, npv_drag_plane_normal);
	if (t > 0) {
		v3 hit = v3_add(ray_o, v3_scale(ray_d, t));
		npv_drag_start_proj = v3_dot(v3_sub(hit, npv_drag_plane_point), axis_d);
	} else {
		npv_drag_start_proj = 0;
	}
}

static void npv_update_translate(v3 ray_o, v3 ray_d, NPV_Shape* s)
{
	float t = npv_ray_plane(ray_o, ray_d, npv_drag_plane_point, npv_drag_plane_normal);
	if (t <= 0) return;
	v3 hit = v3_add(ray_o, v3_scale(ray_d, t));
	v3 axis_d = npv_axis_dirs[npv_drag_axis];
	float proj = v3_dot(v3_sub(hit, npv_drag_plane_point), axis_d);
	float delta = proj - npv_drag_start_proj;
	s->pos = v3_add(npv_drag_start_pos, v3_scale(axis_d, delta));
}

static void npv_begin_rotate(v3 ray_o, v3 ray_d, int axis, NPV_Shape* s)
{
	npv_dragging = 1;
	npv_drag_axis = axis;
	npv_drag_start_rot = s->rot;
	v3 axis_d = npv_axis_dirs[axis];
	float t = npv_ray_plane(ray_o, ray_d, s->pos, axis_d);
	if (t > 0) {
		v3 hit = v3_add(ray_o, v3_scale(ray_d, t));
		v3 v = v3_sub(hit, s->pos);
		v3 pa, pb;
		npv_perp_axes(axis_d, &pa, &pb);
		npv_drag_start_angle = atan2f(v3_dot(v, pb), v3_dot(v, pa));
	} else {
		npv_drag_start_angle = 0;
	}
	npv_drag_plane_point = s->pos;
}

static void npv_update_rotate(v3 ray_o, v3 ray_d, NPV_Shape* s)
{
	v3 axis_d = npv_axis_dirs[npv_drag_axis];
	float t = npv_ray_plane(ray_o, ray_d, npv_drag_plane_point, axis_d);
	if (t <= 0) return;
	v3 hit = v3_add(ray_o, v3_scale(ray_d, t));
	v3 v = v3_sub(hit, npv_drag_plane_point);
	v3 pa, pb;
	npv_perp_axes(axis_d, &pa, &pb);
	float angle = atan2f(v3_dot(v, pb), v3_dot(v, pa));
	float delta = angle - npv_drag_start_angle;
	quat dq = npv_quat_axis_angle(axis_d, delta);
	s->rot = quat_mul(dq, npv_drag_start_rot);
}

// ----------------------------------------------------------------------------
// Gizmo rendering.

static void npv_draw_arrow(v3 from, v3 to, v3 color)
{
	render_debug_line(from, to, color);
	v3 dir = v3_norm(v3_sub(to, from));
	v3 pa, pb;
	npv_perp_axes(dir, &pa, &pb);
	float tip = 0.08f;
	v3 base = v3_add(to, v3_scale(dir, -tip * 3.0f));
	render_debug_line(to, v3_add(base, v3_scale(pa, tip)), color);
	render_debug_line(to, v3_add(base, v3_scale(pa, -tip)), color);
	render_debug_line(to, v3_add(base, v3_scale(pb, tip)), color);
	render_debug_line(to, v3_add(base, v3_scale(pb, -tip)), color);
}

static void npv_draw_circle(v3 center, v3 axis, float radius, v3 color)
{
	const int segs = 48;
	v3 pa, pb;
	npv_perp_axes(axis, &pa, &pb);
	v3 prev = v3_add(center, v3_scale(pa, radius));
	for (int i = 1; i <= segs; i++) {
		float a = 2.0f * 3.14159265f * (float)i / (float)segs;
		v3 pt = v3_add(center, v3_add(v3_scale(pa, cosf(a) * radius), v3_scale(pb, sinf(a) * radius)));
		render_debug_line(prev, pt, color);
		prev = pt;
	}
}

static void npv_draw_gizmo(NPV_Shape* s, int hover_axis)
{
	v3 c = s->pos;
	if (npv_gizmo_mode == 0) {
		// Translation: 3 axis arrows.
		for (int i = 0; i < 3; i++) {
			v3 tip = v3_add(c, v3_scale(npv_axis_dirs[i], npv_gizmo_size));
			v3 col = (i == hover_axis || (npv_dragging && i == npv_drag_axis)) ? npv_axis_colors_bright[i] : npv_axis_colors[i];
			npv_draw_arrow(c, tip, col);
		}
	} else {
		// Rotation: 3 circles.
		for (int i = 0; i < 3; i++) {
			v3 col = (i == hover_axis || (npv_dragging && i == npv_drag_axis)) ? npv_axis_colors_bright[i] : npv_axis_colors[i];
			npv_draw_circle(c, npv_axis_dirs[i], npv_gizmo_size, col);
		}
	}
}

// ----------------------------------------------------------------------------
// Manifold rendering.

static void npv_draw_manifold()
{
	if (!npv_colliding) return;
	for (int i = 0; i < npv_manifold.count; i++) {
		Contact* ct = &npv_manifold.contacts[i];
		v3 p = ct->point;
		v3 n = ct->normal;

		// Contact point: yellow cross.
		float s = 0.06f;
		v3 yellow = V3(1.0f, 1.0f, 0.0f);
		render_debug_line(V3(p.x-s, p.y, p.z), V3(p.x+s, p.y, p.z), yellow);
		render_debug_line(V3(p.x, p.y-s, p.z), V3(p.x, p.y+s, p.z), yellow);
		render_debug_line(V3(p.x, p.y, p.z-s), V3(p.x, p.y, p.z+s), yellow);

		// Normal arrow (green, 0.4 length).
		v3 tip = v3_add(p, v3_scale(n, 0.4f));
		render_debug_line(p, tip, V3(0.0f, 1.0f, 0.5f));

		// Penetration depth line (red, into shape B along normal).
		if (ct->penetration > 0) {
			v3 pen_end = v3_add(p, v3_scale(n, -ct->penetration));
			render_debug_line(p, pen_end, V3(1.0f, 0.3f, 0.3f));
		}
	}
}

// ----------------------------------------------------------------------------
// Initialization.

static void npv_init_shape(NPV_Shape* s, int idx)
{
	*s = (NPV_Shape){0};
	s->kind = (idx == 0) ? NPV_SPHERE : NPV_BOX;
	s->pos = (idx == 0) ? V3(-1.0f, 1.0f, 0.0f) : V3(1.0f, 1.0f, 0.0f);
	s->rot = quat_identity();
	s->radius = 0.5f;
	s->half_height = 0.5f;
	s->half_extents = V3(0.5f, 0.5f, 0.5f);
	s->rand_verts = 12;
	s->rand_seed = 12345 + (uint32_t)idx * 99991;
	s->color = (idx == 0) ? V3(0.85f, 0.35f, 0.3f) : V3(0.3f, 0.45f, 0.85f);
	// Allocate mesh slot.
	assert(r_mesh_count < MAX_MESH_TYPES);
	s->mesh_slot = r_mesh_count++;
	r_mesh_ready[s->mesh_slot] = 0;
	s->mesh_valid = 0;
}

static void npv_init()
{
	for (int i = 0; i < 2; i++) {
		npv_init_shape(&npv_shapes[i], i);
		if (npv_shapes[i].kind >= NPV_HULL_TETRA) npv_rebuild_hull(&npv_shapes[i]);
		npv_rebuild_mesh(&npv_shapes[i]);
	}
	npv_selected = 0;
	npv_gizmo_mode = 0;
	npv_drag_axis = -1;
	npv_dragging = 0;
	npv_colliding = 0;
	npv_manifold = (Manifold){0};
	npv_initialized = 1;
}

// ----------------------------------------------------------------------------
// Shape UI panel for one shape.

static int npv_shape_ui(NPV_Shape* s, const char* label, int shape_idx)
{
	int changed = 0;
	ImGui_PushIDInt(shape_idx);

	int old_kind = s->kind;
	if (ImGui_Combo("Type", &s->kind, s_npv_shape_names)) {
		if (s->kind != old_kind) {
			npv_free_hull(s);
			if (s->kind >= NPV_HULL_TETRA) npv_rebuild_hull(s);
			npv_rebuild_mesh(s);
			changed = 1;
		}
	}

	ImGui_DragFloat3Ex("Position", &s->pos.x, 0.05f, -50.0f, 50.0f, "%.2f", 0);

	// Rotation as euler degrees.
	v3 euler_deg = v3_scale(npv_quat_to_euler(s->rot), 180.0f / 3.14159265f);
	if (ImGui_DragFloat3Ex("Rotation", &euler_deg.x, 1.0f, -360.0f, 360.0f, "%.1f deg", 0)) {
		v3 euler_rad = v3_scale(euler_deg, 3.14159265f / 180.0f);
		s->rot = npv_euler_to_quat(euler_rad.x, euler_rad.y, euler_rad.z);
	}

	// Shape-specific params.
	switch (s->kind) {
	case NPV_SPHERE:
		if (ImGui_DragFloat3Ex("Radius", &s->radius, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		break;
	case NPV_CAPSULE:
		if (ImGui_DragFloat3Ex("Radius##cap", &s->radius, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		if (ImGui_DragFloat3Ex("Half Height", &s->half_height, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		break;
	case NPV_BOX:
		if (ImGui_DragFloat3Ex("Half Extents", &s->half_extents.x, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		break;
	case NPV_CYLINDER:
		if (ImGui_DragFloat3Ex("Radius##cyl", &s->radius, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		if (ImGui_DragFloat3Ex("Half Height##cyl", &s->half_height, 0.01f, 0.05f, 10.0f, "%.2f", 0)) changed = 1;
		break;
	default: // hull types
		if (s->kind == NPV_HULL_RANDOM) {
			if (ImGui_SliderInt("Vertices", &s->rand_verts, 4, 128)) {
				npv_rebuild_hull(s);
				npv_rebuild_mesh(s);
				changed = 1;
			}
			if (ImGui_Button("Randomize")) {
				s->rand_seed += 7919;
				npv_rebuild_hull(s);
				npv_rebuild_mesh(s);
				changed = 1;
			}
		}
		if (s->hull) ImGui_Text("Verts: %d  Faces: %d", s->hull->vert_count, s->hull->face_count);
		break;
	}

	if (changed && (s->kind == NPV_CAPSULE || s->kind == NPV_CYLINDER)) npv_rebuild_mesh(s);

	ImGui_PopID();
	return changed;
}

// ----------------------------------------------------------------------------
// Update: input + UI.

static void npv_update()
{
	if (!npv_initialized) npv_init();

	ImGuiIO* io = ImGui_GetIO();

	// Camera input (same as main app).
	if (!io->WantCaptureMouse) {
		// Scroll zoom always.
		if (io->MouseWheel != 0.0f) cam_zoom(io->MouseWheel);

		// Middle mouse pan.
		if (io->MouseDown[2]) cam_pan(io->MouseDelta.x, io->MouseDelta.y);

		// Left mouse: gizmo drag or orbit.
		if (io->MouseDown[0]) {
			v3 ray_o, ray_d;
			screen_to_ray(io->MousePos.x, io->MousePos.y, &ray_o, &ray_d);

			if (!npv_dragging && !io->MouseDownDuration[0]) {
				// Just clicked -- try gizmo pick.
				NPV_Shape* sel = &npv_shapes[npv_selected];
				int axis = -1;
				if (npv_gizmo_mode == 0) {
					axis = npv_pick_translate_axis(ray_o, ray_d, sel->pos, 0.2f);
				} else {
					axis = npv_pick_rotate_axis(ray_o, ray_d, sel->pos, 0.2f);
				}
				if (axis >= 0) {
					if (npv_gizmo_mode == 0) npv_begin_translate(ray_o, ray_d, axis, sel);
					else npv_begin_rotate(ray_o, ray_d, axis, sel);
				} else {
					// Try shape pick (select shape under cursor).
					for (int i = 0; i < 2; i++) {
						float bounding = fmaxf(fmaxf(npv_shapes[i].half_extents.x, npv_shapes[i].half_extents.y), npv_shapes[i].half_extents.z);
						if (bounding < npv_shapes[i].radius) bounding = npv_shapes[i].radius;
						if (bounding < npv_shapes[i].half_height + npv_shapes[i].radius) bounding = npv_shapes[i].half_height + npv_shapes[i].radius;
						if (bounding < 0.3f) bounding = 0.3f;
						float t = ray_sphere_simple(ray_o, ray_d, npv_shapes[i].pos, bounding);
						if (t >= 0) { npv_selected = i; break; }
					}
				}
			}

			if (npv_dragging) {
				NPV_Shape* sel = &npv_shapes[npv_selected];
				if (npv_gizmo_mode == 0) npv_update_translate(ray_o, ray_d, sel);
				else npv_update_rotate(ray_o, ray_d, sel);
			} else if (io->MouseDownDuration[0] > 0) {
				// Not dragging gizmo -- orbit camera.
				cam_orbit(io->MouseDelta.x, io->MouseDelta.y);
			}
		} else {
			if (npv_dragging) { npv_dragging = 0; npv_drag_axis = -1; }
		}
	} else {
		if (npv_dragging) { npv_dragging = 0; npv_drag_axis = -1; }
	}

	// Keyboard: W=translate, E=rotate.
	if (!io->WantCaptureKeyboard) {
		if (ImGui_IsKeyPressedEx(ImGuiKey_W, false)) npv_gizmo_mode = 0;
		if (ImGui_IsKeyPressedEx(ImGuiKey_E, false)) npv_gizmo_mode = 1;
	}

	// Run collision.
	npv_run_collide();

	// --- UI Panel ---
	ImGui_Begin("NP Viz", NULL, 0);

	if (ImGui_Button("Back to Simulation")) g_npv_mode = 0;

	// Gizmo mode.
	ImGui_SeparatorText("Gizmo");
	if (ImGui_RadioButton("Translate (W)", npv_gizmo_mode == 0)) npv_gizmo_mode = 0;
	ImGui_SameLine();
	if (ImGui_RadioButton("Rotate (E)", npv_gizmo_mode == 1)) npv_gizmo_mode = 1;

	// Shape selector tabs.
	ImGui_SeparatorText("Shapes");
	if (ImGui_RadioButton("Shape A", npv_selected == 0)) npv_selected = 0;
	ImGui_SameLine();
	if (ImGui_RadioButton("Shape B", npv_selected == 1)) npv_selected = 1;
	ImGui_SameLine();
	if (ImGui_Button("Reset")) {
		npv_free_hull(&npv_shapes[npv_selected]);
		npv_init_shape(&npv_shapes[npv_selected], npv_selected);
		// Reclaim the mesh slot (init_shape allocates a new one, but we want to reuse).
		// Undo the allocation: decrement mesh count and use old slot.
		r_mesh_count--;
		npv_shapes[npv_selected].mesh_slot = (npv_selected == 0) ? npv_shapes[0].mesh_slot : npv_shapes[1].mesh_slot;
		if (npv_shapes[npv_selected].kind >= NPV_HULL_TETRA) npv_rebuild_hull(&npv_shapes[npv_selected]);
		npv_rebuild_mesh(&npv_shapes[npv_selected]);
	}

	ImGui_Spacing();
	npv_shape_ui(&npv_shapes[npv_selected], npv_selected == 0 ? "A" : "B", npv_selected);

	// Manifold info.
	ImGui_SeparatorText("Manifold");
	if (npv_colliding) {
		ImGui_TextColored((ImVec4){0.3f, 1.0f, 0.5f, 1.0f}, "Colliding: %d contact%s", npv_manifold.count, npv_manifold.count != 1 ? "s" : "");
		for (int i = 0; i < npv_manifold.count; i++) {
			Contact* ct = &npv_manifold.contacts[i];
			ImGui_Text("  #%d pos=(%.3f, %.3f, %.3f)", i, ct->point.x, ct->point.y, ct->point.z);
			ImGui_Text("     n=(%.3f, %.3f, %.3f) pen=%.4f", ct->normal.x, ct->normal.y, ct->normal.z, ct->penetration);
			ImGui_Text("     feature=0x%08X", ct->feature_id);
		}
	} else {
		ImGui_TextColored((ImVec4){1.0f, 0.4f, 0.4f, 1.0f}, "Not colliding");
	}

	ImGui_End();
}

// ----------------------------------------------------------------------------
// Draw: render shapes + gizmo + manifold.

static void npv_draw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	render_draw_bg(V3(0.30f, 0.33f, 0.38f), V3(0.12f, 0.12f, 0.14f));

	float aspect = (float)g_width / (float)g_height;
	mat4 proj = mat4_perspective(1.0f, aspect, 0.1f, 500.0f);
	mat4 view = cam_view_matrix();
	mat4 vp = mul(proj, view);

	render_set_shadows(0);
	render_set_no_depth_write(1); // always translucent
	render_begin(vp);

	// Ground grid.
	{
		v3 grid_color = V3(0.25f, 0.25f, 0.25f);
		float extent = 5.0f;
		float step = 1.0f;
		for (float x = -extent; x <= extent; x += step) {
			render_debug_line(V3(x, 0, -extent), V3(x, 0, extent), grid_color);
			render_debug_line(V3(-extent, 0, x), V3(extent, 0, x), grid_color);
		}
		// Axis indicators on ground.
		render_debug_line(V3(0, 0, 0), V3(1, 0, 0), V3(0.5f, 0.15f, 0.15f));
		render_debug_line(V3(0, 0, 0), V3(0, 0, 1), V3(0.15f, 0.15f, 0.5f));
	}

	// Render shapes (translucent).
	for (int i = 0; i < 2; i++) {
		NPV_Shape* s = &npv_shapes[i];
		int mt = npv_mesh_type(s);
		v3 sc = npv_mesh_scale(s);
		float opacity = 0.35f;
		v3 color = (i == npv_selected) ? v3_add(v3_scale(s->color, 0.7f), V3(0.3f, 0.3f, 0.3f)) : s->color;
		render_push(mt, mat4_trs(s->pos, s->rot, sc), color, opacity);
	}

	// Draw wireframe outline for selected shape (double render at slightly larger scale).
	{
		NPV_Shape* sel = &npv_shapes[npv_selected];
		int mt = npv_mesh_type(sel);
		v3 sc = v3_add(npv_mesh_scale(sel), V3(0.005f, 0.005f, 0.005f));
		render_push(mt, mat4_trs(sel->pos, sel->rot, sc), V3(1, 1, 1), 0.08f);
	}

	// Draw gizmo for selected shape.
	{
		int hover_axis = -1;
		ImGuiIO* io = ImGui_GetIO();
		if (!io->WantCaptureMouse && !npv_dragging) {
			v3 ray_o, ray_d;
			screen_to_ray(io->MousePos.x, io->MousePos.y, &ray_o, &ray_d);
			NPV_Shape* sel = &npv_shapes[npv_selected];
			if (npv_gizmo_mode == 0) hover_axis = npv_pick_translate_axis(ray_o, ray_d, sel->pos, 0.2f);
			else hover_axis = npv_pick_rotate_axis(ray_o, ray_d, sel->pos, 0.2f);
		}
		npv_draw_gizmo(&npv_shapes[npv_selected], hover_axis);
	}

	// Draw manifold.
	npv_draw_manifold();

	render_end();
}
