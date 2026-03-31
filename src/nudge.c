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
	v3 inv_inertia_local; // diagonal of local-space inverse inertia tensor
} BodyHot;

typedef struct WorldInternal
{
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts; // flat contact array from last step
} WorldInternal;

#include "gjk.c"
#include "quickhull.c"
#include "collision.c"

// -----------------------------------------------------------------------------
// Inertia tensor helpers.

// Multiply world-space inverse inertia tensor by a vector.
// I_world_inv * v = R * diag(inv_i) * R^T * v
static v3 inv_inertia_mul(quat rot, v3 inv_i, v3 v)
{
	v3 local = rotate(inv(rot), v);
	return rotate(rot, V3(local.x * inv_i.x, local.y * inv_i.y, local.z * inv_i.z));
}

// Compute diagonal inertia tensor for a shape (in local principal axes).
static v3 shape_inertia(ShapeInternal* s, float mass)
{
	if (mass <= 0.0f) return V3(0, 0, 0);

	switch (s->type) {
	case SHAPE_SPHERE: {
		float i = 0.4f * mass * s->sphere.radius * s->sphere.radius;
		return V3(i, i, i);
	}
	case SHAPE_BOX: {
		v3 h = s->box.half_extents;
		float x2 = 4.0f*h.x*h.x, y2 = 4.0f*h.y*h.y, z2 = 4.0f*h.z*h.z;
		return V3(mass*(y2+z2)/12.0f, mass*(x2+z2)/12.0f, mass*(x2+y2)/12.0f);
	}
	case SHAPE_CAPSULE: {
		float r = s->capsule.radius, hh = s->capsule.half_height;
		float r2 = r*r, hh2 = hh*hh;
		// Volume-weighted mass split: cylinder vs sphere (two hemispheres)
		float v_cyl = 2.0f * hh * r2; // pi cancels in ratio
		float v_sph = (4.0f/3.0f) * r2 * r;
		float v_tot = v_cyl + v_sph;
		float mc = mass * v_cyl / v_tot;
		float ms = mass * v_sph / v_tot;
		// Axial (Y): cylinder + sphere
		float iy = mc * r2 / 2.0f + ms * 2.0f * r2 / 5.0f;
		// Transverse (X,Z): cylinder + two hemispheres via parallel axis
		float ix_cyl = mc * (r2/4.0f + hh2/3.0f);
		float d = hh + 3.0f*r/8.0f; // hemisphere CoM offset from capsule center
		float ix_hemi = (83.0f/320.0f) * (ms/2.0f) * r2 + (ms/2.0f) * d * d;
		float ix = ix_cyl + 2.0f * ix_hemi;
		return V3(ix, iy, ix);
	}
	case SHAPE_HULL: {
		// Approximate via scaled AABB
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		float lo_x = 1e18f, hi_x = -1e18f;
		float lo_y = 1e18f, hi_y = -1e18f;
		float lo_z = 1e18f, hi_z = -1e18f;
		for (int i = 0; i < hull->vert_count; i++) {
			float x = hull->verts[i].x*sc.x, y = hull->verts[i].y*sc.y, z = hull->verts[i].z*sc.z;
			if (x < lo_x) lo_x = x; if (x > hi_x) hi_x = x;
			if (y < lo_y) lo_y = y; if (y > hi_y) hi_y = y;
			if (z < lo_z) lo_z = z; if (z > hi_z) hi_z = z;
		}
		float sx = hi_x-lo_x, sy = hi_y-lo_y, sz = hi_z-lo_z;
		return V3(mass*(sy*sy+sz*sz)/12.0f, mass*(sx*sx+sz*sz)/12.0f, mass*(sx*sx+sy*sy)/12.0f);
	}
	}
	return V3(0, 0, 0);
}

static v3 inertia_to_inv(v3 inertia)
{
	return V3(
		inertia.x > 0.0f ? 1.0f / inertia.x : 0.0f,
		inertia.y > 0.0f ? 1.0f / inertia.y : 0.0f,
		inertia.z > 0.0f ? 1.0f / inertia.z : 0.0f);
}

// Gyroscopic torque solver (single Newton-Raphson step).
// Corrects angular velocity for gyroscopic precession effects that explicit
// Euler integration misses. Without this, spinning bodies gain energy.
// Reference: Catto, GDC 2015 slide 76.
static v3 solve_gyroscopic(quat q, v3 inv_i, v3 omega, float h)
{
	// Body-space inertia diagonal (invert the inverse).
	v3 ib = V3(
		inv_i.x > 0.0f ? 1.0f / inv_i.x : 0.0f,
		inv_i.y > 0.0f ? 1.0f / inv_i.y : 0.0f,
		inv_i.z > 0.0f ? 1.0f / inv_i.z : 0.0f);

	// Convert to body coordinates.
	v3 wb = rotate(inv(q), omega);

	// Ib * wb (diagonal multiply).
	v3 iw = V3(ib.x * wb.x, ib.y * wb.y, ib.z * wb.z);

	// Residual: f = h * cross(wb, Ib * wb).
	v3 f = scale(cross(wb, iw), h);

	// Jacobian: J = Ib + h * (skew(wb) * Ib - skew(Ib * wb)).
	// skew(a) * diag(d) column k = cross(a, e_k * d_k).
	// skew(a) has row i = [-a cross with basis vectors].
	// J_ij = Ib_ij + h * (skew(wb) * Ib - skew(iw))_ij
	// Since Ib is diagonal, skew(wb)*Ib column j = cross(wb, e_j) * ib_j.
	// skew(iw) is just the skew-symmetric matrix of iw.
	float j[9]; // row-major 3x3
	// skew(wb) * Ib:  column j = cross(wb, e_j * ib_j)
	//   col0 = cross(wb, (ib.x, 0, 0)) = (0, wb.z*ib.x, -wb.y*ib.x)
	//   col1 = cross(wb, (0, ib.y, 0)) = (-wb.z*ib.y, 0, wb.x*ib.y)
	//   col2 = cross(wb, (0, 0, ib.z)) = (wb.y*ib.z, -wb.x*ib.z, 0)
	// skew(iw):
	//   row0 = (0, -iw.z, iw.y)
	//   row1 = (iw.z, 0, -iw.x)
	//   row2 = (-iw.y, iw.x, 0)
	// D = skew(wb)*Ib - skew(iw):
	float d00 = 0            - 0;
	float d01 = -wb.z * ib.y - (-iw.z);
	float d02 = wb.y * ib.z  - iw.y;
	float d10 = wb.z * ib.x  - iw.z;
	float d11 = 0            - 0;
	float d12 = -wb.x * ib.z - (-iw.x);
	float d20 = -wb.y * ib.x - (-iw.y);
	float d21 = wb.x * ib.y  - iw.x;
	float d22 = 0            - 0;

	// J = Ib + h * D  (Ib is diagonal, add to diagonal of D*h).
	j[0] = ib.x + h * d00;  j[1] = h * d01;         j[2] = h * d02;
	j[3] = h * d10;         j[4] = ib.y + h * d11;  j[5] = h * d12;
	j[6] = h * d20;         j[7] = h * d21;         j[8] = ib.z + h * d22;

	// Solve J * x = f via Cramer's rule.
	float det =
		j[0] * (j[4]*j[8] - j[5]*j[7]) -
		j[1] * (j[3]*j[8] - j[5]*j[6]) +
		j[2] * (j[3]*j[7] - j[4]*j[6]);
	if (fabsf(det) < 1e-12f) return omega;
	float inv_det = 1.0f / det;

	v3 x;
	x.x = ((j[4]*j[8]-j[5]*j[7])*f.x - (j[1]*j[8]-j[2]*j[7])*f.y + (j[1]*j[5]-j[2]*j[4])*f.z) * inv_det;
	x.y = (-(j[3]*j[8]-j[5]*j[6])*f.x + (j[0]*j[8]-j[2]*j[6])*f.y - (j[0]*j[5]-j[2]*j[3])*f.z) * inv_det;
	x.z = ((j[3]*j[7]-j[4]*j[6])*f.x - (j[0]*j[7]-j[1]*j[6])*f.y + (j[0]*j[4]-j[1]*j[3])*f.z) * inv_det;

	wb = sub(wb, x);

	// Back to world coordinates.
	return rotate(q, wb);
}

// Volume of a shape (for mass distribution across compound bodies).
static float shape_volume(ShapeInternal* s)
{
	const float PI = 3.14159265f;
	switch (s->type) {
	case SHAPE_SPHERE: {
		float r = s->sphere.radius;
		return (4.0f/3.0f) * PI * r * r * r;
	}
	case SHAPE_CAPSULE: {
		float r = s->capsule.radius, h = s->capsule.half_height;
		return PI * r * r * (2.0f * h + (4.0f/3.0f) * r);
	}
	case SHAPE_BOX: {
		v3 e = s->box.half_extents;
		return 8.0f * e.x * e.y * e.z;
	}
	case SHAPE_HULL: {
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		float lo_x = 1e18f, hi_x = -1e18f;
		float lo_y = 1e18f, hi_y = -1e18f;
		float lo_z = 1e18f, hi_z = -1e18f;
		for (int i = 0; i < hull->vert_count; i++) {
			float x = hull->verts[i].x*sc.x, y = hull->verts[i].y*sc.y, z = hull->verts[i].z*sc.z;
			if (x < lo_x) lo_x = x; if (x > hi_x) hi_x = x;
			if (y < lo_y) lo_y = y; if (y > hi_y) hi_y = y;
			if (z < lo_z) lo_z = z; if (z > hi_z) hi_z = z;
		}
		return (hi_x-lo_x) * (hi_y-lo_y) * (hi_z-lo_z);
	}
	}
	return 0.0f;
}

// Recompute body inertia from all shapes. Mass is distributed by volume ratio,
// each shape's inertia is shifted to body origin via parallel axis theorem.
static void recompute_body_inertia(WorldInternal* w, int idx)
{
	float mass = w->body_cold[idx].mass;
	if (mass <= 0.0f) {
		w->body_hot[idx].inv_inertia_local = V3(0, 0, 0);
		return;
	}

	ShapeInternal* shapes = w->body_cold[idx].shapes;
	int n = asize(shapes);

	float total_vol = 0.0f;
	for (int i = 0; i < n; i++)
		total_vol += shape_volume(&shapes[i]);
	if (total_vol <= 0.0f) {
		w->body_hot[idx].inv_inertia_local = V3(0, 0, 0);
		return;
	}

	v3 total = V3(0, 0, 0);
	for (int i = 0; i < n; i++) {
		float sm = mass * shape_volume(&shapes[i]) / total_vol;
		v3 li = shape_inertia(&shapes[i], sm);

		// Parallel axis theorem: I_body += I_local + m*(|d|^2*I - d*d^T)
		// For diagonal tensor: I_x += m*(dy^2+dz^2), etc.
		v3 d = shapes[i].local_pos;
		total.x += li.x + sm * (d.y*d.y + d.z*d.z);
		total.y += li.y + sm * (d.x*d.x + d.z*d.z);
		total.z += li.z + sm * (d.x*d.x + d.y*d.y);
	}

	w->body_hot[idx].inv_inertia_local = inertia_to_inv(total);
}

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
	afree(w->debug_contacts);
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

	// Store contacts for debug visualization
	aclear(w->debug_contacts);
	for (int i = 0; i < asize(manifolds); i++)
		for (int c = 0; c < manifolds[i].m.count; c++)
			apush(w->debug_contacts, manifolds[i].m.contacts[c]);

	// Resolve contacts (impulse with angular response)
	for (int i = 0; i < asize(manifolds); i++) {
		InternalManifold* im = &manifolds[i];
		BodyHot* a = &w->body_hot[im->body_a];
		BodyHot* b = &w->body_hot[im->body_b];
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		if (inv_mass_sum == 0.0f) continue;

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			v3 ra = sub(ct->point, a->position);
			v3 rb = sub(ct->point, b->position);

			// Velocity at contact point (linear + angular)
			v3 vel_a = add(a->velocity, cross(a->angular_velocity, ra));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, rb));
			float vn = dot(sub(vel_b, vel_a), ct->normal);
			if (vn > 0.0f) continue;

			// Effective inverse mass including angular contribution
			v3 ra_x_n = cross(ra, ct->normal);
			v3 rb_x_n = cross(rb, ct->normal);
			v3 ang_a = inv_inertia_mul(a->rotation, a->inv_inertia_local, ra_x_n);
			v3 ang_b = inv_inertia_mul(b->rotation, b->inv_inertia_local, rb_x_n);
			float eff_mass_inv = inv_mass_sum
				+ dot(cross(ang_a, ra), ct->normal)
				+ dot(cross(ang_b, rb), ct->normal);

			float j = -vn / eff_mass_inv;
			v3 impulse = scale(ct->normal, j);

			// Linear impulse
			a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
			b->velocity = add(b->velocity, scale(impulse, b->inv_mass));

			// Angular impulse
			a->angular_velocity = sub(a->angular_velocity,
				inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(ra, impulse)));
			b->angular_velocity = add(b->angular_velocity,
				inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(rb, impulse)));

			// Positional correction (Baumgarte, linear only)
			float slop = 0.01f;
			float percent = 0.4f;
			float corr = fmaxf(ct->penetration - slop, 0.0f) * percent / inv_mass_sum;
			v3 correction = scale(ct->normal, corr);
			a->position = sub(a->position, scale(correction, a->inv_mass));
			b->position = add(b->position, scale(correction, b->inv_mass));
		}
	}
	afree(manifolds);

	// Integrate positions and rotations
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		h->position = add(h->position, scale(h->velocity, dt));

		// Gyroscopic torque correction before angular integration.
		h->angular_velocity = solve_gyroscopic(
			h->rotation, h->inv_inertia_local, h->angular_velocity, dt);

		// Quaternion integration: q += 0.5 * dt * (0,omega) * q
		v3 w = h->angular_velocity;
		quat spin = { w.x, w.y, w.z, 0.0f };
		quat dq = mul(spin, h->rotation);
		h->rotation.x += 0.5f * dt * dq.x;
		h->rotation.y += 0.5f * dt * dq.y;
		h->rotation.z += 0.5f * dt * dq.z;
		h->rotation.w += 0.5f * dt * dq.w;
		// Renormalize
		float ql = sqrtf(h->rotation.x*h->rotation.x + h->rotation.y*h->rotation.y
			+ h->rotation.z*h->rotation.z + h->rotation.w*h->rotation.w);
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql;
		h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
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
	recompute_body_inertia(w, idx);
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

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}
