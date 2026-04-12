// inertia.c -- mass and inertia tensor computation

// Multiply world-space inverse inertia tensor by a vector.
// I_world_inv * v = R * diag(inv_i) * R^T * v
static v3 inv_inertia_mul(quat rot, v3 inv_i, v3 v)
{
	v3 local = rotate(inv(rot), v);
	return rotate(rot, V3(local.x * inv_i.x, local.y * inv_i.y, local.z * inv_i.z));
}

// Precompute world-space inverse inertia as symmetric 3x3 matrix: I_w = R * diag(inv_i) * R^T.
// Reads inv_inertia_local and rotation from BodyState, writes iw_diag/iw_off to BodyHot.
static void body_compute_inv_inertia_world(BodyHot* h, BodyState* s)
{
	float a = s->inv_inertia_local.x, b = s->inv_inertia_local.y, c = s->inv_inertia_local.z;
	// Uniform inertia (cubes, spheres): I_w = a*I, rotation cancels out.
	if (a == b && b == c) { h->iw_diag = V3(a, a, a); h->iw_off = V3(0, 0, 0); return; }
	quat q = s->rotation;
	float xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
	float xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
	float wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;
	float r00 = 1-2*(yy+zz), r01 = 2*(xy-wz), r02 = 2*(xz+wy);
	float r10 = 2*(xy+wz), r11 = 1-2*(xx+zz), r12 = 2*(yz-wx);
	float r20 = 2*(xz-wy), r21 = 2*(yz+wx), r22 = 1-2*(xx+yy);
	h->iw_diag = V3(a*r00*r00 + b*r01*r01 + c*r02*r02, a*r10*r10 + b*r11*r11 + c*r12*r12, a*r20*r20 + b*r21*r21 + c*r22*r22);
	h->iw_off = V3(a*r00*r10 + b*r01*r11 + c*r02*r12, a*r00*r20 + b*r01*r21 + c*r02*r22, a*r10*r20 + b*r11*r21 + c*r12*r22);
}

// Multiply precomputed world-space inverse inertia by a vector.
static inline v3 inv_inertia_world_mul(BodyHot* h, v3 v)
{
	return V3(h->iw_diag.x*v.x + h->iw_off.x*v.y + h->iw_off.y*v.z, h->iw_off.x*v.x + h->iw_diag.y*v.y + h->iw_off.z*v.z, h->iw_off.y*v.x + h->iw_off.z*v.y + h->iw_diag.z*v.z);
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
	case SHAPE_CYLINDER: {
		// Solid cylinder along Y: I_y = m*r^2/2, I_x=I_z = m*(3r^2 + H^2)/12
		float r = s->cylinder.radius, hh = s->cylinder.half_height;
		float r2 = r*r;
		float H = 2.0f * hh;
		float iy = 0.5f * mass * r2;
		float ix = mass * (3.0f * r2 + H*H) / 12.0f;
		return V3(ix, iy, ix);
	}
	}
	return V3(0, 0, 0);
}

static v3 inertia_to_inv(v3 inertia) { return rcp(inertia); }

// Gyroscopic torque solver (single Newton-Raphson step).
// Corrects angular velocity for gyroscopic precession effects that explicit
// Euler integration misses. Without this, spinning bodies gain energy.
// Reference: Catto, GDC 2015 slide 76.
static v3 solve_gyroscopic(quat q, v3 inv_i, v3 omega, float h)
{
	v3 ib = rcp(inv_i);
	v3 wb = rotate(inv(q), omega);
	v3 iw = hmul(ib, wb);

	// Residual: f = h * cross(wb, Ib * wb)
	v3 f = scale(cross(wb, iw), h);

	// Jacobian: J = Ib + h * (skew(wb) * Ib - skew(Ib * wb))
	m3x3 Ib = diag(ib);
	m3x3 J = add(Ib, scale(sub(mul(skew(wb), Ib), skew(iw)), h));

	// Single Newton-Raphson update
	wb = sub(wb, solve(J, f));

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
	case SHAPE_CYLINDER: {
		float r = s->cylinder.radius, h = s->cylinder.half_height;
		return PI * r * r * (2.0f * h);
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
		w->body_state[idx].inv_inertia_local = V3(0, 0, 0);
		return;
	}

	ShapeInternal* shapes = w->body_cold[idx].shapes;
	int n = asize(shapes);

	float total_vol = 0.0f;
	for (int i = 0; i < n; i++) {
		total_vol += shape_volume(&shapes[i]);
	}
	if (total_vol <= 0.0f) {
		w->body_state[idx].inv_inertia_local = V3(0, 0, 0);
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

	w->body_state[idx].inv_inertia_local = inertia_to_inv(total);
}
