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
		// Volume-weighted mass split: cylindrical body vs sphere (two hemispheres)
		float v_cyl = 2.0f * hh * r2; // pi cancels in ratio
		float v_sph = (4.0f/3.0f) * r2 * r;
		float v_tot = v_cyl + v_sph;
		float mc = mass * v_cyl / v_tot;
		float ms = mass * v_sph / v_tot;
		// Axial (Y): cylindrical body + sphere
		float iy = mc * r2 / 2.0f + ms * 2.0f * r2 / 5.0f;
		// Transverse (X,Z): cylindrical body + two hemispheres via parallel axis
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
	}
	return 0.0f;
}

// Recompute body inertia from all shapes. Mass is distributed by volume ratio,
// each shape's inertia is rotated by its local_rot into the body frame, then
// shifted to the body origin via parallel axis theorem. Off-diagonal moments
// are dropped (BodyState.inv_inertia_local is v3 diagonal) -- exact for 90/180
// deg aligned children, an approximation for arbitrary rotations.
static void recompute_body_inertia(WorldInternal* w, int idx)
{
	float mass = w->body_cold[idx].mass;
	if (mass <= 0.0f) {
		body_inv_inertia_local(w, idx) = V3(0, 0, 0);
		return;
	}

	ShapeInternal* shapes = w->body_cold[idx].shapes;
	int n = asize(shapes);

	float total_vol = 0.0f;
	for (int i = 0; i < n; i++) {
		total_vol += shape_volume(&shapes[i]);
	}
	if (total_vol <= 0.0f) {
		body_inv_inertia_local(w, idx) = V3(0, 0, 0);
		return;
	}

	v3 total = V3(0, 0, 0);
	for (int i = 0; i < n; i++) {
		float sm = mass * shape_volume(&shapes[i]) / total_vol;
		v3 li = shape_inertia(&shapes[i], sm); // diagonal in child local frame

		// Rotate child inertia into body frame: I_b = R * diag(li) * R^T.
		// Only the diagonal of I_b is kept. For a rotation with basis vectors
		// (ex, ey, ez), the i-th body-frame diagonal entry is
		//   (I_b)_ii = li.x*ex_i^2 + li.y*ey_i^2 + li.z*ez_i^2
		// which reduces to li when local_rot = identity.
		v3 ex = rotate(shapes[i].local_rot, V3(1, 0, 0));
		v3 ey = rotate(shapes[i].local_rot, V3(0, 1, 0));
		v3 ez = rotate(shapes[i].local_rot, V3(0, 0, 1));
		v3 body_diag = V3(li.x*ex.x*ex.x + li.y*ey.x*ey.x + li.z*ez.x*ez.x,
		                  li.x*ex.y*ex.y + li.y*ey.y*ey.y + li.z*ez.y*ez.y,
		                  li.x*ex.z*ex.z + li.y*ey.z*ey.z + li.z*ez.z*ez.z);

		// Parallel axis theorem (diagonal of m*(|d|^2*I - d*d^T)).
		v3 d = shapes[i].local_pos;
		total.x += body_diag.x + sm * (d.y*d.y + d.z*d.z);
		total.y += body_diag.y + sm * (d.x*d.x + d.z*d.z);
		total.z += body_diag.z + sm * (d.x*d.x + d.y*d.y);
	}

	body_inv_inertia_local(w, idx) = inertia_to_inv(total);
}

// Fill r_min / r_max on a ShapeInternal (inscribed and bounding sphere radii
// about the shape's local origin). Includes rounding_radius on box/hull.
// Mesh and heightfield are static and skipped (r_min = r_max = 0).
static void shape_compute_bounds(ShapeInternal* s)
{
	float rr = s->rounding_radius;
	switch (s->type) {
	case SHAPE_SPHERE:
		s->r_min = s->sphere.radius;
		s->r_max = s->sphere.radius;
		break;
	case SHAPE_CAPSULE:
		s->r_min = s->capsule.radius;
		s->r_max = s->capsule.half_height + s->capsule.radius;
		break;
	case SHAPE_BOX: {
		v3 h = s->box.half_extents;
		float hmin = h.x < h.y ? h.x : h.y;
		if (h.z < hmin) hmin = h.z;
		s->r_min = hmin + rr;
		s->r_max = sqrtf(h.x*h.x + h.y*h.y + h.z*h.z) + rr;
		break;
	}
	case SHAPE_HULL: {
		const Hull* hull = s->hull.hull;
		v3 sc = s->hull.scale;
		// r_max: scan vertices, measure distance from centroid.
		// Centroid lives in the hull's own local frame; we scale both the
		// centroid reference point and the vertex so the measurement is
		// taken in scaled-shape space.
		v3 c = V3(hull->centroid.x * sc.x, hull->centroid.y * sc.y, hull->centroid.z * sc.z);
		float rmax2 = 0.0f;
		for (int i = 0; i < hull->vert_count; i++) {
			v3 v = V3(hull->verts[i].x * sc.x, hull->verts[i].y * sc.y, hull->verts[i].z * sc.z);
			v3 d = V3(v.x - c.x, v.y - c.y, v.z - c.z);
			float m2 = d.x*d.x + d.y*d.y + d.z*d.z;
			if (m2 > rmax2) rmax2 = m2;
		}
		// r_min: min distance from centroid to any face plane. Planes are
		// in unscaled hull space; with non-uniform scale we conservatively
		// use the minimum axis scale as a lower bound on inscribed radius.
		// Uniform scale case is exact.
		float smin = sc.x < sc.y ? sc.x : sc.y; if (sc.z < smin) smin = sc.z;
		float rmin = 1e18f;
		for (int i = 0; i < hull->face_count; i++) {
			// plane: dot(n, p) = offset, for any p on the face. Distance from
			// centroid = dot(n, centroid) - offset (signed; inside is negative
			// for outward-pointing normals). We want the absolute.
			float d = hull->planes[i].normal.x * hull->centroid.x
			        + hull->planes[i].normal.y * hull->centroid.y
			        + hull->planes[i].normal.z * hull->centroid.z
			        - hull->planes[i].offset;
			float ad = d < 0.0f ? -d : d;
			if (ad < rmin) rmin = ad;
		}
		if (rmin >= 1e18f) rmin = 0.0f;
		s->r_min = rmin * smin + rr;
		s->r_max = sqrtf(rmax2) + rr;
		break;
	}
	case SHAPE_MESH:
	case SHAPE_HEIGHTFIELD:
		s->r_min = 0.0f;
		s->r_max = 0.0f;
		break;
	}
}

// CCD speed test (Catto GDC 2025). A body is "slow" this step iff the maximum
// possible displacement of any surface point within dt fits inside the body's
// inscribed sphere:
//   (|v_lin| + r_max * |omega|) * dt < r_min
// Slow bodies use speculative contacts only. Fast bodies get enrolled for TOI
// advancement against static world geometry.
static int body_is_slow_this_step(BodyHot* h, BodyCold* c, float dt)
{
	if (h->inv_mass == 0.0f) return 1; // static/kinematic never TOI'd
	if (c->r_max_body <= 0.0f) return 1; // no dynamic shapes (mesh/heightfield bodies)
	// Zero r_min is VALID for compound bodies whose COM lies outside any
	// child shape's inscribed sphere. Such bodies always fail the speed
	// test and always enroll for TOI — which is the correct fallback
	// because no speculative guarantee can hold when the body has no
	// interior sphere around the COM.
	float v = sqrtf(h->velocity.x*h->velocity.x + h->velocity.y*h->velocity.y + h->velocity.z*h->velocity.z);
	float w = sqrtf(h->angular_velocity.x*h->angular_velocity.x + h->angular_velocity.y*h->angular_velocity.y + h->angular_velocity.z*h->angular_velocity.z);
	float disp = (v + c->r_max_body * w) * dt;
	return disp < c->r_min_body;
}

// Rebuild WorldInternal.fast_body_list for the upcoming step. Iterates every
// live dynamic body and enrolls those failing the speed test. Cleared and
// repopulated each step. The list is consumed by the post-solve TOI pass.
static void classify_fast_bodies(WorldInternal* w, float dt)
{
	if (w->fast_body_list)    asetlen(w->fast_body_list, 0);
	if (w->fast_body_pre_pos) asetlen(w->fast_body_pre_pos, 0);
	if (w->fast_body_pre_rot) asetlen(w->fast_body_pre_rot, 0);
	int n = asize(w->body_cold);
	for (int i = 0; i < n; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		if (body_is_slow_this_step(&w->body_hot[i], &w->body_cold[i], dt)) continue;
		apush(w->fast_body_list,    i);
		apush(w->fast_body_pre_pos, w->body_state[i].position);
		apush(w->fast_body_pre_rot, w->body_state[i].rotation);
	}
}

// Recompute body-level r_min_body / r_max_body from current shapes.
// Shapes with r_max == 0 (mesh/heightfield) are skipped. Per Catto's speed
// test (GDC 2025), for compound bodies we use the conservative aggregate:
//   r_max_body = max over shapes of (|local_pos| + r_max)
//   r_min_body = min over shapes of max(0, r_min - |local_pos|)
// Sets both to zero when the body has no CCD-eligible shapes.
static void recompute_body_bounds(WorldInternal* w, int idx)
{
	ShapeInternal* shapes = w->body_cold[idx].shapes;
	int n = asize(shapes);
	float r_max_body = 0.0f;
	float r_min_body = 1e18f;
	int any = 0;
	for (int i = 0; i < n; i++) {
		ShapeInternal* s = &shapes[i];
		if (s->r_max <= 0.0f) continue;
		float off = sqrtf(s->local_pos.x*s->local_pos.x + s->local_pos.y*s->local_pos.y + s->local_pos.z*s->local_pos.z);
		float rmax_i = off + s->r_max;
		float rmin_i = s->r_min - off;
		if (rmin_i < 0.0f) rmin_i = 0.0f;
		if (rmax_i > r_max_body) r_max_body = rmax_i;
		if (rmin_i < r_min_body) r_min_body = rmin_i;
		any = 1;
	}
	if (!any) { r_min_body = 0.0f; r_max_body = 0.0f; }
	w->body_cold[idx].r_min_body = r_min_body;
	w->body_cold[idx].r_max_body = r_max_body;
}
