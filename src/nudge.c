// See LICENSE for licensing info.
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
	float friction;
	float restitution;
} BodyHot;

typedef struct WarmManifold WarmManifold; // forward decl for warm cache

typedef struct WorldInternal
{
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts; // flat contact array from last step
	CK_MAP(WarmManifold) warm_cache;     // contact impulse cache for warm starting
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
// Sequential impulses constraint solver.

#define SOLVER_VELOCITY_ITERS      10
#define SOLVER_POSITION_ITERS      4
#define SOLVER_BAUMGARTE           0.2f
#define SOLVER_SLOP                0.01f
#define SOLVER_RESTITUTION_THRESH  1.0f
#define SOLVER_POS_BAUMGARTE       0.2f
#define SOLVER_POS_MAX_CORRECTION  0.2f  // max position correction per step

typedef struct SolverContact
{
	v3 r_a;              // contact offset from body A
	v3 r_b;              // contact offset from body B
	v3 normal;           // cached contact normal
	v3 tangent1;         // friction direction 1
	v3 tangent2;         // friction direction 2
	float eff_mass_n;    // 1/K for normal row
	float eff_mass_t1;   // 1/K for tangent1 row
	float eff_mass_t2;   // 1/K for tangent2 row
	float bias;          // velocity bias for penetration recovery
	float bounce;        // restitution velocity bias
	float lambda_n;      // accumulated normal impulse (>= 0)
	float lambda_t1;     // accumulated tangent1 impulse
	float lambda_t2;     // accumulated tangent2 impulse
	float penetration;   // cached for position correction pass
} SolverContact;

typedef struct SolverManifold
{
	int body_a, body_b;
	int contact_start;
	int contact_count;
	float friction;
} SolverManifold;

// Warm starting: cached impulses from previous frame, keyed by body pair.
typedef struct WarmContact
{
	v3 point;            // world-space position for matching
	float lambda_n;
	float lambda_t1;
	float lambda_t2;
} WarmContact;

struct WarmManifold
{
	WarmContact contacts[MAX_CONTACTS];
	int count;
};

static uint64_t body_pair_key(int a, int b)
{
	uint32_t lo = a < b ? a : b;
	uint32_t hi = a < b ? b : a;
	return ((uint64_t)lo << 32) | (uint64_t)hi;
}

static void contact_tangent_basis(v3 n, v3* t1, v3* t2)
{
	if (fabsf(n.x) >= 0.57735f)
		*t1 = norm(V3(n.y, -n.x, 0.0f));
	else
		*t1 = norm(V3(0.0f, n.z, -n.y));
	*t2 = cross(n, *t1);
}

static float compute_effective_mass(
	BodyHot* a, BodyHot* b, float inv_mass_sum,
	v3 r_a, v3 r_b, v3 dir)
{
	v3 ra_x_d = cross(r_a, dir);
	v3 rb_x_d = cross(r_b, dir);
	float k = inv_mass_sum
		+ dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, ra_x_d), r_a), dir)
		+ dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, rb_x_d), r_b), dir);
	return k > 1e-12f ? 1.0f / k : 0.0f;
}

// Apply an impulse (linear + angular) to a body pair.
static void apply_impulse(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, v3 impulse)
{
	a->velocity = sub(a->velocity, scale(impulse, a->inv_mass));
	b->velocity = add(b->velocity, scale(impulse, b->inv_mass));
	a->angular_velocity = sub(a->angular_velocity,
		inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity,
		inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(r_b, impulse)));
}

// Match a new contact to the nearest cached contact. Returns index or -1.
static int warm_match(WarmManifold* wm, v3 point, float threshold2)
{
	int best = -1;
	float best_d2 = threshold2;
	for (int i = 0; i < wm->count; i++) {
		float d2 = len2(sub(point, wm->contacts[i].point));
		if (d2 < best_d2) { best_d2 = d2; best = i; }
	}
	return best;
}

static void solver_pre_solve(
	WorldInternal* w,
	InternalManifold* manifolds, int manifold_count,
	SolverManifold** out_sm, SolverContact** out_sc,
	float dt)
{
	CK_DYNA SolverManifold* sm = NULL;
	CK_DYNA SolverContact*  sc = NULL;
	float inv_dt = dt > 0.0f ? 1.0f / dt : 0.0f;

	for (int i = 0; i < manifold_count; i++) {
		InternalManifold* im = &manifolds[i];
		BodyHot* a = &w->body_hot[im->body_a];
		BodyHot* b = &w->body_hot[im->body_b];
		float inv_mass_sum = a->inv_mass + b->inv_mass;
		if (inv_mass_sum == 0.0f) continue;

		float mu = sqrtf(a->friction * b->friction);
		float rest = a->restitution > b->restitution ? a->restitution : b->restitution;

		// Look up warm starting data
		uint64_t key = body_pair_key(im->body_a, im->body_b);
		WarmManifold* wm = map_get_ptr(w->warm_cache, key);

		SolverManifold smf = {
			.body_a = im->body_a,
			.body_b = im->body_b,
			.contact_start = asize(sc),
			.contact_count = im->m.count,
			.friction = mu,
		};

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			SolverContact s = {0};

			s.r_a = sub(ct->point, a->position);
			s.r_b = sub(ct->point, b->position);
			s.normal = ct->normal;
			s.penetration = ct->penetration;
			contact_tangent_basis(ct->normal, &s.tangent1, &s.tangent2);

			s.eff_mass_n  = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.normal);
			s.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent1);
			s.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent2);

			s.bias = -SOLVER_BAUMGARTE * inv_dt * fmaxf(ct->penetration - SOLVER_SLOP, 0.0f);

			v3 vel_a = add(a->velocity, cross(a->angular_velocity, s.r_a));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, s.r_b));
			float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
			s.bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;

			// Warm start: match to cached contact by proximity
			if (wm) {
				int match = warm_match(wm, ct->point, 0.05f * 0.05f);
				if (match >= 0) {
					s.lambda_n  = wm->contacts[match].lambda_n;
					s.lambda_t1 = wm->contacts[match].lambda_t1;
					s.lambda_t2 = wm->contacts[match].lambda_t2;
				}
			}

			apush(sc, s);
		}

		apush(sm, smf);
	}

	// Apply warm start impulses
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (s->lambda_n == 0.0f && s->lambda_t1 == 0.0f && s->lambda_t2 == 0.0f)
				continue;
			v3 P = add(add(
				scale(s->normal, s->lambda_n),
				scale(s->tangent1, s->lambda_t1)),
				scale(s->tangent2, s->lambda_t2));
			apply_impulse(a, b, s->r_a, s->r_b, P);
		}
	}

	*out_sm = sm;
	*out_sc = sc;
}

static void solver_iterate(
	WorldInternal* w,
	SolverManifold* sm, int sm_count,
	SolverContact* sc)
{
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];

			// Relative velocity at contact point
			v3 dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));

			// --- Normal constraint ---
			float vn = dot(dv, s->normal);
			float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce));
			float old_n = s->lambda_n;
			s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
			float delta_n = s->lambda_n - old_n;
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, delta_n));

			// --- Friction tangent 1 ---
			dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt1 = dot(dv, s->tangent1);
			float lambda_t1 = s->eff_mass_t1 * (-vt1);
			float max_f = m->friction * s->lambda_n;
			float old_t1 = s->lambda_t1;
			s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + lambda_t1, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

			// --- Friction tangent 2 ---
			dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt2 = dot(dv, s->tangent2);
			float lambda_t2 = s->eff_mass_t2 * (-vt2);
			float old_t2 = s->lambda_t2;
			s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + lambda_t2, max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
		}
	}
}

// NGS position correction: directly fix remaining penetration after velocity solve
// and position integration. Operates on positions only, no velocity modification.
static void solver_position_correct(
	WorldInternal* w,
	SolverManifold* sm, int sm_count,
	SolverContact* sc)
{
	for (int iter = 0; iter < SOLVER_POSITION_ITERS; iter++) {
		for (int i = 0; i < sm_count; i++) {
			SolverManifold* m = &sm[i];
			BodyHot* a = &w->body_hot[m->body_a];
			BodyHot* b = &w->body_hot[m->body_b];
			float inv_mass_sum = a->inv_mass + b->inv_mass;

			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];

				// Recompute separation from current positions
				v3 r_a = sub(add(a->position, rotate(a->rotation,
					rotate(inv(a->rotation), s->r_a))), a->position);
				v3 r_b = sub(add(b->position, rotate(b->rotation,
					rotate(inv(b->rotation), s->r_b))), b->position);
				v3 p_a = add(a->position, r_a);
				v3 p_b = add(b->position, r_b);
				float separation = dot(sub(p_b, p_a), s->normal) - s->penetration;

				float C = fminf(0.0f, separation + SOLVER_SLOP);
				if (C >= 0.0f) continue;

				float correction = -SOLVER_POS_BAUMGARTE * C;
				if (correction > SOLVER_POS_MAX_CORRECTION)
					correction = SOLVER_POS_MAX_CORRECTION;

				// Effective mass for position correction (linear only for speed)
				float k = inv_mass_sum;
				float delta = correction / k;

				v3 P = scale(s->normal, delta);
				a->position = sub(a->position, scale(P, a->inv_mass));
				b->position = add(b->position, scale(P, b->inv_mass));
			}
		}
	}
}

// Store accumulated impulses into warm cache for next frame.
static void solver_post_solve(
	WorldInternal* w,
	SolverManifold* sm, int sm_count,
	SolverContact* sc,
	InternalManifold* manifolds, int manifold_count)
{
	map_clear(w->warm_cache);

	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		uint64_t key = body_pair_key(m->body_a, m->body_b);

		WarmManifold wm = {0};
		wm.count = m->contact_count;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			v3 world_pt = scale(add(
				add(w->body_hot[m->body_a].position, s->r_a),
				add(w->body_hot[m->body_b].position, s->r_b)), 0.5f);
			wm.contacts[ci] = (WarmContact){
				.point = world_pt,
				.lambda_n = s->lambda_n,
				.lambda_t1 = s->lambda_t1,
				.lambda_t2 = s->lambda_t2,
			};
		}

		map_set(w->warm_cache, key, wm);
	}

	afree(sm);
	afree(sc);
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
	map_free(w->warm_cache);
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

	// Sequential impulses solver
	int manifold_count = asize(manifolds);
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, dt);

	for (int iter = 0; iter < SOLVER_VELOCITY_ITERS; iter++)
		solver_iterate(w, sm, asize(sm), sc);

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

	// Position correction pass (after integration)
	solver_position_correct(w, sm, asize(sm), sc);

	// Store warm starting data and free solver arrays
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	afree(manifolds);
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
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
		.friction = fric,
		.restitution = params.restitution,
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
