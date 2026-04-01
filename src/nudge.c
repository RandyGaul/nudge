// See LICENSE for licensing info.
// nudge.c -- physics world implementation

// Portable float validity: rejects NaN, inf, and extreme magnitudes.
// NaN fails the self-equality test; inf and huge values fail the bound test.
static inline int float_valid(float f) { return f == f && f > -1e18f && f < 1e18f; }
static inline int v3_is_valid(v3 v) { return float_valid(v.x) && float_valid(v.y) && float_valid(v.z); }
static inline int quat_is_valid(quat q) { return float_valid(q.x) && float_valid(q.y) && float_valid(q.z) && float_valid(q.w); }
#define is_valid(x) _Generic((x), float: float_valid, v3: v3_is_valid, quat: quat_is_valid)(x)

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

// Joint persistent storage (handle-based, parallel arrays like bodies).
typedef enum JointType { JOINT_BALL_SOCKET, JOINT_DISTANCE } JointType;

typedef struct JointInternal
{
	JointType type;
	int body_a, body_b; // array indices (resolved from handles at creation)
	union {
		struct { v3 local_a, local_b; SpringParams spring; } ball_socket;
		struct { v3 local_a, local_b; float rest_length; SpringParams spring; } distance;
	};
	// Warm starting: accumulated impulses persisted across frames.
	union {
		v3 warm_lambda3;   // ball socket (3 DOF)
		float warm_lambda1; // distance (1 DOF)
	};
} JointInternal;

typedef struct WorldInternal
{
	v3 gravity;
	CK_DYNA BodyCold*    body_cold;
	CK_DYNA BodyHot*     body_hot;
	CK_DYNA uint32_t*    body_gen;
	CK_DYNA int*         body_free;
	CK_DYNA Contact*     debug_contacts;
	CK_MAP(WarmManifold) warm_cache;
	// Joints
	CK_DYNA JointInternal* joints;
	CK_DYNA uint32_t*      joint_gen;
	CK_DYNA int*           joint_free;
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
#define SOLVER_MAX_LINEAR_VEL      500.0f
#define SOLVER_MAX_ANGULAR_VEL     100.0f

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
// Joint solver types.

typedef struct SolverBallSocket
{
	int body_a, body_b;
	v3 r_a, r_b;
	float eff_mass[6]; // symmetric 3x3 inverse (xx,xy,xz,yy,yz,zz)
	v3 bias;
	float softness;
	v3 lambda;
	int joint_idx;
} SolverBallSocket;

typedef struct SolverDistance
{
	int body_a, body_b;
	v3 r_a, r_b;
	v3 axis;
	float eff_mass;
	float bias;
	float softness;
	float lambda;
	int joint_idx;
} SolverDistance;

// Constraint ref for graph coloring dispatch.
enum { CTYPE_CONTACT, CTYPE_BALL_SOCKET, CTYPE_DISTANCE };

typedef struct ConstraintRef
{
	uint8_t type;
	uint8_t color;
	int index;
	int body_a, body_b;
} ConstraintRef;

// Symmetric 3x3 stored as 6 floats: [xx, xy, xz, yy, yz, zz].
static void sym3x3_inverse(const float* in, float* out)
{
	float a = in[0], b = in[1], c = in[2];
	float d = in[3], e = in[4], f = in[5];
	float det = a*(d*f - e*e) - b*(b*f - c*e) + c*(b*e - c*d);
	float inv_det = det != 0.0f ? 1.0f / det : 0.0f;
	out[0] = (d*f - e*e) * inv_det;
	out[1] = (c*e - b*f) * inv_det;
	out[2] = (b*e - c*d) * inv_det;
	out[3] = (a*f - c*c) * inv_det;
	out[4] = (b*c - a*e) * inv_det;
	out[5] = (a*d - b*b) * inv_det;
}

static v3 sym3x3_mul_v3(const float* m, v3 v)
{
	return V3(
		m[0]*v.x + m[1]*v.y + m[2]*v.z,
		m[1]*v.x + m[3]*v.y + m[4]*v.z,
		m[2]*v.x + m[4]*v.y + m[5]*v.z);
}

// Convert spring params to solver coefficients (bepu approach).
static void spring_compute(SpringParams sp, float dt,
	float* pos_to_vel, float* softness)
{
	if (sp.frequency <= 0.0f) {
		// Rigid constraint: Baumgarte stabilization with fractional correction.
		*pos_to_vel = dt > 0.0f ? SOLVER_BAUMGARTE / dt : 0.0f;
		*softness = 0.0f;
		return;
	}
	float omega = 2.0f * 3.14159265f * sp.frequency;
	float d = 2.0f * sp.damping_ratio * omega;
	float k = omega * omega;
	float hd = dt * d, hk = dt * k, hhk = dt * hk;
	float denom = hd + hhk;
	if (denom < 1e-12f) { *pos_to_vel = 0; *softness = 0; return; }
	*softness = 1.0f / denom;
	*pos_to_vel = hk * (*softness);
}

// Build ball socket effective mass (symmetric 3x3).
// K = (inv_mass_a + inv_mass_b)*I + skew(r_a)^T * I_a^-1 * skew(r_a) + same for B.
static void ball_socket_eff_mass(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b,
	float softness, float* out)
{
	float inv_m = a->inv_mass + b->inv_mass;
	// Skew-sandwich: skew(r)^T * I^-1 * skew(r) for diagonal inertia.
	// For body A: I_a^-1 is diagonal in local space, r_a is world space.
	// Transform r_a to local, compute product, transform back.
	// Or equivalently use the explicit formula.
	float K[6] = { inv_m, 0, 0, inv_m, 0, inv_m }; // start with linear term

	// Add angular contribution for body A
	v3 ia = a->inv_inertia_local;
	if (ia.x > 0 || ia.y > 0 || ia.z > 0) {
		// columns of world-space inverse inertia: R * diag(ia) * R^T
		// skew(r)^T * I_world^-1 * skew(r) computed column by column
		v3 e0 = inv_inertia_mul(a->rotation, ia, V3(0, -r_a.z, r_a.y));
		v3 e1 = inv_inertia_mul(a->rotation, ia, V3(r_a.z, 0, -r_a.x));
		v3 e2 = inv_inertia_mul(a->rotation, ia, V3(-r_a.y, r_a.x, 0));
		// skew(r)^T * [e0 e1 e2] -- but skew(r)^T has rows [0 -rz ry; rz 0 -rx; -ry rx 0]
		// Result[i][j] = dot(skew_row_i(r), e_j)
		// Row 0 of skew^T: (0, -rz, ry)
		K[0] += -r_a.z*e0.y + r_a.y*e0.z; // xx
		K[1] += -r_a.z*e1.y + r_a.y*e1.z; // xy
		K[2] += -r_a.z*e2.y + r_a.y*e2.z; // xz
		K[3] +=  r_a.z*e1.x - r_a.x*e1.z; // yy
		K[4] +=  r_a.z*e2.x - r_a.x*e2.z; // yz
		K[5] += -r_a.y*e2.x + r_a.x*e2.y; // zz
	}

	v3 ib = b->inv_inertia_local;
	if (ib.x > 0 || ib.y > 0 || ib.z > 0) {
		v3 e0 = inv_inertia_mul(b->rotation, ib, V3(0, -r_b.z, r_b.y));
		v3 e1 = inv_inertia_mul(b->rotation, ib, V3(r_b.z, 0, -r_b.x));
		v3 e2 = inv_inertia_mul(b->rotation, ib, V3(-r_b.y, r_b.x, 0));
		K[0] += -r_b.z*e0.y + r_b.y*e0.z;
		K[1] += -r_b.z*e1.y + r_b.y*e1.z;
		K[2] += -r_b.z*e2.y + r_b.y*e2.z;
		K[3] +=  r_b.z*e1.x - r_b.x*e1.z;
		K[4] +=  r_b.z*e2.x - r_b.x*e2.z;
		K[5] += -r_b.y*e2.x + r_b.x*e2.y;
	}

	// Add softness to diagonal: K_eff = K + gamma * I
	K[0] += softness; K[3] += softness; K[5] += softness;
	sym3x3_inverse(K, out);
}

// Pre-solve joints: build solver arrays from persistent JointInternal data.
static void joints_pre_solve(WorldInternal* w, float dt,
	SolverBallSocket** out_bs, SolverDistance** out_dist)
{
	CK_DYNA SolverBallSocket* bs = NULL;
	CK_DYNA SolverDistance* dist = NULL;
	int joint_count = asize(w->joints);

	for (int i = 0; i < joint_count; i++) {
		if (!split_alive(w->joint_gen, i)) continue;
		JointInternal* j = &w->joints[i];
		BodyHot* a = &w->body_hot[j->body_a];
		BodyHot* b = &w->body_hot[j->body_b];

		if (j->type == JOINT_BALL_SOCKET) {
			SolverBallSocket s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->ball_socket.local_a);
			s.r_b = rotate(b->rotation, j->ball_socket.local_b);

			float ptv, soft;
			spring_compute(j->ball_socket.spring, dt, &ptv, &soft);
			s.softness = soft;
			ball_socket_eff_mass(a, b, s.r_a, s.r_b, soft, s.eff_mass);

			// Position error: world anchor B - world anchor A
			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			s.bias = scale(sub(anchor_b, anchor_a), ptv);

			// Warm start from persistent storage
			s.lambda = j->warm_lambda3;

			apush(bs, s);
		} else if (j->type == JOINT_DISTANCE) {
			SolverDistance s = {0};
			s.body_a = j->body_a;
			s.body_b = j->body_b;
			s.joint_idx = i;
			s.r_a = rotate(a->rotation, j->distance.local_a);
			s.r_b = rotate(b->rotation, j->distance.local_b);

			v3 anchor_a = add(a->position, s.r_a);
			v3 anchor_b = add(b->position, s.r_b);
			v3 delta = sub(anchor_b, anchor_a);
			float dist_val = len(delta);
			s.axis = dist_val > 1e-6f ? scale(delta, 1.0f / dist_val) : V3(1, 0, 0);

			float ptv, soft;
			spring_compute(j->distance.spring, dt, &ptv, &soft);
			s.softness = soft;

			float inv_mass_sum = a->inv_mass + b->inv_mass;
			float k = inv_mass_sum
				+ dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local,
					cross(s.r_a, s.axis)), s.r_a), s.axis)
				+ dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local,
					cross(s.r_b, s.axis)), s.r_b), s.axis);
			k += soft;
			s.eff_mass = k > 1e-12f ? 1.0f / k : 0.0f;

			float error = dist_val - j->distance.rest_length;
			s.bias = -ptv * error;

			s.lambda = j->warm_lambda1;

			apush(dist, s);
		}
	}

	*out_bs = bs;
	*out_dist = dist;
}

// Apply warm start impulses for joints.
static void joints_warm_start(WorldInternal* w,
	SolverBallSocket* bs, int bs_count,
	SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++) {
		SolverBallSocket* s = &bs[i];
		if (s->lambda.x == 0 && s->lambda.y == 0 && s->lambda.z == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, s->lambda);
	}
	for (int i = 0; i < dist_count; i++) {
		SolverDistance* s = &dist[i];
		if (s->lambda == 0) continue;
		apply_impulse(&w->body_hot[s->body_a], &w->body_hot[s->body_b],
			s->r_a, s->r_b, scale(s->axis, s->lambda));
	}
}

static void solve_ball_socket(WorldInternal* w, SolverBallSocket* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	v3 rhs = sub(neg(add(dv, s->bias)), scale(s->lambda, s->softness));
	v3 impulse = sym3x3_mul_v3(s->eff_mass, rhs);
	s->lambda = add(s->lambda, impulse);
	apply_impulse(a, b, s->r_a, s->r_b, impulse);
}

static void solve_distance(WorldInternal* w, SolverDistance* s)
{
	BodyHot* a = &w->body_hot[s->body_a];
	BodyHot* b = &w->body_hot[s->body_b];
	v3 dv = sub(
		add(b->velocity, cross(b->angular_velocity, s->r_b)),
		add(a->velocity, cross(a->angular_velocity, s->r_a)));
	float cdot = dot(dv, s->axis);
	float lambda = s->eff_mass * (-cdot + s->bias - s->softness * s->lambda);
	s->lambda += lambda;
	apply_impulse(a, b, s->r_a, s->r_b, scale(s->axis, lambda));
}

// Store joint accumulated impulses back to persistent storage.
static void joints_post_solve(WorldInternal* w,
	SolverBallSocket* bs, int bs_count,
	SolverDistance* dist, int dist_count)
{
	for (int i = 0; i < bs_count; i++)
		w->joints[bs[i].joint_idx].warm_lambda3 = bs[i].lambda;
	for (int i = 0; i < dist_count; i++)
		w->joints[dist[i].joint_idx].warm_lambda1 = dist[i].lambda;
	afree(bs);
	afree(dist);
}

// -----------------------------------------------------------------------------
// Graph coloring: assign colors to constraints so no two in same color share a body.
// Uses uint64_t bitmask per body (max 64 colors).

static void color_constraints(ConstraintRef* refs, int count, int body_count,
	int* out_batch_starts, int* out_color_count)
{
	uint64_t* body_colors = (uint64_t*)CK_ALLOC(body_count * sizeof(uint64_t));
	memset(body_colors, 0, body_count * sizeof(uint64_t));
	int max_color = 0;

	for (int i = 0; i < count; i++) {
		uint64_t used = body_colors[refs[i].body_a] | body_colors[refs[i].body_b];
		uint64_t avail = ~used;
		int color = 0;
		if (avail) {
			// Find lowest set bit (MSVC: _BitScanForward64, GCC: __builtin_ctzll)
#ifdef _MSC_VER
			unsigned long idx;
			_BitScanForward64(&idx, avail);
			color = (int)idx;
#else
			color = __builtin_ctzll(avail);
#endif
		}
		assert(color < 64);
		refs[i].color = (uint8_t)color;
		uint64_t bit = 1ULL << color;
		body_colors[refs[i].body_a] |= bit;
		body_colors[refs[i].body_b] |= bit;
		if (color > max_color) max_color = color;
	}

	// Counting sort by color
	int color_count = max_color + 1;
	int counts[64] = {0};
	for (int i = 0; i < count; i++) counts[refs[i].color]++;

	int offsets[64];
	offsets[0] = 0;
	for (int c = 1; c < color_count; c++) offsets[c] = offsets[c-1] + counts[c-1];

	// Record batch starts before sorting
	for (int c = 0; c < color_count; c++) out_batch_starts[c] = offsets[c];
	out_batch_starts[color_count] = count; // sentinel

	ConstraintRef* sorted = (ConstraintRef*)CK_ALLOC(count * sizeof(ConstraintRef));
	for (int i = 0; i < count; i++)
		sorted[offsets[refs[i].color]++] = refs[i];
	memcpy(refs, sorted, count * sizeof(ConstraintRef));

	CK_FREE(sorted);
	CK_FREE(body_colors);
	*out_color_count = color_count;
}

// Dispatch a single constraint solve by type.
static void solve_constraint(WorldInternal* w, ConstraintRef* ref,
	SolverManifold* sm, SolverContact* sc,
	SolverBallSocket* bs, SolverDistance* dist)
{
	switch (ref->type) {
	case CTYPE_CONTACT: {
		SolverManifold* m = &sm[ref->index];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			v3 dv = sub(
				add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));

			float vn = dot(dv, s->normal);
			float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce));
			float old_n = s->lambda_n;
			s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n - old_n));

			dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float max_f = m->friction * s->lambda_n;
			float vt1 = dot(dv, s->tangent1);
			float old_t1 = s->lambda_t1;
			s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + s->eff_mass_t1*(-vt1), max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

			dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)),
				add(a->velocity, cross(a->angular_velocity, s->r_a)));
			float vt2 = dot(dv, s->tangent2);
			float old_t2 = s->lambda_t2;
			s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + s->eff_mass_t2*(-vt2), max_f));
			apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
		}
		break;
	}
	case CTYPE_BALL_SOCKET: solve_ball_socket(w, &bs[ref->index]); break;
	case CTYPE_DISTANCE:    solve_distance(w, &dist[ref->index]); break;
	}
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
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
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

	// Pre-solve contacts
	int manifold_count = asize(manifolds);
	SolverManifold* sm = NULL;
	SolverContact*  sc = NULL;
	solver_pre_solve(w, manifolds, manifold_count, &sm, &sc, dt);

	// Pre-solve joints
	SolverBallSocket* sol_bs = NULL;
	SolverDistance*    sol_dist = NULL;
	joints_pre_solve(w, dt, &sol_bs, &sol_dist);
	joints_warm_start(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	// Build constraint refs for graph coloring
	CK_DYNA ConstraintRef* crefs = NULL;
	int sm_count = asize(sm);
	for (int i = 0; i < sm_count; i++) {
		ConstraintRef r = { .type = CTYPE_CONTACT, .index = i,
			.body_a = sm[i].body_a, .body_b = sm[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_bs); i++) {
		ConstraintRef r = { .type = CTYPE_BALL_SOCKET, .index = i,
			.body_a = sol_bs[i].body_a, .body_b = sol_bs[i].body_b };
		apush(crefs, r);
	}
	for (int i = 0; i < asize(sol_dist); i++) {
		ConstraintRef r = { .type = CTYPE_DISTANCE, .index = i,
			.body_a = sol_dist[i].body_a, .body_b = sol_dist[i].body_b };
		apush(crefs, r);
	}

	// Graph color and sort
	int cref_count = asize(crefs);
	int batch_starts[65] = {0};
	int color_count = 0;
	if (cref_count > 0)
		color_constraints(crefs, cref_count, count, batch_starts, &color_count);

	// Velocity iterations with colored batches
	for (int iter = 0; iter < SOLVER_VELOCITY_ITERS; iter++)
		for (int c = 0; c < color_count; c++)
			for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
				solve_constraint(w, &crefs[i], sm, sc, sol_bs, sol_dist);

	afree(crefs);

	// Integrate positions and rotations
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;

		// Clamp velocities to prevent solver divergence from blowing up.
		float lv2 = len2(h->velocity);
		if (lv2 > SOLVER_MAX_LINEAR_VEL * SOLVER_MAX_LINEAR_VEL)
			h->velocity = scale(h->velocity, SOLVER_MAX_LINEAR_VEL / sqrtf(lv2));
		float av2 = len2(h->angular_velocity);
		if (av2 > SOLVER_MAX_ANGULAR_VEL * SOLVER_MAX_ANGULAR_VEL)
			h->angular_velocity = scale(h->angular_velocity, SOLVER_MAX_ANGULAR_VEL / sqrtf(av2));

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
		if (ql < 1e-15f) ql = 1.0f;  // prevent div-by-zero if quat collapsed
		float inv_ql = 1.0f / ql;
		h->rotation.x *= inv_ql; h->rotation.y *= inv_ql;
		h->rotation.z *= inv_ql; h->rotation.w *= inv_ql;
	}

	// Position correction pass (contacts only, after integration)
	solver_position_correct(w, sm, asize(sm), sc);

	// Store warm starting data and free solver arrays
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));
	afree(manifolds);
}

// -----------------------------------------------------------------------------
// Body.

Body create_body(World world, BodyParams params)
{
	assert(is_valid(params.position) && "create_body: position is NaN/inf");
	assert(is_valid(params.rotation) && "create_body: rotation is NaN/inf");
	assert(is_valid(params.mass) && params.mass >= 0.0f && "create_body: mass must be >= 0 and finite");

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
	assert(is_valid(params.local_pos) && "body_add_shape: local_pos is NaN/inf");

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

// -----------------------------------------------------------------------------
// Joints.

Joint create_ball_socket(World world, BallSocketParams params)
{
	assert(is_valid(params.local_offset_a) && "create_ball_socket: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_ball_socket: local_offset_b is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int idx;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	// Grow joint arrays manually (no split_add -- joints don't need hot/cold split)
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1); // odd = alive
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_BALL_SOCKET,
		.body_a = ba, .body_b = bb,
		.ball_socket = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.spring = params.spring,
		},
	};
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

Joint create_distance(World world, DistanceParams params)
{
	assert(is_valid(params.local_offset_a) && "create_distance: local_offset_a is NaN/inf");
	assert(is_valid(params.local_offset_b) && "create_distance: local_offset_b is NaN/inf");
	assert(is_valid(params.rest_length) && "create_distance: rest_length is NaN/inf");

	WorldInternal* w = (WorldInternal*)world.id;
	int ba = handle_index(params.body_a);
	int bb = handle_index(params.body_b);
	assert(split_valid(w->body_gen, params.body_a));
	assert(split_valid(w->body_gen, params.body_b));

	int idx;
	if (asize(w->joint_free) > 0) {
		idx = apop(w->joint_free);
		w->joint_gen[idx]++;
	} else {
		idx = asize(w->joints);
		JointInternal zero = {0};
		apush(w->joints, zero);
		apush(w->joint_gen, 1);
	}

	// Auto-compute rest length if not specified
	float rest = params.rest_length;
	if (rest <= 0.0f) {
		BodyHot* a = &w->body_hot[ba];
		BodyHot* b = &w->body_hot[bb];
		v3 wa = add(a->position, rotate(a->rotation, params.local_offset_a));
		v3 wb = add(b->position, rotate(b->rotation, params.local_offset_b));
		rest = len(sub(wb, wa));
	}

	w->joints[idx] = (JointInternal){
		.type = JOINT_DISTANCE,
		.body_a = ba, .body_b = bb,
		.distance = {
			.local_a = params.local_offset_a,
			.local_b = params.local_offset_b,
			.rest_length = rest,
			.spring = params.spring,
		},
	};
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	memset(&w->joints[idx], 0, sizeof(JointInternal));
	w->joint_gen[idx]++; // even = dead
	apush(w->joint_free, idx);
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}
