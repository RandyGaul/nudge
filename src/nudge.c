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
	int bvh_leaf;     // -1 if not in BVH
	int island_id;    // -1 = no island (static or unconnected)
	int island_prev;  // prev body in island body list, -1 = head
	int island_next;  // next body in island body list, -1 = tail
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
	float linear_damping;
	float angular_damping;
	float sleep_time; // accumulated seconds below velocity threshold
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
	// Island linked list fields.
	int island_id;    // -1 = none
	int island_prev;  // -1 = head
	int island_next;  // -1 = tail
} JointInternal;

typedef struct BVHTree BVHTree; // forward decl, defined in bvh.c

// Island: group of connected bodies that can sleep/wake together.
typedef struct Island
{
	int head_body, tail_body, body_count;
	int head_joint, tail_joint, joint_count;
	int constraint_remove_count;
	int awake; // 1 = awake, 0 = sleeping
} Island;

#define SLEEP_VEL_THRESHOLD   0.01f  // squared velocity magnitude threshold
#define SLEEP_TIME_THRESHOLD  0.5f   // seconds of stillness before sleep
#define LINEAR_SLOP           0.01f  // contact margin: keeps contacts alive near zero-penetration

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
	// Broadphase
	BroadphaseType broadphase_type;
	BVHTree* bvh_static;
	BVHTree* bvh_dynamic;
	// Islands
	CK_DYNA Island*     islands;
	CK_DYNA uint32_t*   island_gen;
	CK_DYNA int*        island_free;
	CK_MAP(uint8_t)     prev_touching; // body_pair_key -> 1 for pairs touching last frame
	int sleep_enabled; // 1 = island sleep active (default)
	FrictionModel friction_model;
	int velocity_iters;
	int position_iters;
	float contact_hertz;
	float contact_damping_ratio;
	float max_push_velocity;
	int sub_steps;
} WorldInternal;

#include "gjk.c"
#include "gjk_gino.c"
#include "quickhull.c"
#include "bvh.c"
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
#define SOLVER_BAUMGARTE           0.2f   // used by joints only; contacts use NGS
#define SOLVER_SLOP                0.005f // NGS position correction dead zone
#define SOLVER_RESTITUTION_THRESH  1.0f
#define SOLVER_POS_BAUMGARTE       0.2f
#define SOLVER_POS_MAX_CORRECTION  0.2f   // max position correction per step
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
	float softness;      // soft constraint regularization (0 = rigid/NGS)
	float penetration;   // cached for position correction pass
	uint32_t feature_id; // geometric feature key for warm starting
} SolverContact;

typedef struct SolverManifold
{
	int body_a, body_b;
	int contact_start;
	int contact_count;
	float friction;
	// Manifold-level patch friction data (FRICTION_PATCH only)
	v3 centroid_r_a;
	v3 centroid_r_b;
	v3 tangent1, tangent2;
	v3 normal;
	float eff_mass_t1, eff_mass_t2;
	float eff_mass_twist;
	float lambda_t1, lambda_t2;
	float lambda_twist;
	float patch_area;
	float patch_radius;
} SolverManifold;

// Warm starting: cached impulses from previous frame, keyed by body pair.
typedef struct WarmContact
{
	uint32_t feature_id; // geometric feature key for matching
	v3 r_a;              // body-A-relative position for spatial fallback matching
	float lambda_n;
	float lambda_t1;
	float lambda_t2;
} WarmContact;

struct WarmManifold
{
	WarmContact contacts[MAX_CONTACTS];
	int count;
	int stale; // 0 = updated this frame, incremented each frame not touched, evicted at >1
	// Manifold-level friction warm data (FRICTION_PATCH)
	float manifold_lambda_t1;
	float manifold_lambda_t2;
	float manifold_lambda_twist;
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

#define PATCH_MIN_AREA 0.001f

// Estimate contact patch area from manifold points projected onto contact plane.
static float estimate_patch_area(Contact* contacts, int count)
{
	if (count < 3) return PATCH_MIN_AREA;
	// Fan triangulation from contacts[0]
	float area = 0.0f;
	for (int i = 1; i < count - 1; i++)
		area += 0.5f * len(cross(sub(contacts[i].point, contacts[0].point), sub(contacts[i + 1].point, contacts[0].point)));
	return area > PATCH_MIN_AREA ? area : PATCH_MIN_AREA;
}

static float compute_effective_mass(BodyHot* a, BodyHot* b, float inv_mass_sum, v3 r_a, v3 r_b, v3 dir)
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
	a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(r_a, impulse)));
	b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(r_b, impulse)));
}

// Match a new contact to a cached contact by feature ID. Returns index or -1.
static int warm_match(WarmManifold* wm, uint32_t feature_id)
{
	for (int i = 0; i < wm->count; i++)
		if (wm->contacts[i].feature_id == feature_id) return i;
	return -1;
}

static void solver_pre_solve(WorldInternal* w, InternalManifold* manifolds, int manifold_count, SolverManifold** out_sm, SolverContact** out_sc, float dt)
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

		int patch_mode = (w->friction_model == FRICTION_PATCH);

		for (int c = 0; c < im->m.count; c++) {
			Contact* ct = &im->m.contacts[c];
			SolverContact s = {0};

			s.r_a = sub(ct->point, a->position);
			s.r_b = sub(ct->point, b->position);
			s.normal = ct->normal;
			s.penetration = ct->penetration;
			s.feature_id = ct->feature_id;

			s.eff_mass_n = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.normal);

			// Soft contact constraint (Box2D b2MakeSoft)
			{
				float hertz = w->contact_hertz;
				if (a->inv_mass == 0.0f || b->inv_mass == 0.0f) hertz *= 2.0f;
				float omega = 2.0f * 3.14159265f * hertz;
				float d = 2.0f * w->contact_damping_ratio * omega;
				float k = omega * omega;
				float hd = dt * d, hhk = dt * dt * k;
				float denom = hd + hhk;
				if (denom > 1e-12f) {
					s.softness = 1.0f / denom;
					float bias_rate = dt * k * s.softness;
					float K = s.eff_mass_n > 0.0f ? 1.0f / s.eff_mass_n : 0.0f;
					s.eff_mass_n = (K + s.softness) > 1e-12f ? 1.0f / (K + s.softness) : 0.0f;
					float pen = ct->penetration - SOLVER_SLOP;
					s.bias = pen > 0.0f ? -bias_rate * pen : 0.0f;
					if (s.bias < -w->max_push_velocity) s.bias = -w->max_push_velocity;
				}
			}

			if (!patch_mode) {
				contact_tangent_basis(ct->normal, &s.tangent1, &s.tangent2);
				s.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent1);
				s.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, s.r_a, s.r_b, s.tangent2);
			}

			v3 vel_a = add(a->velocity, cross(a->angular_velocity, s.r_a));
			v3 vel_b = add(b->velocity, cross(b->angular_velocity, s.r_b));
			float vn_rel = dot(sub(vel_b, vel_a), ct->normal);
			s.bounce = (-vn_rel > SOLVER_RESTITUTION_THRESH) ? rest * vn_rel : 0.0f;

			apush(sc, s);
		}

		// Patch friction: compute centroid, patch area, tangent basis at manifold level.
		if (patch_mode) {
			v3 centroid_a = V3(0, 0, 0), centroid_b = V3(0, 0, 0);
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				centroid_a = add(centroid_a, s->r_a);
				centroid_b = add(centroid_b, s->r_b);
			}
			float inv_n = 1.0f / (float)smf.contact_count;
			smf.centroid_r_a = scale(centroid_a, inv_n);
			smf.centroid_r_b = scale(centroid_b, inv_n);

			// Use first contact normal for tangent basis (all share same normal in a manifold)
			v3 n = sc[smf.contact_start].normal;
			smf.normal = n;
			contact_tangent_basis(n, &smf.tangent1, &smf.tangent2);
			smf.eff_mass_t1 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent1);
			smf.eff_mass_t2 = compute_effective_mass(a, b, inv_mass_sum, smf.centroid_r_a, smf.centroid_r_b, smf.tangent2);
			smf.patch_area = estimate_patch_area(im->m.contacts, im->m.count);
			smf.patch_radius = 0.6667f * sqrtf(smf.patch_area * (1.0f / 3.14159265f));

			// Torsional friction effective mass: 1 / (n^T * I_a_inv * n + n^T * I_b_inv * n)
			float k_twist = dot(inv_inertia_mul(a->rotation, a->inv_inertia_local, n), n) + dot(inv_inertia_mul(b->rotation, b->inv_inertia_local, n), n);
			smf.eff_mass_twist = k_twist > 1e-12f ? 1.0f / k_twist : 0.0f;
		}

		// Warm start: match new contacts to cached contacts.
		// Pass 1: exact feature ID match.
		// Pass 2: spatial fallback -- closest unmatched old contact by r_a distance.
		// Pass 3: redistribute any remaining unmatched old impulse evenly.
		if (wm) {
			int old_matched[MAX_CONTACTS] = {0};
			int new_matched[MAX_CONTACTS] = {0};
			// Pass 1: feature ID
			for (int c = 0; c < smf.contact_count; c++) {
				SolverContact* s = &sc[smf.contact_start + c];
				if (s->feature_id == 0) continue;
				int match = warm_match(wm, s->feature_id);
				if (match >= 0) {
					s->lambda_n = wm->contacts[match].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[match].lambda_t1;
						s->lambda_t2 = wm->contacts[match].lambda_t2;
					}
					old_matched[match] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 2: spatial fallback for unmatched contacts
			float spatial_tol2 = 0.01f;
			for (int c = 0; c < smf.contact_count; c++) {
				if (new_matched[c]) continue;
				SolverContact* s = &sc[smf.contact_start + c];
				float best_d2 = spatial_tol2;
				int best = -1;
				for (int j = 0; j < wm->count; j++) {
					if (old_matched[j]) continue;
					float d2 = len2(sub(s->r_a, wm->contacts[j].r_a));
					if (d2 < best_d2) { best_d2 = d2; best = j; }
				}
				if (best >= 0) {
					s->lambda_n = wm->contacts[best].lambda_n;
					if (!patch_mode) {
						s->lambda_t1 = wm->contacts[best].lambda_t1;
						s->lambda_t2 = wm->contacts[best].lambda_t2;
					}
					old_matched[best] = 1;
					new_matched[c] = 1;
				}
			}
			// Pass 3: redistribute remaining unmatched old impulse
			int new_unmatched = 0;
			for (int c = 0; c < smf.contact_count; c++)
				if (!new_matched[c]) new_unmatched++;
			if (new_unmatched > 0) {
				float leftover_n = 0, leftover_t1 = 0, leftover_t2 = 0;
				for (int j = 0; j < wm->count; j++) {
					if (!old_matched[j]) {
						leftover_n += wm->contacts[j].lambda_n;
						if (!patch_mode) {
							leftover_t1 += wm->contacts[j].lambda_t1;
							leftover_t2 += wm->contacts[j].lambda_t2;
						}
					}
				}
				float share = 1.0f / (float)new_unmatched;
				for (int c = 0; c < smf.contact_count; c++) {
					if (new_matched[c]) continue;
					SolverContact* s = &sc[smf.contact_start + c];
					s->lambda_n = leftover_n * share;
					if (!patch_mode) {
						s->lambda_t1 = leftover_t1 * share;
						s->lambda_t2 = leftover_t2 * share;
					}
				}
			}
			// Warm start manifold-level friction
			if (patch_mode) {
				smf.lambda_t1 = wm->manifold_lambda_t1;
				smf.lambda_t2 = wm->manifold_lambda_t2;
				smf.lambda_twist = wm->manifold_lambda_twist;
			}
		}

		// Speculative contacts (negative penetration, kept alive by margin) must
		// not carry warm-started impulse -- they exist for cache continuity only.
		for (int c = 0; c < smf.contact_count; c++) {
			SolverContact* s = &sc[smf.contact_start + c];
			if (s->penetration < 0.0f) { s->lambda_n = 0.0f; s->lambda_t1 = 0.0f; s->lambda_t2 = 0.0f; }
		}

		apush(sm, smf);
	}

	// Apply warm start impulses
	int patch_warm = (w->friction_model == FRICTION_PATCH);
	for (int i = 0; i < asize(sm); i++) {
		SolverManifold* m = &sm[i];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			if (patch_warm) {
				if (s->lambda_n == 0.0f) continue;
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n));
			} else {
				if (s->lambda_n == 0.0f && s->lambda_t1 == 0.0f && s->lambda_t2 == 0.0f)
					continue;
				v3 P = add(add(scale(s->normal, s->lambda_n), scale(s->tangent1, s->lambda_t1)), scale(s->tangent2, s->lambda_t2));
				apply_impulse(a, b, s->r_a, s->r_b, P);
			}
		}
		// Warm start manifold-level friction
		if (patch_warm && (m->lambda_t1 != 0.0f || m->lambda_t2 != 0.0f)) {
			v3 P = add(scale(m->tangent1, m->lambda_t1), scale(m->tangent2, m->lambda_t2));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, P);
		}
		// Warm start torsional friction (pure angular impulse along normal)
		if (patch_warm && m->lambda_twist != 0.0f) {
			v3 twist_impulse = scale(m->normal, m->lambda_twist);
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, twist_impulse));
		}
	}

	*out_sm = sm;
	*out_sc = sc;
}

static void solver_iterate(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
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
static void solver_position_correct(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc)
{
	for (int iter = 0; iter < w->position_iters; iter++) {
		for (int i = 0; i < sm_count; i++) {
			SolverManifold* m = &sm[i];
			BodyHot* a = &w->body_hot[m->body_a];
			BodyHot* b = &w->body_hot[m->body_b];
			float inv_mass_sum = a->inv_mass + b->inv_mass;

			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];

				// Recompute separation from current positions
				v3 r_a = sub(add(a->position, rotate(a->rotation, rotate(inv(a->rotation), s->r_a))), a->position);
				v3 r_b = sub(add(b->position, rotate(b->rotation, rotate(inv(b->rotation), s->r_b))), b->position);
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

// Persistent warm cache: update active pairs, age stale pairs, evict old stale.
// BEPU-style freshness: entries survive one extra frame so warm data isn't lost
// when narrowphase misses a pair for a single frame (FP noise).
static void solver_post_solve(WorldInternal* w, SolverManifold* sm, int sm_count, SolverContact* sc, InternalManifold* manifolds, int manifold_count)
{
	// Update active pairs (per sub-step)
	for (int i = 0; i < sm_count; i++) {
		SolverManifold* m = &sm[i];
		uint64_t key = body_pair_key(m->body_a, m->body_b);

		WarmManifold wm = {0};
		wm.count = m->contact_count;
		wm.stale = 0;
		for (int ci = 0; ci < m->contact_count; ci++) {
			SolverContact* s = &sc[m->contact_start + ci];
			wm.contacts[ci] = (WarmContact){
				.feature_id = s->feature_id,
				.r_a = s->r_a,
				.lambda_n = s->lambda_n,
				.lambda_t1 = s->lambda_t1,
				.lambda_t2 = s->lambda_t2,
			};
		}
		if (w->friction_model == FRICTION_PATCH) {
			wm.manifold_lambda_t1 = m->lambda_t1;
			wm.manifold_lambda_t2 = m->lambda_t2;
			wm.manifold_lambda_twist = m->lambda_twist;
		}

		map_set(w->warm_cache, key, wm);
	}

	afree(sm);
	afree(sc);
}

// Age and evict stale warm cache entries. Called once per frame (not per sub-step).
static void warm_cache_age_and_evict(WorldInternal* w)
{
	for (int i = 0; i < map_size(w->warm_cache); i++)
		w->warm_cache[i].stale++;
	int i = 0;
	while (i < map_size(w->warm_cache)) {
		if (w->warm_cache[i].stale > 1)
			map_del(w->warm_cache, map_keys(w->warm_cache)[i]);
		else
			i++;
	}
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
static void spring_compute(SpringParams sp, float dt, float* pos_to_vel, float* softness)
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
static void ball_socket_eff_mass(BodyHot* a, BodyHot* b, v3 r_a, v3 r_b, float softness, float* out)
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
static void joints_pre_solve(WorldInternal* w, float dt, SolverBallSocket** out_bs, SolverDistance** out_dist)
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
			float k = inv_mass_sum + dot(cross(inv_inertia_mul(a->rotation, a->inv_inertia_local, cross(s.r_a, s.axis)), s.r_a), s.axis) + dot(cross(inv_inertia_mul(b->rotation, b->inv_inertia_local, cross(s.r_b, s.axis)), s.r_b), s.axis);
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
static void joints_warm_start(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
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
static void joints_post_solve(WorldInternal* w, SolverBallSocket* bs, int bs_count, SolverDistance* dist, int dist_count)
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

static void color_constraints(ConstraintRef* refs, int count, int body_count, int* out_batch_starts, int* out_color_count)
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
static void solve_constraint(WorldInternal* w, ConstraintRef* ref, SolverManifold* sm, SolverContact* sc, SolverBallSocket* bs, SolverDistance* dist)
{
	switch (ref->type) {
	case CTYPE_CONTACT: {
		SolverManifold* m = &sm[ref->index];
		BodyHot* a = &w->body_hot[m->body_a];
		BodyHot* b = &w->body_hot[m->body_b];

		if (w->friction_model == FRICTION_PATCH) {
			// Solve all normal constraints first (decoupled from friction)
			float total_lambda_n = 0.0f;
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float vn = dot(dv, s->normal);
				float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
				float old_n = s->lambda_n;
				s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n - old_n));
				total_lambda_n += s->lambda_n;
			}

			// Manifold-level 2D friction at centroid, clamped by aggregate normal force
			float max_f = m->friction * total_lambda_n;
			v3 dv = sub(add(b->velocity, cross(b->angular_velocity, m->centroid_r_b)), add(a->velocity, cross(a->angular_velocity, m->centroid_r_a)));
			float vt1 = dot(dv, m->tangent1);
			float old_t1 = m->lambda_t1;
			m->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + m->eff_mass_t1 * (-vt1), max_f));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, scale(m->tangent1, m->lambda_t1 - old_t1));

			dv = sub(add(b->velocity, cross(b->angular_velocity, m->centroid_r_b)), add(a->velocity, cross(a->angular_velocity, m->centroid_r_a)));
			float vt2 = dot(dv, m->tangent2);
			float old_t2 = m->lambda_t2;
			m->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + m->eff_mass_t2 * (-vt2), max_f));
			apply_impulse(a, b, m->centroid_r_a, m->centroid_r_b, scale(m->tangent2, m->lambda_t2 - old_t2));

			// Torsional friction: resist spin around contact normal
			float max_twist = m->friction * total_lambda_n * m->patch_radius;
			float w_rel = dot(sub(b->angular_velocity, a->angular_velocity), m->normal);
			float lambda_tw = m->eff_mass_twist * (-w_rel);
			float old_tw = m->lambda_twist;
			m->lambda_twist = fmaxf(-max_twist, fminf(old_tw + lambda_tw, max_twist));
			float delta_tw = m->lambda_twist - old_tw;
			v3 twist_impulse = scale(m->normal, delta_tw);
			a->angular_velocity = sub(a->angular_velocity, inv_inertia_mul(a->rotation, a->inv_inertia_local, twist_impulse));
			b->angular_velocity = add(b->angular_velocity, inv_inertia_mul(b->rotation, b->inv_inertia_local, twist_impulse));
		} else {
			// Per-point Coulomb friction (original)
			for (int ci = 0; ci < m->contact_count; ci++) {
				SolverContact* s = &sc[m->contact_start + ci];
				v3 dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));

				float vn = dot(dv, s->normal);
				float lambda_n = s->eff_mass_n * (-(vn + s->bias + s->bounce) - s->softness * s->lambda_n);
				float old_n = s->lambda_n;
				s->lambda_n = fmaxf(old_n + lambda_n, 0.0f);
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->normal, s->lambda_n - old_n));

				dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float max_f = m->friction * s->lambda_n;
				float vt1 = dot(dv, s->tangent1);
				float old_t1 = s->lambda_t1;
				s->lambda_t1 = fmaxf(-max_f, fminf(old_t1 + s->eff_mass_t1*(-vt1), max_f));
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent1, s->lambda_t1 - old_t1));

				dv = sub(add(b->velocity, cross(b->angular_velocity, s->r_b)), add(a->velocity, cross(a->angular_velocity, s->r_a)));
				float vt2 = dot(dv, s->tangent2);
				float old_t2 = s->lambda_t2;
				s->lambda_t2 = fmaxf(-max_f, fminf(old_t2 + s->eff_mass_t2*(-vt2), max_f));
				apply_impulse(a, b, s->r_a, s->r_b, scale(s->tangent2, s->lambda_t2 - old_t2));
			}
		}
		break;
	}
	case CTYPE_BALL_SOCKET: solve_ball_socket(w, &bs[ref->index]); break;
	case CTYPE_DISTANCE:    solve_distance(w, &dist[ref->index]); break;
	}
}

// -----------------------------------------------------------------------------
// Island management.

static int island_create(WorldInternal* w)
{
	int idx;
	if (asize(w->island_free) > 0) {
		idx = apop(w->island_free);
		w->island_gen[idx]++;
	} else {
		idx = asize(w->islands);
		Island zero = {0};
		apush(w->islands, zero);
		apush(w->island_gen, 0);
	}
	w->islands[idx] = (Island){
		.head_body = -1, .tail_body = -1, .body_count = 0,
		.head_joint = -1, .tail_joint = -1, .joint_count = 0,
		.constraint_remove_count = 0, .awake = 1,
	};
	w->island_gen[idx] |= 1; // odd = alive
	return idx;
}

static void island_destroy(WorldInternal* w, int id)
{
	w->island_gen[id]++;
	if (w->island_gen[id] & 1) w->island_gen[id]++; // ensure even = dead
	apush(w->island_free, id);
}

static int island_alive(WorldInternal* w, int id)
{
	return id >= 0 && id < asize(w->islands) && (w->island_gen[id] & 1);
}

static void island_add_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	c->island_id = island_id;
	c->island_prev = isl->tail_body;
	c->island_next = -1;
	if (isl->tail_body >= 0)
		w->body_cold[isl->tail_body].island_next = body_idx;
	else
		isl->head_body = body_idx;
	isl->tail_body = body_idx;
	isl->body_count++;
}

static void island_remove_body(WorldInternal* w, int island_id, int body_idx)
{
	Island* isl = &w->islands[island_id];
	BodyCold* c = &w->body_cold[body_idx];
	if (c->island_prev >= 0)
		w->body_cold[c->island_prev].island_next = c->island_next;
	else
		isl->head_body = c->island_next;
	if (c->island_next >= 0)
		w->body_cold[c->island_next].island_prev = c->island_prev;
	else
		isl->tail_body = c->island_prev;
	c->island_id = -1;
	c->island_prev = -1;
	c->island_next = -1;
	isl->body_count--;
}

static void island_add_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	j->island_id = island_id;
	j->island_prev = isl->tail_joint;
	j->island_next = -1;
	if (isl->tail_joint >= 0)
		w->joints[isl->tail_joint].island_next = joint_idx;
	else
		isl->head_joint = joint_idx;
	isl->tail_joint = joint_idx;
	isl->joint_count++;
}

static void island_remove_joint(WorldInternal* w, int island_id, int joint_idx)
{
	Island* isl = &w->islands[island_id];
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_prev >= 0)
		w->joints[j->island_prev].island_next = j->island_next;
	else
		isl->head_joint = j->island_next;
	if (j->island_next >= 0)
		w->joints[j->island_next].island_prev = j->island_prev;
	else
		isl->tail_joint = j->island_prev;
	j->island_id = -1;
	j->island_prev = -1;
	j->island_next = -1;
	isl->joint_count--;
}

// Merge two islands. Keep the bigger one, walk the smaller one remapping bodies/joints.
static int island_merge(WorldInternal* w, int id_a, int id_b)
{
	if (id_a == id_b) return id_a;
	Island* a = &w->islands[id_a];
	Island* b = &w->islands[id_b];
	// Keep the bigger island
	if (a->body_count + a->joint_count < b->body_count + b->joint_count) {
		int tmp = id_a; id_a = id_b; id_b = tmp;
		a = &w->islands[id_a]; b = &w->islands[id_b];
	}
	// Remap smaller's bodies to bigger
	int bi = b->head_body;
	while (bi >= 0) {
		int next = w->body_cold[bi].island_next;
		w->body_cold[bi].island_id = id_a;
		bi = next;
	}
	// Splice body lists: append b's list to a's tail
	if (b->head_body >= 0) {
		if (a->tail_body >= 0) {
			w->body_cold[a->tail_body].island_next = b->head_body;
			w->body_cold[b->head_body].island_prev = a->tail_body;
		} else {
			a->head_body = b->head_body;
		}
		a->tail_body = b->tail_body;
		a->body_count += b->body_count;
	}
	// Remap smaller's joints
	int ji = b->head_joint;
	while (ji >= 0) {
		int next = w->joints[ji].island_next;
		w->joints[ji].island_id = id_a;
		ji = next;
	}
	// Splice joint lists
	if (b->head_joint >= 0) {
		if (a->tail_joint >= 0) {
			w->joints[a->tail_joint].island_next = b->head_joint;
			w->joints[b->head_joint].island_prev = a->tail_joint;
		} else {
			a->head_joint = b->head_joint;
		}
		a->tail_joint = b->tail_joint;
		a->joint_count += b->joint_count;
	}
	a->constraint_remove_count += b->constraint_remove_count;
	// If either was awake, result is awake
	if (b->awake) a->awake = 1;
	island_destroy(w, id_b);
	return id_a;
}

static void island_wake(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 1;
	int bi = isl->head_body;
	while (bi >= 0) {
		w->body_hot[bi].sleep_time = 0.0f;
		bi = w->body_cold[bi].island_next;
	}
}

static void island_sleep(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];
	isl->awake = 0;
	int bi = isl->head_body;
	while (bi >= 0) {
		w->body_hot[bi].velocity = V3(0, 0, 0);
		w->body_hot[bi].angular_velocity = V3(0, 0, 0);
		bi = w->body_cold[bi].island_next;
	}
}

// Link a newly created joint to islands: merge/create as needed, wake if sleeping.
static void link_joint_to_islands(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	int ba = j->body_a, bb = j->body_b;
	int is_static_a = w->body_hot[ba].inv_mass == 0.0f;
	int is_static_b = w->body_hot[bb].inv_mass == 0.0f;
	// Static bodies never get islands
	if (is_static_a && is_static_b) return;
	int isl_a = is_static_a ? -1 : w->body_cold[ba].island_id;
	int isl_b = is_static_b ? -1 : w->body_cold[bb].island_id;
	int target;
	if (isl_a >= 0 && isl_b >= 0) {
		target = island_merge(w, isl_a, isl_b);
	} else if (isl_a >= 0) {
		target = isl_a;
	} else if (isl_b >= 0) {
		target = isl_b;
	} else {
		target = island_create(w);
	}
	// Add dynamic bodies that aren't in the target island yet
	if (!is_static_a && w->body_cold[ba].island_id != target)
		island_add_body(w, target, ba);
	if (!is_static_b && w->body_cold[bb].island_id != target)
		island_add_body(w, target, bb);
	island_add_joint(w, target, joint_idx);
	// Wake if sleeping
	if (!w->islands[target].awake)
		island_wake(w, target);
}

// Unlink a joint from its island before destruction.
static void unlink_joint_from_island(WorldInternal* w, int joint_idx)
{
	JointInternal* j = &w->joints[joint_idx];
	if (j->island_id < 0) return;
	int isl = j->island_id;
	island_remove_joint(w, isl, joint_idx);
	w->islands[isl].constraint_remove_count++;
}

// Update islands from this frame's contact manifolds.
// Merge/create islands for touching dynamic body pairs.
// Detect lost contacts for island splitting.
static void islands_update_contacts(WorldInternal* w, InternalManifold* manifolds, int manifold_count)
{
	CK_MAP(uint8_t) curr_touching = NULL;

	for (int i = 0; i < manifold_count; i++) {
		int a = manifolds[i].body_a, b = manifolds[i].body_b;
		int is_static_a = w->body_hot[a].inv_mass == 0.0f;
		int is_static_b = w->body_hot[b].inv_mass == 0.0f;
		// Skip static-static (shouldn't happen, but guard)
		if (is_static_a && is_static_b) continue;

		uint64_t key = body_pair_key(a, b);
		map_set(curr_touching, key, (uint8_t)1);

		// Static bodies never get island_id
		int isl_a = is_static_a ? -1 : w->body_cold[a].island_id;
		int isl_b = is_static_b ? -1 : w->body_cold[b].island_id;

		int target;
		if (isl_a >= 0 && isl_b >= 0) {
			target = island_merge(w, isl_a, isl_b);
		} else if (isl_a >= 0) {
			target = isl_a;
		} else if (isl_b >= 0) {
			target = isl_b;
		} else {
			target = island_create(w);
		}
		int added_a = 0, added_b = 0;
		if (!is_static_a && w->body_cold[a].island_id != target) {
			island_add_body(w, target, a);
			added_a = 1;
		}
		if (!is_static_b && w->body_cold[b].island_id != target) {
			island_add_body(w, target, b);
			added_b = 1;
		}
		// Wake sleeping island only if a newly added or awake body touches it.
		// Sleeping-vs-static contacts should NOT wake the island.
		if (!w->islands[target].awake) {
			int should_wake = added_a || added_b; // new body joined
			if (!should_wake) {
				// Check if either side was in an awake island before merge
				if (!is_static_a && isl_a >= 0 && isl_a != target) should_wake = 1; // merged from awake
				if (!is_static_b && isl_b >= 0 && isl_b != target) should_wake = 1;
			}
			if (should_wake) island_wake(w, target);
		}
	}

	// Detect lost contacts: iterate prev_touching, find keys not in curr_touching
	if (w->prev_touching) {
		// Walk prev_touching map entries -- use ckit map iteration
		int prev_cap = asize(w->prev_touching);
		for (int i = 0; i < prev_cap; i++) {
			// ckit maps store keys in a parallel array accessible via the header
			// We need to iterate differently -- scan body pairs
		}
		// Simpler approach: we don't have ckit map iteration exposed easily.
		// Instead, mark constraint_remove_count when a contact pair vanishes.
		// We can detect this by checking prev entries against curr.
		// For now, just use a brute-force approach: track via body pair keys.
	}

	// Swap: prev_touching = curr_touching
	map_free(w->prev_touching);
	w->prev_touching = curr_touching;
}

static void island_try_split(WorldInternal* w, int island_id); // forward decl

// Evaluate sleep for all awake islands after post-solve.
static void islands_evaluate_sleep(WorldInternal* w, float dt)
{
	int island_count = asize(w->islands);
	for (int i = 0; i < island_count; i++) {
		if (!(w->island_gen[i] & 1)) continue; // dead slot
		Island* isl = &w->islands[i];
		if (!isl->awake) continue;

		int all_sleepy = 1;
		int bi = isl->head_body;
		while (bi >= 0) {
			BodyHot* h = &w->body_hot[bi];
			float v2 = len2(h->velocity) + len2(h->angular_velocity);
			if (v2 > SLEEP_VEL_THRESHOLD) {
				h->sleep_time = 0.0f;
				all_sleepy = 0;
			} else {
				h->sleep_time += dt;
				if (h->sleep_time < SLEEP_TIME_THRESHOLD)
					all_sleepy = 0;
			}
			bi = w->body_cold[bi].island_next;
		}

		if (all_sleepy) {
			if (isl->constraint_remove_count > 0)
				island_try_split(w, i);
			else
				island_sleep(w, i);
		}
		isl->constraint_remove_count = 0;
	}
}

// DFS island split when contacts have been lost.
static void island_try_split(WorldInternal* w, int island_id)
{
	Island* isl = &w->islands[island_id];

	// Collect all body indices from island
	CK_DYNA int* bodies = NULL;
	int bi = isl->head_body;
	while (bi >= 0) {
		apush(bodies, bi);
		bi = w->body_cold[bi].island_next;
	}
	int n = asize(bodies);
	if (n <= 1) { afree(bodies); island_sleep(w, island_id); return; }

	// Visited array and component assignment
	CK_DYNA int* component = NULL;
	afit(component, n);
	asetlen(component, n);
	for (int i = 0; i < n; i++) component[i] = -1;

	// Map body_idx -> local index for fast lookup
	CK_MAP(int) body_to_local = NULL;
	for (int i = 0; i < n; i++)
		map_set(body_to_local, (uint64_t)bodies[i], i);

	int num_components = 0;
	CK_DYNA int* stack = NULL;

	for (int start = 0; start < n; start++) {
		if (component[start] >= 0) continue;
		int comp_id = num_components++;

		// DFS from start
		aclear(stack);
		apush(stack, start);
		component[start] = comp_id;

		while (asize(stack) > 0) {
			int cur_local = apop(stack);
			int cur_body = bodies[cur_local];

			// Neighbors from joints in this island
			int ji = isl->head_joint;
			while (ji >= 0) {
				JointInternal* j = &w->joints[ji];
				int other = -1;
				if (j->body_a == cur_body) other = j->body_b;
				else if (j->body_b == cur_body) other = j->body_a;
				if (other >= 0) {
					int* loc = map_get_ptr(body_to_local, (uint64_t)other);
					if (loc && component[*loc] < 0) {
						component[*loc] = comp_id;
						apush(stack, *loc);
					}
				}
				ji = w->joints[ji].island_next;
			}

			// Neighbors from current contacts (prev_touching map)
			for (int j = 0; j < n; j++) {
				if (component[j] >= 0) continue;
				uint64_t key = body_pair_key(cur_body, bodies[j]);
				uint8_t* val = map_get_ptr(w->prev_touching, key);
				if (val) {
					component[j] = comp_id;
					apush(stack, j);
				}
			}
		}
	}

	if (num_components <= 1) {
		// Still one island, just sleep it
		island_sleep(w, island_id);
	} else {
		// First component keeps the original island. Remove all bodies/joints,
		// then re-add per component.
		// Detach all bodies from the original island
		for (int i = 0; i < n; i++) {
			w->body_cold[bodies[i]].island_id = -1;
			w->body_cold[bodies[i]].island_prev = -1;
			w->body_cold[bodies[i]].island_next = -1;
		}
		// Detach all joints
		CK_DYNA int* joint_list = NULL;
		int ji = isl->head_joint;
		while (ji >= 0) {
			apush(joint_list, ji);
			int next = w->joints[ji].island_next;
			w->joints[ji].island_id = -1;
			w->joints[ji].island_prev = -1;
			w->joints[ji].island_next = -1;
			ji = next;
		}
		// Reset the original island
		*isl = (Island){ .head_body = -1, .tail_body = -1, .head_joint = -1, .tail_joint = -1, .awake = 1 };

		// Create islands per component (reuse original for component 0)
		CK_DYNA int* comp_islands = NULL;
		afit(comp_islands, num_components);
		asetlen(comp_islands, num_components);
		comp_islands[0] = island_id;
		for (int c = 1; c < num_components; c++)
			comp_islands[c] = island_create(w);

		// Assign bodies to their component's island
		for (int i = 0; i < n; i++)
			island_add_body(w, comp_islands[component[i]], bodies[i]);

		// Assign joints to the island of body_a (both endpoints should be in same component)
		for (int i = 0; i < asize(joint_list); i++) {
			int jidx = joint_list[i];
			int ba = w->joints[jidx].body_a;
			int* loc = map_get_ptr(body_to_local, (uint64_t)ba);
			int comp = loc ? component[*loc] : 0;
			island_add_joint(w, comp_islands[comp], jidx);
		}

		// Evaluate sleep for each sub-island
		for (int c = 0; c < num_components; c++) {
			int cisl = comp_islands[c];
			int all_sleepy = 1;
			int b = w->islands[cisl].head_body;
			while (b >= 0) {
				if (w->body_hot[b].sleep_time < SLEEP_TIME_THRESHOLD) { all_sleepy = 0; break; }
				b = w->body_cold[b].island_next;
			}
			if (all_sleepy) island_sleep(w, cisl);
		}

		afree(comp_islands);
		afree(joint_list);
	}

	afree(stack);
	map_free(body_to_local);
	afree(component);
	afree(bodies);
}

// -----------------------------------------------------------------------------
// World.

World create_world(WorldParams params)
{
	WorldInternal* w = CK_ALLOC(sizeof(WorldInternal));
	memset(w, 0, sizeof(*w));
	w->gravity = params.gravity;
	w->broadphase_type = params.broadphase;
	w->friction_model = params.friction_model;
	w->sleep_enabled = 1;
	w->velocity_iters = params.velocity_iters > 0 ? params.velocity_iters : SOLVER_VELOCITY_ITERS;
	w->position_iters = params.position_iters > 0 ? params.position_iters : SOLVER_POSITION_ITERS;
	w->contact_hertz = params.contact_hertz > 0.0f ? params.contact_hertz : 20.0f;
	w->contact_damping_ratio = params.contact_damping_ratio > 0.0f ? params.contact_damping_ratio : 3.0f;
	w->max_push_velocity = params.max_push_velocity > 0.0f ? params.max_push_velocity : 3.0f;
	w->sub_steps = params.sub_steps > 0 ? params.sub_steps : 4;
	w->bvh_static = CK_ALLOC(sizeof(BVHTree));
	w->bvh_dynamic = CK_ALLOC(sizeof(BVHTree));
	bvh_init(w->bvh_static);
	bvh_init(w->bvh_dynamic);
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
	bvh_free(w->bvh_static); CK_FREE(w->bvh_static);
	bvh_free(w->bvh_dynamic); CK_FREE(w->bvh_dynamic);
	split_free(w->body_cold, w->body_hot, w->body_gen, w->body_free);
	afree(w->joints); afree(w->joint_gen); afree(w->joint_free);
	afree(w->islands); afree(w->island_gen); afree(w->island_free);
	map_free(w->prev_touching);
	CK_FREE(w);
}

static void world_sub_step(WorldInternal* w, float dt)
{
	int count = asize(w->body_hot);

	// Integrate velocities (skip sleeping bodies)
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;
		h->velocity = add(h->velocity, scale(w->gravity, dt));
		// Pade approximation of exponential damping: v *= 1 / (1 + c*dt)
		if (h->linear_damping > 0.0f)
			h->velocity = scale(h->velocity, 1.0f / (1.0f + h->linear_damping * dt));
		if (h->angular_damping > 0.0f)
			h->angular_velocity = scale(h->angular_velocity, 1.0f / (1.0f + h->angular_damping * dt));
	}

	// Collision detection
	CK_DYNA InternalManifold* manifolds = NULL;
	broadphase_and_collide(w, &manifolds);

	// Update island connectivity from contacts
	islands_update_contacts(w, manifolds, asize(manifolds));

	// Store contacts for debug visualization (last sub-step wins)
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
	for (int iter = 0; iter < w->velocity_iters; iter++)
		for (int c = 0; c < color_count; c++)
			for (int i = batch_starts[c]; i < batch_starts[c + 1]; i++)
				solve_constraint(w, &crefs[i], sm, sc, sol_bs, sol_dist);

	afree(crefs);

	// Integrate positions and rotations (skip sleeping bodies)
	for (int i = 0; i < count; i++) {
		if (!split_alive(w->body_gen, i)) continue;
		BodyHot* h = &w->body_hot[i];
		if (h->inv_mass == 0.0f) continue;
		int isl = w->body_cold[i].island_id;
		if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake) continue;

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
		v3 ww = h->angular_velocity;
		quat spin = { ww.x, ww.y, ww.z, 0.0f };
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
	// Skip NGS when soft contacts handle position correction via bias
	if (w->contact_hertz <= 0.0f)
		solver_position_correct(w, sm, asize(sm), sc);

	// Store warm starting data and free solver arrays
	solver_post_solve(w, sm, asize(sm), sc, manifolds, manifold_count);
	joints_post_solve(w, sol_bs, asize(sol_bs), sol_dist, asize(sol_dist));

	afree(manifolds);
}

void world_step(World world, float dt)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int n_sub = w->sub_steps;
	float sub_dt = dt / (float)n_sub;

	// Age warm cache once per frame (not per sub-step)
	warm_cache_age_and_evict(w);

	for (int sub = 0; sub < n_sub; sub++)
		world_sub_step(w, sub_dt);

	// Evaluate sleep once per frame using full dt
	if (w->sleep_enabled) islands_evaluate_sleep(w, dt);
}

void world_set_friction_model(World world, FrictionModel model)
{
	WorldInternal* w = (WorldInternal*)world.id;
	w->friction_model = model;
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
		.bvh_leaf = -1,
		.island_id = -1,
		.island_prev = -1,
		.island_next = -1,
	};
	float fric = params.friction;
	if (fric == 0.0f) fric = 0.5f; // default for all bodies
	float ang_damp = params.angular_damping;
	if (ang_damp == 0.0f) ang_damp = 0.03f; // default: 3%/s (BEPU-style)
	w->body_hot[idx] = (BodyHot){
		.position = params.position,
		.rotation = params.rotation,
		.inv_mass = params.mass > 0.0f ? 1.0f / params.mass : 0.0f,
		.friction = fric,
		.restitution = params.restitution,
		.linear_damping = params.linear_damping,
		.angular_damping = ang_damp,
	};

	return split_handle(Body, w->body_gen, idx);
}

void destroy_body(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	// Remove from island
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl)) {
		// Remove all joints connected to this body
		int ji = w->islands[isl].head_joint;
		while (ji >= 0) {
			int next = w->joints[ji].island_next;
			if (w->joints[ji].body_a == idx || w->joints[ji].body_b == idx) {
				island_remove_joint(w, isl, ji);
				w->islands[isl].constraint_remove_count++;
			}
			ji = next;
		}
		island_remove_body(w, isl, idx);
		w->islands[isl].constraint_remove_count++;
	}
	if (w->body_cold[idx].bvh_leaf >= 0) {
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		bvh_remove(tree, w->body_cold[idx].bvh_leaf);
	}
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

	// Insert into BVH on first shape add.
	if (w->broadphase_type == BROADPHASE_BVH && asize(w->body_cold[idx].shapes) == 1) {
		AABB box = aabb_expand(body_aabb(&w->body_hot[idx], &w->body_cold[idx]), BVH_AABB_MARGIN);
		BVHTree* tree = w->body_hot[idx].inv_mass == 0.0f ? w->bvh_static : w->bvh_dynamic;
		w->body_cold[idx].bvh_leaf = bvh_insert(tree, idx, box);
	}
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

void body_wake(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_velocity(World world, Body body, v3 vel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].velocity = vel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

void body_set_angular_velocity(World world, Body body, v3 avel)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	w->body_hot[idx].angular_velocity = avel;
	int isl = w->body_cold[idx].island_id;
	if (isl >= 0 && island_alive(w, isl) && !w->islands[isl].awake)
		island_wake(w, isl);
}

int body_is_asleep(World world, Body body)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(body);
	assert(split_valid(w->body_gen, body));
	int isl = w->body_cold[idx].island_id;
	if (isl < 0 || !island_alive(w, isl)) return 0;
	return !w->islands[isl].awake;
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
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
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
		.island_id = -1, .island_prev = -1, .island_next = -1,
	};
	link_joint_to_islands(w, idx);
	return (Joint){ handle_make(idx, w->joint_gen[idx]) };
}

void destroy_joint(World world, Joint joint)
{
	WorldInternal* w = (WorldInternal*)world.id;
	int idx = handle_index(joint);
	assert(w->joint_gen[idx] == handle_gen(joint));
	unlink_joint_from_island(w, idx);
	memset(&w->joints[idx], 0, sizeof(JointInternal));
	w->joint_gen[idx]++; // even = dead
	apush(w->joint_free, idx);
}

static void bvh_debug_walk(BVHTree* t, int ni, int depth, BVHDebugFn fn, void* user)
{
	BVHNode* n = &t->nodes[ni];
	for (int s = 0; s < 2; s++) {
		BVHChild* c = bvh_child(n, s);
		if (bvh_child_is_empty(c)) continue;
		fn(c->min, c->max, depth, bvh_child_is_leaf(c), user);
		if (bvh_child_is_internal(c)) bvh_debug_walk(t, c->index, depth + 1, fn, user);
	}
}

void world_debug_bvh(World world, BVHDebugFn fn, void* user)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (w->bvh_dynamic->root >= 0) bvh_debug_walk(w->bvh_dynamic, w->bvh_dynamic->root, 0, fn, user);
	if (w->bvh_static->root >= 0) bvh_debug_walk(w->bvh_static, w->bvh_static->root, 0, fn, user);
}

int world_get_contacts(World world, const Contact** out)
{
	WorldInternal* w = (WorldInternal*)world.id;
	*out = w->debug_contacts;
	return asize(w->debug_contacts);
}
