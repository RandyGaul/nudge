// contacts.c -- ContactSummary build + sort + dedupe + per-body listener dispatch.
//
// Runs once per world_step immediately after narrowphase. Input is the flat
// manifolds[] array (already compacted so every entry has count > 0). Output:
//   - w->contact_summaries: sorted by (a.id, b.id), canonical a.id < b.id,
//     deduped so multi-triangle mesh pairs collapse into one summary.
//   - body listener callbacks fire once per body that owns contacts this step,
//     with a contiguous view where the listener's body is always in the `a`
//     slot (normal flipped, material/sub fields swapped as needed).
//
// One-manifold-per-contact-pair is the common case; only mesh pairs (and a
// future multi-shape body pair case) produce multiple manifolds per pair,
// and the merge step reduces them.

// Decode the triangle index (0-based) from the internal InternalManifold.sub_id,
// which packs ((sa+1) << 16) | (tri_idx + 1) for convex-vs-mesh pairs. Returns
// 0 if the body is not a mesh side.
static inline uint32_t decode_tri_index_plus_one(const InternalManifold* im, int body_is_mesh_side)
{
	if (!body_is_mesh_side) return 0;
	return im->sub_id & 0xFFFFu;
}

// Resolve the material id for the mesh side of a pair. tri_plus_one is the
// triangle index + 1 as encoded in the sub_id low half (0 means no mesh).
static uint8_t resolve_mesh_material(WorldInternal* w, int body_idx, uint32_t tri_plus_one)
{
	uint8_t default_id = w->body_cold[body_idx].material_id;
	if (tri_plus_one == 0) return default_id;
	ShapeInternal* shapes = w->body_cold[body_idx].shapes;
	if (!shapes || shapes[0].type != SHAPE_MESH) return default_id;
	const TriMesh* mesh = shapes[0].mesh.mesh;
	if (!mesh->material_ids) return default_id;
	return mesh->material_ids[tri_plus_one - 1];
}

// Reduce a Manifold to a single ContactSummary row. Canonicalises a.id < b.id.
// Fills material_a/b and sub_a/b with defaults matching side-of-body (mesh is
// always B in an InternalManifold; sub_id carries tri_idx + 1).
static void contact_row_from_manifold(WorldInternal* w, const InternalManifold* im, ContactSummary* out)
{
	const Manifold* m = &im->m;
	int deepest = 0;
	float max_pen = m->contacts[0].penetration;
	v3 centroid = m->contacts[0].point;
	for (int i = 1; i < m->count; i++) {
		centroid = add(centroid, m->contacts[i].point);
		if (m->contacts[i].penetration > max_pen) { max_pen = m->contacts[i].penetration; deepest = i; }
	}
	float inv = 1.0f / (float)m->count;
	centroid = scale(centroid, inv);
	float rad2 = 0;
	for (int i = 0; i < m->count; i++) {
		float d2 = v3_len2(sub(m->contacts[i].point, centroid));
		if (d2 > rad2) rad2 = d2;
	}
	Body ha = split_handle(Body, w->body_gen, im->body_a);
	Body hb = split_handle(Body, w->body_gen, im->body_b);
	v3 normal = m->contacts[deepest].normal;
	// Mesh is always B in the internal manifold (narrowphase canonicalises
	// mesh-to-B before emit). sub_b carries the triangle index + 1 in its
	// low 16 bits; sub_a is 0 for non-mesh sides.
	int b_is_mesh = im->sub_id != 0 && w->body_cold[im->body_b].shapes
	              && w->body_cold[im->body_b].shapes[0].type == SHAPE_MESH;
	uint32_t sub_a = 0;
	uint32_t sub_b = decode_tri_index_plus_one(im, b_is_mesh);
	uint8_t mat_a = w->body_cold[im->body_a].material_id;
	uint8_t mat_b = resolve_mesh_material(w, im->body_b, sub_b);
	if (ha.id > hb.id) {
		Body t = ha; ha = hb; hb = t;
		normal = neg(normal);
		uint32_t ts = sub_a; sub_a = sub_b; sub_b = ts;
		uint8_t tm = mat_a; mat_a = mat_b; mat_b = tm;
	}
	out->a = ha;
	out->b = hb;
	out->normal = normal;
	out->point = centroid;
	out->radius = sqrtf(rad2);
	out->depth = max_pen;
	out->material_a = mat_a;
	out->material_b = mat_b;
	out->_pad = 0;
	out->sub_a = sub_a;
	out->sub_b = sub_b;
}

// Merge src into dst (already has one pair's data). Both must share the same
// (a, b). Used when a body pair produced multiple manifolds (mesh, future
// multi-shape). The merged patch covers the union: centroid is the weighted
// mean, radius encloses both patches, depth is the max, normal is taken from
// whichever input was deeper (single-normal reduction).
static void contact_row_merge(ContactSummary* dst, const ContactSummary* src)
{
	v3 d = sub(src->point, dst->point);
	float dlen = v3_len(d);
	// Weighted centroid (by depth, floored to epsilon so degenerate zero-depth
	// pairs still contribute). Approximation -- we've already lost contact
	// counts by this point.
	float wa = dst->depth > 1e-6f ? dst->depth : 1e-6f;
	float wb = src->depth > 1e-6f ? src->depth : 1e-6f;
	float inv = 1.0f / (wa + wb);
	v3 new_centroid = scale(add(scale(dst->point, wa), scale(src->point, wb)), inv);
	// Radius bound: enclose both old patch centres under the new centre.
	float ra = dst->radius + v3_len(sub(dst->point, new_centroid));
	float rb = src->radius + v3_len(sub(src->point, new_centroid));
	(void)dlen;
	dst->point = new_centroid;
	dst->radius = ra > rb ? ra : rb;
	if (src->depth > dst->depth) {
		dst->depth = src->depth;
		dst->normal = src->normal;
		dst->material_a = src->material_a;
		dst->material_b = src->material_b;
		dst->sub_a = src->sub_a;
		dst->sub_b = src->sub_b;
	}
}

static int contact_summary_cmp(const void* a, const void* b)
{
	const ContactSummary* x = (const ContactSummary*)a;
	const ContactSummary* y = (const ContactSummary*)b;
	if (x->a.id < y->a.id) return -1;
	if (x->a.id > y->a.id) return  1;
	if (x->b.id < y->b.id) return -1;
	if (x->b.id > y->b.id) return  1;
	return 0;
}

static void contacts_build(WorldInternal* w, InternalManifold* manifolds, int mcount)
{
	aclear(w->contact_summaries);
	if (mcount == 0) return;
	afit(w->contact_summaries, mcount);
	asetlen(w->contact_summaries, mcount);
	for (int i = 0; i < mcount; i++) {
		contact_row_from_manifold(w, &manifolds[i], &w->contact_summaries[i]);
	}
	qsort(w->contact_summaries, (size_t)mcount, sizeof(ContactSummary), contact_summary_cmp);
	// Dedupe: merge consecutive rows with matching (a, b).
	int write = 0;
	for (int read = 0; read < mcount; read++) {
		if (write > 0 && w->contact_summaries[write - 1].a.id == w->contact_summaries[read].a.id
		              && w->contact_summaries[write - 1].b.id == w->contact_summaries[read].b.id) {
			contact_row_merge(&w->contact_summaries[write - 1], &w->contact_summaries[read]);
		} else {
			if (write != read) w->contact_summaries[write] = w->contact_summaries[read];
			write++;
		}
	}
	asetlen(w->contact_summaries, write);
}

// Flip a summary so `self` is in the `a` slot. Assumes summary already has
// `self` in either slot. Idempotent if already oriented.
static inline void contact_row_orient_for(ContactSummary* row, Body self)
{
	if (row->a.id == self.id) return;
	Body t = row->a; row->a = row->b; row->b = t;
	row->normal = neg(row->normal);
	uint32_t ts = row->sub_a; row->sub_a = row->sub_b; row->sub_b = ts;
	uint8_t tm = row->material_a; row->material_a = row->material_b; row->material_b = tm;
}

static void contacts_dispatch_listeners(WorldInternal* w)
{
	if (!w->body_listeners || map_size(w->body_listeners) == 0) return;
	int n = asize(w->contact_summaries);
	if (n == 0) return;
	// For each listening body, gather its summaries into the scratch buffer
	// with self in the `a` slot, then fire once.
	//
	// Two-pass O(listeners * pairs). Typical scenes have few listeners and
	// bounded contacts; a reverse index is only worth adding if listener
	// count grows large.
	map_each(w->body_listeners, mi) {
		int body_idx = (int)map_key(w->body_listeners, mi);
		if (!split_alive(w->body_gen, body_idx)) continue;
		BodyListener entry = w->body_listeners[mi];
		if (!entry.fn) continue;
		Body self = split_handle(Body, w->body_gen, body_idx);
		aclear(w->listener_scratch);
		for (int i = 0; i < n; i++) {
			const ContactSummary* s = &w->contact_summaries[i];
			if (s->a.id == self.id || s->b.id == self.id) {
				ContactSummary row = *s;
				contact_row_orient_for(&row, self);
				apush(w->listener_scratch, row);
			}
		}
		int count = asize(w->listener_scratch);
		if (count > 0) entry.fn(self, w->listener_scratch, count, entry.ud);
	}
}

void body_set_contact_listener(World world, Body body, BodyContactListener fn, void* ud)
{
	WorldInternal* w = (WorldInternal*)world.id;
	assert(split_valid(w->body_gen, body));
	int idx = handle_index(body);
	if (fn == NULL) {
		if (w->body_listeners) map_del(w->body_listeners, (uint64_t)idx);
		return;
	}
	BodyListener entry = { .fn = fn, .ud = ud };
	map_set(w->body_listeners, (uint64_t)idx, entry);
}

const ContactSummary* world_contact_summaries(World world, int* out_count)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (out_count) *out_count = asize(w->contact_summaries);
	return w->contact_summaries;
}
