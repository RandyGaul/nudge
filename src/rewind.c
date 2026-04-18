// See LICENSE for licensing info.
// rewind.c -- deterministic world state snapshots.
//
// v1: full snapshots; topology-change invalidation.
// v2: mutation replay (shapes deep-cloned; ring survives create/destroy).
// v3: lossless delta compression. Per frame, we store ONLY bodies whose
//     (BodyHot, BodyState, BodyCold, shapes[]) bytes changed since the
//     previous capture. Then the dirty payload is XOR'd against the prior
//     frame's same bytes and zero-run encoded; this trims near-zero-delta
//     bodies (e.g. bodies that moved an epsilon between substeps) down to
//     a handful of bytes.
//
// Keyframes (first frame, or after body count changes) store the raw
// payload uncompressed. Delta frames store XOR-RLE. On eviction, the
// next-oldest frame is promoted to a keyframe by folding in the evictee's
// data for any bodies the delta didn't cover.
//
// Determinism: every encoding step is lossless. Restore reconstructs bytes
// identical to the captured state, so stepping from a restored snapshot
// produces the same output as the original run.
//
// Captured but not compressed: joints, islands, warm_cache, prev_touching,
// joint_pairs, body_gen, body_free — stored full per frame.
//
// NOT captured: BVH trees (rebuilt on restore), LDL caches (zeroed on
// restore), EPA manifold cache (use SAT backend), perf timers, worker
// arenas.

typedef struct RewindIsland
{
	int head_body, tail_body, body_count;
	int head_joint, tail_joint, joint_count;
	int constraint_remove_count;
	int awake;
} RewindIsland;

// Per-dirty-body payload: hot+state+cold packed into a contiguous byte
// block that's XOR-RLE encoded in delta frames. Shapes travel alongside as
// a deep-owned CK_DYNA per body (shapes are compared byte-wise but not
// XOR-RLE encoded; they're typically tiny and rarely change).
#define REWIND_PAYLOAD_SIZE (sizeof(BodyHot) + sizeof(BodyState) + sizeof(BodyCold))

typedef struct RewindFrame
{
	uint64_t frame_id;
	int sim_frame;

	int n_bodies;
	int is_keyframe;   // 1 = payload stored raw; 0 = XOR-RLE against prior frame
	int n_dirty;
	int* dirty_indices;          // [n_dirty] body indices, sorted ascending

	// Dirty bodies' (hot, state, cold) packed and encoded. Keyframe: raw bytes
	// of REWIND_PAYLOAD_SIZE × n_dirty. Delta frame: variable-length RLE
	// stream per body, indexed by dirty_payload_offsets.
	uint8_t* dirty_payload;
	int dirty_payload_size;
	int* dirty_payload_offsets;  // [n_dirty+1] for delta frames; NULL for keyframe (stride is implicit)

	// Deep-owned per-dirty-body shapes arrays (NULL if that body has no shapes).
	ShapeInternal** dirty_shapes; // [n_dirty]

	uint32_t* body_gen;           // [n_bodies] full copy
	CK_DYNA int* body_free;

	int n_joints;
	JointInternal* joints;
	JointHot*      joint_hot;
	uint32_t*      joint_gen;
	CK_DYNA int*   joint_free;

	int n_islands;
	RewindIsland*  islands;
	uint32_t*      island_gen;
	CK_DYNA int*   island_free;

	int n_warm;
	uint64_t*      warm_keys;
	WarmManifold*  warm_vals;

	int n_prev_touching;
	uint64_t*      prev_touching_keys;
	uint8_t*       prev_touching_vals;

	int n_joint_pairs;
	uint64_t*      joint_pairs_keys;
	uint8_t*       joint_pairs_vals;

	int joint_pairs_version;
	int ldl_topo_version;

	// Sensors: captured as a full list per frame (sensors are few and rarely
	// change, so dirty compression doesn't buy much). Each entry owns a
	// deep-cloned shapes array.
	int n_sensors;
	SensorInternal* sensors;   // [n_sensors]; .shapes is a CK_DYNA owned here
	uint32_t*       sensor_gen;
	CK_DYNA int*    sensor_free;

	// Material palette (256 entries). Stored copy-on-write: keyframes always
	// carry a full copy; delta frames only store when the palette changed vs
	// the rolling baseline. Restore walks backward to the nearest frame with
	// a non-NULL pointer. NULL = inherit from earlier frame. Per-body
	// material_id rides in BodyCold and is already covered by the dirty path.
	Material*      materials;
} RewindFrame;

typedef struct RewindBuffer
{
	int max_frames;
	int auto_capture;
	uint64_t next_frame_id;
	CK_DYNA RewindFrame* frames;
	int head;
	int count;

	// Baseline: full world state as of the previous capture. Used to detect
	// dirty bodies and as the XOR reference when encoding delta frames.
	int baseline_n_bodies;
	BodyHot*   baseline_hot;
	BodyState* baseline_state;
	BodyCold*  baseline_cold;
	ShapeInternal** baseline_shapes; // [n_bodies] deep-owned
	// Material palette baseline (last-captured snapshot). A frame stores a
	// palette copy only when this differs from the world's current palette.
	Material baseline_materials[256];
} RewindBuffer;

static void rewind_frame_free(RewindFrame* f)
{
	CK_FREE(f->dirty_indices);         f->dirty_indices = NULL;
	CK_FREE(f->dirty_payload);         f->dirty_payload = NULL;
	CK_FREE(f->dirty_payload_offsets); f->dirty_payload_offsets = NULL;
	f->dirty_payload_size = 0;
	if (f->dirty_shapes) {
		for (int k = 0; k < f->n_dirty; k++) afree(f->dirty_shapes[k]);
		CK_FREE(f->dirty_shapes);
		f->dirty_shapes = NULL;
	}
	f->n_dirty = 0;

	CK_FREE(f->body_gen);   f->body_gen = NULL;
	afree(f->body_free);    f->body_free = NULL;

	CK_FREE(f->joints);     f->joints = NULL;
	CK_FREE(f->joint_hot);  f->joint_hot = NULL;
	CK_FREE(f->joint_gen);  f->joint_gen = NULL;
	afree(f->joint_free);   f->joint_free = NULL;

	CK_FREE(f->islands);    f->islands = NULL;
	CK_FREE(f->island_gen); f->island_gen = NULL;
	afree(f->island_free);  f->island_free = NULL;

	CK_FREE(f->warm_keys);  f->warm_keys = NULL;
	CK_FREE(f->warm_vals);  f->warm_vals = NULL;

	CK_FREE(f->prev_touching_keys); f->prev_touching_keys = NULL;
	CK_FREE(f->prev_touching_vals); f->prev_touching_vals = NULL;

	CK_FREE(f->joint_pairs_keys);   f->joint_pairs_keys = NULL;
	CK_FREE(f->joint_pairs_vals);   f->joint_pairs_vals = NULL;

	if (f->sensors) {
		for (int i = 0; i < f->n_sensors; i++) afree(f->sensors[i].shapes);
		CK_FREE(f->sensors);
		f->sensors = NULL;
	}
	CK_FREE(f->sensor_gen); f->sensor_gen = NULL;
	afree(f->sensor_free);  f->sensor_free = NULL;
	f->n_sensors = 0;

	CK_FREE(f->materials); f->materials = NULL;
}

static void* rewind_clone_array(const void* src, size_t n, size_t elem_size)
{
	if (n == 0) return NULL;
	void* dst = CK_ALLOC(n * elem_size);
	memcpy(dst, src, n * elem_size);
	return dst;
}

static int* rewind_clone_dyna_int(int* src)
{
	int n = asize(src);
	int* dst = NULL;
	if (n > 0) { afit(dst, n); asetlen(dst, n); memcpy(dst, src, (size_t)n * sizeof(int)); }
	return dst;
}

static ShapeInternal* rewind_clone_shapes(ShapeInternal* src)
{
	int n = asize(src);
	if (n == 0) return NULL;
	ShapeInternal* dst = NULL;
	afit(dst, n); asetlen(dst, n);
	memcpy(dst, src, (size_t)n * sizeof(ShapeInternal));
	return dst;
}

static int rewind_shapes_equal(const ShapeInternal* a, const ShapeInternal* b)
{
	int na = asize((ShapeInternal*)a);
	int nb = asize((ShapeInternal*)b);
	if (na != nb) return 0;
	if (na == 0) return 1;
	return memcmp(a, b, (size_t)na * sizeof(ShapeInternal)) == 0;
}

// ----- XOR-RLE encode / decode --------------------------------------------
//
// Delta payload format, per body:
//   repeat until input consumed:
//     zero_run: 1 byte, 0..255 zero bytes
//     data_run: 1 byte, 0..255 nonzero bytes (may exceed 255 across tuples)
//     data_run bytes of XOR'd data
// Terminator: payload length is known (from offsets), so no sentinel needed.
//
// Encoder packs a fixed REWIND_PAYLOAD_SIZE-byte input (body hot+state+cold
// XOR'd with baseline). Decoder reverses.

static void rewind_pack_body(uint8_t* out, const BodyHot* h, const BodyState* s, const BodyCold* c)
{
	memcpy(out + 0,                                              h, sizeof(BodyHot));
	memcpy(out + sizeof(BodyHot),                                s, sizeof(BodyState));
	memcpy(out + sizeof(BodyHot) + sizeof(BodyState),            c, sizeof(BodyCold));
}

static void rewind_unpack_body(const uint8_t* in, BodyHot* h, BodyState* s, BodyCold* c)
{
	memcpy(h, in + 0,                                              sizeof(BodyHot));
	memcpy(s, in + sizeof(BodyHot),                                sizeof(BodyState));
	memcpy(c, in + sizeof(BodyHot) + sizeof(BodyState),            sizeof(BodyCold));
}

// Encode `len`-byte XOR'd buffer into a dynamically-grown CK_DYNA byte array.
static void rewind_rle_encode(uint8_t** out, const uint8_t* in, int len)
{
	int pos = 0;
	while (pos < len) {
		// Count leading zeros (cap at 255).
		int zr = 0;
		while (pos < len && zr < 255 && in[pos] == 0) { zr++; pos++; }
		// Count leading non-zeros (cap at 255).
		int dr = 0;
		while (pos + dr < len && dr < 255 && in[pos + dr] != 0) dr++;
		apush(*out, (uint8_t)zr);
		apush(*out, (uint8_t)dr);
		for (int i = 0; i < dr; i++) apush(*out, in[pos + i]);
		pos += dr;
	}
}

static void rewind_rle_decode(uint8_t* out_xor, int out_len, const uint8_t* in, int in_len)
{
	int op = 0, ip = 0;
	memset(out_xor, 0, (size_t)out_len);
	while (ip < in_len) {
		int zr = in[ip++];
		int dr = in[ip++];
		op += zr;  // zeros already there from memset
		for (int i = 0; i < dr && op < out_len; i++) out_xor[op++] = in[ip++];
	}
}

// ----- Capture -------------------------------------------------------------

static int rewind_body_differs(const BodyHot* h, const BodyState* s, const BodyCold* c,
                               ShapeInternal* shapes,
                               const BodyHot* bh, const BodyState* bs, const BodyCold* bc,
                               ShapeInternal* bshapes)
{
	if (memcmp(h, bh, sizeof(BodyHot)) != 0)   return 1;
	if (memcmp(s, bs, sizeof(BodyState)) != 0) return 1;
	if (memcmp(c, bc, sizeof(BodyCold)) != 0)  return 1;
	if (!rewind_shapes_equal(shapes, bshapes)) return 1;
	return 0;
}

static void rewind_capture_into(RewindBuffer* rb, RewindFrame* f, WorldInternal* w)
{
	rewind_frame_free(f);

	f->frame_id = ++rb->next_frame_id;
	f->sim_frame = w->frame;

	int nb = asize(w->body_hot);
	f->n_bodies = nb;

	// Keyframe when there's no predecessor or body-count changed.
	int is_keyframe = (rb->count == 0) || (nb != rb->baseline_n_bodies);
	f->is_keyframe = is_keyframe;

	// Determine dirty bodies. For keyframes, all bodies are "dirty" (raw
	// payload storage). For deltas, any field differing from baseline qualifies.
	// We use a snapshot-version of BodyCold with .shapes zeroed so that the
	// live heap pointer never leaks into baseline/frame state.
	int dirty_count = 0;
	int* dirty_tmp = (int*)CK_ALLOC((size_t)(nb > 0 ? nb : 1) * sizeof(int));
	for (int i = 0; i < nb; i++) {
		if (is_keyframe) {
			dirty_tmp[dirty_count++] = i;
			continue;
		}
		BodyCold cur_cold = w->body_cold[i];
		cur_cold.shapes = NULL;
		if (rewind_body_differs(&w->body_hot[i], &w->body_state[i], &cur_cold, w->body_cold[i].shapes,
		                        &rb->baseline_hot[i], &rb->baseline_state[i], &rb->baseline_cold[i], rb->baseline_shapes[i])) {
			dirty_tmp[dirty_count++] = i;
		}
	}

	f->n_dirty = dirty_count;
	if (dirty_count > 0) {
		f->dirty_indices = (int*)CK_ALLOC((size_t)dirty_count * sizeof(int));
		f->dirty_shapes  = (ShapeInternal**)CK_ALLOC((size_t)dirty_count * sizeof(ShapeInternal*));
		memcpy(f->dirty_indices, dirty_tmp, (size_t)dirty_count * sizeof(int));

		uint8_t scratch[REWIND_PAYLOAD_SIZE];
		uint8_t xorbuf[REWIND_PAYLOAD_SIZE];

		if (is_keyframe) {
			// Raw payload: dirty_payload = [payload 0, payload 1, ..., payload n_dirty-1]
			f->dirty_payload_size = dirty_count * (int)REWIND_PAYLOAD_SIZE;
			f->dirty_payload = (uint8_t*)CK_ALLOC((size_t)f->dirty_payload_size);
			f->dirty_payload_offsets = NULL;  // implicit stride REWIND_PAYLOAD_SIZE
			for (int k = 0; k < dirty_count; k++) {
				int i = dirty_tmp[k];
				BodyCold cold_snap = w->body_cold[i];
				cold_snap.shapes = NULL;
				rewind_pack_body(f->dirty_payload + k * REWIND_PAYLOAD_SIZE,
				                 &w->body_hot[i], &w->body_state[i], &cold_snap);
				f->dirty_shapes[k] = rewind_clone_shapes(w->body_cold[i].shapes);
			}
		} else {
			// Delta: XOR-RLE each body's payload against baseline.
			CK_DYNA uint8_t* stream = NULL;
			int* offsets = (int*)CK_ALLOC((size_t)(dirty_count + 1) * sizeof(int));
			offsets[0] = 0;
			for (int k = 0; k < dirty_count; k++) {
				int i = dirty_tmp[k];
				BodyCold cold_snap = w->body_cold[i];
				cold_snap.shapes = NULL;
				rewind_pack_body(scratch, &w->body_hot[i], &w->body_state[i], &cold_snap);
				uint8_t baseline_buf[REWIND_PAYLOAD_SIZE];
				rewind_pack_body(baseline_buf, &rb->baseline_hot[i], &rb->baseline_state[i], &rb->baseline_cold[i]);
				for (size_t b = 0; b < REWIND_PAYLOAD_SIZE; b++) xorbuf[b] = scratch[b] ^ baseline_buf[b];
				rewind_rle_encode(&stream, xorbuf, (int)REWIND_PAYLOAD_SIZE);
				offsets[k + 1] = asize(stream);
				f->dirty_shapes[k] = rewind_clone_shapes(w->body_cold[i].shapes);
			}
			f->dirty_payload_size = asize(stream);
			if (f->dirty_payload_size > 0) {
				f->dirty_payload = (uint8_t*)CK_ALLOC((size_t)f->dirty_payload_size);
				memcpy(f->dirty_payload, stream, (size_t)f->dirty_payload_size);
			}
			afree(stream);
			f->dirty_payload_offsets = offsets;
		}
	}
	CK_FREE(dirty_tmp);

	// Body gen / free — full copy, cheap.
	f->body_gen  = rewind_clone_array(w->body_gen,  nb, sizeof(uint32_t));
	f->body_free = rewind_clone_dyna_int(w->body_free);

	int nj = asize(w->joints);
	f->n_joints = nj;
	f->joints     = rewind_clone_array(w->joints,    nj, sizeof(JointInternal));
	f->joint_hot  = rewind_clone_array(w->joint_hot, nj, sizeof(JointHot));
	f->joint_gen  = rewind_clone_array(w->joint_gen, nj, sizeof(uint32_t));
	f->joint_free = rewind_clone_dyna_int(w->joint_free);

	int ni = asize(w->islands);
	f->n_islands = ni;
	if (ni > 0) {
		f->islands = (RewindIsland*)CK_ALLOC((size_t)ni * sizeof(RewindIsland));
		for (int i = 0; i < ni; i++) {
			f->islands[i].head_body             = w->islands[i].head_body;
			f->islands[i].tail_body             = w->islands[i].tail_body;
			f->islands[i].body_count            = w->islands[i].body_count;
			f->islands[i].head_joint            = w->islands[i].head_joint;
			f->islands[i].tail_joint            = w->islands[i].tail_joint;
			f->islands[i].joint_count           = w->islands[i].joint_count;
			f->islands[i].constraint_remove_count = w->islands[i].constraint_remove_count;
			f->islands[i].awake                 = w->islands[i].awake;
		}
	}
	f->island_gen  = rewind_clone_array(w->island_gen, ni, sizeof(uint32_t));
	f->island_free = rewind_clone_dyna_int(w->island_free);

	int nw = map_size(w->warm_cache);
	f->n_warm = nw;
	if (nw > 0) {
		f->warm_keys = (uint64_t*)CK_ALLOC((size_t)nw * sizeof(uint64_t));
		f->warm_vals = (WarmManifold*)CK_ALLOC((size_t)nw * sizeof(WarmManifold));
		uint64_t* keys = map_keys(w->warm_cache);
		for (int i = 0; i < nw; i++) {
			f->warm_keys[i] = keys[i];
			f->warm_vals[i] = w->warm_cache[i];
		}
	}

	int npt = map_size(w->prev_touching);
	f->n_prev_touching = npt;
	if (npt > 0) {
		f->prev_touching_keys = (uint64_t*)CK_ALLOC((size_t)npt * sizeof(uint64_t));
		f->prev_touching_vals = (uint8_t*)CK_ALLOC((size_t)npt * sizeof(uint8_t));
		uint64_t* keys = map_keys(w->prev_touching);
		for (int i = 0; i < npt; i++) {
			f->prev_touching_keys[i] = keys[i];
			f->prev_touching_vals[i] = w->prev_touching[i];
		}
	}

	int njp = map_size(w->joint_pairs);
	f->n_joint_pairs = njp;
	if (njp > 0) {
		f->joint_pairs_keys = (uint64_t*)CK_ALLOC((size_t)njp * sizeof(uint64_t));
		f->joint_pairs_vals = (uint8_t*)CK_ALLOC((size_t)njp * sizeof(uint8_t));
		uint64_t* keys = map_keys(w->joint_pairs);
		for (int i = 0; i < njp; i++) {
			f->joint_pairs_keys[i] = keys[i];
			f->joint_pairs_vals[i] = w->joint_pairs[i];
		}
	}

	f->joint_pairs_version = w->joint_pairs_version;
	f->ldl_topo_version    = w->ldl_topo_version;

	// Sensors: full-copy per frame. Shapes are deep-cloned so destroy_sensor
	// between captures doesn't invalidate prior snapshots.
	int ns = asize(w->sensors);
	f->n_sensors = ns;
	if (ns > 0) {
		f->sensors = (SensorInternal*)CK_ALLOC((size_t)ns * sizeof(SensorInternal));
		for (int i = 0; i < ns; i++) {
			f->sensors[i] = w->sensors[i];
			f->sensors[i].shapes = NULL;
			int nshapes = asize(w->sensors[i].shapes);
			if (nshapes > 0) {
				ShapeInternal* dst = NULL;
				afit(dst, nshapes); asetlen(dst, nshapes);
				memcpy(dst, w->sensors[i].shapes, (size_t)nshapes * sizeof(ShapeInternal));
				f->sensors[i].shapes = dst;
			}
		}
	}
	f->sensor_gen  = rewind_clone_array(w->sensor_gen,  ns, sizeof(uint32_t));
	f->sensor_free = rewind_clone_dyna_int(w->sensor_free);

	// Material palette: keyframes always carry a copy so restore never has
	// to walk past the ring's oldest frame. Delta frames only store when
	// the palette changed since the last capture (detected against baseline).
	int materials_changed = memcmp(w->materials, rb->baseline_materials, sizeof(w->materials)) != 0;
	if (is_keyframe || materials_changed) {
		f->materials = (Material*)CK_ALLOC(sizeof(w->materials));
		memcpy(f->materials, w->materials, sizeof(w->materials));
	}
	if (materials_changed) memcpy(rb->baseline_materials, w->materials, sizeof(w->materials));

	// Update baseline to match current world state.
	if (nb != rb->baseline_n_bodies) {
		CK_FREE(rb->baseline_hot);   rb->baseline_hot = NULL;
		CK_FREE(rb->baseline_state); rb->baseline_state = NULL;
		CK_FREE(rb->baseline_cold);  rb->baseline_cold = NULL;
		if (rb->baseline_shapes) {
			for (int i = 0; i < rb->baseline_n_bodies; i++) afree(rb->baseline_shapes[i]);
			CK_FREE(rb->baseline_shapes);
			rb->baseline_shapes = NULL;
		}
		if (nb > 0) {
			rb->baseline_hot    = (BodyHot*)CK_ALLOC((size_t)nb * sizeof(BodyHot));
			rb->baseline_state  = (BodyState*)CK_ALLOC((size_t)nb * sizeof(BodyState));
			rb->baseline_cold   = (BodyCold*)CK_ALLOC((size_t)nb * sizeof(BodyCold));
			rb->baseline_shapes = (ShapeInternal**)CK_ALLOC((size_t)nb * sizeof(ShapeInternal*));
			for (int i = 0; i < nb; i++) rb->baseline_shapes[i] = NULL;
		}
		rb->baseline_n_bodies = nb;
	}
	if (nb > 0) {
		memcpy(rb->baseline_hot,   w->body_hot,   (size_t)nb * sizeof(BodyHot));
		memcpy(rb->baseline_state, w->body_state, (size_t)nb * sizeof(BodyState));
		for (int i = 0; i < nb; i++) {
			rb->baseline_cold[i] = w->body_cold[i];
			rb->baseline_cold[i].shapes = NULL;
			// Refresh baseline shape copy only when it actually changed, to
			// avoid churning allocations on every capture.
			if (!rewind_shapes_equal(w->body_cold[i].shapes, rb->baseline_shapes[i])) {
				afree(rb->baseline_shapes[i]);
				rb->baseline_shapes[i] = rewind_clone_shapes(w->body_cold[i].shapes);
			}
		}
	}
}

// ----- Decode helper: read body k's (hot, state, cold) from frame ---------
//
// For keyframe: raw read from payload offset k * REWIND_PAYLOAD_SIZE.
// For delta:    RLE-decode against `baseline_payload` then XOR to reconstruct.
static void rewind_decode_body(const RewindFrame* f, int k,
                               const uint8_t* baseline_payload,
                               BodyHot* out_hot, BodyState* out_state, BodyCold* out_cold)
{
	uint8_t buf[REWIND_PAYLOAD_SIZE];
	if (f->is_keyframe) {
		memcpy(buf, f->dirty_payload + (size_t)k * REWIND_PAYLOAD_SIZE, REWIND_PAYLOAD_SIZE);
	} else {
		int start = f->dirty_payload_offsets[k];
		int end   = f->dirty_payload_offsets[k + 1];
		uint8_t xorbuf[REWIND_PAYLOAD_SIZE];
		rewind_rle_decode(xorbuf, (int)REWIND_PAYLOAD_SIZE, f->dirty_payload + start, end - start);
		for (size_t b = 0; b < REWIND_PAYLOAD_SIZE; b++) buf[b] = xorbuf[b] ^ baseline_payload[b];
	}
	rewind_unpack_body(buf, out_hot, out_state, out_cold);
}

// ----- Eviction + promotion ------------------------------------------------

// Promote a delta frame L to a keyframe by filling in bodies that L didn't
// carry, using the predecessor K (which is a keyframe). After promotion L
// stores raw payload covering all bodies; K can be safely evicted.
static void rewind_promote_to_keyframe(RewindFrame* L, const RewindFrame* K)
{
	if (L->is_keyframe) return;
	int nb = L->n_bodies;

	// Materialize L's full (hot, state, cold, shapes) into scratch using K
	// (keyframe) and L's own deltas.
	BodyHot*   full_hot   = (BodyHot*)CK_ALLOC((size_t)nb * sizeof(BodyHot));
	BodyState* full_state = (BodyState*)CK_ALLOC((size_t)nb * sizeof(BodyState));
	BodyCold*  full_cold  = (BodyCold*)CK_ALLOC((size_t)nb * sizeof(BodyCold));
	ShapeInternal** full_shapes = (ShapeInternal**)CK_ALLOC((size_t)nb * sizeof(ShapeInternal*));
	for (int i = 0; i < nb; i++) full_shapes[i] = NULL;

	// Step 1: fill from K. K is keyframe, so K.dirty_indices covers [0, K->n_bodies).
	for (int k = 0; k < K->n_dirty; k++) {
		int i = K->dirty_indices[k];
		if (i < 0 || i >= nb) continue;
		rewind_unpack_body(K->dirty_payload + (size_t)k * REWIND_PAYLOAD_SIZE,
		                   &full_hot[i], &full_state[i], &full_cold[i]);
		full_shapes[i] = rewind_clone_shapes(K->dirty_shapes[k]);
	}

	// Step 2: overlay L's deltas. Baseline for L's XOR is K's corresponding body data.
	uint8_t baseline_buf[REWIND_PAYLOAD_SIZE];
	for (int k = 0; k < L->n_dirty; k++) {
		int i = L->dirty_indices[k];
		rewind_pack_body(baseline_buf, &full_hot[i], &full_state[i], &full_cold[i]);
		rewind_decode_body(L, k, baseline_buf, &full_hot[i], &full_state[i], &full_cold[i]);
		afree(full_shapes[i]);
		full_shapes[i] = rewind_clone_shapes(L->dirty_shapes[k]);
	}

	// Step 3: rewrite L as a keyframe with raw payload for all bodies.
	CK_FREE(L->dirty_indices);
	CK_FREE(L->dirty_payload);
	CK_FREE(L->dirty_payload_offsets);
	if (L->dirty_shapes) {
		for (int k = 0; k < L->n_dirty; k++) afree(L->dirty_shapes[k]);
		CK_FREE(L->dirty_shapes);
	}

	L->is_keyframe = 1;
	L->n_dirty = nb;
	L->dirty_payload_offsets = NULL;
	L->dirty_payload_size = nb * (int)REWIND_PAYLOAD_SIZE;
	L->dirty_indices = (int*)CK_ALLOC((size_t)nb * sizeof(int));
	L->dirty_payload = (uint8_t*)CK_ALLOC((size_t)L->dirty_payload_size);
	L->dirty_shapes  = (ShapeInternal**)CK_ALLOC((size_t)nb * sizeof(ShapeInternal*));
	for (int i = 0; i < nb; i++) {
		L->dirty_indices[i] = i;
		rewind_pack_body(L->dirty_payload + (size_t)i * REWIND_PAYLOAD_SIZE,
		                 &full_hot[i], &full_state[i], &full_cold[i]);
		L->dirty_shapes[i] = full_shapes[i];  // transfer ownership
	}

	// Keyframes must always carry a palette so restore never walks past them.
	// Inherit from K when L didn't change materials on its own frame.
	if (!L->materials && K->materials) {
		L->materials = (Material*)CK_ALLOC(256 * sizeof(Material));
		memcpy(L->materials, K->materials, 256 * sizeof(Material));
	}

	CK_FREE(full_hot);
	CK_FREE(full_state);
	CK_FREE(full_cold);
	CK_FREE(full_shapes);
}

static void rewind_evict_oldest(RewindBuffer* rb)
{
	if (rb->count <= 0) return;
	int oldest = rb->head;
	if (rb->count >= 2) {
		int next = (rb->head + 1) % rb->max_frames;
		rewind_promote_to_keyframe(&rb->frames[next], &rb->frames[oldest]);
	}
	rewind_frame_free(&rb->frames[oldest]);
	rb->head = (rb->head + 1) % rb->max_frames;
	rb->count--;
}

// ----- Restore -------------------------------------------------------------

#define rewind_resize(arr, n) do { \
	if ((n) > 0) { afit((arr), (n)); asetlen((arr), (n)); } \
	else if ((arr)) { asetlen((arr), 0); } \
} while(0)

static void rewind_restore_from(WorldInternal* w, RewindBuffer* rb, int target_slot)
{
	const RewindFrame* tf = &rb->frames[target_slot];
	int nb = tf->n_bodies;

	// Walk head-to-target, reconstructing (hot, state, cold) + shape aliases
	// into scratch. scratch_shapes[i] aliases frame-owned memory; we deep-copy
	// at apply time.
	BodyHot*   scratch_hot = NULL;
	BodyState* scratch_state = NULL;
	BodyCold*  scratch_cold = NULL;
	ShapeInternal** scratch_shapes_alias = NULL;
	int scratch_n = 0;

	for (int i = 0; i < rb->count; i++) {
		int slot = (rb->head + i) % rb->max_frames;
		const RewindFrame* f = &rb->frames[slot];

		if (f->n_bodies != scratch_n) {
			BodyHot*   new_hot   = (f->n_bodies > 0) ? (BodyHot*)CK_ALLOC((size_t)f->n_bodies * sizeof(BodyHot)) : NULL;
			BodyState* new_state = (f->n_bodies > 0) ? (BodyState*)CK_ALLOC((size_t)f->n_bodies * sizeof(BodyState)) : NULL;
			BodyCold*  new_cold  = (f->n_bodies > 0) ? (BodyCold*)CK_ALLOC((size_t)f->n_bodies * sizeof(BodyCold)) : NULL;
			ShapeInternal** new_sh = (f->n_bodies > 0) ? (ShapeInternal**)CK_ALLOC((size_t)f->n_bodies * sizeof(ShapeInternal*)) : NULL;
			int carry = scratch_n < f->n_bodies ? scratch_n : f->n_bodies;
			for (int j = 0; j < carry; j++) {
				new_hot[j]   = scratch_hot[j];
				new_state[j] = scratch_state[j];
				new_cold[j]  = scratch_cold[j];
				new_sh[j]    = scratch_shapes_alias[j];
			}
			for (int j = carry; j < f->n_bodies; j++) new_sh[j] = NULL;
			CK_FREE(scratch_hot);
			CK_FREE(scratch_state);
			CK_FREE(scratch_cold);
			CK_FREE(scratch_shapes_alias);
			scratch_hot = new_hot;
			scratch_state = new_state;
			scratch_cold = new_cold;
			scratch_shapes_alias = new_sh;
			scratch_n = f->n_bodies;
		}

		for (int k = 0; k < f->n_dirty; k++) {
			int bi = f->dirty_indices[k];
			uint8_t baseline_buf[REWIND_PAYLOAD_SIZE];
			rewind_pack_body(baseline_buf, &scratch_hot[bi], &scratch_state[bi], &scratch_cold[bi]);
			rewind_decode_body(f, k, baseline_buf, &scratch_hot[bi], &scratch_state[bi], &scratch_cold[bi]);
			scratch_shapes_alias[bi] = f->dirty_shapes[k];
		}

		if (slot == target_slot) break;
	}

	// Free live body shapes BEFORE resizing, so trimmed-off slots release
	// their owned heap arrays.
	int nb_live = asize(w->body_cold);
	for (int i = 0; i < nb_live; i++) {
		afree(w->body_cold[i].shapes);
		w->body_cold[i].shapes = NULL;
	}

	rewind_resize(w->body_hot,   nb);
	rewind_resize(w->body_state, nb);
	rewind_resize(w->body_cold,  nb);
	rewind_resize(w->body_gen,   nb);

	if (nb > 0) {
		memcpy(w->body_hot,   scratch_hot,   (size_t)nb * sizeof(BodyHot));
		memcpy(w->body_state, scratch_state, (size_t)nb * sizeof(BodyState));
		memcpy(w->body_cold,  scratch_cold,  (size_t)nb * sizeof(BodyCold));
		memcpy(w->body_gen,   tf->body_gen,  (size_t)nb * sizeof(uint32_t));
	}
	for (int i = 0; i < nb; i++) {
		w->body_cold[i].shapes = rewind_clone_shapes(scratch_shapes_alias[i]);
	}

	if (w->body_free) aclear(w->body_free);
	int nbf = asize(tf->body_free);
	for (int i = 0; i < nbf; i++) apush(w->body_free, tf->body_free[i]);

	int nj = tf->n_joints;
	rewind_resize(w->joints,    nj);
	rewind_resize(w->joint_hot, nj);
	rewind_resize(w->joint_gen, nj);
	if (nj > 0) {
		memcpy(w->joints,    tf->joints,    (size_t)nj * sizeof(JointInternal));
		memcpy(w->joint_hot, tf->joint_hot, (size_t)nj * sizeof(JointHot));
		memcpy(w->joint_gen, tf->joint_gen, (size_t)nj * sizeof(uint32_t));
	}
	if (w->joint_free) aclear(w->joint_free);
	int njf = asize(tf->joint_free);
	for (int i = 0; i < njf; i++) apush(w->joint_free, tf->joint_free[i]);

	int ni = tf->n_islands;
	int ni_live = asize(w->islands);
	for (int i = 0; i < ni_live; i++) ldl_cache_free(&w->islands[i].ldl);
	rewind_resize(w->islands,    ni);
	rewind_resize(w->island_gen, ni);
	for (int i = 0; i < ni; i++) {
		w->islands[i] = (Island){0};
		w->islands[i].head_body             = tf->islands[i].head_body;
		w->islands[i].tail_body             = tf->islands[i].tail_body;
		w->islands[i].body_count            = tf->islands[i].body_count;
		w->islands[i].head_joint            = tf->islands[i].head_joint;
		w->islands[i].tail_joint            = tf->islands[i].tail_joint;
		w->islands[i].joint_count           = tf->islands[i].joint_count;
		w->islands[i].constraint_remove_count = tf->islands[i].constraint_remove_count;
		w->islands[i].awake                 = tf->islands[i].awake;
	}
	if (ni > 0) memcpy(w->island_gen, tf->island_gen, (size_t)ni * sizeof(uint32_t));
	if (w->island_free) aclear(w->island_free);
	int nif = asize(tf->island_free);
	for (int i = 0; i < nif; i++) apush(w->island_free, tf->island_free[i]);

	map_clear(w->warm_cache);
	for (int i = 0; i < tf->n_warm; i++)
		map_set(w->warm_cache, tf->warm_keys[i], tf->warm_vals[i]);

	map_clear(w->prev_touching);
	for (int i = 0; i < tf->n_prev_touching; i++)
		map_set(w->prev_touching, tf->prev_touching_keys[i], tf->prev_touching_vals[i]);

	map_clear(w->joint_pairs);
	for (int i = 0; i < tf->n_joint_pairs; i++)
		map_set(w->joint_pairs, tf->joint_pairs_keys[i], tf->joint_pairs_vals[i]);

	w->joint_pairs_version = tf->joint_pairs_version;
	w->ldl_topo_version    = tf->ldl_topo_version;
	w->frame               = tf->sim_frame;

	// Sensors: free live sensor shape arrays, resize, reinstall deep-cloned
	// shapes from the snapshot (so the snapshot remains reusable).
	int ns_live = asize(w->sensors);
	for (int i = 0; i < ns_live; i++) { afree(w->sensors[i].shapes); w->sensors[i].shapes = NULL; }
	rewind_resize(w->sensors,    tf->n_sensors);
	rewind_resize(w->sensor_gen, tf->n_sensors);
	for (int i = 0; i < tf->n_sensors; i++) {
		w->sensors[i] = tf->sensors[i];
		w->sensors[i].shapes = NULL;
		int nshapes = asize(tf->sensors[i].shapes);
		if (nshapes > 0) {
			ShapeInternal* dst = NULL;
			afit(dst, nshapes); asetlen(dst, nshapes);
			memcpy(dst, tf->sensors[i].shapes, (size_t)nshapes * sizeof(ShapeInternal));
			w->sensors[i].shapes = dst;
		}
	}
	if (tf->n_sensors > 0 && tf->sensor_gen)
		memcpy(w->sensor_gen, tf->sensor_gen, (size_t)tf->n_sensors * sizeof(uint32_t));
	if (w->sensor_free) aclear(w->sensor_free);
	int nsf = asize(tf->sensor_free);
	for (int i = 0; i < nsf; i++) apush(w->sensor_free, tf->sensor_free[i]);

	// Material palette: walk backward from target to the nearest frame with
	// a captured palette. Keyframes always carry one, so the scan is bounded
	// by ring depth; in practice the target is usually the first non-NULL.
	for (int i = 0; i < rb->count; i++) {
		int step = (target_slot - i + rb->max_frames) % rb->max_frames;
		if (rb->frames[step].materials) {
			memcpy(w->materials, rb->frames[step].materials, sizeof(w->materials));
			break;
		}
	}
	memcpy(rb->baseline_materials, w->materials, sizeof(w->materials));

	if (w->broadphase_type == BROADPHASE_BVH) {
		bvh_free(w->bvh_static);   bvh_init(w->bvh_static);
		bvh_free(w->bvh_dynamic);  bvh_init(w->bvh_dynamic);
		bvh_free(w->bvh_sleeping); bvh_init(w->bvh_sleeping);
		for (int i = 0; i < nb; i++) {
			if (!split_alive(w->body_gen, i)) { w->body_cold[i].bvh_leaf = -1; continue; }
			if (asize(w->body_cold[i].shapes) == 0) { w->body_cold[i].bvh_leaf = -1; continue; }
			BVH_Tree* tree;
			if (body_inv_mass(w, i) == 0.0f) tree = w->bvh_static;
			else {
				int isl = w->body_cold[i].island_id;
				int sleeping = (isl >= 0 && (w->island_gen[isl] & 1) && !w->islands[isl].awake);
				tree = sleeping ? w->bvh_sleeping : w->bvh_dynamic;
			}
			AABB box = aabb_expand(body_aabb(&w->body_state[i], &w->body_cold[i]), BVH_AABB_MARGIN);
			w->body_cold[i].bvh_leaf = bvh_insert(tree, i, box);
		}
	} else {
		for (int i = 0; i < nb; i++) w->body_cold[i].bvh_leaf = -1;
	}

	// Sync baseline to restored state so next capture diffs correctly.
	if (nb != rb->baseline_n_bodies) {
		CK_FREE(rb->baseline_hot);   rb->baseline_hot = NULL;
		CK_FREE(rb->baseline_state); rb->baseline_state = NULL;
		CK_FREE(rb->baseline_cold);  rb->baseline_cold = NULL;
		if (rb->baseline_shapes) {
			for (int i = 0; i < rb->baseline_n_bodies; i++) afree(rb->baseline_shapes[i]);
			CK_FREE(rb->baseline_shapes);
			rb->baseline_shapes = NULL;
		}
		if (nb > 0) {
			rb->baseline_hot    = (BodyHot*)CK_ALLOC((size_t)nb * sizeof(BodyHot));
			rb->baseline_state  = (BodyState*)CK_ALLOC((size_t)nb * sizeof(BodyState));
			rb->baseline_cold   = (BodyCold*)CK_ALLOC((size_t)nb * sizeof(BodyCold));
			rb->baseline_shapes = (ShapeInternal**)CK_ALLOC((size_t)nb * sizeof(ShapeInternal*));
			for (int i = 0; i < nb; i++) rb->baseline_shapes[i] = NULL;
		}
		rb->baseline_n_bodies = nb;
	}
	if (nb > 0) {
		memcpy(rb->baseline_hot,   w->body_hot,   (size_t)nb * sizeof(BodyHot));
		memcpy(rb->baseline_state, w->body_state, (size_t)nb * sizeof(BodyState));
		for (int i = 0; i < nb; i++) {
			rb->baseline_cold[i] = w->body_cold[i];
			rb->baseline_cold[i].shapes = NULL;
			if (!rewind_shapes_equal(w->body_cold[i].shapes, rb->baseline_shapes[i])) {
				afree(rb->baseline_shapes[i]);
				rb->baseline_shapes[i] = rewind_clone_shapes(w->body_cold[i].shapes);
			}
		}
	}

	CK_FREE(scratch_hot);
	CK_FREE(scratch_state);
	CK_FREE(scratch_cold);
	CK_FREE(scratch_shapes_alias);
}

// ----- Ring helpers --------------------------------------------------------

static int rewind_ring_find(const RewindBuffer* rb, uint64_t frame_id)
{
	for (int i = 0; i < rb->count; i++) {
		int slot = (rb->head + i) % rb->max_frames;
		if (rb->frames[slot].frame_id == frame_id) return slot;
	}
	return -1;
}

static void rewind_buffer_clear(RewindBuffer* rb)
{
	for (int i = 0; i < rb->count; i++) {
		int slot = (rb->head + i) % rb->max_frames;
		rewind_frame_free(&rb->frames[slot]);
	}
	rb->head = 0;
	rb->count = 0;
}

// ----- Public API ----------------------------------------------------------

void world_rewind_init(World world, RewindParams params)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (w->rewind) { world_rewind_shutdown(world); }
	if (params.max_frames <= 0) return;

	RewindBuffer* rb = (RewindBuffer*)CK_ALLOC(sizeof(RewindBuffer));
	memset(rb, 0, sizeof(*rb));
	rb->max_frames   = params.max_frames;
	rb->auto_capture = params.auto_capture;
	afit(rb->frames, rb->max_frames);
	asetlen(rb->frames, rb->max_frames);
	memset(rb->frames, 0, (size_t)rb->max_frames * sizeof(RewindFrame));
	w->rewind = rb;
}

void world_rewind_shutdown(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	if (!w->rewind) return;
	RewindBuffer* rb = w->rewind;
	rewind_buffer_clear(rb);
	afree(rb->frames);
	CK_FREE(rb->baseline_hot);
	CK_FREE(rb->baseline_state);
	CK_FREE(rb->baseline_cold);
	if (rb->baseline_shapes) {
		for (int i = 0; i < rb->baseline_n_bodies; i++) afree(rb->baseline_shapes[i]);
		CK_FREE(rb->baseline_shapes);
	}
	CK_FREE(rb);
	w->rewind = NULL;
}

uint64_t world_rewind_capture(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	RewindBuffer* rb = w->rewind;
	if (!rb || rb->max_frames <= 0) return 0;

	int slot;
	if (rb->count < rb->max_frames) {
		slot = (rb->head + rb->count) % rb->max_frames;
		rb->count++;
	} else {
		rewind_evict_oldest(rb);
		slot = (rb->head + rb->count) % rb->max_frames;
		rb->count++;
	}
	rewind_capture_into(rb, &rb->frames[slot], w);
	return rb->frames[slot].frame_id;
}

int world_rewind_to_frame(World world, uint64_t frame_id)
{
	WorldInternal* w = (WorldInternal*)world.id;
	RewindBuffer* rb = w->rewind;
	if (!rb) return 0;
	int slot = rewind_ring_find(rb, frame_id);
	if (slot < 0) return 0;
	rewind_restore_from(w, rb, slot);
	return 1;
}

int world_rewind_by_steps(World world, int n)
{
	WorldInternal* w = (WorldInternal*)world.id;
	RewindBuffer* rb = w->rewind;
	if (!rb || rb->count == 0) return 0;
	if (n < 0 || n >= rb->count) return 0;
	int slot = (rb->head + rb->count - 1 - n) % rb->max_frames;
	rewind_restore_from(w, rb, slot);
	return 1;
}

int world_rewind_frames_available(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	return w->rewind ? w->rewind->count : 0;
}

size_t world_rewind_memory_used(World world)
{
	WorldInternal* w = (WorldInternal*)world.id;
	RewindBuffer* rb = w->rewind;
	if (!rb) return 0;
	size_t total = sizeof(RewindBuffer) + (size_t)rb->max_frames * sizeof(RewindFrame);
	total += (size_t)rb->baseline_n_bodies * (sizeof(BodyHot) + sizeof(BodyState) + sizeof(BodyCold) + sizeof(ShapeInternal*));
	if (rb->baseline_shapes) {
		for (int i = 0; i < rb->baseline_n_bodies; i++)
			total += (size_t)asize(rb->baseline_shapes[i]) * sizeof(ShapeInternal);
	}
	for (int i = 0; i < rb->count; i++) {
		int slot = (rb->head + i) % rb->max_frames;
		const RewindFrame* f = &rb->frames[slot];
		total += (size_t)f->n_dirty * (sizeof(int) + sizeof(ShapeInternal*));
		total += (size_t)f->dirty_payload_size;
		if (!f->is_keyframe && f->dirty_payload_offsets)
			total += (size_t)(f->n_dirty + 1) * sizeof(int);
		if (f->dirty_shapes) {
			for (int k = 0; k < f->n_dirty; k++)
				total += (size_t)asize(f->dirty_shapes[k]) * sizeof(ShapeInternal);
		}
		total += (size_t)f->n_bodies * sizeof(uint32_t);
		total += (size_t)asize(f->body_free) * sizeof(int);
		total += (size_t)f->n_joints  * (sizeof(JointInternal) + sizeof(JointHot) + sizeof(uint32_t));
		total += (size_t)asize(f->joint_free) * sizeof(int);
		total += (size_t)f->n_islands * (sizeof(RewindIsland) + sizeof(uint32_t));
		total += (size_t)asize(f->island_free) * sizeof(int);
		total += (size_t)f->n_warm          * (sizeof(uint64_t) + sizeof(WarmManifold));
		total += (size_t)f->n_prev_touching * (sizeof(uint64_t) + sizeof(uint8_t));
		total += (size_t)f->n_joint_pairs   * (sizeof(uint64_t) + sizeof(uint8_t));
		if (f->materials) total += 256 * sizeof(Material);
	}
	return total;
}
