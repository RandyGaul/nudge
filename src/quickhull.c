// See LICENSE for licensing info.
// quickhull.c -- 3D Quickhull convex hull algorithm.
//
// Builds a convex hull from a point cloud using conflict lists,
// recursive horizon computation, and two-pass non-convex face merging.
// Degenerate faces (triangle collapse, vertex skip) are handled
// on-the-fly during the merge splice to maintain manifold topology.
//
// Produces a Hull with polygonal faces, half-edge mesh, and face planes.
//
// Edge convention: each half-edge stores its ORIGIN vertex (= tail).
// The head of edge e is edges[edges[e].next].origin.

#define QH_INVALID    (~0)
#define QH_VISIBLE    1
#define QH_NON_CONVEX 2
#define QH_DELETED    3
#define QH_EDGE_DELETED 0x01

// Debug state: set qh_debug=1 to enable runtime logging.
static int        qh_debug;
static const v3*  qh_dbg_points;
static int        qh_dbg_count;

#define QH_DEBUG(...) do { if (qh_debug) { __VA_ARGS__; } } while(0)

#define QH_ASSERT(cond, msg) do { \
	if (!(cond)) { \
		fprintf(stderr, "quickhull: %s (%s:%d)\n", msg, __FILE__, __LINE__); \
		fprintf(stderr, "  input: %d points\n  v3 pts[] = {\n", qh_dbg_count); \
		for (int _i = 0; _i < qh_dbg_count; _i++) \
			fprintf(stderr, "    {%.9ef, %.9ef, %.9ef},\n", \
				qh_dbg_points[_i].x, qh_dbg_points[_i].y, qh_dbg_points[_i].z); \
		fprintf(stderr, "  };\n"); \
		exit(-1); \
	} \
} while(0)

// -----------------------------------------------------------------------------
// Internal types.

// Vertex data in SoA layout for SIMD-friendly access.
// Positions in aligned float arrays, conflict links in separate int arrays.
typedef struct QH_Verts
{
	float* x;   // aligned position arrays
	float* y;
	float* z;
	int* cnext; // conflict list next (circular doubly-linked)
	int* cprev; // conflict list prev
	int count;
	int cap;
} QH_Verts;

// Edge data in SoA layout. Topology walks only touch enext[]/etwin[]
// (4 bytes/step) instead of striding through 24-byte structs.
typedef struct QH_Edges
{
	int* enext;
	int* eprev;
	int* etwin;
	int* eorigin;
	int* eface;
	int count;
	int cap;
	int free_head; // free list head (linked via enext)
} QH_Edges;

// Hot/cold face split: hot fields (scanned by qh_next_conflict + merge)
// are packed together; cold fields (topology, metadata) accessed only
// during cone creation, merge splicing, and output.
typedef struct QH_Face
{
	// --- hot: conflict scanning + merge distance checks ---
	HullPlane plane;     // 16 bytes
	v3 centroid;         // 12 bytes
	float maxoutside;    // 4 bytes
	int conflict_head;   // 4 bytes
	int mark;            // 4 bytes
	// --- cold: topology + metadata ---
	int edge;            // 4 bytes
	int next, prev;      // 8 bytes
	int conflict_slot;   // 4 bytes
	int num_verts;       // 4 bytes
	float area;          // 4 bytes
} QH_Face;               // 64 bytes total, hot fields in first 40

typedef struct QH_State
{
	QH_Verts verts;
	QH_Edges edges;
	CK_DYNA QH_Face*   faces;
	CK_DYNA int*       conflict_faces; // compact list of faces with non-empty conflict lists
	int face_free;
	v3 interior;
	float epsilon;
} QH_State;

typedef struct QH_FaceList
{
	int first, last;
} QH_FaceList;

// -----------------------------------------------------------------------------
// SoA edge management.

static void qh_edges_reserve(QH_Edges* E, int n)
{
	if (n <= E->cap) return;
	int newcap = E->cap ? E->cap * 2 : 16;
	while (newcap < n) newcap *= 2;
	int* nn = (int*)realloc(E->enext, newcap * sizeof(int));
	int* np = (int*)realloc(E->eprev, newcap * sizeof(int));
	int* nt = (int*)realloc(E->etwin, newcap * sizeof(int));
	int* no = (int*)realloc(E->eorigin, newcap * sizeof(int));
	int* nf = (int*)realloc(E->eface, newcap * sizeof(int));
	E->enext = nn; E->eprev = np; E->etwin = nt; E->eorigin = no; E->eface = nf;
	E->cap = newcap;
}

static void qh_edges_free(QH_Edges* E)
{
	free(E->enext); free(E->eprev); free(E->etwin); free(E->eorigin); free(E->eface);
	*E = (QH_Edges){0};
	E->free_head = QH_INVALID;
}

static int qh_alloc_edge(QH_State* s)
{
	QH_Edges* E = &s->edges;
	if (E->free_head != QH_INVALID) {
		int idx = E->free_head;
		E->free_head = E->enext[idx];
		return idx;
	}
	if (E->count >= E->cap) qh_edges_reserve(E, E->count + 1);
	return E->count++;
}

static int qh_alloc_face(QH_State* s)
{
	if (s->face_free != QH_INVALID) {
		int idx = s->face_free;
		s->face_free = s->faces[idx].next;
		s->faces[idx].conflict_head = QH_INVALID;
		s->faces[idx].conflict_slot = -1;
		s->faces[idx].mark = QH_VISIBLE;
		s->faces[idx].maxoutside = s->epsilon;
		return idx; // caller sets edge, next, prev, plane, centroid, num_verts, area
	}
	QH_Face f = {0};
	f.conflict_head = QH_INVALID;
	f.conflict_slot = -1;
	f.mark = QH_VISIBLE;
	f.maxoutside = s->epsilon;
	apush(s->faces, f);
	return asize(s->faces) - 1;
}

// SoA vertex management.
static void qh_verts_reserve(QH_Verts* v, int n)
{
	if (n <= v->cap) return;
	int newcap = v->cap ? v->cap * 2 : 16;
	while (newcap < n) newcap *= 2;
	int align = 16;
	float* nx = (float*)CK_ALLOC_ALIGNED(newcap * sizeof(float), align);
	float* ny = (float*)CK_ALLOC_ALIGNED(newcap * sizeof(float), align);
	float* nz = (float*)CK_ALLOC_ALIGNED(newcap * sizeof(float), align);
	int* ncnext = (int*)CK_ALLOC_ALIGNED(newcap * sizeof(int), align);
	int* ncprev = (int*)CK_ALLOC_ALIGNED(newcap * sizeof(int), align);
	if (v->count > 0) {
		memcpy(nx, v->x, v->count * sizeof(float));
		memcpy(ny, v->y, v->count * sizeof(float));
		memcpy(nz, v->z, v->count * sizeof(float));
		memcpy(ncnext, v->cnext, v->count * sizeof(int));
		memcpy(ncprev, v->cprev, v->count * sizeof(int));
	}
	if (v->x) { CK_FREE_ALIGNED(v->x); CK_FREE_ALIGNED(v->y); CK_FREE_ALIGNED(v->z); CK_FREE_ALIGNED(v->cnext); CK_FREE_ALIGNED(v->cprev); }
	v->x = nx; v->y = ny; v->z = nz; v->cnext = ncnext; v->cprev = ncprev; v->cap = newcap;
}

static int qh_verts_push(QH_Verts* v, v3 pos)
{
	if (v->count >= v->cap) qh_verts_reserve(v, v->count + 1);
	int idx = v->count++;
	v->x[idx] = pos.x; v->y[idx] = pos.y; v->z[idx] = pos.z;
	v->cnext[idx] = QH_INVALID; v->cprev[idx] = QH_INVALID;
	return idx;
}

static v3 qh_vert_pos(QH_Verts* v, int i) { return V3(v->x[i], v->y[i], v->z[i]); }

static void qh_verts_free(QH_Verts* v)
{
	if (v->x) { CK_FREE_ALIGNED(v->x); CK_FREE_ALIGNED(v->y); CK_FREE_ALIGNED(v->z); CK_FREE_ALIGNED(v->cnext); CK_FREE_ALIGNED(v->cprev); }
	*v = (QH_Verts){0};
}

// Conflict face list: compact array of face indices with non-empty conflict lists.
static void qh_cfl_add(QH_State* s, int fi)
{
	if (s->faces[fi].conflict_slot >= 0) return;
	s->faces[fi].conflict_slot = asize(s->conflict_faces);
	apush(s->conflict_faces, fi);
}

static void qh_cfl_remove(QH_State* s, int fi)
{
	int slot = s->faces[fi].conflict_slot;
	if (slot < 0) return;
	int last = asize(s->conflict_faces) - 1;
	if (slot != last) {
		int moved = s->conflict_faces[last];
		s->conflict_faces[slot] = moved;
		s->faces[moved].conflict_slot = slot;
	}
	asetlen(s->conflict_faces, last);
	s->faces[fi].conflict_slot = -1;
}

// Forward declaration (used by qh_conflict_add before full definition).
static float qh_plane_dist(HullPlane p, v3 pt);

// -----------------------------------------------------------------------------
// Conflict list (circular doubly-linked via vertex next/prev).

static void qh_conflict_add(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Verts* V = &s->verts;
	if (f->conflict_head != QH_INVALID) {
		int tail = V->cprev[f->conflict_head];
		V->cnext[vi] = f->conflict_head;
		V->cprev[vi] = tail;
		V->cprev[f->conflict_head] = vi;
		V->cnext[tail] = vi;
	} else {
		V->cnext[vi] = vi;
		V->cprev[vi] = vi;
		qh_cfl_add(s, fi);
	}
	f->conflict_head = vi;
	float d = qh_plane_dist(f->plane, qh_vert_pos(V, vi));
	if (d > f->maxoutside) f->maxoutside = d;
}

static void qh_conflict_remove(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Verts* V = &s->verts;
	if (V->cnext[vi] == vi) {
		f->conflict_head = QH_INVALID;
		qh_cfl_remove(s, fi);
	} else {
		V->cprev[V->cnext[vi]] = V->cprev[vi];
		V->cnext[V->cprev[vi]] = V->cnext[vi];
		if (f->conflict_head == vi) f->conflict_head = V->cnext[vi];
	}
	V->cnext[vi] = QH_INVALID;
	V->cprev[vi] = QH_INVALID;
}

static int qh_conflict_remove_all(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	int head = f->conflict_head;
	if (head == QH_INVALID) return QH_INVALID;
	int last = s->verts.cprev[head];
	s->verts.cnext[last] = QH_INVALID;
	f->conflict_head = QH_INVALID;
	qh_cfl_remove(s, fi);
	return head;
}

// -----------------------------------------------------------------------------
// Geometry helpers.

static float qh_plane_dist(HullPlane p, v3 pt)
{
	return dot(p.normal, pt) - p.offset;
}

// Compute face normal, centroid, area via Newell method.
// Uses scalar SoA access to avoid v3/simd4f lane extraction overhead.
static void qh_recompute_face(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	QH_Verts* V = &s->verts;
	float nx = 0, ny = 0, nz = 0, cx = 0, cy = 0, cz = 0;
	int count = 0, e = f->edge;
	do {
		int vi = s->edges.eorigin[e], vn = s->edges.eorigin[s->edges.enext[e]];
		float curx = V->x[vi], cury = V->y[vi], curz = V->z[vi];
		float nxtx = V->x[vn], nxty = V->y[vn], nxtz = V->z[vn];
		nx += (cury - nxty) * (curz + nxtz);
		ny += (curz - nxtz) * (curx + nxtx);
		nz += (curx - nxtx) * (cury + nxty);
		cx += curx; cy += cury; cz += curz;
		count++;
		if (count > 1000) { QH_ASSERT(0, "infinite face loop in qh_recompute_face"); }
		e = s->edges.enext[e];
	} while (e != f->edge);

	float a = sqrtf(nx*nx + ny*ny + nz*nz);
	if (a > 0) { float inv = 1.0f / a; nx *= inv; ny *= inv; nz *= inv; }
	float inv_c = 1.0f / count;
	cx *= inv_c; cy *= inv_c; cz *= inv_c;
	v3 normal = V3(nx, ny, nz), centroid = V3(cx, cy, cz);
	f->plane = (HullPlane){ normal, dot(normal, centroid) };
	f->centroid = centroid;
	f->num_verts = count;
	f->area = a;
}

static float qh_face_dist(QH_State* s, int fi, v3 pt)
{
	return qh_plane_dist(s->faces[fi].plane, pt);
}

// Distance of the adjacent face's centroid to this edge's face plane.
static float qh_opp_face_dist(QH_State* s, int ei)
{
	int fi = s->edges.eface[ei];
	int ofi = s->edges.eface[s->edges.etwin[ei]];
	return qh_face_dist(s, fi, s->faces[ofi].centroid);
}

static int qh_face_vert_count(QH_State* s, int fi)
{
	int e = s->faces[fi].edge, start = e, n = 0;
	do { n++; e = s->edges.enext[e]; } while (e != start);
	return n;
}

static int qh_opp_face(QH_State* s, int ei) { return s->edges.eface[s->edges.etwin[ei]]; }
static int qh_edge_head(QH_State* s, int ei) { return s->edges.eorigin[s->edges.enext[ei]]; }
static int qh_edge_tail(QH_State* s, int ei) { return s->edges.eorigin[ei]; }

// Get edge at index i from face's he0 (supports negative indices).
static int qh_get_edge(QH_State* s, int fi, int i)
{
	int e = s->faces[fi].edge;
	while (i > 0) { e = s->edges.enext[e]; i--; }
	while (i < 0) { e = s->edges.eprev[e]; i++; }
	return e;
}

static void qh_set_opposite(QH_State* s, int a, int b)
{
	s->edges.etwin[a] = b;
	s->edges.etwin[b] = a;
}

// Create a triangle face with edges v0->v1->v2->v0.
static int qh_create_triangle(QH_State* s, int v0, int v1, int v2)
{
	int fi = qh_alloc_face(s);
	int e0 = qh_alloc_edge(s), e1 = qh_alloc_edge(s), e2 = qh_alloc_edge(s);
	QH_Edges* E = &s->edges;
	E->enext[e0] = e1; E->eprev[e0] = e2; E->etwin[e0] = QH_INVALID; E->eorigin[e0] = v0; E->eface[e0] = fi;
	E->enext[e1] = e2; E->eprev[e1] = e0; E->etwin[e1] = QH_INVALID; E->eorigin[e1] = v1; E->eface[e1] = fi;
	E->enext[e2] = e0; E->eprev[e2] = e1; E->etwin[e2] = QH_INVALID; E->eorigin[e2] = v2; E->eface[e2] = fi;
	s->faces[fi].edge = e0;
	s->faces[fi].next = fi;
	s->faces[fi].prev = fi;
	// Fast triangle normal via cross product (avoids Newell loop).
	v3 p0 = qh_vert_pos(&s->verts, v0), p1 = qh_vert_pos(&s->verts, v1), p2 = qh_vert_pos(&s->verts, v2);
	v3 n = cross(sub(p1, p0), sub(p2, p0));
	float a = len(n);
	if (a > 0) n = scale(n, 1.0f / a);
	v3 c = scale(add(add(p0, p1), p2), 1.0f / 3.0f);
	s->faces[fi].plane = (HullPlane){ n, dot(n, c) };
	s->faces[fi].centroid = c;
	s->faces[fi].num_verts = 3;
	s->faces[fi].area = a;
	return fi;
}

// -----------------------------------------------------------------------------
// Build initial tetrahedron from 4 non-coplanar extremal points.

static int qh_build_simplex(QH_State* s, int nv)
{
	// Find extremal vertices on each axis.
	int maxV[3], minV[3];
	for (int i = 0; i < 3; i++) { maxV[i] = minV[i] = 0; }
	v3 mx = qh_vert_pos(&s->verts, 0), mn = qh_vert_pos(&s->verts, 0);
	for (int i = 1; i < nv; i++) {
		v3 p = qh_vert_pos(&s->verts, i);
		if (p.x > mx.x) { mx.x = p.x; maxV[0] = i; } else if (p.x < mn.x) { mn.x = p.x; minV[0] = i; }
		if (p.y > mx.y) { mx.y = p.y; maxV[1] = i; } else if (p.y < mn.y) { mn.y = p.y; minV[1] = i; }
		if (p.z > mx.z) { mx.z = p.z; maxV[2] = i; } else if (p.z < mn.z) { mn.z = p.z; minV[2] = i; }
	}

	// Pick axis with greatest spread.
	float best = 0; int imax = 0;
	for (int i = 0; i < 3; i++) {
		v3 a = qh_vert_pos(&s->verts, maxV[i]), b = qh_vert_pos(&s->verts, minV[i]);
		float d = (i == 0) ? a.x - b.x : (i == 1) ? a.y - b.y : a.z - b.z;
		if (d > best) { best = d; imax = i; }
	}
	if (best <= s->epsilon) return 0;

	int vtx[4];
	vtx[0] = maxV[imax]; vtx[1] = minV[imax];

	// Find vertex farthest from line vtx[0]-vtx[1].
	v3 u01 = norm(sub(qh_vert_pos(&s->verts, vtx[1]), qh_vert_pos(&s->verts, vtx[0])));
	float maxSqr = 0; v3 nrml = V3(0, 0, 0); vtx[2] = -1;
	for (int i = 0; i < nv; i++) {
		v3 xp = cross(u01, sub(qh_vert_pos(&s->verts, i), qh_vert_pos(&s->verts, vtx[0])));
		float ls = len2(xp);
		if (ls > maxSqr && i != vtx[0] && i != vtx[1]) { maxSqr = ls; vtx[2] = i; nrml = xp; }
	}
	if (vtx[2] < 0 || sqrtf(maxSqr) <= 100 * s->epsilon) return 0;
	nrml = norm(nrml);
	nrml = norm(sub(nrml, scale(u01, dot(nrml, u01)))); // orthogonalize

	// Find vertex farthest from plane through vtx[0..2].
	float d0 = dot(qh_vert_pos(&s->verts, vtx[2]), nrml);
	float maxDist = 0; vtx[3] = -1;
	for (int i = 0; i < nv; i++) {
		float d = fabsf(dot(qh_vert_pos(&s->verts, i), nrml) - d0);
		if (d > maxDist && i != vtx[0] && i != vtx[1] && i != vtx[2]) { maxDist = d; vtx[3] = i; }
	}
	if (vtx[3] < 0 || maxDist <= 100 * s->epsilon) return 0;

	// Build 4 triangle faces with correct winding.
	// With origin-vertex convention, edge indices within a triangle: e0, e1, e2.
	int tris[4];
	if (dot(qh_vert_pos(&s->verts, vtx[3]), nrml) - d0 < 0) {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[1], vtx[2]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[1], vtx[0]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[2], vtx[1]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[0], vtx[2]);
		for (int i = 0; i < 3; i++) {
			int k = (i + 1) % 3;
			qh_set_opposite(s, qh_get_edge(s, tris[i+1], 0), qh_get_edge(s, tris[k+1], 2));
			int our_k = (k == 0) ? 2 : (k == 1) ? 0 : 1;
			qh_set_opposite(s, qh_get_edge(s, tris[i+1], 1), qh_get_edge(s, tris[0], our_k));
		}
	} else {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[2], vtx[1]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[0], vtx[1]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[1], vtx[2]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[2], vtx[0]);
		for (int i = 0; i < 3; i++) {
			int k = (i + 1) % 3;
			qh_set_opposite(s, qh_get_edge(s, tris[i+1], 2), qh_get_edge(s, tris[k+1], 0));
			int j = (3 - i) % 3, our_j = (j == 0) ? 2 : (j == 1) ? 0 : 1;
			qh_set_opposite(s, qh_get_edge(s, tris[i+1], 1), qh_get_edge(s, tris[0], our_j));
		}
	}

	s->interior = scale(add(add(qh_vert_pos(&s->verts, vtx[0]), qh_vert_pos(&s->verts, vtx[1])), add(qh_vert_pos(&s->verts, vtx[2]), qh_vert_pos(&s->verts, vtx[3]))), 0.25f);

	// Assign conflict vertices: SIMD 4-wide dot product against all 4 initial faces.
	{
		HullPlane tp[4]; for (int k = 0; k < 4; k++) tp[k] = s->faces[tris[k]].plane;
		simd4f tnx = simd_set(tp[0].normal.x, tp[1].normal.x, tp[2].normal.x, tp[3].normal.x);
		simd4f tny = simd_set(tp[0].normal.y, tp[1].normal.y, tp[2].normal.y, tp[3].normal.y);
		simd4f tnz = simd_set(tp[0].normal.z, tp[1].normal.z, tp[2].normal.z, tp[3].normal.z);
		simd4f toff = simd_set(tp[0].offset, tp[1].offset, tp[2].offset, tp[3].offset);
		simd4f veps = simd_set1(s->epsilon);
		for (int i = 0; i < nv; i++) {
			if (i == vtx[0] || i == vtx[1] || i == vtx[2] || i == vtx[3]) continue;
			v3 p = qh_vert_pos(&s->verts, i);
			simd4f d = simd_sub(simd_add(simd_add(simd_mul(tnx, simd_set1(p.x)), simd_mul(tny, simd_set1(p.y))), simd_mul(tnz, simd_set1(p.z))), toff);
			// Find lane with max distance above epsilon.
			simd4f above = simd_cmpgt(d, veps);
			int mask = simd_movemask(above);
			if (mask) {
				float ds[4]; simd_store(ds, d);
				float bd = s->epsilon; int bf = -1;
				for (int k = 0; k < 4; k++) if (ds[k] > bd) { bd = ds[k]; bf = tris[k]; }
				if (bf >= 0) qh_conflict_add(s, bf, i);
			}
		}
	}
	return 1;
}

// -----------------------------------------------------------------------------
// Find next conflict vertex (globally furthest from any face).

static int qh_next_conflict(QH_State* s, int* out_face)
{
	int bv = QH_INVALID, bf = QH_INVALID; float bd = -1e18f;
	QH_Face* faces = s->faces;
	QH_Verts* V = &s->verts;
	int ncf = asize(s->conflict_faces);
	for (int ci_idx = 0; ci_idx < ncf; ci_idx++) {
		int fi = s->conflict_faces[ci_idx];
		if (faces[fi].maxoutside <= bd) continue;
		HullPlane p = faces[fi].plane;
		int head = faces[fi].conflict_head;
		int ci = head;
		do {
			float d = p.normal.x * V->x[ci] + p.normal.y * V->y[ci] + p.normal.z * V->z[ci] - p.offset;
			if (d > bd) { bd = d; bv = ci; bf = fi; }
			ci = V->cnext[ci];
		} while (ci != head);
	}
	*out_face = bf;
	return bv;
}

// -----------------------------------------------------------------------------
// Move conflict vertices from a deleted face to unclaimed or an absorbing face.

static void qh_delete_face_points(QH_State* s, int fi, int absorb, int* unclaimed)
{
	int vlist = qh_conflict_remove_all(s, fi);
	if (vlist == QH_INVALID) return;
	if (absorb == QH_INVALID) {
		// Fast path: push all vertices to unclaimed (no distance check).
		int v = vlist;
		while (s->verts.cnext[v] != QH_INVALID) v = s->verts.cnext[v];
		s->verts.cnext[v] = *unclaimed;
		*unclaimed = vlist;
		return;
	}
	float eps = s->epsilon;
	HullPlane ap = s->faces[absorb].plane;
	QH_Verts* V = &s->verts;
	int v = vlist;
	while (v != QH_INVALID) {
		int nxt = V->cnext[v];
		if (ap.normal.x * V->x[v] + ap.normal.y * V->y[v] + ap.normal.z * V->z[v] - ap.offset > eps) {
			qh_conflict_add(s, absorb, v);
		} else {
			s->verts.cnext[v] = *unclaimed;
			*unclaimed = v;
		}
		v = nxt;
	}
}

static void qh_face_list_add(QH_State* s, QH_FaceList* list, int fi);

// Fused horizon DFS + cone creation: find horizon edges and immediately
// build cone triangles during DFS traversal. Eliminates the intermediate
// horizon array and fuses two passes into one.

typedef struct { int edge0; int edge; } QH_HFrame;

static void qh_horizon_and_cone(QH_State* s, int eye, int eye_face, QH_FaceList* nf, int* unclaimed)
{
	QH_HFrame hstk_buf[64];
	int hstk_n = 0;
	#define HSTK_PUSH(val) do { assert(hstk_n < 64); hstk_buf[hstk_n++] = (val); } while(0)
	#define HSTK_POP() (hstk_n--)
	#define HSTK_TOP() (&hstk_buf[hstk_n - 1])
	#define HSTK_EMPTY() (hstk_n == 0)

	nf->first = nf->last = QH_INVALID;
	int prev = QH_INVALID, begin = QH_INVALID;
	v3 eye_pos = qh_vert_pos(&s->verts, eye);

	qh_delete_face_points(s, eye_face, QH_INVALID, unclaimed);
	s->faces[eye_face].mark = QH_DELETED;
	int e0 = s->faces[eye_face].edge;
	HSTK_PUSH(((QH_HFrame){ e0, e0 }));

	float eps = s->epsilon;
	float eye_x = eye_pos.x, eye_y = eye_pos.y, eye_z = eye_pos.z;
	QH_Edges* E = &s->edges;
	while (!HSTK_EMPTY()) {
		QH_HFrame* f = HSTK_TOP();
		int twin = E->etwin[f->edge];
		int ofi = E->eface[twin];
		if (s->faces[ofi].mark == QH_VISIBLE) {
			HullPlane op = s->faces[ofi].plane;
			if (op.normal.x * eye_x + op.normal.y * eye_y + op.normal.z * eye_z - op.offset > eps) {
				qh_delete_face_points(s, ofi, QH_INVALID, unclaimed);
				s->faces[ofi].mark = QH_DELETED;
				HSTK_PUSH(((QH_HFrame){ twin, E->enext[twin] }));
				continue;
			} else {
				// Horizon edge found: create cone triangle inline.
				int horizon_edge = f->edge;
				int tail = E->eorigin[horizon_edge];
				int head = E->eorigin[E->enext[horizon_edge]];
				int fi = qh_alloc_face(s);
				int ce0 = qh_alloc_edge(s), ce1 = qh_alloc_edge(s), ce2 = qh_alloc_edge(s);
				E = &s->edges; // re-fetch after potential realloc in qh_alloc_edge
				E->enext[ce0] = ce1; E->eprev[ce0] = ce2; E->etwin[ce0] = QH_INVALID; E->eorigin[ce0] = eye; E->eface[ce0] = fi;
				E->enext[ce1] = ce2; E->eprev[ce1] = ce0; E->etwin[ce1] = QH_INVALID; E->eorigin[ce1] = tail; E->eface[ce1] = fi;
				E->enext[ce2] = ce0; E->eprev[ce2] = ce1; E->etwin[ce2] = QH_INVALID; E->eorigin[ce2] = head; E->eface[ce2] = fi;
				s->faces[fi].edge = ce0;
				s->faces[fi].next = fi;
				s->faces[fi].prev = fi;
				v3 p1 = qh_vert_pos(&s->verts, tail), p2 = qh_vert_pos(&s->verts, head);
				v3 n = cross(sub(p1, eye_pos), sub(p2, eye_pos));
				float a = len(n);
				if (a > 0) n = scale(n, 1.0f / a);
				v3 c = scale(add(add(eye_pos, p1), p2), 1.0f / 3.0f);
				s->faces[fi].plane = (HullPlane){ n, dot(n, c) };
				s->faces[fi].centroid = c;
				s->faces[fi].num_verts = 3;
				s->faces[fi].area = a;
				qh_set_opposite(s, ce1, E->etwin[horizon_edge]);
				if (prev != QH_INVALID) qh_set_opposite(s, ce0, prev);
				else begin = ce0;
				qh_face_list_add(s, nf, fi);
				prev = ce2;
			}
		}
		E = &s->edges; // re-fetch after potential realloc
		f = HSTK_TOP();
		f->edge = E->enext[f->edge];
		if (f->edge == f->edge0) HSTK_POP();
	}
	if (begin != QH_INVALID) qh_set_opposite(s, begin, prev);
	#undef HSTK_PUSH
	#undef HSTK_POP
	#undef HSTK_TOP
	#undef HSTK_EMPTY
}

// -----------------------------------------------------------------------------
// Create a triangle face from the eye vertex to a horizon edge,
// twin bottom edge with horizon's old twin, return side edge (head->eye).

static int qh_add_adjoining_face(QH_State* s, int eye, int horizon_edge)
{
	int tail = qh_edge_tail(s, horizon_edge);
	int head = qh_edge_head(s, horizon_edge);
	int fi = qh_create_triangle(s, eye, tail, head);
	// e1 (tail->head) twins with the old opposite of the horizon edge.
	qh_set_opposite(s, qh_get_edge(s, fi, 1), s->edges.etwin[horizon_edge]);
	// Return e2 (head->eye) as the side edge.
	return qh_get_edge(s, fi, 2);
}

// -----------------------------------------------------------------------------
// Splice two edges that should be consecutive after removing a shared boundary.
// Handles degenerate case: if both edges' twins point to the same face, that
// face has become a "fin". Triangle fins are deleted; polygon fins lose a vertex.
// Returns index of discarded face, or QH_INVALID.
//
// Because we store origin explicitly per-edge, the surviving edge must absorb
// the removed edge's origin when splicing.

static int qh_connect_half_edges(QH_State* s, int ep, int en)
{
	int opp_p = qh_opp_face(s, ep);
	int opp_n = qh_opp_face(s, en);
	int this_face = s->edges.eface[en];

	if (opp_p == opp_n) {
		// Degenerate: both edges border the same adjacent face.
		int fin_face = opp_n;
		int discarded = QH_INVALID;
		int new_twin;
		QH_DEBUG(fprintf(stderr, "[qh] connect_degenerate: ep=%d en=%d fin=%d nverts=%d this=%d\n",
			ep, en, fin_face, qh_face_vert_count(s, fin_face), this_face));

		if (ep == s->faces[this_face].edge)
			s->faces[this_face].edge = en;

		if (qh_face_vert_count(s, fin_face) == 3) {
			// Triangle collapse: delete the fin face entirely.
			// Find the third edge (not twin of ep or en) whose twin survives.
			int ht = s->edges.etwin[en];
			int hp = s->edges.etwin[ep];
			int e3 = s->edges.enext[ht];
			if (e3 == hp) e3 = s->edges.eprev[ht];
			new_twin = s->edges.etwin[e3];
			s->faces[fin_face].mark = QH_DELETED;
			discarded = fin_face;
		} else {
			// Polygon fin should have been pre-merged by the caller.
			// If we still get here, use normal linkage as fallback.
			s->edges.enext[ep] = en;
			s->edges.eprev[en] = ep;
			return QH_INVALID;
		}

		// Remove ep from this face's loop; en absorbs ep's origin.
		s->edges.eorigin[en] = s->edges.eorigin[ep];
		int pp = s->edges.eprev[ep];
		QH_DEBUG(fprintf(stderr, "[qh]   splice out ep=%d: pp=%d->en=%d (en.next=%d)\n",
			ep, pp, en, s->edges.enext[en]));
		s->edges.eprev[en] = pp;
		s->edges.enext[pp] = en;

		qh_set_opposite(s, en, new_twin);
		qh_recompute_face(s, s->edges.eface[new_twin]);
		return discarded;
	} else {
		// Normal case: just link consecutively.
		s->edges.enext[ep] = en;
		s->edges.eprev[en] = ep;
		return QH_INVALID;
	}
}

// -----------------------------------------------------------------------------
// Merge the adjacent face (across a shared edge) into this face.

static int qh_merge_adjacent_face(QH_State* s, int hedge_adj, int discarded[], int* unclaimed)
{
	int ofi = qh_opp_face(s, hedge_adj);
	int tfi = s->edges.eface[hedge_adj];
	int nd = 0;
	discarded[nd++] = ofi;
	QH_DEBUG(fprintf(stderr, "[qh] merge: face %d (mark=%d) absorbed into %d (mark=%d) via edge %d\n",
		ofi, s->faces[ofi].mark, tfi, s->faces[tfi].mark, hedge_adj));
	s->faces[ofi].mark = QH_DELETED;

	int ho = s->edges.etwin[hedge_adj];
	int ap = s->edges.eprev[hedge_adj], an = s->edges.enext[hedge_adj];
	int op = s->edges.eprev[ho], on = s->edges.enext[ho];

	// Walk past multiply-shared edges.
	while (qh_opp_face(s, ap) == ofi) { ap = s->edges.eprev[ap]; on = s->edges.enext[on]; }
	while (qh_opp_face(s, an) == ofi) { op = s->edges.eprev[op]; an = s->edges.enext[an]; }

	// Reassign the absorbed face's edges to this face.
	for (int e = on; e != s->edges.enext[op]; e = s->edges.enext[e]) s->edges.eface[e] = tfi;

	if (hedge_adj == s->faces[tfi].edge) s->faces[tfi].edge = an;

	int df = qh_connect_half_edges(s, op, an);
	if (df != QH_INVALID) discarded[nd++] = df;
	df = qh_connect_half_edges(s, ap, on);
	if (df != QH_INVALID) discarded[nd++] = df;

	// Inherit maxoutside from absorbed face.
	if (s->faces[ofi].maxoutside > s->faces[tfi].maxoutside)
		s->faces[tfi].maxoutside = s->faces[ofi].maxoutside;

	// Reassign all edges on the surviving face's loop to tfi, and remove
	// self-edges (both half-edges on the same face) left by degenerate merges.
	{
		int e = s->faces[tfi].edge, start = e, n = 0;
		do {
			n++;
			QH_ASSERT(n < 1000, "face loop corrupted after merge");
			s->edges.eface[e] = tfi;
			e = s->edges.enext[e];
		} while (e != start);

		// Remove stale boundary edges: edges whose twins are on deleted faces
		// are remnants of the shared boundary that cascade splicing didn't clean.
		{
			int cleaned = 1;
			while (cleaned) {
				cleaned = 0;
				e = s->faces[tfi].edge;
				start = e;
				do {
					int tw = s->edges.etwin[e];
					int twf = s->edges.eface[tw];
					if (s->faces[twf].mark == QH_DELETED || twf == tfi) {
						// Stale or self-edge: splice out.
						int p = s->edges.eprev[e];
						int nx = s->edges.enext[e];
						s->edges.enext[p] = nx;
						s->edges.eprev[nx] = p;
						s->edges.eorigin[nx] = s->edges.eorigin[e];
						if (s->faces[tfi].edge == e)
							s->faces[tfi].edge = nx;
						cleaned = 1;
						break;
					}
					e = s->edges.enext[e];
				} while (e != start);
			}
		}

		QH_DEBUG({
			e = s->faces[tfi].edge; start = e; n = 0;
			do { n++; e = s->edges.enext[e]; } while (e != start);
			fprintf(stderr, "[qh]   face %d loop (%d edges):", tfi, n);
			e = start;
			do { fprintf(stderr, " %d", e); e = s->edges.enext[e]; } while (e != start);
			fprintf(stderr, "\n");
		});
	}

	qh_recompute_face(s, tfi);
	return nd;
}

// -----------------------------------------------------------------------------
// Scan face edges for non-convex neighbors and merge one pair. Returns 1 if merged.

#define QH_MERGE_LARGE 1
#define QH_MERGE_ANY   2

static int qh_do_adjacent_merge(QH_State* s, int fi, int type, int* unclaimed)
{
	if (s->faces[fi].mark != QH_VISIBLE) return 0;
	int hedge = s->faces[fi].edge;
	int convex = 1;
	float tol = s->epsilon;
	HullPlane p_fi = s->faces[fi].plane;
	do {
		int ofi = qh_opp_face(s, hedge);
		if (ofi == fi || s->faces[ofi].mark == QH_DELETED) {
			hedge = s->edges.enext[hedge]; continue;
		}

		// Inline opp_face_dist: distance of ofi's centroid to fi's plane, and vice versa.
		float d_fwd = dot(p_fi.normal, s->faces[ofi].centroid) - p_fi.offset;

		int merge = 0;
		if (type == QH_MERGE_ANY) {
			if (d_fwd > -tol) merge = 1;
			else { float d_rev = dot(s->faces[ofi].plane.normal, s->faces[fi].centroid) - s->faces[ofi].plane.offset; if (d_rev > -tol) merge = 1; }
		} else {
			if (s->faces[fi].area > s->faces[ofi].area) {
				if (d_fwd > -tol) merge = 1;
				else { float d_rev = dot(s->faces[ofi].plane.normal, s->faces[fi].centroid) - s->faces[ofi].plane.offset; if (d_rev > -tol) convex = 0; }
			} else {
				float d_rev = dot(s->faces[ofi].plane.normal, s->faces[fi].centroid) - s->faces[ofi].plane.offset;
				if (d_rev > -tol) merge = 1;
				else if (d_fwd > -tol) convex = 0;
			}
		}
		if (merge) {
			QH_DEBUG(fprintf(stderr, "[qh] do_adjacent_merge: fi=%d hedge=%d ofi=%d hedge.face=%d\n",
				fi, hedge, ofi, s->edges.eface[hedge]));
			// Check if the merge will hit a polygon fin degenerate.
			// If so, skip it to avoid infinite loops.
			int will_poly_fin = 0;
			{
				int ho = s->edges.etwin[hedge];
				int ap = s->edges.eprev[hedge], an = s->edges.enext[hedge];
				int op = s->edges.eprev[ho], on = s->edges.enext[ho];
				while (qh_opp_face(s, ap) == ofi) { ap = s->edges.eprev[ap]; on = s->edges.enext[on]; }
				while (qh_opp_face(s, an) == ofi) { op = s->edges.eprev[op]; an = s->edges.enext[an]; }
				if ((qh_opp_face(s, op) == qh_opp_face(s, an) && qh_face_vert_count(s, qh_opp_face(s, op)) >= 4) ||
				    (qh_opp_face(s, ap) == qh_opp_face(s, on) && qh_face_vert_count(s, qh_opp_face(s, ap)) >= 4))
					will_poly_fin = 1;
			}
			if (will_poly_fin) {
				hedge = s->edges.enext[hedge];
				continue;
			}
			int disc[3];
			int nd = qh_merge_adjacent_face(s, hedge, disc, unclaimed);
			for (int i = 0; i < nd; i++) qh_delete_face_points(s, disc[i], fi, unclaimed);
			return 1;
		}
		hedge = s->edges.enext[hedge];
	} while (hedge != s->faces[fi].edge);
	if (!convex) s->faces[fi].mark = QH_NON_CONVEX;
	return 0;
}

// -----------------------------------------------------------------------------
// Build cone of triangle faces from eye vertex to horizon, link adjacent sides.

static void qh_face_list_add(QH_State* s, QH_FaceList* list, int fi)
{
	if (list->first == QH_INVALID) { list->first = fi; } else { s->faces[list->last].next = fi; }
	s->faces[fi].next = QH_INVALID;
	list->last = fi;
}

// Reassign orphaned conflict vertices from the unclaimed list to new faces.
static void qh_resolve_unclaimed(QH_State* s, QH_FaceList* nf, int* unclaimed)
{
	QH_Verts* V = &s->verts;
	float eps = s->epsilon, eps1000 = 1000 * eps;
	int v = *unclaimed;
	while (v != QH_INVALID) {
		int nxt = V->cnext[v];
		float vx = V->x[v], vy = V->y[v], vz = V->z[v];
		float bd = eps; int bf = QH_INVALID;
		for (int fi = nf->first; fi != QH_INVALID; fi = s->faces[fi].next) {
			if (s->faces[fi].mark == QH_VISIBLE) {
				HullPlane p = s->faces[fi].plane;
				float d = p.normal.x * vx + p.normal.y * vy + p.normal.z * vz - p.offset;
				if (d > bd) { bd = d; bf = fi; }
				if (bd > eps1000) break;
			}
		}
		if (bf == QH_INVALID) {
			float closest = -1e18f;
			for (int fi = 0; fi < asize(s->faces); fi++) {
				if (s->faces[fi].mark == QH_DELETED) continue;
				HullPlane p = s->faces[fi].plane;
				float d = p.normal.x * vx + p.normal.y * vy + p.normal.z * vz - p.offset;
				if (d > closest) { closest = d; bf = fi; }
			}
			if (closest < -100.0f * eps) bf = QH_INVALID;
		}
		if (bf != QH_INVALID) qh_conflict_add(s, bf, v);
		v = nxt;
	}
	*unclaimed = QH_INVALID;
}

// Debug: validate mesh topology of all live faces.
static void qh_validate_mesh(QH_State* s, const char* ctx)
{
	for (int fi = 0; fi < asize(s->faces); fi++) {
		if (s->faces[fi].mark != QH_VISIBLE) continue;
		int e = s->faces[fi].edge, start = e, n = 0;
		do {
			if (n > 500) {
				fprintf(stderr, "qh_validate: infinite loop on face %d at %s\n", fi, ctx);
				exit(-1);
			}
			// Check twin reciprocity.
			int tw = s->edges.etwin[e];
			if (s->edges.etwin[tw] != e) {
				fprintf(stderr, "qh_validate: edge %d twin=%d but twin's twin=%d at %s\n",
					e, tw, s->edges.etwin[tw], ctx);
				exit(-1);
			}
			// Check twin's face is alive.
			int twf = s->edges.eface[tw];
			if (s->faces[twf].mark != QH_VISIBLE) {
				fprintf(stderr, "qh_validate: edge %d (face %d) twin=%d on deleted face %d at %s\n",
					e, fi, tw, twf, ctx);
				exit(-1);
			}
			// Check edge belongs to this face.
			if (s->edges.eface[e] != fi) {
				fprintf(stderr, "qh_validate: edge %d face=%d expected %d at %s\n",
					e, s->edges.eface[e], fi, ctx);
				exit(-1);
			}
			n++;
			e = s->edges.enext[e];
		} while (e != start);
	}

	// Euler check: V - E/2 + F = 2 on the live mesh.
	{
		CK_DYNA int* vremap = NULL;
		afit(vremap, s->verts.count);
		for (int i = 0; i < s->verts.count; i++) apush(vremap, 0);
		int nf = 0, ne = 0;
		for (int fi = 0; fi < asize(s->faces); fi++) {
			if (s->faces[fi].mark != QH_VISIBLE) continue;
			nf++;
			int e = s->faces[fi].edge, start = e;
			do {
				ne++;
				vremap[s->edges.eorigin[e]] = 1;
				e = s->edges.enext[e];
			} while (e != start);
		}
		int nv = 0;
		for (int i = 0; i < s->verts.count; i++) nv += vremap[i];
		afree(vremap);
		if (nv - ne/2 + nf != 2) {
			fprintf(stderr, "qh_validate: Euler FAIL V=%d E=%d F=%d (V-E/2+F=%d) at %s\n",
				nv, ne, nf, nv - ne/2 + nf, ctx);
			exit(-1);
		}
	}
}

// -----------------------------------------------------------------------------
// Process one conflict vertex: compute horizon, build cone, merge, reassign.

static void qh_add_point(QH_State* s, int eye, int eye_face)
{
	int unclaimed = QH_INVALID;

	qh_conflict_remove(s, eye_face, eye);

	QH_FaceList nf;
	qh_horizon_and_cone(s, eye, eye_face, &nf, &unclaimed);

	// Two-pass merge: first larger-face-biased, then mutual non-convexity.
	// Each successful merge reduces the live face count by 1, so the total
	// number of merges across all faces is bounded by the initial face count.
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_VISIBLE)
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_LARGE, &unclaimed));
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_NON_CONVEX) {
			s->faces[fi].mark = QH_VISIBLE;
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_ANY, &unclaimed));
		}
	// Global topology fixup: cascade merges + degenerate connect can create
	// face pairs sharing >1 edge, violating manifold topology. Sweep only
	// new faces (and faces touching them) for multi-adjacent pairs.
	// Detect duplicates with a generation-stamped marker array.
	{
		// Ensure marker array is large enough.
		int nfaces = asize(s->faces);
		static int* qh_face_seen = NULL;
		static int qh_face_seen_cap = 0;
		static int qh_face_gen = 0;
		if (nfaces > qh_face_seen_cap) {
			free(qh_face_seen);
			qh_face_seen_cap = nfaces * 2;
			qh_face_seen = (int*)calloc(qh_face_seen_cap, sizeof(int));
			qh_face_gen = 0;
		}
		int fixed = 1;
		while (fixed) {
			fixed = 0;
			for (int fi = nf.first; fi != QH_INVALID && !fixed; fi = s->faces[fi].next) {
				if (s->faces[fi].mark != QH_VISIBLE) continue;
				qh_face_gen++;
				if (qh_face_gen == 0) { memset(qh_face_seen, 0, qh_face_seen_cap * sizeof(int)); qh_face_gen = 1; }
				int hedge = s->faces[fi].edge;
				do {
					int ofi = qh_opp_face(s, hedge);
					if (ofi == fi || s->faces[ofi].mark != QH_VISIBLE) {
						hedge = s->edges.enext[hedge];
						continue;
					}
					if (qh_face_seen[ofi] == qh_face_gen) {
						// Multi-adjacent: ofi seen twice.
						int small = (s->faces[fi].num_verts <= s->faces[ofi].num_verts) ? fi : ofi;
						int large = (small == fi) ? ofi : fi;
						int me = s->faces[large].edge;
						while (qh_opp_face(s, me) != small) me = s->edges.enext[me];
						int disc[3];
						int nd = qh_merge_adjacent_face(s, me, disc, &unclaimed);
						for (int k = 0; k < nd; k++)
							qh_delete_face_points(s, disc[k], large, &unclaimed);
						nfaces = asize(s->faces);
						if (nfaces > qh_face_seen_cap) {
							free(qh_face_seen);
							qh_face_seen_cap = nfaces * 2;
							qh_face_seen = (int*)calloc(qh_face_seen_cap, sizeof(int));
							qh_face_gen = 0;
						}
						fixed = 1;
						break;
					}
					qh_face_seen[ofi] = qh_face_gen;
					hedge = s->edges.enext[hedge];
				} while (hedge != s->faces[fi].edge);
			}
		}
	}

	qh_resolve_unclaimed(s, &nf, &unclaimed);
}

// Shared Newell plane computation: walks a half-edge face loop starting from `start_edge`,
// reading vertex positions from SoA float arrays. Used by both quickhull output and
// hull_face_extension_build for bitwise-deterministic planes.
static HullPlane hull_newell_plane(const uint16_t* edge_next, const uint16_t* edge_origin, const float* vx, const float* vy, const float* vz, int start_edge, v3 hull_centroid)
{
	float fnx = 0, fny = 0, fnz = 0, fcx = 0, fcy = 0, fcz = 0; int cnt = 0;
	int e = start_edge;
	do {
		int vi = edge_origin[e], vn = edge_origin[edge_next[e]];
		float curx = vx[vi], cury = vy[vi], curz = vz[vi], nxtx = vx[vn], nxty = vy[vn], nxtz = vz[vn];
		fnx += (cury - nxty) * (curz + nxtz); fny += (curz - nxtz) * (curx + nxtx); fnz += (curx - nxtx) * (cury + nxty);
		fcx += curx; fcy += cury; fcz += curz; cnt++;
		e = edge_next[e];
	} while (e != start_edge);
	float a = sqrtf(fnx*fnx + fny*fny + fnz*fnz);
	if (a > 0) { float inv = 1.0f / a; fnx *= inv; fny *= inv; fnz *= inv; }
	float inv_c = 1.0f / cnt; fcx *= inv_c; fcy *= inv_c; fcz *= inv_c;
	v3 fn = V3(fnx, fny, fnz), fc = V3(fcx, fcy, fcz);
	if (dot(fn, sub(fc, hull_centroid)) < 0) fn = neg(fn);
	return (HullPlane){ fn, dot(fn, fc) };
}

// -----------------------------------------------------------------------------
// Convert internal state to output Hull with uint16_t-indexed arrays.

static Hull* qh_build_output(QH_State* s, const v3* all_points, int all_count)
{
	CK_DYNA int* live = NULL;
	for (int i = 0; i < asize(s->faces); i++)
		if (s->faces[i].mark == QH_VISIBLE) apush(live, i);

	// Validate face loops before output -- detect broken topology from merges.
	QH_DEBUG({
		int total_edges = s->edges.count;
		for (int i = 0; i < asize(live); i++) {
			int fi = live[i];
			int e = s->faces[fi].edge, start = e, steps = 0;
			do {
				steps++;
				if (steps > total_edges) {
					fprintf(stderr, "qh_build_output: face %d (live[%d]) edge loop did not close after %d steps (start edge %d, cur edge %d)\n",
						fi, i, steps, start, e);
					assert(0 && "qh_build_output: broken face edge loop");
				}
				e = s->edges.enext[e];
			} while (e != start);
		}
	});

	// Count verts/edges first, then allocate hull arrays directly.
	CK_DYNA int* vremap = NULL;
	afit_set(vremap, s->verts.count);
	memset(vremap, -1, s->verts.count * sizeof(int));
	CK_DYNA int* eremap = NULL;
	afit_set(eremap, s->edges.count);
	memset(eremap, -1, s->edges.count * sizeof(int));
	int vc = 0, ec = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do {
			int vi = s->edges.eorigin[e];
			if (vremap[vi] < 0) vremap[vi] = vc++;
			eremap[e] = ec++;
			e = s->edges.enext[e];
		} while (e != start);
	}

	// Allocate hull and its arrays directly (no temp arrays).
	int nlive = asize(live);
	Hull* h = CK_ALLOC(sizeof(Hull));
	h->vert_count = vc; h->edge_count = ec; h->face_count = nlive;
	h->epsilon = s->epsilon;
	v3* vcp = CK_ALLOC(sizeof(v3) * vc); h->verts = vcp;
	uint16_t* et = CK_ALLOC(sizeof(uint16_t) * ec); h->edge_twin = et;
	uint16_t* en = CK_ALLOC(sizeof(uint16_t) * ec); h->edge_next = en;
	uint16_t* eo = CK_ALLOC(sizeof(uint16_t) * ec); h->edge_origin = eo;
	uint16_t* ef = CK_ALLOC(sizeof(uint16_t) * ec); h->edge_face = ef;
	HullFace* fcp = CK_ALLOC(sizeof(HullFace) * nlive); h->faces = fcp;
	HullPlane* pcp = CK_ALLOC(sizeof(HullPlane) * nlive); h->planes = pcp;

	// Fill vertex positions (second pass over vremap to fill only assigned slots).
	memset(vremap, -1, s->verts.count * sizeof(int));
	vc = 0;
	for (int i = 0; i < nlive; i++) {
		int e = s->faces[live[i]].edge, start = e;
		do {
			int vi = s->edges.eorigin[e];
			if (vremap[vi] < 0) { vremap[vi] = vc; vcp[vc++] = qh_vert_pos(&s->verts, vi); }
			e = s->edges.enext[e];
		} while (e != start);
	}

	// Build SoA vertex arrays for SIMD support queries.
	{
		int padded = (vc + 3) & ~3; // pad to multiple of 4 for SIMD
		float* soa = (float*)CK_ALLOC_ALIGNED(padded * 3 * sizeof(float), 16);
		float* sx = soa, *sy = soa + padded, *sz = soa + padded * 2;
		for (int i = 0; i < vc; i++) { sx[i] = vcp[i].x; sy[i] = vcp[i].y; sz[i] = vcp[i].z; }
		for (int i = vc; i < padded; i++) { sx[i] = sx[0]; sy[i] = sy[0]; sz[i] = sz[0]; } // pad with v0
		h->soa_verts = soa;
	}

	// Fill edges and faces directly.
	for (int i = 0; i < nlive; i++) {
		int e = s->faces[live[i]].edge, start = e;
		do {
			int o = eremap[e];
			en[o] = (uint16_t)eremap[s->edges.enext[e]];
			et[o] = (uint16_t)eremap[s->edges.etwin[e]];
			eo[o] = (uint16_t)vremap[s->edges.eorigin[e]];
			ef[o] = (uint16_t)i;
			e = s->edges.enext[e];
		} while (e != start);
		fcp[i] = (HullFace){ .edge = (uint16_t)eremap[s->faces[live[i]].edge] };
		pcp[i] = s->faces[live[i]].plane;
	}

	v3 centroid = V3(0, 0, 0);
	for (int i = 0; i < h->vert_count; i++) centroid = add(centroid, vcp[i]);
	centroid = scale(centroid, 1.0f / h->vert_count);
	h->centroid = centroid;

	// Fix inward-facing normals: degenerate merges can produce Newell normals
	// that point toward the hull interior. Flip them before widening.
	for (int i = 0; i < h->face_count; i++) {
		v3 fc = s->faces[live[i]].centroid;
		if (dot(pcp[i].normal, sub(fc, centroid)) < 0) {
			pcp[i].normal = scale(pcp[i].normal, -1.0f);
			pcp[i].offset = -pcp[i].offset;
		}
	}

	// Post-build plane widening: widen each plane so ALL original input points
	// lie on or behind the plane. Must use original points (not welded subset)
	// since welded-away points may lie outside the welded hull.
	// SoA transpose for SIMD-friendly inner loop.
	h->maxoutside = 0;
	int fc = h->face_count;
	int soa_bytes = fc * 4 * sizeof(float);
	float* nx = (soa_bytes <= 4096) ? (float*)_alloca(soa_bytes) : (float*)CK_ALLOC(soa_bytes);
	float* ny = nx + fc, *nz = ny + fc, *max_ds = nz + fc;
	for (int i = 0; i < fc; i++) { nx[i] = pcp[i].normal.x; ny[i] = pcp[i].normal.y; nz[i] = pcp[i].normal.z; max_ds[i] = pcp[i].offset; }
	int fc4 = fc & ~3;
	for (int pi = 0; pi < all_count; pi++) {
		simd4f vpx = simd_set1(all_points[pi].x), vpy = simd_set1(all_points[pi].y), vpz = simd_set1(all_points[pi].z);
		int i = 0;
		for (; i < fc4; i += 4) {
			simd4f d = simd_add(simd_add(simd_mul(simd_load(nx + i), vpx), simd_mul(simd_load(ny + i), vpy)), simd_mul(simd_load(nz + i), vpz));
			simd_store(max_ds + i, simd_max(simd_load(max_ds + i), d));
		}
		for (; i < fc; i++) {
			float d = nx[i] * all_points[pi].x + ny[i] * all_points[pi].y + nz[i] * all_points[pi].z;
			if (d > max_ds[i]) max_ds[i] = d;
		}
	}
	for (int i = 0; i < fc; i++) {
		float widen = max_ds[i] - pcp[i].offset;
		if (widen > h->maxoutside) h->maxoutside = widen;
		pcp[i].offset = max_ds[i];
	}
	if (soa_bytes > 4096) free(nx);

	afree(live); afree(vremap); afree(eremap);
	return h;
}

// -----------------------------------------------------------------------------
// Public API.

Hull* quickhull(const v3* points, int count)
{
	qh_dbg_points = points;
	qh_dbg_count = count;
	QH_ASSERT(count >= 4, "quickhull needs at least 4 points");

	QH_State state = {0};
	state.edges.free_head = QH_INVALID;
	state.face_free = QH_INVALID;

	// Compute epsilon = 3 * (max|x| + max|y| + max|z|) * FLT_EPSILON.
	float mx = 0, my = 0, mz = 0;
	for (int i = 0; i < count; i++) {
		float ax = fabsf(points[i].x), ay = fabsf(points[i].y), az = fabsf(points[i].z);
		if (ax > mx) mx = ax; if (ay > my) my = ay; if (az > mz) mz = az;
	}
	state.epsilon = 3.0f * (mx + my + mz) * FLT_EPSILON;

	// Weld near-duplicate vertices to prevent degenerate triangles that corrupt
	// topology during merge. Weld distance is epsilon (the floating-point
	// precision floor for this input extent).
	float weld_dist2 = state.epsilon * state.epsilon;
	qh_verts_reserve(&state.verts, count);
	simd4f vwd2 = simd_set1(weld_dist2);
	for (int i = 0; i < count; i++) {
		v3 p = points[i];
		simd4f vpx = simd_set1(p.x), vpy = simd_set1(p.y), vpz = simd_set1(p.z);
		int dup = 0, n = state.verts.count, n4 = n & ~3;
		int j = 0;
		for (; j < n4; j += 4) {
			simd4f dx = simd_sub(vpx, simd_load(state.verts.x + j));
			simd4f dy = simd_sub(vpy, simd_load(state.verts.y + j));
			simd4f dz = simd_sub(vpz, simd_load(state.verts.z + j));
			simd4f d2 = simd_add(simd_add(simd_mul(dx, dx), simd_mul(dy, dy)), simd_mul(dz, dz));
			if (simd_movemask(simd_cmple(d2, vwd2))) { dup = 1; break; }
		}
		for (; !dup && j < n; j++) {
			if (len2(sub(p, qh_vert_pos(&state.verts, j))) <= weld_dist2) { dup = 1; break; }
		}
		if (!dup) qh_verts_push(&state.verts, p);
	}
	int welded_count = state.verts.count;

	if (welded_count < 4) {
		qh_verts_free(&state.verts); qh_edges_free(&state.edges); afree(state.faces); afree(state.conflict_faces);
		return NULL;
	}

	if (!qh_build_simplex(&state, welded_count)) {
		qh_verts_free(&state.verts); qh_edges_free(&state.edges); afree(state.faces); afree(state.conflict_faces);
		return NULL;
	}

	for (int iter = 0; iter < welded_count * 4; iter++) {
		int fi; int vi = qh_next_conflict(&state, &fi);
		if (vi == QH_INVALID) break;
		qh_add_point(&state, vi, fi);
		QH_DEBUG({
			char buf[64];
			snprintf(buf, sizeof(buf), "after adding vertex %d (iter %d)", vi, iter);
			qh_validate_mesh(&state, buf);
		});
	}

	Hull* result = qh_build_output(&state, points, count);
	qh_verts_free(&state.verts); qh_edges_free(&state.edges); afree(state.faces); afree(state.conflict_faces);
	return result;
}
// Build CSR vertex adjacency from half-edge SoA arrays.
// For each vertex, collects all neighbor vertices reachable via edges.
// vert_edge[v] = first edge originating from v (or -1).
// Walk: twin of that edge gives neighbor, then edge_next[twin] gives next
// edge around v, repeat until back to start.
static int hull_build_csr(const Hull* hull, int max_verts, void* offsets_out, int offset_size, void* neighbors_out, int neighbor_size, int* out_neighbor_total)
{
	int nv = hull->vert_count;
	if (nv > max_verts) return -1;

	// Build vert_edge: first edge index for each vertex.
	int vert_edge[1024];
	assert(nv <= 1024);
	for (int i = 0; i < nv; i++) vert_edge[i] = -1;
	for (int i = 0; i < hull->edge_count; i++)
		if (vert_edge[hull->edge_origin[i]] < 0) vert_edge[hull->edge_origin[i]] = i;

	// Count neighbors per vertex, then fill CSR.
	int total = 0;
	for (int v = 0; v < nv; v++) {
		if (offset_size == 1) ((uint8_t*)offsets_out)[v] = (uint8_t)total;
		else ((uint16_t*)offsets_out)[v] = (uint16_t)total;

		int e = vert_edge[v];
		if (e < 0) continue;
		int start = e;
		do {
			int nb = hull->edge_origin[hull->edge_twin[e]];
			if (neighbor_size == 1) ((uint8_t*)neighbors_out)[total] = (uint8_t)nb;
			else ((uint16_t*)neighbors_out)[total] = (uint16_t)nb;
			total++;
			e = hull->edge_next[hull->edge_twin[e]];
		} while (e != start);
	}
	if (offset_size == 1) ((uint8_t*)offsets_out)[nv] = (uint8_t)total;
	else ((uint16_t*)offsets_out)[nv] = (uint16_t)total;
	*out_neighbor_total = total;
	return 0;
}

// -----------------------------------------------------------------------------
// Compact hull converters.
// CompactHull32: fixed-size (<=32 verts), suitable for stack/inline storage.
// CompactHull: heap-allocated, supports arbitrary vert counts.

// Convert a full Hull to a CompactHull32 (fixed-size, <=32 verts).
int compact_hull32_from_hull(CompactHull32* out, const Hull* hull)
{
	if (hull->vert_count > COMPACT_HULL32_MAX_VERTS) return -1;
	*out = (CompactHull32){0};
	out->vert_count = (uint8_t)hull->vert_count;
	out->centroid = hull->centroid;

	// Copy vertex positions into SoA.
	for (int i = 0; i < hull->vert_count; i++) {
		out->verts_x[i] = hull->verts[i].x;
		out->verts_y[i] = hull->verts[i].y;
		out->verts_z[i] = hull->verts[i].z;
	}

	int total = 0;
	if (hull_build_csr(hull, COMPACT_HULL32_MAX_VERTS, out->offsets, 1, out->neighbors, 1, &total) < 0) return -1;
	if (total > COMPACT_HULL32_MAX_NEIGHBORS) return -1;
	out->neighbor_total = (uint8_t)total;

	// Copy face planes for bitwise-deterministic reconstruction.
	if (hull->face_count > COMPACT_HULL32_MAX_FACES) return -1;
	out->face_count = (uint8_t)hull->face_count;
	memcpy(out->planes, hull->planes, hull->face_count * sizeof(HullPlane));
	return 0;
}

// Convert a full Hull to a CompactHull (heap-allocated, arbitrary vert count).
int compact_hull_from_hull(CompactHull* out, const Hull* hull)
{
	*out = (CompactHull){0};
	int nv = hull->vert_count;
	out->vert_count = (uint16_t)nv;
	out->centroid = hull->centroid;

	out->offsets = (uint16_t*)CK_ALLOC((nv + 1) * sizeof(uint16_t));
	// Allocate neighbors conservatively (edge_count = total half-edges = 2 * undirected edges).
	out->neighbors = (uint16_t*)CK_ALLOC(hull->edge_count * sizeof(uint16_t));

	int padded = (nv + 3) & ~3;
	out->verts_x = (float*)CK_ALLOC_ALIGNED(padded * sizeof(float), 16);
	out->verts_y = (float*)CK_ALLOC_ALIGNED(padded * sizeof(float), 16);
	out->verts_z = (float*)CK_ALLOC_ALIGNED(padded * sizeof(float), 16);
	for (int i = 0; i < nv; i++) { out->verts_x[i] = hull->verts[i].x; out->verts_y[i] = hull->verts[i].y; out->verts_z[i] = hull->verts[i].z; }
	for (int i = nv; i < padded; i++) { out->verts_x[i] = out->verts_x[0]; out->verts_y[i] = out->verts_y[0]; out->verts_z[i] = out->verts_z[0]; }

	int total = 0;
	hull_build_csr(hull, 65535, out->offsets, 2, out->neighbors, 2, &total);
	out->neighbor_total = (uint16_t)total;
	// Planes left NULL -- caller uses compact_hull_attach_planes() on demand.
	return 0;
}

void compact_hull_attach_planes(CompactHull* ch, const Hull* hull)
{
	if (ch->planes) CK_FREE(ch->planes);
	ch->face_count = (uint16_t)hull->face_count;
	ch->planes = (HullPlane*)CK_ALLOC(hull->face_count * sizeof(HullPlane));
	memcpy(ch->planes, hull->planes, hull->face_count * sizeof(HullPlane));
}

void compact_hull_free(CompactHull* ch)
{
	if (!ch) return;
	CK_FREE(ch->offsets);
	CK_FREE(ch->neighbors);
	if (ch->planes) CK_FREE(ch->planes);
	CK_FREE_ALIGNED(ch->verts_x);
	CK_FREE_ALIGNED(ch->verts_y);
	CK_FREE_ALIGNED(ch->verts_z);
	*ch = (CompactHull){0};
}

// Validate compact hull round-trip: build face extension, check Euler formula,
// twin reciprocity, plane identity, and face loop closure.
int compact_hull_validate_roundtrip(const CompactHull* ch)
{
	HullFaceExtension ext;
	if (hull_face_extension_build(&ext, ch) != 0) {
		fprintf(stderr, "compact_hull_validate: extension build failed\n");
		return -1;
	}
	int ok = 1;

	// Euler: V - E/2 + F = 2.
	if (ch->vert_count - ext.edge_count / 2 + ext.face_count != 2) {
		fprintf(stderr, "compact_hull_validate: Euler FAIL V=%d E=%d F=%d\n", ch->vert_count, ext.edge_count, ext.face_count);
		ok = 0;
	}

	// Twin reciprocity.
	for (int e = 0; e < ext.edge_count && ok; e++) {
		if (ext.edge_twin[ext.edge_twin[e]] != e) {
			fprintf(stderr, "compact_hull_validate: twin reciprocity FAIL at edge %d\n", e);
			ok = 0;
		}
	}

	// If planes are attached, check bitwise identity.
	if (ch->planes && ch->face_count > 0) {
		if (ext.face_count != ch->face_count) {
			fprintf(stderr, "compact_hull_validate: face count mismatch ext=%d stored=%d\n", ext.face_count, ch->face_count);
			ok = 0;
		}
		for (int f = 0; f < ext.face_count && ok; f++) {
			int found = 0;
			for (int g = 0; g < ch->face_count; g++) {
				if (memcmp(&ext.planes[f], &ch->planes[g], sizeof(HullPlane)) == 0) { found = 1; break; }
			}
			if (!found) {
				fprintf(stderr, "compact_hull_validate: plane %d has no bitwise match\n", f);
				ok = 0;
			}
		}
	}

	// Face loops close.
	for (int f = 0; f < ext.face_count && ok; f++) {
		int e = ext.faces[f].edge, n = 0;
		do { n++; e = ext.edge_next[e]; if (n > 1000) break; } while (e != ext.faces[f].edge);
		if (n > 1000 || n < 3) {
			fprintf(stderr, "compact_hull_validate: face %d loop FAIL (%d edges)\n", f, n);
			ok = 0;
		}
	}

	hull_face_extension_free(&ext);
	return ok ? 0 : -1;
}

// -----------------------------------------------------------------------------
// Face extension: reconstruct half-edge mesh from CompactHull CSR adjacency.
//
// Algorithm: each directed edge (u,v) in the CSR becomes a half-edge.
// The "next" half-edge around a face is found by: for edge (u,v), the next
// edge on the same face starts at v and goes to the neighbor of v that comes
// AFTER u in v's neighbor list (cyclically). This reconstructs the face loops.
//
// Uses hull_newell_plane() (defined in quickhull.c) for plane recomputation.
int hull_face_extension_build(HullFaceExtension* out, const CompactHull* ch)
{
	*out = (HullFaceExtension){0};
	int nv = ch->vert_count;
	int ne = ch->neighbor_total; // total half-edges

	out->edge_count = (uint16_t)ne;
	out->edge_twin = (uint16_t*)CK_ALLOC(ne * sizeof(uint16_t));
	out->edge_next = (uint16_t*)CK_ALLOC(ne * sizeof(uint16_t));
	out->edge_origin = (uint16_t*)CK_ALLOC(ne * sizeof(uint16_t));
	out->edge_face = (uint16_t*)CK_ALLOC(ne * sizeof(uint16_t));

	// Each CSR entry neighbors[i] at offset i is half-edge i: origin = v where offsets[v] <= i < offsets[v+1].
	// Set edge origins.
	for (int v = 0; v < nv; v++)
		for (int i = ch->offsets[v]; i < ch->offsets[v+1]; i++)
			out->edge_origin[i] = (uint16_t)v;

	// Build twin map: edge i goes from origin[i] to neighbors[i].
	// Its twin goes from neighbors[i] back to origin[i].
	// Find twin by scanning the neighbor list of neighbors[i] for origin[i].
	for (int i = 0; i < ne; i++) {
		int u = out->edge_origin[i];
		int v = ch->neighbors[i];
		int twin = -1;
		for (int j = ch->offsets[v]; j < ch->offsets[v+1]; j++) {
			if (ch->neighbors[j] == u) { twin = j; break; }
		}
		assert(twin >= 0);
		out->edge_twin[i] = (uint16_t)twin;
	}

	// Build next: for edge (u,v), next edge is (v, w) where w is the neighbor
	// of v that comes BEFORE u in v's cyclic neighbor list.
	// This gives the CCW face traversal matching quickhull's winding.
	for (int i = 0; i < ne; i++) {
		int v = ch->neighbors[i]; // head of this edge
		int u = out->edge_origin[i]; // tail of this edge
		// Find u in v's neighbor list, take the PREVIOUS entry cyclically.
		int deg = ch->offsets[v+1] - ch->offsets[v];
		int base = ch->offsets[v];
		for (int k = 0; k < deg; k++) {
			if (ch->neighbors[base + k] == u) {
				int next_k = (k + 1 < deg) ? k + 1 : 0;
				// next edge goes from v to w: find that edge index.
				out->edge_next[i] = (uint16_t)(base + next_k);
				(void)next_k;
				break;
			}
		}
	}

	// Extract faces: walk edge loops, assign face indices.
	// Each unvisited edge starts a new face.
	for (int i = 0; i < ne; i++) out->edge_face[i] = 0xFFFF;
	int face_count = 0;
	CK_DYNA HullFace* faces_arr = NULL;
	CK_DYNA HullPlane* planes_arr = NULL;
	for (int i = 0; i < ne; i++) {
		if (out->edge_face[i] != 0xFFFF) continue;
		int fi = face_count++;

		// First pass: assign face index and find lowest-origin vertex for normalized start.
		int best_start = i;
		uint16_t best_origin = out->edge_origin[i];
		int e = i;
		do {
			out->edge_face[e] = (uint16_t)fi;
			if (out->edge_origin[e] < best_origin) { best_origin = out->edge_origin[e]; best_start = e; }
			e = out->edge_next[e];
		} while (e != i);
		apush(faces_arr, ((HullFace){ .edge = (uint16_t)best_start }));

		if (ch->planes && ch->face_count > 0) {
			// Planes attached: match by Newell normal direction for bitwise determinism.
			HullPlane recomputed = hull_newell_plane(out->edge_next, out->edge_origin, ch->verts_x, ch->verts_y, ch->verts_z, best_start, ch->centroid);
			float best_dot = -2; int best_match = fi;
			for (int g = 0; g < ch->face_count; g++) {
				float d = dot(recomputed.normal, ch->planes[g].normal);
				if (d > best_dot) { best_dot = d; best_match = g; }
			}
			apush(planes_arr, ch->planes[best_match]);
		} else {
			// No planes attached: recompute from Newell method.
			apush(planes_arr, hull_newell_plane(out->edge_next, out->edge_origin, ch->verts_x, ch->verts_y, ch->verts_z, best_start, ch->centroid));
		}
	}

	out->face_count = (uint16_t)face_count;
	out->faces = faces_arr;
	out->planes = planes_arr;
	return 0;
}

void hull_face_extension_free(HullFaceExtension* ext)
{
	if (!ext) return;
	CK_FREE(ext->edge_twin);
	CK_FREE(ext->edge_next);
	CK_FREE(ext->edge_origin);
	CK_FREE(ext->edge_face);
	afree(ext->faces);
	afree(ext->planes);
	*ext = (HullFaceExtension){0};
}

// Build a (non-owning) Hull view from a CompactHull + its face extension.
Hull hull_from_compact(const CompactHull* ch, const HullFaceExtension* ext)
{
	int padded = (ch->vert_count + 3) & ~3;
	Hull h = {0};
	h.centroid = ch->centroid;
	// Build AoS verts from SoA (stack or temp alloc -- caller must manage lifetime).
	// For now, return a Hull that points to NULL verts (caller fills or uses SoA directly).
	h.verts = NULL; // caller should set this if AoS access is needed
	h.soa_verts = ch->verts_x; // contiguous x,y,z (padded layout matches)
	h.edge_twin = ext->edge_twin;
	h.edge_next = ext->edge_next;
	h.edge_origin = ext->edge_origin;
	h.edge_face = ext->edge_face;
	h.faces = ext->faces;
	h.planes = ext->planes;
	h.vert_count = ch->vert_count;
	h.edge_count = ext->edge_count;
	h.face_count = ext->face_count;
	return h;
}
