// See LICENSE for licensing info.
#include <float.h>
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
// We use ORIGIN (tail) vertex convention for half-edges.
// for a triangle created with createTriangle(v0, v1, v2).

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
typedef struct QH_Verts {
	float* x;   // aligned position arrays
	float* y;
	float* z;
	int* cnext; // conflict list next (circular doubly-linked)
	int* cprev; // conflict list prev
	int count;
	int cap;
} QH_Verts;

typedef struct QH_Edge {
	int next, prev, twin;
	int origin, face, mark;
} QH_Edge;

// Hot/cold face split: hot fields (scanned by qh_next_conflict + merge)
// are packed together; cold fields (topology, metadata) accessed only
// during cone creation, merge splicing, and output.
typedef struct QH_Face {
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

typedef struct QH_State {
	QH_Verts verts;
	CK_DYNA QH_Edge*   edges;
	CK_DYNA QH_Face*   faces;
	CK_DYNA int*       conflict_faces; // compact list of faces with non-empty conflict lists
	int edge_free, face_free;
	v3 interior;
	float epsilon;
} QH_State;

typedef struct QH_FaceList {
	int first, last;
} QH_FaceList;


// -----------------------------------------------------------------------------
// Free list management.

static int qh_alloc_edge(QH_State* s)
{
	if (s->edge_free != QH_INVALID) {
		int idx = s->edge_free;
		s->edge_free = s->edges[idx].next;
		return idx; // caller writes all fields
	}
	QH_Edge e = {0};
	apush(s->edges, e);
	return asize(s->edges) - 1;
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
static void qh_recompute_face(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	v3 normal = V3(0,0,0), centroid = V3(0,0,0);
	int count = 0, e = f->edge;
	do {
		v3 cur = qh_vert_pos(&s->verts, s->edges[e].origin);
		v3 nxt = qh_vert_pos(&s->verts, s->edges[s->edges[e].next].origin);
		normal.x += (cur.y - nxt.y) * (cur.z + nxt.z);
		normal.y += (cur.z - nxt.z) * (cur.x + nxt.x);
		normal.z += (cur.x - nxt.x) * (cur.y + nxt.y);
		centroid = add(centroid, cur);
		count++;
		if (count > 1000) { QH_ASSERT(0, "infinite face loop in qh_recompute_face"); }
		e = s->edges[e].next;
	} while (e != f->edge);

	float a = len(normal);
	if (a > 0) normal = scale(normal, 1.0f / a);
	centroid = scale(centroid, 1.0f / count);
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
	int fi = s->edges[ei].face;
	int ofi = s->edges[s->edges[ei].twin].face;
	return qh_face_dist(s, fi, s->faces[ofi].centroid);
}

static int qh_face_vert_count(QH_State* s, int fi)
{
	int e = s->faces[fi].edge, start = e, n = 0;
	do { n++; e = s->edges[e].next; } while (e != start);
	return n;
}

static int qh_opp_face(QH_State* s, int ei) { return s->edges[s->edges[ei].twin].face; }
static int qh_edge_head(QH_State* s, int ei) { return s->edges[s->edges[ei].next].origin; }
static int qh_edge_tail(QH_State* s, int ei) { return s->edges[ei].origin; }

// Get edge at index i from face's he0 (supports negative indices).
static int qh_get_edge(QH_State* s, int fi, int i)
{
	int e = s->faces[fi].edge;
	while (i > 0) { e = s->edges[e].next; i--; }
	while (i < 0) { e = s->edges[e].prev; i++; }
	return e;
}

static void qh_set_opposite(QH_State* s, int a, int b)
{
	s->edges[a].twin = b;
	s->edges[b].twin = a;
}

// Create a triangle face with edges v0->v1->v2->v0.
static int qh_create_triangle(QH_State* s, int v0, int v1, int v2)
{
	int fi = qh_alloc_face(s);
	int e0 = qh_alloc_edge(s), e1 = qh_alloc_edge(s), e2 = qh_alloc_edge(s);
	s->edges[e0] = (QH_Edge){ .next=e1, .prev=e2, .twin=QH_INVALID, .origin=v0, .face=fi };
	s->edges[e1] = (QH_Edge){ .next=e2, .prev=e0, .twin=QH_INVALID, .origin=v1, .face=fi };
	s->edges[e2] = (QH_Edge){ .next=e0, .prev=e1, .twin=QH_INVALID, .origin=v2, .face=fi };
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
		float d = (i==0) ? a.x-b.x : (i==1) ? a.y-b.y : a.z-b.z;
		if (d > best) { best = d; imax = i; }
	}
	if (best <= s->epsilon) return 0;

	int vtx[4];
	vtx[0] = maxV[imax]; vtx[1] = minV[imax];

	// Find vertex farthest from line vtx[0]-vtx[1].
	v3 u01 = norm(sub(qh_vert_pos(&s->verts, vtx[1]), qh_vert_pos(&s->verts, vtx[0])));
	float maxSqr = 0; v3 nrml = V3(0,0,0); vtx[2] = -1;
	for (int i = 0; i < nv; i++) {
		v3 xp = cross(u01, sub(qh_vert_pos(&s->verts, i), qh_vert_pos(&s->verts, vtx[0])));
		float ls = len2(xp);
		if (ls > maxSqr && i != vtx[0] && i != vtx[1]) { maxSqr = ls; vtx[2] = i; nrml = xp; }
	}
	if (vtx[2] < 0 || sqrtf(maxSqr) <= 100*s->epsilon) return 0;
	nrml = norm(nrml);
	nrml = norm(sub(nrml, scale(u01, dot(nrml, u01)))); // orthogonalize

	// Find vertex farthest from plane through vtx[0..2].
	float d0 = dot(qh_vert_pos(&s->verts, vtx[2]), nrml);
	float maxDist = 0; vtx[3] = -1;
	for (int i = 0; i < nv; i++) {
		float d = fabsf(dot(qh_vert_pos(&s->verts, i), nrml) - d0);
		if (d > maxDist && i != vtx[0] && i != vtx[1] && i != vtx[2]) { maxDist = d; vtx[3] = i; }
	}
	if (vtx[3] < 0 || maxDist <= 100*s->epsilon) return 0;

	// Build 4 triangle faces with correct winding.
	// With origin-vertex convention, edge indices within a triangle: e0, e1, e2.
	int tris[4];
	if (dot(qh_vert_pos(&s->verts, vtx[3]), nrml) - d0 < 0) {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[1], vtx[2]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[1], vtx[0]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[2], vtx[1]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[0], vtx[2]);
		for (int i = 0; i < 3; i++) {
			int k = (i+1)%3;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],0), qh_get_edge(s,tris[k+1],2));
			int our_k = (k==0)?2:(k==1)?0:1;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],1), qh_get_edge(s,tris[0],our_k));
		}
	} else {
		tris[0] = qh_create_triangle(s, vtx[0], vtx[2], vtx[1]);
		tris[1] = qh_create_triangle(s, vtx[3], vtx[0], vtx[1]);
		tris[2] = qh_create_triangle(s, vtx[3], vtx[1], vtx[2]);
		tris[3] = qh_create_triangle(s, vtx[3], vtx[2], vtx[0]);
		for (int i = 0; i < 3; i++) {
			int k = (i+1)%3;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],2), qh_get_edge(s,tris[k+1],0));
			int j = (3-i)%3, our_j = (j==0)?2:(j==1)?0:1;
			qh_set_opposite(s, qh_get_edge(s,tris[i+1],1), qh_get_edge(s,tris[0],our_j));
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
			if (i==vtx[0]||i==vtx[1]||i==vtx[2]||i==vtx[3]) continue;
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
	int v = vlist;
	while (v != QH_INVALID) {
		int nxt = s->verts.cnext[v];
		if (dot(ap.normal, qh_vert_pos(&s->verts, v)) - ap.offset > eps) {
			qh_conflict_add(s, absorb, v);
		} else {
			s->verts.cnext[v] = *unclaimed;
			*unclaimed = v;
		}
		v = nxt;
	}
}

// Iterative DFS from a visible face, marking visible faces DELETED.
// After a child pops, the parent's edge still points across the child face
// which is now DELETED, so the re-check naturally skips it before advancing.

typedef struct { int edge0; int edge; } QH_HFrame;

static void qh_calculate_horizon(QH_State* s, v3 eye, int edge0_init, int fi_init, CK_DYNA int** horizon, int* unclaimed)
{
	QH_HFrame hstk_buf[64];
	int hstk_n = 0;
	#define HSTK_PUSH(val) do { assert(hstk_n < 64); hstk_buf[hstk_n++] = (val); } while(0)
	#define HSTK_POP() (hstk_n--)
	#define HSTK_TOP() (&hstk_buf[hstk_n - 1])
	#define HSTK_EMPTY() (hstk_n == 0)

	qh_delete_face_points(s, fi_init, QH_INVALID, unclaimed);
	s->faces[fi_init].mark = QH_DELETED;
	int e0 = (edge0_init == QH_INVALID) ? s->faces[fi_init].edge : edge0_init;
	int estart = (edge0_init == QH_INVALID) ? e0 : s->edges[e0].next;
	HSTK_PUSH(((QH_HFrame){ e0, estart }));

	float eps = s->epsilon;
	while (!HSTK_EMPTY()) {
		QH_HFrame* f = HSTK_TOP();
		int ofi = qh_opp_face(s, f->edge);
		if (s->faces[ofi].mark == QH_VISIBLE) {
			if (dot(s->faces[ofi].plane.normal, eye) - s->faces[ofi].plane.offset > eps) {
				int child_e0 = s->edges[f->edge].twin;
				qh_delete_face_points(s, ofi, QH_INVALID, unclaimed);
				s->faces[ofi].mark = QH_DELETED;
				HSTK_PUSH(((QH_HFrame){ child_e0, s->edges[child_e0].next }));
				continue;
			} else {
				apush(*horizon, f->edge);
			}
		}
		f = HSTK_TOP();
		f->edge = s->edges[f->edge].next;
		if (f->edge == f->edge0) HSTK_POP();
	}
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
	qh_set_opposite(s, qh_get_edge(s, fi, 1), s->edges[horizon_edge].twin);
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
	int this_face = s->edges[en].face;

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
			int ht = s->edges[en].twin;
			int hp = s->edges[ep].twin;
			int e3 = s->edges[ht].next;
			if (e3 == hp) e3 = s->edges[ht].prev;
			new_twin = s->edges[e3].twin;
			s->faces[fin_face].mark = QH_DELETED;
			discarded = fin_face;
		} else {
			// Polygon fin should have been pre-merged by the caller.
			// If we still get here, use normal linkage as fallback.
			s->edges[ep].next = en;
			s->edges[en].prev = ep;
			return QH_INVALID;
		}

		// Remove ep from this face's loop; en absorbs ep's origin.
		s->edges[en].origin = s->edges[ep].origin;
		int pp = s->edges[ep].prev;
		QH_DEBUG(fprintf(stderr, "[qh]   splice out ep=%d: pp=%d->en=%d (en.next=%d)\n",
			ep, pp, en, s->edges[en].next));
		s->edges[en].prev = pp;
		s->edges[pp].next = en;

		qh_set_opposite(s, en, new_twin);
		qh_recompute_face(s, s->edges[new_twin].face);
		return discarded;
	} else {
		// Normal case: just link consecutively.
		s->edges[ep].next = en;
		s->edges[en].prev = ep;
		return QH_INVALID;
	}
}

// -----------------------------------------------------------------------------
// Merge the adjacent face (across a shared edge) into this face.

static int qh_merge_adjacent_face(QH_State* s, int hedge_adj, int discarded[], int* unclaimed)
{
	int ofi = qh_opp_face(s, hedge_adj);
	int tfi = s->edges[hedge_adj].face;
	int nd = 0;
	discarded[nd++] = ofi;
	QH_DEBUG(fprintf(stderr, "[qh] merge: face %d (mark=%d) absorbed into %d (mark=%d) via edge %d\n",
		ofi, s->faces[ofi].mark, tfi, s->faces[tfi].mark, hedge_adj));
	s->faces[ofi].mark = QH_DELETED;

	int ho = s->edges[hedge_adj].twin;
	int ap = s->edges[hedge_adj].prev, an = s->edges[hedge_adj].next;
	int op = s->edges[ho].prev, on = s->edges[ho].next;

	// Walk past multiply-shared edges.
	while (qh_opp_face(s, ap) == ofi) { ap = s->edges[ap].prev; on = s->edges[on].next; }
	while (qh_opp_face(s, an) == ofi) { op = s->edges[op].prev; an = s->edges[an].next; }

	// Reassign the absorbed face's edges to this face.
	for (int e = on; e != s->edges[op].next; e = s->edges[e].next) s->edges[e].face = tfi;

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
			s->edges[e].face = tfi;
			e = s->edges[e].next;
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
					int tw = s->edges[e].twin;
					int twf = s->edges[tw].face;
					if (s->faces[twf].mark == QH_DELETED || twf == tfi) {
						// Stale or self-edge: splice out.
						int p = s->edges[e].prev;
						int nx = s->edges[e].next;
						s->edges[p].next = nx;
						s->edges[nx].prev = p;
						s->edges[nx].origin = s->edges[e].origin;
						if (s->faces[tfi].edge == e)
							s->faces[tfi].edge = nx;
						cleaned = 1;
						break;
					}
					e = s->edges[e].next;
				} while (e != start);
			}
		}


		QH_DEBUG({
			e = s->faces[tfi].edge; start = e; n = 0;
			do { n++; e = s->edges[e].next; } while (e != start);
			fprintf(stderr, "[qh]   face %d loop (%d edges):", tfi, n);
			e = start;
			do { fprintf(stderr, " %d", e); e = s->edges[e].next; } while (e != start);
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
			hedge = s->edges[hedge].next; continue;
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
				fi, hedge, ofi, s->edges[hedge].face));
			// Check if the merge will hit a polygon fin degenerate.
			// If so, skip it to avoid infinite loops.
			int will_poly_fin = 0;
			{
				int ho = s->edges[hedge].twin;
				int ap = s->edges[hedge].prev, an = s->edges[hedge].next;
				int op = s->edges[ho].prev, on = s->edges[ho].next;
				while (qh_opp_face(s, ap) == ofi) { ap = s->edges[ap].prev; on = s->edges[on].next; }
				while (qh_opp_face(s, an) == ofi) { op = s->edges[op].prev; an = s->edges[an].next; }
				if ((qh_opp_face(s, op) == qh_opp_face(s, an) && qh_face_vert_count(s, qh_opp_face(s, op)) >= 4) ||
				    (qh_opp_face(s, ap) == qh_opp_face(s, on) && qh_face_vert_count(s, qh_opp_face(s, ap)) >= 4))
					will_poly_fin = 1;
			}
			if (will_poly_fin) {
				hedge = s->edges[hedge].next;
				continue;
			}
			int disc[3];
			int nd = qh_merge_adjacent_face(s, hedge, disc, unclaimed);
			for (int i = 0; i < nd; i++) qh_delete_face_points(s, disc[i], fi, unclaimed);
			return 1;
		}
		hedge = s->edges[hedge].next;
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

// Batch-allocate n contiguous edges at the end of the array.
static int qh_alloc_edges_contiguous(QH_State* s, int n)
{
	int base = asize(s->edges);
	afit(s->edges, base + n);
	asetlen(s->edges, base + n);
	memset(s->edges + base, 0, n * sizeof(QH_Edge));
	return base;
}

static void qh_add_new_faces(QH_State* s, QH_FaceList* nf, int eye, CK_DYNA int* horizon, int nh)
{
	nf->first = nf->last = QH_INVALID;
	int prev = QH_INVALID, begin = QH_INVALID;
	v3 eye_pos = qh_vert_pos(&s->verts, eye);
	// Batch-allocate all cone edges contiguously for cache-friendly writes.
	int edge_base = qh_alloc_edges_contiguous(s, nh * 3);
	for (int i = 0; i < nh; i++) {
		int horizon_edge = horizon[i];
		int tail = s->edges[horizon_edge].origin;
		int head = s->edges[s->edges[horizon_edge].next].origin;
		int fi = qh_alloc_face(s);
		int e0 = edge_base + i*3, e1 = e0+1, e2 = e0+2;
		s->edges[e0] = (QH_Edge){ .next=e1, .prev=e2, .twin=QH_INVALID, .origin=eye, .face=fi };
		s->edges[e1] = (QH_Edge){ .next=e2, .prev=e0, .twin=QH_INVALID, .origin=tail, .face=fi };
		s->edges[e2] = (QH_Edge){ .next=e0, .prev=e1, .twin=QH_INVALID, .origin=head, .face=fi };
		s->faces[fi].edge = e0;
		s->faces[fi].next = fi;
		s->faces[fi].prev = fi;
		v3 p0 = eye_pos, p1 = qh_vert_pos(&s->verts, tail), p2 = qh_vert_pos(&s->verts, head);
		v3 n = cross(sub(p1, p0), sub(p2, p0));
		float a = len(n);
		if (a > 0) n = scale(n, 1.0f / a);
		v3 c = scale(add(add(p0, p1), p2), 1.0f / 3.0f);
		s->faces[fi].plane = (HullPlane){ n, dot(n, c) };
		s->faces[fi].centroid = c;
		s->faces[fi].num_verts = 3;
		s->faces[fi].area = a;
		// Twin e1 with old twin of horizon edge.
		qh_set_opposite(s, e1, s->edges[horizon_edge].twin);
		// Link side edges between consecutive cone triangles.
		// e2.next == e0 (eye->tail), which is the inward side edge.
		if (prev != QH_INVALID) qh_set_opposite(s, e0, prev);
		else begin = e0;
		qh_face_list_add(s, nf, fi);
		prev = e2;
	}
	qh_set_opposite(s, begin, prev);
}

// Reassign orphaned conflict vertices from the unclaimed list to new faces.
static void qh_resolve_unclaimed(QH_State* s, QH_FaceList* nf, int* unclaimed)
{
	int v = *unclaimed;
	while (v != QH_INVALID) {
		int nxt = s->verts.cnext[v];
		float bd = s->epsilon; int bf = QH_INVALID;
		// First try new faces (most likely home for orphaned points).
		for (int fi = nf->first; fi != QH_INVALID; fi = s->faces[fi].next) {
			if (s->faces[fi].mark == QH_VISIBLE) {
				float d = qh_face_dist(s, fi, qh_vert_pos(&s->verts, v));
				if (d > bd) { bd = d; bf = fi; }
				if (bd > 1000*s->epsilon) break;
			}
		}
		// Fallback: if no new face claims the point, scan ALL live faces
		// with a relaxed threshold. A point may be numerically inside the
		// current hull but actually belong on it. Use the least-negative
		// distance (closest to being outside) as a last resort.
		if (bf == QH_INVALID) {
			float closest = -1e18f;
			for (int fi = 0; fi < asize(s->faces); fi++) {
				if (s->faces[fi].mark == QH_DELETED) continue;
				float d = qh_face_dist(s, fi, qh_vert_pos(&s->verts, v));
				if (d > closest) { closest = d; bf = fi; }
			}
			// Only assign if the point is reasonably close to some face plane.
			// If it's deeply inside (> 100x epsilon below all planes), it's truly interior.
			if (closest < -100.0f * s->epsilon) bf = QH_INVALID;
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
			int tw = s->edges[e].twin;
			if (s->edges[tw].twin != e) {
				fprintf(stderr, "qh_validate: edge %d twin=%d but twin's twin=%d at %s\n",
					e, tw, s->edges[tw].twin, ctx);
				exit(-1);
			}
			// Check twin's face is alive.
			int twf = s->edges[tw].face;
			if (s->faces[twf].mark != QH_VISIBLE) {
				fprintf(stderr, "qh_validate: edge %d (face %d) twin=%d on deleted face %d at %s\n",
					e, fi, tw, twf, ctx);
				exit(-1);
			}
			// Check edge belongs to this face.
			if (s->edges[e].face != fi) {
				fprintf(stderr, "qh_validate: edge %d face=%d expected %d at %s\n",
					e, s->edges[e].face, fi, ctx);
				exit(-1);
			}
			n++;
			e = s->edges[e].next;
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
				vremap[s->edges[e].origin] = 1;
				e = s->edges[e].next;
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
	CK_DYNA int* horizon = NULL;
	int unclaimed = QH_INVALID;

	qh_conflict_remove(s, eye_face, eye);
	qh_calculate_horizon(s, qh_vert_pos(&s->verts, eye), QH_INVALID, eye_face, &horizon, &unclaimed);

	QH_FaceList nf;
	qh_add_new_faces(s, &nf, eye, horizon, asize(horizon));

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
						hedge = s->edges[hedge].next;
						continue;
					}
					if (qh_face_seen[ofi] == qh_face_gen) {
						// Multi-adjacent: ofi seen twice.
						int small = (s->faces[fi].num_verts <= s->faces[ofi].num_verts) ? fi : ofi;
						int large = (small == fi) ? ofi : fi;
						int me = s->faces[large].edge;
						while (qh_opp_face(s, me) != small) me = s->edges[me].next;
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
					hedge = s->edges[hedge].next;
				} while (hedge != s->faces[fi].edge);
			}
		}
	}

	qh_resolve_unclaimed(s, &nf, &unclaimed);
	afree(horizon);
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
		int total_edges = asize(s->edges);
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
				e = s->edges[e].next;
			} while (e != start);
		}
	});

	// Vertex + edge remap in a single pass over live face edges.
	CK_DYNA int* vremap = NULL;
	afit_set(vremap, s->verts.count);
	memset(vremap, -1, s->verts.count * sizeof(int));
	CK_DYNA int* eremap = NULL;
	afit_set(eremap, asize(s->edges));
	memset(eremap, -1, asize(s->edges) * sizeof(int));
	CK_DYNA v3* ov = NULL; afit(ov, s->verts.count); int vc = 0;
	CK_DYNA HalfEdge* oe = NULL; afit(oe, asize(s->edges)); int ec = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do {
			int vi = s->edges[e].origin;
			if (vremap[vi] < 0) { vremap[vi] = vc++; apush(ov, qh_vert_pos(&s->verts, vi)); }
			eremap[e] = ec++; HalfEdge he = {0}; apush(oe, he);
			e = s->edges[e].next;
		} while (e != start);
	}
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { int o=eremap[e];
			oe[o].next=(uint16_t)eremap[s->edges[e].next];
			oe[o].twin=(uint16_t)eremap[s->edges[e].twin];
			oe[o].origin=(uint16_t)vremap[s->edges[e].origin];
			oe[o].face=(uint16_t)i;
			e=s->edges[e].next;
		} while (e!=start);
	}

	CK_DYNA HullFace* of = NULL;
	CK_DYNA HullPlane* op = NULL;
	for (int i = 0; i < asize(live); i++) {
		HullFace hf = { .edge=(uint16_t)eremap[s->faces[live[i]].edge] };
		apush(of, hf); apush(op, s->faces[live[i]].plane);
	}

	v3 centroid = V3(0,0,0);
	for (int i = 0; i < vc; i++) centroid = add(centroid, ov[i]);
	centroid = scale(centroid, 1.0f / vc);

	Hull* h = CK_ALLOC(sizeof(Hull));
	h->centroid = centroid; h->vert_count = vc; h->edge_count = ec; h->face_count = asize(live);
	h->epsilon = s->epsilon;

	v3* vcp = CK_ALLOC(sizeof(v3)*vc);       memcpy(vcp, ov, sizeof(v3)*vc);       h->verts = vcp;
	HalfEdge* ecp = CK_ALLOC(sizeof(HalfEdge)*ec); memcpy(ecp, oe, sizeof(HalfEdge)*ec); h->edges = ecp;
	HullFace* fcp = CK_ALLOC(sizeof(HullFace)*asize(live)); memcpy(fcp, of, sizeof(HullFace)*asize(live)); h->faces = fcp;
	HullPlane* pcp = CK_ALLOC(sizeof(HullPlane)*asize(live)); memcpy(pcp, op, sizeof(HullPlane)*asize(live)); h->planes = pcp;

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
			float d = nx[i]*all_points[pi].x + ny[i]*all_points[pi].y + nz[i]*all_points[pi].z;
			if (d > max_ds[i]) max_ds[i] = d;
		}
	}
	for (int i = 0; i < fc; i++) {
		float widen = max_ds[i] - pcp[i].offset;
		if (widen > h->maxoutside) h->maxoutside = widen;
		pcp[i].offset = max_ds[i];
	}
	if (soa_bytes > 4096) free(nx);

	afree(live); afree(vremap); afree(ov); afree(eremap); afree(oe); afree(of); afree(op);
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
	state.edge_free = QH_INVALID;
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
		qh_verts_free(&state.verts); afree(state.edges); afree(state.faces); afree(state.conflict_faces);
		return NULL;
	}

	if (!qh_build_simplex(&state, welded_count)) {
		qh_verts_free(&state.verts); afree(state.edges); afree(state.faces); afree(state.conflict_faces);
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
	qh_verts_free(&state.verts); afree(state.edges); afree(state.faces); afree(state.conflict_faces);
	return result;
}
