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

typedef struct QH_Vertex {
	v3 pos;
	int conflict_next;
	int conflict_prev;
} QH_Vertex;

typedef struct QH_Edge {
	int next, prev, twin;
	int origin, face, mark;
} QH_Edge;

typedef struct QH_Face {
	int edge;
	int next, prev;
	int conflict_head;
	int mark;
	int num_verts;
	float area;
	float maxoutside; // max distance any point was ever seen outside this face
	HullPlane plane;
} QH_Face;

typedef struct QH_State {
	CK_DYNA QH_Vertex* verts;
	CK_DYNA QH_Edge*   edges;
	CK_DYNA QH_Face*   faces;
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
		s->edges[idx] = (QH_Edge){0};
		return idx;
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
		s->faces[idx] = (QH_Face){0};
		s->faces[idx].conflict_head = QH_INVALID;
		s->faces[idx].mark = QH_VISIBLE;
		s->faces[idx].maxoutside = s->epsilon;
		return idx;
	}
	QH_Face f = {0};
	f.conflict_head = QH_INVALID;
	f.mark = QH_VISIBLE;
	f.maxoutside = s->epsilon;
	apush(s->faces, f);
	return asize(s->faces) - 1;
}

// Forward declaration (used by qh_conflict_add before full definition).
static float qh_plane_dist(HullPlane p, v3 pt);

// -----------------------------------------------------------------------------
// Conflict list (circular doubly-linked via vertex next/prev).

static void qh_conflict_add(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Vertex* v = &s->verts[vi];
	if (f->conflict_head != QH_INVALID) {
		QH_Vertex* head = &s->verts[f->conflict_head];
		int tail = head->conflict_prev;
		v->conflict_next = f->conflict_head;
		v->conflict_prev = tail;
		head->conflict_prev = vi;
		s->verts[tail].conflict_next = vi;
	} else {
		v->conflict_next = vi;
		v->conflict_prev = vi;
	}
	f->conflict_head = vi;
	// Track max distance any point was seen outside this face.
	float d = qh_plane_dist(f->plane, v->pos);
	if (d > f->maxoutside) f->maxoutside = d;
}

static void qh_conflict_remove(QH_State* s, int fi, int vi)
{
	QH_Face* f = &s->faces[fi];
	QH_Vertex* v = &s->verts[vi];
	if (v->conflict_next == vi) {
		f->conflict_head = QH_INVALID;
	} else {
		s->verts[v->conflict_next].conflict_prev = v->conflict_prev;
		s->verts[v->conflict_prev].conflict_next = v->conflict_next;
		if (f->conflict_head == vi) f->conflict_head = v->conflict_next;
	}
	v->conflict_next = QH_INVALID;
	v->conflict_prev = QH_INVALID;
}

// Remove all conflict vertices. Returns singly-linked list (via conflict_next),
// terminated by QH_INVALID.
static int qh_conflict_remove_all(QH_State* s, int fi)
{
	QH_Face* f = &s->faces[fi];
	int head = f->conflict_head;
	if (head == QH_INVALID) return QH_INVALID;
	int last = s->verts[head].conflict_prev;
	s->verts[last].conflict_next = QH_INVALID;
	f->conflict_head = QH_INVALID;
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
		v3 cur = s->verts[s->edges[e].origin].pos;
		v3 nxt = s->verts[s->edges[s->edges[e].next].origin].pos;
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
	f->num_verts = count;
	f->area = a;
}

static float qh_face_dist(QH_State* s, int fi, v3 pt)
{
	return qh_plane_dist(s->faces[fi].plane, pt);
}

static v3 qh_face_centroid(QH_State* s, int fi)
{
	int e = s->faces[fi].edge, start = e;
	v3 c = V3(0,0,0); int n = 0;
	do { c = add(c, s->verts[s->edges[e].origin].pos); n++; e = s->edges[e].next; } while (e != start);
	return scale(c, 1.0f / n);
}

// Distance of the adjacent face's centroid to this edge's face plane.
static float qh_opp_face_dist(QH_State* s, int ei)
{
	int fi = s->edges[ei].face;
	int ofi = s->edges[s->edges[ei].twin].face;
	return qh_face_dist(s, fi, qh_face_centroid(s, ofi));
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
	qh_recompute_face(s, fi);
	return fi;
}

// -----------------------------------------------------------------------------
// Build initial tetrahedron from 4 non-coplanar extremal points.

static int qh_build_simplex(QH_State* s, int nv)
{
	// Find extremal vertices on each axis.
	int maxV[3], minV[3];
	for (int i = 0; i < 3; i++) { maxV[i] = minV[i] = 0; }
	v3 mx = s->verts[0].pos, mn = s->verts[0].pos;
	for (int i = 1; i < nv; i++) {
		v3 p = s->verts[i].pos;
		if (p.x > mx.x) { mx.x = p.x; maxV[0] = i; } else if (p.x < mn.x) { mn.x = p.x; minV[0] = i; }
		if (p.y > mx.y) { mx.y = p.y; maxV[1] = i; } else if (p.y < mn.y) { mn.y = p.y; minV[1] = i; }
		if (p.z > mx.z) { mx.z = p.z; maxV[2] = i; } else if (p.z < mn.z) { mn.z = p.z; minV[2] = i; }
	}

	// Pick axis with greatest spread.
	float best = 0; int imax = 0;
	for (int i = 0; i < 3; i++) {
		v3 a = s->verts[maxV[i]].pos, b = s->verts[minV[i]].pos;
		float d = (i==0) ? a.x-b.x : (i==1) ? a.y-b.y : a.z-b.z;
		if (d > best) { best = d; imax = i; }
	}
	if (best <= s->epsilon) return 0;

	int vtx[4];
	vtx[0] = maxV[imax]; vtx[1] = minV[imax];

	// Find vertex farthest from line vtx[0]-vtx[1].
	v3 u01 = norm(sub(s->verts[vtx[1]].pos, s->verts[vtx[0]].pos));
	float maxSqr = 0; v3 nrml = V3(0,0,0); vtx[2] = -1;
	for (int i = 0; i < nv; i++) {
		v3 xp = cross(u01, sub(s->verts[i].pos, s->verts[vtx[0]].pos));
		float ls = len2(xp);
		if (ls > maxSqr && i != vtx[0] && i != vtx[1]) { maxSqr = ls; vtx[2] = i; nrml = xp; }
	}
	if (vtx[2] < 0 || sqrtf(maxSqr) <= 100*s->epsilon) return 0;
	nrml = norm(nrml);
	nrml = norm(sub(nrml, scale(u01, dot(nrml, u01)))); // orthogonalize

	// Find vertex farthest from plane through vtx[0..2].
	float d0 = dot(s->verts[vtx[2]].pos, nrml);
	float maxDist = 0; vtx[3] = -1;
	for (int i = 0; i < nv; i++) {
		float d = fabsf(dot(s->verts[i].pos, nrml) - d0);
		if (d > maxDist && i != vtx[0] && i != vtx[1] && i != vtx[2]) { maxDist = d; vtx[3] = i; }
	}
	if (vtx[3] < 0 || maxDist <= 100*s->epsilon) return 0;

	// Build 4 triangle faces with correct winding.
	// With origin-vertex convention, edge indices within a triangle: e0, e1, e2.
	int tris[4];
	if (dot(s->verts[vtx[3]].pos, nrml) - d0 < 0) {
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

	s->interior = scale(add(add(s->verts[vtx[0]].pos, s->verts[vtx[1]].pos),
	                        add(s->verts[vtx[2]].pos, s->verts[vtx[3]].pos)), 0.25f);

	// Assign conflict vertices: each point goes to the face it's furthest outside of.
	for (int i = 0; i < nv; i++) {
		if (i==vtx[0]||i==vtx[1]||i==vtx[2]||i==vtx[3]) {
			continue;
		}
		float bd = s->epsilon; int bf = -1;
		for (int k = 0; k < 4; k++) {
			float d = qh_face_dist(s, tris[k], s->verts[i].pos);
			if (d > bd) { bd = d; bf = tris[k]; }
		}
		if (bf >= 0) qh_conflict_add(s, bf, i);
	}
	return 1;
}

// -----------------------------------------------------------------------------
// Find next conflict vertex (globally furthest from any face).

static int qh_next_conflict(QH_State* s, int* out_face)
{
	int bv = QH_INVALID, bf = QH_INVALID; float bd = -1e18f;
	for (int fi = 0; fi < asize(s->faces); fi++) {
		if (s->faces[fi].mark != QH_VISIBLE) continue;
		int head = s->faces[fi].conflict_head;
		if (head == QH_INVALID) continue;
		int ci = head;
		do {
			float d = qh_face_dist(s, fi, s->verts[ci].pos);
			if (d > bd) { bd = d; bv = ci; bf = fi; }
			ci = s->verts[ci].conflict_next;
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
	int v = vlist;
	while (v != QH_INVALID) {
		int nxt = s->verts[v].conflict_next;
		if (absorb != QH_INVALID && qh_face_dist(s, absorb, s->verts[v].pos) > s->epsilon) {
			qh_conflict_add(s, absorb, v);
		} else {
			s->verts[v].conflict_next = *unclaimed;
			*unclaimed = v;
		}
		v = nxt;
	}
}

// -----------------------------------------------------------------------------
// Recursive DFS from a visible face, marking visible faces DELETED.

static void qh_calculate_horizon(QH_State* s, v3 eye, int edge0, int fi,
                                  CK_DYNA int** horizon, int* unclaimed)
{
	qh_delete_face_points(s, fi, QH_INVALID, unclaimed);
	s->faces[fi].mark = QH_DELETED;

	int edge;
	if (edge0 == QH_INVALID) { edge0 = s->faces[fi].edge; edge = edge0; }
	else { edge = s->edges[edge0].next; }

	do {
		int ofi = qh_opp_face(s, edge);
		if (s->faces[ofi].mark == QH_VISIBLE) {
			if (qh_face_dist(s, ofi, eye) > s->epsilon)
				qh_calculate_horizon(s, eye, s->edges[edge].twin, ofi, horizon, unclaimed);
			else
				apush(*horizon, edge);
		}
		edge = s->edges[edge].next;
	} while (edge != edge0);
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
			// Polygon: remove one edge from the fin face's loop.
			int ht = s->edges[en].twin;
			new_twin = s->edges[ht].next;
			if (s->faces[fin_face].edge == s->edges[new_twin].prev)
				s->faces[fin_face].edge = new_twin;
			// new_twin absorbs the removed edge's origin.
			int removed = s->edges[new_twin].prev;
			s->edges[new_twin].origin = s->edges[removed].origin;
			int skip = s->edges[removed].prev;
			s->edges[new_twin].prev = skip;
			s->edges[skip].next = new_twin;
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
	do {
		int ofi = qh_opp_face(s, hedge);
		if (ofi == fi || s->faces[ofi].mark == QH_DELETED) { hedge = s->edges[hedge].next; continue; }

		float tol = s->epsilon;

		int merge = 0;
		if (type == QH_MERGE_ANY) {
			if (qh_opp_face_dist(s, hedge) > -tol ||
			    qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) merge = 1;
		} else {
			if (s->faces[fi].area > s->faces[ofi].area) {
				if (qh_opp_face_dist(s, hedge) > -tol) merge = 1;
				else if (qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) convex = 0;
			} else {
				if (qh_opp_face_dist(s, s->edges[hedge].twin) > -tol) merge = 1;
				else if (qh_opp_face_dist(s, hedge) > -tol) convex = 0;
			}
		}
		if (merge) {
			QH_DEBUG(fprintf(stderr, "[qh] do_adjacent_merge: fi=%d hedge=%d ofi=%d hedge.face=%d\n",
				fi, hedge, ofi, s->edges[hedge].face));
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

static void qh_add_new_faces(QH_State* s, QH_FaceList* nf, int eye,
                              CK_DYNA int* horizon, int nh)
{
	nf->first = nf->last = QH_INVALID;
	int prev = QH_INVALID, begin = QH_INVALID;
	for (int i = 0; i < nh; i++) {
		int side = qh_add_adjoining_face(s, eye, horizon[i]);
		if (prev != QH_INVALID) qh_set_opposite(s, s->edges[side].next, prev);
		else begin = side;
		qh_face_list_add(s, nf, s->edges[side].face);
		prev = side;
	}
	qh_set_opposite(s, s->edges[begin].next, prev);
}

// Reassign orphaned conflict vertices from the unclaimed list to new faces.
static void qh_resolve_unclaimed(QH_State* s, QH_FaceList* nf, int* unclaimed)
{
	int v = *unclaimed;
	while (v != QH_INVALID) {
		int nxt = s->verts[v].conflict_next;
		float bd = s->epsilon; int bf = QH_INVALID;
		// First try new faces (most likely home for orphaned points).
		for (int fi = nf->first; fi != QH_INVALID; fi = s->faces[fi].next) {
			if (s->faces[fi].mark == QH_VISIBLE) {
				float d = qh_face_dist(s, fi, s->verts[v].pos);
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
				float d = qh_face_dist(s, fi, s->verts[v].pos);
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
}

// -----------------------------------------------------------------------------
// Process one conflict vertex: compute horizon, build cone, merge, reassign.

static void qh_add_point(QH_State* s, int eye, int eye_face)
{
	CK_DYNA int* horizon = NULL;
	int unclaimed = QH_INVALID;

	qh_conflict_remove(s, eye_face, eye);
	qh_calculate_horizon(s, s->verts[eye].pos, QH_INVALID, eye_face, &horizon, &unclaimed);

	QH_FaceList nf;
	qh_add_new_faces(s, &nf, eye, horizon, asize(horizon));

	// Two-pass merge: first larger-face-biased, then mutual non-convexity.
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_VISIBLE)
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_LARGE, &unclaimed));
	for (int fi = nf.first; fi != QH_INVALID; fi = s->faces[fi].next)
		if (s->faces[fi].mark == QH_NON_CONVEX) {
			s->faces[fi].mark = QH_VISIBLE;
			while (qh_do_adjacent_merge(s, fi, QH_MERGE_ANY, &unclaimed));
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

	// Vertex remap.
	CK_DYNA int* vremap = NULL;
	afit(vremap, asize(s->verts));
	for (int i = 0; i < asize(s->verts); i++) apush(vremap, -1);
	CK_DYNA v3* ov = NULL; int vc = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { int vi = s->edges[e].origin;
			if (vremap[vi]<0) { vremap[vi]=vc++; apush(ov, s->verts[vi].pos); }
			e = s->edges[e].next;
		} while (e != start);
	}

	// Edge remap.
	CK_DYNA int* eremap = NULL;
	afit(eremap, asize(s->edges));
	for (int i = 0; i < asize(s->edges); i++) apush(eremap, -1);
	CK_DYNA HalfEdge* oe = NULL; int ec = 0;
	for (int i = 0; i < asize(live); i++) {
		int e = s->faces[live[i]].edge, start = e;
		do { eremap[e]=ec++; HalfEdge he={0}; apush(oe, he); e=s->edges[e].next; } while (e!=start);
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
		v3 fc = V3(0,0,0); int cnt = 0;
		int start = fcp[i].edge, e = start;
		do { fc = add(fc, vcp[ecp[e].origin]); cnt++; e = ecp[e].next; } while (e != start);
		fc = scale(fc, 1.0f / cnt);
		if (dot(pcp[i].normal, sub(fc, centroid)) < 0) {
			pcp[i].normal = scale(pcp[i].normal, -1.0f);
			pcp[i].offset = -pcp[i].offset;
		}
	}

	// Post-build plane widening: widen each plane so ALL original input points
	// lie on or behind the plane. Must use original points (not welded subset)
	// since welded-away points may lie outside the welded hull.
	h->maxoutside = 0;
	for (int i = 0; i < h->face_count; i++) {
		float max_d = pcp[i].offset;
		for (int pi = 0; pi < all_count; pi++) {
			float d = dot(pcp[i].normal, all_points[pi]);
			if (d > max_d) max_d = d;
		}
		float widen = max_d - pcp[i].offset;
		if (widen > h->maxoutside) h->maxoutside = widen;
		pcp[i].offset = max_d;
	}

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
	for (int i = 0; i < count; i++) {
		v3 p = points[i];
		int dup = 0;
		for (int j = 0; j < asize(state.verts); j++) {
			if (len2(sub(p, state.verts[j].pos)) <= weld_dist2) { dup = 1; break; }
		}
		if (!dup) {
			QH_Vertex v = { .pos=p, .conflict_next=QH_INVALID, .conflict_prev=QH_INVALID };
			apush(state.verts, v);
		}
	}
	int welded_count = asize(state.verts);

	if (welded_count < 4) {
		afree(state.verts); afree(state.edges); afree(state.faces);
		return NULL;
	}

	if (!qh_build_simplex(&state, welded_count)) {
		afree(state.verts); afree(state.edges); afree(state.faces);
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
	afree(state.verts); afree(state.edges); afree(state.faces);
	return result;
}
