// viewer.c -- Remote debug viewer for nudge physics engine.
// Reads engine memory via ReadProcessMemory. Type layout received from engine on connect.
// Zero hardcoded struct layouts -- all driven by reflected type tables from the engine.
//
// Build: cmake target nudge_viewer
// Usage: nudge_viewer [host] [port]

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment(lib, "ws2_32.lib")
#include <windows.h>
#else
#error "ReadProcessMemory viewer is Windows-only"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

// ============================================================================
// ReadProcessMemory helpers
// ============================================================================

static HANDLE g_proc;
static DWORD g_pid;
static uint64_t g_world_ptr;

static int rpm_read(uint64_t addr, void *buf, size_t size)
{
	SIZE_T n;
	return ReadProcessMemory(g_proc, (LPCVOID)addr, buf, size, &n) && n == size;
}

static int rpm_read_ptr(uint64_t addr, uint64_t *out)
{
	*out = 0;
	return rpm_read(addr, out, 8);
}

static int rpm_ptr_valid(uint64_t addr)
{
	if (addr == 0 || addr < 0x10000 || addr > 0x00007FFFFFFFFFFF) return 0;
	MEMORY_BASIC_INFORMATION mbi;
	if (!VirtualQueryEx(g_proc, (LPCVOID)addr, &mbi, sizeof(mbi))) return 0;
	return mbi.State == MEM_COMMIT && !(mbi.Protect & (PAGE_NOACCESS | PAGE_GUARD));
}

// CKit array header (must match ckit.h CK_ArrayHeader).
#define CK_ARRAY_COOKIE 0x59525241
typedef struct { int size, capacity, is_static, alignment; uint64_t alloc_base; uint32_t cookie; } CKArrayHdr;

static int rpm_array_count(uint64_t data_ptr)
{
	if (!data_ptr || !rpm_ptr_valid(data_ptr)) return 0;
	CKArrayHdr hdr;
	if (!rpm_read(data_ptr - sizeof(CKArrayHdr), &hdr, sizeof(hdr))) return 0;
	if (hdr.cookie != CK_ARRAY_COOKIE || hdr.size < 0 || hdr.size > 10000000) return -1;
	return hdr.size;
}

// CKit map header.
typedef struct { uint32_t cookie; int val_size, size, capacity, slot_count, slot_capacity, _pad[2]; } CKMapHdr;

static int rpm_map_count(uint64_t items_ptr)
{
	if (!items_ptr || !rpm_ptr_valid(items_ptr)) return 0;
	CKMapHdr hdr;
	if (!rpm_read(items_ptr - sizeof(CKMapHdr), &hdr, sizeof(hdr))) return 0;
	return (hdr.size >= 0 && hdr.size < 10000000) ? hdr.size : -1;
}

// ============================================================================
// Dynamic type system (parsed from engine on connect)
// ============================================================================

typedef enum FieldKind {
	FK_INT, FK_UINT, FK_FLOAT, FK_DOUBLE, FK_BOOL,
	FK_U8, FK_U16, FK_U64,
	FK_V3, FK_QUAT,
	FK_ENUM, FK_PTR, FK_ARRAY, FK_MAP, FK_STRUCT, FK_PTR_RAW,
	FK_COUNT
} FieldKind;

static const char *fk_names[] = { "INT","UINT","FLOAT","DOUBLE","BOOL","U8","U16","U64","V3","QUAT","ENUM","PTR","ARRAY","MAP","STRUCT","PTR_RAW" };

typedef struct EnumValue { char name[64]; int value; } EnumValue;
typedef struct EnumTable { char name[64]; EnumValue *values; int count; } EnumTable;

typedef struct FieldDesc {
	char name[64];
	FieldKind kind;
	int offset, size;
	char nested_name[64];  // type name for PTR/ARRAY/MAP/STRUCT
	char enum_name[64];    // enum table name for ENUM
} FieldDesc;

typedef struct TypeDesc {
	char name[64];
	int size;
	FieldDesc *fields;
	int field_count;
} TypeDesc;

static TypeDesc *g_types;
static int g_type_count;
static EnumTable *g_enums;
static int g_enum_count;

static TypeDesc *find_type(const char *name)
{
	for (int i = 0; i < g_type_count; i++)
		if (strcmp(g_types[i].name, name) == 0) return &g_types[i];
	return NULL;
}

static EnumTable *find_enum(const char *name)
{
	for (int i = 0; i < g_enum_count; i++)
		if (strcmp(g_enums[i].name, name) == 0) return &g_enums[i];
	return NULL;
}

static const char *enum_value_name(EnumTable *et, int val)
{
	if (!et) return "?";
	for (int i = 0; i < et->count; i++)
		if (et->values[i].value == val) return et->values[i].name;
	return "?";
}

static FieldDesc *find_field(TypeDesc *t, const char *name)
{
	for (int i = 0; i < t->field_count; i++)
		if (strcmp(t->fields[i].name, name) == 0) return &t->fields[i];
	return NULL;
}

static FieldKind parse_field_kind(const char *s)
{
	for (int i = 0; i < FK_COUNT; i++)
		if (strcmp(fk_names[i], s) == 0) return (FieldKind)i;
	return FK_INT;
}

// ============================================================================
// Parse type data from engine (received on connect)
// ============================================================================

static int parse_types(const char *data)
{
	const char *p = data;
	if (sscanf(p, "TYPES %d %d", &g_type_count, &g_enum_count) != 2) return 0;
	p = strchr(p, '\n'); if (!p) return 0; p++;

	g_types = calloc(g_type_count, sizeof(TypeDesc));
	for (int t = 0; t < g_type_count; t++) {
		TypeDesc *td = &g_types[t];
		int sz;
		sscanf(p, "TYPE %63s %d", td->name, &sz);
		td->size = sz;
		p = strchr(p, '\n'); if (!p) return 0; p++;

		// Count fields
		const char *scan = p;
		int fc = 0;
		while (strncmp(scan, "END", 3) != 0) { fc++; scan = strchr(scan, '\n') + 1; }
		td->field_count = fc;
		td->fields = calloc(fc, sizeof(FieldDesc));

		for (int f = 0; f < fc; f++) {
			FieldDesc *fd = &td->fields[f];
			char kind_str[32] = {0};
			char extra[64] = {0};
			int n = sscanf(p, "FIELD %63s %31s %d %d %63s", fd->name, kind_str, &fd->offset, &fd->size, extra);
			fd->kind = parse_field_kind(kind_str);
			if (n >= 5) {
				if (fd->kind == FK_ENUM) strncpy(fd->enum_name, extra, 63);
				else strncpy(fd->nested_name, extra, 63);
			}
			p = strchr(p, '\n'); if (!p) return 0; p++;
		}
		p = strchr(p, '\n'); if (!p) return 0; p++; // END
	}

	g_enums = calloc(g_enum_count, sizeof(EnumTable));
	for (int e = 0; e < g_enum_count; e++) {
		EnumTable *et = &g_enums[e];
		int cnt;
		sscanf(p, "ENUM %63s %d", et->name, &cnt);
		et->count = cnt;
		et->values = calloc(cnt, sizeof(EnumValue));
		p = strchr(p, '\n'); if (!p) return 0; p++;
		for (int v = 0; v < cnt; v++) {
			sscanf(p, "VALUE %63s %d", et->values[v].name, &et->values[v].value);
			p = strchr(p, '\n'); if (!p) return 0; p++;
		}
		p = strchr(p, '\n'); if (!p) return 0; p++; // END
	}
	return 1;
}

// ============================================================================
// Generic text printer (reads remote memory, formats from TypeDesc)
// ============================================================================

static void print_field_value(FieldDesc *f, const void *base, int indent)
{
	const char *p = (const char *)base + f->offset;
	for (int i = 0; i < indent; i++) printf("  ");

	switch (f->kind) {
	case FK_INT:    printf("%s: %d\n", f->name, *(int32_t *)p); break;
	case FK_UINT:   printf("%s: %u\n", f->name, *(uint32_t *)p); break;
	case FK_FLOAT:  printf("%s: %.4f\n", f->name, *(float *)p); break;
	case FK_DOUBLE: printf("%s: %.6f\n", f->name, *(double *)p); break;
	case FK_BOOL:   printf("%s: %s\n", f->name, (f->size == 1 ? *(uint8_t *)p : *(int32_t *)p) ? "true" : "false"); break;
	case FK_U64:    printf("%s: 0x%llX\n", f->name, *(uint64_t *)p); break;
	case FK_V3:     { float *v = (float *)p; printf("%s: (%.3f, %.3f, %.3f)\n", f->name, v[0], v[1], v[2]); break; }
	case FK_QUAT:   { float *v = (float *)p; printf("%s: (%.3f, %.3f, %.3f, %.3f)\n", f->name, v[0], v[1], v[2], v[3]); break; }
	case FK_ENUM:   { int val = *(int32_t *)p; EnumTable *et = find_enum(f->enum_name); printf("%s: %d (%s)\n", f->name, val, enum_value_name(et, val)); break; }
	case FK_PTR: case FK_PTR_RAW: printf("%s: 0x%llX\n", f->name, *(uint64_t *)p); break;
	case FK_ARRAY:  { uint64_t ap = *(uint64_t *)p; printf("%s: [%d] @0x%llX\n", f->name, ap ? rpm_array_count(ap) : 0, ap); break; }
	case FK_MAP:    { uint64_t mp = *(uint64_t *)p; printf("%s: {%d} @0x%llX\n", f->name, mp ? rpm_map_count(mp) : 0, mp); break; }
	case FK_STRUCT: {
		TypeDesc *nt = find_type(f->nested_name);
		if (nt) {
			printf("%s:\n", f->name);
			for (int i = 0; i < nt->field_count; i++)
				print_field_value(&nt->fields[i], p, indent + 1);
		} else printf("%s: (unknown type %s)\n", f->name, f->nested_name);
		break;
	}
	default: printf("%s: ?\n", f->name); break;
	}
}

static void print_struct(TypeDesc *t, const void *data, int indent)
{
	for (int i = 0; i < t->field_count; i++)
		print_field_value(&t->fields[i], data, indent);
}

// ============================================================================
// Path resolution
// ============================================================================

typedef struct Resolved { uint64_t addr; TypeDesc *type; int is_array; int count; int ok; } Resolved;

static Resolved resolve(uint64_t base, TypeDesc *type, const char *path)
{
	Resolved r = { base, type, 0, 0, 1 };
	char seg[128];
	const char *p = path;
	while (*p == '/') p++;
	if (!*p) return r;

	while (*p && r.ok) {
		int si = 0;
		while (*p && *p != '/' && si < 127) seg[si++] = *p++;
		seg[si] = '\0';
		while (*p == '/') p++;

		// Numeric index into array?
		char *endp;
		long idx = strtol(seg, &endp, 10);
		if (*endp == '\0' && r.is_array) {
			if (idx < 0 || idx >= r.count) { printf("ERR index %ld out of range [0, %d)\n", idx, r.count); r.ok = 0; return r; }
			r.addr += (uint64_t)idx * r.type->size;
			r.is_array = 0;
			continue;
		}

		if (r.is_array) { printf("ERR expected index, got '%s'\n", seg); r.ok = 0; return r; }

		// Read struct to find field value
		void *local = malloc(r.type->size);
		if (!rpm_read(r.addr, local, r.type->size)) { printf("ERR read failed at 0x%llX\n", r.addr); free(local); r.ok = 0; return r; }

		FieldDesc *f = find_field(r.type, seg);
		if (!f) {
			printf("ERR no field '%s' in %s (fields:", seg, r.type->name);
			for (int i = 0; i < r.type->field_count; i++) printf(" %s", r.type->fields[i].name);
			printf(")\n");
			free(local); r.ok = 0; return r;
		}

		uint64_t val = 0;
		memcpy(&val, (char *)local + f->offset, f->size > 8 ? 8 : f->size);
		free(local);

		switch (f->kind) {
		case FK_STRUCT:
			r.addr += f->offset;
			r.type = find_type(f->nested_name);
			if (!r.type) { printf("ERR unknown type %s\n", f->nested_name); r.ok = 0; }
			break;
		case FK_PTR:
			if (!val || !rpm_ptr_valid(val)) { printf("ERR null/invalid ptr %s\n", f->name); r.ok = 0; return r; }
			r.addr = val;
			r.type = find_type(f->nested_name);
			if (!r.type) { printf("ERR unknown type %s\n", f->nested_name); r.ok = 0; }
			break;
		case FK_ARRAY: case FK_MAP: {
			if (!val) { printf("ERR null %s\n", f->name); r.ok = 0; return r; }
			int cnt = (f->kind == FK_ARRAY) ? rpm_array_count(val) : rpm_map_count(val);
			if (cnt < 0) { printf("ERR corrupt header for %s\n", f->name); r.ok = 0; return r; }
			r.addr = val; r.is_array = 1; r.count = cnt;
			r.type = find_type(f->nested_name);
			if (!r.type) { printf("ERR unknown type %s\n", f->nested_name); r.ok = 0; }
			break;
		}
		default:
			if (!*p) { r.addr += f->offset; return r; } // leaf
			printf("ERR can't drill into scalar '%s'\n", f->name); r.ok = 0; return r;
		}
	}
	return r;
}

// ============================================================================
// Snapshot for diff
// ============================================================================

typedef struct Snapshot { uint64_t addr; TypeDesc *type; void *data; int count; int valid; } Snapshot;
static Snapshot g_snap;

static void snap_save(Resolved r)
{
	free(g_snap.data);
	g_snap = (Snapshot){0};
	if (!r.ok) return;
	int count = r.is_array ? r.count : 1;
	size_t bytes = (size_t)count * r.type->size;
	g_snap.data = malloc(bytes);
	if (!rpm_read(r.addr, g_snap.data, bytes)) { free(g_snap.data); g_snap.data = NULL; return; }
	g_snap.addr = r.addr;
	g_snap.type = r.type;
	g_snap.count = count;
	g_snap.valid = 1;
}

// ============================================================================
// Commands
// ============================================================================

static SOCKET g_sock = INVALID_SOCKET;

static void driver_cmd(const char *cmd)
{
	if (g_sock == INVALID_SOCKET) { printf("ERR not connected\n"); return; }
	send(g_sock, cmd, (int)strlen(cmd), 0);
	send(g_sock, "\n", 1, 0);
	fd_set fds; FD_ZERO(&fds); FD_SET(g_sock, &fds);
	struct timeval tv = { 2, 0 };
	if (select((int)(g_sock + 1), &fds, NULL, NULL, &tv) > 0) {
		char buf[4096]; int n = recv(g_sock, buf, sizeof(buf) - 1, 0);
		if (n > 0) { buf[n] = '\0'; printf("%s", buf); }
	}
}

static void refresh_world_ptr()
{
	send(g_sock, "info\n", 5, 0);
	fd_set fds; FD_ZERO(&fds); FD_SET(g_sock, &fds);
	struct timeval tv = { 2, 0 };
	if (select((int)(g_sock + 1), &fds, NULL, NULL, &tv) > 0) {
		char buf[512] = {0}; int n = recv(g_sock, buf, sizeof(buf) - 1, 0);
		if (n > 0) {
			buf[n] = '\0';
			char *wp = strstr(buf, "world=0x");
			if (wp) { unsigned long long w; if (sscanf(wp, "world=0x%llX", &w) == 1) g_world_ptr = (uint64_t)w; }
		}
	}
}

static void cmd_get(const char *path)
{
	TypeDesc *root = find_type("WorldInternal");
	if (!root) { printf("ERR no WorldInternal type\n"); return; }
	Resolved r = resolve(g_world_ptr, root, path);
	if (!r.ok) return;
	if (r.is_array) { printf("[%d elements of %s]\n", r.count, r.type->name); return; }
	void *buf = malloc(r.type->size);
	if (!rpm_read(r.addr, buf, r.type->size)) { printf("ERR read failed\n"); free(buf); return; }
	print_struct(r.type, buf, 0);
	free(buf);
}

static void cmd_table(const char *args)
{
	char path[256] = {0};
	const char *rest = args;
	int pi = 0;
	while (*rest && *rest != ' ' && pi < 255) path[pi++] = *rest++;
	while (*rest == ' ') rest++;

	TypeDesc *root = find_type("WorldInternal");
	Resolved r = resolve(g_world_ptr, root, path);
	if (!r.ok || !r.is_array) { if (r.ok) printf("ERR not an array\n"); return; }

	// Select fields
	FieldDesc *show[32]; int show_n = 0;
	if (*rest) {
		char fn[64];
		while (*rest && show_n < 32) {
			int fi = 0;
			while (*rest && *rest != ' ' && fi < 63) fn[fi++] = *rest++;
			fn[fi] = '\0'; while (*rest == ' ') rest++;
			FieldDesc *f = find_field(r.type, fn);
			if (f) show[show_n++] = f;
			else printf("(skip: %s)\n", fn);
		}
	} else {
		for (int i = 0; i < r.type->field_count && show_n < 8; i++) {
			FieldKind k = r.type->fields[i].kind;
			if (k <= FK_ENUM || k == FK_V3 || k == FK_QUAT) show[show_n++] = &r.type->fields[i];
		}
	}
	if (!show_n) { printf("ERR no fields\n"); return; }

	int limit = r.count > 200 ? 200 : r.count;
	void *data = malloc((size_t)limit * r.type->size);
	if (!rpm_read(r.addr, data, (size_t)limit * r.type->size)) { printf("ERR read failed\n"); free(data); return; }

	printf("%6s", "#");
	for (int c = 0; c < show_n; c++) {
		if (show[c]->kind == FK_V3) printf("  %-22s", show[c]->name);
		else printf("  %14s", show[c]->name);
	}
	printf("\n");

	for (int row = 0; row < limit; row++) {
		const char *elem = (const char *)data + (size_t)row * r.type->size;
		printf("%6d", row);
		for (int c = 0; c < show_n; c++) {
			FieldDesc *f = show[c];
			const char *p = elem + f->offset;
			switch (f->kind) {
			case FK_INT:    printf("  %14d", *(int32_t *)p); break;
			case FK_UINT:   printf("  %14u", *(uint32_t *)p); break;
			case FK_FLOAT:  printf("  %14.3f", *(float *)p); break;
			case FK_DOUBLE: printf("  %14.6f", *(double *)p); break;
			case FK_BOOL:   printf("  %14s", (f->size == 1 ? *(uint8_t *)p : *(int32_t *)p) ? "true" : "false"); break;
			case FK_ENUM:   { EnumTable *et = find_enum(f->enum_name); printf("  %14s", enum_value_name(et, *(int32_t *)p)); break; }
			case FK_V3:     { float *v = (float *)p; printf("  (%6.2f,%6.2f,%6.2f)", v[0], v[1], v[2]); break; }
			case FK_QUAT:   { float *v = (float *)p; printf("  (%5.2f,%5.2f,%5.2f,%5.2f)", v[0], v[1], v[2], v[3]); break; }
			default:        printf("  %14s", "..."); break;
			}
		}
		printf("\n");
	}
	if (r.count > limit) printf("  ... (%d more)\n", r.count - limit);
	printf("(%d rows)\n", r.count);
	free(data);
}

static void cmd_summary()
{
	TypeDesc *root = find_type("WorldInternal");
	if (!root) { printf("ERR no WorldInternal\n"); return; }
	void *w = malloc(root->size);
	if (!rpm_read(g_world_ptr, w, root->size)) { printf("ERR read failed\n"); free(w); return; }

	FieldDesc *f_frame = find_field(root, "frame");
	FieldDesc *f_gravity = find_field(root, "gravity");
	FieldDesc *f_body_hot = find_field(root, "body_hot");
	FieldDesc *f_contacts = find_field(root, "debug_contacts");
	FieldDesc *f_warm = find_field(root, "warm_cache");
	FieldDesc *f_islands = find_field(root, "islands");
	FieldDesc *f_perf = find_field(root, "perf");

	if (f_frame) printf("frame: %d\n", *(int32_t *)((char *)w + f_frame->offset));
	if (f_gravity) { float *g = (float *)((char *)w + f_gravity->offset); printf("gravity: (%.1f, %.1f, %.1f)\n", g[0], g[1], g[2]); }
	if (f_body_hot) { uint64_t p = *(uint64_t *)((char *)w + f_body_hot->offset); printf("bodies: %d\n", p ? rpm_array_count(p) : 0); }
	if (f_contacts) { uint64_t p = *(uint64_t *)((char *)w + f_contacts->offset); printf("contacts: %d\n", p ? rpm_array_count(p) : 0); }
	if (f_warm) { uint64_t p = *(uint64_t *)((char *)w + f_warm->offset); printf("warm cache: %d\n", p ? rpm_map_count(p) : 0); }

	if (f_islands) {
		uint64_t ip = *(uint64_t *)((char *)w + f_islands->offset);
		TypeDesc *island_type = find_type("Island");
		FieldDesc *f_ig = find_field(root, "island_gen");
		if (ip && island_type && f_ig) {
			int n = rpm_array_count(ip);
			uint64_t gp = *(uint64_t *)((char *)w + f_ig->offset);
			if (n > 0 && gp) {
				void *islands = malloc(n * island_type->size);
				uint32_t *gens = malloc(n * 4);
				rpm_read(ip, islands, n * island_type->size);
				rpm_read(gp, gens, n * 4);
				FieldDesc *f_awake = find_field(island_type, "awake");
				int total = 0, awake = 0;
				for (int i = 0; i < n; i++) {
					if (gens[i] & 1) {
						total++;
						if (f_awake && *(int32_t *)((char *)islands + i * island_type->size + f_awake->offset)) awake++;
					}
				}
				printf("islands: %d (%d awake, %d sleeping)\n", total, awake, total - awake);
				free(islands); free(gens);
			}
		}
	}

	// Perf
	if (f_perf) {
		TypeDesc *pt = find_type("PerfTimers");
		if (pt) {
			FieldDesc *f_total = find_field(pt, "total");
			if (f_total) {
				double total = *(double *)((char *)w + f_perf->offset + f_total->offset);
				printf("perf: %.0f us\n", total * 1e6);
			}
		}
	}
	free(w);
}

static void cmd_snap(const char *path)
{
	TypeDesc *root = find_type("WorldInternal");
	Resolved r;
	if (*path) {
		r = resolve(g_world_ptr, root, path);
	} else {
		r = (Resolved){ g_world_ptr, root, 0, 0, 1 };
	}
	if (!r.ok) return;
	snap_save(r);
	printf("OK snapshot saved (%s, %d elements, %d bytes)\n", r.type->name, r.is_array ? r.count : 1, (r.is_array ? r.count : 1) * r.type->size);
}

static void cmd_diff(const char *path)
{
	if (!g_snap.valid) { printf("ERR no snapshot (use 'snap [path]' first)\n"); return; }

	// Re-read current state at the same address
	int count = g_snap.count;
	size_t bytes = (size_t)count * g_snap.type->size;
	void *cur = malloc(bytes);
	if (!rpm_read(g_snap.addr, cur, bytes)) { printf("ERR read failed\n"); free(cur); return; }

	TypeDesc *t = g_snap.type;
	int diffs = 0;
	for (int elem = 0; elem < count; elem++) {
		const char *old_p = (const char *)g_snap.data + elem * t->size;
		const char *new_p = (const char *)cur + elem * t->size;
		for (int fi = 0; fi < t->field_count; fi++) {
			FieldDesc *f = &t->fields[fi];
			if (f->kind > FK_QUAT) continue; // only diff scalars/vectors
			const char *ov = old_p + f->offset;
			const char *nv = new_p + f->offset;
			if (memcmp(ov, nv, f->size) == 0) continue;
			if (count > 1) printf("[%d].", elem);
			switch (f->kind) {
			case FK_INT:   printf("%s: %d -> %d\n", f->name, *(int32_t *)ov, *(int32_t *)nv); break;
			case FK_UINT:  printf("%s: %u -> %u\n", f->name, *(uint32_t *)ov, *(uint32_t *)nv); break;
			case FK_FLOAT: printf("%s: %.4f -> %.4f (delta %.4f)\n", f->name, *(float *)ov, *(float *)nv, *(float *)nv - *(float *)ov); break;
			case FK_DOUBLE:printf("%s: %.6f -> %.6f\n", f->name, *(double *)ov, *(double *)nv); break;
			case FK_V3:    { float *a = (float *)ov, *b = (float *)nv; printf("%s: (%.3f,%.3f,%.3f) -> (%.3f,%.3f,%.3f)\n", f->name, a[0],a[1],a[2], b[0],b[1],b[2]); break; }
			case FK_QUAT:  { float *a = (float *)ov, *b = (float *)nv; printf("%s: (%.3f,%.3f,%.3f,%.3f) -> (%.3f,%.3f,%.3f,%.3f)\n", f->name, a[0],a[1],a[2],a[3], b[0],b[1],b[2],b[3]); break; }
			default: printf("%s: changed\n", f->name); break;
			}
			diffs++;
		}
	}
	if (diffs == 0) printf("(no changes)\n");
	else printf("(%d changes)\n", diffs);
	free(cur);
}

static void cmd_filter(const char *args)
{
	char path[256] = {0}, field_name[64] = {0}, op[4] = {0}, val_str[64] = {0};
	sscanf(args, "%255s %63s %3s %63s", path, field_name, op, val_str);
	if (!*path || !*field_name || !*op || !*val_str) {
		printf("usage: filter <path> <field> <op> <value>\n  ops: > < >= <= == !=\n");
		return;
	}

	TypeDesc *root = find_type("WorldInternal");
	Resolved r = resolve(g_world_ptr, root, path);
	if (!r.ok || !r.is_array) { if (r.ok) printf("ERR not an array\n"); return; }

	FieldDesc *f = find_field(r.type, field_name);
	if (!f) { printf("ERR no field '%s' in %s\n", field_name, r.type->name); return; }

	float threshold = (float)atof(val_str);
	int limit = r.count > 10000 ? 10000 : r.count;
	void *data = malloc((size_t)limit * r.type->size);
	if (!rpm_read(r.addr, data, (size_t)limit * r.type->size)) { printf("ERR read\n"); free(data); return; }

	int matches = 0;
	for (int i = 0; i < limit; i++) {
		const char *elem = (const char *)data + (size_t)i * r.type->size;
		float val = 0;
		switch (f->kind) {
		case FK_INT:   val = (float)*(int32_t *)(elem + f->offset); break;
		case FK_UINT:  val = (float)*(uint32_t *)(elem + f->offset); break;
		case FK_FLOAT: val = *(float *)(elem + f->offset); break;
		case FK_DOUBLE:val = (float)*(double *)(elem + f->offset); break;
		default: continue;
		}

		int match = 0;
		if (strcmp(op, ">") == 0) match = val > threshold;
		else if (strcmp(op, "<") == 0) match = val < threshold;
		else if (strcmp(op, ">=") == 0) match = val >= threshold;
		else if (strcmp(op, "<=") == 0) match = val <= threshold;
		else if (strcmp(op, "==") == 0) match = fabsf(val - threshold) < 0.0001f;
		else if (strcmp(op, "!=") == 0) match = fabsf(val - threshold) > 0.0001f;

		if (match) {
			if (matches == 0) {
				printf("%6s  %14s  (all fields follow)\n", "#", field_name);
			}
			printf("%6d  %14.4f  ", i, val);
			// Print a few other scalar fields for context
			int shown = 0;
			for (int fi = 0; fi < r.type->field_count && shown < 4; fi++) {
				FieldDesc *ff = &r.type->fields[fi];
				if (ff == f) continue;
				if (ff->kind == FK_V3) { float *v = (float *)(elem + ff->offset); printf("%s=(%.1f,%.1f,%.1f) ", ff->name, v[0],v[1],v[2]); shown++; }
				else if (ff->kind == FK_FLOAT) { printf("%s=%.3f ", ff->name, *(float *)(elem + ff->offset)); shown++; }
				else if (ff->kind == FK_INT) { printf("%s=%d ", ff->name, *(int32_t *)(elem + ff->offset)); shown++; }
			}
			printf("\n");
			matches++;
		}
	}
	printf("(%d matches out of %d)\n", matches, r.count);
	free(data);
}

// Show all solver contacts for a specific body index.
static void cmd_contacts(const char *args)
{
	int body_idx = atoi(args);
	TypeDesc *root = find_type("WorldInternal");
	TypeDesc *sm_type = find_type("SolverManifold");
	TypeDesc *sc_type = find_type("SolverContact");
	if (!root || !sm_type || !sc_type) { printf("ERR missing type info\n"); return; }

	// Read solver manifold array pointer from world
	FieldDesc *f_sm = find_field(root, "dbg_solver_manifolds");
	FieldDesc *f_sc = find_field(root, "dbg_solver_contacts");
	if (!f_sm || !f_sc) { printf("ERR no solver snapshot fields (engine too old?)\n"); return; }

	void *world = malloc(root->size);
	if (!rpm_read(g_world_ptr, world, root->size)) { printf("ERR read world\n"); free(world); return; }

	uint64_t sm_ptr = *(uint64_t *)((char *)world + f_sm->offset);
	uint64_t sc_ptr = *(uint64_t *)((char *)world + f_sc->offset);
	free(world);

	if (!sm_ptr || !sc_ptr) { printf("ERR solver arrays null (step first?)\n"); return; }
	int sm_count = rpm_array_count(sm_ptr);
	if (sm_count <= 0) { printf("ERR no solver manifolds\n"); return; }

	// Read all manifolds
	void *manifolds = malloc(sm_count * sm_type->size);
	if (!rpm_read(sm_ptr, manifolds, sm_count * sm_type->size)) { printf("ERR read manifolds\n"); free(manifolds); return; }

	FieldDesc *f_ba = find_field(sm_type, "body_a");
	FieldDesc *f_bb = find_field(sm_type, "body_b");
	FieldDesc *f_cs = find_field(sm_type, "contact_start");
	FieldDesc *f_cc = find_field(sm_type, "contact_count");
	FieldDesc *f_fric = find_field(sm_type, "friction");
	FieldDesc *f_norm = find_field(sm_type, "normal");
	FieldDesc *f_lt1 = find_field(sm_type, "lambda_t1");
	FieldDesc *f_lt2 = find_field(sm_type, "lambda_t2");

	// Read all contacts
	int sc_count = rpm_array_count(sc_ptr);
	void *contacts = NULL;
	if (sc_count > 0) {
		contacts = malloc(sc_count * sc_type->size);
		rpm_read(sc_ptr, contacts, sc_count * sc_type->size);
	}

	FieldDesc *f_cn = find_field(sc_type, "normal");
	FieldDesc *f_pen = find_field(sc_type, "penetration");
	FieldDesc *f_lam = find_field(sc_type, "lambda_n");
	FieldDesc *f_em = find_field(sc_type, "eff_mass_n");
	FieldDesc *f_fid = find_field(sc_type, "feature_id");
	FieldDesc *f_soft = find_field(sc_type, "softness");

	int found = 0;
	for (int i = 0; i < sm_count; i++) {
		char *m = (char *)manifolds + i * sm_type->size;
		int ba = *(int32_t *)(m + f_ba->offset);
		int bb = *(int32_t *)(m + f_bb->offset);
		if (ba != body_idx && bb != body_idx) continue;

		int cs = *(int32_t *)(m + f_cs->offset);
		int cc = *(int32_t *)(m + f_cc->offset);
		float *norm = (float *)(m + f_norm->offset);
		float fric = *(float *)(m + f_fric->offset);
		float lt1 = f_lt1 ? *(float *)(m + f_lt1->offset) : 0;
		float lt2 = f_lt2 ? *(float *)(m + f_lt2->offset) : 0;

		printf("manifold %d: body_a=%d body_b=%d contacts=%d friction=%.3f normal=(%.3f,%.3f,%.3f) lambda_t=(%.3f,%.3f)\n",
			i, ba, bb, cc, fric, norm[0], norm[1], norm[2], lt1, lt2);

		if (contacts && f_cn && f_pen && f_lam) {
			for (int ci = 0; ci < cc; ci++) {
				int idx = cs + ci;
				if (idx >= sc_count) break;
				char *c = (char *)contacts + idx * sc_type->size;
				float *cn = (float *)(c + f_cn->offset);
				float pen = *(float *)(c + f_pen->offset);
				float lam = *(float *)(c + f_lam->offset);
				float em = f_em ? *(float *)(c + f_em->offset) : 0;
				uint32_t fid = f_fid ? *(uint32_t *)(c + f_fid->offset) : 0;
				float soft = f_soft ? *(float *)(c + f_soft->offset) : 0;
				printf("  [%d] normal=(%.3f,%.3f,%.3f) pen=%.4f lambda=%.4f eff_mass=%.4f softness=%.4f feature=0x%08X\n",
					ci, cn[0], cn[1], cn[2], pen, lam, em, soft, fid);
			}
		}
		found++;
	}
	if (!found) printf("(no manifolds for body %d)\n", body_idx);
	free(manifolds);
	free(contacts);
}

// Show warm cache entries with their body-pair keys.
static void cmd_warm(const char *args)
{
	TypeDesc *root = find_type("WorldInternal");
	TypeDesc *wm_type = find_type("WarmManifold");
	if (!root || !wm_type) { printf("ERR missing type info\n"); return; }

	void *world = malloc(root->size);
	if (!rpm_read(g_world_ptr, world, root->size)) { printf("ERR read world\n"); free(world); return; }

	FieldDesc *f_wc = find_field(root, "warm_cache");
	uint64_t map_ptr = *(uint64_t *)((char *)world + f_wc->offset);
	free(world);

	if (!map_ptr) { printf("ERR null warm cache\n"); return; }

	// Read map header
	CKMapHdr hdr;
	if (!rpm_read(map_ptr - sizeof(CKMapHdr), &hdr, sizeof(hdr))) { printf("ERR read map header\n"); return; }
	int count = hdr.size;
	if (count <= 0) { printf("(empty warm cache)\n"); return; }

	// Read items (WarmManifold array at map_ptr)
	void *items = malloc(count * wm_type->size);
	rpm_read(map_ptr, items, count * wm_type->size);

	// Read keys (uint64_t array at computed offset)
	// keys_offset = sizeof(CKMapHeader) + ALIGN8(capacity * val_size)
	size_t keys_off = sizeof(CKMapHdr) + ((hdr.capacity * hdr.val_size + 7) & ~(size_t)7);
	uint64_t keys_addr = (map_ptr - sizeof(CKMapHdr)) + keys_off;
	uint64_t *keys = malloc(count * 8);
	rpm_read(keys_addr, keys, count * 8);

	// Decode body pair from key (key = low32 | high32 << 32, where low < high)
	int filter_body = -1;
	if (*args) filter_body = atoi(args);

	FieldDesc *f_cnt = find_field(wm_type, "count");
	FieldDesc *f_stale = find_field(wm_type, "stale");
	FieldDesc *f_sat = find_field(wm_type, "sat_axis");
	FieldDesc *f_lt1 = find_field(wm_type, "manifold_lambda_t1");

	int shown = 0;
	for (int i = 0; i < count; i++) {
		int ba = (int)(keys[i] & 0xFFFFFFFF);
		int bb = (int)(keys[i] >> 32);
		if (filter_body >= 0 && ba != filter_body && bb != filter_body) continue;

		char *wm = (char *)items + i * wm_type->size;
		int cnt = f_cnt ? *(int32_t *)(wm + f_cnt->offset) : 0;
		int stale = f_stale ? *(int32_t *)(wm + f_stale->offset) : 0;
		int sat = f_sat ? *(int32_t *)(wm + f_sat->offset) : 0;
		float lt1 = f_lt1 ? *(float *)(wm + f_lt1->offset) : 0;

		printf("[%d] bodies=(%d,%d) contacts=%d stale=%d sat_axis=%d lambda_t1=%.4f\n",
			i, ba, bb, cnt, stale, sat, lt1);
		shown++;
	}
	if (!shown) printf("(no entries%s)\n", filter_body >= 0 ? " for that body" : "");
	printf("(%d shown, %d total)\n", shown, count);

	free(items);
	free(keys);
}

static void cmd_raw(const char *args)
{
	unsigned long long addr; int len = 64;
	if (sscanf(args, "0x%llX %d", &addr, &len) < 1 && sscanf(args, "%llu %d", &addr, &len) < 1) {
		printf("usage: raw <addr> [len]\n"); return;
	}
	if (len > 4096) len = 4096;
	uint8_t *buf = malloc(len);
	if (!rpm_read((uint64_t)addr, buf, len)) { printf("ERR read failed\n"); free(buf); return; }
	for (int i = 0; i < len; i += 16) {
		printf("%016llX:", addr + i);
		for (int j = 0; j < 16 && i + j < len; j++) printf(" %02X", buf[i + j]);
		printf("  ");
		for (int j = 0; j < 16 && i + j < len; j++) { char c = buf[i + j]; printf("%c", (c >= 32 && c < 127) ? c : '.'); }
		printf("\n");
	}
	free(buf);
}

// Replay a script file: sends commands to engine with pacing.
// Lines starting with # are comments. "wait <ms>" adds a delay.
static void cmd_replay(const char *filename)
{
	FILE *f = fopen(filename, "r");
	if (!f) { printf("ERR can't open %s\n", filename); return; }
	printf("Replaying %s...\n", filename);

	char line[512];
	while (fgets(line, sizeof(line), f)) {
		int len = (int)strlen(line);
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
		if (!len || line[0] == '#') continue;

		// "wait <ms>" — local delay for pacing
		if (strncmp(line, "wait ", 5) == 0) {
			int ms = atoi(line + 5);
			if (ms > 0) { printf("  (wait %dms)\n", ms); Sleep(ms); }
			continue;
		}

		printf("  > %s\n", line);
		driver_cmd(line);

		// After scene/restart, refresh world pointer
		if (strncmp(line, "scene ", 6) == 0 || strncmp(line, "restart", 7) == 0)
			refresh_world_ptr();

		// Small delay between commands so engine processes each one
		Sleep(50);
	}
	fclose(f);
	printf("Replay complete.\n");
}

// Check whether all dynamic bodies have settled (speed below threshold).
// Steps the engine N frames, checking each frame. Reports pass/fail + offenders.
static void cmd_stable(const char *args)
{
	float threshold = 0.05f;
	int frames = 60;
	sscanf(args, "%f %d", &threshold, &frames);

	TypeDesc *root = find_type("WorldInternal");
	TypeDesc *hot_type = find_type("BodyHot");
	TypeDesc *gen_type = find_type("WorldInternal"); // for body_gen
	if (!root || !hot_type) { printf("ERR missing types\n"); return; }

	FieldDesc *f_hot = find_field(root, "body_hot");
	FieldDesc *f_gen = find_field(root, "body_gen");
	if (!f_hot || !f_gen) { printf("ERR missing fields\n"); return; }

	int consecutive_stable = 0;
	int worst_body = -1;
	float worst_speed = 0;

	for (int frame = 0; frame < frames; frame++) {
		// Step one frame
		driver_cmd("step");

		// Read world to get array pointers
		void *world = malloc(root->size);
		if (!rpm_read(g_world_ptr, world, root->size)) { printf("ERR read\n"); free(world); return; }
		uint64_t hot_ptr = *(uint64_t *)((char *)world + f_hot->offset);
		uint64_t gen_ptr = *(uint64_t *)((char *)world + f_gen->offset);
		free(world);

		if (!hot_ptr || !gen_ptr) { printf("ERR null arrays\n"); return; }
		int count = rpm_array_count(hot_ptr);
		if (count <= 0) { printf("ERR no bodies\n"); return; }

		void *hots = malloc(count * hot_type->size);
		uint32_t *gens = malloc(count * 4);
		rpm_read(hot_ptr, hots, count * hot_type->size);
		rpm_read(gen_ptr, gens, count * 4);

		FieldDesc *f_vel = find_field(hot_type, "velocity");
		FieldDesc *f_angvel = find_field(hot_type, "angular_velocity");
		FieldDesc *f_inv = find_field(hot_type, "inv_mass");

		int all_stable = 1;
		worst_body = -1;
		worst_speed = 0;
		for (int i = 0; i < count; i++) {
			if (!(gens[i] & 1)) continue;
			char *h = (char *)hots + i * hot_type->size;
			float im = *(float *)(h + f_inv->offset);
			if (im == 0) continue; // skip static
			float *v = (float *)(h + f_vel->offset);
			float *av = (float *)(h + f_angvel->offset);
			float speed = sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
			float aspeed = sqrtf(av[0]*av[0] + av[1]*av[1] + av[2]*av[2]);
			// Use linear speed for stability check (angular can be large during smooth rolling)
			float vy = fabsf(v[1]);
			if (vy > threshold) {
				all_stable = 0;
				if (vy > worst_speed) { worst_speed = vy; worst_body = i; }
			}
		}
		free(hots);
		free(gens);

		if (all_stable) {
			consecutive_stable++;
			if (consecutive_stable >= 10) {
				printf("STABLE: all bodies settled after %d frames (%d consecutive stable)\n", frame + 1, consecutive_stable);
				return;
			}
		} else {
			consecutive_stable = 0;
		}
	}

	printf("UNSTABLE: after %d frames, body %d still has |vy|=%.4f (threshold=%.4f)\n", frames, worst_body, worst_speed, threshold);
}

static void cmd_help()
{
	printf(
		"Viewer commands (RPM memory inspection):\n"
		"  get <path>                          print struct at path\n"
		"  table <path> [fields...]            array as table\n"
		"  summary                             scene overview\n"
		"  contacts <body_idx>                 all solver contacts for a body\n"
		"  warm [body_idx]                     warm cache entries (optionally filter by body)\n"
		"  filter <path> <field> <op> <val>    find matches (ops: > < >= <= == !=)\n"
		"  snap [path]                         save snapshot for diff\n"
		"  diff                                compare current state with snapshot\n"
		"  replay <file>                       play back a repro script\n"
		"  stable [threshold] [frames]         check if all bodies settled\n"
		"  raw <addr> [len]                    hex dump\n"
		"\n"
		"Driver commands (sent to engine):\n"
		"  pause / step [n]                    pause, step frames\n"
		"  scene <name> / scenes / restart     scene control\n"
		"  push <idx> <f> [r]                  impulse at point\n"
		"  drag <idx> <local> <target>         begin mouse-like drag\n"
		"  dragto <target>                     update drag position\n"
		"  release                             end drag\n"
		"  highlight <idx> <color|r g b>       tint a body\n"
		"  unhighlight <idx|all>               remove tint\n"
		"  label <text>                        show overlay text\n"
		"  slow <factor>                       time scale\n"
		"\n"
		"Paths walk from WorldInternal root:\n"
		"  body_state/0  body_hot/5/velocity  islands/8  bvh_dynamic/nodes/0\n"
		"  dbg_solver_manifolds/0  dbg_solver_contacts/0  warm_cache/0\n"
	);
}

// ============================================================================
// Network + main loop
// ============================================================================

static int net_connect(const char *host, int port)
{
	WSADATA wsa; WSAStartup(MAKEWORD(2, 2), &wsa);
	struct addrinfo hints = {0}, *res;
	hints.ai_family = AF_INET; hints.ai_socktype = SOCK_STREAM;
	char ps[16]; snprintf(ps, sizeof(ps), "%d", port);
	if (getaddrinfo(host, ps, &hints, &res) != 0) return 0;
	g_sock = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
	if (g_sock == INVALID_SOCKET) { freeaddrinfo(res); return 0; }
	if (connect(g_sock, res->ai_addr, (int)res->ai_addrlen) == SOCKET_ERROR) {
		closesocket(g_sock); g_sock = INVALID_SOCKET; freeaddrinfo(res); return 0;
	}
	freeaddrinfo(res);
	u_long nb = 1; ioctlsocket(g_sock, FIONBIO, &nb);
	return 1;
}

// Receive all available data with a timeout. Returns malloc'd string (caller frees).
static char *recv_all(int timeout_ms)
{
	char *buf = NULL;
	int len = 0, cap = 0;
	int deadline = timeout_ms;
	while (deadline > 0) {
		fd_set fds; FD_ZERO(&fds); FD_SET(g_sock, &fds);
		struct timeval tv = { 0, 50000 }; // 50ms poll
		if (select((int)(g_sock + 1), &fds, NULL, NULL, &tv) > 0) {
			if (len + 4096 > cap) { cap = len + 8192; buf = realloc(buf, cap); }
			int n = recv(g_sock, buf + len, 4096, 0);
			if (n > 0) { len += n; deadline = 200; } // reset short deadline after data
			else if (n == 0) break;
		}
		deadline -= 50;
	}
	if (buf) buf[len] = '\0';
	return buf;
}

int main(int argc, char *argv[])
{
	const char *host = "localhost";
	int port = 9999;
	if (argc > 1) host = argv[1];
	if (argc > 2) port = atoi(argv[2]);

	printf("Connecting to %s:%d...\n", host, port);
	if (!net_connect(host, port)) { fprintf(stderr, "Failed to connect\n"); return 1; }

	// Receive banner + type data
	char *init_data = recv_all(3000);
	if (!init_data) { fprintf(stderr, "No data from engine\n"); return 1; }

	// Parse banner line
	unsigned long pid; unsigned long long world;
	if (sscanf(init_data, "nudge pid=%lu world=0x%llX", &pid, &world) != 2) {
		fprintf(stderr, "Bad banner: %.80s\n", init_data); free(init_data); return 1;
	}
	g_pid = (DWORD)pid;
	g_world_ptr = (uint64_t)world;

	// Parse types (starts after first newline)
	char *types_start = strchr(init_data, '\n');
	if (!types_start || !parse_types(types_start + 1)) {
		fprintf(stderr, "Failed to parse type data\n"); free(init_data); return 1;
	}
	free(init_data);

	g_proc = OpenProcess(PROCESS_VM_READ | PROCESS_QUERY_INFORMATION, FALSE, g_pid);
	if (!g_proc) { fprintf(stderr, "OpenProcess failed for PID %lu\n", g_pid); return 1; }

	printf("Connected: PID=%lu world=0x%llX (%d types, %d enums)\n\n", g_pid, g_world_ptr, g_type_count, g_enum_count);

	char line[1024];
	for (;;) {
		printf("> ");
		fflush(stdout);
		if (!fgets(line, sizeof(line), stdin)) break;
		int len = (int)strlen(line);
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
		if (!len) continue;
		if (strcmp(line, "quit") == 0 || strcmp(line, "exit") == 0) break;

		char cmd[32]; const char *rest = line;
		int ci = 0;
		while (*rest && *rest != ' ' && ci < 31) cmd[ci++] = *rest++;
		cmd[ci] = '\0'; while (*rest == ' ') rest++;

		if (strcmp(cmd, "get") == 0) cmd_get(rest);
		else if (strcmp(cmd, "table") == 0) cmd_table(rest);
		else if (strcmp(cmd, "summary") == 0) cmd_summary();
		else if (strcmp(cmd, "contacts") == 0) cmd_contacts(rest);
		else if (strcmp(cmd, "warm") == 0) cmd_warm(rest);
		else if (strcmp(cmd, "filter") == 0) cmd_filter(rest);
		else if (strcmp(cmd, "snap") == 0) cmd_snap(rest);
		else if (strcmp(cmd, "diff") == 0) cmd_diff(rest);
		else if (strcmp(cmd, "raw") == 0) cmd_raw(rest);
		else if (strcmp(cmd, "help") == 0) cmd_help();
		else if (strcmp(cmd, "replay") == 0) cmd_replay(rest);
		else if (strcmp(cmd, "stable") == 0) cmd_stable(rest);
		else if (strcmp(cmd, "pause") == 0 || strcmp(cmd, "step") == 0 || strcmp(cmd, "run") == 0 ||
		         strcmp(cmd, "scene") == 0 || strcmp(cmd, "scenes") == 0 ||
		         strcmp(cmd, "restart") == 0 || strcmp(cmd, "info") == 0 ||
		         strcmp(cmd, "push") == 0 || strcmp(cmd, "playrecording") == 0 ||
		         strcmp(cmd, "drag") == 0 ||
		         strcmp(cmd, "dragto") == 0 || strcmp(cmd, "release") == 0 ||
		         strcmp(cmd, "highlight") == 0 ||
		         strcmp(cmd, "unhighlight") == 0 || strcmp(cmd, "label") == 0 ||
		         strcmp(cmd, "slow") == 0) {
			driver_cmd(line);
			if (strcmp(cmd, "scene") == 0 || strcmp(cmd, "restart") == 0)
				refresh_world_ptr();
		}
		else printf("ERR unknown command '%s'\n", cmd);
	}

	CloseHandle(g_proc);
	closesocket(g_sock);
	WSACleanup();
	return 0;
}
