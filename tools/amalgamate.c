// tools/amalgamate.c -- generate single-file header from nudge source tree.
//
// Usage: amalgamate [-o file] [-x name]... headers... -- impl_files...
//
// Recursively inlines #include "..." from seed files into a single-file header.
// Headers go before #ifdef NUDGE_IMPLEMENTATION, impl files go inside it.
// -x keeps a filename as a literal #include instead of inlining it.
// Already-inlined files are skipped across both phases (dedup).
//
// Example (nudge):
//   amalgamate -o nudge.h src/nudge.h -- src/nudge_amalg.c

#define CKIT_IMPLEMENTATION
#include "../src/ckit.h"
#include <stdio.h>
#include <string.h>

static CK_DYNA char** g_done = NULL;
static CK_DYNA char** g_excl = NULL;

static int list_has(CK_DYNA char** list, const char* s)
{
	for (int i = 0; i < asize(list); i++)
		if (strcmp(list[i], s) == 0) return 1;
	return 0;
}

static void process(FILE* out, const char* path)
{
	if (list_has(g_done, path)) return;

	int n = (int)strlen(path);
	char* saved = CK_ALLOC(n + 1);
	memcpy(saved, path, n + 1);
	apush(g_done, saved);

	FILE* f = fopen(path, "r");
	if (!f) { fprintf(stderr, "amalgamate: %s: not found\n", path); return; }

	// Directory of this file for resolving relative includes.
	char dir[512] = ".";
	const char* sl = strrchr(path, '/');
	if (!sl) sl = strrchr(path, '\\');
	if (sl) { int len = (int)(sl - path); memcpy(dir, path, len); dir[len] = 0; }

	char line[4096];
	while (fgets(line, sizeof(line), f)) {
		char inc[256];
		if (sscanf(line, " #include \"%255[^\"]\"", inc) == 1) {
			if (list_has(g_excl, inc)) { fputs(line, out); continue; }
			char resolved[512];
			snprintf(resolved, sizeof(resolved), "%s/%s", dir, inc);
			process(out, resolved);
		} else {
			fputs(line, out);
		}
	}
	fclose(f);
}

static void emit_preamble(FILE* out)
{
	fprintf(out,
		"/*\n"
		"   ------------------------------------------------------------------------------\n"
		"      Licensing information can be found at the end of the file.\n"
		"   ------------------------------------------------------------------------------\n"
		"\n"
		"   nudge.h\n"
		"\n"
		"   To create implementation (the function definitions)\n"
		"      #define NUDGE_IMPLEMENTATION\n"
		"   in *one* C/CPP file (translation unit) that includes this file\n"
		"\n"
		"\n"
		"   SUMMARY\n"
		"\n"
		"      Nudge is a lightweight 3D rigid-body physics engine for games and\n"
		"      prototyping. Written in C23 with zero external dependencies (beyond the\n"
		"      C standard library). All bodies, shapes, joints, and the world are created\n"
		"      and controlled through simple opaque-handle APIs.\n"
		"\n"
		"\n"
		"   FEATURES\n"
		"\n"
		"      - Shape types: sphere, capsule, box, convex hull\n"
		"      - Collision detection: GJK distance, SAT with Gauss map pruning,\n"
		"        analytical sphere/capsule pairs, Sutherland-Hodgman contact clipping\n"
		"      - Multiple solver backends:\n"
		"          SOLVER_SOFT_STEP  -- soft contacts with sub-step relaxation (default)\n"
		"          SOLVER_SI         -- sequential impulse with NGS position correction\n"
		"          SOLVER_BLOCK      -- direct LCP enumeration for 2-4 contact normals\n"
		"          SOLVER_AVBD       -- augmented vertex block descent (position-level)\n"
		"      - Joints: ball-socket, distance (both rigid or spring-damper)\n"
		"      - BVH broadphase with 64-byte cache-line nodes and incremental SAH\n"
		"        refinement, plus sleep-aware dirty tracking\n"
		"      - Incremental island-based sleeping\n"
		"      - Warm starting with contact feature IDs\n"
		"      - Hot/cold data split for cache-friendly solver iteration\n"
		"\n"
		"\n"
		"   QUICK START\n"
		"\n"
		"      World world = create_world((WorldParams){ .gravity = {0, -9.81f, 0} });\n"
		"\n"
		"      Body floor = create_body(world, (BodyParams){ .mass = 0 });\n"
		"      body_add_shape(world, floor, (ShapeParams){\n"
		"          .type = SHAPE_BOX, .box.half_extents = {10, 0.5f, 10}\n"
		"      });\n"
		"\n"
		"      Body ball = create_body(world, (BodyParams){ .mass = 1, .position = {0,5,0} });\n"
		"      body_add_shape(world, ball, (ShapeParams){\n"
		"          .type = SHAPE_SPHERE, .sphere.radius = 0.5f\n"
		"      });\n"
		"\n"
		"      while (running) {\n"
		"          world_step(world, 1.0f / 60.0f);\n"
		"          v3 pos = body_get_position(world, ball);\n"
		"      }\n"
		"      destroy_world(world);\n"
		"\n"
		"\n"
		"   CREDITS\n"
		"\n"
		"      BEPUphysics v2 (Ross Nordby) -- broadphase, soft constraints, memory layout\n"
		"      Box2D (Erin Catto) -- GJK, sequential impulse, contact feature IDs\n"
		"      Dirk Gregorius -- SAT with Gauss map, contact clipping, block solver\n"
		"      Chris Giles -- AVBD solver reference\n"
		"*/\n"
		"\n"
	);
}

int main(int argc, char** argv)
{
	const char* outpath = NULL;
	CK_DYNA char** hdrs = NULL;
	CK_DYNA char** impl = NULL;
	int phase = 0;

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) outpath = argv[++i];
		else if (strcmp(argv[i], "-x") == 0 && i + 1 < argc) { char* x = argv[++i]; apush(g_excl, x); }
		else if (strcmp(argv[i], "--") == 0) phase = 1;
		else if (!phase) { char* h = argv[i]; apush(hdrs, h); }
		else { char* im = argv[i]; apush(impl, im); }
	}

	FILE* out = outpath ? fopen(outpath, "w") : stdout;
	if (!out) { fprintf(stderr, "amalgamate: can't write %s\n", outpath); return 1; }

	emit_preamble(out);

	fprintf(out, "#ifndef NUDGE_H\n#define NUDGE_H\n\n");
	fprintf(out, "#include <stdint.h>\n");
	fprintf(out, "#include <stddef.h>\n");
	fprintf(out, "#include <math.h>\n");
	fprintf(out, "#include <string.h>\n");
	fprintf(out, "#include <assert.h>\n\n");

	for (int i = 0; i < asize(hdrs); i++) process(out, hdrs[i]);

	fprintf(out, "\n// end of header -- implementation follows\n");
	fprintf(out, "#ifdef NUDGE_IMPLEMENTATION\n\n");
	for (int i = 0; i < asize(impl); i++) process(out, impl[i]);
	fprintf(out, "\n#endif // NUDGE_IMPLEMENTATION\n");
	fprintf(out, "\n/*\n"
		"   ------------------------------------------------------------------------------\n"
		"   This software is available under the public domain (Unlicense).\n"
		"   See LICENSE for details.\n"
		"   ------------------------------------------------------------------------------\n"
		"*/\n");
	fprintf(out, "#endif // NUDGE_H\n");

	if (outpath) fclose(out);
	return 0;
}
