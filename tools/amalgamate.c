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
//   amalgamate -o nudge.h src/ckit.h src/nudge.h -- src/nudge.c

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

	fprintf(out, "\n// ---- %s ----\n\n", path);

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

	fprintf(out, "// nudge.h -- single-file 3D physics library (generated)\n");
	fprintf(out, "// https://github.com/RandyGaul/nudge\n");
	fprintf(out, "//\n");
	fprintf(out, "// Do this in *one* C file:\n");
	fprintf(out, "//   #define NUDGE_IMPLEMENTATION\n");
	fprintf(out, "//   #include \"nudge.h\"\n\n");
	fprintf(out, "#ifndef NUDGE_SINGLE_FILE_H\n#define NUDGE_SINGLE_FILE_H\n\n");
	fprintf(out, "#include <stdint.h>\n");
	fprintf(out, "#include <stddef.h>\n");
	fprintf(out, "#include <math.h>\n");
	fprintf(out, "#include <string.h>\n");
	fprintf(out, "#include <assert.h>\n\n");

	for (int i = 0; i < asize(hdrs); i++) process(out, hdrs[i]);

	fprintf(out, "\n#ifdef NUDGE_IMPLEMENTATION\n");
	for (int i = 0; i < asize(impl); i++) process(out, impl[i]);
	fprintf(out, "\n#endif // NUDGE_IMPLEMENTATION\n");
	fprintf(out, "#endif // NUDGE_SINGLE_FILE_H\n");

	if (outpath) fclose(out);
	return 0;
}
