# Nudge -- 3D Physics Engine

## What this is
- Rigid body physics engine focused on contact constraints
- Shape types: box, convex hull, sphere, capsule, triangle mesh, cylinder
- Written in C, no C++

## Commits
- Do NOT add `Co-Authored-By` trailers to commit messages.

## Non-negotiables
- **Locality of code**: All behavioral code for a feature lives in one file. Helpers, init, solve steps -- everything in one place. Shared struct definitions live centrally in a common header. Do not split behavioral code across multiple files.
- **C23, unity build**: all sources are `#include`d from `src/main.c`. New files must be added to that include list. Headers must ONLY be included from `main.c`, never from other headers.
- **Pointer/API posture**:
  - No defensive NULL checks unless NULL is explicitly a valid input for that API. Let invalid lifetimes crash; asserts are fine.
  - Memory/pointer bugs must be root-caused, not masked. Crashes help find bugs early.
- **Style**:
  - Tabs for indentation; spaces for alignment after initial tabs.
  - Do not write `f(void)`; use `f()`.
  - Struct/enum definitions: opening brace on its own line.
  - No emojis or extended characters in code, comments, or output.
  - Follow existing patterns/macros once established.
  - Prefer `_Generic` macros over type-prefixed functions (e.g. `is_valid(x)` not `v3_is_valid(x)`).
  - Struct naming: use underscore to separate acronym prefixes from words (e.g. `BVH_Tree` not `BVHTree`).
  - Keep expressions on one line. Do not break assignments, returns, or function call arguments across multiple lines with continuation indentation. Long one-liners are fine. Exception: matrix math blocks where rows map to visual rows (e.g. Cramer's rule, symmetric matrix multiply) should stay multi-line.
  - All `for` loops require braces unless the entire loop (header + body) fits on a single line.
  - Function signatures on a single line. No wrapping parameters across lines.
- **Handle-based API**: public types are opaque handles, e.g. `typedef struct Body { uint64_t id; } Body;`. Internally cast/lookup as needed. Distinct struct types per concept (Body, World, etc.).
- **No rsqrt hacks**: Never use `_mm_rsqrt_ss`, fast inverse sqrt, or Newton-Raphson refinement for normalization. Use `1.0f / sqrtf(x)`. The "optimization" saves nothing on modern CPUs and introduces precision loss.
- **Scope**:
  - Only dependency is `ckit.h` (vendored). Standard C library and platform APIs only.
  - Prefer minimal patches and local changes.
  - Do not add features, abstractions, or error handling beyond what was asked.

## Architecture
- **Hot/cold data split**: `split_store.h` provides parallel cold[]/hot[]/gen[]/freelist[] arrays. Hot data (position, velocity) iterated by solver. Cold data (mass, shapes, topology) accessed on demand. Generational handles via odd/even gen scheme.
- **Naming**: `create_X`/`destroy_X` for lifecycle. `obj_verb` for everything else.
- **Collision detection**: N^2 broadphase (placeholder). Narrowphase dispatch by shape pair (upper-triangle, auto-swap).
  - Hull SAT with Gauss map pruning (Gregorius GDC 2013) for hull-hull
  - GJK distance (Catto GDC 2010) with 3D tetrahedron Voronoi regions
  - Dedicated analytical code for sphere/capsule pairs
  - Pattern from lm engine: GJK for shallow hits, SAT/face-search for deep hits
- Contact constraint: Baumgarte positional correction + velocity-level impulse (placeholder)
- **Hull**: half-edge mesh (`HalfEdge`, `HullFace`, `HullPlane`). Edges stored in twin pairs at adjacent indices. Unit box hull scaled by half_extents at query time.

## Shape types (implemented)
- **Sphere**: center + radius
- **Capsule**: segment along local Y (half_height) + radius
- **Box**: unit box hull scaled by half_extents

## Shape types (planned)
- **Convex hull**: arbitrary vertex + face/edge topology
- **Triangle mesh**: static/concave geometry
- **Cylinder**: segment + radius, flat caps

## Build & platform
- CMake with C23 (MSVC `/std:clatest /experimental:c11atomics`)
- SDL3 via FetchContent for window + OpenGL 3.3 core context
- Dear ImGui (docking, v1.92.6) with dear_bindings C API, vendored in `libraries/imgui/`
- App implements `init()`, `update()`, `draw()`
- Rendering: instanced static meshes, premultiplied alpha blending, configurable directional light

## Key files
- `src/main.c` -- unity build root, platform layer, app entry
- `src/nudge.h` -- public API (handles, params, shape types)
- `src/nudge.c` -- world/body lifecycle, solver step, includes collision.c and gjk.c
- `src/collision.c` -- broadphase, SAT, all narrowphase routines, hull data
- `src/gjk.c` -- GJK distance algorithm
- `src/split_store.h` -- hot/cold parallel array abstraction
- `src/vmath.h` -- v3, quat, mat4 math
- `src/render.c` -- GL 3.3 instanced renderer
- `src/ckit.h` -- vendored data structures (dynamic arrays, strings, maps)

## Debugging
- Use the `win-debugger` MCP for debugging crashes and inspecting runtime state. Launch the Debug build under cdb, set breakpoints on nudge functions (module name is `nudge` for the app, `nudge_tests` for tests), inspect locals/structs, and walk the call stack.
- Prefer the debugger over printf debugging -- it gives full struct inspection, call stacks, conditional breakpoints, and catches access violations at the exact instruction.
- Symbol path: `C:/git/nudge/build/Debug`. Force-load symbols with `.reload /f <module>.exe`.

## Reference code
- `ref/` -- Gregorius GDC 2013 SAT reference (rnHull.h, rnSAT.cpp)
- User's lm engine at `C:\Users\randy\Downloads\bitbucket_things\...\src\lm` -- collision dispatch, GJK, shape routines

## Testbed
- `testbed/` -- cross-engine comparison harness. Runs Jolt, Bepu, and nudge together in a single process so identical scenes can be stepped side-by-side for correctness comparison and basic perf captures.
- Written in C# (not C) specifically so Bepu (managed) is first-class; Jolt and nudge are exposed via native DLLs (`testbed/native/jolt_dll.cpp`, `testbed/native/nudge_dll.c`) and P/Invoked from the C# harness.
- Layout: `testbed/src/Testbed.{Common,Bepu,Jolt,Nudge,Visual}` projects + `Testbed.sln`. Native side built via `testbed/native/CMakeLists.txt`.
- The C# language choice for the testbed is deliberate and does not relax the "C, no C++" rule for the engine itself.
