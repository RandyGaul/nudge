# cdb.exe practice findings

Notes from actual cdb sessions against `nudge_tests.exe` Debug build. Every claim below was verified by running the command and observing output (no guessing).

## Environment

- `cdb.exe`: `C:/Program Files (x86)/Windows Kits/10/Debuggers/x64/cdb.exe`, v10.0.26100.7705.
- Target: `build/Debug/nudge_tests.exe` (unity build — all physics symbols under module `nudge_tests`).
- Bash/MSYS shell. cdb output has lots of startup noise (module loads, initial int-3). Filter with `.echo MARKER` + `sed -n '/MARKER/,/quit:/p'`.

## First five minutes — the essentials

These four lines, run at session start, fix more problems than any other combo:

```
.lines -e             ; enable source-line info (DEFAULT IS OFF — huge gotcha)
l+t                   ; step by source line (not assembly)
l+s                   ; display source at each prompt
bu <module>!<func>    ; pending breakpoint, survives if module not yet loaded
```

Without `.lines -e`, `p` and `t` step by *machine instruction*. You'll execute 15 `p` to cross one C statement. This trips up everybody on day 1.

Without `l+s`, you get only the assembly at each stop. With it, every step shows a source-line header like `>  507:   w->frame++;` above the disassembly. Enormous readability win.

## One-shot scripted invocation

```bash
cdb -c ".lines -e; l+t; l+s; bu nudge_tests!world_step 0n106; g; dv; q" \
    ./build/Debug/nudge_tests.exe
```

- `-c "<semicolons>"` is the startup command string. Wrap in double quotes — the shell eats unquoted `;`.
- End with `q` to exit. Without `q`, cdb drops to interactive prompt after running.
- Always put a `.echo AT_SOMETHING` marker before the data you want, then post-filter stdout with sed/grep. Saves you from ~200 lines of module-load noise.

## Breakpoints

| Form | Example | Use |
|------|---------|-----|
| `bp` | `bp nudge_tests!world_step` | Resolved NOW. Fails silently if module isn't loaded yet. |
| `bu` | `bu nudge_tests!world_step` | Unresolved / pending. Survives module-load timing. **Default choice.** |
| `bp <func> 0n<N>` | `bu nudge_tests!world_step 0n106` | Pass-count: skip N-1 hits, stop on the Nth. `0n` prefixes decimal. |
| Conditional (MASM) | `bu <func> "j (<MASM expr>) '';'gc'"` | `j` = branch; `''` = stop (fall through to prompt); `gc` = go-from-conditional = continue. |
| Conditional (C++) | `bu <func> "j (@@c++(<C++ expr>)) '';'gc'"` | MASM doesn't have `&&`, `||`, `->`, etc. `@@c++(...)` gives you C++. |

### Verified conditional recipe: stop at frame N

```
bu nudge_tests!world_step \
   "j (@@c++(@rcx != 0 && ((WorldInternal*)@rcx)->frame == 105)) '.echo TARGET';'gc'"
```

- `@rcx` = first integer-class arg under x64 Windows calling convention. For `world_step(World world, float dt)`, `World` is a 1-field struct `{ u64 id; }` passed by value in `@rcx`. So `@rcx` directly holds `world.id`, which is a `WorldInternal*`.
- The `@rcx != 0` guard is necessary: `bu` fires the first time before cdb has valid args, and accessing `((WorldInternal*)0)->frame` is a segfault → cdb treats exception as "true" → bp stops incorrectly on iteration 0. Guard fixes this.
- `.echo TARGET` is for scriptability — lets you `sed -n '/TARGET/,/quit:/p'` downstream.

### Pass-count bp — the caveat

`bu nudge_tests!world_step 0n106` stops on the 106th call from *process start*. In a test binary that runs many tests, that count spans every `world_step` across every earlier test — the 106th hit may be in `test_capsule_settles_on_floor`, not your bench. **Always verify** by inspecting `frame` / `sub_steps` / etc. at the stop. If it's wrong, scope down by breaking first on a bench-specific entry function, then counting from there.

## Inspection — the commands that earn their keep

### `??` auto-expands whatever you point it at

`?? w` on a pointer prints type + offsets + values for every field, recursively for nested pointers:

```
0:000> ?? w
struct WorldInternal * 0x00000192`541b11c0
   +0x000 frame            : 0n106
   +0x010 gravity          : v3
   +0x020 body_cold        : 0x00000192`541b1410 BodyCold
   +0x028 body_hot         : 0x00000192`541b48a0 BodyHot
   ...
   +0x110 contact_hertz    : 60
   +0x118 max_push_velocity : 3
   +0x11c sub_steps        : 0n4
```

Way more useful than `dt` — it follows pointers and prints values without a separate command per field. Works even on uninit pointers (`0xcccccccc`), in which case it shows the type layout with `????` values.

### `??` on a direct expression prints the value

```
?? ((WorldInternal*)@rcx)->frame         -> int 0n105
?? ((WorldInternal*)@rcx)->sub_steps     -> int 0n4
?? ((WorldInternal*)@rcx)->velocity_iters -> int 0n10
```

This is what you want at Phase-3 breakpoints when locals aren't yet assigned. Work through the register, not the still-uninitialized local.

### Locals at function entry are garbage

At the first instruction of a function, all stack locals = `0xCCCCCCCC` (MSVC Debug uninit-canary). `dv /V /t` will show the variable names and offsets but the values are lies.

Fix: after hitting the bp, do `p 6` (or however many source lines to get past the prologue assignments) before inspecting `?? w`. Or, as above, bypass with `?? ((WorldInternal*)@rcx)->field` which reads via the register.

### Stack trace with source paths

```
kp
```

With `.lines -e` set, each frame shows full file:line. Example:
```
Child-SP          RetAddr               Call Site
00000026`9d3fcdd8 00007ff7`f7a19505     nudge_tests!world_step(struct World world = struct World, float dt = -107374176) [C:\git\nudge\src\nudge.c @ 502]
00000026`9d3fcde0 00007ff7`f7af88e9     nudge_tests!test_capsule_settles_on_floor(void)+0x545 [C:\git\nudge\src\tests.c @ 347]
```

Note: first arg `World world = struct World` shows the *type* but not the `.id` field — cdb's struct-by-value display is shallow in stack-trace context. Use `??` for contents.

### Source-line stepping that actually works

After `.lines -e; l+t; l+s`, each `p` step prints:

```
>  507:     w->frame++;
nudge_tests!world_step+0x75:
00007ff7`f7939a05 488b442458      mov     rax,qword ptr [rsp+58h] ss:000000fb`bf5fb8f8=0000028fea437780
```

So you get: source-line arrow + C statement, then the asm about to execute. Perfect for "read the code while watching registers update".

### x64 Windows calling convention — quick reference

| Arg index | Integer/pointer | Float |
|----------:|-----------------|-------|
| 1 | `@rcx` | `@xmm0` |
| 2 | `@rdx` | `@xmm1` |
| 3 | `@r8`  | `@xmm2` |
| 4 | `@r9`  | `@xmm3` |
| 5+ | `[rsp+0x28 + 8*N]` | same |

For `narrowphase_pair(WorldInternal* w, int i, int j, InternalManifold** manifolds)`:
`@rcx` = w, `@edx` = i, `@r8d` = j, `@r9` = manifolds.

For `world_step(World world, float dt)`:
`@rcx` = world.id, `@xmm1` = dt.

### Useful commands summary

| Goal | Command | Notes |
|------|---------|-------|
| Show vars in current frame | `dv /V /t` | `/V` addresses, `/t` types |
| Evaluate any C expression | `?? expr` | Single most useful command |
| Dump struct by address | `dt ModuleName!StructType <addr>` | Explicit type path |
| Stack with params + lines | `kp` | Source file:line if `.lines -e` |
| Select a frame | `.frame /c <n>` | `/c` also switches register context |
| Step over (source) | `p` | Requires `.lines -e` + `l+t` |
| Step into | `t` | Descends into calls |
| Run to function end | `gu` | Faster than many `p` through body |
| Enumerate symbols | `x <module>!<pattern>` | Wildcards work; some patterns return empty though |
| Set symbol server | `.symfix; .reload` | Falls back to MS symbol server |

## DBG_BREAK + cdb combo (the production workflow)

`DBG_BREAK(name, world)` is defined in `src/debug_server.c:999`:

```c
#define DBG_BREAK(name, world) do { if (dbg_break_match(name)) dbg_break_wait_impl((name), __FILE__, __LINE__, (world)); } while (0)
```

`dbg_break_wait_impl` publishes globals and spins in a `Sleep(5) ; debug_server_poll()` loop until a viewer sends `continue`. **That spin loop is a perfect cdb attach point.**

Workflow:

1. Run test with `--debug --break=<pattern>` so DBG_BREAK arms at call sites matching `<pattern>`.
2. Test hits a `DBG_BREAK("tunnel:first", world)` site → publishes globals, sleeps.
3. Attach cdb to the running process (`cdb -p <pid>` or `cdb -pn nudge_tests.exe`).
4. In cdb: inspect anything via `??`. Use `g_dbg_world`, `g_dbg_break_name`, `g_dbg_break_file`, `g_dbg_break_line` to know where you are.
5. To resume without a viewer: `ed nudge_tests!g_dbg_break_resume 1; g`. (Writes 1 to the flag, then continues.)
6. Next `DBG_BREAK` site pauses again.

This beats both pure-cdb and pure-viewer for live inspection, because you get cdb's **local variables and stack control** (which the viewer can't see via RPM) PLUS the pause is scheduled from inside your source (no flaky frame-counter breakpoint).

## Verified gotchas

| Symptom | Cause | Fix |
|---------|-------|-----|
| `p` takes 20 hits to cross one C line | Line info disabled | `.lines -e` at session start |
| Conditional bp stops every time | MASM expression parse error → treated as true | Use `@@c++(...)`; test predicate manually first |
| `Numeric expression missing from '& ...'` | MASM has no `&&` or `||` | Switch to `@@c++(...)` for boolean expressions |
| `Memory access error at '...'` | Dereffed an uninitialized local | Use `@rcx` or `p` past prologue first |
| Silent bp miss | Used `bp` before module loaded | Use `bu` for pending |
| All locals `0xCCCCCCCC` | At function entry, pre-assignments | `p` source-step past prologue |
| Pass-count bp stops in wrong test | Counter is global across tests | Scope via a bench-specific bp first |
| `*** WARNING: Unable to verify checksum` | No Authenticode sig on exe | Ignore — harmless |
| `x nudge_tests!*pattern*` returns empty | Not all function names match through all wildcards | Try without `*...*`, or use a more specific substring |
| `-c "cmd1; cmd2"` runs only first | Shell ate the `;` | Always double-quote the full `-c` arg |

## Recipes I actually use

### Snapshot at a target frame (non-interactive)

```bash
cdb -c ".lines -e; l+t; l+s;\
 bu nudge_tests!world_step \"j (@@c++(@rcx != 0 && ((WorldInternal*)@rcx)->frame == 105)) '.echo AT_BP';'gc'\";\
 g; ?? ((WorldInternal*)@rcx)->frame; ?? ((WorldInternal*)@rcx)->sub_steps; kp; q" \
 ./build/Debug/nudge_tests.exe 2>&1 | sed -n '/AT_BP/,/quit:/p'
```

Output:
```
AT_BP
int 0n105
int 0n4
Child-SP          RetAddr               Call Site
...world_step [C:\git\nudge\src\nudge.c @ 502]
...test_capsule_settles_on_floor+0x545 [C:\git\nudge\src\tests.c @ 347]
```

### Step-through the bad frame

Once `AT_BP` fires (at world_step entry for frame-before-bug):

```
p                                  ; advance one source line
?? w                               ; inspect w once assigned
bu nudge_tests!collide_cylinder_hull   ; add a bp inside the suspicious function
g                                  ; continue into frame 106's narrowphase
?? a.center; ?? b.hull->face_count ; spot-inspect first args of collide_cylinder_hull
```

### Body-specific break on narrowphase

For `narrowphase_pair(WorldInternal* w, int i, int j, InternalManifold** manifolds)`:

```
bu nudge_tests!narrowphase_pair "j (@@c++(@edx == 45 || @r8d == 45)) '';'gc'"
```

Stops only on pairs where body 45 is either A or B.

### Dumping a Manifold with its Contact array

With some `Manifold* m` in scope:

```
?? m->count
?? m->contacts[0]
?? m->contacts[1]
```

Or for all at once:
```
dt nudge_tests!Contact @@c++(&m->contacts[0]) -a4
```
(`-a4` = array of 4. Requires Contact to be a known type in the PDB — it is for this project.)

## What I won't reach for

- **`dt` as primary inspection.** `??` is almost always better — follows pointers, shows values.
- **Bare `bp` on non-loaded symbols.** Always `bu`.
- **Interactive stepping through solver loops.** PGS iterates over thousands of rows; stepping each is hopeless. Set a bp inside the row body with a predicate matching the target constraint.
- **MASM conditionals with multiple terms.** Just use `@@c++(...)`.

## Summary: cdb's role in the physics-debug skill

- **When you know the function and want state:** use `??` after one source-step past prologue. 30 seconds.
- **When you need to catch a specific frame:** conditional bp with `@@c++(... && frame == N)`. Verified working above.
- **When you need to watch control flow:** `l+t; l+s; p` loop. Source + asm side-by-side is productive.
- **When you want deep reflected state with zero rebuild:** attach to a process paused at `DBG_BREAK` via the viewer's inspection commands. cdb is for the locals/stack that RPM can't reach.
