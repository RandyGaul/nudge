---
name: physics-debug
description: Debug physics bugs (collision, solver, joints, integration) using deterministic repro + cdb step-through. Activates when user reports bodies exploding, tunneling, not settling, false contacts, or wrong collision behavior — especially after speculative fixes failed.
version: 1.0.0
---

# Physics Debug — Deterministic Repro + cdb Step-Through

Debug physics bugs by watching the code execute one step at a time with a debugger attached, from the frame BEFORE the bug manifests. Not by speculating and patching symptoms.

## When to activate

User says things like:
- "bodies are exploding / falling through / tunneling / jittering / not settling"
- "collision isn't working right / objects fly off the mesh"
- "find the collision bug / debug physics / step through the solver"
- "why is X happening at frame Y"

And **especially** when speculative fixes haven't worked.

## When NOT to activate

- General algorithmic questions ("how does SAT work?") — that's research, not debug.
- Build errors, syntax errors, non-physics issues.
- Symptoms reproducible only in interactive scenes with mouse input (those need `remote-debug` in app mode, sometimes with `playrecording`).

## Trust the solver. Suspect its inputs.

Contact constraints in PGS are the most-stable component in the engine:

- **Simple**: one inequality row per contact, single lambda clamped ≥ 0 with accumulated-impulse clamping, friction cone capped by the normal lambda. That's it.
- **Constantly tested**: every sub-step of every scene exercises thousands of rows across many shape pairs; if a correctness bug lived here, every bench would already be on fire.
- **Globally convergent**: contact rows don't need to live in the same manifold to communicate — they share the body's velocity state, and accumulated lambdas converge across the whole island.

Bugs almost never live *inside* the solver. Blowups, launches, large angular velocities, tunneling, jitter — PGS is faithfully integrating the inputs it was handed. The actual bug is virtually always upstream:

- **Narrowphase** emitting bad contacts: wrong normal direction, inflated penetration, false-positive axis from a degenerate hull, duplicate contacts at the same world point, non-symmetric sampling on a rigid body.
- **Warm cache** returning stale lambdas, matching a new contact against a dead feature ID, or handing out a dangling pointer after a map resize.
- **Broadphase / integration** missing a frame of contact so the pair re-contacts at large penetration the next frame, or writing a NaN/garbage position the solver then tries to correct.

**This is why the skill's core move is to break one physics step *before* the visible anomaly.** You want to watch the bad data get written — not inspect the wreckage after the solver has faithfully propagated it. By the time PGS finishes the anomalous step, the symptom is everywhere and the causal chain is gone. Single-stepping from the prior frame through narrowphase → pre-solve → solver-input buffers catches the exact line where a normal flips sign, a pen goes negative, or a stale lambda leaks in.

**Corollary:** never "fix" a blowup by turning `beta`, `max_push_velocity`, `velocity_iters`, `sub_steps`, `contact_hertz`, or `contact_damping_ratio`. The solver is correct; a knob-tune that makes *this* symptom smaller is masking the upstream bug, which will resurface under different inputs.

## Debugging styles — pick the right tool

Three styles are available. Prefer them in this order unless the table below
says otherwise.

| Style | What it is | When it wins |
|-------|-----------|-------------|
| **`remote-debug` pause points** (`DBG_BREAK`) | TCP-attached viewer inspects live memory via RPM while the test/app is paused at a named breakpoint. Zero rebuild to inspect anything new. | Default. Almost all test-side bugs; any time you'd be tempted to add "just one more printf". |
| **cdb / `win-debugger` step-through** | Full debugger with call stack, step-in/step-out, registers, arbitrary expression evaluation, memory writes. | You need to watch logic *inside* a single function execute line-by-line, see the call stack leading to a bad state, or mutate locals to test a hypothesis. |
| **printf / conditional logging** | Targeted `printf` or CSV dump compiled into the test. | You need a time series across many frames, aggregate statistics, or a trace to diff between runs. |

### Why pause points beat printf for most investigations

- **No rebuild to change what you look at.** Printf commits you to one query; inspecting a different field means edit → rebuild → re-run. A pause point lets you fire `get body_state 2`, `table body_hot`, `contacts 1`, `get ldl_debug_info`, and `warm` in any order against the exact same paused state.
- **Full reflected struct access, not a hand-picked view.** `summary` / `get` / `table` traverse the entire reflected type graph. Printfs show only what you remembered to dump — and usually miss the field that actually mattered.
- **One run covers many hypotheses.** Each `continue` advances to the next break. No repeated startup cost per experiment.
- **State is not transcribed.** Printf output is a text snapshot of what a human formatted at commit time. RPM reads the bytes — no chance of formatting bugs misleading you.
- **Works on release code paths.** You don't need to leave formatting code in a hot loop; the pause point is a predicate check only.

### Why printf still wins sometimes

- **Temporal patterns.** "When does body B first tunnel?" needs every frame, not one paused frame. Printf writes a line per frame; pausing 600 times is tedious. Use printf when the bug's signature is a *trend over time*.
- **Hot-loop statistics.** Histograms, running averages, counters — build them in code, dump a summary line.
- **Non-interactive CI / comparisons.** Printf output is diff-able between runs. Pause points require a live process + operator.
- **Inspecting state you can't reach via reflection.** Stack-local variables, intermediate expressions in narrowphase SIMD code, control-flow decisions inside a function — printf (or a debugger) gets there; RPM does not.
- **You already know exactly what one value to watch.** The printf latency is lower than standing up a viewer for a single `body_state[B].position.y < −0.3` gate (but `DBG_BREAK` gated on the same predicate is usually just as cheap).

### Why cdb step-through is still needed

- **"Why did we end up in this branch?"** Reflection sees *state*, cdb sees *control flow + call stack*.
- **Breaking one frame before the bug** to watch the BAD frame execute top to bottom. Pause points pause at a known named site; cdb can stop at any instruction given the right breakpoint condition.
- **Mutating state to test a hypothesis.** "If I force `best.index1 = -1` here, does the bug vanish?" — only a real debugger does that.
- **Bugs in code that never reaches a pause point** (e.g. infinite loops, early-exit returns, crashes inside SIMD intrinsics).

A typical session uses all three: pause point to scope the investigation, cdb to step through the bad frame, printf to snapshot a time series across the fix-verification matrix.

### What each tool can actually see

The three styles are **not interchangeable** — they inspect different memory regions. Keep this table in mind when the question is "which tool do I need for *this* state?"

| State you want to inspect | Viewer (RPM) | cdb / win-debugger | printf |
|---------------------------|:------------:|:-------------------:|:------:|
| Reflected struct fields on `WorldInternal`, bodies, manifolds, joints, warm cache, LDL debug info, NP_DebugSnapshot | ✅ | ✅ | ✅ (if you printed them) |
| Arbitrary memory reads at any address | ✅ | ✅ | ❌ |
| **Function-local variables** (e.g. `best`, `face_a`, `edge_q` inside `collide_hull_hull_ex`) | ❌ | ✅ | ✅ (if you printed) |
| **Callstack** leading to current location | ❌ | ✅ | ❌ |
| **Intermediate expressions** inside a function (mid-SIMD, pre-return) | ❌ | ✅ | ✅ (if you printed) |
| **Control-flow decisions** ("did we take this branch?") | ❌ (inferred from state) | ✅ | ✅ |
| **Mutating state** to test a hypothesis ("what if this were 0?") | ❌ | ✅ | ❌ |
| Multiple hypotheses against same paused moment — no rebuild | ✅ | ✅ | ❌ (rebuild each) |
| **Time-series across N frames** | ❌ (one moment at a time) | ❌ (painful) | ✅ |
| **Diff-able output between runs** | ❌ | ❌ | ✅ |
| **CI / headless** | ❌ (needs operator) | ❌ (needs operator) | ✅ |

**Biggest trap:** thinking RPM will show you locals. It won't. The reflected type graph traverses structs reachable from `WorldInternal*`; function locals live on the thread stack at addresses RPM can't easily resolve. If the value you want is a local in a narrowphase or solver function, you need cdb — or you need to publish that local into reflected storage first (see "Deep narrowphase debugging" below).

### Deep narrowphase debugging: cdb + viewer on the same paused process

DBG_BREAK parks the calling thread in `Sleep(5)` inside `dbg_break_wait_impl`. The callstack above that sleep is still live, so you can attach cdb to the already-paused process and inspect everything cdb normally can — *plus* the viewer keeps giving you RPM-level state in parallel.

Workflow for a deep NP bug (e.g. "SAT picks the wrong axis for this box-hull pair"):

1. **Place a break point inside** (or immediately after) the narrowphase call:
   ```c
   // inside collide_hull_hull_ex_dbg, right before the successful return
   DBG_BREAK("np:hull-hull:wrong-axis", *(World*)&w_handle);
   ```
   Or — since bodies travel by index, not by pointer — use an existing entry point in the dispatch (`np_hull_hull` etc.) and read `c->body_a`/`c->body_b`.

2. **Run** with debug enabled and filtered:
   ```bash
   ./build/Debug/nudge_tests.exe --debug --break=np:* --bench-trimesh
   ```
   Stderr shows `[dbg] break "np:hull-hull:wrong-axis" at ... (world=0x...)` when the filter matches.

3. **Attach cdb to the paused PID** (grab the PID from the viewer's connect banner or the stderr line):
   ```
   cdb -p <pid>
   ~*k                              # every thread's callstack — find the one parked in dbg_break_wait_impl
   ~<nn>s                            # switch to that thread
   .frame 2                          # step up out of Sleep/dbg_break_wait_impl into your NP frame
   dv                                # list locals: face_a, face_b, edge_q, ref_face, flip, ...
   ?? face_a.separation              # evaluate any C expression
   dt FaceQuery @@c++(&face_a)       # dump a struct
   ```

4. **Run the viewer concurrently** on its own socket for the full world state:
   ```
   > summary
   > get body_state 42
   > get np_debug        # <-- reflected SAT intermediates, see below
   > contacts 42
   ```

5. **Detach cdb** (`q` or `.detach`), then `continue` from the viewer to resume the test.

Rule of thumb: when the bug is "which axis did SAT pick and why", cdb wins. When it is "what contact manifold did we hand the solver for body X", the viewer wins. A deep narrowphase session usually needs both.

### `NP_DebugSnapshot`: reflected SAT intermediates (no cdb required)

To close the "I need to see `face_a.separation` mid-NP" gap without cdb, narrowphase writes a reflected snapshot to `WorldInternal.np_debug` when enabled. The viewer reads it via `get np_debug`.

**Turn it on** (from the viewer):
```
> np-debug 1 42 -1          # capture any pair touching body 42; -1 = any for the other side
> step 1                    # (app mode) or `continue` past a DBG_BREAK in test mode
> get np_debug
```

**What it contains** (per hull-hull pair):
- `body_a`, `body_b`, `shape_a_kind`, `shape_b_kind`, `frame`
- Full SAT output: `face_a_index` + `face_a_sep`, `face_b_index` + `face_b_sep`, `edge_a_index` + `edge_b_index` + `edge_sep`
- `winning_axis`: `0` = face_a, `1` = face_b, `2` = edge, `-1` = separated / no contact produced
- Final manifold: `contact_normal`, `contact_count`, `contact_points[0..3]`, `contact_pens[0..3]`

**When it suffices (no cdb needed):**
- Which SAT axis won, and by how much over the runners-up.
- Whether the final normal is plausible (does `contact_normal` point from A to B? Is it unit-length?).
- Whether contact count and penetrations look sane.
- Comparing two runs: capture on frame 53 vs frame 54 to see the transition.

**When you still need cdb:**
- SAT picked the right axis but the *support vertex* it found is wrong — that's inside `hull_support` / `sat_eval_face_ex`, not in the snapshot.
- Hill-climb got stuck at a local minimum — you need the `sat_hint` path and intermediate face indices, which are locals.
- Contact reduction collapsed two points that shouldn't have collapsed — `reduce_contacts` loop internals aren't in the snapshot.
- Anything in GJK / EPA / capsule code — the snapshot only covers hull-hull SAT today (the worst offender for local-state-heavy bugs). Extend it if a new hot spot appears.

**Caveat:** single-slot snapshot. Under parallel narrowphase (`--threads >1`) the slot is raced; force `--threads 1` when using it.

**Zero-cost when off:** the dispatch agents take one predicate branch (`w->np_debug_enabled`) and skip. Safe to ship.

## Core methodology

1. **Deterministic repro first.** Same input + same binary = same numbers every run.
2. **Minimize** to one body, one frame, one anomaly, whenever possible. You want a minimal repro to isolate variables.
3. **Inspect before stepping.** Land a `DBG_BREAK` pause point at/just-before the anomaly, attach the viewer, and read the full paused state via RPM (`summary`, `get`, `table`, `contacts`, `warm`). This is the new default — it replaces most exploratory printfs.
4. **Attach cdb one physics step BEFORE the bug** when you also need to watch *control flow* inside the bad frame. Post-mortem state is already corrupted; watch the bad step unfold.
5. **Step-through with a hypothesis.** Every breakpoint has an expected/actual comparison.
6. **Fix the root cause** in the suspected module. Minimal change, try to isolate variables. You may make surgical changes to gather information to continue root-cause analysis.
7. **Verify with varied inputs**, not just the one repro, but verify with similar setups, variational inputs.

## The five phases

### Phase 1 — Deterministic repro

**Rule:** You cannot debug what you cannot reproduce exactly. "Sometimes fails" = undebuggable.

- Use a bench / test binary, not an interactive scene. Fixed script, fixed initial state, no mouse.
- For nudge: `nudge_tests.exe --bench-trimesh`, `--bench-ragdoll`, etc. Or write a one-off minimal test.
- Confirm determinism: run the same command 2× and diff the output. If it differs, fix determinism first.

```bash
for run in 1 2; do
  ./build/Release/nudge_tests.exe --bench-trimesh > run-$run.txt 2>&1
done
diff run-1.txt run-2.txt   # must be empty
```

**Minimize** until the repro is the smallest configuration that still shows the bug. Bodies: 100→50→10→1. Shapes: mixed→single type. Mesh: simplify until bug disappears, then back off one change.

**Output:** a repro that fires the bug at a specific (body B, frame F), with pre-frame state recorded.

### Phase 2 — Instrument to find the moment

Pinpoint (body, frame) programmatically. **Default to `DBG_BREAK` pause points.**
Reach for printf or conditional CSV dumps only when you need a trend across
many frames.

#### 2a. Default: guarded pause point (preferred)

Insert a `DBG_BREAK` conditional on the anomaly. The first time the predicate
trips, the test thread blocks and waits for a viewer. Then you can inspect
*anything* about the paused state — no rebuild, no pre-committed print format.

```c
// In the test loop, body index `bi`, frame `f`
if (wi->body_state[bi].position.y < -0.3f) {
    DBG_BREAK("tunnel:first", *(World*)&w);   // first offender caught
}
```

Run the test with the debugger server up:
```bash
./build/Debug/nudge_tests.exe --debug --break=tunnel:* --unit
# or for a specific test-only binary that has bench flags:
./build/Debug/nudge_tests.exe --debug --break=tunnel:* --bench-trimesh
```

Attach the viewer from another shell:
```bash
./build/Debug/nudge_viewer.exe localhost 9999
> where                 # file:line of the pause
> summary               # frame, body count, contacts, islands
> get body_state        # positions + rotations of every body
> get body_hot 85       # velocities of the target body
> contacts 85           # its solver contacts
> table body_hot        # all velocities side-by-side
> continue              # resume until next match (or test end)
```

This gives you the full pre-frame picture without committing to any
particular printf format ahead of time.

#### 2b. Fallback: conditional printf / CSV (trend over time)

When the bug's signature is "the value drifts over N frames" or you want to
diff two runs, printf wins. Keep it minimal — one event, one line, stop.

```c
// In the bench loop
if (first_blowup_frame < 0 && wi->body_state[bi].position.y < -0.3f) {
    first_blowup_frame = f;
    first_blowup_body  = bi;
    // snapshot prev_pos/prev_vel for the frame BEFORE, then set the flag so
    // this block doesn't fire again. Avoid flooding stdout — you want one hit.
}
```

Iterate the threshold until you catch the earliest manifestation — a body
penetrating the mesh by 0.01 is a better target than one already at y=-100.

**Output of this phase:** concrete state 1 frame BEFORE the bug, e.g. *"frame
53, body 85 at (−0.73, 0.02, 4.54) vel (−1.42, −8.05, 0.07); frame 54 it will
be at (−0.76, −0.12, 4.54) vel (−1.49, −8.17, 0.08)."* Record it either by
reading `get`/`table` at the pause or by printf — whichever got you here.

### Phase 3 — Attach cdb one frame BEFORE

**Critical insight:** if you break AT the bad frame, state is already corrupted. You can see the wreckage but not the cause. Break on the PREVIOUS frame and watch THAT frame execute.

#### Build Debug

```bash
cmake --build build --config Debug --target nudge_tests
```

#### Choose your entry path

Three ways to land cdb inside the bad frame, each with tradeoffs:

| Entry path | Good for | Tradeoff |
|-----------|----------|----------|
| **Attach to a DBG_BREAK'd process** (preferred if you've already scoped in Phase 2) | Your Phase-2 pause point fired on the right (body, frame). cdb shows full callstack + locals at that exact instant. No conditional-bp gymnastics. | Need to grab the PID, attach another tool, detach before `continue`. See "Deep narrowphase debugging" above for the recipe. |
| **Launch under cdb with conditional bp** | No DBG_BREAK yet; you want to hit a specific function the first time a frame counter matches. | cdb has to evaluate the condition on every hit; slow if the bp is in a hot path. |
| **Win-debugger MCP** (preferred over shelling out when tooling is available) | Scripted debug loops — lets you chain `launch` / `set_breakpoint` / `go` / `locals` / `eval` from the agent layer. | Only useful when you have the MCP available. |

Fallback: direct cdb invocation
```bash
"/c/Program Files (x86)/Windows Kits/10/Debuggers/x64/cdb.exe" \
  ./build/Debug/nudge_tests.exe --bench-trimesh
```

Symbols: symbol path is `C:/git/nudge/build/Debug`. Force-load with `.reload /f nudge_tests.exe`.

Module names: `nudge` (app) or `nudge_tests` (test binary).

#### First five minutes in cdb — mandatory session setup

These four things fix more problems than anything else. Run them at the start of every session:

```
.lines -e          ; enable source-line info — DEFAULT IS OFF
l+t                ; step by source line (not assembly)
l+s                ; show source at every prompt
bu <mod>!<func>    ; pending breakpoint, survives module-load timing
```

Without `.lines -e`, `p` and `t` step one *assembly instruction* at a time — you'll execute 15 `p` per C statement and not know why. With `l+s`, each step prints the C line (`>  507:   w->frame++;`) above the asm. Without `bu` (vs `bp`), your breakpoint silently misses if the module isn't fully loaded.

#### Target a specific frame: verified conditional bp recipe

```
bu nudge_tests!world_step "j (@@c++(@rcx != 0 && ((WorldInternal*)@rcx)->frame == 0n105)) '';'gc'"
```

- `@rcx` holds `world.id` at entry — `World` is a 1-field struct `{ u64 id; }` passed by value in RCX. So `((WorldInternal*)@rcx)` is the world pointer *before* the function body runs (when locals are still `0xCCCCCCCC`).
- `@rcx != 0` guard is **required**: cdb evaluates the predicate on an early dummy hit; dereffing NULL throws an exception which cdb treats as "true" and stops spuriously. The guard skips that.
- `@@c++(...)` for `&&` — MASM's expression parser doesn't support `&&`/`||`. Anything with compound boolean logic must go through C++ mode.
- `0n105` = decimal 105. cdb's default base is hex.
- `'';'gc'` = "when true, stop (empty command = fall through to prompt); when false, `gc` = continue from conditional".

Stop on frame 105 means the bad frame (106) is about to execute — watch it unfold with `p`.

#### Target a specific body inside a hot function

Don't break at `world_step` — you'll stop 900 times. Narrow to the suspected function:
- Narrowphase: `collide_cylinder_hull`, `collide_hull_hull_ex`, `sat_query_edges`
- Constraint setup: `pre_solve_manifold`
- Solver iteration: `solver_pgs_iterate_island`
- Integration: `integrate_positions`, `integrate_velocities`

Match on the arg register for x64 Windows calling convention. For `narrowphase_pair(WorldInternal* w, int i, int j, ...)`, `i` is `@edx`, `j` is `@r8d`:

```
bu nudge_tests!narrowphase_pair "j (@@c++(@edx == 0n45 || @r8d == 0n45)) '';'gc'"
```

Arg register cheat sheet (Windows x64): 1=`@rcx`, 2=`@rdx`, 3=`@r8`, 4=`@r9`; floats in `@xmm0..@xmm3`. 32-bit int args use the low-32 aliases `@ecx/@edx/@r8d/@r9d`.

#### Pass-count bp — simpler, but scope-fragile

```
bu nudge_tests!world_step 0n106
```

Stops on the 106th call from process start. **Fragile in a test binary:** the counter spans every earlier test's `world_step` calls, so the 106th hit may be in an unrelated test. Always verify at the stop (`?? ((WorldInternal*)@rcx)->frame` and a look at `kp`) before trusting it. Scope by first breaking on a bench-specific entry function when the test binary runs many tests before yours.

### Phase 4 — Step-through with a hypothesis

Never step aimlessly. Before every `p`/`t`, write down the hypothesis:

> "By the time we exit `sat_query_edges`, `best.index1` should be −1 (no edge-edge contact) because this is a flat triangle hull. If it's ≥ 0, that's the bug."

Then step, inspect, compare.

#### Inspection commands (cdb / win-debugger)

| Goal | cdb command | win-debugger tool |
|------|-------------|-------------------|
| List locals | `dv` | `locals` |
| Dump struct | `dt Cylinder @@c++(&cyl)` | `eval &cyl` + struct type |
| Evaluate C expression | `?? body_state[85].position.y` | `eval` |
| Stack trace | `k` / `kn` / `kv` | `stack` |
| Step over one line | `p` | `step_over` |
| Step into | `t` | `step_into` |
| Continue | `g` | `go` |

#### Things to watch for

| Symptom | What to check |
|---------|--------------|
| Body launched upward with huge impulse | `lambda_n` on each contact — anything >0.5 for a resting body is suspect |
| Body penetrates but no correction | Number of manifolds emitted for this body; contact `normal` direction |
| Wrong contact direction | `contact_n` vs `tri_n`; sign of `dot(cyl.center - tri_center, tri_n)` |
| Stale warm lambda | `im->warm != NULL` for mesh manifolds (should be NULL if mesh warm-start is off) |
| SAT picks wrong axis | `get np_debug` in the viewer → compare `face_a_sep`, `face_b_sep`, `edge_sep` + `winning_axis`. No cdb needed. |
| NaN appearing | Any computed value is `1.#QNAN` — trace backward to the divide/sqrt |
| Duplicate contacts at same point | Multiple manifolds emitting the same (point, normal) — collapse into 1 effective constraint |
| Feature ID mismatch | `feature_id` across frames: stable for resting, flickering for unstable contact |

**The bug lives at the moment expected diverges from actual.**

When you find divergence:
1. Note the exact file:line where the wrong value is produced.
2. Back up the stack (`k`) to find WHY that code path was chosen.
3. State the hypothesis in one sentence: *"on line X, variable Y is A but should be B because Z."*
4. Design ONE mechanical change that fixes Z without touching anything else.

### Phase 5 — Fix + varied-input verify

Minimal fix. 1–5 lines typical. If you're writing 50+ lines, you're fixing a symptom.

**Verify:**
1. Original repro is fixed.
2. 10+ varied configurations still pass (different body counts, shape combinations, mesh topologies, initial conditions).
3. Unit tests all pass.
4. A few related bench runs to catch collateral regressions.

Commit message includes: bug description, root cause, fix description, before/after numbers.

## Technique that paid off: impulse ledger vs. velocity delta

When you think the solver is the suspect (rare, but it does happen at the
SIMD/scatter/batching layer rather than the math), the fastest signal is to
**compare what the solver says it applied against what the body actually
experienced**. Two prints, side by side:

1. At end of each PGS iter: accumulated `lambda_n` for every contact on the
   suspect body.
2. At the same moment: the body's velocity *in body_hot* (not in the SIMD
   batch register).

If the ledger and the velocity agree — PGS did what it thought it did, so
move upstream to narrowphase / warm cache / broadphase as usual.

If they diverge — **the computed impulse is being written somewhere that
gets discarded.** In our case, `lambda_n` climbed linearly to 219 N·s over
8 iterations (substep 0) while `body_hot[30].velocity.y` changed by less
than 0.5. That single number pair — "solver is confident it pushed 219,
body says 0.5" — pointed the investigation straight at the writeback path.
One look at the SIMD scatter (`SCATTER_V3` with duplicate indices) and the
coloring fallback that produced those duplicates, and the bug was named.

**Why this works:** the solver's own books are the cheapest oracle you
have. It's keeping them correctly, one conservation-of-momentum calculation
at a time, so they reflect what SHOULD have happened. If reality disagrees,
the gap between ledger and effect is a measurable physical quantity
(impulse in N·s) that localizes to a specific substep and body. You don't
have to guess which phase is wrong — the conservation law does it for you.

**When to reach for this:** the body's post-step state doesn't match any
sane application of the inputs you fed it. Huge lambda + small velocity
change → scatter / batching / aliasing. Small lambda + large velocity
change → something OUTSIDE the solver is writing velocity (gravity-bypass,
a direct velocity set in user code, a corrupted warm cache entry
read as velocity).

**Cheap to set up:** one `fprintf` in the scatter site, one in the
integrator — both gated by a single global flag so they're off for the
rest of the run. Pair with `g_b30_trace = (wf >= 55 && wf <= 58)` to scope
it to the frames of interest. Total instrumentation: ~20 lines.

## Cognitive bias guards

| Bias | Antidote |
|------|---------|
| Confirmation bias | Actively search for evidence AGAINST your hypothesis |
| Anchoring on first theory | Have at least 2 alternative hypotheses before stepping |
| Sunk cost | 3 failed experiments in same file → switch files |
| Availability heuristic | Familiar pattern ≠ actual cause |
| Band-aid urge | When tempted to tweak a solver parameter — STOP. Narrowphase first |

## Anti-patterns

| Don't | Why |
|-------|-----|
| Tune solver params (`beta`, iters, hertz) as the fix | Almost never the root cause; masks the real bug |
| Increase substeps to fix tunneling | Temporal-aliasing band-aid; doubles perf cost; bug returns at higher velocities |
| Add speculative contacts | Creates ghost contacts with negative penetration; new bugs |
| Declare victory on one test | Varied-input verification is non-negotiable |
| Step through without a hypothesis | Random stepping wastes hours — you can't recognize "wrong" without "expected" |
| Break AT the bad frame | State already corrupted. Break ONE frame before. |
| Fix the symptom in one place, ignore the pattern | If `(x - eps)` is wrong in one narrowphase function, check all the others |
| Compare across nondeterministic runs | Same binary + same script must give identical numbers |
| Add 10 printfs and rebuild 10 times | Use one `DBG_BREAK` + viewer. All inspection costs zero rebuilds. |
| Leave large printf diagnostic blocks in committed code | Pause points (`DBG_BREAK`) are zero-cost when `--debug` isn't passed; printfs bloat output and rot |
| Accept "lambda sum is huge but velocity barely changes" as normal | Non-conservation of momentum is never a rounding error — it's a writeback path discarding impulses (SIMD scatter race, forgotten scatter on early-out, wrong body index). See "impulse ledger vs velocity delta" above. |

## cdb quick reference

### Session setup (always)

```
.lines -e                       Enable source-line info — DEFAULT IS OFF
l+t                             Step by source line (not asm)
l+s                             Display source at each prompt
.symfix; .reload                Fix up symbol path + reload (harmless if local)
```

### Breakpoints

```
bu <module>!<func>              Pending / unresolved — default choice
bp <module>!<func>              Resolved now (fails if module not loaded)
bu <func> 0n<N>                 Pass count: stop on the Nth hit
bu <func> "j (@@c++(<expr>)) '';'gc'"   Conditional; C++ expression mode
bl / bc / bd / be               List / clear / disable / enable
```

### Running

```
g                               Go until breakpoint / exit
gc                              Go from a conditional bp (used inside `j`)
gu                              Go until current function returns
p [n]                           Step over one source line (n = count)
t                               Step into
pc                              Step to next call
```

### Inspection

```
?? <expr>                       Evaluate any C/C++ expression — default tool
?? <ptr>                        Auto-expand struct: offsets + types + values
kp                              Stack trace with params + file:line (need `.lines -e`)
k / kn / kv                     Variants: plain, numbered, FPO info
dv /V /t                        Locals: addresses and types
dt <Type> @@c++(&var)           Dump struct by explicit type (fallback for ??)
r                               Registers
.frame /c <N>                   Switch stack frame and register context
x <module>!<pattern>            Enumerate symbols (wildcards, sometimes empty)
```

### Verified recipes

**Stop at frame N, then snapshot state** (the "one frame before" workflow):

```bash
cdb -c ".lines -e; l+t; l+s;\
 bu nudge_tests!world_step \"j (@@c++(@rcx != 0 && ((WorldInternal*)@rcx)->frame == 0n105)) '.echo AT_BP';'gc'\";\
 g;\
 ?? ((WorldInternal*)@rcx)->frame;\
 ?? ((WorldInternal*)@rcx)->sub_steps;\
 kp;\
 q" \
 ./build/Debug/nudge_tests.exe 2>&1 | sed -n '/AT_BP/,/quit:/p'
```

Wrap the full `-c` value in double quotes — the shell eats unquoted `;`. Escape inner quotes. Always add a `.echo MARKER` so you can `sed -n '/MARKER/,/quit:/p'` past ~200 lines of module-load noise.

**Stop on a specific body index in narrowphase:**

```
bu nudge_tests!narrowphase_pair "j (@@c++(@edx == 0n45 || @r8d == 0n45)) '';'gc'"
```

**Once stopped inside a function, dump all contacts:**

```
?? m->count
dt nudge_tests!Contact @@c++(&m->contacts[0]) -a4
```

### Top traps — verified

| Symptom | Cause | Fix |
|---------|-------|-----|
| `p` takes 15 hits to cross one C line | Line info disabled | `.lines -e` at session start |
| Locals show `0xCCCCCCCC` / garbage | Stopped at function entry, pre-assignments | Either `p` past prologue, or read via arg register: `?? ((WorldInternal*)@rcx)->field` |
| Conditional bp stops every time | MASM expression error → treated as true | Verify predicate with `?? <expr>` first, then use `@@c++(...)` |
| `Numeric expression missing from '& ...'` | MASM has no `&&`/`||` | Switch the whole predicate to `@@c++(...)` |
| `Memory access error at 'w->frame'` | Dereffed still-uninit local `w` | Use `@rcx` + cast at entry, or `p` past the assignment |
| Conditional bp fires on the very first hit before any frame count matters | cdb probed the predicate with a dummy/NULL context | Guard with `@rcx != 0` as the first conjunct |
| Silent bp miss | `bp` used before module loaded | Use `bu` |
| Pass-count bp stops in the wrong test | Counter is global across tests in the binary | Add a second bp on the bench entry function and count from there |
| `*** WARNING: Unable to verify checksum` | No Authenticode sig on the exe | Ignore |
| `x nudge_tests!*pattern*` returns empty | Not every name matches every wildcard form | Try a more specific substring or drop the outer `*` |
| `-c "cmd1; cmd2"` only runs `cmd1` | Shell ate the `;` | Always double-quote the `-c` arg |

### Why `??` beats `dt`

`?? w` on a pointer auto-expands the full type — offsets, types, and dereferenced values, recursively for nested pointers — in one command. `dt` requires the explicit module+type path and won't chase values through pointers without extra syntax. Reach for `dt` only when you have an address + type but no expression path (e.g. an allocation dump).

### DBG_BREAK + cdb attach (production workflow)

`DBG_BREAK(name, world)` in `debug_server.c:999` publishes globals and spin-sleeps in `dbg_break_wait_impl` until a viewer sends `continue`. That spin loop is a perfect cdb attach point:

1. Run `./build/Debug/nudge_tests.exe --debug --break=<pattern>`.
2. When it pauses, the stderr line prints `[dbg] break "<name>" at <file>:<line> (world=0x...)`.
3. Attach: `cdb -p <pid>` (or `cdb -pn nudge_tests.exe`).
4. `~*k` → find the thread parked in `dbg_break_wait_impl`, `~<n>s` to switch, `.frame 2` to step out into the caller.
5. Inspect anything with `??`, `dv`, `kp` — cdb sees locals the viewer's RPM cannot.
6. To resume without the viewer: `ed nudge_tests!g_dbg_break_resume 1; g`.

This beats either tool alone: cdb gives locals + stack; the viewer (on its own socket) gives reflected world state. Use both at once.

## Companion modules

- `autoresearch:debug` — breadth-first hypothesis testing when you don't know WHERE the bug is.
- `remote-debug` — TCP viewer + driver.
  - **App mode** (`nudge.exe`): mouse-interaction-triggered bugs, scene repros, `playrecording`.
  - **Test mode** (`nudge_tests.exe --debug --break=<pattern>`): pause points (`DBG_BREAK(name, world)`) inside unit tests. Use this instead of printf whenever you would otherwise add exploratory prints.
- `win-debugger` MCP — tool-level cdb access. Prefer this over shelling out to cdb.exe.

Decision flow:
1. Still forming hypotheses / don't know the file? → `autoresearch:debug`.
2. Know the scene / test and need to *inspect* state at a specific moment? → `remote-debug` test-mode pause point (default).
3. Need to *watch code execute line-by-line* to understand control flow or mutate locals? → `physics-debug` + `win-debugger` / cdb.
4. Need a time-series trend across many frames, or diff-able output for CI? → printf/CSV.

`physics-debug` (this skill) applies once you have a deterministic repro and
need to WATCH the bad frame execute. It consumes the output of the pause-point
inspection (or printf trace) to decide exactly where to break in cdb.
