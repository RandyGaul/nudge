---
name: remote-debug
description: Autonomous physics debugging against nudge.exe (app/scene-based) and nudge_tests.exe (unit-test pause points). Scientific method with deterministic repros, mechanical metrics, varied-input verification, and bias guards. Activate when the user asks to debug a physics bug end-to-end and wants the agent to drive from repro to verified fix.
version: 1.0.0
---

# Remote Debug -- Autonomous Physics Debugging

End-to-end loop: find a deterministic repro, isolate one variable at a time, measure against expectations, propose a minimal fix, and re-verify across many varied inputs. The agent drives; the human watches.

## When to activate

- User points at a bug and asks Claude to find + fix it.
- User says things like "debug this", "why is X happening", "find the repro", or invokes `/remote-debug <bug>`.
- A unit test is failing and the user wants pause-point inspection instead of printf spelunking.

## When NOT to activate

- The user wants a feature, refactor, or design discussion.
- The user is asking how an algorithm works — that is research, not debugging.
- You already have a verified repro and root cause and just need to implement the fix.

## Two modes

The viewer (`nudge_viewer.exe`) talks to either target over TCP and inspects memory via `ReadProcessMemory`. Type tables are sent by the engine on connect, so all command syntax is identical -- the difference is the simulation host.

| Mode | Target | When to use |
|------|--------|-------------|
| **App mode** | `nudge.exe` | Scene-based bugs: scene loads, interactive push/drag, rendering-tied repros, mouse-recording playback. |
| **Test mode** | `nudge_tests.exe --debug --break=<pattern>` | Unit-test bugs: a specific test fails and you want to pause at a known point and inspect state -- an alternative to printf debugging. |

Both modes support the same inspection commands (`summary`, `get`, `table`, `contacts`, `warm`, `filter`, `snap`, `diff`, `stable`, `raw`). Test mode additionally supports pause points (`DBG_BREAK(name, world)`) and the control commands `continue` / `where` / `break-filter`.

## Core principles

1. **Find the repro FIRST.** Can't debug what you can't reproduce. Don't investigate until you have a deterministic repro that shows the bug.
2. **Use the user's actual interaction.** If the user triggers it with the mouse, use `playrecording` to replay their mouse recording. Don't substitute with `push` -- it's a different interaction pattern.
3. **Determinism is required for comparison.** Same script + same binary = same result. Fresh engine per run. Non-deterministic comparisons are INVALID.
4. **ONE variable per experiment.** Never change two things between comparisons.
5. **The fix is usually small.** If you're writing 50+ lines, you're probably fixing a symptom. The actual fix for this engine's capsule bug was 2 lines.
6. **Band-aids create new bugs.** Speculative contacts, warm cache thresholds, substep increases -- these mask the real bug and introduce new problems.
7. **Listen to the user.** When they say "narrowphase", look at the narrowphase. When they say "it's not the solver", stop modifying the solver.
8. **Verify with VARIED inputs.** One passing test is not proof. Run 15+ varied scenarios. Failure to find ANY repro after the fix = confirmed.

## Phase 0: Find the repro

**This is the MOST IMPORTANT phase. Do NOT skip or rush it.**

**Step 1: Try the user's exact interaction.**
- If the user has a mouse recording: `playrecording` command
- If the user describes an interaction: replicate with `drag`/`push`
- If the user says "knock capsules over": knock them over at various angles and forces

**Step 2: Vary inputs until you find a failure.**
- Try different push forces: gentle (1), moderate (3), hard (5)
- Try different angles: sideways, diagonal, upward+sideways
- Try different application points: top, middle, offset
- Try different scenes: single body, multiple bodies, mixed shapes
- The bug may only appear at specific angular velocities or tilt angles

**Step 3: Confirm the repro is deterministic.**
```bash
# Run the SAME commands 3 times with fresh engines
for run in 1 2 3; do
  ./build/Debug/nudge.exe &
  sleep 2
  echo -e "scene X\npause\nstep N\npush ...\nstep M\nget body_hot/K\nexit" | timeout T ./build/Debug/nudge_viewer.exe localhost 9999
  taskkill //F //IM nudge.exe
  sleep 1
done
# All 3 outputs must be IDENTICAL
```

**Step 4: If you can't find a repro with `push`, use `playrecording`.**
The user's mouse constraint creates sustained force over many frames -- fundamentally different from an instant impulse. Add commands to replay user recordings rather than trying to approximate the interaction.

**A repro that settles after a few seconds is NOT a repro of "persistent erratic bouncing."** Keep varying inputs until you find one that shows the EXACT symptom the user describes.

## Phase 1: Isolate -- one variable at a time

**Disable everything, then add back one at a time:**
- N^2 broadphase (eliminates BVH)
- Disable sleep, warm start, incremental NP
- Single body + ground (eliminates multi-body)

**For EACH toggle, run the SAME deterministic repro and compare.**
If the result changes when toggling X, X is involved.
If the result is identical, X is not involved -- move on.

## Phase 2: Measure

Record at the frame of anomaly:
- pen depth, lambda, velocity, feature_id, contact count
- Expected pen = velocity * dt / substeps
- Expected lambda for resting ≈ 0.08
- If measured >> expected, the discrepancy IS the bug

**Track the full contact cycle:** contact frame -> gap frames -> re-contact.
Healthy contacts persist many frames. Single-frame contacts are the red flag.

## Phase 3: Hypothesize -- in the right file

**Priority:**
1. Narrowphase early-exit conditions (hard boundaries without slop margin)
2. Broadphase pair filtering (tight AABB check discarding near-contact pairs)
3. Contact filtering between NP and solver
4. Solver convergence (only after 1-3 are ruled out)

**Key lesson learned:** The capsule bounce bug was a 2-line fix in the narrowphase early exit -- `r.distance > a.radius` needed `+ LINEAR_SLOP`. Not in the solver, not in warm start, not in the broadphase.

**If 3 experiments in the same file don't help, you're in the wrong file.**

## Phase 4: Test -- one atomic experiment

- Add ONE printf, rebuild, run repro, read output
- Revert the printf after reading
- ONE experiment per iteration

## Phase 5: Fix -- minimal

- Fix the MEASURED root cause
- **DO NOT** tweak solver parameters, thresholds, substep counts, or margins
- **DO NOT** add speculative contacts or change eviction thresholds
- The fix should be in the narrowphase or broadphase code path, not the solver

## Phase 6: Verify -- varied inputs (critical)

**Running the original repro is NOT ENOUGH. You must search for a new repro with varied inputs and FAIL to find one.**

Run at minimum 15 tests covering:
```
# Push variations (6+ tests)
push (1,0,0), (2,0,0), (3,0,0), (0,0,2), (1,0,1), (2,1,1)

# Push with torque (3+ tests)
push with offset application points

# Drag interactions (3+ tests)
gentle drag, moderate drag, diagonal drag

# Multi-body (2+ tests)
bowling pins, capsule-capsule collision

# Regression (2+ tests)
sphere, box -- must still work

# Test suite
./build/Debug/nudge_tests.exe -- ALL tests pass
```

**For each test: check |vy| < 0.01 after 10 seconds (600 frames).**

**Only declare FIXED when:** Zero failures across all varied tests + all unit tests pass.

## Progress logging

Every 3 iterations:
```
=== Debug Progress (iteration N) ===
Hypotheses: X tested (Y confirmed, Z disproven)
Repro status: [found / searching / confirmed deterministic]
Variables ruled out: [list]
Current file under investigation: [path]
```

## Key lessons (from capsule bounce debugging)

| Lesson | Detail |
|--------|--------|
| The fix was 2 lines | `r.distance > a.radius` -> `r.distance > a.radius + LINEAR_SLOP` in two functions |
| Band-aids created new bugs | Speculative contacts caused negative penetrations and energy injection |
| The user's recording was the real repro | `push` couldn't reproduce the exact mouse-drag dynamics |
| 10+ iterations were wasted in the solver | The bug was in the narrowphase all along |
| N^2 vs BVH comparison was not deterministic | Different frame timing invalidated the comparison |
| "Settled after 10s" is not "no bug" | The user sees erratic bouncing in the FIRST few seconds |
| Contact toggling at tangent boundary | GJK distance oscillating around `radius` caused 1-frame contacts |

## Anti-patterns

| Don't | Why |
|-------|-----|
| Substitute `push` for mouse interaction | Different dynamics -- use `playrecording` |
| Change multiple variables per experiment | Can't attribute results |
| Declare victory on one test | Must fail to find repro across 15+ varied inputs |
| Spend >3 iterations in wrong file | Switch to the file the user pointed at |
| Fix symptoms (big lambda) | Find WHY the contact is intermittent |
| Add speculative contacts | Creates ghost contacts with negative penetration |
| Increase substeps as a fix | Masks temporal aliasing, doubles perf cost |
| Compare non-deterministic runs | Use fresh engine per run, same script |

## Commands

Run `help` in the viewer for all viewer commands (RPM inspection).
Run `help` as a driver command for all engine commands (TCP control).

Key commands for debugging:
- `summary` -- scene overview
- `contacts <body>` -- solver contacts for a body
- `stable [threshold] [frames]` -- mechanical pass/fail stability check
- `push <idx> <f> [r]` -- impulse at body (app mode only)
- `playrecording` -- replay user's mouse recording (app mode only)
- `npviz` + `npvset <A|B> <kind> <pos> <rot> <params>` -- load shapes into NP Viz for visual narrowphase debugging (app mode only)
- `continue` / `c` -- resume from a pause point (test mode)
- `where` -- show current pause point name + file:line (test mode)
- `break-filter [pattern]` -- view or change the break filter at runtime (test mode)
- `np-debug <0|1> [body_a] [body_b]` -- toggle narrowphase snapshot capture + filter (either body index, `-1` = any). After enabling, step one frame (`step 1` in app mode, or `continue` past a DBG_BREAK in test mode), then `get np_debug` dumps the SAT intermediates and emitted contacts for the first matching hull-hull pair. Force `--threads 1` when capturing -- the slot is single-entry and races under parallel NP.

## Test mode -- pause points instead of printf

Use this when a unit test is misbehaving and you need to inspect live engine state mid-test without recompiling with printfs.

### Adding a pause point

Place `DBG_BREAK("name", world)` anywhere in a test, inside the loop or right before an assertion. The macro is defined by `src/debug_server.c` (included into the test binary) and is a zero-cost predicate check unless the binary is started with `--debug`. Example:

```c
for (int f = 0; f < 60; f++) {
    world_step(w, 1.0f / 60.0f);
    if (f == 30) DBG_BREAK("my_test:frame30", w);   // pauses here when armed
}
```

Name pattern: pick something greppable and scope-prefixed. `my_test:frame30` beats `break1`. The viewer glob-matches against this.

### Running the test under the debugger

```bash
# Start the test binary with --debug. --break arms one or more named points.
# Without --break, DBG_BREAK is a no-op. --break=* catches every pause point.
powershell Start-Process -FilePath 'build/Debug/nudge_tests.exe' \
    -ArgumentList '--break=my_test:*', '--unit' -PassThru

# Or run it directly in foreground and attach a viewer from another shell.
./build/Debug/nudge_tests.exe --debug --break=my_test:* --unit
```

When a matching `DBG_BREAK` is reached the test thread blocks. The binary prints `[dbg] break "<name>" at <file:line> (world=0x...)` to stderr and keeps servicing TCP until a viewer sends `continue`.

### Driving the viewer

```bash
# Interactive
./build/Debug/nudge_viewer.exe localhost 9999

# Scripted (pipe commands via stdin)
(echo "where"; echo "summary"; echo "get body_state"; echo "continue"; echo "quit") \
    | ./build/Debug/nudge_viewer.exe localhost 9999
```

Typical inspection loop at a break:
```
> where
OK break="my_test:frame30" at src/tests_weld_bridge_unit.c:83 world=0x...
> summary                    # frame, body count, contact count, islands
> get body_hot 1             # velocities / angular velocities of body 1
> get body_state 2           # position + rotation of body 2
> table body_hot             # spreadsheet view of all bodies
> contacts 1                 # solver contacts touching body 1
> continue                   # resume; fires again at the next matching break
```

### Flow between breaks

After `continue` the test resumes and runs until the next matching pause point (or the test ends). The viewer prompt returns immediately; send `info` or `where` to check whether you're at a new break, still running, or the binary has exited. Tighten the filter at runtime with `break-filter <pattern>` if you only care about later breaks in the same run.

### Disabling without editing

Leaving `DBG_BREAK(...)` calls in checked-in tests is fine -- they cost one predicate check when `--debug` isn't passed. You don't need to `#ifdef` them out. Skipping past a compiled break without editing the source: set the filter so the name doesn't match (`--break=somethingElse`) or omit `--break` entirely.

## Companion modules

- `physics-debug` -- deeper methodology reference (5-phase breakdown, break-one-frame-before rule, cdb step-through recipes, impulse-ledger-vs-velocity check, bias guards). Reach for it when the bug is inside a function body and reflection alone isn't enough.
- `win-debugger` MCP -- tool-level cdb access for function locals and callstacks.

## Now begin

Find the repro first. Use `playrecording` if the user has a recording. Vary inputs. One variable at a time. For unit-test bugs, prefer `DBG_BREAK` + viewer over adding printfs.
