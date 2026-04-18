---
name: remote-debug
description: Autonomous physics debugging against nudge.exe or nudge_tests.exe. Find a deterministic repro, isolate, measure, fix, and verify with varied inputs. Invokes the remote-debug skill.
argument-hint: "<bug description> [--iterations N]"
---

# /remote-debug

Invoke the `remote-debug` skill to debug the physics bug described in $ARGUMENTS.

## Parse arguments

- Bug description from $ARGUMENTS
- `--iterations N` -- max debug cycles (default 10)

## Rules

1. Read `.claude/skills/remote-debug/SKILL.md` first and follow its phases in order.
2. **Phase 0 is blocking.** Do NOT investigate until you have a deterministic repro (3 identical runs, fresh engine each time).
3. If the user has a mouse recording, use `playrecording` -- do NOT substitute `push`.
4. One variable per experiment. Revert every diagnostic change before the next one.
5. Prefer `DBG_BREAK` + viewer (test mode) or scripted driver commands (app mode) over exploratory printfs.
6. Fix the measured root cause in narrowphase / broadphase first. Do not tune solver parameters as the fix.
7. Verify with 15+ varied inputs. Failure to find any repro after the fix = confirmed.

## Don't

- Start investigating before a deterministic repro exists.
- Change multiple variables between comparisons.
- Declare victory on the one repro that triggered the bug.
- Tune `beta`, substeps, iterations, hertz, or eviction thresholds as "the fix".
- Spend >3 iterations in a file the evidence doesn't point at.

Now execute. Begin with Phase 0.
