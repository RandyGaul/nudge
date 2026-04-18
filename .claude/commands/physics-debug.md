---
name: physics-debug
description: Debug a physics bug via deterministic repro + cdb step-through, one frame before the anomaly. Invokes the physics-debug skill.
argument-hint: "<bug description>"
---

# /physics-debug

Invoke the `physics-debug` skill to debug the physics bug described in $ARGUMENTS.

## Rules

1. Read `.claude/skills/physics-debug/SKILL.md` first.
2. Follow its 5 phases in order. Do NOT skip Phase 1 (deterministic repro) — that is the blocking prerequisite.
3. In Phase 2, prefer a `DBG_BREAK(name, world)` pause point + `remote-debug` viewer over speculative printfs. Printf only when you need a time series or diff-able output.
4. Prefer the `win-debugger` MCP tools over shelling out to `cdb.exe` once you're in Phase 3.
5. Break ONE frame before the anomaly, not at it.
6. Fix the root cause, not the symptom. 1–5 line fixes. No solver-parameter tuning.
7. Verify with at least 10 varied repros, not just the one that triggered the bug.

## Don't

- Start investigating before you have a deterministic repro that fires 3 identical runs.
- Add 10 exploratory printfs and rebuild 10 times — put a `DBG_BREAK` in and use the viewer.
- Turn solver knobs (beta, substeps, iterations, hertz) as "the fix".
- Declare done on one passing test.
- Step through without a written hypothesis.

Now execute. Begin with Phase 1.
