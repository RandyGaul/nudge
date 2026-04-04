# Stability Optimization Log

Tracking benchmark results for box stack stability improvements.
Benchmark: `bench_box_stack(N)` -- N unit boxes stacked, 1000 frames at 60Hz, sleep disabled.
Metrics: cumulative motion (sum of per-frame position deltas), max drift from ideal rest positions.

## Baseline (before optimizations)

Solver: 10 velocity iters, 4 position iters (NGS), Baumgarte 0.2, slop 0.005

| Height | Motion    | Max Drift |
|--------|-----------|-----------|
| 5      | 2.944401  | 0.0336    |
| 10     | 78.549147 | 14.3540   |
| 20     | 255.383398| 23.9731   |

Stack of 5 settles well. Stack of 10+ diverges badly -- boxes slide/fall apart.

## Phase 2a: Runtime iteration counts

No behavioral change (same defaults 10/4). Infrastructure for later tuning.

## Phase 2b: Soft contact constraints

Box2D b2MakeSoft spring-damper model. hertz=30, damping_ratio=10, max_push_vel=3.0.
Static body contacts use 2x hertz. NGS position correction skipped.

| Height | Motion    | Max Drift | vs Baseline |
|--------|-----------|-----------|-------------|
| 5      | 2.133138  | 0.0287    | -28% motion, -15% drift |
| 10     | 64.113295 | 11.6528   | -18% motion, -19% drift |
| 20     | 276.434069| 20.7889   | +8% motion, -13% drift  |

Modest improvement alone. Sub-stepping expected to unlock the full benefit.

## Phase 2c: Sub-stepping

Split world_step into N sub-steps at dt/N. Full collision + solve pipeline per sub-step.
Warm cache aging once per frame. Sleep evaluation once per frame.

Default sub_steps=4, hertz=30, damping=10:

| Height | Motion    | Max Drift | vs Baseline |
|--------|-----------|-----------|-------------|
| 5      | 0.352612  | 0.0258    | -88% motion, -23% drift |
| 10     | 1.745535  | 0.0536    | -98% motion, -99.6% drift |
| 20     | 9.889315  | 0.1396    | -96% motion, -99.4% drift |

Sub-stepping is the dominant stabilizer.

## Phase 2d: Angular NGS

Skipped -- NGS is disabled when soft contacts are active (the new default).

## Phase 2e: Parameter tuning

### Hertz sweep (sub=4, vel=10, height=20)
| Hertz | Motion   | Max Drift |
|-------|----------|-----------|
| 10    | 281.6    | 25.9      |
| 20    | 6.06     | 0.185     |
| 30    | 9.89     | 0.140     |
| 60    | 9.53     | 0.138     |
| 120   | 313.0    | 35.9      |

Sweet spot: hertz=20. Below 10 or above 120 causes instability.

### Damping sweep (sub=4, vel=10, hertz=20, height=20)
| Damping | Motion   | Max Drift |
|---------|----------|-----------|
| 1.0     | 4.91     | 0.152     |
| 3.0     | 3.29     | 0.143     |
| 5.0     | 5.00     | 0.165     |
| 10.0    | 6.06     | 0.185     |
| 15.0    | 317.5    | 31.8      |

Sweet spot: damping=3. Range 1-5 all work well.

### Final tuned defaults

`sub_steps=4, velocity_iters=10, contact_hertz=20, contact_damping_ratio=3, max_push_velocity=3.0`

| Height | Motion    | Max Drift | vs Baseline |
|--------|-----------|-----------|-------------|
| 5      | 0.134847  | 0.0276    | -95% motion, -18% drift |
| 10     | 0.563865  | 0.0577    | -99.3% motion, -99.6% drift |
| 20     | 3.286324  | 0.1427    | -98.7% motion, -99.4% drift |

### Stable height limits by config
| Config         | Stable Height (drift < 0.5) |
|----------------|-----------------------------|
| sub=4, vel=10  | ~22                         |
| sub=8, vel=10  | ~28                         |
| sub=16, vel=10 | ~30+                        |

### Scaling bottleneck

Stacks above ~25 boxes topple laterally due to PGS solver convergence limits.
Graph coloring groups contacts into 2 colors for a linear chain; force propagation
through N contacts requires ~N/2 velocity iterations per sub-step. Warm starting
helps but cannot fully compensate for deep chains.

Reaching 1000 boxes would require either:
- A block/direct solver for constraint chains
- Constraint ordering (bottom-to-top) instead of graph coloring
- A fundamentally different solver (e.g., XPBD, TGS)

These are structural changes beyond parameter tuning.

## Local sub-stepping (current)

Replaced full-pipeline sub-stepping with local sub-stepping: collision detection runs
once per frame, sub-steps only re-solve velocity constraints and integrate positions.
Much cheaper (no repeated broadphase/narrowphase).

Rigid joint Baumgarte bias applied on first sub-step only (error is stale after that).
Soft constraint feedback (-softness * lambda_n) handles staleness naturally.

Tuning: higher hertz works better with local sub-stepping (stiffer contacts converge
faster when reusing stale contact data). hz=60 outperforms hz=20 here.

### Defaults

`sub_steps=4, velocity_iters=10, contact_hertz=60, contact_damping_ratio=3, max_push_velocity=3.0`

| Height | Motion    | Max Drift | vs Baseline |
|--------|-----------|-----------|-------------|
| 5      | 0.397561  | 0.0270    | -87% motion, -20% drift |
| 10     | 1.769827  | 0.0534    | -98% motion, -99.6% drift |
| 20     | 4.711212  | 0.1139    | -98% motion, -99.5% drift |
| 25     | 12.656948 | 0.1767    | stable      |

Stable height: ~25 boxes (sub=4). More sub-steps DON'T help here (stale contacts
accumulate error). With hz=90 damp=1, stack-30 reaches drift=0.14 but at cost of
higher motion.
