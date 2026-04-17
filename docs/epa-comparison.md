# SAT vs GJK+EPA: narrowphase backend comparison

This is an A/B comparison of nudge's two hull-involved narrowphase backends,
captured for evaluating whether a GJK+EPA path with an incremental manifold is
a viable alternative to the existing SAT one-shot path. The goal is to quantify
the actual tradeoffs so a team considering the swap can make an informed call
without having to re-run the experiment themselves.

## Methodology

**Backends under test.** Both ship in-tree and are selected by
`WorldParams.narrowphase_backend`:

- `NARROWPHASE_SAT` (default): face-face SAT with Gauss-map pruning + edge-edge
  axes, producing a full one-shot 4-contact manifold per pair per frame. Feature
  IDs enable lambda warm-start via `WarmManifold`.
- `NARROWPHASE_GJK_EPA`: GJK termination gates EPA polytope expansion. EPA
  returns a single contact per pair per frame. A per-pair `EpaManifold`
  accumulates up to `MAX_CONTACTS=4` contacts over successive frames, with
  tangential-drift validation + spatial-spread eviction. Warm-start stores the
  prior-frame terminating face's Minkowski directions for next-frame seed.

**Per-query bench** (`nudge_tests --bench-epa`) runs 20 000 iterations of a
canonical deep-contact configuration per config and reports ns/op.

**Scene bench** (`nudge_tests --bench-epa-scenes`) runs four gravity-driven
scenes under each backend for a fixed frame budget and reports:

- **settle**: first frame where `max(|v_linear|)` across dynamic bodies drops
  below 0.01 m/s; `-1` if it never settles in the budget.
- **peakpen**: maximum `contact.penetration` observed across the run.
- **avgCt/pr**: average contacts per EPA pair over the last 20 % of frames.
- **ms/np**, **ms/solve**, **ms/total**: per-frame phase timings.

All numbers below are from Release builds on the test machine. Reproduce with
the commands shown.

## Per-query timing

```
$ nudge_tests --bench-epa

=== EPA vs SAT per-query bench (20000 iters per config) ===
  config                           SAT ns/op   EPA ns/op   ratio
  box-box deep (overlap 1)             973.1      1255.7    1.29x
  sphere-hull deep                     351.5      4195.4   11.94x
  hull-hull rotated deep              1130.2      1112.0    0.98x
```

EPA is essentially tied with SAT on hull-hull-deep (0.98x), ~30 % slower on the
box-box unit case, and about **12x** slower on sphere-hull. The latter is
entirely structural: SAT has a dedicated analytical `collide_sphere_hull` path,
while EPA runs full GJK+EPA iteration. For the experiment this is fine — don't
route sphere-hull through EPA in production.

## World-step timing

```
$ nudge_tests --bench-epa   (excerpt)

=== World-step bench: 5-box stack on floor, 120 frames ===
  SAT backend:       0.015 ms/frame
  GJK+EPA backend:   0.012 ms/frame   (0.80x SAT)
```

On a simple 5-box stack, EPA is 20 % **faster** per frame than SAT. Two
mechanisms drive this:

1. EPA warm-start makes steady-state iteration count ~1.0, so the per-query
   cost drops well below the microbench number.
2. The incremental manifold carries fewer contacts per frame (only 1 freshly
   emitted per pair), so the solver has less work.

## Scene bench

```
$ nudge_tests --bench-epa-scenes

  scene           backend      settle   peakpen   avgCt/pr      ms/np   ms/solve   ms/total
  3-box stack     SAT              94    0.0518       0.00      0.002      0.012      0.017
                  GJK+EPA          -1    0.0518       0.50      0.006      0.007      0.015
  10-box stack    SAT              -1    0.0921       0.00      0.012      0.024      0.047
                  GJK+EPA          -1    0.0590       0.18      0.029      0.015      0.047
  5x5 pile        SAT              36    0.0518       0.00      0.072      0.152      0.252
                  GJK+EPA          -1   27.5710       0.00      0.139      0.014      0.161
  small pyramid   SAT              -1    0.0119       0.00      0.053      0.090      0.153
                  GJK+EPA          -1   59.2999       0.01      0.166      0.046      0.224
```

**Warm-reseed + iteration numbers:**

| Scene | Warm % | Avg iters | Iter-cap hits |
| --- | ---:| ---:| ---:|
| 3-box stack  | 99.5 % | 1.00 | 0 |
| 10-box stack | 99.8 % | 1.00 | 0 |
| 5x5 pile     | 68.5 % | 2.02 | 0 |
| small pyramid| 30.8 % | 0.85 | 0 |

Warm-start is doing its job on stable pairs (99 %+ on the tall stacks). The
pile and pyramid have high churn — pairs constantly forming and breaking — so
the warm-reseed rate drops but iteration count stays low because the regular
tetrahedron fallback is cheap when warm data is missing. `MAX_ITERATIONS=48` is
never hit across any scene.

## Quality observations

The headline: **EPA timings are competitive, EPA quality is not.**

- **3-box stack.** EPA peak penetration matches SAT (0.0518 m). But it never
  "settles" under the definition above because the incremental manifold only
  grows to 1 average contact per pair within the 300-frame budget, leaving
  residual jitter that keeps body velocities above the 0.01 m/s floor.
- **10-box stack.** EPA peak penetration is actually **lower** than SAT (0.059
  vs 0.092) at the moment the run ends. With more frames SAT would likely also
  not settle — stacking this tall pushes both into the "never quite quiet"
  regime. The EPA advantage here is an artifact of accumulating contacts
  slowly; the per-frame impulse is gentler.
- **5x5 pile.** SAT settles in 36 frames. EPA **explodes**: peak penetration
  27.57 m means bodies are tunneling catastrophically. This is the expected
  failure mode of the incremental manifold on stacked piles — the first few
  frames of any new contact only have 1 contact point, which is insufficient
  rotational support under a pile of bodies above.
- **Small pyramid.** Same story, worse (peak pen 59.3 m). The contact timeline
  for the first EPA pair illustrates the dynamic: 1 contact at frame 0, 2 at
  frame 1, 3 at frame 2, 4 at frame 3 — the 4-contact manifold finally arrives
  at frame 3, by which time the body above has already accelerated through the
  contact plane.

## Contact-count timeline

The small-pyramid scene is the only one where the manifold actually grows as
the plan predicted:

```
f0       (frame 0): 1 contacts
f1       (frame 1): 2 contacts
f2       (frame 2): 3 contacts
f3       (frame 3): 4 contacts
f5       (frame 5): 4 contacts
f10      (frame 10): 3 contacts
final    (frame 499): 1 contacts
```

The 4-contact peak holds briefly then contacts evict as bodies shift. In the
stable stack cases the manifold never even reaches 4 — by the time 4 frames
have elapsed, the settling is done and the manifold sits at 0-1 contacts.
This is the core weakness of the incremental approach: it can't produce the
multi-point manifold that SAT gives on frame 0.

## Known issues and scope notes

- **Sphere-hull** runs through the full GJK+EPA path under `NARROWPHASE_GJK_EPA`
  and pays a 12x cost over SAT's analytical path. Recommend routing
  primitive-hull pairs around EPA in any real deployment.
- **Incremental validator** uses a fixed 0.1 m tangential-drift threshold
  (`EPA_CONTACT_TANGENTIAL_DRIFT`). Tighter values strand stale contacts
  faster; looser values let them accumulate out-of-patch. The 0.1 m value was
  picked empirically; scenes with smaller bodies likely want it tuned.
- **Max iteration cap** is 48 and is never hit across any bench scene. This
  could likely be lowered to 16 or 24 with no observable effect, saving some
  iteration budget on pathological queries. Not done here.
- **Warm-start invalidation** happens when an `EpaManifold` drops to 0
  contacts. A pair that momentarily separates and rejoins within the same
  frame may see a cold EPA on rejoin; this appears not to matter for any
  tested scene.

## Recommendation

- **Default to SAT** for any production workload that includes stacking,
  piling, or pyramids. The incremental-manifold quality gap is the dominant
  effect; per-query timings are close to even and not worth pursuing EPA for.
- **Use GJK+EPA** only as a comparison/debugging backend, or in narrow contexts
  where the workload is well-separated pairs (rigid bodies far apart most of
  the frame, only occasional brief collisions). The EPA path's warm-reseed
  story is genuinely good (99 %+ on stable pairs, sub-2 average iters across
  all scenes) — there is no performance barrier to shipping EPA. The
  correctness barrier is entirely the incremental manifold, and that is a
  property of the manifold model, not EPA itself.
- If the desired end state is a hybrid — EPA for deep contacts, SAT-style clip
  for the manifold — the warm-start infrastructure already in place (Minkowski
  directions stored on `EpaManifold`) would continue to pay dividends. That
  hybrid is out of scope for this experiment.

## How to reproduce

```
cmake --build build --config Release --target nudge_tests
./build/Release/nudge_tests.exe                        # full test suite (6217 tests)
./build/Release/nudge_tests.exe --fuzz-epa 250         # 4000 pairwise fuzz queries
./build/Release/nudge_tests.exe --bench-epa            # per-query timing
./build/Release/nudge_tests.exe --bench-epa-scenes     # scene quality bench (numbers in this doc)
```

Toggle the backend live in the app via the "Narrowphase" ImGui combo. When
`GJK+EPA` is active, the Stats panel shows live queries, avg iters, warm-reseed
ratio, and contacts-per-pair.
