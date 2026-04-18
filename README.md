nudge
=====

Small, handle-based 3D rigid body physics engine in C. Tiny API surface, drop
a couple of headers + one translation unit into your project and go.

```c
#include "nudge.h"

World world = create_world((WorldParams){ .gravity = {0, -9.81f, 0} });

Body floor = create_body(world, (BodyParams){ .mass = 0 }); // mass 0 = static
body_add_shape(world, floor, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = {10, 0.5, 10} });

Body ball = create_body(world, (BodyParams){ .mass = 1.0f, .position = {0, 5, 0} });
body_add_shape(world, ball, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 0.5f });

while (running) {
    world_step(world, 1.0f / 60.0f);
    v3 pos = body_get_position(world, ball);
}
destroy_world(world);
```


What you get
------------

### Shapes

Sphere, capsule, box, convex hull, cylinder, and static triangle mesh. Build
convex hulls from arbitrary point clouds with `quickhull`, or use the built-in
unit box. One body can hold several child shapes with local offsets and
rotations -- compound colliders come for free.

```c
// Compound: three shapes on one body.
Body robot = create_body(world, (BodyParams){ .mass = 3.0f, .position = {0, 2, 0} });
body_add_shape(world, robot, (ShapeParams){
    .type = SHAPE_BOX, .box.half_extents = {0.3f, 0.5f, 0.2f }
});
body_add_shape(world, robot, (ShapeParams){
    .type = SHAPE_CAPSULE,
    .local_pos = {0.4f, 0.2f, 0},
    .capsule.half_height = 0.3f, .capsule.radius = 0.1f
});
body_add_shape(world, robot, (ShapeParams){
    .type = SHAPE_SPHERE,
    .local_pos = {0, 0.9f, 0},
    .sphere.radius = 0.25f
});

// Build a convex hull from a point cloud.
v3 pts[8] = { ... };
Hull* h = quickhull(pts, 8);
Body rock = create_body(world, (BodyParams){ .mass = 2.0f, .position = {5, 3, 0} });
body_add_shape(world, rock, (ShapeParams){
    .type = SHAPE_HULL, .hull.hull = h, .hull.scale = {1, 1, 1}
});

// Static triangle mesh as level geometry.
TriMesh* terrain = trimesh_create(verts, vert_count, indices, tri_count);
Body ground = create_body(world, (BodyParams){ .mass = 0 });
body_add_shape(world, ground, (ShapeParams){ .type = SHAPE_MESH, .mesh.mesh = terrain });
```


### Joints

Ball-socket, distance, hinge, fixed, prismatic, angular motor, twist limit,
cone limit, and swing-twist (ragdoll). Limits and motors on hinge /
prismatic / distance. Optional soft-spring mode (frequency + damping ratio)
on any joint.

```c
// Hinge with angular limits and a motor.
Joint elbow = create_hinge(world, (HingeParams){
    .body_a = upper_arm, .body_b = forearm,
    .local_offset_a = {0, -0.5f, 0}, .local_offset_b = {0, 0.5f, 0},
    .local_axis_a = {1, 0, 0}, .local_axis_b = {1, 0, 0},
});
joint_set_hinge_limits(world, elbow, -1.2f, 2.5f);  // radians
joint_set_hinge_motor(world, elbow, 3.0f, 5.0f);    // speed, max impulse

// Soft distance (rope-spring) with frequency / damping ratio.
create_distance(world, (DistanceParams){
    .body_a = ceiling, .body_b = payload,
    .local_offset_a = {0, 0, 0}, .local_offset_b = {0, 0.5f, 0},
    .rest_length = 2.0f,
    .spring = { .frequency = 3.0f, .damping_ratio = 0.2f }
});
```


### Simulation

Rigid body stacking, friction, restitution, linear / angular damping. Sleep
so stacks and piles stop burning CPU once at rest. Rolling friction for
balls / cylinders that should lose spin on contact. Gyroscopic integration
for tops, cylinders on edge. Picks up free threads when available;
single-threaded works too.

```c
Body die = create_body(world, (BodyParams){
    .mass = 0.8f, .position = {0, 5, 0},
    .friction = 0.6f, .restitution = 0.15f,
    .rolling_friction = 0.04f,      // spin decays once rolling on ground
    .linear_damping = 0.0f, .angular_damping = 0.05f,
});
body_add_shape(world, die, (ShapeParams){ .type = SHAPE_BOX, .box.half_extents = {0.25f, 0.25f, 0.25f} });
```


### Contact reporting

`world_contact_summaries()` returns the current step's collisions sorted by
body pair, deduplicated, with one normal + patch point + radius + depth per
pair. For per-body callbacks, `body_set_contact_listener` fires once per
step with an array of contacts where `self` is always in the `a` slot --
the engine flips normals and swaps fields for you.

```c
void on_player_hit(Body self, const ContactSummary* pairs, int n, void* ud) {
    for (int i = 0; i < n; i++) {
        v3 n = pairs[i].normal;        // points away from self
        float depth = pairs[i].depth;
        uint8_t other_mat = pairs[i].material_b;
        if (depth > 0.1f && pairs[i].normal.y > 0.8f)
            play_landing_sound(other_mat);
    }
}

body_set_contact_listener(world, player, on_player_hit, &game);

// Or poll the sorted buffer after world_step:
int n = 0;
const ContactSummary* s = world_contact_summaries(world, &n);
for (int i = 0; i < n; i++) { /* ... */ }
```


### Materials

World-level palette of 256 entries (friction / restitution / user data).
Bodies carry a default material id, and triangle meshes can carry per-
triangle ids. Contact summaries report the resolved id for each side, so
you can pick audio clips, VFX, or whatever else from a single lookup.

```c
world_set_material(world, 1, (Material){ .friction = 0.9f, .restitution = 0.05f, .user_data = MAT_WOOD });
world_set_material(world, 2, (Material){ .friction = 0.2f, .restitution = 0.30f, .user_data = MAT_ICE });

body_set_material_id(world, crate, 1);

// Per-triangle materials on a trimesh (e.g. grass vs rock patches):
uint8_t tri_mats[tri_count] = { /* one id per triangle */ };
trimesh_set_material_ids(terrain, tri_mats);
```


### Queries

Raycasts and AABB queries today. Sensors are overlapping shapes that never
enter the physics solver -- useful for trigger volumes, damage zones, AI
perception.

```c
RayHit hit;
if (world_raycast(world, V3(0, 10, 0), V3(0, -1, 0), 50.0f, &hit))
    printf("hit body %llu at y=%f\n", hit.body.id, hit.point.y);

Body overlaps[64];
int n = world_query_aabb(world, V3(-5, 0, -5), V3(5, 3, 5), overlaps, 64);

Sensor zone = create_sensor(world, (SensorParams){ .position = {10, 1, 0} });
sensor_add_shape(world, zone, (ShapeParams){ .type = SHAPE_SPHERE, .sphere.radius = 3.0f });
Body bodies_inside[32];
int count = sensor_query(world, zone, bodies_inside, 32);
```


### Persistence

Snapshot save/load is versioned binary, DEFLATE-compressed. Hulls and
triangle meshes are referenced by name so your asset loader rebinds them on
load. Rewind is a ring buffer of deterministic snapshots you can jump back
to; delta-compressed, so resting piles cost almost nothing.

```c
// Save / load.
hull_set_name(h, "rock_medium");
trimesh_set_name(terrain, "level1_terrain");
world_register_hull(world, h);
world_register_mesh(world, terrain);
world_save_snapshot(world, "save.nudge");

World w2 = create_world(wp);
world_register_hull(w2, h);
world_register_mesh(w2, terrain);
world_load_snapshot_into(w2, "save.nudge");

// Rewind ring: capture every step, jump back.
world_rewind_init(world, (RewindParams){ .max_frames = 120, .auto_capture = 1 });
uint64_t checkpoint = world_rewind_capture(world);
// ... run some frames ...
world_rewind_to_frame(world, checkpoint);     // bit-identical to the capture
world_rewind_by_steps(world, 30);             // or: jump back 30 steps
```


### Debugging

A tiny TCP debug server is built into the engine and a remote viewer
(`tools/viewer.c`) attaches to a running game and walks its memory directly
-- no engine changes needed to inspect bodies, manifolds, contacts,
islands, BVH, or the warm cache. There's also a cross-engine testbed that
runs the same scene in nudge, Bepu, and Jolt side-by-side so you can diff
behavior or compare perf.


### Language / packaging

- Plain C API, easy to bind from other languages.
- Handle-based, foot-gun resistant: stale `Body` / `Joint` / `Sensor`
  handles return inert / zero rather than crashing.
- Header is `extern "C"` and self-contained; `#include "nudge.h"` from a
  C++ TU works with no special flags.
- Shared-library builds: define `NUDGE_BUILD_DLL` when building as a
  shared library and `NUDGE_USE_DLL` when linking against it. Or
  `#define NUDGE_API ...` before including to override entirely.


What's missing
--------------

- Continuous collision (CCD) / swept shape casts.
- Shape-cast and overlap-shape world queries (raycast and AABB query
  work today).
- Heightfield (use triangle mesh for now).
- Character controller, vehicles, soft bodies.


Building
--------

```sh
cmake -B build
cmake --build build
```

C23 compiler required. The engine is a unity build rooted at `src/main.c`;
tests live in `src/test_main.c` and build as `nudge_tests`.


License
-------

Public domain.
