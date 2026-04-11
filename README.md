# Nudge
Lightweight 3D physics engine for small games or prototyping

## Preview
```c
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

## Features/Vision
- Pure C API -- highly portable and trivially bindable to other languages
- Really small API, nearly impossible to mess up
- Box/Sphere/Hull/Capsule/Mesh shape types
- Variety of common constraint/joints supported
- Unity build

## What's done?
Still WIP, but plenty of the core is usable.
- Rigid bodies stacking, colliding, and resting against each other
- Type-safe handle API -- stale `Body`/`Joint` handles are inert instead of crashing, foot-gun proof by design
- Shape types: sphere, capsule, box, convex hull, cylinder
- Build convex hulls from arbitrary point clouds
- Joints: ball-socket, distance, hinge, fixed, prismatic
- Joint limits, motors, and optional soft springs on any joint
- Direct solver for equality joints -- ropes, chains, and springs stay stable without jitter or stretch
- Friction, restitution, linear/angular damping per body
- Competitive performance via SIMD and multithreading
- Sleeping so big piles of resting bodies stop burning CPU
- Raycasts and AABB queries
- Demo app with a handful of scenes (pyramid, stacks, friction, mass ratio, ...)
- Side-by-side testbed that runs the same scene in nudge, Bepu, and Jolt for comparisons

## Ideas/Todos
- Triangle mesh shape (static/concave geometry)
- Swept physics step and continuous collision detection
- Single-file-header distribution
- Serialize/inspect tool (maybe remote visualizer??)
- Experiments! Different solvers, new constraint designs/formulations, new broadphase exploration, or CCD techniques
- Rollback

## Building
```
cmake -B build
cmake --build build
```
C23 compiler required.

## License

Nah. Public domain.
