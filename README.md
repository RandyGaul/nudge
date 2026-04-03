# Nudge
Lightweight 3D physics engine for small games or prototyping

## Preview
```c
World world = create_world(256); // max 256 bodies
Body floor = create_body(world, (BodyParams){ .mode = STATIC });
add_box(world, floor, (BoxParams){ .half_extents = {10, 0.5, 10} });

Body ball = create_body(world, (BodyParams){ .mass = 1.0, .position = {0, 5, 0} });
add_sphere(world, ball, (SphereParams){ .radius = 0.5 });

while (running) {
    world_step(world, 1.0f / 60.0f);
    // read back body positions for rendering
}
```

## Features/Vision
- Pure C for strong perf and cross-platformness
- Really small API, nearly impossible to mess up
- Box/Sphere/Hull/Capsule/Mesh shape types
- Variety of common constraint/joints supported
- Unity build (todo: package as single-file-header script)

## What's done?
This is in early WIP state.
- Basic world step
- Add rigid bodies + shapes
- Discrete collision detection (no sweep/time of impact/continuous)
- Basic joints (rod/ball + socket)
- Barely working sleep + broadphase

## Ideas/Todos
- Incremental sub-island sleeping
- Mesh collisions
- Full suite of common constraint/joint implementations
- Swept physics step and collision detection
- Serialize/inspect tool (maybe remote visualizer??)
- Experiments! Different solvers, new constraint designs/formulations, new broadphase exploration, or CCD techniques

## Building

```
cmake -B build
cmake --build build
Requires SDL3 (fetched automatically) for the demo app. The physics engine itself has zero dependencies.
```

## License

Nah. Public domain.
