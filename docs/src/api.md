# API summary

The package exposes a compact, high-level API. This page lists the primary
types and functions; consult docstrings in `src/` for details.

Core types and constructors

- `FrontField(values, mesh, location)` — attach per-vertex/edge/face fields.
- `FrontState(mesh; t=0.0)` — front geometry with attached fields and time.
- `FrontEquation(; front, terms, integrator, redistributor, callback...)` — top-level equation object.

Motion terms

- `AdvectionTerm(u)` — prescribed vector velocity (constant or function).
- `NormalMotionTerm(Vn)` — scalar normal speed (positive means inward for CCW curves).
- `CurvatureMotionTerm(β)` — curvature-driven normal motion.

Integrators

- `ForwardEuler`, `RK2`, `RK3` — explicit time-stepping integrators.

Redistributors / Remeshing

- `NoRedistribution` — do nothing.
- `CurveEqualArcRedistributor` — equal-arc redistribution on curves.
- `AdaptiveCurveRemesher` — fixed-topology adaptive remesher with corner protection.
- `SurfaceTangentialRedistributor` — tangential smoothing on surfaces.
- `ExperimentalSurfaceRemesher` — experimental fixed-topology surface quality remesher.

Core functions

- `integrate!(eq, T; dt, callback...)` — advance solution to time `T`.
- `current_state(eq)` / `current_time(eq)` — query running equation.
- `refresh_geometry!(state)` — recompute geometry (normals, curvature, areas).
- `redistribute!(state, redistributor)` / `remesh!(...)` — redistribute or remesh.
- `transfer_vertex_field!(...)` — transfer fields between meshes (`:piecewise_linear`, `:nearest_vertex`, `:barycentric`).
- `add_field!(state, name, values; location=:vertex)` — attach a field.

Diagnostics & metrics

- `front_enclosed_measure`, `front_measure`, `front_centroid`, `edge_length_spread`
- `surface_edge_length_stats`, `surface_triangle_area_stats`, `surface_triangle_angle_stats`
- `surface_aspect_ratio_stats`, `surface_degenerate_fraction`, `surface_quality_summary`
- `symmetric_hausdorff_curve`, `l2_distance_curve`, `relative_area_error`

Callbacks

- `EveryNSteps(n, callback)` — call `callback(state, t, step)` every `n` steps.
- `TimeIntervalCallback(dt, callback)` — call `callback` at wall-clock time intervals.

For function signatures and examples, inspect the docstrings in the
corresponding `src/*.jl` files (the module exports a compact set of names).
