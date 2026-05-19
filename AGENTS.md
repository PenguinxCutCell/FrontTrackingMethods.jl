# AGENTS.md

Guidance for AI coding agents working in **FrontTrackingMethods.jl**.

## What this package is

A Julia package (v0.2.0) for evolving explicit front-tracking meshes in time.
Operates on Lagrangian polygonal curves (`CurveMesh`), triangulated surfaces
(`SurfaceMesh`), and minimal 1-D point fronts (`PointFront1D`) built on top of
[FrontIntrinsicOps.jl](https://github.com/PenguinxCutCell/FrontIntrinsicOps.jl).

The user-facing API mirrors `LevelSetMethods.jl`: define terms → choose an
integrator → build a `FrontEquation` → call `integrate!`. Static geometry,
normals, curvature, and DEC operators live in `FrontIntrinsicOps`; this package
contributes the *time-evolution* layer.

## Repository layout

```
src/
  FrontTrackingMethods.jl   module + exports
  front_types.jl            abstract types
  frontfield.jl             FrontField (vertex/edge/face data)
  front1d.jl                PointFront1D support
  frontstate.jl             FrontState (mesh + fields + time)
  geometry_refresh.jl       refresh_geometry!
  frontterms.jl             AdvectionTerm, NormalMotionTerm, CurvatureMotionTerm
  timestepping.jl           ForwardEuler, RK2, RK3 + CFL adaptive Δt
  redistribution.jl         redistribute! / repair_front!
  remeshing.jl              AdaptiveCurveRemesher, ExperimentalSurfaceRemesher
  transfer.jl               transfer_vertex_field! (piecewise linear, barycentric)
  frontequation.jl          FrontEquation + integrate!
  diagnostics.jl            edge / area / volume / quality metrics
  benchmark_geometries.jl   circles, spheres, Zalesak disk/sphere, open arcs
  benchmark_fields.jl       rigid_*, rider_kothe_*, deformation16_2d, serpentine_2d, enright_3d
  compare_fronts.jl         Hausdorff / L² metrics, area/volume/centroid errors
  callbacks.jl              EveryNSteps, TimeIntervalCallback, compose_callbacks
  utils.jl                  internal helpers
  plotting_stubs.jl         names exported here, implementations live in ext/MakieExt.jl
ext/MakieExt.jl             weak Makie extension (plot_state, snapshot, record_evolution!)
test/                       runtests.jl + ~20 topical test files
examples/                   runnable demo scripts (translate/rotate/curvature/vortex/...)
benchmark/                  headless convergence + 3-D robustness studies
docs/                       Documenter site source
```

Include order in [src/FrontTrackingMethods.jl](src/FrontTrackingMethods.jl) matters —
types must be defined before consumers. If you add a file, place it in the
existing include block accordingly.

## Dependencies

- `FrontIntrinsicOps` (`= 0.1.0`) — all static geometry, normals, curvature, DEC.
- `LinearAlgebra`, `StaticArrays`, `Statistics` — standard.
- `Makie` is a **weak dep**; plotting goes through `ext/MakieExt.jl`. Do not
  add `using Makie` outside the extension. The base `src/plotting_stubs.jl`
  declares the function names so they can be exported; the extension provides
  the methods.
- Test target only adds `Test`. Do not add new test deps without updating
  `Project.toml`'s `[extras]` / `[targets]`.

## Public API (exported)

Types: `FrontField`, `FrontState`, `FrontEquation`, `AdvectionTerm`,
`NormalMotionTerm`, `CurvatureMotionTerm`, `ForwardEuler`, `RK2`, `RK3`,
`NoRedistribution`, `CurveEqualArcRedistributor`, `SurfaceTangentialRedistributor`,
`AdaptiveCurveRemesher`, `ExperimentalSurfaceRemesher`, `EveryNSteps`,
`TimeIntervalCallback`.

Core: `integrate!`, `current_state`, `current_time`, `refresh_geometry!`,
`redistribute!`, `remesh!`, `repair_front!`, `transfer_vertex_field!`,
`transfer_fields!`, `compute_cfl`, `add_field!`, `get_field`,
`vertex_coordinates`, `set_vertex_coordinates!`, `location`, `mesh`.

Diagnostics + benchmarks + comparison metrics + Makie plotting names are also
exported — see [src/FrontTrackingMethods.jl:78-178](src/FrontTrackingMethods.jl#L78-L178)
for the authoritative list.

## Sign / orientation conventions (do not violate)

**Curves (2-D, CCW closed):** `vertex_normals` point **inward**;
`signed_curvature > 0` for convex CCW (κ = 1/R on a circle); `NormalMotionTerm(Vn>0)`
shrinks the curve; `CurvatureMotionTerm(β>0)` is curve-shortening flow.

**Surfaces (3-D, outward-oriented):** `vertex_normals` point **outward**;
`mean_curvature_normal` points inward on convex surfaces;
`NormalMotionTerm(Vn>0)` expands the surface; `CurvatureMotionTerm(β>0)` is
mean-curvature flow (sphere shrinks).

Exact reference solution for curvature flow: `R(t) = √(R₀² − 2βt)` in both 2-D
and 3-D. The test suite asserts this — if you "fix" a sign and break it,
you are wrong.

## What is in scope vs. not (v0.2)

In scope: fixed-connectivity vertex motion; explicit RK; curve & surface
redistribution / remeshing without topology change; field transfer; benchmark
geometries / velocity fields / comparison metrics.

**Not in scope** (do not add unprompted): topology change (splits / collapses /
flips), merge / break / contact handling, cut-cell bulk coupling, open fronts
with boundary conditions, semi-implicit curvature, surfactant / surface PDEs.
Topology change is planned for v0.3 — see [TODO.md](TODO.md).

## Testing

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

`test/runtests.jl` preserves the original v0.1 blocks inline and then
`include`s ~20 topical files. Makie extension tests are conditionally included
only if `Makie` and `CairoMakie` are installed in the test environment.

When adding features, prefer adding a topical `test/test_*.jl` and `include`ing
it from `runtests.jl` rather than appending to the inline blocks. Keep
tolerances tight where physics gives an exact answer (rigid motion, circle
under curve-shortening): existing tests use ~1e-8 for advection and ~2% for
curvature-flow radii.

## Benchmarks and examples

- `examples/*.jl` — small, runnable demos. Use `examples/example_utils.jl`
  helpers when present rather than duplicating setup.
- `benchmark/run_all.jl` — entry point for advection convergence studies.
- `benchmark/run_{sphere_rigid,zalesak_sphere,enright3d_remeshing}_study.jl` —
  headless 3-D studies, not in CI.
- Output of benchmark / example scripts goes to `benchmark/output/` and
  `examples/output/` respectively; do not commit large artefacts there.

## Conventions for changes

- Match the existing flat-module style; no submodules.
- New geometry / velocity-field constructors go in `benchmark_geometries.jl` /
  `benchmark_fields.jl` and must be exported from
  [src/FrontTrackingMethods.jl](src/FrontTrackingMethods.jl).
- Time integrators implement the interface in `timestepping.jl`. New
  integrators must support `cfl=` and adaptive `dt` selection.
- Terms (`<:AbstractFrontTerm`) contribute a vertex velocity; the
  `FrontEquation` sums them. Add a new term by extending the dispatch in
  `frontterms.jl`.
- `PointFront1D` supports only `AdvectionTerm` and `NormalMotionTerm` (no
  curvature, no remeshing). Don't add curvature dispatches for it without a
  matching definition upstream in `FrontIntrinsicOps`.
- Update `README.md`'s feature matrix when adding user-visible features.

## Quick start (canonical pattern)

```julia
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays

mesh = make_circle_benchmark_curve(center=SVector(0.0,0.0), R=1.0, N=128)
eq   = FrontEquation(;
    terms      = (AdvectionTerm(SVector(1.0, 0.0)),),
    front      = mesh,
    t          = 0.0,
    integrator = RK2(cfl=0.8),
)
integrate!(eq, 1.0)
state = current_state(eq)
```

## Known rough edges

- Surface remeshing (`ExperimentalSurfaceRemesher`,
  `SurfaceTangentialRedistributor`) is tuning-sensitive on severe 3-D
  deformations (Enright). See [TODO.md](TODO.md) — redistribution behaviour in
  the Enright test is flagged for revision.
- `transfer_vertex_field!` `:barycentric` on surfaces beats `:nearest` on
  smooth fields and preserves constants exactly; do not regress this property.
