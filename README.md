# FrontTrackingMethods.jl

A Julia package for evolving explicit front-tracking meshes in time.
Operates on Lagrangian polygonal curves and triangulated surfaces
built on top of [FrontIntrinsicOps.jl](https://github.com/PenguinxCutCell/FrontIntrinsicOps.jl).

## Optional Makie support

Plotting and animation helpers are provided via a weak `Makie` extension.
Load a Makie backend (for example `CairoMakie`) to enable:

- `plot_state`, `plot_equation`
- `snapshot`
- `record_evolution!`

```julia
using CairoMakie
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays

mesh = make_circle_benchmark_curve(center=SVector(0.0, 0.0), R=1.0, N=128)
eq = FrontEquation(; terms=(AdvectionTerm(SVector(1.0, 0.0)),), front=mesh)

snapshot(eq, "initial.png")
record_evolution!(eq, "front.mp4", range(0.0, stop=1.0, length=120))
snapshot(eq, "final.png")
```

## Quick start

```julia
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays

# 1. Build a closed polygonal circle
N = 64; R = 1.0
pts   = [SVector(R*cos(2π*(i-1)/N), R*sin(2π*(i-1)/N)) for i in 1:N]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
mesh  = CurveMesh(pts, edges)

# 2. Define motion terms
u_adv = SVector(1.0, 0.5)                 # constant advection velocity
terms = (AdvectionTerm(u_adv),)

# 3. Build equation and integrate
eq = FrontEquation(; terms=terms, front=mesh, t=0.0, integrator=RK2())
integrate!(eq, 1.0)

# 4. Inspect result
state = current_state(eq)
println("centroid: ", sum(state.mesh.points) / length(state.mesh.points))
```

## v0.3 Feature Matrix

| Feature | Status | Notes |
|---------|--------|-------|
| Closed 2-D curves (`CurveMesh`) | ✅ Stable | |
| Closed 3-D surfaces (`SurfaceMesh`) | ✅ Stable | |
| `AdvectionTerm` (prescribed velocity) | ✅ Stable | constant, function, or `FrontField` |
| `NormalMotionTerm` (scalar normal speed) | ✅ Stable | |
| `CurvatureMotionTerm` (mean-curvature flow) | ✅ Stable | explicit only |
| `ForwardEuler` / `RK2` / `RK3` integrators | ✅ Stable | |
| CFL-based adaptive time-stepping | ✅ Stable | |
| `CurveEqualArcRedistributor` | ✅ Stable | |
| `AdaptiveCurveRemesher` | ✅ Stable | v0.2 new; corner protection |
| `SurfaceTangentialRedistributor` | ⚠️ Experimental | useful for modest deformation |
| `ExperimentalSurfaceRemesher` | ⚠️ Experimental | fixed-topology tangential + edge-length regularization |
| Field transfer – piecewise linear (curves) | ✅ Stable | v0.2 improved |
| Field transfer – barycentric (surfaces) | ✅ Stable | validated vs nearest-vertex on smooth fields |
| Field transfer – nearest-vertex fallback | ✅ Stable | |
| Benchmark geometry constructors | ✅ Stable | v0.2 new |
| Benchmark velocity fields | ✅ Stable | v0.2 new |
| Front-to-front comparison metrics | ✅ Stable | v0.2 new |
| Lightweight callbacks (`EveryNSteps`) | ✅ Stable | v0.2 new |
| Multi-component states (`MultiFrontState`) | ✅ New in v0.3 | additive to `FrontState` |
| Geometry-driven 2-D topology change | ✅ Experimental / usable | local Cartesian reconstruction |
| Geometry-driven 3-D topology change | ⚠️ Prototype | coarse whole-component reconstruction |
| Film-model-based coalescence / rupture | ❌ Not in v0.3 | future work |
| Topology change / remeshing | ⚠️ Partial in v0.3 | remeshing remains separate from event detection |
| Contact / collision handling | ❌ Not in v0.2 | |
| Cut-cell bulk coupling | ❌ Out of scope | |

## v0.2 Benchmark Gallery

The following standard benchmarks are implemented as reusable setups
in `src/benchmark_geometries.jl` and `src/benchmark_fields.jl`:

| Benchmark | Type | Field | Key metric |
|-----------|------|-------|------------|
| Circle rigid rotation | 2-D | `rigid_rotation_2d` | Hausdorff to initial |
| Zalesak disk rotation | 2-D | `rigid_rotation_2d` | area drift + corner preservation |
| Rider–Kothe vortex | 2-D | `rider_kothe_single_vortex` | shape recovery with redistribution |
| 16-vortex deformation | 2-D | `deformation16_2d` | min edge length / shape recovery |
| Serpentine deformation | 2-D | `serpentine_2d` | L² and Hausdorff to initial |
| Sphere translation | 3-D | `rigid_translation_velocity` | centroid error, volume drift |
| Sphere rotation | 3-D | `rigid_rotation_3d` | shape, volume preservation |
| Zalesak sphere | 3-D | `rigid_rotation_3d` | volume drift + qualitative slot |
| Enright deformation | 3-D | `enright_3d` | volume drift, mesh quality |

Headless 3-D study scripts for stronger quantitative evidence (outside default CI):

- `benchmark/run_sphere_rigid_study.jl`
- `benchmark/run_zalesak_sphere_study.jl`
- `benchmark/run_enright3d_remeshing_study.jl`

Examples are in the `examples/` directory.  Advection benchmark drivers and
error scripts are in `benchmark/` (`julia --project=. benchmark/run_all.jl`).

## What is now strong / still experimental / missing

| Category | Status |
|----------|--------|
| 2-D rigid motion tests | **Strong** – tight tolerances, fully automated |
| 2-D reversible deformation tests | **Strong** – good with redistribution |
| 3-D rigid motion tests | **Strong** – centroid/area/volume/quality checks in CI |
| Mild 3-D prescribed deformations | **Usable** – robust on coarse/moderate meshes |
| Severe 3-D deformation with remeshing | **Experimental** – benchmarked, tuning-sensitive |
| Surface remeshing quality | **Experimental** – fixed-topology conservative updates |
| Barycentric surface transfer | **Strong** – constants preserved; smooth fields beat nearest |
| 2-D topology change | **Usable (experimental)** – geometry-driven local reconstruction |
| 3-D topology change | **Prototype** – coarse smoke/regression support |
| Film drainage / rupture physics | **Missing** – intentionally out of v0.3 scope |

## Topology change in v0.3

v0.3 keeps explicit front tracking as the primary representation and introduces
an additive geometry-driven topology pipeline:

- Normal evolution remains explicit front tracking on `CurveMesh` / `SurfaceMesh`.
- Imminent events are detected geometrically.
- A local Cartesian patch is rasterized with an indicator field.
- Marching-squares-style (2-D) / coarse marching-cubes-style (3-D prototype)
  reconstruction creates post-event components.
- Reconstructed components replace old components in `MultiFrontState`.

Current scope and limits:

- 2-D merge/split is the main validated target for v0.3.
- 3-D support is intentionally prototype-grade and conservative.
- Physical film-drainage / rupture models are **not** included in v0.3.

## Relationship to LevelSetMethods.jl

| LevelSetMethods concept | FrontTrackingMethods concept |
|-------------------------|------------------------------|
| `CartesianGrid` | `CurveMesh` / `SurfaceMesh` (from FrontIntrinsicOps) |
| `MeshField` | `FrontField` attached to vertices / edges / faces |
| `LevelSet` | `FrontState` (front geometry + fields + time) |
| `AdvectionTerm` | `AdvectionTerm` (vertex advection) |
| `NormalMotionTerm` | `NormalMotionTerm` (normal speed × normal) |
| `CurvatureTerm` | `CurvatureMotionTerm` (β × mean-curvature-normal) |
| `reinitialize!` | `redistribute!` / `remesh!` |
| `extend_along_normals!` | `transfer_vertex_field!` / `transfer_fields!` |
| `LevelSetEquation` | `FrontEquation` |
| `integrate!` | `integrate!` |

**Same high-level philosophy**: define terms → define integrator → build equation → call `integrate!`.

**Different representation**: Eulerian level-set field on a Cartesian grid
vs. Lagrangian vertex positions on an explicit polygonal / triangulated mesh.

**Reinitialization replaced by redistribution**: in level-set methods,
reinitialization restores the signed-distance property.  In front tracking,
redistribution restores mesh quality (equal arc-length / tangential smoothing)
without changing the physical front shape.

**Velocity extension replaced by field transfer**: in level-set methods,
the velocity field is extended off the interface via `extend_along_normals!`.
In front tracking, after redistribution moves vertices, attached fields are
transferred to the new vertex positions via `transfer_vertex_field!`.

## Sign conventions

### CurveMesh (2-D curves)
- `vertex_normals`: **inward** unit normals for CCW-oriented closed curves.
- `signed_curvature`: **positive** for convex CCW curves (κ = 1/R for a circle).
- `NormalMotionTerm(Vn)`: `Vn > 0` → vertices move **inward** → circle **shrinks**.
- `CurvatureMotionTerm(β)`: `β > 0` → curve-shortening flow → circle **shrinks**.
  - Exact: R(t) = √(R₀² - 2βt).

### SurfaceMesh (3-D surfaces)
- `vertex_normals`: **outward** unit normals for closed outward-oriented surfaces.
- `mean_curvature_normal`: points **inward** for convex outward-oriented surfaces.
- `NormalMotionTerm(Vn)`: `Vn > 0` → vertices move **outward** → sphere **expands**.
- `CurvatureMotionTerm(β)`: `β > 0` → mean-curvature flow → sphere **shrinks**.
  - Exact: R(t) = √(R₀² - 2βt).

## Package structure

```
src/
  FrontTrackingMethods.jl      # module + exports
  front_types.jl               # abstract types
  frontfield.jl                # FrontField
  frontstate.jl                # FrontState
  geometry_refresh.jl          # refresh_geometry!
  frontterms.jl                # AdvectionTerm, NormalMotionTerm, CurvatureMotionTerm
  timestepping.jl              # ForwardEuler, RK2, RK3
  frontequation.jl             # FrontEquation + integrate!
  redistribution.jl            # redistribute! / repair_front!
  remeshing.jl                 # AdaptiveCurveRemesher, ExperimentalSurfaceRemesher (v0.2)
  transfer.jl                  # transfer_vertex_field! (v0.2: piecewise linear + barycentric)
  diagnostics.jl               # quality metrics
  benchmark_geometries.jl      # Zalesak disk/sphere, circle, sphere constructors (v0.2)
  benchmark_fields.jl          # rigid/vortex/deformation velocity fields (v0.2)
  compare_fronts.jl            # Hausdorff / L² front-to-front metrics (v0.2)
  callbacks.jl                 # EveryNSteps, TimeIntervalCallback (v0.2)
  utils.jl                     # internal helpers
```

## Dependencies
- `FrontIntrinsicOps.jl` – all static geometry, normals, curvature, DEC operators
- `LinearAlgebra` – dot, norm, cross
- `StaticArrays` – SVector for vertex coordinates

## Non-goals (v0.2)
- Topology change (edge splits / collapses / flips)
- Merge / break / contact handling
- Bulk cut-cell coupling
- Open fronts with boundary conditions
- Semi-implicit curvature schemes
- Surfactant / surface PDEs
