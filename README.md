# FrontTrackingMethods.jl

A Julia package for evolving explicit front-tracking meshes in time.
Operates on Lagrangian polygonal curves and triangulated surfaces
built on top of [FrontIntrinsicOps.jl](https://github.com/PenguinxCutCell/FrontIntrinsicOps.jl).

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

## v0.2 Feature Matrix

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
| `ExperimentalSurfaceRemesher` | ⚠️ Experimental | v0.2 new; edge flips + smoothing |
| Field transfer – piecewise linear (curves) | ✅ Stable | v0.2 improved |
| Field transfer – barycentric (surfaces) | ✅ Stable | v0.2 new |
| Field transfer – nearest-vertex fallback | ✅ Stable | |
| Benchmark geometry constructors | ✅ Stable | v0.2 new |
| Benchmark velocity fields | ✅ Stable | v0.2 new |
| Front-to-front comparison metrics | ✅ Stable | v0.2 new |
| Lightweight callbacks (`EveryNSteps`) | ✅ Stable | v0.2 new |
| Topology change / remeshing | ❌ Not in v0.2 | planned for v0.3 |
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

Examples are in the `examples/` directory.  Long convergence scripts are in
`benchmarks/`.

## What is now strong / still experimental / missing

| Category | Status |
|----------|--------|
| 2-D rigid motion tests | **Strong** – tight tolerances, fully automated |
| 2-D reversible deformation tests | **Strong** – good with redistribution |
| 3-D rigid motion tests | **Good** – works on coarse meshes |
| 3-D deformation tests | **Moderate** – needs remeshing; coarse only in CI |
| Surface remeshing quality | **Experimental** – works but may change API |
| Barycentric surface transfer | **Good** – materially better than nearest-vertex |
| Topology change | **Missing** – out of v0.2 scope |

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
