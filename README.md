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

## Features

| Feature | Status |
|---------|--------|
| Closed 2-D curves (CurveMesh) | ✅ Implemented |
| Closed 3-D surfaces (SurfaceMesh) | ✅ Implemented |
| `AdvectionTerm` (prescribed vector velocity) | ✅ Implemented |
| `NormalMotionTerm` (prescribed scalar normal speed) | ✅ Implemented |
| `CurvatureMotionTerm` (mean-curvature / curve-shortening flow) | ✅ Implemented |
| `ForwardEuler` / `RK2` / `RK3` integrators | ✅ Implemented |
| CFL-based adaptive time-stepping | ✅ Implemented |
| `CurveEqualArcRedistributor` | ✅ Implemented |
| `SurfaceTangentialRedistributor` | ⚠️ Experimental |
| Field transfer on curves | ✅ Implemented |
| Field transfer on surfaces | ⚠️ Nearest-vertex only |
| Topology change / remeshing | ❌ Not in v0.1 |
| Contact / collision handling | ❌ Not in v0.1 |
| Cut-cell bulk coupling | ❌ Out of scope |

## Relationship to LevelSetMethods.jl

| LevelSetMethods concept | FrontTrackingMethods concept |
|-------------------------|------------------------------|
| `CartesianGrid` | `CurveMesh` / `SurfaceMesh` (from FrontIntrinsicOps) |
| `MeshField` | `FrontField` attached to vertices / edges / faces |
| `LevelSet` | `FrontState` (front geometry + fields + time) |
| `AdvectionTerm` | `AdvectionTerm` (vertex advection) |
| `NormalMotionTerm` | `NormalMotionTerm` (normal speed × normal) |
| `CurvatureTerm` | `CurvatureMotionTerm` (β × mean-curvature-normal) |
| `reinitialize!` | `redistribute!` / `repair_front!` |
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
  FrontTrackingMethods.jl   # module + exports
  front_types.jl            # abstract types
  frontfield.jl             # FrontField
  frontstate.jl             # FrontState
  geometry_refresh.jl       # refresh_geometry!
  frontterms.jl             # AdvectionTerm, NormalMotionTerm, CurvatureMotionTerm
  timestepping.jl           # ForwardEuler, RK2, RK3
  frontequation.jl          # FrontEquation + integrate!
  redistribution.jl         # redistribute! / repair_front!
  transfer.jl               # transfer_vertex_field!
  diagnostics.jl            # quality metrics
  utils.jl                  # internal helpers
```

## Dependencies
- `FrontIntrinsicOps.jl` – all static geometry, normals, curvature, DEC operators
- `LinearAlgebra` – dot, norm, cross
- `StaticArrays` – SVector for vertex coordinates

## Non-goals (v0.1)
- Topology change (edge splits / collapses / flips)
- Merge / break / contact handling
- Bulk cut-cell coupling
- Open fronts with boundary conditions
