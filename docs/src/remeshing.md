# Redistribution and Remeshing

This page explains the difference between physical motion and numerical
redistribution / remeshing, and documents the redistribution and remeshing
strategies available in `FrontTrackingMethods.jl` v0.2.

---

## Physical motion vs. numerical redistribution

**Physical motion** is governed by the PDE:
```
dx/dt = u(x, t)
```
where `u` is the prescribed velocity field.  The integrators (`RK2`, `RK3`,
`ForwardEuler`) advance the marker positions according to this equation.

**Redistribution / remeshing** is a purely numerical operation that
repositions or restructures the front mesh without affecting the physical
interface.  Its purpose is to maintain good mesh quality so that the physical
motion can be integrated accurately.

> **Key principle:** redistribution must not be confused with topology change
> or interface diffusion.  A good redistributor preserves the interface
> geometry while improving marker distribution.

Redistribution is triggered by calling `redistribute!(state, redistributor)`
or by supplying a `redistributor` to `FrontEquation`.

---

## Curve redistributors

### `CurveEqualArcRedistributor`

Redistributes markers uniformly along the arc length of the current curve.
Vertex count is preserved.

```julia
r = CurveEqualArcRedistributor(; every=1)
redistribute!(state, r)
```

**Properties:**
- preserves closure
- preserves vertex count
- area change is typically < 1%
- does not protect corners

**Recommended for:** smooth interfaces where uniform spacing is sufficient.

### `AdaptiveCurveRemesher` (v0.2)

A fixed-topology adaptive remesher for curves.

```julia
r = AdaptiveCurveRemesher(;
    iterations      = 5,      # number of smoothing passes
    protect_corners = true,   # protect sharp-corner vertices from smoothing
    corner_angle    = π/4,    # turning angle threshold for corner detection
)
redistribute!(state, r)
```

**Properties:**
- fixed topology (no vertex insertion / deletion by default)
- optional corner protection for Zalesak-like geometries
- preserves closure
- area change is typically < 2% with default settings

**Recommended for:** Zalesak disk, serpentine flow, any interface with sharp
geometric features that should be preserved.

---

## Surface redistributors

### `SurfaceTangentialRedistributor`

Laplacian tangential smoothing for surface meshes.  Vertices are moved in the
tangent plane of the local surface; the normal component is projected back
onto the surface.

```julia
r = SurfaceTangentialRedistributor(; iterations=3, strength=0.5)
redistribute!(state, r)
```

**Status:** usable for mild/moderate 3-D deformation in v0.2.  
**Properties:**
- no topology change
- volume change is typically < 5%
- strength in [0, 1]; lower values are more conservative

**Recommended for:** sphere translation/rotation and modest deformation tests.

### `ExperimentalSurfaceRemesher` *(experimental)*

A conservative fixed-topology surface remesher for robustness studies.
It combines:

1. tangential centroid relaxation,
2. edge-length regularization toward `[lmin, lmax]`,
3. local area-based move rejection to avoid degeneracy,
4. optional mild global volume correction.

```julia
r = ExperimentalSurfaceRemesher(;
    iterations=3,
    strength=0.25,
    lmin=nothing,
    lmax=nothing,
    min_area=1e-14,
    volume_correction=true,
    volume_relaxation=0.15,
)
redistribute!(state, r)
```

**Status:** experimental in v0.2.  API may change.  
**Properties:**
- no topology change
- deterministic conservative updates
- avoids accepting local moves that create near-degenerate triangles
- useful for benchmark robustness, not yet a full production remeshing framework

**Threshold/tuning notes:**
- `lmin`, `lmax`: if omitted, inferred from current mean edge length.
- `min_area`: larger values are safer but may under-correct poor meshes.
- `max_tangential_step`: limits per-pass vertex displacement.
- `volume_correction=true` is only applied on closed surfaces.

**Recommended for:** Enright 3-D deformation benchmark on longer runs.

---

## Fixed-topology scope

v0.2 remeshing is **fixed topology**: no edge insertion, edge collapse, edge
flip, or hole creation. Vertex/face connectivity is preserved for surface
redistributors and remeshers.

Topology-changing remeshing (e.g., vertex insertion for long edges) is planned
for v0.3 and is **not** part of this release.

---

## Controlling redistribution frequency

Redistribution can be controlled in two ways:

### Per-step via `redistribute_every` (CurveEqualArcRedistributor)

```julia
r = CurveEqualArcRedistributor(; every=5)  # redistribute every 5 steps
```

### Via `EveryNSteps` callback

```julia
cb = EveryNSteps(5, (state, t, step) -> redistribute!(state, my_remesher))
integrate!(eq, T; dt=dt, callback=cb)
```

This pattern is useful when you want to redistribute less frequently than every
step, or when redistribution should depend on mesh quality.

---

## `remesh!` alias

`remesh!` is an alias for `redistribute!`:

```julia
remesh!(state, r)   # same as redistribute!(state, r)
```
