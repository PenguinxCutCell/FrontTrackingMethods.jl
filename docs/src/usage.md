# Quick start / Usage

This quick start demonstrates the typical workflow: construct a front,
define motion terms, assemble a `FrontEquation`, and call `integrate!`.

```julia
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays

# Build a circle mesh
N = 64; R = 1.0
pts   = [SVector(R*cos(2π*(i-1)/N), R*sin(2π*(i-1)/N)) for i in 1:N]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
mesh  = CurveMesh(pts, edges)

# Define motion: constant advection
u_adv = SVector(0.2, -0.1)
terms = (AdvectionTerm(u_adv),)

# Build equation and integrate for T=1.0
eq = FrontEquation(; terms=terms, front=mesh, t=0.0, integrator=RK2())
integrate!(eq, 1.0)

# Inspect final state
state = current_state(eq)
println("centroid: ", front_centroid(state))
```

Redistribution and remeshing are purely numerical: call `redistribute!`
or pass a `redistributor` to `FrontEquation` or use callbacks such as
`EveryNSteps`.

## Surface field transfer

When redistribution/remeshing moves vertices, transfer vertex fields with
`transfer_vertex_field!`:

```julia
transfer_vertex_field!(dst, oldmesh, srcvals, newmesh; method=:barycentric)
```

- `:nearest_vertex` is the simplest fallback.
- `:barycentric` is generally preferred on surfaces for smooth fields and
	preserves constants exactly.
