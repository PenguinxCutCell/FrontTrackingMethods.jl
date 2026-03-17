# Plotting and animations

`FrontTrackingMethods.jl` enables plotting through a Makie weak extension.
That means plotting APIs become available when a Makie backend is loaded,
for example `CairoMakie`.

```julia
using CairoMakie
using FrontIntrinsicOps
using FrontTrackingMethods
```

## Static plots

```julia
using StaticArrays

mesh = make_circle_benchmark_curve(center=SVector(0.0, 0.0), R=1.0, N=128)
state = FrontState(mesh)

set_makie_theme!()
fig, ax, p = plot_state(state; show_vertices=true, title="initial front")
save("state.png", fig)
```

You can also plot equation state directly:

```julia
eq = FrontEquation(; terms=(AdvectionTerm(SVector(1.0, 0.0)),), front=mesh)
fig, ax, p = plot_equation(eq)
save("equation.png", fig)
```

## Snapshots and recording

High-level helpers are provided by the extension:

- `snapshot(obj, filename)` for `CurveMesh`, `SurfaceMesh`, `FrontState`, or `FrontEquation`
- `record_evolution!(eq, filename, times)` to integrate and record a movie

```julia
times = range(0.0, stop=1.0, length=120)
snapshot(eq, "initial.png")
record_evolution!(eq, "front.mp4", times)
snapshot(eq, "final.png")
```

## Notes

- Plotting remains optional: core simulation APIs work without Makie.
- `CairoMakie` is only needed in examples/docs/tests; the package weak-dep is `Makie`.
- For low-level mesh plotting (`CurveMesh`, `SurfaceMesh`, normals, wireframes),
  see `FrontIntrinsicOps` plotting helpers.
