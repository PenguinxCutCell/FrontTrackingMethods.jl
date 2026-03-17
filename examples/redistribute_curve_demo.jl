# redistribute_curve_demo.jl – Start with clustered markers, show improvement.
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays, LinearAlgebra
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

N = 32; R = 1.0
# Cluster 24 vertices in the first quarter arc
angles = vcat(LinRange(0, π/2, 24), LinRange(π/2, 2π, 10)[2:end])[1:N]
pts   = [SVector(R*cos(θ), R*sin(θ)) for θ in angles]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
mesh  = CurveMesh(pts, edges)
state = FrontState(mesh)

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("redistribute_curve_demo")
snapshot(state, joinpath(outdir, "initial.png"); title="redistribution: before", show_vertices=true)

println("Before redistribution:")
println("  min edge = ", round(min_edge_length(state), digits=5))
println("  max edge = ", round(max_edge_length(state), digits=5))
println("  spread   = ", round(edge_length_spread(state), digits=3))

redistribute!(state, CurveEqualArcRedistributor())

println("After redistribution:")
println("  min edge = ", round(min_edge_length(state), digits=5))
println("  max edge = ", round(max_edge_length(state), digits=5))
println("  spread   = ", round(edge_length_spread(state), digits=3))

snapshot(state, joinpath(outdir, "final.png"); title="redistribution: after", show_vertices=true)
println("Saved plotting outputs to: $outdir")
