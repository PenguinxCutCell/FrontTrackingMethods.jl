# normal_advection_poisson_redistribution.jl
#
# Demonstrates the split workflow:
#   1. projected normal advection of an ellipse in a rigid-rotation field
#   2. pure Poisson tangential redistribution on a static circle with clustered markers
#
# Outputs:
#   examples/output/normal_advection_poisson_redistribution/ellipse_normal_advection.gif
#   examples/output/normal_advection_poisson_redistribution/circle_poisson_redistribution.gif

using CairoMakie
using FrontIntrinsicOps
using FrontTrackingMethods
using LinearAlgebra
using Printf
using StaticArrays
using Statistics

include("example_utils.jl")
using .ExampleUtils

FrontTrackingMethods.set_makie_theme!()

const OUTDIR = ExampleUtils.output_dir_for("normal_advection_poisson_redistribution")

function make_ellipse_curve(;
    center::SVector{2,Float64}=SVector(0.5, 0.75),
    a::Float64=0.30,
    b::Float64=0.15,
    N::Int=96,
)
    pts = [
        center + SVector{2,Float64}(a * cos(2π * k / N), b * sin(2π * k / N))
        for k in 0:N-1
    ]
    edges = [SVector{2,Int}(i, mod1(i + 1, N)) for i in 1:N]
    return CurveMesh{Float64}(pts, edges)
end

function make_clustered_circle(;
    center::SVector{2,Float64}=SVector(0.0, 0.0),
    R::Float64=1.0,
    N::Int=72,
    cluster_strength::Float64=0.65,
)
    # Smooth monotone reparametrization: derivative is proportional to
    # 1 + cluster_strength*cos(theta), so markers concentrate near theta = pi.
    angles = [2π * k / N + cluster_strength * sin(2π * k / N) for k in 0:N-1]
    pts = [center + R * SVector{2,Float64}(cos(θ), sin(θ)) for θ in angles]
    edges = [SVector{2,Int}(i, mod1(i + 1, N)) for i in 1:N]
    return CurveMesh{Float64}(pts, edges)
end

closed_points(mesh::CurveMesh) = begin
    pts = [Point2f(p[1], p[2]) for p in mesh.points]
    push!(pts, first(pts))
    pts
end

vertex_points(mesh::CurveMesh) = [Point2f(p[1], p[2]) for p in mesh.points]

function edge_uniformity(mesh::CurveMesh)
    geom = compute_geometry(mesh)
    return std(geom.edge_lengths) / mean(geom.edge_lengths)
end

function record_curve_gif!(
    filename::AbstractString,
    get_mesh::Function,
    step!::Function;
    nframes::Int,
    title::AbstractString,
    xlims,
    ylims,
    reference::Union{Nothing,CurveMesh}=nothing,
    color=:royalblue,
)
    mesh0 = get_mesh()
    fig = Figure(size=(720, 640))
    ax = Axis(fig[1, 1]; title=title, xlabel="x", ylabel="y", aspect=DataAspect())
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)

    if reference !== nothing
        lines!(ax, closed_points(reference); color=(:gray45, 0.35), linewidth=2, linestyle=:dash)
    end

    line_obs = Observable(closed_points(mesh0))
    marker_obs = Observable(vertex_points(mesh0))
    lines!(ax, line_obs; color=color, linewidth=3)
    scatter!(ax, marker_obs; color=:black, markersize=7)

    label = Label(fig[2, 1], ""; tellwidth=false)

    record(fig, filename, 1:nframes; framerate=24) do frame
        frame > 1 && step!(frame)
        mesh = get_mesh()
        line_obs[] = closed_points(mesh)
        marker_obs[] = vertex_points(mesh)
        label.text[] = @sprintf("frame %03d/%03d    edge uniformity std/mean = %.4f",
                                frame, nframes, edge_uniformity(mesh))
    end
    return filename
end

function run_ellipse_normal_advection_example()
    center_ellipse = SVector(0.5, 0.75)
    center_rot = SVector(0.5, 0.5)
    mesh0 = make_ellipse_curve(center=center_ellipse, a=0.30, b=0.15, N=96)
    velocity = rigid_rotation_2d(; center=center_rot, omega=2π)

    eq = FrontEquation(;
        terms = ProjectedAdvectionTerm(velocity),
        front = deepcopy(mesh0),
        integrator = RK2(),
        redistributor = PoissonTangentialRedistributor(;
            iterations=5,
            pseudo_dt=0.05,
            omega=1.0,
            every=1,
            max_step_fraction=0.20,
        ),
    )

    tf = 1.0
    nframes = 80
    times = collect(range(0.0, tf; length=nframes))
    dt = 0.001

    outfile = joinpath(OUTDIR, "ellipse_normal_advection.gif")
    record_curve_gif!(
        outfile,
        () -> current_state(eq).mesh,
        frame -> integrate!(eq, times[frame]; dt=dt),
        nframes=nframes,
        title="Projected normal advection + Poisson redistribution",
        xlims=(0.00, 0.95),
        ylims=(0.00, 1.05),
        reference=mesh0,
        color=:dodgerblue3,
    )

    state = current_state(eq)
    @printf("[ellipse] final time %.3f, edge uniformity %.4f, area %.6f\n",
            current_time(eq), edge_uniformity(state.mesh), front_enclosed_measure(state; correction=:arc))
    println("Saved ellipse GIF: $outfile")
    return outfile
end

function run_circle_redistribution_example()
    mesh0 = make_clustered_circle(N=72, cluster_strength=0.65)
    state = FrontState(mesh0)
    redist = PoissonTangentialRedistributor(;
        iterations=1,
        pseudo_dt=0.18,
        omega=1.0,
        max_step_fraction=0.20,
    )

    outfile = joinpath(OUTDIR, "circle_poisson_redistribution.gif")
    record_curve_gif!(
        outfile,
        () -> state.mesh,
        frame -> redistribute!(state, redist),
        nframes=80,
        title="Static circle: Poisson tangential marker redistribution",
        xlims=(-1.2, 1.2),
        ylims=(-1.2, 1.2),
        reference=mesh0,
        color=:seagreen4,
    )

    @printf("[circle] edge uniformity %.4f -> %.4f, area drift %.3e\n",
            edge_uniformity(mesh0),
            edge_uniformity(state.mesh),
            abs(enclosed_measure(state.mesh; correction=:arc) - π) / π)
    println("Saved circle GIF: $outfile")
    return outfile
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Writing outputs to: $OUTDIR")
    run_ellipse_normal_advection_example()
    run_circle_redistribution_example()
end
