using Test
using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# runtests.jl – Test dispatcher for FrontTrackingMethods v0.2
#
# Test files
# ----------
#   test_utils.jl                  – shared helpers (geometry constructors, assertions)
#   test_frontfield.jl             – FrontField API
#   test_frontstate.jl             – FrontState API
#   test_advection_terms.jl        – AdvectionTerm (pure translation)
#   test_normal_motion.jl          – NormalMotionTerm
#   test_curvature_motion.jl       – CurvatureMotionTerm
#   test_curve_redistribution.jl   – CurveEqualArcRedistributor, AdaptiveCurveRemesher
#   test_surface_redistribution.jl – SurfaceTangentialRedistributor, ExperimentalSurfaceRemesher
#   test_surface_quality_metrics.jl – reusable surface quality diagnostics
#   test_transfer_curve.jl         – field transfer on curves
#   test_transfer_surface.jl       – field transfer on surfaces
#   test_benchmark_rigid_2d.jl     – rigid 2-D translation / rotation
#   test_benchmark_open_curve_2d.jl – open-curve topology/advection benchmark
#   test_benchmark_rigid_3d.jl     – rigid 3-D translation / rotation
#   test_benchmark_zalesak_2d.jl   – Zalesak disk one rotation
#   test_benchmark_zalesak_sphere_3d.jl – slotted-sphere rigid-rotation benchmark
#   test_benchmark_vortex_2d.jl    – Rider–Kothe reversed vortex
#   test_benchmark_serpentine_2d.jl – serpentine deformation smoke test
#   test_benchmark_enright_3d.jl   – Enright 3-D deformation smoke test
# ─────────────────────────────────────────────────────────────────────────────

include("test_utils.jl")
using .TestUtils

# ── Backward-compatible helpers (preserved from v0.1 test suite) ─────────────
# These are kept here so the original test blocks below continue to work.

function make_circle(R=1.0, N=64)
    T = Float64
    angles = [2π * (i-1) / N for i in 1:N]
    pts   = [SVector{2,T}(R * cos(θ), R * sin(θ)) for θ in angles]
    edges = [SVector{2,Int}(i, mod1(i+1, N)) for i in 1:N]
    return CurveMesh(pts, edges)
end

# ─────────────────────────────────────────────────────────────────────────────
# Original v0.1 test blocks (preserved)
# ─────────────────────────────────────────────────────────────────────────────

@testset "Constructor / API" begin
    mesh = make_circle(1.0, 32)
    state = FrontState(mesh; t=0.0)
    @test current_time(state) == 0.0
    @test current_state(state) === state
    @test length(state.mesh.points) == 32

    # FrontField
    vals = ones(32)
    ff = FrontField(vals, mesh, :vertex)
    @test location(ff) == :vertex
    @test length(ff) == 32
    @test all(Base.values(ff) .== 1.0)

    # add_field! / get_field
    add_field!(state, :test_field, ff)
    f2 = get_field(state, :test_field)
    @test f2 === ff

    # FrontEquation
    u = SVector(1.0, 0.0)
    eq = FrontEquation(;
        terms=(AdvectionTerm(u),),
        front=mesh,
        t=0.0,
        integrator=RK2()
    )
    @test current_time(eq) == 0.0
    @test current_state(eq) isa FrontState
end

@testset "Pure translation – circle" begin
    R  = 1.0
    N  = 64
    Ux = 0.7
    Uy = -0.3
    tf = 0.5

    # Exact centroid after advection
    exact_center = SVector(Ux * tf, Uy * tf)

    for IntegT in (ForwardEuler(cfl=0.5), RK2(cfl=0.8), RK3(cfl=1.0))
        mesh = make_circle(R, N)
        u    = SVector(Ux, Uy)
        eq   = FrontEquation(;
            terms=(AdvectionTerm(u),),
            front=mesh, t=0.0, integrator=IntegT
        )
        integrate!(eq, tf; dt=0.01)

        pts   = eq.state.mesh.points
        center = sum(pts) / length(pts)
        @test norm(center - exact_center) < 1e-8

        # Radius should be preserved
        radii = [norm(p - center) for p in pts]
        @test maximum(abs.(radii .- R)) < 1e-8

        # Enclosed area should be preserved (polygon area ≈ πR², within polygon discretization error)
        A = front_enclosed_measure(eq.state)
        @test abs(A - π*R^2) / (π*R^2) < 0.005
    end
end

@testset "Constant normal motion – circle (inward, shrink)" begin
    R0 = 1.0
    Vn = 0.3
    tf = 0.5
    N  = 128
    dt = 1e-3

    mesh = make_circle(R0, N)
    eq   = FrontEquation(;
        terms=(NormalMotionTerm(Vn),),
        front=mesh, t=0.0, integrator=RK2(cfl=0.8)
    )
    integrate!(eq, tf; dt=dt)

    pts    = eq.state.mesh.points
    center = sum(pts) / length(pts)
    radii  = [norm(p - center) for p in pts]
    R_mean = sum(radii) / length(radii)

    @test norm(center) < 1e-3

    R_exact = R0 - Vn * tf
    @test abs(R_mean - R_exact) / R_exact < 0.02
end

@testset "Curve-shortening flow on circle" begin
    R0   = 1.0
    beta = 0.1
    tf   = 2.0
    N    = 64
    dt   = 1e-3

    mesh = make_circle(R0, N)
    eq   = FrontEquation(;
        terms=(CurvatureMotionTerm(beta),),
        front=mesh, t=0.0, integrator=RK2(cfl=0.8)
    )
    integrate!(eq, tf; dt=dt)

    pts    = eq.state.mesh.points
    center = sum(pts) / length(pts)
    radii  = [norm(p - center) for p in pts]
    R_mean = sum(radii) / length(radii)

    R_exact = sqrt(R0^2 - 2*beta*tf)
    @test abs(R_mean - R_exact) / R_exact < 0.02
end

@testset "Equal-arclength redistribution" begin
    N = 32
    R = 1.0
    angles_bad = vcat(
        LinRange(0, π/2, 24),
        LinRange(π/2, 2π, 10)[2:end]
    )
    pts  = [SVector{2,Float64}(R*cos(θ), R*sin(θ)) for θ in angles_bad[1:N]]
    edges = [SVector{2,Int}(i, mod1(i+1, N)) for i in 1:N]
    mesh = CurveMesh(pts, edges)

    state = FrontState(mesh)
    A_before = front_enclosed_measure(state)
    spread_before = edge_length_spread(state)

    r = CurveEqualArcRedistributor()
    redistribute!(state, r)

    A_after  = front_enclosed_measure(state)
    spread_after = edge_length_spread(state)

    @test length(state.mesh.points) == N
    @test abs(A_after - A_before) / A_before < 0.02
    @test spread_after < spread_before
end

@testset "Field transfer – constant field preserved" begin
    N = 32
    R = 1.0
    mesh_old = make_circle(R, N)
    mesh_new = make_circle(R * 1.0001, N)

    oldvals = fill(3.14, N)
    newvals = similar(oldvals)
    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:piecewise_linear)
    @test all(abs.(newvals .- 3.14) .< 1e-10)
end

@testset "Diagnostics" begin
    mesh  = make_circle(1.0, 64)
    state = FrontState(mesh)

    h_min = min_edge_length(state)
    h_max = max_edge_length(state)
    h_mea = mean_edge_length(state)
    @test h_min > 0
    @test h_max >= h_min
    @test h_mea >= h_min && h_mea <= h_max

    c = front_centroid(state)
    @test norm(c) < 1e-10

    A = front_enclosed_measure(state)
    @test abs(A - π) / π < 0.01

    @test check_front_validity(state; warn=false)
end

@testset "Sphere translation" begin
    mesh0 = make_sphere_benchmark_surface(
        center = SVector(0.0, 0.0, 0.0),
        R      = 1.0,
        refinement = 2,
    )
    u  = SVector(1.0, 0.5, -0.3)
    tf = 0.5
    eq = FrontEquation(;
        terms=(AdvectionTerm(u),),
        front=mesh0, t=0.0, integrator=RK2()
    )
    integrate!(eq, tf; dt=0.05)

    pts    = eq.state.mesh.points
    center = sum(pts) / length(pts)
    exact  = SVector(1.0, 0.5, -0.3) .* tf
    @test norm(center - exact) < 1e-6
end

# ─────────────────────────────────────────────────────────────────────────────
# v0.2 test files
# ─────────────────────────────────────────────────────────────────────────────

include("test_frontfield.jl")
include("test_frontstate.jl")
include("test_advection_terms.jl")
include("test_normal_motion.jl")
include("test_curvature_motion.jl")
include("test_curve_redistribution.jl")
include("test_surface_redistribution.jl")
include("test_surface_quality_metrics.jl")
include("test_transfer_curve.jl")
include("test_transfer_surface.jl")
include("test_benchmark_rigid_2d.jl")
include("test_benchmark_open_curve_2d.jl")
include("test_benchmark_rigid_3d.jl")
include("test_benchmark_zalesak_2d.jl")
include("test_benchmark_zalesak_sphere_3d.jl")
include("test_benchmark_vortex_2d.jl")
include("test_benchmark_serpentine_2d.jl")
include("test_benchmark_enright_3d.jl")

if Base.find_package("Makie") !== nothing && Base.find_package("CairoMakie") !== nothing
    include("test_makie_ext.jl")
else
    @info "Skipping Makie extension tests (Makie/CairoMakie not available in test env)."
end

println("All FrontTrackingMethods tests passed.")

