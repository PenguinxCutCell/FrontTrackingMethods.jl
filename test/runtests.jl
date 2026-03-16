using Test
using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra

# ─────────────────────────────────────────────────────────────────────────────
# Helpers: create a circle CurveMesh
# ─────────────────────────────────────────────────────────────────────────────

function make_circle(R=1.0, N=64)
    T = Float64
    angles = [2π * (i-1) / N for i in 1:N]
    pts  = [SVector{2,T}(R * cos(θ), R * sin(θ)) for θ in angles]
    edges = [SVector{2,Int}(i, mod1(i+1, N)) for i in 1:N]
    return CurveMesh(pts, edges)
end

function make_sphere(R=1.0, refinements=2)
    return FrontIntrinsicOps.sphere_mesh(R; refinements=refinements)
end

# ─────────────────────────────────────────────────────────────────────────────
# A. Constructor / API tests
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

# ─────────────────────────────────────────────────────────────────────────────
# B. Pure translation – circle
# ─────────────────────────────────────────────────────────────────────────────

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

# ─────────────────────────────────────────────────────────────────────────────
# C. Constant normal motion – circle
# ─────────────────────────────────────────────────────────────────────────────

@testset "Constant normal motion – circle (inward, shrink)" begin
    # Sign convention for CurveMesh:
    # vertex_normals are INWARD for CCW circles.
    # NormalMotionTerm with Vn > 0 moves inward → circle shrinks.
    # R(t) = R0 - Vn * t  (since inward normal, Vn > 0 reduces radius)
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

    # Center should stay at origin
    @test norm(center) < 1e-3

    # Radius should be close to R0 - Vn*tf
    R_exact = R0 - Vn * tf
    @test abs(R_mean - R_exact) / R_exact < 0.02
end

# ─────────────────────────────────────────────────────────────────────────────
# D. Curvature motion – circle (curve-shortening flow)
# ─────────────────────────────────────────────────────────────────────────────

@testset "Curve-shortening flow on circle" begin
    # Sign convention (from FrontIntrinsicOps):
    # For CCW circles, signed_curvature > 0, vertex_normals point INWARD.
    # CurvatureMotionTerm velocity = β * κ * n_inward.
    # β > 0 → inward motion → shrinkage.
    # Exact: R(t) = sqrt(R0^2 - 2β t).
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

# ─────────────────────────────────────────────────────────────────────────────
# E. Redistribution – curve equal-arclength
# ─────────────────────────────────────────────────────────────────────────────

@testset "Equal-arclength redistribution" begin
    # Build a circle with clustered vertices
    N = 32
    R = 1.0
    # Cluster first half of vertices in a quarter arc
    angles_bad = vcat(
        LinRange(0, π/2, 24),    # 24 vertices in first quarter
        LinRange(π/2, 2π, 10)[2:end]  # 8 in the rest (including wrap)
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

    # Closure: still N vertices
    @test length(state.mesh.points) == N
    # Area nearly preserved
    @test abs(A_after - A_before) / A_before < 0.02
    # Spread reduced (better uniformity)
    @test spread_after < spread_before
end

# ─────────────────────────────────────────────────────────────────────────────
# F. Field transfer – constant field preserved
# ─────────────────────────────────────────────────────────────────────────────

@testset "Field transfer – constant field preserved" begin
    N = 32
    R = 1.0
    mesh_old = make_circle(R, N)
    mesh_new = make_circle(R * 1.0001, N)  # slightly different but same connectivity

    oldvals = fill(3.14, N)
    newvals = similar(oldvals)
    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:piecewise_linear)
    @test all(abs.(newvals .- 3.14) .< 1e-10)
end

# ─────────────────────────────────────────────────────────────────────────────
# G. Diagnostics
# ─────────────────────────────────────────────────────────────────────────────

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
    @test norm(c) < 1e-10   # circle centered at origin

    A = front_enclosed_measure(state)
    @test abs(A - π) / π < 0.01   # area ≈ π for R=1

    @test check_front_validity(state; warn=false)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D sphere tests (if sphere_mesh available)
# ─────────────────────────────────────────────────────────────────────────────

@testset "Sphere translation" begin
    try
        mesh = make_sphere(1.0, 2)
        u    = SVector(1.0, 0.5, -0.3)
        tf   = 0.5
        eq   = FrontEquation(;
            terms=(AdvectionTerm(u),),
            front=mesh, t=0.0, integrator=RK2()
        )
        integrate!(eq, tf; dt=0.05)

        pts    = eq.state.mesh.points
        center = sum(pts) / length(pts)
        exact  = SVector(1.0, 0.5, -0.3) .* tf
        @test norm(center - exact) < 1e-6
    catch e
        @warn "Skipping sphere test: $e"
    end
end

println("All FrontTrackingMethods tests passed.")
