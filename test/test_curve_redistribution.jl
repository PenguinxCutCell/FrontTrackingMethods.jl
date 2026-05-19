# test_curve_redistribution.jl – Tests for curve redistribution.

using Test
using StaticArrays
using LinearAlgebra
using Statistics
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Equal-arclength redistribution" begin
    # After redistribution, all edges should have equal length
    mesh = make_circle_curve(R=1.0, N=64)

    # Manually perturb vertices to create unequal spacing
    pts = copy(mesh.points)
    for i in eachindex(pts)
        θ = 2π * (i-1) / length(pts)
        r = 1.0 + 0.2 * sin(3θ)   # perturb radially
        pts[i] = SVector(r*cos(θ), r*sin(θ))
    end
    perturbed = CurveMesh{Float64}(pts, mesh.edges)
    state = FrontState(perturbed)

    r = CurveEqualArcRedistributor(; every=1)
    redistribute!(state, r)

    geom  = state.geom
    ls    = geom.edge_lengths
    @test maximum(ls) / minimum(ls) < 1.01   # nearly equal spacing
    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)
end

@testset "Adaptive curve remesher – corner protection" begin
    # Zalesak disk should preserve slot corners after remeshing
    mesh  = make_zalesak_disk_curve(N_arc=128, N_slot=32)
    state = FrontState(mesh)

    r = AdaptiveCurveRemesher(; iterations=5, protect_corners=true)
    redistribute!(state, r)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)

    # Area should be approximately preserved (within 5%)
    A_before = curve_area(mesh)
    A_after  = curve_area(state.mesh)
    @test abs(A_after - A_before) / abs(A_before) < 0.05
end

@testset "AdaptiveCurveRemesher – circle (no corners)" begin
    mesh  = make_circle_curve(R=1.0, N=64)
    state = FrontState(mesh)

    r = AdaptiveCurveRemesher(; iterations=10, protect_corners=false)
    redistribute!(state, r)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)

    # Area approximately preserved
    A_before = curve_area(mesh)
    A_after  = curve_area(state.mesh)
    @test abs(A_after - A_before) / abs(A_before) < 0.01
end

@testset "PoissonTangentialRedistributor – clustered curve spacing improves" begin
    mesh = make_circle_curve(R=1.0, N=64)
    pts = copy(mesh.points)
    for i in eachindex(pts)
        θ = 2π * (i - 1) / length(pts)
        pts[i] = SVector((1.0 + 0.15 * sin(3θ)) * cos(θ),
                         (1.0 + 0.15 * sin(3θ)) * sin(θ))
    end
    state = FrontState(CurveMesh{Float64}(pts, mesh.edges))
    A0 = front_enclosed_measure(state)
    spread0 = edge_length_spread(state)

    r = PoissonTangentialRedistributor(; iterations=5, pseudo_dt=0.1, omega=1.0)
    redistribute!(state, r)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)
    @test edge_length_spread(state) < spread0
    @test abs(front_enclosed_measure(state) - A0) / abs(A0) < 0.03
end

@testset "PoissonTangentialRedistributor – open curve rejected" begin
    pts = [SVector{2,Float64}(0.0, 0.0), SVector{2,Float64}(1.0, 0.0), SVector{2,Float64}(2.0, 0.0)]
    edges = [SVector{2,Int}(1, 2), SVector{2,Int}(2, 3)]
    state = FrontState(CurveMesh{Float64}(pts, edges))
    @test_throws ErrorException redistribute!(state, PoissonTangentialRedistributor())
end

@testset "Curve diagnostics – arc correction passthrough" begin
    coarse = make_circle_curve(R=1.0, N=12)
    fine = make_circle_curve(R=1.0, N=512)
    state = FrontState(coarse)

    @test abs(front_measure(state; correction=:arc) - 2π) <
          abs(front_measure(state) - 2π)
    @test abs(front_enclosed_measure(state; correction=:arc) - π) <
          abs(front_enclosed_measure(state) - π)
    @test abs(perimeter_error(coarse, fine; correction=:arc)) <
          abs(perimeter_error(coarse, fine))
    @test abs(relative_area_error(coarse, fine; correction=:arc)) <
          abs(relative_area_error(coarse, fine))
end

@testset "PoissonTangentialRedistributor – every interval honored" begin
    mesh = make_circle_curve(R=1.0, N=48)
    pts = copy(mesh.points)
    for i in eachindex(pts)
        θ = 2π * (i - 1) / length(pts)
        pts[i] = SVector((1.0 + 0.10 * sin(4θ)) * cos(θ),
                         (1.0 + 0.10 * sin(4θ)) * sin(θ))
    end
    state = FrontState(CurveMesh{Float64}(pts, mesh.edges))
    spread0 = edge_length_spread(state)

    eq = FrontEquation(;
        terms = NormalMotionTerm(0.0),
        front = state,
        integrator = ForwardEuler(),
        redistributor = PoissonTangentialRedistributor(; iterations=3, pseudo_dt=0.1, every=2),
    )
    integrate!(eq, 0.01; dt=0.01)
    @test edge_length_spread(current_state(eq)) ≈ spread0

    integrate!(eq, 0.02; dt=0.01)
    @test edge_length_spread(current_state(eq)) < spread0
end

@testset "remesh! alias" begin
    mesh  = make_circle_curve(R=1.0, N=32)
    state = FrontState(mesh)
    r     = CurveEqualArcRedistributor()
    @test remesh!(state, r) === state
end
