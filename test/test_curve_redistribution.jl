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

@testset "remesh! alias" begin
    mesh  = make_circle_curve(R=1.0, N=32)
    state = FrontState(mesh)
    r     = CurveEqualArcRedistributor()
    @test remesh!(state, r) === state
end
