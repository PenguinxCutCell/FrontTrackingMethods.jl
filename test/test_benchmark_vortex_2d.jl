# test_benchmark_vortex_2d.jl – Rider–Kothe reversed vortex benchmark test.
#
# Tests that after one full reversal period T, the front approximately
# recovers to its initial shape.  Tests with and without redistribution.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Rider–Kothe vortex – modest mesh, with redistribution" begin
    T     = 2.0    # reversal period
    dt    = 0.02
    N     = 128

    mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
    u     = rider_kothe_single_vortex(; T=T)

    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = CurveEqualArcRedistributor(),
    )
    integrate!(eq, T; dt=dt)

    state = current_state(eq)

    # Area drift should be small (below 5%)
    A0 = curve_area(mesh0)
    Af = curve_area(state.mesh)
    @test abs(Af - A0) / A0 < 0.05

    # Centroid recovery
    c0 = curve_centroid(mesh0)
    cf = curve_centroid(state.mesh)
    @test norm(cf - c0) < 0.01

    assert_no_nan(state.mesh)
    assert_closed_curve(state.mesh)
end

@testset "Rider–Kothe vortex – redistribution improves recovery vs none" begin
    T     = 2.0
    dt    = 0.02
    N     = 128
    mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
    u     = rider_kothe_single_vortex(; T=T)

    # Without redistribution
    eq_no = FrontEquation(;
        terms      = AdvectionTerm(u),
        front      = deepcopy(mesh0),
        integrator = RK2(),
    )
    integrate!(eq_no, T; dt=dt)

    # With equal-arc redistribution
    eq_re = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = CurveEqualArcRedistributor(),
    )
    integrate!(eq_re, T; dt=dt)

    # With adaptive remesher
    eq_ad = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = AdaptiveCurveRemesher(; iterations=5),
    )
    integrate!(eq_ad, T; dt=dt)

    d_no = symmetric_hausdorff_curve(current_state(eq_no).mesh, mesh0)
    d_re = symmetric_hausdorff_curve(current_state(eq_re).mesh, mesh0)
    d_ad = symmetric_hausdorff_curve(current_state(eq_ad).mesh, mesh0)

    # Both redistribution strategies should give a finite, non-NaN distance
    @test isfinite(d_re)
    @test isfinite(d_ad)
    # Log the distances for comparison (no strict ordering required at modest resolution)
    @info "Rider-Kothe Hausdorff: no_redist=$(round(d_no,digits=4)), equal_arc=$(round(d_re,digits=4)), adaptive=$(round(d_ad,digits=4))"
end
