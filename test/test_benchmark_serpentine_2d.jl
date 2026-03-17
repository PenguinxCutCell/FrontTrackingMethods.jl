# test_benchmark_serpentine_2d.jl – Serpentine deformation benchmark test.
#
# Smoke test: verify the front completes Tmax without catastrophic collapse.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Serpentine deformation – smoke test (no collapse)" begin
    Tmax = 1.0   # modest time (full Tmax=3 reserved for long convergence scripts)
    dt   = 0.02
    N    = 128

    mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
    u     = serpentine_2d(; Tmax=Tmax)

    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = AdaptiveCurveRemesher(; iterations=5, protect_corners=false),
    )
    integrate!(eq, Tmax; dt=dt)

    state = current_state(eq)

    # Area drift within 10% (generous for severe shear)
    A0 = curve_area(mesh0)
    Af = curve_area(state.mesh)
    @test abs(Af - A0) / abs(A0) < 0.10

    # No NaN
    assert_no_nan(state.mesh)
    assert_closed_curve(state.mesh)

    # Min edge length stays positive
    stats = edge_length_stats(state.mesh)
    @test stats.min > 0
end
