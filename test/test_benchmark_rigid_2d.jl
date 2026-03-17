# test_benchmark_rigid_2d.jl – Rigid-motion 2-D benchmark tests.
#
# Tests:
# 1. Rigid translation of a circle (shape, area, centroid preservation)
# 2. Rigid rotation of a circle (shape, area, centroid preservation)

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

include("test_utils.jl")

@testset "Rigid translation – circle benchmark" begin
    mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=128)
    u     = rigid_translation_velocity(SVector(0.1, 0.05))
    tf    = 1.0
    dt    = 0.01

    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, tf; dt=dt)

    state = current_state(eq)
    c0    = curve_centroid(mesh0)
    cf    = curve_centroid(state.mesh)

    # Centroid moves by u * tf
    @test cf[1] ≈ c0[1] + 0.1 * tf  atol=1e-10
    @test cf[2] ≈ c0[2] + 0.05 * tf atol=1e-10

    # Area preserved
    A0 = curve_area(mesh0)
    Af = curve_area(state.mesh)
    @test abs(Af - A0) / A0 < 1e-6

    assert_no_nan(state.mesh)
    assert_closed_curve(state.mesh)
end

@testset "Rigid rotation – circle benchmark (one full period)" begin
    center = SVector(0.5, 0.5)
    mesh0  = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=128)
    omega  = 2π    # one revolution per unit time
    T_rev  = 1.0   # period
    dt     = 0.01

    u  = rigid_rotation_2d(; center=center, omega=omega)
    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2(),
                         redistributor=CurveEqualArcRedistributor())
    integrate!(eq, T_rev; dt=dt)

    state  = current_state(eq)

    # Centroid should have returned to original position
    c0 = curve_centroid(mesh0)
    cf = curve_centroid(state.mesh)
    @test norm(cf - c0) < 2e-3

    # Area preserved
    A0 = curve_area(mesh0)
    Af = curve_area(state.mesh)
    @test abs(Af - A0) / A0 < 2e-3

    # Shape recovery: Hausdorff distance to initial front should be small
    d_hausdorff = symmetric_hausdorff_curve(state.mesh, mesh0)
    @test d_hausdorff < 3e-3

    assert_no_nan(state.mesh)
end
