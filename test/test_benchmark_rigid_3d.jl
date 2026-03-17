# test_benchmark_rigid_3d.jl – Rigid-motion 3-D benchmark tests.
#
# Tests:
# 1. Rigid translation of a sphere
# 2. Rigid rotation of a sphere (one full period)

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

include("test_utils.jl")

@testset "Rigid translation – sphere benchmark" begin
    mesh0 = make_sphere_benchmark_surface(center=SVector(0.35, 0.35, 0.35), R=0.15,
                                          refinement=2)
    u  = rigid_translation_velocity(SVector(0.1, 0.0, 0.0))
    tf = 1.0
    dt = 0.05

    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, tf; dt=dt)

    state = current_state(eq)
    c0    = surface_centroid(mesh0)
    cf    = surface_centroid(state.mesh)

    @test cf[1] ≈ c0[1] + 0.1 * tf  atol=1e-8
    @test cf[2] ≈ c0[2]              atol=1e-8
    @test cf[3] ≈ c0[3]              atol=1e-8

    V0 = surface_volume(mesh0)
    Vf = surface_volume(state.mesh)
    @test abs(Vf - V0) / V0 < 1e-6

    assert_no_nan(state.mesh)
    assert_surface_is_reasonable(state.mesh)
end

@testset "Rigid rotation – sphere benchmark (one full period)" begin
    center = SVector(0.5, 0.5, 0.5)
    mesh0  = make_sphere_benchmark_surface(center=SVector(0.5, 0.75, 0.5), R=0.15,
                                           refinement=2)
    omega  = 2π
    T_rev  = 1.0
    dt     = 0.05

    u  = rigid_rotation_3d(; center=center, axis=:z, omega=omega)
    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, T_rev; dt=dt)

    state = current_state(eq)

    # Centroid returns
    c0 = surface_centroid(mesh0)
    cf = surface_centroid(state.mesh)
    @test norm(cf - c0) < 0.05

    # Volume preserved
    V0 = surface_volume(mesh0)
    Vf = surface_volume(state.mesh)
    @test abs(Vf - V0) / V0 < 0.10

    # Approximate Hausdorff distance to initial
    d_hausdorff = symmetric_hausdorff_surface(state.mesh, mesh0)
    @test d_hausdorff < 0.05

    assert_no_nan(state.mesh)
end
