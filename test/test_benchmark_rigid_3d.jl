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

    @test cf[1] ≈ c0[1] + 0.1 * tf  atol=1e-9
    @test cf[2] ≈ c0[2]              atol=1e-9
    @test cf[3] ≈ c0[3]              atol=1e-9

    V0 = surface_volume(mesh0)
    Vf = surface_volume(state.mesh)
    @test abs(Vf - V0) / V0 < 5e-7

    A0 = FrontIntrinsicOps.measure(mesh0, FrontIntrinsicOps.compute_geometry(mesh0))
    Af = FrontIntrinsicOps.measure(state.mesh, FrontIntrinsicOps.compute_geometry(state.mesh))
    @test abs(Af - A0) / A0 < 2e-6

    q = surface_quality_summary(state.mesh; degenerate_atol=1e-14)
    @test q.angle.min_angle > deg2rad(8.0)
    @test q.degenerate_fraction == 0.0

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
    eq_none = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq_none, T_rev; dt=dt)
    state_none = current_state(eq_none)

    eq_remesh = FrontEquation(;
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=ExperimentalSurfaceRemesher(; iterations=2, strength=0.2, volume_correction=true),
    )
    integrate!(eq_remesh, T_rev; dt=dt)
    state_remesh = current_state(eq_remesh)

    # Centroid returns
    c0 = surface_centroid(mesh0)
    cf_none = surface_centroid(state_none.mesh)
    cf_remesh = surface_centroid(state_remesh.mesh)
    @test norm(cf_none - c0) < 0.03
    @test norm(cf_remesh - c0) < 0.04

    # Volume preserved
    V0 = surface_volume(mesh0)
    Vf_none = surface_volume(state_none.mesh)
    Vf_remesh = surface_volume(state_remesh.mesh)
    @test abs(Vf_none - V0) / V0 < 0.06
    @test abs(Vf_remesh - V0) / V0 < 0.08

    # Approximate Hausdorff distance to initial
    d_none = symmetric_hausdorff_surface(state_none.mesh, mesh0)
    d_remesh = symmetric_hausdorff_surface(state_remesh.mesh, mesh0)
    @test d_none < 0.05
    @test d_remesh <= 1.35 * d_none + 1e-8

    q_none = surface_quality_summary(state_none.mesh; degenerate_atol=1e-12)
    q_remesh = surface_quality_summary(state_remesh.mesh; degenerate_atol=1e-12)
    @test q_none.angle.min_angle > deg2rad(2.5)
    @test q_remesh.angle.min_angle > deg2rad(3.0)
    @test q_none.degenerate_fraction <= 0.01
    @test q_remesh.degenerate_fraction <= 0.01

    assert_no_nan(state_none.mesh)
    assert_no_nan(state_remesh.mesh)
end
