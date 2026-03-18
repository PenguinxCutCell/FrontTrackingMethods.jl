# test_surface_redistribution.jl – Tests for surface redistribution.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "SurfaceTangentialRedistributor" begin
    mesh  = make_sphere_surface(R=1.0, refinement=2)
    state = FrontState(mesh)
    q0 = surface_quality_summary(state.mesh)

    r = SurfaceTangentialRedistributor(; iterations=3, strength=0.5)
    redistribute!(state, r)
    q1 = surface_quality_summary(state.mesh)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)

    # Volume should be approximately preserved (within 5%)
    V_before = surface_volume(mesh)
    V_after  = surface_volume(state.mesh)
    @test abs(V_after - V_before) / abs(V_before) < 0.05

    @test q1.degenerate_fraction <= 0.01
    @test q1.angle.min_angle >= 0.8 * q0.angle.min_angle
end

@testset "ExperimentalSurfaceRemesher" begin
    mesh  = make_sphere_surface(R=1.0, refinement=2)
    state = FrontState(mesh)
    q0 = surface_quality_summary(state.mesh)
    V0 = surface_volume(mesh)

    r = ExperimentalSurfaceRemesher(; iterations=3, strength=0.25, volume_correction=true)
    redistribute!(state, r)
    q1 = surface_quality_summary(state.mesh)
    V1 = surface_volume(state.mesh)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)
    @test abs(V1 - V0) / abs(V0) < 0.05
    @test q1.degenerate_fraction <= 0.01
    @test q1.aspect.max <= 1.25 * q0.aspect.max
    @test q1.angle.min_angle >= 0.75 * q0.angle.min_angle
end
