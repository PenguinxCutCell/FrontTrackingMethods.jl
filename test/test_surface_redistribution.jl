# test_surface_redistribution.jl – Tests for surface redistribution.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "SurfaceTangentialRedistributor" begin
    mesh  = make_sphere_surface(R=1.0, refinement=2)
    state = FrontState(mesh)

    r = SurfaceTangentialRedistributor(; iterations=3, strength=0.5)
    redistribute!(state, r)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)

    # Volume should be approximately preserved (within 5%)
    V_before = surface_volume(mesh)
    V_after  = surface_volume(state.mesh)
    @test abs(V_after - V_before) / abs(V_before) < 0.05
end

@testset "ExperimentalSurfaceRemesher" begin
    mesh  = make_sphere_surface(R=1.0, refinement=2)
    state = FrontState(mesh)

    r = ExperimentalSurfaceRemesher(; iterations=3, strength=0.3)
    redistribute!(state, r)

    @test is_closed(state.mesh)
    assert_no_nan(state.mesh)
end
