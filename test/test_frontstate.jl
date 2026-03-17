# test_frontstate.jl – Tests for FrontState API.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "FrontState – construction and fields" begin
    mesh  = make_circle_curve(R=1.0, N=32)
    state = FrontState(mesh; t=0.5)

    @test current_time(state) ≈ 0.5
    @test state.mesh === mesh

    # Geometry is available
    @test state.geom isa FrontIntrinsicOps.CurveGeometry

    # Add a field
    vals = ones(Float64, length(mesh.points))
    ff   = FrontField(vals, mesh, :vertex)
    add_field!(state, :phi, ff)

    ff2 = get_field(state, :phi)
    @test ff2.values == vals

    # Check that location lookup works
    @test location(ff) === :vertex
    @test mesh === FrontTrackingMethods.mesh(ff)
end

@testset "FrontState – geometry refresh" begin
    mesh  = make_circle_curve(R=1.0, N=64)
    state = FrontState(mesh)

    # Move all vertices outward
    new_pts = [p * 2.0 for p in mesh.points]
    set_vertex_coordinates!(state, new_pts)
    refresh_geometry!(state)

    # Edge lengths should have doubled
    geom_orig = FrontIntrinsicOps.compute_geometry(mesh)
    geom_new  = state.geom

    @test sum(geom_new.edge_lengths) ≈ 2 * sum(geom_orig.edge_lengths) rtol=1e-10
end
