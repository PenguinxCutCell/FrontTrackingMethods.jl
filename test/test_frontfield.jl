# test_frontfield.jl – Tests for FrontField construction and API.

using Test
using StaticArrays
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "FrontField – Constructor / API" begin
    mesh = make_circle_curve(R=1.0, N=32)
    N    = length(mesh.points)

    # Scalar field construction
    vals = ones(Float64, N)
    ff   = FrontField(vals, mesh, :vertex)
    @test ff.location === :vertex
    @test length(ff.values) == N

    # State construction
    state = FrontState(mesh; t=0.0)
    @test current_time(state) == 0.0
    @test length(state.mesh.points) == N

    # Add and retrieve field
    add_field!(state, :phi, ff)
    ff2 = get_field(state, :phi)
    @test ff2.values == vals

    # Missing field throws
    @test_throws ErrorException get_field(state, :missing_field)

    # Vertex coordinates
    pts = vertex_coordinates(state)
    @test length(pts) == N

    # set_vertex_coordinates!
    new_pts = [p + SVector(0.1, 0.0) for p in pts]
    set_vertex_coordinates!(state, new_pts)
    pts2 = vertex_coordinates(state)
    @test pts2[1][1] ≈ pts[1][1] + 0.1

    # FrontEquation constructors
    u  = rigid_translation_velocity(SVector(1.0, 0.0))
    eq = FrontEquation(; terms=AdvectionTerm(u), front=mesh)
    @test current_time(eq) == 0.0
    @test current_state(eq) === eq.state
end
