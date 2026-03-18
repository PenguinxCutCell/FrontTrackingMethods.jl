using Test
using StaticArrays
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Topology state – multi-component basics" begin
    fronts = make_two_circles_merge_setup(
        c1=SVector(0.42, 0.5),
        c2=SVector(0.58, 0.5),
        R=0.10,
        N=96,
    )
    state = MultiFrontState(fronts; t=0.25)

    @test current_time(state) ≈ 0.25
    @test ncomponents(state) == 2
    @test length(all_meshes(state)) == 2
    @test component_mesh(state, 1) isa FrontIntrinsicOps.CurveMesh
    @test component_geom(state, 2) isa FrontIntrinsicOps.CurveGeometry
    @test component_fields(state, 1) isa Dict{Symbol,Any}
    @test length(collect(eachcomponent(state))) == 2
    @test length(map_components(c -> length(c.mesh.points), state)) == 2

    # selective geometry refresh
    p0 = component_mesh(state, 1).points[1]
    c1 = component(state, 1)
    c1.mesh = CurveMesh([p * 1.01 for p in c1.mesh.points], c1.mesh.edges)
    refresh_geometry!(state; which=1)
    @test component_mesh(state, 1).points[1] != p0

    # compatibility conversions
    single = FrontState(fronts[1]; t=0.5)
    ms = MultiFrontState(single)
    @test ncomponents(ms) == 1
    st_back = FrontState(ms)
    @test st_back.mesh isa FrontIntrinsicOps.CurveMesh
    @test current_time(st_back) ≈ current_time(ms)
end
