using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology field transfer 2D – constant scalar" begin
    fronts = make_two_circles_merge_setup(
        c1=SVector(0.47, 0.5),
        c2=SVector(0.53, 0.5),
        R=0.08,
        N=96,
    )
    state = MultiFrontState(fronts)

    for comp in eachcomponent(state)
        vals = fill(3.25, length(comp.mesh.points))
        comp.fields[:phi] = FrontField(vals, comp.mesh, :vertex)
    end

    handler = LocalCartesianTopologyHandler(
        d_merge=2.5,
        d_split=0.8,
        patch_h_factor=0.45,
        reconstruct_scope=:whole_component,
        preserve_fields=true,
    )

    report = handle_topology_change!(state, handler)
    @test report.changed
    @test ncomponents(state) == 1

    comp = component(state, 1)
    @test haskey(comp.fields, :phi)
    fld = comp.fields[:phi]
    @test fld isa FrontField
    @test maximum(abs.(fld.values .- 3.25)) < 1e-10
end
