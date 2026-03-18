using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology merge 3D – prototype smoke" begin
    fronts = make_two_spheres_merge_setup(
        c1=SVector(0.44, 0.5, 0.5),
        c2=SVector(0.56, 0.5, 0.5),
        R=0.10,
        refinement=1,
    )
    state = MultiFrontState(fronts)
    V0 = front_enclosed_measure(state)

    handler = LocalCartesianTopologyHandler(
        d_merge=2.2,
        d_split=0.2,
        patch_h_factor=0.9,
        reconstruct_scope=:whole_component,
        max_components_created=8,
    )

    report = handle_topology_change!(state, handler)
    @test report.changed
    @test report.event_type == :merge
    @test ncomponents(state) == 1

    V1 = front_enclosed_measure(state)
    @test abs(V1 - V0) / max(abs(V0), eps(Float64)) < 0.85

    comp = component(state, 1)
    @test !any(any(!isfinite, p) for p in comp.mesh.points)
    @test !isempty(comp.mesh.faces)
end
