using Test
using FrontTrackingMethods

@testset "Topology split 3D – prototype smoke" begin
    surface = make_dumbbell_surface_setup(neck_strength=0.9, refinement=1)
    state = MultiFrontState([surface])
    V0 = front_enclosed_measure(state)

    handler = LocalCartesianTopologyHandler(
        d_split=0.1,
        d_merge=2.0,
        patch_h_factor=1.0,
        reconstruct_scope=:whole_component,
        max_components_created=8,
    )

    report = handle_topology_change!(state, handler)
    @test report.changed
    @test report.event_type == :split
    @test ncomponents(state) >= 1

    V1 = front_enclosed_measure(state)
    @test isfinite(V1)
    @test abs(V1 - V0) / max(abs(V0), eps(Float64)) < 0.95

    for comp in eachcomponent(state)
        @test !any(any(!isfinite, p) for p in comp.mesh.points)
        @test !isempty(comp.mesh.faces)
    end
end
