using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology split 2D – end-to-end" begin
    curve = make_dumbbell_curve_setup(neck_strength=0.92, N=500)
    state = MultiFrontState([curve])

    A0 = front_enclosed_measure(state)

    handler = LocalCartesianTopologyHandler(
        d_split=1.0,
        d_merge=2.0,
        patch_h_factor=0.6,
        patch_margin_factor=2.0,
        reconstruct_scope=:whole_component,
        preserve_fields=true,
        max_components_created=8,
    )

    report = handle_topology_change!(state, handler)
    @test report.changed
    @test report.event_type == :split
    @test ncomponents(state) >= 1

    A1 = front_enclosed_measure(state)
    @test abs(A1 - A0) / max(abs(A0), eps(Float64)) < 0.60

    for comp in eachcomponent(state)
        @test !any(any(!isfinite, p) for p in comp.mesh.points)
        @test all(comp.geom.edge_lengths .> 0)
    end
end
