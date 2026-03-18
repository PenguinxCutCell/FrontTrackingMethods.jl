using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology merge 2D – end-to-end" begin
    fronts = make_two_circles_merge_setup(
        c1=SVector(0.47, 0.5),
        c2=SVector(0.53, 0.5),
        R=0.08,
        N=96,
    )
    state = MultiFrontState(fronts)

    A0 = front_enclosed_measure(state)

    handler = LocalCartesianTopologyHandler(
        d_merge=2.5,
        d_split=0.8,
        patch_h_factor=0.45,
        patch_margin_factor=2.0,
        reconstruct_scope=:whole_component,
        preserve_fields=true,
    )

    report = handle_topology_change!(state, handler)
    @test report.changed
    @test report.event_type == :merge
    @test ncomponents(state) == 1

    A1 = front_enclosed_measure(state)
    @test abs(A1 - A0) / max(abs(A0), eps(Float64)) < 0.40
end
