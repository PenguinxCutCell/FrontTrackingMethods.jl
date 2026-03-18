using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology detection 3D – merge/split candidates" begin
    handler = LocalCartesianTopologyHandler(
        d_merge=2.2,
        d_split=0.1,
        reconstruct_scope=:whole_component,
    )

    merge_fronts = make_two_spheres_merge_setup(
        c1=SVector(0.44, 0.5, 0.5),
        c2=SVector(0.56, 0.5, 0.5),
        R=0.10,
        refinement=1,
    )
    merge_state = MultiFrontState(merge_fronts)
    merge_candidates = find_topology_candidates(merge_state, handler)
    @test any(c.event_type == :merge for c in merge_candidates)

    dumbbell = make_dumbbell_surface_setup(neck_strength=0.9, refinement=1)
    split_state = MultiFrontState([dumbbell])
    split_candidates = find_topology_candidates(split_state, handler)
    @test any(c.event_type == :split for c in split_candidates)
end
