using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology detection 2D – merge/split candidates" begin
    handler = LocalCartesianTopologyHandler(
        d_merge=2.0,
        d_split=1.2,
        reconstruct_scope=:local_patch,
    )

    # merge candidate: near circles
    near_fronts = make_two_circles_merge_setup(
        c1=SVector(0.45, 0.5),
        c2=SVector(0.55, 0.5),
        R=0.08,
        N=96,
    )
    near_state = MultiFrontState(near_fronts)
    near_candidates = find_topology_candidates(near_state, handler)
    @test any(c.event_type == :merge for c in near_candidates)

    # no merge candidate: far circles
    far_fronts = make_two_circles_merge_setup(
        c1=SVector(0.2, 0.5),
        c2=SVector(0.8, 0.5),
        R=0.08,
        N=96,
    )
    far_state = MultiFrontState(far_fronts)
    far_candidates = find_topology_candidates(far_state, handler)
    @test !any(c.event_type == :merge for c in far_candidates)

    # split candidate: thin dumbbell
    dumbbell = make_dumbbell_curve_setup(neck_strength=0.82, N=320)
    split_state = MultiFrontState([dumbbell])
    split_candidates = find_topology_candidates(split_state, handler)
    @test any(c.event_type == :split for c in split_candidates)

    # not too early: thick peanut
    peanut = make_peanut_curve_setup(neck_strength=0.15, N=320)
    peanut_state = MultiFrontState([peanut])
    peanut_candidates = find_topology_candidates(peanut_state, handler)
    @test !any(c.event_type == :split for c in peanut_candidates)
end
