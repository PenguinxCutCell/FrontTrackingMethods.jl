using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology patch 2D – extraction and rasterization" begin
    handler = LocalCartesianTopologyHandler(d_merge=2.0, patch_h_factor=0.5)
    fronts = make_two_circles_merge_setup(
        c1=SVector(0.45, 0.5),
        c2=SVector(0.55, 0.5),
        R=0.08,
        N=96,
    )
    state = MultiFrontState(fronts)

    candidates = find_topology_candidates(state, handler)
    @test !isempty(candidates)

    cand = select_topology_candidate(candidates)
    patch = extract_event_patch(state, cand, handler)

    @test length(patch.component_ids) >= 1
    @test length(patch.patch_grid.dims) == 2
    nx, ny = patch.patch_grid.dims
    @test nx >= 2 && ny >= 2

    χ = zeros(Float64, patch.patch_grid.dims)
    rasterize_indicator!(χ, patch, patch.local_mesh_data)
    occ = sum(χ .> 0.5)
    @test occ > 0
    @test occ < length(χ)
end
