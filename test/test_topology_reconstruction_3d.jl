using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology reconstruction 3D – voxel prototype" begin
    grid = make_patch_grid((SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0)), 0.1)
    patch = EventPatch(Int[], nothing, nothing, Dict{Symbol,Any}(), Any[], grid)

    χ = zeros(Float64, grid.dims)
    nx, ny, nz = grid.dims
    for k in 1:nz, j in 1:ny, i in 1:nx
        c = grid.centers[i, j, k]
        b1 = abs(c[1] - 0.3) < 0.12 && abs(c[2] - 0.5) < 0.12 && abs(c[3] - 0.5) < 0.12
        b2 = abs(c[1] - 0.7) < 0.12 && abs(c[2] - 0.5) < 0.12 && abs(c[3] - 0.5) < 0.12
        χ[i, j, k] = (b1 || b2) ? 1.0 : 0.0
    end

    surfaces = reconstruct_surface_from_patch(χ, patch)
    @test length(surfaces) >= 2

    for surf in surfaces
        @test !isempty(surf.faces)
        @test !any(any(!isfinite, p) for p in surf.points)
    end
end
