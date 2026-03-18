using Test
using StaticArrays
using FrontTrackingMethods

@testset "Topology reconstruction 2D – loop extraction" begin
    # Synthetic two-blob occupancy should reconstruct to at least two loops.
    grid = make_patch_grid((SVector(0.0, 0.0), SVector(1.0, 1.0)), 0.05)
    patch = EventPatch(Int[], nothing, nothing, Dict{Symbol,Any}(), Any[], grid)

    χ = zeros(Float64, grid.dims)
    nx, ny = grid.dims
    for j in 1:ny, i in 1:nx
        c = grid.centers[i, j]
        in1 = (c[1] - 0.3)^2 + (c[2] - 0.5)^2 < 0.12^2
        in2 = (c[1] - 0.7)^2 + (c[2] - 0.5)^2 < 0.12^2
        χ[i, j] = (in1 || in2) ? 1.0 : 0.0
    end

    curves = reconstruct_curve_from_patch(χ, patch)
    @test length(curves) >= 2

    for curve in curves
        @test length(curve.points) >= 4
        @test !any(any(!isfinite, p) for p in curve.points)
    end
end
