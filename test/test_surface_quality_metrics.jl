# test_surface_quality_metrics.jl – Surface quality diagnostic tests.

using Test
using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays

@testset "Surface quality diagnostics" begin
    mesh = make_sphere_benchmark_surface(center=SVector(0.0, 0.0, 0.0), R=1.0, refinement=2)
    geom = compute_geometry(mesh)

    edge = surface_edge_length_stats(mesh)
    @test edge.min > 0
    @test edge.max >= edge.min
    @test edge.mean >= edge.min
    @test edge.mean <= edge.max
    @test edge.std >= 0
    @test edge.ratio >= 1

    area = surface_triangle_area_stats(mesh)
    @test area.min > 0
    @test area.max >= area.min
    @test area.mean >= area.min
    @test area.mean <= area.max

    ang = surface_triangle_angle_stats(mesh)
    @test ang.min_angle > 0
    @test ang.max_angle < π
    @test ang.mean_min_angle > 0

    asp = surface_aspect_ratio_stats(mesh)
    @test asp.min >= 1.0 - 1e-12
    @test asp.max >= asp.min
    @test asp.mean >= asp.min

    degen = surface_degenerate_fraction(mesh; atol=1e-12)
    @test degen == 0.0

    q = surface_quality_summary(mesh)
    @test q.edge.ratio == edge.ratio
    @test q.area.mean == area.mean
    @test q.degenerate_fraction == degen

    nc = surface_normal_consistency(mesh, geom)
    @test nc.inconsistent_fraction <= 0.05
    @test nc.mean_alignment > 0.5
    @test nc.min_alignment > -0.5
end
