# test_transfer_surface.jl – Tests for field transfer on surfaces.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

include("test_utils.jl")

@testset "Field transfer – constant field preserved (surface, barycentric)" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    mesh_new = make_sphere_surface(R=1.0, refinement=2)

    oldvals = ones(Float64, length(mesh_old.points))
    newvals = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:barycentric)
    @test all(v -> isapprox(v, 1.0; atol=1e-10), newvals)
end

@testset "Field transfer – constant field preserved (surface, nearest_vertex)" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    mesh_new = make_sphere_surface(R=1.0, refinement=2)

    oldvals = fill(42.0, length(mesh_old.points))
    newvals = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:nearest_vertex)
    @test all(v -> isapprox(v, 42.0; atol=1e-10), newvals)
end

@testset "Field transfer – smooth z field (surface, barycentric vs nearest)" begin
    # Transfer f(x,y,z) = z from a coarse sphere to a finer sphere
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    mesh_new = make_sphere_surface(R=1.0, refinement=3)

    oldvals = [p[3] for p in mesh_old.points]
    newvals_bary = zeros(Float64, length(mesh_new.points))
    newvals_near = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals_bary, mesh_old, oldvals, mesh_new; method=:barycentric)
    transfer_vertex_field!(newvals_near, mesh_old, oldvals, mesh_new; method=:nearest_vertex)

    # Both methods should roughly match the analytic value f=z on unit sphere
    expected = [p[3] for p in mesh_new.points]

    err_bary = sqrt(sum((newvals_bary .- expected).^2) / length(expected))
    err_near = sqrt(sum((newvals_near .- expected).^2) / length(expected))

    # Barycentric should be more accurate than nearest-vertex
    @test err_bary <= err_near + 1e-8
    @test err_bary < 0.1   # reasonable accuracy
end
