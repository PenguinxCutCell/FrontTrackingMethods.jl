# test_transfer_curve.jl – Tests for field transfer on curves.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Field transfer – constant field preserved (curve)" begin
    mesh_old = make_circle_curve(R=1.0, N=32)
    mesh_new = make_circle_curve(R=1.0, N=48)   # different resolution

    oldvals = ones(Float64, length(mesh_old.points))
    newvals = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:piecewise_linear)
    @test all(v -> isapprox(v, 1.0; atol=1e-10), newvals)
end

@testset "Field transfer – nearest-vertex fallback (curve)" begin
    mesh_old = make_circle_curve(R=1.0, N=32)
    mesh_new = make_circle_curve(R=1.0, N=16)

    oldvals = Float64.(1:length(mesh_old.points))
    newvals = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:nearest_vertex)
    @test all(v -> v >= 1.0 && v <= length(mesh_old.points), newvals)
end

@testset "Field transfer – piecewise linear smooth field (curve)" begin
    # Transfer sin(θ) from a 64-vertex circle to a 32-vertex circle
    N_old = 64
    N_new = 32
    mesh_old = make_circle_curve(R=1.0, N=N_old)
    mesh_new = make_circle_curve(R=1.0, N=N_new)

    oldvals = [sin(2π * k / N_old) for k in 0:N_old-1]
    newvals = zeros(Float64, N_new)

    transfer_vertex_field!(newvals, mesh_old, oldvals, mesh_new; method=:piecewise_linear)

    # Compare to analytic values at new vertex positions
    for (i, p) in enumerate(mesh_new.points)
        θ = atan(p[2], p[1])
        expected = sin(θ)
        @test isapprox(newvals[i], expected; atol=0.05)
    end
end
