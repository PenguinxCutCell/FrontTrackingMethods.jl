# test_benchmark_open_curve_2d.jl – Open-curve 2-D benchmark tests.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Open arc benchmark constructor" begin
    center = SVector(0.5, 0.5)
    R = 0.2
    θ_start = π / 6
    θ_end = 5π / 6
    N = 65

    mesh = make_open_arc_benchmark_curve(; center=center, R=R, θ_start=θ_start, θ_end=θ_end, N=N)
    state = FrontState(mesh)

    @test !is_closed(mesh)
    @test length(mesh.points) == N
    @test length(mesh.edges) == N - 1
    @test curve_vertex_order(mesh) == collect(1:N)
    @test front_measure(state) ≈ R * abs(θ_end - θ_start) rtol=5e-4
    assert_no_nan(mesh)
end

@testset "Rigid translation – open arc benchmark" begin
    mesh0 = make_open_arc_benchmark_curve(; center=SVector(0.5, 0.5), R=0.2, N=96)
    shift = SVector(0.1, -0.05)
    u     = rigid_translation_velocity(shift)
    tf    = 0.8
    dt    = 0.01

    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, tf; dt=dt)

    state = current_state(eq)
    exact_shift = shift * tf

    @test !is_closed(state.mesh)
    @test length(state.mesh.edges) == length(mesh0.edges)
    @test maximum(norm.(state.mesh.points .- (mesh0.points .+ Ref(exact_shift)))) < 1e-10
    @test front_measure(state) ≈ front_measure(FrontState(mesh0)) rtol=1e-10
    @test curve_centroid(state.mesh) ≈ curve_centroid(mesh0) + exact_shift atol=1e-10
    assert_no_nan(state.mesh)
end