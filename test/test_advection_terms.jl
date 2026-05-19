# test_advection_terms.jl – Tests for AdvectionTerm (pure translation).

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Pure translation – circle" begin
    R  = 1.0
    N  = 64
    T0 = Float64
    mesh0 = make_circle_curve(R=R, N=N)

    # Translate by (1, 0) for time 1.0
    ux = 1.0
    u  = rigid_translation_velocity(SVector{2,T0}(ux, 0.0))
    eq = FrontEquation(; terms=AdvectionTerm(u), front=mesh0, integrator=RK2())
    integrate!(eq, 1.0; dt=0.01)

    state = current_state(eq)
    pts   = state.mesh.points
    orig  = mesh0.points

    # Every vertex should have moved by (1, 0)
    for (p, p0) in zip(pts, orig)
        @test p[1] ≈ p0[1] + ux   atol=1e-12
        @test p[2] ≈ p0[2]        atol=1e-12
    end

    # Enclosed area should be preserved
    A_new = curve_area(state.mesh)
    A_ref = π * R^2
    @test A_new ≈ A_ref  rtol=1e-2

    # Centroid shift
    c_new = curve_centroid(state.mesh)
    c_old = curve_centroid(mesh0)
    @test c_new[1] ≈ c_old[1] + ux  atol=1e-12
    @test c_new[2] ≈ c_old[2]       atol=1e-12
end

@testset "Pure translation – RK3" begin
    mesh0 = make_circle_curve(R=0.5, N=32)
    u     = rigid_translation_velocity(SVector(0.0, 2.0))
    eq    = FrontEquation(; terms=AdvectionTerm(u), front=mesh0, integrator=RK3())
    integrate!(eq, 1.0; dt=0.05)

    pts  = current_state(eq).mesh.points
    orig = mesh0.points
    for (p, p0) in zip(pts, orig)
        @test p[1] ≈ p0[1]        atol=1e-12
        @test p[2] ≈ p0[2] + 2.0  atol=1e-12
    end
end

@testset "Pure translation – ForwardEuler" begin
    mesh0 = make_circle_curve(R=0.5, N=32)
    u     = rigid_translation_velocity(SVector(1.0, 1.0))
    eq    = FrontEquation(; terms=AdvectionTerm(u), front=mesh0, integrator=ForwardEuler())
    integrate!(eq, 1.0; dt=0.001)

    pts  = current_state(eq).mesh.points
    orig = mesh0.points
    for (p, p0) in zip(pts, orig)
        @test p[1] ≈ p0[1] + 1.0  atol=1e-10
        @test p[2] ≈ p0[2] + 1.0  atol=1e-10
    end
end

@testset "ProjectedAdvectionTerm – rigid rotation has zero normal velocity" begin
    mesh0 = make_circle_curve(R=1.0, N=64, center=SVector(0.0, 0.0))
    u = (x, t, state) -> SVector(-x[2], x[1])
    eq = FrontEquation(; terms=ProjectedAdvectionTerm(u), front=mesh0, integrator=RK2())
    integrate!(eq, 0.5; dt=0.05)

    pts = current_state(eq).mesh.points
    for (p, p0) in zip(pts, mesh0.points)
        @test norm(p - p0) < 1e-12
    end
end

@testset "ProjectedAdvectionTerm – FrontField and constant velocity modes" begin
    mesh0 = make_circle_curve(R=1.0, N=32)
    state = FrontState(mesh0)
    normals = state.geom.vertex_normals

    field = FrontField(copy(normals), mesh0, :vertex)
    V = [zero(SVector{2,Float64}) for _ in mesh0.points]
    FrontTrackingMethods.accumulate_term!(V, ProjectedAdvectionTerm(field), state, 0.0)
    for i in eachindex(V)
        @test norm(V[i] - normals[i]) < 1e-12
    end

    fill!(V, zero(SVector{2,Float64}))
    FrontTrackingMethods.accumulate_term!(V, ProjectedAdvectionTerm(SVector(1.0, 0.0)), state, 0.0)
    for i in eachindex(V)
        @test abs(dot(V[i], normals[i]) - dot(SVector(1.0, 0.0), normals[i])) < 1e-12
        @test norm(V[i] - dot(SVector(1.0, 0.0), normals[i]) * normals[i]) < 1e-12
    end
end

@testset "ProjectedAdvectionTerm – PointFront1D rejected" begin
    front = PointFront1D([0.0])
    eq = FrontEquation(; terms=ProjectedAdvectionTerm(1.0), front=front)
    @test_throws ErrorException integrate!(eq, 0.1; dt=0.01)
end
