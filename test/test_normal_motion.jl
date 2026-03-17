# test_normal_motion.jl – Tests for NormalMotionTerm (inward/outward motion).

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Constant normal motion – circle (inward, shrink)" begin
    R0    = 1.0
    speed = 0.1   # inward speed (Vn > 0 means inward for CCW circles)
    dt    = 0.001
    tf    = 0.3
    mesh0 = make_circle_curve(R=R0, N=128)

    vn    = (x, t, state) -> speed   # positive = inward = shrink
    eq    = FrontEquation(; terms=NormalMotionTerm(vn), front=mesh0, integrator=RK2())
    integrate!(eq, tf; dt=dt)

    # Radius should have decreased by speed * tf
    R_expected = R0 - speed * tf
    pts = current_state(eq).mesh.points
    center = sum(pts) / length(pts)
    radii = [norm(p - center) for p in pts]
    @test all(r -> isapprox(r, R_expected; rtol=0.02), radii)
    @test minimum(radii) > 0
end

@testset "Outward normal motion – circle expands" begin
    R0    = 0.5
    speed = 0.1
    dt    = 0.01
    tf    = 0.5
    mesh0 = make_circle_curve(R=R0, N=64)

    vn = (x, t, state) -> -speed   # negative = outward = expand
    eq = FrontEquation(; terms=NormalMotionTerm(vn), front=mesh0, integrator=RK2())
    integrate!(eq, tf; dt=dt)

    pts = current_state(eq).mesh.points
    center = sum(pts) / length(pts)
    radii = [norm(p - center) for p in pts]
    R_expected = R0 + speed * tf
    @test all(r -> isapprox(r, R_expected; rtol=0.02), radii)
end
