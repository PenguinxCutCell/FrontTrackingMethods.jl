# test_curvature_motion.jl – Tests for CurvatureMotionTerm.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Curve-shortening flow on circle" begin
    # A circle is a stationary solution of curve-shortening flow:
    # κ = 1/R, so normal speed = -κ → radius shrinks as dR/dt = -1/R
    # This is equivalent to a speed = 1/R inward.
    R0  = 1.0
    N   = 128
    dt  = 0.001
    tf  = 0.1
    mesh0 = make_circle_curve(R=R0, N=N)
    eq    = FrontEquation(; terms=CurvatureMotionTerm(1.0), front=mesh0, integrator=RK2())
    integrate!(eq, tf; dt=dt)

    # The circle should have shrunk. Exact ODE: R(t) = sqrt(R0^2 - 2t)
    R_expected = sqrt(R0^2 - 2 * tf)
    pts   = current_state(eq).mesh.points
    radii = [norm(p) for p in pts]
    @test mean(radii) ≈ R_expected  rtol=0.02
    @test std(radii) / mean(radii)  < 0.01   # still approximately circular
end
