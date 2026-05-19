# diffeq_integrator_demo.jl
#
# Demonstrates the optional OrdinaryDiffEq.jl weak extension:
#   using OrdinaryDiffEq
#   integrator = DiffEqIntegrator(Tsit5(); abstol=1e-11, reltol=1e-11)
#
# The example rotates an off-center circle for one full period and compares
# a coarse fixed-step RK2 run with a Tsit5-backed DiffEqIntegrator run using
# the same outer step size.

using FrontIntrinsicOps
using FrontTrackingMethods
using LinearAlgebra
using OrdinaryDiffEq
using Printf
using StaticArrays

center_disk = SVector(0.5, 0.75)
center_rot = SVector(0.5, 0.5)
R = 0.15
N = 96
tf = 1.0
dt = 0.125

mesh0 = make_circle_benchmark_curve(center=center_disk, R=R, N=N)
velocity = rigid_rotation_2d(; center=center_rot, omega=2π)
state0 = FrontState(mesh0)
A0 = front_enclosed_measure(state0; correction=:arc)
c0 = front_centroid(state0)

function run_case(label, integrator)
    eq = FrontEquation(;
        terms=AdvectionTerm(velocity),
        front=deepcopy(mesh0),
        t=0.0,
        integrator=integrator,
    )
    elapsed = @elapsed integrate!(eq, tf; dt=dt)

    state = current_state(eq)
    area = front_enclosed_measure(state; correction=:arc)
    centroid = front_centroid(state)
    hausdorff = symmetric_hausdorff_curve(state.mesh, mesh0)

    @printf("\n[%s]\n", label)
    @printf("  outer steps:      %d\n", eq.step)
    @printf("  elapsed:          %.4e s\n", elapsed)
    @printf("  area rel. error:  %.4e\n", abs(area - A0) / A0)
    @printf("  centroid drift:   %.4e\n", norm(centroid - c0))
    @printf("  Hausdorff error:  %.4e\n", hausdorff)
    return eq
end

println("=" ^ 64)
println("OrdinaryDiffEq weak-extension demo")
println("Rigid rotation of an off-center circle for one period")
@printf("N = %d, outer dt = %.3f, tf = %.3f\n", N, dt, tf)
println("=" ^ 64)

run_case("RK2 fixed-stage integrator", RK2())
run_case("DiffEqIntegrator(Tsit5)", DiffEqIntegrator(Tsit5(); abstol=1e-11, reltol=1e-11))

println("\nDone.")
