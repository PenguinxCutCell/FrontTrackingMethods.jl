# translate_circle.jl – Constant advection of a circle.
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays, LinearAlgebra

# Setup
N  = 64;  R = 1.0
pts   = [SVector(R*cos(2π*(i-1)/N), R*sin(2π*(i-1)/N)) for i in 1:N]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
mesh  = CurveMesh(pts, edges)

u_adv = SVector(1.0, 0.5)
eq = FrontEquation(; terms=(AdvectionTerm(u_adv),), front=mesh, integrator=RK2())

println("t = 0.0  centroid = ", round.(sum(mesh.points)/N, digits=4))
integrate!(eq, 1.0; dt=0.05)
pts_final = eq.state.mesh.points
c = sum(pts_final) / N
println("t = 1.0  centroid = ", round.(c, digits=4), "  exact = (1.0, 0.5)")
println("Area drift = ", abs(front_enclosed_measure(eq.state) - π*R^2) / (π*R^2))
