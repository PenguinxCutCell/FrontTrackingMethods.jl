# expand_circle_normal_speed.jl – Circle with constant normal speed.
#
# Sign convention for CurveMesh:
#   NormalMotionTerm(Vn) with Vn < 0 expands (inward normal → outward motion).
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays, LinearAlgebra

N = 64; R0 = 1.0; Vn = -0.5   # negative → expand
pts   = [SVector(R0*cos(2π*(i-1)/N), R0*sin(2π*(i-1)/N)) for i in 1:N]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
mesh  = CurveMesh(pts, edges)

tf = 0.5
eq = FrontEquation(; terms=(NormalMotionTerm(Vn),), front=mesh, integrator=RK2())
integrate!(eq, tf; dt=1e-3)

pts_f  = eq.state.mesh.points
center = sum(pts_f) / N
R_num  = sum(norm(p - center) for p in pts_f) / N
R_exact = R0 - Vn * tf   # inward normal, Vn<0 so R increases
println("R_numeric = $(round(R_num, digits=4)),  R_exact = $(round(R_exact, digits=4))")
