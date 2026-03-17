# curve_shortening_circle.jl – Circle under curve-shortening flow.
#
# CurvatureMotionTerm(β) with β > 0 shrinks a CCW circle.
# Exact: R(t) = sqrt(R0^2 - 2β t).
using FrontIntrinsicOps, FrontTrackingMethods, StaticArrays, LinearAlgebra

N = 64; R0 = 1.0; beta = 0.1
pts   = [SVector(R0*cos(2π*(i-1)/N), R0*sin(2π*(i-1)/N)) for i in 1:N]
edges = [SVector(i, mod1(i+1,N)) for i in 1:N]
global meshh  = CurveMesh(pts, edges)


make_circle(R, N) = (
    pts = [SVector(R*cos(2π*(i-1)/N), R*sin(2π*(i-1)/N)) for i in 1:N];
    CurveMesh(pts, [SVector(i, mod1(i+1,N)) for i in 1:N])
)


for tf in [0.5, 1.0, 2.0, 4.0]
    eq = FrontEquation(; terms=(CurvatureMotionTerm(beta),), front=meshh, integrator=RK2())
    integrate!(eq, tf; dt=1e-3)
    pts_f  = eq.state.mesh.points
    center = sum(pts_f) / N
    R_num  = sum(norm(p - center) for p in pts_f) / N
    R_exact = sqrt(max(0.0, R0^2 - 2*beta*tf))
    println("t=$tf  R_num=$(round(R_num,digits=4))  R_exact=$(round(R_exact,digits=4)) Error=$(round(abs(R_num - R_exact)/R_exact * 100, digits=2))%")
    global meshh = make_circle(R0, N)   # reset for each test
end
