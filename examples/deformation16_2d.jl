# deformation16_2d.jl – Strong 16-vortex deformation benchmark.
#
# Setup:
#   - domain [0,1]²
#   - initial circle centered at (0.5, 0.75), radius 0.15
#   - time-reversed deformation field
#   - final time T = 2.0

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf

T     = 2.0
dt    = 0.02
N     = 256

mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
u     = deformation16_2d(; T=T)

A0 = front_enclosed_measure(FrontState(mesh0))
c0 = front_centroid(FrontState(mesh0))

println("=" ^ 60)
println("16-vortex strong deformation benchmark")
println("  T = $T, dt = $dt, N = $N")
println("  Initial area = $(round(A0, sigdigits=6))")
println("=" ^ 60)

function run_and_report(label, redistributor)
    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = redistributor,
    )
    integrate!(eq, T; dt=dt)

    state = current_state(eq)
    Af    = front_enclosed_measure(state)
    cf    = front_centroid(state)
    dH    = symmetric_hausdorff_curve(state.mesh, mesh0)
    geom  = state.geom
    ls    = geom.edge_lengths

    println("\n[$label]")
    @printf "  Area drift:       %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / A0
    @printf "  Centroid drift:    %.4e\n" norm(cf - c0)
    @printf "  Hausdorff to init: %.4e\n" dH
    @printf "  Min edge length:   %.4e\n" minimum(ls)
end

run_and_report("No redistribution",  nothing)
run_and_report("Equal-arclength",    CurveEqualArcRedistributor())
run_and_report("Adaptive remesher",  AdaptiveCurveRemesher(; iterations=10))

println("\nDone.")
