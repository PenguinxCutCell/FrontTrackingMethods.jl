# rider_kothe_single_vortex.jl – Rider–Kothe reversed single-vortex benchmark.
#
# This is the canonical reversible deformation test.
# A stream-function-based solenoidal velocity stretches the interface into a
# long thin filament; the velocity is then reversed to recover the initial shape.
#
# Setup:
#   - unit-square domain [0,1]²
#   - initial circle centered at (0.5, 0.75), radius 0.15
#   - reversal period T = 2.0
#
# Comparison:
#   A. no redistribution
#   B. equal-arclength redistribution
#   C. adaptive curve remeshing with corner protection
#
# Metrics:
#   - relative area drift at t=T
#   - centroid recovery error
#   - front-to-front Hausdorff distance to initial

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Parameters ──────────────────────────────────────────────────────────────
T     = 2.0
dt    = 0.02
N     = 256

mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
u     = rider_kothe_single_vortex(; T=T)

state0 = FrontState(mesh0)
A0 = front_enclosed_measure(state0)
c0 = front_centroid(state0)

println("=" ^ 60)
println("Rider–Kothe reversed single-vortex benchmark")
println("  T = $T, dt = $dt, N = $N")
println("  Initial area = $(round(A0, sigdigits=6))")
println("=" ^ 60)

# ─── Helper ─────────────────────────────────────────────────────────────────
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

    println("\n[$label]")
    @printf "  Area drift:       %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / A0
    @printf "  Centroid drift:    %.4e\n" norm(cf - c0)
    @printf "  Hausdorff to init: %.4e\n" dH
end

run_and_report("A. No redistribution",        nothing)
run_and_report("B. Equal-arclength",          CurveEqualArcRedistributor())
run_and_report("C. Adaptive remesher",        AdaptiveCurveRemesher(; iterations=5))

println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("rider_kothe_single_vortex")

eq_anim = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = AdaptiveCurveRemesher(; iterations=5),
)
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="rider-kothe: initial")
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T, 180);
    title="rider-kothe")
snapshot(eq_anim, joinpath(outdir, "final.png"); title="rider-kothe: final")
println("Saved plotting outputs to: $outdir")
