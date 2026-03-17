# serpentine_2d.jl – Serpentine deformation benchmark.
#
# This is a severe shear-stretch test.  The divergence-free serpentine flow
# creates long thin filaments that must be tracked without catastrophic marker
# collapse.
#
# Setup:
#   - domain [0,1]²
#   - initial circle centered at (0.5, 0.75), radius 0.15
#   - time-reversed field with reversal at Tmax/2
#   - final time Tmax = 3.0
#
# Metrics:
#   - relative area drift
#   - front-to-front L² and Hausdorff distances to initial
#   - minimum edge length over time (collapse indicator)

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Parameters ──────────────────────────────────────────────────────────────
Tmax  = 3.0
dt    = 0.01
N     = 256

mesh0 = make_circle_benchmark_curve(center=SVector(0.5, 0.75), R=0.15, N=N)
u     = serpentine_2d(; Tmax=Tmax)

state0 = FrontState(mesh0)
A0 = front_enclosed_measure(state0)
c0 = front_centroid(state0)

println("=" ^ 60)
println("Serpentine deformation benchmark")
println("  Tmax = $Tmax, dt = $dt, N = $N")
println("  Initial area = $(round(A0, sigdigits=6))")
println("=" ^ 60)

# Track minimum edge length over time
min_edge_history = Float64[]
step_cb = EveryNSteps(10, (state, t, step) -> begin
    geom = state.geom
    push!(min_edge_history, minimum(geom.edge_lengths))
end)

eq = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = AdaptiveCurveRemesher(; iterations=5, protect_corners=false),
)
integrate!(eq, Tmax; dt=dt, callback=step_cb)

state = current_state(eq)
Af    = front_enclosed_measure(state)
cf    = front_centroid(state)
dH    = symmetric_hausdorff_curve(state.mesh, mesh0)
dL2   = l2_distance_curve(state.mesh, mesh0)

println("\nFinal diagnostics (t = $Tmax):")
@printf "  Area drift:       %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / A0
@printf "  Centroid drift:    %.4e\n" norm(cf - c0)
@printf "  Hausdorff to init: %.4e\n" dH
@printf "  L² to init:        %.4e\n" dL2
if !isempty(min_edge_history)
    @printf "  Min edge (min over time): %.4e\n" minimum(min_edge_history)
end
println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("serpentine_2d")

eq_anim = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = AdaptiveCurveRemesher(; iterations=5, protect_corners=false),
)
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="serpentine: initial")
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(Tmax, 180);
    title="serpentine")
snapshot(eq_anim, joinpath(outdir, "final.png"); title="serpentine: final")
println("Saved plotting outputs to: $outdir")
