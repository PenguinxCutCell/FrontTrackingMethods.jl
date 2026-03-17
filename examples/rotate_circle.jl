# rotate_circle.jl – Rigid rotation benchmark for a 2-D circle.
#
# Benchmark: rigid rotation of a circle about the domain center.
#
# Setup:
#   - circle centered at (0.5, 0.75) with radius 0.15
#   - rotation center (0.5, 0.5), angular speed ω = 2π (one revolution per unit time)
#   - run for one full period T = 1.0
#
# Expected behavior:
#   - centroid, area, and shape all return to initial values
#
# Diagnostics printed to stdout:
#   - area drift after one period
#   - centroid drift
#   - Hausdorff distance to initial front
#   - edge-length stats
#
# Comparison:
#   - no redistribution
#   - equal-arclength redistribution
#   - adaptive remeshing

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Setup ───────────────────────────────────────────────────────────────────
center_disk = SVector(0.5, 0.75)
center_rot  = SVector(0.5, 0.5)
R           = 0.15
N           = 256
omega       = 2π
T_rev       = 1.0
dt          = 0.005

mesh0 = make_circle_benchmark_curve(center=center_disk, R=R, N=N)
u     = rigid_rotation_2d(; center=center_rot, omega=omega)

A0 = front_enclosed_measure(FrontState(mesh0))
c0 = front_centroid(FrontState(mesh0))

println("=" ^ 60)
println("Rigid rotation – circle benchmark")
println("  R = $R, N = $N, dt = $dt, T = $T_rev")
println("  Initial area = $(round(A0, sigdigits=6))")
println("=" ^ 60)

# ─── Helper to run and report ────────────────────────────────────────────────
function run_and_report(label, redistributor)
    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = redistributor,
    )
    integrate!(eq, T_rev; dt=dt)

    state = current_state(eq)
    Af    = front_enclosed_measure(state)
    cf    = front_centroid(state)
    dH    = symmetric_hausdorff_curve(state.mesh, mesh0)

    println("\n[$label]")
    @printf "  Area drift:       %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / A0
    @printf "  Centroid drift:    %.4e\n" norm(cf - c0)
    @printf "  Hausdorff to init: %.4e\n" dH
    h_stats = edge_length_stats(state.mesh)
    @printf "  Edge lengths:  min=%.4e  max=%.4e  spread=%.2f\n" h_stats.min h_stats.max h_stats.spread
    return state
end

edge_length_stats(m::CurveMesh) = begin
    geom = FrontIntrinsicOps.compute_geometry(m)
    ls = geom.edge_lengths
    (min=minimum(ls), max=maximum(ls), mean=sum(ls)/length(ls),
     spread=maximum(ls)/max(minimum(ls), eps()))
end

run_and_report("No redistribution",        nothing)
run_and_report("Equal-arclength",          CurveEqualArcRedistributor())
run_and_report("Adaptive remesher",        AdaptiveCurveRemesher(; iterations=5))

println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("rotate_circle")

eq_anim = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = CurveEqualArcRedistributor(; every=5),
)
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="rotate_circle: initial")
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T_rev, 140);
    title="rotate_circle")
snapshot(eq_anim, joinpath(outdir, "final.png"); title="rotate_circle: final")
println("Saved plotting outputs to: $outdir")
