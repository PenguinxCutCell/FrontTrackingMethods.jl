# zalesak_disk_rotation.jl – Zalesak slotted-disk rigid-rotation benchmark.
#
# This is the canonical 2-D benchmark for interface-tracking methods.
# The sharp slot corners are preserved as long as the method handles
# them without excessive smoothing.
#
# Setup (following Zalesak 1979):
#   - unit-square domain [0,1]²
#   - slotted disk centered at (0.5, 0.75), radius 0.15
#   - slot width 0.05, slot depth 0.25
#   - rigid-body rotation about (0.5, 0.5), ω = 2π (one revolution per unit time)
#
# Metrics:
#   - final area drift
#   - centroid drift
#   - front-to-front Hausdorff distance to initial
#   - edge-length stats

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Parameters ──────────────────────────────────────────────────────────────
center_disk = SVector(0.5, 0.75)
center_rot  = SVector(0.5, 0.5)
omega       = 2π
T_rev       = 1.0
dt          = 0.005

mesh0 = make_zalesak_disk_curve(;
    center     = center_disk,
    R          = 0.15,
    slot_width = 0.05,
    slot_depth = 0.25,
    N_arc      = 256,
    N_slot     = 64,
)
u = rigid_rotation_2d(; center=center_rot, omega=omega)

state0 = FrontState(mesh0)
A0 = front_enclosed_measure(state0)
c0 = front_centroid(state0)

println("=" ^ 60)
println("Zalesak disk – rigid rotation benchmark")
println("  dt = $dt, T = $T_rev, N = $(length(mesh0.points))")
println("  Initial area = $(round(A0, sigdigits=6))")
println("=" ^ 60)

# ─── Run ─────────────────────────────────────────────────────────────────────
eq = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = CurveEqualArcRedistributor(; every=5),
)
integrate!(eq, T_rev; dt=dt)

state = current_state(eq)
Af    = front_enclosed_measure(state)
cf    = front_centroid(state)
dH    = symmetric_hausdorff_curve(state.mesh, mesh0)

geom  = state.geom
ls    = geom.edge_lengths

println("\nFinal diagnostics (after one full revolution):")
@printf "  Area drift:       %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / A0
@printf "  Centroid drift:    %.4e\n" norm(cf - c0)
@printf "  Hausdorff to init: %.4e\n" dH
@printf "  Edge lengths:  min=%.4e  max=%.4e\n" minimum(ls) maximum(ls)
println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("zalesak_disk_rotation")

eq_anim = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = CurveEqualArcRedistributor(; every=5),
)
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="zalesak: initial")
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T_rev, 160);
    title="zalesak")
snapshot(eq_anim, joinpath(outdir, "final.png"); title="zalesak: final")
println("Saved plotting outputs to: $outdir")
