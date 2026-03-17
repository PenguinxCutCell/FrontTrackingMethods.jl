# translate_sphere.jl – Rigid translation of a 3-D sphere.
#
# Sanity check: rigid translation must preserve centroid shift,
# volume, and shape (Hausdorff distance to translated initial sphere).

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Parameters ──────────────────────────────────────────────────────────────
center0    = SVector(0.35, 0.35, 0.35)
R          = 0.15
refinement = 3
u_vec      = SVector(0.1, 0.05, -0.05)
tf         = 1.0
dt         = 0.05

mesh0 = make_sphere_benchmark_surface(center=center0, R=R, refinement=refinement)
u     = rigid_translation_velocity(u_vec)

state0 = FrontState(mesh0)
V0 = front_enclosed_measure(state0)
c0 = front_centroid(state0)

println("=" ^ 60)
println("Rigid translation – sphere")
println("  R = $R, refinement = $refinement, N_verts = $(length(mesh0.points))")
println("  u = $u_vec, T = $tf")
println("  Initial volume = $(round(V0, sigdigits=6))")
println("=" ^ 60)

eq = FrontEquation(;
    terms      = AdvectionTerm(u),
    front      = deepcopy(mesh0),
    integrator = RK2(),
)
integrate!(eq, tf; dt=dt)

state = current_state(eq)
Vf    = front_enclosed_measure(state)
cf    = front_centroid(state)

# Expected centroid after translation
c_exact = c0 + u_vec * tf
dH = symmetric_hausdorff_surface(state.mesh, mesh0)

println("\nFinal diagnostics (t = $tf):")
@printf "  Centroid exact:    %s\n" string(round.(c_exact, sigdigits=4))
@printf "  Centroid actual:   %s\n" string(round.(cf, sigdigits=4))
@printf "  Centroid error:    %.4e\n" norm(cf - c_exact)
@printf "  Volume drift:      %+.4e  (relative: %.4e)\n" (Vf - V0) abs(Vf - V0) / V0
@printf "  Hausdorff to init: %.4e\n" dH
println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("translate_sphere")

eq_anim = FrontEquation(;
    terms      = AdvectionTerm(u),
    front      = deepcopy(mesh0),
    integrator = RK2(),
)
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="translate_sphere: initial", wireframe=true)
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(tf, 90);
    title="translate_sphere", wireframe=true)
snapshot(eq_anim, joinpath(outdir, "final.png"); title="translate_sphere: final", wireframe=true)
println("Saved plotting outputs to: $outdir")
