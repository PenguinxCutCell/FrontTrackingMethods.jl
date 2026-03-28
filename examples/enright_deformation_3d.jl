# enright_deformation_3d.jl – Enright 3-D deformation benchmark.
#
# One of the most demanding surface-tracking tests.
# The divergence-free velocity stretches a sphere into a thin sheet,
# which is then reversed.
#
# Setup (Enright et al. 2002):
#   - unit-cube domain [0,1]³
#   - initial sphere centered at (0.35, 0.35, 0.35), radius 0.15
#   - reversal at T/2, with full cycle T = 3.0
#
# Comparison:
#   A. no remeshing
#   B. tangential redistribution every 5 steps
#   C. experimental surface remesher every 5 steps
#
# Metrics:
#   - final volume drift
#   - surface area drift
#   - approximate symmetric Hausdorff to initial
#   - minimum triangle quality over time

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

# ─── Parameters ──────────────────────────────────────────────────────────────
T_end      = 3.0   # full cycle
dt         = 0.01
refinement = 3
redistribute_every = 5

mesh0 = make_sphere_benchmark_surface(
    center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=refinement)
u = enright_3d(; T=T_end)

state0 = FrontState(mesh0)
V0 = front_enclosed_measure(state0)
A0 = front_measure(state0)
c0 = front_centroid(state0)

println("=" ^ 60)
println("Enright 3-D deformation benchmark")
println("  T = $T_end, dt = $dt, refinement = $refinement")
println("  N_verts = $(length(mesh0.points)),  N_faces = $(length(mesh0.faces))")
println("  Initial volume = $(round(V0, sigdigits=5))")
println("  Initial surface area = $(round(A0, sigdigits=5))")
println("=" ^ 60)

function run_and_report(label, redistributor)
    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK3(),
        # Keep equation redistribution off; we control frequency via callbacks.
        redistributor = nothing,
    )
    min_angle_hist = Float64[]
    degen_hist = Float64[]
    cb_quality = EveryNSteps(5, (state, t, step) -> begin
        q = surface_quality_summary(state.mesh; degenerate_atol=1e-12)
        push!(min_angle_hist, q.angle.min_angle)
        push!(degen_hist, q.degenerate_fraction)
        nothing
    end)
    cb_redistribute = EveryNSteps(redistribute_every, (state, t, step) -> begin
        redistributor === nothing || redistribute!(state, redistributor)
        nothing
    end)
    cb = compose_callbacks(cb_quality, cb_redistribute)

    integrate!(eq, T_end; dt=dt, callback=cb)

    state = current_state(eq)
    Vf    = front_enclosed_measure(state)
    Af    = front_measure(state)
    cf    = front_centroid(state)
    dH    = symmetric_hausdorff_surface(state.mesh, mesh0)
    qf    = surface_quality_summary(state.mesh; degenerate_atol=1e-12)

    println("\n[$label]")
    @printf "  Volume drift:       %+.4e  (relative: %.4e)\n" (Vf - V0) abs(Vf - V0) / max(V0, eps())
    @printf "  Surface area drift: %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / max(A0, eps())
    @printf "  Centroid drift:      %.4e\n" norm(cf - c0)
    @printf "  Hausdorff to init:   %.4e\n" dH
    @printf "  Final min angle:     %.4e rad\n" qf.angle.min_angle
    @printf "  Worst min angle:     %.4e rad\n" (isempty(min_angle_hist) ? qf.angle.min_angle : minimum(min_angle_hist))
    @printf "  Worst degen frac:    %.4e\n" (isempty(degen_hist) ? qf.degenerate_fraction : maximum(degen_hist))
    println("  Remeshing stats:     splits/collapses/flips not tracked in v0.2 fixed-topology path")
end

run_and_report("A. No remeshing",                  nothing)
run_and_report("B. Tangential redistribution",     SurfaceTangentialRedistributor(; iterations=3))
run_and_report("C. Experimental surface remesher", ExperimentalSurfaceRemesher(; iterations=3))

# ─── Reversibility Analysis ──────────────────────────────────────────────────
println("\n" * "=" ^ 60)
println("Reversibility Analysis")
println("=" ^ 60)

# Test reversibility with best method
eq_rev = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK2(),
    redistributor = nothing,  # no remeshing for reversibility test (best-case scenario)
)

# Integrate forward then backward
times_for_rev = range(0.0, T_end, length=61)
half_idx = div(length(times_for_rev), 2) + 1
half_time = times_for_rev[half_idx]

# Integrate to half time
integrate!(eq_rev, half_time; dt=dt)
state_half = deepcopy(current_state(eq_rev))
V_half = front_enclosed_measure(state_half)
A_half = front_measure(state_half)
@printf "At t = T/2 = %.2f:\n" half_time
@printf "  Volume:       %.6e\n" V_half
@printf "  Surface area: %.6e\n" A_half

# Integrate to end
integrate!(eq_rev, T_end; dt=dt)
state_end = current_state(eq_rev)
V_end = front_enclosed_measure(state_end)
A_end = front_measure(state_end)
dH_end = symmetric_hausdorff_surface(state_end.mesh, mesh0)

@printf "\nAt t = T = %.2f:\n" T_end
@printf "  Volume:              %.6e (initial: %.6e, drift: %+.4e)\n" V_end V0 (V_end - V0)
@printf "  Surface area:        %.6e (initial: %.6e, drift: %+.4e)\n" A_end A0 (A_end - A0)
@printf "  Hausdorff to init:   %.6e\n" dH_end
@printf "  Expected reversion: Should be CLOSE to initial state!\n"
if dH_end < 0.01
    println("  ✓ GOOD: Sphere appears to revert properly (HausdorffOk < 0.01)")
else
    println("  ✗ WARNING: Large Hausdorff distance suggests poor reversibility")
    println("    Possible causes:")
    println("    - dt too large (try dt = 0.01 or smaller)")
    println("    - Insufficient remeshing (try more iterations)")
    println("    - Solver error accumulation (try higher-order integrator)")
end

println("\nDone.")

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("enright_deformation_3d")

# ─── Visualization Setup ─────────────────────────────────────────────────────
# Camera parameters for better 3D structure visualization
# Adjust azimuth (0-360°) and elevation (-90 to 90°) for different views:
#   - azimuth=0, elevation=0:    looking along +x axis
#   - azimuth=45, elevation=30:  isometric view (good default)
#   - azimuth=0, elevation=90:   top-down view
#   - azimuth=90, elevation=0:   looking along +y axis
camera_params = Dict(
    :azimuth => 95.0,      # rotate around z-axis (degrees)
    :elevation => 3.0,    # tilt up from xy-plane (degrees)
)

eq_anim = FrontEquation(;
    terms        = AdvectionTerm(u),
    front        = deepcopy(mesh0),
    integrator   = RK3(),
    redistributor = nothing,  # no remeshing for animation (best-case scenario)
)

xlims = (0.0, 1.0)
ylims = (0.0, 1.0)
zlims = (0.0, 1.0)

snapshot(eq_anim, joinpath(outdir, "initial.png"); title="enright: initial", wireframe=true, xlims=xlims, ylims=ylims, zlims=zlims, camera_params...)
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T_end, 120);
    title="enright", wireframe=true, xlims=xlims, ylims=ylims, zlims=zlims, camera_params...)
snapshot(eq_anim, joinpath(outdir, "final.png"); title="enright: final", wireframe=true, xlims=xlims, ylims=ylims, zlims=zlims, camera_params...)
println("Saved plotting outputs to: $outdir")
