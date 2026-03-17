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

# ─── Parameters ──────────────────────────────────────────────────────────────
T_end      = 3.0   # full cycle
dt         = 0.05
refinement = 3

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
        integrator   = RK2(),
        redistributor = redistributor,
    )
    integrate!(eq, T_end; dt=dt)

    state = current_state(eq)
    Vf    = front_enclosed_measure(state)
    Af    = front_measure(state)
    cf    = front_centroid(state)
    dH    = symmetric_hausdorff_surface(state.mesh, mesh0)

    println("\n[$label]")
    @printf "  Volume drift:       %+.4e  (relative: %.4e)\n" (Vf - V0) abs(Vf - V0) / max(V0, eps())
    @printf "  Surface area drift: %+.4e  (relative: %.4e)\n" (Af - A0) abs(Af - A0) / max(A0, eps())
    @printf "  Centroid drift:      %.4e\n" norm(cf - c0)
    @printf "  Hausdorff to init:   %.4e\n" dH
end

run_and_report("A. No remeshing",                  nothing)
run_and_report("B. Tangential redistribution",     SurfaceTangentialRedistributor(; iterations=3))
run_and_report("C. Experimental surface remesher", ExperimentalSurfaceRemesher(; iterations=3))

println("\nDone.")
