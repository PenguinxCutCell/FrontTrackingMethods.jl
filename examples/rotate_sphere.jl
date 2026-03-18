# rotate_sphere.jl – Rigid rotation of a 3-D sphere.

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

center_rot = SVector(0.5, 0.5, 0.5)
center0    = SVector(0.5, 0.75, 0.5)
R          = 0.15
refinement = 3
omega      = 2π
T_rev      = 1.0
dt         = 0.05

mesh0 = make_sphere_benchmark_surface(center=center0, R=R, refinement=refinement)
u = rigid_rotation_3d(; center=center_rot, axis=:z, omega=omega)

eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
integrate!(eq, T_rev; dt=dt)
meshf = current_state(eq).mesh

c0 = sum(mesh0.points) / length(mesh0.points)
cf = sum(meshf.points) / length(meshf.points)
V0 = enclosed_measure(mesh0)
Vf = enclosed_measure(meshf)
A0 = FrontIntrinsicOps.measure(mesh0, FrontIntrinsicOps.compute_geometry(mesh0))
Af = FrontIntrinsicOps.measure(meshf, FrontIntrinsicOps.compute_geometry(meshf))
dH = symmetric_hausdorff_surface(meshf, mesh0)
qf = surface_quality_summary(meshf; degenerate_atol=1e-12)

println("=" ^ 60)
println("Rigid rotation – sphere")
println("  refinement = $refinement")
@printf "  Centroid return error: %.4e\n" norm(cf - c0)
@printf "  Volume drift:          %.4e\n" abs(Vf - V0) / max(abs(V0), eps())
@printf "  Surface-area drift:    %.4e\n" abs(Af - A0) / max(abs(A0), eps())
@printf "  Hausdorff to init:     %.4e\n" dH
@printf "  Min angle (final):     %.4e rad\n" qf.angle.min_angle
println("=" ^ 60)

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("rotate_sphere")

eq_anim = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="rotate_sphere: initial", wireframe=true)
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T_rev, 90);
    title="rotate_sphere", wireframe=true)
snapshot(eq_anim, joinpath(outdir, "final.png"); title="rotate_sphere: final", wireframe=true)
println("Saved plotting outputs to: $outdir")
