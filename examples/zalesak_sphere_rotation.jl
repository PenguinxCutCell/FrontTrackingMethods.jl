# zalesak_sphere_rotation.jl – 3-D slotted-sphere rigid-rotation benchmark.

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

center_shape = SVector(0.5, 0.75, 0.5)
center_rot   = SVector(0.5, 0.5, 0.5)
R            = 0.15
slot_width   = 0.05
slot_depth   = 0.125
refinement   = 2
T_rev        = 1.0
dt           = 0.05

mesh0 = make_zalesak_sphere_surface(
    center=center_shape,
    R=R,
    slot_width=slot_width,
    slot_depth=slot_depth,
    refinement=refinement,
)

u = rigid_rotation_3d(; center=center_rot, axis=:z, omega=2π)

eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
integrate!(eq, T_rev; dt=dt)
meshf = current_state(eq).mesh

c0 = sum(mesh0.points) / length(mesh0.points)
cf = sum(meshf.points) / length(meshf.points)
A0 = FrontIntrinsicOps.measure(mesh0, FrontIntrinsicOps.compute_geometry(mesh0))
Af = FrontIntrinsicOps.measure(meshf, FrontIntrinsicOps.compute_geometry(meshf))
dH = symmetric_hausdorff_surface(meshf, mesh0)
qf = surface_quality_summary(meshf; degenerate_atol=1e-12)

slot_pts = [p for p in mesh0.points if abs(p[1] - center_shape[1]) <= 0.03 && p[3] <= center_shape[3] - 0.08]
slot_proxy = isempty(slot_pts) ? NaN : maximum(minimum(norm(p - q) for q in meshf.points) for p in slot_pts)

println("=" ^ 60)
println("Zalesak sphere rigid rotation (coarse)")
println("  refinement = $refinement, closed = $(is_closed(mesh0))")
@printf "  Centroid drift:      %.4e\n" norm(cf - c0)
@printf "  Surface area drift:  %.4e\n" abs(Af - A0) / max(abs(A0), eps())
@printf "  Hausdorff to init:   %.4e\n" dH
@printf "  Slot proxy error:    %.4e\n" slot_proxy
@printf "  Min angle (final):   %.4e rad\n" qf.angle.min_angle
@printf "  Degenerate fraction: %.4e\n" qf.degenerate_fraction
println("=" ^ 60)

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("zalesak_sphere_rotation")

eq_anim = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
snapshot(eq_anim, joinpath(outdir, "initial.png"); title="zalesak_sphere: initial", wireframe=true)
record_evolution!(eq_anim, joinpath(outdir, "animation.mp4"), default_times(T_rev, 100);
    title="zalesak_sphere", wireframe=true)
snapshot(eq_anim, joinpath(outdir, "final.png"); title="zalesak_sphere: final", wireframe=true)
println("Saved plotting outputs to: $outdir")
