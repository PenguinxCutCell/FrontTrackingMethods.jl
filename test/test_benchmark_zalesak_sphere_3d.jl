# test_benchmark_zalesak_sphere_3d.jl – Slotted-sphere rigid-rotation benchmark.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Zalesak sphere – rigid rotation (coarse quantitative)" begin
    center_shape = SVector(0.5, 0.75, 0.5)
    center_rot   = SVector(0.5, 0.5, 0.5)
    omega        = 2π
    T_rev        = 1.0
    dt           = 0.05

    mesh0 = make_zalesak_sphere_surface(
        center=center_shape,
        R=0.15,
        slot_width=0.05,
        slot_depth=0.125,
        refinement=2,
    )

    u = rigid_rotation_3d(; center=center_rot, axis=:z, omega=omega)
    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, T_rev; dt=dt)
    meshf = current_state(eq).mesh

    assert_no_nan(meshf)

    c0 = sum(mesh0.points) / length(mesh0.points)
    cf = sum(meshf.points) / length(meshf.points)
    @test norm(cf - c0) < 0.04

    A0 = FrontIntrinsicOps.measure(mesh0, FrontIntrinsicOps.compute_geometry(mesh0))
    Af = FrontIntrinsicOps.measure(meshf, FrontIntrinsicOps.compute_geometry(meshf))
    @test abs(Af - A0) / A0 < 0.06

    dH = symmetric_hausdorff_surface(meshf, mesh0)
    @test dH < 0.08

    q = surface_quality_summary(meshf; degenerate_atol=1e-12)
    @test q.degenerate_fraction <= 0.01
    @test q.angle.min_angle > deg2rad(8.0)

    # Slot-preservation proxy: points near initial slot should remain geometrically close.
    slot_pts = [p for p in mesh0.points if abs(p[1] - center_shape[1]) <= 0.03 && p[3] <= center_shape[3] - 0.08]
    @test !isempty(slot_pts)
    max_slot_dist = maximum(minimum(norm(p - q) for q in meshf.points) for p in slot_pts)
    @test max_slot_dist < 0.09
end
