# test_benchmark_zalesak_2d.jl – Zalesak disk rotation benchmark test.
#
# Verifies that a slotted disk completes one rigid rotation with:
# - acceptable area drift
# - centroid returning to initial position
# - front-to-front distance within tolerance

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Zalesak disk – one full rotation" begin
    center_disk = SVector(0.5, 0.75)
    center_rot  = SVector(0.5, 0.5)
    omega       = 2π
    T_rev       = 1.0
    dt          = 0.01

    mesh0 = make_zalesak_disk_curve(;
        center     = center_disk,
        R          = 0.15,
        slot_width = 0.05,
        slot_depth = 0.25,
        N_arc      = 128,
        N_slot     = 32,
    )

    u  = rigid_rotation_2d(; center=center_rot, omega=omega)
    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = CurveEqualArcRedistributor(),
    )
    integrate!(eq, T_rev; dt=dt)

    state = current_state(eq)

    # Area should be preserved (within 2%)
    A0 = curve_area(mesh0)
    Af = curve_area(state.mesh)
    @test abs(Af - A0) / abs(A0) < 0.02

    # Centroid should return to initial position
    c0 = curve_centroid(mesh0)
    cf = curve_centroid(state.mesh)
    @test norm(cf - c0) < 5e-3

    # Front-to-front recovery
    d_hausdorff = symmetric_hausdorff_curve(state.mesh, mesh0)
    @test d_hausdorff < 0.02

    assert_no_nan(state.mesh)
    assert_closed_curve(state.mesh)
end
