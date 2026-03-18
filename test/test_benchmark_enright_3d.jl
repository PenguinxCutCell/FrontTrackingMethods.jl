# test_benchmark_enright_3d.jl – Enright 3-D deformation benchmark test.
#
# Smoke test on a coarse sphere:
# - verify no NaN
# - volume drift below threshold
# - surface quality does not collapse

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "Enright 3-D deformation – robustness (coarse mesh)" begin
    # Partial cycle for CI runtime while preserving severe deformation signal.
    T_end = 0.75
    dt    = 0.05

    mesh0 = make_sphere_benchmark_surface(
        center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=2)
    u = enright_3d(; T=3.0)

    eq_none = FrontEquation(;
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
    )
    integrate!(eq_none, T_end; dt=dt)
    state_none = current_state(eq_none)

    eq_remesh = FrontEquation(;
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=ExperimentalSurfaceRemesher(; iterations=2, strength=0.25, volume_correction=true),
    )
    integrate!(eq_remesh, T_end; dt=dt)
    state_remesh = current_state(eq_remesh)

    assert_no_nan(state_none.mesh)
    assert_no_nan(state_remesh.mesh)

    V0 = surface_volume(mesh0)
    V_none = surface_volume(state_none.mesh)
    V_remesh = surface_volume(state_remesh.mesh)
    @test abs(V_none - V0) / abs(V0) < 0.40
    @test abs(V_remesh - V0) / abs(V0) < 0.30

    d_none = symmetric_hausdorff_surface(state_none.mesh, mesh0)
    d_remesh = symmetric_hausdorff_surface(state_remesh.mesh, mesh0)
    @test isfinite(d_none)
    @test isfinite(d_remesh)
    @test d_remesh <= 1.15 * d_none + 1e-8

    q_none = surface_quality_summary(state_none.mesh; degenerate_atol=1e-12)
    q_remesh = surface_quality_summary(state_remesh.mesh; degenerate_atol=1e-12)
    @test q_none.degenerate_fraction <= 0.05
    @test q_remesh.degenerate_fraction <= 0.05
    @test q_remesh.angle.min_angle + 1e-12 >= q_none.angle.min_angle * 0.85

    assert_surface_is_reasonable(state_none.mesh; min_area_tol=1e-12)
    assert_surface_is_reasonable(state_remesh.mesh; min_area_tol=1e-12)
end
