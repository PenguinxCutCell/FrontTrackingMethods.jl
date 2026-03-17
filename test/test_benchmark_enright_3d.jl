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

include("test_utils.jl")

@testset "Enright 3-D deformation – smoke test (coarse mesh)" begin
    # Use T=0.5 (partial cycle) to keep CI runtime manageable
    T_end = 0.5
    dt    = 0.05

    mesh0 = make_sphere_benchmark_surface(
        center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=2)
    u = enright_3d(; T=3.0)

    eq = FrontEquation(;
        terms        = AdvectionTerm(u),
        front        = deepcopy(mesh0),
        integrator   = RK2(),
        redistributor = SurfaceTangentialRedistributor(; iterations=3),
    )
    integrate!(eq, T_end; dt=dt)

    state = current_state(eq)

    # No NaN / Inf
    assert_no_nan(state.mesh)

    # Volume drift within 10% (generous for severe deformation)
    V0 = surface_volume(mesh0)
    Vf = surface_volume(state.mesh)
    @test abs(Vf - V0) / abs(V0) < 0.10

    # Surface should still be reasonable
    assert_surface_is_reasonable(state.mesh; min_area_tol=1e-12)
end
