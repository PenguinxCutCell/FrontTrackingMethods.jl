using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Printf

include("error_metrics.jl")
include("velocity_interpolation.jl")
include("visualization.jl")

function run_zalesak_disk_rotation(; N_arc::Int=256, face_grid::Int=128, save_visuals::Bool=true)
    center = SVector(0.5, 0.75)
    Tfinal = 1.0
    cfl = π / 16

    mesh0 = make_zalesak_disk_curve(
        center=center,
        R=0.15,
        slot_width=0.05,
        slot_depth=0.25,
        N_arc=N_arc,
        N_slot=max(16, N_arc ÷ 4),
    )

    h = mean_edge_size(mesh0)
    umax = π
    dt = cfl * h / umax

    u_analytic = rigid_rotation_2d(center=SVector(0.5, 0.5), omega=2π)
    u = face_bilinear_marker_velocity(u_analytic; Nx=face_grid, Ny=face_grid)
    eq = FrontEquation(
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=CurveEqualArcRedistributor(every=5),
    )
    integrate!(eq, Tfinal; dt=dt)

    meshf = current_state(eq).mesh
    errs = benchmark_errors(mesh0, meshf)

    @printf("[zalesak_disk_rotation] N_arc=%d h=%.5e dt=%.5e CFL=%.5f\n", N_arc, h, dt, cfl)
    @printf("  E_area = %.6e\n", errs.E_area)
    @printf("  E_sym = %.6e\n", errs.E_sym)

    if save_visuals
        eq_anim = FrontEquation(
            terms=AdvectionTerm(u),
            front=deepcopy(mesh0),
            integrator=RK2(),
            redistributor=CurveEqualArcRedistributor(every=5),
        )
        save_case_visuals!("zalesak_disk_rotation", eq_anim, Tfinal)
    end

    return (; N_arc=N_arc, h=h, dt=dt, CFL=cfl, E_area=errs.E_area, E_sym=errs.E_sym)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_zalesak_disk_rotation()
end
