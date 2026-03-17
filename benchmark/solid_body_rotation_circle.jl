using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Printf

include("error_metrics.jl")
include("velocity_interpolation.jl")
include("visualization.jl")

function run_solid_body_rotation_circle(; N::Int=256, face_grid::Int=128, save_visuals::Bool=true)
    R = 0.15
    c0 = SVector(0.5, 0.75)
    center_rot = SVector(0.5, 0.5)
    ω = 2π
    Tfinal = 1.0
    cfl = π / 16

    mesh0 = make_circle_benchmark_curve(center=c0, R=R, N=N)
    h = mean_edge_size(mesh0)
    umax = π
    dt = cfl * h / umax

    u_analytic = rigid_rotation_2d(center=center_rot, omega=ω)
    u = face_bilinear_marker_velocity(u_analytic; Nx=face_grid, Ny=face_grid)
    eq = FrontEquation(
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=CurveEqualArcRedistributor(every=5),
    )
    integrate!(eq, Tfinal; dt=dt)

    meshf = current_state(eq).mesh
    errs = benchmark_errors(mesh0, meshf; circle_center=c0, circle_radius=R)

    @printf("[solid_body_rotation_circle] N=%d h=%.5e dt=%.5e CFL=%.5f\n", N, h, dt, cfl)
    @printf("  E_area = %.6e\n", errs.E_area)
    @printf("  E_shape = %.6e\n", errs.E_shape)
    @printf("  E_sym = %.6e\n", errs.E_sym)

    if save_visuals
        eq_anim = FrontEquation(
            terms=AdvectionTerm(u),
            front=deepcopy(mesh0),
            integrator=RK2(),
            redistributor=CurveEqualArcRedistributor(every=5),
        )
        save_case_visuals!("solid_body_rotation_circle", eq_anim, Tfinal)
    end

    return (; N=N, h=h, dt=dt, CFL=cfl, errs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_solid_body_rotation_circle()
end
