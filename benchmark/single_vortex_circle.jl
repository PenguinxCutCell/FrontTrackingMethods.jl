using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Printf

include("error_metrics.jl")
include("velocity_interpolation.jl")
include("visualization.jl")

function run_single_vortex_circle(; N::Int=256, Tperiod::Float64=8.0, face_grid::Int=128, save_visuals::Bool=true,
    xlims=(0.0, 1.0), ylims=(0.0, 1.0), zlims=nothing)
    R = 0.15
    c0 = SVector(0.5, 0.75)
    cfl = 0.125

    mesh0 = make_circle_benchmark_curve(center=c0, R=R, N=N)
    h = mean_edge_size(mesh0)

    umax = 1 / π
    dt = cfl * h / umax

    u_analytic = rider_kothe_single_vortex(T=Tperiod)
    u = face_bilinear_marker_velocity(u_analytic; Nx=face_grid, Ny=face_grid)
    eq = FrontEquation(
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=AdaptiveCurveRemesher(iterations=3, tangential_smooth=0.25),
    )
    integrate!(eq, Tperiod; dt=dt)

    meshf = current_state(eq).mesh
    errs = benchmark_errors(mesh0, meshf; circle_center=c0, circle_radius=R)

    @printf("[single_vortex_circle] N=%d h=%.5e dt=%.5e CFL=%.3f T=%.3f\n", N, h, dt, cfl, Tperiod)
    @printf("  E_area = %.6e\n", errs.E_area)
    @printf("  E_shape = %.6e\n", errs.E_shape)
    @printf("  E_sym = %.6e\n", errs.E_sym)

    if save_visuals
        eq_anim = FrontEquation(
            terms=AdvectionTerm(u),
            front=deepcopy(mesh0),
            integrator=RK2(),
            redistributor=AdaptiveCurveRemesher(iterations=3, tangential_smooth=0.25),
        )
        save_case_visuals!("single_vortex_circle", eq_anim, Tperiod;
            xlims=xlims,
            ylims=ylims,
            zlims=zlims,
        )
    end

    return (; N=N, h=h, dt=dt, CFL=cfl, Tperiod=Tperiod, errs...)
end

    run_single_vortex_circle()

