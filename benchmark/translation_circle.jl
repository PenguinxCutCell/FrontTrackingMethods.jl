using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Printf

include("error_metrics.jl")
include("velocity_interpolation.jl")
include("visualization.jl")

function run_translation_circle(; N::Int=256, face_grid::Int=128, save_visuals::Bool=true)
    R = 0.15
    c0 = SVector(0.25, 0.75)
    Tfinal = 1.0
    cfl = 0.125

    mesh0 = make_circle_benchmark_curve(center=c0, R=R, N=N)
    h = mean_edge_size(mesh0)

    u0 = SVector{2,Float64}(1.0, -1.0)
    u_analytic = (x, t, state) -> (t <= 0.5Tfinal ? u0 : -u0)
    u = face_bilinear_marker_velocity(u_analytic; Nx=face_grid, Ny=face_grid)
    dt = cfl * h / abs(u0[1])

    eq = FrontEquation(
        terms=AdvectionTerm(u),
        front=deepcopy(mesh0),
        integrator=RK2(),
        redistributor=CurveEqualArcRedistributor(every=5),
    )
    integrate!(eq, Tfinal; dt=dt)

    meshf = current_state(eq).mesh
    errs = benchmark_errors(mesh0, meshf; circle_center=c0, circle_radius=R)

    @printf("[translation_circle] N=%d h=%.5e dt=%.5e CFL=%.3f\n", N, h, dt, cfl)
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
        save_case_visuals!("translation_circle", eq_anim, Tfinal)
    end

    return (; N=N, h=h, dt=dt, CFL=cfl, errs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_translation_circle()
end
