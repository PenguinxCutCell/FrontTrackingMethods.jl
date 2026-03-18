# oscillating_drop_2d.jl – 2-D oscillating drop with animation + diagnostic plot.

using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra
using Statistics
using Printf
using CairoMakie
include("example_utils.jl")
using .ExampleUtils

center = SVector(0.5, 0.5)
R0     = 0.18
N      = 256
A      = 0.6
omega  = 8π
tf     = 2.0
dt     = 0.002

mesh0 = make_circle_benchmark_curve(; center=center, R=R0, N=N)

@inline function mode2_cos(dx, dy)
    r2 = dx * dx + dy * dy
    r2 > eps(Float64) ? (dx * dx - dy * dy) / r2 : 0.0
end

vn = (x, t, state) -> begin
    c = front_centroid(state)
    dx = x[1] - c[1]
    dy = x[2] - c[2]
    -A * mode2_cos(dx, dy) * sin(omega * t)
end

eq = FrontEquation(;
    terms         = (NormalMotionTerm(vn),),
    front         = deepcopy(mesh0),
    integrator    = RK2(),
    redistributor = CurveEqualArcRedistributor(; every=5),
)

ts = Float64[]
rbar = Float64[]
areas = Float64[]

function save_diag!(state, t)
    c = front_centroid(state)
    push!(ts, Float64(t))
    push!(rbar, mean(norm(p - c) for p in state.mesh.points))
    push!(areas, front_enclosed_measure(state))
    return nothing
end

save_diag!(eq.state, 0.0)
integrate!(eq, tf; dt=dt, callback=TimeIntervalCallback(tf / 220, (state, t, step) -> save_diag!(state, t)))

statef = current_state(eq)
meshf  = statef.mesh

A0 = front_enclosed_measure(FrontState(mesh0))
Af = front_enclosed_measure(statef)
c0 = front_centroid(FrontState(mesh0))
cf = front_centroid(statef)
dH = symmetric_hausdorff_curve(meshf, mesh0)

println("=" ^ 60)
println("Oscillating drop – 2-D")
println("  R0 = $R0, N = $N, A = $A, omega = $omega, tf = $tf, dt = $dt")
@printf "  Final centroid drift: %.4e\n" norm(cf - c0)
@printf "  Final area drift:     %.4e\n" abs(Af - A0) / max(abs(A0), eps())
@printf "  Hausdorff to initial: %.4e\n" dH
println("=" ^ 60)

function closed_xy(mesh::CurveMesh)
    xs = [p[1] for p in mesh.points]
    ys = [p[2] for p in mesh.points]
    push!(xs, xs[1])
    push!(ys, ys[1])
    return xs, ys
end

FrontTrackingMethods.set_makie_theme!()
outdir = ExampleUtils.output_dir_for("oscillating_drop_2d")

eq_anim = FrontEquation(;
    terms         = (NormalMotionTerm(vn),),
    front         = deepcopy(mesh0),
    integrator    = RK2(),
    redistributor = CurveEqualArcRedistributor(; every=5),
)

snapshot(eq_anim, joinpath(outdir, "initial.png"); title="oscillating_drop_2d: initial")
fig_anim, ax_anim, _ = plot_equation(eq_anim; title="t = 0.0")
times = default_times(tf, 200)
record(fig_anim, joinpath(outdir, "animation.mp4"), times) do ttarget
    integrate!(eq_anim, ttarget; dt=dt)
    plot_state(current_state(eq_anim);
        figure=fig_anim,
        axis=ax_anim,
        clear_axis=true,
        title="t = $(round(current_time(eq_anim), digits=4))")
end
snapshot(eq_anim, joinpath(outdir, "final.png"); title="oscillating_drop_2d: final")

x0, y0 = closed_xy(mesh0)
xf, yf = closed_xy(meshf)

fig = Figure(size=(1100, 450))
ax1 = Axis(fig[1, 1]; title="Oscillating drop contour", xlabel="x", ylabel="y", aspect=DataAspect())
lines!(ax1, x0, y0; linewidth=2.0, label="initial")
lines!(ax1, xf, yf; linewidth=2.0, label="final")
axislegend(ax1; position=:rb)

ax2 = Axis(fig[1, 2]; title="Mean radius history", xlabel="t", ylabel="R̄(t)", aspect=DataAspect())
lines!(ax2, ts, rbar; linewidth=2.0)

save(joinpath(outdir, "plot.png"), fig)
println("Saved plotting outputs to: $outdir")
