using Test
using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using Makie
using CairoMakie

@testset "FrontTrackingMethods MakieExt smoke" begin
    FrontTrackingMethods.set_makie_theme!()
    @test FrontTrackingMethods.makie_theme() isa Makie.Theme

    mesh = make_circle_benchmark_curve(center=SVector(0.0, 0.0), R=1.0, N=64)
    state = FrontState(mesh)
    eq = FrontEquation(; terms=(AdvectionTerm(SVector(0.1, 0.0)),), front=deepcopy(mesh), integrator=RK2())

    fig1, _, _ = plot_state(state; show_vertices=true, title="state", xlims=(-1.5, 1.5), ylims=(-1.5, 1.5))
    @test fig1 isa Makie.Figure

    fig2, _, _ = plot_equation(eq)
    @test fig2 isa Makie.Figure

    png_eq = tempname() * ".png"
    out = snapshot(eq, png_eq)
    @test out == png_eq
    @test isfile(png_eq)

    # animation helper (GIF is generally robust in CI)
    gif = tempname() * ".gif"
    times = range(current_time(eq) + 0.01, stop=current_time(eq) + 0.05, length=5)
    outgif = record_evolution!(eq, gif, times; xlims=(-1.5, 1.5), ylims=(-1.5, 1.5))
    @test outgif == gif
    @test isfile(gif)

    # 3-D surface smoke
    mesh3 = make_sphere_benchmark_surface(center=SVector(0.0, 0.0, 0.0), R=1.0, refinement=1)
    state3 = FrontState(mesh3)
    fig3, _, _ = plot_state(state3; wireframe=true, xlims=(-1.5, 1.5), ylims=(-1.5, 1.5), zlims=(-1.5, 1.5))
    @test fig3 isa Makie.Figure
    png3 = tempname() * ".png"
    save(png3, fig3)
    @test isfile(png3)
end
