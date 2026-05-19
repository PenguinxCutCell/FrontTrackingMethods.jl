@testset "DiffEqIntegrator extension" begin
    if Base.find_package("OrdinaryDiffEq") === nothing
        @test_skip "OrdinaryDiffEq not available in this environment"
    else
        @eval using OrdinaryDiffEq

        R = 1.0
        N = 64
        u = SVector(0.7, -0.3)
        tf = 0.5

        mesh = make_circle(R, N)
        eq = FrontEquation(;
            terms=(AdvectionTerm(u),),
            front=mesh,
            t=0.0,
            integrator=DiffEqIntegrator(Tsit5(); cfl=0.8, abstol=1e-11, reltol=1e-11),
        )
        integrate!(eq, tf; dt=0.05)

        pts = eq.state.mesh.points
        center = sum(pts) / length(pts)
        @test norm(center - u * tf) < 1e-10

        radii = [norm(p - center) for p in pts]
        @test maximum(abs.(radii .- R)) < 1e-10
    end
end
