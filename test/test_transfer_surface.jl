# test_transfer_surface.jl – Tests for field transfer on surfaces.

using Test
using StaticArrays
using LinearAlgebra
using FrontTrackingMethods
using FrontIntrinsicOps

_rmse(a, b) = sqrt(sum((a .- b).^2) / max(length(a), 1))

function _field_errors(mesh_old, mesh_new, f)
    oldvals = [f(p) for p in mesh_old.points]
    expected = [f(p) for p in mesh_new.points]

    newvals_bary = zeros(Float64, length(mesh_new.points))
    newvals_near = zeros(Float64, length(mesh_new.points))
    transfer_vertex_field!(newvals_bary, mesh_old, oldvals, mesh_new; method=:barycentric)
    transfer_vertex_field!(newvals_near, mesh_old, oldvals, mesh_new; method=:nearest_vertex)
    return _rmse(newvals_bary, expected), _rmse(newvals_near, expected)
end

@testset "Field transfer – constants preserved (surface)" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    mesh_new = make_sphere_surface(R=1.0, refinement=2)

    oldvals = fill(42.0, length(mesh_old.points))
    newvals_b = zeros(Float64, length(mesh_new.points))
    newvals_n = zeros(Float64, length(mesh_new.points))

    transfer_vertex_field!(newvals_b, mesh_old, oldvals, mesh_new; method=:barycentric)
    transfer_vertex_field!(newvals_n, mesh_old, oldvals, mesh_new; method=:nearest_vertex)
    @test all(v -> isapprox(v, 42.0; atol=1e-10), newvals_b)
    @test all(v -> isapprox(v, 42.0; atol=1e-10), newvals_n)
end

@testset "Field transfer – smooth fields: barycentric vs nearest" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    mesh_new = make_sphere_surface(R=1.0, refinement=3)

    fields = [
        p -> p[1] + 2p[2] - p[3],
        p -> p[1]^2 - p[2] * p[3],
        p -> p[3] * (3p[3]^2 - 1),
    ]

    materially_better = 0
    for f in fields
        err_bary, err_near = _field_errors(mesh_old, mesh_new, f)
        @test err_bary <= err_near + 1e-10
        if err_near > 1e-12 && err_bary <= 0.95 * err_near
            materially_better += 1
        end
    end
    @test materially_better >= 2
end

@testset "Field transfer – after mild surface remeshing" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)
    state_new = FrontState(deepcopy(mesh_old))
    redistribute!(state_new, ExperimentalSurfaceRemesher(; iterations=2, strength=0.25, volume_correction=true))
    mesh_new = state_new.mesh

    err_bary, err_near = _field_errors(mesh_old, mesh_new, p -> p[1] + 2p[2] - p[3])
    @test err_bary <= err_near + 1e-10
end

@testset "Field transfer – after rigid sphere rotation" begin
    mesh_old = make_sphere_surface(R=1.0, refinement=2)

    θ = 0.35
    Rz = @SMatrix [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
    rot_pts = [Rz * p for p in mesh_old.points]
    mesh_new = SurfaceMesh(rot_pts, mesh_old.faces)

    err_bary, err_near = _field_errors(mesh_old, mesh_new, p -> p[1]^2 - p[2] * p[3])
    @test err_bary <= err_near + 1e-10
end

@testset "Field transfer – coarse Enright snapshot" begin
    mesh0 = make_sphere_benchmark_surface(center=SVector(0.35, 0.35, 0.35), R=0.15, refinement=2)
    u = enright_3d(; T=3.0)
    eq = FrontEquation(; terms=AdvectionTerm(u), front=deepcopy(mesh0), integrator=RK2())
    integrate!(eq, 0.25; dt=0.05)
    mesh_new = current_state(eq).mesh

    err_bary, err_near = _field_errors(mesh0, mesh_new, p -> p[1] + 2p[2] - p[3])
    @test isfinite(err_bary)
    @test isfinite(err_near)
    @test err_bary <= err_near + 1e-8
end
