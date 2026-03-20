# test_pointfront1d.jl – Minimal PointFront1D integration path tests.

using Test
using FrontTrackingMethods
using FrontIntrinsicOps

@testset "PointFront1D – state construction and geometry refresh" begin
    f1 = PointFront1D([0.3], true)
    s1 = FrontState(f1)
    @test s1.mesh == f1
    @test s1.geom.marker_positions == [0.3]
    @test s1.geom.vertex_normals == [1.0]

    f2 = PointFront1D([0.2, 0.8], true)
    s2 = FrontState(f2)
    @test s2.mesh == f2
    @test s2.geom.marker_positions == [0.2, 0.8]
    @test s2.geom.vertex_normals == [-1.0, 1.0]
    @test marker_gap(s2) ≈ 0.6

    set_vertex_coordinates!(s2, [0.1, 0.9])
    refresh_geometry!(s2)
    @test s2.geom.marker_positions == [0.1, 0.9]
end

@testset "PointFront1D – one-marker advection" begin
    x0 = 0.3
    u  = 2.0
    dt = 0.1
    tf = dt

    for integ in (ForwardEuler(), RK2(), RK3())
        front = PointFront1D([x0], true)
        eq = FrontEquation(; terms=(AdvectionTerm(u),), front=front, integrator=integ)
        integrate!(eq, tf; dt=dt)
        @test current_state(eq).mesh.x[1] ≈ x0 + u * tf atol=1e-12
    end
end

@testset "PointFront1D – one-marker normal motion sign follows normals" begin
    dt = 0.2
    tf = dt
    vn = 1.0

    front_right = PointFront1D([0.3], true)   # normal +1
    eq_right = FrontEquation(; terms=(NormalMotionTerm(vn),), front=front_right, integrator=RK2())
    integrate!(eq_right, tf; dt=dt)
    @test current_state(eq_right).mesh.x[1] ≈ 0.3 + vn * tf atol=1e-12

    front_left = PointFront1D([0.3], false)   # normal -1
    eq_left = FrontEquation(; terms=(NormalMotionTerm(vn),), front=front_left, integrator=RK2())
    integrate!(eq_left, tf; dt=dt)
    @test current_state(eq_left).mesh.x[1] ≈ 0.3 - vn * tf atol=1e-12
end

@testset "PointFront1D – two-marker normal motion (interval inside)" begin
    front = PointFront1D([0.2, 0.8], true)   # normals [-1, +1]
    eq = FrontEquation(; terms=(NormalMotionTerm(0.5),), front=front, integrator=RK2())
    integrate!(eq, 0.2; dt=0.2)

    x = current_state(eq).mesh.x
    @test x[1] < 0.2
    @test x[2] > 0.8
    @test (x[2] - x[1]) > 0.6
end

@testset "PointFront1D – two-marker normal motion (interval outside)" begin
    front = PointFront1D([0.2, 0.8], false)  # normals [+1, -1]
    eq = FrontEquation(; terms=(NormalMotionTerm(0.5),), front=front, integrator=RK2())
    integrate!(eq, 0.2; dt=0.2)

    x = current_state(eq).mesh.x
    @test x[1] > 0.2
    @test x[2] < 0.8
    @test (x[2] - x[1]) < 0.6
end

@testset "PointFront1D – two-marker advection preserves width" begin
    front = PointFront1D([0.2, 0.8], true)
    eq = FrontEquation(; terms=(AdvectionTerm(1.5),), front=front, integrator=RK3())
    integrate!(eq, 0.25; dt=0.25)
    x = current_state(eq).mesh.x
    @test x[1] ≈ 0.2 + 1.5 * 0.25 atol=1e-12
    @test x[2] ≈ 0.8 + 1.5 * 0.25 atol=1e-12
    @test (x[2] - x[1]) ≈ 0.6 atol=1e-12
end

@testset "PointFront1D – repair and invariants" begin
    s1 = FrontState(PointFront1D([0.3], true))
    x_before = copy(s1.mesh.x)
    repair_front!(s1)
    @test s1.mesh.x == x_before

    s2 = FrontState(PointFront1D([0.2, 0.8], true))
    repair_front!(s2; xminsep=0.5)
    @test s2.mesh.x[1] < s2.mesh.x[2]

    @test_throws ErrorException repair_front!(s2; xminsep=0.7)

    s3 = FrontState(PointFront1D([0.2, 0.8], true))
    repair_front!(s3; xminsep=0.7, clamp_small_gap=true)
    @test s3.mesh.x[2] - s3.mesh.x[1] ≈ 0.7 atol=1e-12

    # Crossed state during update must fail and never be accepted.
    eq_cross = FrontEquation(;
        terms=(NormalMotionTerm(2.0),),
        front=PointFront1D([0.2, 0.8], false),   # normals [+1, -1]
        integrator=ForwardEuler(),
    )
    @test_throws ArgumentError integrate!(eq_cross, 0.2; dt=0.2)
end

@testset "PointFront1D – curvature term rejected" begin
    eq = FrontEquation(;
        terms=(CurvatureMotionTerm(1.0),),
        front=PointFront1D([0.3], true),
        integrator=RK2(),
    )
    @test_throws ErrorException integrate!(eq, 0.1; dt=0.1)
end
