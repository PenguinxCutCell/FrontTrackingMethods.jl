# test_utils.jl – Shared test helpers for FrontTrackingMethods test suite.
#
# Provides geometry constructors and assertion helpers used across multiple
# test files.

module TestUtils

using StaticArrays
using LinearAlgebra
using Test
using FrontIntrinsicOps

export make_circle_curve, make_ellipse_curve, make_sphere_surface
export curve_centroid, surface_centroid, curve_area, surface_volume
export edge_length_stats, triangle_quality_stats
export assert_no_nan, assert_closed_curve, assert_positive_orientation_curve
export assert_surface_is_reasonable

# ─────────────────────────────────────────────────────────────────────────────
# Geometry constructors
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_circle_curve(; R=1.0, N=128, center=SVector(0.0,0.0)) -> CurveMesh

Create a regular N-gon approximation of a circle.
"""
function make_circle_curve(;
    R      :: Float64            = 1.0,
    N      :: Int                = 128,
    center :: SVector{2,Float64} = SVector(0.0, 0.0),
) :: CurveMesh{Float64}
    T = Float64
    pts   = [center + SVector{2,T}(R*cos(2T(π)*k/N), R*sin(2T(π)*k/N)) for k in 0:N-1]
    edges = [SVector{2,Int}(k, mod1(k+1, N)) for k in 1:N]
    return CurveMesh{T}(pts, edges)
end

"""
    make_ellipse_curve(; a=1.0, b=0.5, N=128, center=SVector(0.0,0.0)) -> CurveMesh
"""
function make_ellipse_curve(;
    a      :: Float64            = 1.0,
    b      :: Float64            = 0.5,
    N      :: Int                = 128,
    center :: SVector{2,Float64} = SVector(0.0, 0.0),
) :: CurveMesh{Float64}
    T = Float64
    pts   = [center + SVector{2,T}(a*cos(2T(π)*k/N), b*sin(2T(π)*k/N)) for k in 0:N-1]
    edges = [SVector{2,Int}(k, mod1(k+1, N)) for k in 1:N]
    return CurveMesh{T}(pts, edges)
end

"""
    make_sphere_surface(; R=1.0, refinement=2, center=SVector(0.0,0.0,0.0)) -> SurfaceMesh
"""
function make_sphere_surface(;
    R          :: Float64            = 1.0,
    refinement :: Int                = 2,
    center     :: SVector{3,Float64} = SVector(0.0, 0.0, 0.0),
) :: SurfaceMesh{Float64}
    mesh = FrontIntrinsicOps.generate_icosphere(R, refinement)
    if !iszero(center)
        new_pts = [p + center for p in mesh.points]
        return SurfaceMesh{Float64}(new_pts, mesh.faces)
    end
    return mesh
end

# ─────────────────────────────────────────────────────────────────────────────
# Shape metrics on meshes
# ─────────────────────────────────────────────────────────────────────────────

"""
    curve_centroid(mesh::CurveMesh) -> SVector{2}

Return the vertex centroid of the curve.
"""
curve_centroid(mesh::CurveMesh) = sum(mesh.points) / length(mesh.points)

"""
    surface_centroid(mesh::SurfaceMesh) -> SVector{3}

Return the vertex centroid of the surface.
"""
surface_centroid(mesh::SurfaceMesh) = sum(mesh.points) / length(mesh.points)

"""
    curve_area(mesh::CurveMesh) -> Float64

Return the enclosed area of the curve using FrontIntrinsicOps.
"""
curve_area(mesh::CurveMesh) = FrontIntrinsicOps.enclosed_measure(mesh)

"""
    surface_volume(mesh::SurfaceMesh) -> Float64

Return the enclosed volume of the surface using FrontIntrinsicOps.
"""
surface_volume(mesh::SurfaceMesh) = FrontIntrinsicOps.enclosed_measure(mesh)

"""
    edge_length_stats(mesh) -> (min, max, mean, spread)

Return (min, max, mean, spread) of edge lengths.
"""
function edge_length_stats(mesh)
    geom = FrontIntrinsicOps.compute_geometry(mesh)
    ls   = geom.edge_lengths
    mn   = minimum(ls)
    mx   = maximum(ls)
    me   = sum(ls) / length(ls)
    sp   = mx / max(mn, eps(Float64))
    return (min=mn, max=mx, mean=me, spread=sp)
end

"""
    triangle_quality_stats(mesh::SurfaceMesh) -> (min_angle, mean_aspect)

Compute minimum angle and mean aspect ratio of triangles.
"""
function triangle_quality_stats(mesh::SurfaceMesh)
    pts    = mesh.points
    faces  = mesh.faces
    T      = Float64
    min_angle  = T(π)
    sum_aspect = zero(T)
    n = length(faces)
    n == 0 && return (min_angle=T(π), mean_aspect=one(T))

    for f in faces
        p1 = pts[f[1]]
        p2 = pts[f[2]]
        p3 = pts[f[3]]
        e1 = norm(p2 - p1)
        e2 = norm(p3 - p2)
        e3 = norm(p1 - p3)
        a  = acos(clamp(dot(p2-p1, p3-p1) / (max(e1*norm(p3-p1), eps(T))), -one(T), one(T)))
        b  = acos(clamp(dot(p1-p2, p3-p2) / (max(e1*e2, eps(T))), -one(T), one(T)))
        c  = T(π) - a - b
        min_angle  = min(min_angle, a, b, c)
        # Aspect ratio: circumradius / (2 * inradius)
        s = (e1 + e2 + e3) / 2
        area = sqrt(max(s*(s-e1)*(s-e2)*(s-e3), zero(T)))
        if area > eps(T)
            r_in   = area / s
            r_circ = e1 * e2 * e3 / (4 * area)
            sum_aspect += r_circ / (2 * r_in)
        else
            sum_aspect += T(Inf)
        end
    end
    return (min_angle=min_angle, mean_aspect=sum_aspect / n)
end

# ─────────────────────────────────────────────────────────────────────────────
# Assertion helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    assert_no_nan(mesh)

Assert that no vertex has NaN or Inf coordinates.
"""
function assert_no_nan(mesh)
    for p in mesh.points
        @test all(isfinite, p)
    end
end

"""
    assert_closed_curve(mesh::CurveMesh)

Assert that the curve edges form a closed loop.
"""
function assert_closed_curve(mesh::CurveMesh)
    @test is_closed(mesh)
end

"""
    assert_positive_orientation_curve(mesh::CurveMesh)

Assert that the curve has positive (CCW) orientation.
"""
function assert_positive_orientation_curve(mesh::CurveMesh)
    A = curve_area(mesh)
    @test A > 0
end

"""
    assert_surface_is_reasonable(mesh::SurfaceMesh; min_area_tol=1e-10)

Assert that the surface has no degenerate triangles and that the enclosed
volume is positive.
"""
function assert_surface_is_reasonable(mesh::SurfaceMesh; min_area_tol::Float64=1e-10)
    geom = FrontIntrinsicOps.compute_geometry(mesh)
    for fa in geom.face_areas
        @test fa >= min_area_tol
    end
    V = surface_volume(mesh)
    @test V > 0
end

end  # module TestUtils
