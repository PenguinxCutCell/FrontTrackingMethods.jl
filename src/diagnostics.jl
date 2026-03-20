# diagnostics.jl – Front geometry quality metrics and diagnostics.
#
# These helpers are useful for monitoring simulation quality, checking
# convergence, and detecting degenerate configurations.

"""
    min_marker_position(front::PointFront1D) -> Real
"""
min_marker_position(front::PointFront1D) = minimum(front.x)
min_marker_position(state::FrontState{<:PointFront1D}) = min_marker_position(state.mesh)

"""
    max_marker_position(front::PointFront1D) -> Real
"""
max_marker_position(front::PointFront1D) = maximum(front.x)
max_marker_position(state::FrontState{<:PointFront1D}) = max_marker_position(state.mesh)

"""
    marker_gap(front::PointFront1D) -> Real

Return `x[2] - x[1]` for two-marker fronts, otherwise `0`.
"""
marker_gap(front::PointFront1D) = length(front.x) == 2 ? front.x[2] - front.x[1] : zero(eltype(front.x))
marker_gap(state::FrontState{<:PointFront1D}) = marker_gap(state.mesh)

"""
    is_valid_front(front::PointFront1D) -> Bool

Check PointFront1D invariants.
"""
function is_valid_front(front::PointFront1D)
    return try
        FrontIntrinsicOps.check(front)
    catch
        false
    end
end

is_valid_front(state::FrontState{<:PointFront1D}) = is_valid_front(state.mesh)

"""
    min_edge_length(state::FrontState) -> Float64

Return the length of the shortest edge in the current front mesh.
"""
function min_edge_length(state::FrontState)
    return minimum(state.geom.edge_lengths)
end
min_edge_length(state::FrontState{<:PointFront1D}) = marker_gap(state)

"""
    max_edge_length(state::FrontState) -> Float64

Return the length of the longest edge in the current front mesh.
"""
function max_edge_length(state::FrontState)
    return maximum(state.geom.edge_lengths)
end
max_edge_length(state::FrontState{<:PointFront1D}) = marker_gap(state)

"""
    mean_edge_length(state::FrontState) -> Float64

Return the mean edge length in the current front mesh.
"""
function mean_edge_length(state::FrontState)
    ls = state.geom.edge_lengths
    return sum(ls) / length(ls)
end
mean_edge_length(state::FrontState{<:PointFront1D}) = marker_gap(state)

"""
    edge_length_spread(state::FrontState) -> Float64

Return the ratio max_edge_length / min_edge_length.
A well-distributed front has this close to 1.
"""
function edge_length_spread(state::FrontState)
    ls = state.geom.edge_lengths
    return maximum(ls) / max(minimum(ls), eps(eltype(ls)))
end
edge_length_spread(state::FrontState{<:PointFront1D}) = 1.0

"""
    front_enclosed_measure(state::FrontState) -> Float64

Compute the enclosed measure (area for curves, volume for surfaces)
using FrontIntrinsicOps.

Returns the shoelace area for a closed `CurveMesh` or the divergence-
theorem volume for a closed `SurfaceMesh`.
"""
function front_enclosed_measure(state::FrontState)
    return enclosed_measure(state.mesh)
end

"""
    front_measure(state::FrontState) -> Float64

Return the total intrinsic measure of the front (arc-length for curves,
surface area for surfaces).
"""
function front_measure(state::FrontState)
    return measure(state.mesh, state.geom)
end

"""
    front_centroid(state::FrontState) -> SVector

Compute the centroid of the front vertices (simple average of all vertex
positions).
"""
function front_centroid(state::FrontState)
    pts = state.mesh.points
    isempty(pts) && error("front_centroid: empty mesh.")
    return sum(pts) / length(pts)
end

"""
    normal_tangential_decomposition(v::SVector{2}, n::SVector{2})
        -> (v_n, v_t)

Decompose a 2-D vector `v` into its normal (scalar) and tangential (SVector)
components with respect to normal `n`.
"""
function normal_tangential_decomposition(v::SVector{2,T}, n::SVector{2,T}) where {T}
    vn = dot(v, n)
    vt = v - vn * n
    return vn, vt
end

"""
    normal_tangential_decomposition(v::SVector{3}, n::SVector{3})
        -> (v_n, v_t)

Decompose a 3-D vector `v` into its normal (scalar) and tangential (SVector)
components with respect to normal `n`.
"""
function normal_tangential_decomposition(v::SVector{3,T}, n::SVector{3,T}) where {T}
    vn = dot(v, n)
    vt = v - vn * n
    return vn, vt
end

"""
    check_front_validity(state::FrontState; warn=true) -> Bool

Perform basic sanity checks on the front:
- No NaN or Inf in vertex positions.
- No zero-length edges.
- Closure (if the mesh is expected to be closed).

Returns `true` if all checks pass.  If `warn=true`, prints warnings instead
of throwing errors.
"""
function check_front_validity(state::FrontState; warn::Bool=true)
    ok    = true
    pts   = state.mesh.points
    geom  = state.geom

    # NaN / Inf check
    for (vi, p) in enumerate(pts)
        if any(!isfinite, p)
            msg = "check_front_validity: vertex $vi has non-finite position: $p"
            warn ? (@warn msg) : error(msg)
            ok = false
        end
    end

    # Zero-length edge check
    for (ei, len) in enumerate(geom.edge_lengths)
        if len < eps(eltype(geom.edge_lengths))
            msg = "check_front_validity: edge $ei has near-zero length: $len"
            warn ? (@warn msg) : error(msg)
            ok = false
        end
    end

    return ok
end

function check_front_validity(state::FrontState{<:PointFront1D}; warn::Bool=true)
    x = state.mesh.x
    ok = true
    if any(!isfinite, x)
        msg = "check_front_validity: PointFront1D has non-finite marker coordinates."
        warn ? (@warn msg) : error(msg)
        ok = false
    end
    if length(x) == 2 && !(x[1] < x[2])
        msg = "check_front_validity: PointFront1D requires strict ordering x[1] < x[2]."
        warn ? (@warn msg) : error(msg)
        ok = false
    end
    return ok
end

"""
    surface_edge_length_stats(mesh::SurfaceMesh)

Return edge-length statistics for a surface mesh as
`(min, max, mean, std, ratio)` where `ratio = max / min`.
"""
function surface_edge_length_stats(mesh::SurfaceMesh)
    topo = build_topology(mesh)
    pts  = mesh.points
    T    = eltype(eltype(pts))

    lengths = T[]
    sizehint!(lengths, length(topo.edges))
    for e in topo.edges
        push!(lengths, norm(pts[e[2]] - pts[e[1]]))
    end
    isempty(lengths) && return (min=zero(T), max=zero(T), mean=zero(T), std=zero(T), ratio=one(T))

    mn = minimum(lengths)
    mx = maximum(lengths)
    μ  = sum(lengths) / length(lengths)
    var = sum((l - μ)^2 for l in lengths) / length(lengths)
    σ  = sqrt(var)
    return (min=mn, max=mx, mean=μ, std=σ, ratio=mx / max(mn, eps(T)))
end

"""
    surface_triangle_area_stats(mesh::SurfaceMesh)

Return triangle-area statistics for a surface mesh as `(min, max, mean)`.
"""
function surface_triangle_area_stats(mesh::SurfaceMesh)
    pts = mesh.points
    T   = eltype(eltype(pts))

    if isempty(mesh.faces)
        return (min=zero(T), max=zero(T), mean=zero(T))
    end

    amin = T(Inf)
    amax = zero(T)
    asum = zero(T)
    for f in mesh.faces
        p1, p2, p3 = pts[f[1]], pts[f[2]], pts[f[3]]
        A = norm((p2 - p1) × (p3 - p1)) / 2
        amin = min(amin, A)
        amax = max(amax, A)
        asum += A
    end
    return (min=amin, max=amax, mean=asum / length(mesh.faces))
end

"""
    surface_triangle_angle_stats(mesh::SurfaceMesh)

Return triangle-angle statistics as `(min_angle, max_angle, mean_min_angle)` in radians.
`mean_min_angle` is the mean over triangles of each triangle's minimum interior angle.
"""
function surface_triangle_angle_stats(mesh::SurfaceMesh)
    pts = mesh.points
    T   = eltype(eltype(pts))

    if isempty(mesh.faces)
        return (min_angle=zero(T), max_angle=zero(T), mean_min_angle=zero(T))
    end

    minang = T(Inf)
    maxang = zero(T)
    sum_minang = zero(T)

    for f in mesh.faces
        p1, p2, p3 = pts[f[1]], pts[f[2]], pts[f[3]]
        e12 = norm(p2 - p1)
        e23 = norm(p3 - p2)
        e31 = norm(p1 - p3)

        c1 = clamp(dot(p2 - p1, p3 - p1) / max(e12 * e31, eps(T)), -one(T), one(T))
        c2 = clamp(dot(p1 - p2, p3 - p2) / max(e12 * e23, eps(T)), -one(T), one(T))
        a1 = acos(c1)
        a2 = acos(c2)
        a3 = T(π) - a1 - a2

        tri_min = min(a1, a2, a3)
        tri_max = max(a1, a2, a3)
        minang = min(minang, tri_min)
        maxang = max(maxang, tri_max)
        sum_minang += tri_min
    end

    return (min_angle=minang, max_angle=maxang, mean_min_angle=sum_minang / length(mesh.faces))
end

"""
    surface_aspect_ratio_stats(mesh::SurfaceMesh)

Return triangle aspect-ratio statistics `(min, max, mean)` using the stable metric

`aspect = (lmax * perimeter) / (4 * sqrt(3) * area)`.

This metric is dimensionless and equals 1 for an equilateral triangle.
"""
function surface_aspect_ratio_stats(mesh::SurfaceMesh)
    pts = mesh.points
    T   = eltype(eltype(pts))

    if isempty(mesh.faces)
        return (min=one(T), max=one(T), mean=one(T))
    end

    amin = T(Inf)
    amax = zero(T)
    asum = zero(T)
    c_eq = T(4 * sqrt(3.0))

    for f in mesh.faces
        p1, p2, p3 = pts[f[1]], pts[f[2]], pts[f[3]]
        l1 = norm(p2 - p1)
        l2 = norm(p3 - p2)
        l3 = norm(p1 - p3)
        lmax = max(l1, l2, l3)
        perim = l1 + l2 + l3
        area = norm((p2 - p1) × (p3 - p1)) / 2
        aspect = area > eps(T) ? (lmax * perim) / (c_eq * area) : T(Inf)
        amin = min(amin, aspect)
        amax = max(amax, aspect)
        asum += aspect
    end

    return (min=amin, max=amax, mean=asum / length(mesh.faces))
end

"""
    surface_degenerate_fraction(mesh::SurfaceMesh; atol=1e-12)

Return the fraction of triangles with area `<= atol`.
"""
function surface_degenerate_fraction(mesh::SurfaceMesh; atol::Real=1e-12)
    isempty(mesh.faces) && return 0.0
    pts = mesh.points
    bad = 0
    for f in mesh.faces
        p1, p2, p3 = pts[f[1]], pts[f[2]], pts[f[3]]
        A = norm((p2 - p1) × (p3 - p1)) / 2
        A <= atol && (bad += 1)
    end
    return bad / length(mesh.faces)
end

"""
    surface_quality_summary(mesh::SurfaceMesh; degenerate_atol=1e-12)

Return a compact named-tuple summary of surface quality metrics:
edge, triangle area, angle, aspect ratio, and degenerate fraction.
"""
function surface_quality_summary(mesh::SurfaceMesh; degenerate_atol::Real=1e-12)
    edge   = surface_edge_length_stats(mesh)
    area   = surface_triangle_area_stats(mesh)
    angle  = surface_triangle_angle_stats(mesh)
    aspect = surface_aspect_ratio_stats(mesh)
    degen  = surface_degenerate_fraction(mesh; atol=degenerate_atol)
    return (
        edge=edge,
        area=area,
        angle=angle,
        aspect=aspect,
        degenerate_fraction=degen,
    )
end

"""
    surface_normal_consistency(mesh::SurfaceMesh, geom::SurfaceGeometry)

Return `(inconsistent_fraction, mean_alignment, min_alignment)` based on
alignment of each face normal with the average of the corresponding
vertex normals.
"""
function surface_normal_consistency(mesh::SurfaceMesh, geom::SurfaceGeometry)
    isempty(mesh.faces) && return (inconsistent_fraction=0.0, mean_alignment=1.0, min_alignment=1.0)

    pts = mesh.points
    vn  = geom.vertex_normals
    bad = 0
    min_align = 1.0
    sum_align = 0.0

    for f in mesh.faces
        p1, p2, p3 = pts[f[1]], pts[f[2]], pts[f[3]]
        n = (p2 - p1) × (p3 - p1)
        nn = norm(n)
        nn <= eps(eltype(eltype(pts))) && continue
        nf = n / nn

        vavg = vn[f[1]] + vn[f[2]] + vn[f[3]]
        vv = norm(vavg)
        vv <= eps(eltype(eltype(pts))) && continue
        vdir = vavg / vv

        a = dot(nf, vdir)
        a < 0 && (bad += 1)
        min_align = min(min_align, Float64(a))
        sum_align += Float64(a)
    end

    nfaces = max(length(mesh.faces), 1)
    return (
        inconsistent_fraction=bad / nfaces,
        mean_alignment=sum_align / nfaces,
        min_alignment=min_align,
    )
end
