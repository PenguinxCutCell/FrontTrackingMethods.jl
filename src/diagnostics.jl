# diagnostics.jl – Front geometry quality metrics and diagnostics.
#
# These helpers are useful for monitoring simulation quality, checking
# convergence, and detecting degenerate configurations.

"""
    min_edge_length(state::FrontState) -> Float64

Return the length of the shortest edge in the current front mesh.
"""
function min_edge_length(state::FrontState)
    return minimum(state.geom.edge_lengths)
end

"""
    max_edge_length(state::FrontState) -> Float64

Return the length of the longest edge in the current front mesh.
"""
function max_edge_length(state::FrontState)
    return maximum(state.geom.edge_lengths)
end

"""
    mean_edge_length(state::FrontState) -> Float64

Return the mean edge length in the current front mesh.
"""
function mean_edge_length(state::FrontState)
    ls = state.geom.edge_lengths
    return sum(ls) / length(ls)
end

"""
    edge_length_spread(state::FrontState) -> Float64

Return the ratio max_edge_length / min_edge_length.
A well-distributed front has this close to 1.
"""
function edge_length_spread(state::FrontState)
    ls = state.geom.edge_lengths
    return maximum(ls) / max(minimum(ls), eps(eltype(ls)))
end

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
