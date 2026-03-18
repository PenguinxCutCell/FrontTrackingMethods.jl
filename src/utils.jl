# utils.jl – Internal helpers.

"""
    front_spacing(state::FrontState) -> Float64

Return the minimum edge length in the current front mesh.
Used as the reference length scale for CFL estimates.
"""
function front_spacing(state::FrontState)
    geom = state.geom
    lengths = if geom isa CurveGeometry
        geom.edge_lengths
    elseif geom isa SurfaceGeometry
        geom.edge_lengths
    else
        error("front_spacing: unsupported geometry type $(typeof(geom))")
    end
    isempty(lengths) && return one(eltype(lengths))
    h = minimum(lengths)
    return max(h, eps(eltype(lengths)))
end

"""
    _zeros_like_points(mesh) -> Vector{SVector}

Allocate a zero-filled velocity buffer with the same SVector type as mesh.points.
"""
function _zeros_like_points(mesh::CurveMesh{T}) where {T}
    return [zero(SVector{2,T}) for _ in eachindex(mesh.points)]
end

function _zeros_like_points(mesh::SurfaceMesh{T}) where {T}
    return [zero(SVector{3,T}) for _ in eachindex(mesh.points)]
end

"""
    _advance_coordinates(mesh, x0, V, dt) -> new_points

Compute new point positions: p_a = x0[a] + dt * V[a].
"""
function _advance_coordinates(mesh, x0::Vector{S}, V::Vector{S}, dt::Real) where {S}
    T = eltype(S)
    dt_T = T(dt)
    return [x0[a] + dt_T * V[a] for a in eachindex(x0)]
end

"""
    _replace_mesh(mesh::CurveMesh, new_points) -> CurveMesh

Create a new CurveMesh with `new_points` but the same connectivity.
"""
_replace_mesh(mesh::CurveMesh, new_points) = CurveMesh(new_points, mesh.edges)

"""
    _replace_mesh(mesh::SurfaceMesh, new_points) -> SurfaceMesh

Create a new SurfaceMesh with `new_points` but the same connectivity.
"""
_replace_mesh(mesh::SurfaceMesh, new_points) = SurfaceMesh(new_points, mesh.faces)
