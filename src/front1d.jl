# front1d.jl – Minimal helpers for PointFront1D support.

struct PointFront1DGeometry{T<:AbstractFloat}
    marker_positions :: Vector{T}
    vertex_normals   :: Vector{T}
    edge_lengths     :: Vector{T}
end

function compute_geometry(front::PointFront1D{T}) where {T<:Real}
    U = float(T)
    x = U.(front.x)
    nrm = FrontIntrinsicOps.interface_normals(front)

    edge_lengths = if length(x) == 2
        U[x[2] - x[1]]
    else
        U[]
    end

    return PointFront1DGeometry{U}(x, nrm, edge_lengths)
end
