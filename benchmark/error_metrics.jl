using FrontTrackingMethods
using FrontIntrinsicOps
using StaticArrays
using LinearAlgebra

function mean_edge_size(mesh::CurveMesh)
    geom = compute_geometry(mesh)
    return sum(geom.edge_lengths) / length(geom.edge_lengths)
end

function area_error(mesh_final::CurveMesh, mesh_init::CurveMesh)
    A0 = enclosed_measure(mesh_init)
    AT = enclosed_measure(mesh_final)
    return abs(AT - A0) / max(abs(A0), eps(Float64))
end

function circle_shape_error(mesh::CurveMesh; center::SVector{2,Float64}, R::Float64)
    return maximum(abs(norm(p - center) - R) for p in mesh.points)
end

function _ordered_polygon(mesh::CurveMesh)
    ord = curve_vertex_order(mesh)
    return mesh.points[ord]
end

function _point_in_polygon(p::SVector{2,Float64}, poly::Vector{SVector{2,Float64}})
    inside = false
    n = length(poly)
    j = n
    @inbounds for i in 1:n
        xi, yi = poly[i][1], poly[i][2]
        xj, yj = poly[j][1], poly[j][2]
        intersects = ((yi > p[2]) != (yj > p[2])) &&
                     (p[1] < (xj - xi) * (p[2] - yi) / (yj - yi + eps(Float64)) + xi)
        intersects && (inside = !inside)
        j = i
    end
    return inside
end

function symmetric_difference_area(
    meshA::CurveMesh,
    meshB::CurveMesh;
    xmin::Float64=0.0,
    xmax::Float64=1.0,
    ymin::Float64=0.0,
    ymax::Float64=1.0,
    nx::Int=600,
    ny::Int=600,
)
    polyA = _ordered_polygon(meshA)
    polyB = _ordered_polygon(meshB)

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    cellA = dx * dy

    xor_count = 0
    for j in 1:ny
        y = ymin + (j - 0.5) * dy
        for i in 1:nx
            x = xmin + (i - 0.5) * dx
            p = SVector{2,Float64}(x, y)
            inA = _point_in_polygon(p, polyA)
            inB = _point_in_polygon(p, polyB)
            (inA ⊻ inB) && (xor_count += 1)
        end
    end
    return xor_count * cellA
end

function benchmark_errors(mesh_init::CurveMesh, mesh_final::CurveMesh;
                          circle_center::Union{Nothing,SVector{2,Float64}}=nothing,
                          circle_radius::Union{Nothing,Float64}=nothing,
                          sym_grid::Int=600)
    E_area = area_error(mesh_final, mesh_init)
    E_sym  = symmetric_difference_area(mesh_init, mesh_final; nx=sym_grid, ny=sym_grid)
    E_shape = if circle_center === nothing || circle_radius === nothing
        NaN
    else
        circle_shape_error(mesh_final; center=circle_center, R=circle_radius)
    end
    return (E_area=E_area, E_shape=E_shape, E_sym=E_sym)
end
