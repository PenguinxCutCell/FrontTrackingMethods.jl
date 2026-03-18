# topology_rasterize.jl – Rasterize local front data onto Cartesian indicator patches.

@inline function _point_in_polygon_2d(p::SVector{2,Float64}, mesh::CurveMesh)
    x, y = p
    inside = false
    for e in mesh.edges
        p1 = mesh.points[e[1]]
        p2 = mesh.points[e[2]]
        x1, y1 = p1[1], p1[2]
        x2, y2 = p2[1], p2[2]

        intersects = ((y1 > y) != (y2 > y)) &&
                     (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)
        intersects && (inside = !inside)
    end
    return inside
end

@inline function _ray_intersects_triangle_x(p::SVector{3,Float64}, a, b, c)
    dir = SVector(1.0, 0.0, 0.0)
    epsv = 1e-12
    e1 = b - a
    e2 = c - a
    h = dir × e2
    det = dot(e1, h)
    abs(det) < epsv && return false

    inv_det = 1.0 / det
    s = p - a
    u = inv_det * dot(s, h)
    (u < -epsv || u > 1 + epsv) && return false

    q = s × e1
    v = inv_det * dot(dir, q)
    (v < -epsv || (u + v) > 1 + epsv) && return false

    t = inv_det * dot(e2, q)
    return t > epsv
end

function rasterize_indicator!(χ::AbstractArray{<:Real}, patch::EventPatch, local_front_data=patch.local_mesh_data)
    grid = patch.patch_grid
    D = length(grid.dims)

    if D == 2
        nx, ny = grid.dims
        @assert size(χ) == (nx, ny)
        curves = [m for m in local_front_data if m isa CurveMesh]

        for j in 1:ny, i in 1:nx
            c = grid.centers[i, j]
            inside = false
            for curve in curves
                if _point_in_polygon_2d(c, curve)
                    inside = true
                    break
                end
            end
            χ[i, j] = inside ? one(eltype(χ)) : zero(eltype(χ))
        end
    elseif D == 3
        nx, ny, nz = grid.dims
        @assert size(χ) == (nx, ny, nz)
        surfaces = [m for m in local_front_data if m isa SurfaceMesh]

        for k in 1:nz, j in 1:ny, i in 1:nx
            p = grid.centers[i, j, k]
            inside = false
            for surf in surfaces
                nint = 0
                for f in surf.faces
                    a = surf.points[f[1]]
                    b = surf.points[f[2]]
                    c = surf.points[f[3]]
                    _ray_intersects_triangle_x(p, a, b, c) && (nint += 1)
                end
                if isodd(nint)
                    inside = true
                    break
                end
            end
            χ[i, j, k] = inside ? one(eltype(χ)) : zero(eltype(χ))
        end
    else
        error("rasterize_indicator!: only 2D/3D supported.")
    end

    return χ
end
