# topology_reconstruct.jl – Reconstruct front meshes from patch indicators.

@inline _ptkey2(p::SVector{2,Float64}) = (round(Int, p[1] * 1e10), round(Int, p[2] * 1e10))
@inline _ptkey3(p::SVector{3,Float64}) = (round(Int, p[1] * 1e10), round(Int, p[2] * 1e10), round(Int, p[3] * 1e10))

function _clean_loop_points(loop::Vector{SVector{2,Float64}}, h::Float64)
    isempty(loop) && return loop
    tol = max(1e-12, 1e-6 * h)
    out = SVector{2,Float64}[]
    push!(out, loop[1])
    for i in 2:length(loop)
        norm(loop[i] - out[end]) > tol && push!(out, loop[i])
    end
    if length(out) >= 2 && norm(out[1] - out[end]) <= tol
        pop!(out)
    end
    return out
end

function _boundary_loops_from_occupancy_2d(χ, grid::CartesianPatch{T,2}) where {T}
    nx, ny = grid.dims
    lo = grid.xmin
    h = grid.h

    PointI2 = NTuple{2,Int}
    EdgeI2 = Tuple{PointI2,PointI2}

    edge_store = Dict{EdgeI2,EdgeI2}()

    canon(a, b) = a <= b ? (a, b) : (b, a)

    function toggle_edge!(a, b)
        key = canon(a, b)
        if haskey(edge_store, key)
            delete!(edge_store, key)
        else
            edge_store[key] = (a, b)
        end
    end

    for j in 1:ny, i in 1:nx
        χ[i, j] > 0.5 || continue

        p00 = (i - 1, j - 1)
        p10 = (i,     j - 1)
        p11 = (i,     j)
        p01 = (i - 1, j)

        toggle_edge!(p00, p10)
        toggle_edge!(p10, p11)
        toggle_edge!(p11, p01)
        toggle_edge!(p01, p00)
    end

    adjacency = Dict{PointI2,Vector{PointI2}}()
    for (_, (a, b)) in edge_store
        push!(get!(adjacency, a, PointI2[]), b)
        push!(get!(adjacency, b, PointI2[]), a)
    end

    unused = Set{EdgeI2}(keys(edge_store))
    loops = Vector{Vector{SVector{2,Float64}}}()

    while !isempty(unused)
        ekey = first(unused)
        start, nxt = ekey
        path = PointI2[start, nxt]
        delete!(unused, canon(start, nxt))

        prev = start
        cur = nxt
        while cur != start
            nbrs = get(adjacency, cur, PointI2[])
            isempty(nbrs) && break

            candidate = nothing
            for nb in nbrs
                nb == prev && continue
                ckey = canon(cur, nb)
                if ckey in unused
                    candidate = nb
                    break
                end
            end

            candidate === nothing && break
            nb = candidate::PointI2
            delete!(unused, canon(cur, nb))
            push!(path, nb)
            prev, cur = cur, nb
        end

        if length(path) >= 4
            pts = [SVector(lo[1] + p[1] * h, lo[2] + p[2] * h) for p in path]
            push!(loops, pts)
        end
    end

    return loops
end

function reconstruct_curve_from_patch(χ::AbstractMatrix{<:Real}, patch::EventPatch)
    grid = patch.patch_grid
    loops = _boundary_loops_from_occupancy_2d(χ, grid)
    curves = CurveMesh[]

    for loop in loops
        loop = _clean_loop_points(loop, grid.h)
        length(loop) >= 4 || continue

        area2 = zero(Float64)
        n = length(loop)
        for i in 1:n
            p = loop[i]
            q = loop[mod1(i + 1, n)]
            area2 += p[1] * q[2] - q[1] * p[2]
        end
        abs(area2) < 2 * grid.h^2 && continue

        pts = copy(loop)
        if area2 < 0
            reverse!(pts)
        end
        m = length(pts)
        edges = [SVector{2,Int}(i, mod1(i + 1, m)) for i in 1:m]
        push!(curves, CurveMesh{Float64}(pts, edges))
    end

    return curves
end

function _occupied_components_3d(χ::AbstractArray{<:Real,3})
    nx, ny, nz = size(χ)
    visited = falses(nx, ny, nz)
    comps = Vector{Vector{NTuple{3,Int}}}()

    dirs = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

    for k in 1:nz, j in 1:ny, i in 1:nx
        if visited[i,j,k] || χ[i,j,k] <= 0.5
            continue
        end
        stack = [(i,j,k)]
        visited[i,j,k] = true
        comp = NTuple{3,Int}[]

        while !isempty(stack)
            c = pop!(stack)
            push!(comp, c)
            ci, cj, ck = c
            for (di,dj,dk) in dirs
                ni, nj, nk = ci + di, cj + dj, ck + dk
                if 1 <= ni <= nx && 1 <= nj <= ny && 1 <= nk <= nz && !visited[ni,nj,nk] && χ[ni,nj,nk] > 0.5
                    visited[ni,nj,nk] = true
                    push!(stack, (ni,nj,nk))
                end
            end
        end

        push!(comps, comp)
    end

    return comps
end

function _voxel_surface_mesh(comp_cells::Vector{NTuple{3,Int}}, χ::AbstractArray{<:Real,3}, grid::CartesianPatch{T,3}) where {T}
    nx, ny, nz = size(χ)
    lo = grid.xmin
    h = grid.h

    cellset = Set(comp_cells)
    vmap = Dict{NTuple{3,Int},Int}()
    pts = SVector{3,Float64}[]
    faces = SVector{3,Int}[]

    function vindex(ii, jj, kk)
        key = (ii, jj, kk)
        if haskey(vmap, key)
            return vmap[key]
        end
        p = SVector(lo[1] + (ii - 1) * h,
                    lo[2] + (jj - 1) * h,
                    lo[3] + (kk - 1) * h)
        push!(pts, p)
        idx = length(pts)
        vmap[key] = idx
        return idx
    end

    function add_quad(a,b,c,d)
        push!(faces, SVector{3,Int}(a,b,c))
        push!(faces, SVector{3,Int}(a,c,d))
    end

    dirs = (
        (1,0,0), (-1,0,0),
        (0,1,0), (0,-1,0),
        (0,0,1), (0,0,-1),
    )

    for (i,j,k) in comp_cells
        for (di,dj,dk) in dirs
            ni, nj, nk = i + di, j + dj, k + dk
            exposed = !(1 <= ni <= nx && 1 <= nj <= ny && 1 <= nk <= nz && ((ni,nj,nk) in cellset))
            exposed || continue

            if (di,dj,dk) == (1,0,0)
                v1 = vindex(i+1,j,  k)
                v2 = vindex(i+1,j+1,k)
                v3 = vindex(i+1,j+1,k+1)
                v4 = vindex(i+1,j,  k+1)
            elseif (di,dj,dk) == (-1,0,0)
                v1 = vindex(i,j,  k)
                v2 = vindex(i,j,  k+1)
                v3 = vindex(i,j+1,k+1)
                v4 = vindex(i,j+1,k)
            elseif (di,dj,dk) == (0,1,0)
                v1 = vindex(i,  j+1,k)
                v2 = vindex(i,  j+1,k+1)
                v3 = vindex(i+1,j+1,k+1)
                v4 = vindex(i+1,j+1,k)
            elseif (di,dj,dk) == (0,-1,0)
                v1 = vindex(i,  j,k)
                v2 = vindex(i+1,j,k)
                v3 = vindex(i+1,j,k+1)
                v4 = vindex(i,  j,k+1)
            elseif (di,dj,dk) == (0,0,1)
                v1 = vindex(i,  j,  k+1)
                v2 = vindex(i+1,j,  k+1)
                v3 = vindex(i+1,j+1,k+1)
                v4 = vindex(i,  j+1,k+1)
            else
                v1 = vindex(i,  j,  k)
                v2 = vindex(i,  j+1,k)
                v3 = vindex(i+1,j+1,k)
                v4 = vindex(i+1,j,  k)
            end
            add_quad(v1,v2,v3,v4)
        end
    end

    return SurfaceMesh{Float64}(pts, faces)
end

function reconstruct_surface_from_patch(χ::AbstractArray{<:Real,3}, patch::EventPatch)
    grid = patch.patch_grid
    comps = _occupied_components_3d(χ)
    surfaces = SurfaceMesh[]

    for comp in comps
        length(comp) < 2 && continue
        mesh = _voxel_surface_mesh(comp, χ, grid)
        isempty(mesh.faces) && continue
        push!(surfaces, mesh)
    end

    return surfaces
end
