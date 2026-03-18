# topology_patch.jl – Event localization and Cartesian patch construction.

struct CartesianPatch{T,D}
    xmin::SVector{D,T}
    xmax::SVector{D,T}
    h::T
    dims::NTuple{D,Int}
    centers
end

struct EventPatch
    component_ids::Vector{Int}
    bbox
    margin_bbox
    anchor_data
    local_mesh_data
    patch_grid
end

function _bbox_union(boxes)
    isempty(boxes) && error("_bbox_union: empty box list.")
    mins = collect(first(boxes)[1])
    maxs = collect(first(boxes)[2])
    D = length(mins)
    for b in boxes
        lo, hi = b
        for d in 1:D
            mins[d] = min(mins[d], lo[d])
            maxs[d] = max(maxs[d], hi[d])
        end
    end
    return (SVector{D,Float64}(mins), SVector{D,Float64}(maxs))
end

function _expand_bbox(bbox, margin::Real)
    lo, hi = bbox
    D = length(lo)
    margin_v = SVector{D,Float64}(ntuple(_ -> Float64(margin), D))
    return (lo .- margin_v, hi .+ margin_v)
end

function make_patch_grid(bbox, h::Real)
    lo, hi = bbox
    D = length(lo)
    hh = max(Float64(h), eps(Float64))

    ext = hi .- lo
    dims = ntuple(d -> max(2, ceil(Int, ext[d] / hh)), D)

    if D == 2
        nx, ny = dims
        centers = Matrix{SVector{2,Float64}}(undef, nx, ny)
        for j in 1:ny, i in 1:nx
            centers[i, j] = SVector(lo[1] + (i - 0.5) * hh,
                                    lo[2] + (j - 0.5) * hh)
        end
    elseif D == 3
        nx, ny, nz = dims
        centers = Array{SVector{3,Float64},3}(undef, nx, ny, nz)
        for k in 1:nz, j in 1:ny, i in 1:nx
            centers[i, j, k] = SVector(lo[1] + (i - 0.5) * hh,
                                       lo[2] + (j - 0.5) * hh,
                                       lo[3] + (k - 0.5) * hh)
        end
    else
        error("make_patch_grid: only 2D/3D patches are supported.")
    end

    return CartesianPatch{Float64,D}(lo, hi, hh, dims, centers)
end

patch_cell_centers(grid::CartesianPatch) = grid.centers

function _component_bbox(comp::FrontComponentState)
    pts = comp.mesh.points
    D = length(first(pts))
    mins = collect(first(pts))
    maxs = collect(first(pts))
    for p in pts
        for d in 1:D
            mins[d] = min(mins[d], p[d])
            maxs[d] = max(maxs[d], p[d])
        end
    end
    return (SVector{D,Float64}(mins), SVector{D,Float64}(maxs))
end

function extract_event_patch(state::MultiFrontState,
                             candidate::TopologyCandidate,
                             handler::LocalCartesianTopologyHandler)
    comps = [component(state, cid) for cid in candidate.component_ids]
    boxes = [_component_bbox(comp) for comp in comps]
    base_bbox = if handler.reconstruct_scope == :whole_component
        _bbox_union(boxes)
    else
        candidate.bbox === nothing ? _bbox_union(boxes) : candidate.bbox
    end

    margin = handler.patch_margin_factor * max(candidate.h_local, eps(Float64))
    margin_bbox = _expand_bbox(base_bbox, margin)

    h_patch = handler.patch_h_factor * max(candidate.h_local, eps(Float64))
    grid = make_patch_grid(margin_bbox, h_patch)

    anchor_data = Dict{Symbol,Any}(
        :scope => handler.reconstruct_scope,
        :mode => (handler.reconstruct_scope == :local_patch ? :local_patch : :whole_component),
    )

    local_mesh_data = [comp.mesh for comp in comps]

    return EventPatch(candidate.component_ids, base_bbox, margin_bbox,
                      anchor_data, local_mesh_data, grid)
end

function extract_event_patch(state::FrontState,
                             candidate::TopologyCandidate,
                             handler::LocalCartesianTopologyHandler)
    return extract_event_patch(MultiFrontState(state), candidate, handler)
end
