# topology_detection.jl – Conservative geometric topology-event candidate detection.

@inline function _bbox_from_points(points)
    D = length(first(points))
    mins = collect(points[1])
    maxs = collect(points[1])
    @inbounds for p in points
        for d in 1:D
            mins[d] = min(mins[d], p[d])
            maxs[d] = max(maxs[d], p[d])
        end
    end
    return (SVector{D,Float64}(mins), SVector{D,Float64}(maxs))
end

@inline _pair_score(distance::Real, h_local::Real) = h_local / max(distance, eps(Float64))

@inline function _is_local_neighbor(i::Int, j::Int, n::Int, span::Int, closed::Bool)
    d = abs(i - j)
    if d <= span
        return true
    end
    if closed
        return (n - d) <= span
    end
    return false
end

function local_length_scale(component::FrontComponentState, _feature_ids...)
    lengths = component.geom.edge_lengths
    isempty(lengths) && return 1.0
    return max(sum(lengths) / length(lengths), eps(Float64))
end

function _closest_vertex_pair(points_a, points_b)
    dmin = Inf
    pair = (0, 0)
    @inbounds for i in eachindex(points_a)
        pa = points_a[i]
        for j in eachindex(points_b)
            d = norm(pa - points_b[j])
            if d < dmin
                dmin = d
                pair = (i, j)
            end
        end
    end
    return dmin, pair
end

function detect_imminent_merge_2d(componentA::FrontComponentState,
                                  componentB::FrontComponentState,
                                  handler::LocalCartesianTopologyHandler)
    meshA = componentA.mesh
    meshB = componentB.mesh
    (meshA isa CurveMesh && meshB isa CurveMesh) || return nothing

    dmin, pair = _closest_vertex_pair(meshA.points, meshB.points)
    hloc = min(local_length_scale(componentA), local_length_scale(componentB))
    trigger = dmin < handler.d_merge * hloc
    if !trigger
        return nothing
    end
    pa = meshA.points[pair[1]]
    pb = meshB.points[pair[2]]
    bbox = _bbox_from_points((pa, pb))
    return (distance=dmin,
            h_local=hloc,
            feature_pairs=pair,
            bbox=bbox,
            score=_pair_score(dmin, hloc))
end

function detect_imminent_self_merge_2d(component::FrontComponentState,
                                       handler::LocalCartesianTopologyHandler)
    mesh = component.mesh
    mesh isa CurveMesh || return nothing
    n = length(mesh.points)
    n < 6 && return nothing

    closed = is_closed(mesh)
    span = max(3, round(Int, 0.02 * n))
    dmin = Inf
    best = (0, 0)
    @inbounds for i in 1:n
        for j in i+1:n
            if _is_local_neighbor(i, j, n, span, closed)
                continue
            end
            d = norm(mesh.points[i] - mesh.points[j])
            if d < dmin
                dmin = d
                best = (i, j)
            end
        end
    end

    hloc = local_length_scale(component)
    dmin < handler.d_merge * hloc || return nothing

    ni = component.geom.vertex_normals[best[1]]
    nj = component.geom.vertex_normals[best[2]]
    dot(ni, nj) < -0.2 || return nothing

    pa = mesh.points[best[1]]
    pb = mesh.points[best[2]]
    bbox = _bbox_from_points((pa, pb))
    return (distance=dmin,
            h_local=hloc,
            feature_pairs=best,
            bbox=bbox,
            score=_pair_score(dmin, hloc))
end

function detect_imminent_split_2d(component::FrontComponentState,
                                  handler::LocalCartesianTopologyHandler)
    mesh = component.mesh
    mesh isa CurveMesh || return nothing
    n = length(mesh.points)
    n < 6 && return nothing

    closed = is_closed(mesh)
    span = max(3, round(Int, 0.02 * n))
    dmin = Inf
    best = (0, 0)
    @inbounds for i in 1:n
        for j in i+1:n
            if _is_local_neighbor(i, j, n, span, closed)
                continue
            end
            d = norm(mesh.points[i] - mesh.points[j])
            if d < dmin
                dmin = d
                best = (i, j)
            end
        end
    end

    hloc = local_length_scale(component)
    dmin < handler.d_split * hloc || return nothing

    ni = component.geom.vertex_normals[best[1]]
    nj = component.geom.vertex_normals[best[2]]
    dot(ni, nj) < 0.8 || return nothing

    pa = mesh.points[best[1]]
    pb = mesh.points[best[2]]
    bbox = _bbox_from_points((pa, pb))
    return (distance=dmin,
            h_local=hloc,
            feature_pairs=best,
            bbox=bbox,
            score=_pair_score(dmin, hloc))
end

function detect_imminent_merge_3d(componentA::FrontComponentState,
                                  componentB::FrontComponentState,
                                  handler::LocalCartesianTopologyHandler)
    meshA = componentA.mesh
    meshB = componentB.mesh
    (meshA isa SurfaceMesh && meshB isa SurfaceMesh) || return nothing

    dmin, pair = _closest_vertex_pair(meshA.points, meshB.points)
    hloc = min(local_length_scale(componentA), local_length_scale(componentB))
    dmin < handler.d_merge * hloc || return nothing

    pa = meshA.points[pair[1]]
    pb = meshB.points[pair[2]]
    bbox = _bbox_from_points((pa, pb))
    return (distance=dmin,
            h_local=hloc,
            feature_pairs=pair,
            bbox=bbox,
            score=_pair_score(dmin, hloc))
end

function detect_imminent_self_merge_3d(component::FrontComponentState,
                                       handler::LocalCartesianTopologyHandler)
    mesh = component.mesh
    mesh isa SurfaceMesh || return nothing
    n = length(mesh.points)
    n < 10 && return nothing

    dmin = Inf
    best = (0, 0)
    @inbounds for i in 1:n
        for j in i+3:n
            d = norm(mesh.points[i] - mesh.points[j])
            if d < dmin
                dmin = d
                best = (i, j)
            end
        end
    end
    hloc = local_length_scale(component)
    dmin < handler.d_merge * hloc || return nothing

    pa = mesh.points[best[1]]
    pb = mesh.points[best[2]]
    bbox = _bbox_from_points((pa, pb))
    return (distance=dmin,
            h_local=hloc,
            feature_pairs=best,
            bbox=bbox,
            score=_pair_score(dmin, hloc))
end

function detect_imminent_split_3d(component::FrontComponentState,
                                  handler::LocalCartesianTopologyHandler)
    mesh = component.mesh
    mesh isa SurfaceMesh || return nothing
    n = length(mesh.points)
    n < 10 && return nothing

    center = sum(mesh.points) / n
    radii = [norm(p - center) for p in mesh.points]
    rmin = minimum(radii)
    rmean = sum(radii) / n
    hloc = local_length_scale(component)

    neck = rmean - rmin
    neck > handler.d_split * hloc || return nothing

    i = argmin(radii)
    p = mesh.points[i]
    bbox = _bbox_from_points((center, p))
    return (distance=rmin,
            h_local=hloc,
            feature_pairs=(i,),
            bbox=bbox,
            score=neck / max(hloc, eps(Float64)))
end

function find_topology_candidates(state::FrontState, handler::LocalCartesianTopologyHandler)
    return find_topology_candidates(MultiFrontState(state), handler)
end

function find_topology_candidates(state::MultiFrontState, handler::LocalCartesianTopologyHandler)
    candidates = TopologyCandidate[]
    n = ncomponents(state)

    for i in 1:n-1
        ci = component(state, i)
        for j in i+1:n
            cj = component(state, j)
            res = if ci.mesh isa CurveMesh && cj.mesh isa CurveMesh
                detect_imminent_merge_2d(ci, cj, handler)
            elseif ci.mesh isa SurfaceMesh && cj.mesh isa SurfaceMesh
                detect_imminent_merge_3d(ci, cj, handler)
            else
                nothing
            end
            if res !== nothing
                push!(candidates, TopologyCandidate(
                    :merge,
                    [i, j],
                    res.feature_pairs,
                    Float64(res.distance),
                    Float64(res.h_local),
                    Float64(res.score),
                    res.bbox,
                ))
            end
        end
    end

    for i in 1:n
        comp = component(state, i)
        split_res = comp.mesh isa CurveMesh ?
            detect_imminent_split_2d(comp, handler) :
            detect_imminent_split_3d(comp, handler)
        if split_res !== nothing
            push!(candidates, TopologyCandidate(
                :split,
                [i],
                split_res.feature_pairs,
                Float64(split_res.distance),
                Float64(split_res.h_local),
                Float64(split_res.score),
                split_res.bbox,
            ))
        end
    end

    return candidates
end

select_topology_candidate(candidates::AbstractVector{TopologyCandidate}) =
    isempty(candidates) ? nothing : candidates[argmax(getfield.(candidates, :score))]
