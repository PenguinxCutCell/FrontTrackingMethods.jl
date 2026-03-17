# compare_fronts.jl – Front-to-front comparison metrics.
#
# Provides approximate geometric distance and shape-integral diagnostics
# that do NOT depend on a background grid.
#
# All functions operate directly on FrontState or CurveMesh / SurfaceMesh.
#
# 2-D curve metrics
# -----------------
#   nearest_distance_curve_to_curve(meshA, meshB)
#   symmetric_hausdorff_curve(meshA, meshB)
#   l2_distance_curve(meshA, meshB; nsamples)
#
# 3-D surface metrics
# -------------------
#   nearest_distance_surface_to_surface(meshA, meshB)
#   symmetric_hausdorff_surface(meshA, meshB)
#   l2_distance_surface(meshA, meshB; nsamples)
#
# Shape-integral diagnostics
# --------------------------
#   relative_area_error(mesh, mesh_ref)
#   relative_volume_error(mesh, mesh_ref)
#   centroid_error(mesh, mesh_ref)
#   perimeter_error(mesh, mesh_ref)
#   surface_area_error(mesh, mesh_ref)

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

# Sample N points uniformly along a closed curve
function _sample_curve(mesh::CurveMesh, nsamples::Int)
    pts   = mesh.points
    edges = mesh.edges
    N     = length(edges)
    T     = eltype(eltype(pts))

    # Build cumulative arc-lengths
    s = Vector{T}(undef, N + 1)
    s[1] = zero(T)
    for k in 1:N
        e = edges[k]
        s[k+1] = s[k] + norm(pts[e[2]] - pts[e[1]])
    end
    L = s[N+1]
    L < eps(T) && return pts

    samples = Vector{eltype(pts)}(undef, nsamples)
    for i in 1:nsamples
        t_i = L * (i - 1) / nsamples
        seg = searchsortedfirst(s, t_i, 1, N, Base.Order.Forward)
        seg = clamp(seg, 2, N + 1)
        k   = seg - 1
        e   = edges[k]
        ds  = s[k+1] - s[k]
        if ds < eps(T)
            samples[i] = pts[e[1]]
        else
            α = (t_i - s[k]) / ds
            samples[i] = pts[e[1]] + T(α) * (pts[e[2]] - pts[e[1]])
        end
    end
    return samples
end

# Sample points from a surface (triangle centroids + vertices)
function _sample_surface(mesh::SurfaceMesh, nsamples::Int)
    pts   = mesh.points
    faces = mesh.faces
    nv    = length(pts)
    nf    = length(faces)
    T     = eltype(eltype(pts))

    # Use vertices + face centroids
    all_pts = Vector{eltype(pts)}(undef, nv + nf)
    for (i, p) in enumerate(pts)
        all_pts[i] = p
    end
    for (i, f) in enumerate(faces)
        all_pts[nv + i] = (pts[f[1]] + pts[f[2]] + pts[f[3]]) / 3
    end

    # Subsample if too many
    if length(all_pts) <= nsamples
        return all_pts
    else
        step = max(1, div(length(all_pts), nsamples))
        return all_pts[1:step:end]
    end
end

# Distance from one point to its nearest neighbor in a set
function _min_dist_to_set(p, pts)
    best = Inf
    for q in pts
        d = norm(p - q)
        d < best && (best = d)
    end
    return best
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D curve metrics
# ─────────────────────────────────────────────────────────────────────────────

"""
    nearest_distance_curve_to_curve(meshA::CurveMesh, meshB::CurveMesh)
        -> Float64

Compute the one-sided nearest-point distance from each vertex of `meshA`
to the nearest point in `meshB` (sampled densely), and return the mean.
"""
function nearest_distance_curve_to_curve(
    meshA :: CurveMesh,
    meshB :: CurveMesh;
    nsamples :: Int = 2048,
)
    samplesB = _sample_curve(meshB, nsamples)
    dists = [_min_dist_to_set(p, samplesB) for p in meshA.points]
    return sum(dists) / length(dists)
end

"""
    symmetric_hausdorff_curve(meshA::CurveMesh, meshB::CurveMesh)
        -> Float64

Compute the approximate symmetric Hausdorff distance between two curves.
Both curves are sampled and the maximum of both one-sided max-distances
is returned.
"""
function symmetric_hausdorff_curve(
    meshA :: CurveMesh,
    meshB :: CurveMesh;
    nsamples :: Int = 2048,
)
    samplesA = _sample_curve(meshA, nsamples)
    samplesB = _sample_curve(meshB, nsamples)
    dAB = maximum(_min_dist_to_set(p, samplesB) for p in samplesA)
    dBA = maximum(_min_dist_to_set(p, samplesA) for p in samplesB)
    return max(dAB, dBA)
end

"""
    l2_distance_curve(meshA::CurveMesh, meshB::CurveMesh; nsamples=1024)
        -> Float64

Approximate L² distance between two curves: sample `nsamples` points
on `meshA`, find nearest point in `meshB`, and compute the RMS of distances.
"""
function l2_distance_curve(
    meshA    :: CurveMesh,
    meshB    :: CurveMesh;
    nsamples :: Int = 1024,
)
    samplesA = _sample_curve(meshA, nsamples)
    samplesB = _sample_curve(meshB, nsamples)
    d2 = sum((_min_dist_to_set(p, samplesB))^2 for p in samplesA) / length(samplesA)
    return sqrt(d2)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D surface metrics
# ─────────────────────────────────────────────────────────────────────────────

"""
    nearest_distance_surface_to_surface(meshA::SurfaceMesh, meshB::SurfaceMesh)
        -> Float64

Compute the mean one-sided nearest-point distance from each vertex of
`meshA` to the nearest sample point of `meshB`.
"""
function nearest_distance_surface_to_surface(
    meshA    :: SurfaceMesh,
    meshB    :: SurfaceMesh;
    nsamples :: Int = 4096,
)
    samplesB = _sample_surface(meshB, nsamples)
    dists = [_min_dist_to_set(p, samplesB) for p in meshA.points]
    return sum(dists) / length(dists)
end

"""
    symmetric_hausdorff_surface(meshA::SurfaceMesh, meshB::SurfaceMesh)
        -> Float64

Approximate symmetric Hausdorff distance between two surfaces.
"""
function symmetric_hausdorff_surface(
    meshA    :: SurfaceMesh,
    meshB    :: SurfaceMesh;
    nsamples :: Int = 4096,
)
    samplesA = _sample_surface(meshA, nsamples)
    samplesB = _sample_surface(meshB, nsamples)
    dAB = maximum(_min_dist_to_set(p, samplesB) for p in samplesA)
    dBA = maximum(_min_dist_to_set(p, samplesA) for p in samplesB)
    return max(dAB, dBA)
end

"""
    l2_distance_surface(meshA::SurfaceMesh, meshB::SurfaceMesh; nsamples=2048)
        -> Float64

Approximate L² distance between two surfaces.
"""
function l2_distance_surface(
    meshA    :: SurfaceMesh,
    meshB    :: SurfaceMesh;
    nsamples :: Int = 2048,
)
    samplesA = _sample_surface(meshA, nsamples)
    samplesB = _sample_surface(meshB, nsamples)
    d2 = sum((_min_dist_to_set(p, samplesB))^2 for p in samplesA) / length(samplesA)
    return sqrt(d2)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shape-integral diagnostics
# ─────────────────────────────────────────────────────────────────────────────

"""
    relative_area_error(mesh::CurveMesh, mesh_ref::CurveMesh) -> Float64

Relative difference in enclosed area between `mesh` and `mesh_ref`.
Returns `(area - area_ref) / area_ref`.
"""
function relative_area_error(mesh::CurveMesh, mesh_ref::CurveMesh)
    A     = enclosed_measure(mesh)
    A_ref = enclosed_measure(mesh_ref)
    abs(A_ref) < eps(Float64) && return abs(A)
    return (A - A_ref) / A_ref
end

"""
    relative_volume_error(mesh::SurfaceMesh, mesh_ref::SurfaceMesh) -> Float64

Relative difference in enclosed volume between `mesh` and `mesh_ref`.
Returns `(volume - volume_ref) / volume_ref`.
"""
function relative_volume_error(mesh::SurfaceMesh, mesh_ref::SurfaceMesh)
    V     = enclosed_measure(mesh)
    V_ref = enclosed_measure(mesh_ref)
    abs(V_ref) < eps(Float64) && return abs(V)
    return (V - V_ref) / V_ref
end

"""
    centroid_error(meshA, meshB) -> Float64

Euclidean distance between the vertex centroids of `meshA` and `meshB`.
"""
function centroid_error(meshA, meshB)
    cA = sum(meshA.points) / length(meshA.points)
    cB = sum(meshB.points) / length(meshB.points)
    return norm(cA - cB)
end

"""
    perimeter_error(mesh::CurveMesh, mesh_ref::CurveMesh) -> Float64

Relative difference in arc-length (perimeter) between `mesh` and `mesh_ref`.
Returns `(L - L_ref) / L_ref`.
"""
function perimeter_error(mesh::CurveMesh, mesh_ref::CurveMesh)
    geom     = compute_geometry(mesh)
    geom_ref = compute_geometry(mesh_ref)
    L     = sum(geom.edge_lengths)
    L_ref = sum(geom_ref.edge_lengths)
    abs(L_ref) < eps(Float64) && return abs(L)
    return (L - L_ref) / L_ref
end

"""
    surface_area_error(mesh::SurfaceMesh, mesh_ref::SurfaceMesh) -> Float64

Relative difference in surface area between `mesh` and `mesh_ref`.
Returns `(A - A_ref) / A_ref`.
"""
function surface_area_error(mesh::SurfaceMesh, mesh_ref::SurfaceMesh)
    geom     = compute_geometry(mesh)
    geom_ref = compute_geometry(mesh_ref)
    A     = measure(mesh, geom)
    A_ref = measure(mesh_ref, geom_ref)
    abs(A_ref) < eps(Float64) && return abs(A)
    return (A - A_ref) / A_ref
end
