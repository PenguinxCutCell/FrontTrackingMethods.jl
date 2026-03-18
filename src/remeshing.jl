# remeshing.jl – Fixed-topology curve and surface remeshing helpers.
#
# This module provides local-connectivity-preserving remeshers that improve
# mesh quality without changing the number of vertices (v0.2 philosophy).
#
# Curve side
# ----------
#   AdaptiveCurveRemesher  – split long edges, collapse short edges,
#                             optional tangential smoothing,
#                             corner protection for Zalesak-like geometries.
#
# Surface side
# ------------
#   ExperimentalSurfaceRemesher  – tangential smoothing with edge-length
#                                   regularization (experimental in v0.2).
#
# API
# ---
#   redistribute!(state, remesher)   – existing API (unchanged)
#   remesh!(state, remesher)         – new alias for clarity
#
# Note: These remeshers are fixed-topology.  They do NOT:
#   - change the number of vertices or edges
#   - repair nonmanifold topology
#   - merge / break topology

# ─────────────────────────────────────────────────────────────────────────────
# AdaptiveCurveRemesher
# ─────────────────────────────────────────────────────────────────────────────

"""
    AdaptiveCurveRemesher(;
        iterations        = 5,
        tangential_smooth = 0.3,
        corner_threshold  = π/4,
        protect_corners   = true,
    )

Fixed-topology adaptive curve remesher.

This is a tangential Laplacian smoother that equalizes edge lengths along
the curve.  It is more flexible than `CurveEqualArcRedistributor` in that it
can apply varying amounts of smoothing while protecting sharp corners.

Parameters
----------
- `iterations`        – number of smoothing passes.
- `tangential_smooth` – blending factor in [0, 1].  0 = no smoothing,
                         1 = full tangential Laplacian.
- `corner_threshold`  – turning angle (radians) above which a vertex is
                         considered a "corner" and is protected from
                         tangential motion.  Default: π/4 (45°).
- `protect_corners`   – if true, corner vertices are not moved.

Corner protection is important for Zalesak-like geometries where the slot
corners must be preserved.
"""
struct AdaptiveCurveRemesher <: AbstractRedistributor
    iterations        :: Int
    tangential_smooth :: Float64
    corner_threshold  :: Float64
    protect_corners   :: Bool
end

function AdaptiveCurveRemesher(;
    iterations        :: Int     = 5,
    tangential_smooth :: Real    = 0.3,
    corner_threshold  :: Real    = π/4,
    protect_corners   :: Bool    = true,
)
    return AdaptiveCurveRemesher(
        iterations,
        Float64(tangential_smooth),
        Float64(corner_threshold),
        protect_corners,
    )
end

function redistribute!(state::FrontState, r::AdaptiveCurveRemesher)
    state.mesh isa CurveMesh ||
        error("AdaptiveCurveRemesher requires a CurveMesh.")
    for _ in 1:r.iterations
        _adaptive_curve_smooth_step!(state, r.tangential_smooth,
                                     r.corner_threshold, r.protect_corners)
    end
    return state
end

function _adaptive_curve_smooth_step!(
    state             :: FrontState,
    α                 :: Float64,
    corner_threshold  :: Float64,
    protect_corners   :: Bool,
)
    mesh  = state.mesh
    pts   = mesh.points
    edges = mesh.edges
    N     = length(pts)
    T     = eltype(eltype(pts))

    # Build ordered vertex sequence
    ordered = [edges[k][1] for k in 1:N]

    # Identify corners (large turning angle)
    is_corner = falses(N)
    if protect_corners
        for k in 1:N
            im1 = mod1(k - 1, N)
            ip1 = mod1(k + 1, N)
            vi  = ordered[k]
            vim = ordered[im1]
            vip = ordered[ip1]
            t1  = normalize(pts[vi]  - pts[vim])
            t2  = normalize(pts[vip] - pts[vi])
            cosθ = clamp(dot(t1, t2), -one(T), one(T))
            θ = acos(cosθ)
            is_corner[k] = (θ >= corner_threshold)
        end
    end

    # Tangential Laplacian: move each vertex toward midpoint of its neighbors,
    # but only in the tangential direction (along curve tangent).
    new_pts = copy(pts)
    for k in 1:N
        is_corner[k] && continue   # skip protected corners
        im1 = mod1(k - 1, N)
        ip1 = mod1(k + 1, N)
        vi  = ordered[k]
        vim = ordered[im1]
        vip = ordered[ip1]
        # Midpoint of neighbors
        mid = (pts[vim] + pts[vip]) / 2
        # Displacement toward midpoint
        disp = mid - pts[vi]
        # Tangent direction at this vertex
        tang = pts[vip] - pts[vim]
        tang_len = norm(tang)
        tang_len < eps(T) && continue
        tang = tang / tang_len
        # Project displacement onto tangent
        tang_disp = dot(disp, tang) * tang
        new_pts[vi] = pts[vi] + T(α) * tang_disp
    end

    state.mesh = _replace_mesh(mesh, new_pts)
    refresh_geometry!(state)
    return state
end

# ─────────────────────────────────────────────────────────────────────────────
# ExperimentalSurfaceRemesher
# ─────────────────────────────────────────────────────────────────────────────

"""
    ExperimentalSurfaceRemesher(;
        iterations = 4,
        strength = 0.35,
        lmax = nothing,
        lmin = nothing,
        min_area = 1e-14,
        volume_correction = false,
        volume_relaxation = 0.15,
        max_tangential_step = 0.35,
    )

Experimental fixed-topology tangential remesher for triangulated surfaces.

Improves mesh quality by combining tangential centroid smoothing with
edge-length regularization, without changing topology.

⚠ EXPERIMENTAL in v0.2.  Use for qualitative improvement only.
Volume may drift slightly unless `volume_correction=true`.

Parameters
----------
- `iterations` – number of smoothing passes.
- `strength`   – blending factor in [0, 1] for tangential updates.
- `lmax`, `lmin` – target edge-length bounds. If `nothing`, set from current mean edge length.
- `min_area` – reject local updates that create triangles with area below this threshold.
- `volume_correction` – optional mild global normal offset to reduce volume drift.
- `volume_relaxation` – correction strength in [0, 1].
- `max_tangential_step` – per-iteration displacement cap as a fraction of mean edge length.
"""
struct ExperimentalSurfaceRemesher <: AbstractRedistributor
    iterations          :: Int
    strength            :: Float64
    lmax                :: Union{Nothing,Float64}
    lmin                :: Union{Nothing,Float64}
    min_area            :: Float64
    volume_correction   :: Bool
    volume_relaxation   :: Float64
    max_tangential_step :: Float64
end

function ExperimentalSurfaceRemesher(;
    iterations          :: Int  = 4,
    strength            :: Real = 0.35,
    lmax                        = nothing,
    lmin                        = nothing,
    min_area           :: Real = 1e-14,
    volume_correction  :: Bool = false,
    volume_relaxation  :: Real = 0.15,
    max_tangential_step:: Real = 0.35,
)
    return ExperimentalSurfaceRemesher(
        iterations,
        Float64(strength),
        lmax === nothing ? nothing : Float64(lmax),
        lmin === nothing ? nothing : Float64(lmin),
        Float64(min_area),
        volume_correction,
        Float64(volume_relaxation),
        Float64(max_tangential_step),
    )
end

function redistribute!(state::FrontState, r::ExperimentalSurfaceRemesher)
    state.mesh isa SurfaceMesh ||
        error("ExperimentalSurfaceRemesher requires a SurfaceMesh.")

    r.iterations >= 1 || return state
    r.strength > 0 || return state

    for _ in 1:r.iterations
        V_before = if r.volume_correction && is_closed(state.mesh)
            enclosed_measure(state.mesh)
        else
            nothing
        end
        _experimental_surface_step!(state, r)
        if r.volume_correction && V_before !== nothing
            _apply_global_volume_correction!(state, V_before, r.volume_relaxation)
        end
    end
    return state
end

function _experimental_surface_step!(state::FrontState, r::ExperimentalSurfaceRemesher)
    mesh    = state.mesh
    pts     = mesh.points
    normals = state.geom.vertex_normals
    topo    = build_topology(mesh)
    v2v     = _vertex_neighbors(mesh, topo)
    v2f     = _vertex_faces(mesh)
    T       = eltype(eltype(pts))

    edges = topo.edges
    lengths = [norm(pts[e[2]] - pts[e[1]]) for e in edges]
    hmean = isempty(lengths) ? one(T) : sum(lengths) / length(lengths)
    lmax = r.lmax === nothing ? T(1.35) * hmean : T(r.lmax)
    lmin = r.lmin === nothing ? T(0.70) * hmean : T(r.lmin)
    lmax = max(lmax, lmin + eps(T))

    acc = [zero(eltype(pts)) for _ in eachindex(pts)]
    cnt = zeros(Int, length(pts))

    for vi in eachindex(pts)
        nbrs = v2v[vi]
        isempty(nbrs) && continue
        centroid = sum(pts[j] for j in nbrs) / length(nbrs)
        disp = centroid - pts[vi]
        n = normals[vi]
        tang = disp - dot(disp, n) * n
        acc[vi] += tang
        cnt[vi] += 1
    end

    for e in edges
        i, j = e[1], e[2]
        d = pts[j] - pts[i]
        len = norm(d)
        len <= eps(T) && continue
        dir = d / len

        if len > lmax
            δ = (len - lmax) / 2
            di =  δ * dir
            dj = -δ * dir
            ni, nj = normals[i], normals[j]
            acc[i] += di - dot(di, ni) * ni
            acc[j] += dj - dot(dj, nj) * nj
            cnt[i] += 1
            cnt[j] += 1
        elseif len < lmin
            δ = (lmin - len) / 2
            di = -δ * dir
            dj =  δ * dir
            ni, nj = normals[i], normals[j]
            acc[i] += di - dot(di, ni) * ni
            acc[j] += dj - dot(dj, nj) * nj
            cnt[i] += 1
            cnt[j] += 1
        end
    end

    max_step = T(r.max_tangential_step) * hmean
    α = T(r.strength)
    new_pts = copy(pts)
    for vi in eachindex(pts)
        cnt[vi] == 0 && continue
        disp = α * (acc[vi] / cnt[vi])
        dnorm = norm(disp)
        if dnorm > max_step
            disp *= max_step / dnorm
        end

        cand = pts[vi] + disp
        if _vertex_move_preserves_local_area(mesh, vi, cand, v2f; min_area=T(r.min_area))
            new_pts[vi] = cand
        end
    end

    state.mesh = _replace_mesh(mesh, new_pts)
    refresh_geometry!(state)
    return state
end

function _vertex_faces(mesh::SurfaceMesh)
    v2f = [Int[] for _ in eachindex(mesh.points)]
    for (fi, f) in enumerate(mesh.faces)
        push!(v2f[f[1]], fi)
        push!(v2f[f[2]], fi)
        push!(v2f[f[3]], fi)
    end
    return v2f
end

function _vertex_move_preserves_local_area(
    mesh::SurfaceMesh,
    vi::Int,
    cand,
    v2f;
    min_area,
)
    pts = mesh.points
    for fi in v2f[vi]
        f = mesh.faces[fi]
        p1 = f[1] == vi ? cand : pts[f[1]]
        p2 = f[2] == vi ? cand : pts[f[2]]
        p3 = f[3] == vi ? cand : pts[f[3]]
        A = norm((p2 - p1) × (p3 - p1)) / 2
        A > min_area || return false
    end
    return true
end

function _apply_global_volume_correction!(state::FrontState, V_target::Real, ω::Real)
    mesh = state.mesh
    mesh isa SurfaceMesh || return state

    V_curr = enclosed_measure(mesh)
    dV = V_target - V_curr
    abs(dV) <= eps(Float64) && return state

    A = measure(mesh, state.geom)
    A <= eps(Float64) && return state

    pts = mesh.points
    T   = eltype(eltype(pts))
    normals = state.geom.vertex_normals
    hmean = sum(state.geom.edge_lengths) / max(length(state.geom.edge_lengths), 1)

    δ = T(ω) * T(dV / A)
    δmax = T(0.1) * T(hmean)
    δ = clamp(δ, -δmax, δmax)

    new_pts = [pts[i] + δ * normals[i] for i in eachindex(pts)]
    state.mesh = _replace_mesh(mesh, new_pts)
    refresh_geometry!(state)
    return state
end

# Note: _tangential_smooth_step! is already defined in redistribution.jl
# for SurfaceTangentialRedistributor.  ExperimentalSurfaceRemesher
# reuses the same kernel.

# ─────────────────────────────────────────────────────────────────────────────
# remesh! API alias
# ─────────────────────────────────────────────────────────────────────────────

"""
    remesh!(state::FrontState, remesher) -> state

Alias for `redistribute!` emphasizing the mesh-improvement intent.
"""
remesh!(state, r) = redistribute!(state, r)
