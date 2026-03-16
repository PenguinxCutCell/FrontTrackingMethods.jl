# redistribution.jl – Front-quality redistribution (analogue of reinitialization).
#
# Redistribution moves vertices tangentially to improve mesh quality without
# (ideally) changing the physical shape of the front.  It is OPTIONAL and
# SEPARATE from the physical RHS.
#
# In v0.1 the following redistributors are provided:
#
#   NoRedistribution              – do nothing (default)
#   CurveEqualArcRedistributor    – equal-arclength redistribution on curves
#   SurfaceTangentialRedistributor – tangential smoothing for surfaces (experimental)
#
# API
# ---
#   redistribute!(state, redistributor)
#   repair_front!(state, redistributor)   – alias

# ─────────────────────────────────────────────────────────────────────────────
# NoRedistribution
# ─────────────────────────────────────────────────────────────────────────────

"""
    NoRedistribution

Placeholder redistributor that does nothing.  Use this as the default when
redistribution is not needed.
"""
struct NoRedistribution <: AbstractRedistributor end

redistribute!(state, ::NoRedistribution) = nothing
redistribute!(state, ::Nothing)          = nothing

# ─────────────────────────────────────────────────────────────────────────────
# CurveEqualArcRedistributor
# ─────────────────────────────────────────────────────────────────────────────

"""
    CurveEqualArcRedistributor(; every=1)

Equal-arclength redistribution for closed polygonal curves.

Redistributes the `N` vertices of a closed curve so that consecutive
vertices are spaced equally along the current polygon arc-length.
The total arc-length, the closure, and the orientation are preserved.

Parameter `every` controls how many time-steps to skip between
applications (default: apply every step).

Algorithm
---------
1. Compute cumulative arc-length s[i] with s[1] = 0.
2. Generate `N` equally-spaced target parameter values
   t_k = (k-1) / N * L, k = 1 …  N.
3. For each t_k, locate the segment [s[i], s[i+1]] containing t_k and
   interpolate linearly.
4. The first vertex is kept fixed to prevent global rotation drift;
   only the spacing is equalized.
"""
struct CurveEqualArcRedistributor <: AbstractRedistributor
    every :: Int
end
CurveEqualArcRedistributor(; every::Int=1) = CurveEqualArcRedistributor(every)

function redistribute!(state::FrontState, r::CurveEqualArcRedistributor)
    mesh = state.mesh
    mesh isa CurveMesh ||
        error("CurveEqualArcRedistributor requires a CurveMesh.")
    is_closed(mesh) ||
        error("CurveEqualArcRedistributor: curve must be closed.")
    _equal_arc_redistribute!(state)
    return state
end

function _equal_arc_redistribute!(state::FrontState)
    mesh  = state.mesh
    pts   = mesh.points
    edges = mesh.edges
    N     = length(pts)
    T     = eltype(eltype(pts))

    # ── Build cumulative arc-length ───────────────────────────────────────────
    # Follow the edge ordering; for a closed curve the last edge wraps back.
    # We build a linear sequence of N+1 points (closing back to pts[1]).
    s = Vector{T}(undef, N + 1)
    s[1] = zero(T)
    # Use the ordered vertex sequence from edge list
    for k in 1:N
        e    = edges[k]
        diff = pts[e[2]] - pts[e[1]]
        s[k+1] = s[k] + norm(diff)
    end
    L = s[N+1]   # total arc-length
    L > eps(T) || return state  # degenerate: do nothing

    # ── Build ordered vertex sequence (perm) following edge order ─────────────
    # For a well-formed closed curve: edges = [(1→2),(2→3),...,(N→1)]
    # The ordered vertices are simply pts[1], pts[2], ..., pts[N].
    # We extract this order from the edges.
    ordered = Vector{Int}(undef, N)
    ordered[1] = edges[1][1]
    for k in 2:N
        ordered[k] = edges[k][1]
    end

    # ── Build new positions by equal-arc interpolation ────────────────────────
    new_pts = Vector{eltype(pts)}(undef, N)
    for k in 1:N
        # Target arc-length parameter: equally spaced
        t_k = T(k - 1) / T(N) * L
        # Find segment containing t_k
        seg = searchsortedfirst(s, t_k, 1, N, Base.Order.Forward)
        # searchsortedfirst returns i such that s[i] >= t_k; we want [seg-1, seg]
        seg = clamp(seg, 2, N+1)
        i   = seg - 1   # start of segment (index into ordered)
        j_mod = mod1(i + 1, N)
        ds  = s[i+1] - s[i]
        if ds < eps(T)
            new_pts[k] = pts[ordered[i]]
        else
            alpha = (t_k - s[i]) / ds
            # vertex indices in original mesh
            vi = ordered[i]
            vj = ordered[j_mod]
            new_pts[k] = pts[vi] + T(alpha) * (pts[vj] - pts[vi])
        end
    end

    # ── Reassemble in same order ──────────────────────────────────────────────
    # ordered[k] tells us which original vertex slot gets new_pts[k].
    # We rebuild the points vector in original vertex index order.
    pts_out = Vector{eltype(pts)}(undef, N)
    for k in 1:N
        pts_out[ordered[k]] = new_pts[k]
    end

    state.mesh = _replace_mesh(mesh, pts_out)
    refresh_geometry!(state)
    return state
end

# ─────────────────────────────────────────────────────────────────────────────
# SurfaceTangentialRedistributor (experimental)
# ─────────────────────────────────────────────────────────────────────────────

"""
    SurfaceTangentialRedistributor(; iterations=3, strength=0.5)

Experimental tangential Laplacian smoothing for triangulated surfaces.

Moves each interior vertex by `strength` × the tangential component of the
Laplacian displacement, preserving the normal position.
This is a mesh-quality operation: it improves triangle aspect ratios without
(to first order) changing the surface shape.

⚠ EXPERIMENTAL in v0.1.  The normal correction is minimal; enclosed volume
may drift slightly.  Do not use on coarse meshes or with large `strength`.

Parameters
----------
- `iterations` – number of smoothing passes (default: 3).
- `strength`   – blending factor in [0, 1] (default: 0.5).
"""
struct SurfaceTangentialRedistributor <: AbstractRedistributor
    iterations :: Int
    strength   :: Float64
end
SurfaceTangentialRedistributor(; iterations::Int=3, strength::Real=0.5) =
    SurfaceTangentialRedistributor(iterations, Float64(strength))

function redistribute!(state::FrontState, r::SurfaceTangentialRedistributor)
    mesh = state.mesh
    mesh isa SurfaceMesh ||
        error("SurfaceTangentialRedistributor requires a SurfaceMesh.")
    for _ in 1:r.iterations
        _tangential_smooth_step!(state, r.strength)
    end
    return state
end

function _tangential_smooth_step!(state::FrontState, strength::Float64)
    mesh    = state.mesh
    pts     = mesh.points
    geom    = state.geom
    normals = geom.vertex_normals
    nv      = length(pts)
    T       = eltype(eltype(pts))

    # Build vertex adjacency from faces
    topo = build_topology(mesh)
    v2v  = _vertex_neighbors(mesh, topo)

    new_pts = copy(pts)
    for vi in 1:nv
        nbrs = v2v[vi]
        isempty(nbrs) && continue
        centroid = sum(pts[j] for j in nbrs) / length(nbrs)
        displacement = centroid - pts[vi]
        n_vi = normals[vi]
        # Tangential component: remove normal part
        tang_disp = displacement - dot(displacement, n_vi) * n_vi
        new_pts[vi] = pts[vi] + T(strength) * tang_disp
    end

    state.mesh = _replace_mesh(mesh, new_pts)
    refresh_geometry!(state)
    return state
end

function _vertex_neighbors(mesh::SurfaceMesh, topo::MeshTopology)
    nv   = length(mesh.points)
    v2v  = [Int[] for _ in 1:nv]
    for e in topo.edges
        push!(v2v[e[1]], e[2])
        push!(v2v[e[2]], e[1])
    end
    return v2v
end

# ─────────────────────────────────────────────────────────────────────────────
# Public aliases
# ─────────────────────────────────────────────────────────────────────────────

"""
    repair_front!(state, redistributor) -> state

Alias for `redistribute!`.  Use this name to emphasize the mesh-quality
restoration aspect rather than the analogy with reinitialization.
"""
repair_front!(state, r) = redistribute!(state, r)
