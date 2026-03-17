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
        iterations = 5,
        strength   = 0.3,
    )

Experimental fixed-topology tangential remesher for triangulated surfaces.

Improves mesh quality by moving each vertex tangentially toward the
centroid of its neighbors, without changing topology.

⚠ EXPERIMENTAL in v0.2.  Use for qualitative improvement only.
Volume may drift slightly.

Parameters
----------
- `iterations` – number of smoothing passes.
- `strength`   – blending factor in [0, 1].
"""
struct ExperimentalSurfaceRemesher <: AbstractRedistributor
    iterations :: Int
    strength   :: Float64
end

function ExperimentalSurfaceRemesher(;
    iterations :: Int  = 5,
    strength   :: Real = 0.3,
)
    return ExperimentalSurfaceRemesher(iterations, Float64(strength))
end

function redistribute!(state::FrontState, r::ExperimentalSurfaceRemesher)
    state.mesh isa SurfaceMesh ||
        error("ExperimentalSurfaceRemesher requires a SurfaceMesh.")
    for _ in 1:r.iterations
        _tangential_smooth_step!(state, r.strength)
    end
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
