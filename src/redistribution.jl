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
redistribution_interval(r::CurveEqualArcRedistributor) = max(r.every, 1)

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
# PoissonTangentialRedistributor
# ─────────────────────────────────────────────────────────────────────────────

"""
    PoissonTangentialRedistributor(;
        omega=1.0,
        pseudo_dt=0.2,
        iterations=1,
        every=1,
        pin_vertex=1,
        max_step_fraction=0.25,
    )

Poisson-potential tangential redistribution for closed curves and surfaces.

The redistributor solves a pinned scalar Poisson problem on the current front.
For curves it uses the finite-volume/DEC form `K * psi = M0 * b`, converts the
potential to edge tangential speeds, removes mean drift, and averages edge
speeds back to vertices.  For surfaces it applies the analogous surface
stiffness solve and moves vertices by the projected negative potential gradient.

Parameters
----------
- `omega` – source strength.
- `pseudo_dt` – pseudo-time multiplier for each redistribution iteration.
- `iterations` – number of Poisson redistribution passes per call.
- `every` – apply only every `every` accepted time steps when used in `integrate!`.
- `pin_vertex` – vertex used to remove the constant nullspace from the solve.
- `max_step_fraction` – cap each pseudo-step relative to mean edge length.
"""
struct PoissonTangentialRedistributor <: AbstractRedistributor
    omega             :: Float64
    pseudo_dt         :: Float64
    iterations        :: Int
    every             :: Int
    pin_vertex        :: Int
    max_step_fraction :: Float64
end

function PoissonTangentialRedistributor(;
    omega             :: Real = 1.0,
    pseudo_dt         :: Real = 0.2,
    iterations        :: Int  = 1,
    every             :: Int  = 1,
    pin_vertex        :: Int  = 1,
    max_step_fraction :: Real = 0.25,
)
    iterations >= 0 || throw(ArgumentError("iterations must be nonnegative."))
    every >= 1 || throw(ArgumentError("every must be >= 1."))
    pin_vertex >= 1 || throw(ArgumentError("pin_vertex must be >= 1."))
    max_step_fraction >= 0 || throw(ArgumentError("max_step_fraction must be nonnegative."))
    return PoissonTangentialRedistributor(
        Float64(omega),
        Float64(pseudo_dt),
        iterations,
        every,
        pin_vertex,
        Float64(max_step_fraction),
    )
end

redistribution_interval(r::PoissonTangentialRedistributor) = max(r.every, 1)

function redistribute!(state::FrontState, r::PoissonTangentialRedistributor)
    mesh = state.mesh
    if mesh isa CurveMesh
        is_closed(mesh) ||
            error("PoissonTangentialRedistributor requires a closed CurveMesh.")
        for _ in 1:r.iterations
            _poisson_tangential_curve_step!(state, r)
        end
    elseif mesh isa SurfaceMesh
        is_closed(mesh) ||
            error("PoissonTangentialRedistributor requires a closed SurfaceMesh.")
        for _ in 1:r.iterations
            _poisson_tangential_surface_step!(state, r)
        end
    else
        error("PoissonTangentialRedistributor requires a CurveMesh or SurfaceMesh.")
    end
    return state
end

function _poisson_rhs(weights::AbstractVector{T}, omega::Real) where {T}
    total = sum(weights)
    n = length(weights)
    total > eps(T) || return zeros(T, n)
    target = total / n
    b = Vector{T}(undef, n)
    for i in 1:n
        b[i] = weights[i] > eps(T) ? T(omega) * (target / weights[i] - one(T)) : zero(T)
    end
    mean_b = dot(weights, b) / total
    b .-= mean_b
    return b
end

function _pinned_poisson_solve(L, b::AbstractVector{T}, pin_vertex::Int) where {T}
    n = length(b)
    1 <= pin_vertex <= n || throw(ArgumentError("pin_vertex=$pin_vertex is outside 1:$n."))
    A = copy(L)
    rhs = copy(b)
    A[pin_vertex, :] .= zero(T)
    A[:, pin_vertex] .= zero(T)
    A[pin_vertex, pin_vertex] = one(T)
    rhs[pin_vertex] = zero(T)
    return A \ rhs
end

function _edge_orientation_map(mesh::CurveMesh{T}) where {T}
    edge_map = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for (ei, e) in enumerate(mesh.edges)
        i, j = e[1], e[2]
        edge_map[(i, j)] = (ei, 1)
        edge_map[(j, i)] = (ei, -1)
    end
    return edge_map
end

@inline function _cap_step(disp, max_step)
    dnorm = norm(disp)
    return dnorm > max_step && dnorm > 0 ? disp * (max_step / dnorm) : disp
end

function _poisson_tangential_curve_step!(state::FrontState, r::PoissonTangentialRedistributor)
    mesh = state.mesh
    geom = state.geom
    pts = mesh.points
    n = length(pts)
    n >= 3 || return state
    T = eltype(eltype(pts))

    dec = FrontIntrinsicOps.build_dec(mesh, geom)
    b = _poisson_rhs(geom.vertex_dual_lengths, r.omega)
    K = dec.d0' * dec.star1 * dec.d0
    rhs = geom.vertex_dual_lengths .* b
    ψ = _pinned_poisson_solve(K, rhs, min(r.pin_vertex, n))

    order = FrontIntrinsicOps.curve_vertex_order(mesh)
    nloc = length(order)
    nloc == n || error("PoissonTangentialRedistributor currently expects one closed curve component.")
    edge_map = _edge_orientation_map(mesh)

    α_by_mesh_edge = Vector{T}(undef, length(mesh.edges))
    for (ei, e) in enumerate(mesh.edges)
        i, j = e[1], e[2]
        ℓ = geom.edge_lengths[ei]
        α_by_mesh_edge[ei] = ℓ > eps(T) ? -(ψ[j] - ψ[i]) / ℓ : zero(T)
    end
    mean_α = dot(α_by_mesh_edge, geom.edge_lengths) / sum(geom.edge_lengths)
    α_by_mesh_edge .-= mean_α

    α_edge = Vector{T}(undef, nloc)
    ℓ_edge = Vector{T}(undef, nloc)
    t_edge = Vector{eltype(pts)}(undef, nloc)
    for k in 1:nloc
        i = order[k]
        j = order[mod1(k + 1, nloc)]
        entry = get(edge_map, (i, j), nothing)
        entry === nothing && error("PoissonTangentialRedistributor: missing curve edge ($i, $j).")
        ei, sign = entry
        d = pts[j] - pts[i]
        ℓ = norm(d)
        ℓ > eps(T) || error("PoissonTangentialRedistributor: zero-length curve edge ($i, $j).")
        α_edge[k] = T(sign) * α_by_mesh_edge[ei]
        ℓ_edge[k] = ℓ
        t_edge[k] = d / ℓ
    end

    new_pts = copy(pts)
    dt = T(r.pseudo_dt)

    for k in 1:nloc
        i  = order[k]

        t = t_edge[mod1(k - 1, nloc)] + t_edge[k]
        tnorm = norm(t)
        tnorm > eps(T) || continue
        tangent = t / tnorm

        α_vertex = T(0.5) * (α_edge[mod1(k - 1, nloc)] + α_edge[k])
        dual_length = T(0.5) * (ℓ_edge[mod1(k - 1, nloc)] + ℓ_edge[k])
        max_step = T(r.max_step_fraction) * dual_length
        disp = dt * α_vertex * tangent
        new_pts[i] = pts[i] + _cap_step(disp, max_step)
    end

    state.mesh = _replace_mesh(mesh, new_pts)
    refresh_geometry!(state)
    return state
end

function _poisson_tangential_surface_step!(state::FrontState, r::PoissonTangentialRedistributor)
    mesh = state.mesh
    geom = state.geom
    pts = mesh.points
    n = length(pts)
    n >= 3 || return state
    T = eltype(eltype(pts))

    dec = FrontIntrinsicOps.build_dec(mesh, geom)
    b = _poisson_rhs(geom.vertex_dual_areas, r.omega)
    K = dec.d0' * dec.star1 * dec.d0
    rhs = geom.vertex_dual_areas .* b
    ψ = _pinned_poisson_solve(K, rhs, min(r.pin_vertex, n))

    topo = build_topology(mesh)
    acc = [zero(eltype(pts)) for _ in 1:n]
    cnt = zeros(Int, n)
    for e in topo.edges
        i, j = e[1], e[2]
        d = pts[j] - pts[i]
        l2 = dot(d, d)
        l2 > eps(T) || continue
        edge_grad = ((ψ[j] - ψ[i]) / l2) * d
        acc[i] -= edge_grad
        acc[j] -= edge_grad
        cnt[i] += 1
        cnt[j] += 1
    end

    mean_len = isempty(geom.edge_lengths) ? one(T) : sum(geom.edge_lengths) / length(geom.edge_lengths)
    max_step = T(r.max_step_fraction) * mean_len
    dt = T(r.pseudo_dt)
    normals = geom.vertex_normals
    new_pts = copy(pts)
    for i in 1:n
        cnt[i] == 0 && continue
        disp = dt * (acc[i] / cnt[i])
        normal = normals[i]
        tang_disp = disp - dot(disp, normal) * normal
        new_pts[i] = pts[i] + _cap_step(tang_disp, max_step)
    end

    state.mesh = _replace_mesh(mesh, new_pts)
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
# PointFront1D validity repair (not geometric redistribution)
# ─────────────────────────────────────────────────────────────────────────────

"""
    repair_front!(state::FrontState{<:PointFront1D}; xminsep=0.0, clamp_small_gap=false)

Validate and (optionally) minimally repair a `PointFront1D` state.

- one marker: finite check only.
- two markers: finite check + strict ordering `x[1] < x[2]`.
- if `xminsep > 0` and gap is too small, either throw or clamp to `xminsep`
  when `clamp_small_gap=true`.
"""
function repair_front!(
    state::FrontState{<:PointFront1D};
    xminsep::Real=0.0,
    clamp_small_gap::Bool=false,
)
    mesh = state.mesh
    x = copy(mesh.x)
    all(isfinite, x) || error("repair_front!(PointFront1D): marker coordinates must be finite.")

    if length(x) == 2
        x[1] < x[2] ||
            error("repair_front!(PointFront1D): crossed/invalid markers (need x[1] < x[2], got $(x[1]) and $(x[2])).")
        if xminsep > 0
            gap = x[2] - x[1]
            if gap < xminsep
                if clamp_small_gap
                    x[2] = x[1] + xminsep
                else
                    error("repair_front!(PointFront1D): marker gap $gap is below xminsep=$xminsep.")
                end
            end
        end
    end

    state.mesh = PointFront1D(x, mesh.interval_is_inside)
    refresh_geometry!(state; rebuild_dec=false)
    return state
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
