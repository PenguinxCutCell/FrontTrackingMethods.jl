# transfer.jl – Field transfer between old and new fronts.
#
# After redistribution moves vertices, any vertex-attached fields may need
# to be transferred from the old vertex positions to the new ones.
# This is the front-tracking analogue of velocity extension.
#
# v0.1 methods
# ------------
#   :piecewise_linear  (default for curves)  – linear interpolation on arc-length
#   :nearest_vertex    (fallback)            – nearest-vertex projection
#
# v0.2 methods
# ------------
#   :barycentric  (surfaces)  – barycentric interpolation on nearest triangle
#
# Constant fields are preserved exactly by :piecewise_linear and :barycentric.

"""
    transfer_vertex_field!(newvals, oldmesh, oldvals, newmesh;
                           method=:piecewise_linear) -> newvals

Transfer a vertex scalar field `oldvals` defined on `oldmesh` to vertex
positions of `newmesh`, writing results into `newvals`.

`method` for curves:
- `:piecewise_linear` – interpolate on piecewise-linear arc-length (default).
- `:nearest_vertex`   – copy the value of the nearest old vertex.

`method` for surfaces:
- `:barycentric`     – barycentric interpolation on the nearest triangle (default).
- `:nearest_vertex`  – copy the value of the nearest old vertex (fallback).

For constant fields, `:piecewise_linear` and `:barycentric` preserve the
constant exactly.
"""
function transfer_vertex_field!(
        newvals   :: AbstractVector,
        oldmesh   :: CurveMesh{T},
        oldvals   :: AbstractVector,
        newmesh   :: CurveMesh{T};
        method    :: Symbol = :piecewise_linear,
) :: AbstractVector where {T}

    if method === :nearest_vertex
        return _transfer_nearest!(newvals, oldmesh, oldvals, newmesh)
    elseif method === :piecewise_linear
        return _transfer_arc_linear!(newvals, oldmesh, oldvals, newmesh)
    else
        error("transfer_vertex_field!: unknown method $(repr(method)). " *
              "Use :piecewise_linear or :nearest_vertex.")
    end
end

function transfer_vertex_field!(
        newvals   :: AbstractVector,
        oldmesh   :: SurfaceMesh{T},
        oldvals   :: AbstractVector,
        newmesh   :: SurfaceMesh{T};
        method    :: Symbol = :barycentric,
) :: AbstractVector where {T}

    if method === :nearest_vertex
        return _transfer_nearest!(newvals, oldmesh, oldvals, newmesh)
    elseif method === :barycentric
        return _transfer_barycentric!(newvals, oldmesh, oldvals, newmesh)
    else
        error("transfer_vertex_field! on SurfaceMesh: supported methods are " *
              ":barycentric (default) and :nearest_vertex.")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Nearest-vertex transfer (works for any mesh type)
# ─────────────────────────────────────────────────────────────────────────────

function _transfer_nearest!(newvals, oldmesh, oldvals, newmesh)
    old_pts = oldmesh.points
    new_pts = newmesh.points
    for (bi, p) in enumerate(new_pts)
        best_ai   = 1
        best_dist = norm(p - old_pts[1])
        for ai in 2:length(old_pts)
            d = norm(p - old_pts[ai])
            if d < best_dist
                best_dist = d
                best_ai   = ai
            end
        end
        newvals[bi] = oldvals[best_ai]
    end
    return newvals
end

# ─────────────────────────────────────────────────────────────────────────────
# Piecewise-linear arc-length interpolation (curves only)
# ─────────────────────────────────────────────────────────────────────────────

function _transfer_arc_linear!(newvals, oldmesh::CurveMesh{T}, oldvals, newmesh::CurveMesh{T}) where {T}
    old_pts  = oldmesh.points
    new_pts  = newmesh.points
    old_edges = oldmesh.edges
    N_old    = length(old_pts)

    # Build cumulative arc-length along old curve
    s_old = Vector{T}(undef, N_old + 1)
    s_old[1] = zero(T)
    # ordered vertex sequence from edges
    ordered = [old_edges[k][1] for k in 1:N_old]
    for k in 1:N_old
        vi = ordered[k]
        vj = ordered[mod1(k+1, N_old)]
        s_old[k+1] = s_old[k] + norm(old_pts[vj] - old_pts[vi])
    end
    L_old = s_old[N_old+1]
    L_old < eps(T) && (fill!(newvals, oldvals[1]); return newvals)

    # New arc-length (only for closed curves: approximate by matching total length)
    new_edges = newmesh.edges
    N_new     = length(new_pts)
    ordered_new = [new_edges[k][1] for k in 1:N_new]
    s_new = Vector{T}(undef, N_new + 1)
    s_new[1] = zero(T)
    for k in 1:N_new
        vi = ordered_new[k]
        vj = ordered_new[mod1(k+1, N_new)]
        s_new[k+1] = s_new[k] + norm(new_pts[vj] - new_pts[vi])
    end
    L_new = s_new[N_new+1]

    for bi in 1:N_new
        # Map new arc position to old arc parameter space
        t_new = s_new[bi] / max(L_new, eps(T))
        t_old = t_new * L_old   # target position on old curve

        # Find segment in old arc-length
        seg = searchsortedfirst(s_old, t_old, 1, N_old, Base.Order.Forward)
        seg = clamp(seg, 2, N_old + 1)
        i   = seg - 1
        ds  = s_old[i+1] - s_old[i]
        if ds < eps(T)
            newvals[ordered_new[bi]] = oldvals[ordered[i]]
        else
            alpha = (t_old - s_old[i]) / ds
            vi    = ordered[i]
            vj    = ordered[mod1(i+1, N_old)]
            val_i = oldvals[vi]
            val_j = oldvals[vj]
            newvals[ordered_new[bi]] = val_i + T(alpha) * (val_j - val_i)
        end
    end
    return newvals
end

# ─────────────────────────────────────────────────────────────────────────────
# transfer_fields! – transfer all FrontFields in a state after redistribution
# ─────────────────────────────────────────────────────────────────────────────

"""
    transfer_fields!(state::FrontState, oldmesh; method=:piecewise_linear)

Transfer all vertex fields in `state.fields` from `oldmesh` to the current
`state.mesh`.  Non-vertex fields are left unchanged.

`oldmesh` must have the same connectivity as `state.mesh`.
"""
function transfer_fields!(state::FrontState, oldmesh; method::Symbol=:piecewise_linear)
    for (name, field) in state.fields
        if field isa FrontField && field.location === :vertex
            newvals = similar(field.values)
            transfer_vertex_field!(newvals, oldmesh, field.values, state.mesh;
                                   method=method)
            state.fields[name] = FrontField(newvals, state.mesh, :vertex)
        end
    end
    return state
end

# ─────────────────────────────────────────────────────────────────────────────
# Barycentric transfer for surfaces (v0.2)
# ─────────────────────────────────────────────────────────────────────────────

"""
    _transfer_barycentric!(newvals, oldmesh, oldvals, newmesh)

Transfer a vertex field from `oldmesh` to `newmesh` using barycentric
interpolation on the closest triangle.

For each new vertex, we find the closest triangle in the old mesh and
interpolate the field using barycentric coordinates of the projected point.
Falls back to nearest-vertex if projection fails.

Constant fields are preserved exactly.
"""
function _transfer_barycentric!(
    newvals  :: AbstractVector,
    oldmesh  :: SurfaceMesh{T},
    oldvals  :: AbstractVector,
    newmesh  :: SurfaceMesh{T},
) where {T}
    old_pts  = oldmesh.points
    old_faces = oldmesh.faces
    new_pts  = newmesh.points

    for (bi, p) in enumerate(new_pts)
        best_dist  = Inf
        best_val   = oldvals[1]

        for f in old_faces
            p1 = old_pts[f[1]]
            p2 = old_pts[f[2]]
            p3 = old_pts[f[3]]

            # Project p onto the plane of triangle (p1, p2, p3)
            e1   = p2 - p1
            e2   = p3 - p1
            n    = e1 × e2
            n_len = norm(n)
            n_len < eps(T) && continue   # degenerate triangle

            # Closest point on triangle to p
            proj, λ1, λ2, λ3 = _closest_point_triangle(p, p1, p2, p3, e1, e2, n, n_len)

            d = norm(p - proj)
            if d < best_dist
                best_dist = d
                best_val  = λ1 * oldvals[f[1]] + λ2 * oldvals[f[2]] + λ3 * oldvals[f[3]]
            end
        end

        newvals[bi] = best_val
    end
    return newvals
end

# Project point p onto triangle (p1, p2, p3) and return the closest point
# plus its barycentric coordinates (λ1, λ2, λ3) clamped to [0,1].
@inline function _closest_point_triangle(
    p    :: SVector{3,T},
    p1   :: SVector{3,T},
    p2   :: SVector{3,T},
    p3   :: SVector{3,T},
    e1   :: SVector{3,T},    # p2 - p1
    e2   :: SVector{3,T},    # p3 - p1
    n    :: SVector{3,T},    # e1 × e2
    n_len :: T,
) where {T}
    # Project p onto the triangle plane
    d    = dot(n, p - p1) / (n_len * n_len)
    proj = p - d * n

    # Barycentric coordinates of proj w.r.t. triangle
    r = proj - p1
    d11 = dot(e1, e1)
    d12 = dot(e1, e2)
    d22 = dot(e2, e2)
    dr1 = dot(r,  e1)
    dr2 = dot(r,  e2)
    denom = d11 * d22 - d12 * d12
    if abs(denom) < eps(T)
        return p1, one(T), zero(T), zero(T)
    end
    λ2 = (d22 * dr1 - d12 * dr2) / denom
    λ3 = (d11 * dr2 - d12 * dr1) / denom
    λ1 = one(T) - λ2 - λ3

    # Clamp to triangle
    λ1 = max(zero(T), λ1)
    λ2 = max(zero(T), λ2)
    λ3 = max(zero(T), λ3)
    s  = λ1 + λ2 + λ3
    if s > eps(T)
        λ1 /= s; λ2 /= s; λ3 /= s
    else
        λ1 = one(T); λ2 = zero(T); λ3 = zero(T)
    end

    closest = λ1 * p1 + λ2 * p2 + λ3 * p3
    return closest, λ1, λ2, λ3
end
