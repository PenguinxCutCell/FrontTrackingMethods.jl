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
# Constant fields are preserved exactly by :piecewise_linear interpolation.

"""
    transfer_vertex_field!(newvals, oldmesh, oldvals, newmesh;
                           method=:piecewise_linear) -> newvals

Transfer a vertex scalar field `oldvals` defined on `oldmesh` to vertex
positions of `newmesh`, writing results into `newvals`.

`method`:
- `:piecewise_linear` – interpolate on piecewise-linear parameter for curves.
- `:nearest_vertex`   – copy the value of the nearest old vertex.

For constant fields, `:piecewise_linear` preserves the constant exactly.
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
        method    :: Symbol = :nearest_vertex,
) :: AbstractVector where {T}

    if method === :nearest_vertex
        return _transfer_nearest!(newvals, oldmesh, oldvals, newmesh)
    else
        error("transfer_vertex_field! on SurfaceMesh: only :nearest_vertex is " *
              "supported in v0.1.")
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
