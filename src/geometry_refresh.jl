# geometry_refresh.jl – Rebuild geometric / DEC data after vertex motion.
#
# After vertex coordinates change (by time integration, redistribution, etc.),
# all geometry-dependent quantities must be recomputed.  This module provides
# the single public entry point `refresh_geometry!`.
#
# Design
# ------
# * For `CurveMesh`: `compute_geometry` is cheap (O(N) per vertex) and already
#   includes `signed_curvature`.  No DEC rebuild is needed unless explicitly
#   requested.
# * For `SurfaceMesh`: `compute_geometry` gives normals and areas; the DEC
#   (cotan Laplacian) must also be rebuilt whenever curvature-driven motion is
#   active, because the cotan weights change with vertex positions.
# * Connectivity (edges, face-edge incidence) is NOT rebuilt here because
#   v0.1 has fixed connectivity.

"""
    refresh_geometry!(state::FrontState; rebuild_dec=true)

Recompute the geometric and (optionally) DEC data for the current mesh.

After this call, `state.geom` contains up-to-date geometry.  If
`state.dec !== nothing` and `rebuild_dec=true`, `state.dec` is also
rebuilt (required for curvature-driven motion on surfaces).

Steps
-----
1. Call `FrontIntrinsicOps.compute_geometry(state.mesh)` → new `geom`.
2. If `rebuild_dec && state.dec !== nothing`:
   a. Call `FrontIntrinsicOps.build_dec(state.mesh, geom)` → new `dec`.
   b. If `SurfaceMesh`: call `FrontIntrinsicOps.compute_curvature(...)` to
      fill `mean_curvature_normal` into `geom`.
"""
function refresh_geometry!(state::FrontState; rebuild_dec::Bool=true)
    mesh = state.mesh
    state.geom = compute_geometry(mesh)
    if rebuild_dec && state.dec !== nothing
        new_dec = FrontIntrinsicOps.build_dec(mesh, state.geom)
        state.dec = new_dec
        if mesh isa SurfaceMesh
            state.geom = FrontIntrinsicOps.compute_curvature(mesh, state.geom, new_dec)
        end
    end
    return state
end

function _refresh_component_geometry!(component::FrontComponentState;
                                      rebuild_dec::Bool=true,
                                      build_dec::Bool=false)
    mesh = component.mesh
    component.geom = compute_geometry(mesh)

    if build_dec && component.dec === nothing
        component.dec = FrontIntrinsicOps.build_dec(mesh, component.geom)
        if mesh isa SurfaceMesh
            component.geom = FrontIntrinsicOps.compute_curvature(mesh, component.geom, component.dec)
        end
    elseif rebuild_dec && component.dec !== nothing
        component.dec = FrontIntrinsicOps.build_dec(mesh, component.geom)
        if mesh isa SurfaceMesh
            component.geom = FrontIntrinsicOps.compute_curvature(mesh, component.geom, component.dec)
        end
    end

    empty!(component.cache)
    return component
end

function refresh_geometry!(state::MultiFrontState; build_dec::Bool=false, which=:all)
    idxs = if which === :all
        eachindex(state.components)
    elseif which isa Integer
        (Int(which),)
    elseif which isa AbstractVector{<:Integer} || which isa Tuple
        Int.(collect(which))
    else
        error("refresh_geometry!(MultiFrontState): `which` must be :all, Int, or collection of Int.")
    end

    for i in idxs
        _refresh_component_geometry!(state.components[i]; rebuild_dec=true, build_dec=build_dec)
    end
    empty!(state.cache)
    return state
end
