"""
    FrontTrackingMethods

A package for evolving explicit front-tracking meshes in time.

Mirrors the user workflow of `LevelSetMethods.jl`, but operates on Lagrangian
polygonal curves and triangulated surfaces built on top of `FrontIntrinsicOps.jl`.

v0.1 supports:
- Closed 2-D polygonal curves (`CurveMesh`) and closed 3-D surfaces (`SurfaceMesh`).
- Prescribed vector advection (`AdvectionTerm`), prescribed normal motion
  (`NormalMotionTerm`), and curvature-driven motion (`CurvatureMotionTerm`).
- Fixed-connectivity vertex motion (no remeshing).
- Forward-Euler, RK2, and RK3 time integrators.
- Equal-arclength redistribution on curves (`CurveEqualArcRedistributor`).
- Experimental tangential redistribution on surfaces (`SurfaceTangentialRedistributor`).
- Piecewise-linear field transfer on curves.

See also: `FrontIntrinsicOps` for static geometry and DEC operators.
"""
module FrontTrackingMethods

using LinearAlgebra
using StaticArrays
using FrontIntrinsicOps

# Bring key FrontIntrinsicOps names into scope
const CurveMesh      = FrontIntrinsicOps.CurveMesh
const SurfaceMesh    = FrontIntrinsicOps.SurfaceMesh
const CurveGeometry  = FrontIntrinsicOps.CurveGeometry
const SurfaceGeometry = FrontIntrinsicOps.SurfaceGeometry
const CurveDEC       = FrontIntrinsicOps.CurveDEC
const SurfaceDEC     = FrontIntrinsicOps.SurfaceDEC
const MeshTopology   = FrontIntrinsicOps.MeshTopology

# Convenience aliases for FrontIntrinsicOps functions used across the package
const compute_geometry = FrontIntrinsicOps.compute_geometry
const build_topology   = FrontIntrinsicOps.build_topology
const is_closed        = FrontIntrinsicOps.is_closed
const measure          = FrontIntrinsicOps.measure
const enclosed_measure = FrontIntrinsicOps.enclosed_measure

# ── Source files (order matters: types first) ─────────────────────────────────
include("front_types.jl")
include("frontfield.jl")
include("frontstate.jl")
include("utils.jl")
include("geometry_refresh.jl")
include("frontterms.jl")
include("timestepping.jl")
include("redistribution.jl")
include("transfer.jl")
include("frontequation.jl")
include("diagnostics.jl")

# ── Public API ────────────────────────────────────────────────────────────────

export
    # Types
    FrontField,
    FrontState,
    FrontEquation,
    AdvectionTerm,
    NormalMotionTerm,
    CurvatureMotionTerm,
    ForwardEuler,
    RK2,
    RK3,
    NoRedistribution,
    CurveEqualArcRedistributor,
    SurfaceTangentialRedistributor,

    # Core functions
    integrate!,
    current_state,
    current_time,
    refresh_geometry!,
    redistribute!,
    repair_front!,
    transfer_vertex_field!,
    transfer_fields!,
    compute_cfl,
    add_field!,
    get_field,
    vertex_coordinates,
    set_vertex_coordinates!,
    location,
    mesh,

    # Diagnostics
    min_edge_length,
    max_edge_length,
    mean_edge_length,
    edge_length_spread,
    front_enclosed_measure,
    front_measure,
    front_centroid,
    front_spacing,
    check_front_validity,
    normal_tangential_decomposition

end # module FrontTrackingMethods
