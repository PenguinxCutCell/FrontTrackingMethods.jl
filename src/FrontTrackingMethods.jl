"""
    FrontTrackingMethods

A package for evolving explicit front-tracking meshes in time.

Mirrors the user workflow of `LevelSetMethods.jl`, but operates on Lagrangian
polygonal curves and triangulated surfaces built on top of `FrontIntrinsicOps.jl`.

v0.2 supports:
- Closed 2-D polygonal curves (`CurveMesh`) and closed 3-D surfaces (`SurfaceMesh`).
- Prescribed vector advection (`AdvectionTerm`), prescribed normal motion
  (`NormalMotionTerm`), and curvature-driven motion (`CurvatureMotionTerm`).
- Fixed-connectivity vertex motion (no topology change).
- Forward-Euler, RK2, and RK3 time integrators.
- Equal-arclength redistribution on curves (`CurveEqualArcRedistributor`).
- Adaptive tangential curve remeshing with corner protection (`AdaptiveCurveRemesher`).
- Experimental tangential redistribution on surfaces (`SurfaceTangentialRedistributor`,
  `ExperimentalSurfaceRemesher`).
- Piecewise-linear field transfer on curves.
- Barycentric field transfer on surfaces (`transfer_vertex_field!` with `:barycentric`).
- Standard benchmark geometries (`make_circle_benchmark_curve`, `make_zalesak_disk_curve`,
  `make_sphere_benchmark_surface`, `make_zalesak_sphere_surface`).
- Standard benchmark velocity fields (`rigid_translation_velocity`,
  `rigid_rotation_2d`, `rigid_rotation_3d`, `rider_kothe_single_vortex`,
  `deformation16_2d`, `serpentine_2d`, `enright_3d`).
- Front-to-front comparison metrics (Hausdorff, L², area/volume errors).
- Lightweight callbacks (`EveryNSteps`, `TimeIntervalCallback`).

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
include("remeshing.jl")
include("transfer.jl")
include("frontequation.jl")
include("diagnostics.jl")
include("benchmark_geometries.jl")
include("benchmark_fields.jl")
include("compare_fronts.jl")
include("callbacks.jl")

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
    AdaptiveCurveRemesher,
    ExperimentalSurfaceRemesher,
    EveryNSteps,
    TimeIntervalCallback,

    # Core functions
    integrate!,
    current_state,
    current_time,
    refresh_geometry!,
    redistribute!,
    remesh!,
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
    normal_tangential_decomposition,

    # Benchmark geometry constructors (v0.2)
    make_circle_benchmark_curve,
    make_zalesak_disk_curve,
    make_sphere_benchmark_surface,
    make_zalesak_sphere_surface,

    # Benchmark velocity fields (v0.2)
    rigid_translation_velocity,
    rigid_rotation_2d,
    rigid_rotation_3d,
    rider_kothe_single_vortex,
    deformation16_2d,
    serpentine_2d,
    enright_3d,

    # Front comparison metrics (v0.2)
    nearest_distance_curve_to_curve,
    symmetric_hausdorff_curve,
    l2_distance_curve,
    nearest_distance_surface_to_surface,
    symmetric_hausdorff_surface,
    l2_distance_surface,
    relative_area_error,
    relative_volume_error,
    centroid_error,
    perimeter_error,
    surface_area_error,

    # Callbacks (v0.2)
    compose_callbacks

end # module FrontTrackingMethods
