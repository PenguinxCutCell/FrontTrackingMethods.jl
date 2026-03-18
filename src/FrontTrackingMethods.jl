"""
    FrontTrackingMethods

A package for evolving explicit front-tracking meshes in time.

Mirrors the user workflow of `LevelSetMethods.jl`, but operates on Lagrangian
polygonal curves and triangulated surfaces built on top of `FrontIntrinsicOps.jl`.

v0.2 supports:
- Closed 2-D polygonal curves (`CurveMesh`) and closed 3-D surfaces (`SurfaceMesh`).
- Limited open 2-D curve support for topology/geometry and prescribed-advection cases.
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
- Standard benchmark geometries (`make_circle_benchmark_curve`, `make_open_arc_benchmark_curve`,
  `make_zalesak_disk_curve`, `make_sphere_benchmark_surface`, `make_zalesak_sphere_surface`).
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
include("topology_state.jl")
include("utils.jl")
include("geometry_refresh.jl")
include("frontterms.jl")
include("timestepping.jl")
include("redistribution.jl")
include("remeshing.jl")
include("transfer.jl")
include("topology_change.jl")
include("topology_detection.jl")
include("topology_patch.jl")
include("topology_rasterize.jl")
include("topology_reconstruct.jl")
include("topology_replace.jl")
include("topology_transfer.jl")
include("frontequation.jl")
include("diagnostics.jl")
include("benchmark_geometries.jl")
include("benchmark_fields.jl")
include("compare_fronts.jl")
include("callbacks.jl")
include("plotting_stubs.jl")

# ── Public API ────────────────────────────────────────────────────────────────

export
    # Types
    AbstractFrontState,
    FrontField,
    FrontState,
    FrontComponentState,
    MultiFrontState,
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
    AbstractTopologyHandler,
    NoTopologyChange,
    LocalCartesianTopologyHandler,
    TopologyEventReport,
    TopologyCandidate,
    EventPatch,
    CartesianPatch,
    EveryNSteps,
    TimeIntervalCallback,

    # Core functions
    integrate!,
    current_state,
    current_time,
    ncomponents,
    component,
    component_mesh,
    component_geom,
    component_fields,
    all_meshes,
    eachcomponent,
    map_components,
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
    handle_topology_change!,
    find_topology_candidates,
    select_topology_candidate,
    extract_event_patch,
    make_patch_grid,
    patch_cell_centers,
    rasterize_indicator!,
    reconstruct_curve_from_patch,
    reconstruct_surface_from_patch,
    replace_components!,
    replace_curve_patch!,
    transfer_fields_after_topology_change!,
    local_length_scale,
    detect_imminent_merge_2d,
    detect_imminent_self_merge_2d,
    detect_imminent_split_2d,
    detect_imminent_merge_3d,
    detect_imminent_self_merge_3d,
    detect_imminent_split_3d,
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
    surface_edge_length_stats,
    surface_triangle_area_stats,
    surface_triangle_angle_stats,
    surface_aspect_ratio_stats,
    surface_quality_summary,
    surface_degenerate_fraction,
    surface_normal_consistency,

    # Benchmark geometry constructors (v0.2)
    make_circle_benchmark_curve,
    make_open_arc_benchmark_curve,
    make_zalesak_disk_curve,
    make_sphere_benchmark_surface,
    make_zalesak_sphere_surface,
    make_two_circles_merge_setup,
    make_dumbbell_curve_setup,
    make_peanut_curve_setup,
    make_two_spheres_merge_setup,
    make_dumbbell_surface_setup,

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
    compose_callbacks,

    # Optional Makie plotting API (implemented in ext/MakieExt.jl)
    makie_theme,
    set_makie_theme!,
    plot_front,
    plot_state,
    plot_equation,
    animate_equation!,
    record_evolution!,
    snapshot

end # module FrontTrackingMethods
