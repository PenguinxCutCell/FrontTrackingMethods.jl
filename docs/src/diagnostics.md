# Diagnostics

This page defines the front-tracking diagnostic metrics available in
`FrontTrackingMethods.jl` v0.2 and explains why each is necessary.

---

## Why volume preservation alone is not enough

A common mistake in interface-tracking is to report only volume / area
preservation.  Volume preservation is necessary but not sufficient:

- A highly distorted front that has lost its physical shape can still have
  the correct volume if marker points have been redistributed aggressively.
- Volume-only metrics do not penalize shape diffusion, loss of sharp corners,
  filament tearing, or numerical dissipation of the front geometry.

For a complete assessment of an interface-tracking method, you need:

1. **Volume / area conservation** — global, easy to compute, necessary.
2. **Centroid position error** — detects systematic drift or spurious translation.
3. **Front-to-front geometric error** — detects shape diffusion and filament loss.
4. **Mesh quality metrics** — detect marker collapse / degenerate elements.

---

## Shape-integral diagnostics

### `relative_area_error(mesh, mesh_ref)`

Relative signed area difference:
```
(area(mesh) - area(mesh_ref)) / area(mesh_ref)
```
For curves.  Positive means the mesh has grown; negative means it has shrunk.

### `relative_volume_error(mesh, mesh_ref)`

Same as above for surfaces (enclosed volume).

### `centroid_error(mesh, mesh_ref)`

Euclidean distance between the vertex centroids of the two meshes.
```
‖centroid(mesh) - centroid(mesh_ref)‖₂
```

### `perimeter_error(mesh, mesh_ref)`

Relative difference in total arc length (curve perimeter).

### `surface_area_error(mesh, mesh_ref)`

Relative difference in total surface area.

---

## Front-to-front geometric error

These metrics measure how well two front meshes agree geometrically, without
assuming vertex-to-vertex correspondence.  They are essential for reversible
benchmark tests where the initial and final fronts should coincide.

### `nearest_distance_curve_to_curve(meshA, meshB)`

For each vertex of `meshA`, find the minimum distance to any segment of
`meshB`.  Returns the maximum such distance (one-sided Hausdorff).

### `symmetric_hausdorff_curve(meshA, meshB)`

```
max(nearest_distance_curve_to_curve(A,B),
    nearest_distance_curve_to_curve(B,A))
```
The approximate symmetric Hausdorff distance between two curves.

### `l2_distance_curve(meshA, meshB; nsamples=1000)`

An approximate L² distance between two curves computed by sampling `nsamples`
uniformly-spaced points on `meshA` and summing squared nearest distances to
`meshB`.

### 3-D equivalents

`nearest_distance_surface_to_surface`, `symmetric_hausdorff_surface`, and
`l2_distance_surface` provide the same functionality for surface meshes.

---

## Mesh quality metrics

### `edge_length_spread(state)` / `edge_length_stats(mesh)`

The ratio `max_edge_length / min_edge_length`.  A spread close to 1 means
uniform spacing; a large spread indicates clustering that will degrade
numerical accuracy or cause stability issues.

### `min_edge_length(state)` / `max_edge_length(state)` / `mean_edge_length(state)`

Basic edge-length statistics.

### `triangle_quality_stats(mesh::SurfaceMesh)`

Returns minimum triangle angle and mean aspect ratio.  Low minimum angle
(< 5°) or high aspect ratio (> 20) indicate degenerate elements that will
cause integration errors.

### `surface_edge_length_stats(mesh)`

Returns `(min, max, mean, std, ratio)` for surface edge lengths, with
`ratio = max/min`.

### `surface_triangle_area_stats(mesh)`

Returns `(min, max, mean)` for triangle areas.

### `surface_triangle_angle_stats(mesh)`

Returns `(min_angle, max_angle, mean_min_angle)` in radians.

### `surface_aspect_ratio_stats(mesh)`

Returns `(min, max, mean)` using

```
aspect = (lmax * perimeter) / (4 * sqrt(3) * area)
```

for each triangle (1 for equilateral).

### `surface_degenerate_fraction(mesh; atol=...)`

Fraction of triangles with area `<= atol`.

### `surface_quality_summary(mesh)`

Compact named-tuple aggregation of edge, area, angle, aspect, and degenerate
fraction metrics.

### `surface_normal_consistency(mesh, geom)`

Optional orientation-consistency diagnostic based on alignment between face
normals and averaged vertex normals.

### `check_front_validity(state; warn=true)`

Returns `true` if the front passes basic sanity checks:
- no NaN / Inf vertices
- correct number of edges / faces
- non-zero area / volume

---

## Why front-to-front geometric error matters

Consider the Rider–Kothe reversed vortex benchmark.  A method may preserve
area nearly perfectly while completely destroying the filament shape:

- **No redistribution:** the interface stretches into a filament and markers
  cluster near the tips.  Volume is preserved because markers cannot cross.
  But the filament becomes unresolved and the Hausdorff distance to the
  initial front is large.

- **With redistribution:** markers are spread uniformly along the filament.
  The shape is better resolved and the Hausdorff distance is materially
  smaller.

This is why `symmetric_hausdorff_curve` is the primary shape metric in
the reversible benchmark tests, not just area drift.
