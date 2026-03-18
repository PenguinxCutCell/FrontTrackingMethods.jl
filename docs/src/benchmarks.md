# Benchmarks

This page describes the standard front-tracking benchmarks implemented in
`FrontTrackingMethods.jl` v0.2.  Each benchmark is available through the
reusable helpers in `src/benchmark_geometries.jl` and
`src/benchmark_fields.jl`, which are called by the test suite and by the
example scripts in `examples/`.

---

## 2-D Rigid Motions

### Circle rotation (rigid body)

**Geometry:** circle centered at `(0.5, 0.75)`, radius `0.15`.  
**Field:** `rigid_rotation_2d(; center=(0.5,0.5), omega=2π)` — one full revolution per unit time.  
**Final time:** `T = 1.0`.  
**What it stresses:** time-integration accuracy, area preservation of the polygon.  
**Primary metrics:** centroid return error, relative area drift, Hausdorff distance to initial.  
**Exact reversibility:** yes — all errors should be due to time-integration truncation only.

### Circle translation

**Field:** `rigid_translation_velocity(u)`.  
**What it stresses:** pure transport accuracy (all errors should be machine precision).  
**Primary metrics:** centroid error, area error.

---

## 2-D Zalesak disk rotation

**Reference:** Zalesak (1979), "Fully Multidimensional Flux-Corrected Transport".

**Geometry:** slotted disk centered at `(0.5, 0.75)`, radius `0.15`, slot width `0.05`, slot depth `0.25`. Built by `make_zalesak_disk_curve`.  
**Field:** rigid-body rotation about `(0.5, 0.5)`, `ω = 2π`.  
**Final time:** `T = 1.0`.  
**What it stresses:** preservation of sharp slot corners, area conservation under redistribution.  
**Primary metrics:** area drift, centroid drift, Hausdorff distance to initial front.  
**Expected behavior:** the slot opening should remain clearly recognizable after one revolution.

---

## 2-D Rider–Kothe Reversed Single Vortex

**Reference:** Rider & Kothe (1998), "Reconstructing Volume Tracking".

**Geometry:** circle centered at `(0.5, 0.75)`, radius `0.15`.  
**Field:** `rider_kothe_single_vortex(; T=2.0)` — stream-function-based solenoidal field, reversed at `T/2`.  
**Final time:** `T = 2.0`.  
**What it stresses:** interface stretching and filament formation; reversibility under time reversal; redistribution quality.  
**Primary metrics:** shape recovery (Hausdorff / L² distance to initial), area drift.  
**Exact reversibility:** yes in the exact case; numerical errors are non-zero.  
**Note:** redistribution/remeshing should materially improve final-time shape recovery.

---

## 2-D Strong Deformation (16-Vortex)

**Geometry:** circle centered at `(0.5, 0.75)`, radius `0.15`.  
**Field:** `deformation16_2d(; T=2.0)` — periodic stream-function field with 16 vortices, reversed at `T/2`.  
**What it stresses:** severe shear-stretch, robustness of redistribution under strong deformation.  
**Primary metrics:** minimum edge length (collapse indicator), shape recovery error.

---

## 2-D Serpentine Deformation

**Reference:** LeVeque (1996), "High-resolution conservative algorithms for advection in incompressible flow".

**Geometry:** circle centered at `(0.5, 0.75)`, radius `0.15`.  
**Field:** `serpentine_2d(; Tmax=3.0)` — divergence-free serpentine flow, reversed at `Tmax/2`.  
**Final time:** `Tmax = 3.0`.  
**What it stresses:** long-filament formation, severe shear, marker redistribution under extreme stretch.  
**Primary metrics:** L² and Hausdorff shape recovery, relative area drift, minimum edge length.  
**Note:** this is the most demanding 2-D reversible test.  Redistribution is mandatory to prevent catastrophic marker collapse.

---

## 3-D Rigid Motions

### Sphere translation

**Geometry:** sphere centered at `(0.35, 0.35, 0.35)`, radius `0.15`.  
**Field:** `rigid_translation_velocity(u)`.  
**Primary metrics:** centroid error (should be near machine precision), volume drift.

### Sphere rotation

**Geometry:** sphere centered at `(0.5, 0.75, 0.5)`, radius `0.15`.  
**Field:** `rigid_rotation_3d(; center=(0.5,0.5,0.5), axis=:z, omega=2π)`.  
**Final time:** `T = 1.0`.  
**Primary metrics:** centroid return error, volume drift, Hausdorff to initial.

---

## 3-D Zalesak Sphere Rotation

**Geometry:** slotted sphere built by `make_zalesak_sphere_surface`.  
**Field:** rigid-body rotation.  
**What it stresses:** fixed-topology surface tracking under rigid motion, slot-corner preservation on a triangulated mesh.  
**Primary metrics (v0.2):** surface-area drift, centroid drift, symmetric Hausdorff to initial, slot-region proxy distance.  
**Important limitation:** the current slotted-sphere constructor is approximate and not guaranteed closed at coarse resolution, so enclosed-volume drift is not the primary metric here.

---

## 3-D Enright Deformation

**Reference:** Enright et al. (2002), "A Hybrid Particle Level Set Method for Improved Interface Tracking".

**Geometry:** sphere centered at `(0.35, 0.35, 0.35)`, radius `0.15`, in a unit cube.  
**Field:** `enright_3d(; T=3.0)` — divergence-free prescribed velocity, reversed at `T/2`.  
**Final time:** `T = 3.0`.  
**What it stresses:** severe surface deformation into thin sheets; volume preservation; surface mesh quality degradation.  
**Primary metrics:**
- final volume drift
- surface area drift
- approximate symmetric Hausdorff distance to initial sphere
- quality over time (minimum triangle angle, degenerate-triangle fraction)
**Which metric to check first:** for Enright, check symmetric Hausdorff (shape recovery) first, then volume drift, then quality-over-time.  
**Note:** surface remeshing is often helpful but remains experimental; benchmark studies should compare remeshing-enabled and no-remeshing runs.
