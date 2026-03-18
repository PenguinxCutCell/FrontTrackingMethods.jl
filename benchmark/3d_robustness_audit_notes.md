# 3-D robustness audit notes (pre-change baseline)

Scope audited before edits:
- `src/remeshing.jl` (`ExperimentalSurfaceRemesher`)
- `src/redistribution.jl` (`SurfaceTangentialRedistributor`)
- `src/transfer.jl` (surface `:barycentric` transfer)
- `src/diagnostics.jl`
- examples: `examples/translate_sphere.jl`, `examples/enright_deformation_3d.jl`
- tests: `test_surface_redistribution.jl`, `test_transfer_surface.jl`, `test_benchmark_rigid_3d.jl`, `test_benchmark_enright_3d.jl`

## Weak points found

1. **Surface remesher behavior was too limited**
   - `ExperimentalSurfaceRemesher` only called `_tangential_smooth_step!` repeatedly.
   - No explicit edge-length thresholds (`lmin`, `lmax`), no local degeneracy rejection, no quality-based trigger logic.
   - No optional volume correction.

2. **Surface diagnostics were underpowered for 3-D regression**
   - Existing diagnostics focused on generic edge-length helpers.
   - Missing reusable surface summaries for edge/area/angle/aspect/degeneracy.
   - This made it hard to assert mesh-quality trends in CI and benchmark scripts.

3. **Surface transfer validation was narrow**
   - Existing tests covered constants and one smooth field (`f=z`) only.
   - No transfer checks after remeshing, rigid rotation, or deformation snapshots.
   - No multi-field evidence that barycentric transfer is consistently at least as good as nearest.

4. **3-D benchmark tests were mostly smoke-level**
   - Rigid 3-D tests had coarse tolerances and no quality assertions.
   - Enright 3-D test asserted only no-NaN + broad volume tolerance + generic reasonableness.
   - No comparative assertion for remeshing-enabled vs no-remeshing runs.

5. **Zalesak sphere specifics**
   - `make_zalesak_sphere_surface` is intentionally approximate and **not closed** (`is_closed(mesh) == false` at coarse refinement).
   - Therefore enclosed-volume metrics are not valid for this geometry in current form.
   - Practical metrics should focus on area drift, centroid drift, front-distance, and slot-preservation proxies.

## Baseline observations from quick probes

- `make_zalesak_sphere_surface(refinement=2)`:
  - `closed=false`
  - no degenerate triangles at `atol=1e-12`
- Rigid rotation on coarse Zalesak sphere (no remeshing):
  - centroid drift around `2.66e-2`
  - relative area drift around `3.26e-2`
  - symmetric Hausdorff around `3.37e-2`

These baseline values were used to set coarse but meaningful CI thresholds.
