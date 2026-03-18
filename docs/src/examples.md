# Examples

Short descriptions of example scripts included in the repository (`examples/`).

- `curve_shortening_circle.jl` — curve-shortening (mean-curvature) flow on a circle.
- `deformation16_2d.jl` — 2-D 16-vortex deformation benchmark setup.
- `enright_deformation_3d.jl` — Enright 3-D deformation smoke test on a coarse sphere.
- `expand_circle_normal_speed.jl` — outward normal motion example.
- `redistribute_curve_demo.jl` — demonstrates curve redistribution strategies.
- `rider_kothe_single_vortex.jl` — Rider–Kothe reversible vortex benchmark.
- `rotate_circle.jl` — circle rigid rotation demo.
- `rotate_sphere.jl` — rigid rotation demo for a sphere.
- `serpentine_2d.jl` — serpentine deformation benchmark setup.
- `translate_circle.jl` — translation demo for a circle.
- `translate_sphere.jl` — translation demo for a sphere.
- `topology_merge_two_circles_2d.jl` — 2-D merge event via local Cartesian reconstruction.
- `topology_split_dumbbell_2d.jl` — 2-D split event on a dumbbell-like curve.
- `topology_merge_two_spheres_3d.jl` — coarse 3-D merge prototype smoke example.
- `topology_split_dumbbell_3d.jl` — coarse 3-D split prototype smoke example.
- `zalesak_disk_rotation.jl` — Zalesak disk rotation demo.
- `zalesak_sphere_rotation.jl` — coarse 3-D slotted-sphere rigid-rotation benchmark.

To run an example interactively:

```bash
julia --project=. examples/rotate_circle.jl
```

Examples that evolve fronts now use `CairoMakie` and write outputs under:

- `examples/output/<example_name>/initial.png`
- `examples/output/<example_name>/final.png`
- `examples/output/<example_name>/animation.mp4`

Examples are intentionally small; for longer convergence or benchmarking
scripts see the `test/` harness and the `benchmark_*.jl` helpers in `src/`.
