# Codex instruction for v0.3 topology-change work

You are working on `PenguinxCutCell/FrontTrackingMethods.jl`.

## Goal

Implement v0.3 topology change support for explicit front tracking.

## Important scope decision

Do NOT try to build a fully general, production-grade 3-D breakup/coalescence framework in one PR.

v0.3 should deliver:
1. a clean architecture for topology-aware states and equations,
2. a robust 2-D topology-change pipeline,
3. a first coarse 3-D topology-change prototype,
4. clear benchmark cases and diagnostics,
5. honest docs about what is strong vs experimental.

## Core design choice

Use a hybrid local-reconstruction strategy:

- normal evolution remains explicit front tracking on `CurveMesh` / `SurfaceMesh`,
- imminent topology events are detected geometrically,
- the event is localized to a component or local patch,
- the interface in that region is rasterized to an auxiliary Cartesian field,
- marching squares (2-D) or marching cubes (3-D prototype) reconstructs the post-event interface,
- the reconstructed component(s) replace the old one(s),
- time stepping continues explicitly.

Do NOT start with general direct mesh-surgery-only topology change.
Do NOT start with a physical film-drainage coalescence model.
Those can come later.

## Non-goals for v0.3

- no full physical film drainage / rupture model
- no rebound-vs-coalescence collision physics
- no non-manifold or triple-line topology handling
- no contact-angle/contact-line topology changes
- no cut-cell bulk coupling
- no multiphysics field coupling across topology change
- no perfect 3-D robustness at first release

## High-level release target

v0.3 should mean:

- 2-D:
  - two loops can merge into one,
  - one loop can pinch off into two,
  - local reconstruction is robust enough for standard tests,
  - attached scalar fields survive the event reasonably.

- 3-D:
  - first prototype for merge/split using reconstruction,
  - coarse benchmark demonstrations only,
  - clearly marked experimental.

## Part A — architecture changes

### A1. Replace single-front-only state logic with multi-component state logic

The current `FrontState` stores exactly one `CurveMesh` or `SurfaceMesh`.
That is no longer enough once one component can split into two or two can merge into one.

Add a new state layer instead of breaking everything at once.

Implement:

```julia
abstract type AbstractFrontState end

mutable struct FrontComponentState{M,G,D}
    mesh::M
    geom::G
    dec::D
    fields::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end

mutable struct MultiFrontState <: AbstractFrontState
    components::Vector{FrontComponentState}
    t::Float64
    global_fields::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end
```

Requirements:

- one connected curve/surface = one `FrontComponentState`
- `MultiFrontState` owns the simulation time
- topology change modifies the `components` vector
- existing fixed-topology cases should still work cleanly

Important compatibility rule:

- do NOT delete `FrontState` immediately
- keep it for v0.3 as a convenience wrapper or legacy single-component state
- add conversion helpers:
  - `MultiFrontState(state::FrontState)`
  - `FrontState(multistate::MultiFrontState)` only when `length(components)==1`

### A2. Add state-level accessors

Add:

- `current_time(state::MultiFrontState)`
- `ncomponents(state)`
- `component(state, i)`
- `component_mesh(state, i)`
- `component_geom(state, i)`
- `component_fields(state, i)`
- `all_meshes(state)`

Add iteration helpers:

- `eachcomponent(state)`
- `map_components(f, state)`

### A3. Geometry refresh on multi-component states

Implement:

- `refresh_geometry!(state::MultiFrontState; build_dec=false, which=:all)`

Requirements:

- refresh all components or only selected component indices
- preserve component ordering unless explicitly changed
- keep cache invalidation simple and safe

### A4. Equation compatibility

Extend `FrontEquation` so it can evolve:

- `FrontState`
- `MultiFrontState`

Strategy:

- keep `FrontEquation` user API as stable as possible
- allow `front=` to accept either a single mesh or a vector of meshes
- if vector of meshes is provided, build a `MultiFrontState`

Recommended API:

```julia
eq = FrontEquation(; terms=..., front=mesh)
eq = FrontEquation(; terms=..., front=[mesh1, mesh2])
eq = FrontEquation(; terms=..., state=multistate)
```

Do not redesign the whole integrator interface.

## Part B — topology handler abstraction

### B1. Add topology-handler types

Create a new file:

- `src/topology_change.jl`

Define:

```julia
abstract type AbstractTopologyHandler end

struct NoTopologyChange <: AbstractTopologyHandler end

struct LocalCartesianTopologyHandler <: AbstractTopologyHandler
    d_trigger::Float64
    d_merge::Float64
    d_split::Float64
    patch_h_factor::Float64
    patch_margin_factor::Float64
    reconstruct_scope::Symbol   # :local_patch or :whole_component
    use_indicator::Bool
    preserve_fields::Bool
    max_components_created::Int
end
```

Defaults:

- conservative, stable, easy to reason about
- choose `:local_patch` in 2-D by default
- choose `:whole_component` or `:event_box` in 3-D prototype if that is more robust

### B2. Hook topology handling into the time integrator

Add a topology-change stage to `integrate!`.

At a high level:

1. advance one explicit step/stage as usual,
2. refresh geometry,
3. call `handle_topology_change!(state, handler)` if handler is not `NoTopologyChange`,
4. refresh geometry again if topology changed,
5. continue.

Add a small event report type:

```julia
struct TopologyEventReport
    changed::Bool
    event_type::Symbol      # :none, :merge, :split, :reconnect, :ambiguous
    component_ids::Vector{Int}
    patch_bbox
    n_old_components::Int
    n_new_components::Int
    notes::String
end
```

Return this report from topology handling.
Store the most recent report in `state.cache` or `eq.buffers`.

Do not overbuild a full event system.
Keep it small.

## Part C — geometric detection of imminent topology change

### C1. Detection philosophy

Implement purely geometric detection first.

Do NOT try to decide coalescence from film physics in v0.3.
The first detector only answers:

- are two non-neighboring front regions close enough and approaching enough that a merge is imminent?
- is a neck thin enough that a split is imminent?
- is the current front already self-intersecting or about to self-intersect?

### C2. 2-D detectors

Implement:

- `detect_imminent_merge_2d(componentA, componentB, handler)`
- `detect_imminent_self_merge_2d(component, handler)`
- `detect_imminent_split_2d(component, handler)`

Use geometric criteria based on:

1. segment-segment distance between non-adjacent edges,
2. vertex-segment distance where useful,
3. approximate local thickness,
4. predicted crossing if velocities are available,
5. normal alignment.

Recommended criteria:

- gap distance `d < η * h_loc`
- normals approximately opposite for film/neck situations
- relative normal velocity negative for approaching sheets if available
- exclude directly adjacent edges/vertices from self-detection

Define a helper:

- `local_length_scale(component, feature_ids...)`

### C3. 3-D detectors

Implement prototype versions:

- `detect_imminent_merge_3d(componentA, componentB, handler)`
- `detect_imminent_self_merge_3d(component, handler)`
- `detect_imminent_split_3d(component, handler)`

Use:

- triangle-triangle near-contact
- edge-edge near-contact
- vertex-face near-contact
- neck-thickness heuristics
- optional bounding-volume acceleration

Important:

- start conservative
- false negatives are better than many false positives
- if robust 3-D self-contact detection is difficult, start with:
  - component-component merge detection,
  - coarse neck-thickness split detection,
  - self-intersection smoke detection only

### C4. Candidate event structure

Define:

```julia
struct TopologyCandidate
    event_type::Symbol         # :merge or :split
    component_ids::Vector{Int}
    feature_pairs::Any
    distance::Float64
    h_local::Float64
    score::Float64
    bbox
end
```

Create:

- `find_topology_candidates(state, handler)`
- `select_topology_candidate(candidates)`

At first:

- pick the single strongest candidate per step
- do not resolve simultaneous multiple events in one pass
- after one event, refresh and continue

## Part D — event localization

### D1. Patch extraction API

Implement:

- `extract_event_patch(state, candidate, handler)`

Patch object:

```julia
struct EventPatch
    component_ids::Vector{Int}
    bbox
    margin_bbox
    anchor_data
    local_mesh_data
    patch_grid
end
```

The patch should contain:

- involved components
- a bounding box around the candidate features
- a safety margin
- the subset of front elements in the patch
- anchor information for reconnecting to untouched regions if using local-patch mode

### D2. Scope modes

Support two modes:

1. `:local_patch`
   - reconstruct only a local portion of a component
   - harder to stitch
   - cheaper
   - preferred target for 2-D

2. `:whole_component`
   - reconstruct the full involved component(s)
   - easier logic
   - more expensive
   - acceptable for first 3-D prototype

Recommended v0.3 plan:

- 2-D: support `:local_patch` and `:whole_component`
- 3-D: support `:whole_component` first, add true local patch later only if easy

### D3. Anchor strategy for local patch mode

For 2-D local patch reconstruction:

- identify curve vertices/edges entering and leaving the patch box
- freeze a small outer ring of front data as anchors
- reconstruct inside, then reconnect to those anchors

Do NOT attempt fancy C1 stitching at first.
Piecewise-linear reconnection is acceptable if robust.

## Part E — auxiliary Cartesian reconstruction

### E1. Patch grid

Implement a small Cartesian patch grid type:

```julia
struct CartesianPatch{T}
    xmin::SVector
    xmax::SVector
    h::T
    dims::NTuple
    centers
end
```

Add:

- `make_patch_grid(bbox, h)`
- `patch_cell_centers(grid)`

Choose patch resolution from local mesh scale:

- `h_patch = handler.patch_h_factor * h_local`
- clamp to reasonable min/max if needed

### E2. Rasterization target

Support at least one scalar representation on the patch:

- indicator field `χ`
  Optionally later:
- signed distance `φ`

Recommended v0.3 default:

- indicator / occupancy reconstruction first

Implement:

- `rasterize_indicator!(χ, patch, local_front_data)`

2-D:

- use winding/ray-crossing or polygon fill logic

3-D:

- use inside/outside tests from oriented triangles
- if exact robust sign tests are difficult, use a practical ray-casting method with tolerances

Do not overengineer exact computational geometry initially.
Prefer robust-enough, well-tested rasterization.

### E3. Reconstruction

Implement:

- `reconstruct_curve_from_patch(χ, patch)` using marching squares
- `reconstruct_surface_from_patch(χ, patch)` using marching cubes

Requirements:

- extract one or more reconstructed connected components
- orient each component consistently
- remove tiny spurious components below a threshold
- optionally smooth one pass if necessary, but avoid destroying volume/area too much

Output should be:

- vector of new `CurveMesh` or `SurfaceMesh` objects

Do NOT rely on a global reconstruction of the whole domain in v0.3.
Keep reconstruction local to the event scope.

### E4. Post-reconstruction cleanup

After extracting reconstructed components:

- run a quality cleanup / redistribution / remeshing pass
- ensure no duplicate vertices
- ensure no zero-length edges / zero-area triangles
- orient components consistently
- refresh geometry

For 2-D:

- ensure closed curves
- ensure orientation convention is restored

For 3-D:

- ensure triangles are consistently oriented
- ensure closedness if the benchmark expects closed surfaces
- remove tiny disconnected artifacts if below a documented threshold

## Part F — replacing components in the state

### F1. Whole-component replacement

Implement:

- `replace_components!(state::MultiFrontState, old_ids, new_components)`

Behavior:

- remove old components
- insert new reconstructed components
- preserve other unaffected components
- keep time unchanged
- refresh geometry for inserted components

### F2. Local-patch replacement for curves

Implement a 2-D-specific path:

- `replace_curve_patch!(component, patch, reconstructed_curves, anchors)`

Recommended first target:

- support local patch replacement only for curves
- support whole-component replacement for surfaces in v0.3

### F3. Field transfer across topology events

When replacing components, attached fields must be transferred.

Implement:

- `transfer_fields_after_topology_change!(old_components, new_components; method=...)`

For v0.3:

- support vertex-based scalar fields first
- use closest-point or barycentric transfer where applicable
- preserve constant fields exactly if possible
- document limitations clearly

If field transfer becomes difficult in ambiguous multi-component splits:

- keep scalar fields only
- error clearly for unsupported field types

## Part G — tests and benchmark cases

### G1. 2-D benchmark geometries for topology change

Add reusable geometry constructors:

- `make_two_circles_merge_setup(...)`
- `make_dumbbell_curve_setup(...)`
- `make_peanut_curve_setup(...)`
- `make_figure8_like_near_contact_setup(...)` if useful

These should return meshes or vectors of meshes suitable for `MultiFrontState`.

### G2. 2-D test cases

Add a dedicated test file set:

- `test_topology_state.jl`
- `test_topology_detection_2d.jl`
- `test_topology_patch_2d.jl`
- `test_topology_reconstruction_2d.jl`
- `test_topology_merge_2d.jl`
- `test_topology_split_2d.jl`
- `test_topology_field_transfer_2d.jl`

Required tests:

1. Multi-component state
   - create state with two circles
   - component count correct
   - accessors work
   - geometry refresh works

2. Merge detection
   - two circles moving together trigger a merge candidate
   - far-apart circles do not trigger

3. Split detection
   - dumbbell curve with thin neck triggers split candidate
   - thick peanut shape does not trigger too early

4. Patch rasterization
   - indicator field is correct on a small synthetic patch

5. Marching-squares reconstruction
   - reconstruct simple closed loop
   - reconstruct split into two loops
   - reconstruct merge into one loop

6. End-to-end 2-D merge
   - two circles advect into each other
   - after event, one component remains
   - area is approximately conserved
   - no NaNs
   - shape is reasonable

7. End-to-end 2-D split
   - dumbbell / peanut under prescribed inward motion or curvature-like shrinking
   - after event, two components remain
   - total area approximately conserved
   - no degenerate edges

8. Field transfer
   - constant field preserved
   - smooth scalar field transferred reasonably through event

### G3. 3-D benchmark geometries

Add:

- `make_two_spheres_merge_setup(...)`
- `make_dumbbell_surface_setup(...)`

Keep them modest and well documented.
Prefer robust coarse geometries over fancy CAD-quality shapes.

### G4. 3-D prototype tests

Add:

- `test_topology_detection_3d.jl`
- `test_topology_reconstruction_3d.jl`
- `test_topology_merge_3d.jl`
- `test_topology_split_3d.jl`

Keep them coarse and forgiving.

Required tests:

1. two spheres approaching trigger a merge candidate
2. whole-component reconstruction on a patch produces one merged surface
3. dumbbell surface split prototype produces two surfaces
4. no NaNs / no catastrophic invalid geometry
5. volume change remains bounded, with loose tolerances

Important:

- 3-D tests should be smoke/regression tests in v0.3
- stronger studies belong in `benchmark/`, not strict CI

## Part H — examples

Add examples:

- `examples/topology_merge_two_circles_2d.jl`
- `examples/topology_split_dumbbell_2d.jl`
- `examples/topology_merge_two_spheres_3d.jl`
- `examples/topology_split_dumbbell_3d.jl`

Each example should:

1. build the topology-capable equation/state
2. choose a `LocalCartesianTopologyHandler`
3. evolve to a target time
4. save initial/final snapshots
5. record an animation if Makie is available
6. print:
   - number of components over time
   - event type
   - area/volume drift
   - quality metrics

Use output paths like:

- `output/<example_name>/initial.png`
- `output/<example_name>/final.png`
- `output/<example_name>/animation.mp4`
- `output/<example_name>/metrics.csv`

## Part I — benchmark scripts

Add:

- `benchmark/run_topology_merge_2d_study.jl`
- `benchmark/run_topology_split_2d_study.jl`
- `benchmark/run_topology_merge_3d_study.jl`
- `benchmark/run_topology_split_3d_study.jl`

Each script should sweep:

- patch resolution factor
- trigger distance
- reconstruction scope
- remeshing on/off after event

Output:

- area/volume drift
- front-to-front distance if a reference exists
- runtime
- number of components
- event count
- quality metrics

Do not make benchmark scripts depend on plotting.

## Part J — docs and README

### J1. Update feature matrix

Mark:

- `MultiFrontState` / multi-component states: ✅ new in v0.3
- 2-D topology change: ✅ experimental/usable
- 3-D topology change: ⚠️ prototype only
- film-model-based coalescence: ❌ not in v0.3

### J2. Add a dedicated docs page

Add `docs/src/topology_change.md`

Explain:

1. why the state model changed
2. geometric detection vs physical coalescence models
3. local patch vs whole-component reconstruction
4. marching squares / marching cubes reconstruction
5. current limitations

Be explicit:

- v0.3 topology change is geometry-driven, not film-physics-driven
- 2-D is the main validated target
- 3-D is first-generation and experimental

### J3. README wording

Add a short section:
“Topology change in v0.3”

State clearly:

- explicit front tracking remains the main representation
- topology events are handled by local Cartesian reconstruction
- 2-D merge/split supported
- 3-D coarse prototype supported
- physical film drainage / rupture models are future work

## Part K — implementation details and practical constraints

### K1. Keep the topology handler additive

Do not break existing no-topology workflows.
Default behavior should remain:

- `topology_handler = NoTopologyChange()`

### K2. Keep 2-D stronger than 3-D

Spend most implementation effort on:

- 2-D detection
- 2-D patching
- 2-D marching-squares reconstruction
- 2-D field transfer across events

3-D only needs to prove the architecture works.

### K3. Prefer conservative triggering

Use thresholds that reduce false positives.
If uncertain:

- do not trigger topology change
- continue explicit evolution

### K4. Handle one event at a time

Do not try to process many simultaneous events in the first version.
Algorithm:

- detect candidates
- pick best candidate
- reconstruct
- refresh state
- continue

### K5. Keep ambiguous cases explicit

If a reconstruction is ambiguous:

- emit `TopologyEventReport(changed=false, event_type=:ambiguous, ...)`
- do not silently corrupt the mesh

### K6. Keep remeshing separate from topology logic

After topology change, call existing remeshing / redistribution tools.
Do not bury mesh-quality restoration inside the detector itself.

## Part L — suggested file changes

Create / modify at least:

src/
- FrontTrackingMethods.jl
- frontstate.jl
- frontequation.jl
- topology_change.jl
- topology_state.jl
- topology_detection.jl
- topology_patch.jl
- topology_rasterize.jl
- topology_reconstruct.jl
- topology_replace.jl
- topology_transfer.jl
- benchmark_geometries.jl   # add merge/split setups
- diagnostics.jl            # add per-component metrics

test/
- runtests.jl
- test_topology_state.jl
- test_topology_detection_2d.jl
- test_topology_patch_2d.jl
- test_topology_reconstruction_2d.jl
- test_topology_merge_2d.jl
- test_topology_split_2d.jl
- test_topology_field_transfer_2d.jl
- test_topology_detection_3d.jl
- test_topology_reconstruction_3d.jl
- test_topology_merge_3d.jl
- test_topology_split_3d.jl

examples/
- topology_merge_two_circles_2d.jl
- topology_split_dumbbell_2d.jl
- topology_merge_two_spheres_3d.jl
- topology_split_dumbbell_3d.jl

benchmark/
- run_topology_merge_2d_study.jl
- run_topology_split_2d_study.jl
- run_topology_merge_3d_study.jl
- run_topology_split_3d_study.jl

docs/src/
- topology_change.md

## Part M — acceptance criteria

This work is done only if all of the following hold:

1. `MultiFrontState` exists and supports multiple components cleanly.
2. Existing single-component workflows still run without topology handling.
3. A topology handler API exists and defaults to `NoTopologyChange`.
4. 2-D merge detection works on two circles approaching each other.
5. 2-D split detection works on a dumbbell/peanut necking case.
6. 2-D marching-squares reconstruction works and is used end-to-end.
7. End-to-end 2-D merge and split examples/tests pass.
8. Constant scalar fields survive topology change reasonably.
9. A first 3-D prototype exists for merge/split using whole-component reconstruction.
10. README/docs clearly state that 3-D topology change is still experimental.
11. No physical film-drainage model is falsely implied to exist.
12. Existing v0.2 tests remain green.

## Suggested implementation order

1. add `MultiFrontState` and compatibility layer
2. add topology handler types
3. add 2-D detection
4. add 2-D patch extraction
5. add 2-D rasterization + marching squares reconstruction
6. add 2-D replacement + field transfer
7. add 2-D end-to-end tests/examples
8. add 3-D whole-component reconstruction prototype
9. add 3-D smoke tests/examples
10. update docs/README

## Style requirements

- keep the public API additive
- prefer readable code over generic abstraction
- keep topology logic isolated in dedicated files
- keep 2-D robust before trying to make 3-D fancy
- be explicit and honest about limitations

## Strategic note

Film drainage and rupture physics are intentionally out of v0.3 scope.
v0.3 should focus on a robust geometry-driven topology-change engine first; 
sub-grid coalescence/rupture physics can be added in later releases.
