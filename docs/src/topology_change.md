# Topology change (v0.3)

v0.3 introduces an additive topology-change pipeline for explicit front tracking.
The design goal is to keep the core explicit representation while enabling
robust, geometry-driven merge/split handling in 2-D and a first coarse 3-D
prototype.

## Why the state model changed

Topology events can change component count (one splits to two, two merge to one).
A single-mesh state is no longer sufficient for those transitions.

v0.3 therefore adds:

- `FrontComponentState` for one connected component,
- `MultiFrontState` for a vector of components plus global time/fields,
- compatibility helpers to/from legacy single-component `FrontState`.

This keeps existing fixed-topology workflows working while enabling topology-aware
workflows in a new state layer.

## Detection philosophy: geometry first

v0.3 topology detection is geometric and conservative:

- near-contact and self-contact proximity checks,
- neck-thickness style split heuristics,
- one-event-at-a-time candidate selection.

Importantly, v0.3 does **not** include film-drainage / rupture physics. The
current detector answers geometric imminence, not physical coalescence onset.

## Local patch vs whole-component reconstruction

The topology handler supports reconstruction scopes:

- `:local_patch` (preferred target for 2-D):
  reconstruct only the local event neighborhood.
- `:whole_component` (default practical path for 3-D prototype):
  reconstruct all involved components for simpler, more robust logic.

## Auxiliary Cartesian reconstruction

For a selected event candidate:

1. Extract an event patch (bounding box + margin).
2. Build a Cartesian patch grid from local mesh scale.
3. Rasterize an indicator field (`χ`) on the patch.
4. Reconstruct interface components from the patch indicator.
5. Replace old components in state and refresh geometry.

Current implementation emphasis:

- 2-D: robust occupancy-based reconstruction of closed loops,
- 3-D: coarse prototype reconstruction suitable for smoke/regression tests.

## Field transfer across events

v0.3 transfers attached scalar vertex fields through topology replacement using
conservative nearest-source sampling across old components.

Current limitations:

- scalar vertex fields are the supported target,
- advanced tensor/vector and high-order conservative transfer are future work.

## Current limitations and scope

- 2-D topology change is the main validated target in v0.3.
- 3-D topology change is first-generation and explicitly experimental.
- No film-drainage / rupture sub-grid model in v0.3.
- No non-manifold or triple-line topology handling.
- No contact-line/contact-angle topology physics.

Use this feature as a geometry engine foundation; higher-fidelity physical
coalescence criteria are intended for later releases.
