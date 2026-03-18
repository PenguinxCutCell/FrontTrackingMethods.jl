# topology_change.jl – Topology-change handler API and integration hooks.

abstract type AbstractTopologyHandler end

struct NoTopologyChange <: AbstractTopologyHandler end

struct LocalCartesianTopologyHandler <: AbstractTopologyHandler
    d_trigger::Float64
    d_merge::Float64
    d_split::Float64
    patch_h_factor::Float64
    patch_margin_factor::Float64
    reconstruct_scope::Symbol
    use_indicator::Bool
    preserve_fields::Bool
    max_components_created::Int
end

function LocalCartesianTopologyHandler(;
    d_trigger::Real=1.0,
    d_merge::Real=1.2,
    d_split::Real=1.0,
    patch_h_factor::Real=0.75,
    patch_margin_factor::Real=2.0,
    reconstruct_scope::Symbol=:local_patch,
    use_indicator::Bool=true,
    preserve_fields::Bool=true,
    max_components_created::Int=4,
)
    return LocalCartesianTopologyHandler(
        Float64(d_trigger),
        Float64(d_merge),
        Float64(d_split),
        Float64(patch_h_factor),
        Float64(patch_margin_factor),
        reconstruct_scope,
        use_indicator,
        preserve_fields,
        max_components_created,
    )
end

struct TopologyEventReport
    changed::Bool
    event_type::Symbol
    component_ids::Vector{Int}
    patch_bbox
    n_old_components::Int
    n_new_components::Int
    notes::String
end

struct TopologyCandidate
    event_type::Symbol
    component_ids::Vector{Int}
    feature_pairs::Any
    distance::Float64
    h_local::Float64
    score::Float64
    bbox
end

function _none_topology_report(state::AbstractFrontState; notes::String="")
    n = state isa MultiFrontState ? ncomponents(state) : 1
    return TopologyEventReport(false, :none, Int[], nothing, n, n, notes)
end

handle_topology_change!(state::FrontState, ::NoTopologyChange) = _none_topology_report(state)
handle_topology_change!(state::MultiFrontState, ::NoTopologyChange) = _none_topology_report(state)

function handle_topology_change!(state::FrontState, handler::LocalCartesianTopologyHandler)
    return TopologyEventReport(
        false,
        :ambiguous,
        [1],
        nothing,
        1,
        1,
        "Topology-change reconstruction requires MultiFrontState. Construct FrontEquation with front=[mesh] or state=MultiFrontState(...).",
    )
end

function handle_topology_change!(state::MultiFrontState, handler::LocalCartesianTopologyHandler)
    candidates = find_topology_candidates(state, handler)
    selected = select_topology_candidate(candidates)

    if selected === nothing
        return _none_topology_report(state)
    end

    patch = extract_event_patch(state, selected, handler)
    D = length(patch.patch_grid.dims)

    reconstructed = if D == 2
        χ = zeros(Float64, patch.patch_grid.dims)
        rasterize_indicator!(χ, patch, patch.local_mesh_data)
        reconstruct_curve_from_patch(χ, patch)
    elseif D == 3
        χ = zeros(Float64, patch.patch_grid.dims)
        rasterize_indicator!(χ, patch, patch.local_mesh_data)
        reconstruct_surface_from_patch(χ, patch)
    else
        Any[]
    end

    if isempty(reconstructed)
        return TopologyEventReport(
            false,
            :ambiguous,
            selected.component_ids,
            patch.margin_bbox,
            ncomponents(state),
            ncomponents(state),
            "Reconstruction produced no components; topology update skipped.",
        )
    end

    if length(reconstructed) > handler.max_components_created
        return TopologyEventReport(
            false,
            :ambiguous,
            selected.component_ids,
            patch.margin_bbox,
            ncomponents(state),
            ncomponents(state),
            "Reconstruction produced too many components; topology update skipped.",
        )
    end

    old_ids = selected.component_ids
    old_components = [component(state, i) for i in old_ids]
    new_components = FrontComponentState[_to_component_state(m) for m in reconstructed]

    if handler.preserve_fields
        transfer_fields_after_topology_change!(old_components, new_components)
    end

    n_old = ncomponents(state)
    replace_components!(state, old_ids, new_components)

    for comp in eachcomponent(state)
        if comp.mesh isa CurveMesh
            cstate = _component_front_state(comp, state.t)
            redistribute!(cstate, CurveEqualArcRedistributor())
            _sync_component_from_state!(comp, cstate)
        end
    end
    refresh_geometry!(state; which=:all)

    n_new = ncomponents(state)

    return TopologyEventReport(
        true,
        selected.event_type,
        old_ids,
        patch.margin_bbox,
        n_old,
        n_new,
        "Topology updated via local Cartesian reconstruction.",
    )
end
