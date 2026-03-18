# topology_replace.jl – Replace state components after topology reconstruction.

function _to_component_state(obj)
    if obj isa FrontComponentState
        return obj
    elseif obj isa CurveMesh
        return FrontComponentState(obj)
    elseif obj isa SurfaceMesh
        return FrontComponentState(obj)
    else
        error("replace_components!: unsupported component type $(typeof(obj)).")
    end
end

function replace_components!(state::MultiFrontState, old_ids, new_components)
    isempty(old_ids) && return state
    old_set = Set(Int.(old_ids))

    old_components = [state.components[i] for i in sort(collect(old_set))]

    kept = FrontComponentState[]
    for (i, comp) in enumerate(state.components)
        if !(i in old_set)
            push!(kept, comp)
        end
    end

    insert_pos = minimum(old_set)
    insert_pos = clamp(insert_pos, 1, length(kept) + 1)

    new_states = [_to_component_state(c) for c in new_components]

    merged = FrontComponentState[]
    append!(merged, kept[1:insert_pos-1])
    append!(merged, new_states)
    append!(merged, kept[insert_pos:end])

    state.components = merged
    refresh_geometry!(state; which=:all)

    state.cache[:last_replaced_old_components] = old_components
    state.cache[:last_replaced_new_components] = new_states

    return state
end

function replace_curve_patch!(component::FrontComponentState,
                              patch::EventPatch,
                              reconstructed_curves,
                              anchors)
    isempty(reconstructed_curves) && return component
    # v0.3 first pass: use first reconstructed curve as full-component replacement.
    component.mesh = reconstructed_curves[1]
    _refresh_component_geometry!(component)
    component.cache[:patch_replace] = Dict(:patch => patch, :anchors => anchors)
    return component
end
