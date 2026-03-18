# topology_transfer.jl – Transfer scalar vertex fields across topology changes.

function _nearest_scalar_from_old_components(p, old_components::AbstractVector{<:FrontComponentState}, field_name::Symbol)
    best_dist = Inf
    best_val = nothing
    for comp in old_components
        haskey(comp.fields, field_name) || continue
        fld = comp.fields[field_name]
        fld isa FrontField || continue
        fld.location === :vertex || continue
        vals = fld.values
        eltype(vals) <: Number || error("transfer_fields_after_topology_change!: only scalar vertex fields are supported.")

        for (i, q) in enumerate(comp.mesh.points)
            d = norm(p - q)
            if d < best_dist
                best_dist = d
                best_val = vals[i]
            end
        end
    end
    return best_val
end

function transfer_fields_after_topology_change!(old_components::AbstractVector{<:FrontComponentState},
                                                new_components::AbstractVector{<:FrontComponentState};
                                                method::Symbol=:auto)
    if isempty(old_components) || isempty(new_components)
        return new_components
    end

    field_names = Set{Symbol}()
    for comp in old_components
        for (name, fld) in comp.fields
            if fld isa FrontField && fld.location === :vertex
                push!(field_names, name)
            end
        end
    end

    for new_comp in new_components
        for name in field_names
            vals = Vector{Float64}(undef, length(new_comp.mesh.points))
            for (i, p) in enumerate(new_comp.mesh.points)
                v = _nearest_scalar_from_old_components(p, old_components, name)
                if v === nothing
                    vals[i] = 0.0
                else
                    vals[i] = Float64(v)
                end
            end
            new_comp.fields[name] = FrontField(vals, new_comp.mesh, :vertex)
        end
    end

    return new_components
end
