# topology_state.jl – Multi-component front state for topology-aware workflows.

"""
    FrontComponentState

State container for one connected front component.
"""
mutable struct FrontComponentState{M,G,D}
    mesh::M
    geom::G
    dec::D
    fields::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end

"""
    MultiFrontState <: AbstractFrontState

Topology-aware front state containing one or more connected components.
`MultiFrontState` owns simulation time.
"""
mutable struct MultiFrontState <: AbstractFrontState
    components::Vector{FrontComponentState}
    t::Float64
    global_fields::Dict{Symbol,Any}
    cache::Dict{Symbol,Any}
end

function FrontComponentState(mesh::CurveMesh; fields=Dict{Symbol,Any}(), build_dec::Bool=false)
    geom = compute_geometry(mesh)
    dec = build_dec ? FrontIntrinsicOps.build_dec(mesh, geom) : nothing
    return FrontComponentState{typeof(mesh),typeof(geom),typeof(dec)}(
        mesh,
        geom,
        dec,
        Dict{Symbol,Any}(fields),
        Dict{Symbol,Any}(),
    )
end

function FrontComponentState(mesh::SurfaceMesh; fields=Dict{Symbol,Any}(), build_dec::Bool=false)
    geom = compute_geometry(mesh)
    if build_dec
        dec = FrontIntrinsicOps.build_dec(mesh, geom)
        geom = FrontIntrinsicOps.compute_curvature(mesh, geom, dec)
    else
        dec = nothing
    end
    return FrontComponentState{typeof(mesh),typeof(geom),typeof(dec)}(
        mesh,
        geom,
        dec,
        Dict{Symbol,Any}(fields),
        Dict{Symbol,Any}(),
    )
end

function MultiFrontState(components::AbstractVector{<:FrontComponentState};
                         t::Real=0.0,
                         global_fields=Dict{Symbol,Any}(),
                         cache=Dict{Symbol,Any}())
    return MultiFrontState(
    FrontComponentState[components...],
        Float64(t),
        Dict{Symbol,Any}(global_fields),
        Dict{Symbol,Any}(cache),
    )
end

function MultiFrontState(fronts::AbstractVector{<:Union{CurveMesh,SurfaceMesh}};
                         t::Real=0.0,
                         build_dec::Bool=false)
    components = FrontComponentState[]
    sizehint!(components, length(fronts))
    for mesh in fronts
        push!(components, FrontComponentState(mesh; build_dec=build_dec))
    end
    return MultiFrontState(components; t=t)
end

MultiFrontState(state::FrontState) = MultiFrontState([
    FrontComponentState(state.mesh, state.geom, state.dec,
                        Dict{Symbol,Any}(state.fields), Dict{Symbol,Any}(state.cache))
]; t=state.t)

function FrontState(multistate::MultiFrontState; build_dec::Bool=false)
    ncomponents(multistate) == 1 ||
        error("FrontState(multistate): conversion requires exactly one component.")
    comp = component(multistate, 1)
    st = if comp.mesh isa CurveMesh
        FrontState(comp.mesh; t=multistate.t, fields=comp.fields, build_dec=build_dec || comp.dec !== nothing)
    elseif comp.mesh isa SurfaceMesh
        FrontState(comp.mesh; t=multistate.t, fields=comp.fields, build_dec=build_dec || comp.dec !== nothing)
    else
        error("FrontState(multistate): unsupported mesh type $(typeof(comp.mesh)).")
    end
    return st
end

current_state(state::MultiFrontState) = state
current_time(state::MultiFrontState) = state.t
ncomponents(state::MultiFrontState) = length(state.components)
component(state::MultiFrontState, i::Integer) = state.components[i]
component_mesh(state::MultiFrontState, i::Integer) = state.components[i].mesh
component_geom(state::MultiFrontState, i::Integer) = state.components[i].geom
component_fields(state::MultiFrontState, i::Integer) = state.components[i].fields
all_meshes(state::MultiFrontState) = [comp.mesh for comp in state.components]
vertex_coordinates(state::MultiFrontState) = [copy(comp.mesh.points) for comp in eachcomponent(state)]

eachcomponent(state::MultiFrontState) = state.components
map_components(f, state::MultiFrontState) = map(f, state.components)

function add_field!(state::MultiFrontState, name::Symbol, value)
    state.global_fields[name] = value
    return state
end

function get_field(state::MultiFrontState, name::Symbol)
    haskey(state.global_fields, name) ||
        error("MultiFrontState: no global field named $(repr(name)).")
    return state.global_fields[name]
end

function Base.show(io::IO, state::MultiFrontState)
    print(io, "MultiFrontState: ", ncomponents(state), " components, t = ", state.t)
end

# Internal helper used by equations/terms to reuse existing FrontState machinery.
_component_front_state(comp::FrontComponentState, t::Real) =
    FrontState{typeof(comp.mesh),typeof(comp.geom),typeof(comp.dec)}(
        comp.mesh,
        comp.geom,
        comp.dec,
        Float64(t),
        comp.fields,
        comp.cache,
    )

_sync_component_from_state!(comp::FrontComponentState, state::FrontState) = begin
    comp.mesh = state.mesh
    comp.geom = state.geom
    comp.dec = state.dec
    comp.fields = state.fields
    comp.cache = state.cache
    comp
end
