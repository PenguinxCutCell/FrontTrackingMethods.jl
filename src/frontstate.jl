# frontstate.jl – Mutable state object for front tracking.
#
# FrontState wraps a front mesh together with its current geometry,
# optional DEC operators, current time, attached field dictionary, and
# an internal cache for low-allocation stepping.
#
# Because CurveMesh and SurfaceMesh are immutable, updating vertex
# coordinates is done by replacing the mesh field with a new instance
# that shares the same connectivity.

"""
    FrontState

Mutable state object for a front-tracking simulation.

Fields
------
- `mesh   :: Union{CurveMesh,SurfaceMesh}` – current front mesh.
- `geom   :: Union{CurveGeometry,SurfaceGeometry}` – current geometry.
- `dec    :: Union{CurveDEC,SurfaceDEC,Nothing}` – optional DEC operators.
- `t      :: Float64` – current simulation time.
- `fields :: Dict{Symbol,Any}` – user-attached named fields (see `FrontField`).
- `cache  :: Dict{Symbol,Any}` – internal preallocated arrays.

Constructors
------------
    FrontState(mesh; t=0.0, fields=Dict(), build_dec=false)

Builds a `FrontState` from a `CurveMesh` or `SurfaceMesh`.
If `build_dec=true`, the DEC operators are assembled (required for
`CurvatureMotionTerm` on surfaces).
"""
mutable struct FrontState{M,G,D}
    mesh   :: M
    geom   :: G
    dec    :: D
    t      :: Float64
    fields :: Dict{Symbol,Any}
    cache  :: Dict{Symbol,Any}
end

# ── Constructor ───────────────────────────────────────────────────────────────

function FrontState(mesh::CurveMesh; t::Real=0.0, fields=Dict{Symbol,Any}(),
                    build_dec::Bool=false)
    geom = compute_geometry(mesh)
    dec  = build_dec ? FrontIntrinsicOps.build_dec(mesh, geom) : nothing
    cache = Dict{Symbol,Any}()
    return FrontState{typeof(mesh), typeof(geom), typeof(dec)}(
        mesh, geom, dec, Float64(t), Dict{Symbol,Any}(fields), cache)
end

function FrontState(mesh::SurfaceMesh; t::Real=0.0, fields=Dict{Symbol,Any}(),
                    build_dec::Bool=false)
    geom = compute_geometry(mesh)
    if build_dec
        dec  = FrontIntrinsicOps.build_dec(mesh, geom)
        geom = FrontIntrinsicOps.compute_curvature(mesh, geom, dec)
    else
        dec = nothing
    end
    cache = Dict{Symbol,Any}()
    return FrontState{typeof(mesh), typeof(geom), typeof(dec)}(
        mesh, geom, dec, Float64(t), Dict{Symbol,Any}(fields), cache)
end

# ── Accessors ─────────────────────────────────────────────────────────────────

"""
    current_state(state::FrontState) -> FrontState

Return the state itself (for API symmetry with `FrontEquation`).
"""
current_state(state::FrontState) = state

"""
    current_time(state::FrontState) -> Float64

Return the current simulation time.
"""
current_time(state::FrontState) = state.t

"""
    vertex_coordinates(state::FrontState) -> Vector{SVector}

Return a copy of the current vertex coordinate array.
"""
vertex_coordinates(state::FrontState) = copy(state.mesh.points)

"""
    set_vertex_coordinates!(state::FrontState, new_points)

Replace the vertex coordinates in `state` with `new_points`.
Creates a new mesh with the same connectivity but updated points.
Geometry is NOT automatically refreshed; call `refresh_geometry!` afterwards.
"""
function set_vertex_coordinates!(state::FrontState, new_points)
    mesh = state.mesh
    if mesh isa CurveMesh
        state.mesh = CurveMesh(convert(typeof(mesh.points), new_points), mesh.edges)
    elseif mesh isa SurfaceMesh
        state.mesh = SurfaceMesh(convert(typeof(mesh.points), new_points), mesh.faces)
    else
        error("set_vertex_coordinates!: unsupported mesh type $(typeof(mesh))")
    end
    return state
end

# ── add_field! ────────────────────────────────────────────────────────────────

"""
    add_field!(state::FrontState, name::Symbol, field::FrontField)

Attach a named field to the state.
"""
function add_field!(state::FrontState, name::Symbol, field)
    state.fields[name] = field
    return state
end

"""
    get_field(state::FrontState, name::Symbol) -> FrontField

Retrieve a named field from the state. Throws if not found.
"""
function get_field(state::FrontState, name::Symbol)
    haskey(state.fields, name) ||
        error("FrontState: no field named $(repr(name)).")
    return state.fields[name]
end

# ── show ──────────────────────────────────────────────────────────────────────

function Base.show(io::IO, state::FrontState)
    nv   = length(state.mesh.points)
    mtyp = typeof(state.mesh)
    print(io, "FrontState: $mtyp with $nv vertices, t = $(state.t)")
    if !isempty(state.fields)
        print(io, ", fields: $(collect(keys(state.fields)))")
    end
end
