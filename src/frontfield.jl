# frontfield.jl – Front-attached data containers.
#
# A FrontField stores a scalar or vector field defined at vertices, edges,
# or faces of a CurveMesh or SurfaceMesh.  In v0.1 only vertex fields are
# fully supported.

"""
    FrontField{V,M,L}

A field attached to a front mesh.

Fields
------
- `values   :: V` – data array; typically `Vector{T}` or `Vector{SVector{d,T}}`.
- `mesh     :: M` – the `CurveMesh` or `SurfaceMesh` this field lives on.
- `location :: Symbol` – `:vertex`, `:edge`, or `:face`.

Only `:vertex` fields are fully supported in v0.1.

Constructors
------------
    FrontField(values, mesh, location=:vertex)

Checks that `length(values)` is consistent with the location:
- `:vertex` → `length(mesh.points)` (or `mesh.points` for `CurveMesh`/`SurfaceMesh`)
- `:edge`   → number of edges
- `:face`   → number of faces (`length(mesh.faces)`)
"""
struct FrontField{V,M,L}
    values   :: V
    mesh     :: M
    location :: Symbol

    function FrontField{V,M,L}(values::V, mesh::M, location::Symbol) where {V,M,L}
        _check_frontfield_size(values, mesh, location)
        return new{V,M,L}(values, mesh, location)
    end
end

function FrontField(values::V, mesh::M, location::Symbol=:vertex) where {V,M}
    L = location
    return FrontField{V,M,typeof(L)}(values, mesh, location)
end

# ── size consistency check ────────────────────────────────────────────────────

function _frontfield_expected_length(mesh, location::Symbol)
    if location === :vertex
        return length(mesh.points)
    elseif location === :edge
        return _n_edges(mesh)
    elseif location === :face
        return _n_faces(mesh)
    else
        error("FrontField: unsupported location $(repr(location)). Use :vertex, :edge, or :face.")
    end
end

_n_edges(mesh::CurveMesh)    = length(mesh.edges)
_n_edges(mesh::SurfaceMesh)  = length(build_topology(mesh).edges)
_n_faces(mesh::CurveMesh)    = length(mesh.edges)   # edges serve as "faces" for curves
_n_faces(mesh::SurfaceMesh)  = length(mesh.faces)

function _check_frontfield_size(values, mesh, location::Symbol)
    n_expected = _frontfield_expected_length(mesh, location)
    n_actual   = length(values)
    n_actual == n_expected ||
        error("FrontField: size mismatch for location=$(repr(location)). " *
              "Expected $n_expected values, got $n_actual.")
    return nothing
end

# ── accessors ─────────────────────────────────────────────────────────────────

"""Return the mesh of a `FrontField`."""
mesh(f::FrontField)     = f.mesh

"""Return the value array of a `FrontField`."""
Base.values(f::FrontField) = f.values

"""Return the location tag (`:vertex`, `:edge`, or `:face`) of a `FrontField`."""
location(f::FrontField) = f.location

Base.length(f::FrontField) = length(f.values)
Base.getindex(f::FrontField, i) = f.values[i]
Base.setindex!(f::FrontField, v, i) = (f.values[i] = v)

function Base.show(io::IO, f::FrontField)
    T_vals = eltype(f.values)
    print(io, "FrontField{$T_vals} on $(typeof(f.mesh)) at $(f.location) " *
              "($(length(f.values)) values)")
end
