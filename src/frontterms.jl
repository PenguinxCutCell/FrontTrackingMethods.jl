# frontterms.jl – Motion terms for front evolution.
#
# Three physical motion terms mirror the LevelSetMethods API:
#
#   AdvectionTerm        – prescribed vector velocity  u(x, t)
#   NormalMotionTerm     – prescribed scalar normal speed  Vn(x, t)
#   CurvatureMotionTerm  – coefficient * mean-curvature-normal vector
#
# General vertex evolution:
#   dx_a/dt = u_adv(x_a, t) + Vn(x_a, t) * n_a + β(x_a, t) * Hn_a
#
# Sign conventions (from FrontIntrinsicOps)
# ------------------------------------------
# CurveMesh:
#   `geom.vertex_normals[a]` – unit INWARD normal for CCW curves.
#   `geom.signed_curvature[a]` – positive for convex CCW curves.
#   CurvatureMotionTerm velocity = β * κ * n_inward.
#   → for β > 0 a CCW circle SHRINKS (curve-shortening flow).
#   Exact radius evolution: R(t) = sqrt(R₀² - 2 β t).
#
# SurfaceMesh:
#   `geom.vertex_normals[a]` – unit OUTWARD normal.
#   `geom.mean_curvature_normal[a]` – INWARD pointing for convex outward-
#     oriented surface; ‖Hn‖ → 1/R for a sphere of radius R.
#   CurvatureMotionTerm velocity = β * Hn.
#   → for β > 0 a sphere SHRINKS (mean-curvature flow).
#   Exact radius evolution: R(t) = sqrt(R₀² - 2 β t).
#
# For pure normal motion (NormalMotionTerm):
#   `geom.vertex_normals[a]` is used as the normal direction.
#   For CCW curves: inward normal → Vn > 0 shrinks, Vn < 0 expands.
#   For outward-oriented surfaces: outward normal → Vn > 0 expands, Vn < 0 shrinks.

# ─────────────────────────────────────────────────────────────────────────────
# Helper: evaluate a velocity/speed/coefficient that may be a function,
# a FrontField, a constant vector/scalar, or a Number.
# ─────────────────────────────────────────────────────────────────────────────

@inline function _eval_velocity(vel::Function, x, t, state)
    return vel(x, t, state)
end
@inline function _eval_velocity(vel::FrontField, x, t, state)
    # Find vertex index by position lookup – used only for function-like API
    # For now, accept a FrontField and index externally; see accumulate_term!
    error("_eval_velocity: use FrontField via accumulate_term! with index.")
end
@inline function _eval_velocity(vel, x, t, state)
    return vel  # constant vector
end

@inline function _eval_speed(sp::Function, x, t, state)
    return sp(x, t, state)
end
@inline function _eval_speed(sp::FrontField, x, t, state)
    error("_eval_speed: use FrontField via accumulate_term! with index.")
end
@inline function _eval_speed(sp::Number, x, t, state)
    return sp
end
@inline function _eval_speed(sp, x, t, state)
    return sp
end

@inline function _eval_coeff(coeff::Function, x, t, state)
    return coeff(x, t, state)
end
@inline function _eval_coeff(coeff::Number, x, t, state)
    return coeff
end
@inline function _eval_coeff(coeff, x, t, state)
    return coeff
end

# ─────────────────────────────────────────────────────────────────────────────
# AdvectionTerm
# ─────────────────────────────────────────────────────────────────────────────

"""
    AdvectionTerm(velocity[, update_func])

Prescribed vector advection: each vertex moves as dx/dt = u(x, t).

`velocity` may be:
- A function `(x, t, state) -> SVector` evaluated at each vertex.
- A `FrontField` of vectors at vertices.
- A constant vector (same for all vertices).

`update_func`, if provided, is called as `update_func(velocity, state, t)`
at the beginning of each RHS evaluation stage; it can modify `velocity`
in-place.
"""
struct AdvectionTerm{V,F} <: AbstractFrontTerm
    velocity    :: V
    update_func :: F
end

AdvectionTerm(vel) = AdvectionTerm(vel, (v, s, t) -> nothing)
AdvectionTerm(vel, f) = AdvectionTerm{typeof(vel),typeof(f)}(vel, f)

velocity(term::AdvectionTerm)     = term.velocity
update_func(term::AdvectionTerm)  = term.update_func

Base.show(io::IO, ::AdvectionTerm) = print(io, "AdvectionTerm (dx/dt = u(x,t))")

function update_term!(term::AdvectionTerm, state, t)
    term.update_func(term.velocity, state, t)
    return nothing
end

function accumulate_term!(V, term::AdvectionTerm{<:Function}, state, t)
    pts = state.mesh.points
    for a in eachindex(pts)
        V[a] = V[a] + term.velocity(pts[a], t, state)
    end
    return V
end

function accumulate_term!(V, term::AdvectionTerm{<:FrontField}, state, t)
    fv = term.velocity.values
    for a in eachindex(V)
        V[a] = V[a] + fv[a]
    end
    return V
end

function accumulate_term!(V, term::AdvectionTerm, state, t)
    # Constant velocity
    vel = term.velocity
    for a in eachindex(V)
        V[a] = V[a] + vel
    end
    return V
end

function compute_cfl(term::AdvectionTerm{<:Function}, state, t)
    h = front_spacing(state)
    pts = state.mesh.points
    max_speed = maximum(norm(term.velocity(pts[a], t, state)) for a in eachindex(pts))
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

function compute_cfl(term::AdvectionTerm{<:FrontField}, state, t)
    h = front_spacing(state)
    fv = term.velocity.values
    max_speed = maximum(norm(fv[a]) for a in eachindex(fv))
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

function compute_cfl(term::AdvectionTerm, state, t)
    h = front_spacing(state)
    max_speed = norm(term.velocity)
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

# ─────────────────────────────────────────────────────────────────────────────
# NormalMotionTerm
# ─────────────────────────────────────────────────────────────────────────────

"""
    NormalMotionTerm(speed[, update_func])

Prescribed normal speed: dx/dt = Vn(x, t) * n_a, where n_a is the vertex
normal from `FrontIntrinsicOps`.

Sign convention:
- CurveMesh: n_a is the INWARD normal for CCW curves. Vn > 0 shrinks; Vn < 0 expands.
- SurfaceMesh: n_a is the OUTWARD normal. Vn > 0 expands; Vn < 0 shrinks.

`speed` may be:
- A function `(x, t, state) -> Real`.
- A scalar `FrontField`.
- A constant scalar.

`update_func`, if provided, is called as `update_func(speed, state, t)`.
"""
struct NormalMotionTerm{V,F} <: AbstractFrontTerm
    speed       :: V
    update_func :: F
end

NormalMotionTerm(sp) = NormalMotionTerm(sp, (v, s, t) -> nothing)
NormalMotionTerm(sp, f) = NormalMotionTerm{typeof(sp),typeof(f)}(sp, f)

speed(term::NormalMotionTerm)       = term.speed
update_func(term::NormalMotionTerm) = term.update_func

Base.show(io::IO, ::NormalMotionTerm) = print(io, "NormalMotionTerm (dx/dt = Vn * n)")

function update_term!(term::NormalMotionTerm, state, t)
    term.update_func(term.speed, state, t)
    return nothing
end

function accumulate_term!(V, term::NormalMotionTerm{<:Function}, state, t)
    pts     = state.mesh.points
    normals = state.geom.vertex_normals
    for a in eachindex(pts)
        vn = term.speed(pts[a], t, state)
        V[a] = V[a] + vn * normals[a]
    end
    return V
end

function accumulate_term!(V, term::NormalMotionTerm{<:FrontField}, state, t)
    normals = state.geom.vertex_normals
    sv = term.speed.values
    for a in eachindex(V)
        V[a] = V[a] + sv[a] * normals[a]
    end
    return V
end

function accumulate_term!(V, term::NormalMotionTerm, state, t)
    normals = state.geom.vertex_normals
    vn = term.speed
    for a in eachindex(V)
        V[a] = V[a] + vn * normals[a]
    end
    return V
end

function compute_cfl(term::NormalMotionTerm{<:Function}, state, t)
    h = front_spacing(state)
    pts = state.mesh.points
    max_speed = maximum(abs(term.speed(pts[a], t, state)) for a in eachindex(pts))
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

function compute_cfl(term::NormalMotionTerm{<:FrontField}, state, t)
    h = front_spacing(state)
    sv = term.speed.values
    max_speed = maximum(abs(sv[a]) for a in eachindex(sv))
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

function compute_cfl(term::NormalMotionTerm, state, t)
    h = front_spacing(state)
    max_speed = abs(term.speed)
    max_speed < eps(typeof(h)) && return Inf
    return h / max_speed
end

# ─────────────────────────────────────────────────────────────────────────────
# CurvatureMotionTerm
# ─────────────────────────────────────────────────────────────────────────────

"""
    CurvatureMotionTerm(coefficient[, update_func])

Curvature-driven motion: dx/dt = β(x, t) * Hn_a, where Hn_a is the
mean-curvature-normal vector from `FrontIntrinsicOps`.

Sign conventions (from FrontIntrinsicOps)
------------------------------------------
- CurveMesh: Hn_a = κ_a * n_a, where κ_a is the discrete signed curvature
  (positive for CCW curves) and n_a is the inward unit normal.
  For β > 0, a CCW circle shrinks (curve-shortening flow).
  Exact: R(t) = sqrt(R₀² - 2 β t).

- SurfaceMesh: Hn_a is the discrete mean-curvature-normal vector; it points
  INWARD for convex outward-oriented surfaces, with ‖Hn‖ → 1/R.
  For β > 0, a sphere shrinks (mean-curvature flow).
  Exact: R(t) = sqrt(R₀² - 2 β t).

`coefficient` may be:
- A function `(x, t, state) -> Real`.
- A constant scalar.

For surface curvature, `FrontState` must be built with `build_dec=true`.
"""
struct CurvatureMotionTerm{C,F} <: AbstractFrontTerm
    coefficient :: C
    update_func :: F
end

CurvatureMotionTerm(coeff) = CurvatureMotionTerm(coeff, (c, s, t) -> nothing)
CurvatureMotionTerm(coeff, f) = CurvatureMotionTerm{typeof(coeff),typeof(f)}(coeff, f)

coefficient(term::CurvatureMotionTerm)  = term.coefficient
update_func(term::CurvatureMotionTerm)  = term.update_func

Base.show(io::IO, ::CurvatureMotionTerm) = print(io, "CurvatureMotionTerm (dx/dt = β Hn)")

function update_term!(term::CurvatureMotionTerm, state, t)
    term.update_func(term.coefficient, state, t)
    return nothing
end

function accumulate_term!(V, term::CurvatureMotionTerm, state, t)
    mesh = state.mesh
    if mesh isa CurveMesh
        kappa   = state.geom.signed_curvature
        normals = state.geom.vertex_normals
        for a in eachindex(mesh.points)
            beta  = _eval_coeff(term.coefficient, mesh.points[a], t, state)
            V[a]  = V[a] + beta * kappa[a] * normals[a]
        end
    elseif mesh isa SurfaceMesh
        Hn = state.geom.mean_curvature_normal
        isempty(Hn) &&
            error("CurvatureMotionTerm on SurfaceMesh requires `build_dec=true` " *
                  "when constructing FrontState (or FrontEquation).")
        for a in eachindex(mesh.points)
            beta = _eval_coeff(term.coefficient, mesh.points[a], t, state)
            V[a] = V[a] + beta * Hn[a]
        end
    else
        error("CurvatureMotionTerm: unsupported mesh type $(typeof(mesh))")
    end
    return V
end

function compute_cfl(term::CurvatureMotionTerm, state, t)
    h    = front_spacing(state)
    mesh = state.mesh
    pts  = mesh.points
    if mesh isa CurveMesh
        max_abs_beta = maximum(abs(_eval_coeff(term.coefficient, pts[a], t, state))
                               for a in eachindex(pts))
    else
        max_abs_beta = maximum(abs(_eval_coeff(term.coefficient, pts[a], t, state))
                               for a in eachindex(pts))
    end
    max_abs_beta < eps(typeof(h)) && return Inf
    return h^2 / (2 * max_abs_beta)
end

# ─────────────────────────────────────────────────────────────────────────────
# Aggregate CFL over a tuple of terms
# ─────────────────────────────────────────────────────────────────────────────

"""
    compute_cfl(terms, state, t) -> Float64

Compute the minimum stable time-step from all terms.
Returns `Inf` if all terms have zero effect.
"""
function compute_cfl(terms::Tuple, state, t)
    dt = minimum(compute_cfl(term, state, t) for term in terms)
    isfinite(dt) && dt > 0 || error("compute_cfl: got invalid dt = $dt")
    return dt
end

# ─────────────────────────────────────────────────────────────────────────────
# Generic update_term! / accumulate_term! fallback
# ─────────────────────────────────────────────────────────────────────────────

"""
    update_term!(term, state, t)

Prepare a term before RHS evaluation (e.g. update a time-dependent field).
Default is a no-op.
"""
update_term!(::AbstractFrontTerm, state, t) = nothing

"""
    compute_rhs!(V, terms, state, t) -> V

Assemble the total vertex velocity `V` as the sum of all term contributions.
`V` is a `Vector` of `SVector`s, preallocated to match the number of vertices.
"""
function compute_rhs!(V, terms, state, t)
    fill!(V, zero(eltype(V)))
    for term in terms
        update_term!(term, state, t)
        accumulate_term!(V, term, state, t)
    end
    return V
end
