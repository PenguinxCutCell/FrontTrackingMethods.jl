# frontequation.jl – FrontEquation object and integrate! driver.
#
# FrontEquation is the top-level user-facing object, analogous to
# LevelSetMethods.LevelSetEquation.  It stores the motion terms, the time
# integrator, the front state, and optional redistribution / field-transfer
# strategies.
#
# Usage:
#   eq = FrontEquation(; terms=(AdvectionTerm(u),), front=mesh, t=0.0)
#   integrate!(eq, 1.0)

# ─────────────────────────────────────────────────────────────────────────────
# FrontEquation
# ─────────────────────────────────────────────────────────────────────────────

"""
    FrontEquation

Top-level object for a front-tracking simulation.

Fields
------
- `terms        :: Tuple`  – tuple of `AbstractFrontTerm` objects.
- `integrator   :: TimeIntegrator` – time integration scheme.
- `state        :: FrontState`    – current state (mesh, geometry, fields).
- `redistributor`                 – `AbstractRedistributor` or `nothing`.
- `transfer`                      – `AbstractFieldTransfer` or `:auto` or `nothing`.
- `step         :: Int`           – current step count.
- `buffers      :: NamedTuple`    – preallocated velocity buffers.

Constructors
------------
    FrontEquation(; terms, front, t=0.0, integrator=RK2(),
                    redistributor=nothing, transfer=nothing,
                    fields=nothing, build_dec=false)

where `front` is a `CurveMesh`, `SurfaceMesh`, `PointFront1D`, or an existing `FrontState`.
`build_dec=true` is required for `CurvatureMotionTerm` on surfaces.
"""
mutable struct FrontEquation
    terms         :: Tuple
    integrator    :: TimeIntegrator
    state         :: FrontState
    redistributor
    transfer
    step          :: Int
    buffers       :: NamedTuple
end

function FrontEquation(;
    terms,
    front,
    t        :: Real                = 0.0,
    integrator :: TimeIntegrator    = RK2(),
    redistributor                   = nothing,
    transfer                        = nothing,
    fields                          = nothing,
    build_dec :: Bool               = false,
)
    # Normalize terms to a tuple
    if terms isa AbstractFrontTerm
        terms = (terms,)
    else
        terms = Tuple(terms)
    end
    all(isa(t, AbstractFrontTerm) for t in terms) ||
        error("FrontEquation: all entries in `terms` must be AbstractFrontTerm.")

    # Build state
    state = if front isa FrontState
        front.t = Float64(t)
        front
    elseif is_supported_front(front)
        st = FrontState(front; t=t, build_dec=build_dec)
        # attach user-provided fields
        if fields !== nothing
            for (name, fld) in pairs(fields)
                add_field!(st, Symbol(name), fld)
            end
        end
        st
    else
        error("FrontEquation: `front` must be a CurveMesh, SurfaceMesh, PointFront1D, or FrontState.")
    end

    # Preallocate velocity buffer
    vel_buf = _zeros_like_points(state.mesh)
    buffers = (vel = vel_buf,)

    return FrontEquation(terms, integrator, state, redistributor, transfer, 0, buffers)
end

# ── Accessors ─────────────────────────────────────────────────────────────────

"""
    current_state(eq::FrontEquation) -> FrontState

Return the current state of the front equation.
"""
current_state(eq::FrontEquation) = eq.state

"""
    current_time(eq::FrontEquation) -> Float64

Return the current simulation time.
"""
current_time(eq::FrontEquation) = eq.state.t

# ── show ──────────────────────────────────────────────────────────────────────

function Base.show(io::IO, eq::FrontEquation)
    print(io, "FrontEquation:\n")
    print(io, "  terms: ")
    join(io, (repr(t) for t in eq.terms), " + ")
    print(io, "\n  integrator: $(eq.integrator)")
    print(io, "\n  state: $(eq.state)")
    print(io, "\n  step: $(eq.step)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Stage helpers (ForwardEuler / RK2 / RK3)
# ─────────────────────────────────────────────────────────────────────────────

# Return a fresh velocity buffer of the right size/type for `state`.
function _get_vel_buf(eq::FrontEquation, state::FrontState)
    buf = eq.buffers.vel
    if length(buf) == nmarkers(state.mesh)
        return buf
    else
        # Re-allocate if front size changed (edge case; normally fixed)
        return _zeros_like_points(state.mesh)
    end
end

# Advance the state's vertex coordinates by dt * V, refresh geometry.
function _apply_velocity!(state::FrontState, x0, V, dt::Real)
    new_pts = _advance_coordinates(state.mesh, x0, V, dt)
    state.mesh = _replace_mesh(state.mesh, new_pts)
    refresh_geometry!(state)
    return state
end

function _step_fe!(eq::FrontEquation, state::FrontState, t::Float64, dt::Float64)
    V  = _get_vel_buf(eq, state)
    x0 = vertex_coordinates(state)
    compute_rhs!(V, eq.terms, state, t)
    _apply_velocity!(state, x0, V, dt)
end

function _step_rk2!(eq::FrontEquation, state::FrontState, t::Float64, dt::Float64)
    V   = _get_vel_buf(eq, state)
    x0  = vertex_coordinates(state)

    # Stage 1: k1 = f(x0, t)
    compute_rhs!(V, eq.terms, state, t)
    k1 = copy(V)

    # Intermediate: x1 = x0 + dt * k1
    _apply_velocity!(state, x0, k1, dt)

    # Stage 2: k2 = f(x1, t + dt)
    compute_rhs!(V, eq.terms, state, t + dt)
    k2 = copy(V)

    # Final: x_new = x0 + (dt/2) * (k1 + k2)
    dt_half = dt / 2
    new_pts = [x0[a] + dt_half * (k1[a] + k2[a]) for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts)
    refresh_geometry!(state)
end

function _step_rk3!(eq::FrontEquation, state::FrontState, t::Float64, dt::Float64)
    V  = _get_vel_buf(eq, state)
    x0 = vertex_coordinates(state)

    # Stage 1: k1 = f(x0, t)
    compute_rhs!(V, eq.terms, state, t)
    k1 = copy(V)

    # x1 = x0 + dt * k1
    _apply_velocity!(state, x0, k1, dt)
    x1 = vertex_coordinates(state)

    # Stage 2: k2 = f(x1, t + dt)
    compute_rhs!(V, eq.terms, state, t + dt)
    k2 = copy(V)

    # x2 = (3/4)*x0 + (1/4)*x1 + (1/4)*dt*k2
    new_pts2 = [0.75 * x0[a] + 0.25 * x1[a] + 0.25 * dt * k2[a]
                for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts2)
    refresh_geometry!(state)

    # Stage 3: k3 = f(x2, t + dt/2)
    compute_rhs!(V, eq.terms, state, t + dt/2)
    k3 = copy(V)

    # x_new = (1/3)*x0 + (2/3)*x2 + (2/3)*dt*k3
    x2 = vertex_coordinates(state)
    new_pts = [x0[a] / 3 + (2 / 3) * x2[a] + (2 / 3) * dt * k3[a]
               for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts)
    refresh_geometry!(state)
end

# ─────────────────────────────────────────────────────────────────────────────
# integrate!
# ─────────────────────────────────────────────────────────────────────────────

"""
    integrate!(eq::FrontEquation, tf;
               dt=nothing, callback=nothing, max_steps=10^6) -> FrontEquation

Evolve the front equation from `current_time(eq)` to `tf`.

Keyword arguments
-----------------
- `dt`         – fixed time step; if `nothing` (default), use CFL-based dt.
- `callback`   – function `f(state, t, step)` called after each accepted step.
- `max_steps`  – maximum number of steps before stopping (safety limit).

After each step:
1. Geometry is refreshed (done inside the integrator stage).
2. If `eq.redistributor !== nothing`, `redistribute!(state, redistributor)` is called.
3. If a `callback` is provided, it is called with the current state, time, and step count.

Returns `eq` (mutated in-place).
"""
function integrate!(
    eq         :: FrontEquation,
    tf         :: Real;
    dt         = nothing,
    callback              = nothing,
    max_steps  :: Int     = 10^6,
)
    state = eq.state
    t     = state.t
    tf_f  = Float64(tf)

    t >= tf_f && return eq

    step = eq.step

    while t < tf_f && step < max_steps
        # ── time-step selection ──────────────────────────────────────────────
        dt_use = if dt !== nothing
            Float64(dt)
        else
            dt_cfl = compute_cfl(eq.terms, state, t)
            dt_cfl * cfl(eq.integrator)
        end
        # Clamp to not overshoot tf
        dt_use = min(dt_use, tf_f - t)
        dt_use > 0 || break

        # ── take one step ────────────────────────────────────────────────────
        integ = eq.integrator
        if integ isa ForwardEuler
            _step_fe!(eq, state, t, dt_use)
        elseif integ isa RK2
            _step_rk2!(eq, state, t, dt_use)
        elseif integ isa RK3
            _step_rk3!(eq, state, t, dt_use)
        else
            error("integrate!: unsupported integrator $(typeof(integ))")
        end

        t      += dt_use
        state.t = t
        step   += 1

        # ── optional redistribution ─────────────────────────────────────────
        if eq.redistributor !== nothing
            redistribute!(state, eq.redistributor)
        end

        # ── callback ────────────────────────────────────────────────────────
        callback !== nothing && callback(state, t, step)
    end

    eq.step = step
    return eq
end
