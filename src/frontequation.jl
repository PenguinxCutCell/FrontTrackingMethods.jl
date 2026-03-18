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

where `front` is a `CurveMesh`, `SurfaceMesh`, or an existing `FrontState`.
`build_dec=true` is required for `CurvatureMotionTerm` on surfaces.
"""
mutable struct FrontEquation
    terms         :: Tuple
    integrator    :: TimeIntegrator
    state         :: AbstractFrontState
    topology_handler
    redistributor
    transfer
    step          :: Int
    buffers       :: NamedTuple
end

function FrontEquation(;
    terms,
    front                           = nothing,
    state                           = nothing,
    t        :: Real                = 0.0,
    integrator :: TimeIntegrator    = RK2(),
    topology_handler                = NoTopologyChange(),
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
    (state === nothing || front === nothing) ||
        error("FrontEquation: provide either `front` or `state`, not both.")

    built_state = if state !== nothing
        state isa AbstractFrontState ||
            error("FrontEquation: `state` must be a FrontState or MultiFrontState.")
        state.t = Float64(t)
        state
    elseif front isa FrontState || front isa MultiFrontState
        front.t = Float64(t)
        front
    elseif front isa CurveMesh || front isa SurfaceMesh
        st = FrontState(front; t=t, build_dec=build_dec)
        if fields !== nothing
            for (name, fld) in pairs(fields)
                add_field!(st, Symbol(name), fld)
            end
        end
        st
    elseif front isa AbstractVector{<:Union{CurveMesh,SurfaceMesh}}
        mst = MultiFrontState(front; t=t, build_dec=build_dec)
        if fields !== nothing
            for (name, fld) in pairs(fields)
                mst.global_fields[Symbol(name)] = fld
            end
        end
        mst
    else
        error("FrontEquation: `front` must be mesh/mesh-vector/state or provide `state=`.")
    end

    if !(topology_handler isa NoTopologyChange) && built_state isa FrontState
        built_state = MultiFrontState(built_state)
    end

    # Preallocate velocity buffer
    buffers = if built_state isa FrontState
        (vel = _zeros_like_points(built_state.mesh), vel_components = Vector{Any}())
    else
        vc = [_zeros_like_points(comp.mesh) for comp in eachcomponent(built_state)]
        (vel = nothing, vel_components = vc)
    end

    return FrontEquation(terms, integrator, built_state, topology_handler,
                         redistributor, transfer, 0, buffers)
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
    if length(buf) == length(state.mesh.points)
        return buf
    else
        # Re-allocate if front size changed (edge case; normally fixed)
        return _zeros_like_points(state.mesh)
    end
end

function _get_vel_buf(eq::FrontEquation, state::MultiFrontState, i::Int)
    bufs = eq.buffers.vel_components
    if length(bufs) < i
        resize!(bufs, i)
        bufs[i] = _zeros_like_points(component_mesh(state, i))
        eq.buffers = (vel = eq.buffers.vel, vel_components = bufs)
    end
    buf = bufs[i]
    if length(buf) != length(component_mesh(state, i).points)
        bufs[i] = _zeros_like_points(component_mesh(state, i))
        eq.buffers = (vel = eq.buffers.vel, vel_components = bufs)
    end
    return eq.buffers.vel_components[i]
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

function _step_fe!(eq::FrontEquation, state::MultiFrontState, t::Float64, dt::Float64)
    for (i, comp) in enumerate(eachcomponent(state))
        cstate = _component_front_state(comp, t)
        V  = _get_vel_buf(eq, state, i)
        x0 = vertex_coordinates(cstate)
        compute_rhs!(V, eq.terms, cstate, t)
        _apply_velocity!(cstate, x0, V, dt)
        _sync_component_from_state!(comp, cstate)
    end
    return state
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
    T = eltype(eltype(x0))
    half_dt = T(dt) / 2
    k_avg = [half_dt * (k1[a] + k2[a]) / T(dt) * T(dt) for a in eachindex(k1)]
    # equivalently: x0 + (dt/2)*(k1+k2)
    new_pts = [x0[a] + T(dt)/2 * (k1[a] + k2[a]) for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts)
    refresh_geometry!(state)
end

function _step_rk2!(eq::FrontEquation, state::MultiFrontState, t::Float64, dt::Float64)
    for (i, comp) in enumerate(eachcomponent(state))
        cstate = _component_front_state(comp, t)
        V   = _get_vel_buf(eq, state, i)
        x0  = vertex_coordinates(cstate)

        compute_rhs!(V, eq.terms, cstate, t)
        k1 = copy(V)

        _apply_velocity!(cstate, x0, k1, dt)

        compute_rhs!(V, eq.terms, cstate, t + dt)
        k2 = copy(V)

        T = eltype(eltype(x0))
        new_pts = [x0[a] + T(dt)/2 * (k1[a] + k2[a]) for a in eachindex(x0)]
        cstate.mesh = _replace_mesh(cstate.mesh, new_pts)
        refresh_geometry!(cstate)

        _sync_component_from_state!(comp, cstate)
    end
    return state
end

function _step_rk3!(eq::FrontEquation, state::FrontState, t::Float64, dt::Float64)
    V  = _get_vel_buf(eq, state)
    x0 = vertex_coordinates(state)
    T  = eltype(eltype(x0))
    dt_T = T(dt)

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
    new_pts2 = [T(3)/4 * x0[a] + T(1)/4 * x1[a] + T(1)/4 * dt_T * k2[a]
                for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts2)
    refresh_geometry!(state)

    # Stage 3: k3 = f(x2, t + dt/2)
    compute_rhs!(V, eq.terms, state, t + dt/2)
    k3 = copy(V)

    # x_new = (1/3)*x0 + (2/3)*x2 + (2/3)*dt*k3
    x2 = vertex_coordinates(state)
    new_pts = [T(1)/3 * x0[a] + T(2)/3 * x2[a] + T(2)/3 * dt_T * k3[a]
               for a in eachindex(x0)]
    state.mesh = _replace_mesh(state.mesh, new_pts)
    refresh_geometry!(state)
end

function _step_rk3!(eq::FrontEquation, state::MultiFrontState, t::Float64, dt::Float64)
    for (i, comp) in enumerate(eachcomponent(state))
        cstate = _component_front_state(comp, t)
        V  = _get_vel_buf(eq, state, i)
        x0 = vertex_coordinates(cstate)
        T  = eltype(eltype(x0))
        dt_T = T(dt)

        compute_rhs!(V, eq.terms, cstate, t)
        k1 = copy(V)

        _apply_velocity!(cstate, x0, k1, dt)
        x1 = vertex_coordinates(cstate)

        compute_rhs!(V, eq.terms, cstate, t + dt)
        k2 = copy(V)

        new_pts2 = [T(3)/4 * x0[a] + T(1)/4 * x1[a] + T(1)/4 * dt_T * k2[a]
                    for a in eachindex(x0)]
        cstate.mesh = _replace_mesh(cstate.mesh, new_pts2)
        refresh_geometry!(cstate)

        compute_rhs!(V, eq.terms, cstate, t + dt/2)
        k3 = copy(V)

        x2 = vertex_coordinates(cstate)
        new_pts = [T(1)/3 * x0[a] + T(2)/3 * x2[a] + T(2)/3 * dt_T * k3[a]
                   for a in eachindex(x0)]
        cstate.mesh = _replace_mesh(cstate.mesh, new_pts)
        refresh_geometry!(cstate)

        _sync_component_from_state!(comp, cstate)
    end
    return state
end

function _compute_cfl(eq::FrontEquation, state::FrontState, t::Float64)
    return compute_cfl(eq.terms, state, t)
end

function _compute_cfl(eq::FrontEquation, state::MultiFrontState, t::Float64)
    ncomponents(state) > 0 || return Inf
    return minimum(compute_cfl(eq.terms, _component_front_state(comp, t), t) for comp in eachcomponent(state))
end

function _maybe_redistribute!(eq::FrontEquation, state::FrontState)
    if eq.redistributor !== nothing
        redistribute!(state, eq.redistributor)
    end
    return state
end

function _maybe_redistribute!(eq::FrontEquation, state::MultiFrontState)
    if eq.redistributor !== nothing
        for comp in eachcomponent(state)
            cstate = _component_front_state(comp, state.t)
            redistribute!(cstate, eq.redistributor)
            _sync_component_from_state!(comp, cstate)
        end
    end
    return state
end

function _handle_topology_stage!(eq::FrontEquation, state::AbstractFrontState)
    report = handle_topology_change!(state, eq.topology_handler)
    state.cache[:last_topology_report] = report
    if report.changed
        if state isa FrontState
            refresh_geometry!(state)
        else
            refresh_geometry!(state; which=:all)
        end
    end
    return report
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
            dt_cfl = _compute_cfl(eq, state, t)
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

        # ── topology change stage ──────────────────────────────────────────
        _handle_topology_stage!(eq, state)

        # ── optional redistribution ─────────────────────────────────────────
        _maybe_redistribute!(eq, state)

        # ── callback ────────────────────────────────────────────────────────
        callback !== nothing && callback(state, t, step)
    end

    eq.step = step
    return eq
end
