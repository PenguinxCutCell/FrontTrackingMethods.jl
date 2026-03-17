# callbacks.jl – Lightweight callback support for front-tracking simulations.
#
# Callbacks are small callables that are invoked during integrate!.
# They are useful for:
#   - saving diagnostics
#   - triggering redistribution/remeshing every N steps
#   - stopping when mesh quality becomes unacceptable
#   - writing CSV snapshots
#
# Two lightweight callback types are provided:
#   EveryNSteps(n, f)             – call f every n time steps
#   TimeIntervalCallback(Δt, f)   – call f approximately every Δt in time
#
# Both have the signature: f(state, t, step)
#
# Usage with integrate!:
#   cb = EveryNSteps(10, (state, t, step) -> println("step $step, t=$t"))
#   integrate!(eq, tf; callback=cb)
#
# Or combine multiple callbacks with compose_callbacks:
#   cb = compose_callbacks(cb1, cb2, cb3)

# ─────────────────────────────────────────────────────────────────────────────
# EveryNSteps
# ─────────────────────────────────────────────────────────────────────────────

"""
    EveryNSteps(n::Int, f)

A callback that calls `f(state, t, step)` every `n` time steps.

Example:
    cb = EveryNSteps(10, (state, t, step) -> println("step=\$step t=\$t"))
"""
struct EveryNSteps
    n :: Int
    f :: Any
end

function (cb::EveryNSteps)(state, t, step)
    if mod(step, cb.n) == 0
        cb.f(state, t, step)
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# TimeIntervalCallback
# ─────────────────────────────────────────────────────────────────────────────

"""
    TimeIntervalCallback(Δtout::Float64, f)

A callback that calls `f(state, t, step)` approximately every `Δtout` in
simulation time.

The callback fires when `t` has advanced by at least `Δtout` since the last
fire.

Example:
    cb = TimeIntervalCallback(0.1, (state, t, step) -> println("t=\$t"))
"""
mutable struct TimeIntervalCallback
    Δtout     :: Float64
    f         :: Any
    t_last    :: Float64
end

TimeIntervalCallback(Δtout::Real, f) = TimeIntervalCallback(Float64(Δtout), f, -Inf)

function (cb::TimeIntervalCallback)(state, t, step)
    if t >= cb.t_last + cb.Δtout - eps(Float64)
        cb.f(state, t, step)
        cb.t_last = t
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# compose_callbacks
# ─────────────────────────────────────────────────────────────────────────────

"""
    compose_callbacks(callbacks...) -> function

Combine multiple callbacks into a single callable.
Each callback is called in order with the same `(state, t, step)` arguments.

Example:
    cb = compose_callbacks(
        EveryNSteps(10, print_diagnostics),
        TimeIntervalCallback(0.5, save_snapshot),
    )
"""
function compose_callbacks(callbacks...)
    return (state, t, step) -> begin
        for cb in callbacks
            cb(state, t, step)
        end
        return nothing
    end
end
