# timestepping.jl – Time integration schemes.
#
# Mirrors the LevelSetMethods naming.  Each integrator carries a CFL number
# that is multiplied by the term-based stable time-step estimate.

"""
    ForwardEuler(; cfl=0.5)

First-order explicit (Euler) time integrator.

`cfl` is the safety factor (in (0, 1]) applied to the CFL-based step estimate.
"""
struct ForwardEuler <: TimeIntegrator
    cfl :: Float64
end
ForwardEuler(; cfl::Real=0.5) = ForwardEuler(Float64(cfl))

"""
    RK2(; cfl=0.8)

Second-order explicit Runge-Kutta (Heun's method).

`cfl` is the safety factor applied to the CFL-based step estimate.
"""
struct RK2 <: TimeIntegrator
    cfl :: Float64
end
RK2(; cfl::Real=0.8) = RK2(Float64(cfl))

"""
    RK3(; cfl=1.0)

Third-order Shu-Osher TVD Runge-Kutta.

`cfl` is the safety factor applied to the CFL-based step estimate.
"""
struct RK3 <: TimeIntegrator
    cfl :: Float64
end
RK3(; cfl::Real=1.0) = RK3(Float64(cfl))

"""
    DiffEqIntegrator(algorithm; cfl=0.8, kwargs...)

Adaptive OrdinaryDiffEq.jl-backed time integrator.

`algorithm` is an OrdinaryDiffEq algorithm object, for example `Tsit5()`.
Extra keyword arguments are forwarded to `OrdinaryDiffEq.solve`, so users can
set tolerances such as `abstol` and `reltol`.

This is a weak-extension integrator. Load OrdinaryDiffEq before using it:

```julia
using OrdinaryDiffEq
integrator = DiffEqIntegrator(Tsit5(); abstol=1e-10, reltol=1e-10)
```
"""
struct DiffEqIntegrator{A,K} <: TimeIntegrator
    algorithm :: A
    cfl       :: Float64
    kwargs    :: K
end

function DiffEqIntegrator(algorithm; cfl::Real=0.8, kwargs...)
    return DiffEqIntegrator(algorithm, Float64(cfl), (; kwargs...))
end

"""
    cfl(integrator::TimeIntegrator) -> Float64

Return the CFL safety factor of the integrator.
"""
cfl(integ::ForwardEuler) = integ.cfl
cfl(integ::RK2)          = integ.cfl
cfl(integ::RK3)          = integ.cfl
cfl(integ::DiffEqIntegrator) = integ.cfl

Base.show(io::IO, ::ForwardEuler) = print(io, "ForwardEuler")
Base.show(io::IO, ::RK2)          = print(io, "RK2 (Heun)")
Base.show(io::IO, ::RK3)          = print(io, "RK3 (Shu-Osher TVD)")
Base.show(io::IO, integ::DiffEqIntegrator) =
    print(io, "DiffEqIntegrator ($(typeof(integ.algorithm)))")
