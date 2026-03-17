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
    cfl(integrator::TimeIntegrator) -> Float64

Return the CFL safety factor of the integrator.
"""
cfl(integ::ForwardEuler) = integ.cfl
cfl(integ::RK2)          = integ.cfl
cfl(integ::RK3)          = integ.cfl

Base.show(io::IO, ::ForwardEuler) = print(io, "ForwardEuler")
Base.show(io::IO, ::RK2)          = print(io, "RK2 (Heun)")
Base.show(io::IO, ::RK3)          = print(io, "RK3 (Shu-Osher TVD)")
