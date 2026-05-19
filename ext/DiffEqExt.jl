module DiffEqExt

using OrdinaryDiffEq
using StaticArrays

import FrontIntrinsicOps as FIO
import FrontTrackingMethods as FTM

function _flatten_mesh(mesh::FIO.CurveMesh{T}) where {T}
    z = Vector{T}(undef, 2 * length(mesh.points))
    for (i, p) in pairs(mesh.points)
        j = 2i - 1
        z[j] = p[1]
        z[j + 1] = p[2]
    end
    return z
end

function _flatten_mesh(mesh::FIO.SurfaceMesh{T}) where {T}
    z = Vector{T}(undef, 3 * length(mesh.points))
    for (i, p) in pairs(mesh.points)
        j = 3i - 2
        z[j] = p[1]
        z[j + 1] = p[2]
        z[j + 2] = p[3]
    end
    return z
end

_flatten_mesh(mesh::FIO.PointFront1D{T}) where {T} = T.(mesh.x)

function _flatten_velocity!(dz, V::AbstractVector{<:SVector{2}})
    for (i, v) in pairs(V)
        j = 2i - 1
        dz[j] = v[1]
        dz[j + 1] = v[2]
    end
    return dz
end

function _flatten_velocity!(dz, V::AbstractVector{<:SVector{3}})
    for (i, v) in pairs(V)
        j = 3i - 2
        dz[j] = v[1]
        dz[j + 1] = v[2]
        dz[j + 2] = v[3]
    end
    return dz
end

function _flatten_velocity!(dz, V::AbstractVector{<:Real})
    for i in eachindex(V)
        dz[i] = V[i]
    end
    return dz
end

function _points_from_vector(mesh::FIO.CurveMesh{T}, z) where {T}
    n = length(mesh.points)
    length(z) == 2n || error("DiffEqIntegrator: state vector length mismatch for CurveMesh.")
    return [SVector{2,T}(T(z[2i - 1]), T(z[2i])) for i in 1:n]
end

function _points_from_vector(mesh::FIO.SurfaceMesh{T}, z) where {T}
    n = length(mesh.points)
    length(z) == 3n || error("DiffEqIntegrator: state vector length mismatch for SurfaceMesh.")
    return [SVector{3,T}(T(z[3i - 2]), T(z[3i - 1]), T(z[3i])) for i in 1:n]
end

function _points_from_vector(mesh::FIO.PointFront1D{T}, z) where {T}
    length(z) == length(mesh.x) || error("DiffEqIntegrator: state vector length mismatch for PointFront1D.")
    return [T(zi) for zi in z]
end

function _restore_state_from_vector!(state::FTM.FrontState, z, t)
    state.mesh = FTM._replace_mesh(state.mesh, _points_from_vector(state.mesh, z))
    state.t = Float64(t)
    FTM.refresh_geometry!(state)
    return state
end

function FTM._step_integrator!(
    eq::FTM.FrontEquation,
    state::FTM.FrontState,
    t::Float64,
    dt::Float64,
    integ::FTM.DiffEqIntegrator,
)
    if dt <= 10 * eps(max(abs(t), abs(t + dt), 1.0))
        state.t = t + dt
        return state
    end

    z0 = _flatten_mesh(state.mesh)
    V = FTM._zeros_like_points(state.mesh)
    params = (eq=eq, state=state, V=V)

    function rhs!(dz, z, p, t_eval)
        _restore_state_from_vector!(p.state, z, t_eval)
        FTM.compute_rhs!(p.V, p.eq.terms, p.state, t_eval)
        _flatten_velocity!(dz, p.V)
        return nothing
    end

    tspan = (t, t + dt)
    prob = OrdinaryDiffEq.ODEProblem(rhs!, z0, tspan, params)
    solve_kwargs = merge((save_everystep=false,), integ.kwargs)
    sol = OrdinaryDiffEq.solve(prob, integ.algorithm; solve_kwargs...)
    if hasproperty(sol, :retcode) && !occursin("Success", string(sol.retcode))
        error("DiffEqIntegrator solve failed with retcode $(sol.retcode).")
    end
    _restore_state_from_vector!(state, sol.u[end], t + dt)
    return state
end

end # module DiffEqExt
