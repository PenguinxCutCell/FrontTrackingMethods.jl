# benchmark_fields.jl – Reusable benchmark velocity fields.
#
# Provides standard benchmark velocity fields for front-tracking tests.
# Each function returns a callable `(x, t, state) -> velocity_vector`.
#
# 2-D fields
# ----------
#   rigid_translation_velocity  – constant translation
#   rigid_rotation_2d           – solid-body rotation in 2-D
#   rider_kothe_single_vortex   – Rider–Kothe reversed vortex (Rider & Kothe 1998)
#   deformation16_2d            – 16-vortex strong deformation
#   serpentine_2d               – serpentine shear deformation
#
# 3-D fields
# ----------
#   rigid_rotation_3d           – solid-body rotation about an axis in 3-D
#   enright_3d                  – Enright deformation benchmark (Enright et al. 2002)
#
# Stream-function helpers
# -----------------------
#   _stream2vel_2d(ψx, ψy) -> (u, v) = (ψy, -ψx)
#
# References
# ----------
# - Rider & Kothe (1998). Reconstructing Volume Tracking. J. Comput. Phys.
# - Enright, Fedkiw, Ferziger & Mitchell (2002). A Hybrid Particle Level Set
#   Method for Improved Interface Capturing. J. Comput. Phys.
# - Leveque (1996). High-resolution conservative algorithms for advection in
#   incompressible flow. SIAM J. Numer. Anal.
# - Barth & Sethian (1998).

# ─────────────────────────────────────────────────────────────────────────────
# Helper: stream function → velocity
# ─────────────────────────────────────────────────────────────────────────────

"""
    _stream2vel_2d(dψ_dy, dψ_dx) -> (u, v)

Convert stream-function partial derivatives to 2-D velocity:
    u =  ∂ψ/∂y
    v = -∂ψ/∂x
"""
@inline _stream2vel_2d(dψ_dy, dψ_dx) = (dψ_dy, -dψ_dx)

# ─────────────────────────────────────────────────────────────────────────────
# Rigid translation
# ─────────────────────────────────────────────────────────────────────────────

"""
    rigid_translation_velocity(u::SVector) -> (x, t, state) -> SVector

Return a callable that applies constant rigid-body translation velocity `u`
everywhere.

Domain: any
Reversible: no (net displacement u*T after time T)
"""
function rigid_translation_velocity(u::SVector)
    return (x, t, state) -> u
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: rigid rotation
# ─────────────────────────────────────────────────────────────────────────────

"""
    rigid_rotation_2d(; center=SVector(0.5,0.5), omega=2π)
        -> (x, t, state) -> SVector{2,Float64}

Return a callable for 2-D solid-body rotation about `center` with angular
velocity `omega` (rad/time).

Velocity field:
    u(x, y) = omega * (-(y - cy), (x - cx))

Domain: ℝ²
Period: T = 2π / |omega|
Reversibility: exact after one period T
"""
function rigid_rotation_2d(;
    center :: SVector{2,Float64} = SVector(0.5, 0.5),
    omega  :: Float64            = 2π,
)
    cx, cy = center[1], center[2]
    return (x, t, state) -> begin
        dx = x[1] - cx
        dy = x[2] - cy
        SVector{2,Float64}(-omega * dy, omega * dx)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D: rigid rotation
# ─────────────────────────────────────────────────────────────────────────────

"""
    rigid_rotation_3d(; center=SVector(0.5,0.5,0.5), axis=:z, omega=2π)
        -> (x, t, state) -> SVector{3,Float64}

Return a callable for 3-D solid-body rotation about an axis through `center`.

`axis` can be `:x`, `:y`, or `:z`, or an `SVector{3}` unit vector.

Velocity field:
    u(x) = omega * (axis_vec × (x - center))
"""
function rigid_rotation_3d(;
    center :: SVector{3,Float64} = SVector(0.5, 0.5, 0.5),
    axis                         = :z,
    omega  :: Float64            = 2π,
)
    axis_vec = if axis === :x
        SVector{3,Float64}(1.0, 0.0, 0.0)
    elseif axis === :y
        SVector{3,Float64}(0.0, 1.0, 0.0)
    elseif axis === :z
        SVector{3,Float64}(0.0, 0.0, 1.0)
    elseif axis isa SVector{3}
        normalize(axis)
    else
        error("rigid_rotation_3d: axis must be :x, :y, :z, or an SVector{3}")
    end
    return (x, t, state) -> begin
        r = x - center
        omega * (axis_vec × r)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: Rider–Kothe reversed single vortex
# ─────────────────────────────────────────────────────────────────────────────

"""
    rider_kothe_single_vortex(; T=2.0)
        -> (x, t, state) -> SVector{2,Float64}

Rider–Kothe reversible single-vortex field on [0,1]².

Stream function:
    ψ(x, y, t) = (1/π) sin²(πx) sin²(πy) cos(πt/T)

Velocity:
    u =  ∂ψ/∂y = (2/π) sin²(πx) sin(πy) cos(πy) cos(πt/T)
    v = -∂ψ/∂x = -(2/π) sin(πx) cos(πx) sin²(πy) cos(πt/T)

Domain: [0,1]²
Final time: T (one full reversal cycle)
Initial geometry: circle at (0.5, 0.75) with R = 0.15

Reversibility:
    At t = T/2, cos(πt/T) = 0 → velocity = 0, maximum deformation.
    At t = T,   cos(πT/T) = cos(π) = -1 → reversed flow returns to IC.

Reference: Rider & Kothe (1998), Figure 3.
"""
function rider_kothe_single_vortex(; T::Float64=2.0)
    return (x, t, state) -> begin
        xv, yv = x[1], x[2]
        amp  = cos(π * t / T)
        u    = (2.0/π) * sin(π*xv)^2 * sin(π*yv) * cos(π*yv) * amp
        v    = -(2.0/π) * sin(π*xv) * cos(π*xv) * sin(π*yv)^2 * amp
        SVector{2,Float64}(u, v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: 16-vortex / strong deformation
# ─────────────────────────────────────────────────────────────────────────────

"""
    deformation16_2d(; T=2.0)
        -> (x, t, state) -> SVector{2,Float64}

Strong deformation field on [0,1]² with cos(πt/T) time modulation.

Stream function (4-vortex version of Rider–Kothe):
    ψ(x, y, t) = (1/(4π)) sin⁴(πx) sin⁴(πy) cos(πt/T)

Velocity:
    u =  ∂ψ/∂y = (1/π) sin⁴(πx) sin³(πy) cos(πy) cos(πt/T)
    v = -∂ψ/∂x = -(1/π) sin³(πx) cos(πx) sin⁴(πy) cos(πt/T)

Domain: [0,1]²
Final time: T
Reversibility: exact at time T.

Reference: LeVeque (1996) strong deformation test.
"""
function deformation16_2d(; T::Float64=2.0)
    return (x, t, state) -> begin
        xv, yv = x[1], x[2]
        amp  = cos(π * t / T)
        s4x  = sin(π*xv)^4
        s3y  = sin(π*yv)^3
        csy  = cos(π*yv)
        s3x  = sin(π*xv)^3
        csx  = cos(π*xv)
        s4y  = sin(π*yv)^4
        u    = (1.0/π) * s4x * s3y * csy * amp
        v    = -(1.0/π) * s3x * csx * s4y * amp
        SVector{2,Float64}(u, v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: serpentine deformation
# ─────────────────────────────────────────────────────────────────────────────

"""
    serpentine_2d(; Tmax=3.0)
        -> (x, t, state) -> SVector{2,Float64}

Serpentine (long-filament shear) deformation field on [0,1]².

Stream function:
    ψ(x, y, t) = (1/π) sin(2πx) sin²(πy) g(t)
    g(t) = cos(πt/Tmax)

Velocity:
    u =  ∂ψ/∂y = (2/π) sin(2πx) sin(πy) cos(πy) g(t)
            = (1/π) sin(2πx) sin(2πy) g(t)
    v = -∂ψ/∂x = -2 cos(2πx) sin²(πy) g(t)

Domain: [0,1]²
Final time: Tmax
Initial geometry: circle at (0.5, 0.75) with R = 0.15

Note: This is a severe shear/stretch test that can create long thin filaments.
Remeshing is essential for this benchmark.

Reference: Leveque (1996), serpentine test; Barth & Sethian (1998).
"""
function serpentine_2d(; Tmax::Float64=3.0)
    return (x, t, state) -> begin
        xv, yv = x[1], x[2]
        g    = cos(π * t / Tmax)
        u    = (1.0/π) * sin(2π*xv) * sin(2π*yv) * g
        v    = -2.0 * cos(2π*xv) * sin(π*yv)^2 * g
        SVector{2,Float64}(u, v)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D: Enright deformation
# ─────────────────────────────────────────────────────────────────────────────

"""
    enright_3d(; T=3.0)
        -> (x, t, state) -> SVector{3,Float64}

Enright 3-D deformation benchmark on [0,1]³.

Stream function (vector potential formulation):
    u(x,y,z,t) =  2 sin²(πx) sin(2πy) sin(2πz) cos(πt/T)
    v(x,y,z,t) = -sin(2πx) sin²(πy) sin(2πz) cos(πt/T)
    w(x,y,z,t) = -sin(2πx) sin(2πy) sin²(πz) cos(πt/T)

This is a divergence-free field.

Domain: [0,1]³
Final time: T
Initial geometry: sphere at (0.35, 0.35, 0.35) with R = 0.15

Reversibility: exact at time T.

Reference: Enright, Fedkiw, Ferziger & Mitchell (2002), Section 4.
"""
function enright_3d(; T::Float64=3.0)
    return (x, t, state) -> begin
        xv, yv, zv = x[1], x[2], x[3]
        amp = cos(π * t / T)
        u   = 2.0 * sin(π*xv)^2 * sin(2π*yv) * sin(2π*zv) * amp
        v   = -sin(2π*xv) * sin(π*yv)^2 * sin(2π*zv) * amp
        w   = -sin(2π*xv) * sin(2π*yv) * sin(π*zv)^2 * amp
        SVector{3,Float64}(u, v, w)
    end
end
