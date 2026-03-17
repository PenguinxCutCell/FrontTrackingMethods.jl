using StaticArrays

function sample_velocity_on_faces(
    vel_analytic;
    t::Float64,
    state=nothing,
    xmin::Float64=0.0,
    xmax::Float64=1.0,
    ymin::Float64=0.0,
    ymax::Float64=1.0,
    Nx::Int=128,
    Ny::Int=128,
)
    Nx >= 2 || throw(ArgumentError("Nx must be ≥ 2, got $Nx"))
    Ny >= 2 || throw(ArgumentError("Ny must be ≥ 2, got $Ny"))

    dx = (xmax - xmin) / Nx
    dy = (ymax - ymin) / Ny

    uface = Matrix{Float64}(undef, Nx + 1, Ny)
    for i in 1:(Nx + 1)
        x = xmin + (i - 1) * dx
        for j in 1:Ny
            y = ymin + (j - 0.5) * dy
            uface[i, j] = vel_analytic(SVector{2,Float64}(x, y), t, state)[1]
        end
    end

    vface = Matrix{Float64}(undef, Nx, Ny + 1)
    for i in 1:Nx
        x = xmin + (i - 0.5) * dx
        for j in 1:(Ny + 1)
            y = ymin + (j - 1) * dy
            vface[i, j] = vel_analytic(SVector{2,Float64}(x, y), t, state)[2]
        end
    end

    return (
        uface=uface,
        vface=vface,
        xmin=xmin,
        xmax=xmax,
        ymin=ymin,
        ymax=ymax,
        dx=dx,
        dy=dy,
        Nx=Nx,
        Ny=Ny,
    )
end

@inline function _bilinear(A::AbstractMatrix{<:Real}, i0::Int, i1::Int, j0::Int, j1::Int, ξ::Float64, η::Float64)
    return (1 - ξ) * (1 - η) * A[i0, j0] +
           ξ * (1 - η) * A[i1, j0] +
           (1 - ξ) * η * A[i0, j1] +
           ξ * η * A[i1, j1]
end

@inline function _clamp_pair(k::Int, n::Int)
    if k < 1
        return 1, 1
    elseif k >= n
        return n, n
    else
        return k, k + 1
    end
end

function _interp_u(faces, x::Float64, y::Float64)
    gx = (x - faces.xmin) / faces.dx + 1.0
    gy = (y - (faces.ymin + 0.5 * faces.dy)) / faces.dy + 1.0

    i0raw = floor(Int, gx)
    j0raw = floor(Int, gy)
    i0, i1 = _clamp_pair(i0raw, faces.Nx + 1)
    j0, j1 = _clamp_pair(j0raw, faces.Ny)

    ξ = i0 == i1 ? 0.0 : clamp(gx - i0, 0.0, 1.0)
    η = j0 == j1 ? 0.0 : clamp(gy - j0, 0.0, 1.0)
    return _bilinear(faces.uface, i0, i1, j0, j1, ξ, η)
end

function _interp_v(faces, x::Float64, y::Float64)
    gx = (x - (faces.xmin + 0.5 * faces.dx)) / faces.dx + 1.0
    gy = (y - faces.ymin) / faces.dy + 1.0

    i0raw = floor(Int, gx)
    j0raw = floor(Int, gy)
    i0, i1 = _clamp_pair(i0raw, faces.Nx)
    j0, j1 = _clamp_pair(j0raw, faces.Ny + 1)

    ξ = i0 == i1 ? 0.0 : clamp(gx - i0, 0.0, 1.0)
    η = j0 == j1 ? 0.0 : clamp(gy - j0, 0.0, 1.0)
    return _bilinear(faces.vface, i0, i1, j0, j1, ξ, η)
end

function face_bilinear_marker_velocity(
    vel_analytic;
    xmin::Float64=0.0,
    xmax::Float64=1.0,
    ymin::Float64=0.0,
    ymax::Float64=1.0,
    Nx::Int=128,
    Ny::Int=128,
)
    has_faces = Ref(false)
    last_t = Ref(typemin(Float64))
    faces_ref = Ref{Any}(nothing)

    function u_interp(x, t, state)
        tt = Float64(t)
        if !has_faces[] || !isapprox(tt, last_t[]; atol=0.0, rtol=0.0)
            faces_ref[] = sample_velocity_on_faces(
                vel_analytic;
                t=tt,
                state=state,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                Nx=Nx,
                Ny=Ny,
            )
            last_t[] = tt
            has_faces[] = true
        end

        xx = clamp(Float64(x[1]), xmin, xmax)
        yy = clamp(Float64(x[2]), ymin, ymax)
        faces = faces_ref[]
        return SVector{2,Float64}(_interp_u(faces, xx, yy), _interp_v(faces, xx, yy))
    end

    return u_interp
end
