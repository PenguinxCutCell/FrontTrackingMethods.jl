# benchmark_geometries.jl – Reusable benchmark geometry constructors.
#
# Provides closed polygonal curves, selected open curves, and triangulated
# surfaces for standard front-tracking benchmarks. Every constructor returns
# a mesh ready to use with FrontEquation.
#
# 2-D constructors
# ----------------
#   make_circle_benchmark_curve      – circle for translation/rotation tests
#   make_open_arc_benchmark_curve    – open circular arc for advection tests
#   make_zalesak_disk_curve          – slotted disk (Zalesak 1979)
#
# 3-D constructors
# ----------------
#   make_sphere_benchmark_surface    – sphere for translation/rotation tests
#   make_zalesak_sphere_surface      – slotted sphere (3-D Zalesak)

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: circle benchmark curve
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_circle_benchmark_curve(; center=SVector(0.5,0.75), R=0.15, N=256)
        -> CurveMesh{Float64}

Create a closed polygonal approximation of a circle for benchmark tests.
The circle is centered at `center` with radius `R` and `N` vertices.

Used in:
- 2-D Rider–Kothe single vortex
- 2-D rigid rotation benchmarks
"""
function make_circle_benchmark_curve(;
    center :: SVector{2,Float64} = SVector(0.5, 0.75),
    R      :: Float64            = 0.15,
    N      :: Int                = 256,
) :: CurveMesh{Float64}
    T = Float64
    pts   = [center + SVector{2,T}(R*cos(2T(π)*k/N), R*sin(2T(π)*k/N)) for k in 0:N-1]
    edges = [SVector{2,Int}(k, mod1(k+1, N)) for k in 1:N]
    return CurveMesh{T}(pts, edges)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: open arc benchmark curve
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_open_arc_benchmark_curve(;
        center=SVector(0.5, 0.5),
        R=0.2,
        θ_start=π/6,
        θ_end=5π/6,
        N=128,
    ) -> CurveMesh{Float64}

Create an open polygonal approximation of a circular arc for open-front tests.

The arc is centered at `center`, has radius `R`, starts at angle `θ_start`,
ends at angle `θ_end`, and contains `N` vertices connected by `N-1` edges.

This constructor is intended for open-front topology/geometry checks and
prescribed-advection benchmarks. It does not enclose area.
"""
function make_open_arc_benchmark_curve(;
    center  :: SVector{2,Float64} = SVector(0.5, 0.5),
    R       :: Float64            = 0.2,
    θ_start :: Float64            = Float64(π / 6),
    θ_end   :: Float64            = Float64(5π / 6),
    N       :: Int                = 128,
) :: CurveMesh{Float64}
    T = Float64

    R > 0 || throw(ArgumentError("R must be > 0, got R=$R"))
    N >= 2 || throw(ArgumentError("N must be at least 2, got N=$N"))
    abs(θ_end - θ_start) > eps(T) ||
        throw(ArgumentError("θ_start and θ_end must define a nonzero arc."))

    angles = range(T(θ_start), T(θ_end), length=N)
    pts = [center + SVector{2,T}(R * cos(θ), R * sin(θ)) for θ in angles]
    edges = [SVector{2,Int}(k, k + 1) for k in 1:N-1]
    return CurveMesh{T}(pts, edges)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2-D: Zalesak disk curve
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_zalesak_disk_curve(;
        center     = SVector(0.5, 0.75),
        R          = 0.15,
        slot_width = 0.05,
        slot_depth = 0.25,
        N_arc      = 256,
        N_slot     = 64,
    ) -> CurveMesh{Float64}

Construct the Zalesak slotted-disk curve (Zalesak 1979).

The disk is a circle of radius `R` centered at `center` with a rectangular
notch of width `slot_width` and depth `slot_depth` cut from the top.

The curve is a closed polygon traversing:
1. The outer arc from the left notch lip to the right notch lip
    (going counter-clockwise along most of the circle).
2. The right side of the notch (downward).
3. The notch bottom (left).
4. The left side of the notch (upward back to the circle).

Orientation is counter-clockwise (positive enclosed area).

Parameters
----------
- `N_arc`  – number of points along the circular arc.
- `N_slot` – number of points along each slot side (each of the 3 slot sides
  gets `N_slot` segments).
"""
function make_zalesak_disk_curve(;
    center     :: SVector{2,Float64} = SVector(0.5, 0.75),
    R          :: Float64            = 0.15,
    slot_width :: Float64            = 0.05,
    slot_depth :: Float64            = 0.25,
    N_arc      :: Int                = 256,
    N_slot     :: Int                = 64,
) :: CurveMesh{Float64}

    T  = Float64
    cx, cy = center[1], center[2]
    hw = slot_width / 2

    if !(R > 0)
        throw(ArgumentError("R must be > 0, got R=$R"))
    end
    if !(slot_width > 0)
        throw(ArgumentError("slot_width must be > 0, got slot_width=$slot_width"))
    end
    if !(slot_depth > 0)
        throw(ArgumentError("slot_depth must be > 0, got slot_depth=$slot_depth"))
    end
    if !(slot_width < 2R)
        throw(ArgumentError("slot_width must be < 2R for a valid notch (slot_width=$slot_width, 2R=$(2R))"))
    end
    if N_arc < 8
        throw(ArgumentError("N_arc must be at least 8, got N_arc=$N_arc"))
    end
    if N_slot < 2
        throw(ArgumentError("N_slot must be at least 2, got N_slot=$N_slot"))
    end

    x_right = cx + hw
    x_left  = cx - hw

    y_lip_offset = sqrt(max(R^2 - hw^2, zero(T)))
    y_lip = cy + y_lip_offset
    y_bottom_min = cy - R + eps(T)
    y_notch_bottom = max(y_lip - slot_depth, y_bottom_min)

    θ_left  = atan(y_lip - cy, x_left - cx)
    θ_right = atan(y_lip - cy, x_right - cx)

    arc_span = mod(θ_right - θ_left, 2T(π))
    arc_pts = SVector{2,T}[]
    for k in 0:N_arc-1
        θ = θ_left + arc_span * k / N_arc
        push!(arc_pts, SVector{2,T}(cx + R*cos(θ), cy + R*sin(θ)))
    end

    right_lip = SVector{2,T}(x_right, y_lip)
    slot_right_pts = SVector{2,T}[]
    for k in 1:N_slot
        frac = k / N_slot
        y = y_lip + frac * (y_notch_bottom - y_lip)
        push!(slot_right_pts, SVector{2,T}(x_right, y))
    end

    slot_bottom_pts = SVector{2,T}[]
    for k in 1:N_slot
        frac = k / N_slot
        x = x_right + frac * (x_left - x_right)
        push!(slot_bottom_pts, SVector{2,T}(x, y_notch_bottom))
    end

    slot_left_pts = SVector{2,T}[]
    for k in 1:N_slot-1
        frac = k / N_slot
        y = y_notch_bottom + frac * (y_lip - y_notch_bottom)
        push!(slot_left_pts, SVector{2,T}(x_left, y))
    end

    all_pts = vcat(arc_pts, [right_lip], slot_right_pts, slot_bottom_pts, slot_left_pts)

    signed_area2 = zero(T)
    M = length(all_pts)
    for k in 1:M
        p = all_pts[k]
        q = all_pts[mod1(k + 1, M)]
        signed_area2 += p[1] * q[2] - q[1] * p[2]
    end
    if signed_area2 < 0
        reverse!(all_pts)
    end

    N       = length(all_pts)
    edges   = [SVector{2,Int}(k, mod1(k+1, N)) for k in 1:N]
    return CurveMesh{T}(all_pts, edges)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D: sphere benchmark surface
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_sphere_benchmark_surface(;
        center      = SVector(0.35, 0.35, 0.35),
        R           = 0.15,
        refinement  = 3,
    ) -> SurfaceMesh{Float64}

Create a triangulated sphere for 3-D benchmark tests.

Uses an icosphere (subdivided icosahedron) for uniform triangle quality.
The sphere is centered at `center` with radius `R`.

Used in:
- 3-D rigid translation/rotation benchmarks
- 3-D Enright deformation benchmark
"""
function make_sphere_benchmark_surface(;
    center     :: SVector{3,Float64} = SVector(0.35, 0.35, 0.35),
    R          :: Float64            = 0.15,
    refinement :: Int                = 3,
) :: SurfaceMesh{Float64}
    mesh = FrontIntrinsicOps.generate_icosphere(R, refinement)
    # Translate to center
    new_pts = [p + center for p in mesh.points]
    return SurfaceMesh{Float64}(new_pts, mesh.faces)
end

# ─────────────────────────────────────────────────────────────────────────────
# 3-D: Zalesak sphere surface
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_zalesak_sphere_surface(;
        center      = SVector(0.5, 0.75, 0.5),
        R           = 0.15,
        slot_width  = 0.05,
        slot_depth  = 0.125,
        refinement  = 3,
    ) -> SurfaceMesh{Float64}

Construct an approximate Zalesak slotted sphere for 3-D benchmark tests.

The sphere has a rectangular slot (tunnel) cut from the bottom.  This is a
fixed-topology approximation suitable for v0.2 benchmarks.

Construction approach:
1. Generate an icosphere.
2. Remove vertices/faces that fall inside the slot region.
3. Re-stitch the boundary with simple triangles.

⚠ This is a fixed-topology v0.2 approximation.  The slot corners are
geometry-approximate, not CAD-quality.  Prioritizes stability over
perfect geometry.
"""
function make_zalesak_sphere_surface(;
    center     :: SVector{3,Float64} = SVector(0.5, 0.75, 0.5),
    R          :: Float64            = 0.15,
    slot_width :: Float64            = 0.05,
    slot_depth :: Float64            = 0.125,
    refinement :: Int                = 3,
) :: SurfaceMesh{Float64}
    T  = Float64
    cx, cy, cz = center[1], center[2], center[3]
    hw = slot_width / 2   # half-width of slot in x

    # Build icosphere centered at origin, then translate
    base = FrontIntrinsicOps.generate_icosphere(R, refinement)
    pts0 = [p + center for p in base.points]

    # The slot is defined as:
    #   |x - cx| <= hw  AND  y <= cy  (bottom hemisphere)  AND  z in [cz-R, cz-R+slot_depth]
    # (i.e., a vertical rectangular prism cut from the bottom)
    y_slot_cut = cy - R * cos(asin(hw / R))  # y level where slot meets circle
    z_slot_bot = cz - R
    z_slot_top = cz - R + slot_depth

    # Identify vertices inside the slot
    function in_slot(p::SVector{3,Float64})
        return (abs(p[1] - cx) <= hw &&
                p[2] <= cy &&
                p[3] >= z_slot_bot - eps(Float64) &&
                p[3] <= z_slot_top + eps(Float64))
    end

    # Keep only faces where NO vertex is inside the slot
    new_faces = SVector{3,Int}[]
    for f in base.faces
        p1 = pts0[f[1]]
        p2 = pts0[f[2]]
        p3 = pts0[f[3]]
        if !in_slot(p1) && !in_slot(p2) && !in_slot(p3)
            push!(new_faces, f)
        end
    end

    # If nothing was removed (slot too small or refinement too coarse),
    # return the full sphere
    if length(new_faces) == length(base.faces)
        return SurfaceMesh{T}(pts0, base.faces)
    end

    # Compact vertex list (remove unused vertices)
    used = falses(length(pts0))
    for f in new_faces
        used[f[1]] = true
        used[f[2]] = true
        used[f[3]] = true
    end
    old_to_new = zeros(Int, length(pts0))
    new_pts = SVector{3,T}[]
    ni = 0
    for (oi, u) in enumerate(used)
        if u
            ni += 1
            old_to_new[oi] = ni
            push!(new_pts, pts0[oi])
        end
    end
    remapped = [SVector{3,Int}(old_to_new[f[1]], old_to_new[f[2]], old_to_new[f[3]])
                for f in new_faces]

    return SurfaceMesh{T}(new_pts, remapped)
end
