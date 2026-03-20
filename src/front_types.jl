# front_types.jl – Core abstract types for FrontTrackingMethods.

"""
    AbstractFrontTerm

Abstract supertype for all motion terms in a front evolution equation.
Concrete subtypes: `AdvectionTerm`, `NormalMotionTerm`, `CurvatureMotionTerm`.
"""
abstract type AbstractFrontTerm end

"""
    TimeIntegrator

Abstract supertype for time integration schemes.
Concrete subtypes: `ForwardEuler`, `RK2`, `RK3`.
"""
abstract type TimeIntegrator end

"""
    AbstractRedistributor

Abstract supertype for front-quality redistribution strategies.
Redistribution is the front-tracking analogue of reinitialization:
it moves vertices tangentially (and in v0.1 never changes connectivity)
to improve mesh quality without changing the physical front shape.

Concrete subtypes: `NoRedistribution`, `CurveEqualArcRedistributor`,
`SurfaceTangentialRedistributor`.
"""
abstract type AbstractRedistributor end

"""
    AbstractFieldTransfer

Abstract supertype for field-transfer strategies used after redistribution
moves vertices.  Field transfer is the front-tracking analogue of velocity
extension.

Concrete subtypes: `ClosestPointTransfer`.
"""
abstract type AbstractFieldTransfer end

# ── Internal front-kind helpers ───────────────────────────────────────────────

is_supported_front(::CurveMesh)   = true
is_supported_front(::SurfaceMesh) = true
is_supported_front(::PointFront1D) = true
is_supported_front(::Any)         = false

ambient_dimension(::CurveMesh)    = 2
ambient_dimension(::SurfaceMesh)  = 3
ambient_dimension(::PointFront1D) = 1

nmarkers(front::CurveMesh)    = length(front.points)
nmarkers(front::SurfaceMesh)  = length(front.points)
nmarkers(front::PointFront1D) = length(front.x)

front_markers(front::CurveMesh)    = front.points
front_markers(front::SurfaceMesh)  = front.points
front_markers(front::PointFront1D) = front.x
