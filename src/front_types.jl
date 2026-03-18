# front_types.jl – Core abstract types for FrontTrackingMethods.

"""
    AbstractFrontState

Abstract supertype for front simulation states.
Concrete state types include `FrontState` (single component) and
`MultiFrontState` (multiple connected components).
"""
abstract type AbstractFrontState end

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
