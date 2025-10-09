"""
    SpecialUnitaryGroup{T}

The special orthogonal group ``$(_math(:SU))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`GeneralUnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) with determinant one.

# Constructor
    SpecialUnitaryGroup(n::Int; kwargs...)

Generate special unitary group ``$(_math(:SU))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialUnitaryGroup{T} = LieGroup{
    ManifoldsBase.ℂ,
    MatrixMultiplicationGroupOperation,
    GeneralUnitaryMatrices{ℂ, T, DeterminantOneMatrixType},
}

function SpecialUnitaryGroup(n::Int; kwargs...)
    GU = GeneralUnitaryMatrices(n, ℂ, DeterminantOneMatrixType; kwargs...)
    return SpecialUnitaryGroup{typeof(GU).parameters[2]}(
        GU, MatrixMultiplicationGroupOperation()
    )
end

function Base.show(io::IO, G::SpecialUnitaryGroup)
    size = get_parameter(G.manifold.size)[1]
    return print(io, "SpecialUnitaryGroup($(size))")
end
