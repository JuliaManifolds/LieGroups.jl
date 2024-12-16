"""
   UnitaryGroup{T}

The special orthogonal group ``$(_math(:U))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`UnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) with absolute value of the determinant equal to one.

# Constructor
   UnitaryGroup(n; kwargs...)

Generate unitary group ``$(_math(:U))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const UnitaryGroup{T} = LieGroup{
    ManifoldsBase.ℂ,
    MatrixMultiplicationGroupOperation,
    Manifolds.UnitaryMatrices{T,ManifoldsBase.ℂ},
}

function UnitaryGroup(n; kwargs...)
    U = Manifolds.GeneralUnitaryMatrices(
        n, ManifoldsBase.ℂ, Manifolds.AbsoluteDeterminantOneMatrices; kwargs...
    )
    return UnitaryGroup{typeof(U).parameters[1]}(U, MatrixMultiplicationGroupOperation())
end

function Base.show(io::IO, G::UnitaryGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "UnitaryGroup($(size))")
end
