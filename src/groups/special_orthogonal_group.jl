"""
    SpecialOrthogonalGroup{T}

The special orthogonal group ``$(_math(:SO))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`Rotations`](@extref `Manifolds.Rotations`).

# Constructor
    SpecialOrthogonalGroup(n; kwargs...)

Generate  special orthogonal group ``$(_math(:SO))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialOrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ,MatrixMultiplicationGroupOperation,Manifolds.Rotations{T}
}

function SpecialOrthogonalGroup(n; kwargs...)
    R = Manifolds.Rotations(n...; kwargs...)
    return SpecialOrthogonalGroup{typeof(R).parameters[1]}(
        R, MatrixMultiplicationGroupOperation()
    )
end

function Base.show(io::IO, G::SpecialOrthogonalGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "SpecialOrthogonalGroup($(size))")
end
