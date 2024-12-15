"""
    OrthogonalGroup{T}

The orthogonal group ``$(_math(:O))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`).

# Constructor
    OrthogonalGroup(n; kwargs...)

Generate  orthogonal group ``$(_math(:O))(n)``.
All keyword arguments in `kwargs...` are passed on to [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) as well.
"""
const OrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ,MatrixMultiplicationGroupOperation,Manifolds.OrthogonalMatrices{T}
}

function OrthogonalGroup(n; kwargs...)
    R = Manifolds.OrthogonalMatrices(n; kwargs...)
    return OrthogonalGroup{typeof(R).parameters[1]}(R, MatrixMultiplicationGroupOperation())
end

function Base.show(io::IO, G::OrthogonalGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "OrthogonalGroup($(size))")
end
