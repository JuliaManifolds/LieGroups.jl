
@doc """
    GeneralLinearGroup{𝔽,T}

The general linear group ``\operatorname{GL}(n)`` is the set of all invertible matrices in ``𝔽^{n×n}``
equipped with the [`MatrixMultiplicationGroupOperation`](@ref) as the group operation.

# Constructor

    GeneralLinearGroup(n::Int; kwargs...)

Generate the general linear group  group on ``𝔽^{n×n}``.
All keyword arguments in `kwargs...` are passed on to [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`).
"""
const GeneralLinearGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,Manifolds.InvertibleMatrices{𝔽,T}
}

function GeneralLinearGroup(n::Int...; kwargs...)
    Im = Manifolds.InvertibleMatrices(n...; kwargs...)
    return TranslationGroup{typeof(Im).parameters[[2, 1]]...}(
        Im, MatrixMultiplicationGroupOperation()
    )
end

function Base.show(io::IO, G::GeneralLinearGroup{𝔽}) where {𝔽}
    n = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "GeneralLinearGroup($n; field=$(𝔽))")
end
