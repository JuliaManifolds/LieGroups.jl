"""
    TranslationGroup{𝔽,T}

The Lie group consisting of the [`AdditionGroupOperation`](@ref) on some
[`Euclidean`](@extref `Manifolds.Euclidean`) space.

# Constructor
    TranslationGroup(n₁,...,nᵢ; kwargs...)

Generate the translation group on ``𝔽^{n₁,…,nᵢ}`` = `Euclidean(n₁,...,nᵢ; field=𝔽)`,
which is isomorphic to the group itself. all keyword arguments in `kwargs...`
are passed on to [`Euclidean`](@extref `Manifolds.Euclidean`) as well
"""
const TranslationGroup{𝔽,T} = LieGroup{𝔽,AdditionGroupOperation,Manifolds.Euclidean{T,𝔽}}

function TranslationGroup(n::Int...; kwargs...)
    Rn = Manifolds.Euclidean(n...; kwargs...)
    return TranslationGroup{typeof(Rn).parameters[[2, 1]]...}(Rn, AdditionGroupOperation())
end

function Base.show(io::IO, G::TranslationGroup)
    return print(io, "TranslationGroup($(G.manifold), $(G.op))")
end
