@doc """
    AdditiveGroup <: LieGroup{bbF   }

The group ``(𝕂^n, +)`` is a [`LieGroup`](@ref) using the [`AdditiveGroupOperation`](@ref)
and serves as the prototype for the [`AdditiveGroupOperation`](@ref).

# Constructor

    AdditiveGroup(n₁,n₂,...,nᵢ; field=ℝ, kwargs...)

Generate the additive group. All parameters are used to initialise the internal
[`Euclidean`](@extref `Manifolds.Euclidean`) manifold.
"""
const AdditiveGroup{𝔽} = LieGroup{𝔽,Euclidean{T,𝔽} where T,AdditiveGroupOperation} where {𝔽}

function AdditiveGroup(n::Vararg{Int,I}; field::AbstractNumbers=ℝ, kwargs...) where {I}
    return AdditiveGroup(Euclidean(n; field=field, kwargs...), AdditiveGroupOperation())
end
