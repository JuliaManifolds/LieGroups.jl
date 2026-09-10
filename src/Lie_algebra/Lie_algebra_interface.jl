"""
    LieAlgebra{𝔽, G} <: AbstractManifold{𝔽}

Represent the Lie algebra ``$(_math(:𝔤))``, that is a ``𝔽`` vector space with an associated
[`lie_bracket`](@ref) ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The Lie algebras considered here are those related to a [`AbstractLieGroup`](@ref) ``$(_math(:G))``,
namely the tangent space ``T_{$(_math(:e))}$(_math(:G))`` at the [`Identity`](@ref),
this is internally just a `const` of the corresponding $(_link(:TangentSpace)).

!!! note "Convention for representing tangent vectors in the Lie algebra"
    A vector field ``$(_tex(:Cal, "X")): $(_math(:G)) → T$(_math(:G))``, ``X(g) ∈ T_g$(_math(:G))``
    is called a left-invariant vector field if it satisfies

    ```math
    $(_tex(:Cal, "X"))(λ_g(h)) = Dλ_g(h)[$(_tex(:Cal, "X"))(h)], $(_tex(:quad))$(_tex(:text, "for all"))$(_tex(:quad)) g, h ∈ $(_math(:G)),
    ```

    where ``λ_g: $(_math(:G)) → $(_math(:G))`` is the left multiplication by ``g``.
    Hence ``$(_tex(:Cal, "X"))`` is determined already when ``X ∈ $(_math(:𝔤))`` is given,
    since ``$(_tex(:Cal, "X"))(g) = Dλ_g(e)[X]``, cf [HilgertNeeb:2012; Definition 9.1.7](@cite).

    Throughout `LieGroups.jl`, we use this left-invariant convention to store tangent vectors at points on a Lie group as elements of the corresponding Lie algebra.

# Constructor

    LieAlgebra(G::AbstractLieGroup)

Return the Lie Algebra belonging to the [`AbstractLieGroup`](@ref) `G`.
"""
const LieAlgebra{𝔽, O <: AbstractGroupOperation, G <: AbstractLieGroup{<:Any, O}} = ManifoldsBase.Fiber{
    𝔽, ManifoldsBase.TangentSpaceType, G, Identity{O},
}

function LieAlgebra(G::AbstractLieGroup{𝔽, O}) where {𝔽, O <: AbstractGroupOperation}
    return LieAlgebra{𝔽, O, typeof(G)}(G, Identity(G), ManifoldsBase.TangentSpaceType())
end

"""
    base_manifold(𝔤::LieAlgebra)

Return the [`base_manifold`](@extref `ManifoldsBase.base_manifold`) the
[`AbstractLieGroup`](@ref) of the given [`LieAlgebra`](@ref) is based on.
"""
ManifoldsBase.base_manifold(𝔤::LieAlgebra) = base_manifold(base_lie_group(𝔤))

"""
    base_lie_group(𝔤::LieAlgebra)

Return the [`base_lie_group`](@ref) of the given [`LieAlgebra`](@ref) belongs to.
"""
base_lie_group(𝔤::LieAlgebra) = 𝔤.manifold

_doc_get_coordinates = """
    get_coordinates(𝔤::LieAlgebra, X::T, B::AbstractBasis)
    get_coordinates!(𝔤::LieAlgebra, c, X::T, B::AbstractBasis)

Return the vector of coordinates to the decomposition of `X` with respect to an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`)
of the [`LieAlgebra`](@ref) `𝔤`.
The operation can be performed in-place of `c`.

By default this function requires that [`identity_element`](@ref)`(G, T)` is available and calls
the corresponding [`get_coordinates`](@extref ManifoldsBase :jl:function:`ManifoldsBase.get_coordinates`) function
of the Riemannian manifold the Lie group is build on.

The inverse operation is [`get_vector`](@ref).

See also [`vee`](@ref).
"""

@doc "$(_doc_get_coordinates)"
function ManifoldsBase.get_coordinates(
        𝔤::LieAlgebra, X, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
    )
    return ManifoldsBase._get_coordinates(𝔤, X, B)
end
# Mimic the levels from ManifoldsBase just without the base point p
function ManifoldsBase._get_coordinates(𝔤::LieAlgebra, X::T, B::AbstractBasis) where {T}
    G = base_lie_group(𝔤)
    return get_coordinates(base_manifold(G), identity_element(G, T), X, B)
end
@doc "$(_doc_get_coordinates)"
function ManifoldsBase.get_coordinates!(
        𝔤::LieAlgebra, c, X, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
    )
    return ManifoldsBase._get_coordinates!(𝔤, c, X, B)
end
function ManifoldsBase._get_coordinates!(𝔤::LieAlgebra, c, X::T, B::AbstractBasis) where {T}
    G = 𝔤.manifold
    return ManifoldsBase.get_coordinates!(base_manifold(G), c, identity_element(G, T), X, B)
end
function ManifoldsBase._get_coordinates(
        𝔤::LieAlgebra, X, B::DefaultLieAlgebraOrthogonalBasis
    )
    return get_coordinates_lie(𝔤, X, B)
end
function ManifoldsBase._get_coordinates!(
        𝔤::LieAlgebra, c, X, B::DefaultLieAlgebraOrthogonalBasis
    )
    get_coordinates_lie!(𝔤, c, X, B)
    return c
end

# the hat/vee variant
function get_coordinates_lie(𝔤::LieAlgebra, X, B::DefaultLieAlgebraOrthogonalBasis)
    c = allocate_result(𝔤, get_coordinates, X, B)
    return get_coordinates_lie!(𝔤, c, X, B)
end
function get_coordinates_lie!(
        𝔤::LieAlgebra, c, X::T, ::DefaultLieAlgebraOrthogonalBasis{𝔽}
    ) where {T, 𝔽}
    # Provide a default fallback
    G = base_lie_group(𝔤)
    return get_coordinates!(
        base_manifold(G), c, identity_element(G, T), X, DefaultOrthogonalBasis(𝔽)
    )
end

_doc_get_vector = """
    get_vector(G::AbstractLieGroup, c, B::AbstractBasis; kwargs...)
    get_vector(𝔤::LieAlgebra, c, B::AbstractBasis; kwargs...)
    get_vector!(G::AbstractLieGroup, X::T, c, B::AbstractBasis; kwargs...)
    get_vector!(𝔤::LieAlgebra, X::T, c, B::AbstractBasis; kwargs...)

Return the vector corresponding to a set of coefficients in an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`)
of the [`LieAlgebra`](@ref) `𝔤`.
Since all tangent vectors are assumed to be represented in the Lie algebra,
both signatures are equivalent.
The operation can be performed in-place of a tangent vector `X` of type `::T`.

By default this function requires [`identity_element`](@ref)`(G)` and calls
the corresponding [`get_vector`](@extref ManifoldsBase :jl:function:`ManifoldsBase.get_vectors`) function
of the Riemannian manifold the Lie group is build on.

The inverse operation is [`get_coordinates`](@ref).

# Keyword arguments

* `tangent_vector_type` specify the tangent vector type to use for the allocating variants.

See also [`hat`](@ref)
"""

@doc "$(_doc_get_vector)"
function ManifoldsBase.get_vector(
        𝔤::LieAlgebra,
        c,
        B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis();
        tangent_vector_type = nothing,
        kwargs...,
    )
    return ManifoldsBase._get_vector(𝔤, c, B, tangent_vector_type)
end
# Overwrite layer 2 since we do not have a base point and as well if a basis is provided and if we get nothing
# (define for all basis when moving this to Base)
@inline function ManifoldsBase._get_vector(
        𝔤::LieAlgebra, c, B::DefaultLieAlgebraOrthogonalBasis, ::Nothing
    )
    return get_vector_lie(𝔤, c, B)
end
@inline function ManifoldsBase._get_vector(
        𝔤::LieAlgebra, c, B::DefaultLieAlgebraOrthogonalBasis, T::Type
    )
    return get_vector_lie(𝔤, c, B, T)
end
@inline function ManifoldsBase._get_vector(𝔤::LieAlgebra, c, B::AbstractBasis, ::Nothing)
    G = base_lie_group(𝔤)
    return get_vector(base_manifold(G), identity_element(G), c, B)
end
@inline function ManifoldsBase._get_vector(𝔤::LieAlgebra, c, B::AbstractBasis, T::Type)
    G = base_lie_group(𝔤)
    return get_vector(base_manifold(G), identity_element(G, T), c, B)
end

@doc "$(_doc_get_vector)"
function ManifoldsBase.get_vector!(
        𝔤::LieAlgebra, X, c, B::ManifoldsBase.AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
    )
    return ManifoldsBase._get_vector!(𝔤, X, c, B)
end

function ManifoldsBase._get_vector!(
        𝔤::LieAlgebra, X::T, c, B::DefaultLieAlgebraOrthogonalBasis
    ) where {T}
    return get_vector_lie!(𝔤, X, c, B)
end
function ManifoldsBase._get_vector!(
        𝔤::LieAlgebra, X::T, c, B::ManifoldsBase.AbstractBasis
    ) where {T}
    G = base_lie_group(𝔤)
    return ManifoldsBase.get_vector!(base_manifold(G), X, identity_element(G, T), c, B)
end

@inline function get_vector_lie(𝔤::LieAlgebra, c, B::DefaultLieAlgebraOrthogonalBasis)
    X = zero_vector(𝔤)
    return get_vector_lie!(𝔤, X, c, B)
end
@inline function get_vector_lie(
        𝔤::LieAlgebra, c, B::DefaultLieAlgebraOrthogonalBasis, T::Type
    )
    X = zero_vector(𝔤, T)
    return get_vector_lie!(𝔤, X, c, B)
end
@inline function get_vector_lie!(
        𝔤::LieAlgebra, X::T, c, ::DefaultLieAlgebraOrthogonalBasis{𝔽}
    ) where {T, 𝔽}
    G = base_lie_group(𝔤)
    return get_vector!(
        base_manifold(G), X, identity_element(G, T), c, DefaultOrthogonalBasis(𝔽)
    )
end

_doc_hat = """
    hat(G::LieAlgebra, c)
    hat(G::LieAlgebra, c, T::Type)
    hat!(G::LieAlgebra, X::T, c)

Compute the hat map ``(⋅)^{\\wedge}: $(_tex(:Cal, "V")) → 𝔤`` that maps a vector of coordinates ``$(_tex(:vec, "c")) ∈ $(_tex(:Cal, "V"))``,
to a tangent vector ``X ∈ $(_math(:𝔤))``.
The coefficients are given with respect to a specific basis to a tangent vector in the Lie algebra

```math
X = $(_tex(:sum))_{i∈$(_tex(:Cal, "I"))} c_iB_i,
```

where ``$(_tex(:Set, "B_i"))_{i∈$(_tex(:Cal, "I"))}`` is a basis of the Lie algebra
and ``$(_tex(:Cal, "I"))`` a corresponding index set, which is usually ``$(_tex(:Cal, "I"))=$(_tex(:Set, raw"1,\ldots,n"))``.
Then ``$(_tex(:Cal, "V")) = ℝ^n``.

For the allocating variant, you can specify the type `T` of the tangent vector to obtain,
in case there are different representations. The first signature produces the default representation.

The computation can be performed in-place of `X`. The inverse of `hat` is [`vee`](@ref).
Technically, `hat` is a specific case of [`get_vector`](@ref) and is implemented using the
[`DefaultLieAlgebraOrthogonalBasis`](@ref).
"""

@doc "$(_doc_hat)"
function ManifoldsBase.hat(𝔤::LieAlgebra, c)
    return get_vector(𝔤, c, DefaultLieAlgebraOrthogonalBasis())
end

@doc "$(_doc_hat)"
function ManifoldsBase.hat(𝔤::LieAlgebra, c, T::Type)
    return get_vector(𝔤, c, DefaultLieAlgebraOrthogonalBasis(); tangent_vector_type = T)
end

@doc "$(_doc_hat)"
function ManifoldsBase.hat!(𝔤::LieAlgebra, X, c)
    get_vector!(𝔤, X, c, DefaultLieAlgebraOrthogonalBasis())
    return X
end

"""
    inner(𝔤::LieAlgebra, X, Y)

Compute the inner product ``⟨⋅,⋅⟩: $(_math(:𝔤))×$(_math(:𝔤)) → ℝ`` on the Lie algebra ``𝔤``.
By default this uses the inner product on the underlying manifold of the
[`AbstractLieGroup`](@ref) of `𝔤` at the [`identity_element`](@ref)`(G)`.
Note that this method allocates an appropriate identity element.
"""
ManifoldsBase.inner(𝔤::LieAlgebra, X, Y)

"""
    is_point(𝔤::LieAlgebra, X; kwargs...)

Check whether `X` is a valid point on the Lie Algebra `𝔤`.
This falls back to checking whether `X` is a valid point on the tangent space
at the [`identity_element`](@ref)`(G)` on the [`base_manifold`](@ref)`(G)`
on the [`AbstractLieGroup`](@ref) of `𝔤`
"""
function ManifoldsBase.is_point(𝔤::LieAlgebra, X::T, b::Bool = false; kwargs...) where {T}
    return ManifoldsBase.is_vector(
        base_lie_group(𝔤), identity_element(base_lie_group(𝔤), T), X, b; kwargs...
    )
end

_doc_lie_bracket = """
    lie_bracket(𝔤::LieAlgebra, X, Y)
    lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The computation can be done in-place of `Z`.
"""
function lie_bracket end
@doc "$(_doc_lie_bracket)"
function lie_bracket(𝔤::LieAlgebra, X, Y)
    Z = ManifoldsBase.allocate_result(𝔤, lie_bracket, X, Y)
    return lie_bracket!(𝔤, Z, X, Y)
end

function lie_bracket! end
@doc "$(_doc_lie_bracket)"
lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

function LinearAlgebra.norm(𝔤::LieAlgebra, X)
    return sqrt(real(inner(𝔤, X, X)))
end
# Avoid an ambiguity
function LinearAlgebra.norm(𝔤::LA, X::Real) where {LA <: LieAlgebra}
    return LinearAlgebra.norm(base_manifold(𝔤), identity_element(base_lie_group(𝔤)), X)
end

function ManifoldsBase.project!(𝔤::LieAlgebra, Y, X)
    return ManifoldsBase.project!(
        base_manifold(𝔤), Y, identity_element(base_lie_group(𝔤)), X
    )
end
function ManifoldsBase.project!(𝔤::LieAlgebra, W, X, V)
    return ManifoldsBase.project!(
        base_manifold(𝔤), W, identity_element(base_lie_group(𝔤)), V
    )
end

_doc_rand_algebra = """
    rand(::LieGroup; vector_at=nothing, σ=1.0, kwargs...)
    rand(::LieAlgebra; σ=1.0, kwargs...)
    rand!(::AbstractLieGroup, gX; vector_at=nothing, kwargs...)
    rand!(::LieAlgebra, X; σ=1.0, kwargs...)

Compute a random point or tangent vector on a Lie group.

For points this just means to generate a random point on the
underlying manifold itself.

For tangent vectors, an element in the Lie Algebra is generated,
see also [`rand(::LieAlgebra; kwargs...)`](@ref)
"""

@doc "$(_doc_rand_algebra)"
Random.rand(::LieAlgebra; kwargs...)

function Random.rand(𝔤::LieAlgebra, T::Type; vector_at = nothing, kwargs...)
    X = allocate_on(𝔤, TangentSpaceType(), T)
    G = base_lie_group(𝔤)
    rand!(𝔤, X; vector_at = identity_element(G), kwargs...)
    return X
end
function Random.rand(𝔤::LieAlgebra, d::Integer; vectpr_at = nothing, kwargs...)
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return [rand(M; vector_at = identity_element(G), kwargs...) for _ in 1:d]
end
function Random.rand(rng::AbstractRNG, 𝔤::LieAlgebra, T::Type; vector_at = nothing, kwargs...)
    X = allocate_on(base_lie_group(𝔤), TangentSpaceType(), T)
    G = base_lie_group(𝔤)
    # Here we also have to turn T into a point type P for the identity.
    rand!(rng, 𝔤, X; vector_at = identity_element(G), kwargs...)
    return X
end

@doc "$(_doc_rand_algebra)"
Random.rand!(::LieAlgebra, X; kwargs...)

function Base.show(io::IO, 𝔤::LieAlgebra)
    return print(io, "LieAlgebra($(base_lie_group(𝔤)))")
end

# Overwrite the longer version for tangent spaces
function Base.show(io::IO, ::MIME"text/plain", 𝔤::LieAlgebra)
    return print(io, "The Lie algebra of the Lie Group $(base_lie_group(𝔤))")
end

_doc_vee = """
    vee(𝔤::LieAlgebra, X)
    vee!(𝔤::LieAlgebra, c, X)

Compute the vee map ``(⋅){\\vee}: $(_math(:𝔤)) → $(_tex(:Cal, "V"))`` that maps a tangent vector `X`
from the [`LieAlgebra`](@ref) ``$(_math(:𝔤))`` to its coordinates with respect to the [`DefaultLieAlgebraOrthogonalBasis`](@ref) basis in the Lie algebra

```math
X = $(_tex(:sum))_{i∈$(_tex(:Cal, "I"))} c_iB_i,
```

where ``$(_tex(:Set, "B_i"))_{i∈$(_tex(:Cal, "I"))}`` is a basis of the Lie algebra
and ``$(_tex(:Cal, "I"))`` a corresponding index set, which is usually ``$(_tex(:Cal, "I"))=$(_tex(:Set, raw"1,\ldots,n"))``.
Then ``$(_tex(:Cal, "V")) = ℝ^n``

The computation can be performed in-place of `c`. The inverse of `vee` is [`hat`](@ref).
Technically, `vee` is a specific case of [`get_coordinates`](@ref) and is implemented using
the [`DefaultLieAlgebraOrthogonalBasis`](@ref).
"""

@doc "$(_doc_vee)"
function ManifoldsBase.vee(𝔤::LieAlgebra, X)
    return get_coordinates(𝔤, X, DefaultLieAlgebraOrthogonalBasis())
end

@doc "$(_doc_vee)"
function ManifoldsBase.vee!(𝔤::LieAlgebra, c, X)
    get_coordinates!(𝔤, c, X, DefaultLieAlgebraOrthogonalBasis())
    return c
end

"""
    zero_vector(𝔤::LieAlgebra)
    zero_vector(𝔤::LieAlgebra, T::Type)
    zero_vector!(𝔤::LieAlgebra, X::T)

Generate a $(_link(:zero_vector)) of type `T` in the [`LieAlgebra`](@ref) ``𝔤`` of
the [`AbstractLieGroup`](@ref) `G`.
By default this calls `zero_vector` on the manifold of `G` at the `identity_element(G,T)`

For the allocating variant the type `T` of the zero vector can be specified.
"""
ManifoldsBase.zero_vector(G::LieGroup{𝔽, <:O}, T::Type) where {𝔽, O <: AbstractGroupOperation}

function ManifoldsBase.zero_vector(𝔤::LieAlgebra, T::Type)
    G = base_lie_group(𝔤) # access manifold twice -> pass to manifold directly
    return ManifoldsBase.zero_vector(
        base_manifold(G), identity_element(G, point_type(G, T))
    )
end
function ManifoldsBase.zero_vector(𝔤::LieAlgebra)
    G = base_lie_group(𝔤) # access manifold twice -> pass to manifold directly
    return ManifoldsBase.zero_vector(base_manifold(G), identity_element(G))
end
function ManifoldsBase.zero_vector!(𝔤::LieAlgebra, X::T) where {T}
    G = base_lie_group(𝔤) # access manifold twice -> pass to manifold directly
    return ManifoldsBase.zero_vector!(base_manifold(G), X, identity_element(G, T))
end

#
#
# A fallback macro for types that merely wrap the actual data
"""
    default_lie_algebra_fallbacks(TG, TF, Op, TV, Xfield::Symbol)

Introduce default fallbacks for all basic functions on Lie algebras, for Lie group of type
`TG` with number system `TF`, an group operation `Op`, tangent vectors of type `TV`, with
forwarding to fields `Xfield` and tangent vector functions
"""
macro default_lie_algebra_fallbacks(TG, TF, Op, TV, Xfield::Symbol)
    block = quote
        function ManifoldsBase.check_size(𝔤::LieAlgebra{$TF, <:$Op, <:$TG}, X::$TV; kwargs...)
            return ManifoldsBase.check_size(𝔤, X.$Xfield; kwargs...)
        end
        function ManifoldsBase.is_point(𝔤::LieAlgebra{$TF, <:$Op, <:$TG}, X::$TV; kwargs...)
            return ManifoldsBase.is_point(𝔤, X.$Xfield; kwargs...)
        end
        function LinearAlgebra.norm(𝔤::LieAlgebra{$TF, <:$Op, <:$TG}, X::$TV)
            return LinearAlgebra.norm(𝔤, X.$Xfield)
        end
        function ManifoldsBase.isapprox(
                𝔤::LieAlgebra{$TF, <:$Op, <:$TG}, X::$TV, Y::$TV; kwargs...
            )
            return ManifoldsBase.isapprox(𝔤, X.$Xfield, Y.$Xfield; kwargs...)
        end
    end
    return esc(block)
end

#
#
# allocation helpers

function ManifoldsBase.allocate_result(
        𝔤::LieAlgebra, f::typeof(ManifoldsBase.get_coordinates), X, basis::AbstractBasis{𝔽}
    ) where {𝔽}
    X_ = internal_value(X)
    T = ManifoldsBase.coordinate_eltype(𝔤, X_, 𝔽)
    return ManifoldsBase.allocate_coordinates(𝔤, X_, T, number_of_coordinates(𝔤, basis))
end
