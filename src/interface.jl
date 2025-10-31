#
#
# Generic types for the interface of a Lie group
@doc """
    AbstractGroupOperation

Represent a type of group operation for a [`AbstractLieGroup`](@ref) ``$(_math(:G))``, that is a
smooth binary operation ``$(_math(:∘)) : $(_math(:G)) × $(_math(:G)) → $(_math(:G))``
on elements of a Lie group ``$(_math(:G))``.
"""
abstract type AbstractGroupOperation end

@doc raw"""
    BaseManifoldInverseRetraction{IRM<:AbstractInverseRetractionMethod} <: AbstractInverseRetractionMethod

Compute an inverse retraction by using the inverse retraction of type `IRM` on the base manifold of
a [`LieGroup`](@ref).

# Constructor

    BaseManifoldInverseRetraction(irm::AbstractInverseRetractionMethod)

Generate the inverse retraction with inverse retraction `rm` to use on the base manifold.
"""
struct BaseManifoldInverseRetraction{IRM <: AbstractInverseRetractionMethod} <:
    AbstractInverseRetractionMethod
    inverse_retraction::IRM
end

@doc raw"""
    BaseManifoldRetraction{RM<:AbstractRetractionMethod} <: AbstractRetractionMethod

Compute a retraction by using the retraction of type `RM` on the base manifold of
a [`LieGroup`](@ref).

# Constructor

    BaseManifoldRetraction(rm::AbstractRetractionMethod)

Generate the retraction with retraction `rm` to use on the base manifold.
"""
struct BaseManifoldRetraction{RM <: AbstractRetractionMethod} <: AbstractRetractionMethod
    retraction::RM
end

"""
    DefaultLieAlgebraOrthogonalBasis{𝔽} <: ManifoldsBase.AbstractOrthogonalBasis{𝔽,ManifoldsBase.TangentSpaceType}

Specify an orthogonal basis for a Lie algebra.
This is used as the default within [`hat`](@ref) and [`vee`](@ref).

If not specifically overwritten/implemented for a Lie group, the [`DefaultOrthogonalBasis`](@extref `ManifoldsBase.DefaultOrthogonalBasis`)
at the [`identity_element`](@ref) on the [`base_manifold](@ref base_manifold(::AbstractLieGroup)) acts as a fallback.

!!! note
    In order to implement the corresponding [`get_coordinates`](@ref) and [`get_vector`](@ref) functions,
    define `get_coordinates_lie(::AbstractLieAlgebra, X, B)` and `get_vector_lie(::AbstractLieAlgebra, X, B)`, resp.
"""
struct DefaultLieAlgebraOrthogonalBasis{𝔽} <:
    ManifoldsBase.AbstractOrthogonalBasis{𝔽, ManifoldsBase.TangentSpaceType} end
function DefaultLieAlgebraOrthogonalBasis(𝔽::ManifoldsBase.AbstractNumbers = ℝ)
    return DefaultLieAlgebraOrthogonalBasis{𝔽}()
end

"""
    AbstractLieGroup{𝔽, O<:AbstractGroupOperation, M<:AbstractManifold{𝔽}} <: AbstractManifold{𝔽}

An abstract type to represent Lie groups. For most cases it should suffice to “combine”
an $(_link(:AbstractManifold)) with an [`AbstractGroupOperation`](@ref), see [`LieGroup`](@ref).
"""
abstract type AbstractLieGroup{𝔽, O <: AbstractGroupOperation, M <: AbstractManifold{𝔽}} <:
AbstractManifold{𝔽} end

"""
    LieGroup{𝔽, O<:AbstractGroupOperation, M<:AbstractManifold{𝔽}} <:  AbstractLieGroup{𝔽, O, M}

Represent a Lie Group ``$(_math(:G))``.

A *Lie Group* ``$(_math(:G))`` is a group endowed with the structure of a manifold such that the
group operations ``$(_math(:∘)): $(_math(:G))×$(_math(:G)) → $(_math(:G))``, see [`compose`](@ref)
and the inverse operation ``⋅^{-1}: $(_math(:G)) → $(_math(:G))``, see [`inv`](@ref) are smooth,
see for example [HilgertNeeb:2012; Definition 9.1.1](@cite).

Lie groups are named after the Norwegian mathematician [Marius Sophus Lie](https://en.wikipedia.org/wiki/Sophus_Lie) (1842–1899).

# Fields

* `manifold`: an $(_link(:AbstractManifold)) ``$(_math(:M))``
* `op`: an [`AbstractGroupOperation`](@ref) ``$(_math(:∘))`` on that manifold

# Constructor

    LieGroup(M::AbstractManifold, op::AbstractGroupOperation)

Generate a Lie group based on a manifold `M` and a group operation `op`, where vectors by default are stored in the Lie Algebra.
"""
struct LieGroup{𝔽, O <: AbstractGroupOperation, M <: ManifoldsBase.AbstractManifold{𝔽}} <:
    AbstractLieGroup{𝔽, O, M}
    manifold::M
    op::O
end

@doc """
    Identity{O<:AbstractGroupOperation}

Represent the group identity element ``e ∈ $(_math(:G))`` on an [`AbstractLieGroup`](@ref) ``$(_math(:G))``
with [`AbstractGroupOperation`](@ref) of type `O`.

Similar to the philosophy that points are agnostic of their group at hand, the identity
does not store the group ``$(_math(:G))`` it belongs to. However it depends on the type of the [`AbstractGroupOperation`](@ref) used.

See also [`identity_element`](@ref) on how to obtain the corresponding [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`) or array representation.

# Constructors

    Identity(::AbstractLieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation}
    Identity(o::AbstractGroupOperation)
    Identity(::Type{AbstractGroupOperation})

create the identity of the corresponding subtype `O<:`[`AbstractGroupOperation`](@ref)
"""
struct Identity{O <: AbstractGroupOperation} end

Identity(::AbstractLieGroup{𝔽, O}) where {𝔽, O <: AbstractGroupOperation} = Identity{O}()
Identity(::O) where {O <: AbstractGroupOperation} = Identity(O)
Identity(::Type{O}) where {O <: AbstractGroupOperation} = Identity{O}()

"""
    AbstractLieGroupPoint <: ManifoldsBase.AbstractManifoldPoint end

An abstract type for a point on an [`AbstractLieGroup`](@ref).
While an points and tangent vectors are usually kept untyped for flexibility,
it might be necessary to distinguish different types of points, for example

* for complicated representations that require a struct
* semantic verification
* when there exist different representations

By sub-typing the [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`), this follows the same idea as in $(_link(:ManifoldsBase)).
"""
abstract type AbstractLieGroupPoint <: ManifoldsBase.AbstractManifoldPoint end

"""
    AbstractLieAlgebraTangentVector <: ManifoldsBase.AbstractTangentVector

An abstract type for a tangent vector represented in a [`LieAlgebra`](@ref).

While an tangent vectors are usually kept untyped for flexibility,
it might be necessary to distinguish different types of points, for example

* for complicated representations that require a struct
@ semantic verification
* when there exist different representations

By sub-typing the [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`),
this follows the same idea as in $(_link(:ManifoldsBase)).
"""
abstract type AbstractLieAlgebraTangentVector <: ManifoldsBase.AbstractTangentVector end

#
#
# --- Functions ---

_doc_adjoint = """
    adjoint(G::AbstractLieGroup, g, X)
    adjoint!(G::AbstractLieGroup, Y, g, X)

Compute the adjoint ``$(_math(:Ad))(g): $(_math(:𝔤)) → $(_math(:𝔤))``, which is defined as
the differential [`diff_conjugate`](@ref) of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``
evaluated at the [`Identity`](@ref) ``h=$(_math(:e))``.
The operation can be performed in-place of `Y`.

```math
  $(_math(:Ad))(g)[X] = $(_math(:d)) c_g($(_math(:e))) [X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```

see [HilgertNeeb:2012; Section 9.2.3](@cite).

On matrix Lie groups the adjoint reads ``$(_math(:Ad))(g)[X] = gXg^{-1}``.
"""

@doc "$(_doc_adjoint)"
function Base.adjoint(G::AbstractLieGroup, g, X)
    Y = ManifoldsBase.allocate_result(G, adjoint, g, X)
    return adjoint!(G, Y, g, X)
end

function adjoint! end
@doc "$(_doc_adjoint)"
function adjoint!(G::AbstractLieGroup, Y, g, X)
    diff_conjugate!(G, Y, g, Identity(G), X)
    return Y
end

@doc """
    base_manifold(G::AbstractLieGroup)

Return the manifold stored within the [`AbstractLieGroup`](@ref) `G`.
"""
ManifoldsBase.base_manifold(G::AbstractLieGroup)
ManifoldsBase.base_manifold(G::LieGroup) = G.manifold

# Since we dispatch per point here, identity is already checked on the `is_point` level.
function ManifoldsBase.check_point(
        G::AbstractLieGroup{𝔽, O}, g; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return ManifoldsBase.check_point(base_manifold(G), g; kwargs...)
end

function ManifoldsBase.check_vector(G::AbstractLieGroup, g::P, X; kwargs...) where {P}
    return ManifoldsBase.check_vector(
        base_manifold(G), identity_element(G, P), X; kwargs...
    )
end

# compose g ∘ h
_doc_compose = """
    compose(G::AbstractLieGroup, g, h)
    compose!(G::AbstractLieGroup, k, g, h)

Perform the group operation ``g $(_math(:∘)) h`` for two ``g, h ∈ $(_math(:G))``
on the [`AbstractLieGroup`](@ref) `G`. This can also be done in-place of `h`.

!!! info
    This function also handles the case where `g` or/and `h` are the [`Identity`](@ref)`(G)`.
    Since this would lead to ambiguities when implementing a new group operations,
    this function calls `_compose` and `_compose!`, respectively, which is meant for the actual computation of
    group operations on (non-[`Identity`](@ref)` but maybe its numerical representation) elements.
"""
@doc "$(_doc_compose)"
compose(G::AbstractLieGroup, g, h) = _compose(G, g, h)
compose(::AbstractLieGroup{𝔽, O}, g::Identity{O}, h) where {𝔽, O <: AbstractGroupOperation} = h
compose(::AbstractLieGroup{𝔽, O}, g, h::Identity{O}) where {𝔽, O <: AbstractGroupOperation} = g
function compose(
        ::AbstractLieGroup{𝔽, O}, g::Identity{O}, h::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return g
end

function _compose(G::AbstractLieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose, g, h)
    return _compose!(G, k, g, h)
end

function compose! end

@doc "$(_doc_compose)"
compose!(G::AbstractLieGroup, k, g, h) = _compose!(G, k, g, h)
function compose!(
        G::AbstractLieGroup{𝔽, O}, k, ::Identity{O}, h
    ) where {𝔽, O <: AbstractGroupOperation}
    return copyto!(G, k, h)
end
function compose!(
        G::AbstractLieGroup{𝔽, O}, k, g, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return copyto!(G, k, g)
end
function compose!(
        G::AbstractLieGroup{𝔽, O}, k, ::Identity{O}, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return identity_element!(G, k)
end
function compose!(
        ::AbstractLieGroup{𝔽, O}, k::Identity{O}, ::Identity{O}, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return k
end

function _compose! end

_doc_conjugate = """
    conjugate(G::AbstractLieGroup, g, h)
    conjugate!(G::AbstractLieGroup, k, g, h)

Compute the conjugation map ``c_g: $(_math(:G)) → $(_math(:G))`` given by ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_conjugate)"
function conjugate(G::AbstractLieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, conjugate, h, g)
    return conjugate!(G, k, g, h)
end

function conjugate! end
@doc "$(_doc_conjugate)"
function conjugate!(G::AbstractLieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, k) # `h∘k` in-place of k
    compose!(G, k, g, k) # `g∘k` in-place of k
    return k
end

ManifoldsBase.copyto!(G::AbstractLieGroup, h, g) = copyto!(base_manifold(G), h, g)
function ManifoldsBase.copyto!(
        G::AbstractLieGroup{𝔽, O}, h::P, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation, P}
    return identity_element!(G, h)
end
function ManifoldsBase.copyto!(
        ::AbstractLieGroup{𝔽, O}, h::Identity{O}, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return h
end
function ManifoldsBase.copyto!(
        G::AbstractLieGroup{𝔽, O}, h::Identity{O}, g
    ) where {𝔽, O <: AbstractGroupOperation}
    (is_identity(G, g)) && return h
    throw(
        DomainError(
            g,
            "copyto! into the identity element of $G ($h) is not defined for a non-identity element g ($g)",
        ),
    )
end

function ManifoldsBase.default_basis(
        ::AbstractLieGroup, ::Type{T}; field::AbstractNumbers = ℝ
    ) where {T}
    return DefaultLieAlgebraOrthogonalBasis(field)
end
function ManifoldsBase.default_basis(::AbstractLieGroup; field::AbstractNumbers = ℝ)
    return DefaultLieAlgebraOrthogonalBasis(field)
end

_doc_diff_conjugate = """
    diff_conjugate(G::AbstractLieGroup, g, h, X)
    diff_conjugate!(G::AbstractLieGroup, Y, g, h, X)

Compute the differential of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``
on the [`AbstractLieGroup`](@ref) `G`. The operation can be performed in-place of `Y`.

```math
  $(_math(:d))(c_g(h))[X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```
"""
@doc "$(_doc_diff_conjugate)"
function diff_conjugate(G::AbstractLieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_conjugate, g, h, X)
    return diff_conjugate!(G, Y, g, h, X)
end

function diff_conjugate! end
@doc "$(_doc_diff_conjugate)"
diff_conjugate!(::AbstractLieGroup, Y, g, h, X)

_doc_diff_inv = """
    diff_inv(G::AbstractLieGroup, g, X)
    diff_inv!(G::AbstractLieGroup, Y, g, X)

Compute the differential of the function ``ι_{$(_math(:G))}(g) = g^{-1}``, where
``$(_math(:d))ι_{$(_math(:G))}(g): $(_math(:𝔤)) → $(_math(:𝔤))``.
This can be done in-place of `Y`.
Note that we represent tangent vectors in the Lie algebra ``𝔤``.

For example on matrix manifolds this means, we use ``X ∈ 𝔤`` and hence ``W = gX ∈ T_g$(_math(:G))``.
The (classical) differential ``$(_math(:D))ι_{$(_math(:G))}(g): T_g$(_math(:G)) → T_{g^{-1}}$(_math(:G))`` reads

```math
  $(_math(:D))ι_{$(_math(:G))}(g)[W] = -g^{-1}Wg^{-1} = -Xg^{-1} = -g^{-1}(gXg^{-1}) = -g^{-1}$(_math(:Ad))(g)[X] = V ∈ T_{g^{-1}}$(_math(:G)),
```

see e.g. [Giles:2008](@cite). To bring this back to the Lie algebra, we Write ``V = g^{-1}Y ∈ T_{g^{-1}}$(_math(:G))``
for some ``Y ∈ 𝔤`` and obtain

```math
  $(_math(:d)) ι_{$(_math(:G))}(g)[X] = -$(_math(:Ad))(g)[X] ∈ 𝔤,
```

where we use ``$(_math(:d))`` to denote the differential in the Lie algebra.
"""

@doc "$_doc_diff_inv"
function diff_inv(G::AbstractLieGroup, g, X)
    Y = allocate_result(G, diff_inv, g, X)
    return diff_inv!(G, Y, g, X)
end

function diff_inv! end
@doc "$_doc_diff_inv"
diff_inv!(G::AbstractLieGroup, Y, g, X)

_doc_diff_left_compose = """
    diff_left_compose(G::AbstractLieGroup, g, h, X)
    diff_left_compose!(G::AbstractLieGroup, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its first (left) argument `g`.

Another interpretation is to consider a function where we do a fixed multiplication from the right with `h`.
i..e. the right group multiplication function ``ρ_h(g) = g$(_math(:∘))h`` (where the _right_ refers to the fixed argument ``h``).

In this notation, this function computes the differential ``$(_math(:d))ρ_h: 𝔤 → 𝔤``.

For example on matrix Lie groups this means that for ``X ∈ 𝔤`` we can start with ``W = gX ∈ T_g$(_math(:G))``
and compute the (classical) differential ``$(_math(:D))ρ_h(g): T_g$(_math(:G)) → T_{g$(_math(:∘))h}$(_math(:G))``.

It reads

```math
  $(_math(:D))ρ_h(g)[W] = Wh = gXh = V ∈ T_{g$(_math(:∘))h}$(_math(:G)).
```

To obtain the Lie group differential ``$(_math(:d))ρ_h(g)`` we have to “pull back” ``V``
from the tangent space ``T_{g$(_math(:∘))h}$(_math(:G))`` to the Lie algebra ``𝔤``.
We use the same identification, that we can write ``V = ghY ∈ T_{g$(_math(:∘))h}$(_math(:G))``.
This means in practice that with ``V = gXh = gh(h^{-1}Xh)`` differential reads

```math
$(_math(:d)) ρ_h(g)[X] = h^{-1}Xh = $(_math(:Ad))(h^{-1})[X] ∈ 𝔤,
```

where ``$(_math(:Ad))`` denotes the [`adjoint`](@ref).
"""
@doc "$(_doc_diff_left_compose)"
function diff_left_compose(G::AbstractLieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_left_compose, X, g, h)
    return diff_left_compose!(G, Y, g, h, X)
end
function diff_left_compose(G::AbstractLieGroup{𝔽, O}, g::Identity{O}, h, X) where {𝔽, O <: AbstractGroupOperation}
    Y = ManifoldsBase.allocate_result(G, diff_left_compose, X, h)
    return diff_left_compose!(G, Y, g, h, X)
end

function diff_left_compose! end
@doc "$(_doc_diff_left_compose)"
diff_left_compose!(::AbstractLieGroup, Y, g, h, X)

_doc_diff_right_compose = """
    diff_right_compose(G::AbstractLieGroup, g, h, X)
    diff_right_compose!(G::AbstractLieGroup, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its second (right) argument `h`.

Another interpretation is to consider a function where we do a fixed multiplication from the left with `g`.
i..e. the left group multiplication function ``λ_g(h) = g$(_math(:∘))h`` (where the _left_ refers to the fixed argument ``g``.).

In this notation, this function ``$(_math(:d))λ_g: 𝔤 → 𝔤``.

For example on matrix Lie groups this means that for ``X ∈ 𝔤`` we can start with ``W = hX ∈ T_h$(_math(:G))``
and compute the (classical) differential ``$(_math(:D))λ_g(h): T_h$(_math(:G)) → T_{g$(_math(:∘))h}$(_math(:G))``.

It reads

```math
  $(_math(:D))λ_g(h)[W] = gW = ghX ∈ T_{g$(_math(:∘))h}$(_math(:G)).
```

To obtain the Lie group differential ``$(_math(:d))λ_g(h)`` we have to multiply the result with ``(gh)^{-1}``
from the left and move from ``W`` to ``X``. Then the differential just simplifies to the identity. It reads

```math
$(_math(:d)) λ_g(h)[X] = X ∈ 𝔤.
```
"""
@doc "$(_doc_diff_right_compose)"
function diff_right_compose(G::AbstractLieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_right_compose, X, g, h)
    return diff_right_compose!(G, Y, g, h, X)
end

function diff_right_compose(G::AbstractLieGroup{𝔽, O}, g::Identity{O}, h, X) where {𝔽, O <: AbstractGroupOperation}
    Y = ManifoldsBase.allocate_result(G, diff_right_compose, X, h)
    return diff_right_compose!(G, Y, g, h, X)
end

function diff_right_compose! end
@doc "$(_doc_diff_right_compose)"
diff_right_compose!(::AbstractLieGroup, Y, g, h, X)

_doc_exp = """
    exp(G::AbstractLieGroup, g, X)
    exp!(G::AbstractLieGroup, h, g, X)

Compute the Lie group exponential map for ``g∈$(_math(:G))`` and ``X∈$(_math(:𝔤))``,
where ``$(_math(:𝔤))`` denotes the [`LieAlgebra`](@ref) of ``$(_math(:G))``.
It is given by

```math
$(_tex(:exp))_g X = g$(_math(:∘))$(_tex(:exp))_{$(_math(:G))}(X)
```

where `X` can be scaled by `t`, the computation can be performed in-place of `h`,
and ``$(_tex(:exp))_{$(_math(:G))}`` denotes the  [Lie group exponential function](@ref exp(::AbstractLieGroup, ::Identity, :Any)).

If `g` is the [`Identity`](@ref) the [Lie group exponential function](@ref exp(::AbstractLieGroup, ::Identity, :Any)) ``$(_tex(:exp))_{$(_math(:G))}`` is computed directly.
Implementing the Lie group exponential function introduces a default implementation with the formula above.

!!! note
    The Lie group exponential map is usually different from the exponential map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) exponential map, use `exp(`[`base_manifold`](@ref)`(G), g, X)`.
"""

@doc "$_doc_exp"
ManifoldsBase.exp(G::AbstractLieGroup, ::Any, ::Any)

@doc "$_doc_exp"
function ManifoldsBase.exp!(G::AbstractLieGroup, h, g, X)
    if Base.mightalias(g, h)
        compose!(G, h, g, exp(G, X))
        return h
    else
        exp!(G, h, X)
        compose!(G, h, g, h)
        return h
    end
end

_doc_exponential = """
    exp(G::AbstractLieGroup, X::T)
    exp!(G::AbstractLieGroup, g, X)

Compute the (Lie group) exponential function

```math
$(_tex(:exp))_{$(_math(:G))}: $(_math(:𝔤)) → $(_math(:G)),$(_tex(:qquad)) $(_tex(:exp))_{$(_math(:G))}(X) = γ_X(1),
```

where ``γ_X`` is the unique solution of the initial value problem

```math
γ(0) = $(_math(:e)), $(_tex(:quad)) γ'(s) = γ(s)$(_math(:act))X.
```

See also [HilgertNeeb:2012; Definition 9.2.2](@cite).
On matrix Lie groups this is the same as the [matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential).

The computation can be performed in-place of `g`.

!!! info "Naming convention"
    There are at least two different objects usually called “exponential” that need to be distinguished

    * the [(Riemannian) exponential map](@extref `Base.exp-Tuple{AbstractManifold, Any, Any}`) `exp(M, p, X)` from $(_link(:ManifoldsBase)).
      This can be accessed here using `exp(base_manifold(G), p, X)`
    * the exponential map for a (left/right/bi-invariant) Cartan-Schouten (pseudo-)metric `exp(G, g, X)`, which we use as a default within this package
    * the (matrix/Lie group) exponential function `exp(G, g)` which agrees with the previous one for `g` being the identity there.
"""

@doc "$(_doc_exponential)"
function ManifoldsBase.exp(G::AbstractLieGroup, X)
    g = allocate_result(G, exp, X)
    exp!(G, g, X)
    return g
end

@doc "$(_doc_exponential)"
ManifoldsBase.exp!(G::AbstractLieGroup, ::Any, ::Any)

function ManifoldsBase.get_coordinates(
        G::AbstractLieGroup, g, X, B::AbstractBasis{<:Any, TangentSpaceType}
    )
    return get_coordinates(LieAlgebra(G), X, B)
end
function ManifoldsBase.get_coordinates!(
        G::AbstractLieGroup, c, g, X, B::AbstractBasis{<:Any, TangentSpaceType}
    )
    get_coordinates!(LieAlgebra(G), c, X, B)
    return c
end

function ManifoldsBase.get_vector(
        G::AbstractLieGroup, g, c, B::AbstractBasis{<:Any, TangentSpaceType}
    )
    return get_vector(
        LieAlgebra(G),
        c,
        B;
        tangent_vector_type = ManifoldsBase.tangent_vector_type(G, typeof(g)),
    )
end
function ManifoldsBase.get_vector!(
        G::AbstractLieGroup, X, g, c, B::AbstractBasis{<:Any, TangentSpaceType}
    )
    return get_vector!(LieAlgebra(G), X, c, B)
end

_doc_identity_element = """
    identity_element(G::AbstractLieGroup)
    identity_element(G::AbstractLieGroup, T)
    identity_element!(G::AbstractLieGroup, e::T)

Return a point representation of the [`Identity`](@ref) on the [`AbstractLieGroup`](@ref) `G`.
By default this representation is the default array or number representation.
If there exist several representations, the type `T` can be used to distinguish between them,
and it should be provided for both the [`AbstractLieGroupPoint`](@ref) as well as the [`AbstractLieAlgebraTangentVector`](@ref)
if they differ, since maybe only one of these types might be available for the second signature.

It returns the corresponding default representation of ``e`` as a point on `G`.
This can be performed in-place of `e`.
"""
# `function identity_element end`
@doc "$(_doc_identity_element)"
function identity_element(G::AbstractLieGroup)
    e = ManifoldsBase.allocate_result(G, identity_element)
    return identity_element!(G, e)
end

function identity_element! end
@doc "$(_doc_identity_element)"
identity_element!(G::AbstractLieGroup, e)

ManifoldsBase.inner(G::AbstractLieGroup, g, X, Y) = inner(LieAlgebra(G), X, Y)

_doc_inv = """
    inv(G::AbstractLieGroup, g)
    inv!(G::AbstractLieGroup, h, g)

Compute the inverse group element ``g^{-1}`` with respect to the [`AbstractGroupOperation`](@ref) ``$(_math(:∘))``
on the [`AbstractLieGroup`](@ref) ``$(_math(:G))``,
that is, return the unique element ``h=g^{-1}`` such that ``h$(_math(:∘))g=$(_math(:e))``, where ``$(_math(:e))`` denotes the [`Identity`](@ref).

This can be done in-place of `h`, without side effects, that is you can do `inv!(G, g, g)`.

!!! info
    This function also handles the case where `g` is the [`Identity`](@ref)`(G)`.
    Since this would lead to ambiguities when implementing a new group operations,
    this function calls `_inv` and `_inv!`, respectively, which is meant for the actual computation of
    group operations on (non-[`Identity`](@ref)` but maybe its numerical representation) elements.
"""

@doc "$_doc_inv"
Base.inv(G::AbstractLieGroup, g) = _inv(G, g)
function _inv(G::AbstractLieGroup, g)
    h = allocate_result(G, inv, g)
    return inv!(G, h, g)  # while we could go to _inv! as well, someone might just have done inv!
end

function inv! end
@doc "$_doc_inv"
Manifolds.inv!(G::AbstractLieGroup, h, g) = _inv!(G, h, g)

function Base.inv(
        ::AbstractLieGroup{𝔽, O}, e::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return e
end

function Manifolds.inv!(
        G::AbstractLieGroup{𝔽, O}, g, ::Identity{O}
    ) where {𝔽, O <: AbstractGroupOperation}
    return identity_element!(G, g)
end

function _inv! end

_doc_inv_left_compose = """
    inv_left_compose(G::AbstractLieGroup, g, h)
    inv_left_compose!(G::AbstractLieGroup, k, g, h)

Compute the inverse of the left group operation ``λ_g(h) = g$(_math(:∘))h``,
on the [`AbstractLieGroup`](@ref) `G`, that is, compute ``λ_g^{-1}(h) = g^{-1}$(_math(:∘))h``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_left_compose)"
function inv_left_compose(G::AbstractLieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, inv_left_compose, g, h)
    return inv_left_compose!(G, k, g, h)
end

function inv_left_compose! end
@doc "$(_doc_inv_left_compose)"
function inv_left_compose!(G::AbstractLieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, k, h) # compose `k∘h` in-place of k
    return k
end

_doc_inv_right_compose = """
    inv_right_compose(G::AbstractLieGroup, h, g)
    inv_right_compose!(G::AbstractLieGroup, k, h, g)

Compute the inverse of the right group operation ``ρ_g(h) = h$(_math(:∘))g``,
on the [`AbstractLieGroup`](@ref) `G`, that is compute ``ρ_g^{-1}(h) = h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_right_compose)"
function inv_right_compose(G::AbstractLieGroup, h, g)
    k = ManifoldsBase.allocate_result(G, inv_right_compose, h, g)
    return inv_right_compose!(G, k, h, g)
end

function inv_right_compose! end
@doc "$(_doc_inv_right_compose)"
function inv_right_compose!(G::AbstractLieGroup, k, h, g)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, k) # compose `h∘k` in-place of k
    return k
end

"""
    inverse_retract(G::AbstractLieGroup, g, h, m::BaseManifoldInverseRetraction)

Compute the inverse retraction of `g` and `h` on the [`AbstractLieGroup`](@ref) `G`
by using an inverse retraction on the underlying manifold and pulling the result back to the Lie algebra.
"""
ManifoldsBase.inverse_retract(G::AbstractLieGroup, g, h, m::BaseManifoldInverseRetraction)

# Layer 3
function ManifoldsBase._inverse_retract!(
        G::AbstractLieGroup, X, g, h, m::BaseManifoldInverseRetraction
    )
    return inverse_retract_base_manifold!(G, X, g, h, m)
end
function ManifoldsBase._inverse_retract(
        G::AbstractLieGroup, g, h, m::BaseManifoldInverseRetraction
    )
    return inverse_retract_base_manifold(G, g, h, m)
end
function inverse_retract_base_manifold!(
        G::AbstractLieGroup, X, g, h, m::BaseManifoldInverseRetraction
    )
    inverse_retract!(base_manifold(G), X, g, h, m.inverse_retraction)
    # X is in TgM so we still ave to pull it back to TeM using
    # the left group opp diff.
    pull_back_tangent!(G, X, g, X)
    return X
end
function inverse_retract_base_manifold(
        G::AbstractLieGroup, g, h, m::BaseManifoldInverseRetraction
    )
    X = inverse_retract(base_manifold(G), g, h, m.inverse_retraction)
    # X is in TgM so we still ave to pull it back to TeM using
    # the left group opp diff.
    return pull_back_tangent(G, g, X)
end

function is_identity end
@doc """
    is_identity(G::AbstractLieGroup, q; kwargs...)

Check whether `q` is the identity on the [`AbstractLieGroup`](@ref) ``$(_math(:G))``.
This means it is either the [`Identity`](@ref)`{O}` with the respect to the corresponding
[`AbstractGroupOperation`](@ref) `O`, or (approximately) the correct point representation.

# See also

[`identity_element`](@ref), [`identity_element!`](@ref)
"""
is_identity(G::AbstractLieGroup, q)

# Declare as “fallback” for types

function is_identity(
        G::AbstractLieGroup{𝔽, O}, h::P; kwargs...
    ) where {𝔽, P, O <: AbstractGroupOperation}
    return ManifoldsBase.isapprox(G, identity_element(G, P), h; kwargs...)
end
function is_identity(
        ::AbstractLieGroup{𝔽, O}, ::Identity{O}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return true
end
# any other identity than the fitting one
function is_identity(
        G::AbstractLieGroup{𝔽, <:AbstractGroupOperation},
        h::Identity{<:AbstractGroupOperation};
        kwargs...,
    ) where {𝔽}
    return false
end

"""
    is_point(G::AbstractLieGroup, g; kwargs...)

Check whether `g` is a valid point on the Lie Group `G`.
This falls back to checking whether `g` is a valid point on the [`base_manifold`](@ref)`G`.
unless `g` is an [`Identity`](@ref). Then, it is checked whether it is the
identity element corresponding to `G`.
"""
ManifoldsBase.is_point(G::AbstractLieGroup, g; kwargs...)

# resolve identity already here, everything else passes down to checks.

function ManifoldsBase.is_point(
        G::AbstractLieGroup{𝔽, O}, e::Identity{O}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return true
end
function ManifoldsBase.is_point(
        G::AbstractLieGroup, e::Identity; error::Symbol = :none, kwargs...
    )
    s = """
    The provided point $e is not the Identity on $G.
    Expected an Identity corresponding to $(G.op).
    """
    (error === :error) && throw(DomainError(s))
    (error === :info) && @info s
    (error === :warn) && @warn s
    return false
end

function ManifoldsBase.is_vector(
        G::AbstractLieGroup{𝔽, O}, ::Identity{O}, X; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return ManifoldsBase.is_point(LieAlgebra(G), X; kwargs...)
end

"""
    isapprox(M::AbstractLieGroup, g, h; kwargs...)

Check if points `g` and `h` from [`AbstractLieGroup`](@ref) are approximately equal.
this function calls the corresponding $(_link(:isapprox)) on the $(_link(:AbstractManifold))
after handling the cases where one or more
of the points are the [`Identity`](@ref).
All keyword argments are passed to this function as well.
"""
ManifoldsBase.isapprox(G::AbstractLieGroup, g, h; kwargs...) =
    isapprox(base_manifold(G), g, h; kwargs...)
function ManifoldsBase.isapprox(
        G::AbstractLieGroup{𝔽, O}, g::Identity{O}, h; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return ManifoldsBase.isapprox(G, identity_element(G, typeof(h)), h; kwargs...)
end
function ManifoldsBase.isapprox(
        G::AbstractLieGroup{𝔽, O}, g, h::Identity{O}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return ManifoldsBase.isapprox(G, g, identity_element(G, typeof(g)); kwargs...)
end
function ManifoldsBase.isapprox(
        G::AbstractLieGroup{𝔽, O}, g::Identity{O}, h::Identity{O}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return true
end
function ManifoldsBase.isapprox(
        G::AbstractLieGroup{𝔽, O}, g::Identity{O}, h::Identity{O2}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation, O2 <: AbstractGroupOperation}
    return false
end
function ManifoldsBase.isapprox(
        G::AbstractLieGroup{𝔽, O}, g::Identity{O2}, h::Identity{O}; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation, O2 <: AbstractGroupOperation}
    return false
end

_doc_jacobian_conjugate = """
    jacobian_conjugate(G::AbstractLieGroup, g, h; basis::AbstractBasis=DefaultLieAlgebraOrthogonalBasis(); X=zero_vector(LieAlgebar(G)))
    jacobian_conjugate!(G::AbstractLieGroup, J, g, h; basis::AbstractBasis=DefaultLieAlgebraOrthogonalBasis(); X=zero_vector(LieAlgebar(G)))

Compute the Jacobian of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``,
with respect to an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`) of the [`LieAlgebra`](@ref).

A default is implemented using [`diff_conjugate`](@ref) ``$(_math(:d))(c_g(h))[X]``:
the ``j``th column of of the Jacobian matrix ``J`` are given by the coefficients of
the tangent vector ``$(_math(:d))(c_g(h))[X_j]`` with respect to the basis ``B``,
where ``X_j`` is the ``j``th basis vector of ``B``.

!!! note
    For the case that `h` is the [`Identity`](@ref) and the relation of ``$(_math(:d))(c_g(h))[X]``
    to the [`adjoint`](@ref) ``$(_math(:Ad))(g)``, the Jacobian then sometimes called “adjoint matrix”,
    e.g. in [SolaDerayAtchuthan:2021](@cite), when choosing as a basis the
    [`DefaultLieAlgebraOrthogonalBasis`](@ref)`()` that is used for [`hat`](@ref) and [`vee`](@ref).

## Keyword arguments

* `X=zero_vector(LieAlgebra(G))` pass an interims memory to store the Lie algebra tangent vector in.
"""
@doc "$(_doc_jacobian_conjugate)"
function jacobian_conjugate(
        G::AbstractLieGroup, g, h, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
    )
    J = ManifoldsBase.allocate_result(G, jacobian_conjugate, g, h, B)
    return jacobian_conjugate!(G, J, g, h, B)
end

function jacobian_conjugate! end
@doc "$(_doc_jacobian_conjugate)"
function jacobian_conjugate!(
        G::AbstractLieGroup,
        J,
        g,
        h,
        B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis();
        X = zero_vector(LieAlgebra(G)),
    )
    n = number_of_coordinates(base_manifold(G), B)
    𝔤 = LieAlgebra(G)
    c = zeros(eltype(J), n)
    for j in 1:n
        c .= 0
        c[j] = 1
        get_vector!(𝔤, X, c, B)         # store ``X_j`` in X
        diff_conjugate!(G, X, g, h, X) # compute the differential in-place
        get_coordinates!(𝔤, view(J, :, j), X, B)   # compute its coordinates in B and store the result in J
    end
    return J
end

_doc_jac_exp = """
    jacobian_exp(G::AbstractLieGroup, g, X, b)
    jacobian_exp!(G::AbstractLieGroup, J, g, X, b)

Compute the Jacobian of the [`exp`](@ref) ``$(_tex(:exp))_g(X)`` with respect to
an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`) of the [`LieAlgebra`](@ref).
"""

"$(_doc_jac_exp)"
function jacobian_exp(
        G::AbstractLieGroup, g, X, B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
    )
    J = ManifoldsBase.allocate_result(G, jacobian_exp, g, X, B)
    return jacobian_exp!(G, J, g, X, B)
end

function jacobian_exp! end
@doc "$(_doc_jac_exp)"
jacobian_exp!(
    G::AbstractLieGroup, J, g, X,
    B::AbstractBasis = DefaultLieAlgebraOrthogonalBasis()
)

_doc_log = """
    log(G::AbstractLieGroup, g, h)
    log!(G::AbstractLieGroup, X, g, h)

Compute the Lie group logarithmic map ``$(_tex(:log))_g: $(_math(:G)) → $(_math(:𝔤))``,
where ``$(_math(:𝔤))`` denotes the [`LieAlgebra`](@ref) of ``$(_math(:G))``.
It is given by

```math
$(_tex(:log))_g h = $(_tex(:log))_{$(_math(:G))}(g^{-1}$(_math(:∘))h)
```

where ``$(_tex(:log))_{$(_math(:G))}`` denotes the [Lie group logarithmic function](@ref log(::AbstractLieGroup, :Any))
The computation can be performed in-place of `X`.


!!! info "Naming convention"
    There are at least two different objects usually called “logarithm” that need to be distinguished
    * the [(Riemannian) logarithmic map](@extref `Base.log-Tuple{AbstractManifold, Any, Any}`) `log(M, p, X)` from $(_link(:ManifoldsBase))
    * the exponential map for a (left/right/bi-invariant) Cartan-Schouten (pseudo-)metric `exp(G, g, X)`, which we use as a default within this package
    * the (matrix/Lie group) exponential function `exp(G, g)` which agrees with the previous one for `g` being the identity there.
"""

@doc "$_doc_log"
function ManifoldsBase.log(G::AbstractLieGroup, g, h)
    X = allocate_result(G, log, g, h)
    log!(G, X, g, h)
    return X
end

@doc "$_doc_log"
function ManifoldsBase.log!(G::AbstractLieGroup, X, g, h)
    log!(G, X, compose(G, inv(G, g), h))
    return h
end

_doc_log = """
    log(G::AbstractLieGroup, g, h)
    log(G::AbstractLieGroup, g)
    log(G::AbstractLieGroup, g::Identity, T)
    log!(G::AbstractLieGroup, X::T, g)

Compute the (Lie group) logarithmic function ``$(_tex(:log))_{$(_math(:G))}: $(_math(:G)) → $(_math(:𝔤))``,
which is the inverse of the [Lie group exponential function](@ref exp(::AbstractLieGroup, :Any)).
For the allocating variant, you can specify the type `T`, when the point argument is the identity and hence does not provide the representation used.
The computation can be performed in-place of `X::T`, which then determines the type.

!!! info "Naming convention"
    There are at least two different objects usually called “logarithm” that need to be distinguished

    * the [(Riemannian) logarithm](@extref `Base.log-Tuple{AbstractManifold, Any, Any}`) map `log(M, p, q)` from $(_link(:ManifoldsBase)). This can be accessed here using `log(base_manifold(G), p, q)`.
    * the logarithmic map for a (left/right/bi-invariant) Cartan-Schouten (pseudo-)metric `log(G, g, h)`, which we use as a default within this package
    * the (matrix/Lie group) logarithm function `log(G, h)` which agrees with the previous one for `g` being the identity there.
"""

@doc "$(_doc_log)"
function ManifoldsBase.log(G::AbstractLieGroup, g)
    X = allocate_result(G, log, g)
    log!(G, X, g)
    return X
end
function ManifoldsBase.log(
        G::AbstractLieGroup{𝔽, Op}, e::Identity{Op}
    ) where {𝔽, Op <: AbstractGroupOperation}
    return zero_vector(LieAlgebra(G))
end
function ManifoldsBase.log(
        G::AbstractLieGroup{𝔽, Op}, ::Identity{Op}, T::Type
    ) where {𝔽, Op <: AbstractGroupOperation}
    return zero_vector(LieAlgebra(G), T)
end

@doc "$(_doc_log)"
ManifoldsBase.log!(G::AbstractLieGroup, ::Any, ::Any)

function ManifoldsBase.log!(
        G::AbstractLieGroup{𝔽, Op}, X, e::Identity{Op}
    ) where {𝔽, Op <: AbstractGroupOperation}
    return zero_vector!(LieAlgebra(G), X)
end

ManifoldsBase.manifold_dimension(G::AbstractLieGroup) = manifold_dimension(base_manifold(G))

ManifoldsBase.norm(G::AbstractLieGroup, g, X) = norm(LieAlgebra(G), X)

_doc_rand = """
    rand(::AbstractLieGroup; vector_at=nothing, σ::Real=1.0, kwargs...)
    rand(::AbstractLieGroup, PT::Type; vector_at=nothing, σ::Real=1.0, kwargs...)
    rand!(::LieAlgebra, T::Type; σ::Real=1.0, kwargs...)
    rand!(::AbstractLieGroup, gX::PT; vector_at=nothing, σ::Real=1.0, kwargs...)
    rand!(::LieAlgebra, X::T; σ::Real=1.0, kwargs...)

Compute a random point or tangent vector on a Lie group.

For points this just means to generate a random point on the
underlying manifold itself.

For tangent vectors, an element in the Lie Algebra is generated,
see also [`rand(::LieAlgebra; kwargs...)`](@ref)

For both cases, you can provide the type ``T`` for the tangent vector and/or point ``PT``,
if you want to generate a random point in a certain representation.
For the in-place variants the type is inferred from `pX` and `X`, respectively.
"""

function ManifoldsBase.project!(G::AbstractLieGroup, g, p)
    return ManifoldsBase.project!(base_manifold(G), g, p)
end

# TODO: Move to ManifoldsBase at some point
@doc """
    point_type(G::AbstractLieGroup, tangent_vector_type::Type)

Change `tangent_vector_type` that is a type of tangent vector type on Lie group `G`
to its matching type for representing points.

By default both these types are assumed to be identical.
"""
point_type(::AbstractLieGroup, tangent_vector_type::Type) = tangent_vector_type

_doc_pull_back_t = """
    pull_back_tangent(G::AnstractLieGroup, g, X; kwargs...)
    pull_back_tangent!(G::AbstractLiegroup, Y, g, X; kwargs...)

Given a tangent vector `X` on the tangent space at `g` interpreted as the one on the manifold,
this function pulls it back to the [`LieAlgebra`](@ref).

By default this function falls back to calling [`diff_left_compose`](@ref), but compared to
that function, this function also takes care about the change of representation.

For example if a default representation for tangent vectors on a manifold is also the Lie algebra,
then this function simplifies to the identity. This is for example the case for the
[`SpecialOrthogonalGroup`](@ref) and its [`Rotations`](@extref `Manifolds.Rotations`) manifold.

# Keyword argument
* `e = identity_element(G, typeof(g))` – if you have a memory available to store an identity point in,
  you can pass that memory here.
"""

@doc "$(_doc_pull_back_t)"
function pull_back_tangent(G::AbstractLieGroup, g, X; e = identity_element(G, typeof(g)))
    Y = zero_vector(LieAlgebra(G), typeof(X))
    return pull_back_tangent!(G, Y, g, X)
end

function pull_back_tangent! end
@doc "$(_doc_pull_back_t)"
pull_back_tangent!(G::AbstractLieGroup, Y, g, X; e = identity_element(G, typeof(g)))

_doc_push_fwd_t = """
    push_forward_tangent(G::AnstractLieGroup, g, X)
    push_forward_tangent!(G::AbstractLiegroup, Y, g, X)

Given a Lie algebra vector `X` on the [`LieAlgebra`](@ref), this function pushes
the vector forward to the tangent space at `g` interpreted as the one on the manifold.

By default this function falls back to calling [`diff_left_compose`](@ref), but compared to
that function, this function also takes care about the change of representation.

For example if a default representation for tangent vectors on a manifold is also the Lie algebra,
then this function simplifies to the identity. This is for example the case for the
[`SpecialOrthogonalGroup`](@ref) and its [`Rotations`](@extref `Manifolds.Rotations`) manifold.

# Keyword argument
* `e = identity_element(G, typeof(g))` – if you have a memory available to store an identity point in,
  you can pass that memory here.
"""

@doc "$(_doc_push_fwd_t)"
function push_forward_tangent(G::AbstractLieGroup, g, X; e = identity_element(G, typeof(g)))
    Y = zero_vector(LieAlgebra(G), typeof(X))
    return push_forward_tangent!(G, Y, g, X; e = e)
end

function push_forward_tangent! end
@doc "$(_doc_push_fwd_t)"
push_forward_tangent!(G::AbstractLieGroup, Y, g, X; e = identity_element(G, typeof(g)))

@doc "$(_doc_rand)"
Random.rand(::AbstractLieGroup; kwargs...)

# New in LIeGroups, maybe move to ManifoldsBase at some point
@doc "$(_doc_rand)"
Random.rand(G::AbstractLieGroup, T::Type; vector_at = nothing, kwargs...)

function Random.rand(G::AbstractLieGroup, T::Type, d::Integer; kwargs...)
    return [rand(G, T; kwargs...) for _ in 1:d]
end
function Random.rand(rng::AbstractRNG, G::AbstractLieGroup, T::Type, d::Integer; kwargs...)
    return [rand(rng, G, T; kwargs...) for _ in 1:d]
end
function Random.rand(G::AbstractLieGroup, d::Integer; kwargs...)
    return [rand(G; kwargs...) for _ in 1:d]
end
function Random.rand(G::AbstractLieGroup, T::Type; vector_at = nothing, kwargs...)
    if vector_at === nothing
        gX = allocate_on(G, T)
    else
        gX = allocate_on(G, TangentSpaceType(), T)
    end
    rand!(G, gX; vector_at = vector_at, kwargs...)
    return gX
end
function Random.rand(
        rng::AbstractRNG, M::AbstractLieGroup, T::Type; vector_at = nothing, kwargs...
    )
    if vector_at === nothing
        gX = allocate_on(M, T)
    else
        gX = allocate_on(M, TangentSpaceType(), T)
    end
    rand!(rng, M, gX; vector_at = vector_at, kwargs...)
    return gX
end

@doc "$(_doc_rand)"
function Random.rand!(G::AbstractLieGroup, pX; kwargs...)
    return rand!(Random.default_rng(), G, pX; kwargs...)
end

function Random.rand!(
        rng::AbstractRNG, G::AbstractLieGroup, pX::T; vector_at = nothing, kwargs...
    ) where {T}
    M = base_manifold(G)
    return if vector_at === nothing # for points -> pass to manifold
        rand!(rng, M, pX; kwargs...)
    else # for tangent vectors -> materialize identity, pass to tangent space there.
        rand!(rng, M, pX; vector_at = identity_element(G, T), kwargs...)
    end
end

"""
    retract(G::AbstractLieGroup, g, h, m::BaseManifoldRetraction)

Compute the retraction of `g` and `X` on the [`AbstractLieGroup`](@ref) `G`
by pushing `X` forward to the tangent space at `g` and using a retraction on the underlying manifold.
"""
ManifoldsBase.retract(::LieGroup, p, X, m::BaseManifoldRetraction)

# Layer 2
function ManifoldsBase._retract!(G::AbstractLieGroup, h, g, X, m::BaseManifoldRetraction)
    return retract_base_manifold!(G, h, g, X, m)
end
function ManifoldsBase._retract(G::AbstractLieGroup, g, X, m::BaseManifoldRetraction)
    return retract_base_manifold(G, g, X, m)
end
function retract_base_manifold!(G, h, g, X, m::BaseManifoldRetraction)
    # X is in TeM so we first push it to TpM using
    Y = push_forward_tangent(G, g, X)
    # now we can use the retraction on the base manifold
    retract!(base_manifold(G), h, g, Y, m.retraction)
    return h
end
function retract_base_manifold(G, g, X, m::BaseManifoldRetraction)
    # X is in TeM so we first push it to TpM using
    Y = push_forward_tangent(G, g, X)
    # now we can use the retraction on the base manifold
    return retract(base_manifold(G), g, Y, m.retraction)
end

function ManifoldsBase.representation_size(G::AbstractLieGroup)
    return representation_size(base_manifold(G))
end

function Base.show(io::IO, G::LieGroup)
    return print(io, "LieGroup($(base_manifold(G)), $(G.op))")
end

"""
    BaseManifoldVectorTransportMethod{VTM<:AbstractVectorTransportMethod} <:
        AbstractVectorTransportMethod

Compute a vector transport by using the transport of type `VTM` on the base manifold of
a [`LieGroup`](@ref).

# Constructor

    BaseManifoldVectorTransportMethod(vtm::AbstractVectorTransportMethod)

Generate the vector transport with transport `vtm` to use on the base manifold.
"""
struct BaseManifoldVectorTransportMethod{VTM <: AbstractVectorTransportMethod} <:
    AbstractVectorTransportMethod
    vector_transport_method::VTM
end

"""
    vector_transport_to(G::AbstractLieGroup, g, X, h, m::BaseManifoldVectorTransportMethod)

Compute the vector transport of a Lie algebra `X` from  `g` to `h` using a vector transport
on the underlying manifold. This is done by pushing `X` forward to the tangent space at `g`,
then performing the vector transport on the base manifold, and finally pulling the resulting
tangent vector back to the Lie algebra.

This method merely exists for experimental reasons, since the parallel transport on Lie groups,
where all tangent vectors are represented in the Lie algebra is the identity.
Hence any of the methods performed here are more costly than plain parallel transport.
"""
ManifoldsBase.vector_transport_to(
    G::AbstractLieGroup, g, X, h, m::BaseManifoldVectorTransportMethod
)

function ManifoldsBase._vector_transport_to!(
        G::AbstractLieGroup, Y, g, X, h, m::BaseManifoldVectorTransportMethod
    )
    return _vector_transport_to_basemanifold!(G, Y, g, X, h, m)
end
function ManifoldsBase._vector_transport_to(
        G::AbstractLieGroup, g, X, h, m::BaseManifoldVectorTransportMethod
    )
    return _vector_transport_to_basemanifold(G, g, X, h, m)
end

function _vector_transport_to_basemanifold!(
        G::AbstractLieGroup, Y, g, X, h, m::BaseManifoldVectorTransportMethod
    )
    # (a) we have to push forward X from TeG to TgG
    # we can do this in-place of Y
    push_forward_tangent!(G, Y, g, X)
    # then we do the vector transport purely in place of Y
    vector_transport_to!(base_manifold(G), Y, g, Y, h, m.vector_transport_method)
    # now Y is in ThM so we still ave to pull it back to TeM using
    # the left group opp diff.
    pull_back_tangent!(G, Y, g, X)
    return Y
end
function _vector_transport_to_basemanifold(
        G::AbstractLieGroup, g, X, h, m::BaseManifoldVectorTransportMethod
    )
    # (a) we have to push forward X from TeG to TgG
    Y = push_forward_tangent(G, g, X)
    # then we do the vector transport
    Y = vector_transport_to(base_manifold(G), g, Y, h, m.vector_transport_method)
    # now Y is in ThM so we still ave to pull it back to TeM using
    # the left group opp diff.
    return pull_back_tangent(G, g, X)
end

function ManifoldsBase.zero_vector(G::AbstractLieGroup{𝔽, O}, ::Identity{O}) where {𝔽, O <: AbstractGroupOperation}
    return zero_vector(LieAlgebra(G))
end

#
# Allocation hints - mainly pass-through, especially for power manifolds
function ManifoldsBase.allocate_on(G::AbstractLieGroup, T::Type{<:AbstractArray})
    return ManifoldsBase.allocate_on(base_manifold(G), T)
end

function ManifoldsBase.allocate_result(
        G::AbstractLieGroup,
        f::Union{typeof(compose), typeof(inv), typeof(conjugate), typeof(exp)},
        args...,
    )
    return ManifoldsBase.allocate_result(base_manifold(G), ManifoldsBase.exp, args...)
end
function ManifoldsBase.allocate_result(G::LieGroup, f::typeof(jacobian_conjugate), g, h, B)
    n = number_of_coordinates(G.manifold, B)
    return zeros(float(number_eltype(g)), n, n)
end
function ManifoldsBase.allocate_result(G::LieGroup, f::typeof(jacobian_exp), g, X, B)
    n = number_of_coordinates(G.manifold, B)
    return zeros(float(number_eltype(g)), n, n)
end
function ManifoldsBase.allocate_result(G::AbstractLieGroup, f::typeof(log), args...)
    return ManifoldsBase.allocate_result(base_manifold(G), f, args...)
end
function ManifoldsBase.allocate_result(
        G::AbstractLieGroup, f::Union{typeof(rand), typeof(identity_element)}
    )
    # both get a type allocated like rand
    return ManifoldsBase.allocate_result(base_manifold(G), rand)
end

#
#
# A fallback macro for types that merely wrap the actual data
"""
    default_lie_group_fallbacks(TG, TF, TP, TV, pfield::Symbol, Xfield::Symbol, groupOp)

Introduce default fallbacks for all basic functions on Lie groups, for Lie group of type
`TG` with group operation `Op`, points of type `TP`, tangent vectors of type `TV`, with
forwarding to fields `pfield` and `Xfield` for point and tangent vector functions,
respectively.
"""
macro default_lie_group_fallbacks(TG, Op, TP, TV, gfield::Symbol, Xfield::Symbol)
    block = quote
        function ManifoldsBase.allocate_result(::$TG, ::typeof(adjoint), ::$TP, X::$TV)
            return $TV(allocate(X.$Xfield))
        end
        function ManifoldsBase.allocate_result(::$TG, ::typeof(compose), g::$TP, ::$TP)
            return $TP(allocate(g.$gfield))
        end
        function ManifoldsBase.allocate_result(::$TG, ::typeof(exp), X::$TV)
            return $TP(allocate(X.$Xfield))
        end
        function ManifoldsBase.allocate_result(::$TG, ::typeof(inv), g::$TP)
            return $TP(allocate(g.$gfield))
        end
        function ManifoldsBase.allocate_result(::$TG, ::typeof(log), g::$TP)
            return $TV(allocate(g.$gfield))
        end

        function LieGroups.adjoint!(G::$TG, Y::$TV, g::$TP, X::$TV)
            LieGroups.adjoint!(G, Y.$Xfield, g.$gfield, X.$Xfield)
            return Y
        end

        # Could probably even be moved to MainfoldsBase?
        function ManifoldsBase.check_size(G::$TG, g::$TP; kwargs...)
            return ManifoldsBase.check_size(G, g.$gfield; kwargs...)
        end
        function ManifoldsBase.check_size(G::$TG, g::$TP, X::$TV; kwargs...)
            return ManifoldsBase.check_size(G, g.$gfield, X.$Xfield; kwargs...)
        end

        function LieGroups._compose!(G::$TG, k::$TP, g::$TP, h::$TP)
            LieGroups._compose!(G, k.$gfield, g.$gfield, h.$gfield)
            return k
        end
        function LieGroups.conjugate!(G::$TG, k::$TP, g::$TP, h::$TP)
            LieGroups.conjugate!(G, k.$gfield, g.$gfield, h.$gfield)
            return k
        end

        function LieGroups.diff_inv!(G::$TG, Y::$TV, g::$TP, X::$TV)
            LieGroups.diff_inv!(G, Y.$Xfield, g.$gfield, X.$Xfield)
            return Y
        end
        function LieGroups.diff_left_compose!(G::$TG, Y::$TV, g::$TP, h::$TP, X::$TV)
            LieGroups.diff_left_compose!(G, Y.$Xfield, g.$gfield, h.$gfield, X.$Xfield)
            return Y
        end
        function LieGroups.diff_right_compose!(G::$TG, Y::$TV, g::$TP, h::$TP, X::$TV)
            LieGroups.diff_right_compose!(G, Y.$Xfield, g.$gfield, h.$gfield, X.$Xfield)
            return Y
        end

        function LieGroups.exp!(G::$TG, g::$TP, X::$TV)
            LieGroups.exp!(G, g.$gfield, X.$Xfield)
            return g
        end

        function identity_element!(G::$TG, g::$TP)
            identity_element!(G, g.$gfield)
            return g
        end
        function LieGroups._inv!(G::$TG, h::$TP, g::$TP)
            LieGroups._inv!(G, h.$gfield, g.$gfield)
            return h
        end

        function LieGroups.is_identity(G::$TG, g::$TP; kwargs...)
            return LieGroups.is_identity(G, g.$gfield; kwargs...)
        end
        function ManifoldsBase.isapprox(
                G::$TG, e::Identity{<:$Op}, X::$TV, Y::$TV; kwargs...
            )
            return ManifoldsBase.isapprox(G, e, X.$Xfield, Y.$Xfield; kwargs...)
        end
        function LieGroups.log!(G::$TG, X::$TV, g::$TP)
            LieGroups.log!(G, X.$Xfield, g.$gfield)
            return X
        end
        function LieGroups.log!(G::$TG, X::$TV, e::Identity{<:$Op})
            LieGroups.log!(G, X.$Xfield, e)
            return X
        end
    end
    return esc(block)
end
