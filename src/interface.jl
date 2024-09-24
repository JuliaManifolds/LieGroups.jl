@doc """
    AbstractGroupAction

Represent a type of Group action for a [`LieGroup`](@ref) ``$(_math(:G))``
Abstract type for smooth binary operations ``$(_math(:op))`` on elements of a Lie group ``$(_math(:G))``:

```math
$(_math(:op)) : $(_math(:G)) × $(_math(:G)) → $(_math(:G))
```
"""
abstract type AbstractGroupOperation end

"""
    LieGroup{𝔽,M<:AbstractManifold{𝔽}, O<:AbstractGroupOperation}

Represent a Lie Group ``$(_math(:G))``.

A *Lie Group* ``$(_math(:G))`` is a group endowed with the structure of a manifold such that the
group operations ``$(_math(:op)): $(_math(:G))×$(_math(:G)) → $(_math(:G))``, see [`compose`](@ref)
and the inverse operation ``⋅^{-1}: $(_math(:G)) → $(_math(:G))``, see [`inv`](@ref) are smooth.

See for example [HilgertNeeb:2012; Definition 9.1.1](@cite).

# Fields

* `manifold`: an $(_link(:AbstractManifold)) ``$(_math(:M))``
* `op`: an [`AbstractGroupOperation`](@ref) ``$(_math(:op))`` on that manifold

# Constructor

    LieGroup(M::AbstractManifold, op::AbstractGroupOperation; vectors=LeftInvariantRepresentation())

Generate a Lie group based on a manifold `M` and a group operation `op`, where vectors by default are stored in the Lie Algebra.
"""
struct LieGroup{𝔽,M<:ManifoldsBase.AbstractManifold{𝔽},O<:AbstractGroupOperation} <:
       ManifoldsBase.AbstractManifold{𝔽}
    manifold::M
    op::O
end

@doc """
    Identity{O<:AbstractGroupOperation}

Represent the group identity element ``e ∈ $(_math(:G))`` on a [`LieGroup`](@ref) ``$(_math(:G))``
with [`AbstractGroupOperation`](@ref) of type `O`.

Similar to the philosophy that points are agnostic of their group at hand, the identity
does not store the group ``$(_math(:G))`` it belongs to. However it depends on the type of the [`AbstractGroupOperation`](@ref) used.

See also [`identity_element`](@ref) on how to obtain the corresponding [`AbstractManifoldPoint`](@extref `ManifoldsBase.AbstractManifoldPoint`) or array representation.

# Constructors

    Identity(::LieGroup{𝔽,M,O}) where {𝔽,M,O<:AbstractGroupOperation}
    Identity(o::AbstractGroupOperation)
    Identity(::Type{AbstractGroupOperation})

create the identity of the corresponding subtype `O<:`[`AbstractGroupOperation`](@ref)
"""
struct Identity{O<:AbstractGroupOperation} end

ManifoldsBase.@trait_function Identity(M::ManifoldsBase.AbstractDecoratorManifold)
Identity(::LieGroup{𝔽,M,O}) where {𝔽,M,O<:AbstractGroupOperation} = Identity{O}()
Identity(::O) where {O<:AbstractGroupOperation} = Identity(O)
Identity(::Type{O}) where {O<:AbstractGroupOperation} = Identity{O}()

"""
    LieAlgebra{𝔽, G} <: AbstractManifold{𝔽}

Represent the Lie Algebra ``$(_math(:𝔤))``, that is a ``𝔽``vector space with an associated
[`Lie_bracket`](@ref) ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfils

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity holds ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The Lie algebras considered here are those related to a [`LieGroup`](@ref) ``$(_math(:G))``,
namely the tangent space the tangent space ``T_{$(_math(:e))}$(_math(:G))`` at the [`Identity`](@ref),
this is internally just a `const` of the corresponding $(_link(:TangentSpace)).

# Constructor

    LieAlgebra(G::LieGroup)

Return the Lie Algebra belonging to the [`LieGroup`](@ref) `G`.
"""
const LieAlgebra{𝔽,G,I} = ManifoldsBase.Fiber{
    𝔽,ManifoldsBase.TangentSpaceType,G,I
} where {𝔽,G<:LieGroup{𝔽},I<:Identity}

function LieAlgebra(G::LieGroup{𝔽}) where {𝔽}
    return LieAlgebra{𝔽,G,typeof(Identity(G))}(
        G, Identity(G), ManifoldsBase.TangentSpaceType()
    )
end

#
#
# Traits for LIe groups
"""
    AbstractInvarianceTrait <: AbstractTrait

A common supertype for anz [`AbstractTrait`](@extref `ManifoldsBase.AbstractTrait`) related to metric invariance
"""
abstract type AbstractInvarianceTrait <: ManifoldsBase.AbstractTrait end

"""
    HasLeftInvariantMetric <: AbstractInvarianceTrait

Specify that the defaut metric on a [`LieGroup`](@ref) ``$(_math(:G))`` is left-invariant.
"""
struct HasLeftInvariantMetric <: AbstractInvarianceTrait end

"""
    HasRightInvariantMetric <: AbstractInvarianceTrait

Specify that the defaut metric on a [`LieGroup`](@ref) ``$(_math(:G))`` is right-invariant.
"""
struct HasRightInvariantMetric <: AbstractInvarianceTrait end

"""
    HasBiinvariantMetric <: AbstractInvarianceTrait

Specify that the defaut metric on a [`LieGroup`](@ref) ``$(_math(:G))`` is bi-invariant.
"""
struct HasBiinvariantMetric <: AbstractInvarianceTrait end
function parent_trait(::HasBiinvariantMetric)
    return ManifoldsBase.TraitList(HasLeftInvariantMetric(), HasRightInvariantMetric())
end

#
#
# --- Functions ---
# compose g_1 ∘ g_2 and its frieds
_doc_compose = """
    compose(G::LieGroup, g, h)
    compose!(G::LieGroup, k, g, h)

Perform the group oepration ``g $(_math(:op)) h`` for two ``g, h ∈ $(_math(:G))``
on the [`LieGroup`](@ref) `G`. This can also be done in-place of `h`.
"""
@doc "$(_doc_compose)"
function compose(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose, g, h)
    return compose!(G, k, g, h)
end

function compose! end
@doc "$(_doc_compose)"
compose!(::LieGroup, k, g, h)

_doc_compose_diff_left = """
    compose_diff_left(G::LieGroup, g, h, X)
    compose_diff_left!(G::LieGroup, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:op))h``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dλ_g(h)[X]``, ``X ∈ 𝔤``.
This can be done in-place of `Y`.
"""
@doc "$(_doc_compose_diff_left)"
function compose_diff_left(G::LieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, compose_diff_left, g, h, X)
    return compose_diff_left!(G, Y, g, h, X)
end

function compose_diff_left! end
@doc "$(_doc_compose_diff_left)"
compose_diff_left!(::LieGroup, Y, g, h, X)

_doc_compose_diff_right = """
    compose_diff_right(G::LieGroup, h, g, X)
    compose_diff_right!(G::LieGroup, Y, h, g, X)

Compute the differential of the right group multiplication ``ρ_g(h) = h$(_math(:op))g``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dρ_g(h)[X]``, ``X ∈ 𝔤``
This can be done in-place of `Y`.
"""
@doc "$(_doc_compose_diff_right)"
function compose_diff_right(G::LieGroup, h, g, X)
    Y = ManifoldsBase.allocate_result(G, compose_diff_right, h, g, X)
    return compose_diff_right!(G, Y, h, g, X)
end

function compose_diff_right! end
@doc "$(_doc_compose_diff_right)"
compose_diff_right!(::LieGroup, h, g1, g2)

# ---
_doc_compose_inv_left = """
    compose_inv_left(G::LieGroup, g, h)
    compose_inv_left!(G::LieGroup, k, g, h)

Compute the inverse of the left group multiplication ``λ_g(h) = g$(_math(:op))h``,
on the [`LieGroup`](@ref) `G`, that is, compute ``λ_g^{-1}(h) = g^{-1}$(_math(:op))h``.
This can be done in-place of `k`.
"""
@doc "$(_doc_compose_inv_left)"
function compose_inv_left(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose_inv_left, g, h)
    return compose_inv_left!(G, k, g, h)
end

function compose_inv_left! end
@doc "$(_doc_compose)"
function compose_inv_left!(::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, k, h) # kh inplace of k
    return k
end

_doc_compose_inv_right = """
    compose_inv_right(G::LieGroup, h, g)
    compose_inv_right!(G::LieGroup, k, h, g)

Compute the inverse of the right group multiplication ``ρ_g(h) = h$(_math(:op))g``,
on the [`LieGroup`](@ref) `G`, that is Compute ``ρ_g^{-1}(h) = h$(_math(:op))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_compose_inv_right)"
function compose_inv_right(G::LieGroup, h, g)
    k = ManifoldsBase.allocate_result(G, compose_inv_right, h, g)
    return compose_inv_right!(G, k, h, g)
end

function compose_inv_right! end
@doc "$(_doc_compose_inv_right)"
function compose_inv_right!(::LieGroup, k, h, g)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, h) # hk inplace of k
    return k
end

_doc_conjugate = """
    conjugate(G::LieGroup, g, h)
    conjugate!(G::LieGroup, k, g, h)

Compute the conjugation map ``c_g: $(_math(:G)) → $(_math(:G))`` given by ``c_g(h) = g$(_math(:op))h$(_math(:op))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_conjugate)"
function conjugate(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose_inv_right, h, g)
    return conjugate!(G, k, g, h)
end

function conjugate! end
@doc "$(_doc_conjugate)"
function conjugate!(::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, h) # hk inplace of k
    compose!(G, k, g, k) # gk inplace of k
    return k
end

_doc_exp = """
    exp(G::LieGroup, g, X, t::Number=1)
    exp!(G::LieGroup, h, g, X, t::Number=1)

Compute the Lie group exponential map given by

```math
$(_tex(:exp))_g X = g$(_math(:op))$(_tex(:exp))_{$(_math(:G))}(X)
```

where `X` can be scaled by `t`, the computation can be performed in-place of `h`,
and ``$(_tex(:exp))_{$(_math(:G))}`` denotes the  [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)).

!!! note
    If `g` is the [`Identity`](@ref) the [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)) is computed directly.
    Implementing the Lie group exponential function introduces a default implementation for this function.

!!! note
    The Lie group exponential map is usually different from the exponential map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) exponential map, use `exp(`[`manifold`](@ref)`G, g, X)`.
"""

@doc "$_doc_exp"
function ManifoldsBase.exp(G::LieGroup, g, X, t::Number=1)
    h = allocate_result(G, exp, g)
    exp!(G, h, g, X, t)
    return h
end

@doc "$_doc_exp"
function ManifoldsBase.exp!(G::LieGroup, h, g, X, t::Number=1)
    exp!(G, h, Identity(G), X, t)
    compose!(G, h, g, h)
    return h
end

_doc_exp_id = """
    exp(G::LieGroup, e::Identity, X)
    exp!(G::LieGroup, h, e::Identity, X)

Compute the (Lie group) exponential function

```math
$(_tex(:exp))_{$(_math(:G))}: $(_math(:𝔤)) → $(_math(:G)),$(_tex(:qquad)) $(_tex(:exp))_{$(_math(:G))}(X) = γ_X(1),
```

where ``γ_X`` is the unique solution of the initial value problem

```math
γ(0) = $(_math(:e)), $(_tex(:quad)) γ'(t) = γ(t)$(_math(:act))X.
```

The computation can be performed in-place of `h`.

See also [HilgertNeeb:2012; Definition 9.2.2](@cite).
"""

@doc "$(_doc_exp_id)"
function ManifoldsBase.exp(G::LieGroup, e::Identity, X, t::Number=1)
    h = identity_element(G) # allocate
    exp!(G, h, e, X, t)
    return h
end

@doc "$(_doc_exp_id)"
ManifoldsBase.exp!(::LieGroup, h, ::Identity, X, ::Number=1)

function is_identity end
@doc """
    is_identity(G::LieGroup, q; kwargs)

Check whether `q` is the identity on the [`LieGroup`](@ref) ``$(_math(:G))``.
This means it is either the [`Identity`](@ref)`{O}` with the respect to the corresponding
[`AbstractGroupOperation`](@ref) `O`, or (approximately) the correct point representation.

# See also

[`identity_element`](@ref), [`identity_element!`](@ref)
"""
is_identity(G::AbstractDecoratorManifold, q)

_doc_identity_element = """
    identity_element(G::LieGroup)
    identity_element!(G::LieGroup, g)

Return a point representation of the [`Identity`](@ref) on the [`LieGroup`](@ref) `G`.
By default this representation is the default array or number representation.
It should return the corresponding default representation of ``e`` as a point on `G` if
points are not represented by arrays.
This can be performed in-place of `g`.
"""
function identity_element end
@doc "$(_doc_identity_element)"
function identity_element(G::AbstractDecoratorManifold)
    g = ManifoldsBase.allocate_result(G, identity_element)
    return identity_element!(G, g)
end
function identity_element! end
@doc "$(_doc_identity_element)"
identity_element!(G::AbstractDecoratorManifold, g)

_doc_inv = """
    inv(G, g)
    inv!(G, h, g)

Compute the inverse group element ``g^{-1}`` with respect to the [`AbstractGroupOperation`](@ref) ``$(_math(:op))``
on the [`LieGroup`](@ref) ``$(_math(:G))``,
that is, return the unique element ``h=g^{-1}`` such that ``h$(_math(:op))g=$(_math(:e))``, where ``$(_math(:e))`` denotes the [`Identity`](@ref).

This can be done in-place of `h`, without side effects, that is you can do `inv!(G, g, g)`.
"""

@doc "$_doc_inv"
function Base.inv(::LieGroup, g)
    h = allocate_result(G, inv, g)
    return inv!(G, h, g)
end

function inv! end
@doc "$_doc_inv"
inv!(G::LieGroup, h, g)

_doc_inv_diff = """
    inv_diff(G, g, X)
    inv_diff!(G, Y, g, X)

Compute the differential of the function ``ι_{$(_math(:G))}(p) = p^-1``, where
``Dι_{$(_math(:G))}(p): $(_math(:𝔤)) → $(_math(:𝔤)).
This can be done in-place of `Y`.
"""

@doc "$_doc_inv_diff"
function inv_diff(::LieGroup, g, X)
    Y = allocate_result(G, inv_diff, g, X)
    return inv_diff!(G, Y, g, X)
end

function inv_diff! end
@doc "$_doc_inv_diff"
inv_diff!(G::LieGroup, Y, g, X)

_doc_Lie_bracket = """
    Lie_bracket!(𝔤::LieAlgebra, X, Y)
    Lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfils

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity holds ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The computation can be done in-place of `Z`.
"""
function Lie_bracket end
@doc "$(_doc_Lie_bracket)"
function Lie_bracket(𝔤::LieAlgebra, X, Y)
    Z = ManifoldsBase.allocate_result(𝔤, Lie_bracket, X, Y)
    return Lie_bracket!(𝔤, Z, X, Y)
end

function Lie_bracket! end
@doc "$(_doc_Lie_bracket)"
Lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

_doc_log = """
    log(G::LieGroup, g, h)
    log!(G::LieGroup, X, g, h)

Compute the Lie group logarithmic map

```math
$(_tex(:log))_g h = $(_math(:op))$(_tex(:log))_{$(_math(:G))}(g^{-1}$(_math(:op))h)
```

where ``$(_tex(:log))_{$(_math(:G))}`` denotes the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any))
The computation can be performed in-place of `X`.

!!! note
    If `g` is the [`Identity`](@ref) the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any)) is computed directly.
    Implementing the Lie group logarithmic function introduces a default implementation for this function.

!!! note
    The Lie group logarithmic map is usually different from the logarithmic map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) logarithmic map, use `log(`[`manifold`](@ref)`G, g, h)`.
"""

@doc "$_doc_log"
function ManifoldsBase.log(G::LieGroup, g, h)
    X = allocate_result(G, log, g, h)
    log!(G, X, g, h)
    return X
end

@doc "$_doc_log"
function ManifoldsBase.log!(G::LieGroup, X, g, h)
    log!(G, X, Identity(G), compose(G, inv(G, g), h))
    return h
end

_doc_log_id = """
    log(G::LieGroup, e::Identity, g)
    log!(G::LieGroup, X, e::Identity, g)

Compute the (Lie group) logarithmic function ``$(_tex(:log))_{$(_math(:G))}: $(_math(:G)) → $(_math(:𝔤))``,
which is the inverse of the [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any))
The compuation can be performed in-place of `X`.
"""

@doc "$(_doc_log_id)"
function ManifoldsBase.log(G::LieGroup, e::Identity, g)
    X = allocate_result(G, log, g)
    log!(G, X, e, g)
    return X
end

@doc "$(_doc_log_id)"
ManifoldsBase.log!(::LieGroup, X, ::Identity, g)

@doc """
    base_manifold(G::LieGroup)

Return the manifold stored within the [`LieGroup`](@ref) `G`.
"""
Manifolds.base_manifold(G::LieGroup) = G.manifold
