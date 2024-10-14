#
#
# Generic types for the interface of a Lie group
@doc """
    AbstractGroupOperation

Represent a type of group operation for a [`LieGroup`](@ref) ``$(_math(:G))``, that is a
smooth binary operation ``$(_math(:∘)) : $(_math(:G)) × $(_math(:G)) → $(_math(:G))``
on elements of a Lie group ``$(_math(:G))``.
"""
abstract type AbstractGroupOperation end

"""
    LieGroup{𝔽, O<:AbstractGroupOperation, M<:AbstractManifold{𝔽}}

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
struct LieGroup{𝔽,O<:AbstractGroupOperation,M<:ManifoldsBase.AbstractManifold{𝔽}} <:
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

    Identity(::LieGroup{𝔽,O}) where {𝔽,M,O<:AbstractGroupOperation}
    Identity(o::AbstractGroupOperation)
    Identity(::Type{AbstractGroupOperation})

create the identity of the corresponding subtype `O<:`[`AbstractGroupOperation`](@ref)
"""
struct Identity{O<:AbstractGroupOperation} end

Identity(::LieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation} = Identity{O}()
Identity(::O) where {O<:AbstractGroupOperation} = Identity(O)
Identity(::Type{O}) where {O<:AbstractGroupOperation} = Identity{O}()

"""
    LieAlgebra{𝔽, G} <: AbstractManifold{𝔽}

Represent the Lie Algebra ``$(_math(:𝔤))``, that is a ``𝔽``vector space with an associated
[`lie_bracket`](@ref) ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

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
    return LieAlgebra{𝔽,typeof(G),typeof(Identity(G))}(
        G, Identity(G), ManifoldsBase.TangentSpaceType()
    )
end

#
#
# --- Functions ---

_doc_adjoint = """
    adjoint(G::LieGroup, g X)
    adjoint!(G::LieGroup, Y, g, X)

Compute the adjoint ``$(_math(:Ad))(g): $(_math(:𝔤)) → $(_math(:𝔤))``, which is defined as
the differential [`diff_conjugate`](@ref) of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``
evaluated at the [`Identity`](@ref) ``h=$(_math(:e))``.
The operation can be performed in-place of `Y`.

```math
  $(_math(:Ad))(g)[X] = D c_g($(_math(:e))) [X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```

see [HilgertNeeb:2012; Section 9.2.3](@cite).

On matrix Lie groups the adjoint reads ``$(_math(:Ad))(g)[X] = g$(_math(:∘))X$(_math(:∘))g^{-1}``.
"""

@doc "$(_doc_adjoint)"
function Base.adjoint(G::LieGroup, g, X)
    Y = ManifoldsBase.allocate_result(G, adjoint, g, X)
    return adjoint!(G, Y, g, X)
end

function adjoint! end
@doc "$(_doc_adjoint)"
function adjoint!(::LieGroup, Y, g, X)
    diff_conjugate!(G, Y, g, Identity(G), X)
    return Y
end

#
# Allocation hints
function ManifoldsBase.allocate_result(G::LieGroup, f::typeof(identity_element))
    apf = ManifoldsBase.allocation_promotion_function(G, f, ())
    rs = ManifoldsBase.representation_size(G)
    return ManifoldsBase.allocate_result_array(G, f, apf(Float64), rs)
end

@doc """
    base_manifold(G::LieGroup)

Return the manifold stored within the [`LieGroup`](@ref) `G`.
"""
Manifolds.base_manifold(G::LieGroup) = G.manifold

function ManifoldsBase.check_point(
    G::LieGroup{𝔽,O}, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return nothing
end
function ManifoldsBase.check_point(
    G::LieGroup{𝔽,O}, e::Identity{O2}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,O2}
    return DomainError(
        e,
        """
        The provided point $e is not the Identity on $G.
        Expected an Identity corresponding to $(G.op).
        """,
    )
end

# compose g ∘ h
_doc_compose = """
    compose(G::LieGroup, g, h)
    compose!(G::LieGroup, k, g, h)

Perform the group oepration ``g $(_math(:∘)) h`` for two ``g, h ∈ $(_math(:G))``
on the [`LieGroup`](@ref) `G`. This can also be done in-place of `h`.

!!! info
    This function also handles the case where `g` or/and `h` are the [`Identity`](@ref)`(G)`.
    Since this would lead to ambiguities when implementing a new group operations,
    this function calls `_compose` and `compose!`, respectively, which is meant for the actual computation of
    group operations on (non-[`Identity`](@ref)` but maybe its numerical representation) elements.
"""
@doc "$(_doc_compose)"
compose(G::LieGroup, g, h) = _compose(G, g, h)
compose(::LieGroup{𝔽,O}, g::Identity{O}, h) where {𝔽,O<:AbstractGroupOperation} = h
compose(::LieGroup{𝔽,O}, g, h::Identity{O}) where {𝔽,O<:AbstractGroupOperation} = g
function compose(
    ::LieGroup{𝔽,O}, g::Identity{O}, h::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return g
end

function _compose(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, compose, g, h)
    return _compose!(G, k, g, h)
end

function compose! end

@doc "$(_doc_compose)"
compose!(G::LieGroup, k, g, h) = _compose!(G, k, g, h)
function compose!(G::LieGroup{𝔽,O}, k, ::Identity{O}, h) where {𝔽,O<:AbstractGroupOperation}
    return copyto!(G, k, h)
end
function compose!(G::LieGroup{𝔽,O}, k, g, ::Identity{O}) where {𝔽,O<:AbstractGroupOperation}
    return copyto!(G, k, g)
end
function compose!(
    G::LieGroup{𝔽,O}, k, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return identity_element!(G, k)
end
function compose!(
    ::LieGroup{𝔽,O}, k::Identity{O}, ::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return k
end

function _compose! end

_doc_conjugate = """
    conjugate(G::LieGroup, g, h)
    conjugate!(G::LieGroup, k, g, h)

Compute the conjugation map ``c_g: $(_math(:G)) → $(_math(:G))`` given by ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_conjugate)"
function conjugate(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, inv_right_compose, h, g)
    return conjugate!(G, k, g, h)
end

function conjugate! end
@doc "$(_doc_conjugate)"
function conjugate!(::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, h) # `h∘k` in-place of k
    compose!(G, k, g, k) # `g∘k` in-place of k
    return k
end

ManifoldsBase.copyto!(G::LieGroup, h, g) = ManifoldsBase.copyto!(G.manifold, h, g)
function ManifoldsBase.copyto!(G::LieGroup{𝔽,O}, h, g::Identity{O}) where {𝔽,O}
    return ManifoldsBase.copyto!(G.manifold, h, identity_element(g))
end
function ManifoldsBase.copyto!(G::LieGroup{𝔽,O}, h::Identity{O}, g::Identity{O}) where {𝔽,O}
    return ManifoldsBase.copyto!(G.manifold, h, identity_element(g))
end

_doc_diff_conjugate = """
    diff_conjugate(G::LieGroup, g, h, X)
    diff_conjugate!(G::LieGroup, Y, g, h, X)

Compute the differential of the [`conjugate`](@ref) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}``,
which can be performed in-place of `Y`.

```math
  D(c_g(h))[X], $(_tex(:qquad)) X ∈ $(_math(:𝔤)).
```
"""
@doc "$(_doc_diff_conjugate)"
function diff_conjugate(G::LieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_conjugate, g, h, X)
    return diff_conjugate!(G, Y, g, h, X)
end

_doc_diff_left_compose = """
    diff_left_compose(G::LieGroup, g, h, X)
    diff_left_compose!(G::LieGroup, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dλ_g(h)[X]``, ``X ∈ 𝔤``.
This can be done in-place of `Y`.
"""
@doc "$(_doc_diff_left_compose)"
function diff_left_compose(G::LieGroup, g, h, X)
    Y = ManifoldsBase.allocate_result(G, diff_left_compose, g, h, X)
    return diff_left_compose!(G, Y, g, h, X)
end

function diff_left_compose! end
@doc "$(_doc_diff_left_compose)"
diff_left_compose!(::LieGroup, Y, g, h, X)

_doc_diff_right_compose = """
    diff_right_compose(G::LieGroup, h, g, X)
    diff_right_compose!(G::LieGroup, Y, h, g, X)

Compute the differential of the right group multiplication ``ρ_g(h) = h$(_math(:∘))g``,
on the [`LieGroup`](@ref) `G`, that is Compute ``Dρ_g(h)[X]``, ``X ∈ 𝔤``
This can be done in-place of `Y`.
"""
@doc "$(_doc_diff_right_compose)"
function diff_right_compose(G::LieGroup, h, g, X)
    Y = ManifoldsBase.allocate_result(G, diff_right_compose, h, g, X)
    return diff_right_compose!(G, Y, h, g, X)
end

function diff_right_compose! end
@doc "$(_doc_diff_right_compose)"
diff_right_compose!(::LieGroup, h, g1, g2)

function diff_conjugate! end
@doc "$(_doc_diff_conjugate)"
diff_conjugate!(::LieGroup, Y, g, h, X)

_doc_diff_inv = """
    diff_inv(G::LieGroup, g, X)
    diff_inv!(G::LieGroup, Y, g, X)

Compute the differential of the function ``ι_{$(_math(:G))}(p) = p^-1``, where
``Dι_{$(_math(:G))}(p): $(_math(:𝔤)) → $(_math(:𝔤))``.
This can be done in-place of `Y`.
"""

@doc "$_doc_diff_inv"
function diff_inv(::LieGroup, g, X)
    Y = allocate_result(G, diff_inv, g, X)
    return diff_inv!(G, Y, g, X)
end

function diff_inv! end
@doc "$_doc_diff_inv"
diff_inv!(G::LieGroup, Y, g, X)

_doc_exp = """
    exp(G::LieGroup, g, X, t::Number=1)
    exp!(G::LieGroup, h, g, X, t::Number=1)

Compute the Lie group exponential map given by

```math
$(_tex(:exp))_g X = g$(_math(:∘))$(_tex(:exp))_{$(_math(:G))}(X)
```

where `X` can be scaled by `t`, the computation can be performed in-place of `h`,
and ``$(_tex(:exp))_{$(_math(:G))}`` denotes the  [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)).

!!! note
    If `g` is the [`Identity`](@ref) the [Lie group exponential function](@ref exp(::LieGroup, ::Identity, :Any)) is computed directly.
    Implementing the Lie group exponential function introduces a default implementation for this function.

!!! note
    The Lie group exponential map is usually different from the exponential map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) exponential map, use `exp(`[`base_manifold`](@ref)`G, g, X)`.
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
    h = identity_element(G)
    exp!(G, h, e, X, t)
    return h
end

@doc "$(_doc_exp_id)"
ManifoldsBase.exp!(::LieGroup, h, ::Identity, X, ::Number=1)

_doc_inv_left_compose = """
    inv_left_compose(G::LieGroup, g, h)
    inv_left_compose!(G::LieGroup, k, g, h)

Compute the inverse of the left group multiplication ``λ_g(h) = g$(_math(:∘))h``,
on the [`LieGroup`](@ref) `G`, that is, compute ``λ_g^{-1}(h) = g^{-1}$(_math(:∘))h``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_left_compose)"
function inv_left_compose(G::LieGroup, g, h)
    k = ManifoldsBase.allocate_result(G, inv_left_compose, g, h)
    return inv_left_compose!(G, k, g, h)
end

function inv_left_compose! end
@doc "$(_doc_compose)"
function inv_left_compose!(::LieGroup, k, g, h)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, k, h) # compose `k∘h` in-place of k
    return k
end

_doc_inv_right_compose = """
    inv_right_compose(G::LieGroup, h, g)
    inv_right_compose!(G::LieGroup, k, h, g)

Compute the inverse of the right group multiplication ``ρ_g(h) = h$(_math(:∘))g``,
on the [`LieGroup`](@ref) `G`, that is Compute ``ρ_g^{-1}(h) = h$(_math(:∘))g^{-1}``.
This can be done in-place of `k`.
"""
@doc "$(_doc_inv_right_compose)"
function inv_right_compose(G::LieGroup, h, g)
    k = ManifoldsBase.allocate_result(G, inv_right_compose, h, g)
    return inv_right_compose!(G, k, h, g)
end

function inv_right_compose! end
@doc "$(_doc_inv_right_compose)"
function inv_right_compose!(::LieGroup, k, h, g)
    inv!(G, k, g) # g^{-1} in-place of k
    compose!(G, k, h, k) # compose `h∘k` in-place of k
    return k
end

_doc_identity_element = """
    identity_element(G::LieGroup)
    identity_element!(G::LieGroup, g)

Return a point representation of the [`Identity`](@ref) on the [`LieGroup`](@ref) `G`.
By default this representation is the default array or number representation.
It should return the corresponding default representation of ``e`` as a point on `G` if
points are not represented by arrays.
This can be performed in-place of `g`.
"""
# `function identity_element end`
@doc "$(_doc_identity_element)"
function identity_element(G::LieGroup)
    g = ManifoldsBase.allocate_result(G, identity_element)
    return identity_element!(G, g)
end

function identity_element! end
@doc "$(_doc_identity_element)"
identity_element!(G::LieGroup, g)

_doc_inv = """
    inv(G::LieGroup, g)
    inv!(G::LieGroup, h, g)

Compute the inverse group element ``g^{-1}`` with respect to the [`AbstractGroupOperation`](@ref) ``$(_math(:∘))``
on the [`LieGroup`](@ref) ``$(_math(:G))``,
that is, return the unique element ``h=g^{-1}`` such that ``h$(_math(:∘))g=$(_math(:e))``, where ``$(_math(:e))`` denotes the [`Identity`](@ref).

This can be done in-place of `h`, without side effects, that is you can do `inv!(G, g, g)`.
"""

@doc "$_doc_inv"
function Base.inv(G::LieGroup, g)
    h = allocate_result(G, inv, g)
    return inv!(G, h, g)
end

function inv! end
@doc "$_doc_inv"
inv!(G::LieGroup, h, g)

function is_identity end
@doc """
    is_identity(G::LieGroup, q; kwargs)

Check whether `q` is the identity on the [`LieGroup`](@ref) ``$(_math(:G))``.
This means it is either the [`Identity`](@ref)`{O}` with the respect to the corresponding
[`AbstractGroupOperation`](@ref) `O`, or (approximately) the correct point representation.

# See also

[`identity_element`](@ref), [`identity_element!`](@ref)
"""
is_identity(G::LieGroup, q)

"""
    is_point(𝔤::LieAlgebra, X; kwargs...)

Check whether `X` is a valid point on the Lie Algebra `𝔤`.
This falls back to checking whether `X` is a valid point on the tangent space
at the [`identity_element``](@ref)`]`(G)` on `G.manifold` on the [`LieGroup`](@ref)
of `G`
"""
function ManifoldsBase.is_point(𝔤::LieAlgebra, g; kwargs...)
    # the manifold stored in the Fiber / Lie Algebra is the Lie Group G
    G = 𝔤.manifold
    e = identity_element(G)
    return ManifoldsBase.is_vector(G.manifold, e, X; kwargs...)
end

"""
    is_point(G::LieGroup, g; kwargs...)

Check whether `g` is a valid point on the Lie Group `G`.
This falls back to checking whether `g` is a valid point on `G.manifold`,
unless `g` is an [`Identity`](@ref). Then, it is checked whether it is the
idenity element corresponding to `G`.
"""
ManifoldsBase.is_point(G::LieGroup, g; kwargs...) =
    ManifoldsBase.is_point(G.manifold, g; kwargs...)

# Any Identity: pass to check point
function ManifoldsBase.is_point(G::LieGroup, e::Identity; kwargs...)
    return ManifoldsBase.is_point(G, e; kwargs...)
end

_doc_is_vector = """
    is_vector(G::LieGroup, X; kwargs...)
    is_vector(G::LieGroup{𝔽,O}, e::Indentity{O}, X; kwargs...)

Check whether `X` is a tangent vector, that is an element of the |`LieAlgebra`](@ref)
of `G`.
The first variant calls [`is_point`](@extref ManifoldsBase.is_point) on the [`LieAlgebra`](@ref) `𝔤` of `G`.
The second variant calls [`is_vector`](@extref ManifoldsBase.is_vector) on the $(_link(:AbstractManifold)) at the [`identity_element`](@ref).

All kwyword arguments are passed on to the corresponding call
"""

@doc "$(_doc_is_vector)"
ManifoldsBase.is_vector(G::LieGroup, X; kwargs...) =
    ManifoldsBase.is_point(LieAlgebra(G), X; kwargs...)

@doc "$(_doc_is_vector)"
function ManifoldsBase.is_vector(G::LieGroup{𝔽,O}, e::Identity{O}, X; kwargs...) where {𝔽,O}
    return ManifoldsBase.is_vector(G.manifold, identity_element(G), X; kwargs...)
end

"""
    isapprox(M::LieGroup, g, h; kwargs...)

Check if points `g` and `h` from [`LieGroup`](@ref) are approximately equal.
this function calls the corresponding $(_link(:isapprox)) on the $(_link(:AbstractManifold))
after handling the cases where one or more
of the points are the [`Identity`](@ref).
All keyword argments are passed to this function as well.
"""
ManifoldsBase.isapprox(G::LieGroup, g, h; kwargs...) = isapprox(G.manifold, g, h; kwargs...)
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g::Identity{<:O}, h; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.isapprox(G.manifold, identity_element(G), h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.isapprox(G.manifold, g, identity_element(G); kwargs...)
end
function ManifoldsBase.isapprox(
    G::LieGroup{𝔽,O}, g::Identity{O}, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return true
end

_doc_lie_bracket = """
    lie_bracket!(𝔤::LieAlgebra, X, Y)
    lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity holds ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

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

_doc_log = """
    log(G::LieGroup, g, h)
    log!(G::LieGroup, X, g, h)

Compute the Lie group logarithmic map

```math
$(_tex(:log))_g h = $(_math(:∘))$(_tex(:log))_{$(_math(:G))}(g^{-1}$(_math(:∘))h)
```

where ``$(_tex(:log))_{$(_math(:G))}`` denotes the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any))
The computation can be performed in-place of `X`.

!!! note
    If `g` is the [`Identity`](@ref) the [Lie group logarithmic function](@ref log(::LieGroup, ::Identity, :Any)) is computed directly.
    Implementing the Lie group logarithmic function introduces a default implementation for this function.

!!! note
    The Lie group logarithmic map is usually different from the logarithmic map with respect
    to the metric of the underlying Riemannian manifold ``$(_math(:M))``.
    To access the (Riemannian) logarithmic map, use `log(`[`base_manifold`](@ref)`G, g, h)`.
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

function ManifoldsBase.representation_size(G::LieGroup)
    return ManifoldsBase.representation_size(G.manifold)
end

Base.show(io::IO, G::LieGroup) = print(io, "LieGroup($(G.manifold), $(G.op))")
