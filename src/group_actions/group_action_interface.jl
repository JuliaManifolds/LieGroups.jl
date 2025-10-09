#
#
# Group Action related types and functions

@doc """
    AbstractGroupActionType

An abstract supertype for group action types, which are used within a [`GroupAction`](@ref).
"""
abstract type AbstractGroupActionType end

@doc """
    AbstractLeftGroupActionType <: AbstractGroupActionType

A type representing a (smooth) group action ``σ: $(_math(:G)) × $(_math(:M)) → $(_math(:M))``
of a [`AbstractLieGroup`](@ref) ``$(_math(:G))`` acting (from the left) on an $(_link(:AbstractManifold)) ``$(_math(:M))``.
with the following properties

1. ``σ($(_math(:e)), p) = p`` holds for all ``p ∈ $(_math(:M))``
2. ``σ(g, σ(h, p)) = σ(g$(_math(:∘))h, p)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``


where ``$(_math(:∘))`` denotes the group operation of the [`AbstractLieGroup`](@ref) ``$(_math(:G))``.
See also [HilgertNeeb:2012; Definition 9.1.11](@cite).

The type of action can be seen a bit better when writing the action as a family ``σ_g(p)``:
we obtain from the second property as

```math
  σ_g(σ_h(p)) = σ_{gh}(p)
```

and see that ``g`` appears on the left.

When writing about general group actions, the symbol ``α`` is often used. The order of arguments then follows the same as the one of the left action.
Most often we use the index notation ``α_g(p)``.
The definition of functions also follows this notation, i.e. we use e.g. `apply(A, g, p)`

One notable example of a left action is the inverse of an action of [`AbstractRightGroupActionType`](@ref) ``τ``,
which is given by ``σ_g = (τ_g)^{-1} = τ_{g^{-1}}``.
We obtain

$(_note(:RightInverseActionIsLeft))
"""
abstract type AbstractLeftGroupActionType <: AbstractGroupActionType end

@doc """
    AbstractRightGroupActionType <: AbstractGroupActionType

A type representing a (smooth) group action ``τ: $(_math(:M)) × $(_math(:G)) → $(_math(:M))``
of a [`AbstractLieGroup`](@ref) ``$(_math(:G))`` acting (from the right) on an $(_link(:AbstractManifold)) ``$(_math(:M))``.
with the following properties

1. ``τ(p, $(_math(:e))) = p`` holds for all ``p ∈ $(_math(:M))``
2. ``τ(τ(p, g), h) = τ(τ(p, g), h)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``

where ``$(_math(:∘))`` denotes the group operation of the [`AbstractLieGroup`](@ref) ``$(_math(:G))``.
See also [HilgertNeeb:2012; Remark 9.1.12](@cite).

The type of action can be seen a bit better when writing the action as a family ``τ_g(p)``:
we obtain from the second property as

```math
  τ_g(τ_h(p)) = τ_{hg}(p)
```

and see that ``g`` appears on the right.

When writing about general group actions, the symbol ``α`` is often used.
In that case the order of arguments follows either the one from the left action, but most often we use the index notation.
The definition of functions also follows this notation, i.e. we use e.g. `apply(A, g, p)`

One notable example of a right action is the inverse of an action of [`AbstractLeftGroupActionType`](@ref) ``σ``,
which is given by ``τ_g = (σ_g)^{-1} = σ_{g^{-1}}``.
We obtain

$(_note(:LeftInverseActionIsRight))
"""
abstract type AbstractRightGroupActionType <: AbstractGroupActionType end

"""
    AbstractActionActsOnType

An abstract type representing what an action acts on,
Most notably these are the [`ActionActsOnLeft`](@ref) and [`ActionActsOnRight`](@ref),
see their documentations for more details.

For its practical use see the [`SemidirectProductGroupOperation`](@ref).
"""
abstract type AbstractActionActsOnType end

"""
    ActionActsOnLeft <: AbstractActionActsOnType

An [`AbstractActionActsOnType`](@ref) representing that an action acts on the left.

This is meant in the following way: Given a [`GroupAction`](@ref) ``α: $(_math(:G)) × $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``
where a Lie group ``$(_math(:G))`` acts on another Lie group ``$(_tex(:Cal, "H"))`` with an arbitrary action.

Then, e.g. within the definition of the [`LeftSemidirectProductGroupOperation`](@ref) or [`RightSemidirectProductGroupOperation`](@ref),
we have two choices where the group action ``α`` acts on, namely:
Let ``g ∈ $(_math(:G))`` and ``h_1,h_2 ∈ $(_tex(:Cal, "H"))`` be given, then this type represents the variant

```math
α_g(h_1) ⋅ h_2,
```
where ``⋅`` denotes the group operation on ``$(_tex(:Cal, "H"))``.
The `Left` in the name of this type refers to the fact that the action is applied to the left element ``h_1``.

Note that this is independent of both the type of action (left or right) and whether the semidirect product is a left or a right semidirect one.

For its practical use see the [`SemidirectProductGroupOperation`](@ref).
"""
struct ActionActsOnLeft <: AbstractActionActsOnType end

"""
    ActionActsOnRight <: AbstractActionActsOnType

An [`AbstractActionActsOnType`](@ref) representing that an action acts on the right.

This is meant in the following way: Given a [`GroupAction`](@ref) ``α: $(_math(:G)) × $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``
where a Lie group ``$(_math(:G))`` acts on another Lie group ``$(_tex(:Cal, "H"))`` with an arbitrary action.

Then, e.g. within the definition of the [`LeftSemidirectProductGroupOperation`](@ref) or [`RightSemidirectProductGroupOperation`](@ref),
we have two choices where the group action ``α`` acts on, namely:
Let ``g ∈ $(_math(:G))`` and ``h_1,h_2 ∈ $(_tex(:Cal, "H"))`` be given, then this type represents the variant

```math
h_1 ⋅ α_g(h_2),
```
where ``⋅`` denotes the group operation on ``$(_tex(:Cal, "H"))``.
The `Right` in the name of this type refers to the fact that the action is applied to the right element ``h_2``.

Note that this is independent of both the type of action (left or right) and whether the semidirect product is a left or a right semidirect one.

For its practical use see the [`SemidirectProductGroupOperation`](@ref).
"""
struct ActionActsOnRight <: AbstractActionActsOnType end

"""
    GroupAction{T<:GroupActionType, L<:LieGroup, M<:AbstractManifold}

Specify a group action of [`AbstractGroupActionType`](@ref) `T` of a [`AbstractLieGroup`](@ref) `G`
acting on an $(_link(:AbstractManifold)) `M`.

Let ``$(_math(:M))`` be a $(_link(:AbstractManifold)) and ``$(_math(:G))`` be a [`AbstractLieGroup`](@ref)
with group operation ``$(_math(:∘))``.

A (smooth) action of the group ``$(_math(:G))`` on the manifold ``$(_math(:M))`` is a map

```math
α: $(_math(:G)) × $(_math(:M)) → $(_math(:M))
```

with the properties

**Identity.** ``α($(_math(:e)), p) = p`` holds for all ``p ∈ $(_math(:M))``

**Compatibility.**
If ``α`` is a [``AbstractLeftGroupActionType``](@ref) we usually denote it by ``σ`` and the compatibility reads

```math
σ_g(σ_h(p)) = σ_{g$(_math(:∘))h}(p) $(_tex(:text, " holds for all")) g,h ∈ $(_math(:G)) $(_tex(:text, " and")) p ∈ $(_math(:M))
```

If ``α`` is a [``AbstractRightGroupActionType``](@ref) we usually denote it by ``τ`` and the compatibility reads

```math
τ_g(τ_h(p)) = τ_{h$(_math(:∘))g}(p)`` holds for all ``g,h ∈ $(_math(:G)) $(_tex(:text, " for all")) p ∈ $(_math(:M))
```

_intuitively_ the left/right property of an action specifies on which side the “outer” group actions element ``g``
gets “appended” in the composition.

# Fields

* `type::`[`AbstractGroupActionType`](@ref): The type of the group action.
* `group::`[`AbstractLieGroup`](@ref): The group acting.
* `manifold::`$(_link(:AbstractManifold)): The manifold the group acts upon.

See [HilgertNeeb:2012; Section 9.1.3](@cite) for more details.

# Constructors
    GroupAction(
        group::AbstractLieGroup, manifold::ManifoldsBase.AbstractManifold, type::AbstractGroupActionType
    )
Generate a group action where the type of the action and what it acts on are keyword arguments.
They default to the most common choice, that the [`ActionActsOnRight`](@ref).
"""
struct GroupAction{
        T <: AbstractGroupActionType, L <: AbstractLieGroup, M <: ManifoldsBase.AbstractManifold,
    }
    type::T
    group::L
    manifold::M
end
function GroupAction(
        group::G, manifold::M, type::T
    ) where {G <: AbstractLieGroup, M <: ManifoldsBase.AbstractManifold, T <: AbstractGroupActionType}
    return GroupAction{T, G, M}(type, group, manifold)
end
#
#
# Functions

_doc_apply = """
    apply(A::GroupAction{T}, g, p)
    apply!(A::GroupAction{T}, q, g, p)

Apply the group action induced by ``g ∈ $(_math(:G))`` to ``p ∈ $(_math(:M))``,
where the kind of group action is indicated by the [`AbstractGroupActionType`](@ref) `T`.
This can be performed in-place of `q`.
"""

# function apply end
# un-comment the preceding line and remove this, once GroupManifolds no longer exists in Manifolds.jl
@doc "$(_doc_apply)"
function apply(A::GroupAction, g, p)
    q = allocate_result(base_manifold(A), apply, p)
    apply!(A, q, g, p)
    return q
end

# Define `function apply! end` here as well
# un-comment (remove this comment) when removing this function from Manifolds.jl
@doc "$(_doc_apply)"
apply!(A::GroupAction, q, g, p)

function base_lie_group end
@doc """
    base_lie_group(A::GroupAction)

Return the [`AbstractLieGroup`](@ref) of the [`GroupAction`](@ref)
specifying the action.
"""
base_lie_group(A::GroupAction) = A.group

@doc """
    base_manifold(A::GroupAction)

Return the $(_link(:AbstractManifold)) the group action acts upon.
"""
ManifoldsBase.base_manifold(A::GroupAction) = A.manifold

function default_left_action end
"""
    default_left_action(G::AbstractLieGroup, M::AbstractManifold)

Return the default left action for a Lie group `G` acting on a manifold `M`.
"""
default_left_action(N::AbstractLieGroup, M::AbstractManifold)

function default_right_action end
"""
    default_right_action(G::AbstractLieGroup, M::AbstractManifold)

Return the default right action for a Lie group `G` acting on a manifold `M`.
"""
default_right_action(N::AbstractLieGroup, M::AbstractManifold)

_doc_diff_apply = """
    diff_apply(A::GroupAction, g, p, X)
    diff_apply!(A::GroupAction, Y, g, p, X)

Compute the differential ``$(_math(:D))_p α_g(p): T_p$(_math(:M)) → T_{σ_g(p)}$(_math(:M))``.
"""

function diff_apply end
@doc "$(_doc_diff_apply)"
function diff_apply(A::GroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), diff_apply, p, X)
    diff_apply!(A, Y, g, p, X)
    return Y
end

function diff_apply! end
@doc "$(_doc_diff_apply)"
diff_apply!(A::GroupAction, q, g, p)

_doc_diff_group_apply = raw"""
    diff_group_apply(A::GroupAction, g, p, X)
    diff_group_apply!(A::GroupAction, Y, g, p, X)


Compute the value of differential of action [`GroupAction`](@ref) `A` on vector `X`,
where element `g` is acting on `p`, with respect to the group element.

Let ``\mathcal G`` be the group acting on manifold ``\mathcal M`` by the action `A`.
The action is of element ``g ∈ \mathcal G`` on a point ``p ∈ \mathcal M``.
The differential transforms vector `X` from the tangent space at `a ∈ \mathcal G`,
``X ∈ T_a \mathcal G`` into a tangent space of the manifold ``\mathcal M``.
When action on element `p` is written as ``\mathrm{d}τ^p``, with the specified left or right
convention, the differential transforms vectors

````math
(\mathrm{d}τ^p) : T_{a} \mathcal G → T_{τ_a p} \mathcal M
````

# See also

[`apply`](@ref), [`diff_apply`](@ref)
"""

function diff_group_apply end
@doc "$(_doc_diff_group_apply)"
function diff_group_apply(A::GroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), diff_group_apply, g, X)
    diff_group_apply!(A, Y, g, p, X)
    return Y
end
function diff_group_apply(A::GroupAction, g::Identity, p, X)
    Y = allocate_on(A.manifold, TangentSpaceType())
    diff_group_apply!(A, Y, g, p, X)
    return Y
end

function diff_group_apply! end
@doc "$(_doc_diff_group_apply)"
diff_group_apply!(A::GroupAction, q, g, p)

@doc """
     inv(A::GroupAction{T})

Return the [`GroupAction`](@ref) representing the inverse of an [`GroupAction`](@ref) of [`AbstractGroupActionType`](@ref) `T`.
This is usually done by returning the group action with the inverse type of `T`.
"""
Base.inv(A::GroupAction) = GroupAction(A.group, A.manifold, inv(A.type))

"""
    inv(::AbstractGroupActionType)

return the inverse group operation action, that is, use the type representing the
inverse operation.
"""
Base.inv(::AbstractGroupActionType)

function Base.show(io::IO, A::GroupAction)
    return print(io, "GroupAction($(A.group), $(A.manifold), $(A.type))")
end

function switch end

@doc """
     switch(A::GroupAction{T})

Return the group operation action representing the similar [`GroupAction`](@ref) of [`AbstractGroupActionType`](@ref) `T`
but acting from the other side. It switches left to right and vice versa.
This is done by returning the group action with the “switched” type of `T`.
"""
switch(A::GroupAction) = GroupAction(A.group, A.manifold, switch(A.type))

@doc """
    switch(T::AbstractGroupActionType)

Return the object representing an [`AbstractGroupActionType`](@ref) related to a group operation action that switched the side,
that is it turns a left action type into its corresponding right action type.
"""
switch(::AbstractGroupActionType)
