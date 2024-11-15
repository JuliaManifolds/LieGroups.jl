#
#
# Group Action related types and functions

@doc """
    AbstractGroupActionType

An abstract supertype for group action types.
"""
abstract type AbstractGroupActionType end

_note_action_argument_order = """
!!! note
    In function definitions we follow the idea of the family of actions and use the order `(M, g, p)` in function signatures.
"""

@doc """
    AbstractLeftGroupActionType <: AbstractGroupActionType

A type representing a (smooth) group action ``σ: $(_math(:G)) × $(_math(:M)) → $(_math(:M))``
of a [`LieGroup`](@ref) ``$(_math(:G))`` acting (from the left) on an $(_link(:AbstractManifold)) ``$(_math(:M))``.
with the following properties

1. ``σ($(_math(:e)), p) = p`` holds for all ``p ∈ $(_math(:M))``
2. ``σ(g, σ(h, p)) = σ(g$(_math(:∘))h, p)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``


where ``$(_math(:∘))`` denotes the group operation of the [`LieGroup`](@ref) ``$(_math(:G))``.
See also [HilgertNeeb:2012; Definition 9.1.11](@cite).

The type of action can be seen a bit better when writing the action as a family ``σ_g(p)``:
we obtain from the second property as

```math
  σ_g(σ_h(p)) = σ_{gh}(p)
```

and see that ``g`` appears on the left.

To emphasize the side the group operation is acting from, we sometimes write ``σ^{$(_tex(:rm,"L"))}``.
If the action is clear from context we write ``σ(g, p) = g $(_math(:act)) p``.

One notable example of a left action is the inverse of an action of [`AbstractRightGroupActionType`](@ref) ``σ^{$(_tex(:rm,"R"))}``,
which is given by ``τ_g = (σ^{$(_tex(:rm,"R"))}_g)^{-1} = σ^{$(_tex(:rm,"R"))}_{g^{-1}}``.
We obtain

$(_note(:RightInverseActionIsLeft))

$(_note_action_argument_order)
"""
abstract type AbstractLeftGroupActionType <: AbstractGroupActionType end

@doc """
    AbstractRightGroupActionType <: AbstractGroupActionType

A type representing a (smooth) group action ``σ: $(_math(:M)) × $(_math(:G)) → $(_math(:M))``
of a [`LieGroup`](@ref) ``$(_math(:G))`` acting (from the right) on an $(_link(:AbstractManifold)) ``$(_math(:M))``.
with the following properties

1. ``σ(p, $(_math(:e))) = p`` holds for all ``p ∈ $(_math(:M))``
2. ``σ(σ(p, g), h) = σ(p, g$(_math(:∘))h)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``

where ``$(_math(:∘))`` denotes the group operation of the [`LieGroup`](@ref) ``$(_math(:G))``.
See also [HilgertNeeb:2012; Remark 9.1.12](@cite).

The type of action can be seen a bit better when writing the action as a family ``σ_g(p)``:
we obtain from the second property as

```math
  σ_g(σ_h(p)) = σ_{hg}(p)
```

and see that ``g`` appears on the right.

To emphasize the side the group operation is acting from, we sometimes write ``σ^{$(_tex(:rm,"R"))}``.
If the action is clear from context we write ``σ(p, g) = p $(_math(:act)) g``.

One notable example of a right action is the inverse of an action of  [`AbstractLeftGroupActionType`](@ref) ``σ^{$(_tex(:rm,"L"))}``,
which is given by ``τ_g = (σ^{$(_tex(:rm,"L"))}_g)^{-1} = σ^{$(_tex(:rm,"L"))}_{g^{-1}}``.
We obtain

$(_note(:LeftInverseActionIsRight))

$(_note_action_argument_order)
"""
abstract type AbstractRightGroupActionType <: AbstractGroupActionType end

"""
    LieGroupOperationAction{T<:AbstractLeftGroupActionType,G<:LieGroup} <: AbstractGroupAction{T,G,G}

A group action of [`AbstractGroupActionType`](@ref) `T` of a [`LieGroup`](@ref) of type `L`
acting on an $(_link(:AbstractManifold)) of type `M`.

# Fields

* `type::T`: The type of the group action.
* `group::L`: The group acting.
* `manifold::M`: The manifold the group acts upon.

See [HilgertNeeb:2012; Section 9.1.3](@cite) for more details.

"""
struct GroupAction{T<:AbstractGroupActionType, L<:LieGroup, M<:Manifold}
    type::T
    group::L
    manifold::M
end



function base_lie_group end
@doc """
    base_lie_group(A::AbstractGroupAction)

Return the [`LieGroup`](@ref) of the [`AbstractGroupAction`](@ref)
specifying the action.
"""
base_lie_group(A::GroupAction) = A.group

@doc """
    base_manifold(A::AbstractGroupAction)

Return the $(_link(:AbstractManifold)) the group action acts upon.
"""
ManifoldsBase.base_manifold(A::GroupOperationAction) = A.manifold

#
#
# Functions

_doc_apply = """
    apply(A::AbstractGroupAction{T, L, M}, g, p)
    apply!(A::AbstractGroupAction{T, L, M}, q, g, p)

Apply the group action induced by ``g ∈ $(_math(:G))`` to ``p ∈ $(_math(:M))``,
where the kind of group action is indicated by the [`AbstractGroupActionType`](@ref) `T`.
This can be perfomed in-place of `q`.
"""

# function apply end
# un-comment the preceding line and remove this, once GroupManifolds no longer exists in Manifolds.jl
@doc "$(_doc_apply)"
function apply(A::AbstractGroupAction, g, p)
    q = allocate_result(base_manifold(A), apply, g, p)
    apply!(A, q, g, p)
    return q
end

# Define `function apply! end` here as well
# un-comment (remove this comment) when removing this function from Manifolds.jl
@doc "$(_doc_apply)"
apply!(A::AbstractGroupAction, q, g, p)

_doc_diff_apply = """
    diff_apply(A::AbstractGroupAction{T, L, M}, g, p, X)
    diff_apply!(A::AbstractGroupAction{T, L, M}, Y, g, p, X)

Compute the differential ``D_p σ_g(p): T_p$(_math(:M)) → T_{σ_g(p)}$(_math(:M))``,
where for a left group action we have ``σ_g(p) = σ(g,p)``, for a right action ``σ_g(p) = σ(p, g)``.
"""

function diff_apply end
@doc "$(_doc_diff_apply)"
function diff_apply(A::AbstractGroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), apply_diff_group, p, g, X)
    diff_apply!(A, Y, g, p, X)
    return Y
end

function diff_apply! end
@doc "$(_doc_diff_apply)"
diff_apply!(A::AbstractGroupAction, q, g, p)

_doc_diff_group_apply = """
    diff_group_apply(A::AbstractGroupAction{T, L, M}, g, p, X)
    diff_group_apply!(A::AbstractGroupAction{T, L, M}, Y, g, p, X)

Compute the differential ``D_g σ_g(p): $(_math(:𝔤)) → $(_math(:𝔤))``,
where we use the short hand notation ``σ_p(g) = σ(g,p)`` for a left action,
and for a right action ``σ_p(g) = σ(p, g)``.
"""

function diff_group_apply end
@doc "$(_doc_diff_group_apply)"
function diff_group_apply(A::AbstractGroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), apply, g, p, X)
    diff_group_apply!(A, Y, g, p, X)
    return Y
end

function diff_group_apply! end
@doc "$(_doc_diff_group_apply)"
diff_group_apply!(A::AbstractGroupAction, q, g, p)

@doc """
     inv(A::AbstractGroupAction{T})

Return the tuple representing the inverse of an [`AbstractGroupAction`](@ref) of [`AbstractGroupActionType`](@ref) `T`.
This is usually done by returning the group action with the inverse type of `T`.
"""
Base.inv(::AbstractGroupAction)

@doc """
     inv(T::AbstractGroupActionType)

Return the type representing the inverse of an [`AbstractGroupActionType`](@ref).
"""
Base.inv(::AbstractGroupActionType)

function switch end

@doc """
     switch(A::AbstractGroupAction{T})

Return the group operation action representing the similar [`AbstractGroupAction`](@ref) of [`AbstractGroupActionType`](@ref) `T`
but acting from the other side. It switches left to right and vice versa.
This is done by returning the group action with the “switched” type of `T`.
"""
switch(::AbstractGroupAction)

@doc """
    switch(T::AbstractGroupActionType)

Return the object representing an [`AbstractGroupActionType`](@ref) related to a group operation action that switched the side,
that is it turns a left action type into its corresponding right action type.
"""
switch(::AbstractGroupActionType)
