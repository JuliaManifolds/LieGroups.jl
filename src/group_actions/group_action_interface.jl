#
#
# Group Action related types and functions

@doc """
    AbstractGroupActionType

An abstract supertype for group action types, which are used within a [`GroupAction`](@ref).
"""
abstract type AbstractGroupActionType end

_note_action_argument_order = """
!!! note
    In function definitions we follow the idea of the family of actions and use the order `(M, g, p)` in function signatures.
"""

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

To emphasize the side the group operation is acting from, we sometimes write ``σ^{$(_tex(:rm, "L"))}``.
If the action is clear from context we write ``σ(g, p) = g $(_math(:act)) p``.

One notable example of a left action is the inverse of an action of [`AbstractRightGroupActionType`](@ref) ``σ^{$(_tex(:rm, "R"))}``,
which is given by ``τ_g = (σ^{$(_tex(:rm, "R"))}_g)^{-1} = σ^{$(_tex(:rm, "R"))}_{g^{-1}}``.
We obtain

$(_note(:RightInverseActionIsLeft))

$(_note_action_argument_order)
"""
abstract type AbstractLeftGroupActionType <: AbstractGroupActionType end

@doc """
    AbstractRightGroupActionType <: AbstractGroupActionType

A type representing a (smooth) group action ``σ: $(_math(:M)) × $(_math(:G)) → $(_math(:M))``
of a [`AbstractLieGroup`](@ref) ``$(_math(:G))`` acting (from the right) on an $(_link(:AbstractManifold)) ``$(_math(:M))``.
with the following properties

1. ``σ(p, $(_math(:e))) = p`` holds for all ``p ∈ $(_math(:M))``
2. ``σ(σ(p, g), h) = σ(p, g$(_math(:∘))h)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``

where ``$(_math(:∘))`` denotes the group operation of the [`AbstractLieGroup`](@ref) ``$(_math(:G))``.
See also [HilgertNeeb:2012; Remark 9.1.12](@cite).

The type of action can be seen a bit better when writing the action as a family ``σ_g(p)``:
we obtain from the second property as

```math
  σ_g(σ_h(p)) = σ_{hg}(p)
```

and see that ``g`` appears on the right.

To emphasize the side the group operation is acting from, we sometimes write ``σ^{$(_tex(:rm, "R"))}``.
If the action is clear from context we write ``σ(p, g) = p $(_math(:act)) g``.

One notable example of a right action is the inverse of an action of  [`AbstractLeftGroupActionType`](@ref) ``σ^{$(_tex(:rm, "L"))}``,
which is given by ``τ_g = (σ^{$(_tex(:rm, "L"))}_g)^{-1} = σ^{$(_tex(:rm, "L"))}_{g^{-1}}``.
We obtain

$(_note(:LeftInverseActionIsRight))

$(_note_action_argument_order)
"""
abstract type AbstractRightGroupActionType <: AbstractGroupActionType end

"""
    GroupAction{T<:GroupActionType, L<:LieGroup, M<:AbstractManifold}

Specify a group action of [`AbstractGroupActionType`](@ref) `T` of a [`AbstractLieGroup`](@ref) `G` acting on `M`.

Let ``$(_math(:M))`` be a $(_link(:AbstractManifold)) and ``$(_math(:G))`` be a [`AbstractLieGroup`](@ref)
with group operation ``$(_math(:∘))``.

A (smooth) action of the group ``$(_math(:G))`` on the manifold ``$(_math(:M))`` is a map

```math
σ: $(_math(:G)) × $(_math(:M)) → $(_math(:M))
```

with the properties

* ``σ($(_math(:e)), p) = p`` holds for all ``p ∈ $(_math(:M))``
* ``σ(g, σ(h, p)) = σ(g$(_math(:∘))h, p)`` holds for all ``g,h ∈ $(_math(:G))``, ``p ∈ $(_math(:M))``


# Fields

* `type::T`: The type of the group action.
* `group::L`: The group acting.
* `manifold::M`: The manifold the group acts upon.

See [HilgertNeeb:2012; Section 9.1.3](@cite) for more details.

"""
struct GroupAction{
        T <: AbstractGroupActionType, L <: AbstractLieGroup, M <: ManifoldsBase.AbstractManifold,
    }
    type::T
    group::L
    manifold::M
end

#
#
# Functions

_doc_apply = """
    apply(A::GroupAction{T, L, M}, g, p)
    apply!(A::GroupAction{T, L, M}, q, g, p)

Apply the group action induced by ``g ∈ $(_math(:G))`` to ``p ∈ $(_math(:M))``,
where the kind of group action is indicated by the [`AbstractGroupActionType`](@ref) `T`.
This can be performed in-place of `q`.
"""

# function apply end
# un-comment the preceding line and remove this, once GroupManifolds no longer exists in Manifolds.jl
@doc "$(_doc_apply)"
function apply(A::GroupAction, g, p)
    q = allocate_result(base_manifold(A), apply, g, p)
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
    diff_apply(A::GroupAction{T, L, M}, g, p, X)
    diff_apply!(A::GroupAction{T, L, M}, Y, g, p, X)

Compute the differential ``D_p σ_g(p): T_p$(_math(:M)) → T_{σ_g(p)}$(_math(:M))``,
where for a left group action we have ``σ_g(p) = σ(g,p)``, for a right action ``σ_g(p) = σ(p, g)``.
"""

function diff_apply end
@doc "$(_doc_diff_apply)"
function diff_apply(A::GroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), diff_apply, p, g, X)
    diff_apply!(A, Y, g, p, X)
    return Y
end

function diff_apply! end
@doc "$(_doc_diff_apply)"
diff_apply!(A::GroupAction, q, g, p)

_doc_diff_group_apply = """
    diff_group_apply(A::GroupAction{T, L, M}, g, p, X)
    diff_group_apply!(A::GroupAction{T, L, M}, Y, g, p, X)

Compute the differential ``D_g σ_g(p): $(_math(:𝔤)) → $(_math(:𝔤))``,
where we use the short hand notation ``σ_p(g) = σ(g,p)`` for a left action,
and for a right action ``σ_p(g) = σ(p, g)``.
"""

function diff_group_apply end
@doc "$(_doc_diff_group_apply)"
function diff_group_apply(A::GroupAction, g, p, X)
    Y = allocate_result(base_manifold(A), apply, g, p, X)
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
Base.inv(A::GroupAction) = GroupAction(inv(A.type), A.group, A.manifold)

"""
    inv(::AbstractGroupActionType)

return the inverse group operation action, that is, use the type representing the
inverse operation.
"""
Base.inv(::AbstractGroupActionType)

function Base.show(io::IO, A::GroupAction)
    return print(io, "GroupAction($(A.type), $(A.group), $(A.manifold))")
end

function switch end

@doc """
     switch(A::GroupAction{T})

Return the group operation action representing the similar [`GroupAction`](@ref) of [`AbstractGroupActionType`](@ref) `T`
but acting from the other side. It switches left to right and vice versa.
This is done by returning the group action with the “switched” type of `T`.
"""
switch(A::GroupAction) = GroupAction(switch(A.type), A.group, A.manifold)

@doc """
    switch(T::AbstractGroupActionType)

Return the object representing an [`AbstractGroupActionType`](@ref) related to a group operation action that switched the side,
that is it turns a left action type into its corresponding right action type.
"""
switch(::AbstractGroupActionType)
