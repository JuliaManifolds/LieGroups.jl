"""
    LeftGroupOperation <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
from the left, that is

```math
σ_h(g) = σ(h,g) = h$(_math(:∘))g
```

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperation`](@ref).
"""
struct LeftGroupOperation <: AbstractLeftGroupActionType end

"""
    RightGroupOperation <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
gfrom the right

```math
σ_h(g) = σ(h,g) = g$(_math(:∘))h
```

for its inverse ``(σ_h)^{-1}`` see [`InverseRightGroupOperation`](@ref).
"""
struct RightGroupOperation <: AbstractLeftGroupActionType end

"""
    InverseLeftGroupOperation <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`LeftGroupOperation`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = h^{-1}$(_math(:∘))g
```

Note that while in naming it is the inverse of the left action, it's
properties yield that is is an [`AbstractRightGroupActionType`](@ref), since

$(_note(:LeftInverseActionIsRight))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperation`](@ref).

!!! note
    Some literature also calls this by itself _the_ right group operation action.
"""
struct InverseLeftGroupOperation <: AbstractRightGroupActionType end

"""
    InverseRightGroupOperation <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`RightGroupOperation`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = g$(_math(:∘))h^{-1}
```

Note that while in naming it is the inverse of the right action, it's
properties yield that is is an [`AbstractLeftGroupActionType`](@ref), since

$(_note(:RightInverseActionIsLeft))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperation`](@ref).
"""
struct InverseRightGroupOperation <: AbstractLeftGroupActionType end

"""
    GroupOperationAction{T<:AbstractLeftGroupActionType,G<:LieGroup} <: AbstractGroupAction{T,G,G}

The [`AbstractGroupAction`](@ref).

"""
struct GroupOperationAction{T<:AbstractGroupActionType,G<:LieGroup} <:
       AbstractGroupAction{T,G,G}
    type::T
    group::G
end

base_lie_group(A::GroupOperationAction) = A.group

ManifoldsBase.base_manifold(A::GroupOperationAction) = A.group

_doc_apply_groupop = """
    apply(A::GroupOperationAction, g, h)
    apply!(A::GroupOperationAction, k, g, h)

apply the stored group operation action, using [`compose`](@ref) ``$(_math(:∘))``.
this can be done in-place of `k`.
"""

@doc "$(_doc_apply_groupop)"
apply(A::GroupOperationAction, g, h)

@doc "$(_doc_apply_groupop)"
apply!(A::GroupOperationAction, k, g, h)

function apply!(A::GroupOperationAction{LeftGroupOperation}, k, g, h)
    return compose!(A.group, k, g, h) #apply/compose g from left
end
function apply!(A::GroupOperationAction{RightGroupOperation}, k, g, h)
    return compose!(A.group, k, h, g) #apply/compose g from right
end
function apply!(A::GroupOperationAction{InverseLeftGroupOperation}, k, g, h)
    return inv_left_compose!(A.group, k, g, h) #apply/compose inv(g) from left
end
function apply!(A::GroupOperationAction{InverseRightGroupOperation}, k, g, h)
    return inv_right_compose!(A.group, k, g, h) #apply/compose inv(g) from right
end

_doc_diff_apply_groupop = """

    diff_apply(A::GroupOperationAction, g, p, X)

For the group operation action ``σ_g(p)``, compute the differential
``D_p σ_g(p): T_p$(_math(:G)) → T_{σ_g(p)}$(_math(:M))``, that is

* for the [`LeftGroupOperation`](@ref) this calls [`diff_right_compose`](@ref)`(G, g, p, X)`, since here ``σ_g(p) = g$(_math(:∘))p``
* for the [`RightGroupOperation`](@ref) this calls [`diff_left_compose`](@ref)`(G, p, g, X)`, since here ``σ_g(p) = p$(_math(:∘))g``
* for the [`InverseLeftGroupOperation`](@ref) this calls [`diff_right_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = g^{-1}$(_math(:∘))p``
* for the [`InverseRightGroupOperation`](@ref) this calls [`diff_left_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = p$(_math(:∘))g^{-1}``
"""
@doc "$(_doc_diff_apply_groupop)"
diff_apply(A::GroupOperationAction, g, p, X)

@doc "$(_doc_diff_apply_groupop)"
diff_apply!(A::GroupOperationAction, Y, g, p, X)

function diff_apply!(A::GroupOperationAction{LeftGroupOperation}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, g, p, X)
end
function diff_apply!(A::GroupOperationAction{RightGroupOperation}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, g, X)
end
function diff_apply!(A::GroupOperationAction{InverseLeftGroupOperation}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, inv(A.group, g), p, X)
end
function diff_apply!(A::GroupOperationAction{InverseRightGroupOperation}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, inv(A.group, g), X)
end

_doc_diff_group_apply_groupop = """

    diff_group_apply(A::GroupOperationAction, g, p, X)

Compute the differential ``D_g σ_g(p): $(_math(:𝔤)) → $(_math(:𝔤))`` of a group operation action,
that is

* for the [`LeftGroupOperation`](@ref) this calls [`diff_left_compose`](@ref)`(G, g, p, X)`, since here ``σ_g(p) = g$(_math(:∘))p``
* for the [`RightGroupOperation`](@ref) this calls [`diff_right_compose`](@ref)`(G, p, g, X)`, since here ``σ_g(p) = p$(_math(:∘))g``
* for the [`InverseLeftGroupOperation`](@ref) this calls [`diff_left_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = g^{-1}$(_math(:∘))p`` together with [`diff_inv`](@ref)
* for the [`InverseRightGroupOperation`](@ref) this calls [`diff_right_compose`](@ref) with ``g^{-1}``, since here ``σ_g(p) = p$(_math(:∘))g^{-1}`` together with [`diff_inv`](@ref)
"""

@doc "$(_doc_diff_apply_groupop)"
diff_group_apply(A::GroupOperationAction, g, p, X)

@doc "$(_doc_diff_apply_groupop)"
diff_group_apply!(A::GroupOperationAction, Y, g, p, X)

function diff_group_apply!(A::GroupOperationAction{LeftGroupOperation}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, g, p, X)
end
function diff_group_apply!(A::GroupOperationAction{RightGroupOperation}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, p, g, X)
end
function diff_group_apply!(A::GroupOperationAction{InverseLeftGroupOperation}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_left_compose!(A.group, Y, inv(A.group, g), p, Y)
end
function diff_group_apply!(A::GroupOperationAction{InverseRightGroupOperation}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_right_compose!(A.group, Y, p, inv(A.group, g), Y)
end

"""
    inv(::GroupOperationAction)

return the inverse group operation action, that is, use the type representing the
inverse operation.
"""
Base.inv(A::GroupOperationAction) = GroupOperationAction(inv(A.type), A.group)

"""
    inv(::LeftGroupOperation)

Return the inverse of the [`LeftGroupOperation`](@ref), that is the [`InverseLeftGroupOperation`](@ref).
"""
Base.inv(::LeftGroupOperation) = InverseLeftGroupOperation()
"""
    inv(::RightGroupOperation)

Return the inverse of the [`RightGroupOperation`](@ref), that is the [`InverseRightGroupOperation`](@ref).
"""
Base.inv(::RightGroupOperation) = InverseRightGroupOperation()
"""
    inv(::InverseLeftGroupOperation)

Return the inverse of the [`InverseLeftGroupOperation`](@ref), that is the [`LeftGroupOperation`](@ref).
"""
Base.inv(::InverseLeftGroupOperation) = LeftGroupOperation()
"""
    inv(::InverseRightGroupOperation)

Return the inverse of the [`InverseRightGroupOperation`](@ref), that is the [`RightGroupOperation`](@ref).
"""
Base.inv(::InverseRightGroupOperation) = RightGroupOperation()

function Base.show(io::IO, A::GroupOperationAction)
    return print(io, "GroupOperationAction($(A.type),$(A.group))")
end

"""
    switch(::LeftGroupOperation)

Return the [`RightGroupOperation`](@ref), that is,
turns ``σ_g = g$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g``
"""
switch(::LeftGroupOperation) = RightGroupOperation()

"""
    switch(::RightGroupOperation)

Return the [`LeftGroupOperation`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g`` into ``τ_g(h) = g$(_math(:∘))h``
"""
switch(::RightGroupOperation) = LeftGroupOperation()

"""
    switch(::InverseLeftGroupOperation)

Return the [`InverseRightGroupOperation`](@ref), that is,
turns ``σ_g = g^{-1}$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g^{-1}``
"""
switch(::InverseLeftGroupOperation) = InverseRightGroupOperation()

"""
    switch(::InverseRightGroupOperation)

Return the [`InverseLeftGroupOperation`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g^{-1}`` into ``τ_g(h) = g^{-1}$(_math(:∘))h``
"""
switch(::InverseRightGroupOperation) = InverseLeftGroupOperation()
