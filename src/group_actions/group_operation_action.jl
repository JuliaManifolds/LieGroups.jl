"""
    LeftGroupOperationAction <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
from the left, that is

```math
σ_h(g) = σ(h,g) = h$(_math(:∘))g
```

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationAction`](@ref).
"""
struct LeftGroupOperationAction <: AbstractLeftGroupActionType end

"""
    RightGroupOperationAction <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
from the right.

```math
σ_h(g) = σ(h,g) = g$(_math(:∘))h
```

for its inverse ``(σ_h)^{-1}`` see [`InverseRightGroupOperationAction`](@ref).
"""
struct RightGroupOperationAction <: AbstractLeftGroupActionType end

"""
    InverseLeftGroupOperationAction <: AbstractRightGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`LeftGroupOperationAction`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = h^{-1}$(_math(:∘))g
```

Note that while in naming it is the inverse of the left action, it's
properties yield that is is an [`AbstractRightGroupActionType`](@ref), since

$(_note(:LeftInverseActionIsRight))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationAction`](@ref).

!!! note
    Some literature also calls this by itself _the_ right group operation action.
"""
struct InverseLeftGroupOperationAction <: AbstractRightGroupActionType end

"""
    InverseRightGroupOperationAction <: AbstractLeftGroupActionType

A type for the [`AbstractLeftGroupActionType`](@ref) when acting on the group itself
given by the inverse of a [`RightGroupOperationAction`](@ref) ``σ_h`` as

```math
τ_h(g) $(_tex(:def)) σ_h^{-1}(g) = σ(h^{-1},g) = g$(_math(:∘))h^{-1}
```

Note that while in naming it is the inverse of the right action, it's
properties yield that is is an [`AbstractLeftGroupActionType`](@ref), since

$(_note(:RightInverseActionIsLeft))

for its inverse ``(σ_h)^{-1}`` see [`InverseLeftGroupOperationAction`](@ref).
"""
struct InverseRightGroupOperationAction <: AbstractLeftGroupActionType end

"""
    GroupOperationAction(action::AbstractGroupActionType, group::LieGroup)

Return a [`GroupAction`](@ref) for an [`AbstractGroupActionType`](@ref) `action`
representing the group operation as an action of the group on itself.
"""
function GroupOperationAction(action::AbstractGroupActionType, G::LieGroup)
    return GroupAction(action, G, G)
end
function Base.show(
    io::IO, GA::GroupAction{A,G,G}
) where {A<:AbstractGroupActionType,G<:LieGroup}
    return print(io, "GroupOperationAction($(GA.type), $(GA.group))")
end

function apply!(A::GroupAction{LeftGroupOperationAction}, k, g, h)
    return compose!(A.group, k, g, h) #apply/compose g from left
end
function apply!(A::GroupAction{RightGroupOperationAction}, k, g, h)
    return compose!(A.group, k, h, g) #apply/compose g from right
end
function apply!(A::GroupAction{InverseLeftGroupOperationAction}, k, g, h)
    return inv_left_compose!(A.group, k, g, h) #apply/compose inv(g) from left
end
function apply!(A::GroupAction{InverseRightGroupOperationAction}, k, g, h)
    return inv_right_compose!(A.group, k, g, h) #apply/compose inv(g) from right
end

function diff_apply!(A::GroupAction{LeftGroupOperationAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, g, p, X)
end
function diff_apply!(A::GroupAction{RightGroupOperationAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, g, X)
end
function diff_apply!(A::GroupAction{InverseLeftGroupOperationAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, inv(A.group, g), p, X)
end
function diff_apply!(A::GroupAction{InverseRightGroupOperationAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, p, inv(A.group, g), X)
end

function diff_group_apply!(A::GroupAction{LeftGroupOperationAction}, Y, g, p, X)
    return diff_left_compose!(A.group, Y, g, p, X)
end
function diff_group_apply!(A::GroupAction{RightGroupOperationAction}, Y, g, p, X)
    return diff_right_compose!(A.group, Y, p, g, X)
end
function diff_group_apply!(A::GroupAction{InverseLeftGroupOperationAction}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_left_compose!(A.group, Y, inv(A.group, g), p, Y)
end
function diff_group_apply!(A::GroupAction{InverseRightGroupOperationAction}, Y, g, p, X)
    diff_inv!(A.group, Y, g, X)
    return diff_right_compose!(A.group, Y, p, inv(A.group, g), Y)
end

"""
    inv(::LeftGroupOperationAction)

Return the inverse of the [`LeftGroupOperationAction`](@ref), that is the [`InverseLeftGroupOperationAction`](@ref).
"""
Base.inv(::LeftGroupOperationAction) = InverseLeftGroupOperationAction()
"""
    inv(::RightGroupOperationAction)

Return the inverse of the [`RightGroupOperationAction`](@ref), that is the [`InverseRightGroupOperationAction`](@ref).
"""
Base.inv(::RightGroupOperationAction) = InverseRightGroupOperationAction()
"""
    inv(::InverseLeftGroupOperationAction)

Return the inverse of the [`InverseLeftGroupOperationAction`](@ref), that is the [`LeftGroupOperationAction`](@ref).
"""
Base.inv(::InverseLeftGroupOperationAction) = LeftGroupOperationAction()
"""
    inv(::InverseRightGroupOperationAction)

Return the inverse of the [`InverseRightGroupOperationAction`](@ref), that is the [`RightGroupOperationAction`](@ref).
"""
Base.inv(::InverseRightGroupOperationAction) = RightGroupOperationAction()

"""
    switch(::LeftGroupOperationAction)

Return the [`RightGroupOperationAction`](@ref), that is,
turns ``σ_g = g$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g``
"""
switch(::LeftGroupOperationAction) = RightGroupOperationAction()

"""
    switch(::RightGroupOperationAction)

Return the [`LeftGroupOperationAction`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g`` into ``τ_g(h) = g$(_math(:∘))h``
"""
switch(::RightGroupOperationAction) = LeftGroupOperationAction()

"""
    switch(::InverseLeftGroupOperationAction)

Return the [`InverseRightGroupOperationAction`](@ref), that is,
turns ``σ_g = g^{-1}$(_math(:∘))h`` into ``τ_g(h) = h$(_math(:∘))g^{-1}``
"""
switch(::InverseLeftGroupOperationAction) = InverseRightGroupOperationAction()

"""
    switch(::InverseRightGroupOperationAction)

Return the [`InverseLeftGroupOperationAction`](@ref), that is,
turns ``σ_g = h$(_math(:∘))g^{-1}`` into ``τ_g(h) = g^{-1}$(_math(:∘))h``
"""
switch(::InverseRightGroupOperationAction) = InverseLeftGroupOperationAction()
