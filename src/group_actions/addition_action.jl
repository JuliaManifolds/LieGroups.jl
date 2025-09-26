"""
    AdditionGroupAction <: AbstractLeftGroupActionType

Specify that in a [`GroupAction`](@ref) with Lie group ``$(_tex(:Cal, "G"))``
and manifold ``$(_tex(:Cal, "M"))``, where the group action is given by addition

Given an element ``g ∈ $(_tex(:Cal, "G"))`` and a point ``p ∈ $(_tex(:Cal, "M"))``,
the family of actions ``σ_g: $(_tex(:Cal, "M")) → $(_tex(:Cal, "M"))`` is given by

```math
σ_g(p) = g + p
```

where ``+`` denotes the vector or matrix addition.
"""
struct AdditionGroupAction <: AbstractLeftGroupActionType end

function apply!(::GroupAction{AdditionGroupAction}, q, g, p)
    return q .= g .+ p
end

function diff_apply!(::GroupAction{AdditionGroupAction}, Y, g, p, X)
    return Y .= X
end

function diff_group_apply!(::GroupAction{AdditionGroupAction}, Y, g, p, X)
    return Y .= X
end
