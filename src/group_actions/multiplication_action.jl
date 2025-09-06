"""
    LeftMultiplicationGroupAction <: AbstractLeftGroupActionType

Specify that in a [`GroupAction`](@ref) with Lie group ``$(_tex(:Cal, "G"))``
and manifold ``$(_tex(:Cal, "M"))``, the group action is given by multiplication
from the left:

Given an element ``g ∈ $(_tex(:Cal, "G"))`` and a point ``p ∈ $(_tex(:Cal, "M"))``,
the family of actions ``σ_g: $(_tex(:Cal, "M")) → $(_tex(:Cal, "M"))`` is given by

```math
σ_g(p) = g*p
```

where ``*`` denotes the matrix(-vector) multiplication.
"""
struct LeftMultiplicationGroupAction <: AbstractLeftGroupActionType end

function apply!(::GroupAction{LeftMultiplicationGroupAction}, q, g, p)
    return Base.mightalias(q, p) ? q .= g * p : mul!(q, g, p)
end
#=
# I (kellertuer)  am not sure about these – they should be tested and checked once we
# have a use case for them.

function diff_apply!(::GroupAction{LeftMultiplicationGroupAction}, Y, g, p, X)
    return Base.mightalias(Y, X) ? Y .= g * X : mul!(Y, g, X)
end

function diff_group_apply!(::GroupAction{LeftMultiplicationGroupAction}, Y, g, p, X)
    return Base.mightalias(Y, X) ? Y .= X * p : mul!(Y, X, p)
end
=#