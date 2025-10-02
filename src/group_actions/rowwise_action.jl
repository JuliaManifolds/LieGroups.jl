"""
    RowwiseGroupAction{A<:AbstractGroupActionType} <: AbstractGroupActionType

A type for an action that applies a group action row-wise on a manifold, where a row-wise
interpretation makes sense, e.g. a matrix-manifold.

# Fields
* `action::A`: The group action to be applied row-wise.

# Constructor
    RowwiseGroupAction(action::AbstractGroupActionType)
"""
struct RowwiseGroupAction{A <: AbstractGroupActionType} <: AbstractGroupActionType
    action::A
end

function apply(::GroupAction{RowwiseGroupAction{LeftMultiplicationGroupAction}}, g, p)
    return (g * p')'
end
function apply!(a::GroupAction{RowwiseGroupAction{LeftMultiplicationGroupAction}}, q, g, p)
    b = GroupAction(a.type.action, a.group, a.manifold)
    map((qrow, prow) -> apply!(b, qrow, g, prow), eachrow(q), eachrow(p))
    return q
end
