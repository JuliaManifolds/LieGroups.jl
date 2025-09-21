"""
    ColumnwiseGroupAction{A<:AbstractGroupActionType} <: AbstractGroupActionType

A type for an action that applies a group action column-wise on a manifold, where a column-wise
interpretation makes sense, e.g. a matrix-manifold.

# Fields
* `action::A`: The group action to be applied column-wise.

# Constructor
    ColumnwiseGroupAction(action::AbstractGroupActionType)
"""
struct ColumnwiseGroupAction{A <: AbstractGroupActionType} <: AbstractGroupActionType
    action::A
end

function apply(a::GroupAction{ColumnwiseGroupAction{<:LeftMultiplicationGroupAction}}, g, p)
    return g * p
end
function apply!(a::GroupAction{ColumnwiseGroupAction}, q, g, p)
    b = GroupAction(a.action, a.group, a.manifold)
    return map((qcol, pcol) -> apply!(b, qcol, g, pcol), eachcol(q), eachcol(p))
end
