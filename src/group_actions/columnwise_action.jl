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

function apply(::GroupAction{ColumnwiseGroupAction{LeftMultiplicationGroupAction}}, g, p)
    return g * p
end
function apply!(a::GroupAction{ColumnwiseGroupAction{LeftMultiplicationGroupAction}}, q, g, p)
    b = GroupAction(a.type.action, a.group, a.manifold)
    map((qcol, pcol) -> apply!(b, qcol, g, pcol), eachcol(q), eachcol(p))
    return q
end

function diff_group_apply!(
        ::GroupAction{ColumnwiseGroupAction{LeftMultiplicationGroupAction}},
        Y,
        ::Identity{MatrixMultiplicationGroupOperation},
        p,
        X,
    )
    return Base.mightalias(Y, X) ? Y .= X * p : mul!(Y, X, p)
end
