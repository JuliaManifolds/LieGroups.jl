"""
    ColumnwiseAction{A<:AbstractGroupActionType} <: AbstractGroupActionType

A type for an action that applies a group action column-wise on a manifold, where a column-wise
interpretation makes sense, e.g. a matrix-manifold.

# Fields
* `action::A`: The group action to be applied column-wise.

# Constructor
    ColumnwiseAction(action::AbstractGroupActionType)
"""
struct ColumnwiseAction{A <: AbstractGroupActionType} <: AbstractGroupActionType
    action::A
end
