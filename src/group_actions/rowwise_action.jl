"""
    RowwiseAction{A<:AbstractGroupActionType} <: AbstractGroupActionType

A type for an action that applies a group action row-wise on a manifold, where a row-wise
interpretation makes sense, e.g. a matrix-manifold.

# Fields
* `action::A`: The group action to be applied row-wise.

# Constructor
    RowwiseAction(action::AbstractGroupActionType)
"""
struct RowwiseAction{A <: AbstractGroupActionType} <: AbstractGroupActionType
    action::A
end
