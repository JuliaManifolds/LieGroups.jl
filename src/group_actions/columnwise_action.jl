"""
    ColumnwiseAction{A<:AbstractGroupActionType} <: AbstractGroupActionType

A type for an action that applies a group action column-wise on a manifold, where a column-wise
interpretation makes sense, e.g. a matrix-manifold
"""
struct ColumnwiseAction{A <: AbstractGroupActionType} <: AbstractGroupActionType
    action::A
end
ColumnwiseAction(action::A) where {A <: AbstractGroupActionType} = ColumnwiseAction{A}(action)
