"""
   LinearRepresentation{T<:AbstractGroupActionType} <: AbstractGroupActionType

Indicate that a certain group action type has a linear representation.

# Fields
* ` type` – the [`AbstractGroupActionType`](@ref) that has a linear representation
"""
struct LinearRepresentation{T<:AbstractGroupActionType} <: AbstractGroupActionType
    type::T
end

function diff_apply!(
    A::GroupAction{<:LinearRepresentation{<:AbstractGroupActionType},L,M}, Y, g, p, X
)
    return apply!(A, Y, g, X)
end
