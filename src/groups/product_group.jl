#
#
# Together with a product group operation, that is a tuple of operations,
# One for each factor, we define the product operation as acting element wise.

"""
    ProductOperation{O} <: AbstractGroupOperation

A struct do model a tuple of group operations, one for each factor of a product group,
that together forms a new group operation.

# Constructor

    ProductOperation(O...)
"""
struct ProductOperation{OTM<:Tuple} <: AbstractGroupOperation
    operations::OTM
end
ProductOperation(operations::AbstractGroupOperation...) = ProductOperation(operations)
