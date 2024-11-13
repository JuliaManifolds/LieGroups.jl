#
#
# Together with a product group operation, that is a tuple of operations,
# One for each factor, we define the product operation as acting element wise.

"""
    ProductGroupOperation{O} <: AbstractGroupOperation

A struct do model a tuple of group operations, one for each factor of a product group,
that together forms a new group operation.

# Constructor

    ProductGroupOperation(O...)
    ×(O...) = ProductGroupOperation(O...)
"""
struct ProductGroupOperation{OTM<:Tuple} <: AbstractGroupOperation
    operations::OTM
end
function ProductGroupOperation(operations::AbstractGroupOperation...)
    return ProductGroupOperation(operations)
end

@doc raw"""
    cross(O1, O2)
    O1 × O2
    O1 × O2 × O3 × ...

Return the [`ProductGroupOperation`](@ref) For two [AbstractGroupOperation`](@ref) `O1` and `O2`,
where for the case that one of them is a [`ProductGroupOperation`](@ref) itself,
the other is either prepended (if `O1` is a product) or appenden (if `O2` is).
If both are product operations, they are combined into one, keeping the order of operations.

For the case that more than two are concatenated with `×` this is iterated.
"""
cross(::AbstractGroupOperation...)
function LinearAlgebra.cross(O1::AbstractGroupOperation, O2::AbstractGroupOperation)
    return ProductGroupOperation(O1, O2)
end
function LinearAlgebra.cross(P::ProductGroupOperation, O::AbstractGroupOperation)
    return ProductGroupOperation(P.operations..., O)
end
function LinearAlgebra.cross(O::AbstractGroupOperation, P::ProductGroupOperation)
    return ProductGroupOperation(O, P.operations...)
end
function LinearAlgebra.cross(P1::ProductGroupOperation, P2::ProductGroupOperation)
    return ProductGroupOperation(P1.operations..., P2.operations...)
end

"""
    ProductLieGroup(G, H, ...)

Return the [`LieGroup`](@ref) of the product of Lie groups `G` and `H`.

Alternatively, the short hand `G × H` can be used.
"""
function ProductLieGroup(G::LieGroup, H::LieGroup)
    return LieGroup(G.manifold × H.manifold, G.op × H.op)
end

@doc raw"""
    cross(G, H)
    G × H
    G1 × G2 × G3 × ...

Return the [`ProductLieGroup`](@ref) For two [`LieGroups`](@ref) `G` and `H`,
where for the case that one of them is a [`ProductLieGroup`](@ref) itself,
the other is either prepended (if `H` is a product) or appenden (if `G` is).
If both are product Lie groups, they are combined into one, keeping the order of operations.

For the case that more than two are concatenated with `×` this is iterated.
"""
cross(::LieGroup...)
function LinearAlgebra.cross(G::LieGroup, H::LieGroup)
    return ProductLieGroup(G, H)
end

function Base.show(
    io::IO, G::LieGroup{𝔽,<:ProductGroupOperation,<:ManifoldsBase.ProductManifold}
) where {𝔽}
    M = G.manifold.manifolds
    ops = G.op.operations
    return print(io, "ProductLieGroup($(join(M, " × ")), $(join(ops, " × ")))")
end
