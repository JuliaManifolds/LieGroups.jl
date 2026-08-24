#
#
# Together with a product group operation, that is a tuple of operations,
# One for each factor, we define the product operation as acting element wise.

"""
    AbstractProductGroupOperation <: AbstractGroupOperation

An abstract type to model group operations on a product manifold
"""
abstract type AbstractProductGroupOperation <: AbstractGroupOperation end

"""
    ProductGroupOperation{O<:NTuple{N,AbstractGroupOperation} where N} <: AbstractProductGroupOperation

A struct do model a tuple of group operations, one for each factor of a product group,
that together forms a new group operation.

Access to the single operations can be done by `pgo[i]`.

# Constructor

    ProductGroupOperation(o::AbstractGroupOperation...)
    ×(o::AbstractGroupOperation...) = ProductGroupOperation(o...)
"""
struct ProductGroupOperation{OTM <: NTuple{N, AbstractGroupOperation} where {N}} <:
    AbstractProductGroupOperation
    operations::OTM
end
function ProductGroupOperation(operations::AbstractGroupOperation...)
    return ProductGroupOperation(operations)
end

@doc raw"""
    cross(O1::AbstractGroupOperation, O2::AbstractGroupOperation)
    O1 × O2
    O1 × O2 × O3 × ...

Return the [`ProductGroupOperation`](@ref) For two [AbstractGroupOperation`](@ref) `O1` and `O2`,
where for the case that one of them is a [`ProductGroupOperation`](@ref) itself,
the other is either prepended (if `O1` is a product) or appended (if `O2` is).
If both are product operations, they are combined into one, keeping the order of operations.

For the case that more than two are concatenated with `×` this is iterated.
"""
LinearAlgebra.cross(::AbstractGroupOperation...)
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
    ProductLieGroup(G1, G2, ..., Gn)

Return the [`LieGroup`](@ref) of the product of Lie groups `G1`, `G2`, up to `Gn`
and all following Lie groups.
This can be considered as a vector of Lie group, where the vector is always
of the same length as the number of provided Lie Groups.

If none of the Lie groups are product Lie groups themselves,
this is equivalent to `G1 × G2 × ... × Gn`.

For an example illustrating the differences see [`x`](@ref LinearAlgebra.cross(::LieGroup...)`(::LieGroup...)`.
"""
function ProductLieGroup(lie_groups::LieGroup...)
    return LieGroup(
        ProductManifold([L.manifold for L in lie_groups]...),
        ProductGroupOperation([L.op for L in lie_groups]...)
    )
end
# Do not “wrap twice”
function ProductLieGroup(G::LieGroup{𝔽, Op, M}) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    return G
end


function ManifoldsBase.submanifold_components(
        ::LieGroup{𝔽, Op, M}, op::ProductGroupOperation
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    return op.operations
end
function ManifoldsBase.submanifold_components(
        PrG::LieGroup{𝔽, Op, M}, ::Identity{Op}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    return map(Identity, PrG.op.operations)
end

function _compose!(
        PrG::LieGroup{𝔽, Op, M}, k, g, h
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    map(
        compose!,
        map(LieGroup, PrG.manifold.manifolds, PrG.op.operations),
        submanifold_components(PrG.manifold, k),
        submanifold_components(PrG.manifold, g),
        submanifold_components(PrG.manifold, h),
    )
    return k
end

function ManifoldsBase.check_size(
        PrG::LieGroup{𝔽, Op, M}, g
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, g)
end
function ManifoldsBase.check_size(
        ::LieGroup{𝔽, Op, M}, ::Identity
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    return nothing
end
function ManifoldsBase.check_size(
        PrG::LieGroup{𝔽, Op, M}, g, X
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, g, X)
end
function ManifoldsBase.check_size(
        PrG::LieGroup{𝔽, Op, M}, ::Identity, X
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    return ManifoldsBase.check_size(PrG.manifold, identity_element(PrG, typeof(X)), X)
end

function conjugate!(
        PrG::LieGroup{𝔽, Op, M}, k, g, h
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        conjugate!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, k),
        submanifold_components(PrM, g),
        submanifold_components(PrM, h),
    )
    return k
end

@doc raw"""
    cross(G::LieGroup, H::LieGroup)
    G × H
    G1 × G2 × G3 × ...

Return the [`ProductLieGroup`](@ref) For two [`LieGroups`](@ref) `G` and `H`.

For the case that one of them is a product Lie group already,
the other one is appended or prepended, depending on which one is the product;
they are joined into one large product if both are product Lie groups.

In order to build a [`ProductLieGroup`](@ref) that does not “splat” its arguments,
or in other words to obtain “nested” products,
use [`ProductLieGroup`](@ref)`(G1, G2, G3, ...)`.

# Example.

For
```
G1 = TranslationGroup(2)
G2 = SpecialOrthogonalGroup(2)
G3 = GeneralLinearGroup(2)
```

We can have one large product Lie group

```
G = G1 × G2 × G3 # or equivalently ProductLieGroup(G1, G2, G3)
```

and alternatively generate a product of a Lie group with an existing product using

```
H = ProductLieGroup(G1, G2 × G3)
```

!!! note "Technical detail"
    Since for the first, single Lie group, the order should be irrelevant, it means
    in practice that `×` behaves slightly different than `ProductLieGroup` in that it “splats” its arguments.
    `G` is equivalent to calling
    `ProductLieGroup(G1, G2) × G3` or ` G1 × ProductLieGroup(G2, G3)`.
    Both, as `G` would consist of vectors of length 3.
    These are different from both
    `ProductLieGroup(ProductLieGroup(G1, G2), G3)` and `ProductLieGroup(G1, ProductLieGroup(G2, G3))`,
    which are both vectors of length 2, where the first has a vector of length 2 in its first component,
    the second such a vector in its second component.
"""
LinearAlgebra.cross(::LieGroup...)
function LinearAlgebra.cross(G::LieGroup, H::LieGroup)
    return ProductLieGroup(G, H)
end
function LinearAlgebra.cross(
        G::LieGroup{𝔽, Op, M}, H::LieGroup
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    return ProductLieGroup(map(LieGroup, G.manifold.manifolds, G.op.operations)..., H)
end
function LinearAlgebra.cross(
        G::LieGroup, H::LieGroup{𝔽, Op, M}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    return ProductLieGroup(G, map(LieGroup, H.manifold.manifolds, H.op.operations)...)
end
function LinearAlgebra.cross(
        G::LieGroup{𝔽1, Op1, M1}, H::LieGroup{𝔽2, Op2, M2}
    ) where {𝔽1, Op1 <: ProductGroupOperation, M1 <: ProductManifold, 𝔽2, Op2 <: ProductGroupOperation, M2 <: ProductManifold}
    return ProductLieGroup(
        map(LieGroup, G.manifold.manifolds, G.op.operations)...,
        map(LieGroup, H.manifold.manifolds, H.op.operations)...
    )
end

function diff_conjugate!(
        PrG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        diff_conjugate!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_inv!(
        PrG::LieGroup{𝔽, Op, M}, Y, g, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        diff_inv!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_left_compose!(
        PrG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        diff_left_compose!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function diff_right_compose!(
        PrG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        diff_right_compose!,
        LieGroup.(PrM.manifolds, PrG.op.operations),
        submanifold_components(PrG, Y),
        submanifold_components(PrG, g),
        submanifold_components(PrG, h),
        submanifold_components(PrG, X),
    )
    return Y
end

function ManifoldsBase.exp!(
        PrG::LieGroup{𝔽, Op, M}, h, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        (M, h, X) -> exp!(M, h, X), # introduce a function with “hard coded” t
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, X),
    )
    return h
end

function ManifoldsBase.exp!(
        PrG::LieGroup{𝔽, Op, M}, h, g, X
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        (M, h, g, X) -> exp!(M, h, g, X), # introduce a function with “hard coded” t
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, g),
        submanifold_components(PrM, X),
    )
    return h
end

_doc_jacobian_exp_product = raw"""
    jacobian_exp(G::ProductLieGroup, X, ::DefaultLieAlgebraOrthogonalBasis)
    jacobian_exp!(G::ProductLieGroup, J, X, ::DefaultLieAlgebraOrthogonalBasis)

Compute the Jacobian of the Lie group exponential in a basis of the Lie algebra on a
[`ProductLieGroup`](@ref).

Since on a (direct) product ``G = G_1 × ⋯ × G_n`` the exponential map acts componentwise
and the basis of the [`LieAlgebra`](@ref) is the concatenation of the bases of the factors,
the Jacobian is block-diagonal with the [`jacobian_exp`](@ref) of the factors on its
diagonal,

````math
J = \begin{pmatrix} J_1 & & \\ & \ddots & \\ & & J_n \end{pmatrix},
\qquad J_i = \text{jacobian\_exp}(G_i, X_i).
````
"""

@doc "$(_doc_jacobian_exp_product)"
jacobian_exp(
    ::LieGroup{𝔽, <:ProductGroupOperation, <:ProductManifold}, X, ::AbstractBasis
) where {𝔽}

@doc "$(_doc_jacobian_exp_product)"
function jacobian_exp!(
        PrG::LieGroup{𝔽, Op, M}, J, X, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    fill!(J, 0)
    foreach(
        (Gi, Xi, dr) -> jacobian_exp!(Gi, view(J, dr, dr), Xi, B),
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        dim_ranges,
    )
    return J
end

function get_vector_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, c, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, submanifold_components(PrM, X), Prc)
    map(ts) do t
        return get_vector_lie!(t..., B)
    end
    return X
end

@inline Base.getindex(pgo::ProductGroupOperation, i::Integer) = pgo.operations[i]
@inline Base.getindex(pgo::ProductGroupOperation, ::Colon) = pgo.operations

function identity_element!(
        PrG::LieGroup{𝔽, Op, M}, e
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        identity_element!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, e),
    )
    return e
end

function ManifoldsBase.inner(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, Y
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold # The product Lie group
    PrM = PrG.manifold # The product manifold
    return sum(
        map(
            inner,
            LieAlgebra.(map(LieGroup, PrM.manifolds, PrG.op.operations)),
            submanifold_components(PrM, X),
            submanifold_components(PrM, Y),
        ),
    )
end

function _inv!(
        PrG::LieGroup{𝔽, Op, M}, h, g
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        submanifold_components(PrM, g),
    )
    return h
end
function Manifolds.inv!(
        PrG::LieGroup{𝔽, Op, M}, h, ::Identity{Op}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, h),
        map(Identity, PrG.op.operations),
    )
    return h
end

function lie_bracket!(
        PrA::LieAlgebra{𝔽, Op, <:LieGroup{𝔽, Op, M}}, Z, X, Y
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrA.manifold.manifold
    map(
        lie_bracket!,
        map(LieAlgebra, LieGroup.(PrM.manifolds, PrA.manifold.op.operations)),
        submanifold_components(PrM, Z),
        submanifold_components(PrM, X),
        submanifold_components(PrM, Y),
    )
    return Z
end

function ManifoldsBase.log!(
        PrG::LieGroup{𝔽, Op, M}, X, ::Identity{Op}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    zero_vector!(LieAlgebra(PrG), X)
    return X
end
function ManifoldsBase.log!(
        PrG::LieGroup{𝔽, Op, M}, X, ::Identity{Op}, ::Identity{Op}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    zero_vector!(LieAlgebra(PrG), X)
    return X
end
function ManifoldsBase.log!(
        PrG::LieGroup{𝔽, Op, M}, X, g, h
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        log!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        submanifold_components(PrM, g),
        submanifold_components(PrM, h),
    )
    return X
end

function ManifoldsBase.log!(
        PrG::LieGroup{𝔽, Op, M}, X, h
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = PrG.manifold
    map(
        log!,
        map(LieGroup, PrM.manifolds, PrG.op.operations),
        submanifold_components(PrM, X),
        submanifold_components(PrM, h),
    )
    return X
end

function Base.show(
        io::IO, G::LieGroup{𝔽, Op, M}
    ) where {𝔽, Op <: ProductGroupOperation, M <: ProductManifold}
    PrM = G.manifold.manifolds
    ops = G.op.operations
    return print(io, "ProductLieGroup($(join(PrM, " × ")), $(join(ops, " × ")))")
end

function get_coordinates_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, c, X, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: AbstractProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc, submanifold_components(PrM, X))
    map(ts) do t
        return get_coordinates_lie!(t..., B)
    end
    return c
end
