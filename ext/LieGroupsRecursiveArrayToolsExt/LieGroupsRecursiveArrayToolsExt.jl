module LieGroupsRecursiveArrayToolsExt

using LieGroups
# As long as Manifolds.jl exports these, we have to specify them here specifically
using LieGroups: identity_element, identity_element!
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using Manifolds
using ManifoldsBase
using ManifoldsBase: base_manifold, submanifold_components, submanifold_component

include("special_euclidean_group_RAT_ext.jl")

function LieGroups.identity_element(
        G::LieGroup{𝔽, <:LieGroups.AbstractProductGroupOperation}, ::Type{ArrayPartition}
    ) where {𝔽}
    Gs = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return ArrayPartition(map(identity_element, Gs)...)
end
function LieGroups.identity_element(
        G::LieGroup{𝔽, <:LieGroups.AbstractProductGroupOperation}, ::Type{<:ArrayPartition{T, U}}
    ) where {𝔽, T, U <: Tuple}
    Gs = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return ArrayPartition(map(identity_element, Gs, U.parameters)...)
end

function LieGroups.get_vector_lie(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}},
        c,
        B::DefaultLieAlgebraOrthogonalBasis,
        ::Type{ArrayPartition},
    ) where {𝔽, Op <: LieGroups.AbstractProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = map(LieAlgebra, map(LieGroup, PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc)
    parts = map(ts) do t
        return LieGroups.get_vector_lie(t..., B)
    end
    return ArrayPartition(parts...)
end

function LieGroups.get_vector_lie(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}},
        c,
        B::DefaultLieAlgebraOrthogonalBasis,
        ::Type{<:ArrayPartition{T, U}},
    ) where {𝔽, Op <: LieGroups.AbstractProductGroupOperation, M <: ProductManifold, T, U <: Tuple}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = map(LieAlgebra, map(LieGroup, PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc, tuple(U.parameters...))
    parts = map(ts) do t
        return LieGroups.get_vector_lie(t[1], t[2], B, t[3])
    end
    return ArrayPartition(parts...)
end
end
