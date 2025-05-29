module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: base_manifold

include("special_euclidean_group_RAT_ext.jl")

function LieGroups.identity_element(
    G::LieGroup{𝔽,<:LieGroups.AbstractProductGroupOperation}, ::Type{ArrayPartition}
) where {𝔽}
    Gs = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return ArrayPartition(map(identity_element, Gs)...)
end
function LieGroups.identity_element(
    G::LieGroup{𝔽,<:LieGroups.AbstractProductGroupOperation}, ::Type{<:ArrayPartition{T,U}}
) where {𝔽,T,U<:Tuple}
    Gs = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return ArrayPartition(map(identity_element, Gs, U.parameters)...)
end

end
