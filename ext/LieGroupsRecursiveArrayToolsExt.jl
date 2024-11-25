module LieGroupsRecursiveArrayToolsExt

using LieGroups, ManifoldsBase, RecursiveArrayTools

using LieGroups: identity_element, LieGroup

function ManifoldsBase.allocate_result(
    G::LieGroup{𝔽,Op,M}, ::typeof(identity_element)
) where {𝔽,Op,<:ManifoldsBase.ProductManifold}
    M = base_manifold(G)
    Ls = LieGroup.(M.manifolds, G.op.operations)
    ps = ManifoldsBase.allocate_result.(Ls, Ref(identity_element))
    return ArrayPartition(ps...)
end
end
