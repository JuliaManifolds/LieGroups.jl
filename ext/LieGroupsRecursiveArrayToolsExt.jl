module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
# Implement SE(n) also on an Array Partition

function identity_element(G::SpecialEuclideanGroup, T::Type{ArrayPartition})
    # Allocate for the inner manifold (back to default)
    e = ManifoldsBase.allocate_result(G, identity_element)
    return identity_element!(G, e)
end
# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; v=1) = nothing

function ManifoldsBase.submanifold_component(
    G::LieGroups.SpecialEuclideanGroup, g::ArrayPartition, ::Val{I}
) where {I}
    # pass down to manifold by default
    return ManifoldsBase.submanifold_component(G.manifold, g, I)
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 1)
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Translation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Translation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 1)
end

function ManifoldsBase.zero_vector(
    G::LieGroups.LeftSpecialEuclideanGroup,
    e::Identity{LieGroups.SpecialEuclideanOperation},
    ::Union{ComponentsLieAlgebraTVector,ArrayPartition},
)
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(n, n), zeros(n))
end
function ManifoldsBase.zero_vector(
    G::LieGroups.RightSpecialEuclideanGroup,
    e::Identity{LieGroups.SpecialEuclideanOperation},
    ::Union{ComponentsLieAlgebraTVector,ArrayPartition},
)
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(n), zeros(n, n))
end

# TODO: Implement?
# The following three should also work due to the
# AbstractProductGroup s implementations
# * `check_point`
# * `check_size`
# * `check_vector`
end
