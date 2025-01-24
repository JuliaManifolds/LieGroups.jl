module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
# Implement SE(n) also on an Array Partition

function _value(v::Union{SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector})
    return v.value
end
_value(v::ArrayPartition) = v

function LieGroups.identity_element(
    G::LieGroups.LeftSpecialEuclideanGroup, ::Type{<:ArrayPartition}
)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(SOn), identity_element(Tn))
end
function LieGroups.identity_element(
    G::LieGroups.RightSpecialEuclideanGroup, ::Type{<:ArrayPartition}
)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(Tn), identity_element(SOn))
end
# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; v=1) = nothing

function ManifoldsBase.submanifold_component(
    G::LieGroups.SpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{I},
) where {I}
    # pass down to manifold by default
    return ManifoldsBase.submanifold_component(G.manifold, _value(g), I)
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Rotation},
)
    return ManifoldsBase.submanifold_component(G.manifold, _value(g), 1)
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(G.manifold, _value(g), 2)
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup, g::ArrayPartition, ::Val{:Rotation}
)
    return ManifoldsBase.submanifold_component(G.manifold, g, 2)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(G.manifold, _value(g), 1)
end

function ManifoldsBase.zero_vector(
    G::LieAlgebra{
        𝔽,LieGroups.SpecialEuclideanGroupOperation,LieGroups.SpecialEuclideanGroup
    },
    ::ArrayPartition,
) where {𝔽}
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(n), zeros(n, n))
end
function ManifoldsBase.zero_vector(
    G::LieAlgebra{
        𝔽,LieGroups.SpecialEuclideanGroupOperation,LieGroups.SpecialEuclideanGroup
    },
    ::SpecialEuclideanProductTangentVector,
) where {𝔽}
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return SpecialEuclideanProductTangentVector(ArrayPartition(zeros(n), zeros(n, n)))
end

# TODO: Implement?
# The following three should also work due to the
# AbstractProductGroup s implementations
# * `check_point`
# * `check_size`
# * `check_vector`
end
