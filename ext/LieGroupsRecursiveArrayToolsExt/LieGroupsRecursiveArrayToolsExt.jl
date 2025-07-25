module LieGroupsRecursiveArrayToolsExt

using LieGroups
using RecursiveArrayTools: ArrayPartition
using LinearAlgebra
using ManifoldsBase
using ManifoldsBase: base_manifold, submanifold_components

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

function LieGroups.get_vector_lie(
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}},
    c,
    B::DefaultLieAlgebraOrthogonalBasis,
    ::Type{ArrayPartition},
) where {𝔽,Op<:LieGroups.AbstractProductGroupOperation,M<:ProductManifold}
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
    Pr𝔤::LieAlgebra{𝔽,Op,LieGroup{𝔽,Op,M}},
    c,
    B::DefaultLieAlgebraOrthogonalBasis,
    ::Type{<:ArrayPartition{T,U}},
) where {𝔽,Op<:LieGroups.AbstractProductGroupOperation,M<:ProductManifold,T,U<:Tuple}
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

#
#
# Conversions on SE(n)

"""
TODO: Document and check whether we can have nicer accessors for R and t.
"""
function Base.convert(
    ::Type{<:SpecialEuclideanMatrixPoint}, g::SpecialEuclideanProductPoint
)
    return SpecialEuclideanMatrixPoint(convert(AbstractMatrix, g))
end
function Base.convert(
    ::Type{LinearAlgebra.AbstractMatrix},
    # g has matrix first, vector second, so it is a left semidirect product Lie group point
    g::SpecialEuclideanProductPoint{
        <:ArrayPartition{T,Tuple{<:AbstractVector{T},<:AbstractMatrix{T}}}
    },
) where {T}
    n = length(g.value.x[2])
    A = zeros(T, n + 1, n + 1)
    A[1:n, end] = g.value.x[1] # translation part
    A[1:n, 1:n] = g.value.x[2] # rotation part
    A[end, end] = 1.0 # last entry is always 1
    return A
end
function Base.convert(
    ::Type{LinearAlgebra.AbstractMatrix},
    # g has matrix first, vector second, so it is a left semidirect product Lie group point
    g::SpecialEuclideanProductPoint{
        <:ArrayPartition{T,Tuple{<:AbstractMatrix{T},<:AbstractVector{T}}}
    },
) where {T}
    n = length(g.value.x[2])
    A = zeros(T, n + 1, n + 1)
    A[1:n, end] = g.value.x[2] # translation part
    A[1:n, 1:n] = g.value.x[1] # rotation part
    A[end, end] = 1.0 # last entry is always 1
    return SpecialEuclideanMatrixPoint(A)
end

"""
TODO: Document and check whether we can have nicer accessors for R and t.
"""
function Base.convert(
    ::Type{<:SpecialEuclideanMatrixTangentVector{A}},
    g::SpecialEuclideanProductTangentVector,
) where {A}
    return SpecialEuclideanMatrixTangentVector(convert(AbstractMatrix, g))
end
function Base.convert(
    ::Type{LinearAlgebra.AbstractMatrix},
    # g has matrix first, vector second, so it is a left semidirect product Lie group point
    g::SpecialEuclideanProductTangentVector{
        <:ArrayPartition{T,Tuple{<:AbstractVector{T},<:AbstractMatrix{T}}}
    },
) where {T}
    n = length(g.value.x[2])
    A = zeros(T, n + 1, n + 1)
    A[1:n, end] = g.value.x[1] # translation part
    A[1:n, 1:n] = g.value.x[2] # rotation part
    A[end, end] = 0.0 # last entry is always 0
    return A
end
function Base.convert(
    ::Type{LinearAlgebra.AbstractMatrix},
    # g has matrix first, vector second, so it is a left semidirect product Lie group point
    g::SpecialEuclideanProductTangentVector{
        <:ArrayPartition{T,Tuple{<:AbstractMatrix{T},<:AbstractVector{T}}}
    },
) where {T}
    n = length(g.value.x[2])
    A = zeros(T, n + 1, n + 1)
    A[1:n, end] = g.value.x[2] # translation part
    A[1:n, 1:n] = g.value.x[1] # rotation part
    A[end, end] = 0.0 # last entry is always 0
    return SpecialEuclideanMatrixPoint(A)
end
end
