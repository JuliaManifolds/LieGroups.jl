
# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; v=1) = nothing

#
# Conversions
#
Base.convert(::Type{<:ArrayPartition}, g::SpecialEuclideanProductPoint) = g.value
function Base.convert(::Type{SpecialEuclideanProductPoint}, p::ArrayPartition)
    return SpecialEuclideanProductPoint(p)
end
Base.convert(::Type{<:ArrayPartition}, X::SpecialEuclideanProductTangentVector) = X.value
function Base.convert(::Type{SpecialEuclideanProductTangentVector}, X::ArrayPartition)
    return SpecialEuclideanProductTangentVector(X)
end
# convert between both representations
function Base.convert(
    ::Type{<:AbstractMatrix}, g::SpecialEuclideanProductPoint{<:ArrayPartition}
)
    n = size(g.value.x[1])[1]
    A = zeros((n + 1), (n + 1))
    A[(n + 1), (n + 1)] = 1.0
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(g.value.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= g.value.x[1]
        A[1:n, (n + 1)] .= g.value.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= g.value.x[2]
        A[1:n, (n + 1)] .= g.value.x[1]
    end
    return A
end
function Base.convert(::Type{<:SpecialEuclideanProductPoint}, g::AbstractMatrix)
    return SpecialEuclideanProductPoint(
        convert(ArrayPartition, SpecialEuclideanMatrixPoint(g))
    )
end
function Base.convert(::Type{<:SpecialEuclideanProductTangentVector}, X::AbstractMatrix)
    return SpecialEuclideanProductTangentVector(
        convert(ArrayPartition, SpecialEuclideanMatrixTangentVector(X))
    )
end

function Base.convert(
    ::Type{<:ArrayPartition}, g::SpecialEuclideanMatrixPoint{<:AbstractMatrix}
)
    A = g.value
    n = size(A)[1] - 1
    return ArrayPartition(A[1:n, 1:n], A[1:n, (n + 1)])
end
function Base.convert(
    ::Type{<:AbstractMatrix}, X::SpecialEuclideanProductTangentVector{<:ArrayPartition}
)
    n = size(X.value.x[1])[1]
    A = zeros(n + 1, n + 1)
    A[n + 1, n + 1] = 0.0
    # We do not know whether g.value is (R,t) or (t,R) so we have to depend on sizes
    s = size(X.value.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= X.value.x[1]
        A[1:n, n + 1] .= X.value.x[2]
    else # format (t,R)
        A[1:n, 1:n] .= X.value.x[2]
        A[1:n, n + 1] .= X.value.x[1]
    end
    return A
end
# The reverse is hence not unique, but we can still cast to the special matrix point
function Base.convert(::Type{T}, g::ArrayPartition) where {T<:SpecialEuclideanMatrixPoint}
    return T(convert(AbstractMatrix, SpecialEuclideanProductPoint(g)))
end
function Base.convert(
    ::Type{T}, g::ArrayPartition
) where {T<:SpecialEuclideanMatrixTangentVector}
    return T(convert(AbstractMatrix, SpecialEuclideanProductTangentVector(g)))
end
function Base.convert(
    ::M, g::SpecialEuclideanProductPoint{<:ArrayPartition}
) where {T,M<:AbstractMatrix{T}}
    v = g.value.x
    n = size(v[1])[1]
    A = zeros(T, n + 1, n + 1)
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(g.x[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= v[1]
        A[1:n, n + 1] .= v[2]
    else # format (t,R)
        A[1:n, 1:n] .= v[2]
        A[1:n, n + 1] .= v[1]
    end
    return A
end
# convert between both representation explicitly
# the inverse is again not unique since both (R,t) and (t,R) are possible results.
function Base.convert(
    ::Type{<:SpecialEuclideanMatrixPoint}, g::SpecialEuclideanProductPoint{<:ArrayPartition}
)
    return SpecialEuclideanMatrixPoint(convert(AbstractMatrix, g))
end
function Base.convert(
    ::Type{<:SpecialEuclideanProductPoint}, g::SpecialEuclideanMatrixPoint{<:AbstractMatrix}
)
    return SpecialEuclideanProductPoint(convert(ArrayPartition, g))
end
function Base.convert(
    ::Type{SpecialEuclideanMatrixTangentVector},
    g::SpecialEuclideanProductTangentVector{<:ArrayPartition},
)
    return SpecialEuclideanMatrixTangentVector(convert(AbstractMatrix, g))
end
function Base.convert(
    ::Type{<:ArrayPartition}, g::SpecialEuclideanMatrixTangentVector{<:AbstractMatrix}
)
    A = g.value
    n = size(A)[1] - 1
    return ArrayPartition(A[1:n, 1:n], A[1:n, (n + 1)])
end
function Base.convert(
    ::Type{<:SpecialEuclideanProductTangentVector},
    g::SpecialEuclideanMatrixTangentVector{<:AbstractMatrix},
)
    return SpecialEuclideanProductTangentVector(convert(ArrayPartition, g))
end
#
# Functions specialised from the interface
#

function ManifoldsBase.exp!(
    G::SpecialEuclideanGroup{<:ManifoldsBase.TypeParameter{Tuple{2}}},
    g::ArrayPartition,
    X::ArrayPartition,
)
    # This is up to the dispatch types a nearly copy of the matrix case, since
    # but we can skip to initialise the constant areas
    LieGroups._exp_SE2!(G, g, X)
    return g
end

function ManifoldsBase.exp!(
    G::SpecialEuclideanGroup{<:ManifoldsBase.TypeParameter{Tuple{3}}},
    g::ArrayPartition,
    X::ArrayPartition,
)
    # This is up to the dispatch types a nearly copy of the matrix case, since
    # but we can skip to initialise the constant areas
    LieGroups._exp_SE3!(G, g, X)
    return g
end

function LieGroups.identity_element(
    G::LieGroups.LeftSpecialEuclideanGroup, ::Type{<:ArrayPartition{T,Tuple{A,B}}}
) where {T,A,B}
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(SOn, A), identity_element(Tn, B))
end
function LieGroups.identity_element(
    G::LieGroups.RightSpecialEuclideanGroup, ::Type{<:ArrayPartition{T,Tuple{A,B}}}
) where {T,A,B}
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    return ArrayPartition(identity_element(Tn, A), identity_element(SOn, B))
end
function LieGroups.identity_element(
    G::SpecialEuclideanGroup, ::Type{<:SpecialEuclideanProductPoint{A}}
) where {A<:ArrayPartition}
    return SpecialEuclideanProductPoint(identity_element(G, A))
end
function LieGroups.identity_element!(G::LieGroups.SpecialEuclideanGroup, g::ArrayPartition)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    identity_element!(SOn, ManifoldsBase.submanifold_component(G, g, :Rotation))
    identity_element!(Tn, ManifoldsBase.submanifold_component(G, g, :Translation))
    return g
end

function LieGroups.inv!(G::SpecialEuclideanGroup, h::ArrayPartition, g::ArrayPartition)
    LieGroups._inv_SE!(G, h, g)
    return h
end

function ManifoldsBase.log!(
    G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    X::ArrayPartition,
    g::ArrayPartition,
)
    LieGroups._log_SE2!(G, X, g)
    return X
end

function ManifoldsBase.log!(
    G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    X::ArrayPartition,
    g::ArrayPartition,
)
    LieGroups._log_SE3!(G, X, g)
    return X
end

function ManifoldsBase.submanifold_component(
    G::SE,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{I},
) where {I<:Int,SE<:SpecialEuclideanGroup}
    # pass down to manifold by default
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), I
    )
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Rotation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 1
    )
end
function ManifoldsBase.submanifold_component(
    G::LieGroups.LeftSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 2
    )
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Rotation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 2
    )
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::LieGroups.RightSpecialEuclideanGroup,
    g::Union{
        ArrayPartition,SpecialEuclideanProductPoint,SpecialEuclideanProductTangentVector
    },
    ::Val{:Translation},
)
    return ManifoldsBase.submanifold_component(
        G.manifold, ManifoldsBase.internal_value(g), 1
    )
end

function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{
        𝔽,
        <:LieGroups.LeftSpecialEuclideanGroupOperation,
        <:LieGroups.LeftSpecialEuclideanGroup,
    },
    ::ArrayPartition{T},
) where {𝔽,T}
    G = 𝔤.manifold
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(T, zeros(n, n), zeros(n))
end
function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{
        𝔽,
        <:LieGroups.RightSpecialEuclideanGroupOperation,
        <:LieGroups.RightSpecialEuclideanGroup,
    },
    ::ArrayPartition{T},
) where {𝔽,T}
    G = 𝔤.manifold
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(T, zeros(n), zeros(n, n))
end
function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{
        𝔽,<:LieGroups.SpecialEuclideanGroupOperation,<:LieGroups.SpecialEuclideanGroup
    },
    ::Type{LieGroups.SpecialEuclideanProductTangentVector{AP}},
) where {𝔽,AP<:ArrayPartition}
    return SpecialEuclideanProductTangentVector(zero_vector(𝔤, AP))
end
