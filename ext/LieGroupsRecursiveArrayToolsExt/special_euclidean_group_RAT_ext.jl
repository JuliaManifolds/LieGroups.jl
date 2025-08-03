# disable affine check
LieGroups._check_matrix_affine(::ArrayPartition, ::Int; kwargs...) = nothing

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
        ::Type{AbstractMatrix}, g::SpecialEuclideanProductPoint{<:ArrayPartition{T}}
    ) where {T}
    v = g.value.x
    n = size(v[1])[1]
    A = zeros(T, n + 1, n + 1)
    A[(n + 1), (n + 1)] = one(T)
    # We do not know whether g.value ir (R,t) or (t,R) so we have to depend on sizes
    s = size(v[1])
    if length(s) == 2 # (R,t)
        A[1:n, 1:n] .= v[1]
        A[1:n, n + 1] .= v[2]
    else # format (t,R)
        A[1:n, 1:n] .= v[2]
        A[1:n, n + 1] .= v[1]
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
        ::Type{AbstractMatrix}, X::SpecialEuclideanProductTangentVector{<:ArrayPartition{T}}
    ) where {T}
    n = size(X.value.x[1])[1]
    A = zeros(T, n + 1, n + 1)
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
function Base.convert(::Type{SpecialEuclideanMatrixTangentVector}, g::ArrayPartition)
    return SpecialEuclideanMatrixTangentVector(
        convert(AbstractMatrix, SpecialEuclideanProductTangentVector(g))
    )
end
function Base.convert(::Type{SpecialEuclideanMatrixPoint}, g::ArrayPartition)
    return SpecialEuclideanMatrixPoint(
        convert(AbstractMatrix, SpecialEuclideanProductPoint(g))
    )
end
# convert between both representation explicitly
# the inverse is again not unique since both (R,t) and (t,R) are possible results.
function Base.convert(
        ::Type{<:SpecialEuclideanMatrixPoint}, g::SpecialEuclideanProductPoint{<:ArrayPartition}
    )
    return SpecialEuclideanMatrixPoint(convert(AbstractMatrix, g))
end
function Base.convert(
        ::Type{SpecialEuclideanProductPoint}, g::SpecialEuclideanMatrixPoint{<:AbstractMatrix}
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
        G::SpecialEuclideanGroup, ::Type{<:SpecialEuclideanProductPoint{A}}
    ) where {A <: ArrayPartition}
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

function LinearAlgebra.norm(
        𝔤::LieAlgebra{
            ℝ, <:LieGroups.SpecialEuclideanGroupOperation, <:LieGroups.SpecialEuclideanGroup,
        },
        X::ArrayPartition,
    )
    G = LieGroups.base_lie_group(𝔤)
    SOn, Tn = LieGroups._SOn_and_Tn(G)
    n1 = LinearAlgebra.norm(
        LieGroups.LieAlgebra(SOn), ManifoldsBase.submanifold_component(𝔤, X, :Rotation)
    )
    n2 = LinearAlgebra.norm(
        LieGroups.LieAlgebra(Tn), ManifoldsBase.submanifold_component(𝔤, X, :Translation)
    )
    return LinearAlgebra.norm([n1, n2])
end

_doc_lie_bracket_SEn_RAT = """
    lie_bracket(𝔰𝔢::LieAlgebra{ℝ, SpecialEuclideanGroupOperation, SpecialEuclideanGroup}, X::ArrayPartition, Y::ArrayPartition)
    lie_bracket!(𝔰𝔢::LieAlgebra{ℝ, SpecialEuclideanGroupOperation, SpecialEuclideanGroup}, Z::ArrayPartition, X::ArrayPartition, Y::ArrayPartition)

Calculate the Lie bracket between elements `X` and `Y` of the Lie algebra of the [`SpecialEuclideanGroup`](@ref).
For the representation as a matrix and a vector, cf. [`SpecialEuclideanProductTangentVector`](@ref) or a `ArrayPartition`
every Lie algebra element is represented as a pair ``X = (X_{$(LieGroups._tex(:text, "R"))}, X_$(LieGroups._tex(:text, "t")))``
or a rotation matrix and a translation vector, respectively.

Then the formula for the Lie bracket is given by
```math
[X, Y] = [(X_{$(LieGroups._tex(:text, "R"))}, X_{$(LieGroups._tex(:text, "t"))}), (Y_{$(LieGroups._tex(:text, "R"))}, Y_{$(LieGroups._tex(:text, "t"))})]
= (X_{$(LieGroups._tex(:text, "R"))} * Y_{$(LieGroups._tex(:text, "R"))} - Y_{$(LieGroups._tex(:text, "R"))} * X_{$(LieGroups._tex(:text, "R"))}, X_{$(LieGroups._tex(:text, "R"))} * Y_{$(LieGroups._tex(:text, "t"))} - Y_{$(LieGroups._tex(:text, "R"))} * X_{$(LieGroups._tex(:text, "t"))}),
```

where for the right semidirect product variant, the order of the pair is switched.
"""

"$(_doc_lie_bracket_SEn_RAT)"
LieGroups.lie_bracket(
    𝔤::LieGroups.LieAlgebra{
        ℝ, <:LieGroups.SpecialEuclideanGroupOperation, <:LieGroups.SpecialEuclideanGroup,
    },
    X::Union{<:ArrayPartition, <:SpecialEuclideanProductTangentVector},
    Y::Union{<:ArrayPartition, <:SpecialEuclideanProductTangentVector},
)

"$(_doc_lie_bracket_SEn_RAT)"
function LieGroups.lie_bracket!(
        𝔤::LieGroups.LieAlgebra{
            ℝ, <:LieGroups.SpecialEuclideanGroupOperation, <:LieGroups.SpecialEuclideanGroup,
        },
        Z::Union{<:ArrayPartition, <:SpecialEuclideanProductTangentVector},
        X::Union{<:ArrayPartition, <:SpecialEuclideanProductTangentVector},
        Y::Union{<:ArrayPartition, <:SpecialEuclideanProductTangentVector},
    )
    G = LieGroups.base_lie_group(𝔤)
    SOn, _ = LieGroups._SOn_and_Tn(G)
    X_t = submanifold_component(LieAlgebra(G), X, Val(:Translation))
    X_R = submanifold_component(LieAlgebra(G), X, Val(:Rotation))
    Y_t = submanifold_component(LieAlgebra(G), Y, Val(:Translation))
    Y_R = submanifold_component(LieAlgebra(G), Y, Val(:Rotation))
    Z_t = submanifold_component(LieAlgebra(G), Z, Val(:Translation))
    Z_R = submanifold_component(LieAlgebra(G), Z, Val(:Rotation))
    LieGroups.lie_bracket!(LieAlgebra(SOn), Z_R, X_R, Y_R)
    Z_t .= X_R * Y_t .- Y_R * X_t
    return Z
end

function ManifoldsBase.submanifold_component(
        G::LieGroups.LeftSpecialEuclideanGroup,
        g::Union{ArrayPartition, SpecialEuclideanProductPoint},
        ::Val{:Rotation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(G), ManifoldsBase.internal_value(g), Val(1)
    )
end
function ManifoldsBase.submanifold_component(
        G::LieGroups.LeftSpecialEuclideanGroup,
        g::Union{ArrayPartition, SpecialEuclideanProductPoint},
        ::Val{:Translation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(G), ManifoldsBase.internal_value(g), Val(2)
    )
end

function ManifoldsBase.submanifold_component(
        𝔤::LieGroups.LieAlgebra{
            ℝ,
            <:LieGroups.LeftSpecialEuclideanGroupOperation,
            <:LieGroups.LeftSpecialEuclideanGroup,
        },
        X::Union{ArrayPartition, SpecialEuclideanProductTangentVector},
        ::Val{:Rotation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(𝔤), ManifoldsBase.internal_value(X), Val(1)
    )
end
function ManifoldsBase.submanifold_component(
        𝔤::LieGroups.LieAlgebra{
            ℝ,
            <:LieGroups.LeftSpecialEuclideanGroupOperation,
            <:LieGroups.LeftSpecialEuclideanGroup,
        },
        X::Union{ArrayPartition, SpecialEuclideanProductTangentVector},
        ::Val{:Translation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(𝔤), ManifoldsBase.internal_value(X), Val(2)
    )
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::LieGroups.RightSpecialEuclideanGroup,
        g::Union{ArrayPartition, SpecialEuclideanProductPoint},
        ::Val{:Rotation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(G), ManifoldsBase.internal_value(g), Val(2)
    )
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::LieGroups.RightSpecialEuclideanGroup,
        g::Union{ArrayPartition, SpecialEuclideanProductPoint},
        ::Val{:Translation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(G), ManifoldsBase.internal_value(g), Val(1)
    )
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        𝔤::LieGroups.LieAlgebra{
            ℝ,
            <:LieGroups.RightSpecialEuclideanGroupOperation,
            <:LieGroups.RightSpecialEuclideanGroup,
        },
        X::Union{ArrayPartition, SpecialEuclideanProductTangentVector},
        ::Val{:Rotation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(𝔤), ManifoldsBase.internal_value(X), Val(2)
    )
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        𝔤::LieGroups.LieAlgebra{
            ℝ,
            <:LieGroups.RightSpecialEuclideanGroupOperation,
            <:LieGroups.RightSpecialEuclideanGroup,
        },
        X::Union{ArrayPartition, SpecialEuclideanProductTangentVector},
        ::Val{:Translation},
    )
    return ManifoldsBase.submanifold_component(
        base_manifold(𝔤), ManifoldsBase.internal_value(X), Val(1)
    )
end

function ManifoldsBase.zero_vector(
        𝔤::LieAlgebra{
            ManifoldsBase.ℝ,
            <:LieGroups.LeftSpecialEuclideanGroupOperation,
            <:LieGroups.LeftSpecialEuclideanGroup,
        },
        ::Type{<:ArrayPartition{T}},
    ) where {T}
    G = 𝔤.manifold
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(T, n, n), zeros(T, n))
end
function ManifoldsBase.zero_vector(
        𝔤::LieAlgebra{
            𝔽,
            <:LieGroups.RightSpecialEuclideanGroupOperation,
            <:LieGroups.RightSpecialEuclideanGroup,
        },
        ::Type{<:ArrayPartition{T}},
    ) where {𝔽, T}
    G = 𝔤.manifold
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return ArrayPartition(zeros(T, n), zeros(T, n, n))
end
function ManifoldsBase.zero_vector(
        𝔤::LieAlgebra{
            𝔽, <:LieGroups.SpecialEuclideanGroupOperation, <:LieGroups.SpecialEuclideanGroup,
        },
        ::Type{LieGroups.SpecialEuclideanProductTangentVector{AP}},
    ) where {𝔽, AP <: ArrayPartition}
    return SpecialEuclideanProductTangentVector(zero_vector(𝔤, AP))
end
