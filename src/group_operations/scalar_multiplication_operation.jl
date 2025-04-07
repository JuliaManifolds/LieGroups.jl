"""
    AbelianMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a scalar multiplication.
"""
struct AbelianMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

function _compose(
    ::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * h
end

function _compose(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    p::AbstractArray{<:Any,0},
    q::AbstractArray{<:Any,0},
) where {𝔽}
    return map((pp, qq) -> compose(G, pp, qq), p, q)
end

function _compose(
    ::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    g::Number,
    h::AbstractArray{<:Any,0},
) where {𝔽}
    return g .* h
end

function _compose(
    ::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::Number,
) where {𝔽}
    return g .* h
end

function _compose!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, k, g, h) where {𝔽}
    return copyto!(k, compose(G, g, h))
end

function conjugate(::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, h) where {𝔽}
    return g
end

function conjugate!(
    ::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, k::AbstractArray{<:Any,0}, g, h
) where {𝔽}
    return copyto!(k, g)
end

function diff_conjugate(
    ::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, ::Any, ::Any, X::Number
) where {𝔽}
    return X
end

function diff_conjugate(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map(XX -> diff_conjugate(G, g, h, XX), X)
end

function diff_conjugate!(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_conjugate(G, g, h, X))
end

function diff_inv!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, -X)
end

function diff_left_compose(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> diff_left_compose(G, gg, h, XX), g, X)
end

function diff_left_compose!(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_left_compose(G, g, h, X))
end

function diff_right_compose(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::Any,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> diff_right_compose(G, gg, h, XX), g, X)
end

function diff_right_compose!(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_right_compose(G, g, h, X))
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, X::AbstractArray{<:Any,0}
) where {𝔽}
    return map(XX -> exp(G, XX), X)
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, g::Number, X::Number
) where {𝔽}
    return g * exp(G, X)
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> gg * exp(G, XX), g, X)
end

function ManifoldsBase.exp!(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, g, X
) where {𝔽}
    return copyto!(g, exp(G, X))
end

function ManifoldsBase.exp!(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, h, g, X
) where {𝔽}
    return copyto!(h, exp(G, g, X))
end

Base.inv(::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g::Number) where {𝔽} = inv(g)
function Base.inv(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g::AbstractArray{<:Any,0}
) where {𝔽}
    return map(gg -> inv(G, gg), g)
end

function inv!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, h, g) where {𝔽}
    copyto!(h, inv(G, g))
    return h
end

function inv!(
    G::LieGroup{𝔽,O}, g, ::Identity{O}
) where {𝔽,O<:AbelianMultiplicationGroupOperation}
    return identity_element!(G, g)
end

function inv_left_compose(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return inv(g) * h
end

function inv_left_compose(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::AbstractArray{<:Any,0},
) where {𝔽}
    return map(\, g, h)
end

function inv_left_compose!(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, inv_left_compose(G, g, h))
end

function inv_right_compose(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * inv(h)
end

function inv_right_compose(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::AbstractArray{<:Any,0},
) where {𝔽}
    return map(/, g, h)
end

function inv_right_compose!(
    G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, inv_right_compose(G, g, h))
end

_doc_identity_element_scalar_mult = """
    identity_element(G::LieGroup{𝔽,AbelianMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`AbelianMultiplicationGroupOperation`](@ref) is the one-element.
"""

@doc "$(_doc_identity_element_scalar_mult)"
identity_element(::LieGroup{𝔽,AbelianMultiplicationGroupOperation}) where {𝔽} = 1.0

function identity_element(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, ::Type{T}
) where {𝔽,T<:Union{Number,AbstractArray{0,<:Number}}}
    return one(T)
end
function identity_element(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, ::Type{Array{T,0}}
) where {𝔽,T<:Number}
    return fill(one(T))
end
function identity_element(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, e::Number
) where {𝔽}
    return one(e)
end

@doc "$(_doc_identity_element_scalar_mult)"
identity_element!(::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, e) where {𝔽}

function identity_element!(
    ::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, e::AbstractArray{<:Number,0}
) where {𝔽}
    return fill!(e, 1)
end

function lie_bracket(
    ::LieAlgebra{𝔽,AbelianMultiplicationGroupOperation}, X::Any, ::Any
) where {𝔽}
    return zero(X)
end

function lie_bracket!(
    ::LieAlgebra{𝔽,AbelianMultiplicationGroupOperation}, Z, X, Y
) where {𝔽}
    return copyto!(Z, zero(Z))
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, g, h
) where {𝔽}
    return log(G, compose(G, inv(G, g), h))
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    ::Identity{AbelianMultiplicationGroupOperation},
) where {𝔽}
    return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    ::Identity{AbelianMultiplicationGroupOperation},
    T::Type,
) where {𝔽}
    return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(
    G::LieGroup{𝔽,AbelianMultiplicationGroupOperation},
    X,
    ::Identity{AbelianMultiplicationGroupOperation},
) where {𝔽}
    return zero_vector!(LieAlgebra(G), X)
end
