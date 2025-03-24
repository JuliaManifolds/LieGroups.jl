"""
    ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a scalar multiplication.
"""
struct ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

function compose(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * h
end

function compose(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    p::AbstractArray{<:Any,0},
    q::AbstractArray{<:Any,0},
) where {𝔽}
    return map((pp, qq) -> compose(G, pp, qq), p, q)
end

function compose(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    p::AbstractArray{<:Any,0},
    q::Number,
) where {𝔽}
    return map(pp -> compose(G, pp, q), p)
end

function compose(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    p::Number,
    q::AbstractArray{<:Any,0},
) where {𝔽}
    return map(qq -> compose(G, p, qq), q)
end

function _compose!(G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k, g, h) where {𝔽}
    return copyto!(k, compose(G, g, h))
end

function conjugate(::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g, h) where {𝔽}
    return g
end

function conjugate!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k::AbstractArray{<:Any,0}, g, h
) where {𝔽}
    return copyto!(k, g)
end

function diff_conjugate(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, ::Any, ::Any, X::Number
) where {𝔽}
    return X
end

function diff_conjugate(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map(XX -> diff_conjugate(G, g, h, XX), X)
end

function diff_conjugate!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_conjugate(G, g, h, X))
end

diff_inv(::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g, X) where {𝔽} = -X

function diff_inv!(G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, Y, g, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, -X)
end

function diff_left_compose(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::Any,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> diff_left_compose(G, gg, h, XX), g, X)
end

function diff_left_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_left_compose(G, g, h, X))
end

function diff_right_compose(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::Any,
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> diff_right_compose(G, gg, h, XX), g, X)
end

function diff_right_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_right_compose(G, g, h, X))
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, X::AbstractArray{<:Any,0}
) where {𝔽}
    return map(XX -> exp(G, XX), X)
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::Number, X::Number
) where {𝔽}
    return g * exp(G, X)
end

function ManifoldsBase.exp(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    X::AbstractArray{<:Any,0},
) where {𝔽}
    return map((gg, XX) -> gg * exp(G, XX), g, X)
end

function ManifoldsBase.exp!(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g, X
) where {𝔽}
    return copyto!(g, exp(G, X))
end

function ManifoldsBase.exp!(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, h, g, X
) where {𝔽}
    return copyto!(h, exp(G, g, X))
end

#diff_right_compose(G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g, h, X) where {𝔽} = diff_left_compose(G, g, h, X)

Base.inv(::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::Number) where {𝔽} = inv(g)
function Base.inv(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::AbstractArray{<:Any,0}
) where {𝔽}
    return map(gg -> inv(G, gg), g)
end

function inv!(G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, h, g) where {𝔽}
    copyto!(h, inv(G, g))
    return h
end

function inv!(
    G::LieGroup{𝔽,O}, g, ::Identity{O}
) where {𝔽,O<:ScalarMultiplicationGroupOperation}
    return identity_element!(G, g)
end

function inv_left_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return inv(g) * h
end

function inv_left_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::AbstractArray{<:Any,0},
) where {𝔽}
    return map(\, g, h)
end

function inv_left_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, inv_left_compose(G, g, h))
end

function inv_right_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * inv(h)
end

function inv_right_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::AbstractArray{<:Any,0},
    h::AbstractArray{<:Any,0},
) where {𝔽}
    return map(/, g, h)
end

function inv_right_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, inv_right_compose(G, g, h))
end

_doc_identity_element_scalar_mult = """
    identity_element(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`ScalarMultiplicationGroupOperation`](@ref) is the one-element.
"""

@doc "$(_doc_identity_element_scalar_mult)"
identity_element(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}) where {𝔽} = 1.0

function identity_element(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e::Union{<:Number,Type{<:Number}}
) where {𝔽}
    return one(e)
end

@doc "$(_doc_identity_element_scalar_mult)"
identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}

function identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}
    return fill!(e, 1.0)
end

function lie_bracket(
    ::LieAlgebra{𝔽,ScalarMultiplicationGroupOperation}, X::Any, ::Any
) where {𝔽}
    return zero(X)
end

function lie_bracket!(
    ::LieAlgebra{𝔽,ScalarMultiplicationGroupOperation}, Z, ::Any, ::Any
) where {𝔽}
    return copyto!(Z, zero(Z))
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::AbstractArray{<:Any,0}
) where {𝔽}
    return map(gg -> log(G, gg), g)
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g, h
) where {𝔽}
    return log(G, compose(G, inv(G, g), h))
end

function ManifoldsBase.log!(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, X, g
) where {𝔽}
    return copyto!(X, log(G, g))
end

function ManifoldsBase.log(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    ::Identity{ScalarMultiplicationGroupOperation},
) where {𝔽}
    return zero_vector(LieAlgebra(G))
end
function ManifoldsBase.log(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    ::Identity{ScalarMultiplicationGroupOperation},
    T::Type,
) where {𝔽}
    return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(
    G::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    X,
    ::Identity{ScalarMultiplicationGroupOperation},
) where {𝔽}
    return zero_vector!(LieAlgebra(G), X)
end
