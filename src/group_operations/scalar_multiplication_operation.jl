"""
    ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a scalar multiplication.
"""
struct ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

function get_number_type(x::Number)
    return typeof(x)
end

function get_number_type(x::Array{<:Number,0})
    return typeof(x[])
end

function get_number_type(x::Ref{<:Number})
    return typeof(x[])
end

get_num(x::Number) = identity(x)
function get_num(x::Array{<:Number,0})
    return x[]
end
function get_num(x::Ref{<:Number})
    return x[]
end

function compose(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * h
end
function compose(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    return get_num(g) * get_num(h)
end

function compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{T},Array{T,0}},
    g::T1,
    h::T2,
) where {𝔽,T1<:Number,T2<:Number,T<:Number}
    k[] = ((k === g || k === h) ? copy(g * h) : g * h)
    return k
end

function compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{T},Array{T,0}},
    g::Union{Ref{T1},Array{T1,0}},
    h::T2,
) where {𝔽,T1<:Number,T2<:Number,T<:Number}
    k[] = ((k === g || k === h) ? copy(g[] * h) : g[] * h)
    return k
end

function compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{T},Array{T,0}},
    g::T1,
    h::Union{Ref{T2},Array{T2,0}},
) where {𝔽,T1<:Number,T2<:Number,T<:Number}
    k[] = ((k === g || k === h) ? copy(g * h[]) : g * h[])
    return k
end

function compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{T},Array{T,0}},
    g::Union{Ref{T1},Array{T1,0}},
    h::Union{Ref{T2},Array{T2,0}},
) where {𝔽,T1<:Number,T2<:Number,T<:Number}
    k[] = ((k === g || k === h) ? copy(g[] * h[]) : g[] * h[])
    return k
end

function conjugate(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * h * inv(g)
end

function conjugate(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    return get_num(g) * get_num(h) * inv(get_num(g))
end

_doc_exp_scalar_mult = """
    exp(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e::Identity{ScalarMultiplicationGroupOperation}, X, t::Number=1)
    exp!(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g, e::Identity{ScalarMultiplicationGroupOperation}, X, t::Number=1)

Compute the Lie group exponential on a [`LieGroup`](@ref) with a [`ScalarMultiplicationGroupOperation`](@ref),
which simplifies to the [ordinary exponential](https://en.wikipedia.org/wiki/Matrix_exponential).

This can be computed in-place of `g`.
"""

@doc "$(_doc_exp_scalar_mult)"
ManifoldsBase.exp(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    X::Union{<:Number,Ref{<:Number},<:Array{<:Number,0}},
    t::Number=1,
) where {𝔽} = fill(exp(t * (X isa Number ? X : X[])))

@doc "$(_doc_exp_scalar_mult)"
function ManifoldsBase.exp!(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::Union{Ref{<:Number},<:Array{<:Number,0}},
    X::Union{<:Number,Ref{<:Number},<:Array{<:Number,0}},
    t::Number=1,
) where {𝔽}
    g[] = exp(t * (X isa Number ? X : X[]))
    return g
end

Base.inv(::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, g::Number) where {𝔽} = inv(g)
function Base.inv(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    g::Union{Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    return inv(g[])
end

function inv!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    h::Union{Ref{<:Number},Array{<:Number,0}},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    h[] = inv(get_num(h))
    return h
end

function inv_left_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return inv(g) * h
end
function inv_left_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    return inv(get_num(g)) * get_num(h)
end

function inv_left_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{<:Number},Array{<:Number,0}},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    k[] = inv_left_compose(G, g, h)
    return k
end

function inv_right_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, g::Number, h::Number
) where {𝔽}
    return g * inv(h)
end
function inv_right_compose(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    return get_num(g) * inv(get_num(h))
end

function inv_right_compose!(
    G::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation},
    k::Union{Ref{<:Number},Array{<:Number,0}},
    g::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
    h::Union{<:Number,Ref{<:Number},Array{<:Number,0}},
) where {𝔽}
    k[] = inv_right_compose(G, g, h)
    return k
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
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e::T
) where {𝔽,T<:Number}
    return one(e)
end

@doc "$(_doc_identity_element_scalar_mult)"
identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}
    return fill!(e, 1)
end
