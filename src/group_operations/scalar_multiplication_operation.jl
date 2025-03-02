"""
    ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a scalar multiplication.
"""
struct ScalarMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

function _compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k::Union{Ref{T}, Array{T,0}}, g::T1, h::T2
) where {𝔽, T1<:Number, T2<:Number, T<:Number}
    k[] = ((k === g || k === h) ? copy(g * h) : g * h)
    return k
end

function _compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k::Union{Ref{T}, Array{T,0}}, g::Union{Ref{T1}, Array{T1,0}}, h::T2
) where {𝔽, T1<:Number, T2<:Number, T<:Number}
    k[] = ((k === g || k === h) ? copy(g[] * h) : g[] * h)
    return k
end

function _compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k::Union{Ref{T}, Array{T,0}}, g::T1, h::Union{Ref{T2}, Array{T2,0}}
) where {𝔽, T1<:Number, T2<:Number, T<:Number}
    k[] = ((k === g || k === h) ? copy(g * h[]) : g * h[])
    return k
end

function _compose!(
    ::LieGroup{𝔽,<:ScalarMultiplicationGroupOperation}, k::Union{Ref{T}, Array{T,0}}, g::Union{Ref{T1}, Array{T1,0}}, h::Union{Ref{T2}, Array{T2,0}}
) where {𝔽, T1<:Number, T2<:Number, T<:Number}
    k[] = ((k === g || k === h) ? copy(g[] * h[]) : g[] * h[])
    return k
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
    X::Union{<:Number, Ref{<:Number}, <:Array{<:Number, 0}},
    t::Number=1,
) where {𝔽} = fill(exp(t * (X isa Number ? X : X[])))

@doc "$(_doc_exp_scalar_mult)"
function ManifoldsBase.exp!(
    ::LieGroup{𝔽,ScalarMultiplicationGroupOperation},
    g::Union{Ref{<:Number}, <:Array{<:Number, 0}},
    X::Union{<:Number, Ref{<:Number}, <:Array{<:Number, 0}},
    t::Number=1,
) where {𝔽}
    g[] = exp(t * (X isa Number ? X : X[]))
    return g
end

_doc_identity_element_scalar_mult = """
    identity_element(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`ScalarMultiplicationGroupOperation`](@ref) is the one-element.
"""

@doc "$(_doc_identity_element_scalar_mult)"
identity_element(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}) where {𝔽} = 1.0

@doc "$(_doc_identity_element_scalar_mult)"
identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(::LieGroup{𝔽,ScalarMultiplicationGroupOperation}, e) where {𝔽}
    return e[] = 1.0
end