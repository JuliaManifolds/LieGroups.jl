_doc_sym_rem = """
    sym_rem(x,[T=π])

Compute symmetric remainder of `x` with respect to the interall 2*`T`, i.e.
`(x+T)%2T`, where the default for `T` is ``π``
"""
@doc "$(_doc_sym_rem)"
function sym_rem(x::N, T=π) where {N<:Number}
    return (x ≈ T ? convert(N, -T) : rem(x, convert(N, 2 * T), RoundNearest))
end
sym_rem(x, T=π) = map(sym_rem, x, Ref(T))

_compose(::RealCircleGroup, p::Number, q::Number) = sym_rem(p + q)
_compose(G::RealCircleGroup, p::AbstractArray{<:Any,0}, q::AbstractArray{<:Any,0}) = map((pp, qq) -> compose(G, pp, qq), p, q)

_compose!(G::RealCircleGroup, x, p, q) = copyto!(x, compose(G, p, q))

conjugate(::RealCircleGroup, g, h) = g
conjugate!(::RealCircleGroup, k, g, ::Any) = copyto!(k, g)

diff_conjugate(::RealCircleGroup, g, h, X::Number) = X

diff_inv(::RealCircleGroup, g, X) = -X
diff_inv(G::RealCircleGroup, Y, g, X) = copyto!(LieAlgebra(G), Y, -X)

diff_left_compose(::RealCircleGroup, g, h, X::Number) = X
diff_left_compose(G::RealCircleGroup, g, h, X) = map((gg, XX) -> diff_left_compose(G, gg, h, XX), g, X)
diff_left_compose!(G::RealCircleGroup, Y, g, h, X) = copyto!(LieAlgebra(G), Y, diff_left_compose(G, g, h, X))

diff_right_compose(::RealCircleGroup, g, h, X::Number) = X
diff_right_compose(G::RealCircleGroup, g, h, X) = map((gg, XX) -> diff_right_compose(G, gg, h, XX), g, X)
diff_right_compose!(G::RealCircleGroup, Y, g, h, X) = copyto!(LieAlgebra(G), Y, diff_right_compose(G, g, h, X))

ManifoldsBase.exp(::RealCircleGroup, X::Number) = sym_rem(X)
ManifoldsBase.exp(G::RealCircleGroup, X::AbstractArray{<:Any,0}) = map(XX-> exp(G, XX), X)
ManifoldsBase.exp(::RealCircleGroup, g::Number, X::Number) = sym_rem(g + X)
ManifoldsBase.exp(G::RealCircleGroup, g::AbstractArray{<:Any,0}, X::AbstractArray{<:Any,0}) = map((gg, XX) -> exp(G, gg, XX), g, X)

ManifoldsBase.exp!(G::RealCircleGroup, g, X) = copyto!(g, exp(G, X))
ManifoldsBase.exp!(G::RealCircleGroup, h, g, X) = copyto!(h, exp(G, g, X))

identity_element(::RealCircleGroup) = 0.0
identity_element(::RealCircleGroup, p::Union{<:Number,Type{<:Number}}) = zero(p)
function identity_element(::RealCircleGroup, ::Type{<:SArray{S,T}}) where {S,T}
  @SArray fill(one(T))
end
function identity_element(::RealCircleGroup, ::Type{<:MArray{S,T}}) where {S,T}
  @MArray fill(one(T))
end

Base.inv(::RealCircleGroup, p::Number) = sym_rem(-p)
Base.inv(G::RealCircleGroup, p::AbstractArray{<:Any,0}) = map(pp -> inv(G, pp), p)

inv_left_compose(::RealCircleGroup, g::Number, h::Number) = sym_rem(-g + h)
#inv_left_compose(G::RealCircleGroup, g::AbstractArray{<:Any,0}, h::AbstractArray{<:Any,0}) = map((gg,hh) -> inv_left_compose(G, gg, hh), g, h)
#inv_left_compose!(G::RealCircleGroup, x, g, h) = copyto!(x, inv_left_compose(G, g, h))

inv_right_compose(::RealCircleGroup, g::Number, h::Number) = sym_rem(g - h)
#inv_right_compose(G::RealCircleGroup, g::AbstractArray{<:Any,0}, h::AbstractArray{<:Any,0}) = map((gg,hh) -> inv_right_compose(G, gg, hh), g, h)
#inv_right_compose!(G::RealCircleGroup, x, g, h) = copyto!(x, inv_right_compose(G, g, h))


lie_bracket(::LieAlgebra{ℝ, AdditionGroupOperation, RealCircleGroup}, X::Any, ::Any) = zero(X)
lie_bracket!(::LieAlgebra{ℝ, AdditionGroupOperation, RealCircleGroup}, Z, X, ::Any) = copyto!(Z, zero(X))






ManifoldsBase.log(::RealCircleGroup, g::Number) = g
ManifoldsBase.log(G::RealCircleGroup, g) = map(gg -> log(G, gg), g)
ManifoldsBase.log(G::RealCircleGroup, g, h) = log(G, compose(G, inv(G, g), h))
ManifoldsBase.log!(G::RealCircleGroup, X, g) = copyto!(X, log(G, g))
ManifoldsBase.log!(G::RealCircleGroup, X, g, h) = copyto!(X, log(G, g, h))

function ManifoldsBase.log(
  G::RealCircleGroup,
  ::Identity{AdditionGroupOperation},
)
  return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(
  G::RealCircleGroup,
  ::Identity{AdditionGroupOperation},
  T::Type,
)
  return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(G::RealCircleGroup, X,::Identity{AdditionGroupOperation},)
  return zero_vector!(LieAlgebra(G), X)
end



_doc_exp_real_circ = """
    exp(::RealCircleGroup, X)
    exp!(::RealCircleGroup, g, X)

The Lie group exponential on the [`RealCircleGroup`](@ref) is given by the projection into the equivalence class of its defining relation.

This can be computed in-place of `X`.
"""

function Base.show(io::IO, ::RealCircleGroup)
    return print(io, "RealCircleGroup()")
end