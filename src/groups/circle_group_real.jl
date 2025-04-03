#
#circle group represented in ℝ mod 2π = [-π, π), operation: addition mod 2π
#
function CircleGroup(M::Manifolds.Circle{ℝ})
    return LieGroup{ℝ,AdditionGroupOperation,typeof(M)}(M, AdditionGroupOperation())
end

#construct CircleGroup(Circle(ℝ)) or CircleGroup(Circle(ℂ)) by just using the field as input, ℂ is the default
CircleGroup(𝔽::ManifoldsBase.AbstractNumbers=ℂ) = CircleGroup(Circle(𝔽))

const _RealCircleGroup = LieGroup{ℝ,AdditionGroupOperation,<:Circle{ℝ}}

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

function _compose(::_RealCircleGroup, p::Number, q::Number)
    return sym_rem(p + q)
end
function _compose(G::_RealCircleGroup, p::AbstractArray{<:Any,0}, q::AbstractArray{<:Any,0})
    return map((pp, qq) -> compose(G, pp, qq), p, q)
end
function _compose(::_RealCircleGroup, p::Number, q::AbstractArray{<:Any,0})
    return p .+ q
end
function _compose(::_RealCircleGroup, p::AbstractArray{<:Any,0}, q::Number)
    return p .+ q
end

function _compose!(G::_RealCircleGroup, x, p, q)
    return copyto!(x, compose(G, p, q))
end

conjugate(::_RealCircleGroup, g, h) = g
conjugate!(::_RealCircleGroup, k, g, ::Any) = copyto!(k, g)

diff_conjugate(::_RealCircleGroup, g, h, X::Number) = X

diff_inv(::_RealCircleGroup, g, X) = -X
function diff_inv(G::_RealCircleGroup, Y, g, X)
    return copyto!(LieAlgebra(G), Y, -X)
end

diff_left_compose(::_RealCircleGroup, g, h, X::Number) = X

diff_right_compose(::_RealCircleGroup, g, h, X::Number) = X

ManifoldsBase.exp(::_RealCircleGroup, X::Number) = sym_rem(X)
function ManifoldsBase.exp(G::_RealCircleGroup, X)
    return map(XX -> exp(G, XX), X)
end
function ManifoldsBase.exp(::_RealCircleGroup, g::Number, X::Number)
    return sym_rem(g + X)
end
function ManifoldsBase.exp(G::_RealCircleGroup, g, X)
    return map((gg, XX) -> exp(G, gg, XX), g, X)
end

function ManifoldsBase.exp!(G::_RealCircleGroup, g, X)
    return copyto!(g, exp(G, X))
end
function ManifoldsBase.exp!(G::_RealCircleGroup, h, g, X)
    return copyto!(h, exp(G, g, X))
end

# This can be combined with the functions above once we only have one circle group const
#
function get_coordinates_lie(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,_RealCircleGroup},
    X::T,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G, T), X, DefaultOrthonormalBasis(ℝ))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,_RealCircleGroup},
    c,
    X::T,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G, T), X, DefaultOrthonormalBasis(ℝ))
end

function get_vector_lie(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,_RealCircleGroup},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    T::Type=Float64,
)
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(ℝ))
end
function get_vector_lie!(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,_RealCircleGroup},
    X::T,
    c,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(ℝ))
end

identity_element(::_RealCircleGroup) = 0.0

Base.inv(::_RealCircleGroup, p::Number) = sym_rem(-p)

function inv_left_compose(::_RealCircleGroup, g::Number, h::Number)
    return sym_rem(-g + h)
end

function inv_right_compose(::_RealCircleGroup, g::Number, h::Number)
    return sym_rem(g - h)
end

function ManifoldsBase.isapprox(::_RealCircleGroup, p, X, Y; kwargs...)
    return isapprox(X[], Y[]; kwargs...)
end

function lie_bracket(::LieAlgebra{ℝ,AdditionGroupOperation,_RealCircleGroup}, X::Any, ::Any)
    return zero(X)
end

ManifoldsBase.log(::_RealCircleGroup, g::Number) = g
function ManifoldsBase.log(G::_RealCircleGroup, g)
    return map(gg -> log(G, gg), g)
end
function ManifoldsBase.log(G::_RealCircleGroup, g, h)
    return log(G, compose(G, inv(G, g), h))
end
function ManifoldsBase.log!(G::_RealCircleGroup, X, g)
    return copyto!(X, log(G, g))
end
function ManifoldsBase.log!(G::_RealCircleGroup, X, g, h)
    return copyto!(X, log(G, g, h))
end

function ManifoldsBase.log(G::_RealCircleGroup, ::Identity{AdditionGroupOperation})
    return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(G::_RealCircleGroup, ::Identity{AdditionGroupOperation}, T::Type)
    return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(G::_RealCircleGroup, X, ::Identity{AdditionGroupOperation})
    return zero_vector!(LieAlgebra(G), X)
end

_doc_exp_real_circ = """
    exp(::CircleGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, X)
    exp!(::CircleGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g, X)

The Lie group exponential on the [`CircleGroup`](@ref) represented in ℝ is given by the projection into the equivalence class of its defining relation.

This can be computed in-place of `X`.
"""

function Base.show(io::IO, ::_RealCircleGroup)
    return print(io, "CircleGroup(ℝ)")
end
