#
#circle group represented in ℝ mod 2π = [-π, π), operation: addition mod 2π
#
function CircleGroup(M::Manifolds.Circle{ℝ})
    return CircleGroup{ℝ,AdditionGroupOperation,typeof(M)}(M, AdditionGroupOperation())
end

#construct CircleGroup(Circle(ℝ)) or CircleGroup(Circle(ℂ)) by just using the field as input, ℂ is the default
CircleGroup(𝔽::ManifoldsBase.AbstractNumbers=ℂ) = CircleGroup(Circle(𝔽))

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

function _compose(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, p::Number, q::Number)
    return sym_rem(p + q)
end
function _compose(
    G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}},
    p::AbstractArray{<:Any,0},
    q::AbstractArray{<:Any,0},
)
    return map((pp, qq) -> compose(G, pp, qq), p, q)
end

function _compose!(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, x, p, q)
    return copyto!(x, compose(G, p, q))
end

conjugate(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, h) = g
conjugate!(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, k, g, ::Any) = copyto!(k, g)

diff_conjugate(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, h, X::Number) = X

diff_inv(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, X) = -X
function diff_inv(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, Y, g, X)
    return copyto!(LieAlgebra(G), Y, -X)
end

diff_left_compose(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, h, X::Number) = X

diff_right_compose(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, h, X::Number) = X

ManifoldsBase.exp(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, X::Number) = sym_rem(X)
function ManifoldsBase.exp(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, X)
    return map(XX -> exp(G, XX), X)
end
function ManifoldsBase.exp(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g::Number, X::Number
)
    return sym_rem(g + X)
end
function ManifoldsBase.exp(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, X)
    return map((gg, XX) -> exp(G, gg, XX), g, X)
end

function ManifoldsBase.exp!(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, X)
    return copyto!(g, exp(G, X))
end
function ManifoldsBase.exp!(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, h, g, X)
    return copyto!(h, exp(G, g, X))
end

# This can be combined with the functions above once we only have one circle group const
#
function get_coordinates_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}},
    X,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}},
    c,
    X,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end

function get_vector_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
    T::Type=Float64,
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}},
    X::T,
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
) where {𝔽,Op<:AbstractGroupOperation,T,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

identity_element(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}) = 0.0
function identity_element(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, p::Union{<:Number,Type{<:Number}}
)
    return zero(p)
end
function identity_element(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, ::Type{<:SArray{S,T}}
) where {S,T}
    @SArray fill(one(T))
end
function identity_element(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, ::Type{<:MArray{S,T}}
) where {S,T}
    @MArray fill(one(T))
end

Base.inv(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, p::Number) = sym_rem(-p)

function inv_left_compose(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g::Number, h::Number
)
    return sym_rem(-g + h)
end

function inv_right_compose(
    ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g::Number, h::Number
)
    return sym_rem(g - h)
end

function lie_bracket(
    ::LieAlgebra{ℝ,AdditionGroupOperation,CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}},
    X::Any,
    ::Any,
)
    return zero(X)
end

ManifoldsBase.log(::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g::Number) = g
function ManifoldsBase.log(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g)
    return map(gg -> log(G, gg), g)
end
function ManifoldsBase.log(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, g, h)
    return log(G, compose(G, inv(G, g), h))
end
function ManifoldsBase.log!(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, X, g)
    return copyto!(X, log(G, g))
end
function ManifoldsBase.log!(G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, X, g, h)
    return copyto!(X, log(G, g, h))
end

function ManifoldsBase.log(
    G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}}, ::Identity{AdditionGroupOperation}
)
    return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(
    G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}},
    ::Identity{AdditionGroupOperation},
    T::Type,
)
    return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(
    G::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}},
    X,
    ::Identity{AdditionGroupOperation},
)
    return zero_vector!(LieAlgebra(G), X)
end

_doc_exp_real_circ = """
    exp(::CircleGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, X)
    exp!(::CircleGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g, X)

The Lie group exponential on the [`CircleGroup`](@ref) represented in ℝ is given by the projection into the equivalence class of its defining relation.

This can be computed in-place of `X`.
"""

function Base.show(io::IO, ::CircleGroup{ℝ,AdditionGroupOperation,Circle{ℝ}})
    return print(io, "CircleGroup(ℝ)")
end
