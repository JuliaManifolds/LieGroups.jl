#
#circle group represented in ℝ mod 2π = [-π, π), operation: addition mod 2π
#
function CircleGroup(M::Manifolds.Circle{ℝ})  
    return CircleGroup{ℝ, typeof(M)}(
        M, AdditionGroupOperation()
    )
end
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

_compose(::CircleGroup{ℝ, Circle{ℝ}}, p::Number, q::Number) = sym_rem(p + q)
_compose(G::CircleGroup{ℝ, Circle{ℝ}}, p::AbstractArray{<:Any,0}, q::AbstractArray{<:Any,0}) = map((pp, qq) -> compose(G, pp, qq), p, q)

_compose!(G::CircleGroup{ℝ, Circle{ℝ}}, x, p, q) = copyto!(x, compose(G, p, q))

conjugate(::CircleGroup{ℝ, Circle{ℝ}}, g, h) = g
conjugate!(::CircleGroup{ℝ, Circle{ℝ}}, k, g, ::Any) = copyto!(k, g)

diff_conjugate(::CircleGroup{ℝ, Circle{ℝ}}, g, h, X::Number) = X

diff_inv(::CircleGroup{ℝ, Circle{ℝ}}, g, X) = -X
diff_inv(G::CircleGroup{ℝ, Circle{ℝ}}, Y, g, X) = copyto!(LieAlgebra(G), Y, -X)

diff_left_compose(::CircleGroup{ℝ, Circle{ℝ}}, g, h, X::Number) = X

diff_right_compose(::CircleGroup{ℝ, Circle{ℝ}}, g, h, X::Number) = X

ManifoldsBase.exp(::CircleGroup{ℝ, Circle{ℝ}}, X::Number) = sym_rem(X)
ManifoldsBase.exp(G::CircleGroup{ℝ, Circle{ℝ}}, X) = map(XX-> exp(G, XX), X)
ManifoldsBase.exp(::CircleGroup{ℝ, Circle{ℝ}}, g::Number, X::Number) = sym_rem(g + X)
ManifoldsBase.exp(G::CircleGroup{ℝ, Circle{ℝ}}, g, X) = map((gg, XX) -> exp(G, gg, XX), g, X)

ManifoldsBase.exp!(G::CircleGroup{ℝ, Circle{ℝ}}, g, X) = copyto!(g, exp(G, X))
ManifoldsBase.exp!(G::CircleGroup{ℝ, Circle{ℝ}}, h, g, X) = copyto!(h, exp(G, g, X))


# This can be combined with the functions above once we only have one circle group const
#
function get_coordinates_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ, Circle{ℝ}}}, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ, Circle{ℝ}}}, c, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end

function get_vector_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ, Circle{ℝ}}},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
    T::Type=Float64,
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℝ, Circle{ℝ}}}, X::T, c, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,T,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

identity_element(::CircleGroup{ℝ, Circle{ℝ}}) = 0.0
identity_element(::CircleGroup{ℝ, Circle{ℝ}}, p::Union{<:Number,Type{<:Number}}) = zero(p)
function identity_element(::CircleGroup{ℝ, Circle{ℝ}}, ::Type{<:SArray{S,T}}) where {S,T}
  @SArray fill(one(T))
end
function identity_element(::CircleGroup{ℝ, Circle{ℝ}}, ::Type{<:MArray{S,T}}) where {S,T}
  @MArray fill(one(T))
end

Base.inv(::CircleGroup{ℝ, Circle{ℝ}}, p::Number) = sym_rem(-p)

inv_left_compose(::CircleGroup{ℝ, Circle{ℝ}}, g::Number, h::Number) = sym_rem(-g + h)

inv_right_compose(::CircleGroup{ℝ, Circle{ℝ}}, g::Number, h::Number) = sym_rem(g - h)

lie_bracket(::LieAlgebra{ℝ, AbstractGroupOperation, CircleGroup{ℝ, Circle{ℝ}}}, X::Any, ::Any) = zero(X)







ManifoldsBase.log(::CircleGroup{ℝ, Circle{ℝ}}, g::Number) = g
ManifoldsBase.log(G::CircleGroup{ℝ, Circle{ℝ}}, g) = map(gg -> log(G, gg), g)
ManifoldsBase.log(G::CircleGroup{ℝ, Circle{ℝ}}, g, h) = log(G, compose(G, inv(G, g), h))
ManifoldsBase.log!(G::CircleGroup{ℝ, Circle{ℝ}}, X, g) = copyto!(X, log(G, g))
ManifoldsBase.log!(G::CircleGroup{ℝ, Circle{ℝ}}, X, g, h) = copyto!(X, log(G, g, h))

function ManifoldsBase.log(
  G::CircleGroup{ℝ, Circle{ℝ}},
  ::Identity{AdditionGroupOperation},
)
  return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(
  G::CircleGroup{ℝ, Circle{ℝ}},
  ::Identity{AdditionGroupOperation},
  T::Type,
)
  return zero_vector(LieAlgebra(G), T)
end

function ManifoldsBase.log!(G::CircleGroup{ℝ, Circle{ℝ}}, X,::Identity{AdditionGroupOperation},)
  return zero_vector!(LieAlgebra(G), X)
end



_doc_exp_real_circ = """
    exp(::CircleGroup{ℝ, Circle{ℝ}}, X)
    exp!(::CircleGroup{ℝ, Circle{ℝ}}, g, X)

The Lie group exponential on the [`CircleGroup`](@ref) represented in ℝ is given by the projection into the equivalence class of its defining relation.

This can be computed in-place of `X`.
"""

function Base.show(io::IO, ::CircleGroup{ℝ, Circle{ℝ}})
    return print(io, "CircleGroup(ℝ)")
end