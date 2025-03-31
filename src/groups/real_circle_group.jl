"""
    RealCircleGroup = LieGroup{ℝ, AdditionGroupOperation, Manifolds.Circle{ℝ}}

The real circle group ``𝕊^1`` is the set of points on the unit circle in ``ℝ^2``, represented by its angles  ``x ∈ [-π,π)``.
It is equipped with the group operation of adding angles ``$(_tex(:rm, raw"mod\, ")) 2π`` via [`AdditionGroupOperation`](@ref).

It it is obtained as a quotient space of the real numbers

```math
 𝕊¹ := ℝ / 2πℤ = $(_tex(:SetDef, "[x] ∈ ℝ / 2πℤ", "x ∈ [-π,π)", "big")).
```

It can be identified with the set of complex numbers of absolute value 1, i.e. the one dimensional [`Sphere`](@extref `Manifolds.Sphere`),
together with the group operation of multiplying two complex numbers. For that construction see [`CircleGroup`](@ref).

The (real) circle group is a one dimensional Riemannian manifold and a Lie group.

# Constructor

    RealCircleGroup()

Generate the real circle group.
"""
const RealCircleGroup = LieGroup{ℝ,AdditionGroupOperation,Manifolds.Circle{ℝ}}

function RealCircleGroup()
    circ = Manifolds.Circle(ℝ)
    return RealCircleGroup(circ, AdditionGroupOperation())
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

_compose(::RealCircleGroup, p::Number, q::Number) = sym_rem(p + q)
_compose(G::RealCircleGroup, p::AbstractArray{<:Any,0}, q::AbstractArray{<:Any,0}) = map((pp, qq) -> compose(G, pp, qq), p, q)

_compose!(G::RealCircleGroup, x, p, q) = copyto!(x, compose(G, p, q))

conjugate(::RealCircleGroup, g, h) = g
conjugate!(::RealCircleGroup, k, g, ::Any) = copyto!(k, g)

diff_conjugate(::RealCircleGroup, g, h, X::Number) = X

diff_inv(::RealCircleGroup, g, X) = -X
diff_inv(G::RealCircleGroup, Y, g, X) = copyto!(LieAlgebra(G), Y, -X)

diff_left_compose(::RealCircleGroup, g, h, X::Number) = X

diff_right_compose(::RealCircleGroup, g, h, X::Number) = X

ManifoldsBase.exp(::RealCircleGroup, X::Number) = sym_rem(X)
ManifoldsBase.exp(G::RealCircleGroup, X) = map(XX-> exp(G, XX), X)
ManifoldsBase.exp(::RealCircleGroup, g::Number, X::Number) = sym_rem(g + X)
ManifoldsBase.exp(G::RealCircleGroup, g, X) = map((gg, XX) -> exp(G, gg, XX), g, X)

ManifoldsBase.exp!(G::RealCircleGroup, g, X) = copyto!(g, exp(G, X))
ManifoldsBase.exp!(G::RealCircleGroup, h, g, X) = copyto!(h, exp(G, g, X))


# This can be combined with the functions above once we only have one circle group const
#
function get_coordinates_lie(
    𝔤::LieAlgebra{𝔽,Op,RealCircleGroup}, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{𝔽,Op,RealCircleGroup}, c, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end

function get_vector_lie(
    𝔤::LieAlgebra{𝔽,Op,RealCircleGroup},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
    T::Type=Float64,
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
    𝔤::LieAlgebra{𝔽,Op,RealCircleGroup}, X::T, c, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,T,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

identity_element(::RealCircleGroup) = 0.0
identity_element(::RealCircleGroup, p::Union{<:Number,Type{<:Number}}) = zero(p)
function identity_element(::RealCircleGroup, ::Type{<:SArray{S,T}}) where {S,T}
  @SArray fill(one(T))
end
function identity_element(::RealCircleGroup, ::Type{<:MArray{S,T}}) where {S,T}
  @MArray fill(one(T))
end

Base.inv(::RealCircleGroup, p::Number) = sym_rem(-p)

inv_left_compose(::RealCircleGroup, g::Number, h::Number) = sym_rem(-g + h)

inv_right_compose(::RealCircleGroup, g::Number, h::Number) = sym_rem(g - h)

lie_bracket(::LieAlgebra{ℝ, AdditionGroupOperation, RealCircleGroup}, X::Any, ::Any) = zero(X)







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