#
#circle group represented in ℝ mod 2π = [-π, π), operation: addition mod 2π
#
function CircleGroup(M::Circle{ℝ})
    return LieGroup{ℝ,AdditionGroupOperation,typeof(M)}(M, AdditionGroupOperation())
end
CircleGroup(𝔽::AbstractNumbers=ℂ) = CircleGroup(Circle(𝔽))

const _RealCircleGroup = LieGroup{ℝ,AdditionGroupOperation,<:Circle{ℝ}}

_doc_sym_rem = """
    sym_rem(x,[T=π])

Compute symmetric remainder of `x` with respect to the interall 2*`T`, i.e.
`(x+T)%2T`, where the default for `T` is ``π``.
"""
@doc "$(_doc_sym_rem)"
function sym_rem(x::N, T=π) where {N<:Number}
    return (x ≈ T ? convert(N, -T) : rem(x, convert(N, 2 * T), RoundNearest))
end
function sym_rem(x, T=π)
    return map(sym_rem, x, Ref(T))
end

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

_doc_conjugate_circle_group = """
    conjugate(CircleGroup, g, h)
    conjugate!(CircleGroup, k, g, h)

Compute the conjugation map ``c_g: $(_math(:G)) → $(_math(:G))`` given by
``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = h``.
It simplifies to the identity since the group operation on the circle group is abelian.

This can be computed in-place of `k` if `k` is `mutable`.
"""

@doc "$(_doc_conjugate_circle_group)"
conjugate(::_RealCircleGroup, g, h) = g

@doc "$(_doc_conjugate_circle_group)"
conjugate!(::_RealCircleGroup, k, g, ::Any) = copyto!(k, g)

_doc_diff_conjugate_real_circle = """
    diff_conjugate(G::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g, h, X)
    diff_conjugate!(G::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, Y, g, h, X)

Compute the differential of the [conjugation map](@ref conjugate) ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1}=h``.
On the circle group represented as [part of the real line](@ref circle-group-real), this simplifies to ``D(c_g(h))[X] = X``.

This can be computed in-place of `Y` if `Y` is `mutable`.
"""

@doc "$(_doc_diff_conjugate_real_circle)"
diff_conjugate(::_RealCircleGroup, g, h, X::Number) = X

_doc_diff_inv_real_circle = """
    diff_inv(G::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g, X)
    diff_inv!(G::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, Y, g, X)

Compute the the differential ``Dι_{$(_math(:G))}([g])[X]`` of the inversion ``ι_{$(_math(:G))}([g]) := [g]^{-1} = [-g]`` at ``X ∈ 𝔤``
in the [`LieAlgebra`](@ref) ``𝔤`` of the [real `CircleGroup`](@ref circle-group-real) `G` ``=$(_math(:G))``.

The computation simplifies due to commutativity to

```math
Dι_{$(_math(:G))}([g])[X] = -X.
```

This can be computed in-place of `Y` if `Y` is `mutable`.
"""

@doc "$(_doc_diff_inv_real_circle)"
diff_inv(::_RealCircleGroup, g, X) = -X

@doc "$(_doc_diff_inv_real_circle)"
function diff_inv!(G::_RealCircleGroup, Y, g, X)
    return copyto!(LieAlgebra(G), Y, -X)
end

diff_left_compose(::_RealCircleGroup, g, h, X::Number) = X

diff_right_compose(::_RealCircleGroup, g, h, X::Number) = X

_doc_exp_real_circ = """
    exp(::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, X)
    exp!(::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g, X)

Compute the Lie group exponential of a vector `X` of the [`LieAlgebra`](@ref)
of the circle group, represented as angles in ``[-π, π)``.
In that case, the Lie algebra is the real line and the Lie group exponential of
a vector ``X ∈ ℝ`` is its equivalence class
```math
    $(_tex(:exp))(X) = [X] ∈ $(_tex(:SetDef, "[x] ∈ ℝ / 2πℤ", "x ∈ [-π,π)", "big")).
```

This can be computed in-place of `g`.
"""

@doc "$(_doc_exp_real_circ)"
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

@doc "$(_doc_exp_real_circ)"
function ManifoldsBase.exp!(G::_RealCircleGroup, g, X)
    return copyto!(g, exp(G, X))
end
function ManifoldsBase.exp!(G::_RealCircleGroup, h, g, X)
    return copyto!(h, exp(G, g, X))
end

# This can be combined with the functions above once we only have one circle group const
#
function get_coordinates_lie(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup},
    X::T,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G, T), X, DefaultOrthonormalBasis(ℝ))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup},
    c,
    X::T,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G, T), X, DefaultOrthonormalBasis(ℝ))
end

function get_vector_lie(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    T::Type=Float64,
)
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(ℝ))
end
function get_vector_lie!(
    𝔤::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup},
    X::T,
    c,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
) where {T}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(ℝ))
end

identity_element(::_RealCircleGroup) = 0.0

function ManifoldsBase.inner(
    ::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup}, X, Y
)
    return dot(X, Y)
end

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

function lie_bracket(
    ::LieAlgebra{ℝ,AdditionGroupOperation,<:_RealCircleGroup}, X::Any, ::Any
)
    return zero(X)
end

_doc_log_real_circ = """
    log(::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, g)
    log!(::LieGroup{ℝ, AdditionGroupOperation, Circle{ℝ}}, X, g)

Compute the Lie group logarithm on the [`CircleGroup`](@ref circle-group-real), represented as angles in ``[-π,π)``.
The [`LieAlgebra`](@ref) is the real line and ``$(_tex(:log))`` is given by the identity map.

Formally ``$(_tex(:log))`` promotes an equivalence class ``[X]`` to a representative ``X∈ℝ``.

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_real_circ)"
ManifoldsBase.log(::_RealCircleGroup, g::Number) = g

function ManifoldsBase.log(G::_RealCircleGroup, g)
    return map(gg -> log(G, gg), g)
end

function ManifoldsBase.log(G::_RealCircleGroup, g, h)
    return log(G, compose(G, inv(G, g), h))
end

@doc "$(_doc_log_real_circ)"
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

function Base.show(io::IO, ::_RealCircleGroup)
    return print(io, "CircleGroup(ℝ)")
end
