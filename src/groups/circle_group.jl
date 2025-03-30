"""
    CircleGroup = LieGroup{ℂ, ScalarMultiplicationGroupOperation, Manifolds.Circle{ℂ}}

The complex circle group ``𝕊^1`` is the set of complex numbers ``z ∈ ℂ`` of absolute value ``1``

```math
𝕊¹ := $(_tex(:SetDef, "z ∈ ℂ", "|z| = 1", "big")) = $(_tex(:SetDef, "a + bi ∈ ℂ", "a^2+b^2 = 1", "big")),
```

equipped with the group operation of complex multiplication [`ScalarMultiplicationGroupOperation`](@ref).

It can be identified with the unit circle in ``ℝ^2``, i.e. the one dimensional [`Sphere`](@extref `Manifolds.Sphere`),
together with the group operation of adding the angles of two points on the circle. For that construction see [`RealCircleGroup`](@ref).

The (complex) circle group is a one dimensional Riemannian manifold and a Lie group.

# Constructor

    CircleGroup()

Generate the complex circle group.
"""
const CircleGroup = LieGroup{ℂ,ScalarMultiplicationGroupOperation,Manifolds.Circle{ℂ}}

function CircleGroup()
    circ = Manifolds.Circle(ℂ)
    return CircleGroup(circ, ScalarMultiplicationGroupOperation())
end

function diff_left_compose(::CircleGroup, g::Number, h::Any, X::Number)
    return g * X
end

function diff_right_compose(::CircleGroup, g::Number, h::Any, X::Number)
    return X * g
end

_doc_exp_complex_circ = """
    exp(::CircleGroup, X)
    exp!(::CircleGroup, g, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```
"""

@doc "$(_doc_exp_complex_circ)"
Base.exp(::CircleGroup, X::Number) = exp(X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup, g, X)

function get_coordinates_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup}, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup}, c, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
    T::Type=ComplexF64,
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup}, X::T, c, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,T,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

_doc_log_complex_circ = """
    log(::CircleGroup, g)
    log!(::CircleGroup, X, g)

Compute the Lie group logarithm on the complex [`CircleGroup`](@ref), which coincides with the
ordinary complex logarithm.
"""

identity_element(::CircleGroup) = 1.0 + 0.0im
identity_element(::CircleGroup, T::Union{<:Number,Type{<:Number}}) = one(T)
identity_element(::CircleGroup, ::Type{<:SArray{S,T}}) where {S,T} = @SArray fill(one(T))
identity_element(::CircleGroup, ::Type{<:MArray{S,T}}) where {S,T} = @MArray fill(one(T))

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log(::CircleGroup, g)

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log!(M::CircleGroup, X, g)

function ManifoldsBase.log(::CircleGroup, g::Number)
    return log(g)
end

function Base.show(io::IO, ::CircleGroup)
    return print(io, "CircleGroup()")
end

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


