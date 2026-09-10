#
#circle group represented by complex numbers, operation: complex multiplication
#
function CircleGroup(M::Circle{ℂ})
    return LieGroup{ℂ, AbelianMultiplicationGroupOperation, typeof(M)}(
        M, AbelianMultiplicationGroupOperation()
    )
end

const _ComplexCircleGroup = LieGroup{ℂ, AbelianMultiplicationGroupOperation, <:Circle{ℂ}}

function LieAlgebra(G::_ComplexCircleGroup)
    return LieAlgebra{ℝ, AbelianMultiplicationGroupOperation, typeof(G)}(
        G, Identity(G), ManifoldsBase.TangentSpaceType()
    )
end

_doc_diff_left_compose_complex_circ = """
    diff_left_compose(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, g, h, X)
    diff_left_compose!(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, Y, g, h, X)

Compute the differential of the group operation ``gh`` with respeect to the lefy argument `g`.

The formula reads

```math
   $(_math(:d)) ρ_h(g) = X.
```

since the group operation is abelian. This can be computed in-place of `Y` if `Y` is mutable.
"""

@doc "$(_doc_diff_left_compose_complex_circ)"
function diff_left_compose(::_ComplexCircleGroup, g::Number, h::Any, X::Number)
    return X
end

_doc_diff_right_compose_complex_circ = """
    diff_right_compose(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, g, h, X)
    diff_right_compose!(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its second (right) argument `h`.

Another interpretation is to consider a function where we do a fixed multiplication from the left with `g`.
i..e. the left group multiplication function ``λ_g(h) = g$(_math(:∘))h``.

It simplifies for the [`AbelianMultiplicationGroupOperation`](@ref) to ``$(_math(:d))λ_g(h)[X] = X``.
This can be computed in-place of `Y`.
"""

@doc "$(_doc_diff_right_compose_complex_circ)"
function diff_right_compose(::_ComplexCircleGroup, g::Number, h::Any, X::Number)
    return X
end

_doc_exp_complex_circ = """
    exp(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, X)
    exp!(::LieGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, g, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```
"""

@doc "$(_doc_exp_complex_circ)"
Base.exp(::_ComplexCircleGroup, X::Number) = exp(X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::_ComplexCircleGroup, g, X)

function get_coordinates_lie(
        𝔤::LieAlgebra{ℝ, AbelianMultiplicationGroupOperation, <:_ComplexCircleGroup},
        X::T,
        ::DefaultLieAlgebraOrthogonalBasis{𝔽},
    ) where {T, 𝔽}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G, T), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
        𝔤::LieAlgebra{ℝ, AbelianMultiplicationGroupOperation, <:_ComplexCircleGroup},
        c,
        X,
        ::DefaultLieAlgebraOrthogonalBasis{𝔽},
    ) where {𝔽}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie(
        𝔤::LieAlgebra{ℝ, AbelianMultiplicationGroupOperation, <:_ComplexCircleGroup},
        c,
        ::DefaultLieAlgebraOrthogonalBasis{𝔽},
        T::Type = ComplexF64,
    ) where {𝔽}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
        𝔤::LieAlgebra{ℝ, AbelianMultiplicationGroupOperation, <:_ComplexCircleGroup},
        X::T,
        c,
        ::DefaultLieAlgebraOrthogonalBasis{𝔽},
    ) where {T, 𝔽}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

function identity_element(::_ComplexCircleGroup)
    return 1.0 + 0.0im
end

function ManifoldsBase.isapprox(::_ComplexCircleGroup, p, X, Y; kwargs...)
    return isapprox(X[], Y[]; kwargs...)
end

_doc_log_complex_circ = """
    log(::CircleGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, g)
    log!(::CircleGroup{ℂ, AbelianMultiplicationGroupOperation, Circle{ℂ}}, X, g)

Compute the Lie group logarithm on the complex [`CircleGroup`](@ref), which coincides with the
ordinary complex logarithm.

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log(::_ComplexCircleGroup, g)

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log!(M::_ComplexCircleGroup, X, g)

function ManifoldsBase.log(::_ComplexCircleGroup, g::Number)
    return log(g)
end

function ManifoldsBase.log!(::_ComplexCircleGroup, X, g)
    X[] = log(g[])
    return X
end

function ManifoldsBase.log!(
        G::_ComplexCircleGroup, X, ::Identity{AbelianMultiplicationGroupOperation}
    )
    return zero_vector!(LieAlgebra(G), X)
end

function Base.show(io::IO, ::_ComplexCircleGroup)
    return print(io, "CircleGroup()")
end
