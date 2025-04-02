#
#circle group represented by complex numbers, operation: complex multiplication
#
function CircleGroup(M::Manifolds.Circle{ℂ})
    return CircleGroup{ℂ, ScalarMultiplicationGroupOperation, typeof(M)}(
        M, ScalarMultiplicationGroupOperation()
    )
end 

function diff_left_compose(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g::Number, h::Any, X::Number)
    return g * X
end

function diff_right_compose(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g::Number, h::Any, X::Number)
    return X * g
end

_doc_exp_complex_circ = """
    exp(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, X)
    exp!(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```
"""

@doc "$(_doc_exp_complex_circ)"
Base.exp(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, X::Number) = exp(X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g, X)

function get_coordinates_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}}, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates(M, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_coordinates_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}}, c, X, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_coordinates!(M, c, identity_element(G), X, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}},
    c,
    ::DefaultLieAlgebraOrthogonalBasis{𝔾},
    T::Type=ComplexF64,
) where {𝔽,Op<:AbstractGroupOperation,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector(M, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end
function get_vector_lie!(
    𝔤::LieAlgebra{𝔽,Op,CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}}, X::T, c, ::DefaultLieAlgebraOrthogonalBasis{𝔾}
) where {𝔽,Op<:AbstractGroupOperation,T,𝔾}
    G = base_lie_group(𝔤)
    M = base_manifold(G)
    return get_vector!(M, X, identity_element(G, T), c, DefaultOrthonormalBasis(𝔽))
end

_doc_log_complex_circ = """
    log(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g)
    log!(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, X, g)

Compute the Lie group logarithm on the complex [`CircleGroup`](@ref), which coincides with the
ordinary complex logarithm.
"""

identity_element(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}) = 1.0 + 0.0im
identity_element(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, T::Union{<:Number,Type{<:Number}}) = one(T)
identity_element(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, ::Type{<:SArray{S,T}}) where {S,T} = @SArray fill(one(T))
identity_element(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, ::Type{<:MArray{S,T}}) where {S,T} = @MArray fill(one(T))

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g)

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log!(M::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, X, g)

function ManifoldsBase.log(::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}}, g::Number)
    return log(g)
end

function Base.show(io::IO, ::CircleGroup{ℂ, ScalarMultiplicationGroupOperation, Circle{ℂ}})
    return print(io, "CircleGroup()")
end