#
#circle group represented by complex numbers, operation: complex multiplication
#
function CircleGroup(M::Manifolds.Circle{ℂ})
    return CircleGroup{ℂ, typeof(M)}(
        M, ScalarMultiplicationGroupOperation()
    )
end 

CircleGroup(𝔽::ManifoldsBase.AbstractNumbers=ℂ) = CircleGroup(Circle(𝔽))

function diff_left_compose(::CircleGroup{ℂ, Circle{ℂ}}, g::Number, h::Any, X::Number)
    return g * X
end

function diff_right_compose(::CircleGroup{ℂ, Circle{ℂ}}, g::Number, h::Any, X::Number)
    return X * g
end

_doc_exp_complex_circ = """
    exp(::CircleGroup{ℂ, Circle{ℂ}}, X)
    exp!(::CircleGroup{ℂ, Circle{ℂ}}, g, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```
"""

@doc "$(_doc_exp_complex_circ)"
Base.exp(::CircleGroup{ℂ, Circle{ℂ}}, X::Number) = exp(X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup{ℂ, Circle{ℂ}}, g, X)

_doc_log_complex_circ = """
    log(::CircleGroup{ℂ, Circle{ℂ}}, g)
    log!(::CircleGroup{ℂ, Circle{ℂ}}, X, g)

Compute the Lie group logarithm on the complex [`CircleGroup`](@ref), which coincides with the
ordinary complex logarithm.
"""

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log(::CircleGroup{ℂ, Circle{ℂ}}, g)

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log!(M::CircleGroup{ℂ, Circle{ℂ}}, X, g)

function ManifoldsBase.log(::CircleGroup{ℂ, Circle{ℂ}}, g::Number)
    return log(g)
end

function Base.show(io::IO, ::CircleGroup{ℂ, Circle{ℂ}})
    return print(io, "CircleGroup()")
end