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
const CircleGroup = LieGroup{
    ℂ, ScalarMultiplicationGroupOperation, Manifolds.Circle{ℂ}
}


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
const RealCircleGroup = LieGroup{
    ℝ, AdditionGroupOperation, Manifolds.Circle{ℝ}
}



function CircleGroup()
    circ = Manifolds.Circle(ℂ)
    return CircleGroup(circ, ScalarMultiplicationGroupOperation())
end

function RealCircleGroup()
    circ = Manifolds.Circle(ℝ)
    return RealCircleGroup(circ, AdditionGroupOperation())
end

_doc_exp_real_circ = """
    exp(::RealCircleGroup, ::Identity{AdditionGroupOperation}, X)
    exp!(::RealCircleGroup, g, ::Identity{AdditionGroupOperation}, X)

The Lie group exponential on the [`RealCircleGroup`](@ref) is given by the projection into the equivalence class of its defining relation. 

This can be computed in-place of `X`.
"""

@doc "$(_doc_exp_real_circ)"
exp(::RealCircleGroup, ::Identity{AdditionGroupOperation}, X)

@doc "$(_doc_exp_real_circ)"
exp!(M::RealCircleGroup, g, Id::Identity{AdditionGroupOperation}, X)


_doc_exp_complex_circ = """
    exp(::CircleGroup, ::Identity{ScalarMultiplicationGroupOperation}, X)
    exp!(::CircleGroup, g, ::Identity{ScalarMultiplicationGroupOperation}, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `X`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```

"""


@doc "$(_doc_exp_complex_circ)"
exp(::CircleGroup, ::Identity{ScalarMultiplicationGroupOperation}, X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup, g, Id::Identity{ScalarMultiplicationGroupOperation}, X)



function Base.show(io::IO, ::CircleGroup)
    return print(io, "CircleGroup()")
end
function Base.show(io::IO, ::RealCircleGroup)
    return print(io, "RealCircleGroup()")
end