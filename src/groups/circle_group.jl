function CircleGroup end
"""
    CircleGroup

The circle group ``𝕊^1`` the circle together with composing points on the circle
by either adding angles.
The circle itself is a one dimensional Riemannian manifold.
Hence the Lie algebra is the real line.

The elements of the circle group can be represented in three different ways.


## As complex numbers

Representing the circle as complex numbers of absolute value one, that is

```math
𝕊¹ = $(_tex(:SetDef, "z ∈ ℂ", "|z| = 1", "big")) = $(_tex(:SetDef, "a + b$(_math(:i)) ∈ ℂ", "a^2+b^2 = 1", "big")),
```

where ``$(_math(:i))`` denotes the imaginary unit.
It is equipped with the group operation of complex
multiplication [`AbelianMultiplicationGroupOperation`](@ref).
That operation is given by

```math
(a + b$(_math(:i))) ∘ (c + d$(_math(:i))) := (ac - bd) + (ad + bc)*$(_math(:i)),
```

for complex numbers ``a + b$(_math(:i)), c + d$(_math(:i)) ∈ ℂ``.

## As part of the real line

Elements of the circle group can be represented by the angle
on the unit circle that they correspond to. In that case
the elements are represented by real numbers ``x ∈ [-π,π)`` and the
circle group is identified with a quotient space of the real numbers

```math
 𝕊¹ = ℝ / 2πℤ = $(_tex(:SetDef, "[x] ∈ ℝ / 2πℤ", "x ∈ [-π,π)", "big")).
```

It is equipped with the group operation of adding angles
``$(_tex(:rm, raw"mod\, ")) 2π`` via [`AdditionGroupOperation`](@ref).

## As part of the 2D plane ``ℝ^2``

Elements of the circle group can be represented as two dimensional
real valued vectors ``x ∈ ℝ`` of length 1.
In that case the circle group is identified with the unit circle in ``ℝ^2``,
that is the one dimensional [`Sphere`](@extref `Manifolds.Sphere`).

```math
𝕊^1 = $(_tex(:SetDef, "(x, y) ∈ ℝ^2", "x^2 + y^2 = 1", "big")).
```

It is equipped with the group operation of adding the angles
of two points on the unit circle which corresponds to the complex
multiplication

```math
(x_1, y_1) ∘ (x_2, y_2) := (x_1*x_2 - y_1*y_2, (x_1*y_2 + x_2*y_1)),
```
for real valued vectors ``(x_1, y_1)^$(_tex(:transp)), (x_2, y_2)^$(_tex(:transp)) ∈ ℝ^2`` via [`AbelianMultiplicationGroupOperation`](@ref).

# Constructors

    CircleGroup(Circle(ℂ))
    CircleGroup(ℂ)
    CircleGroup()

Generate the circle group represented as complex numbers.

    CircleGroup(Circle(ℝ))
    CircleGroup(ℝ)

Generate the circle group represented as real valued angles ``x ∈ [-π, π)``.

    CircleGroup(Sphere(1))
    CircleGroup(ℝ^2)

Generate the circle group represented as two dimensional real valued vectors of unit norm.

The default representation is by complex numbers and can be constructed with `CircleGroup()`.
"""
CircleGroup(::Union{ManifoldsBase.AbstractNumbers,Circle,Euclidean,Sphere})

#
#
# Common definitions
