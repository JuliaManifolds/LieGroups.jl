"""
    CircleGroup{𝔽, M}

The circle group ``𝕊^1`` is the multiplicative group of complex numbers 
``z ∈ ℂ`` of absolute value ``1``. 
It is  a one dimensional Riemannian manifold and a Lie group. The Lie algebra is precisely the imaginary axis of the complex plane.
The elements of the circle group can be represented in three different ways.


The first way is to represent the elements of the circle group as complex numbers

```math
𝕊¹ = $(_tex(:SetDef, "z ∈ ℂ", "|z| = 1", "big")) = $(_tex(:SetDef, "a + bi ∈ ℂ", "a^2+b^2 = 1", "big")).
```

It is equipped with the group operation of complex 
multiplication [`ScalarMultiplicationGroupOperation`](@ref). 
That operation is given by

```math
(a + b*im) ∘ (c + d*im) := (ac - bd) + (ad + bc)*im,
```
for complex numbers ``(a + b*im), (c + d*im) ∈ ℂ``.


The second way to represent elements of the circle group is by the angle 
on the unit circle that they correspond to. In that case
the elements are represented by real numbers ``x ∈ [-π,π)`` and the 
circle group is identified with a quotient space of the real numbers

```math
 𝕊¹ = ℝ / 2πℤ = $(_tex(:SetDef, "[x] ∈ ℝ / 2πℤ", "x ∈ [-π,π)", "big")).
```

It is equipped with the group operation of adding angles 
``$(_tex(:rm, raw"mod\, ")) 2π`` via [`AdditionGroupOperation`](@ref).


The third way is to represent elements of the circle group as two dimensional 
real valued vectors. In that case the circle group
is identified with the unit circle in ``ℝ^2``, i.e. the 
one dimensional [`Sphere`](@extref `Manifolds.Sphere`).

```math
𝕊^1 = $(_tex(:SetDef, "(x, y) ∈ ℝ^2", "x^2 + y^2 = 1", "big")).
```

It is equipped with the group operation of adding the angles 
of two points on the unit circle which corresponds to the complex 
multiplication

```math
(x_1, y_1) ∘ (x_2, y_2) := ((x_1*x_2 - y_1*y_2), (x_1*y_2 + x_2*y_1)),
```
for real valued vectors ``(x_1, y_1), (x_2, y_2) ∈ ℂ`` via [`MultiplicationGroupOperation`](@ref).

# Constructor
    	
    CircleGroup(Circle(ℂ))
    CircleGroup(ℂ)
    CircleGroup()

Generate the circle group represented as complex numbers.

    CircleGroup(Circle(ℝ))
    CircleGroup(ℝ)

Generate the circle group represented as real valued angles
``x ∈ [-π, π)``.

    CircleGroup(Sphere(1))
    CircleGroup(ℝ^2)

Generate the circle group represented as two dimensional real valued vectors.

The default representation is by complex numbers and can be constructed with `CircleGroup()`.
"""
const CircleGroup{𝔽, Op, M <: AbstractManifold{𝔽}} = LieGroup{𝔽, Op, M}

#functions for different representations in seperate files