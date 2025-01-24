"""
    CircleGroup = LieGroup{ℂ, ScalarMultiplicationGroupOperation, Manifolds.Circle{ℂ}}


The complex circle ``𝕊^1`` is the manifold represented by all complex-valued points ``z ∈ ℂ`` of norm ``1``:

```math
𝕊¹ := $(_tex(:SetDef, "z ∈ ℂ", "|z| = 1", "big")).
```

The standard complex multiplication (internally via [`ScalarMultiplicationGroupOperation`](@ref)) makes it a Lie group.

# Constructor

    CircleGroup()

Generate the complex circle group.
"""
const CircleGroup = LieGroup{
    ℂ, ScalarMultiplicationGroupOperation, Manifolds.Circle{ℂ}
}


"""
    RealCircleGroup = LieGroup{ℝ, AdditionGroupOperation, Manifolds.Circle{ℝ}}

The real circle ``𝕊^1`` is the manifold represented by all real-valued points  ``x ∈ [-π,π)`` and can therefore be understood as a symmetric system of representatives of ``ℝ$(_tex(:rm, raw"\, mod\, ")) 2πℤ``.
```math
 𝕊¹ :=  [-π,π) = ℝ $(_tex(:rm, raw"\, mod\, ")) 2πℤ

```


Addition ``$(_tex(:rm, raw"mod\, ")) 2π`` via [`AdditionGroupOperation`](@ref) defines its structure as a Lie group.

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
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples)

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