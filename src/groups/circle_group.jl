@doc raw"""
    CircleGroup <: GroupManifold{Circle{ℂ},MultiplicationOperation}


The circle ``𝕊^1`` is a manifold here represented by
real-valued points in ``[-π,π)`` or complex-valued points ``z ∈ ℂ`` of absolute value
``\lvert z\rvert = 1``.
The real circle becomes a Lie-Group by addition of the representing angles modulo π using [`AdditionGroupOperation`](@ref).     
The complex circle is equipped with the group operation of standard complex multiplication ([`ScalarMultiplicationGroupOperation`](@ref)).
"""

const CircleGroup = LieGroup{
    ℂ, ScalarMultiplicationGroupOperation, Manifolds.Circle{ℂ}
}
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
lol

"""

@doc "$(_doc_exp_real_circ)"
exp(::RealCircleGroup, ::Identity{AdditionGroupOperation}, X)

@doc "$(_doc_exp_real_circ)"
exp!(M::RealCircleGroup, g, Id::Identity{AdditionGroupOperation}, X)


_doc_exp_complex_circ = """s
    exp(::CircleGroup{ℂ}, ::Identity{ScalarMultiplicationGroupOperation}, X)
    exp!(::CircleGroup{ℂ}, g, ::Identity{ScalarMultiplicationGroupOperation}, X)

Compute the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples)


```mathS1
$(_tex(:exp)) it = 
```

"""


@doc "$(_doc_exp_complex_circ)"
exp(::CircleGroup, ::Identity{ScalarMultiplicationGroupOperation}, X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup, g, Id::Identity{ScalarMultiplicationGroupOperation}, X)



function Base.show(io::IO, G::CircleGroup)
    return print(io, "CircleGroup()")
end