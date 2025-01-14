@doc raw"""
    CircleGroup <: GroupManifold{Circle{ℂ},MultiplicationOperation}

The circle group is the complex circle ([`Circle(ℂ)`](@ref)) equipped with
the group operation of complex multiplication ([`ScalarMultiplicationGroupOperation`](@ref)).
"""

const CircleGroup = LieGroup{
   ℂ, ScalarMultiplicationGroupOperation, Manifolds.AbstractSphere{ℂ}
}

function CircleGroup()
    circ = Manifolds.Sphere(0, ℂ)
    return CircleGroup(
        circ, ScalarMultiplicationGroupOperation()
        )
end

_doc_exp_circ = """
    exp(::CircleGroup, ::Identity{ScalarMultiplicationGroupOperation}, X)
    exp!(::CircleGroup, g, ::Identity{ScalarMultiplicationGroupOperation}, X)

Compute the Lie group exponential on the [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples)


```math
$(_tex(:exp)) it = 
```

"""





@doc "$(_doc_exp_circ)"
exp(::CircleGroup, ::Identity{ScalarMultiplicationGroupOperation}, X)

@doc "$(_doc_exp_circ)"
exp!(M::CircleGroup, g, Id::Identity{ScalarMultiplicationGroupOperation}, X)



function Base.show(io::IO, G::CircleGroup)
    n = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "CircleGroup()")
end