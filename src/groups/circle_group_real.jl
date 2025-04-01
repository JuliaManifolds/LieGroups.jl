#
#circle group represented in ℝ mod 2π = [-π, π), operation: addition mod 2π
#
function CircleGroup(M::Manifolds.Circle{ℝ})  
    return CircleGroup{ℝ, typeof(M)}(
        M, AdditionGroupOperation()
    )
end

_doc_exp_real_circ = """
    exp(::CircleGroup{ℝ, Circle{ℝ}}, X)
    exp!(::CircleGroup{ℝ, Circle{ℝ}}, g, X)

The Lie group exponential on the [`CircleGroup`](@ref)
represented as angles in [-π,π)
is given by the projection into the equivalence class of its defining relation.

This can be computed in-place of `X`.
"""

@doc raw"""
    sym_rem(x,[T=π])

Compute symmetric remainder of `x` with respect to the interall 2*`T`, i.e.
`(x+T)%2T`, where the default for `T` is ``π``
"""
function sym_rem(x::N, T=π) where {N<:Number}
    return (x ≈ T ? convert(N, -T) : rem(x, convert(N, 2 * T), RoundNearest))
end
sym_rem(x, T=π) = map(sym_rem, x, Ref(T))

compose(::CircleGroup{ℝ, Circle{ℝ}}, p, q) = sym_rem(p + q)
compose(::CircleGroup{ℝ, Circle{ℝ}}, ::Identity{AdditionGroupOperation}, q) = sym_rem(q)
compose(::CircleGroup{ℝ, Circle{ℝ}}, p, ::Identity{AdditionGroupOperation}) = sym_rem(p)
function compose(
    ::CircleGroup{ℝ, Circle{ℝ}},
    e::Identity{AdditionGroupOperation},
    ::Identity{AdditionGroupOperation},
)
    return e
end

compose!(::CircleGroup{ℝ, Circle{ℝ}}, x, p, q) = copyto!(x, sym_rem(p + q))
function compose!(::CircleGroup{ℝ, Circle{ℝ}}, x, ::Identity{AdditionGroupOperation}, q)
    return copyto!(x, sym_rem(q))
end
function compose!(::CircleGroup{ℝ, Circle{ℝ}}, x, p, ::Identity{AdditionGroupOperation})
    return copyto!(x, sym_rem(p))
end
function compose!(
    ::CircleGroup{ℝ, Circle{ℝ}},
    ::Identity{AdditionGroupOperation},
    e::Identity{AdditionGroupOperation},
    ::Identity{AdditionGroupOperation},
)
    return e
end

identity_element(::CircleGroup{ℝ, Circle{ℝ}}) = 0.0
identity_element(::CircleGroup{ℝ, Circle{ℝ}}, p::Union{<:Number,Type{<:Number}}) = zero(p)

Base.inv(G::CircleGroup{ℝ, Circle{ℝ}}, p::Number) = sym_rem(-p)

Base.inv(G::CircleGroup{ℝ, Circle{ℝ}}, p::AbstractArray{<:Any,0}) = map(pp -> inv(G, pp), p)

@doc "$(_doc_exp_real_circ)"
exp(::CircleGroup{ℝ, Circle{ℝ}}, X)

@doc "$(_doc_exp_real_circ)"
exp!(M::CircleGroup{ℝ, Circle{ℝ}}, g, X)

function Base.show(io::IO, ::CircleGroup{ℝ, Circle{ℝ}})
    return print(io, "CircleGroup(ℝ)")
end