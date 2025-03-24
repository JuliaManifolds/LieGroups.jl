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
const CircleGroup = LieGroup{ℂ,ScalarMultiplicationGroupOperation,Manifolds.Circle{ℂ}}

function CircleGroup()
    circ = Manifolds.Circle(ℂ)
    return CircleGroup(circ, ScalarMultiplicationGroupOperation())
end

function diff_left_compose(::CircleGroup, g::Number, h::Any, X::Number)
    return g * X
end

function diff_right_compose(::CircleGroup, g::Number, h::Any, X::Number)
    return X * g
end

_doc_exp_complex_circ = """
    exp(::CircleGroup, X)
    exp!(::CircleGroup, g, X)

Computes the Lie group exponential on the complex [`CircleGroup`](@ref), which coincides with the
[ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples).

The Lie algebra is precisely the imaginary axis of the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) ($(_math(:i))t) = $(_tex(:cos))(t) + $(_math(:i))$(_tex(:sin))(t)
```
"""

@doc "$(_doc_exp_complex_circ)"
Base.exp(::CircleGroup, X::Number) = exp(X)

@doc "$(_doc_exp_complex_circ)"
exp!(M::CircleGroup, g, X)

_doc_log_complex_circ = """
    log(::CircleGroup, g)
    log!(::CircleGroup, X, g)

Compute the Lie group logarithm on the complex [`CircleGroup`](@ref), which coincides with the
ordinary complex logarithm.
"""

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log(::CircleGroup, g)

@doc "$(_doc_log_complex_circ)"
ManifoldsBase.log!(M::CircleGroup, X, g)

function ManifoldsBase.log(::CircleGroup, g::Number)
    return log(g)
end

function Base.show(io::IO, ::CircleGroup)
    return print(io, "CircleGroup()")
end

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
const RealCircleGroup = LieGroup{ℝ,AdditionGroupOperation,Manifolds.Circle{ℝ}}

function RealCircleGroup()
    circ = Manifolds.Circle(ℝ)
    return RealCircleGroup(circ, AdditionGroupOperation())
end

_doc_exp_real_circ = """
    exp(::RealCircleGroup, X)
    exp!(::RealCircleGroup, g, X)

The Lie group exponential on the [`RealCircleGroup`](@ref) is given by the projection into the equivalence class of its defining relation.

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

compose(::RealCircleGroup, p, q) = sym_rem(p + q)
compose(::RealCircleGroup, ::Identity{AdditionGroupOperation}, q) = sym_rem(q)
compose(::RealCircleGroup, p, ::Identity{AdditionGroupOperation}) = sym_rem(p)
function compose(
    ::RealCircleGroup,
    e::Identity{AdditionGroupOperation},
    ::Identity{AdditionGroupOperation},
)
    return e
end

compose!(::RealCircleGroup, x, p, q) = copyto!(x, sym_rem(p + q))
function compose!(::RealCircleGroup, x, ::Identity{AdditionGroupOperation}, q)
    return copyto!(x, sym_rem(q))
end
function compose!(::RealCircleGroup, x, p, ::Identity{AdditionGroupOperation})
    return copyto!(x, sym_rem(p))
end
function compose!(
    ::RealCircleGroup,
    ::Identity{AdditionGroupOperation},
    e::Identity{AdditionGroupOperation},
    ::Identity{AdditionGroupOperation},
)
    return e
end

identity_element(::RealCircleGroup) = 0.0
identity_element(::RealCircleGroup, p::Union{<:Number,Type{<:Number}}) = zero(p)

Base.inv(G::RealCircleGroup, p::Number) = sym_rem(-p)

Base.inv(G::RealCircleGroup, p::AbstractArray{<:Any,0}) = map(pp -> inv(G, pp), p)

@doc "$(_doc_exp_real_circ)"
exp(::RealCircleGroup, X)

@doc "$(_doc_exp_real_circ)"
exp!(M::RealCircleGroup, g, X)

function Base.show(io::IO, ::RealCircleGroup)
    return print(io, "RealCircleGroup()")
end
