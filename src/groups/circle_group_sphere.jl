#
#circle group represented as vectors in ℝ^2, operation: adding angles in the unit circle
#
function CircleGroup(M::Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ})
    return LieGroup{ℝ,AbelianMultiplicationGroupOperation,typeof(M)}(
        M, AbelianMultiplicationGroupOperation()
    )
end

function CircleGroup(::Euclidean{ManifoldsBase.TypeParameter{Tuple{2}},ℝ})
    return CircleGroup(Sphere(1))
end

const _PlanarCircleGroup = LieGroup{
    ℝ,AbelianMultiplicationGroupOperation,<:Sphere{ManifoldsBase.TypeParameter{Tuple{1}},ℝ}
}

#=
function diff_conjugate!(::_PlanarCircleGroup, g, h, X)
    return X
end
=#
_doc_diff_left_compose_mult_planar_circ = """
    diff_left_compose(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, h, X)
    diff_left_compose!(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h`` of the Circle Group,
represented as two dimensional vectors in ``ℝ^2``.

It simplifies for the [`AbelianMultiplicationGroupOperation`](@ref) to ``Dλ_g(h)[X] = g*X``,
where the multiplication corresponds to the complex multiplication after canonical identification
of the real plane with the complex plane.

This can be computed in-place of `Y`.
"""

@doc "$(_doc_diff_left_compose_mult_planar_circ)"
diff_left_compose!(C::_PlanarCircleGroup, Y, g, h, X) = compose!(C, Y, g, X)

@doc "$(_doc_diff_left_compose_mult_planar_circ)"
diff_left_compose(_PlanarCircleGroup, g, h, X)

_doc_diff_right_compose_mult_planar_circ = """
    diff_right_compose(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, h, X)
    diff_right_compose!(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, Y, g, h, X)

Compute the differential of the right group multiplication ``λ_g(h) = g$(_math(:∘))h`` of the Circle Group,
represented as two dimensional vectors in ``ℝ^2``.

It simplifies for the [`AbelianMultiplicationGroupOperation`](@ref) to ``Dλ_g(h)[X] = X*g``,
where the multiplication corresponds to the complex multiplication after canonical identification
of the real plane with the complex plane.
"""

@doc "$(_doc_diff_right_compose_mult_planar_circ)"
diff_right_compose!(C::_PlanarCircleGroup, Y, g, h, X) = compose!(C, Y, X, g)

_doc_exp_planar_circ = """
    exp(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, X)
    exp!(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, X)

Compute the Lie group exponential on the [`CircleGroup`](@ref), represented as two dimensional vectors in the real plane.
It coincides with the [ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples) after
canonical identification of the real plane with the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) $(_tex(:pmatrix, 0, "t")) = $(_tex(:pmatrix, _tex(:cos)*"(t)", _tex(:sin)*"(t)"))
```
"""
@doc "$(_doc_exp_planar_circ)"
ManifoldsBase.exp(::_PlanarCircleGroup, X)

@doc "$(_doc_exp_planar_circ)"
function ManifoldsBase.exp!(::_PlanarCircleGroup, g, X)
    z = exp(X[1] + X[2] * im)
    g[1] = real(z)
    g[2] = imag(z)
    return g
end

function _compose!(::_PlanarCircleGroup, k, g, h)
    a = g[1]
    b = g[2]
    c = h[1]
    d = h[2]
    k[1] = a * c - b * d
    k[2] = a * d + b * c
    return k
end

identity_element(::_PlanarCircleGroup, ::Type{<:AbstractVector}) = [1.0, 0.0]
function identity_element(::_PlanarCircleGroup, ::Type{<:AbstractVector{T}}) where {T}
    return [one(T), zero(T)]
end
function identity_element!(::_PlanarCircleGroup, p::V) where {T,V<:AbstractVector{T}}
    p .= [one(T), zero(T)]
    return p
end

function inv!(::_PlanarCircleGroup, q, p)
    q[1] = p[1]
    q[2] = -p[2]
    return q
end

function inv!(
    G::PG, q, ::Identity{AbelianMultiplicationGroupOperation}
) where {PG<:_PlanarCircleGroup}
    return identity_element!(G, q)
end

function inv_left_compose!(C::_PlanarCircleGroup, k, g, h)
    inv!(C, k, g)
    return compose!(C, k, k, h)
end

function inv_right_compose!(C::_PlanarCircleGroup, k, g, h)
    inv!(C, k, h)
    return compose!(C, k, g, k)
end

_doc_log_planar_circ = """
    log(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g)
    log!(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, X, g)

Compute the Lie group logarithm on the [`CircleGroup`](@ref), represented as two dimensional vectors in the real plane.
It coincides with the ordinary complex logarithm after canonical identification of the real plane with the complex plane.

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_planar_circ)"
ManifoldsBase.log(::_PlanarCircleGroup, g)

@doc "$(_doc_log_planar_circ)"
ManifoldsBase.log!(M::_PlanarCircleGroup, X, g)

function ManifoldsBase.log(::_PlanarCircleGroup, g)
    z = log(g[1] + g[2] * im)
    return [z.re, z.im]
end

function ManifoldsBase.log(
    C::_PlanarCircleGroup, ::Identity{AbelianMultiplicationGroupOperation}
)
    return zero_vector(LieAlgebra(C))
end

function Base.show(io::IO, ::_PlanarCircleGroup)
    return print(io, "CircleGroup(Sphere(1))")
end
