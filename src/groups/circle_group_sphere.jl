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

function diff_conjugate(::_PlanarCircleGroup, g, h, X)
    return X
end

_doc_diff_left_compose_mult_planar_circ = """
    diff_left_compose(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, h, X)
    diff_left_compose!(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, Y, g, h, X)

Compute the differential of the left group multiplication ``λ_g(h) = g$(_math(:∘))h`` of the Circle Group,
represented as two dimensional vectors in ``ℝ^2``.

It simplifies for the [`AbelianMultiplicationGroupOperation`](@ref) to ``Dλ_g(h)[X] = gX``,
where the multiplication corresponds to the complex multiplication after canonical identification
of the real plane with the complex plane.
"""

@doc "$(_doc_diff_left_compose_mult_planar_circ)"
diff_left_compose(C::_PlanarCircleGroup, g, h, X) = compose(C, g, X)

_doc_diff_right_compose_mult_planar_circ = """
    diff_right_compose(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, h, X)
    diff_right_compose!(G::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, Y, g, h, X)

Compute the differential of the right group multiplication ``λ_g(h) = g$(_math(:∘))h`` of the Circle Group,
represented as two dimensional vectors in ``ℝ^2``.

It simplifies for the [`AbelianMultiplicationGroupOperation`](@ref) to ``Dλ_g(h)[X] = Xg``,
where the multiplication corresponds to the complex multiplication after canonical identification
of the real plane with the complex plane.
"""

@doc "$(_doc_diff_right_compose_mult_planar_circ)"
diff_right_compose(C::_PlanarCircleGroup, g, h, X) = compose(C, X, g)

_doc_exp_planar_circ = """
    exp(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, X)
    exp!(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g, X)

Computes the Lie group exponential on the [`CircleGroup`](@ref), represented as two dimensional vectors in the real plane.
It coincides with the [ordinary complex exponential](https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)#Examples) after
canonical identification of the real plane with the complex plane.

This can be computed in-place of `g`.
```math
$(_tex(:exp)) (t) = $(_tex(:pmatrix, _tex(:cos)*"(t)", _tex(:sin)*"(t)"))
```
"""
@doc "$(_doc_exp_planar_circ)"
function Base.exp(::_PlanarCircleGroup, X)
    z = exp(X[1] + X[2] * im)
    return [z.re, z.im]
end

@doc "$(_doc_exp_planar_circ)"
ManifoldsBase.exp!(G::_PlanarCircleGroup, g, X) = copyto!(g, exp(G, X))

function _compose(::_PlanarCircleGroup, p, q)
    a = p[1]
    b = p[2]
    c = q[1]
    d = q[2]
    z = [(a * c - b * d), (a * d + b * c)]
    return z
end

identity_element(::_PlanarCircleGroup) = [1.0, 0.0]
function identity_element!(::_PlanarCircleGroup, p::Array{<:Any,1})
    p = [1.0, 0.0]
    return p
end

Base.inv(::_PlanarCircleGroup, p) = [p[1], -p[2]]
Base.inv(::_PlanarCircleGroup, e::Identity{<:AbelianMultiplicationGroupOperation}) = e

function inv_left_compose(C::_PlanarCircleGroup, g, h)
    g1 = Base.inv(C, g)
    return compose(C, g1, h)
end

function inv_right_compose(C::_PlanarCircleGroup, g, h)
    h1 = Base.inv(C, h)
    return compose(C, g, h1)
end

_doc_log_planar_circ = """
    log(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, g)
    log!(::LieGroup{ℝ, AbelianMultiplicationGroupOperation, Sphere}, X, g)

Compute the Lie group logarithm on the [`CircleGroup`](@ref), , represented as two dimensional vectors in the real plane.
It coincides with the ordinary complex logarithm after canonical identification of the real plane with the complex plane.

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_planar_circ)"
ManifoldsBase.log(::_PlanarCircleGroup, g)

@doc "$(_doc_log_planar_circ)"
ManifoldsBase.log!(M::_PlanarCircleGroup, X, g)

function ManifoldsBase.log(::_PlanarCircleGroup, g)
    z = log(g[1] + g[2]*im)
    return [z.re, z.im]
end

function ManifoldsBase.log(C::_PlanarCircleGroup, ::Identity{AbelianMultiplicationGroupOperation})
    return zero_vector(LieAlgebra(C))
end

function Base.show(io::IO, ::_PlanarCircleGroup)
    return print(io, "CircleGroup(Sphere(1))")
end
