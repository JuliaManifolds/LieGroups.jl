"""
    AbelianMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by an abelian multiplication.
"""
struct AbelianMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

function _compose(
        ::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, g::Number, h::Number
    ) where {𝔽}
    return g * h
end

function _compose(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        p::AbstractArray{<:Any, 0},
        q::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map((pp, qq) -> compose(G, pp, qq), p, q)
end

function _compose(
        ::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        g::Number,
        h::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return g .* h
end

function _compose(
        ::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h::Number,
    ) where {𝔽}
    return g .* h
end

function _compose!(G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, k, g, h) where {𝔽}
    return copyto!(k, compose(G, g, h))
end

function conjugate(::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, g, h) where {𝔽}
    return g
end

function conjugate!(::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, k, g, h) where {𝔽}
    return copyto!(k, g)
end

_doc_diff_conjugate_abelmul = """
    diff_conjugate(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the conjugate ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = ghg^{-1}``,
which simplifies for an [`AbelianMultiplicationGroupOperation`](@ref) to ``$(_math(:d))(c_g(h))[X] = X``.

This can be computed in-place of `Y`.
"""

@doc "$(_doc_diff_conjugate_abelmul)"
function diff_conjugate(
        ::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, g, h, X::Number
    ) where {𝔽}
    return X
end

function diff_conjugate(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h::AbstractArray{<:Any, 0},
        X::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map(XX -> diff_conjugate(G, g, h, XX), X)
end

@doc "$(_doc_diff_conjugate_abelmul)"
function diff_conjugate!(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        Y::AbstractArray{<:Any, 0},
        g::AbstractArray{<:Any, 0},
        h::AbstractArray{<:Any, 0},
        X::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return copyto!(LieAlgebra(G), Y, diff_conjugate(G, g, h, X))
end

_doc_diff_inv_abelmult = """
    diff_inv(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, X)
    diff_inv!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, X)

Compute the value of the differential ``$(_math(:d))ι_{$(_math(:G))}(g)[X]`` of the inversion ``ι_{$(_math(:G))}(g) := g^{-1}`` at ``X ∈ 𝔤``
in the [`LieAlgebra`](@ref) ``𝔤`` of the [`LieGroup`](@ref) `G`.

In the Abelian case, the computation simplifies to

```math
$(_math(:d))ι_{$(_math(:G))}(g)[X] = -gXg^{-1} = -X.
```

This can be computed in-place of `Y` if `Y` is `mutable`.
"""

@doc "$(_doc_diff_inv_abelmult)"
function diff_inv!(G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, Y, g, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, -X)
end

_doc_diff_left_compose_abelmult = """
    diff_left_compose(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, h, X)
    diff_left_compose!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its first (left) argument `g`.

Another interpretation is to consider a function where we do a fixed multiplication from the right with `h`.
i..e. the right group multiplication function ``ρ_h(g) = g$(_math(:∘))h``.

The differential simplifies for an [`AbelianMultiplicationGroupOperation`](@ref) to the identity, i.e.
``$(_math(:d))ρ_h(g)[X] = X``.

This can be computed in-place of `Y` if `Y` is `mutable`.
""" #

@doc "$(_doc_diff_left_compose_abelmult)"
function diff_left_compose(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h,
        X::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map((gg, XX) -> diff_left_compose(G, gg, h, XX), g, X)
end

@doc "$(_doc_diff_left_compose_abelmult)"
function diff_left_compose!(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, Y, g, h, X
    ) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_diff_right_compose_abelmult = """
    diff_right_compose(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, h, X)
    diff_right_compose!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its second (right) argument `h`.

Another interpretation is to consider a function where we do a fixed multiplication from the left with `g`.
i..e. the left group multiplication function ``λ_g(h) = g$(_math(:∘))h``.

The differential simplifies for an [`AbelianMultiplicationGroupOperation`](@ref) to the identity, i.e.
``$(_math(:d))λ_g(h)[X] = X``.
"""

@doc "$(_doc_diff_right_compose_abelmult)"
function diff_right_compose(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h,
        X::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map((gg, XX) -> diff_right_compose(G, gg, h, XX), g, X)
end

@doc "$(_doc_diff_right_compose_abelmult)"
function diff_right_compose!(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, Y, g, h, X
    ) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_exp_abelmult = """
    exp(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, X)
    exp(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, X)
    exp!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, h, X)
    exp!(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, h, g, X)

Compute the Lie group exponential on a [`LieGroup`](@ref) at a point `g` or the [`Identity`](@ref) with an [`AbelianMultiplicationGroupOperation`](@ref).

Due to differences in the representation of some abelian Lie groups, this method wraps a concrete implementation of a specific abelian LieGroup with inputs of type `AbstractArray{<:Any,0}` and supports in-place computation.

This can be computed in-place of `h` if `h` is `mutable`.
"""

@doc "$(_doc_exp_abelmult)"
function ManifoldsBase.exp(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, X::AbstractArray{<:Any, 0}
    ) where {𝔽}
    return map(XX -> exp(G, XX), X)
end

function ManifoldsBase.exp(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, g, X
    ) where {𝔽}
    return compose(G, g, exp(G, X))
end

function ManifoldsBase.exp(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        X::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map((gg, XX) -> gg * exp(G, XX), g, X)
end

@doc "$(_doc_exp_abelmult)"
function ManifoldsBase.exp!(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, h, X
    ) where {𝔽}
    return copyto!(h, exp(G, X))
end

function ManifoldsBase.exp!(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, h, g, X
    ) where {𝔽}
    return copyto!(h, exp(G, g, X))
end

_doc_inv_abelmult = """
    inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g)
    inv!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for the [`AbelianMultiplicationGroupOperation`](@ref) simplifies for a scalar input to the ordinary scalar inverse ``g^{-1}``.

Due to differences in the representation of some abelian Lie groups, this method wraps a concrete implementation of a specific abelian LieGroup with inputs of type `AbstractArray{<:Any,0}` and supports in-place computation.

This can be computed in-place of `h` if `h` is `mutable`.
"""

@doc "$(_doc_inv_abelmult)"
Base.inv(::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, g::Number) where {𝔽} = inv(g)

function Base.inv(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, g::AbstractArray{<:Any, 0}
    ) where {𝔽}
    return map(gg -> inv(G, gg), g)
end

@doc "$(_doc_inv_abelmult)"
inv!(::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, h, g) where {𝔽}

function _inv!(G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, h, g) where {𝔽}
    h .= inv.(g)
    return h
end

function inv_left_compose(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, g::Number, h::Number
    ) where {𝔽}
    return inv(g) * h
end

function inv_left_compose(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map(\, g, h)
end

function inv_left_compose!(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, k, g, h
    ) where {𝔽}
    return copyto!(k, inv_left_compose(G, g, h))
end

function inv_right_compose(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, g::Number, h::Number
    ) where {𝔽}
    return g * inv(h)
end

function inv_right_compose(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        g::AbstractArray{<:Any, 0},
        h::AbstractArray{<:Any, 0},
    ) where {𝔽}
    return map(/, g, h)
end

function inv_right_compose!(
        G::LieGroup{𝔽, <:AbelianMultiplicationGroupOperation}, k, g, h
    ) where {𝔽}
    return copyto!(k, inv_right_compose(G, g, h))
end

_doc_identity_element_abelian_mult = """
    identity_element(G::LieGroup{𝔽,AbelianMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,AbelianMultiplicationGroupOperation}, e)

Return the point representation of the [`Identity`](@ref),
which for an [`AbelianMultiplicationGroupOperation`](@ref) is the one-element.

This can be computed in `e` if `e` is `mutable`.
"""

@doc "$(_doc_identity_element_abelian_mult)"
function identity_element(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, ::Type{T}
    ) where {𝔽, T <: Union{Number, AbstractArray{0, <:Number}}}
    return one(T)
end
function identity_element(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, ::Type{Array{T, 0}}
    ) where {𝔽, T <: Number}
    return fill(one(T))
end

@doc "$(_doc_identity_element_abelian_mult)"
identity_element!(::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, e) where {𝔽}

function identity_element!(
        ::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, e::AbstractArray{<:Number, 0}
    ) where {𝔽}
    return fill!(e, 1)
end

_doc_lie_bracket_abelmult = """
    lie_bracket(::LieAlgebra{𝔽,AbelianMultiplicationGroupOperation}, X, Y)
    lie_bracket!(::LieAlgebra{𝔽,AbelianMultiplicationGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`AbelianMultiplicationGroupOperation`](@ref) yields the zero vector of the [`LieAlgebra`](@ref) due to commutativity.

```math
[X, Y] = XY-YX = 0
```

The computation can be done in-place of `Z` if `Z` is `mutable`.
"""

@doc "$(_doc_lie_bracket_abelmult)"
function lie_bracket(
        ::LieAlgebra{𝔽, O, <:LieGroup{𝔾, O}}, X::Number, Y::Number
    ) where {𝔽, 𝔾, O <: AbelianMultiplicationGroupOperation}
    return zero(X)
end

@doc "$(_doc_lie_bracket_abelmult)"
function lie_bracket!(
        ::LieAlgebra{𝔽, O, <:LieGroup{𝔾, O}}, Z, X, Y
    ) where {𝔽, 𝔾, O <: AbelianMultiplicationGroupOperation}
    return fill!(Z, 0)
end

_doc_log_abelmult = """
    log(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, h)
    log(G::LieGroup{𝔽,<:AbelianMultiplicationGroupOperation}, g, h)
    log!(G::LieGroup{𝔽,<:AbeliantMultiplicationGroupOperation}, X, g)
    log!(G::LieGroup{𝔽,<:AbeliantMultiplicationGroupOperation}, X, g, h)

Compute the Lie group logarithm on a [`LieGroup`](@ref) at a point `g` or the [`Identity`](@ref) with a concrete instance of [`AbelianMultiplicationGroupOperation`](@ref).

Due to differences in the representation of some abelian Lie groups, this method wraps a concrete implementation of a specific abelian LieGroup with inputs of type `AbstractArray{<:Any,0}` and supports in-place computation.

This can be computed in-place of `X` if `X` is `mutable`.
"""

@doc "$(_doc_log_abelmult)"
function ManifoldsBase.log(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, g, h
    ) where {𝔽}
    return log(G, compose(G, inv(G, g), h))
end

function ManifoldsBase.log(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        ::Identity{AbelianMultiplicationGroupOperation},
    ) where {𝔽}
    return zero_vector(LieAlgebra(G))
end

function ManifoldsBase.log(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        ::Identity{AbelianMultiplicationGroupOperation},
        T::Type,
    ) where {𝔽}
    return zero_vector(LieAlgebra(G), T)
end

@doc "$(_doc_log_abelmult)"
function ManifoldsBase.log!(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation}, X, g
    ) where {𝔽}
    return copyto!(X, log(G, g))
end

function ManifoldsBase.log!(
        G::LieGroup{𝔽, AbelianMultiplicationGroupOperation},
        X,
        ::Identity{AbelianMultiplicationGroupOperation},
    ) where {𝔽}
    return zero_vector!(LieAlgebra(G), X)
end
