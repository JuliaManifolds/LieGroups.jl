"""
    AbstractMultiplicationGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `*` being overloaded, for example
`_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, a, b) = a * b`
"""
abstract type AbstractMultiplicationGroupOperation <: AbstractGroupOperation end

"""
    MatrixMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A group operation that is realised by a matrix multiplication.
"""
struct MatrixMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

Base.:*(::Identity{MatrixMultiplicationGroupOperation}, p::Union{AbstractMatrix, Number}) = p
function Base.:*(
        p::Union{AbstractMatrix, Number}, ::Identity{MatrixMultiplicationGroupOperation}
    )
    return p
end
function Base.:*(
        e::Identity{<:AbstractMultiplicationGroupOperation},
        ::Identity{<:AbstractMultiplicationGroupOperation},
    )
    return e
end
function Base.:*(
        ::Identity{<:AbstractMultiplicationGroupOperation}, e::Identity{AdditionGroupOperation}
    )
    return e
end
function Base.:*(
        e::Identity{AdditionGroupOperation}, ::Identity{<:AbstractMultiplicationGroupOperation}
    )
    return e
end

Base.:/(p, ::Identity{<:AbstractMultiplicationGroupOperation}) = p
Base.:/(::Identity{<:AbstractMultiplicationGroupOperation}, p) = inv(p)
function Base.:/(
        e::Identity{<:AbstractMultiplicationGroupOperation},
        ::Identity{<:AbstractMultiplicationGroupOperation},
    )
    return e
end

Base.:\(p, ::Identity{<:AbstractMultiplicationGroupOperation}) = inv(p)
Base.:\(::Identity{<:AbstractMultiplicationGroupOperation}, p) = p
function Base.:\(
        e::Identity{<:AbstractMultiplicationGroupOperation},
        ::Identity{<:AbstractMultiplicationGroupOperation},
    )
    return e
end

_doc_compose_mult = """
    compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h)

Compute the group operation composition of `g` and `h` on an [`LieGroup`](@ref) `G`.

For an [`AbstractMultiplicationGroupOperation`](@ref) this falls back to calling
`g*h`, where `*` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_mult)"
compose(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any, ::Any) where {𝔽}
@doc "$(_doc_compose_mult)"
compose!(
    ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any, ::Any, ::Any
) where {𝔽}

function _compose!(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}
    # perform the multiplication “safe”, that is, even when providing
    # one of the inputs `g,h`` and as output `k`
    (k === g || k === h) ? copyto!(k, g * h) : mul!(k, g, h)
    return k
end

Base.inv(e::Identity{<:AbstractMultiplicationGroupOperation}) = e

LinearAlgebra.det(::Identity{<:AbstractMultiplicationGroupOperation}) = true

_doc_diff_conjugate_mul = """
    diff_conjugate(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the conjugate ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = ghg^{-1}``,
which simplifies for an [`AbstractMultiplicationGroupOperation`](@ref)to ``$(_math(:d))(c_g(h))[X] = gXg^{-1}``.
"""

@doc "$(_doc_diff_conjugate_mul)"
diff_conjugate(
    ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any, ::Any, ::Any
) where {𝔽}

@doc "$(_doc_diff_conjugate_mul)"
function diff_conjugate!(
        G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, h, X
    ) where {𝔽}
    inv_right_compose!(G, Y, X, g) # Y = Xg^{-1}
    compose!(G, Y, g, Y) # Y = gY
    return Y
end

_doc_diff_inv_mult = """
    diff_inv(G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, g, X)
    diff_inv!(G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, X)

Compute the value of differential ``$(_math(:d))ι_{$(_math(:G))}(g)[X]`` of matrix inversion ``ι_{$(_math(:G))}(g) := g^{-1}`` at ``X ∈ 𝔤``
in the [`LieAlgebra`](@ref) ``𝔤`` of the [`LieGroup`](@ref) `G`.

The (classical) differential ``$(_math(:D))ι_{$(_math(:G))}(g): T_g$(_math(:G)) → T_{g^{-1}}$(_math(:G))`` reads

```math
  $(_math(:D))ι_{$(_math(:G))}(g)[W] = -g^{-1}Wg^{-1} = -Xg^{-1} = -g^{-1}(gXg^{-1}) = -g^{-1}$(_math(:Ad))(g)[X] = V ∈ T_{g^{-1}}$(_math(:G)),
```

see e.g. [Giles:2008](@cite). To bring this back to the Lie algebra, we Write ``V = g^{-1}Y ∈ T_{g^{-1}}$(_math(:G))``
for some ``Y ∈ 𝔤`` and obtain

```math
  $(_math(:d)) ι_{$(_math(:G))}(g)[X] = -$(_math(:Ad))(g)[X] ∈ 𝔤,
```

where we use ``$(_math(:d))`` to denote the differential in the Lie algebra.
"""

@doc "$(_doc_diff_inv_mult)"
diff_inv(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, g, X::Number) where {𝔽} = -X

function diff_inv(
        ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation},
        p::AbstractArray{<:Number, 0},
        X::AbstractArray{<:Number, 0},
    ) where {𝔽}
    p_inv = inv(p[])
    return -(p[] * X * p_inv)
end

@doc "$(_doc_diff_inv_mult)"
function diff_inv!(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, p, X) where {𝔽}
    Z = X / p
    mul!(Y, p, Z)
    Y .*= -1
    return Y
end

_doc_diff_left_compose_mult = """
    diff_left_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_left_compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its first (left) argument `g`.

Another interpretation is to consider a function where we do a fixed multiplication from the right with `h`.
i..e. the right group multiplication function ``ρ_h(g) = g$(_math(:∘))h``.

The differential simplifies for an [`AbstractMultiplicationGroupOperation`](@ref) to

```math
$(_math(:d)) ρ_h(g)[X] = h^{-1}Xh = $(_math(:Ad))(h^{-1})[X] ∈ 𝔤,
```

where ``$(_math(:Ad))`` denotes the [`adjoint`](@ref).
"""

@doc "$(_doc_diff_left_compose_mult)"
diff_left_compose(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_left_compose_mult)"
function diff_left_compose!(
        G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, h, X
    ) where {𝔽}
    return copyto!(LieAlgebra(G), Y, (h \ X) * h)
end

_doc_diff_right_compose_mult = """
    diff_right_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_right_compose!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its second (right) argument `h`.

Another interpretation is to consider a function where we do a fixed multiplication from the left with `g`.
i..e. the left group multiplication function ``λ_g(h) = g$(_math(:∘))h``.

It reads for an [`AbstractMultiplicationGroupOperation`](@ref) ``$(_math(:d))λ_g(h)[X] = X``.
"""

@doc "$(_doc_diff_right_compose_mult)"
diff_right_compose(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_right_compose_mult)"
function diff_right_compose!(
        G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, h, X
    ) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_exp_mult = """
    exp(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, X)
    exp!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, X)

Compute the Lie group exponential on a [`LieGroup`](@ref) with an [`AbstractMultiplicationGroupOperation`](@ref),
which simplifies to the [matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential).

This can be computed in-place of `g`.
"""

@doc "$(_doc_exp_mult)"
ManifoldsBase.exp(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any) where {𝔽}

@doc "$(_doc_exp_mult)"
function ManifoldsBase.exp!(
        ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, g, X
    ) where {𝔽}
    copyto!(g, exp(X))
    return g
end

_doc_identity_element_mult = """
    identity_element(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`AbstractMultiplicationGroupOperation`](@ref) is the one-element or identity array.
"""

@doc "$(_doc_identity_element_mult)"
identity_element(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}) where {𝔽}

function identity_element(
        G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Type{<:AbstractArray{T, 2}}
    ) where {𝔽, T}
    (N, M) = representation_size(G.manifold)
    return Matrix{T}(I(N))
end

@doc "$(_doc_identity_element_mult)"
identity_element!(::LieGroup{𝔽, <:MatrixMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(
        ::LieGroup{𝔽, <:MatrixMultiplicationGroupOperation}, e::AbstractMatrix
    ) where {𝔽}
    return copyto!(e, LinearAlgebra.I)
end

function ManifoldsBase.inner(
        G::LieAlgebra{𝔽, <:AbstractMultiplicationGroupOperation}, X, Y
    ) where {𝔽}
    return dot(X, Y)
end

_doc_inv_mult = """
    inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g)
    inv!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for an [`AbstractMultiplicationGroupOperation`](@ref)
simplifies to the multiplicative inverse ``g^{-1}``. This can be done in-place of `h`.
"""

@doc "$(_doc_inv_mult)"
Base.inv(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any) where {𝔽}

@doc "$(_doc_inv_mult)"
inv!(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, h, g) where {𝔽}

function _inv!(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, h, g) where {𝔽}
    copyto!(h, inv(g))
    return h
end

# Compute g^{-1}h more efficient than inverting g
function inv_left_compose!(
        ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, k, g, h
    ) where {𝔽}
    return copyto!(k, g \ h)
end
# Compute g∘h^{-1} more efficient than inverting h
function inv_right_compose!(
        ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, k, g, h
    ) where {𝔽}
    return copyto!(k, g / h)
end

_doc_lie_bracket_mult = """
    lie_bracket(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, X, Y)
    lie_bracket!(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`MatrixMultiplicationGroupOperation`](@ref) yields the
commutator bracket

```math
[X, Y] = XY-YX
```

The computation can be done in-place of `Z`.
"""

@doc "$(_doc_lie_bracket_mult)"
lie_bracket(::LieAlgebra{𝔽, MatrixMultiplicationGroupOperation}, ::Any, ::Any) where {𝔽}

@doc "$(_doc_lie_bracket_mult)"
function lie_bracket!(
        ::LieAlgebra{𝔽, O, <:LieGroup{𝔽, O}}, Z, X, Y
    ) where {𝔽, O <: MatrixMultiplicationGroupOperation}
    Z .= X * Y - Y * X
    return Z
end

_doc_log_mult = """
    log(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g)
    log!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, X, g)

Compute the Lie group logarithm on a [`LieGroup`](@ref) with a concrete instance of [`AbstractMultiplicationGroupOperation`](@ref),
which simplifies to the [(matrix) logarithm](https://en.wikipedia.org/wiki/Logarithm_of_a_matrix).

This can be computed in-place of `X`.
"""

@doc "$(_doc_log_mult)"
ManifoldsBase.log(::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, ::Any) where {𝔽}

function ManifoldsBase.log(
        G::LieGroup{𝔽, O}, ::Identity{O}
    ) where {𝔽, O <: AbstractMultiplicationGroupOperation}
    return zero_vector(LieAlgebra(G))
end
function ManifoldsBase.log(
        G::LieGroup{𝔽, O}, ::Identity{O}, T::Type
    ) where {𝔽, O <: AbstractMultiplicationGroupOperation}
    return zero_vector(LieAlgebra(G), T)
end

@doc "$(_doc_log_mult)"
function ManifoldsBase.log!(
        ::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, X, g
    ) where {𝔽}
    copyto!(X, log(g))
    return X
end
function ManifoldsBase.log!(
        G::LieGroup{𝔽, O}, X, ::Identity{O}
    ) where {𝔽, O <: AbstractMultiplicationGroupOperation}
    zero_vector!(LieAlgebra(G), X)
    return X
end

LinearAlgebra.mul!(q, ::Identity{MatrixMultiplicationGroupOperation}, p) = copyto!(q, p)
function LinearAlgebra.mul!(
        q::AbstractMatrix, p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}
    )
    return copyto!(q, p)
end
function LinearAlgebra.mul!(
        q::Union{AbstractMatrix},
        ::Identity{MatrixMultiplicationGroupOperation},
        ::Identity{MatrixMultiplicationGroupOperation},
    )
    return copyto!(q, I)
end
function LinearAlgebra.mul!(
        q::Identity{MatrixMultiplicationGroupOperation},
        ::Identity{MatrixMultiplicationGroupOperation},
        ::Identity{MatrixMultiplicationGroupOperation},
    )
    return q
end

Base.one(e::Identity{<:AbstractMultiplicationGroupOperation}) = e

function push_forward_tangent!(G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, X; kwargs...) where {𝔽}
    return copyto!(Y, g * X)
end
function pull_back_tangent!(G::LieGroup{𝔽, <:AbstractMultiplicationGroupOperation}, Y, g, X; kwargs...) where {𝔽}
    return copyto!(Y, g \ X)
end
