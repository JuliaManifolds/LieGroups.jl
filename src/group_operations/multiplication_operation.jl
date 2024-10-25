"""
    AbstractMultiplicationGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `*` being overloaded, for example
`_compose(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, a, b) = a * b`
"""
abstract type AbstractMultiplicationGroupOperation <: AbstractGroupOperation end

"""
    AbstractMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A grou poperation that is realised by a matrix multiplication.
"""
struct MatrixMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

Base.:*(e::Identity{<:AbstractMultiplicationGroupOperation}) = e
Base.:*(::Identity{MatrixMultiplicationGroupOperation}, p::AbstractMatrix) = p
Base.:*(p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}) = p
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

Compute the group operation composition of `g` and `h` with respect to
an [`AbstractMultiplicationGroupOperation`](@ref) on `G`, which falls back to calling
`g*h`, where `*` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_mult)"
compose(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_mult)"
compose!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}

function _compose!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}
    # perform the multiplication “safe”, that is, even when providing
    # one of the inputs `g,h`` and as output `k`
    (k === g || k === h) ? copyto!(k, g * h) : mul!(k, g, h)
    return k
end

Base.inv(e::Identity{<:AbstractMultiplicationGroupOperation}) = e

LinearAlgebra.det(::Identity{<:AbstractMultiplicationGroupOperation}) = true
LinearAlgebra.adjoint(e::Identity{<:AbstractMultiplicationGroupOperation}) = e

_doc_diff_conjugate_add = """
    diff_conjugate(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X)

Compute the differential of the conjutage ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = ghg^{-1}``,
which simplifies for an [`AbstractMultiplicationGroupOperation`](@ref) to ``D(c_g(h))[X] = gXg^{-1}``.
"""

@doc "$(_doc_diff_conjugate_add)"
diff_conjugate(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_conjugate_add)"
function diff_conjugate!(
    G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, Y, g, h, X
) where {𝔽}
    inv_right_compose!(G, Y, X, g) # Y = Xg^{-1}
    compose!(G, Y, g, Y) # Y = gY
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_identity_element_mult = """
    identity_element(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for an [`AbstractMultiplicationGroupOperation`](@ref) is the one-element or identity array.
"""

@doc "$(_doc_identity_element_mult)"
identity_element(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}) where {𝔽}

@doc "$(_doc_identity_element_mult)"
identity_element!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(
    ::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, e::AbstractMatrix
) where {𝔽}
    return copyto!(e, LinearAlgebra.I)
end

_doc_inv_mult = """
    inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperationroupOperation}, g)
    inv!(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for an [`AbstractMultiplicationGroupOperation`](@ref)
simplifies to the multiplicative inverse ``g^{-1}``. This can be done in-place of `h`.
"""

@doc "$(_doc_inv_mult)"
Base.inv(G::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, g) where {𝔽}

function Base.inv(
    ::LieGroup{𝔽,O}, e::Identity{O}
) where {𝔽,O<:AbstractMultiplicationGroupOperation}
    return e
end

@doc "$(_doc_inv_mult)"
function inv!(::LieGroup{𝔽,<:AbstractMultiplicationGroupOperation}, h, g) where {𝔽}
    copyto!(h, inv(g))
    return h
end
function inv!(
    G::LieGroup{𝔽,O}, q, ::Identity{O}
) where {𝔽,O<:AbstractMultiplicationGroupOperation}
    return identity_element!(G, q)
end

# Compute g^{-1}h more efficient than inverting g
function inv_left_compose!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, g \ h)
end
# Compute g∘h^{-1} more efficient than inverting h
function inv_right_compose!(
    ::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, k, g, h
) where {𝔽}
    return copyto!(k, g / h)
end

_doc_lie_bracket_mult = """
    lie_bracket(::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, X, Y)
    lie_bracket!(::LieGroup{𝔽,MatrixMultiplicationGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`MatrixMultiplicationGroupOperation`](@ref) yields the
commutator bracket

```math
[X, Y] = XY-YX
```

The computation can be done in-place of `Z`.
"""

@doc "$(_doc_lie_bracket_mult)"
lie_bracket(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, X, Y) where {𝔽}

@doc "$(_doc_lie_bracket_mult)"
function lie_bracket!(::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation}, Z, X, Y) where {𝔽}
    mul!(Z, X, Y)
    mul!(Z, Y, X, -1, true)
    return Z
end

LinearAlgebra.mul!(q, ::Identity{<:AbstractMultiplicationGroupOperation}, p) = copyto!(q, p)
function LinearAlgebra.mul!(
    q::AbstractMatrix, p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}
)
    return copyto!(q, p)
end
function LinearAlgebra.mul!(
    q::AbstractMatrix,
    ::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return copyto!(q, I)
end
function LinearAlgebra.mul!(
    q,
    ::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return copyto!(q, one(q))
end
function LinearAlgebra.mul!(
    q::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
    ::Identity{<:AbstractMultiplicationGroupOperation},
)
    return q
end
Base.one(e::Identity{<:AbstractMultiplicationGroupOperation}) = e
