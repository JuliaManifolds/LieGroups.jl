"""
    AbstractMultiplicationGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `*` being overloaded, for example
`_compose(G::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, a, b) = a * b`
"""
abstract type AbstractMultiplicationGroupOperation <: AbstractGroupOperation end

"""
    AbstractMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation

A grou poperation that is realised by a matrix multiplication.
"""
struct MatrixMultiplicationGroupOperation <: AbstractMultiplicationGroupOperation end

Base.:*(e::Identity{AbstractMultiplicationGroupOperation}) = e
Base.:*(::Identity{MatrixMultiplicationGroupOperation}, p::AbstractMatrix) = p
Base.:*(p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}) = p
function Base.:*(
    e::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return e
end
function Base.:*(
    ::Identity{AbstractMultiplicationGroupOperation}, e::Identity{AdditionGroupOperation}
)
    return e
end
function Base.:*(
    e::Identity{AdditionGroupOperation}, ::Identity{AbstractMultiplicationGroupOperation}
)
    return e
end

Base.:/(p, ::Identity{AbstractMultiplicationGroupOperation}) = p
Base.:/(::Identity{AbstractMultiplicationGroupOperation}, p) = inv(p)
function Base.:/(
    e::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return e
end

Base.:\(p, ::Identity{AbstractMultiplicationGroupOperation}) = inv(p)
Base.:\(::Identity{AbstractMultiplicationGroupOperation}, p) = p
function Base.:\(
    e::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return e
end

Base.inv(e::Identity{AbstractMultiplicationGroupOperation}) = e

LinearAlgebra.det(::Identity{AbstractMultiplicationGroupOperation}) = true
LinearAlgebra.adjoint(e::Identity{AbstractMultiplicationGroupOperation}) = e

_doc_compose_mult = """
    compose(G::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, k, g, h)

Compute the group operation composition of `g` and `h` with respect to
an [`AbstractMultiplicationGroupOperation`](@ref) on `G`, which falls back to calling
`g*h`, where `*` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_mult)"
compose(::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_mult)"
compose!(::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}

function _compose!(::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, k, g, h) where {𝔽}
    k .= g * h
    return k
end

_doc_identity_element_mult = """
    identity_element(G::LieGroup{𝔽,AbstractMultiplicationGroupOperation})
    identity_element!(G::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for the [`AdditionGroupOperation`](@ref) is the one-element or identity array.
"""

@doc "$(_doc_identity_element_mult)"
identity_element(::LieGroup{𝔽,AbstractMultiplicationGroupOperation}) where {𝔽}

@doc "$(_doc_identity_element_mult)"
identity_element!(::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, e) where {𝔽}
function identity_element!(
    ::LieGroup{𝔽,AbstractMultiplicationGroupOperation}, e::AbstractMatrix
) where {𝔽}
    return copyto!(e, LinearAlgebra.I)
end

LinearAlgebra.mul!(q, ::Identity{AbstractMultiplicationGroupOperation}, p) = copyto!(q, p)
function LinearAlgebra.mul!(
    q::AbstractMatrix, p::AbstractMatrix, ::Identity{MatrixMultiplicationGroupOperation}
)
    return copyto!(q, p)
end
function LinearAlgebra.mul!(
    q::AbstractMatrix,
    ::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return copyto!(q, I)
end
function LinearAlgebra.mul!(
    q,
    ::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return copyto!(q, one(q))
end
function LinearAlgebra.mul!(
    q::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
    ::Identity{AbstractMultiplicationGroupOperation},
)
    return q
end
Base.one(e::Identity{AbstractMultiplicationGroupOperation}) = e
