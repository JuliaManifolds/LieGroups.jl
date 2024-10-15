"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`_compose(G::LieGroup{𝔽,AdditionGroupOperation}, a,b) = a+b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

#
#
# Handle interactions of `+` and `-` with the identity element, though they are
# also already handled on the `compose()` level
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(g, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(g, ::Identity{AdditionGroupOperation}) = g

_doc_compose_add = """
    compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h)

Copmute the group operation composition of `g` and `h` with respect to
the [`AdditionGroupOperation`](@ref) on `G`, which falls back to calling
`g+h`, where `+` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_add)"
compose(::LieGroup{𝔽,AdditionGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_add)"
compose!(::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}

function _compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}
    ManifoldsBase.copyto!(G, k, g + h)
    return k
end

_doc_exp_add = """
    exp(G::LieGroup{𝔽,AdditionGroupOperation}, e::Identity{AdditionGroupOperation}, X, t=1)
    exp!(G::LieGroup{𝔽,AdditionGroupOperation}, g, e::Identity{AdditionGroupOperation}, X, t)

Compute the Lie group exponential on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `g`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``g=0+X=X``.
"""

@doc "$(_doc_exp_add)"
Base.exp(
    ::LieGroup{𝔽,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}, X, t
) where {𝔽}

@doc "$(_doc_exp_add)"
function ManifoldsBase.exp!(
    G::LieGroup{𝔽,AdditionGroupOperation},
    g,
    ::Identity{AdditionGroupOperation},
    X,
    t::Number=1,
) where {𝔽}
    return copyto!(G, g, X)
end

function identity_element!(::LieGroup{𝔽,AdditionGroupOperation}, e) where {𝔽}
    return fill!(e, 0)
end

function inv!(G::LieGroup{𝔽,AdditionGroupOperation}, h, g) where {𝔽}
    return copyto!(G, h, -g)
end

_doc_log_add = """
    log(G::LieGroup{𝔽,AdditionGroupOperation}, e::Identity{AdditionGroupOperation}, g)
    log!(G::LieGroup{𝔽,AdditionGroupOperation}, X, e::Identity{AdditionGroupOperation}, g)

Compute the Lie group logarithm on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `X`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``X=g-0=g``.
"""

@doc "$(_doc_log_add)"
ManifoldsBase.log(
    G::LieGroup{𝔽,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}, q
) where {𝔽}

@doc "$(_doc_log_add)"
function ManifoldsBase.log!(
    G::LieGroup{𝔽,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}, g
) where {𝔽}
    return copyto!(G, X, g)
end
