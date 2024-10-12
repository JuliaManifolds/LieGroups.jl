"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`_compose(G::LieGroup{𝔽,AdditionGroupOperation}, a,b) = a+b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

#
#
# Handle interactions of `+` and `-` with the identity element
# still necessary? This is handled in compose now anyways
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(g, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(::Any, ::Identity{AdditionGroupOperation}) = g

Base.:*(e::Identity{AdditionGroupOperation}, p) = e
Base.:*(::Any, e::Identity{AdditionGroupOperation}) = e
Base.:*(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e

function _compose(::LieGroup{𝔽,AdditionGroupOperation}, g, h) where {𝔽}
    return g + h
end

function _compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}
    ManifoldsBase.copyto!(G, k, g + h)
    return k
end

function Base.exp(::LieGroup{𝔽,AdditionGroupOperation}, X) where {𝔽}
    return X
end

function ManifoldsBase.exp!(G::LieGroup{𝔽,AdditionGroupOperation}, g, X) where {𝔽}
    return copyto!(G, g, X)
end

function identity_element!(::LieGroup{𝔽,AdditionGroupOperation}, e) where {𝔽}
    return fill!(e, 0)
end

function Base.inv(::LieGroup{𝔽,AdditionGroupOperation}, g) where {𝔽}
    return -g
end
function inv!(G::LieGroup{𝔽,AdditionGroupOperation}, h, g) where {𝔽}
    return copyto!(G, h, -g)
end

function is_identity(G::LieGroup{𝔽,AdditionGroupOperation}, h; kwargs...) where {𝔽}
    return ManifoldsBase.isapprox(G, Identity{AdditionGroupOperation}(), h; kwargs...)
end
function is_identity(
    G::LieGroup{𝔽,AdditionGroupOperation}, h::Identity{AdditionGroupOperation}; kwargs...
) where {𝔽}
    return true
end
function is_identity(
    G::LieGroup{𝔽,AdditionGroupOperation}, h::Identity; kwargs...
) where {𝔽}
    return true
end

function Base.log(::LieGroup{𝔽,AdditionGroupOperation}, q) where {𝔽}
    return q
end
function Base.log(
    ::LieGroup{𝔽,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}
) where {𝔽}
    return zero_vector(G, identity_element(G))
end
function ManifoldsBase.log!(G::LieGroup{𝔽,AdditionGroupOperation}, X, q) where {𝔽}
    return copyto!(G, X, q)
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}
) where {𝔽}
    return zero_vector!(G, X, identity_element(G))
end
