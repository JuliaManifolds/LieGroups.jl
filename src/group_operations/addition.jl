"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`compose(G::LieGroup{𝔽,M,AdditionGroupOperation}, a,b) = a+b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

const _AdditionLieGroup =
    LieGroup{𝔽,M,AdditionGroupOperation} where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}

#
#
# Handle interactions of `+` and `-` with the identity element
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

function compose(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, g, h
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return g + h
end

function compose!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, k, g, h
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    ManifoldsBase.copyto!(G, k, g + h)
    return k
end

function Base.exp(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, X
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return X
end

function ManifoldsBase.exp!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, g, X
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return copyto!(G, g, X)
end

function identity_element!(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, e
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return fill!(e, 0)
end

function Base.inv(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, g
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return -g
end
function inv!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, h, g
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return copyto!(G, h, -g)
end

function is_identity(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, h; kwargs...
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return ManifoldsBase.isapprox(G, Identity{AdditionGroupOperation}(), h; kwargs...)
end
function is_identity(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, h::Identity{AdditionGroupOperation}; kwargs...
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return true
end
function is_identity(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, h::Identity; kwargs...
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return true
end

function Base.log(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, q
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return q
end
function Base.log(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return zero_vector(G, identity_element(G))
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, X, q
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return copyto!(G, X, q)
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}
) where {𝔽,M<:ManifoldsBase.AbstractManifold{𝔽}}
    return zero_vector!(G, X, identity_element(G))
end
