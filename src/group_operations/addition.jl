"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`compose(G::LieGroup{𝔽,M,AdditionGroupOperation}, a,b) = a+b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

const _AdditionLieGroup = LieGroup{𝔽,M,AdditionGroupOperation} where {𝔽,M}

#
#
# Handle interactions of `+` and `-` with the identity element
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(p, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(p, ::Identity{AdditionGroupOperation}) = g

Base.:*(e::Identity{AdditionGroupOperation}, p) = e
Base.:*(p, e::Identity{AdditionGroupOperation}) = e
Base.:*(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e

compose(::LieGroup{𝔽,M,AdditionGroupOperation}, g, h) where {𝔽,M} = g + h

function compose!(::LieGroup{𝔽,M,AdditionGroupOperation}, k, g, h) where {𝔽,M}
    k .= g .+ h
    return k
end

Base.exp(::LieGroup{𝔽,M,AdditionGroupOperation}, X) where {𝔽,M} = X

function ManifoldsBase.exp!(G::LieGroup{𝔽,M,AdditionGroupOperation}, g, X) where {𝔽,M}
    return copyto!(G, g, X)
end

Base.log(::LieGroup{𝔽,M,AdditionGroupOperation}, q) where {𝔽,M} = q
function Base.log(
    ::LieGroup{𝔽,M,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}
) where {𝔽,M}
    return zero_vector(G, identity_element(G))
end
function ManifoldsBase.log!(G::LieGroup{𝔽,M,AdditionGroupOperation}, X, q) where {𝔽,M}
    return copyto!(G, X, q)
end
function ManifoldsBase.log!(
    G::LieGroup{𝔽,M,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}
) where {𝔽,M}
    return zero_vector!(G, X, identity_element(G))
end
