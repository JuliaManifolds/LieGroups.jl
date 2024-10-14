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
# still necessary? This is handled in compose now anyways - but just to be safe, these can be kept
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(g, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(g, ::Identity{AdditionGroupOperation}) = g

function _compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h) where {𝔽}
    ManifoldsBase.copyto!(G, k, g + h)
    return k
end

"""
"""
Base.exp(
    ::LieGroup{𝔽,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}, X
) where {𝔽}

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

"""
"""
ManifoldsBase.log(
    G::LieGroup{𝔽,AdditionGroupOperation}, ::Identity{AdditionGroupOperation}, q
) where {𝔽}

function ManifoldsBase.log!(
    G::LieGroup{𝔽,AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}, q
) where {𝔽}
    return copyto!(G, X, q)
end
