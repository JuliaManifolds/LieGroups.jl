#
#
# Power Lie groups: work element wise

@doc """
    PowerLieGroup(G::LieGroup, args...; kwargs...)
    PowerLieGroup(M::AbstractManifold, op::AbstractGroupOperation, args...; kwargs...)
    (L::LueGroup)^(n...) = PowerLieGroup(L, n...)

Generate the [`LieGroup`](@ref) of the `n`-th power of a Lie group `G` or manifold `M`.
If passed a Lie group `G`, the group operation is the same as on `G`, but applied elementwise.

Bot the arguments `args...` as well as the keyword arguments `kwargs...` are passed on to
the constructor of the [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`).
This especially includes the `size` of the manifold and allows to specify a [`NestedPowerRepresentation`](@extref `ManifoldsBase.NestedPowerRepresentation`).

"""
PowerLieGroup(::AbstractManifold, args...; kwargs...)

function PowerLieGroup(M::AbstractManifold, op::AbstractGroupOperation, args...; kwargs...)
    pM = Manifolds.PowerManifold(M, args...; kwargs...)
    return LieGroup(pM, op)
end
function PowerLieGroup(G::LieGroup, args...; kwargs...)
    return PowerLieGroup(G.manifold, G.op, args...; kwargs...)
end

Base.:^(G::LieGroup, n...) = PowerLieGroup(G, n...)

function compose!(PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, k, g, h) where {𝔽,Op}
    PM = PG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PG.op)
    for i in get_iterator(PM)
        compose!(
            G,
            ManifoldsBase._write(PM, rep_size, k, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, h, i),
        )
    end
    return k
end

function ManifoldsBase.check_size(
    PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, g
) where {𝔽,Op}
    return ManifoldsBase.check_size(PG.manifold, g)
end
function ManifoldsBase.check_size(
    PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, g, X
) where {𝔽,Op}
    return ManifoldsBase.check_size(PG.manifold, g, X)
end

function exp!(PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, h, g, X, t::Number=1) where {𝔽,Op}
    PM = PG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PG.op)
    for i in get_iterator(PM)
        exp!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return h
end

function identity_element!(PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, e) where {𝔽,Op}
    PM = PG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PG.op)
    for i in get_iterator(PM)
        identity_element!(G, ManifoldsBase._write(PM, rep_size, e, i))
    end
    return e
end

function inv!(PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, h, g) where {𝔽,Op}
    PM = PG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PG.op)
    for i in get_iterator(PM)
        inv!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, g, i),
        )
    end
    return h
end

function log!(PG::LieGroup{𝔽,Op,<:AbstractPowerManifold}, X, g) where {𝔽,Op}
    PM = PG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PG.op)
    for i in get_iterator(PM)
        log!(
            G,
            ManifoldsBase._write(PM, rep_size, X, i),
            ManifoldsBase._read(PM, rep_size, g, i),
        )
    end
    return X
end

function Base.show(
    io::IO, G::LieGroup{𝔽,O,<:ManifoldsBase.AbstractPowerManifold}
) where {𝔽,O<:AbstractGroupOperation}
    M = G.manifold.manifold
    size = Manifolds.get_parameter(G.manifold.size)
    return print(io, "PowerLieGroup($(M), $(G.op), $(join(size, ", ")))")
end
