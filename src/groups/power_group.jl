#
#
# Lie groups over  a power manifold work as large products, where on every factor
# the same group operation is used. Hence we also only need one group operation,
# we can recognize power Lie groups when the manifold is a power manifold.

@doc """
    PowerLieGroup(G::LieGroup, size::Int...; kwargs...)
    PowerLieGroup(M::AbstractManifold, op::AbstractGroupOperation, size:Int...; kwargs...)
    (L::LueGroup)^(n...) = PowerLieGroup(L, n...)

Generate the [`LieGroup`](@ref) of the `n`-th power of a Lie group `G` or manifold `M`.
If passed a Lie group `G`, the group operation is the same as on `G`, but applied elementwise.

The keyword arguments `kwargs...` are passed on to the constructor of the [`PowerManifold`](@exref).
"""
PowerLieGroup(::AbstractManifold, size::Int...; kwargs...)

function PowerLieGroup(
    M::AbstractManifold, op::AbstractGroupOperation, size::Int...; kwargs...
)
    pM = Manifolds.PowerManifold(M, size...; kwargs...)
    return LieGroup(pM, op)
end
function PowerLieGroup(L::LieGroup, size::Int...; kwargs...)
    return PowerLieGroup(L.manifold, L.op, size...; kwargs...)
end

Base.:^(L::LieGroup, n...) = PowerLieGroup(L, n...)

function show(io::IO, G::LieGroup{𝔽,O,<:ManifoldsBase.AbstractPowerManifold}) where {𝔽,O}
    M = G.manifold.manifold
    size = get_parameter(G.manifold.size)
    return print(io, "PowerLieGroup($(M), $(G.op), $(join(size, ", ")))")
end
