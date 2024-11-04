#
#
# Lie groups over  a power manifold work as large products, where on every factor
# the same group operation is used. Hence we also only need one group operation,
# we can recognize power Lie groups when the manifold is a power manifold.

"""
    PowerLieGroup(G::LieGroup, size::Int...; kwargs...)
    PowerLieGroup(M::AbstractManifold, op::AbstractGroupOperation, size:Int...; kwargs...)
    (L::LueGroup)^(n...) = PowerLieGroup(L, n...)

Generate the [`LieGroup`](@ref) of the `n`-th power of a Lie group `G` or manifold `M`.
If passed a Lie group `G`, the group operation is the same as on `G`, but applied elementwise.

The keyword arguments `kwargs...` are passed on to the constructor of the [`PowerManifold`](@exref).
"""
PowerLieGroup(::AbstractManifold, size::Int...; kwargs...)

function PowerLieGroup(L::LieGroup, size::Int...; kwargs...)
    pM = Manifolds.PowerManifold(L.manifold, size; kwargs...)
    return LieGroup(pM, L.op)
end

function PowerLieGroup(
    M::AbstractManifold, op::AbstractGroupOperation, size::Int...; kwargs...
)
    pM = Manifolds.PowerManifold(M, size; kwargs...)
    return LieGroup(pM, L.op)
end

Base.:^(L::LieGroup, n...) = PowerLieGroup(L, n...)
