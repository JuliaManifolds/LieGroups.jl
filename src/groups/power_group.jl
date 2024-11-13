#
#
# Power Lie groups: work element wise

@doc """
    PowerLieGroup(G::LieGroup, size::Int...; kwargs...)
    PowerLieGroup(M::AbstractManifold, op::AbstractGroupOperation, size:Int...; kwargs...)
    (L::LueGroup)^(n...) = PowerLieGroup(L, n...)

Generate the [`LieGroup`](@ref) of the `n`-th power of a Lie group `G` or manifold `M`.
If passed a Lie group `G`, the group operation is the same as on `G`, but applied elementwise.

The keyword arguments `kwargs...` are passed on to the constructor of the [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`).
"""
PowerLieGroup(::AbstractManifold, size::Int...; kwargs...)

function PowerLieGroup(
    M::AbstractManifold, op::AbstractGroupOperation, size::Int...; kwargs...
)
    pM = Manifolds.PowerManifold(M, size...; kwargs...)
    return LieGroup(pM, op)
end
function PowerLieGroup(G::LieGroup, size::Int...; kwargs...)
    return PowerLieGroup(G.manifold, G.op, size...; kwargs...)
end

Base.:^(G::LieGroup, n...) = PowerLieGroup(G, n...)

function Base.show(
    io::IO, G::LieGroup{𝔽,O,<:ManifoldsBase.AbstractPowerManifold}
) where {𝔽,O<:AbstractGroupOperation}
    M = G.manifold.manifold
    size = Manifolds.get_parameter(G.manifold.size)
    return print(io, "PowerLieGroup($(M), $(G.op), $(join(size, ", ")))")
end
