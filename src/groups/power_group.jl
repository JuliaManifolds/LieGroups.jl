#
#
# Power Lie groups: work element wise

@doc """
    PowerGroupOperation{O<:AbstractGroupOperation} <: AbstractGroupOperation

A struct do model a that a certain group operation is applied element-wise on a [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`).

# Constructor

    PowerGroupOperation(o::AbstractGroupOperation)
"""
struct PowerGroupOperation{O <: AbstractGroupOperation} <: AbstractGroupOperation
    op::O
end

@doc """
    PowerLieGroup(G::LieGroup, args...; kwargs...)
    (G::LieGroup)^(n::Integer) = PowerLieGroup(G, n)

Generate the [`LieGroup`](@ref) of the `n`-th power of a Lie group `G` or manifold `M`.
If passed a Lie group `G`, the group operation on the [`PowerLieGroup`](@ref) is the same as on `G`,
but applied elementwise. Internally, the corresponding [`PowerGroupOperation`](@ref) is created.
If you pass a manifold `M`, you have to provide the corresponding [`PowerGroupOperation`](@ref) yourself.

Both the arguments `args...` as well as the keyword arguments `kwargs...` are passed on to
the constructor of the [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`).
This especially includes the `size` of the manifold and allows to specify a [`NestedPowerRepresentation`](@extref `ManifoldsBase.NestedPowerRepresentation`).
"""
PowerLieGroup(::AbstractManifold, args...; kwargs...)

function PowerLieGroup(G::LieGroup, args...; kwargs...)
    M = G.manifold
    pM = ManifoldsBase.PowerManifold(M, args...; kwargs...)
    return LieGroup(pM, PowerGroupOperation(G.op))
end

Base.:^(G::LieGroup, n::Integer) = PowerLieGroup(G, n)

function _compose!(
        PoG::LieGroup{𝔽, Op, M}, k, g, h
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
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
        PoG::LieGroup{𝔽, Op, M}, g
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return ManifoldsBase.check_size(PoG.manifold, g)
end
function ManifoldsBase.check_size(
        ::LieGroup{𝔽, Op, M}, ::Identity{Op}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return nothing
end
function ManifoldsBase.check_size(
        G::LieGroup{𝔽, Op, M}, e::Identity
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return DomainError(
        "The Identity $e is not the identity of the group, expected $(Identity(G.op))."
    )
end
function ManifoldsBase.check_size(
        PoG::LieGroup{𝔽, Op, M}, g, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return ManifoldsBase.check_size(PoG.manifold, g, X)
end

function conjugate!(
        PoG::LieGroup{𝔽, Op, M}, h, g, k
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        conjugate!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, k, i),
        )
    end
    return h
end

function diff_conjugate!(
        PoG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_conjugate!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end
function diff_conjugate!(
        PoG::LieGroup{𝔽, Op, M}, Y, g, ::Identity{Op}, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    e_g = Identity(G)
    for i in ManifoldsBase.get_iterator(PM)
        diff_conjugate!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            e_g,
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end
function diff_inv!(
        PoG::LieGroup{𝔽, Op, M}, Y, g, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_inv!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end

function diff_left_compose!(
        PoG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_left_compose!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end
function diff_left_compose!(
        PoG::LieGroup{𝔽, Op, M}, Y, ::Identity{Op}, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_left_compose!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            PoG.op.op,
            ManifoldsBase._read(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end

function diff_right_compose!(
        PoG::LieGroup{𝔽, Op, M}, Y, g, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_right_compose!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end
function diff_right_compose!(
        PoG::LieGroup{𝔽, Op, M}, Y, ::Identity{Op}, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        diff_right_compose!(
            G,
            ManifoldsBase._write(PM, rep_size, Y, i),
            PoG.op.op,
            ManifoldsBase._read(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return Y
end

function ManifoldsBase.exp!(
        PoG::LieGroup{𝔽, Op, M}, h, g, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        exp!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return h
end

function ManifoldsBase.exp!(
        PoG::LieGroup{𝔽, Op, M}, h, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        exp!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, X, i),
        )
    end
    return h
end

Base.@propagate_inbounds function Base.getindex(
        p::AbstractArray, PoG::LieGroup{𝔽, Op, M}, I::Union{Integer, Colon, AbstractVector}...,
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return get_component(base_manifold(PoG), p, I...)
end

function jacobian_exp!(
        PoG::LieGroup{𝔽, Op, M}, J, X, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    dim = manifold_dimension(PM.manifold)
    G = LieGroup(PM.manifold, PoG.op.op)
    fill!(J, 0)
    v_iter = 1
    for i in ManifoldsBase.get_iterator(PM)
        dr = v_iter:(v_iter + dim - 1)
        jacobian_exp!(G, view(J, dr, dr), ManifoldsBase._read(PM, rep_size, X, i), B)
        v_iter += dim
    end
    return J
end

function ManifoldsBase.hat!(
        Po𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, c
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PoG = Po𝔤.manifold
    PM = PoG.manifold
    rep_size = representation_size(PM.manifold)
    dim = manifold_dimension(PM.manifold)
    v_iter = 1
    𝔤 = LieAlgebra(LieGroup(PM.manifold, PoG.op.op))
    for i in ManifoldsBase.get_iterator(PM)
        hat!(𝔤, ManifoldsBase._write(PM, rep_size, X, i), c[v_iter:(v_iter + dim - 1)])
        v_iter += dim
    end
    return X
end

function LieGroups.identity_element(
        PoG::LieGroup{𝔽, Op, M}, ::Type{Vector{T}}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold, T}
    PM = PoG.manifold
    G = LieGroup(PM.manifold, PoG.op.op)
    return [identity_element(G, T) for _ in ManifoldsBase.get_iterator(PM)]
end

function identity_element!(
        PoG::LieGroup{𝔽, Op, M}, e
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        identity_element!(G, ManifoldsBase._write(PM, rep_size, e, i))
    end
    return e
end

function ManifoldsBase.inner(
        Po𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, Y
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    v = PoG = Po𝔤.manifold
    PM = PoG.manifold
    𝔤 = LieAlgebra(LieGroup(PM.manifold, PoG.op.op))
    rep_size = representation_size(PM.manifold)
    return sum(
        inner(
                𝔤,
                ManifoldsBase._read(PM, rep_size, X, i),
                ManifoldsBase._read(PM, rep_size, Y, i),
            ) for i in ManifoldsBase.get_iterator(PM)
    )
end

function _inv!(
        PoG::LieGroup{𝔽, Op, M}, h, g
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    for i in ManifoldsBase.get_iterator(PM)
        inv!(
            G,
            ManifoldsBase._write(PM, rep_size, h, i),
            ManifoldsBase._read(PM, rep_size, g, i),
        )
    end
    return h
end
function Manifolds.inv!(
        PoG::LieGroup{𝔽, Op, M}, h, ::Identity{Op}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    e_g = Identity(G.op)
    for i in ManifoldsBase.get_iterator(PM)
        inv!(G, ManifoldsBase._write(PM, rep_size, h, i), e_g)
    end
    return h
end

function lie_bracket!(
        PoA::LieAlgebra{𝔽, Op, <:LieGroup{𝔽, Op, M}}, Z, X, Y
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoA.manifold.manifold
    rep_size = representation_size(PM)
    𝔤 = LieAlgebra(LieGroup(PM.manifold, PoA.manifold.op.op))
    for i in ManifoldsBase.get_iterator(PM)
        lie_bracket!(
            𝔤,
            ManifoldsBase._write(PM, rep_size, Z, i),
            ManifoldsBase._read(PM, rep_size, X, i),
            ManifoldsBase._read(PM, rep_size, Y, i),
        )
    end
    return Z
end

function ManifoldsBase.log!(
        PoG::LieGroup{𝔽, Op, M}, X, g
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op)
    e_g = Identity(G.op)
    for i in ManifoldsBase.get_iterator(PM)
        log!(
            G,
            ManifoldsBase._write(PM, rep_size, X, i),
            ManifoldsBase._read(PM, rep_size, g, i),
        )
    end
    return X
end

function ManifoldsBase.log!(
        G::LieGroup{𝔽, Op, M}, X, ::Identity{Op}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return zero_vector!(LieAlgebra(G), X)
end
function ManifoldsBase.log!(
        G::LieGroup{𝔽, Op, M}, X, ::Identity{Op}, ::Identity{Op}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    return zero_vector!(LieAlgebra(G), X)
end
function ManifoldsBase.log!(
        PoG::LieGroup{𝔽, Op, M}, X, g, h
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = PoG.manifold
    rep_size = representation_size(PM)
    G = LieGroup(PM.manifold, PoG.op.op) # generate the single Lie group
    for i in ManifoldsBase.get_iterator(PM)
        log!(
            G,
            ManifoldsBase._write(PM, rep_size, X, i),
            ManifoldsBase._read(PM, rep_size, g, i),
            ManifoldsBase._read(PM, rep_size, h, i),
        )
    end
    return X
end

function Base.show(
        io::IO, G::LieGroup{𝔽, Op, M}
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PM = G.manifold
    POp = G.op
    L = LieGroup(PM.manifold, POp.op)
    size = get_parameter(G.manifold.size)
    return print(io, "PowerLieGroup($L, $(join(size, ", ")))")
end

function ManifoldsBase.vee!(
        Po𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, c, X
    ) where {𝔽, Op <: PowerGroupOperation, M <: ManifoldsBase.AbstractPowerManifold}
    PoG = Po𝔤.manifold
    PM = PoG.manifold
    rep_size = representation_size(PM.manifold)
    dim = manifold_dimension(PM.manifold)
    𝔤 = LieAlgebra(LieGroup(PM.manifold, PoG.op.op))
    v_iter = 1
    for i in ManifoldsBase.get_iterator(PM)
        vee!(𝔤, view(c, v_iter:(v_iter + dim - 1)), ManifoldsBase._read(PM, rep_size, X, i))
        v_iter += dim
    end
    return c
end
