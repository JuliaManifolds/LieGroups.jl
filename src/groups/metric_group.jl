"""
    MetricLieGroup{𝔽,M<:AbstractLieGroup{𝔽},G<:AbstractMetric} <: AbstractDecoratorManifold{𝔽}

Equip an [`AbstractLieGroup`](@ref) `L` explicitly with an [`AbstractMetric`](@extref `ManifoldsBase.AbstractMetric`) `G`.

By default every [`AbstractLieGroup`](@ref) `L` is assumed to implicitly
implement all functions on their [`LieAlgebra`](@ref) with respect to a metric,
that corresponds to a certain “default metric”, either because it is a widely recognized metric
or because the first implementation within `LieGroups.jl` was done with respect to this metric.

This default is usually indicated by checking that [`metric(L)`](@extref `Manifolds`) returns [`DefaultMetric`](@extref `Manifolds.DefaultMetric`).

This decorator type allows to explicitly equip the Lie group with a different metric.
Note that for any functions unrelated to the metric, this new Lie group will simply forward
to the Lie group that is internally stored.
This for example includes group operations like [`compose`], `inv`, or `identity_element`,
or functions that are already otherwise passed on to the inner manifold like [`manifold_dimension`](@extref `ManifoldsBase.manifold_dimension-Tuple{AbstractManifold}`).

# Constructor

    MetricManifold(M, G)

Generate the [`AbstractLieGroup`](https://juliamanifolds.github.io/Manifolds.jl/latest/interface.html#ManifoldsBase.AbstractLieGroup) `M` as a manifold with the `AbstractMetric` `G`.
"""
struct MetricLieGroup{
        𝔽,
        O <: AbstractGroupOperation,
        M <: ManifoldsBase.AbstractManifold{𝔽},
        L <: LieGroup{𝔽, O, M},
        G <: AbstractMetric,
    } <: AbstractLieGroup{𝔽, O, M}
    lie_group::L
    metric::G
end

function ManifoldsBase.submanifold_component(G::MetricLieGroup, g, i::Val)
    return submanifold_component(base_lie_group(G), g, i)
end
function ManifoldsBase.submanifold_component(G::MetricLieGroup, g, i::Integer)
    return submanifold_component(base_lie_group(G), g, i)
end
function ManifoldsBase.submanifold_components(G::MetricLieGroup, g)
    return submanifold_components(G.lie_group, g)
end

#
#
# passthrough for group operations

Identity(G::MetricLieGroup) = Identity(G.lie_group)

adjoint(G::MetricLieGroup, g, X) = adjoint(G.lie_group, g, X)
adjoint!(G::MetricLieGroup, Y, g, X) = adjoint!(G.lie_group, Y, g, X)

base_lie_group(G::MetricLieGroup) = G.lie_group
ManifoldsBase.base_manifold(G::MetricLieGroup) = base_manifold(base_lie_group(G))

_compose(G::MetricLieGroup, g, h) = compose(base_lie_group(G), g, h)
_compose!(G::MetricLieGroup, k, g, h) = compose!(base_lie_group(G), k, g, h)

conjugate(G::MetricLieGroup, g, h) = conjugate(base_lie_group(G), g, h)
conjugate!(G::MetricLieGroup, k, g, h) = conjugate!(base_lie_group(G), k, g, h)

diff_conjugate(G::MetricLieGroup, g, h, X) = diff_conjugate(base_lie_group(G), g, h, X)
diff_conjugate!(G::MetricLieGroup, Y, g, h, X) = diff_conjugate!(base_lie_group(G), Y, g, h, X)

diff_inv(G::MetricLieGroup, g, X) = diff_inv(base_lie_group(G), g, X)
diff_inv!(G::MetricLieGroup, Y, g, X) = diff_inv!(base_lie_group(G), Y, g, X)

diff_left_compose(G::MetricLieGroup, g, h, X) = diff_left_compose(base_lie_group(G), g, h, X)
diff_left_compose!(G::MetricLieGroup, Y, g, h, X) = diff_left_compose!(base_lie_group(G), Y, g, h, X)

diff_right_compose(G::MetricLieGroup, g, h, X) = diff_right_compose(base_lie_group(G), g, h, X)
diff_right_compose!(G::MetricLieGroup, Y, g, h, X) = diff_right_compose!(base_lie_group(G), Y, g, h, X)

function identity_element(G::MetricLieGroup)
    return identity_element(base_lie_group(G))
end
function identity_element(G::MetricLieGroup, T::Type)
    return identity_element(base_lie_group(G), T)
end
function identity_element!(G::MetricLieGroup, e)
    return identity_element!(base_lie_group(G), e)
end


_inv(G::MetricLieGroup, g) = inv(base_lie_group(G), g)
_inv!(G::MetricLieGroup, k, g) = inv!(base_lie_group(G), k, g)

_inv_left_compose(G::MetricLieGroup, h, g) = _inv_left_compose(base_lie_group(G), h, g)
_inv_left_compose!(G::MetricLieGroup, k, h, g) = _inv_left_compose!(base_lie_group(G), k, h, g)

_inv_right_compose(G::MetricLieGroup, h, g) = _inv_right_compose(base_lie_group(G), h, g)
_inv_right_compose!(G::MetricLieGroup, k, h, g) = _inv_right_compose!(base_lie_group(G), k, h, g)

ManifoldsBase.is_point(G::MetricLieGroup, g; kwargs...) = is_point(base_lie_group(G), g; kwargs...)
ManifoldsBase.is_point(G::MetricLieGroup, e::Identity; kwargs...) = is_point(base_lie_group(G), e; kwargs...)

function ManifoldsBase.is_vector(
        G::MetricLieGroup{𝔽, O}, e::Identity{O}, X; kwargs...
    ) where {𝔽, O <: AbstractGroupOperation}
    return is_vector(base_lie_group(G), e, X; kwargs...)
end

function Base.show(io::IO, G::MetricLieGroup)
    return print(io, "MetricLieGroup($(G.lie_group), $(G.metric))")
end
