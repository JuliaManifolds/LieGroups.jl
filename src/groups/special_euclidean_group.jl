# for (R, t)
const LeftSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:LeftSemidirectProductGroupOperation{
        <:MatrixMultiplicationGroupOperation,
        <:AdditionGroupOperation,
        LeftGroupOperationAction,
    },
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Rotations{T},<:Manifolds.Euclidean{T,ℝ}}
    },
}

# for (t, R)
const RightSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:RightSemidirectProductGroupOperation{
        <:AdditionGroupOperation,
        <:MatrixMultiplicationGroupOperation,
        LeftGroupOperationAction,
    },
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Rotations{T},<:Manifolds.Euclidean{T,ℝ}}
    },
}

"""
    SpecialEuclideanGroup{T}

The special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``. is the Lie group consisting of the
[`LeftSemidirectProductGroupOperation`](@ref) of the [`SpecialOrthogonalGroup`](@ref) and the
[`TranslationGroup`](@ref) together with the [`GroupOperationAction`](@ref)`{`[`LeftGroupOperationAction`](@ref)`}`.

To be precise, the group operation is defined on ``$(_math(:SO))(n) ⋉ $(_math(:T))(n)`` as follows:

```math
(g_1, t_1) ⋅ (g_2, t_2) = (g_1$(_math(:∘))g_2, t_1 + g_1$(_math(:⋅))t_2)
```

Analogously you can write this on elements of ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)`` as

```math
(s_1, h_1) ⋅ (s_2, h_2) = (s_1 + h_1$(_math(:⋅))s_2, h_1$(_math(:∘))h_2)
```

# Constructor
    SpecialEuclideanGroup(n; variant=:left, kwargs...)
    SpecialOrthogonalGroup(n; kwargs...) ⋉ TranslationGroup(n; kwargs...)
    TranslationGroup(n; kwargs...) ⋊ SpecialOrthogonalGroup(n; kwargs...)

Generate special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``, where the first
constructor is equivalent to the second.

All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.

The default representation for ``$(_math(:SE))(n)`` is in [affine form](https://en.wikipedia.org/wiki/Affine_group#Matrix_representation)

```math
$(_tex(:pmatrix, "g & t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")),
$(_tex(:qquad)) g ∈ $(_math(:SO))(n), t ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

Alternatively you can use the `ArrayPartition` from [`RecursiveArrayTools.jl`](https://docs.sciml.ai/RecursiveArrayTools/stable/)
to work on ``(g,t)``.

Alternatively you can use ``$(_math(:T))(n) ⋊ $(_math(:SO))(n)`` for `ArrayPartition`s ``(t,g)``.
or set in the first constructor `variant=:right`
"""
const SpecialEuclideanGroup{T} = Union{
    <:LeftSpecialEuclideanGroup{T},<:RightSpecialEuclideanGroup{T}
}

# For matrices, both should be fine, we clarify this in the access of subcomponents

function SpecialEuclideanGroup(n; variant=:left, kwargs...)
    SOn = SpecialOrthogonalGroup(n; kwargs...)
    Tn = TranslationGroup(n; kwargs...)
    variant ∉ [:left, :right] &&
        error("SE(n) requires a  variant ∉ [:left, :right] but you provided $variant")
    return variant === :left ? SOn ⋉ Tn : Tn ⋊ SOn
end

"""
    default_left_action(G::SpecialOrthogonalGroup, ::TranslationGroup)

Return the default left action for the special Euclidean group ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)``,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
default_left_action(::SpecialOrthogonalGroup, ::TranslationGroup) =
    LeftGroupOperationAction()

"""
    default_right_action(::TranslationGroup, G::SpecialOrthogonalGroup)

Return the default right action for the special Euclidean group,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
function default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)
    return LeftGroupOperationAction()
end

function ManifoldsBase.check_point(G::SpecialEuclideanGroup, p; kwargs...)
    n = get_parameter(G.manifold[1])
    errs = DomainError[]
    #
    if !isapprox(p[end, :], [zeros(size(p, 2) - 1)..., 1]; kwargs...)
        push!(
            errs,
            DomainError(
                p[end, :], "The last row of $p is not homogeneous, i.e. of form [0,..,0,1]."
            ),
        )
    end
    # SOn
    errS = check_point(G.manifold[1], get_submanifold_component(G, p, 1); kwargs...)
    !isnothing(errS) && push!(errs, errS)
    # translate part
    errT = check_point(G.manifold[2], get_submanifold_component(G, p, 2); kwargs...)
    !isnothing(errT) && push!(errs, errT)
    if length(errs) > 1
        return ManifoldsBase.CompositeManifoldError(errs)
    end
    return length(errs) == 0 ? nothing : first(errs)
end

function ManifoldsBase.check_vector(G::SpecialEuclideanGroup, X::AbstractMatrix; kwargs...)
    n = get_parameter(G.manifolds[1])
    errs = DomainError[]
    # homogeneous
    if !isapprox(X[end, :], zeros(size(X, 2)); kwargs...)
        push!(
            errs,
            DomainError(
                X[end, :],
                "The last row of $X is not homogeneous. Expected a zero-vector, got $(X[end, :]).",
            ),
        )
    end
    # SO(n)  part
    SOn = LieGroup(G.manifold[1], G.op[1])
    errS = check_vector(SOn, get_submanifold_component(G, X, 1); kwargs...)
    !isnothing(errS) && push!(errs, errS)
    # T(n) part
    Tn = LieGroup(G.manifold[2], G.op[2])
    errT = check_vector(Tn, get_submanifold_component(G, X, 2); kwargs...)
    !isnothing(errT) && push!(errs, errT)
    (length(errs) > 1) && (return CompositeManifoldError(errs))
    return length(errs) == 0 ? nothing : first(errs)
end

function compose!(
    ::SpecialEuclideanGroup, x::AbstractMatrix, p::AbstractMatrix, q::AbstractMatrix
)
    copyto!(x, p * q)
    return x
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{1}
)
    n = get_parameter(G.manifold[1])
    return view(p, 1:n, n + 1)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{2}
)
    n = _get_parameter(G)
    return view(p, 1:n, 1:n)
end

function ManifoldsBase.submanifold_components(G::SpecialOrthogonalGroup, p::AbstractMatrix)
    n = _get_parameter(G)
    @assert size(p) == (n + 1, n + 1)
    @inbounds t = submanifold_component(G, p, Val(1))
    @inbounds R = submanifold_component(G, p, Val(2))
    return (t, R)
end

function Base.show(io::IO, G::SpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[2].size)[1]
    return print(io, "SpecialEuclideanGroup($(size))")
end
