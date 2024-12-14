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
    SpecialEuclideanGroup(n; kwargs...)
    SpecialOrthogonalGroup(n; kwargs...) ⋉ TranslationGroup(n; kwargs...)

Generate special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``, where the first
constructor is equivalent to the second.

Alternatively you  can also use

    TranslationGroup(n; kwargs...) ⋊ SpecialOrthogonalGroup(n; kwargs...)

to define ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)``.

If you prefer to have the order of the elements reversed, i.e. a representation with first
the translation vector and then the rotation matrix, use the third constructor.

All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialEuclideanGroup{T} = LieGroup{
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

"""
    default_left_action(G::SpecialOrthogonalGroup, ::TranslationGroup)

Return the default left action for the special Euclidean group ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)``,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
default_left_action(G::SpecialOrthogonalGroup, ::TranslationGroup) =
    LeftGroupOperationAction()

"""
    default_right_action(::TranslationGroup, G::SpecialOrthogonalGroup)

Return the default right action for the special Euclidean group,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
function default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)
    return LeftGroupOperationAction()
end

function SpecialEuclideanGroup(n; kwargs...)
    SOn = SpecialOrthogonalGroup(n; kwargs...)
    Tn = TranslationGroup(n; kwargs...)
    return SOn ⋉ Tn
end

function Base.show(io::IO, G::SpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[2].size)[1]
    return print(io, "SpecialEuclideanGroup($(size))")
end
