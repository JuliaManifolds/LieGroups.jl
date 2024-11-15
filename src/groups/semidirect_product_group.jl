
#
#
# Semidirect product groups – model semidirect products of rwo Lie groups
#
"""
    SemidirectProductGroupOperation{O1,O2,A} <: AbstractGroupOperation

A struct to model a semidirect Lie group product.

Let ``($(_tex(:Cal, "N")), ⋄)`` and ``($(_tex(:Cal, "H")), ⋆)`` be two Lie groups
with group operations ``⋄`` and ``⋆``, respectively, as well as
a (left) group action
``σ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``, cf [`AbstractLeftGroupActionType`](#ref).
We use here as well use the notation ``σ_h: $(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))`` as a family of maps on ``$(_tex(:Cal, "N"))``

Then we define a group operation ``∘`` on the product manifold $(_tex(:Cal, "N"))×$(_tex(:Cal, "H")) by

```math
    (n_1,h_1) ∘ (n_2,h_2) := (n_1 ⋄ σ_{h_1}(n_2), h_1 ⋆ h_2)
```

and similarly for a right group action, for example ``τ_h = (σ_h)^{-1}``
    a group operation on $(_tex(:Cal, "H"))×$(_tex(:Cal, "N"))

```math
    (h_1,n_1) ∘ (h_2,n_2) := (h_1 ⋆ h_2, τ_{h_2}(n_1) ⋄ n_1).
```

See [HilgertNeeb:2012; Definition 9.2.22](@cite) for more details.

# Constructor

    SemidirectProductGroupOperation(op1::AbstractGroupOperation, op2::AbstractGroupOperation, action)

Create a `SemidirectProductGroupOperation` for [`AbstractGroupOperation`](@ref) `op1=```⋄`` and `pü2=```⋆``
together with an [`AbstractGroupActionType`](@ref) `action` that defaults to the [`LeftGroupOperationAction`](@ref), which yields the corresponding [`GroupOperationAction`](@ref),
that is where the action is given by the group operation of ``$(_tex(:Cal, "H"))`` acting on ``$(_tex(:Cal, "N"))``.
The first shorthand `op1 ⋊ op2` is a short form for this.

For the second form, use for example the [`RightGroupOperationAction`](@ref) to also have the [`GroupOperationAction`](@ref) mentioned.
"""
struct SemidirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: AbstractGroupOperation
    op1::O1
    op2::O2
    action_type::A
    function SemidirectProductGroupOperation(
        op1::O1, op2::O2; action::A=LeftGroupOperationAction()
    ) where {
        O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
    }
        return SemidirectProductGroupOperation{O1,O2,A}(op1, op2, action)
    end
end


# LeftSemidirectProductLieGroup
# RightSemidirectProductLieGroup
# default_left_action(L1,L2)
# default_right_action(L1,L2)
"""
    SemidirectProductLieGroup(N::LieGroup, H::LieGroup; action=LeftGroupOperationAction())

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) = N ⋊ H`` for an [`AbstractLeftGroupActionType`](@ref),
and ``$(_tex(:Cal, "G")) = N ⋉ H`` for an [`AbstractRightGroupActionType`](@ref),

see [`SemidirectProductGroupOperation`](@ref) for the group operation definition as well as [HilgertNeeb:2012; Definition 9.2.22](@cite) for more details.
"""
function SemidirectProductLieGroup(N::LieGroup, H::LieGroup; action=LeftGroupOperationAction())
    return LieGroup(
        N.manifold × H.manifold, SemidirectProductGroupOperation(N.op, H.op; action=action)
    )
end

"""
    L1 ⋊ L2

For two [`LieGroups`](@ref), generate the [`SemidirectProductLieGroup`](@ref)`(L1, L2; action=`[`LeftGroupOperationAction`](@ref)`())`
"""
function ⋊(L1::LieGroup, L2::LieGroup)
    return SemidirectProductLieGroup(L1, L2; action=LeftGroupOperationAction())
end

"""
    L1 ⋉ L2

For two [`LieGroups`](@ref), generate the [`SemidirectProductLieGroup`](@ref)`(L1, L2; action=`[`RightGroupOperationAction`](@ref)`())`
"""
function ⋉(L1::LieGroup, L2::LieGroup)
    return SemidirectProductLieGroup(L1, L2; action=RightGroupOperationAction())
end
