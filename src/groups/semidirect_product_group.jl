"""
    SemidirectProductGroupOperation{
        O1<:AbstractGroupOperation,
        O2<:AbstractGroupOperation,
        A<:AbstractGroupActionType,
        AO <: AbstractActionActsOnType
    } <: AbstractProductGroupOperation

An abstract type for all semidirect product group operations.

Most notably there are the left and right semidirect product group operations,
see [`LeftSemidirectProductGroupOperation`](@ref) and [`RightSemidirectProductGroupOperation`](@ref), respectively.
"""
abstract type SemidirectProductGroupOperation{
    O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
} <: AbstractProductGroupOperation end

"""
    LeftSemidirectProductGroupOperation{O1,O2,A,AO} <: SemidirectProductGroupOperation{O1,O2,A,AO}

A struct to model a left semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.

Then this group operation ``∘`` is defined on the product manifold ``$(_tex(:Cal, "G"))×$(_tex(:Cal, "H"))``
and uses the group operations ``⋆`` in the first component.
The second component depends on the choice of the actual [`AbstractGroupActionType`](@ref) `A`
and what it acts on, i.e. the [`AbstractActionActsOnType`](@ref) `AO`.

The resulting group operations are documented in the corresponding `compose` documentation.

For all four possible cases, we still use the shorthand notation ``$(_tex(:Cal, "G"))``[`⋉`](@ref)``$(_tex(:Cal, "H")) = ($(_tex(:Cal, "G"))×$(_tex(:Cal, "H")),∘)`` when it is clear which variant we refer to.
See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition for more details.

# Constructor

    LeftSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType,
        action_on::AbstractActionActsOnType=ActionActsOnLeft()
    )

# Parameters

* `op1::`[`AbstractGroupOperation`](@ref): The group operation ``⋆`` on ``$(_tex(:Cal, "G"))``
* `op2::`[`AbstractGroupOperation`](@ref): The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `action::`[`AbstractGroupActionType`](@ref): The group action ``α`` of ``$(_tex(:Cal, "G"))`` acting on ``$(_tex(:Cal, "H"))``.
* `action_on::`[`AbstractActionActsOnType`](@ref)`=`[`ActionActsOnLeft`](@ref)`()`: The type of element in ``$(_tex(:Cal, "H"))`` the action is applied to.

!!! note "A note on left/right"
    The “left” in the left semidirect product refers to the side, where the “pure” group operation takes place
    The “left/right” for the action refers to the type of group action used
    The “left/right” to act on refers to the left or right element in the second component, the action is applied to, see e.g. the explanation in [`ActionActsOnLeft`](@ref)
"""
struct LeftSemidirectProductGroupOperation{
        O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
    } <: SemidirectProductGroupOperation{O1, O2, A, AO}
    operations::Tuple{O1, O2}
    action_type::A
    action_on::AO
    function LeftSemidirectProductGroupOperation(
            op1::O1, op2::O2, action::A, action_on::AO = ActionActsOnLeft()
        ) where {
            O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
        }
        return new{O1, O2, A, AO}((op1, op2), action, action_on)
    end
end
@inline Base.getindex(spgo::SemidirectProductGroupOperation, i::Integer) =
    spgo.operations[i]

"""
    RightSemidirectProductGroupOperation{O1,O2,A} <: SemidirectProductGroupOperation{O1,O2,A}

A struct to model a right semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.


Then this group operation ``∘`` is defined on the product manifold ``$(_tex(:Cal, "H"))×$(_tex(:Cal, "G"))``
and uses the group operations ``⋆`` in the second component.
The first component depends on the choice of the actual [`AbstractGroupActionType`](@ref) `A`
and what it acts on, i.e. the [`AbstractActionActsOnType`](@ref) `AO`.

The resulting group operations are documented in the corresponding `compose` documentation.

For all four possible cases, we still use the shorthand notation ``$(_tex(:Cal, "H"))``[`⋊`](@ref)``$(_tex(:Cal, "G")) = ($(_tex(:Cal, "H"))×$(_tex(:Cal, "G")),∘)`` when it is clear which variant we refer to.
See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition for more details.

# Constructor

    RightSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
        action_on::AbstractActionActsOnType=ActionActsOnRight()
    )

# Parameters

* `op1::`[`AbstractGroupOperation`](@ref): The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `op2::`[`AbstractGroupOperation`](@ref): The group operation ``⋆`` on ``$(_tex(:Cal, "G"))``
* `action::`[`AbstractGroupActionType`](@ref): The group action ``α`` of ``$(_tex(:Cal, "G"))`` acting on ``$(_tex(:Cal, "H"))``.
* `action_on::`[`AbstractActionActsOnType`](@ref)`=`[`ActionActsOnRight`](@ref)`()`: The type of element in ``$(_tex(:Cal, "H"))`` the action is applied to.

!!! note "A note on left/right"
    The “right” in the right semidirect product refers to the side, where the “pure” group operation takes place
    The “left/right” for the action refers to the type of group action used
    The “left/right” to act on refers to the left or right element in the second component, the action is applied to, see e.g. the explanation in [`ActionActsOnLeft`](@ref)
"""
struct RightSemidirectProductGroupOperation{
        O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
    } <: SemidirectProductGroupOperation{O1, O2, A, AO}
    operations::Tuple{O1, O2}
    action_type::A
    action_on::AO
    function RightSemidirectProductGroupOperation(
            op1::O1, op2::O2, action::A, action_on::AO = ActionActsOnRight()
        ) where {
            O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
        }
        return new{O1, O2, A, AO}((op1, op2), action, action_on)
    end
end

"""
    LeftSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_left_action(N, H);
        action_on::AbstractActionActsOnType=ActionActsOnLeft()
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` for an [`AbstractGroupActionType`](@ref)
using the [`LeftSemidirectProductGroupOperation`](@ref) as group operation definition.
See [HilgertNeeb:2012; Definition 9.2.22](@cite), second definition, for more details.

The short form [`G ⋉ H`](@ref ⋉(L1::LieGroup, L2::LieGroup)) can be used if the
corresponding [`default_left_action(G,H)`](@ref default_left_action) as well as the [`ActionActsOnLeft`](@ref)
are the ones you want to use.
"""
function LeftSemidirectProductLieGroup(
        G::LieGroup, H::LieGroup, action::AbstractGroupActionType = default_left_action(G, H);
        action_on::AbstractActionActsOnType = ActionActsOnLeft()
    )
    # Use product manifold instead of × to not accidentally splat.
    return LieGroup(
        ProductManifold(G.manifold, H.manifold), LeftSemidirectProductGroupOperation(G.op, H.op, action, action_on)
    )
end

"""
    RightSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_right_action(N,H);
        action_on::AbstractActionActsOnType=ActionActsOnRight()
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` for an [`AbstractGroupActionType`](@ref)
using the [`RightSemidirectProductGroupOperation`](@ref) for the group operation definition.
See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition, for more details.

The short form [`H ⋊ G`](@ref ⋊(L1::LieGroup, L2::LieGroup)) can be used if the
corresponding [`default_right_action`](@ref)`(H,G)` and the [`ActionActsOnRight`](@ref)
are the ones you want to use.
"""
function RightSemidirectProductLieGroup(
        H::LieGroup, G::LieGroup, action::AbstractGroupActionType = default_right_action(H, G);
        action_on::AbstractActionActsOnType = ActionActsOnRight()
    )
    # Use product manifold instead of × to not accidentally splat.
    return LieGroup(
        ProductManifold(H.manifold, G.manifold), RightSemidirectProductGroupOperation(H.op, G.op, action, action_on)
    )
end

"""
    G ⋉ H
    ⋉(G, H)

For two [`LieGroups`](@ref) `G`, `H`, generate the [`LeftSemidirectProductLieGroup`](@ref)`(G, H)`,
where the corresponding [`default_left_action`](@ref)`(G, H)` and [`ActionActsOnRight`](@ref) are used.
"""
function ⋉(G::LieGroup, H::LieGroup)
    return LeftSemidirectProductLieGroup(G, H, default_left_action(G, H); action_on = ActionActsOnRight())
end

"""
    H ⋊ G
    ⋊(H, G)

For two [`LieGroups`](@ref) `H`, `G`, generate the [`RightSemidirectProductLieGroup`](@ref)`(H, G)`,
where the corresponding [`default_right_action`](@ref)`(H, G)` and [`ActionActsOnRight`](@ref) are used.
"""
function ⋊(H::LieGroup, G::LieGroup)
    return RightSemidirectProductLieGroup(H, G, default_right_action(H, G); action_on = ActionActsOnRight())
end

# A small helper to extract the product manifold, both Lie Groups and the action A
# It returns
# PM the product manifold
# G the first Lie group
# H the second Lie group that G acts on
# a the group action
# 1 the index of G in points of the semidirect Lie group/manifold
# 2 the index of H in points of the semidirect Lie group/manifold
function _semidirect_parts(SDPG::LieGroup{𝔽, <:LeftSemidirectProductGroupOperation, <:ProductManifold}) where {𝔽}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    a = GroupAction(G, H, SDPG.op.action_type)
    return PM, G, H, a, 1, 2
end
function _semidirect_parts(SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}) where {𝔽}
    PM = SDPG.manifold
    H, G = map(LieGroup, PM.manifolds, SDPG.op.operations)
    a = GroupAction(G, H, SDPG.op.action_type)
    return PM, G, H, a, 2, 1
end
# A major difference between left and right actions is that for right, we have to invert the action while for left we do not
function _semidirect_maybe_inv(::GroupAction{<:AbstractRightGroupActionType}, G, g)
    return inv(G, g)
end
function _semidirect_maybe_inv(::GroupAction{<:AbstractLeftGroupActionType}, G, g)
    return copy(G, g)
end
# and in in-place
function _semidirect_maybe_inv!(::GroupAction{<:AbstractRightGroupActionType}, G, k, g)
    return inv!(G, k, g)
end
function _semidirect_maybe_inv!(::GroupAction{<:AbstractLeftGroupActionType}, G, k, g)
    return copyto!(G, k, g)
end
#
#
# Functions
# ------------------------------------------------------------------------------------------
# For every function we to the following order of the 8 cases
# 1. Left semidirect, left action, act on left
# 2. Left semidirect, left action, act on right
# 3. Left semidirect, right action, act on left
# 4. Left semidirect, right action, act on right
# 5. Right semidirect, left action, act on left
# 6. Right semidirect, left action, act on right
# 7. Right semidirect, right action, act on left
# 8. Right semidirect, right action, act on right

_doc_semidirect_sub_groups = "Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.
"

# 1. Left semidirect, left action, act on left
# 5. Right semidirect, left action, act on left
"""
    compose(L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnLeft}}, g, h)

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the left.

The [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, σ_{g_2}(h_1) ⋄ h_2 $(_tex(:bigr))).
```

The [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( σ_{g_2}(h_1) ⋄ h_2, g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractLeftGroupActionType`](@ref) and [`ActionActsOnLeft`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnLeft}

# 3. Left semidirect, right action, act on left
# 7. Right semidirect, right action, act on left
"""
    compose(L::LieGroup{𝔽,SemidirectProductGroupOperation{⋄,⋆,<:AbstractRightGroupActionType,ActionActsOnLeft}}, g, h)

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the left.

The [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, τ_{g_2^{-1}}(h_1) ⋄ h_2 $(_tex(:bigr))).
```

The [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( τ_{g_2^{-1}}(h_1) ⋄ h_2, g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractRightGroupActionType`](@ref) and [`ActionActsOnLeft`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnLeft}

# An implementation for 1,3 (no inverse for left) and 5,7 (inverse for right)
function _compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, k, g, h
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnLeft}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    if Base.mightalias(k, g) || Base.mightalias(h, g) || Base.mightalias(k, h) # k/h may not overlap with g
        kH = copy(H, submanifold_component(SDPG, k, Val(h_ind))) # to
    else # if it does not alias, we can just use kG & kH
        kH = submanifold_component(SDPG, k, Val(h_ind))
    end
    if Base.mightalias(k, h) || Base.mightalias(k, g) # h/g may not overlap with k
        kG = copy(G, submanifold_component(SDPG, k, Val(g_ind))) # to
    else # if it does not alias, we can just use kG & kH
        kG = submanifold_component(SDPG, k, Val(g_ind))
    end
    # invert hG for right, copy for left
    # this is inplace if both are not aliased and creates a copy kG otherwise to avoid overwriting hG
    _semidirect_maybe_inv!(a, G, kG, submanifold_component(SDPG, h, Val(g_ind)))
    # a) group action  (first to avoid side effects in g, set kH to σ_{kG}(gH), with the above this avoids aliasing
    apply!(a, kH, kG, submanifold_component(SDPG, g, Val(h_ind)))
    # b) group operation on G
    _compose!(G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, h, Val(g_ind)))
    # c) group operation on H (note that the action happened in the de-aliased kH that hG or its inverse is already in kG
    _compose!(H, submanifold_component(SDPG, k, Val(h_ind)), kH, submanifold_component(SDPG, h, Val(h_ind)))
    return k
end

# 2. Left semidirect, left action, act on right
# 6. Right semidirect, left action, act on right
@doc """
    compose(L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h)

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the right.

The [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by


```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ σ_{g_1}(h_2) $(_tex(:bigr))).
```

The [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( h_1 ⋄ σ_{g_1}(h_2), g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractLeftGroupActionType`](@ref) and [`ActionActsOnRight`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnRight}

# 4 Left semidirect, right action, act on right
# 8. Right semidirect, right action, act on right
@doc """
    compose(L::LieGroup{𝔽,LeftSemidirectProductGroupOperation{⋆,⋄,<:AbstractRightGroupActionType,ActionActsOnRight}}, g, h)

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the right.

The [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by


```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2) $(_tex(:bigr))).
```

The [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( h_1 ⋄ τ_{g_1^{-1}}(h_2),  g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractRightGroupActionType`](@ref) and [`ActionActsOnRight`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnRight}

# a common implementation for 2,4 (left, no inverse) and 6,8 (right, with inverse)
function _compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, k, g, h
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnRight}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    if Base.mightalias(k, g) || Base.mightalias(h, g) # copy kH to avoid overlap aliased effects
        kH = copy(H, submanifold_component(SDPG, k, Val(h_ind))) # to
    else # if it does not alias, we can just use kG & kH
        kH = submanifold_component(SDPG, k, Val(h_ind))
    end
    if Base.mightalias(k, h) || Base.mightalias(k, g) # copy kH to avoid overlap aliased effects for both
        kG = copy(G, submanifold_component(SDPG, k, Val(g_ind))) # to
    else # if it does not alias, we can just use kG & kH
        kG = submanifold_component(SDPG, k, Val(g_ind))
    end
    # invert gG for right, copy for left
    # this is inplace if both are not aliased and creates a copy G otherwise to avoid overwriting hG
    _semidirect_maybe_inv!(a, G, kG, submanifold_component(SDPG, g, Val(g_ind)))
    # a) group action (first to avoid side effects in g, set kH to σ_{gG}(hH) - since we might have inverted, we have to use kG
    apply!(a, kH, kG, submanifold_component(PM, h, h_ind)) #accidentially overwriting hH is fine, we do not need it.
    # b) group operation on G
    _compose!(G, submanifold_component(PM, k, g_ind), submanifold_component(PM, g, g_ind), submanifold_component(PM, h, g_ind))
    # c) group operation on H (note that the action on hH is already in kH
    _compose!(H, submanifold_component(PM, k, h_ind), submanifold_component(PM, g, h_ind), kH)
    return k
end

_semidirect_diff_compose_notation = """
## Notation

* ``Dρ_{g'}(g)[X]`` and ``Dρ_{h'}(h)[Y]`` denote [`diff_left_compose`](@ref) on ``$(_tex(:Cal, "G"))`` and ``$(_tex(:Cal, "H"))``, respectively.
* ``Dλ_{g'}(g)[X]`` and ``Dλ_{h'}(h)[Y]`` denote the [`diff_right_compose`](@ref) on ``$(_tex(:Cal, "G"))`` and  ``$(_tex(:Cal, "H"))``, respectively.
* ``Dσ_g(h)[Y]`` denotes the [`diff_apply`](@ref)
* ``D_{$(_tex(:Cal, "G"))}σ_g(h)[X]`` denotes the [`diff_group_apply`](@ref), i.e., differential of the group action with respect to its base point ``g``.
"""

# diff compose left

# 1. Left semidirect, left action, act on left
# 5. Right semidirect, left action, act on left
"""
    diff_left_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnLeft}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the left argument.
This means that we consider the right group operation man ``ρ_h(g) = g ∘ h`` with respect to ``g``
or the left group operation ``λ_g(h) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the left.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    ρ_{(g_2,h_2)}(g_1,h_1) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, σ_{g_2}(h_1) ⋄ h_2).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D ρ_{(g_2,h_2)}(g_1,h_1)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{g_2}(g_1)[X],
  Dρ_{h_2}(σ_{g_2}(h_1)[Dσ_{g_2}(h_1)[Y]]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    ρ_{(h_2,g_2)}(h_1,g_1) := (h_1,g_1) ∘ (h_2,g_2) = (σ_{g_2}(h_1) ⋄ h_2, g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{h_2}(σ_{g_2}(h_1)[Dσ_{g_2}(h_1)[Y]],
  Dρ_{g_2}(g_1)[X],
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_left_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnLeft}

# 3. Left semidirect, right action, act on left
# 7. Right semidirect, right action, act on left
"""
    diff_left_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the left argument.
This means that we consider the right group operation man ``ρ_h(g) = g ∘ h`` with respect to ``g``
or the left group operation ``λ_g(h) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the left.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    ρ_{(g_2,h_2)}(g_1,h_1) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, τ_{g_2^{-1}}(h_1) ⋄ h_2).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D ρ_{(g_2,h_2)}(g_1,h_1)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{g_2}(g_1)[X],
  Dρ_{h_2}(τ_{g_2^{-1}}(h_1)[Dτ_{g_2^{-1}}(h_1)[Y]]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    ρ_{(h_2,g_2)}(h_1,g_1) := (h_1,g_1) ∘ (h_2,g_2) = (τ_{g_2^{-1}}(h_1) ⋄ h_2, g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{h_2}(τ_{g_2^{-1}}(h_1)[Dτ_{g_2^{-1}}(h_1)[Y]],
  Dρ_{g_2}(g_1)[X]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_left_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnLeft}

# An implementation for 1,3 (no inverse for left) and 5,7 (inverse for right)
# in the computation of the element coming from the group action.
function diff_left_compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, Y, g, h, X
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnLeft}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    sdpg = LieAlgebra(SDPG)
    # We use a naming close to the left variant: g = (gG, gH) = (g_1,h_1)
    YG, YH = submanifold_component(sdpg, Y, Val(g_ind)), submanifold_component(sdpg, Y, Val(h_ind))
    XG, XH = submanifold_component(sdpg, X, Val(g_ind)), submanifold_component(sdpg, X, Val(h_ind))
    gG, gH = submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, g, Val(h_ind))
    hG, hH = submanifold_component(SDPG, h, Val(g_ind)), submanifold_component(SDPG, h, Val(h_ind))
    # For right actions we have to invert hG - this allocates when it has to invert
    hG_mod = _semidirect_maybe_inv(a, G, hG)
    # we need one allocation to compute the action
    σg2h1 = apply(a, hG_mod, gH)
    # Step 1: Compute the argument for the second components diff
    diff_apply!(a, YH, hG_mod, gH, XH)
    # Step 2: Differential of right group compose (argument from 1)
    diff_left_compose!(H, YH, σg2h1, hH, YH)
    # last: the plain diff compose on G
    diff_left_compose!(G, YG, gG, hG, XG)
    return Y
end

# 2. Left semidirect, left action, act on right
# 6. Right semidirect, left action, act on right
"""
    diff_left_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the left argument.
This means that we consider the right group operation man ``ρ_h(g) = g ∘ h`` with respect to ``g``
or the left group operation ``λ_g(h) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the right.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    ρ_{(g_2,h_2)}(g_1,h_1) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, h_1 ⋄ σ_{g_1}(h_2)).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D ρ_{(g_2,h_2)}(g_1,h_1)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{g_2}(g_1)[X],
  Dρ_{σ_{g_1}(h_2)}(h_1)[Y] + Dλ_{h_1}(σ_{g_1}(h_2)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}σ_{g_1}(h_2)[X] $(_tex(:bigr))]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    ρ_{(h_2,g_2)}(h_1,g_1) := (h_1,g_1) ∘ (h_2,g_2) = (h_1 ⋄ σ_{g_2}(h_2), g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{σ_{g_1}(h_2)}(h_1)[Y] + Dλ_{h_1}(σ_{g_1}(h_2)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}σ_{g_1}(h_2)[X] $(_tex(:bigr))],
  Dρ_{g_2}(g_1)[X]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_left_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnRight}

# 4 Left semidirect, right action, act on right
# 8. Right semidirect, right action, act on right
"""
    diff_left_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractRightGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the left argument.
This means that we consider the right group operation man ``ρ_h(g) = g ∘ h`` with respect to ``g``
or the left group operation ``λ_g(h) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the right.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    ρ_{(g_2,h_2)}(g_1,h_1) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2)).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D ρ_{(g_2,h_2)}(g_1,h_1)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{g_2}(g_1)[X],
  Dρ_{τ_{g_1^{-1}}}(h_1)[Y] + Dλ_{h_1}(τ_{g_1^{-1}}(h_2)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}τ_{g_1^{-1}}(h_2)[X] $(_tex(:bigr))]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    ρ_{(h_2,g_2)}(h_1,g_1) := (h_1,g_1) ∘ (h_2,g_2) = (h_1 ⋄ σ_{g_2}(h_2), g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dρ_{τ_{g_1^{-1}}}(h_1)[Y] + Dλ_{h_1}(τ_{g_1^{-1}}(h_2)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}τ_{g_1^{-1}}(h_2)[X] $(_tex(:bigr))],
  Dρ_{g_2}(g_1)[X]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_left_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnRight}

# An implementation for 2,6 (no inverse for left) and 4,8 (inverse for right)
# in the computation of the element coming from the group action.
function diff_left_compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, Y, g, h, X
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnRight}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    sdpg = LieAlgebra(SDPG)
    # We use a naming close to the left variant:
    # g = (gG, gH) = (g_1,h_1) h = (hG, hH) = (g_2,h_2)
    # X = (XG, XH) = (X, Y), Y = (YG, YH) the two components of the result
    YG, YH = submanifold_component(sdpg, Y, Val(g_ind)), submanifold_component(sdpg, Y, Val(h_ind))
    XG, XH = submanifold_component(sdpg, X, Val(g_ind)), submanifold_component(sdpg, X, Val(h_ind))
    gG, gH = submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, g, Val(h_ind))
    hG, hH = submanifold_component(SDPG, h, Val(g_ind)), submanifold_component(SDPG, h, Val(h_ind))
    # avoid aliasing
    _YH = Base.mightalias(YH, XH) ? copy(H, hH, YH) : YH
    # For right actions we have to invert gG - this allocates when it has to invert
    gG_mod = _semidirect_maybe_inv(a, G, gG)
    # one allocation for applying the action
    # we need one allocation to compute the action
    σg1mh2 = apply(a, gG_mod, hH)
    # Step 1: Compute the argument for the second summand on H - use the memory of YH for the result
    diff_group_apply!(a, _YH, gG_mod, hH, XG)
    # Step 2: Differential of right group compose (argument from 1)
    diff_right_compose!(H, _YH, gH, σg1mh2, _YH)
    # Step 3: a second allocation for the other (first) differential, we want to add to that.
    YH .= _YH + diff_left_compose(H, gH, σg1mh2, XH)
    # last: the plain diff compose on G
    diff_left_compose!(G, YG, gG, hG, XG)
    return Y
end

# diff_right_compose

# 1. Left semidirect, left action, act on left
# 5. Right semidirect, left action, act on left
"""
    diff_right_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnLeft}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the right argument.
This means that we consider the left group operation man ``λ_g(h) = g ∘ h`` with respect to ``h``
or the right group operation ``ρ_h(g) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the left.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    λ_{(g_1,h_1)}(g_2,h_2) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, σ_{g_2}(h_1) ⋄ h_2).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D λ_{(g_1,h_1)}((g_2,h_2))$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{σ_{g_2}(h_1)}(h_2)[Y] + Dρ_{h_2}(σ_{g_2}(h_1)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}σ_{g_2}(h_1)[X] $(_tex(:bigr))]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    λ_{(h_1,g_1)}(h_2,g_2) := (h_1,g_1) ∘ (h_2,g_2) = (σ_{g_2}(h_1) ⋄ h_2, g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{σ_{g_2}(h_1)}(h_2)[Y] + Dρ_{h_2}(σ_{g_2}(h_1)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}σ_{g_2}(h_1)[X] $(_tex(:bigr))]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_right_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnLeft}

# 3. Left semidirect, right action, act on left
# 7. Right semidirect, right action, act on left
"""
    diff_right_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the right argument.
This means that we consider the left group operation man ``λ_g(h) = g ∘ h`` with respect to ``h``
or the right group operation ``ρ_h(g) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the left.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    λ_{(g_1,h_1)}((g_2,h_2)) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, τ_{g_2^{-1}}(h_1) ⋄ h_2).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D λ_{(g_1,h_1)}(g_2,h_2)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{τ_{g_2^{-1}}(h_1)}(h_2)[Y] + Dρ_{h_2}(τ_{g_2^{-1}}(h_1)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}τ_{g_2^{-1}}(h_1)[X] $(_tex(:bigr))]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    λ_{(h_1,g_1)}(h_2,g_2) := (h_1,g_1) ∘ (h_2,g_2) = (τ_{g_2^{-1}}(h_1) ⋄ h_2, g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D λ_{(h_1,g_1)}(h_2,g_2)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{τ_{g_2^{-1}}(h_1)}(h_2)[Y] + Dρ_{h_2}(τ_{g_2^{-1}}(h_1)) $(_tex(:bigl))[ D_{$(_tex(:Cal, "G"))}τ_{g_2^{-1}}(h_1)[X] $(_tex(:bigr))],
  Dλ_{g_1}(g_2)[X],
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_right_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnLeft}

# An implementation for 1,3 (no inverse for left) and 5,7 (inverse for right)
# in the computation of the element coming from the group action.
function diff_right_compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, Y, g, h, X
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnLeft}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    sdpg = LieAlgebra(SDPG)
    # We use a naming close to the left variant: g = (gG, gH) = (g_1,h_1)
    YG, YH = submanifold_component(sdpg, Y, Val(g_ind)), submanifold_component(sdpg, Y, Val(h_ind))
    XG, XH = submanifold_component(sdpg, X, Val(g_ind)), submanifold_component(sdpg, X, Val(h_ind))
    gG, gH = submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, g, Val(h_ind))
    hG, hH = submanifold_component(SDPG, h, Val(g_ind)), submanifold_component(SDPG, h, Val(h_ind))

    _YH = Base.mightalias(YH, XH) ? copy(H, hH, YH) : YH
    # For right actions we have to invert hG (g_2) - this allocates when it has to invert
    hG_mod = _semidirect_maybe_inv(a, G, hG)
    # one allocation for applying the action
    # we need one allocation to compute the action (g_2 to h_1)
    σg2mh1 = apply(a, hG_mod, gH)
    # Step 1: Compute the argument for the second summand on H - use the memory of YH for the result
    diff_group_apply!(a, _YH, hG_mod, gH, XG)
    # Step 2: Differential of right group compose (argument from 1)
    diff_left_compose!(H, _YH, σg2mh1, hH, _YH)
    # Step 3: a second allocation for the other (first) differential, we want to add to that.
    YH .= _YH + diff_right_compose(H, σg2mh1, hH, XH)
    # last: the plain diff compose on G
    diff_right_compose!(G, YG, gG, hG, XG)
    return Y
end

# 2. Left semidirect, left action, act on right
# 6. Right semidirect, left action, act on right
"""
    diff_right_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the right argument.
This means that we consider the left group operation man ``λ_g(h) = g ∘ h`` with respect to ``h``
or the right group operation ``ρ_h(g) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the right.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    λ_{(g_1,h_1)}(g_2,h_2) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, h_1 ⋄ σ_{g_1}(h_2)).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D λ_{(g_1,h_1)}(g_2,h_2)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{h_1}(σ_{g_1}(h_2))[Dσ_{g_1}(h_2)[Y]]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    λ_{(h_1,g_1)}(h_2,g_2) := (h_1,g_1) ∘ (h_2,g_2) = (h_1 ⋄ σ_{g_2}(h_2), g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{h_1}(σ_{g_1}(h_2))[Dσ_{g_1}(h_2)[Y]]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_right_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnRight}

# 4 Left semidirect, right action, act on right
# 8. Right semidirect, right action, act on right
"""
    diff_right_compose(
        L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h, X
    )

Compute the differential of the group operation ``⋅∘⋅`` with respect to the right argument.
This means that we consider the left group operation man ``λ_g(h) = g ∘ h`` with respect to ``h``
or the right group operation ``ρ_h(g) = g ∘ h`` with respect to its index ``g``.

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the right.

For the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``$(_tex(:Cal, "G")) ⋉ $(_tex(:Cal, "H"))`` we
have

```math
    λ_{(g_1,h_1)}(g_2,h_2) := (g_1,h_1) ∘ (g_2,h_2) = (g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2)).
```

such that their differential reads for some ``(X, Y)`` from the Lie algebra that

```math
D λ_{(g_1,h_1)}(g_2,h_2)$(_tex(:bigl))[ (X,Y) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{g_1}(g_2)[X],
  Dλ_{h_1}(τ_{g_1^{-1}}(h_2))[Dτ_{g_1^{-1}}(h_2)[Y]]
  $(_tex(:bigr))).
```

For the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on `` $(_tex(:Cal, "H")) ⋊ $(_tex(:Cal, "G"))`` we
have

```math
    ρ_{(h_2,g_2)}(h_1,g_1) := (h_1,g_1) ∘ (h_2,g_2) = (h_1 ⋄ σ_{g_2}(h_2), g_1 ⋆ g_2).
```

such that their differential reads for some ``(Y, X)`` from the Lie algebra that

```math
D ρ_{(h_2,g_2)}(h_1,g_1)$(_tex(:bigl))[ (Y, X) $(_tex(:bigr))]
  = $(_tex(:bigl))(
  Dλ_{h_1}(τ_{g_1^{-1}}(h_2))[Dτ_{g_1^{-1}}(h_2)[Y]],
  Dλ_{g_1}(g_2)[X]
  $(_tex(:bigr))).
```

$(_semidirect_diff_compose_notation)
"""
diff_right_compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, g, h, X
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: ActionActsOnRight}

# An implementation for 2,6 (no inverse for left) and 4,8 (inverse for right)
# in the computation of the element coming from the group action.
function diff_right_compose!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, Y, g, h, X
    ) where {𝔽, O1, O2, A <: AbstractGroupActionType, AO <: ActionActsOnRight}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    sdpg = LieAlgebra(SDPG)
    # We use a naming close to the left variant:
    # g = (gG, gH) = (g_1,h_1) h = (hG, hH) = (g_2,h_2)
    # X = (XG, XH) = (X, Y), Y = (YG, YH) the two components of the result
    YG, YH = submanifold_component(sdpg, Y, Val(g_ind)), submanifold_component(sdpg, Y, Val(h_ind))
    XG, XH = submanifold_component(sdpg, X, Val(g_ind)), submanifold_component(sdpg, X, Val(h_ind))
    gG, gH = submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, g, Val(h_ind))
    hG, hH = submanifold_component(SDPG, h, Val(g_ind)), submanifold_component(SDPG, h, Val(h_ind))
    # For right actions we have to invert gG - this allocates when it has to invert
    gG_mod = _semidirect_maybe_inv(a, G, gG)
    # we need one allocation to compute the action
    σg1h2 = apply(a, gG_mod, hH)
    # Step 1: Compute the argument for the second components diff
    diff_apply!(a, YH, gG_mod, hH, XH)
    # Step 2: Differential of left group compose (argument from 1)
    diff_right_compose!(H, YH, gH, σg1h2, YH)
    # last: the plain diff compose on G
    diff_right_compose!(G, YG, gG, hG, XG)
    return Y
end

# get coordinates

function get_coordinates_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, c, X, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: SemidirectProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, Prc, submanifold_components(Pr𝔤, X))
    map(ts) do t
        return get_coordinates_lie!(t..., B)
    end
    return c
end

# get vector

function get_vector_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, c, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: SemidirectProductGroupOperation, M <: ProductManifold}
    PrG = Pr𝔤.manifold
    PrM = PrG.manifold
    dims = map(manifold_dimension, PrM.manifolds)
    @assert length(c) == sum(dims)
    dim_ranges = ManifoldsBase._get_dim_ranges(dims)
    Prc = map(dr -> (@inbounds view(c, dr)), dim_ranges)
    PrL = LieAlgebra.(LieGroup.(PrM.manifolds, PrG.op.operations))
    ts = ManifoldsBase.ziptuples(PrL, submanifold_components(Pr𝔤, X), Prc)
    map(ts) do t
        return get_vector_lie!(t..., B)
    end
    return X
end

# 1. Left semidirect, left action, act on left
# 2. Left semidirect, left action, act on right
# 5. Right semidirect, left action, act on left
# 6. Right semidirect, left action, act on right
"""
    inv(L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,AbstractActionActsOn}}, g)

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. The inverse is the same for both [`ActionActsOnLeft`](@ref) and [`ActionActsOnRight`](@ref).

The inverse for the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by

```math
(g,h)^{-1} = $(_tex(:bigl))( g^{-1}, σ_{g^{-1}}(h^{-1}) $(_tex(:bigr))).
```

The inverse for the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋊ H`` is given by

```math
(h,g)^{-1} = $(_tex(:bigl))( σ_{g^{-1}} (h^{-1}), g^{-1} $(_tex(:bigr)))
```
"""
inv(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any,
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: AbstractActionActsOnType}

function _inv!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, k, g
    ) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: AbstractActionActsOnType}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    # (a) compute the inverse in G
    inv!(G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(PM, g, Val(g_ind)))
    # (b) compute the inverse in H
    inv!(H, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(PM, g, Val(h_ind)))
    # (c) apply the group action w.r.t. the inverse in G (from a) to the inverse from (b)
    apply!( # Apply the group action with g1^-1 to g2^-1 - works with aliases if apply does
        a,
        submanifold_component(SDPG, k, Val(h_ind)),
        submanifold_component(SDPG, k, Val(g_ind)),
        submanifold_component(SDPG, k, Val(h_ind)),
    )
    return k
end

# 3. Left semidirect, right action, act on left
# 4. Left semidirect, right action, act on right
# 8. Right semidirect, right action, act on right
# 7. Right semidirect, right action, act on left
"""
    inv(L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractRightGroupActionType,AbstractActionActsOn}}, g)

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. The inverse is the same for both [`ActionActsOnLeft`](@ref) and [`ActionActsOnRight`](@ref).

The inverse for the [`LeftSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋉ H`` is given by

```math
(g,h)^{-1} = $(_tex(:bigl))( g^{-1}, τ_{g}(h^{-1}) $(_tex(:bigr))).
```

The inverse for the [`RightSemidirectProductGroupOperation`](@ref) ``$(_math(:∘))`` on ``G ⋊ H`` is given by

```math
(h,g)^{-1} = $(_tex(:bigl))( τ_{g}(h^{-1}), g^{-1} $(_tex(:bigr)))
```
"""
inv(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any,
) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: AbstractActionActsOnType}

function _inv!(
        SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, k, g
    ) where {𝔽, O1, O2, A <: AbstractRightGroupActionType, AO <: AbstractActionActsOnType}
    PM, G, H, a, g_ind, h_ind = _semidirect_parts(SDPG)
    # (a) compute the inverse in H
    inv!(H, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(PM, g, Val(h_ind)))
    # (b) apply the group action w.r.t. g_g in G (not yet inverted in place) to the inverse g_h from (a)) that is already in k_h; apply it in-place of k_h
    apply!( # Apply the group action with g1^-1 to g2^-1
        a,
        submanifold_component(SDPG, k, Val(h_ind)),
        submanifold_component(SDPG, g, Val(g_ind)),
        submanifold_component(SDPG, k, Val(h_ind)),
    )
    # (c) compute the inverse in G
    inv!(G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(PM, g, Val(g_ind)))
    return k
end

function identity_element!(
        SDPG::LieGroup{𝔽, Op, M}, e
    ) where {𝔽, Op <: SemidirectProductGroupOperation, M <: ProductManifold}
    GH = map(LieGroup, SDPG.manifold.manifolds, SDPG.op.operations)
    identity_element!.(GH, submanifold_components(SDPG.manifold, e))
    return e
end

function Base.show(
        io::IO, SDPG::LieGroup{𝔽, <:LeftSemidirectProductGroupOperation, <:ProductManifold}
    ) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = SDPG.op.action_type
    return print(io, "LeftSemidirectProductLieGroup($G, $H, $at)")
end
function Base.show(
        io::IO, SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}
    ) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = SDPG.op.action_type
    return print(io, "RightSemidirectProductLieGroup($G, $H, $at)")
end
