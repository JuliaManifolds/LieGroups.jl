"""
    SemiDirectProductGroupOperation{
        O1<:AbstractGroupOperation,
        O2<:AbstractGroupOperation,
        A<:AbstractGroupActionType,
        AO <: AbstractActionActsOnType
    } <: AbstractProductGroupOperation

An abstract type for all semidirect product group operations.

Most notably there are the left and right semidirect product group operations,
see [`LeftSemidirectProductGroupOperation`](@ref) and [`RightSemidirectProductGroupOperation`](@ref), respectively.
"""
abstract type SemiDirectProductGroupOperation{
    O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType, AO <: AbstractActionActsOnType,
} <: AbstractProductGroupOperation end

"""
    LeftSemidirectProductGroupOperation{O1,O2,A,AO} <: SemiDirectProductGroupOperation{O1,O2,A,AO}

A struct to model a left semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.

Then this group operation ``∘`` is defined on the product manifold ``$(_tex(:Cal, "G"))×$(_tex(:Cal, "H"))``
and uses the group operations ``⋆`` in the first component.
The second component depends on the choice of the actual [``AbstractGroupActionType`](@ref) `A`
and what it acts on, i.e. the [`AbstractActionActsOnType`](@ref) `AO`.

The resulting group operations are documented in the corresponding `compose` documentation.

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
    } <: SemiDirectProductGroupOperation{O1, O2, A, AO}
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
@inline Base.getindex(spgo::SemiDirectProductGroupOperation, i::Integer) =
    spgo.operations[i]

"""
    RightSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a right semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.


Then this group operation ``∘`` is defined on the product manifold ``$(_tex(:Cal, "H"))×$(_tex(:Cal, "G"))``
and uses the group operations ``⋆`` in the second component.
The first component depends on the choice of the actual [``AbstractGroupActionType`](@ref) `A`
and what it acts on, i.e. the [`AbstractActionActsOnType`](@ref) `AO`.

The resulting group operations are documented in the corresponding `compose` documentation.

For both we use the shorthand notation ``H``[`⋊`](@ref)``G = (H×G,∘)``.
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
    } <: SemiDirectProductGroupOperation{O1, O2, A, AO}
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
    L1 ⋉ L2
    ⋉(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`LeftSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_left_action`](@ref)`(L1, L2)` and [`ActionActsOnLeft`](@ref) are used.
"""
function ⋉(L1::LieGroup, L2::LieGroup)
    return LeftSemidirectProductLieGroup(L1, L2, default_left_action(L1, L2); action_on = ActionActsOnLeft())
end

"""
    L1 ⋊ L2
    ⋊(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`RightSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_right_action`](@ref)`(L1, L2)` and [`ActionActsOnRight`](@ref) are used.
"""
function ⋊(L1::LieGroup, L2::LieGroup)
    return RightSemidirectProductLieGroup(L1, L2, default_right_action(L1, L2); action_on = ActionActsOnRight())
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
    a = GroupAction(SDPG.op.action_type, G, H)
    return PM, G, H, a, 1, 2
end
function _semidirect_parts(SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}) where {𝔽}
    PM = SDPG.manifold
    H, G = map(LieGroup, PM.manifolds, SDPG.op.operations)
    a = GroupAction(SDPG.op.action_type, G, H)
    return PM, G, H, a, 2, 1
end
# A major difference betrween left and right actions is that for right, we have to invert the action while for left we do not
function _semidirect_maybe_inv!(a::A, G, k, g) where {A <: AbstractRightGroupActionType}
    return inv!(G, k, g)
end
function _semidirect_maybe_inv!(a::A, G, k, g) where {A <: AbstractLeftGroupActionType}
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
    compose(L::LieGroup{𝔽,LeftSemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnLeft}}, g, h)

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the left.

The group operation ``$(_math(:∘))`` on the [`LeftSemidirectProductGroup`](@ref) ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, σ_{g_2}(h_1) ⋄ h_2 $(_tex(:bigr))).
```

The group operation ``$(_math(:∘))`` on the [`RightSemidirectProductGroup`](@ref) ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( σ_{g_2}(h_1) ⋄ h_2, g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractLeftGroupActionType`](@ref) and [`ActionActsOnLeft`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnLeft}

# 3. Leftsemidirect, right action, act on left
# 7. Rightsemidirect, right action, act on left
"""
    compose(L::LieGroup{𝔽,SemidirectProductGroupOperation{⋄,⋆,<:AbstractRightGroupActionType,ActionActsOnLeft}}, g, h)

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the left.

The group operation ``$(_math(:∘))`` on the [`LeftSemidirectProductGroup`](@ref) ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, τ_{g_2^{-1}}(h_1) ⋄ h_2 $(_tex(:bigr))).
```

The group operation ``$(_math(:∘))`` on the [`RightSemidirectProductGroup`](@ref) ``H ⋊ G`` is given by

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
    _semidirect_maybe_inv!(a, G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(SDPG, h, Val(g_ind)))
    # Note difference from math to code: g=(g_1,h_1), h=(g_2,h_2)=(h1,)
    # a) group action (first to avoid side effects in g, set k_2 to σ_{g_2}(h_1)
    apply!(a, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(SDPG, h, Val(g_ind)), submanifold_component(SDPG, g, Val(k_ind)))
    # b) group operation on G
    _compose!(G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, h, Val(g_ind)))
    # c) group operation on H (not that h_1 is already in k2)
    _compose!(H, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(SDPG, k, Val(k_ind)), submanifold_component(SDPG, h, Val(k_ind)))
    return k
end

# 2. Left semidirect, left action, act on right
# 6. Right semidirect, left action, act on right
@doc """
    compose(L::LieGroup{𝔽,<:SemidirectProductGroupOperation{⋆,⋄,<:AbstractLeftGroupActionType,ActionActsOnRight}}, g, h)

$(_doc_semidirect_sub_groups) Let ``σ`` denote a left group action. It here acts on the right.

The group operation ``$(_math(:∘))`` on the [`LeftSemidirectProductGroup`](@ref) ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ σ_{g_1}(h_2) $(_tex(:bigr))).
```

The group operation ``$(_math(:∘))`` on the [`RightSemidirectProductGroup`](@ref) ``H ⋊ G`` is given by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( h_1 ⋄ σ_{g_1}(h_2), g_1 ⋆ g_2 $(_tex(:bigr))).
```

See also [`AbstractLeftGroupActionType`](@ref) and [`ActionActsOnRight`](@ref).
"""
compose(
    SDPG::LieGroup{𝔽, <:SemidirectProductGroupOperation{O1, O2, A, AO}, <:ProductManifold}, ::Any, ::Any
) where {𝔽, O1, O2, A <: AbstractLeftGroupActionType, AO <: ActionActsOnRight}

# 4 Leftsemidirect, right action, act on right
# 8. Rightsemidirect, right action, act on right
@doc """
    compose(L::LieGroup{𝔽,LeftSemidirectProductGroupOperation{⋆,⋄,<:AbstractRightGroupActionType,ActionActsOnRight}}, g, h)

$(_doc_semidirect_sub_groups) Let ``τ`` denote a right group action. It here acts on the right.
The group operation ``$(_math(:∘))`` on the [`LeftSemidirectProductGroup`](@ref) ``G ⋉ H`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2) $(_tex(:bigr))).
```

The group operation ``$(_math(:∘))`` on the [`RightSemidirectProductGroup`](@ref) ``H ⋊ G`` is given by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2) $(_tex(:bigr))).
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
    # invert for right, copy for left
    _semidirect_maybe_inv!(a, G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(SDPG, g, Val(g_ind)))
    # Note difference from math to code: g=(g_1,h_1), h=(g_2,h_2)=(h1,h2)
    # a) group action (first to avoid side effects in g, set k_2 to σ_{g_1}(h_2)
    apply!(a, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, h, Val(h_ind)))
    # b) group operation on G
    _compose!(G, submanifold_component(SDPG, k, Val(g_ind)), submanifold_component(SDPG, g, Val(g_ind)), submanifold_component(SDPG, h, Val(g_ind)))
    # c) group operation on H (not that h_2 is already in k2)
    _compose!(H, submanifold_component(SDPG, k, Val(h_ind)), submanifold_component(SDPG, g, Val(h_ind)), submanifold_component(SDPG, k, Val(h_ind)))
    return k
end

_doc_LSDP_diff_left_compose = """
    diff_left_compose(
        SDPG::LieGroup{𝔽,LeftSemidirectProductGroupOperation,<:ProductManifold}, g, h, X
    ) where {𝔽}
    diff_left_compose!(
        SDPG::LieGroup{𝔽,LeftSemidirectProductGroupOperation,<:ProductManifold}, Y, g, h, X
    ) where {𝔽}

Compute the differential of the left group operation ``λ_g``, that is ``D_{λ_g}(h)[X]``.
For this case it is given by

TODO Update formula the Diff of Lambda is the one with respect to g.

```math
    D_{λ_g}(h)[X] = $(_tex(:bigl))( D_{λ_{g_1}}(h_1)[X_1], D_{λ_{g_2}}(σ_{g_1}(h_2))$(_tex(:bigl))[ D_{σ_{g_1}}(h_2)[X_2]$(_tex(:bigr))]$(_tex(:bigr)))
```
where ``D_{λ_{g_2}}(σ_{g_1}(h_2))`` is given by [`diff_group_apply`](@ref)`(A, h_2, g_1, X_2)` with ``A`` denotes the [`GroupAction`](@ref) ``σ``.
"""

"$(_doc_LSDP_diff_left_compose)"
diff_left_compose(
    SDPG::LieGroup{𝔽, <:LeftSemidirectProductGroupOperation, <:ProductManifold}, g, h, X
) where {𝔽}

"$(_doc_LSDP_diff_left_compose)"
function diff_left_compose!(
        SDPG::LieGroup{𝔽, <:LeftSemidirectProductGroupOperation, <:ProductManifold}, Y, g, h, X
    ) where {𝔽}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)

    Y1, Y2 = submanifold_components(LieAlgebra(SDPG), Y)
    X1, X2 = submanifold_components(LieAlgebra(SDPG), X)
    g1, g2 = submanifold_components(SDPG, g)
    h1, h2 = submanifold_components(SDPG, h)

    # We have to perform 3 steps applying the group action
    # 1) for the left this is just a diff on that group
    diff_left_compose!(G, Y1, g1, h1, X1)
    # For the second (right) it is diff_compose applied to the diff_apply of the group action
    # where we can do that diff apply already in-place
    diff_group_apply!(A, Y2, g1, h2, X2)
    # and then apply diff compose for the right
    x = copy(G, h2)
    # we need the point on G where we apply to
    apply!(A, x, g1, h2)
    diff_left_compose!(H, Y2, g2, x, Y2)
    return Y
end

_doc_RSDP_diff_left_compose = """
    diff_left_compose(
        SDPG::LieGroup{𝔽,RightSemidirectProductGroupOperation,<:ProductManifold}, g, h, X
    ) where {𝔽}
    diff_left_compose!(
        SDPG::LieGroup{𝔽,RightSemidirectProductGroupOperation,<:ProductManifold}, Y, g, h, X
    ) where {𝔽}

Compute the differential of the left group operation ``λ_g``, that is ``D_{λ_g}(h)[X]``.
For this case it is given by

TODO Update formula the Diff of Lambda is the one with respect to g.

```math
    D_{λ_g}(h)[X] = $(_tex(:bigl))( D_{λ_{g_1}}(σ_{g_2}(h_1))$(_tex(:bigl))[ D_{σ_{g_2}}(h_1)[X_1], D_{λ_{g_2}}(h_2)[X_2]$(_tex(:bigr))]$(_tex(:bigr)))
```
where ``D_{λ_{g_2}}(σ_{g_1}(h_2))`` is given by [`diff_group_apply`](@ref)`(A, h_2, g_1, X_2)` with ``A`` denotes the [`GroupAction`](@ref) ``σ``.
"""

"$(_doc_RSDP_diff_left_compose)"
diff_left_compose(
    SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}, g, h, X
) where {𝔽}

"$(_doc_RSDP_diff_left_compose)"
function diff_left_compose!(
        SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}, Y, g, h, X
    ) where {𝔽}
    PM = SDPG.manifold
    H, G = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)

    Y1, Y2 = submanifold_components(LieAlgebra(SDPG), Y)
    X1, X2 = submanifold_components(LieAlgebra(SDPG), X)
    g1, g2 = submanifold_components(SDPG, g)
    h1, h2 = submanifold_components(SDPG, h)

    # We have to perform 3 steps applying the group action
    # 1) for the right this is just a diff on that group
    diff_left_compose!(G, Y2, g2, h2, X2)
    # For the second (left) it is diff_compose applied to the diff_apply of the group action
    # where we can do that diff apply already in-place
    diff_group_apply!(A, Y1, g2, h1, X1)
    # and then apply diff compose for the right
    x = copy(G, h1)
    # we need the point on G where we apply to
    apply!(A, x, g2, h1)
    diff_left_compose!(H, Y1, g1, x, Y1)
    return Y
end

function get_vector_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, X, c, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: SemiDirectProductGroupOperation, M <: ProductManifold}
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

"""
    inv(SDPG::LieGroup{𝔽,Op,M}, g) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ProductManifold}

Compute the inverse element of an element ``g = (g_1, g_2)`` given by

```math
g^{-1} = (g_1^{-1}, σ_{g_1^{-1}}g_2^{-1}).
```

for the left variant and

```math
g^{-1} = (σ_{g_2^{-1}} g_1^{-1}, g_2^{-1})
```

for the right variant, respectively. See also [HilgertNeeb:2012; Proof of Lemma 2.2.3](@cite).
"""
Base.inv(
    SDPG::LieGroup{𝔽, Op, M}, g
) where {𝔽, Op <: SemiDirectProductGroupOperation, M <: ProductManifold}

function inv!(
        SDPG::LieGroup{𝔽, O, M}, k, g
    ) where {𝔽, O <: LeftSemidirectProductGroupOperation, M <: ProductManifold}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    inv!(G, submanifold_component(SDPG, k, Val(1)), submanifold_component(PM, g, Val(1)))
    inv!(H, submanifold_component(SDPG, k, Val(2)), submanifold_component(PM, g, Val(2)))
    apply!( # Apply the group action with g1^-1 to g2^-1
        A,
        submanifold_component(SDPG, k, Val(2)),
        submanifold_component(SDPG, k, Val(1)),
        submanifold_component(SDPG, k, Val(2)), #TODO make sure apply is safe against side effects
    )
    return k
end
function inv!(
        SDPG::LieGroup{𝔽, O, M}, k, g
    ) where {𝔽, O <: RightSemidirectProductGroupOperation, M <: ProductManifold}
    PM = SDPG.manifold
    H, G = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    inv!(H, submanifold_component(SDPG, k, Val(1)), submanifold_component(PM, g, Val(1)))
    inv!(G, submanifold_component(SDPG, k, Val(2)), submanifold_component(PM, g, Val(2)))
    apply!( # Apply the group action with g2^-1 to g1^-1
        A,
        submanifold_component(SDPG, k, Val(1)),
        submanifold_component(SDPG, k, Val(2)),
        submanifold_component(SDPG, k, Val(1)), #TODO make sure apply is safe against side effects
    )
    return k
end
function inv!(
        SDPG::LieGroup{𝔽, O, M}, k, ::Identity{O}
    ) where {𝔽, O <: LeftSemidirectProductGroupOperation, M <: ProductManifold}
    PrM = SDPG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, SDPG.op.operations),
        submanifold_components(PrM, k),
        map(Identity, SDPG.op.operations),
    )
    return k
end
function inv!(
        SDPG::LieGroup{𝔽, O, M}, k, ::Identity{O}
    ) where {𝔽, O <: RightSemidirectProductGroupOperation, M <: ProductManifold}
    PrM = SDPG.manifold
    map(
        inv!,
        map(LieGroup, PrM.manifolds, SDPG.op.operations),
        submanifold_components(PrM, k),
        map(Identity, SDPG.op.operations),
    )
    return k
end
function identity_element!(
        SDPG::LieGroup{𝔽, Op, M}, e
    ) where {𝔽, Op <: SemiDirectProductGroupOperation, M <: ProductManifold}
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

function get_coordinates_lie!(
        Pr𝔤::LieAlgebra{𝔽, Op, LieGroup{𝔽, Op, M}}, c, X, B::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽, Op <: SemiDirectProductGroupOperation, M <: ProductManifold}
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
