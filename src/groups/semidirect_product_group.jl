
"""
    SemiDirectProductGroupOperation{
        O1<:AbstractGroupOperation,
        O2<:AbstractGroupOperation,
        A<:AbstractGroupActionType
    } <: AbstractProductGroupOperation

An abstract type for all semdirect product group operations
"""
abstract type SemiDirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: AbstractProductGroupOperation end

"""
    LeftSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a semidirect Lie group product.

Let ``($(_tex(:Cal, "N")), ⋄)`` and ``($(_tex(:Cal, "H")), ⋆)`` be two Lie groups
with group operations ``⋄`` and ``⋆``, respectively, as well as a group action
``σ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``, cf [`AbstractLeftGroupActionType`](@ref).

We use here as well use the notation ``σ_h: $(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``
as a family of maps on ``$(_tex(:Cal, "N"))``

Then we define a group operation ``∘`` on the product manifold ``$(_tex(:Cal, "N"))×$(_tex(:Cal, "H"))`` by

```math
    (h_1,n_1) ∘ (h_2,n_2) := (h_1 ⋆ h_2, σ_{h_2}(n_1) ⋄ n_2).
```

See [HilgertNeeb:2012; Definition 9.2.22](@cite), second definition for more details.

# Constructor

    LeftSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `op2::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "N"))``
* `action::AbstractGroupActionType` The group action ``σ`` of ``$(_tex(:Cal, "H"))`` on ``$(_tex(:Cal, "N"))``

"""
struct LeftSemidirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: SemiDirectProductGroupOperation{O1,O2,A}
    operations::Tuple{O1,O2}
    action_type::A
    function LeftSemidirectProductGroupOperation(
        op1::O1, op2::O2, action::A
    ) where {
        O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
    }
        return new{O1,O2,A}((op1, op2), action)
    end
end

"""
    RightSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a right semidirect Lie group product.

Let ``($(_tex(:Cal, "N")), ⋄)`` and ``($(_tex(:Cal, "H")), ⋆)`` be two Lie groups
with group operations ``⋄`` and ``⋆``, respectively, as well as a group action
``σ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``, cf [`AbstractGroupActionType`](#ref).

We use here as well use the notation ``σ_h: $(_tex(:Cal, "N")) → $(_tex(:Cal, "N"))``
as a family of maps on ``$(_tex(:Cal, "N"))``

Then we define a group operation ``∘`` on the product manifold ``$(_tex(:Cal, "N"))×$(_tex(:Cal, "H"))`` by

```math
    (n_1,h_1) ∘ (n_2,h_2) := (n_1 ⋄ σ_{h_1}(n_2), h_1 ⋆ h_2)
```

See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition for more details.

# Constructor

    RightSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "N"))``
* `op2::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `action::AbstractGroupActionType`: The group action ``σ`` of ``$(_tex(:Cal, "H"))`` on ``$(_tex(:Cal, "N"))``

"""
struct RightSemidirectProductGroupOperation{
    O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
} <: SemiDirectProductGroupOperation{O1,O2,A}
    operations::Tuple{O1,O2}
    action_type::A
    function RightSemidirectProductGroupOperation(
        op1::O1, op2::O2, action::A
    ) where {
        O1<:AbstractGroupOperation,O2<:AbstractGroupOperation,A<:AbstractGroupActionType
    }
        return new{O1,O2,A}((op1, op2), action)
    end
end

"""
    LeftSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_left_action(N, H)
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) = N ⋉ H`` for an [`AbstractLeftGroupActionType`](@ref)
using the [`LeftSemidirectProductGroupOperation`](@ref) for the group operation definition
as well as [HilgertNeeb:2012; Definition 9.2.22](@cite), second definition, for more details.

The short form `N `[`⋉`](@ref ⋉(L1::LieGroup, L2::LieGroup))` H` can be used if the
corresponding [`default_left_action`](@ref)`(N,H)` is the one you want to use.
"""
function LeftSemidirectProductLieGroup(
    N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_left_action(N, H)
)
    return LieGroup(
        N.manifold × H.manifold, LeftSemidirectProductGroupOperation(N.op, H.op, action)
    )
end

"""
    RightSemidirectProductLieGroup(
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_right_action(N,H)
    )

Generate the semidirect product Lie Group ``$(_tex(:Cal, "G")) = N ⋊ H`` for an [`AbstractLeftGroupActionType`](@ref)
using the [`RightSemidirectProductGroupOperation`](@ref) for the group operation definition
as well as [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition, for more details.

The short form `N `[`⋊`](@ref ⋊(L1::LieGroup, L2::LieGroup))` H` can be used if the
corresponding [`default_right_action`](@ref)`(N,H)` is the one you want to use.
"""
function RightSemidirectProductLieGroup(
    N::LieGroup, H::LieGroup, action::AbstractGroupActionType=default_right_action(N, H)
)
    return LieGroup(
        N.manifold × H.manifold, RightSemidirectProductGroupOperation(N.op, H.op, action)
    )
end

"""
    L1 ⋉ L2
    ⋉(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`LeftSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_left_action`](@ref)`(L1, L2)` is used.
"""
function ⋉(L1::LieGroup, L2::LieGroup)
    return LeftSemidirectProductLieGroup(L1, L2, default_left_action(L1, L2))
end

"""
    L1 ⋊ L2
    ⋊(L1, L2)

For two [`LieGroups`](@ref) `L1`, `L2`, generate the [`RightSemidirectProductLieGroup`](@ref)`(L1, L2)`,
where the corresponding [`default_right_action`](@ref)`(L1, L2)` is used.
"""
function ⋊(L1::LieGroup, L2::LieGroup)
    return RightSemidirectProductLieGroup(L1, L2, default_right_action(L1, L2))
end

#
#
# Functions

"""
    compose(L::LieGroup{𝔽,LeftSemidirectProductGroupOperation}, g, h)

Compute the group operation $(_math(:∘))``on the semidirect product Lie group ``L = G ⋉ H``,
that is for `g` = ``(g_1,h_1)``, `h` ``= (g_2,h_2)`` with ``g_1,g_2 ∈ G``, ``h_1,h_2 ∈ H``
this computes

```math
    (g_1,h_1) ∘ (g_2,h_2) := (g_1 ⋄ σ_{h_1}(g_2), h_1 ⋆ h_2).
```
"""
compose!(
    SDPG::LieGroup{𝔽,LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold}
) where {𝔽}

function _compose!(
    SDPG::LieGroup{𝔽,<:LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
    k,
    g,
    h,
) where {𝔽}
    PM = SDPG.manifold
    G, H = LieGroup.(PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    # for the first components, just perform the group op
    _compose!(
        G,
        submanifold_component(PM, k, 1),
        submanifold_component(PM, g, 1),
        submanifold_component(PM, h, 1),
    )
    # apply the first element from g to
    apply!(
        A,
        submanifold_component(PM, k, 2),
        submanifold_component(PM, h, 1),
        submanifold_component(PM, g, 2),
    )
    _compose!(
        H,
        submanifold_component(PM, k, 2),
        submanifold_component(PM, k, 2),
        submanifold_component(PM, h, 2),
    )
    return k
end
function _compose!(
    SDPG::LieGroup{
        𝔽,<:RightSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold
    },
    k,
    g,
    h,
) where {𝔽}
    PM = SPDG.manifold
    G, H = LieGroup.(SDPG.manifold, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    # for the first components, just perform the group op
    # apply the first element from g to
    apply!(
        A,
        submanifold_component(PM, k, 1),
        submanifold_component(PM, g, 2),
        submanifold_component(PM, h, 1),
    )
    _compose!(
        G,
        submanifold_component(PM, k, 1),
        submanifold_component(PM, g, 1),
        submanifold_component(PM, k, 1),
    )
    _compose!(
        H,
        submanifold_component(PM, k, 2),
        submanifold_component(PM, g, 2),
        submanifold_component(PM, h, 2),
    )
    return k
end

function identity_element!(
    SDPG::LieGroup{𝔽,Op,M}, e
) where {𝔽,Op<:SemiDirectProductGroupOperation,M<:ManifoldsBase.ProductManifold}
    GH = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    identity_element!.(GH, submanifold_components(SDPG.manifold, e))
    return e
end

function Base.show(
    io::IO,
    SDPG::LieGroup{𝔽,<:LeftSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold},
) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = SDPG.op.action_type
    return print(io, "LeftSemidirectProductLieGroup($G, $H, $at)")
end
function Base.show(
    io::IO,
    RSDL::LieGroup{
        𝔽,<:RightSemidirectProductGroupOperation,<:ManifoldsBase.ProductManifold
    },
) where {𝔽}
    G, H = LieGroup.(SDPG.manifold.manifolds, SDPG.op.operations)
    at = RSDL.op.action_type
    return print(io, "RightSemidirectProductLieGroup($G, $H, $at)")
end
