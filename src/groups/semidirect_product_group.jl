"""
    SemiDirectProductGroupOperation{
        O1<:AbstractGroupOperation,
        O2<:AbstractGroupOperation,
        A<:AbstractGroupActionType
    } <: AbstractProductGroupOperation

An abstract type for all semidirect product group operations.
"""
abstract type SemiDirectProductGroupOperation{
    O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType,
} <: AbstractProductGroupOperation end

"""
    LeftSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a left semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.

We further need a [`GroupAction`](@ref). We have two choices:

For a left group action, we use an [`AbstractLeftGroupActionType`](@ref)
``σ: $(_tex(:Cal, "G"))×$(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``
and use here the notation ``σ_g: $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``, ``g  ∈ G``,
as a family of maps on ``$(_tex(:Cal, "H"))``.

In this case the group operation ``∘`` on the product manifold ``$(_tex(:Cal, "G"))×$(_tex(:Cal, "H"))`` is defined by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ σ_{g_1}(h_2) $(_tex(:bigr))).
```

For a right group action, we use an [`AbstractRightGroupActionType`](@ref)
``τ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "G")) → $(_tex(:Cal, "H"))``
and use here the notation ``τ_g: $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``, ``g  ∈ G``,
as a family of maps on ``$(_tex(:Cal, "H"))``.

In this case the group operation ``∘`` on the product manifold ``$(_tex(:Cal, "G"))×$(_tex(:Cal, "H'"))`` is defined by

```math
    (g_1,h_1) ∘ (g_2,h_2) := $(_tex(:bigl))( g_1 ⋆ g_2, h_1 ⋄ τ_{g_1^{-1}}(h_2) $(_tex(:bigr))).
```

For both we use the shorthand notation ``G``[`⋉`](@ref)``H`` = (G×H,∘)``.
More details can be found in [`HilgertNeeb:2012; Definition 9.2.22](@cite), second definition.

# Constructor

    LeftSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "G"))``
* `op2::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `action::AbstractGroupActionType` The group action ``σ`` or ``τ`` of ``$(_tex(:Cal, "G"))`` acting on ``$(_tex(:Cal, "H"))``

!!! note "A note on left/right"
    The “left” in the left semidirect product refers to the side, where the “pure” group operation takes place
    The “left/right” for the action refers to the type of group action used
    Formally there is a third “left/right” choice: which side the action is applied to.
    Since 2 choices add enough complexity to the interface and the implementation, this third
    choice is for now omitted and we decided for the more common case as an action on the right (``h_2``) as default.
"""
struct LeftSemidirectProductGroupOperation{
        O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType,
    } <: SemiDirectProductGroupOperation{O1, O2, A}
    operations::Tuple{O1, O2}
    action_type::A
    function LeftSemidirectProductGroupOperation(
            op1::O1, op2::O2, action::A
        ) where {
            O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType,
        }
        return new{O1, O2, A}((op1, op2), action)
    end
end
@inline Base.getindex(spgo::SemiDirectProductGroupOperation, i::Integer) =
    spgo.operations[i]

"""
    RightSemidirectProductGroupOperation{O1,O2,A} <: SemiDirectProductGroupOperation{O1,O2,A}

A struct to model a right semidirect Lie group product.

Let ``($(_tex(:Cal, "G")), ⋆)`` and ``($(_tex(:Cal, "H")), ⋄)`` be two Lie groups
with group operations ``⋆`` and ``⋄``, respectively.

We further need a [`GroupAction`](@ref). We have two choices:

For a right group action, we use an [`AbstractRightGroupActionType`](@ref)
``τ: $(_tex(:Cal, "H"))×$(_tex(:Cal, "G")) → $(_tex(:Cal, "H"))``
and use here the notation ``τ_g: $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``, ``g  ∈ G``,
as a family of maps on ``$(_tex(:Cal, "H"))``.

In this case the group operation ``∘`` on the product manifold ``$(_tex(:Cal, "G"))×$(_tex(:Cal, "H"))`` is defined by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( h_1 ⋄ τ_{g_1^{-1}}(h_2), g_1 ⋆ g_2 $(_tex(:bigr))).
```

For a left group action, we use an [`AbstractLeftGroupActionType`](@ref)
``σ: $(_tex(:Cal, "G"))×$(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``
and use here the notation ``σ_g: $(_tex(:Cal, "H")) → $(_tex(:Cal, "H"))``, ``g  ∈ G``,
as a family of maps on ``$(_tex(:Cal, "H"))``.

In this case the group operation ``∘`` on the product manifold ``$(_tex(:Cal, "H"))×$(_tex(:Cal, "G"))`` is defined by

```math
    (h_1,g_1) ∘ (h_2,g_2) := $(_tex(:bigl))( h_1 ⋄ σ_{g_1}(h_2), g_1 ⋆ g_2 $(_tex(:bigr))).
```

For both we use the shorthand notation ``H``[`⋊`](@ref)``G = (H×G,∘)``.
See [HilgertNeeb:2012; Definition 9.2.22](@cite), first definition for more details.

# Constructor

    RightSemidirectProductGroupOperation(
        op1::AbstractGroupOperation,
        op2::AbstractGroupOperation,
        action::AbstractGroupActionType
    )

# Parameters

* `op1::AbstractGroupOperation`: The group operation ``⋄`` on ``$(_tex(:Cal, "H"))``
* `op2::AbstractGroupOperation`: The group operation ``⋆`` on ``$(_tex(:Cal, "G"))``
* `action::AbstractGroupActionType`: The group action ``σ`` or ``τ`` of ``$(_tex(:Cal, "G"))`` acting on ``$(_tex(:Cal, "H"))``

!!! note "A note on left/right"
    The “right” in the right semidirect product refers to the side, where the “pure” group operation takes place.
    The “left/right” for the action refers to the type of group action used.
    Formally there is a third “left/right” choice: which side the action is applied to.
    Since the first two choices add enough complexity to the interface and the implementation, this third
    choice is for now omitted and we decided for the more common case as an action on the right (``h_2``) as default.
"""
struct RightSemidirectProductGroupOperation{
        O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType,
    } <: SemiDirectProductGroupOperation{O1, O2, A}
    operations::Tuple{O1, O2}
    action_type::A
    function RightSemidirectProductGroupOperation(
            op1::O1, op2::O2, action::A
        ) where {
            O1 <: AbstractGroupOperation, O2 <: AbstractGroupOperation, A <: AbstractGroupActionType,
        }
        return new{O1, O2, A}((op1, op2), action)
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
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType = default_left_action(N, H)
    )
    # Use product manifold instead of × to not accidentally splat.
    return LieGroup(
        ProductManifold(N.manifold, H.manifold), LeftSemidirectProductGroupOperation(N.op, H.op, action)
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
        N::LieGroup, H::LieGroup, action::AbstractGroupActionType = default_right_action(N, H)
    )
    # Use product manifold instead of × to not accidentally splat.
    return LieGroup(
        ProductManifold(N.manifold, H.manifold), RightSemidirectProductGroupOperation(N.op, H.op, action)
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

Compute the group operation ``$(_math(:∘))`` on the semidirect product Lie group ``L = G ⋉ H``,
that is for `g` `` = (g_1,h_1)``, `h` ``= (g_2,h_2)`` with ``g_1,g_2 ∈ G``, ``h_1,h_2 ∈ H``
this computes

```math
    (g_1,h_1) ∘ (g_2,h_2) := (g_1 ⋄ g_2, h_1 ⋆ σ_{g_1}(h_2)).
```
where ``∘`` denotes the group operation on ``L``, ``⋄`` and ``⋆`` those on ``G`` and ``H``,
respectively, and ``σ`` is the group action specified by the [`AbstractGroupActionType`](#ref) within the [`LeftSemidirectProductLieGroup`](@ref)  ``L``.
"""
compose(
    SDPG::LieGroup{𝔽, LeftSemidirectProductGroupOperation, <:ProductManifold}, ::Any, ::Any
) where {𝔽}

function _compose!(
        SDPG::LieGroup{𝔽, <:LeftSemidirectProductGroupOperation, <:ProductManifold}, k, g, h
    ) where {𝔽}
    PM = SDPG.manifold
    G, H = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    # We have to perform 3 steps applying the group action
    # 1) x = σ_g[1](h[2]) (a point on H)
    # 2) compose g[1] and h[1] (a point on G, neither input needed afterwards)
    # 3) compose g[2] and x (a point on H, neither input needed afterwards)
    # to avoid to overwrite elements in case k=g or k=h: allocate for result of (1)
    # especially after (1) we still need g[2] (in case k=g)
    x = copy(H, submanifold_component(SDPG, k, Val(2)))
    # (1)
    apply!(
        A, x, submanifold_component(SDPG, g, Val(1)), submanifold_component(SDPG, h, Val(2))
    )
    # (2)
    _compose!(
        G,
        submanifold_component(SDPG, k, Val(1)),
        submanifold_component(SDPG, g, Val(1)),
        submanifold_component(SDPG, h, Val(1)),
    )
    _compose!(
        H, submanifold_component(SDPG, k, Val(2)), submanifold_component(SDPG, g, Val(2)), x
    )
    return k
end
"""
    compose(L::LieGroup{𝔽,RightSemidirectProductGroupOperation}, g, h)

Compute the group operation ``$(_math(:∘))`` on the semidirect product Lie group ``L = G ⋊ H``,
that is for `g` `` = (g_1,h_1)``, `h` ``= (g_2,h_2)`` with ``g_1,g_2 ∈ G``, ``h_1,h_2 ∈ H``
this computes

```math
    (g_1,h_1) ∘ (g_2,h_2) := (g_1 ⋄ σ_{h_1}(g_2), h_1 ⋆ h_2),
```

where ``∘`` denotes the group operation on ``L``, ``⋄`` and ``⋆`` those on ``G`` and ``H``,
respectively, and ``σ`` is the group action specified by the [`AbstractGroupActionType`](#ref) within the [`RightSemidirectProductLieGroup`](@ref) ``L``.
"""
compose(
    SDPG::LieGroup{𝔽, RightSemidirectProductGroupOperation, <:ProductManifold}, ::Any, ::Any
) where {𝔽}
function _compose!(
        SDPG::LieGroup{𝔽, <:RightSemidirectProductGroupOperation, <:ProductManifold}, k, g, h
    ) where {𝔽}
    PM = SDPG.manifold
    H, G = map(LieGroup, PM.manifolds, SDPG.op.operations)
    A = GroupAction(SDPG.op.action_type, G, H)
    # We have to perform 3 steps applying the group action
    # 1) x = σ_g[2](h[1]) (a point on H)
    # 2) compose g[1] and x (a point on H)
    # 3) compose g[2] and h[2] (a point on G)
    # to avoid to overwrite elements in case k=g or k=h: allocate for result of (1)
    # especially after (1) we still need g[1] (in case k=g)
    x = copy(H, submanifold_component(SDPG, k, Val(1)))
    # (1)
    apply!(
        A, x, submanifold_component(SDPG, g, Val(2)), submanifold_component(SDPG, h, Val(1))
    )
    # (2)
    _compose!(
        H, submanifold_component(SDPG, k, Val(1)), submanifold_component(SDPG, g, Val(1)), x
    )
    # (3)
    _compose!(
        G,
        submanifold_component(SDPG, k, Val(2)),
        submanifold_component(SDPG, g, Val(2)),
        submanifold_component(SDPG, h, Val(2)),
    )
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
