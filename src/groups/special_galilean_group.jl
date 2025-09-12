"""
    RotationBoostAction

The group action of the semidirect product of spatial rotations and velocity boosts (SO(n) ⋉ ℝⁿ) on the space of events (position, time).
"""
struct RotationBoostAction <: AbstractLeftGroupActionType end

"""
    LieGroups.apply!(A::GroupAction{RotationBoostAction}, k, g, h)

Apply the action of the rotation-boost semidirect product group (SO(n) ⋉ ℝⁿ) on an event (x, t).
Given group element `g = (R, v)` and event `h = (x, t)`, computes the transformed event `k = (R*x + v*t, t)`.
See [Kelly:2025; section 4.1](@cite).
"""
function LieGroups.apply!(A::GroupAction{RotationBoostAction}, k, g, h)
    G = A.group
    R = submanifold_component(G, g, Val(1))
    v = submanifold_component(G, g, Val(2))

    H = A.manifold
    p = submanifold_component(H, h, Val(1))
    t = submanifold_component(H, h, Val(2))[1]

    k_p = submanifold_component(H, k, Val(1))
    k_p .= v * t + R * p

    k_t = submanifold_component(H, k, Val(2))
    k_t .= t
    return k
end

const LeftSpecialGalileanGroupOperation = LeftSemidirectProductGroupOperation{
    LeftSemidirectProductGroupOperation{
        MatrixMultiplicationGroupOperation, AdditionGroupOperation, LeftMultiplicationGroupAction, ActionActsOnRight,
    },
    ProductGroupOperation{Tuple{AdditionGroupOperation, AdditionGroupOperation}},
    RotationBoostAction,
    ActionActsOnRight,
}

const SpecialGalileanGroup{T} = LieGroup{
    ℝ,
    <:LeftSpecialGalileanGroupOperation,
    <:ProductManifold{
        ℝ,
        Tuple{
            <:ProductManifold{ℝ, Tuple{<:Rotations{T}, <:Euclidean{T, ℝ}}},
            <:ProductManifold{
                ℝ,
                Tuple{
                    <:Euclidean{T, ℝ}, <:Euclidean{ManifoldsBase.TypeParameter{Tuple{1}}, ℝ},
                },
            },
        },
    },
}

"""
    SpecialGalileanGroup(n::Int)

Construct the special Galilean group SGal(n) as a nested semidirect product:
    ``(SO(n) ⋉ ℝⁿ) ⋉ (ℝⁿ × ℝ)``
where ``R ∈ SO(n)`` are spatial rotations, ``v ∈ ℝⁿ`` are velocity boosts, and ``(p, t) ∈ (ℝⁿ × ℝ)`` are the (position, time) events.
The affine representation of the group is given by the matrix:
```math
\\mathrm{SGal}(3) = \\begin{bmatrix}
R & v & p \\\\
0 & 1 & t \\\\
0 & 0 & 1
\\end{bmatrix} \\subset \\mathbb{R}^{5\\times 5}
```
And the ArrayPartition representation as:
``((R, v), (p, t))``

!!! note "Technical Detail"
    The ArrayPartition (default) implementation requires `RecursiveArrayTools.jl` to be loaded.

[Kelly:2025](@cite)
"""
function SpecialGalileanGroup(n::Int)
    G = SpecialOrthogonalGroup(n) ⋉ TranslationGroup(n)
    N = TranslationGroup(n) × TranslationGroup(1)
    return LieGroup(
        ProductManifold(G.manifold, N.manifold),
        LeftSemidirectProductGroupOperation(G.op, N.op, RotationBoostAction(), ActionActsOnRight())
    )
end

#
#
# doc strings

_doc_SGal3_exp = """
    LieGroups.exp(M::SpecialGalileanGroup, X)
    LieGroups.exp!(M::SpecialGalileanGroup, h, X)

Compute the Lie group exponential function on the [`SpecialGalileanGroup`](@ref)`(3)`,
where `X` is an element of the Lie algebra represented as an `ArrayPartition`.

The closed-form expression for the matrix exponential from [Kelly:2025; section 6](@cite) is used.

The computation can be done in-place of `h`.
"""

@doc "$(_doc_SGal3_exp)"
ManifoldsBase.exp!(
    ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, h, X,
)

@doc "$(_doc_SGal3_exp)"
LieGroups.exp(
    ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X
)

_doc_SGal3_log = """
    LieGroups.log(M::SpecialGalileanGroup, g)
    LieGroups.log!(M::SpecialGalileanGroup, X, g)

Compute the Lie group logarithm function on the [`SpecialGalileanGroup`](@ref)`(3)`,
where `g` is a group element represented as an `ArrayPartition`.    

The closed-form expression from [Kelly:2025; section 6](@cite) is used.

The computation can be done in-place of `X`.
"""

@doc "$(_doc_SGal3_log)"
LieGroups.log(
    ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g
)

@doc "$(_doc_SGal3_log)"
ManifoldsBase.log!(
    ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, g,
)
