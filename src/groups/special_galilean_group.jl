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
    The ArrayPartition (default) implementation requires `RecursiveArrayTools.jl` to be loaded. The matrix representation is not implemented yet.

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
where `X` is an element of the Lie algebra.

The closed-form expression for the matrix exponential from [Kelly:2025; section 6](@cite) is used.

```math
\\exp X
=
\\exp{\\begin{bmatrix}
\\boldsymbol{\\phi}^\\wedge & \\nu & \\rho \\\\
0 & 0 & \\iota \\\\
0 & 0 & 0
\\end{bmatrix}}
= \\begin{bmatrix}
C & Dν & Dρ + Eνι \\\\
0 & 1 & ι \\\\
0 & 0 & 1
\\end{bmatrix},
```
where
```math
C = I_3 + \\sin(\\phi)\\, \\mathbf{u}^{\\wedge} + \\bigl(1 - \\cos(\\phi)\\bigr)\\, \\mathbf{u}^{\\wedge}\\mathbf{u}^{\\wedge}, \\\\

D = I_3 + \\frac{1 - \\cos(\\phi)}{\\phi} \\, \\mathbf{u}^{\\wedge}
+ \\frac{\\phi - \\sin(\\phi)}{\\phi} \\, \\mathbf{u}^{\\wedge}\\mathbf{u}^{\\wedge}, \\\\

E = \\tfrac12 I_3
+ \\frac{\\phi - \\sin(\\phi)}{\\phi^2} \\, \\mathbf{u}^{\\wedge}
+ \\frac{\\phi^2 + 2\\cos(\\phi) - 2}{2\\phi^2} \\, \\mathbf{u}^{\\wedge}\\mathbf{u}^{\\wedge}.
```
``\\boldsymbol{\\phi}=\\phi \\mathbf{u}`` is the angle-axis rotation parameterization with
``\\phi = \\|\\boldsymbol{\\phi}\\|`` and ``\\mathbf{u} = \\boldsymbol{\\phi}/\\phi``. 

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
where `g` is a group element.

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

_doc_hat_special_galilean = """
    X = hat(𝔤::LieAlgebra{ℝ,LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, c)
    hat!(𝔤::LieAlgebra{ℝ,LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, X, c)

Compute the hat map ``(⋅)^{\\wedge} : ℝ^{10} → 𝔤`` that turns a vector of coordinates `c`
into a tangent vector in the Lie algebra.

```math
\\begin{bmatrix}
\\rho \\\\ \\nu \\\\ \\phi \\\\ \\iota
\\end{bmatrix}^\\wedge
=
\\begin{bmatrix}
\\phi^\\wedge & \\nu & \\rho \\\\
0 & 0 & \\iota \\\\
0 & 0 & 0
\\end{bmatrix}
\\in \\mathbb{R}^{5\\times 5}
```
The basis is defined in eq 14 of [Kelly:2025](@cite).

This can be computed in-place of `X`.
"""

@doc "$(_doc_hat_special_galilean)"
ManifoldsBase.hat(::typeof(LieAlgebra(SpecialGalileanGroup(3))), c)

@doc "$(_doc_hat_special_galilean)"
ManifoldsBase.hat!(::typeof(LieAlgebra(SpecialGalileanGroup(3))), X, c)

_doc_vee_special_galilean = """
    c = vee(𝔤::LieAlgebra{ℝ,<:LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, X)
    vee!(𝔤::LieAlgebra{ℝ,LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, c, X)

Compute the vee map ``(⋅)^{\\vee}: $(_math(:𝔤)) →  ℝ^{10}`` that maps a tangent vector
from the Lie algebra to a vector of coordinates `c`.

```math
\\begin{bmatrix}
\\phi^\\wedge & \\nu & \\rho \\\\
0 & 0 & \\iota \\\\
0 & 0 & 0
\\end{bmatrix}^\\vee =
\\begin{bmatrix}
\\rho \\\\ \\nu \\\\ \\phi \\\\ \\iota
\\end{bmatrix}
\\in \\mathbb{R}^{10}
```
The basis is defined in eq 14 of [Kelly:2025](@cite).

This can be computed in-place of `c`.
"""

@doc "$(_doc_vee_special_galilean)"
ManifoldsBase.vee(
    ::typeof(LieAlgebra(SpecialGalileanGroup(3))), X
)

@doc "$(_doc_vee_special_galilean)"
ManifoldsBase.vee!(
    ::typeof(LieAlgebra(SpecialGalileanGroup(3))), c, X
)
