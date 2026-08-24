"""
    RotationBoostAction

The group action of the semidirect product of spatial rotations and velocity boosts (``(R, v) \\in SO(n) ⋉ ℝⁿ``) 
on the space of events (position, time) (``(p, t) \\in ℝⁿ × ℝ``).
See [Kelly:2025; section 4.1](@cite) and [`apply!`](@ref).
"""
struct RotationBoostAction <: AbstractLeftGroupActionType end

"""
    LieGroups.apply!(A::GroupAction{RotationBoostAction}, k, g, h)

Apply the action of the rotation-boost semidirect product group (SO(n) ⋉ ℝⁿ) on an event ``(p, t)``.
Given group element ``g = (R, v)`` and event ``h = (p, t)``, computes the transformed event ``k = (Rp + vt, t)``.
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
            <:ProductManifold{ℝ, Tuple{<:Rotations{T}, <:Euclidean{ℝ, T}}},
            <:ProductManifold{
                ℝ,
                Tuple{
                    <:Euclidean{ℝ, T}, <:Euclidean{ℝ, ManifoldsBase.TypeParameter{Tuple{1}}},
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

The group operation ([`compose`](@ref)) is given by:
```math
((R_1, v_1), (p_1, t_1)) \\circ ((R_2, v_2), (p_2, t_2))
= ((R_1 R_2, v_1 + R_1 v_2), (p_1 + v_1 t_2 + R_1 p_2, t_1 + t_2))
``` 
and the identity element ([`identity_element`](@ref)) is ``((I_n, \\mathbf{0}), (\\mathbf{0}, 0))``.

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

_doc_jacobian_exp_SGal3 = raw"""
    jacobian_exp(G::SpecialGalileanGroup, X, ::DefaultLieAlgebraOrthogonalBasis)
    jacobian_exp!(G::SpecialGalileanGroup, J, X, ::DefaultLieAlgebraOrthogonalBasis)

Compute the Jacobian of the Lie group exponential in a basis of the Lie algebra on the
[`SpecialGalileanGroup`](@ref)`(3)`.

The closed form ``\mathbf{J}_ℓ`` of [Kelly:2025; section 8, equations (31)–(36)](@cite) is
used, evaluated at ``-ξ`` to match the convention of [`jacobian_exp`](@ref).
In the coordinate order ``ξ = (ρ, ν, ϕ, ι)`` (see [`hat`](@ref)) it has the block structure

````math
\mathbf{J}_ℓ(ξ) = \begin{pmatrix}
\mathbf{D} & -\mathbf{L}ι & \mathbf{N} & \mathbf{E}ν \\
\mathbf{0} & \mathbf{D} & \mathbf{M} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{D} & \mathbf{0} \\
\mathbf{0} & \mathbf{0} & \mathbf{0} & 1
\end{pmatrix} ∈ ℝ^{10×10},
````

where ``\mathbf{D}`` is the corresponding matrix for ``\mathrm{SO}(3)``, ``\mathbf{E}`` and
``\mathbf{L}`` are given by [Kelly:2025; equations (19) and (32)](@cite),
``\mathbf{M}`` and ``\mathbf{N} = \mathbf{N}_1 - \mathbf{N}_2`` by
[Kelly:2025; equations (33)–(36)](@cite).

For small rotation angles the Jacobian is evaluated by truncating the series
``\sum_{n ≥ 0} \frac{1}{(n+1)!} \operatorname{ad}_ξ^n`` of the adjoint matrix
[Kelly:2025; equation (28)](@cite), which is numerically robust there.
"""

@doc "$(_doc_jacobian_exp_SGal3)"
jacobian_exp(::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, basis = DefaultLieAlgebraOrthogonalBasis())

@doc "$(_doc_jacobian_exp_SGal3)"
jacobian_exp!(::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, J, X, basis = DefaultLieAlgebraOrthogonalBasis())

# Internal function to compute the skew-symmetric matrix as an SMatrix used for performance.
# Can be replaced with hat(SO(3), v) once that works without allocations.
function _skew(v::AbstractVector{T}) where {T <: Real}
    return SMatrix{3, 3, T}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
end

# Jacobian of SGal(3), [Kelly:2025, eqs. (31)-(36)], coordinate order (ρ, ν, ϕ, ι);
# the blocks M and N₁ share one closed form with different arguments, see _sgal3_MN₁
function _jacobian_exp_left_SGal3!(J::AbstractMatrix, ρ, ν, ω, ι)
    φ = sqrt(ω[1]^2 + ω[2]^2 + ω[3]^2)
    fill!(J, 0)
    if φ < 0.15
        # Truncated series of the adjoint matrix [Kelly:2025, eq. (28)] to avoid the
        # catastrophic cancellation of the closed form below, whose N₂ block (eq. (36))
        # divides by φ³ and loses ~3 extra digits per decade as φ → 0. The switch point
        # φ = 0.15 is where the closed form stops being the more accurate branch;
        # so 9 terms reach ~1e-14 at φ = 0.15. Both were chosen by measuring each
        # branch against a BigFloat ground truth.
        mad = zeros(MMatrix{10, 10, Float64})
        Ωx = _skew(ω)
        mad[1:3, 1:3] .= Ωx
        mad[4:6, 4:6] .= Ωx
        mad[7:9, 7:9] .= Ωx
        for i in 1:3
            mad[i, 3 + i] = -ι
        end
        mad[1:3, 7:9] .= _skew(ρ)
        mad[1:3, 10] .= ν
        mad[4:6, 7:9] .= _skew(ν)
        ad = SMatrix(mad)
        term = SMatrix{10, 10, Float64}(LinearAlgebra.I)
        Js = term
        for k in 1:9
            term = term * ad ./ (k + 1) # tₖ = adᵏ/(k+1)!
            Js += term
        end
        J .= Js
        return J
    end
    u = ω ./ φ
    ux = _skew(u)
    ux² = ux * ux
    νx = _skew(ν)
    sφ, cφ = sincos(φ)
    I₃ = Matrix{Float64}(LinearAlgebra.I, 3, 3)
    # [Kelly:2025, eq. (18)]
    D = I₃ .+ ((1 - cφ) / φ) .* ux .+ ((φ - sφ) / φ) .* ux²
    # [Kelly:2025, eq. (19)]
    E = I₃ ./ 2 .+ ((φ - sφ) / φ^2) .* ux .+ ((φ^2 + 2cφ - 2) / (2φ^2)) .* ux²
    # [Kelly:2025, eq. (32)]
    L = I₃ ./ 2 .+ ((sφ - φ * cφ) / φ^2) .* ux .+
        ((φ^2 + 2 - 2φ * sφ - 2cφ) / (2φ^2)) .* ux²
    # [Kelly:2025, eqs. (33) and (35)]
    M = _sgal3_MN₁(νx, ux, ux², φ, sφ, cφ)
    N₁ = _sgal3_MN₁(_skew(ρ), ux, ux², φ, sφ, cφ)
    # [Kelly:2025, eq. (36)]
    N₂ = (
        ((2 - φ * sφ - 2cφ) / φ^3) .* ux .+ ((φ + φ * cφ - 2sφ) / φ^3) .* ux²
    ) * νx .* ι .+
        ((4sφ - φ^2 * sφ - 4φ * cφ) / (2φ^3)) .* (ux * νx * ux) .* ι .+
        ((4 + φ^2 + φ^2 * cφ - 4φ * sφ - 4cφ) / (2φ^3)) .* (ux² * νx * ux) .* ι .+
        νx * (
        ((φ^2 + 2cφ - 2) / (2φ^3)) .* ux .+ ((sφ - φ) / φ^3) .* ux²
    ) .* ι
    N = N₁ .- N₂
    J[1:3, 1:3] .= D
    J[1:3, 4:6] .= (-ι) .* L
    J[1:3, 7:9] .= N
    J[1:3, 10] .= E * ν
    J[4:6, 4:6] .= D
    J[4:6, 7:9] .= M
    J[7:9, 7:9] .= D
    J[10, 10] = 1
    return J
end

# shared closed form of the blocks M [Kelly:2025, eq. (33)] and N₁ [eq. (35)], with
# wx the skew matrix of ν and ρ, respectively
function _sgal3_MN₁(wx, ux, ux², φ, sφ, cφ)
    return ((1 - cφ) / φ^2) .* wx .+
        ((φ - sφ) / φ^2) .* (ux * wx .+ wx * ux) .+
        ((2 - φ * sφ - 2cφ) / φ^2) .* (ux * wx * ux) .+
        ((2φ + φ * cφ - 3sφ) / φ^2) .* (ux * wx * ux²)
end
