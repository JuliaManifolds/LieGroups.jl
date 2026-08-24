using LieGroups: SpecialGalileanGroup, _skew
using StaticArrays
using LinearAlgebra

# Internal function to compute the matrix Q used in the exponential and logarithm maps for the Special Galilean group. (D matrix in Kelly:2025)
function _Q(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return SMatrix{3, 3, T}(I)
    else
        u = θ⃗ / θ
        sθ, cθ = sincos(θ)
        uₓ = _skew(u)
        return SMatrix{3, 3, T}(I) + (1 - cθ) / θ * uₓ + (θ - sθ) / θ * uₓ^2
    end
end

# Internal function to compute the matrix P used in the exponential and logarithm maps for the Special Galilean group. (E matrix in Kelly:2025)
function _P(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return 1 / 2 * SMatrix{3, 3, T}(I)
    else
        u = θ⃗ / θ
        sθ, cθ = sincos(θ)
        uₓ = _skew(u)
        return 1 / 2 * SMatrix{3, 3, T}(I) +
            (θ - sθ) / θ^2 * uₓ +
            (cθ + 1 / 2 * θ^2 - 1) / θ^2 * uₓ^2
    end
end

function LieGroups.exp!(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        h::ArrayPartition,
        X::ArrayPartition,
    )
    θ⃗ₓ = X.x[1].x[1] # ωΔt
    ν = X.x[1].x[2]  # aΔt
    ρ = X.x[2].x[1]  # vΔt

    Δt = X.x[2].x[2][1]

    θ⃗ = SA[θ⃗ₓ[3, 2]; θ⃗ₓ[1, 3]; θ⃗ₓ[2, 1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonalGroup(3)
    exp!(M_SO3, h.x[1].x[1], θ⃗ₓ)
    h.x[1].x[2] .= Q * ν
    h.x[2].x[1] .= Q * ρ + P * ν .* Δt
    h.x[2].x[2] .= Δt

    return h
end

function LieGroups.exp(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        X::ArrayPartition{T}
    ) where {T <: Real}

    θ⃗ₓ = X.x[1].x[1] # ωΔt
    ν = X.x[1].x[2]  # aΔt
    ρ = X.x[2].x[1]  # vΔt
    Δt = X.x[2].x[2][1]

    θ⃗ = SA[θ⃗ₓ[3, 2]; θ⃗ₓ[1, 3]; θ⃗ₓ[2, 1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonalGroup(3)
    h = ArrayPartition(
        ArrayPartition(
            exp(M_SO3, θ⃗ₓ),
            Q * ν
        ),
        ArrayPartition(
            Q * ρ + P * ν * Δt,
            copy(X.x[2].x[2])
        )
    )
    return h
end

function LieGroups.log!(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        X::ArrayPartition,
        g::ArrayPartition,
    )
    ΔR = g.x[1].x[1]
    Δv = g.x[1].x[2]
    Δp = g.x[2].x[1]
    Δt = g.x[2].x[2][1]

    SO3 = SpecialOrthogonalGroup(3)
    log!(SO3, X.x[1].x[1], ΔR) # θ⃗ₓ # FIXME allocates
    θ⃗ = vee(LieAlgebra(SO3), X.x[1].x[1])

    P = _P(θ⃗)
    Q = _Q(θ⃗)
    iQ = inv(Q)

    X.x[1].x[2] .= iQ * Δv # ν aΔt
    X.x[2].x[1] .= iQ * (Δp - P * iQ * Δv * Δt) # ρ vΔt
    X.x[2].x[2] .= Δt
    return X
end

function LieGroups.log(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        g::ArrayPartition
    )
    ΔR = g.x[1].x[1]
    Δv = g.x[1].x[2]
    Δp = g.x[2].x[1]
    Δt = g.x[2].x[2][1]

    SO3 = SpecialOrthogonalGroup(3)
    θ⃗ₓ = log(SO3, ΔR)
    θ⃗ = vee(LieAlgebra(SO3), θ⃗ₓ)

    P = _P(θ⃗)
    Q = _Q(θ⃗)
    iQ = inv(Q)
    return ArrayPartition(
        ArrayPartition(
            θ⃗ₓ,
            iQ * Δv # ν aΔt
        ),
        ArrayPartition(
            iQ * (Δp - P * iQ * Δv * Δt), # ρ vΔt
            copy(g.x[2].x[2]) # Δt
        )
    )
end

function LieGroups.identity_element(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{N}}}, ::Type{<:StaticArray}
    ) where {N}
    return ArrayPartition(
        ArrayPartition(
            SMatrix{3, 3, Float64}(I), # ΔR
            @SVector(zeros(3)),        # Δv
        ),
        ArrayPartition(
            @SVector(zeros(3)),        # Δp
            @SVector([0.0]),           # Δt
        ),
    )
end

function LieGroups.inv(::SpecialGalileanGroup, g::ArrayPartition)
    ΔR = g.x[1].x[1]
    Δv = g.x[1].x[2]
    Δp = g.x[2].x[1]
    Δt = g.x[2].x[2]

    return ArrayPartition(
        ArrayPartition(
            ΔR',                      # ΔR
            -ΔR' * Δv,                # Δv
        ),                            #
        ArrayPartition(               #
            -ΔR' * (Δp - Δv * Δt[1]), # Δp
            -Δt,                      # Δt
        ),
    )
end

function LieGroups.compose(::SpecialGalileanGroup, g::ArrayPartition, h::ArrayPartition)
    ΔR = g.x[1].x[1]
    Δv = g.x[1].x[2]
    Δp = g.x[2].x[1]
    Δt = g.x[2].x[2]

    δR = h.x[1].x[1]
    δv = h.x[1].x[2]
    δp = h.x[2].x[1]
    δt = h.x[2].x[2]

    return ArrayPartition(
        ArrayPartition(
            ΔR * δR,                   # ΔR
            Δv + ΔR * δv,              # Δv
        ),
        ArrayPartition(
            Δp + Δv * δt[1] + ΔR * δp, # Δp
            Δt + δt,                   # Δt
        ),
    )
end

_doc_lie_bracket_SGal3 = """
    lie_bracket(𝔰𝔤𝔞𝔩::LieAlgebra{ℝ,<:LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, X, Y)
    lie_bracket!(𝔰𝔤𝔞𝔩::LieAlgebra{ℝ,<:LeftSpecialGalileanGroupOperation,<:SpecialGalileanGroup}, Z, X, Y)

Compute the Lie bracket ``[X, Y] = XY - YX`` of two tangent vectors `X`, `Y` of the Lie algebra of the
[`SpecialGalileanGroup`](@ref)`(3)`, i.e. the matrix commutator of their ``5×5`` `hat` representations.

In the ``((\\Omega, \\nu), (\\rho, \\iota))`` block form (see [`hat`](@ref)) this reduces to
```math
[X, Y] = \\bigl(
(\\Omega_X \\Omega_Y - \\Omega_Y \\Omega_X,\\ \\Omega_X \\nu_Y - \\Omega_Y \\nu_X),\\
(\\Omega_X \\rho_Y - \\Omega_Y \\rho_X + \\iota_Y \\nu_X - \\iota_X \\nu_Y,\\ 0)
\\bigr).
```
The basis is defined in eq. (14) of [Kelly:2025](@cite).

This can be computed in-place of `Z`.
"""

"$(_doc_lie_bracket_SGal3)"
function LieGroups.lie_bracket(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        X::ArrayPartition,
        Y::ArrayPartition,
    )
    ΩX, νX, ρX, ιX = X.x[1].x[1], X.x[1].x[2], X.x[2].x[1], X.x[2].x[2][1]
    ΩY, νY, ρY, ιY = Y.x[1].x[1], Y.x[1].x[2], Y.x[2].x[1], Y.x[2].x[2][1]
    return ArrayPartition(
        ArrayPartition(
            ΩX * ΩY - ΩY * ΩX,                        # Ω
            ΩX * νY - ΩY * νX,                         # ν
        ),
        ArrayPartition(
            ΩX * ρY - ΩY * ρX + ιY * νX - ιX * νY,     # ρ
            zero(X.x[2].x[2]),                        # ι
        ),
    )
end

"$(_doc_lie_bracket_SGal3)"
function LieGroups.lie_bracket!(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        Z::ArrayPartition,
        X::ArrayPartition,
        Y::ArrayPartition,
    )
    ΩX, νX, ρX, ιX = X.x[1].x[1], X.x[1].x[2], X.x[2].x[1], X.x[2].x[2][1]
    ΩY, νY, ρY, ιY = Y.x[1].x[1], Y.x[1].x[2], Y.x[2].x[1], Y.x[2].x[2][1]
    Z.x[1].x[1] .= ΩX * ΩY .- ΩY * ΩX
    Z.x[1].x[2] .= ΩX * νY .- ΩY * νX
    Z.x[2].x[1] .= ΩX * ρY .- ΩY * ρX .+ ιY .* νX .- ιX .* νY
    Z.x[2].x[2] .= 0
    return Z
end

# Dev NOTE: hat and vee use a different bases order than that of the underlining semidirect + direct product groups,
# therefore, get_vector_lie and get_coordinates_lie are implemented explicitly. see hat/vee docstrings for details.
function LieGroups.get_vector_lie(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        c::SVector{10, T},
        ::DefaultLieAlgebraOrthogonalBasis,
        ::Type{<:ArrayPartition{<:Real}}
    ) where {T <: Real}
    return ArrayPartition(
        ArrayPartition(
            _skew(c[SA[7:9...]]), # θ ωΔt
            c[SA[4:6...]]         # ν aΔt
        ),
        ArrayPartition(
            c[SA[1:3...]],        # ρ vΔt
            c[SA[10]],            # Δt
        )
    )
end

function LieGroups.get_vector_lie(
        𝔤::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        c,
        B::DefaultLieAlgebraOrthogonalBasis,
        ::Type{T}
    ) where {T <: ArrayPartition{<:Real}}
    X = zero_vector(𝔤, T)
    return LieGroups.get_vector_lie!(𝔤, X, c, B)
end

function LieGroups.get_vector_lie!(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        X::ArrayPartition{T},
        c,
        ::DefaultLieAlgebraOrthogonalBasis
    ) where {T <: Real}
    X.x[1].x[1] .= _skew(c[SA[7:9...]]) # θ ωΔt
    X.x[1].x[2] .= c[SA[4:6...]]        # ν aΔt
    X.x[2].x[1] .= c[SA[1:3...]]        # ρ vΔt
    X.x[2].x[2] .= c[10]                # Δt
    return X
end

function LieGroups.get_coordinates_lie(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        X::ArrayPartition{
            T, Tuple{
                ArrayPartition{T, Tuple{SMatrix{3, 3, T, 9}, SVector{3, T}}},
                ArrayPartition{T, Tuple{SVector{3, T}, SVector{1, T}}},
            },
        },
        ::DefaultLieAlgebraOrthogonalBasis
    ) where {T <: Real}
    return SVector{10, T}(
        X.x[2].x[1][1],   # ρ vΔt
        X.x[2].x[1][2],
        X.x[2].x[1][3],
        X.x[1].x[2][1],   # ν aΔt
        X.x[1].x[2][2],
        X.x[1].x[2][3],
        X.x[1].x[1][3, 2], # θ⃗ₓ[3,2]
        X.x[1].x[1][1, 3], # θ⃗ₓ[1,3]
        X.x[1].x[1][2, 1], # θ⃗ₓ[2,1]
        X.x[2].x[2][],    # Δt
    )
end

function LieGroups.get_coordinates_lie!(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        c,
        X::ArrayPartition,
        ::DefaultLieAlgebraOrthogonalBasis
    )
    c[1] = X.x[2].x[1][1]   # ρ vΔt
    c[2] = X.x[2].x[1][2]
    c[3] = X.x[2].x[1][3]
    c[4] = X.x[1].x[2][1]   # ν aΔt
    c[5] = X.x[1].x[2][2]
    c[6] = X.x[1].x[2][3]
    c[7] = X.x[1].x[1][3, 2] # θ⃗ₓ[3,2]
    c[8] = X.x[1].x[1][1, 3] # θ⃗ₓ[1,3]
    c[9] = X.x[1].x[1][2, 1] # θ⃗ₓ[2,1]
    c[10] = X.x[2].x[2][]   # Δt
    return c
end

function LieGroups.jacobian_exp!(
        ::LieGroups.SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        J::AbstractMatrix,
        X::ArrayPartition,
        ::DefaultLieAlgebraOrthogonalBasis,
    )
    Ω = X.x[1].x[1]
    ν = X.x[1].x[2]
    ρ = X.x[2].x[1]
    ι = X.x[2].x[2][]
    ω = [Ω[3, 2], Ω[1, 3], Ω[2, 1]]
    return LieGroups._jacobian_exp_left_SGal3!(J, -ρ, -ν, -ω, -ι)
end
