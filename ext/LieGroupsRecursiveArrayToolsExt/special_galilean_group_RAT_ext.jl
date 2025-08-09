using Rotations: RotationVec
using LieGroups: SpecialGalileanGroup
using StaticArrays
using LinearAlgebra

function _skew(v::AbstractVector{T}) where {T <: Real}
    return SMatrix{3, 3, T}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
end

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
        q::ArrayPartition,
        X::ArrayPartition,
    )
    θ⃗ₓ = X.x[1].x[1] # ωΔt
    ν = X.x[1].x[2]  # aΔt
    ρ = X.x[2].x[1]  # vΔt

    Δt = X.x[2].x[2][1]

    # ωΔt = vee(θ⃗ₓ)
    θ⃗ = SA[θ⃗ₓ[3, 2]; θ⃗ₓ[1, 3]; θ⃗ₓ[2, 1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonalGroup(3)
    exp!(M_SO3, q.x[1].x[1], θ⃗ₓ)
    q.x[1].x[2] .= Q * ν
    q.x[2].x[1] .= Q * ρ + P * ν * Δt #FIXME allocates if X is not a static array
    q.x[2].x[2] .= Δt

    return q
end

function LieGroups.exp(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        X::ArrayPartition{T}
    ) where {T <: Real}

    θ⃗ₓ = X.x[1].x[1] # ωΔt
    ν = X.x[1].x[2]  # aΔt
    ρ = X.x[2].x[1]  # vΔt
    Δt = X.x[2].x[2][1]

    # ωΔt = vee(θ⃗ₓ)
    θ⃗ = SA[θ⃗ₓ[3, 2]; θ⃗ₓ[1, 3]; θ⃗ₓ[2, 1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonalGroup(3)
    q = ArrayPartition(
        ArrayPartition(
            exp(M_SO3, θ⃗ₓ),
            Q * ν
        ),
        ArrayPartition(
            Q * ρ + P * ν * Δt,
            [Δt]
        )
    )
    return q
end

function LieGroups.log!(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        X::ArrayPartition,
        p::ArrayPartition,
    )
    ΔR = p.x[1].x[1]
    Δv = p.x[1].x[2]
    Δp = p.x[2].x[1]
    Δt = p.x[2].x[2][1]

    #FIXME Rotations is not a dependency, find alternative for RotationVec
    Rv = RotationVec(ΔR)
    θ⃗ = SA[Rv.sx, Rv.sy, Rv.sz]

    P = _P(θ⃗)
    Q = _Q(θ⃗)
    iQ = inv(Q)

    log!(SpecialOrthogonalGroup(3), X.x[1].x[1], ΔR) # θ⃗ₓ # FIXME allocates
    X.x[1].x[2] .= iQ * Δv # ν aΔt
    X.x[2].x[1] .= iQ * (Δp - P * iQ * Δv * Δt) # ρ vΔt
    X.x[2].x[2] .= Δt
    return X
end

function LieGroups.log(
        ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        p::ArrayPartition
    )
    ΔR = p.x[1].x[1]
    Δv = p.x[1].x[2]
    Δp = p.x[2].x[1]
    Δt = p.x[2].x[2][1]

    Rv = RotationVec(ΔR)
    θ⃗ = SA[Rv.sx, Rv.sy, Rv.sz]

    P = _P(θ⃗)
    Q = _Q(θ⃗)
    iQ = inv(Q)
    return ArrayPartition(
        ArrayPartition(
            log(SpecialOrthogonalGroup(3), ΔR), # θ⃗ₓ
            iQ * Δv # ν aΔt
        ),
        ArrayPartition(
            iQ * (Δp - P * iQ * Δv * Δt), # ρ vΔt
            [Δt]
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

function LieGroups.inv(::SpecialGalileanGroup, p::ArrayPartition)
    ΔR = p.x[1].x[1]
    Δv = p.x[1].x[2]
    Δp = p.x[2].x[1]
    Δt = p.x[2].x[2]

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

function LieGroups.compose(::SpecialGalileanGroup, p::ArrayPartition, q::ArrayPartition)
    ΔR = p.x[1].x[1]
    Δv = p.x[1].x[2]
    Δp = p.x[2].x[1]
    Δt = p.x[2].x[2]

    δR = q.x[1].x[1]
    δv = q.x[1].x[2]
    δp = q.x[2].x[1]
    δt = q.x[2].x[2]

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

#NOTE The ordering is hat and vee are different that the default in LieGroups.jl
# see eq 14 for the basis used for hat and vee

function LieGroups.get_vector_lie(
        ::typeof(LieAlgebra(SpecialGalileanGroup(3))),
        Xⁱ::SVector{10, T},
        ::DefaultLieAlgebraOrthogonalBasis,
        ::Type{<:ArrayPartition{<:Real}} # = ArrayPartition{T}
    ) where {T <: Real}
    return ArrayPartition(
        ArrayPartition(
            _skew(Xⁱ[SA[7:9...]]), # θ ωΔt
            Xⁱ[SA[4:6...]]         # ν aΔt
        ),
        ArrayPartition(
            Xⁱ[SA[1:3...]],        # ρ vΔt
            Xⁱ[SA[10]],            # Δt
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
        Xⁱ,
        ::DefaultLieAlgebraOrthogonalBasis
    ) where {T <: Real}
    X.x[1].x[1] .= _skew(Xⁱ[SA[7:9...]]) # θ ωΔt
    X.x[1].x[2] .= Xⁱ[SA[4:6...]]        # ν aΔt
    X.x[2].x[1] .= Xⁱ[SA[1:3...]]        # ρ vΔt
    X.x[2].x[2] .= Xⁱ[10]                # Δt
    return X
end
# 𝔤::LieAlgebra, c, B::DefaultLieAlgebraOrthogonalBasis, T::Type

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
        X.x[2].x[1][1],   # ν aΔt
        X.x[2].x[1][2],
        X.x[2].x[1][3],
        X.x[1].x[2][1],   # ρ vΔt
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
        Xⁱ,
        X::ArrayPartition,
        ::DefaultLieAlgebraOrthogonalBasis
    )
    Xⁱ[1] = X.x[2].x[1][1]   # ν aΔt
    Xⁱ[2] = X.x[2].x[1][2]
    Xⁱ[3] = X.x[2].x[1][3]
    Xⁱ[4] = X.x[1].x[2][1]   # ρ vΔt
    Xⁱ[5] = X.x[1].x[2][2]
    Xⁱ[6] = X.x[1].x[2][3]
    Xⁱ[7] = X.x[1].x[1][3, 2] # θ⃗ₓ[3,2]
    Xⁱ[8] = X.x[1].x[1][1, 3] # θ⃗ₓ[1,3]
    Xⁱ[9] = X.x[1].x[1][2, 1] # θ⃗ₓ[2,1]
    Xⁱ[10] = X.x[2].x[2][]    # Δt
    return Xⁱ
end
