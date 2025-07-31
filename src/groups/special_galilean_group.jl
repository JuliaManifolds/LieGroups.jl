# using LieGroups
# using ManifoldsBase
# using ManifoldsBase: submanifold_component
# using LieGroups:
#     ℝ, Rotations, Euclidean, ProductManifold, LeftSemidirectProductGroupOperation
# import LieGroups: default_left_action

# using RecursiveArrayTools: ArrayPartition

# using StaticArrays, LinearAlgebra

using Rotations: RotationVec
using RecursiveArrayTools: ArrayPartition

##
struct RotationBoostAction <: AbstractLeftGroupActionType end

# the action of the rotation-boost semidirect product group on the events (x,t) group
function LieGroups.apply!(A::GroupAction{RotationBoostAction}, k, g, h)
    G = A.group
    R = submanifold_component(G, g, Val(1))
    v = submanifold_component(G, g, Val(2))
    
    H = A.manifold
    p = submanifold_component(H, h, Val(1))
    t = submanifold_component(H, h, Val(2))[1]

    k_p = submanifold_component(H, k, Val(1))
    k_p .= v * t + R * p
    # NOTE k_t unchanged
    return k
end

const LeftSpecialGalileanGroupOperation = LeftSemidirectProductGroupOperation{
    # The group of spatial rotations and velocity boosts SO(n) ⋉ ℝⁿ
    LeftSemidirectProductGroupOperation{
        MatrixMultiplicationGroupOperation,AdditionGroupOperation,LeftGroupOperationAction
    },
    ProductGroupOperation{Tuple{AdditionGroupOperation,AdditionGroupOperation}},
    RotationBoostAction,
}

SpecialGalileanGroup{T} = LieGroup{
    ℝ,
    <:LeftSpecialGalileanGroupOperation,
    <:ProductManifold{
        ℝ,
        Tuple{
            <:ProductManifold{ℝ,Tuple{<:Rotations{T},<:Euclidean{T,ℝ}}},
            <:ProductManifold{
                ℝ,
                Tuple{
                    <:Euclidean{T,ℝ},<:Euclidean{ManifoldsBase.TypeParameter{Tuple{1}},ℝ}
                },
            },
        },
    },
}

# EventsGroup is a Product group of Translation(n) x Time
# see eq 1. in https://arxiv.org/pdf/2312.07555
EventsGroup{T} = LieGroup{
    ℝ,
    ProductGroupOperation{Tuple{AdditionGroupOperation,AdditionGroupOperation}},
    ProductManifold{
        ℝ,Tuple{Euclidean{T,ℝ},Euclidean{ManifoldsBase.TypeParameter{Tuple{1}},ℝ}}
    },
}

function default_left_action(::SpecialEuclideanGroup, ::EventsGroup)
    return RotationBoostAction()
end

"""
    SpecialGalileanGroup

References: 
- https://hal.science/hal-02183498/document
- TODO new reference: https://arxiv.org/pdf/2312.07555
- TODO new reference: https://arxiv.org/pdf/2409.14276

Affine representation 
Δ = [ΔR Δv Δp;
     0   1 Δt;
     0   0  1] ⊂ ℝ⁵ˣ⁵

ArrayPartition representation
Δ = ((ΔR, Δv), (Δp, Δt))
"""
function SpecialGalileanGroup(n::Int)
    return (SpecialOrthogonalGroup(n) ⋉ TranslationGroup(n)) ⋉
           (TranslationGroup(n) × TranslationGroup(1))
end

#TODO don't know if something similar already exists in LieGroups.jl
function _skew(v::AbstractVector{T}) where T<:Real
    return SMatrix{3,3,T}(
            0,
         v[3],
        -v[2],
        -v[3],  
            0,  
         v[1],
         v[2],
        -v[1],
            0
    )
end

function _Q(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return SMatrix{3,3,T}(I)
    else
        u = θ⃗ / θ
        sθ, cθ = sincos(θ)
        uₓ = _skew(u)
        return SMatrix{3,3,T}(I) + (1 - cθ) / θ * uₓ + (θ - sθ) / θ * uₓ^2
    end
end

function _P(θ⃗)
    T = eltype(θ⃗)
    θ = norm(θ⃗)
    if θ ≈ 0
        return 1 / 2 * SMatrix{3,3,T}(I)
    else
        u = θ⃗ / θ
        sθ, cθ = sincos(θ)
        uₓ = _skew(u)
        return 1 / 2 * SMatrix{3,3,T}(I) +
               (θ - sθ) / θ^2 * uₓ +
               (cθ + 1 / 2 * θ^2 - 1) / θ^2 * uₓ^2
    end
end

function LieGroups.exp!(G::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, q::ArrayPartition, X::ArrayPartition)
    
    θ⃗ₓ = X.x[1].x[1] # ωΔt

    ν = X.x[1].x[2]  # aΔt
    ρ = X.x[2].x[1]  # vΔt

    Δt = X.x[2].x[2][1]

    # ωΔt = vee(θ⃗ₓ)
    θ⃗ = SA[θ⃗ₓ[3, 2]; θ⃗ₓ[1, 3]; θ⃗ₓ[2, 1]]

    P = _P(θ⃗)
    Q = _Q(θ⃗)

    M_SO3 = SpecialOrthogonalGroup(3)
    q.x[1].x[1] .= exp(M_SO3, θ⃗ₓ)
    q.x[1].x[2] .= Q * ν
    q.x[2].x[1] .= Q * ρ + P * ν * Δt
    q.x[2].x[2] .= Δt

    return q
end

function LieGroups.log!(
    M::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    X::ArrayPartition,
    p::ArrayPartition
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

    
    X.x[1].x[1] .= log(SpecialOrthogonalGroup(3), ΔR) # θ⃗ₓ
    X.x[1].x[2] .= iQ * Δv # ν aΔt
    X.x[2].x[1] .= iQ * (Δp - P * iQ * Δv * Δt) # ρ vΔt
    X.x[2].x[2] .= Δt
    return X
end

function LieGroups.identity_element(
    ::SpecialGalileanGroup{ManifoldsBase.TypeParameter{Tuple{N}}},
    ::Type{<:StaticArray}
) where {N}
    return ArrayPartition(
        ArrayPartition(
            SMatrix{3,3,Float64}(I), # ΔR
            @SVector(zeros(3)),      # Δv
        ),
        ArrayPartition(
            @SVector(zeros(3)),      # Δp
            @SVector([0.0]),         # Δt
        )
    )
end

function LieGroups.inv(G::SpecialGalileanGroup, p::ArrayPartition)
    ΔR = p.x[1].x[1]
    Δv = p.x[1].x[2]
    Δp = p.x[2].x[1]
    Δt = p.x[2].x[2]

    return ArrayPartition(
        ArrayPartition(
            ΔR',
            -ΔR' * Δv
        ),
        ArrayPartition(
            -ΔR' * (Δp - Δv * Δt[1]),
            -Δt
        )
    )
end

function LieGroups.compose(G::SpecialGalileanGroup, p::ArrayPartition, q::ArrayPartition)
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
            ΔR * δR,
            Δv + ΔR * δv,
        ),
        ArrayPartition(
            Δp + Δv * δt[1] + ΔR * δp,
            Δt + δt
        )
    )
end
