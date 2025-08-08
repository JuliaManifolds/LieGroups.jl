
"""
    RotationBoostAction

The group action of the semidirect product of spatial rotations and velocity boosts (SO(n) ⋉ ℝⁿ) on the space of events (position, time).
"""
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

    k_t = submanifold_component(H, k, Val(2))
    k_t .= t
    return k
end

const LeftSpecialGalileanGroupOperation = LeftSemidirectProductGroupOperation{
    # The group of spatial rotations and velocity boosts SO(n) ⋉ ℝⁿ
    LeftSemidirectProductGroupOperation{
        MatrixMultiplicationGroupOperation, AdditionGroupOperation, LeftGroupOperationAction,
    },
    ProductGroupOperation{Tuple{AdditionGroupOperation, AdditionGroupOperation}},
    RotationBoostAction,
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

# EventsGroup is a Product group of Translation(n) x Time
# see eq 1. in https://arxiv.org/pdf/2312.07555
EventsGroup{T} = LieGroup{
    ℝ,
    ProductGroupOperation{Tuple{AdditionGroupOperation, AdditionGroupOperation}},
    ProductManifold{
        ℝ, Tuple{Euclidean{T, ℝ}, Euclidean{ManifoldsBase.TypeParameter{Tuple{1}}, ℝ}},
    },
}

function default_left_action(::SpecialEuclideanGroup, ::EventsGroup)
    return RotationBoostAction()
end

"""
    SpecialGalileanGroup(n::Int)

Construct the special Galilean group SGal(n) as a nested semidirect product:
    (SO(n) ⋉ ℝⁿ) ⋉ (ℝⁿ × ℝ)
where SO(n) are spatial rotations, ℝⁿ are velocity boosts, and (ℝⁿ × ℝ) are spacetime translations.
Affine representation 
Δ = [ΔR Δv Δp;
     0   1 Δt;
     0   0  1] ⊂ ℝ⁵ˣ⁵

ArrayPartition representation
Δ = ((ΔR, Δv), (Δp, Δt))
See Eq. (SGal3_definition) and Section "The Matrix Representation of SGal(3)".
References: 
- https://hal.science/hal-02183498/document
- https://arxiv.org/pdf/2312.07555
- https://arxiv.org/pdf/2409.14276
"""
function SpecialGalileanGroup(n::Int)
    return (SpecialOrthogonalGroup(n) ⋉ TranslationGroup(n)) ⋉ (TranslationGroup(n) × TranslationGroup(1))
end
