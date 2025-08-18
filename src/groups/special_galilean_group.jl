"""
    RotationBoostAction

The group action of the semidirect product of spatial rotations and velocity boosts (SO(n) ⋉ ℝⁿ) on the space of events (position, time).
"""
struct RotationBoostAction <: AbstractLeftGroupActionType end

# the action of the rotation-boost semidirect product group on the events (x,t) group
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
    # The group of spatial rotations and velocity boosts SO(n) ⋉ ℝⁿ
    LeftSemidirectProductGroupOperation{
        MatrixMultiplicationGroupOperation, AdditionGroupOperation, LeftGroupOperationAction,
    },
    ProductGroupOperation{Tuple{AdditionGroupOperation, AdditionGroupOperation}},
    RotationBoostAction,
}

"""
    SpecialGalileanGroup{T}

The special Galilean group SGal(3), a 10-dimensional Lie group of spacetime transformations preserving spatial distances and absolute time intervals.
It consists of spatial rotations, velocity boosts, and spacetime translations [Kelly:2025](@cite).
"""
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
    EventsGroup{T}

An event is a point in Galilean spacetime, specified by three spatial coordinates and one temporal coordinate and 
denoted by a tuple (x, t) ∈ ℝ³ × ℝ.
See [Kelly:2025; section 4.1](@cite).
"""
EventsGroup{T} = LieGroup{
    ℝ,
    ProductGroupOperation{Tuple{AdditionGroupOperation, AdditionGroupOperation}},
    ProductManifold{
        ℝ, Tuple{Euclidean{T, ℝ}, Euclidean{ManifoldsBase.TypeParameter{Tuple{1}}, ℝ}},
    },
}

"""
    default_left_action(::SpecialEuclideanGroup, ::EventsGroup)

Return the default left group action for SE(n) acting on events (position, time), i.e., the RotationBoostAction.
"""
function default_left_action(::SpecialEuclideanGroup, ::EventsGroup)
    return RotationBoostAction()
end

"""
    SpecialGalileanGroup(n::Int)

Construct the special Galilean group SGal(n) as a nested semidirect product:
    ``(SO(n) ⋉ ℝⁿ) ⋉ (ℝⁿ × ℝ)``
where SO(n) are spatial rotations, ℝⁿ are velocity boosts, and (ℝⁿ × ℝ) are spacetime translations.
The affine representation of the group is given by the matrix:
```math
SGal(3) = [R v p;
           0 1 t;
           0 0 1] ⊂ ℝ⁵ˣ⁵
```
ArrayPartition representation
``SGal(3) = ((R, v), (p, t))``
[Kelly:2025](@cite)
"""
function SpecialGalileanGroup(n::Int)
    return (SpecialOrthogonalGroup(n) ⋉ TranslationGroup(n)) ⋉ (TranslationGroup(n) × TranslationGroup(1))
end
