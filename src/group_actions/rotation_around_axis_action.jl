@doc raw"""
    RotationAroundAxisAction(axis::AbstractVector)

Space of actions of the circle group [`CircleGroup`](@ref) on $ℝ^3$ around given `axis`.
"""
struct RotationAroundAxisAction{TA <: AbstractVector} <: AbstractLeftGroupActionType
    axis::TA
end

function apply(A::GroupAction{<:RotationAroundAxisAction}, θ, p)
    sθ, cθ = sincos(θ)
    apd = dot(A.type.axis, p)
    return p .* cθ .+ cross(A.type.axis, p) .* sθ .+ A.type.axis .* apd .* (1 - cθ)
end

@doc raw"""
    apply(A::GroupAction{<:RotationAroundAxisAction}, θ, p)

Rotate point `p` from [`Euclidean`](@extref `Manifolds.Euclidean`) manifold around axis `A.axis` by angle `θ`.
The formula reads

````math
p_{rot} = (\cos(θ))p + (k×p) \sin(θ) + k (k⋅p) (1-\cos(θ)),
````

where ``k`` is the vector `A.type.axis` and `⋅` is the dot product.
"""
apply(A::GroupAction{<:RotationAroundAxisAction}, θ, p)

function apply(A::GroupAction{<:RotationAroundAxisAction}, θ::AbstractArray, p)
    # this method is here to make sure that θ represented by 1-element vectors works
    return apply(A, θ[], p)
end

function apply!(A::GroupAction{<:RotationAroundAxisAction}, q, θ, p)
    return copyto!(q, apply(A, θ, p))
end
