"""
    OrthogonalGroup{T}

The orthogonal group ``$(_math(:O))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`).

# Constructor
    OrthogonalGroup(n; kwargs...)

Generate  orthogonal group ``$(_math(:O))(n)``.
All keyword arguments in `kwargs...` are passed on to [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) as well.
"""
const OrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ,MatrixMultiplicationGroupOperation,Manifolds.OrthogonalMatrices{T}
}

function OrthogonalGroup(n; kwargs...)
    R = Manifolds.OrthogonalMatrices(n; kwargs...)
    return OrthogonalGroup{typeof(R).parameters[1]}(R, MatrixMultiplicationGroupOperation())
end

#
#
# Generic special cases for O(n) and SO(n)

@doc raw"""
    angles_4d_skew_sym_matrix(A)

The Lie algebra of [`Rotations(4)`](@ref) in ``ℝ^{4×4}``, ``𝔰𝔬(4)``, consists of ``4×4``
skew-symmetric matrices. The unique imaginary components of their eigenvalues are the
angles of the two plane rotations. This function computes these more efficiently than
`eigvals`.

By convention, the returned values are sorted in decreasing order
(corresponding to the same ordering of _angles_ as
[`cos_angles_4d_rotation_matrix`](@extref Manifolds.cos_angles_4d_rotation_matrix)).
"""
function angles_4d_skew_sym_matrix(A)
    @assert size(A) == (4, 4)
    @inbounds begin
        halfb = (A[1, 2]^2 + A[1, 3]^2 + A[2, 3]^2 + A[1, 4]^2 + A[2, 4]^2 + A[3, 4]^2) / 2
        c = (A[1, 2] * A[3, 4] - A[1, 3] * A[2, 4] + A[1, 4] * A[2, 3])^2
    end
    sqrtdisc = sqrt(halfb^2 - c)
    return sqrt(halfb + sqrtdisc), sqrt(halfb - sqrtdisc)
end

function ManifoldsBase.exp!(
    ::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    q,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
    @assert size(X) == (2, 2)
    @inbounds θ = (X[2, 1] - X[1, 2]) / 2
    sinθ, cosθ = sincos(θ)
    @inbounds begin
        q[1, 1] = cosθ
        q[2, 1] = sinθ
        q[1, 2] = -sinθ
        q[2, 2] = cosθ
    end
    return q
end
function ManifoldsBase.exp!(
    ::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{3}}},
    q,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
    θ = norm(X) / sqrt(2)
    if θ ≈ 0
        a = 1 - θ^2 / 6
        b = θ / 2
    else
        a = sin(θ) / θ
        b = (1 - cos(θ)) / θ^2
    end
    copyto!(q, I)
    q .+= a .* X
    mul!(q, X, X, b, true)
    return q
end
function ManifoldsBase.exp!(
    ::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{4}}},
    q,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
    T = eltype(X)
    α, β = angles_4d_skew_sym_matrix(X)
    sinα, cosα = sincos(α)
    sinβ, cosβ = sincos(β)
    α² = α^2
    β² = β^2
    Δ = β² - α²
    if !isapprox(Δ, 0; atol=1e-6)  # Case α > β ≥ 0
        sincα = sinα / α
        sincβ = β == 0 ? one(T) : sinβ / β
        a₀ = (β² * cosα - α² * cosβ) / Δ
        a₁ = (β² * sincα - α² * sincβ) / Δ
        a₂ = (cosα - cosβ) / Δ
        a₃ = (sincα - sincβ) / Δ
    elseif α == 0 # Case α = β = 0
        a₀ = a₁ = one(T)
        a₂ = inv(T(2))
        a₃ = inv(T(6))
    else  # Case α ⪆ β ≥ 0, α ≠ 0
        sincα = sinα / α
        r = β / α
        c = 1 / (1 + r)
        d = α * (α - β) / 2
        if α < 1e-2
            e = evalpoly(α², (inv(T(3)), inv(T(-30)), inv(T(840)), inv(T(-45360))))
        else
            e = (sincα - cosα) / α²
        end
        a₀ = (α * sinα + (1 + r - d) * cosα) * c
        a₁ = ((3 - d) * sincα - (2 - r) * cosα) * c
        a₂ = (sincα - (1 - r) / 2 * cosα) * c
        a₃ = (e + (1 - r) * (e - sincα / 2)) * c
    end

    X² = X * X
    X³ = X² * X
    q = a₀ * LinearAlgebra.I + a₁ .* X .+ a₂ .* X² .+ a₃ .* X³
    return q
end

function ManifoldsBase.log!(
    G::CommonUnitarySubGroups{ℝ,<:Any},
    X::AbstractMatrix,
    ::Identity{MatrixMultiplicationGroupOperation},
    q::AbstractMatrix,
)
    log_safe!(X, q)
    return project!(G, X, Identity(G), X)
end
function ManifoldsBase.log!(
    ::CommonUnitarySubGroups{ℝ,<:Any},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
)
    fill!(X, 0)
    return X
end
function ManifoldsBase.log!(
    G::CommonUnitarySubGroups{
        ℝ,
        ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.ManifoldsBase.TypeParameter{
            Tuple{2}
        },
    },
    X::AbstractMatrix,
    ::Identity{MatrixMultiplicationGroupOperation},
    q::AbstractMatrix,
)
    @assert size(q) == (2, 2)
    @inbounds θ = atan(q[2, 1], q[1, 1])
    return get_vector!(G, X, Identity(G), θ, DefaultOrthogonalBasis())
end
function ManifoldsBase.log!(
    G::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{3}}},
    X::AbstractMatrix,
    e::Identity{MatrixMultiplicationGroupOperation},
    q::AbstractMatrix,
)
    cosθ = (tr(q) - 1) / 2
    if cosθ ≈ -1
        eig = eigen_safe(q)
        ival = findfirst(λ -> isapprox(λ, 1), eig.values)
        inds = SVector{3}(1:3)
        ax = eig.vectors[inds, ival]
        return get_vector!(G, X, e, π * ax, DefaultOrthogonalBasis())
    end
    X .= q ./ usinc_from_cos(cosθ)
    # TODO: Is this necessary or just for numerical stability?
    return project!(G, X, e, X)
end
function ManifoldsBase.log!(
    G::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{4}}},
    X::AbstractMatrix,
    e::Identity{MatrixMultiplicationGroupOperation},
    q::AbstractMatrix,
)
    cosα, cosβ = Manifolds.cos_angles_4d_rotation_matrix(q)
    α = acos(clamp(cosα, -1, 1))
    β = acos(clamp(cosβ, -1, 1))
    if α ≈ 0 && β ≈ π
        A² = Symmetric((q - I) ./ 2)
        P = eigvecs(A²)
        E = similar(q)
        fill!(E, 0)
        @inbounds begin
            E[2, 1] = -β
            E[1, 2] = β
        end
        copyto!(X, P * E * transpose(P))
    else
        det(q) < 0 && throw(
            DomainError(
                "The Lie group logarithm is not defined for $q with a negative determinant ($(det(q)) < 0). Point `q` is in a different connected component of the manifold $G",
            ),
        )
        log_safe!(X, q)
    end
    # TODO: Is this necessary or just for numerical stability?
    return project!(G, X, Identity(G), X)
end

function Base.show(io::IO, G::OrthogonalGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "OrthogonalGroup($(size))")
end
