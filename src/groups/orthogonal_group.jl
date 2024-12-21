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

_doc_exp_O2_id = """
    exp(G::OrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, X)
    exp(G::SpecialOrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, X)
    exp!(G::OrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, g, X)
    exp!(G::SpecialOrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, g, X)

Compute the Lie group exponential function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(2)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(2)``.

Since the Lie algebra of both groups agrees and consist of the set of skew symmetric matrices,
these simplify for the case of ``2×2`` matrices to ``X=$(_tex(:pmatrix, "0 & -α", "α & 0"))``, for some ``α∈ℝ``.

Their exponential is

```math
$(_tex(:exp))_{$(_math(:G))}(X) =  $(_tex(:pmatrix, "$(_tex(:cos))(α) & -$(_tex(:sin))(α)", "$(_tex(:sin))(α) & $(_tex(:cos))(α)")).
```

This result can be computed in-place of `g`.

Note that since ``$(_math(:SO))(2)`` consists of two disjoint connected components and the exponential map is smooth,
the result ``g`` always lies in the connected component of the identity.
"""

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp(
    ::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp!(
    ::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    g,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
function ManifoldsBase.exp(
    G::CommonUnitarySubGroups{ManifoldsBase.ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    e::Identity{MatrixMultiplicationGroupOperation},
    X,
)
    g = ManifoldsBase.allocate_result(G, exp, X)
    exp!(G, g, e, X)
    return g
end

function ManifoldsBase.exp!(
    ::CommonUnitarySubGroups{ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    g,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
    @assert size(X) == (2, 2)
    @assert size(g) == (2, 2)
    @inbounds α = (X[2, 1] - X[1, 2]) / 2
    sinα, cosα = sincos(α)
    @inbounds begin
        g[1, 1] = cosα
        g[2, 1] = sinα
        g[1, 2] = -sinα
        g[2, 2] = cosα
    end
    return g
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

_doc_log_O2_id = """
    log(G::OrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, g)
    log(G::SpecialOrthogonalGroup{TypeParameter{Tuple{2}}}, ::Identity{MatrixMultiplicationGroupOperation}, g)
    log!(G::OrthogonalGroup{TypeParameter{Tuple{2}}}, X, ::Identity{MatrixMultiplicationGroupOperation}, g)
    log!(G::SpecialOrthogonalGroup{TypeParameter{Tuple{2}}}, X, ::Identity{MatrixMultiplicationGroupOperation}, g)

Compute the Lie group logarithm function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(2)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(2)``.

For the two-dimensional case, any rotation matrix ``g`` can be represented as ``$(_tex(:pmatrix, "$(_tex(:cos))(α) & -$(_tex(:sin))(α)", "$(_tex(:sin))(α) & $(_tex(:cos))(α)"))``.
For the [`SpecialOrthogonalGroup`](@ref), ``g`` might also include reflections.

The logarithm is then

```math
$(_tex(:log))_{$(_math(:G))}(g) =  $(_tex(:pmatrix, "0 & α &", "-α & 0")).
```

This result can be computed in-place of `X`

Note the logarithmic map is only locally around the identity uniquely determined.
Especially, since ``$(_math(:SO))(2)`` consists of two disjoint connected components and the exponential map is smooth,
for any ``g`` in the other component, the logarithmic map is defined, but not the inverse of the exponential map.
"""

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log(
    ::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log!(
    ::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

function ManifoldsBase.log(
    G::CommonUnitarySubGroups{ManifoldsBase.ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    e::Identity{MatrixMultiplicationGroupOperation},
    g,
)
    Y = ManifoldsBase.allocate_result(G, log, g)
    log!(G, Y, e, g)
    return Y
end
# Resolve an ambiguity compared to the general matrix multiplication definitions
function Base.log(
    G::CommonUnitarySubGroups{ManifoldsBase.ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    e::Identity{MatrixMultiplicationGroupOperation},
    ::Identity{MatrixMultiplicationGroupOperation},
)
    return zero_vector(G, e)
end

function ManifoldsBase.log!(
    ::CommonUnitarySubGroups{ManifoldsBase.ℝ,ManifoldsBase.TypeParameter{Tuple{2}}},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)
    @assert size(X) == (2, 2)
    @assert size(g) == (2, 2)
    @inbounds α = atan(g[2, 1], g[1, 1])
    @inbounds begin
        X[1, 1] = 0
        X[2, 1] = α
        X[1, 2] = -α
        X[2, 2] = 0
    end
    return X
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
