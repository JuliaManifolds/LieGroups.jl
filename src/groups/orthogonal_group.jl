"""
    OrthogonalGroup{T}

The orthogonal group ``$(_math(:O))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`).

# Constructor
    OrthogonalGroup(n::Int; kwargs...)

Generate  orthogonal group ``$(_math(:O))(n)``.
All keyword arguments in `kwargs...` are passed on to [`OrthogonalMatrices`](@extref `Manifolds.OrthogonalMatrices`) as well.
"""
const OrthogonalGroup{T} = LieGroup{
    ℝ, MatrixMultiplicationGroupOperation, OrthogonalMatrices{T},
}

function OrthogonalGroup(n::Int; kwargs...)
    R = OrthogonalMatrices(n; kwargs...)
    return OrthogonalGroup{typeof(R).parameters[2]}(R, MatrixMultiplicationGroupOperation())
end

const OrthogonalLieAlgebra{T} = LieAlgebra{
    ℝ,
    MatrixMultiplicationGroupOperation,
    LieGroup{ℝ, MatrixMultiplicationGroupOperation, OrthogonalMatrices{T}},
}

#
#
# Generic special cases for O(n) and SO(n)

_doc_exp_O2_id = """
    exp(G, X)
    exp!(G, g, X)

Compute the Lie group exponential function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(2)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(2)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

Since the Lie algebra of both groups agrees and consist of the set of skew symmetric matrices,
these simplify for the case of ``2×2`` matrices to ``X=$(_tex(:pmatrix, "0 & -α", "α & 0"))``, for some ``α∈ℝ``.

Their exponential is

```math
$(_tex(:exp))_{$(_math(:G))}(X) =  $(_tex(:pmatrix, "$(_tex(:cos))(α) & -$(_tex(:sin))(α)", "$(_tex(:sin))(α) & $(_tex(:cos))(α)")).
```

This result can be computed in-place of `g`.

Note that since ``$(_math(:O))(2)`` consists of two disjoint connected components and the exponential map is smooth,
the result ``g`` always lies in the connected component of the identity.
"""

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X)

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g, X)

_doc_exp_O3_id = """
    exp(G, X)
    exp!(G, g, X)

Compute the Lie group exponential function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(3)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(3)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

Since the Lie algebra of both groups agrees and consist of the set of skew symmetric matrices,
the ``3×3`` skew symmetric matrices are of the form

```math
    X = $(_tex(:pmatrix, "0 & -c & b", "c & 0 & -a", "-b & a & 0")),
```
for some ``a, b, c ∈ ℝ``. To compute the exponential, the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Olinde_Rodrigues)
can be used. With ``α = $(_tex(:sqrt, "a^2+b^2+c^2")) = $(_tex(:frac, "1", _tex(:sqrt, "2")))$(_tex(:norm, "X"))``
we obtain for ``α ≠ 0``

```math
$(_tex(:exp))_{$(_math(:G))}(X) = I_3 + $(_tex(:frac, _tex(:sin, "α"), "α"))X + $(_tex(:frac, "(1 - $(_tex(:cos, "α")))", "α^2"))X^2,
```

and $(_tex(:exp))_{$(_math(:G))}(X) = I_3`` otherwise.

This result can be computed in-place of `g`.

Note that since ``$(_math(:SO))(3)`` consists of two disjoint connected components and the exponential map is smooth,
the result ``g`` always lies in the connected component of the identity.
"""

@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X)

@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g, X)

_doc_exp_O4_id = """
    exp(G, e, X)
    exp!(G, e, g, X)

Compute the Lie group exponential function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(4)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(4)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

Similar to the ``3×3`` case, an efficient computation is provided,
adapted from [GallierXu:2002](@cite), [AndricaRohan:2013](@cite) with a few numerical stabilisations.
"""

@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X)

@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g, X)

function ManifoldsBase.exp!(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}}, g, X
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
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}}, g, X
    )
    θ = norm(X) / sqrt(2)
    if θ ≈ 0
        a = 1 - θ^2 / 6
        b = θ / 2
    else
        a = sin(θ) / θ
        b = (1 - cos(θ)) / θ^2
    end
    copyto!(g, LinearAlgebra.I)
    g .+= a .* X
    mul!(g, X, X, b, true)
    return g
end
function ManifoldsBase.exp!(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{4}}}, g, X
    )
    T = eltype(X)
    α, β = angles_4d_skew_sym_matrix(X)
    sinα, cosα = sincos(α)
    sinβ, cosβ = sincos(β)
    α² = α^2
    β² = β^2
    Δ = β² - α²
    if !isapprox(Δ, 0; atol = 1.0e-6)  # Case α > β ≥ 0
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
        if α < 1.0e-2
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
    g .= a₀ * LinearAlgebra.I + a₁ .* X .+ a₂ .* X² .+ a₃ .* X³
    return g
end

function ManifoldsBase.exp(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        p::SMatrix,
        X::SMatrix{2, 2, T},
    ) where {T}
    # Use the default orthogonal basis and not vee
    θ = X[2] * sqrt(T(2))
    sinθ, cosθ = sincos(θ)
    return p * SA[cosθ -sinθ; sinθ cosθ]
end
function ManifoldsBase.exp(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        X::SMatrix{2, 2, T},
    ) where {T}
    θ = X[2] * sqrt(T(2))
    sinθ, cosθ = sincos(θ)
    return SA[cosθ -sinθ; sinθ cosθ]
end
function ManifoldsBase.exp(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}},
        p::SMatrix,
        X::SMatrix,
    )
    return exp(G.manifold, p, X)
end
function ManifoldsBase.exp(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}},
        X::SMatrix,
    )
    return exp(G.manifold, SMatrix{3, 3, eltype(X)}(I), X)
end

function ManifoldsBase.exp_fused(
        M::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        p::SMatrix,
        X::SMatrix,
        t::Real,
    )
    return exp(M, p, t * X)
end

_doc_get_coordinates_On = """
    get_coordinates(𝔤::OrthogonalLieAlgebra, X, ::DefaultLieAlgebraOrthogonalBasis)
    get_coordinates(G::SpecialOrthogonalLieAlgebra, X, ::DefaultLieAlgebraOrthogonalBasis)
    get_coordinates!(G::OrthogonalLieAlgebra, c, X ::DefaultLieAlgebraOrthogonalBasis)
    get_coordinates!(G::SpecialOrthogonalLieAlgebra, c, X ::DefaultLieAlgebraOrthogonalBasis)

Compute the vector of coordinates ``c ∈ ℝ^d`` from the Lie algebra tangent vector ``X ∈ 𝔬(n)``
of the [`OrthogonalGroup`](@ref) `O(n)` in the [`DefaultLieAlgebraOrthogonalBasis`](@ref),
where ``d`` is the dimension of the Lie algebra. This is also the version used in [`vee`](@ref).

For ``O(2)`` there is only one coefficient ``α`` in the basis ``$(_tex(:pmatrix, "0 & -α", "α & 0"))``,
which is returned as ``c = (α)^$(_tex(:transp))``.

A usual basis representation of ``𝔬(3) is given by
```math
    X = $(_tex(:pmatrix, "0 & -γ & β", "γ & 0 & -α", "-β & α & 0")),
```
hence the coordinate vector is ``c = (α, β, γ)^$(_tex(:transp)) ∈ ℝ^3``.

For `n ≥ 4`` the lower triangular part is added to ``c`` row-wise.
"""

@doc "$(_doc_get_coordinates_On)"
get_coordinates(G::OrthogonalLieAlgebra, X, ::DefaultLieAlgebraOrthogonalBasis)

@doc "$(_doc_get_coordinates_On)"
get_coordinates!(G::OrthogonalLieAlgebra, c, X, ::DefaultLieAlgebraOrthogonalBasis)

function get_coordinates_lie!(
        ::CUSA, c, X, ::DefaultLieAlgebraOrthogonalBasis{ℝ}
    ) where {
        CUSA <: CommonUnitarySubAlgebra{ManifoldsBase.ℝ, <:ManifoldsBase.TypeParameter{Tuple{2}}},
    }
    @assert size(X) == (2, 2)
    @assert size(c) == (1,)
    c[1] = X[2, 1]
    return c
end
function get_coordinates_lie!(
        G::CommonUnitarySubAlgebra{ManifoldsBase.ℝ, <:ManifoldsBase.TypeParameter{Tuple{n}}},
        c,
        X,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    ) where {n}
    @assert size(X) == (n, n)
    @assert length(c) == manifold_dimension(G)
    @assert n > 2
    _get_coordinates_lie_On!(c, X)
    return c
end
function get_coordinates_lie!(
        𝔤::CommonUnitarySubAlgebra{ManifoldsBase.ℝ}, c, X, ::DefaultLieAlgebraOrthogonalBasis{ℝ}
    )
    G = 𝔤.manifold
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    @assert length(c) == manifold_dimension(𝔤)
    @assert size(X) == (n, n)
    if n == 2
        c[1] = X[2, 1]
    else
        _get_coordinates_lie_On!(c, X)
    end
    return c
end

function _get_coordinates_lie_On!(c, X)
    n = size(X, 1)
    @inbounds begin
        c[1] = X[3, 2]
        c[2] = X[1, 3]
        c[3] = X[2, 1]
        k = 4
        for i in 4:n, j in 1:(i - 1)
            c[k] = X[i, j]
            k += 1
        end
    end
    return c
end

function ManifoldsBase.get_coordinates(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        p::SMatrix,
        X::SMatrix,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    )
    return SA[X[2]]
end

function ManifoldsBase.get_coordinates(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}},
        p::SMatrix,
        X::SMatrix,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    )
    return SA[X[3, 2], X[1, 3], X[2, 1]]
end

_doc_get_vector_On = """
    get_vector(G::OrthogonalLieAlgebra, e, c, ::DefaultLieAlgebraOrthogonalBasis)
    get_vector(G::SpecialOrthogonalLieAlgebra, e, c, ::DefaultLieAlgebraOrthogonalBasis)
    get_vector!(G::OrthogonalLieAlgebra, X, e, c, ::DefaultLieAlgebraOrthogonalBasis)
    get_vector!(G::SpecialOrthogonalLieAlgebra, X, e, c, ::DefaultLieAlgebraOrthogonalBasis)

Compute the tangent vector ``X ∈ 𝔬(n)`` from a vector of coordinates ``c ∈ ℝ^d``.

Here ``d`` is the dimension of the Lie algebra of the [`OrthogonalGroup`](@ref) `O(n)`
and the coordinates are with respect to the [`DefaultLieAlgebraOrthogonalBasis`](@ref).
This is also the version used in [`hat`](@ref).

For ``O(2)`` there is only one coefficient ````c = (α)^$(_tex(:transp))`` and hence
``X = $(_tex(:pmatrix, "0 & -α", "α & 0"))`` is returned.

For ``n=3`` a usual representtion turns ``c = (α, β, γ)^$(_tex(:transp)) ∈ ℝ^3`` into
```math
    X = $(_tex(:pmatrix, "0 & -γ & β", "γ & 0 & -α", "-β & α & 0")),
```
hence the coordinate vector is .

For `n ≥ 4`` all further coefficients are used to fill up the following rows of the lower
triangular part – which determines the upper triangular part due to skew-symmetry
"""

@doc "$(_doc_get_vector_On)"
get_vector(𝔤::OrthogonalLieAlgebra, X, ::DefaultLieAlgebraOrthogonalBasis)

@doc "$(_doc_get_vector_On)"
get_vector!(𝔤::OrthogonalLieAlgebra, c, X::DefaultLieAlgebraOrthogonalBasis)

function get_vector_lie(
        ::CommonUnitarySubAlgebra{ManifoldsBase.ℝ, <:ManifoldsBase.TypeParameter{Tuple{2}}},
        c,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
        ::Type{<:SMatrix{2, 2, T}},
    ) where {T}
    @assert size(c) == (1,)
    return SMatrix{2, 2, T}(0, c[1], -c[1], 0)
end

function get_vector_lie!(
        ::CommonUnitarySubAlgebra{ManifoldsBase.ℝ, <:ManifoldsBase.TypeParameter{Tuple{2}}},
        X,
        c,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    )
    @assert size(X) == (2, 2)
    @assert size(c) == (1,)
    X[1, 1] = 0.0
    X[2, 1] = c[1]
    X[1, 2] = -c[1]
    X[2, 2] = 0.0
    return X
end
function get_vector_lie!(
        𝔤::CommonUnitarySubAlgebra{ℝ, <:ManifoldsBase.TypeParameter{Tuple{n}}},
        X,
        c,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    ) where {n}
    @assert size(X) == (n, n)
    @assert length(c) == manifold_dimension(𝔤)
    @assert n > 2
    _get_vector_lie_On!(X, c)
    return X
end
function get_vector_lie!(
        𝔤::CommonUnitarySubAlgebra{ManifoldsBase.ℝ}, X, c, ::DefaultLieAlgebraOrthogonalBasis{ℝ}
    )
    G = 𝔤.manifold
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    @assert length(c) == manifold_dimension(𝔤)
    @assert size(X) == (n, n)
    if n == 2
        X[1, 1] = 0.0
        X[2, 1] = c[1]
        X[1, 2] = -c[1]
        X[2, 2] = 0.0
    else
        _get_vector_lie_On!(X, c)
    end
    return X
end

function _get_vector_lie_On!(X, c)
    n = size(X, 1)
    @inbounds begin
        X[1, 1] = 0
        X[2, 2] = 0
        X[3, 3] = 0
        X[1, 2] = -c[3]
        X[2, 1] = c[3]
        X[1, 3] = c[2]
        X[3, 1] = -c[2]
        X[2, 3] = -c[1]
        X[3, 2] = c[1]
        k = 4
        for i in 4:n
            X[i, i] = 0.0
            for j in 1:(i - 1)
                X[i, j] = c[k]
                X[j, i] = -c[k]
                k += 1
            end
        end
    end
    return X
end

function ManifoldsBase.get_vector(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        p::SMatrix,
        Xⁱ,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    )
    return @SMatrix [0 -Xⁱ[]; Xⁱ[] 0]
end
function ManifoldsBase.get_vector(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}},
        p::SMatrix,
        Xⁱ,
        ::DefaultLieAlgebraOrthogonalBasis{ℝ},
    )
    return @SMatrix [0 -Xⁱ[3] Xⁱ[2]; Xⁱ[3] 0 -Xⁱ[1]; -Xⁱ[2] Xⁱ[1] 0]
end

_doc_log_O2_id = """
    log(G, g)
    log!(G, X, g)

Compute the Lie group logarithm function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(2)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(2)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

For the two-dimensional case, any rotation matrix ``g`` can be represented as ``$(_tex(:pmatrix, "$(_tex(:cos))(α) & -$(_tex(:sin))(α)", "$(_tex(:sin))(α) & $(_tex(:cos))(α)"))``.
For the [`SpecialOrthogonalGroup`](@ref), ``g`` might also include reflections.

The logarithm is then

```math
$(_tex(:log))_{$(_math(:G))}(g) =  $(_tex(:pmatrix, "0 & α", "-α & 0")).
```

This result can be computed in-place of `X`

Note the logarithmic map is only locally around the identity uniquely determined.
Especially, since ``$(_math(:O))(2)`` consists of two disjoint connected components and the exponential map is smooth,
for any ``g`` in the other component, the logarithmic map is defined, but not the inverse of the exponential map.
"""

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g)

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X, g)

_doc_log_O3_id = """
    log(G, g)
    log!(G, X, g)

Compute the Lie group logarithm function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(3)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(3)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

Here, ``$(_tex(:exp))_{$(_math(:G))}(X) = g`` is inverted using the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Olinde_Rodrigues)

```math
$(_tex(:exp))_{$(_math(:G))}(X) = I_3 + $(_tex(:frac, _tex(:sin, "α"), "α"))X + $(_tex(:frac, "(1 - $(_tex(:cos, "α")))", "α^2"))X^2,
```

For ``α ∉ $(_tex(:Set, "0, π"))`` we obtain ``X`` from the observation that
```math
$(_tex(:rm, "tr"))(g) = 1 + 2$(_tex(:cos))(α)
$(_tex(:qquad))$(_tex(:text, " and "))$(_tex(:qquad))
$(_tex(:frac, "1", "2"))(g-g^$(_tex(:transp))) = $(_tex(:sin))(α)X.
```

For ``α = 0`` we have ``g = I_3`` and ``X = 0``.

For ``α = π`` we have to solve ``X^2 = $(_tex(:frac, "1", "2"))(g-I_3)``,
where ``X`` is skew-symmetric and hence we have to solve for three unknowns.

```math
$(_tex(:log))_{$(_math(:G))}(g) = X.
```

This result can be computed in-place of `X`

Note the logarithmic map is only locally around the identity uniquely determined.
Especially, since ``$(_math(:O))(3)`` consists of two disjoint connected components and the exponential map is smooth,
for any ``g`` in the other component, the logarithmic map is defined, but not the inverse of the exponential map.
"""

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, g)

_doc_log_O4_id = """
    log(G, g)
    log!(G, X, g)

Compute the Lie group logarithm function on the [`OrthogonalGroup`](@ref) ``$(_math(:O))(4)`` or [`SpecialOrthogonalGroup`](@ref) ``$(_math(:SO))(4)``,
where `e` is the [`Identity`](@ref)`{MatrixMultiplicationGroupOperation}` and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`) for dispatch.

The implementation is based on a generalized variant of the Rodrigues' like formula.
For details, see [GallierXu:2002; Section 3](@cite).

This result can be computed in-place of `X`

Note the logarithmic map is only locally around the identity uniquely determined.
"""

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g)

function ManifoldsBase.log(
        M::CommonUnitarySubGroup{ManifoldsBase.ℝ}, p::SMatrix, q::SMatrix
    )
    return log(M.manifold, p, q)
end

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log!(::OrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X, g)

function ManifoldsBase.log!(
        ::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}}, X, g
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
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{2}}},
        X,
        ::Identity{MatrixMultiplicationGroupOperation},
    )
    return zero_vector!(LieAlgebra(G), X)
end

function ManifoldsBase.log!(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{3}}},
        X::AbstractMatrix,
        q::AbstractMatrix,
    )
    cosθ = (tr(q) - 1) / 2
    if cosθ ≈ -1
        eig = eigen(q)
        ival = findfirst(λ -> isapprox(λ, 1), eig.values)
        inds = SVector{3}(1:3)
        ax = eig.vectors[inds, ival]
        return get_vector!(LieAlgebra(G), X, π * ax, DefaultLieAlgebraOrthogonalBasis())
    end
    X .= q ./ usinc_from_cos(cosθ)
    # project onto 𝔰𝔬(3) for numerical stability
    return project!(LieAlgebra(G), X, X)
end
function ManifoldsBase.log!(
        G::CommonUnitarySubGroup{ManifoldsBase.ℝ, ManifoldsBase.TypeParameter{Tuple{4}}},
        X::AbstractMatrix,
        q::AbstractMatrix,
    )
    cosα, cosβ = cos_angles_4d_rotation_matrix(q)
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
    X .= (X .- X') ./ 2
    return X
end

function Base.show(io::IO, G::OrthogonalGroup)
    size = get_parameter(G.manifold.size)[1]
    return print(io, "OrthogonalGroup($(size))")
end
