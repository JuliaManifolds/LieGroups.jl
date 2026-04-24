# for (g, t)
const LeftSpecialEuclideanGroupOperation = LeftSemidirectProductGroupOperation{
    MatrixMultiplicationGroupOperation, AdditionGroupOperation, LeftMultiplicationGroupAction, ActionActsOnRight,
}

const LeftSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:LeftSpecialEuclideanGroupOperation,
    <:ProductManifold{ℝ, Tuple{<:Rotations{T}, <:Euclidean{ℝ, T}}},
}

# for (t, g)
const RightSpecialEuclideanGroupOperation = RightSemidirectProductGroupOperation{
    AdditionGroupOperation, MatrixMultiplicationGroupOperation, LeftMultiplicationGroupAction, ActionActsOnRight,
}
const RightSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:RightSpecialEuclideanGroupOperation,
    <:ProductManifold{ℝ, Tuple{<:Euclidean{ℝ, T}, <:Rotations{T}}},
}

"""
    SpecialEuclideanGroup{T}

The special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)`` is the Lie group consisting of the
[`LeftSemidirectProductGroupOperation`](@ref) of the [`SpecialOrthogonalGroup`](@ref) and the
[`TranslationGroup`](@ref) together with the [`GroupOperationAction`](@ref)`{`[`LeftGroupOperationAction`](@ref)`}`.

To be precise, the group operation is defined on ``$(_math(:SO))(n) ⋉ $(_math(:T))(n)`` as follows:

```math
(r_1, t_1) ⋅ (r_2, t_2) = (r_1$(_math(:∘))r_2, t_1 + r_1$(_math(:⋅))t_2),
```

where ``r_1,r_2 ∈ $(_math(:SO))(n)`` and ``t_1,t_2 ∈ $(_math(:T))(n)``

Analogously you can write this on elements of ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)`` as

```math
(t_1, r_1) ⋅ (t_2, r_2) = (t_1 + r_1$(_math(:⋅))s_2, r_1$(_math(:∘))r_2)
```

Both these cases can be represented in a single matrix in [affine form](https://en.wikipedia.org/wiki/Affine_group#Matrix_representation)

```math
g = $(_tex(:pmatrix, "r & t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")),
$(_tex(:qquad)) r ∈ $(_math(:SO))(n), t ∈ $(_math(:T))(n),
```

where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

We refer also in general to elements on ``$(_math(:SE))(n)`` as ``g``
and their rotation and translation components as ``r`` and ``t``, respectively.

Note further that in the notation above and in matrix form the default is the [`ActionActsOnRight`](@ref) action.

# Constructor
    SpecialEuclideanGroup(n::Int; variant=:left, kwargs...)
    SpecialOrthogonalGroup(n; kwargs...) ⋉ TranslationGroup(n; kwargs...)
    TranslationGroup(n; kwargs...) ⋊ SpecialOrthogonalGroup(n; kwargs...)

Generate special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``, where the first
constructor is equivalent to the second.

All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.

The default representation for ``$(_math(:SE))(n)`` is the affine form. Alternatively
you can use the `ArrayPartition` from [`RecursiveArrayTools.jl`](https://docs.sciml.ai/RecursiveArrayTools/stable/) to work on ``(r,t)``.
or for ``$(_math(:T))(n) ⋊ $(_math(:SO))(n)`` using the `ArrayPartition`s ``(t,r)``;
which corresponds to setting `variant=:right` in the first constructor.
"""
const SpecialEuclideanGroup{T} = Union{<:LeftSpecialEuclideanGroup{T}, <:RightSpecialEuclideanGroup{T}}

const SpecialEuclideanGroupOperation = Union{<:LeftSpecialEuclideanGroupOperation, <:RightSpecialEuclideanGroupOperation}

"""
    SpecialEuclideanMatrixPoint <: AbstractLieGroupPoint

represent a point on some [`AbstractLieGroup`](@ref) by an [affine matrix](https://en.wikipedia.org/wiki/Affine_group#Matrix_representation).

```math
$(_tex(:pmatrix, "M & v", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")) ∈ ℝ^{(n+1)×(n+1)},
$(_tex(:qquad)) M ∈ ℝ^{n×n}, v ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.
"""
struct SpecialEuclideanMatrixPoint{T} <: AbstractLieGroupPoint
    value::T
end

"""
    SpecialEuclideanMatrixTangentVector <: AbstractLieAlgebraTangentVector

represent a tangent vector on some [`AbstractLieGroup`](@ref) by a matrix of the form

```math
$(_tex(:pmatrix, "M & v", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 0")) ∈ ℝ^{(n+1)×(n+1)},
$(_tex(:qquad)) M ∈ ℝ^{n×n}, v ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

While this tangent vector itself is not an affine matrix itself, it can be used for the Lie algebra of the affine group
"""
struct SpecialEuclideanMatrixTangentVector{T} <: AbstractLieAlgebraTangentVector
    value::T
end

ManifoldsBase.@manifold_element_forwards SpecialEuclideanMatrixPoint value
ManifoldsBase.@manifold_vector_forwards SpecialEuclideanMatrixTangentVector value
ManifoldsBase.@default_manifold_fallbacks LeftSpecialEuclideanGroup SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTangentVector value value
ManifoldsBase.@default_manifold_fallbacks RightSpecialEuclideanGroup SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTangentVector value value

@default_lie_group_fallbacks LeftSpecialEuclideanGroup LeftSpecialEuclideanGroupOperation SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTangentVector value value
@default_lie_group_fallbacks RightSpecialEuclideanGroup RightSpecialEuclideanGroupOperation SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTangentVector value value

@default_lie_algebra_fallbacks LeftSpecialEuclideanGroup ℝ LeftSpecialEuclideanGroupOperation SpecialEuclideanMatrixTangentVector value
@default_lie_algebra_fallbacks RightSpecialEuclideanGroup ℝ RightSpecialEuclideanGroupOperation SpecialEuclideanMatrixTangentVector value

ManifoldsBase.internal_value(semp::SpecialEuclideanMatrixPoint) = semp.value
ManifoldsBase.internal_value(semtv::SpecialEuclideanMatrixTangentVector) = semtv.value

function ManifoldsBase.tangent_vector_type(
        ::SpecialEuclideanGroup, ::Type{SpecialEuclideanMatrixPoint}
    )
    return SpecialEuclideanMatrixTangentVector
end
function ManifoldsBase.tangent_vector_type(
        ::SpecialEuclideanGroup, ::Type{SpecialEuclideanMatrixPoint{T}}
    ) where {T}
    return SpecialEuclideanMatrixTangentVector{T}
end

"""
    SpecialEuclideanProductPoint <: AbstractLieGroupPoint

Represent a point on a Lie group (explicitly) as a point that consists of components
"""
struct SpecialEuclideanProductPoint{T} <: AbstractLieGroupPoint
    value::T
end

"""
    SpecialEuclideanProductTangentVector <: AbstractLieAlgebraTangentVector

Represent a point on a Lie algebra (explicitly) as a tangent vector that consists of components.
"""
struct SpecialEuclideanProductTangentVector{T} <: AbstractLieAlgebraTangentVector
    value::T
end

ManifoldsBase.@manifold_element_forwards SpecialEuclideanProductPoint value
ManifoldsBase.@manifold_vector_forwards SpecialEuclideanProductTangentVector value
ManifoldsBase.@default_manifold_fallbacks LeftSpecialEuclideanGroup SpecialEuclideanProductPoint SpecialEuclideanProductTangentVector value value
ManifoldsBase.@default_manifold_fallbacks RightSpecialEuclideanGroup SpecialEuclideanProductPoint SpecialEuclideanProductTangentVector value value

@default_lie_group_fallbacks LeftSpecialEuclideanGroup LeftSpecialEuclideanGroupOperation SpecialEuclideanProductPoint SpecialEuclideanProductTangentVector value value
@default_lie_group_fallbacks RightSpecialEuclideanGroup RightSpecialEuclideanGroupOperation SpecialEuclideanProductPoint SpecialEuclideanProductTangentVector value value

@default_lie_algebra_fallbacks LeftSpecialEuclideanGroup ℝ LeftSpecialEuclideanGroupOperation SpecialEuclideanProductTangentVector value
@default_lie_algebra_fallbacks RightSpecialEuclideanGroup ℝ RightSpecialEuclideanGroupOperation SpecialEuclideanProductTangentVector value

ManifoldsBase.internal_value(sepp::SpecialEuclideanProductPoint) = sepp.value
ManifoldsBase.internal_value(septv::SpecialEuclideanProductTangentVector) = septv.value

function ManifoldsBase.tangent_vector_type(
        ::SpecialEuclideanGroup, ::Type{SpecialEuclideanProductPoint}
    )
    return SpecialEuclideanProductTangentVector
end
function ManifoldsBase.tangent_vector_type(
        ::SpecialEuclideanGroup, ::Type{SpecialEuclideanProductPoint{T}}
    ) where {T}
    return SpecialEuclideanProductTangentVector{T}
end

# This union we can also use for the matrix case where we do not care

function SpecialEuclideanGroup(n::Int; variant::Symbol = :left, kwargs...)
    SOn = SpecialOrthogonalGroup(n; kwargs...)
    Tn = TranslationGroup(n; kwargs...)
    variant ∉ SA[:left, :right] && error(
        "SE(n) requires a  variant ∉ [:left, :right] but you provided `variant=:$variant`",
    )
    return variant === :left ? SOn ⋉ Tn : Tn ⋊ SOn
end

function _check_matrix_affine(g, n; v = 1, kwargs...)
    if !isapprox(g[end, :], [zeros(size(g, 2) - 1)..., v]; kwargs...)
        return DomainError(g[end, :], "The last row of $g is not of form [0,..,0,$v].")
    end
    return nothing
end
# Order in a unified way
function ManifoldsBase.check_point(G::LeftSpecialEuclideanGroup, g; kwargs...)
    return _check_point(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], g; kwargs...)
end
function ManifoldsBase.check_point(G::RightSpecialEuclideanGroup, g; kwargs...)
    return _check_point(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], g; kwargs...)
end

function _check_point(
        G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, g; kwargs...
    ) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(g, n; v = 1, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SOn
    errS = ManifoldsBase.check_point(
        Rotn, ManifoldsBase.submanifold_component(G, g, :Rotation); kwargs...
    )
    !isnothing(errS) && push!(errs, errS)
    # translate part
    errT = ManifoldsBase.check_point(
        Rn, ManifoldsBase.submanifold_component(G, g, :Translation); kwargs...
    )
    !isnothing(errT) && push!(errs, errT)
    if length(errs) > 1
        return ManifoldsBase.CompositeManifoldError(errs)
    end
    return length(errs) == 0 ? nothing : first(errs)
end

# Order in a unified way with identities as well for resolving ambiguities
function ManifoldsBase.check_vector(G::LeftSpecialEuclideanGroup, g, X; kwargs...)
    return _check_vector(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], g, X; kwargs...)
end
function ManifoldsBase.check_vector(G::RightSpecialEuclideanGroup, g, X; kwargs...)
    return _check_vector(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], g, X; kwargs...)
end
function _check_vector(
        G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, g, X; kwargs...
    ) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(X, n; v = 0, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SO(n)  part
    SOn = LieGroup(Rotn, op1)
    errS = ManifoldsBase.check_vector(
        SOn,
        ManifoldsBase.submanifold_component(G, g, :Rotation),
        ManifoldsBase.submanifold_component(G, X, :Rotation);
        kwargs...,
    )
    !isnothing(errS) && push!(errs, errS)
    # T(n) part
    Tn = LieGroup(Rn, op2)
    errT = ManifoldsBase.check_vector(
        Tn,
        ManifoldsBase.submanifold_component(G, g, :Translation),
        ManifoldsBase.submanifold_component(G, X, :Translation);
        kwargs...,
    )
    !isnothing(errT) && push!(errs, errT)
    (length(errs) > 1) && (return ManifoldsBase.CompositeManifoldError(errs))
    return length(errs) == 0 ? nothing : first(errs)
end
function ManifoldsBase.check_size(
        G::LG, g::AbstractMatrix; kwargs...
    ) where {LG <: SpecialEuclideanGroup}
    n = size(g)
    m = ManifoldsBase.representation_size(G)
    if n != m
        return DomainError(
            n,
            "The point $(g) can not belong to the Lie Group $(G), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
end
function ManifoldsBase.check_size(
        G::LG, g::AbstractMatrix, X::AbstractMatrix; kwargs...
    ) where {LG <: SpecialEuclideanGroup}
    n = size(X)
    m = ManifoldsBase.representation_size(G)
    if n != m
        return DomainError(
            n,
            "The point $(X) can not belong to the Lie Algebra $(LieAlgebra(G)), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
end
function _compose!(
        ::SpecialEuclideanGroup, k::AbstractMatrix, g::AbstractMatrix, h::AbstractMatrix
    )
    copyto!(k, g * h)
    return k
end

Base.convert(::Type{<:AbstractMatrix}, g::SpecialEuclideanMatrixPoint) = g.value
function Base.convert(::Type{<:SpecialEuclideanMatrixPoint}, p::AbstractMatrix)
    return SpecialEuclideanMatrixPoint(p)
end
Base.convert(::Type{<:AbstractMatrix}, X::SpecialEuclideanMatrixTangentVector) = X.value
function Base.convert(::Type{<:SpecialEuclideanMatrixTangentVector}, X::AbstractMatrix)
    return SpecialEuclideanMatrixTangentVector(X)
end

"""
    default_left_action(::SpecialOrthogonalGroup, ::TranslationGroup)

Return the default [`AbstractGroupActionType`](@ref) for the special Euclidean group ``$(_math(:SO))(n) ⋉ $(_math(:T))(n)``,
which is the [`LeftMultiplicationGroupAction`](@ref)
"""
default_left_action(::SpecialOrthogonalGroup, ::TranslationGroup) =
    LeftMultiplicationGroupAction()

"""
    default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)

Return the default [`AbstractGroupActionType`](@ref) for the special Euclidean group ``$(_math(:T))(n) ⋊ $(_math(:SO))(n)``,
which is the [`LeftMultiplicationGroupAction`](@ref)
"""
function default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)
    return LeftMultiplicationGroupAction()
end

@doc raw"""
    diff_left_compose(G::SpecialEuclideanGroup, g, h, X)

Compute the differential of left composition by `h` on [`SpecialEuclideanGroup`](@ref)
for tangent vector `X` at `g`.

Let
```math
h=\begin{pmatrix} R & t\\[4pt] 0 & 1\end{pmatrix},\qquad
X=\begin{pmatrix} X_{\mathrm{R}} & X_{\mathrm{t}}\\[4pt] 0 & 0\end{pmatrix},
```

where ``R\in SO(n)``, ``t\in\mathbb{R}^n``, ``X_{\mathrm{R}}`` is the skew-symmetric rotation block and ``X_{\mathrm{t}}`` the translation column.
Then the differential is the adjoint action by ``h^{-1}``:
```math
\mathrm{D}λ_h(X)=h^{-1}Xh=
\begin{pmatrix}
R^\top X_{\mathrm{R}} R & R^\top\bigl(X_{\mathrm{R}}\,t + X_{\mathrm{t}}\bigr)\\[4pt]
0 & 0
\end{pmatrix}.
```

Component-wise:
```math
Y_{\mathrm{R}} = R^\top X_{\mathrm{R}} R,\qquad Y_{\mathrm{t}} = R^\top\bigl(X_{\mathrm{R}}\,t + X_{\mathrm{t}}\bigr).
```
"""
diff_left_compose(G::SpecialEuclideanGroup, g, h, X)

function diff_left_compose!(G::SpecialEuclideanGroup, Y, g, h, X)
    GA = LieAlgebra(G)
    init_constants!(GA, Y)
    XR = submanifold_component(GA, X, Val(:Rotation))
    Xt = submanifold_component(GA, X, Val(:Translation))
    YR = submanifold_component(GA, Y, Val(:Rotation))
    Yt = submanifold_component(GA, Y, Val(:Translation))
    R = submanifold_component(G, h, Val(:Rotation))
    t = submanifold_component(G, h, Val(:Translation))
    A = R' * XR
    mul!(YR, A, R)
    Yt .= A * t .+ R' * Xt
    return Y
end


@doc raw"""
    diff_right_compose(G::SpecialEuclideanGroup, g, h, X)

Compute the differential of right composition by `h` on [`SpecialEuclideanGroup`](@ref)
for tangent vector `X` at `g`, which is equal to `X`.
"""
diff_right_compose(G::SpecialEuclideanGroup, g, h, X)

function diff_right_compose!(G::SpecialEuclideanGroup, Y, g, h, X)
    copyto!(Y, X)
    return Y
end

_doc_exp_SE2_id = """
    exp(G::SpecialEuclidean, X)
    exp!(G::SpecialEuclidean, g, X)

Compute the Lie group exponential function on the [`SpecialEuclideanGroup`](@ref) `G```=$(_math(:SE))(2)``
using a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{2}}` for dispatch.

The Lie algebra vector ``X = (Y, v) ∈ $(_math(:se))(2)`` consists of a rotation component ``Y ∈ $(_math(:so))(2)``
and a translation component ``v ∈ $(_math(:t))(2)``, so we can use [`vee`](@ref) on ``$(_math(:SO))(2)``
to obtain the angle of rotation ``α`` (or alternatively using ``$(_tex(:sqrt, 2))α = $(_tex(:norm, "Y"))``)

For ``α ≠ 0`` define
```math
U_α = $(_tex(:frac, _tex(:sin) * "α", "α"))I_2 + $(_tex(:frac, "1-$(_tex(:cos))α", "α^2"))Y,
```
and ``U_0 = I_2``, where ``I_2`` is the identity matrix.

Then the result ``g=(R,t)`` is given by
```math
R = $(_tex(:exp))_{$(_math(:SO))(2)}Y ∈ $(_math(:SO))(2),
$(_tex(:quad))
$(_tex(:text, " and "))
$(_tex(:quad))
t = U_αv ∈ $(_math(:T))(2).
```

This result can be computed in-place of `g`.
"""

@doc "$(_doc_exp_SE2_id)"
ManifoldsBase.exp(::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, ::Any)

@doc "$(_doc_exp_SE2_id)"
function ManifoldsBase.exp!(
        G::SpecialEuclideanGroup{<:ManifoldsBase.TypeParameter{Tuple{2}}},
        g::AbstractMatrix,
        X::AbstractMatrix,
    )
    init_constants!(G, g)
    _exp_SE2!(G, g, X)
    return g
end

function _exp_SE2!(G::SpecialEuclideanGroup{<:ManifoldsBase.TypeParameter{Tuple{2}}}, g, X)
    Y = submanifold_component(G, X, :Rotation)
    v = submanifold_component(G, X, :Translation)
    R = submanifold_component(G, g, :Rotation)
    t = submanifold_component(G, g, :Translation)
    α = norm(Y) / sqrt(2) # skew symmetric, so the norm counts everything “twice” in the sqrt.
    SO2, T2 = _SOn_and_Tn(G)
    # (1) use the space of R to compute the U(α)
    if α ≈ 0
        copyto!(T2, t, v) # U(α) is the identity
    else
        copyto!(R, LinearAlgebra.I)
        R .*= (sin(α) / α)
        R .+= (1 - cos(α)) / α^2 .* Y
        copyto!(T2, t, R * v)
    end
    # compute exp(Y) in-place of R
    exp!(SO2, R, Y)
    return g
end
_doc_exp_SE3_id = """
    exp(G::SpecialEuclidean, X)
    exp!(G::SpecialEuclidean, g, X)

Compute the Lie group exponential function on the [`SpecialEuclideanGroup`](@ref) `G```=$(_math(:SE))(3)``
using a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{3}}` for dispatch.

Since ``X = (Y, v) ∈ $(_math(:se))(3)`` consists of a rotation component ``Y ∈ $(_math(:se))(3)`` and a translation component ``v ∈ $(_math(:t))(3)``,
we can use [`vee`](@ref) on ``$(_math(:SO))(3)`` computing the coefficients norm to obtain the angle of rotation ``α``
(or alternatively using ``$(_tex(:sqrt, 2))α = $(_tex(:norm, "Y"))``).

For ``α ≠ 0`` define
```math
U_α = I_3 + $(_tex(:frac, "1-$(_tex(:cos))α", "α^2"))Y + $(_tex(:frac, "α-$(_tex(:sin))α", "α^3"))Y^2,
```
and ``U_0 = I_3``, where ``I_2`` is the identity matrix.

Then the result ``g=(R,t)`` is given by
```math
R = $(_tex(:exp))_{$(_math(:SO))(3)}Y ∈ $(_math(:SO))(3),
$(_tex(:quad))
$(_tex(:text, " and "))
$(_tex(:quad))
t = U_αv ∈ $(_math(:T))(3).
```

This result can be computed in-place of `g`.
"""

@doc "$(_doc_exp_SE3_id)"
ManifoldsBase.exp(::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, ::Any)

@doc "$(_doc_exp_SE3_id)"
function ManifoldsBase.exp!(
        G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        g::AbstractMatrix,
        X::AbstractMatrix,
    )
    init_constants!(G, g)
    _exp_SE3!(G, g, X)
    return g
end

function _exp_SE3!(G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g, X)
    Y = submanifold_component(G, X, :Rotation)
    v = submanifold_component(G, X, :Translation)
    R = submanifold_component(G, g, :Rotation)
    t = submanifold_component(G, g, :Translation)
    α = norm(Y) / sqrt(2) # skew symmetric, so the norm counts everything “twice” in the sqrt.
    SO3, T3 = _SOn_and_Tn(G)
    # (1) use the space of R to compute the U(α)
    if α ≈ 0
        copyto!(T3, t, v) # U(α) is the identity
    else
        copyto!(R, LinearAlgebra.I)
        R .+= ((1 - cos(α)) / α^2) .* Y
        R .+= ((α - sin(α)) / α^3) .* Y^2
        copyto!(T3, t, R * v)
    end
    # compute exp(Y) in-place of R
    exp!(SO3, R, Y)
    return g
end

function ManifoldsBase.exp!(::SpecialEuclideanGroup, g::AbstractMatrix, X::AbstractMatrix)
    copyto!(g, exp(X))
    return g
end

_doc_getindex_SE = """
    g[G::SpecialEuclideanGroup,s]
    getindex(g, G::SpecialEuclideanGroup, s)
    X[𝔤,s]
    getindex(g, 𝔤, s)

Access sub-parts of a [`SpecialEuclideanGroup`](@ref) `G` or its Lie algebra `𝔤`.
where `s` can be an index `:Rotation`, `:Translation` to access the single parts.
Use `:` to access all submanifold components as a unified tuple.
"""

@doc "$(_doc_getindex_SE)"
function Base.getindex(
        g::Union{SpecialEuclideanMatrixPoint, SpecialEuclideanProductPoint, Identity, AbstractMatrix},
        G::SpecialEuclideanGroup,
        s::Union{Symbol, Int},
    )
    return submanifold_component(G, g, s)
end

function Base.getindex(
        g::Union{SpecialEuclideanMatrixPoint, SpecialEuclideanProductPoint, Identity, AbstractMatrix},
        G::SpecialEuclideanGroup,
        ::Colon,
    )
    return submanifold_components(G, g)
end

@doc "$(_doc_getindex_SE)"
function Base.getindex(
        g::Union{SpecialEuclideanMatrixTangentVector, SpecialEuclideanProductTangentVector, AbstractMatrix},
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        s::Union{Symbol, Int},
    )
    return submanifold_component(𝔤, g, s)
end

@doc "$(_doc_getindex_SE)"
function Base.getindex(
        g::Union{SpecialEuclideanMatrixTangentVector, SpecialEuclideanProductTangentVector, AbstractMatrix},
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        ::Colon,
    )
    return submanifold_components(𝔤, g)
end

function identity_element(G::SpecialEuclideanGroup)
    return identity_element(G, AbstractMatrix{Float64})
end
function identity_element(G::SpecialEuclideanGroup, ::Type{<:AbstractMatrix{T}}) where {T}
    q = zeros(T, ManifoldsBase.representation_size(G)...)
    return identity_element!(G, q)
end
function identity_element(G::SpecialEuclideanGroup, ::Type{<:SpecialEuclideanMatrixPoint})
    return SpecialEuclideanMatrixPoint(identity_element(G, AbstractMatrix{Float64}))
end
function identity_element(
        G::SpecialEuclideanGroup, ::Type{<:SpecialEuclideanMatrixPoint{<:AbstractMatrix{T}}}
    ) where {T}
    q = zeros(T, ManifoldsBase.representation_size(G)...)
    identity_element!(G, q)
    return SpecialEuclideanMatrixPoint(q)
end
function identity_element!(::SpecialEuclideanGroup, q::AbstractMatrix)
    copyto!(q, I)
    return q
end

_doc_init_constants = """
    init_constants!(G::SpecialEuclidean, g)
    init_Constants!(𝔰𝔢::LieAlgebra{ℝ, SpecialEuclideanGroupOperation, SpecialEuclidean}, X)

Initialize the constant elements of `g` or `X`.

The matrix representation of ``g∈$(_math(:SE))(n)`` has a last row,
that contains zeros, besides the diagonal element, which is ``g_{n+1,n+1} = 1``.

The matrix representation of ``X∈$(_math(:se))(n)`` has a last row that contains zeros.

this function sets these entries accordingly.

Per default for other representations, this function does not change entries for them.
"""

@doc "$(_doc_init_constants)"
function init_constants!(G::SpecialEuclideanGroup, g::AbstractMatrix)
    n = get_parameter(G.manifold[1].size)[1]
    g[(n + 1), 1:n] .= 0
    g[n + 1, n + 1] = 1
    return g
end
function init_constants!(G::SpecialEuclideanGroup, g::SpecialEuclideanMatrixPoint)
    init_constants!(G, g.value)
    return g
end

@doc "$(_doc_init_constants)"
function init_constants!(
        G::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        X::AbstractMatrix,
    )
    n = get_parameter(G.manifold.manifold[1].size)[1]
    X[(n + 1), :] .= 0
    return X
end
function init_constants!(
        G::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        X::SpecialEuclideanMatrixTangentVector,
    )
    init_constants!(G, X.value)
    return X
end

# default: Do nothing
init_constants!(::AbstractManifold, gX) = gX

#overwrite default inner since here the access is a bit tricky.
function ManifoldsBase.inner(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup}, X, Y
    )
    G = base_lie_group(𝔤)
    SOn, Tn = _SOn_and_Tn(G)
    i1 = inner(
        LieAlgebra(SOn),
        submanifold_component(𝔤, X, :Rotation),
        submanifold_component(𝔤, Y, :Rotation),
    )
    i2 = inner(
        LieAlgebra(Tn),
        submanifold_component(𝔤, X, :Translation),
        submanifold_component(𝔤, Y, :Translation),
    )
    return i1 + i2
end

_doc_inv_SEn = """
    inv(G::SpecialEuclideanGroup, g)
    inv(G::SpecialEuclideanGroup, h, g)

Compute the inverse element of a ``g ∈ $(_math(:SE))(n)`` from the [`SpecialEuclideanGroup`](@ref)`(n)`.

In affine form, ``g = $(_tex(:pmatrix, "r & t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1"))``,
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

The inverse reads

```math
g^{-1} = $(_tex(:pmatrix, "r^{$(_tex(:transp))} & -r^{$(_tex(:transp))}t", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")).
```

This computation can be done in-place of `h`.
"""

@doc "$(_doc_inv_SEn)"
Base.inv(G::SpecialEuclideanGroup, g)

@doc "$(_doc_inv_SEn)"
function _inv!(G::SpecialEuclideanGroup, h::AbstractMatrix, g::AbstractMatrix)
    init_constants!(G, h)
    _inv_SE!(G, h, g)
    return h
end
function _inv_SE!(G::SpecialEuclideanGroup, h, g)
    rg = submanifold_component(G, g, :Rotation)
    tg = submanifold_component(G, g, :Translation)
    rh = submanifold_component(G, h, :Rotation)
    th = submanifold_component(G, h, :Translation)
    copyto!(rh, transpose(rg))
    return copyto!(th, -rh * tg)
end

function ManifoldsBase.isapprox(
        G::SpecialEuclideanGroup, g::AbstractMatrix, h::AbstractMatrix; kwargs...
    )
    return isapprox(g, h; kwargs...)
end

function ManifoldsBase.is_flat(G::SpecialEuclideanGroup)
    size = get_parameter(G.manifold[1].size)[1]
    return size <= 2
end

function is_identity(G::SpecialEuclideanGroup, g::AbstractMatrix; kwargs...)
    return isapprox(g, identity_element(G); kwargs...)
end

_doc_log_SE2_id = """
    log(G::SpecialEuclidean, g)
    log!(G::SpecialEuclidean, X, g)

Compute the Lie group logarithm function on the [`SpecialEuclideanGroup`](@ref) `G```=$(_math(:SE))(2)``,
and `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{2}}` for dispatch.

Since ``g=(R,t) ∈ $(_math(:SE))(2)`` consists of a rotation component ``R ∈ $(_math(:SO))(2)`` and a translation component ``t ∈ $(_math(:T))(2)``,
we first compute ``Y = $(_tex(:log))_{$(_math(:SO))(2)}R``.
Then we can use [`vee`](@ref) on ``$(_math(:SO))(2)`` to obtain the angle of rotation ``α`` (or alternatively using ``$(_tex(:sqrt, 2))α = $(_tex(:norm, "Y"))``)

For ``α ≠ 0`` define
```math
V_α = $(_tex(:frac, "α", "2")) $(_tex(:pmatrix, "$(_tex(:frac, _tex(:sin) * "α", "1-$(_tex(:cos))α")) & 1", "-1 & $(_tex(:frac, _tex(:sin) * "α", "1-$(_tex(:cos))α"))"))
```
and ``V_0 = I_2``, where ``I_2`` is the identity matrix. Note that this is the inverse of ``U_α`` as given in the group exponential

Then the result ``X = (Y, v) ∈ $(_math(:se))(2)`` is given by ``Y ∈ $(_math(:so))(2)`` as computed above and
```math
v = V_αr ∈ $(_math(:T))(2),
```
where ``v`` is computed in-place without setting up ``V_α``

This result can be computed in-place of `g`.
"""

@doc "$(_doc_log_SE2_id)"
ManifoldsBase.log(::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, ::Any)

@doc "$(_doc_log_SE2_id)"
function ManifoldsBase.log!(
        G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
        X::AbstractMatrix,
        g::AbstractMatrix,
    )
    init_constants!(LieAlgebra(G), X)
    _log_SE2!(G, X, g)
    return X
end

function _log_SE2!(G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X, g)
    R = submanifold_component(G, g, :Rotation)
    t = submanifold_component(G, g, :Translation)
    Y = submanifold_component(G, X, :Rotation)
    v = submanifold_component(G, X, :Translation)
    @assert size(v) == (2,)

    SO2, T2 = _SOn_and_Tn(G)
    log!(SO2, Y, R)

    @inbounds θ = Y[2]
    β = θ / 2
    α = θ ≈ 0 ? 1 - β^2 / 3 : β * cot(β)

    @inbounds begin
        v[1] = α * t[1] + β * t[2]
        v[2] = α * t[2] - β * t[1]
    end
    return X
end

_doc_log_SE3_id = """
    log(G::SpecialEuclidean, g)
    log!(G::SpecialEuclidean, X, g)

Compute the Lie group logarithm function on the [`SpecialEuclideanGroup`](@ref) `G```=$(_math(:SE))(3)``,
where `e` is the [`Identity`](@ref) on ``$(_math(:SE))(3)`` `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{3}}` for dispatch.

Since ``g=(R,t) ∈ $(_math(:SE))(3)`` consists of a rotation component ``R ∈ $(_math(:SO))(3)`` and a translation component ``t ∈ $(_math(:T))(2)``,
we first compute ``Y = $(_tex(:log))_{$(_math(:SO))(3)}R``.
Then we can use [`vee`](@ref) on ``$(_math(:SO))(3)`` to obtain the angle of rotation ``α`` (or alternatively using ``$(_tex(:sqrt, 2))α = $(_tex(:norm, "Y"))``)

For ``α ≠ 0`` define
```math
V_α = I_3 - $(_tex(:frac, "1", "2"))Y + β Y^2, $(_tex(:quad))$(_tex(:text, " where ")) β = $(_tex(:frac, "1", "α^2")) - $(_tex(:frac, "1 + $(_tex(:cos))(α)", "2α$(_tex(:sin))(α)"))
```
and ``V_0 = I_3``, where ``I_3`` is the identity matrix. Note that this is the inverse of ``U_α`` as given in the group exponential

Then the result ``X = (Y, v) ∈ $(_math(:se))(3)`` is given by ``Y ∈ $(_math(:so))(3)`` as computed above and
``v = V_α t ∈ $(_math(:t))(3)``.

This result can be computed in-place of `X`.
"""

@doc "$(_doc_log_SE3_id)"
ManifoldsBase.log(::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, ::Any)

@doc "$(_doc_log_SE3_id)"
function ManifoldsBase.log!(
        G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
        X::AbstractMatrix,
        g::AbstractMatrix,
    )
    init_constants!(LieAlgebra(G), X)
    _log_SE3!(G, X, g)
    return X
end
function _log_SE3!(G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, g)
    R = submanifold_component(G, g, :Rotation)
    t = submanifold_component(G, g, :Translation)
    Y = submanifold_component(G, X, :Rotation)
    v = submanifold_component(G, X, :Translation)

    @assert size(Y) == (3, 3)
    @assert size(v) == (3,)

    trR = tr(R)
    cosθ = (trR - 1) / 2
    θ = acos(clamp(cosθ, -1, 1))
    θ² = θ^2
    if θ ≈ 0
        α = 1 / 2 + θ² / 12
        β = 1 / 12 + θ² / 720
    else
        sinθ = sin(θ)
        α = θ / sinθ / 2
        β = 1 / θ² - (1 + cosθ) / 2 / θ / sinθ
    end

    Y .= (R .- transpose(R)) .* α
    Jₗ⁻¹ = I - Y ./ 2 .+ β .* Y^2
    mul!(v, Jₗ⁻¹, t)

    return X
end

function LinearAlgebra.norm(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        X::AbstractMatrix,
    )
    G = base_lie_group(𝔤)
    SOn, Tn = _SOn_and_Tn(G)
    n1 = norm(LieAlgebra(SOn), submanifold_component(𝔤, X, :Rotation))
    n2 = norm(LieAlgebra(Tn), submanifold_component(𝔤, X, :Translation))
    return norm([n1, n2])
end

function Base.:*(
        X::SpecialEuclideanMatrixTangentVector, Y::SpecialEuclideanMatrixTangentVector
    )
    return SpecialEuclideanMatrixTangentVector(X.value * Y.value)
end

_doc_lie_bracket_SEn = """
    lie_bracket(𝔰𝔢::LieAlgebra{ℝ, SpecialEuclideanGroupOperation, SpecialEuclideanGroup}, X, Y)
    lie_bracket!(𝔰𝔢::LieAlgebra{ℝ, SpecialEuclideanGroupOperation, SpecialEuclideanGroup}, Z, X, Y)

Calculate the Lie bracket between elements `X` and `Y` of the Lie algebra of the [`SpecialEuclideanGroup`](@ref).
For the matrix representation, cf. [`SpecialEuclideanMatrixTangentVector`](@ref) or a `<:AbstractMatrix`, the formula reads

```math
[X, Y] = XY-YX
```
"""

@doc "$(_doc_lie_bracket_SEn)"
lie_bracket(
    𝔰𝔢::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
    X::Union{<:AbstractMatrix, SpecialEuclideanMatrixTangentVector},
    Y::Union{<:AbstractMatrix, SpecialEuclideanMatrixTangentVector},
)

@doc "$(_doc_lie_bracket_SEn)"
function lie_bracket!(
        ::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        Z::Union{<:AbstractMatrix, SpecialEuclideanMatrixTangentVector},
        X::Union{<:AbstractMatrix, SpecialEuclideanMatrixTangentVector},
        Y::Union{<:AbstractMatrix, SpecialEuclideanMatrixTangentVector},
    )
    Z .= X * Y - Y * X
    return Z
end

function ManifoldsBase.log!(::SpecialEuclideanGroup, X::AbstractMatrix, g::AbstractMatrix)
    copyto!(X, log(g))
    return X
end

function ManifoldsBase.representation_size(G::SpecialEuclideanGroup)
    s = get_parameter(G.manifold[1].size)[1]
    return (s + 1, s + 1)
end

function Random.rand(
        rng::AbstractRNG, G::SEG, T::Type = Matrix{Float64}; kwargs...
    ) where {SEG <: SpecialEuclideanGroup}
    g = identity_element(G, T)
    return rand!(rng, G, g; kwargs...)
end
function Random.rand(
        G::SEG, T::Type = Matrix{Float64}; kwargs...
    ) where {SEG <: SpecialEuclideanGroup}
    g = identity_element(G, T)
    return rand!(G, g; kwargs...)
end
# this is always with vector_at=nothing
function Random.rand!(
        rng::AbstractRNG,
        G::SEG,
        g::Union{SpecialEuclideanMatrixPoint, SpecialEuclideanProductPoint};
        vector_at = nothing,
        kwargs..., #but ignore vector_at, since this is a point
    ) where {SEG <: SpecialEuclideanGroup}
    Random.rand!(rng, G, g.value; vector_at = nothing, kwargs...)
    return g
end

function Random.rand!(
        rng::AbstractRNG, G::SEG, gX::AbstractMatrix; vector_at = nothing, kwargs...
    ) where {SEG <: SpecialEuclideanGroup}
    SOn, Tn = _SOn_and_Tn(G)
    if vector_at === nothing # for points -> pass to manifold
        init_constants!(G, gX)
        rand!(
            rng,
            SOn,
            submanifold_component(G, gX, :Rotation);
            vector_at = vector_at,
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, gX, :Translation);
            vector_at = vector_at,
            kwargs...,
        )
    else # for tangent vectors -> subset the vector_at as well.
        init_constants!(LieAlgebra(G), gX)
        rand!(
            rng,
            SOn,
            submanifold_component(G, gX, :Rotation);
            vector_at = submanifold_component(G, vector_at, :Rotation),
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, gX, :Translation);
            vector_at = submanifold_component(G, vector_at, :Translation),
            kwargs...,
        )
    end
    return gX
end

function _SOn_and_Tn(G::LeftSpecialEuclideanGroup)
    return map(LieGroup, G.manifold.manifolds, G.op.operations)
end
function _SOn_and_Tn(G::RightSpecialEuclideanGroup)
    A = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return A[2], A[1]
end

function Base.show(io::IO, G::SpecialEuclideanGroup)
    size = get_parameter(G.manifold[1].size)[1]
    return print(io, "SpecialEuclideanGroup($(size))")
end
function Base.show(io::IO, G::RightSpecialEuclideanGroup)
    size = get_parameter(G.manifold[1].size)[1]
    return print(io, "SpecialEuclideanGroup($(size); variant=:right)")
end

#
#
# We overwrite the `submanifold_component` a bit different, since they are not ordered
# as (1) and (2) but more semantically as :Rotation and :Translation, we also use that here

@inline function ManifoldsBase.submanifold_component(G::SpecialEuclideanGroup, g, s::Symbol)
    return ManifoldsBase.submanifold_component(G, g, Val(s))
end
@inline function ManifoldsBase.submanifold_component(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup}, X, s::Symbol
    )
    return ManifoldsBase.submanifold_component(𝔤, X, Val(s))
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::SpecialEuclideanGroup,
        p::Union{AbstractMatrix, SpecialEuclideanMatrixPoint},
        ::Val{:Rotation},
    )
    n = ManifoldsBase.get_parameter(base_manifold(G)[1].size)[1]
    # view to be able to write, internal_value to “unpack” SE Matrices
    return view(ManifoldsBase.internal_value(p), 1:n, 1:n)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        X::Union{AbstractMatrix, SpecialEuclideanMatrixTangentVector},
        ::Val{:Rotation},
    )
    n = ManifoldsBase.get_parameter(base_manifold(𝔤)[1].size)[1]
    return view(ManifoldsBase.internal_value(X), 1:n, 1:n)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::SpecialEuclideanGroup,
        e::Identity{<:SpecialEuclideanGroupOperation},
        ::Val{:Rotation},
    )
    SOn, Tn = _SOn_and_Tn(G)
    return Identity(SOn)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::SpecialEuclideanGroup, p, ::Val{:Translation}
    )
    n = ManifoldsBase.get_parameter(base_manifold(G)[1].size)[1]
    return view(ManifoldsBase.internal_value(p), 1:n, n + 1)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        X,
        ::Val{:Translation},
    )
    n = ManifoldsBase.get_parameter(base_manifold(𝔤)[1].size)[1]
    return view(ManifoldsBase.internal_value(X), 1:n, n + 1)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
        G::SpecialEuclideanGroup,
        e::Identity{<:SpecialEuclideanGroupOperation},
        ::Val{:Translation},
    )
    SOn, Tn = _SOn_and_Tn(G)
    return Identity(Tn)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_components(
        G::LeftSpecialEuclideanGroup, p
    )
    return (
        submanifold_component(G, p, Val(:Rotation)),
        submanifold_component(G, p, Val(:Translation)),
    )
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_components(
        G::RightSpecialEuclideanGroup, p
    )
    return (
        submanifold_component(G, p, Val(:Translation)),
        submanifold_component(G, p, Val(:Rotation)),
    )
end

Base.@propagate_inbounds function ManifoldsBase.submanifold_components(
        𝔤::LieAlgebra{ℝ, <:LeftSpecialEuclideanGroupOperation, <:LeftSpecialEuclideanGroup}, X
    )
    return (
        submanifold_component(𝔤, X, Val(:Rotation)),
        submanifold_component(𝔤, X, Val(:Translation)),
    )
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_components(
        𝔤::LieGroups.LieAlgebra{
            ℝ, <:RightSpecialEuclideanGroupOperation, <:RightSpecialEuclideanGroup,
        },
        X,
    )
    return (
        submanifold_component(𝔤, X, Val(:Translation)),
        submanifold_component(𝔤, X, Val(:Rotation)),
    )
end

function ManifoldsBase.zero_vector(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        ::Type{<:AbstractMatrix{T}} = AbstractMatrix{Float64},
    ) where {T}
    G = 𝔤.manifold
    n = get_parameter(G.manifold[1].size)[1]
    return zeros(T, n + 1, n + 1)
end

function ManifoldsBase.zero_vector(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup},
        ::Type{SpecialEuclideanMatrixTangentVector{Matrix{T}}},
    ) where {T}
    G = 𝔤.manifold
    n = get_parameter(G.manifold[1].size)[1]
    return SpecialEuclideanMatrixTangentVector(zeros(T, n + 1, n + 1))
end
function ManifoldsBase.zero_vector!(
        𝔤::LieAlgebra{ℝ, <:SpecialEuclideanGroupOperation, <:SpecialEuclideanGroup}, X
    )
    G = 𝔤.manifold
    init_constants!(𝔤, X)
    Y = submanifold_component(G, X, :Rotation)
    v = submanifold_component(G, X, :Translation)
    SO3, T3 = _SOn_and_Tn(G)
    zero_vector!(LieAlgebra(SO3), Y)
    zero_vector!(LieAlgebra(T3), v)
    return X
end

#
# Group Action SE(n) on R^n
# on Matrices this is done in homogeneous coordinates - or by cutting the matrix accordingly
#
_doc_apple_SE_Rn = """
    apply(::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, g, p)
    apply!(::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, q, g, p)

Given the Lie group [`SpecialEuclideanGroup`](@ref)  and the [`Euclidean`](@extref `Manifolds.Euclidean`) manifold ``ℝ^n``,
this action performs both the rotation and translation on a vector ``p ∈ ℝ^n``,
that is, for ``g = (R, t) ∈ $(_math(:SE))(n)``, the action is given by

```math
q = σ_g(p) = g ⋅ p = Rp + t,
```

where the name of the action, [`LeftMultiplicationGroupAction`](@ref), indicates that the group element `g` acts on the left of the vector `p`,
and directly yields a multiplication if interpreted in homogeneous coordinates.
"""

@doc "$(_doc_apple_SE_Rn)"
apply(::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, g, p)

@doc "$(_doc_apple_SE_Rn)"
apply!(::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, q, g, p)

function apply!(
        A::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, q, g, p
    )
    Rn = A.manifold
    n = get_parameter(Rn.size)[1]
    mul!(q, g[1:n, 1:n], p)
    q .+= g[1:n, (n + 1)]
    return q
end

_doc_diff_apply_SE_Rn = """
    diff_apply(::GroupAction{LeftMultiplicationGroupAction,SpecialEuclideanGroup,Euclidean}, g, p, X)
    diff_apply!(::GroupAction{LeftMultiplicationGroupAction,SpecialEuclideanGroup,Euclidean}, Y, g, p, X)

Given the Lie group [`SpecialEuclideanGroup`](@ref)  and the [`Euclidean`](@extref `Manifolds.Euclidean`) manifold ``ℝ^n``,
the differential of the group action
this action performs both the rotation and translation on a vector ``p ∈ ℝ^n``,
that is, for ``g = (R, t) ∈ $(_math(:SE))(n)``, the differential is given by


```math
Y = Dσ_g(p)[X] = RX,
```

where the name of the action, [`LeftMultiplicationGroupAction`](@ref), indicates that the group element `g` acts on the left of the vector `p`,
and directly yields a multiplication if interpreted in homogeneous coordinates.
"""

@doc "$(LieGroups._doc_diff_apply_SE_Rn)"
LieGroups.diff_apply(
    ::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, g, p, X
)

@doc "$(_doc_diff_apply_SE_Rn)"
diff_apply!(
    ::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean}, Y, g, p, X
)

function diff_apply!(
        A::GroupAction{LeftMultiplicationGroupAction, <:SpecialEuclideanGroup, <:Euclidean},
        Y, g, p, X
    )
    Rn = A.manifold
    n = get_parameter(Rn.size)[1]
    mul!(Y, g[1:n, 1:n], X)
    return Y
end
