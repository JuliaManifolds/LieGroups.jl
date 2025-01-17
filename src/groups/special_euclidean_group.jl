# for (g, t)
const LeftSpecialEuclideanOperation = LeftSemidirectProductGroupOperation{
    <:MatrixMultiplicationGroupOperation,<:AdditionGroupOperation,LeftGroupOperationAction
}

const LeftSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:LeftSpecialEuclideanOperation,
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Rotations{T},<:Manifolds.Euclidean{T,ℝ}}
    },
}

# for (t, g)
const RightSpecialEuclideanOperation = RightSemidirectProductGroupOperation{
    <:AdditionGroupOperation,<:MatrixMultiplicationGroupOperation,LeftGroupOperationAction
}
const RightSpecialEuclideanGroup{T} = LieGroup{
    ℝ,
    <:RightSpecialEuclideanOperation,
    <:Manifolds.ProductManifold{
        ℝ,Tuple{<:Manifolds.Euclidean{T,ℝ},<:Manifolds.Rotations{T}}
    },
}

"""
    SpecialEuclideanGroup{T}

The special Euclidean group ``$(_math(:SE))(n) = $(_math(:SO))(n) ⋉ $(_math(:T))(n)``. is the Lie group consisting of the
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

# Constructor
    SpecialEuclideanGroup(n; variant=:left, kwargs...)
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
const SpecialEuclideanGroup{T} = Union{
    <:LeftSpecialEuclideanGroup{T},<:RightSpecialEuclideanGroup{T}
}

const SpecialEuclideanOperation = Union{
    <:LeftSemidirectProductGroupOperation{
        <:MatrixMultiplicationGroupOperation,
        <:AdditionGroupOperation,
        LeftGroupOperationAction,
    },
    <:RightSemidirectProductGroupOperation{
        <:AdditionGroupOperation,
        <:MatrixMultiplicationGroupOperation,
        LeftGroupOperationAction,
    },
}

"""
    SpecialEuclideanMatrixPoint <: AbstractLieGroupPoint

represent a point on some [`LieGroup`](@ref) by an [affine matrix](https://en.wikipedia.org/wiki/Affine_group#Matrix_representation).

```math
$(_tex(:pmatrix, "M & v", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 1")) ∈ ℝ^{(n+1)×(n+1)},
$(_tex(:qquad)) M ∈ ℝ^{n×n}, v ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.
"""
struct SpecialEuclideanMatrixPoint{T} <: AbstractLieGroupPoint
    value::T
end

function ManifoldsBase.allocate_on(
    M::SpecialEuclideanGroup, ::Type{SpecialEuclideanMatrixPoint}
)
    return SpecialEuclideanMatrixPoint(Matrix(undef, representation_size(M)...))
end

"""
    SpecialEuclideanMatrixTVector <: AbstractLieGroupPoint

represent a tangent vector on some [`LieGroup`](@ref) by a matrix of the form

```math
$(_tex(:pmatrix, "M & v", "$(_tex(:vec, "0"))_n^{$(_tex(:transp))} & 0")) ∈ ℝ^{(n+1)×(n+1)},
$(_tex(:qquad)) M ∈ ℝ^{n×n}, v ∈ $(_math(:T))(n),
```
where ``$(_tex(:vec, "0"))_n ∈ ℝ^n`` denotes the vector containing zeros.

While this tangent vector itself is not an affine matrix itself, it can be used for the Lie algebra of the affine group
"""
struct SpecialEuclideanMatrixTVector{T} <: AbstractLieAlgebraTVector
    value::T
end

ManifoldsBase.@manifold_element_forwards SpecialEuclideanMatrixPoint value
ManifoldsBase.@manifold_vector_forwards SpecialEuclideanMatrixTVector value
ManifoldsBase.@default_manifold_fallbacks SpecialEuclideanGroup SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTVector value value

@default_lie_group_fallbacks SpecialEuclideanGroup SpecialEuclideanMatrixPoint SpecialEuclideanMatrixTVector value value

"""
    SpecialEuclideanProductPoint <: AbstractLieGroupPoint

represent a point on a Lie group (explicitly) as a point that consists of components
"""
struct SpecialEuclideanProductPoint{T} <: AbstractLieGroupPoint
    value::T
end

"""
    SpecialEuclideanProductTVector <: AbstractLieGroupPoint

represent a point on a Lie algebra (explicitly) as a point that consists of components
"""
struct SpecialEuclideanProductTVector{T} <: AbstractLieAlgebraTVector
    value::T
end

ManifoldsBase.@manifold_element_forwards SpecialEuclideanProductPoint value
ManifoldsBase.@manifold_vector_forwards SpecialEuclideanProductTVector value

# This union we can also use for the matrix case where we do not care

function SpecialEuclideanGroup(n; variant=:left, kwargs...)
    SOn = SpecialOrthogonalGroup(n; kwargs...)
    Tn = TranslationGroup(n; kwargs...)
    variant ∉ [:left, :right] && error(
        "SE(n) requires a  variant ∉ [:left, :right] but you provided `variant=:$variant`",
    )
    return variant === :left ? SOn ⋉ Tn : Tn ⋊ SOn
end

"""
    default_left_action(G::SpecialOrthogonalGroup, ::TranslationGroup)

Return the default left action for the special Euclidean group ``$(_math(:SO))(n) ⋊ $(_math(:T))(n)``,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
default_left_action(::SpecialOrthogonalGroup, ::TranslationGroup) =
    LeftGroupOperationAction()

"""
    default_right_action(::TranslationGroup, G::SpecialOrthogonalGroup)

Return the default right action for the special Euclidean group,
that is the [`GroupOperationAction`](@ref)`(`[`LeftGroupOperationAction`](@ref)`(G.op))`.
"""
function default_right_action(::TranslationGroup, ::SpecialOrthogonalGroup)
    return LeftGroupOperationAction()
end

function _check_matrix_affine(p, n; v=1, kwargs...)
    if !isapprox(p[end, :], [zeros(size(p, 2) - 1)..., v]; kwargs...)
        return nothing
        DomainError(p[end, :], "The last row of $p is not of form [0,..,0,$v].")
    end
    return nothing
end
# Order in a unified way
function ManifoldsBase.check_point(G::LeftSpecialEuclideanGroup, p; kwargs...)
    return _check_point(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], p; kwargs...)
end
function ManifoldsBase.check_point(G::RightSpecialEuclideanGroup, p; kwargs...)
    return _check_point(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], p; kwargs...)
end

function _check_point(
    G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, p; kwargs...
) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(p, n; v=1, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SOn
    errS = ManifoldsBase.check_point(
        Rotn, ManifoldsBase.submanifold_component(G, p, :Rotation); kwargs...
    )
    !isnothing(errS) && push!(errs, errS)
    # translate part
    errT = ManifoldsBase.check_point(
        Rn, ManifoldsBase.submanifold_component(G, p, :Translation); kwargs...
    )
    !isnothing(errT) && push!(errs, errT)
    if length(errs) > 1
        return ManifoldsBase.CompositeManifoldError(errs)
    end
    return length(errs) == 0 ? nothing : first(errs)
end

# Order in a unified way
function ManifoldsBase.check_vector(G::LeftSpecialEuclideanGroup, X; kwargs...)
    return _check_vector(G, G.manifold[1], G.manifold[2], G.op[1], G.op[2], X; kwargs...)
end
function ManifoldsBase.check_vector(G::RightSpecialEuclideanGroup, X; kwargs...)
    return _check_vector(G, G.manifold[2], G.manifold[1], G.op[2], G.op[1], X; kwargs...)
end

function _check_vector(
    G::SpecialEuclideanGroup{T}, Rotn, Rn, op1, op2, X; kwargs...
) where {T}
    errs = DomainError[]
    n = ManifoldsBase.get_parameter(Rotn.size)[1]
    errA = _check_matrix_affine(p, n; v=0, kwargs...)
    !isnothing(errA) && push!(errs, errA)
    # SO(n)  part
    SOn = LieGroup(Rotn, op1)
    errS = ManifoldsBase.check_vector(
        SOn, ManifoldsBase.submanifold_component(G, X, :Rotation); kwargs...
    )
    !isnothing(errS) && push!(errs, errS)
    # T(n) part
    Tn = LieGroup(Rn, G.op2)
    errT = ManifoldsBase.check_vector(
        Tn, ManifoldsBase.submanifold_component(G, X, :Translation); kwargs...
    )
    !isnothing(errT) && push!(errs, errT)
    (length(errs) > 1) && (return ManifoldsBase.CompositeManifoldError(errs))
    return length(errs) == 0 ? nothing : first(errs)
end

function ManifoldsBase.check_size(G::SpecialEuclideanGroup, p::AbstractMatrix)
    n = size(p)
    m = ManifoldsBase.representation_size(G)
    if n != m
        return DomainError(
            n,
            "The point $(p) can not belong to the Lie Group $(M), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
end

function ManifoldsBase.check_size(G::SpecialEuclideanGroup, p, X::AbstractMatrix)
    n = size(X)
    m = ManifoldsBase.representation_size(G)
    if n != m
        return DomainError(
            n,
            "The vector $(X) can not belong to the Lie Algebra of $(M), since its size $(n) is not equal to the manifolds representation size ($(m)).",
        )
    end
end

function _compose!(
    ::SpecialEuclideanGroup, k::AbstractMatrix, g::AbstractMatrix, h::AbstractMatrix
)
    copyto!(k, g * h)
    return k
end

_doc_exp_SE2_id = """
    exp(G::SpecialEuclidean, e, X)
    exp!(G::SpecialEuclideanG, e, g, X)

Compute the Lie group exponential function on the [`SpecialEuclidean`](@ref) `G```=$(_math(:SE))(2)``,
where `e` is the [`Identity`](@ref) on ``$(_math(:SE))(2)`` `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{2}}` for dispatch.

Since ``X = (Ω, v) ∈ $(_math(:se))(2)`` consists of a rotation component ``Ω ∈ $(_math(:se))(2)`` and a translation component ``v ∈ $(_math(:t))(2)``,
we can use [`vee`](@ref) on ``$(_math(:SO))(2)`` to obtain the angle of rotation ``α`` (or alternatively that ``$(_tex(:sqrt,2))α = $(_tex(:norm, "Ω"))``

For ``α ≠ 0`` define
```math
U_α = $(_tex(:frac, _tex(:sin)*"α", "α"))I_2 + $(_tex(:frac, "1-$(_tex(:cos))α", "α^2"))Ω,
```
and ``U_0 = I_2``, where ``I_2`` is the identity matrix.

Then the result ``g=(R,t)`` is given by
```math
R = $(_tex(:exp))_{$(_math(:SO))(2)}Ω ∈ $(_math(:SO))(2),
$(_tex(:quad))
$(_tex(:text, " and "))
$(_tex(:quad))
t = U_αv ∈ $(_math(:T))(2).
```

This result can be computed in-place of `g`.
"""

@doc "$(_doc_exp_SE2_id)"
ManifoldsBase.exp(
    ::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    ::Identity{SpecialEuclideanOperation},
    ::Any,
)

@doc "$(_doc_exp_SE2_id)"
function ManifoldsBase.exp!(
    G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    g,
    ::Identity{SpecialEuclideanOperation},
    X,
)
    Ω = submanifold_component(M, X, :Rotation)
    v = submanifold_component(M, X, :Translation)
    R = submanifold_component(M, g, :Rotation)
    t = submanifold_component(M, g, :Translation)
    α = norm(Ω) / sqrt(2) # skew symmetric, so the norm counts everything “twice” in the sqrt.
    SO2, T2 = _SOn_and_Tn(G)
    # (1) use the space of R to compute the U(α)
    if α ≈ 0
        copyto!(T2, t, v) # U(α) is the identity
    else
        copyto!(R, LinearAlgebra.I)
        R .*= (sin(α) / α)
        R .+= (1 - cos(α)) / α^2 .* Ω
        copyto!(T2, t, R * v)
    end
    # compute exp(Ω) in-place of R
    exp!(SO2, R, Identity(SO2), Ω)
    return g
end

_doc_exp_SE3_id = """
    exp(G::SpecialEuclidean, e, X)
    exp!(G::SpecialEuclideanG, e, g, X)

Compute the Lie group exponential function on the [`SpecialEuclidean`](@ref) `G```=$(_math(:SE))(3)``,
where `e` is the [`Identity`](@ref) on ``$(_math(:SE))(3)`` `G` uses a [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)`{Tuple{3}}` for dispatch.

Since ``X = (Ω, v) ∈ $(_math(:se))(3)`` consists of a rotation component ``Ω ∈ $(_math(:se))(3)`` and a translation component ``v ∈ $(_math(:t))(3)``,
we can use [`vee`](@ref) on ``$(_math(:SO))(3)`` computing the coefficients norm to obtain the angle of rotation ``α``
(or alternatively that ``$(_tex(:sqrt,2))α = $(_tex(:norm, "Ω"))``.

For ``α ≠ 0`` define
```math
U_α = I_3 + $(_tex(:frac, "1-$(_tex(:cos))α", "α^2"))Ω + $(_tex(:frac, "α-$(_tex(:sin))α", "α^3"))Ω^2,
```
and ``U_0 = I_3``, where ``I_2`` is the identity matrix.

Then the result ``g=(R,t)`` is given by
```math
R = $(_tex(:exp))_{$(_math(:SO))(3)}Ω ∈ $(_math(:SO))(3),
$(_tex(:quad))
$(_tex(:text, " and "))
$(_tex(:quad))
t = U_αv ∈ $(_math(:T))(3).
```

This result can be computed in-place of `g`.
"""

@doc "$(_doc_exp_SE3_id)"
ManifoldsBase.exp(
    ::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    ::Identity{SpecialEuclideanOperation},
    ::Any,
)

@doc "$(_doc_exp_SE3_id)"
function ManifoldsBase.exp!(
    G::SpecialEuclideanGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    g,
    ::Identity{SpecialEuclideanOperation},
    X,
)
    Ω = submanifold_component(M, X, :Rotation)
    v = submanifold_component(M, X, :Translation)
    R = submanifold_component(M, g, :Rotation)
    t = submanifold_component(M, g, :Translation)
    α = norm(Ω) / sqrt(2) # skew symmetric, so the norm counts everything “twice” in the sqrt.
    SO3, T3 = _SOn_and_Tn(G)
    # (1) use the space of R to compute the U(α)
    if α ≈ 0
        copyto!(T3, t, v) # U(α) is the identity
    else
        copyto!(R, LinearAlgebra.I)
        R .+= ((1 - cos(α)) / α^2) .* Ω
        R .+= ((α - sin(α)) / α^3) .* Ω^2
        copyto!(T3, t, R * v)
    end
    # compute exp(Ω) in-place of R
    exp!(SO3, R, Identity(SO3), Ω)
    return g
end

function ManifoldsBase.exp!(
    G::SpecialEuclideanGroup,
    g::AbstractMatrix,
    ::Identity{SpecialEuclideanOperation},
    X::AbstractMatrix,
)
    copyto!(g, exp(X))
    return g
end
# Do something similar for the component case?

function identity_element(G::SpecialEuclideanGroup)
    return identity_element(G, AbstractMatrix)
end
function identity_element(G::SpecialEuclideanGroup, ::Type{<:AbstractMatrix})
    q = zeros(ManifoldsBase.representation_size(G)...)
    return identity_element!(G, q)
end
function identity_element(G::SpecialEuclideanGroup, ::Type{SpecialEuclideanMatrixPoint})
    q = zeros(ManifoldsBase.representation_size(G)...)
    identity_element!(G, q)
    return SpecialEuclideanMatrixPoint(q)
end
function identity_element!(::SpecialEuclideanGroup, q::AbstractMatrix)
    copyto!(q, I)
    return q
end
function identity_element!(G::SpecialEuclideanGroup, q::SpecialEuclideanMatrixPoint)
    identity_element!(G, q.value)
    return q
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
function inv!(G::SpecialEuclideanGroup, h, g)
    rg = submanifold_component(G, g, :Rotation)
    tg = submanifold_component(G, g, :Translation)
    rh = submanifold_component(G, h, :Rotation)
    th = submanifold_component(G, h, :Translation)
    copyto!(rh, transpose(rg))
    copyto!(th, -rh * tg)
    return h
end
function inv!(G::SpecialEuclideanGroup, q, ::Identity{<:SpecialEuclideanOperation})
    return identity_element!(G, q)
end

function ManifoldsBase.isapprox(
    G::SpecialEuclideanGroup, g::AbstractMatrix, h::AbstractMatrix; kwargs...
)
    return isapprox(g, h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::SpecialEuclideanGroup,
    g::Identity{SpecialEuclideanOperation},
    h::AbstractMatrix;
    kwargs...,
)
    return isapprox(h, identity_element(G); kwargs...)
end
function ManifoldsBase.isapprox(
    G::SpecialEuclideanGroup,
    g::AbstractMatrix,
    h::Identity{SpecialEuclideanOperation};
    kwargs...,
)
    return isapprox(g, identity_element(G); kwargs...)
end

function is_identity(G::SpecialEuclideanGroup, g::AbstractMatrix; kwargs...)
    return isapprox(g, identity_element(G); kwargs...)
end

function ManifoldsBase.representation_size(G::SpecialEuclideanGroup)
    s = Manifolds.get_parameter(G.manifold[1].size)[1]
    return (s + 1, s + 1)
end

function Random.rand!(
    rng::AbstractRNG,
    G::SpecialEuclideanGroup,
    pX::AbstractMatrix;
    vector_at=nothing,
    kwargs...,
)
    SOn, Tn = _SOn_and_Tn(G)
    if vector_at === nothing # for points -> pass to manifold
        rand!(
            rng,
            SOn,
            submanifold_component(G, pX, :Rotation);
            vector_at=vector_at,
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, pX, :Translation);
            vector_at=vector_at,
            kwargs...,
        )
    else # for tangent vectors -> subset the vector_at as well.
        rand!(
            rng,
            SOn,
            submanifold_component(G, pX, :Rotation);
            vector_at=submanifold_component(G, vector_at, :Rotation),
            kwargs...,
        )
        rand!(
            rng,
            Tn,
            submanifold_component(G, pX, :Translation);
            vector_at=submanifold_component(G, vector_at, :Translation),
            kwargs...,
        )
    end
    return pX
end

function _SOn_and_Tn(G::LeftSpecialEuclideanGroup)
    return map(LieGroup, G.manifold.manifolds, G.op.operations)
end
function _SOn_and_Tn(G::RightSpecialEuclideanGroup)
    A = map(LieGroup, G.manifold.manifolds, G.op.operations)
    return A[2], A[1]
end

function ManifoldsBase.zero_vector(
    G::SpecialEuclidean, e::Identity{SpecialEuclideanOperation}
)
    n = Manifolds.get_parameter(G.manifold[1].size)[1]
    return zeros(n + 1, n + 1)
end

function Base.show(io::IO, G::SpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[1].size)[1]
    return print(io, "SpecialEuclideanGroup($(size))")
end
function Base.show(io::IO, G::RightSpecialEuclideanGroup)
    size = Manifolds.get_parameter(G.manifold[1].size)[1]
    return print(io, "SpecialEuclideanGroup($(size); variant=:right)")
end

#
#
# We overwrite the `submanifold_component` a bit different, since they are not ordered
# as (1) and (2) but more semantically as :Rotation and :Translation, we also use that here

@inline function ManifoldsBase.submanifold_component(G::SpecialEuclideanGroup, g, s::Symbol)
    return ManifoldsBase.submanifold_component(G, g, Val(s))
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{:Rotation}
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p, 1:n, 1:n)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup,
    p::Union{SpecialEuclideanMatrixPoint,SpecialEuclideanMatrixTVector},
    ::Val{:Translation},
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p.value, 1:n, n + 1)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup,
    p::Union{SpecialEuclideanMatrixPoint,SpecialEuclideanMatrixTVector},
    ::Val{:Rotation},
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p.value, 1:n, 1:n)
end
Base.@propagate_inbounds function ManifoldsBase.submanifold_component(
    G::SpecialEuclideanGroup, p, ::Val{:Translation}
)
    n = ManifoldsBase.get_parameter(G.manifold[1].size)[1]
    return view(p, 1:n, n + 1)
end
