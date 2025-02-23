"""
    ValidationLieGroup{L<:AbstractLieGroup} <: AbstractLieGroup

A Lie group to add tests to input parameters and ouptut values of functions defined
for [`LieGroups`](@ref).


Using the `ignore_contexts` keyword allows to specify a single `Symbol` or a vector of `Symbols`
Of which contexts to ignore.

Current contexts are
* `:All`: disable all checks
* `:Point`: checks for points
* `:Algebra`: checks related to the [`LieAlgebra`](@ref)
* `:Output`: checks for output
* `:Input`: checks for input variables

# Fields

* `lie_group::L` the [`AbstractLieGroup`](@ref) to be decorated
* ``

# Constructor
    ValidationLieGroup(L::AbstractLieGroup, check_manifold=true; kwargs...)

Generate the Validation Lie Group for the given [`AbstractLieGroup`](@ref) `L`.
If `check_manifold` is set to `true` the inner manifold is additionally wrapped
in a [`ValidationManifold`](@exref `ManifoldsBase.ValidationManifold`).
All suitable keywords are passed to the constructor of the validation manifold as well.

# Keyword arguments

* `error::Symbol=:error`: specify how errors in the validation should be reported.
  this is passed to [`is_point`](@ref) and [`is_vector`](@ref) as the `error` keyword argument.
  Available values are `:error`, `:warn`, `:info`, and `:none`. Every other value is treated as `:none`.
* `ignore_contexts = Vector{Symbol}()` a vector to indicate which validation contexts should not be performed.
* `ignore_functions=Dict{Function,Union{Symbol,Vector{Symbol}}}()` a dictionary to disable certain contexts within functions.
  The key here is the non-mutating function variant (if it exists). The contexts are the same as in `ignore_contexts`.
"""
struct ValidationLieGroup{
    𝔽,
    O<:AbstractGroupOperation,
    M<:ManifoldsBase.AbstractManifold{𝔽},
    L<:LieGroup{𝔽,O,M},
    D<:Dict{<:Function,<:Union{Symbol,<:AbstractVector{Symbol}}},
    V<:AbstractVector{Symbol},
} <: AbstractLieGroup{𝔽,O,M}
    lie_group::L
    mode::Symbol
    ignore_functions::D
    ignore_contexts::V
end
function ValidationLieGroup(
    L::LieGroup,
    check_manifold::Bool=true;
    error::Symbol=:error,
    ignore_functions::D=Dict{Function,Union{Symbol,<:Vector{Symbol}}}(),
    ignore_contexts::V=Vector{Symbol}(),
) where {
    D<:Dict{<:Function,<:Union{Symbol,<:AbstractVector{Symbol}}},V<:AbstractVector{Symbol}
}
    if check_manifold
        M = ManifoldBase.ValidationManifold(
            L.manifold;
            error=error,
            ignore_functions=ignore_functions,
            ignore_contexts=ignore_contexts,
        )
        _lie_group = LieGroup(L.op, M)
    else
        _lie_group = L
    end
    return ValidationLieGroup(L, error, ignore_functions, ignore_contexts)
end

struct ValidationLieAlgebraTangentVector{T} <: AbstractLieAlgebraTangentVector
    value::T
end
ManifoldsBase.internal_value(X::ValidationLieAlgebraTangentVector) = X.value
ManifoldsBase.@manifold_vector_forwards ValidationLieAlgebraTangentVector value
@default_lie_group_fallbacks ValidationLieGroup AbstractGroupOperation ManifoldsBase.ValidationMPoint ValidationLieAlgebraTangentVector value value

#
#
# An access helper function

"""
    _vLc(M::ValidationLieGroup, f::Function, context::Symbol)
    _vLc(M::ValidationLieGroup, f::Function, context::NTuple{N,Symbol}) where {N}

Return whether a check should be performed within `f` and the `context`(`s`) provided.

This function returns false and hence indicates not to check, when
* (one of the) `context`(`s`) is in the ignore list for `f` within `ignore_functions`
* (one of the) `context`(`s`) is in the `ignore_contexts` list

Otherwise the test is active.

!!! Note
   This function is internal and used very often, so it has a very short name;
    `_vLc` stands for "`ValidationLieGroup` check".
"""
function _vLc end

function _vLc(G::ValidationLieGroup, ::Nothing, context::Symbol)
    # Similarly for the global contexts
    (:All ∈ G.ignore_contexts) && return false
    (context ∈ G.ignore_contexts) && return false
    return true
end
function _vLc(G::ValidationLieGroup, f::Function, context::Symbol)
    if haskey(G.ignore_functions, f)
        # If :All is present -> deactivate
        !_vLc(G.ignore_functions[f], :All) && return false
        # If any of the provided contexts is present -> deactivate
        !_vLc(G.ignore_functions[f], context) && return false
    end
    !_vLc(G, nothing, context) && return false
    return true
end
function _vLc(G::ValidationLieGroup, f, contexts::NTuple{N,Symbol}) where {N}
    for c in contexts
        !_vLc(G, f, c) && return false
    end
    return true
end
# Sub tests: is any of a in b? Then return false – b is from before always a symbol already
# If a and b are symbols, equality is checked
_vLc(a::Symbol, b::Symbol) = !(a === b)
# If a is a vector multiple, then return false if b appears in a
_vLc(a::Union{<:NTuple{N,Symbol} where {N},<:AbstractVector{Symbol}}, b::Symbol) = !(b ∈ a)

"""
    _msg(G::ValidationLieGroup, str; error=:None, within::Union{Nothing,<:Function} = nothing,
    context::Union{NTuple{N,Symbol} where N} = NTuple{0,Symbol}())

issue a message `str` according to the mode `mode` (as `@error`, `@warn`, `@info`).
"""
function _msg(
    G::ValidationLieGroup,
    str;
    error=G.mode,
    within::Union{Nothing,<:Function}=nothing,
    context::Union{NTuple{N,Symbol} where N}=NTuple{0,Symbol}(),
)
    !_vLc(G, within, context) && return nothing
    (error === :error) && (throw(ErrorException(str)))
    (error === :warn) && (@warn str)
    (error === :info) && (@info str)
    return nothing
end
function ManifoldsBase._msg(
    G::ValidationLieGroup,
    err::Union{DomainError,ArgumentError,ErrorException};
    error=G.mode,
    within::Union{Nothing,<:Function}=nothing,
    context::Union{NTuple{N,Symbol} where N}=NTuple{0,Symbol}(),
)
    !_vLc(G, within, context) && return nothing
    (error === :error) && (throw(err))
    (error === :warn) && (@warn "$err")
    (error === :info) && (@info "$err")
    return nothing
end

#
#
# Implement all of the interface but include checks
Identity(VG::ValidationLieGroup) = Identity(VG.lie_group)

function Base.adjoint(G::ValidationLieGroup, g, X; kwargs...)
    is_point(G, g; widthin=adjoint, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=adjoint, contect=(:Input,), kwargs...)
    Y = adjoint(G.lie_group, internal_value(g), internal_value(X))
    is_point(LieAlgebra(G), Y; widthin=adjoint, contect=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function adjoint!(G::ValidationLieGroup, Y, g, X; kwargs...)
    is_point(G, g; widthin=adjoint, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=adjoint, contect=(:Input,), kwargs...)
    adjoint!(G.lie_group, internal_value(Y), internal_value(g), internal_value(X))
    is_point(LieAlgebra(G), Y; widthin=adjoint, contect=(:Output,), kwargs...)
    return Y
end

Manifolds.base_manifold(G::ValidationLieGroup) = base_manifold(G.lie_group)

function _compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; widthin=compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=compose, contect=(:Input,), kwargs...)
    k = compose(G.lie_group, internal_value(g), internal_value(h))
    is_point(G, k; widthin=compose, contect=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function _compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; widthin=compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=compose, contect=(:Input,), kwargs...)
    compose!(G.lie_group, internal_value(k), internal_value(g), internal_value(h))
    is_point(G, k; widthin=compose, contect=(:Output,), kwargs...)
    return k
end

function conjugate(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; widthin=conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=conjugate, contect=(:Input,), kwargs...)
    k = conjugate(G.lie_group, internal_value(g), internal_value(h))
    is_point(G, k; widthin=conjugate, contect=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function conjugate!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; widthin=conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=conjugate, contect=(:Input,), kwargs...)
    conjugate!(G.lie_group, internal_value(k), internal_value(g), internal_value(h))
    is_point(G, k; widthin=conjugate, contect=(:Output,), kwargs...)
    return k
end

function diff_conjugate(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    Y = diff_conjugate(G.lie_group, internal_value(g), internal_value(h), internal_value(X))
    is_point(LieAlgebra(G), Y; widthin=diff_conjugate, contect=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_conjugate!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_conjugate, contect=(:Input,), kwargs...)
    diff_conjugate!(
        G.lie_group,
        internal_value(Y),
        internal_value(g),
        internal_value(h),
        internal_value(X),
    )
    is_point(LieAlgebra(G), Y; widthin=diff_conjugate, contect=(:Output,), kwargs...)
    return Y
end

function diff_inv(G::ValidationLieGroup, g, X; kwargs...)
    is_point(G, g; widthin=diff_inv, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_inv, contect=(:Input,), kwargs...)
    Y = diff_inv(G.lie_group, internal_value(g), internal_value(X))
    is_point(LieAlgebra(G), Y; widthin=diff_inv, contect=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_inv!(G::ValidationLieGroup, Y, h, X; kwargs...)
    is_point(G, g; widthin=diff_inv, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_inv, contect=(:Input,), kwargs...)
    diff_inv!(G.lie_group, internal_value(Y), internal_value(g), internal_value(X))
    is_point(LieAlgebra(G), Y; widthin=diff_inv, contect=(:Output,), kwargs...)
    return Y
end

function diff_left_compose(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    Y = diff_left_compose(
        G.lie_group, internal_value(g), internal_value(h), internal_value(X)
    )
    is_point(LieAlgebra(G), Y; widthin=diff_left_compose, contect=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_left_compose!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_left_compose, contect=(:Input,), kwargs...)
    diff_left_compose!(
        G.lie_group,
        internal_value(Y),
        internal_value(g),
        internal_value(h),
        internal_value(X),
    )
    is_point(LieAlgebra(G), Y; widthin=diff_left_compose, contect=(:Output,), kwargs...)
    return Y
end

function diff_right_compose(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    Y = diff_right_compose(
        G.lie_group, internal_value(g), internal_value(h), internal_value(X)
    )
    is_point(LieAlgebra(G), Y; widthin=diff_right_compose, contect=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_right_compose!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=diff_right_compose, contect=(:Input,), kwargs...)
    diff_right_compose!(
        G.lie_group,
        internal_value(Y),
        internal_value(g),
        internal_value(h),
        internal_value(X),
    )
    is_point(LieAlgebra(G), Y; widthin=diff_right_compose, contect=(:Output,), kwargs...)
    return Y
end

function Base.exp(G::ValidationLieGroup, X; kwargs...)
    is_point(LieAlgebra(G), X; widthin=exp, contect=(:Input,), kwargs...)
    g = exp(G.lie_group, internal_value(X))
    is_point(G, g; widthin=exp, contect=(:Output,), kwargs...)
    return ValidationMPoint(g)
end
function ManifoldsBase.exp!(G::ValidationLieGroup, g, X; kwargs...)
    is_point(LieAlgebra(G), X; widthin=exp, contect=(:Input,), kwargs...)
    exp!(G.lie_group, internal_value(g), internal_value(X))
    is_point(G, g; widthin=exp, contect=(:Output,), kwargs...)
    return g
end

function Base.exp(G::ValidationLieGroup, g, X; kwargs...)
    is_point(G, g; widthin=exp, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=exp, contect=(:Input,), kwargs...)
    h = exp(G.lie_group, internal_value(g), internal_value(X))
    is_point(G, h; widthin=exp, contect=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function ManifoldsBase.exp!(G::ValidationLieGroup, h, g, X; kwargs...)
    is_point(G, g; widthin=exp, contect=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; widthin=exp, contect=(:Input,), kwargs...)
    exp!(G.lie_group, internal_value(h), internal_value(g), internal_value(X))
    is_point(G, h; widthin=exp, contect=(:Output,), kwargs...)
    return g
end

function identity_element(G::ValidationLieGroup; kwargs...)
    g = identity_element(G.lie_group)
    is_point(G, g; widthin=identity_element, contect=(:Output,), kwargs...)
    return ValidationMPoint(g)
end

function identity_element!(G::ValidationLieGroup, g; kwargs...)
    identity_element!(G.lie_group, internal_value(g))
    is_point(G, g; widthin=identity_element, contect=(:Output,), kwargs...)
    return g
end

function Base.inv(G::ValidationLieGroup, g; kwargs...)
    is_point(G, g; widthin=inv, contect=(:Input,), kwargs...)
    h = inv(G.lie_group, internal_value(g))
    is_point(G, h; widthin=inv, contect=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function inv!(G::ValidationLieGroup, h, g; kwargs...)
    is_point(G, g; widthin=inv, contect=(:Input,), kwargs...)
    inv!(G.lie_group, internal_value(h), internal_value(g))
    is_point(G, h; widthin=inv, contect=(:Output,), kwargs...)
    return h
end

function inv_left_compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; widthin=inv_left_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=inv_left_compose, contect=(:Input,), kwargs...)
    k = inv_left_compose(G.lie_group, internal_value(g), internal_value(h))
    is_point(G, k; widthin=inv_left_compose, contect=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function inv_left_compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; widthin=inv_left_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=inv_left_compose, contect=(:Input,), kwargs...)
    inv_left_compose!(G.lie_group, internal_value(k), internal_value(g), internal_value(h))
    is_point(G, k; widthin=inv_left_compose, contect=(:Output,), kwargs...)
    return k
end

function inv_right_compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; widthin=inv_right_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=inv_right_compose, contect=(:Input,), kwargs...)
    k = inv_right_compose(G.lie_group, internal_value(g), internal_value(h))
    is_point(G, k; widthin=inv_right_compose, contect=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function inv_right_compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; widthin=inv_right_compose, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=inv_right_compose, contect=(:Input,), kwargs...)
    inv_right_compose!(G.lie_group, internal_value(k), internal_value(g), internal_value(h))
    is_point(G, k; widthin=inv_right_compose, contect=(:Output,), kwargs...)
    return k
end

is_identity(G::ValidationLieGroup, g) = is_identity(G.lie_group, internal_value(g))

function ManifoldsBase.is_point(G::ValidationLieGroup, g; kwargs...)
    return is_point(G.lie_group, internal_value(g); kwargs...)
end
function ManifoldsBase.is_point(G::ValidationLieGroup, e::Identity; kwargs...)
    return is_point(G.lie_group, e; kwargs...)
end
function ManifoldsBase.is_point(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lieGroup
    return is_point(LieAlgebra(G), internal_value(X); kwargs...)
end

function ManifoldsBase.isapprox(G::ValidationLieGroup, g, h; kwargs...)
    return isapprox(G.lie_group, internal_value(g), internal_value(h), h; kwargs...)
end

function jacobian_conjugate(G::ValidationLieGroup, g, h, args...; kwargs...)
    is_point(G, g; widthin=jacobian_conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=jacobian_conjugate, contect=(:Input,), kwargs...)
    J = jacobian_conjugate(G.lie_group, internal_value(g), internal_value(h), args...)
    return J
end
function jacobian_conjugate!(G::ValidationLieGroup, J, g, h, args...; kwargs...)
    is_point(G, g; widthin=jacobian_conjugate, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=jacobian_conjugate, contect=(:Input,), kwargs...)
    jacobian_conjugate!(G.lie_group, J, internal_value(g), internal_value(h), args...)
    return J
end

function Base.log(G::ValidationLieGroup, g; kwargs...)
    is_point(G, g; widthin=log, contect=(:Input,), kwargs...)
    X = log(G.lie_group, internal_value(g))
    is_point(LieAlgebra(G), X; widthin=log, contect=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function ManifoldsBase.log!(G::ValidationLieGroup, X, g; kwargs...)
    is_point(G, g; widthin=log, contect=(:Input,), kwargs...)
    log!(G.lie_group, internal_value(X), internal_value(g))
    is_point(LieAlgebra(G), X; widthin=log, contect=(:Input,), kwargs...)
    return X
end

function Base.log(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; widthin=log, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=log, contect=(:Input,), kwargs...)
    X = log(G.lie_group, internal_value(g), internal_value(h))
    is_point(LieAlgebra(G), X; widthin=log, contect=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function ManifoldsBase.log!(G::ValidationLieGroup, X, g, h; kwargs...)
    is_point(G, g; widthin=log, contect=(:Input,), kwargs...)
    is_point(G, h; widthin=log, contect=(:Input,), kwargs...)
    log!(G.lie_group, internal_value(X), internal_value(g), internal_value(h))
    is_point(LieAlgebra(G), X; widthin=log, contect=(:Input,), kwargs...)
    return X
end

function Base.rand(G::ValidationLieGroup; vector_at=nothing, kwargs...)
    if vector_at !== nothing
        is_point(G, vector_at; within=rand, context=(:Input,), kwargs...)
    end
    gX = rand(G.manifold; vector_at=vector_at, kwargs...)
    if vector_at !== nothing
        is_point(LieAlgebra(G), gX; within=rand, context=(:Output,), kwargs...)
    else
        is_point(G, gX; within=rand, context=(:Output,), kwargs...)
    end
    return gX
end

ManifoldsBase.manifold_dimension(G::ValidationLieGroup) = manifold_dimension(G.lie_group)

ManifoldsBase.representation_size(G::ValidationLieGroup) = representation_size(G.lie_group)
