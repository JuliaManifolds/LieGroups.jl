"""
    ValidationLieGroup{L<:AbstractLieGroup} <: AbstractLieGroup

A Lie group to add tests to input parameters and output values of functions defined
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
* `mode::Symbol`: The mode to be used for error handling, either `:error` or `:warn`
* `ignore_contexts::AbstractVector{Symbol}`: store contexts to be ignored of validation.
* `ignore_functions::Dict{<:Function,<:Union{Symbol,<:AbstractVector{Symbol}}`:
  store contexts to be ignored with in a function or its mutating variant.

where all but the first field are analogous to the setups of the [`ValidationManifold`](@extref `ManifoldsBase.ValidationManifold`).
We refer to those docs for more examples on their meaning.

# Constructor
    ValidationLieGroup(L::AbstractLieGroup, check_manifold=true; kwargs...)

Generate the Validation Lie Group for the given [`AbstractLieGroup`](@ref) `L`.
If `check_manifold` is set to `true` the inner manifold is additionally wrapped
in a [`ValidationManifold`](@extref `ManifoldsBase.ValidationManifold`).
All suitable keywords are passed to the constructor of the validation manifold as well.

# Keyword arguments

* `error::Symbol=:error`: specify how errors in the validation should be reported.
  this is passed to [`is_point`](@extref `ManifoldsBase.is_point`) and [`is_vector`](@extref `ManifoldsBase.is_vector`) as the `error` keyword argument.
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
        M = ValidationManifold(
            L.manifold;
            error=error,
            ignore_functions=ignore_functions,
            ignore_contexts=ignore_contexts,
        )
        _lie_group = LieGroup(M, L.op)
    else
        _lie_group = L
    end
    return ValidationLieGroup(L, error, ignore_functions, ignore_contexts)
end

struct ValidationLieAlgebraTangentVector{T} <: AbstractLieAlgebraTangentVector
    value::T
end
ValidationLieAlgebraTangentVector(value::ValidationLieAlgebraTangentVector) = value
ManifoldsBase.internal_value(X::ValidationLieAlgebraTangentVector) = X.value
ManifoldsBase.@manifold_vector_forwards ValidationLieAlgebraTangentVector value
@default_lie_group_fallbacks ValidationLieGroup AbstractGroupOperation ManifoldsBase.ValidationMPoint ValidationLieAlgebraTangentVector value value

Identity(VG::ValidationLieGroup) = Identity(VG.lie_group)

# similar to `internal_value` but just for Validation  types
unwrap_validation(v) = v
unwrap_validation(vTV::ValidationLieAlgebraTangentVector) = vTV.value
unwrap_validation(vP::ValidationMPoint) = vP.value

unwrap_validation_type(v) = v
unwrap_validation_type(::Type{<:ValidationLieAlgebraTangentVector{T}}) where {T} = T
unwrap_validation_type(::Type{<:ValidationMPoint{<:P}}) where {P} = P

#
#
# An access helper function

"""
    _vLc(M::ValidationLieGroup, f::Function, context::Symbol)
    _vLc(M::ValidationLieGroup, f::Function, context::NTuple{N,Symbol}) where {N}
    _vLc(M::ValidationLieGroup, ::Nothing, context::NTuple{N,Symbol}) where {N}

Return whether a check should be performed within `f` and the `context`(`s`) provided,
if the second argument is `:Nothing`, only the context is checked

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
function _msg(
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

function Base.adjoint(G::ValidationLieGroup, g, X; kwargs...)
    is_point(G, g; within=adjoint, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=adjoint, context=(:Input,), kwargs...)
    Y = adjoint(G.lie_group, unwrap_validation(g), unwrap_validation(X))
    is_point(LieAlgebra(G), Y; within=adjoint, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function adjoint!(G::ValidationLieGroup, Y, g, X; kwargs...)
    is_point(G, g; within=adjoint, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=adjoint, context=(:Input,), kwargs...)
    adjoint!(G.lie_group, unwrap_validation(Y), unwrap_validation(g), unwrap_validation(X))
    is_point(LieAlgebra(G), Y; within=adjoint, context=(:Output,), kwargs...)
    return Y
end

ManifoldsBase.base_manifold(G::ValidationLieGroup) = base_manifold(G.lie_group)

function _compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; within=compose, context=(:Input,), kwargs...)
    is_point(G, h; within=compose, context=(:Input,), kwargs...)
    k = compose(G.lie_group, unwrap_validation(g), unwrap_validation(h))
    is_point(G, k; within=compose, context=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function _compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; within=compose, context=(:Input,), kwargs...)
    is_point(G, h; within=compose, context=(:Input,), kwargs...)
    compose!(G.lie_group, unwrap_validation(k), unwrap_validation(g), unwrap_validation(h))
    is_point(G, k; within=compose, context=(:Output,), kwargs...)
    return k
end

function conjugate(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; within=conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=conjugate, context=(:Input,), kwargs...)
    k = conjugate(G.lie_group, unwrap_validation(g), unwrap_validation(h))
    is_point(G, k; within=conjugate, context=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function conjugate!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; within=conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=conjugate, context=(:Input,), kwargs...)
    conjugate!(
        G.lie_group, unwrap_validation(k), unwrap_validation(g), unwrap_validation(h)
    )
    is_point(G, k; within=conjugate, context=(:Output,), kwargs...)
    return k
end

function Base.copyto!(G::ValidationLieGroup, h, g)
    copyto!(G.lie_group, unwrap_validation(h), unwrap_validation(g))
    return h
end
function Base.copyto!(
    G::ValidationLieGroup{𝔽,O}, e::Identity{O}, g
) where {𝔽,O<:AbstractGroupOperation}
    copyto!(G.lie_group, e, unwrap_validation(g))
    return e
end
function Base.copyto!(
    G::ValidationLieGroup{𝔽,O}, h, e::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    copyto!(G.lie_group, unwrap_validation(h), e)
    return h
end
function Base.copyto!(
    ::ValidationLieGroup{𝔽,O}, h::Identity{O}, ::Identity{O}
) where {𝔽,O<:AbstractGroupOperation}
    return h
end

function diff_conjugate(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; within=diff_conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_conjugate, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_conjugate, context=(:Input,), kwargs...)
    Y = diff_conjugate(
        G.lie_group, unwrap_validation(g), unwrap_validation(h), unwrap_validation(X)
    )
    is_point(LieAlgebra(G), Y; within=diff_conjugate, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_conjugate!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; within=diff_conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_conjugate, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_conjugate, context=(:Input,), kwargs...)
    diff_conjugate!(
        G.lie_group,
        unwrap_validation(Y),
        unwrap_validation(g),
        unwrap_validation(h),
        unwrap_validation(X),
    )
    is_point(LieAlgebra(G), Y; within=diff_conjugate, context=(:Output,), kwargs...)
    return Y
end

function diff_inv!(G::ValidationLieGroup, Y, g, X; kwargs...)
    is_point(G, g; within=diff_inv, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_inv, context=(:Input,), kwargs...)
    diff_inv!(G.lie_group, unwrap_validation(Y), unwrap_validation(g), unwrap_validation(X))
    is_point(LieAlgebra(G), Y; within=diff_inv, context=(:Output,), kwargs...)
    return Y
end

function diff_left_compose(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; within=diff_left_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_left_compose, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_left_compose, context=(:Input,), kwargs...)
    Y = diff_left_compose(
        G.lie_group, unwrap_validation(g), unwrap_validation(h), unwrap_validation(X)
    )
    is_point(LieAlgebra(G), Y; within=diff_left_compose, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_left_compose!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; within=diff_left_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_left_compose, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_left_compose, context=(:Input,), kwargs...)
    diff_left_compose!(
        G.lie_group,
        unwrap_validation(Y),
        unwrap_validation(g),
        unwrap_validation(h),
        unwrap_validation(X),
    )
    is_point(LieAlgebra(G), Y; within=diff_left_compose, context=(:Output,), kwargs...)
    return Y
end

function diff_right_compose(G::ValidationLieGroup, g, h, X; kwargs...)
    is_point(G, g; within=diff_right_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_right_compose, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_right_compose, context=(:Input,), kwargs...)
    Y = diff_right_compose(
        G.lie_group, unwrap_validation(g), unwrap_validation(h), unwrap_validation(X)
    )
    is_point(LieAlgebra(G), Y; within=diff_right_compose, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Y)
end
function diff_right_compose!(G::ValidationLieGroup, Y, g, h, X; kwargs...)
    is_point(G, g; within=diff_right_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=diff_right_compose, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=diff_right_compose, context=(:Input,), kwargs...)
    diff_right_compose!(
        G.lie_group,
        unwrap_validation(Y),
        unwrap_validation(g),
        unwrap_validation(h),
        unwrap_validation(X),
    )
    is_point(LieAlgebra(G), Y; within=diff_right_compose, context=(:Output,), kwargs...)
    return Y
end

function Base.exp(G::ValidationLieGroup, X; kwargs...)
    is_point(LieAlgebra(G), X; within=exp, context=(:Input,), kwargs...)
    g = exp(G.lie_group, unwrap_validation(X))
    is_point(G, g; within=exp, context=(:Output,), kwargs...)
    return ValidationMPoint(g)
end
function ManifoldsBase.exp!(G::ValidationLieGroup, g, X; kwargs...)
    is_point(LieAlgebra(G), X; within=exp, context=(:Input,), kwargs...)
    exp!(G.lie_group, unwrap_validation(g), unwrap_validation(X))
    is_point(G, g; within=exp, context=(:Output,), kwargs...)
    return g
end

function Base.exp(G::ValidationLieGroup, g, X; kwargs...)
    is_point(G, g; within=exp, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=exp, context=(:Input,), kwargs...)
    h = exp(G.lie_group, unwrap_validation(g), unwrap_validation(X))
    is_point(G, h; within=exp, context=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function ManifoldsBase.exp!(G::ValidationLieGroup, h, g, X; kwargs...)
    is_point(G, g; within=exp, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), X; within=exp, context=(:Input,), kwargs...)
    exp!(G.lie_group, unwrap_validation(h), unwrap_validation(g), unwrap_validation(X))
    is_point(G, h; within=exp, context=(:Output,), kwargs...)
    return g
end

function ManifoldsBase.hat(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, c; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    X = hat(LieAlgebra(G), c)
    is_point(𝔤, X; within=hat, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function ManifoldsBase.hat(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, c, T::Type; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    X = hat(LieAlgebra(G), c, unwrap_validation_type(T))
    is_point(𝔤, X; within=hat, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function ManifoldsBase.hat!(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X, c; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    hat!(LieAlgebra(G), unwrap_validation(X), c)
    is_point(𝔤, X; within=hat, context=(:Output,), kwargs...)
    return X
end

function identity_element(G::ValidationLieGroup, ::Type{T}; kwargs...) where {T}
    g = identity_element(G.lie_group, unwrap_validation_type(T))
    is_point(G, g; within=identity_element, context=(:Output,), kwargs...)
    return ValidationMPoint(g)
end
function identity_element!(G::ValidationLieGroup, g; kwargs...)
    identity_element!(G.lie_group, unwrap_validation(g))
    is_point(G, g; within=identity_element, context=(:Output,), kwargs...)
    return g
end

function Base.inv(G::ValidationLieGroup, g; kwargs...)
    is_point(G, g; within=inv, context=(:Input,), kwargs...)
    h = inv(G.lie_group, unwrap_validation(g))
    is_point(G, h; within=inv, context=(:Output,), kwargs...)
    return ValidationMPoint(h)
end
function Base.inv(
    G::LieGroups.ValidationLieGroup{𝔽,O}, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_identity(G.lie_group, e; kwargs...)
    return inv(G.lie_group, e)
end
function inv!(G::ValidationLieGroup, h, g; kwargs...)
    is_point(G, g; within=inv, context=(:Input,), kwargs...)
    inv!(G.lie_group, unwrap_validation(h), unwrap_validation(g))
    is_point(G, h; within=inv, context=(:Output,), kwargs...)
    return h
end
function inv!(
    G::ValidationLieGroup{𝔽,O}, h, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_identity(G.lie_group, e; kwargs...)
    inv!(G.lie_group, unwrap_validation(h), e)
    is_point(G, h; within=inv, context=(:Output,), kwargs...)
    return h
end

function inv_left_compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; within=inv_left_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=inv_left_compose, context=(:Input,), kwargs...)
    k = inv_left_compose(G.lie_group, unwrap_validation(g), unwrap_validation(h))
    is_point(G, k; within=inv_left_compose, context=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function inv_left_compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; within=inv_left_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=inv_left_compose, context=(:Input,), kwargs...)
    inv_left_compose!(
        G.lie_group, unwrap_validation(k), unwrap_validation(g), unwrap_validation(h)
    )
    is_point(G, k; within=inv_left_compose, context=(:Output,), kwargs...)
    return k
end

function inv_right_compose(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; within=inv_right_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=inv_right_compose, context=(:Input,), kwargs...)
    k = inv_right_compose(G.lie_group, unwrap_validation(g), unwrap_validation(h))
    is_point(G, k; within=inv_right_compose, context=(:Output,), kwargs...)
    return ValidationMPoint(k)
end
function inv_right_compose!(G::ValidationLieGroup, k, g, h; kwargs...)
    is_point(G, g; within=inv_right_compose, context=(:Input,), kwargs...)
    is_point(G, h; within=inv_right_compose, context=(:Input,), kwargs...)
    inv_right_compose!(
        G.lie_group, unwrap_validation(k), unwrap_validation(g), unwrap_validation(h)
    )
    is_point(G, k; within=inv_right_compose, context=(:Output,), kwargs...)
    return k
end

is_identity(G::ValidationLieGroup, g) = is_identity(G.lie_group, unwrap_validation(g))
function is_identity(
    G::ValidationLieGroup{𝔽,O}, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return is_identity(G.lie_group, e)
end
function is_identity(
    G::ValidationLieGroup{𝔽,<:AbstractGroupOperation},
    e::Identity{<:AbstractGroupOperation};
    kwargs...,
) where {𝔽}
    return is_identity(G.lie_group, e)
end
function ManifoldsBase.is_point(
    G::ValidationLieGroup,
    g;
    error::Symbol=G.mode,
    within::Union{Nothing,Function}=nothing,
    context::NTuple{N,Symbol} where {N}=(),
    kwargs...,
)
    !_vLc(G, within, (:Point, context...)) && return true
    return is_point(G.lie_group, unwrap_validation(g); error=error, kwargs...)
end
function ManifoldsBase.is_point(
    G::ValidationLieGroup,
    e::Identity;
    error::Symbol=G.mode,
    within::Union{Nothing,Function}=nothing,
    context::NTuple{N,Symbol} where {N}=(),
    kwargs...,
)
    !_vLc(G, within, (:Point, context...)) && return true
    return is_point(G.lie_group, e; error=error, kwargs...)
end
function ManifoldsBase.is_point(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup},
    X;
    error::Symbol=base_lie_group(𝔤).mode,
    within::Union{Nothing,Function}=nothing,
    context::NTuple{N,Symbol} where {N}=(),
    kwargs...,
) where {𝔽,O<:AbstractGroupOperation}
    vG = base_lie_group(𝔤)
    !_vLc(vG, within, (:Vector, context...)) && return true
    return is_point(LieAlgebra(vG.lie_group), unwrap_validation(X); error=error, kwargs...)
end
function ManifoldsBase.is_vector(
    G::ValidationLieGroup{𝔽,O}, e::Identity{O}, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return ManifoldsBase.is_vector(G.lie_group, e, unwrap_validation(X); kwargs...)
end
function ManifoldsBase.is_vector(G::ValidationLieGroup, g, X; kwargs...)
    return ManifoldsBase.is_vector(
        G.lie_group, unwrap_validation(g), unwrap_validation(X); kwargs...
    )
end

function ManifoldsBase.isapprox(G::ValidationLieGroup, g, h; kwargs...)
    is_point(G, g; within=isapprox, context=(:Input,), kwargs...)
    is_point(G, h; within=isapprox, context=(:Input,), kwargs...)
    return isapprox(G.lie_group, unwrap_validation(g), unwrap_validation(h); kwargs...)
end
function ManifoldsBase.isapprox(
    G::ValidationLieGroup{𝔽,O}, g::Identity{O}, h; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_point(G, h; within=isapprox, context=(:Input,), kwargs...)
    return isapprox(G.lie_group, g, unwrap_validation(h); kwargs...)
end
function ManifoldsBase.isapprox(
    G::ValidationLieGroup{𝔽,O}, g, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_point(G, g; within=isapprox, context=(:Input,), kwargs...)
    return isapprox(G.lie_group, unwrap_validation(g), h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::ValidationLieGroup{𝔽,O}, g::Identity{O}, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    return isapprox(G.lie_group, g, h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::ValidationLieGroup{𝔽,O}, g::Identity{O}, h::Identity{O2}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,O2<:AbstractGroupOperation}
    return isapprox(G.lie_group, g, h; kwargs...)
end
function ManifoldsBase.isapprox(
    G::ValidationLieGroup{𝔽,O}, g::Identity{O2}, h::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,O2<:AbstractGroupOperation}
    return isapprox(G.lie_group, g, h; kwargs...)
end
function ManifoldsBase.isapprox(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X, Y; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    _X = unwrap_validation(X)
    _Y = unwrap_validation(Y)
    is_point(LieAlgebra(G), _X; within=isapprox, context=(:Input,), kwargs...)
    is_point(LieAlgebra(G), _Y; within=isapprox, context=(:Input,), kwargs...)
    return isapprox(LieAlgebra(G), _X, _Y; kwargs...)
end

function jacobian_conjugate(
    G::ValidationLieGroup,
    g,
    h,
    B::AbstractBasis=DefaultLieAlgebraOrthogonalBasis();
    kwargs...,
)
    is_point(G, g; within=jacobian_conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=jacobian_conjugate, context=(:Input,), kwargs...)
    J = jacobian_conjugate(G.lie_group, unwrap_validation(g), unwrap_validation(h))
    return J
end
function jacobian_conjugate!(
    G::ValidationLieGroup,
    J,
    g,
    h,
    B::AbstractBasis=DefaultLieAlgebraOrthogonalBasis();
    kwargs...,
)
    is_point(G, g; within=jacobian_conjugate, context=(:Input,), kwargs...)
    is_point(G, h; within=jacobian_conjugate, context=(:Input,), kwargs...)
    jacobian_conjugate!(G.lie_group, J, unwrap_validation(g), unwrap_validation(h))
    return J
end

function lie_bracket(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X, Y; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    is_point(𝔤, X; within=lie_bracket, context=(:Input,), kwargs...)
    is_point(𝔤, X; within=lie_bracket, context=(:Input,), kwargs...)
    Z = lie_bracket(LieAlgebra(G), unwrap_validation(X), unwrap_validation(Y))
    is_point(𝔤, Z; within=lie_bracket, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(Z)
end
function lie_bracket!(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, Z, X, Y; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    is_point(𝔤, X; within=lie_bracket, context=(:Input,), kwargs...)
    is_point(𝔤, X; within=lie_bracket, context=(:Input,), kwargs...)
    lie_bracket!(
        LieAlgebra(G), unwrap_validation(Z), unwrap_validation(X), unwrap_validation(Y)
    )
    is_point(𝔤, X; within=lie_bracket, context=(:Output,), kwargs...)
    return Z
end

function Base.log(G::ValidationLieGroup, g; kwargs...)
    is_point(G, g; within=log, context=(:Input,), kwargs...)
    X = log(G.lie_group, unwrap_validation(g))
    is_point(LieAlgebra(G), X; within=log, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function Base.log(
    G::ValidationLieGroup{𝔽,O}, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    X = log(G.lie_group, e)
    is_point(LieAlgebra(G), X; within=log, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function Base.log(
    G::ValidationLieGroup{𝔽,O}, e::Identity{O}, ::Type{T}; kwargs...
) where {𝔽,O<:AbstractGroupOperation,T}
    X = log(G.lie_group, e, T)
    is_point(LieAlgebra(G), X; within=log, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function Base.log(
    G::ValidationLieGroup{𝔽,O},
    e::Identity{O},
    ::Type{<:ValidationLieAlgebraTangentVector{T}};
    kwargs...,
) where {𝔽,O<:AbstractGroupOperation,T}
    X = log(G.lie_group, e, T)
    is_point(LieAlgebra(G), X; within=log, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function ManifoldsBase.log!(G::ValidationLieGroup, X, g; kwargs...)
    is_point(G, g; within=log, context=(:Input,), kwargs...)
    log!(G.lie_group, unwrap_validation(X), unwrap_validation(g))
    is_point(LieAlgebra(G), X; within=log, context=(:Input,), kwargs...)
    return X
end
function ManifoldsBase.log!(
    G::ValidationLieGroup{𝔽,O}, X, e::Identity{O}; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    log!(G.lie_group, X, e)
    is_point(LieAlgebra(G), X; within=log, context=(:Output,), kwargs...)
    return ValidationLieAlgebraTangentVector(X)
end
function ManifoldsBase.log!(G::ValidationLieGroup, X, g, h; kwargs...)
    is_point(G, g; within=log, context=(:Input,), kwargs...)
    is_point(G, h; within=log, context=(:Input,), kwargs...)
    log!(G.lie_group, unwrap_validation(X), unwrap_validation(g), unwrap_validation(h))
    is_point(LieAlgebra(G), X; within=log, context=(:Input,), kwargs...)
    return X
end

function LinearAlgebra.norm(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    is_point(𝔤, X; within=norm, context=(:Input,), kwargs...)
    return norm(LieAlgebra(G), unwrap_validation(X))
end
function LinearAlgebra.norm(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X::Real; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    is_point(𝔤, X; within=norm, context=(:Input,), kwargs...)
    return norm(LieAlgebra(G), X)
end

function Base.rand(G::ValidationLieGroup; vector_at=nothing, kwargs...)
    if vector_at !== nothing
        is_point(G, vector_at; within=rand, context=(:Input,), kwargs...)
    end
    gX = rand(G.lie_group; vector_at=vector_at, kwargs...)
    if vector_at !== nothing
        is_point(LieAlgebra(G), gX; within=rand, context=(:Output,), kwargs...)
    else
        is_point(G, gX; within=rand, context=(:Output,), kwargs...)
    end
    return gX
end

function Random.rand!(
    rng::AbstractRNG, G::ValidationLieGroup, gX; vector_at=nothing, kwargs...
)
    if vector_at !== nothing
        is_point(G, vector_at; within=rand, context=(:Input,), kwargs...)
        rand!(
            rng,
            G.lie_group,
            unwrap_validation(gX);
            vector_at=unwrap_validation(vector_at),
            kwargs...,
        )
    else
        rand!(rng, G.lie_group, unwrap_validation(gX); kwargs...)
    end
    if vector_at !== nothing
        is_point(LieAlgebra(G), gX; within=rand, context=(:Output,), kwargs...)
    else
        is_point(G, gX; within=rand, context=(:Output,), kwargs...)
    end
    return gX
end

ManifoldsBase.manifold_dimension(G::ValidationLieGroup) = manifold_dimension(G.lie_group)

ManifoldsBase.representation_size(G::ValidationLieGroup) = representation_size(G.lie_group)

function Base.show(io::IO, G::ValidationLieGroup)
    s = """
    ValidationLieGroup of $(G.lie_group)
        * mode = :$(G.mode)
    """
    G_ig = G.ignore_contexts
    (length(G_ig) > 0) && (s *= "    * ignore_context = $(G_ig)\n")
    G_if = G.ignore_functions
    (length(G_if) > 0) && (s *= "    * ignore_functions = $(G_if)")
    return print(io, s)
end

function ManifoldsBase.vee(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_point(𝔤, X; within=vee, context=(:Input,), kwargs...)
    G = base_lie_group(𝔤).lie_group
    return vee(LieAlgebra(G), unwrap_validation(X))
end
function ManifoldsBase.vee!(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, c, X; kwargs...
) where {𝔽,O<:AbstractGroupOperation}
    is_point(𝔤, X; within=vee, context=(:Input,), kwargs...)
    G = base_lie_group(𝔤).lie_group
    vee!(LieAlgebra(G), unwrap_validation(c), unwrap_validation(X))
    return X
end

function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, T::Type
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    return ValidationLieAlgebraTangentVector(zero_vector(LieAlgebra(G), T))
end
function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, ::Type{ValidationLieAlgebraTangentVector{T}}
) where {𝔽,O<:AbstractGroupOperation,T}
    G = base_lie_group(𝔤).lie_group
    return ValidationLieAlgebraTangentVector(zero_vector(LieAlgebra(G), T))
end
function ManifoldsBase.zero_vector(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}
) where {𝔽,O<:AbstractGroupOperation}
    G = base_lie_group(𝔤).lie_group
    return ValidationLieAlgebraTangentVector(zero_vector(LieAlgebra(G)))
end

function ManifoldsBase.zero_vector!(
    𝔤::LieAlgebra{𝔽,O,<:ValidationLieGroup}, X::T
) where {𝔽,O<:AbstractGroupOperation,T}
    G = base_lie_group(𝔤).lie_group
    T2 = typeof(unwrap_validation(X))
    return ValidationLieAlgebraTangentVector(zero_vector(G, T2))
end
