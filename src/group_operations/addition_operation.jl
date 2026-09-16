"""
    AdditionGroupOperation <: AbstractGroupOperation

A group operation that is realised introducing defaults that fall back
to `+` and `-` being overloaded, for example
`_compose(G::LieGroup{𝔽,AdditionGroupOperation}, a, b) = a + b`
"""
struct AdditionGroupOperation <: AbstractGroupOperation end

#
#
# Handle interactions of `+` and `-` with the identity element, though they are
# also already handled on the `compose()` level
Base.:+(e::Identity{AdditionGroupOperation}) = e
Base.:+(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:+(::Identity{AdditionGroupOperation}, g) = g
Base.:+(g, ::Identity{AdditionGroupOperation}) = g

Base.:-(e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::Identity{AdditionGroupOperation}) = e
Base.:-(::Identity{AdditionGroupOperation}, g) = -g
Base.:-(g, ::Identity{AdditionGroupOperation}) = g

# Identity{AdditionGroupOperation} is a valid point on a Lie algebra and the same as the zero vector.
Base.:+(e::Identity{AdditionGroupOperation}, ::ManifoldsBase.ZeroVector) = e
Base.:+(::ManifoldsBase.ZeroVector, e::Identity{AdditionGroupOperation}) = e
Base.:-(e::Identity{AdditionGroupOperation}, ::ManifoldsBase.ZeroVector) = e
Base.:-(::ManifoldsBase.ZeroVector, e::Identity{AdditionGroupOperation}) = e

_doc_compose_add = """
    compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h)
    compose!(G::LieGroup{𝔽,AdditionGroupOperation}, k, g, h)

Compute the group operation composition of `g` and `h` with respect to
the [`AdditionGroupOperation`](@ref) on `G`, which falls back to calling
`g+h`, where `+` is assumed to be overloaded accordingly.

This can be computed in-place of `k`.
"""

@doc "$(_doc_compose_add)"
compose(::LieGroup{𝔽, AdditionGroupOperation}, g, h) where {𝔽}

@doc "$(_doc_compose_add)"
compose!(::LieGroup{𝔽, AdditionGroupOperation}, k, g, h) where {𝔽}

function _compose!(G::LieGroup{𝔽, AdditionGroupOperation}, k, g, h) where {𝔽}
    k .= g .+ h
    return k
end

_doc_diff_conjugate_add = """
    diff_conjugate(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X)
    diff_conjugate!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X)

Compute the differential of the conjugate ``c_g(h) = g$(_math(:∘))h$(_math(:∘))g^{-1} = g+h-g = h``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``$(_math(:d))(c_g(h))[X] = X``.
"""

@doc "$(_doc_diff_conjugate_add)"
diff_conjugate(G::LieGroup{𝔽, AdditionGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_conjugate_add)"
function diff_conjugate!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_diff_inv_add = """
    diff_inv(G::LieGroup{𝔽,AdditionGroupOperation}, g, X)
    diff_inv!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, X)

Compute the differential of the inverse operation ``ι_{$(_math(:G))}(g) = g^-1 = -g``,
which simplifies for [`AdditionGroupOperation`](@ref) to ``$(_math(:d))ι_{$(_math(:G))}(g)[X] = -X``
"""

@doc "$(_doc_diff_inv_add)"
diff_inv(G::LieGroup{𝔽, AdditionGroupOperation}, g, X) where {𝔽} = -X

@doc "$(_doc_diff_inv_add)"
function diff_inv!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, X) where {𝔽}
    Y .= (-1) .* X
    return Y
end

_doc_diff_left_compose_add = """
    diff_left_compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X)
    diff_left_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X)

Compute the differential of the the group operation ``g+h`` with respect to the left argument `g`.
Here it simplifies for [`AdditionGroupOperation`](@ref) to ``$(_math(:d))ρ_h(g)[X] = X``.
"""

@doc "$(_doc_diff_left_compose_add)"
diff_left_compose(G::LieGroup{𝔽, AdditionGroupOperation}, g, h, X) where {𝔽} = X
diff_left_compose(G::LieGroup{𝔽, AdditionGroupOperation}, ::Identity{AdditionGroupOperation}, h, X) where {𝔽} = X

@doc "$(_doc_diff_left_compose_add)"
function diff_left_compose!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_diff_right_compose_add = """
    diff_right_compose(G::LieGroup{𝔽,AdditionGroupOperation}, g, h, X)
    diff_right_compose!(G::LieGroup{𝔽,AdditionGroupOperation}, Y, g, h, X)

Compute the differential of the group operation ``g$(_math(:∘))h``, on an [`AbstractLieGroup`](@ref) `G`
with respect to its second (right) argument `h`.

Another interpretation is to consider a function where we do a fixed multiplication from the left with `g`.
i.e. the left group multiplication function ``λ_g(h) = g$(_math(:∘))h`` (where the _left_ refers to the fixed argument ``g``.).

For the [`AdditionGroupOperation`](@ref) it reads ``$(_math(:d))λ_g(h)[X] = X``.
"""

@doc "$(_doc_diff_right_compose_add)"
diff_right_compose(::LieGroup{𝔽, AdditionGroupOperation}, g, h, X) where {𝔽}

@doc "$(_doc_diff_right_compose_add)"
function diff_right_compose!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, h, X) where {𝔽}
    return copyto!(LieAlgebra(G), Y, X)
end

_doc_exp_add = """
    exp(G::LieGroup{𝔽,AdditionGroupOperation}, X)
    exp!(G::LieGroup{𝔽,AdditionGroupOperation}, g, X)

Compute the Lie group exponential on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `g`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``g=0+X=X``.
"""

@doc "$(_doc_exp_add)"
ManifoldsBase.exp(::LieGroup{𝔽, AdditionGroupOperation}, X) where {𝔽} = X

@doc "$(_doc_exp_add)"
function ManifoldsBase.exp!(::LieGroup{𝔽, AdditionGroupOperation}, g, X) where {𝔽}
    g .= X
    return g
end

@inline function get_vector_lie(
        ::LieAlgebra{ℝ, AdditionGroupOperation},
        c,
        B::DefaultLieAlgebraOrthogonalBasis{ℝ},
        T::Type{<:SArray},
    )
    return convert(T, c)
end

_doc_jacobian_exp_add = """
    jacobian_exp(G::LieGroup{𝔽,AdditionGroupOperation}, X, ::DefaultLieAlgebraOrthogonalBasis)
    jacobian_exp!(G::LieGroup{𝔽,AdditionGroupOperation}, J, X, ::DefaultLieAlgebraOrthogonalBasis)

Compute the Jacobian of the Lie group exponential in a basis of the Lie algebra on a
[`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).

Since such a group is Abelian and flat, the differential of the exponential map is the
identity, so ``J = I_n`` for every ``X``.
"""

@doc "$(_doc_jacobian_exp_add)"
jacobian_exp(::LieGroup{𝔽, AdditionGroupOperation}, X, basis = DefaultLieAlgebraOrthogonalBasis()) where {𝔽}

@doc "$(_doc_jacobian_exp_add)"
function jacobian_exp!(
        ::LieGroup{𝔽, AdditionGroupOperation}, J::AbstractMatrix, X, ::DefaultLieAlgebraOrthogonalBasis
    ) where {𝔽}
    copyto!(J, LinearAlgebra.I)
    return J
end

_doc_identity_element_add = """
    identity_element(G::LieGroup{𝔽,AdditionGroupOperation})
    identity_element!(G::LieGroup{𝔽,AdditionGroupOperation}, e)

Return the a point representation of the [`Identity`](@ref),
which for the [`AdditionGroupOperation`](@ref) is the zero element or array.
"""

@doc "$(_doc_identity_element_add)"
identity_element(::LieGroup{𝔽, AdditionGroupOperation}) where {𝔽}

function identity_element(
        G::LieGroup{𝔽, AdditionGroupOperation}, ::Type{T}
    ) where {𝔽, T <: AbstractArray}
    return zeros(representation_size(G.manifold))
end
function identity_element(
        ::LieGroup{𝔽, AdditionGroupOperation}, ::Type{T}
    ) where {𝔽, T <: Union{Number, AbstractArray{<:Number, 0}}}
    return zero(T)
end
function identity_element(
        ::LieGroup{𝔽, AdditionGroupOperation}, ::Type{Array{T, 0}}
    ) where {𝔽, T <: Number}
    return fill(zero(T))
end
function identity_element(
        ::LieGroup{𝔽, AdditionGroupOperation}, T::Type{<:StaticArray}
    ) where {𝔽}
    return zero(T)
end

@doc "$(_doc_identity_element_add)"
function identity_element!(::LieGroup{𝔽, AdditionGroupOperation}, e) where {𝔽}
    return fill!(e, 0)
end

function ManifoldsBase.inner(::LieAlgebra{ℝ, AdditionGroupOperation}, X, Y)
    return dot(X, Y)
end

_doc_inv_add = """
    inv(G::LieGroup{𝔽,AdditionGroupOperation}, g)
    inv!(G::LieGroup{𝔽,AdditionGroupOperation}, h, g)

Compute the inverse group element ``g^{-1}``, which for the [`AdditionGroupOperation`](@ref)
simplifies to ``-g``. This can be done in-place of `h`.
"""

@doc "$(_doc_inv_add)"
Base.inv(G::LieGroup{𝔽, AdditionGroupOperation}, g) where {𝔽}

@doc "$(_doc_inv_add)"
inv!(::LieGroup{𝔽, AdditionGroupOperation}, h, g) where {𝔽}
function _inv!(::LieGroup{𝔽, AdditionGroupOperation}, h, g) where {𝔽}
    h .= (-1) .* g
    return h
end

_doc_lie_bracket_add = """
    lie_bracket(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, X, Y)
    lie_bracket!(𝔤::LieAlgebra{𝔽,AdditionGroupOperation}, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))``,
which for the for the [`AdditionGroupOperation`](@ref) simplifies to the
corresponding $(_link(:zero_vector)).
The computation can be done in-place of `Z`.
"""

@doc "$(_doc_lie_bracket_add)"
lie_bracket(𝔤::LieAlgebra{𝔽, AdditionGroupOperation}, X, Y) where {𝔽}

@doc "$(_doc_lie_bracket_add)"
function lie_bracket!(
        𝔤::LieAlgebra{𝔽, O, <:LieGroup{𝔽, O}}, Z, X, Y
    ) where {𝔽, O <: AdditionGroupOperation}
    return zero_vector!(𝔤, Z)
end

_doc_log_add = """
    log(G::LieGroup{𝔽,AdditionGroupOperation}, g)
    log!(G::LieGroup{𝔽,AdditionGroupOperation}, X, g)

Compute the Lie group logarithm on a [`LieGroup`](@ref) with an [`AdditionGroupOperation`](@ref).
This can be computed in-place of `X`.

Since `e` is just the zero-element with respect to the corresponding `+`, the formula reads ``X=g-0=g``.
"""

@doc "$(_doc_log_add)"
ManifoldsBase.log(::LieGroup{𝔽, AdditionGroupOperation}, q) where {𝔽} = q
function ManifoldsBase.log(
        G::LieGroup{𝔽, AdditionGroupOperation}, e::Identity{AdditionGroupOperation}
    ) where {𝔽}
    return zero_vector(LieAlgebra(G))
end
@doc "$(_doc_log_add)"
function ManifoldsBase.log!(G::LieGroup{𝔽, AdditionGroupOperation}, X, g) where {𝔽}
    return copyto!(G, X, g)
end
function ManifoldsBase.log!(
        ::LieGroup{𝔽, AdditionGroupOperation}, X, ::Identity{AdditionGroupOperation}
    ) where {𝔽}
    return fill!(X, 0)
end

function push_forward_tangent!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, X) where {𝔽}
    return copyto!(Y, X)
end
function pull_back_tangent!(G::LieGroup{𝔽, AdditionGroupOperation}, Y, g, X) where {𝔽}
    return copyto!(Y, X)
end
