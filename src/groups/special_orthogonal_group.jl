"""
    SpecialOrthogonalGroup{T}

The special orthogonal group ``$(_math(:SO))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`Rotations`](@extref `Manifolds.Rotations`).

# Constructor
    SpecialOrthogonalGroup(n::Int; kwargs...)

Generate  special orthogonal group ``$(_math(:SO))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialOrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ, MatrixMultiplicationGroupOperation, Rotations{T},
}

function SpecialOrthogonalGroup(n::Int; kwargs...)
    R = Rotations(n; kwargs...)
    return SpecialOrthogonalGroup{typeof(R).parameters[1]}(
        R, MatrixMultiplicationGroupOperation()
    )
end

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X)
@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g, X)

@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X)
@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g, X)

@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X)
@doc "$(_doc_exp_O4_id)"
ManifoldsBase.exp!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g, X)

@doc "$(_doc_get_coordinates_On)"
ManifoldsBase.get_coordinates(
    G::SpecialOrthogonalGroup,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_coordinates_On)"
ManifoldsBase.get_coordinates!(
    G::SpecialOrthogonalGroup,
    c,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_vector_On)"
ManifoldsBase.get_vector(
    G::SpecialOrthogonalGroup,
    ::Identity{MatrixMultiplicationGroupOperation},
    c,
    ::DefaultLieAlgebraOrthogonalBasis,
)

@doc "$(_doc_get_vector_On)"
ManifoldsBase.get_vector!(
    G::SpecialOrthogonalGroup,
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    c,
    ::DefaultLieAlgebraOrthogonalBasis,
)

function identity_element(G::SpecialOrthogonalGroup)
    return identity_element(G, AbstractMatrix{Float64})
end
function identity_element(G::SpecialOrthogonalGroup, ::Type{<:AbstractMatrix{T}}) where {T}
    e = zeros(T, ManifoldsBase.representation_size(G)...)
    return identity_element!(G, e)
end

_inv!(G::SpecialOrthogonalGroup, k, g) = copyto!(G, k, transpose(g))


_doc_jacobian_exp_SO2 = """
    jacobian_exp(M::Rotations{TypeParameter{Tuple{2}}}, g, X, ::DefaultLieAlgebraOrthogonalBasis)
    jacobian_exp!(M::Rotations{TypeParameter{Tuple{2}}}, J, g, X, ::DefaultLieAlgebraOrthogonalBasis)

Compute Jacobian of the Lie group exponential in a basis of the Lie algebra on the [`SpecialOrthogonalGroup`](@ref)`(2)` manifold.

It is equal to matrix ``[1]``, see [SolaDerayAtchuthan:2021](@cite), Appendix A.
"""

@doc "$(_doc_jacobian_exp_SO2)"
jacobian_exp(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, p, X, basis = DefaultLieAlgebraOrthogonalBasis())

@doc "$(_doc_jacobian_exp_SO2)"
function jacobian_exp!(
        ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, J::AbstractMatrix, p, X, ::DefaultLieAlgebraOrthogonalBasis
    )
    J .= 1
    return J
end

_doc_jacobian_exp_SO3 = raw"""
    jacobian_exp(M::Rotations{TypeParameter{Tuple{3}}}, p, X, ::DefaultLieAlgebraOrthogonalBasis)
    jacobian_exp!(M::Rotations{TypeParameter{Tuple{3}}}, J, p, X, ::DefaultLieAlgebraOrthogonalBasis)

Compute Jacobian of the Lie group exponential in a basis of the Lie algebra on the [`Rotations`](@ref)`(3)` manifold. The formula reads

````math
J = 𝕀 + \frac{\cos(θ) - 1}{θ^2} X + \frac{θ - \sin(θ)}{θ^3} X^2,
````

where ``θ`` is the norm of `X`.
It is adapted from [Chirikjian:2012](@cite), Eq. (10.86), to `LieGroups.jl` conventions.
"""

@doc "$(_doc_jacobian_exp_SO3)"
jacobian_exp(M::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, p, X, basis = DefaultLieAlgebraOrthogonalBasis())

@doc "$(_doc_jacobian_exp_SO3)"
function jacobian_exp!(
        M::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, J::AbstractMatrix, p, X, ::DefaultLieAlgebraOrthogonalBasis
    )
    θ = norm(M, p, X) / sqrt(2)
    copyto!(J, I)
    if θ ≉ 0
        a = (cos(θ) - 1) / θ^2
        b = (θ - sin(θ)) / θ^3
        J .+= a .* X .+ b .* (X^2)
    end
    return J
end

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, g)

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}}, X, g)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, g)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}}, X, g)

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, g)

@doc "$(_doc_log_O4_id)"
ManifoldsBase.log!(::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{4}}}, X, g)

#
# Since the manifold also uses the Lie algebra to represent tangent vectors push forward and pull back are the identity.
pull_back_tangent!(::SpecialOrthogonalGroup, Y, g, X; kwargs...) = copyto!(Y, X)
push_forward_tangent!(::SpecialOrthogonalGroup, Y, g, X; kwargs...) = copyto!(Y, X)

function Base.show(io::IO, G::SpecialOrthogonalGroup)
    size = get_parameter(G.manifold.size)[1]
    return print(io, "SpecialOrthogonalGroup($(size))")
end
