
"""
    LieAlgebra{𝔽, G} <: AbstractManifold{𝔽}

Represent the Lie algebra ``$(_math(:𝔤))``, that is a ``𝔽`` vector space with an associated
[`lie_bracket`](@ref) ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The Lie algebras considered here are those related to a [`LieGroup`](@ref) ``$(_math(:G))``,
namely the tangent space ``T_{$(_math(:e))}$(_math(:G))`` at the [`Identity`](@ref),
this is internally just a `const` of the corresponding $(_link(:TangentSpace)).

# Constructor

    LieAlgebra(G::LieGroup)

Return the Lie Algebra belonging to the [`LieGroup`](@ref) `G`.
"""
const LieAlgebra{𝔽,O<:AbstractGroupOperation,G<:LieGroup{𝔽,O}} = ManifoldsBase.Fiber{
    𝔽,ManifoldsBase.TangentSpaceType,G,Identity{O}
}

function LieAlgebra(G::LieGroup{𝔽,O}) where {𝔽,O<:AbstractGroupOperation}
    return LieAlgebra{𝔽,O,typeof(G)}(G, Identity(G), ManifoldsBase.TangentSpaceType())
end

_doc_lie_bracket = """
    lie_bracket!(𝔤::LieAlgebra, X, Y)
    lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

Compute the Lie bracket ``[⋅,⋅]: $(_math(:𝔤))×$(_math(:𝔤)) → $(_math(:𝔤))`` which fulfills

1. ``[X,X] = 0`` for all ``X ∈ $(_math(:𝔤))``
2. The Jacobi identity ``[X, [Y,Z]] = [[X,Y],Z] = [Y, [X,Z]]`` holds for all ``X, Y, Z ∈ $(_math(:𝔤))``.

The computation can be done in-place of `Z`.
"""
function lie_bracket end
@doc "$(_doc_lie_bracket)"
function lie_bracket(𝔤::LieAlgebra, X, Y)
    Z = ManifoldsBase.allocate_result(𝔤, lie_bracket, X, Y)
    return lie_bracket!(𝔤, Z, X, Y)
end

function lie_bracket! end
@doc "$(_doc_lie_bracket)"
lie_bracket!(𝔤::LieAlgebra, Z, X, Y)

"""
    is_point(𝔤::LieAlgebra, X; kwargs...)

Check whether `X` is a valid point on the Lie Algebra `𝔤`.
This falls back to checking whether `X` is a valid point on the tangent space
at the [`identity_element`](@ref)`(G)` on `G.manifold` on the [`LieGroup`](@ref)
of `G`
"""
function ManifoldsBase.is_point(𝔤::LieAlgebra, X; kwargs...)
    # the manifold stored in the Fiber / Lie algebra is the Lie group G
    G = 𝔤.manifold
    e = identity_element(G)
    return ManifoldsBase.is_vector(G.manifold, e, X; kwargs...)
end

LinearAlgebra.norm(𝔤::LieAlgebra, X::Real) = LinearAlgebra.norm(𝔤.manifold, 𝔤.point, X)
function LinearAlgebra.norm(
    G::LieGroup{𝔽,O}, ::Identity{O}, X
) where {𝔽,O<:AbstractGroupOperation}
    return LinearAlgebra.norm(G, identity_element(G), X)
end

function Base.show(io::IO, 𝔤::LieAlgebra)
    return print(io, "LieAlgebra( $(𝔤.manifold) )")
end

function ManifoldsBase.zero_vector(𝔤::LieAlgebra)
    return ManifoldsBase.zero_vector(𝔤.manifold, identity_element(𝔤.manifold))
end

function ManifoldsBase.zero_vector!(𝔤::LieAlgebra, X)
    return ManifoldsBase.zero_vector!(𝔤.manifold, X, identity_element(𝔤.manifold))
end
