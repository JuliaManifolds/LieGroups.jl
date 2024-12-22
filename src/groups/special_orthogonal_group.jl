"""
    SpecialOrthogonalGroup{T}

The special orthogonal group ``$(_math(:SO))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`Rotations`](@extref `Manifolds.Rotations`).

# Constructor
    SpecialOrthogonalGroup(n; kwargs...)

Generate  special orthogonal group ``$(_math(:SO))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const SpecialOrthogonalGroup{T} = LieGroup{
    ManifoldsBase.ℝ,MatrixMultiplicationGroupOperation,Manifolds.Rotations{T}
}

function SpecialOrthogonalGroup(n; kwargs...)
    R = Manifolds.Rotations(n; kwargs...)
    return SpecialOrthogonalGroup{typeof(R).parameters[1]}(
        R, MatrixMultiplicationGroupOperation()
    )
end

@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
@doc "$(_doc_exp_O2_id)"
ManifoldsBase.exp!(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    g,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)

@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)
@doc "$(_doc_exp_O3_id)"
ManifoldsBase.exp!(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    g,
    ::Identity{MatrixMultiplicationGroupOperation},
    X,
)

inv!(G::SpecialOrthogonalGroup, k, g) = copyto!(G, k, transpose(g))
function inv!(
    G::SpecialOrthogonalGroup, q, ::Identity{O}
) where {O<:AbstractMultiplicationGroupOperation}
    return identity_element!(G, q)
end

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

@doc "$(_doc_log_O2_id)"
ManifoldsBase.log!(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{2}}},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

@doc "$(_doc_log_O3_id)"
ManifoldsBase.log!(
    ::SpecialOrthogonalGroup{ManifoldsBase.TypeParameter{Tuple{3}}},
    X,
    ::Identity{MatrixMultiplicationGroupOperation},
    g,
)

function Base.show(io::IO, G::SpecialOrthogonalGroup)
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "SpecialOrthogonalGroup($(size))")
end
