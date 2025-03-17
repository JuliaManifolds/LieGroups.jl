"""
    UnitaryGroup{T}

The special orthogonal group ``$(_math(:U))(n)`` is the Lie group consisting of the [`MatrixMultiplicationGroupOperation`](@ref) on the
manifold of rotations [`UnitaryMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`) with absolute value of the determinant equal to one.

# Constructor
    UnitaryGroup(n::Int, 𝔽::AbstractNumbers=ℂ; kwargs...)

Generate unitary group ``$(_math(:U))(n)``.
All keyword arguments in `kwargs...` are passed on to [`Rotations`](@extref `Manifolds.Rotations`) as well.
"""
const UnitaryGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,Manifolds.UnitaryMatrices{T,𝔽}
}

function UnitaryGroup(n::Int, 𝔽::AbstractNumbers=ManifoldsBase.ℂ; kwargs...)
    U = Manifolds.UnitaryMatrices(n, 𝔽; kwargs...)
    return UnitaryGroup{𝔽,typeof(U).parameters[1]}(U, MatrixMultiplicationGroupOperation())
end

function ManifoldsBase.check_size(
    ::UnitaryGroup{ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, p::Number
)
    return nothing
end
function ManifoldsBase.check_size(
    ::UnitaryGroup{ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, p, X::Number
)
    return nothing
end
function ManifoldsBase.check_size(
    G::UnitaryGroup{ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, p
)
    if size(p) != ()
        return DomainError(
            size(p),
            "The point $(p) can not belong to the manifold $(G), since its size $(size(p)) is not equal to the manifolds representation size for (mutable) Quaternions (()).",
        )
    end
    return nothing
end
function ManifoldsBase.check_size(
    G::UnitaryGroup{ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, p, X
)
    if size(X) != ()
        return DomainError(
            size(X),
            "The point $(p) can not belong to the manifold $(G), since its size $(size(X)) is not equal to the manifolds representation size for (mutable) Quaternions (()).",
        )
    end
    return nothing
end

function Base.exp(
    ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, X::Number
)
    return exp(X)
end
function Base.exp(
    ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}},
    g::Number,
    X::Number,
)
    return g * exp(X)
end

function identity_element(
    ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}
)
    return Quaternions.quat(1.0)
end

Base.inv(::UnitaryGroup, g) = adjoint(g)
inv!(G::UnitaryGroup, h, g) = copyto!(G, h, adjoint(g))
function inv!(::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, h, g)
    h[] = adjoint(g[])
    return h
end

function Base.log(
    ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, g::Number
)
    return log(g)
end
function Base.log(
    G::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}},
    g::Number,
    h::Number,
)
    return log(inv(G, g) * h)
end
function ManifoldsBase.log!(
    ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, X, g
)
    X[] = log(g[])
    return X
end
function ManifoldsBase.log!(
    G::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{1}}}, X, g, h
)
    X[] = log(inv(G, g) * h)
    return nothing
end
#
#
# A common type for all 4 groups: O, SO, SU, U, because they share quite some implementations
# Stored in this, since this is the most generic case
#
"""
    CommonUnitarySubGroup{𝔽,T}

A constant that allows to refer to several subgroups of ``$(_math(:U))(n)`` for
implementations where
* certain subgroups real/complex share a common implementation, e.g. for the same sizes `T` usually via [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)
* certain functions are the same for all sizes `T` as long as the field `𝔽` is the same
"""
const CommonUnitarySubGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,<:Manifolds.GeneralUnitaryMatrices{T,𝔽}
}

#
#
# A common type for all 4 groups: O, SO, SU, U, because they share quite some implementations
# Stored in this, since this is the most generic case
#
"""
    CommonUnitarySubAlgebra{𝔽,T}

A constant that allows to refer to several sub Algebras of ``$(_math(:u))(n)`` for
implementations where
* certain sub algebras real/complex share a common implementation, e.g. for the same sizes `T` usually via [`TypeParameter`](@extref `ManifoldsBase.TypeParameter`)
* certain functions are the same for all sizes `T` as long as the field `𝔽` is the same
"""
const CommonUnitarySubAlgebra{𝔽,T} = LieAlgebra{
    𝔽,MatrixMultiplicationGroupOperation,<:CommonUnitarySubGroup{𝔽,T}
}

function Base.show(
    io::IO, ::UnitaryGroup{ManifoldsBase.ℂ,ManifoldsBase.TypeParameter{Tuple{n}}}
) where {n}
    return print(io, "UnitaryGroup($(n))")
end
function Base.show(io::IO, M::UnitaryGroup{ManifoldsBase.ℂ,Tuple{Int}})
    n = ManifoldsBase.get_parameter(M.manifold.size)[1]
    return print(io, "UnitaryGroup($(n); parameter=:field)")
end
function Base.show(
    io::IO, ::UnitaryGroup{ManifoldsBase.ℍ,ManifoldsBase.TypeParameter{Tuple{n}}}
) where {n}
    return print(io, "UnitaryGroup($(n), ℍ)")
end
function Base.show(io::IO, G::UnitaryGroup{ManifoldsBase.ℍ,Tuple{Int}})
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return print(io, "UnitaryGroup($(n), ℍ; parameter=:field)")
end
