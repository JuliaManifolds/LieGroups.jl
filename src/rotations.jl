# Algebra Interfaces

abstract type AbstractRotationAlgebra{N} <: AbstractLieAlgebra end

dim(::Type{<:AbstractRotationAlgebra{N}}) where {N} = N
dim(::AbstractRotationAlgebra{N}) where {N} = N

dof(::Type{<:AbstractRotationAlgebra{N}}) where {N} = sum(1:(N-1))
dof(::AbstractRotationAlgebra{N}) where {N} = sum(1:(N-1))

struct so{N,V} <: AbstractRotationAlgebra{N}
    θ::V

    function so{N}(x::T) where {N,T<:AbstractVector}
        check_dof(so{N}, length(x))
        return new{N,T}(x)
    end
    
    function so{N}(X::T) where {N,T<:AbstractMatrix}
        check_dim(so{N}, size(X, 1))
        check_skewsymmetric(X)
        return new{N,T}(X)
    end
end

Base.angle(alg::so) = alg.θ

(==)(alg1::so{N}, alg2::so{N}) where {N} = alg1.θ == alg2.θ
Base.isapprox(alg1::so{N}, alg2::so{N}) where {N} = isapprox(alg1.θ, alg2.θ)

identity(alg::so{N,T}) where {N,T<:AbstractVector} =
    so{N}(fill!(similar(alg.θ), 0))

inv(alg::so{N,T}) where {N,T<:AbstractVector} = so{N}(-alg.θ)

(+)(alg1::so{N}, alg2::so{N}) where {N} = so{N}(alg1.θ + alg2.θ)

∧(alg::so{N,T}) where {N,T<:AbstractVector} = so{N}(∧(so{N}, alg.θ))
∧(alg::so{N,T}) where {N,T<:AbstractMatrix} = alg

∨(alg::so{N,T}) where {N,T<:AbstractVector} = alg
∨(alg::so{N,T}) where {N,T<:AbstractMatrix} = so{N}(∨(so{N}, alg.θ))

Base.Vector(alg::so{N,T}) where {N,T<:AbstractVector} = alg.θ
Base.Vector(alg::so{N,T}) where {N,T<:AbstractMatrix} = ∨(so{N}, alg.θ)

Base.Matrix(alg::so{N,T}) where {N,T<:AbstractVector} = ∧(so{N}, alg.θ)
Base.Matrix(alg::so{N,T}) where {N,T<:AbstractMatrix} = alg.θ

function left_jacobian(alg::so{3})
    θ² = sum(abs2, alg.θ)
    θ = √θ²
    W = skewsymmetric(alg.θ)
    M1 = (1. - cos(θ))/(θ²) * W
    M2 = (θ - sin(θ))/(θ² * θ) * W^2
    return I(3) + M1 + M2
end

right_jacobian(alg::so{3}) = left_jacobian(alg)'


# Group Interfaces

abstract type AbstractRotationGroup{N} <: AbstractLieGroup end

dim(::Type{<:AbstractRotationGroup{N}}) where {N} = N
dim(::AbstractRotationGroup{N}) where {N} = N

dof(::Type{<:AbstractRotationGroup{N}}) where {N} = sum(1:(N-1))
dof(::AbstractRotationGroup{N}) where {N} = sum(1:(N-1))

struct SO{N,T} <: AbstractRotationGroup{N}
    R::T

    function SO{N}(R::T) where {N,T<:AbstractMatrix}
        @assert size(R, 1) == N
        return new{N,T}(R)
    end
end

rotation(g::SO) = Array(g.R)

identity(::SO{N}) where {N} = SO{N}(I(N))
identity(::Type{SO{N}}) where {N} = SO{N}(I(N))

inv(g::SO{N}) where {N} = SO{N}(inv(rotation(g)))

(*)(::SO{M}, ::SO{N}) where {M,N} =
    throw(ArgumentError("* operation for SO{$M} and SO{$N} group is not defined."))

(*)(g1::SO{N}, g2::SO{N}) where {N} = SO{N}(rotation(g1) * rotation(g2))

(==)(g1::SO{N}, g2::SO{N}) where {N} = Matrix(g1) == Matrix(g2)
Base.isapprox(g1::SO{N}, g2::SO{N}) where {N} = isapprox(Matrix(g1), Matrix(g2))

Base.Matrix(g::SO) = rotation(g)

Base.show(io::IO, g::SO{N}) where {N} = print(io, "SO{$N}(R=", rotation(g), ")")

⋉(g::SO{N}, x::T) where {N,T<:AbstractVector} = Matrix(g) * x


# Array Interfaces

function ∧(::Type{so{N}}, θ::AbstractVector) where {N}
    check_dof(so{N}, length(θ))
    return skewsymmetric(θ)
end

function ∨(::Type{so{N}}, Θ::AbstractMatrix) where {N}
    check_dim(so{N}, size(Θ, 1))
    if N == 2
        return [Θ[2, 1]]
    elseif N == 3
        return [Θ[3, 2], Θ[1, 3], Θ[2, 1]]
    else
        throw(ArgumentError("not support."))
    end
end


# Maps

Base.exp(alg::so{N,T}) where {N,T<:AbstractMatrix} = SO{N}(exp(alg.θ))
Base.exp(alg::so{N,T}) where {N,T<:AbstractVector} = SO{N}(exp(∧(so{N}, alg.θ)))
Base.log(g::SO{N}) where {N} = so{N}(∨(so{N}, log(rotation(g))))
