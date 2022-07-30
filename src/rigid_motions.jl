# Algebra Interfaces

abstract type SpecialEuclideanAlgebra{N} <: AbstractLieAlgebra end

dim(::Type{<:SpecialEuclideanAlgebra{N}}) where {N} = N
dim(::SpecialEuclideanAlgebra{N}) where {N} = N

dof(::Type{<:SpecialEuclideanAlgebra{N}}) where {N} = sum(1:N)
dof(::SpecialEuclideanAlgebra{N}) where {N} = sum(1:N)

struct se{N,T,S} <: SpecialEuclideanAlgebra{N}
    ρ::S
    θ::T

    function se{N}(ρ::T, θ::S) where {N,T<:AbstractVector,S<:AbstractVector}
        check_dof(so{N}, length(θ))
        check_dim(se{N}, length(ρ))
        return new{N,S,T}(ρ, θ)
    end
    
    function se{N}(ρ::T, θ::S) where {N,T<:AbstractVector,S<:AbstractMatrix}
        check_dim(se{N}, size(θ, 1))
        check_skewsymmetric(θ)
        check_dim(se{N}, length(ρ))
        return new{N,S,T}(ρ, θ)
    end
end

se{N}(x::AbstractVector) where {N} = se{N}(x[1:N], x[N+1:end])
se{N}(X::AbstractMatrix) where {N} = se{N}(X[1:N, end], X[1:N, 1:N])

Base.angle(alg::se{N}) where {N} = alg.θ
translation(alg::se{N}) where {N} = alg.ρ

(==)(alg1::se{N}, alg2::se{N}) where {N} = alg1.ρ == alg2.ρ && alg1.θ == alg2.θ
Base.isapprox(alg1::se{N}, alg2::se{N}) where {N} = alg1.ρ ≈ alg2.ρ && alg1.θ ≈ alg2.θ

identity(alg::se{N,T}) where {N,T<:AbstractVector} =
    se{N}(fill!(similar(alg.ρ), 0), fill!(similar(alg.θ), 0))

inv(alg::se{N,T}) where {N,T<:AbstractVector} = se{N}(-translation(alg), -angle(alg))

(+)(alg1::se{N}, alg2::se{N}) where {N} = se{N}(alg1.ρ + alg2.ρ, angle(alg1) + angle(alg2))

∧(alg::se{N,T}) where {N,T<:AbstractVector} = se{N}(∧(se{N}, translation(alg), angle(alg)))
∧(alg::se{N,T}) where {N,T<:AbstractMatrix} = alg

∨(alg::se{N,T}) where {N,T<:AbstractVector} = alg
∨(alg::se{N,T}) where {N,T<:AbstractMatrix} = se{N}(∨(se{N}, translation(alg), angle(alg)))

Base.show(io::IO, alg::se{N}) where {N} =
    print(io, "se{$N}(ρ=", translation(alg), ", θ=", angle(alg), ")")

function left_jacobian(alg::se{N}) where {N}
    ρ, θ = translation(alg), angle(alg)
    J_l = left_jacobian(so{N}(θ))
    z = fill!(similar(J_l), 0)
    return [J_l Q(ρ, θ);
              z     J_l]
end

right_jacobian(alg::se{N}) where {N} = left_jacobian(se{N}(-translation(alg), -angle(alg)))

function Q(ρ::AbstractVector, θ::AbstractVector)
    θ² = sum(abs2, θ)
    θ_angle = √θ²
    ρ_x, θ_x = skewsymmetric(ρ), skewsymmetric(θ)
    return 0.5*ρ_x +
           (θ_angle - sin(θ_angle)) / θ_angle^3 * (θ_x*ρ_x + ρ_x*θ_x + θ_x*ρ_x*θ_x) +
           -(1 - 0.5*θ_angle^2 - cos(θ_angle)) / θ_angle^4 * (θ_x^2*ρ_x + ρ_x*θ_x^2 - 3θ_x*ρ_x*θ_x) +
           -0.5((1 - 0.5*θ_angle^2 - cos(θ_angle)) / θ_angle^4 - 3(θ_angle - sin(θ_angle) - θ_angle^3/6) / θ_angle^5)*
           (θ_x*ρ_x*θ_x^2 + θ_x^2*ρ_x*θ_x)
end


# Group Interfaces

abstract type SpecialEuclideanGroup{N} <: AbstractLieGroup end

dim(::Type{<:SpecialEuclideanGroup{N}}) where {N} = N
dim(::SpecialEuclideanGroup{N}) where {N} = N

dof(::Type{<:SpecialEuclideanGroup{N}}) where {N} = sum(1:N)
dof(::SpecialEuclideanGroup{N}) where {N} = sum(1:N)

struct SE{N, T} <: SpecialEuclideanGroup{N}
    R
    t

    function SE{N}(R::AbstractMatrix{T}, t::AbstractVector{S}) where {N,T,S}
        @assert size(R, 1) == N
        @assert size(t) == (N, )
        Te = float(promote_type(T, S))
        return new{N, Te}(Te.(R), Te.(t))
    end
end

function SE{N}(A::AbstractMatrix) where {N}
    @assert size(A, 1) == N + 1
    R = A[1:N, 1:N]
    t = A[1:N, end]
    return SE{N}(R, t)
end

rotation(g::SE) = g.R
translation(g::SE) = g.t

identity(::Type{SE{N}}) where {N} = SE{N}(I(N+1))
identity(::SE{N}) where {N} = SE{N}(I(N+1))

function inv(g::SE{N}) where {N}
    R, t = rotation(g), translation(g)
    return SE{N}(R', -R'*t)
end

function (*)(::SE{M}, ::SE{N}) where {M,N}
    throw(ArgumentError("+ operation for SE{$M} and SE{$N} group is not defined."))
end

(*)(g1::SE{N}, g2::SE{N}) where {N} = SE{N}(Matrix(g1) * Matrix(g2))

function jacobian(::typeof(*), g1::SE{N}, g2::SE{N}) where {N}
    R2, t2 = rotation(g2), translation(g2)
    T2 = skewsymmetric(t2)
    z = fill!(similar(R2, N, N), 0)
    J = [R2' -R2'*T2;
           z     R2']
    return J, I(2N)
end

(==)(g1::SE{N}, g2::SE{N}) where {N} = Matrix(g1) == Matrix(g2)
Base.isapprox(g1::SE{N}, g2::SE{N}) where {N} = isapprox(Matrix(g1), Matrix(g2))

function Base.Matrix(g::SE{N}) where {N}
    R, t = rotation(g), translation(g)
    z = fill!(similar(t, 1, N), 0)
    return [R t;
            z 1]
end

function ⋉(g::SE{N}, x::AbstractVector) where {N}
    y = Matrix(g) * [x..., 1]
    return y[1:N]
end

function LinearAlgebra.adjoint(g::SE{N}) where {N}
    R, t = rotation(g), translation(g)
    T = skewsymmetric(t)
    z = fill!(similar(R), 0)
    return [R T*R;
            z   R]
end

jacobian(::typeof(inv), g::SE{N}) where {N} = -adjoint(g)

"""
Jacobian of action wrt `g`
"""
function jacobian(::typeof(⋉), g::SE{N}, x::AbstractVector) where {N}
    R = rotation(g)
    X = skewsymmetric(x)
    return [R -R*X]
end

(⊕)(g::SE{N}, alg::se{N}) where {N} = g * exp(alg)

jacobian(::typeof(⊕), g::SE{N}, alg::se{N}) where {N} =
    jacobian(*, g, exp(alg))[1], right_jacobian(alg)

# Array Interfaces

function ∧(::Type{se{N}}, ρ::AbstractVector{T}, θ::AbstractVector{T}) where {N,T}
    check_dim(se{N}, length(ρ))
    return ρ, ∧(so{N}, θ)
end

∨(::Type{se{N}}, ρ::AbstractVector, Θ::AbstractMatrix) where {N} = ρ, ∨(so{N}, Θ)

Base.show(io::IO, g::SE{N}) where {N} =
    print(io, "SE{$N}(R=", rotation(g), ", t=", translation(g), ")")


# Maps

V(θ::AbstractVector) = left_jacobian(so{3}(θ))

function Base.exp(alg::se{N,T}) where {N,T<:AbstractMatrix}
    ρ, θ = translation(alg), angle(alg)
    ρ, θ = ∨(se{N}, ρ, θ)
    return exp(se{N}(ρ, θ))
end

function Base.exp(alg::se{N,T}) where {N,T<:AbstractVector}
    ρ, θ = translation(alg), angle(alg)
    R = rotation(exp(so{N}(θ)))
    t = V(θ) * ρ
    return SE{N}(R, t)
end

function Base.log(g::SE{N}) where {N}
    R, t = rotation(g), translation(g)
    θ = angle(log(SO{N}(R)))
    ρ = inv(V(θ)) * t
    return se{N}(ρ, θ)
end
