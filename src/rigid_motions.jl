# Algebra Interfaces

abstract type SpecialEuclideanAlgebra <: AbstractLieAlgebra end

struct se{N,V} <: SpecialEuclideanAlgebra
    ρ::V

    function se{N}(x::T) where {N,T<:AbstractVector}
        d = length(x)
        @assert check_dim(se{N}, d)
        return new{N,T}(x)
    end
    
    function se{N}(X::T) where {N,T<:AbstractMatrix}
        @assert size(X, 1) == N + 1
        @assert isskewsymmetric(X[1:N, 1:N])
        return new{N,T}(X)
    end
end

check_dim(::Type{se{N}}, d::Int) where {N} = d == sum(1:N)

(==)(alg1::se{N}, alg2::se{N}) where {N} = alg1.ρ == alg2.ρ
Base.isapprox(alg1::se{N}, alg2::se{N}) where {N} = isapprox(alg1.ρ, alg2.ρ)

identity(alg::se{N,T}) where {N,T<:AbstractVector} =
    se{N}(fill!(similar(alg.ρ), 0))

inv(alg::se{N,T}) where {N,T<:AbstractVector} = se{N}(-alg.ρ)

(+)(alg1::se{N}, alg2::se{N}) where {N} = se{N}(alg1.ρ + alg2.ρ)

function Base.show(io::IO, alg::se{3})
    print(io, "se{3}(x=", alg.ρ[1], ", y=", alg.ρ[2], ", z=", alg.ρ[3])
    print(io, ", θ_x=", alg.ρ[4], ", θ_y=", alg.ρ[5], ", θ_z=", alg.ρ[6], ")")
end

function left_jacobian(alg::se{3})
    ρ, θ = alg.ρ[1:3], alg.ρ[4:end]
    J_l = left_jacobian(so{3}(θ))
    z = fill!(similar(J_l), 0)
    return [J_l Q(ρ, θ);
              z     J_l]
end

right_jacobian(alg::se{3}) = left_jacobian(se{3}(-alg.ρ))

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

abstract type SpecialEuclideanGroup <: AbstractLieGroup end

struct SE{N, T} <: SpecialEuclideanGroup
    A::T

    function SE{N}(A::T) where {N,T<:AbstractMatrix}
        @assert size(A, 1) == N + 1
        return new{N, T}(A)
    end
end

dim(::Type{SE{N}}) where {N} = N
dim(::SE{N}) where {N} = N

dof(::Type{SE{3}}) = 6
dof(::SE{3}) = 6

identity(::Type{SE{N}}) where {N} = SE{N}(I(N+1))
identity(::SE{N}) where {N} = SE{N}(I(N+1))

function inv(g::SE{N}) where {N}
    R = g.A[1:N, 1:N]
    t = g.A[1:N, N+1]
    E = SE_matrix(R', -R'*t)
    return SE{N}(E)
end

function (*)(::SE{M}, ::SE{N}) where {M,N}
    throw(ArgumentError("+ operation for SE{$M} and SE{$N} group is not defined."))
end

(*)(g1::SE{N}, g2::SE{N}) where {N} = SE{N}(g1.A * g2.A)

(==)(g1::SE{N}, g2::SE{N}) where {N} = g1.A == g2.A
Base.isapprox(g1::SE{N}, g2::SE{N}) where {N} = isapprox(g1.A, g2.A)

Base.Matrix(g::SE) = g.A

(⊕)(g::SE{N}, alg::se{N}) where {N} = g * exp(alg)

# function jacobian(::typeof(⊕), g::SE{N}, alg::se{N}) where {N}
#     return , right_jacobian(alg)
# end

# Array Interfaces

function SE_matrix(R::AbstractMatrix{T}, t::AbstractVector{S}) where {T,S}
    n = size(R, 2)
    Te = promote_type(T, S)
    return [              R t;
            zeros(Te, 1, n) 1]
end

function ∧(::Type{se{N}}, alg::AbstractVector{T}) where {N,T}
    d = length(alg)
    @assert check_dim(se{N}, d)

    p = alg[1:N]
    Ω = ∧(so{N}, alg[N+1:end])
    z = zeros(T, 1, N)
    return [Ω p;
            z 0]
end

function ∨(::Type{se{N}}, alg::AbstractMatrix) where {N}
    d = size(alg, 1)
    @assert check_dim(se{N}, sum(1:d-1))
    
    return [alg[1:N, N]..., ∨(so{N}, alg[1:N, 1:N])...]
end

Base.show(io::IO, g::SE{N}) where {N} = print(io, "SE{$N}(A=", g.A, ")")


# Connection between groups and algebra

Base.exp(alg::se{N,T}) where {N,T<:AbstractMatrix} = SE{N}(exp(alg.ρ))
Base.exp(alg::se{N,T}) where {N,T<:AbstractVector} = SE{N}(exp(∧(se{N}, alg.ρ)))
Base.log(g::SE{N}) where {N} = se{N}(∨(se{N}, log(g.A)))
