# Algebra Interfaces

abstract type SpecialEuclideanAlgebra{N} <: AbstractLieAlgebra end

dim(::Type{<:SpecialEuclideanAlgebra{N}}) where {N} = N
dim(::SpecialEuclideanAlgebra{N}) where {N} = N

dof(::Type{<:SpecialEuclideanAlgebra{N}}) where {N} = sum(1:N)
dof(::SpecialEuclideanAlgebra{N}) where {N} = sum(1:N)

struct se{N,V} <: SpecialEuclideanAlgebra{N}
    ρ::V

    function se{N}(x::T) where {N,T<:AbstractVector}
        @assert check_dof(se{N}, length(x))
        return new{N,T}(x)
    end
    
    function se{N}(X::T) where {N,T<:AbstractMatrix}
        @assert size(X, 1) == N + 1
        @assert isskewsymmetric(X[1:N, 1:N])
        return new{N,T}(X)
    end
end

check_dof(::Type{se{N}}, d::Int) where {N} = d == dof(SE{N})

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

abstract type SpecialEuclideanGroup{N} <: AbstractLieGroup end

dim(::Type{<:SpecialEuclideanGroup{N}}) where {N} = N
dim(::SpecialEuclideanGroup{N}) where {N} = N

dof(::Type{<:SpecialEuclideanGroup{N}}) where {N} = sum(1:N)
dof(::SpecialEuclideanGroup{N}) where {N} = sum(1:N)

struct SE{N, T} <: SpecialEuclideanGroup{N}
    A::T

    function SE{N}(A::T) where {N,T<:AbstractMatrix}
        @assert size(A, 1) == N + 1
        return new{N, T}(A)
    end
end

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

function ⋉(g::SE{N}, x::AbstractVector) where {N}
    y = Matrix(g) * [x..., 1]
    return y[1:N]
end

function LinearAlgebra.adjoint(g::SE{N}) where {N}
    R, t = g.A[1:N, 1:N], g.A[1:N, end]
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
    R = g.A[1:N, 1:N]
    X = skewsymmetric(x)
    return [R -R*X]
end

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
    @assert check_dof(se{N}, length(alg))

    p = alg[1:N]
    Ω = ∧(so{N}, alg[N+1:end])
    z = zeros(T, 1, N)
    return [Ω p;
            z 0]
end

function ∨(::Type{se{N}}, alg::AbstractMatrix) where {N}
    d = size(alg, 1) - 1
    @assert check_dof(se{N}, sum(1:d))
    p, R = alg[1:N, N], alg[1:N, 1:N]
    return [p..., ∨(so{N}, R)...]
end

Base.show(io::IO, g::SE{N}) where {N} = print(io, "SE{$N}(A=", g.A, ")")


# Connection between groups and algebra

V(θ::AbstractVector) = left_jacobian(so{3}(θ))

function Base.exp(alg::se{N,T}) where {N,T<:AbstractMatrix}
    alg = se{N}(∨(se{N}, alg.ρ))
    return exp(alg)
end

function Base.exp(alg::se{N,T}) where {N,T<:AbstractVector}
    ρ, θ = alg.ρ[1:N], alg.ρ[N+1:end]
    R = exp(so{N}(θ)).A
    ρ = V(θ) * ρ
    z = fill!(similar(ρ, 1, N), 0)
    return SE{N}(
        [R ρ;
         z 1]
    )
end

function Base.log(g::SE{N}) where {N}
    p, R = g.A[1:N, end], g.A[1:N, 1:N]
    so = log(SO{N}(R))
    θ = so.θ
    p = inv(V(θ)) * p
    return se{N}([p..., θ...])
end
