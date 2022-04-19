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

check_dim(::Type{se{N}}, d::Int) where {N} = d == N*(N-1)/2 + N

(==)(alg1::se{N}, alg2::se{N}) where {N} = alg1.ρ == alg2.ρ
Base.isapprox(alg1::se{N}, alg2::se{N}) where {N} = isapprox(alg1.ρ, alg2.ρ)

identity(alg::se{N,T}) where {N,T<:AbstractVector} =
    se{N}(fill!(similar(alg.ρ), 0))

inv(alg::se{N,T}) where {N,T<:AbstractVector} = se{N}(-alg.ρ)

(+)(alg1::se{N}, alg2::se{N}) where {N} = se{N}(alg1.ρ + alg2.ρ)


# Group Interfaces

abstract type SpecialEuclideanGroup <: AbstractLieGroup end

struct SE{N, T} <: SpecialEuclideanGroup
    A::T

    function SE{N}(A::T) where {N,T<:AbstractMatrix}
        @assert size(A, 1) == N + 1
        return new{N, T}(A)
    end
end

identity(::Type{SE{N}}) where {N} = SE{N}(I(N+1))
identity(::SE{N}) where {N} = SE{N}(I(N+1))

function inv(g::SE{N}) where {N}
    T = eltype(g.A)
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

function ⋉(g::SE{N,T}, x::T) where {N,T<:AbstractVector}
    y = Matrix(g) * [x..., 1]
    return y[1:N]
end


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
    E = SE_matrix(Ω, p)
    return E
end

function ∨(::Type{se{N}}, alg::AbstractMatrix) where {N}
    d = size(alg, 1)
    @assert check_dim(se{N}, d)
    
    return [alg[1:N, N]..., ∨(so{N}, alg)...]
end
