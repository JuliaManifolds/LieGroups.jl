abstract type AbstractRotationGroup end

struct IdentityRotationGroup <: AbstractRotationGroup end

struct SO{N,T<:Real,P<:AbstractRotationGroup} <: AbstractRotationGroup
    pred::P  # predecessor
    θ::T
    axis::Int
end

SO{2}(θ::T, axis::Int) where {T<:Real} = SO{2,T,IdentityRotationGroup}(IdentityRotationGroup(), θ, 0)
SO{2}(θ::T) where {T<:Real} = SO{2}(θ, 0)
SO{3}(θ::T, axis::Int) where {T<:Real} = SO{3,T,IdentityRotationGroup}(IdentityRotationGroup(), θ, axis)

identity(::IdentityRotationGroup) = IdentityRotationGroup()
identity(::SO) = IdentityRotationGroup()

inv(::IdentityRotationGroup) = IdentityRotationGroup()
inv(g::SO{N}) where {N} = inv(g.pred) * SO{N}(-g.θ, g.axis)

(*)(::IdentityRotationGroup, g::SO{N}) where {N} = g
(*)(g::SO{N}, ::IdentityRotationGroup) where {N} = g

function (*)(::SO{M}, ::SO{N}) where {M,N}
    throw(ArgumentError("* operation for SO{$M} and SO{$N} group is not defined."))
end

function (*)(g1::SO{N}, g2::SO{N}) where {N}
    if g1.axis == g2.axis
        θ = g1.θ + g2.θ
        if θ == 0
            return IdentityRotationGroup()
        else
            return SO{N,typeof(θ),typeof(g1.pred)}(g1.pred, θ, g1.axis)
        end
    else
        g = (g2.pred isa IdentityRotationGroup) ? g1 : g1 * g2.pred
        return SO{N,typeof(g2.θ),typeof(g)}(g, g2.θ, g2.axis)
    end
end

(==)(g1::SO{N}, g2::SO{N}) where {N} =
    g1.pred == g2.pred && g1.axis == g2.axis && g1.θ == g2.θ

Base.Matrix(g::SO{2}) = [cos(g.θ) -sin(g.θ); sin(g.θ) cos(g.θ)]

Base.Matrix(g::SO{3}) = Matrix(g.pred, g)
Base.Matrix(g1::SO{3}, g2::SO{3}) = Matrix(IdentityRotationGroup(), g2) * Matrix(g1)

function Base.Matrix(::IdentityRotationGroup, g::SO{3})
    R = zeros(Float64, 3, 3)
    for i in 1:3
        if i == g.axis
            R[i, i] = 1
        else
            R[i, i] = cos(g.θ)
        end
    end

    ax = to0base(g.axis)
    i, j = (ax+1) % 3, (ax+2) % 3
    i, j = to1base(i), to1base(j)

    R[i, j] = -sin(g.θ)
    R[j, i] =  sin(g.θ)

    return R
end

to0base(x) = x - 1
to1base(x) = x + 1

function Base.show(io::IO, g::SO{N}) where N
    if g.pred isa IdentityRotationGroup
        print(io, "SO{$N}(θ=", g.θ, ", axis=", g.axis, ")")
    else
        print(io, g.pred)
        print(io, " ∘ SO{$N}(θ=", g.θ, ", axis=", g.axis, ")")
    end
end
