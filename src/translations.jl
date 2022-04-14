abstract type AbstractTranslationGroup end

struct TranslationGroup{N, T} <: AbstractTranslationGroup
    δ::T
end

function TranslationGroup{N}(δ::T) where {N,T<:AbstractVector}
    return TranslationGroup{N, T}(δ)
end

identity(::Type{TranslationGroup{N}}) where {N} = TranslationGroup{N}(zeros(N))
identity(g::TranslationGroup{N}) where {N} = TranslationGroup{N}(fill!(similar(g.δ), 0))

function (+)(::TranslationGroup{M}, ::TranslationGroup{N}) where {M,N}
    throw(ArgumentError("+ operation for TranslationGroup{$M} and TranslationGroup{$N} group is not defined."))
end

(+)(g1::TranslationGroup{N}, g2::TranslationGroup{N}) where {N} = TranslationGroup{N}(g1.δ + g2.δ)

function (-)(::TranslationGroup{M}, ::TranslationGroup{N}) where {M,N}
    throw(ArgumentError("- operation for TranslationGroup{$M} and TranslationGroup{$N} group is not defined."))
end

(-)(g1::TranslationGroup{N}, g2::TranslationGroup{N}) where {N} = TranslationGroup{N}(g1.δ - g2.δ)

# \rtimes
function ⋊(x::T, g::TranslationGroup{N,T}) where {N,T}
    return x + g.δ
end
