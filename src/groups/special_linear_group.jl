@doc """
    SpecialLinear{𝔽,T}

The special linear group ``$(_math(:SL))(n,𝔽)`` is the group of all invertible matrices
with unit determinant in ``𝔽^{n×n}`` and the [`MatrixMultiplicationGroupOperation`](@ref) as group operation.

The Lie algebra ``$(_math(:sl))(n, 𝔽) = T_e $(_math(:SL))(n,𝔽)`` is the set of all matrices in
``𝔽^{n×n}`` with trace of zero. By default, tangent vectors ``X_p ∈ T_p $(_math(:SL))(n,𝔽)``
for ``p ∈ $(_math(:SL))(n,𝔽)`` are represented with their corresponding Lie algebra vector
``X_e = p^{-1}X_p ∈ 𝔰𝔩(n, 𝔽)``.

# Constructor

    GeneralLinearGroup(n::Int, field=ℝ; kwargs...)

Generate the general linear group  group on ``𝔽^{n×n}``.
All keyword arguments in `kwargs...` are passed on to [`InvertibleMatrices`](@extref `Manifolds.GeneralUnitaryMatrices`)`{T, 𝔽,`[`DeterminantOneMatrices`](@extref `Manifolds.DeterminantOneMatrices`)`}`.
"""
const SpecialLinearGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,DeterminantOneMatrices{𝔽,T}
}

function SpecialLinearGroup(n::Int, field=ManifoldsBase.ℝ; kwargs...)
    M = Manifolds.DeterminantOneMatrices(n, field; kwargs...)
    return SpecialLinearGroup{typeof(M).parameters...}(
        M, MatrixMultiplicationGroupOperation()
    )
end

# TODO: document hat/vee with the corresponding formulae

function get_coordinates_lie!(
    ::LieAlgebra{ℝ,MatrixMultiplicationGroupOperation,<:SpecialLinearGroup},
    c,
    X,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
)
    c .= X[1:(end - 1)]
    return c
end

function get_vector_lie!(
    𝔤::LieAlgebra{ℝ,MatrixMultiplicationGroupOperation,<:SpecialLinearGroup},
    X,
    c,
    ::DefaultLieAlgebraOrthogonalBasis{ℝ},
)
    X[1:(end - 1)] .= c
    X[end] = 0
    X[end] = -tr(X)
    return X
end

function Base.show(
    io::IO, ::SpecialLinearGroup{𝔽,ManifoldsBase.TypeParameter{Tuple{n}}}
) where {𝔽,n}
    return print(io, "SpecialLinearGroup($n, $(𝔽))")
end
function Base.show(io::IO, G::SpecialLinearGroup{𝔽,Tuple{Int}}) where {𝔽}
    M = base_manifold(G)
    n = ManifoldsBase.get_parameter(M.size)[1]
    return print(io, "SpecialLinearGroup($n, $(𝔽); parameter=:field)")
end
