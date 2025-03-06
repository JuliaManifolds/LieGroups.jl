"""
    SymplecticGroup{T}
"""
const SymplecticGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,SymplecticMatrices{T,𝔽}
}

function SymplecticGroup(n, field=ℝ; kwargs...)
    S = Manifolds.SymplecticMatrices(n, field; kwargs...)
    return SymplecticGroup{field,typeof(S).parameters[1]}(
        S, MatrixMultiplicationGroupOperation()
    )
end

function Base.show(io::IO, G::SymplecticGroup{𝔽}) where {𝔽}
    size = Manifolds.get_parameter(G.manifold.size)[1]
    return print(io, "SymplecticGroup($(2*size), $(𝔽))")
end
