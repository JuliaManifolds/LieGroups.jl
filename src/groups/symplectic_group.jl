"""
    SymplecticGroup{T}

On the set of [real symplectic matrices](@extref `Manifolds.SymplecticMatrices`),
of size ``2n×2n`` for some ``n∈ℕ``

```math
$(_math(:Sp))(2n, ℝ) = $(_tex(:SetDef, "p ∈ ℝ^{2n×2n}", "p^$(_tex(:transp))J_{2n}p = J_{2n}", "big"))
```

where ``J_{2n} = $(_tex(:pmatrix, "0_n & I_n", "-I_n & 0_n"))`` denotes the [`SymplecticElement`](@extref `Manifolds.SymplecticElement`).

the corresponding LIe algebra is given by the [`HamiltonianMatrices`](@extref `Manifolds.HamiltonianMatrices`)

```math
$(_math(:sp))(2n, ℝ) = $(_tex(:SetDef, "X ∈ ℝ^{2n×2n}", "X^+ = -X", "big")),
```

where ``⋅^+`` denotes the [`symplectic_inverse`](@extref `Manifolds.symplectic_inverse`).

See [BendokatZimmermann:2021; Section 2](@cite) for more information.
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
