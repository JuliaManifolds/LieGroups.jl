
@doc """
    GeneralLinearGroup{𝔽,T}

The general linear group ``$(_tex(:rm,"GL"))(n)`` is the set of all invertible matrices

```math
$(_tex(:rm,"GL"))(n) = $(_tex(:SetDef, "M ∈ 𝔽^{n×n}", "$(_tex(:rm,"det"))(M) ≠ 0", "big")),
$(_tex(:qquad)) 𝔽 ∈ $(_tex(:Set, "ℝ, ℂ")),
```
equipped with the [`MatrixMultiplicationGroupOperation`](@ref) as the group operation.

The set of invertible matrices is a Riemannian manifold, since it inherits its structure from
the embedding as an open subset of the space of matrices ``ℝ^{n×n}``.

# Constructor

    GeneralLinearGroup(n::Int; field=ℝ, kwargs...)

Generate the general linear group on ``𝔽^{n×n}``.
All keyword arguments in `kwargs...` are passed on to [`InvertibleMatrices`](@extref `Manifolds.InvertibleMatrices`).
"""
const GeneralLinearGroup{𝔽,T} = LieGroup{
    𝔽,MatrixMultiplicationGroupOperation,InvertibleMatrices{𝔽,T}
}

function GeneralLinearGroup(n::Int; field=ManifoldsBase.ℝ, kwargs...)
    Im = InvertibleMatrices(n, field; kwargs...)
    return GeneralLinearGroup{typeof(Im).parameters...}(
        Im, MatrixMultiplicationGroupOperation()
    )
end

_doc_exp_GLn = """
    exp(::GeneralLinearGroup, X)
    exp!(::GeneralLinearGroup, g, X)

Compute the Lie group exponential on the [`GeneralLinearGroup`](@ref), which is given by the
[matrix exponential](https://en.wikipedia.org/wiki/Matrix_exponential)

```math
$(_tex(:exp)) X = $(_tex(:sum))_{k=0}^{∞} $(_tex(:frac, "1", "k!"))X^k
```

see also [HilgertNeeb:2012; Example 9.2.3 (b)](@cite)
"""

@doc "$(_doc_exp_GLn)"
ManifoldsBase.exp(::GeneralLinearGroup, X)

@doc "$(_doc_exp_GLn)"
ManifoldsBase.exp!(::GeneralLinearGroup, g, X)

function ManifoldsBase.inner(
    ::LieAlgebra{𝔽,MatrixMultiplicationGroupOperation,GeneralLinearGroup}, X, Y
) where {𝔽}
    return dot(X, Y)
end

function Base.show(io::IO, G::GeneralLinearGroup{𝔽}) where {𝔽}
    n = ManifoldsBase.get_parameter(G.manifold.size)[1]
    return print(io, "GeneralLinearGroup($n; field=$(𝔽))")
end
