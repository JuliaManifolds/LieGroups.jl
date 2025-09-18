# Notation on Lie groups

In this package,the notation introduced in [Manifolds.jl Notation](https://juliamanifolds.github.io/Manifolds.jl/latest/misc/notation.html) is used with the following additional parts.

| Symbol | Description | Also used | Comment |
|:----:|:--------------- |:----:|:--- |
| ``α`` | a general group action, when it is not specified whether it is a left (``α=σ``) or right (``α=τ``) action. | | |
| ``∘`` | a group operation | | |
| ``c_g:\mathcal G → \mathcal G`` | the conjugation map (with `g`) | | |
| ``Df(p)[X]`` | the differential of a map `f` at point `p` in direction `X` | | |
| ``\mathrm{d}f`` | the differential of a map `f` as a function on the Lie group and the differential on the Lie algebra. | | see also note below |
| ``\mathrm{D}_af`` | the differential of a map `f`. An index is used to indicate a certain parameter. If ``f`` is defined on or maps into the Lie group, this differential indicates the one with respect to tangent spaces  | | see also note below|
| ``\mathrm{e}`` | identity element of a group | | |
| ``\exp_{\mathcal G}(X)`` | The Lie group exponential function | | |
| ``\exp_g(X)`` | The Lie group exponential map (w.r.t. a Cartan Schouten connection) | | |
| ``g, h, k`` | elements on a (Lie) group. Sometimes called points. | ``g_1, g_2, ...`` | |
| ``\mathfrak g`` | a Lie algebra | | |
| ``\mathcal{G}`` | a Lie group | | |
| ``\operatorname{J}_f(p)`` | the Jacobian of a map `f` at point `p` | | sometimes left Jacobian, see note below. |
| ``λ_g: \mathcal G → \mathcal G`` | the left group operation map ``λ_g(h) = g∘h`` | | |
| ``\log_{\mathcal G}(g)`` | The Lie group logarithmic function | | |
| ``\log_g(h)`` | The Lie group logarithmic map (w.r.t. a Cartan Schouten connection) | | |
| ``α: \mathcal M → \mathcal G → \mathcal M`` | a (general) group action | | |
| ``ρ_g: \mathcal G → \mathcal G`` | the right group operation map ``ρ_g(h) = h∘g`` | | |
| ``σ: \mathcal G × \mathcal M → \mathcal M`` | a left group action | | ``σ_g(p)`` to emphasize a fixed group element |
| ``τ: \mathcal G × \mathcal M → \mathcal M`` | a right group action | ``σ_\mathrm{R}`` | ``τ_g(p)`` to emphasize a fixed group element |

## About differentials and Jacobians

For a function defined on a [manifold](@exref `ManifoldsBase.AbstractManifold`) ``f:\mathcal M → \mathcal N``, the differential at a point ``p ∈ \mathcal M`` is a map between the tangent spaces

```
Df(p) : T_p\mathcal M → T_{f(p)}\mathcal N.
```

For the case where ``\mathcal M = \mathcal N = \mathcal G`` is a [`AbstractLieGroup`](@ref), the differential can be expressed in terms of the Lie algebra ``\mathfrak g`` as

```
\mathrm{d}f(g) : \mathfrak g → \mathfrak g,
```

one alternate way to define this differential on the Lie algebra is to consider the differential ``Dg(p)`` of ``g(q) = f(q⋅p)⋅f(p)^{-1}``.

where we use a different notation on purpose. This second notation is the default throughout `LieGroups.jl`.

The Jacobian ``\operatorname{J}_f(p)`` of ``f`` at ``p`` is the matrix representation of the differential with respect to a basis of each of the tangent spaces.
For the default representation ``\mathrm{d}f(g)`` we have to choose a basis of the [Lie algebra](@ref LieAlgebras) ``\mathfrak g``.
Throughout `LieGroups.jl` this is the [`DefaultLieAlgebraOrthogonalBasis`](@ref).

## About left and right Jacobians

This jacobian using ``\mathrm{d}f(g)`` and the [`DefaultLieAlgebraOrthogonalBasis`](@ref) is sometimes called the left Jacobian.
The other choice using on two tangent space bases of ``T_p\mathcal G`` and ``T_{f(p)}\mathcal G``, respectively, is sometimes called the right Jacobian.
The default throughout `LieGroups.jl` is the left Jacobian.
