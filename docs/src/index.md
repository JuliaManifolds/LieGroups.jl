# LieGroups.jl

Welcome to the Documentation of `LieGroups.jl`.

```@meta
CurrentModule = LieGroups
```

```@docs
LieGroups.LieGroups
```

The implemented [Lie groups](https://en.wikipedia.org/wiki/Lie_group) use the interface for manifolds in [`ManifoldsBase.jl`](@extref ManifoldsBase :doc:`index`) together with an [interface for Lie groups](interface/group.md) and [Lie algebras](interface/algebra.md) as well as internally using the manifolds implemented in [`Manifolds.jl`](@extref Manifolds :doc:`index`).

For more general information about the history of and contributions to the package see the [About](about.md) page.

## Getting started

To install the package just type

```julia
using Pkg; Pkg.add("LieGroups")
```

Then you can directly start, for example consider the [`SpecialEuclideanGroup`](@ref)
``\mathrm{SE}(3)`` representing all orientations and places an object can take in ``ℝ^3``.
These are characterised by a ``3×3`` rotation matrix together with a point the object is at.
For example.
having such a point, we can use the Lie group logarithmic function [`log(G::SpecialEuclideanGroup, g)`](@ref)
and the Lie group exponential function [`exp(G::SpecialEuclideanGroup, X)`](@ref)
to create an orientation “half the way” from the origin pose.

The default representation is in [homogeneous coordinates]()

```@example start
using LieGroups
SE3 = SpecialEuclideanGroup(3)
g = 1/sqrt(2) .* [1.0 -1.0 0.0 0.0; 1.0 1.0 0.0 3.0*sqrt(2); 0.0 0.0 sqrt(2) 0.0; 0.0 0.0 0.0 sqrt(2)]
```

Then half that pose is

```@example start
h = exp(SE3, 0.5 .* log(SE3, g))
```

To check, just “perform that movement” twice with the group operation
[`compose`](@ref) of `h` with itself to get `g` back

```@example start
compose(SE3, h, h)
```

for more details see the [get started](tutorials/getstarted.md) tutorial.