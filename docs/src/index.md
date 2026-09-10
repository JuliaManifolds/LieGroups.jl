```@raw html
---
layout: home

hero:
  name: LieGroups.jl
  text: Lie groups and Lie algebras in Julia
  tagline: Work with Lie groups, their Lie algebras, and group actions
  actions:
    - theme: brand
      text: Get started
      link: tutorials/getstarted/index.html
    - theme: alt
      text: List of Lie groups
      link: /groups/index.html
    - theme: alt
      text: Lie group interface
      link: /interface/group/index.html
  image:
    src: /logo.png            # primary image (light themes)
    dark: /logo-dark.png      # variant for dark themes
    alt: LieGroups.jl         # accessibility text

features:
  - icon: 🪶
    title: Lightweight Interface
    details: This package provides a lightweight interface to define Lie groups based on the `LieGroup` combining a `AbstractManifold` from `ManifoldsBase.jl` with a `GroupOperation`. Especially the `exp` and `log` maps are here the Lie group ones.
    link: interface/group.html
  - icon: 🧮
    title: Lie Algebra and group actions
    details: The `LieAlgebra` is generically provided including an . Furthermore `GroupActions` allow Lie groups to act on other manifolds.
    link: interface/algebra.html
  - icon: ⚡️
    title: Efficient
    details: When possible, functions are available working in-place, like `exp!` or `log!` to reduce memory allocations. The neutral element, the `Identity{<:GroupOperation}` provides an allocation free implementation that can “materialise” into the correct actual point on the Lie group when necessary.
    link: /interface/operations.html
  - icon:
        light: /logo.png
        dark: /logo-dark.png
        alt: LieGroups.jl
        wrap: true
    title: Library of Lie groups
    details: This package provides a library of Lie groups. On the one hand there are abstract product- and power- as well as semidirect product Lie groups. On the other hand a library of concrete Lie groups is available as well.
    link: /groups/index.html
  - icon: 📚
    title: Well-documented and -tested
    details: All Lie groups are documented – both their theoretical foundation and all numerical functionality. The theoretical background also refers to further literature. A test suite provides a comprehensive verification for any Lie group, existing and newly written, as well.
    link: test_suite.html
  - icon:
        light: /logo-manifolds.png
        dark: /logo-manifolds-dark.png
        alt: Manifolds.jl
        wrap: true
    title: Manifolds.jl
    details: "The manifolds from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/) build the foundation of the Lie groups implemented here. At the same time every `LieGroup` is also a manifold (with a different connection), so they can for example be used with [Manopt.jl](https://manoptjl.org/stable/).
    "
---
```

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