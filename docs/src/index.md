```@raw html
---
layout: home

hero:
  name: LieGroups.jl
  text: Lie groups and Lie algebras in Julia
  tagline: Work with Lie groups, their Lie algebras, and group actions – on top of the manifolds interface.
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
  - icon: 🎬
    title: Lie groups
    details: A library of Lie groups, from the classical matrix groups to power, product, and semidirect product groups, all with a common interface.
    link: /groups/index.html
  - icon: 🧮
    title: Lie algebras
    details: Every Lie group comes with its Lie algebra, including the Lie bracket, `hat` and `vee`, and the adjoint representation.
    link: /interface/algebra/index.html
  - icon: 🎯
    title: Group actions
    details: Left and right group actions on manifolds are available as their own interface, so groups can act on the data your problem lives on.
    link: /interface/actions/index.html
  - icon:
        light: /logo-manifoldsbase.png
        dark: /logo-manifoldsbase-dark.png
        alt: ManifoldsBase.jl
        wrap: true
    title: ManifoldsBase.jl
    details: "Lie groups here are built on the manifolds interface of [ManifoldsBase.jl](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/), so every Lie group is also a manifold and all functions from that interface are available."
  - icon:
        light: /logo-manifolds.png
        dark: /logo-manifolds-dark.png
        alt: Manifolds.jl
        wrap: true
    title: Manifolds.jl
    details: "The manifolds a Lie group is based on are taken from [Manifolds.jl](https://juliamanifolds.github.io/Manifolds.jl/stable/), which provides a comprehensive library of Riemannian manifolds."
  - icon:
        src: /logo-manopt.png
        alt: Manopt.jl
        wrap: true
    title: Manopt.jl
    details: "Since every Lie group is a manifold, the optimization algorithms in [Manopt.jl](https://manoptjl.org/stable/) can be used to solve optimization problems on Lie groups."
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