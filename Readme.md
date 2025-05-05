<div align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl/main/docs/src/assets/logo_text_readme_dark.png">
      <img alt="Manifolds.jl logo with text on the side" src="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl//main/docs/src/assets/logo_text_readme.png">
    </picture>
</div>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/LieGroups.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/LieGroups.jl/dev/)
[![DOI](https://zenodo.org/badge/481478376.svg)](https://doi.org/10.5281/zenodo.15343362)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![CI](https://github.com/JuliaManifolds/LieGroups.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaManifolds/LieGroups.jl/actions?query=workflow%3ACI+branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaManifolds/LieGroups.jl/graph/badge.svg?token=32odCSyJX5)](https://codecov.io/gh/JuliaManifolds/LieGroups.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package is a rework of the Lie group features of [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) in a unified way into a separate package. It especially puts more focus on the Lie group defaults and handling the corresponding Lie algebra.

## Installation

In Julia you can install this package by typing

```julia
using Pkg; Pkg.add("LieGroups")
```

in the Julia REPL. For a first start, see the [get started tutorial](https://juliamanifolds.github.io/LieGroups.jl/stable/tutorials/getstarted/).

> [!NOTE]
> Since this is a rework of the features from [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/), both `LieGroups.jl` and `Manifolds.jl` 0.10 export a few types of same name, for example `Identity`.
While `LieGroups.jl` depends on `Manifolds.jl`, it is not recommended to load both into the same name space, that is, doing `using Manifolds.jl, LieGroups.jl`, since then these conflicts might lead to unforeseen errors, where you would need to specify the name space to resolve this ambiguity.
> See [transition from Manifolds.jl](https://juliamanifolds.github.io/LieGroups.jl/stable/tutorials/transition/) for a comprehensive list.
> This especially also includes a few different choices in default behaviour that
is different from the [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) one. For purely manifold-based operations, any Lie group still is “build upon” a Riemannian manifold.
