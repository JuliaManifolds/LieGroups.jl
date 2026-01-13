<div align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl/main/docs/src/assets/logo_text_readme_dark.png">
      <img alt="Manifolds.jl logo with text on the side" src="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl//main/docs/src/assets/logo_text_readme.png">
    </picture>
</div>

| **Documentation** | **Source** | **Citation** |
|:----------------- |:---------------------- |:------------ |
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliamanifolds.github.io/LieGroups.jl/stable/) | [![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl) | [![DOI](https://proceedings.juliacon.org/papers/10.21105/jcon.00195/status.svg)](https://doi.org/10.21105/jcon.00195)
 |
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/LieGroups.jl/dev/) | [![CI](https://github.com/JuliaManifolds/LieGroups.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaManifolds/LieGroups.jl/actions?query=workflow%3ACI+branch%3Amain) | [![DOI](https://zenodo.org/badge/481478376.svg)](https://doi.org/10.5281/zenodo.15343362)
| | [![codecov](https://codecov.io/gh/JuliaManifolds/LieGroups.jl/graph/badge.svg?token=32odCSyJX5)](https://codecov.io/gh/JuliaManifolds/LieGroups.jl) | |
| | [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl) | |

This package is a rework of the Lie group features of [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) in a unified way into a separate package. It especially puts more focus on the Lie group defaults and handling the corresponding Lie algebra.

## Installation

In Julia you can install this package by typing

```julia
using Pkg; Pkg.add("LieGroups")
```

in the Julia REPL. For a first start, see the [get started tutorial](https://juliamanifolds.github.io/LieGroups.jl/stable/tutorials/getstarted/).

## Contributing

Contributions are encouraged and appreciated! See [the Contributing page](https://juliamanifolds.github.io/LieGroups.jl/stable/contributing/) in the documentation for further notes, for example the code style.

## Citation

If you use `LieGroups.jl` in your work, please cite the following open access JuliaCon proceedings paper

```biblatex
@article{BergmannBaran:2026,
    Author = {Bergmann, Ronny and Baran, Mateusz},
    Doi = {10.21105/jcon.00195},
    Journal = {Proceedings of the JuliaCon Conferences},
    Volume = {8},
    Number = {79},
    Pages = {195}, Year = {2026},
    Publisher = {The Open Journal},
    Title = {Groups and smooth geometry using LieGroups.jl},
}
```

To refer to a certain version, we recommend to also cite for example

```biblatex
@software{axen_2025_17737448,
  Author = {Axen, Seth D. and Baran, Mateusz and Bergmann, Ronny and Tu, Yueh-Hua and Verdier, Olivier},
  Doi = {10.5281/zenodo.15343362},
  Publisher    = {Zenodo},
  Title        = {LieGroups.jl},
  Year         = {2025},
}
```

for the most recent version or a corresponding version specific DOI, see [the list of all versions](https://zenodo.org/search?q=parent.id%3A15343362&f=allversions%3Atrue&l=list&p=1&s=10&sort=version).
Note that both citations are in [BibLaTeX](https://ctan.org/pkg/biblatex) format.


> [!NOTE]
> This is a rework of the features from [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/).
> See [transition from Manifolds.jl](https://juliamanifolds.github.io/LieGroups.jl/stable/tutorials/transition/) for a comprehensive list how to update your code.
> This especially also includes a few different choices in default behaviour that
is different from the [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) one. For purely manifold-based operations, any Lie group still is “build upon” a Riemannian manifold.