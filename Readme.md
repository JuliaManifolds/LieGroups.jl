<div align="center">
    <picture>
        <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl/main/docs/src/assets/logo_text_readme_dark.png">
      <img alt="Manifolds.jl logo with text on the side" src="https://raw.githubusercontent.com/JuliaManifolds/LieGroups.jl//main/docs/src/assets/logo_text_readme.png">
    </picture>
</div>

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliamanifolds.github.io/Manifolds.jl/latest/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle) [![CI](https://github.com/JuliaManifolds/LieGroups.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaManifolds/LieGroups.jl/actions?query=workflow%3ACI+branch%3Amain) [![codecov.io](http://codecov.io/github/JuliaManifolds/LieGroups.jl/coverage.svg?branch=main)](https://codecov.io/gh/JuliaManifolds/LieGroups.jl/)

This is a package to rework the Lie group features of [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) in a unified way into a separate package.

This especially also includes a few different choices in default behabiour that
is different from the [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/) one. For purely manifold-based operations, any Lie group still is “build upon” a Riemannian manifold.

See [#5](https://github.com/JuliaManifolds/LieGroups.jl/issues/5) for an overview of fatures that are plannned.
Feel free to add to this list with opening further issues.
