# Changelog

All notable Changes to the Julia package `LieGroups.jl` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - unreleased

Everything denoted by “formerly” refers to the previous name in [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/).

### Added

* `compose` and `compose!`
* `compose_diff_left`, `compose_diff_left!`, `compose_diff_right`, `compose_diff_right!` (formerly `translate_diff` with different sides)
* `compose_inv_left` and `compose_inv_right` (formerly only the `inverse_translate` but it had modes that also just did compose, the odes here are `LeftForward` and `RightBackward`)
* `conjugate` and `conjugate!`
* `exp(G::LieGroup, g, X)` and `exp!(G::LieGroup, h, g, X)` (formerly `exp_inv` and `exp_inv!`)
* `exp(G::LieGroup, ::Identity, X)` and `exp!(G::LieGroup, h, ::Identity, X)` (formerly `exp_lie` and `exp_lie!`)
* `Identity`
* `idenity_element` and `identity_element!`
* `is_identity`
* `inv` and `inv!`
* `Lie_bracket` and `Lie_bracket!` (formerly `lie_bracket`)
* `log(G::LieGroup, g, h)` and `log!(G::LieGroup, X, g, h)` (formerly `log_inv` and `log_inv!`)
* `log(G::LieGroup, ::Identity, g)` and `log!(G::LieGroup, X, ::Identity, g)` (formerly `log_lie` and `log_lie!`)
* `base_manifold` to access the manifold within a Lie group
* `LieGroup` (formerly `GroupManifold`) as well as the concrete groups
  * `AdditiveGroup` (formerly `TranslationGroup`)
* `LieGroups.AbstractGroupOperation` as well as its concrete subtypes
  * `AdditiveGroupOperation`
