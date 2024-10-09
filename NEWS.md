# Changelog

All notable Changes to the Julia package `LieGroups.jl` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - unreleased

Everything denoted by “formerly” refers to the previous name in [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/).

### Added

* `LieAlgebra`
* `LieGroup` (formerly `GroupManifold`) as well as the concrete groups
  * `AdditiveGroup` (formerly `TranslationGroup`)
* `AbstractGroupOperation` as well as its concrete subtypes
  * `AdditiveGroupOperation`
* `AbstractGroupActionType` with its 2 specific (new) abstract
  * `AbstractLeftGroupActionType`
  * `AbstractRightGroupActionType`
* For the group operation actions there are now
  * `LeftGroupOperation` (formerly `LeftForwardAction`)
  * `RightGroupOperation` (formerly `RightBackwardAction`)
  * `InverseLeftGroupOperation` (formerly `RightForwardAction`)
  * `InverseRightGroupOperation` (formerly `LeftBackwardAction`)
* `LieGroups.AbstractGroupOperation` as well as its concrete subtypes
  * `AdditiveGroupOperation`
* `apply`and `apply!`
* `base_manifold` to access the manifold within a Lie group
* `compose` and `compose!`
* `conjugate` and `conjugate!`
* `diff_apply`, `diff_apply!`, `diff_group_apply`, and `diff_group_apply!` (formerly `apply_diff_[group][!]`)
* `diff_left_compose`, `diff_left_compose!`, `diff_right_compose`, `diff_right_compose!` (formerly `translate_diff` with different sides)
* `exp(G::LieGroup, g, X)` and `exp!(G::LieGroup, h, g, X)` (formerly `exp_inv` and `exp_inv!`)
* `exp(G::LieGroup, ::Identity, X)` and `exp!(G::LieGroup, h, ::Identity, X)` (formerly `exp_lie` and `exp_lie!`)
* `Identity`
* `idenity_element` and `identity_element!`
* `inv` and `inv!` (`inv(::AbstractGroupAction)` was formerly `switch_direction`)
* `inv_left_compose`, `inv_left_compose!` and `inv_right_compose`, `inv_right_compose!` (these functions correspond to `inverse_translate` with corresponding direction and side)
* `is_identity`
* `Lie_bracket` and `Lie_bracket!` (formerly `lie_bracket`)
* `log(G::LieGroup, g, h)` and `log!(G::LieGroup, X, g, h)` (formerly `log_inv` and `log_inv!`)
* `log(G::LieGroup, ::Identity, g)` and `log!(G::LieGroup, X, ::Identity, g)` (formerly `log_lie` and `log_lie!`)
* `switch` (formerly `switch_side`)

Compared to `Manifolds.jl`
* all `translate` functions are not implemented here, since you can just use `compose`. The differentials are implemented as listed above with respect to both left and right argument of compose
* all `inverse_apply` functions are not implemented here, since it is recommended to use `apply(inv(A), g, p)` as a replacement.

## Old Changelog pre 0.1.0

__Two previous releases where done by Yueh-Hua Tu in 2022 before he was so kind to transfer the development to the JuliaManifolds GitHub organisation.__

All notable changes to this project will be documented in this file.

### [0.0.2]

* fix SE{3} and add jacobian
* define dof and dim for Lie algebra and jacobian of inv
* add action for SE{N}
* add se3_location example