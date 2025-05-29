# Changelog

All notable Changes to the Julia package `LieGroups.jl` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] unreleased

### Added

* `is_flat` for `SpecialEuclideanGroup`
* `inner` and `norm` for `LieAlgebra` to compute the inner product and norm on the Lie algebra.

### Fixed

* `get_vector` on `SpecialEuclideanGroup` with `ArrayPartition` point type.

## [0.1.1] 2025-05-05

### Added

* `identity_element` on `TranslationGroup` supports now `StaticArrays.jl` types.
* introduce `get_vector` in legacy form to work on Lie groups, but they pass on to their Lie algebra.
* adapt to the new `default_basis` from ManifoldsBase.jl 1.1.

### Changed

* the tutorials are now rendered with `quarto` using the [`QuartoNotebookRunner.jl`](https://github.com/PumasAI/QuartoNotebookRunner.jl) and are hence purely julia based.

### Fixed

* `identity_element` on `TranslationGroup` no longer accepts a number as a second argument (it accepts number type instead).

## [0.1.0] 2025-04-22

Everything denoted by “formerly” refers to the previous name in [`Manifolds.jl`](https://juliamanifolds.github.io/Manifolds.jl/stable/).
Several structs have been changed from the pre-release, so these are breaking.

### Added

* `LieAlgebra`
* `LieGroup` (formerly `GroupManifold`) as well as the concrete groups
  * `TranslationGroup`
  * `SpecialEuclideanGroup` (formerly `SpecialEuclidean`) including
    * `SpecialEuclideanMatrixPoint` and `SpecialEuclideanMatrixTangentVector` when representing the points as affine (abstract) matrices
    * `SpecialEuclideanProductPoint` and `SpecialEuclideanProductTangentVector` when representing them in a product structure, that is as an `ArrayPartition` from [`RecursiveArrayTools`](https://github.com/SciML/RecursiveArrayTools.jl).
    * neither of those types is necessary, besides for conversion between both. The product representation differs for the left and right semidirect product, while the affine matrix variant does not.
  * `SpecialOrthogonalGroup` (formerly `SpecialOrthogonal`)
  * `SpecialUnitaryGroup` (formerly `SpecialUnitary`)
  * `OrthogonalGroup` (formerly `Orthogonal`)
  * `UnitaryGroup` (formerly `Unitary`) also for quaternions.
  * `GeneralLinearGroup` (formerly `GeneralLinear`)
  * `HeisenbergGroup`
  * `LeftSemidirectProductLieGroup` (formerly `SemidirectProductGroup`)
  * `⋉` (alias for `LeftSemidirectProductGroupOperation` when a `default_left_action(G,H)` is defined for the two groups)
  * `PowerLieGroup` (formerly `PowerGroup`)
  * `PowerGroupOperation` to internally avoid ambiguities. Since the constructor always expects a Lie group, this is only necessary internally
  * `ProductLieGroup` (formerly `ProductGroup`)
  * `RightSemidirectProductLieGroup`
  * `SpecialLinearGroup` (formerly `SpecialLinear`)
  * `SymplecticGroup`
  * `CircleGroup` now with even three representations: Real line (mod 2π), Complex and plane circle
  * `⋊` (alias for `RightSemidirectProductGroupOperation` when a `default_right_action(G,H)` is defined for the two groups)
  * a `ValidationLieGroup` verifying input and output of all interface functions, similar to the [`ValidationManifold`](https://juliamanifolds.github.io/ManifoldsBase.jl/stable/manifolds/#A-manifold-for-validation) which can also be used internally.
* `AbstractGroupOperation` as well as its concrete subtypes
  * `AdditionGroupOperation` (formerly `AdditionOperation`)
  * `MatrixMultiplicationGroupOperation` (formerly `MultiplicationOperation`)
  * `PowerGroupOperation` (formerly the Lie group was stored inside a power manifold)
  * `ProductGroupOperation` (formerly the Lie groups were stored inside a product manifold)
  * `LeftSemidirectProductGroupOperation` (this was formerly only implicitly stored in the `SemidirectProductGroup`)
  * `RightSemidirectProductGroupOperation`
* `AbstractGroupActionType` with its 2 specific (new) abstract subtypes
  * `AbstractLeftGroupActionType`
  * `AbstractRightGroupActionType`
* For the group operation actions there are now
  * `LeftGroupOperationAction` (formerly `LeftForwardAction`)
  * `RightGroupOperationAction` (formerly `RightBackwardAction`)
  * `InverseLeftGroupOperationAction` (formerly `RightForwardAction`)
  * `InverseRightGroupOperationAction` (formerly `LeftBackwardAction`)
* `DefaultLieAlgebraOrthogonalBasis` (replaces `VeeOrthogonalBasis`, which is still available in `ManifoldsBase.jl`)
* `AbstractLieGroupPoint` and `AbstractLieAlgebraTangentVector` as abstract types to introduce point and Lie algebra tangent vector representations
* `Identity`
* `apply`and `apply!`
* `base_manifold` to access the manifold within a Lie group
* `compose` and `compose!`
* `conjugate` and `conjugate!`
* `diff_apply`, `diff_apply!`, `diff_group_apply`, and `diff_group_apply!` (formerly `apply_diff_[group][!]`)
* `diff_conjugate` and `diff_conjugate!`
* `diff_left_compose`, `diff_left_compose!`, `diff_right_compose`, `diff_right_compose!` (formerly `translate_diff` with different sides)
* `exp(G::LieGroup, g, X)` and `exp!(G::LieGroup, h, g, X)` (formerly `exp_inv` and `exp_inv!`)
* `exp(G::LieGroup, X)` and `exp!(G::LieGroup, h, X)` (formerly `exp_lie` and `exp_lie!`)
* `hat` and `hat!`, with slightly different signatures, since the base point is omitted.
* `identity_element` and `identity_element!`
* `inv` and `inv!` (`inv(::AbstractGroupAction)` was formerly `switch_direction`)
* `inv_left_compose`, `inv_left_compose!` and `inv_right_compose`, `inv_right_compose!` (these functions correspond to `inverse_translate` with corresponding direction and side)
* `is_identity`
* `lie_bracket` and `lie_bracket!`
* `jacobian_conjugate` (formerly `adjoint_matrix`, which is now a special case of this)
* `log(G::LieGroup, g, h)` and `log!(G::LieGroup, X, g, h)` (formerly `log_inv` and `log_inv!`)
* `log(G::LieGroup, ::Identity, g)` and `log!(G::LieGroup, X, ::Identity, g)` (formerly `log_lie` and `log_lie!`)
* `switch` (formerly `switch_side`)
* `vee` and `vee!`, with slightly different signatures, since the base point is omitted.

Compared to `Manifolds.jl`
* all `translate` functions are not implemented here, since you can just use `compose`. The differentials are implemented as listed above with respect to both left and right argument of compose
* all `inverse_apply` functions are not implemented here, since it is recommended to use `apply(inv(A), g, p)` as a replacement.

## [0.0.3] – 2025-02-19


### Added

* Finishes most of the work on the interface for the `LieGroup` type and the new `LieAlgebra` type.
* Finishes a generic implementation of a `SemiDirectProductGroupOperation`
* All details will be detailed in the next release

## Old Changelog pre 0.0.3

__Two previous releases where done by Yueh-Hua Tu in 2022 before he was so kind to transfer the development to the JuliaManifolds GitHub organisation.__

All notable changes to this project will be documented in this file.

### [0.0.2]

* fix SE{3} and add jacobian
* define dof and dim for Lie algebra and jacobian of inv
* add action for SE{N}
* add se3_location example
