# Changelog

All notable Changes to the Julia package `LieGroups.jl` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

* There is now a [JuliaCon proceedings paper](https://doi.org/10.21105/jcon.00195) about `LieGroups.jl`, see the citation in the Readme.md.

## [0.1.9] 2025-11-27

### Changed

* the formerly internal `LieGroupsTestSuite` module inside tests is now available as `LieGroups.Test` as soon as `Test.jl` is loaded.
* refactored the Project.toml to use a `[workspace]`
* moved the test dependencies into a separate `test/Project.toml`

## [0.1.8] 2025-10-31

### Added

* a `MetricLieGroup` as a meta type to be able to change the metric imposed on the Lie algebra and its effects e.g. on the coordinates, inner product and the exponential and logarithmic map.
* `diff_left_compose` and `diff_right_compose` for the special Euclidean group.
* `submanifold_components` now works with `Identity{<:ProductGroupOperation}`.
* `zero_vector(G, Identity(G))` returns the zero vector in default representation.

### Fixed

* fixed the contributing.md to mention runic as the code formatter.
* `diff_left_compose` correctly allocates when `g` is `Identity`.
* `init_constants!` fixed for `SpecialEuclideanMatrixTangentVector` and `SpecialEuclideanMatrixPoint`.

## [0.1.7] 2025-10-28

### Added

* `NestedPowerRepresentation` and `NestedReplacingPowerRepresentation` are now re-exported from `ManifoldsBase.jl`.
* `getindex` access to parts of a point of tangent vector represented by a matrix with `:Rotation` and `:Translation` as indices. For example, when `g` is a point on any variant of the special euclidean group, `g[G, :Translation]` will return the translation part of `g` regardless of whether `G` is the left or right semidirect product even when `g` is a matrix.

### Fixed

* Zenodo metadata.

## [0.1.6] 2025-10-10

### Fixed

* `diff_group_apply` has corrected documentation and works with `g` equal to `Identity`.

## [0.1.5] 2025-10-09

### Changed

* `LieGroups.jl` now requires `ManifoldsBase.jl` v2.0 and `Manifolds.jl` v0.11.

## [0.1.4] 2025-10-02

### Added

* mention `adjoint_matrix` in the transition documentation (#52)
* introduce `jacobian_exp` for the Jacobian of the exponential function.
* Introduce a `AbstractActionActsOnType` to distinguish, what previously was called “side”,
  i.e. whether an action acts on the left (`ActionActsOnLeft`) or right (`ActionActsOnRight`)
* Introduce `_inv` and `_inv!` functions for the inverse operation to work the same way as `_compose` and `_compose!`, respectively.
* introduce a `LeftMultiplicationGroupAction` to represent left multiplication actions on Lie groups. While this is a bit more of a technical
  name, it replaces the old `ComplexPlanarRotation`, `QuaternionRotation`, and `RotationAction`, since in the `GroupAction` this type is coupled with a group and a manifold anyway.
* further methods for the `LeftMultiplicationGroupAction` to cover the previous functionality of e.g.
  * `RotationTranslationActionOnVector`
* Add the Special Galilean Group (#60)

### Fixed

* Fixed an issue where internally the product manifold in a `SemidirectProductLieGroup` was accidentally splashed
* make `×` and `ProductLieGroup` behave the same way as `×` and `ProductManifold` do
* Fixed an issue, where `exp!(G, h, g, X)` would return a wrong result if the input `g`and the output `h`are aliased (#63).
* Fixed issues with the documentation of `diff_left_compose` and `diff_right_compose` that were inconsistent
* fixed `push_forward_tangent` and `pull_back_tangent` which on general Lie groups provided a wrong default.
* fixes an allocation bug for `apply` and `diff_apply` (#52)


## [0.1.3] 2025-08-04

### Added

* introduce `push_forward_tangent` and `pull_back_tangent` to combine the differential of left compose and its inverse to “move” from the Lie algebra to a certain tangent space and back, but also takes care of adapting the representation, for the case where the representation on the manifold is different from the one on the Lie group / Lie algebra.
* introduce a `BaseManifoldRetraction` to be able to use retractions on the underlying manifold also as a retraction on the Lie group, cf. (#43) and (#47). This feature assumes that the representation of points and tangent vectors on Lie group and the underlying manifold are the same (so it doesn't work with special Euclidean group with homogeneous coordinates).

### Changed

* Switch to using [Runic.jl](https://github.com/fredrikekre/Runic.jl) as code formatter

### Fixed

* Fixed a typo, where `within` was misspelled as `widthin` which caused errors in a few places.
* fix `default_basis` for `LieGroup` to return a `DefaultLieAlgebraOrthogonalBasis` also when providing a point type. That way `get_vector` falls back to the manifold when called with a Lie group and a point, though this is mere a historical format and the Lie algebra approach is the recommended one.
* mention `get_coordinates`, `get_vector`, `hat`, and `vee` in the transition documentation since it moved to using the `LieAlgebra` instead of the Lie group and a point.
* Fixed `RightGroupOperationAction` to be a subtype of `AbstractRightGroupActionType`
* Add `lie_bracket` for the SpecialEuclideanGroup.
* For the CircleGroup(ℝ), fixed compose StackOverflowError and a bug where result could be outside [-π,π), see (#62) for detail.

## [0.1.2] 2025-06-24

### Added

* `is_flat` for `SpecialEuclideanGroup`
* `inner` and `norm` for `LieAlgebra` to compute the inner product and norm on the Lie algebra.
* a test suite function for `identity_element`.
* New StaticArrays.jl specializations for multiple functions, including:
  * `exp` and `log` on the orthogonal and special orthogonal group in 2 and 3 dimensions.
  * `get_coordinates` and `get_vector` on the orthogonal and special orthogonal group in 2 and 3 dimensions, for `LieAlgebraOrthogonalBasis`.
* More generic implementation of non-mutating `get_vector_lie` on `AbstractProductGroupOperation` groups.

### Changed

* `identity_element` on `UnitaryGroup(1, ℍ)` now returns by default a 1x1 `Matrix` instead of a number to be consistent with higher-dimensional unitary quaternionic groups. Use `identity_element(UnitaryGroup(1, ℍ), QuaternionF64)` to get a number corresponding to the identity.

### Fixed

* `get_vector` on `SpecialEuclideanGroup` with `ArrayPartition` point type.
* `identity_element` and `zero_vector` are now all using a type as second argument and
  respect this type more thoroughly.
* fixes (#44) (accuracy of `log` on SE(2) and SE(3) for small angles).

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
* Finishes a generic implementation of a `SemidirectProductGroupOperation`
* All details will be detailed in the next release

## Old Changelog pre 0.0.3

__Two previous releases where done by Yueh-Hua Tu in 2022 before he was so kind to transfer the development to the JuliaManifolds GitHub organisation.__

All notable changes to this project will be documented in this file.

### [0.0.2]

* fix SE{3} and add jacobian
* define dof and dim for Lie algebra and jacobian of inv
* add action for SE{N}
* add se3_location example
