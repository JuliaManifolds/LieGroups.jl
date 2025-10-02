# Transition from `GroupManifolds` in `Manifolds.jl`

One predecessor of `LieGroups.jl` are the [`GroupManifold`](@extref `Manifolds.GroupManifold`)s in `Manifolds.jl`.
While this package provides the same features, one reason for a new package is,
that a “restart” offers the opportunity to put the main focus for the functions in this package
really on Lie groups.

This tutorial provides an overview of the necessary changes to your code if you based it on the predecessor.

## Table of function names and its successors

The following table lists all functions related to `GroupManifolds` and their new names
or replacements here in `LieGroups.jl`. In this code `G` always refers to the `GroupManifold`
in the first column and the [`LieGroup`](@ref) in the second.
Lie group elements (points) are always `g,h`,
Lie algebra elements (vectors) always `X, Y`.

New functions and types in this package are only mentioned, if they are worth a comment and if something changed.

The list lists first types, then functions. Within both blocks, the order is alphabetical

| `Manifolds.jl` | `LieGroups.jl` | Comment |
|:---------- |:---------- |:-------------- |
| `AdditionOperation` | [`AdditionGroupOperation`](@ref) | |
| `ColumnwiseMultiplicationAction` | [`ColumnwiseGroupAction`](@ref)`(`[`LeftMultiplicationGroupAction`](@ref)`)` | within a [`GroupAction`](@ref)`(action, group, manifold)` |
| `ColumnwiseSpecialEuclideanAction` |  [`ColumnwiseGroupAction`](@ref)`(`[`LeftMultiplicationGroupAction`](@ref)`)` | within a [`GroupAction`](@ref)`(action, group, manifold)`, where the `group`is a [`SpecialEuclideanGroup`](@ref) |
| `ComplexPlanarRotation` | [`LeftMultiplicationGroupAction`](@ref) | a slightly more general type for all actions that are implemented by (matrix) multiplication |
| `GroupActionSide` | [`AbstractActionActsOnType`](@ref) | Switching to a new, hopefully more descriptive naming. |
| `LeftBackwardAction` | [`AbstractRightGroupActionType`](@ref) and [`ActionActsOnRight`](@ref) | This tuple form has been discontinued. |
| `LeftForwardAction` | [`AbstractLeftGroupActionType`](@ref) and [`ActionActsOnLeft`](@ref) | This tuple form has been discontinued. |
| `LeftSide` | [`ActionActsOnLeft`](@ref) | |
| | [`LieAlgebra`](@ref)`(G)` | new alias to emphasize its manifold- and vector structure as well as for a few dispatch methods. |
| `GroupManifold(M, op)` | [`LieGroup`](@ref)`(M, op)` | |
| `PowerGroup(M)` | [`PowerLieGroup`](@ref)`(G,n)` | The main change is, that the constructor now requires a Lie group to build the power Lie group; This also allows for `G^n`. The other two former constructors for nested and nested-replacing are no longer necessary. `PowerLieGroup` behaves exactly the same as [`PowerManifold`](@extref `ManifoldsBase.PowerManifold`). |
| `ProductGroup(M)` | [`ProductLieGroup`](@ref)`(G, H)` | The main change is, that the constructor now requires two Lie groups to build their product. This also allows for the short hand `G×H` to generate this product. |
| `QuaternionRotation` | [`LeftMultiplicationGroupAction`](@ref) | a slightly more general type for all actions that are implemented by (matrix) multiplication |
| `RightBackwardAction` | [`AbstractRightGroupActionType`](@ref) and [`ActionActsOnRight`](@ref) | This tuple form has been discontinued. |
| `RightForwardAction` | [`AbstractRightGroupActionType`](@ref) and [`ActionActsOnLeft`](@ref) | This tuple form has been discontinued. |
| `RightSide` | [`ActionActsOnRight`](@ref) | |
| `RotationAction` | [`LeftMultiplicationGroupAction`](@ref) | a slightly more general type for all actions that are implemented by (matrix) multiplication |
| `RotationTranslationActionOnVector` | [`LeftMultiplicationGroupAction`](@ref) | e.g. in a [`GroupAction`](@ref) with [`SpecialEuclideanGroup`](@ref) and [`Euclidean`](@extref `Manifolds.Euclidean`)`(n)`. |
| `RowwiseMultiplicationAction` | [`RowwiseGroupAction`](@ref)`(`[`LeftMultiplicationGroupAction`](@ref)`)` | |
| `SemidirectProductGroup(G, H, a)` | [`LeftSemidirectProductLieGroup`](@ref)`(G, H, a)` | While this staid the same, there is now also the [`default_left_action`](@ref)`(G,H)`. When this agrees with `a` you can use the short hand `G⋉H` to generate this semidirect product. Analogously there now also exists the [`RightSemidirectProductLieGroup`](@ref)`(G,H)` with[`default_left_action`](@ref)`(G,H)` that allows for the short cut `G⋊H` |
| `SpecialEuclidean(n)` | [`SpecialEuclideanGroup`](@ref)`(n; variant=:right)` | |
| `SpecialEuclideanInGeneralLinear(n)` | [`SpecialEuclideanGroup`](@ref)`(n; variant=:right)` | `SpecialEuclideanGroup` supports both `ArrayPartition` and matrix representations. Point and tangent vectors can be converted between them using `convert`, which replaces old `project` and `embed` on `SpecialEuclideanInGeneralLinear`, see [`SpecialEuclideanMatrixPoint`](@ref), [`SpecialEuclideanMatrixTangentVector`](@ref), [`SpecialEuclideanProductPoint`](@ref) and [`SpecialEuclideanProductTangentVector`](@ref). |
| `TranslationAction` | [`AdditionGroupAction`](@ref) | this slightly more general name allows to reuse the action in other contexts easier |
| `VeeOrthogonalBasis` | [`DefaultLieAlgebraOrthogonalBasis`](@ref) | |
| `adjoint_action` | [`adjoint`](@ref) | now implemented with a default, when you provide [`diff_conjugate!`](@ref).
| `adjoint_matrix(G, p, b)` | [`jacobian_conjugate`](@ref)`(G, p, e, b)` | `e` is either the [`Identity`](@ref)`(G)` or its [`identity_element`](@ref)`(G)` |
| `apply_diff` | [`diff_apply`](@ref) | modifiers (diff) come first, consistent with [`ManifoldsDiff.jl`](https://juliamanifolds.github.io/ManifoldDiff.jl/stable/) |
| `apply_diff_group` | [`diff_group_apply`](@ref) | modifiers (diff/group) come first, consistent with [`ManifoldsDiff.jl`](https://juliamanifolds.github.io/ManifoldDiff.jl/stable/) |
| | [`conjugate`](@ref), [`diff_conjugate`](@ref) | a new function to model ``c_g: \mathcal G → \mathcal G`` given by ``c_g(h) = g∘h∘g^{-1}`` |
| `differential_exp_argument_lie_approx` | - | Scheduled for update and renaming. Though available in `ManifoldDiff.jl` for `GroupManifolds`, that will move to `differential_exp_argument_approx` instead, since `exp_lie` changed to now just `exp`. |
| `exp(G, g, X)` | `exp(`[`base_manifold`](@ref base_manifold(G::LieGroup))`(G), g, X)` | the previous defaults whenever not agreeing with the Riemannian one can now be accessed on the internal manifold |
| `exp_inv(G, g, X)` | [`exp`](@ref exp(G::LieGroup, g, X))`(G, g, X)`  | the exponential map invariant to the group operation is the default on Lie groups here |
| `exp_lie(G, X)` | [`exp`](@ref exp(G::LieGroup, X))`(G, X)` | the (matrix/Lie group) exponential |
| `get_coordinates(G, p, X)` | [`get_coordinates`](@ref get_coordinates(G::LieAlgebra))`(`[`LieAlgebra`](@ref)`(G), X)` | hat/vee moved to using the new [`LieAlgebra`](@ref). The old format should still work. |
| `get_vector(G, p, c)` | [`get_vector`](@ref get_vector(::LieAlgebra))`(`[`LieAlgebra`](@ref)`(G), c[, T])` | moved to using the new [`LieAlgebra`](@ref). The old format should still work. `T` indicates the type of tangent vector. |
| `hat(G, p, c)` | [`hat`](@ref hat(G::LieAlgebra))`(`[`LieAlgebra`](@ref)`(G), p, c[, T])` | hat/vee moved to using the new [`LieAlgebra`](@ref). The old format should still work. `T` indicates the type of tangent vector. |
| `inner(G, g, X, Y)` | [`inner`](@ref)`(`[`LieAlgebra`](@ref)`(G), X, Y)` | the inner product on the Lie Algebra. The old variant still calls the new one.|
| `inverse_translate(G, g, h, c)` | [`inv_left_compose`](@ref)`(G, g, h)`, [`inv_right_compose`](@ref)`(G, g, h)` | compute ``g^{-1}∘h`` and ``g∘h^{-1}``, resp. |
| `inverse_translate_diff(G, g, h, X, LeftForwardAction())` | - | discontinued, use `diff_left_compose(G, inv(G,g), h)` |
| `inverse_translate_diff(G, g, h, X, RightBackwardAction())` | - | discontinued, use `diff_left_compose(G, h, inv(G,g))` |
| `jacobian_exp_argument(G, g, X, b)` | [`jacobian_exp`](@ref)`(G, g, X, b)` | the Jacobian of the exponential map w.r.t. an [`AbstractBasis`](@extref `ManifoldsBase.AbstractBasis`) of the [`LieAlgebra`](@ref). The old name is resevered for the Riemannian exponential map. |
| `log(G, g, h)` | `log(`[`base_manifold`](@ref base_manifold(G::LieGroup))`(G), g, h)` | you can now access the previous defaults on the internal manifold whenever they do not agree with the invariant one |
| `log_inv(G, g, h)` | [`log`](@ref log(G::LieGroup, g, h))`(G, g, h)` | the logarithmic map invariant to the group operation is the default on Lie groups here |
| `log_lie(G, g)` | [`log`](@ref log(G::LieGroup, g))`(G, g)` | the (matrix/Lie group) logarithm |
| `norm(G, p, X)` | [`norm`](@extref ManifoldsBase :jl:function:`LinearAlgebra.norm`)`(`[`LieAlgebra`](@ref)`(G), X)` | the norm product on the Lie Algebra. The old variant still calls the new one. |
| `switch_direction(A)` | [`inv`](@ref inv(::GroupAction))`(A)` | switches from an action to its inverse action (formerly the direction forward/backward, sometimes even left/right, do not confuse with the side left/right). |
| `switch_side(A)` | - | discontinued |
| `translate(G, g, h)` | [`compose`](@ref)`(G, g, h)` | unified to `compose` |
| `translate_diff(G, g, X, c)` | [`diff_left_compose`](@ref)`(G, g, h, X)`, [`diff_right_compose`](@ref)`(G, g, h, X)` | for compose ``g∘h`` the functions now specify whether the derivative is taken w.r.t. to the left (`g`) or right (`h`) argument |
| `vee(G, p, X)` | [`vee`](@ref vee(G::LieAlgebra, X))`(`[`LieAlgebra`](@ref)`(G), X)` | hat/vee moved to using the new [`LieAlgebra`](@ref). The old format should still work. |

## Further notable changes

1. In general the default for tangent vectors is now to represent them in the [`LieAlgebra`](@ref), which obtains its own name now, though defined as a constant of a certain tangent space.
2. In accordance with point 1., the [`GeneralLinearGroup`](@ref) (formerly `GeneralLinear`) switched to using its Lie algebra to represent tangent vectors.
3. Formerly, both a power manifold of Lie groups as a manifold as well as a Lie group of a power manifold as a Lie group were possible. This is unified to just defining `G^n` as the Lie group on the power manifold with the element-wise group operation.
4. Formerly, product manifolds were stored as a [`ProductManifold`](@extref) of Lie groups and an indicator for the group operation, that the direct product should be used. This is switched to internally only store a [`ProductManifold`](@extref) as well as a (new) [`ProductGroupOperation`](@ref) that specifies one group operation for every factor.
5. The last two points achieve one unified modelling aspect of Lie groups: they are now always a manifold `M` together with a group operation `op`, but a Lie group does not store another Lie group (or product of them) internally with one common type they share, a [`LieGroup`](@ref).
6. Formerly, [`identity_element`](@ref)`(group[, gT])` to provide a materialized version of the [`Identity`](@ref) point on a Lie group accepted either a concrete point `g` on the Lie group _or_ the `Type` `T` of `g`. This is now unified for example with [`zero_vector`](@ref) and only accepts a `Type` as second argument.
7. There are two usual representations of the [`SpecialEuclideanGroup`](@ref), one with the rotation matrix on the left and the translation vector on the right and vice versa. The default is now to have the rotation matrix on the left, while formerly it was on the right. To be precise, [`SpecialEuclideanGroup`](@ref)`(3)` is now equivalent to [`SpecialOrthogonalGroup`](@ref)`(3)` [`⋉`](@ref LeftSemidirectProductGroupOperation) [`TranslationGroup`](@ref)`(3)`. To get the previous default use `SpecialEuclideanGroup(; variant=:right)` which is the same as writing [`TranslationGroup`](@ref)`(3)` [`⋊`](@ref RightSemidirectProductGroupOperation) [`SpecialOrthogonalGroup`](@ref)`(3)`.
