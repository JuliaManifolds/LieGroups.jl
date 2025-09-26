@doc raw"""
LieGroups.jl: Lie groups and Lie algebras in Julia.

The package is named after the Norwegian mathematician [Marius Sophus Lie](https://en.wikipedia.org/wiki/Sophus_Lie) (1842–1899).

* 📚 Documentation: [juliamanifolds.github.io/LieGroups.jl/](https://juliamanifolds.github.io/LieGroups.jl/)
* 📦 Repository: [github.com/JuliaManifolds/LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl)
* 💬 Discussions: [github.com/JuliaManifolds/LieGroups.jl/discussions](https://github.com/JuliaManifolds/LieGroups.jl/discussions)
* 🎯 Issues: [github.com/JuliaManifolds/LieGroups.jl/issues](https://github.com/JuliaManifolds/LieGroups.jl/issues)
"""
module LieGroups

using LinearAlgebra, ManifoldsBase, Quaternions, StaticArrays, Random

#
#
# For the intermediate time be a bit more careful with Manifolds.jl and really only import
# what we need.
using Manifolds:
    Circle,
    DeterminantOneMatrices,
    Euclidean,
    GeneralUnitaryMatrices,
    InvertibleMatrices,
    HeisenbergMatrices,
    OrthogonalMatrices,
    ProductManifold,
    Rotations,
    Sphere,
    SymplecticMatrices,
    UnitaryMatrices
using Manifolds: DeterminantOneMatrixType
using Manifolds: base_manifold

import LinearAlgebra: adjoint, adjoint!
using ManifoldsBase:
    AbstractBasis, AbstractNumbers, RealNumbers, ComplexNumbers, QuaternionNumbers, ℍ, ℝ, ℂ
using ManifoldsBase:
    allocate_result,
    get_parameter,
    internal_value,
    submanifold_component,
    submanifold_components,
    tangent_vector_type
using StaticArrays
#
#
# = Compatibility (and a bit of type piracy for now)
# The following imports are necessary to use Manifolds.jl 0.10 with Lie groups
# The line is removed when the Groups are removed from possibly 0.11
import Manifolds: apply, apply!, compose, identity_element, is_identity
# Both define the following structs, so these for now lead to asking for explicit prefixes
# Manifolds: Identity, TranslationGroup
#
include("documentation_glossary.jl")
include("utils.jl")
include("interface.jl")
include("Lie_algebra/Lie_algebra_interface.jl")
# Generic Operations
include("group_operations/addition_operation.jl")
include("group_operations/multiplication_operation.jl")
include("group_operations/multiplication_operation_abelian.jl")

# Actions
include("group_actions/group_action_interface.jl")
include("group_actions/addition_action.jl")
include("group_actions/group_operation_action.jl")
include("group_actions/multiplication_action.jl")
include("group_actions/rotation_around_axis_action.jl")
# These might use some of the above so we include them second
include("group_actions/columnwise_action.jl")
include("group_actions/rowwise_action.jl")

# Meta Lie groups
include("groups/power_group.jl")
include("groups/product_group.jl")
include("groups/semidirect_product_group.jl")
include("groups/validation_group.jl")

# Lie groups
include("groups/translation_group.jl")
include("groups/general_linear_group.jl")
include("groups/circle_group_complex.jl")
include("groups/circle_group_real.jl")
include("groups/circle_group_sphere.jl")
include("groups/circle_group.jl")
include("groups/heisenberg_group.jl")
include("groups/special_linear_group.jl")
include("groups/symplectic_group.jl")

# includes generic implementations for O(n), U(n), SO(n), SO(n), so we load this first
include("groups/unitary_group.jl")
include("groups/orthogonal_group.jl")
include("groups/special_unitary_group.jl")
include("groups/special_orthogonal_group.jl")

# Products of Groups
include("groups/special_euclidean_group.jl")

export AbstractLieGroup
export LieGroup, LieAlgebra
export PowerLieGroup, ProductLieGroup
export LeftSemidirectProductLieGroup, RightSemidirectProductLieGroup
export ValidationLieGroup
export DefaultLieAlgebraOrthogonalBasis
export ×, ^, ⋉, ⋊
#
#
# Group Operations
export AbstractGroupOperation, Identity
export AdditionGroupOperation
export AbstractMultiplicationGroupOperation
export MatrixMultiplicationGroupOperation
export AbelianMultiplicationGroupOperation
export PowerGroupOperation, ProductGroupOperation
export LeftSemidirectProductGroupOperation, RightSemidirectProductGroupOperation

#
#
# Group Actions
export AbstractGroupActionType
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export AdditionGroupAction
export LeftGroupOperationAction, RightGroupOperationAction
export GroupAction, GroupOperationAction
export LeftMultiplicationGroupAction
export InverseLeftGroupOperationAction, InverseRightGroupOperationAction
export AbstractActionActsOnType, ActionActsOnLeft, ActionActsOnRight
export RotationAroundAxisAction
export ColumnwiseGroupAction, RowwiseGroupAction

#
#
# Specific groups
export CircleGroup
export GeneralLinearGroup
export HeisenbergGroup
export OrthogonalGroup
export SpecialEuclideanGroup, SpecialLinearGroup
export SpecialOrthogonalGroup, SpecialUnitaryGroup
export SymplecticGroup
export TranslationGroup
export UnitaryGroup

# Points and Tangent representations
export AbstractLieGroupPoint, AbstractLieAlgebraTangentVector
export SpecialEuclideanMatrixPoint, SpecialEuclideanMatrixTangentVector
export SpecialEuclideanProductPoint, SpecialEuclideanProductTangentVector
export ValidationMPoint, ValidationLieAlgebraTangentVector

export BaseManifoldInverseRetraction,
    BaseManifoldRetraction, BaseManifoldVectorTransportMethod
# Errors
export CompositeManifoldError
# Functions
export adjoint, adjoint!, apply, apply!
export base_lie_group, base_manifold
export compose, compose!
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export default_left_action,
    default_right_action,
    diff_apply,
    diff_apply!,
    diff_group_apply,
    diff_group_apply!,
    diff_left_compose,
    diff_left_compose!,
    diff_right_compose,
    diff_right_compose!
export get_coordinates, get_coordinates!, get_vector, get_vector!
export hat, hat!
export inv, inv!, inv_left_compose, inv_left_compose!, inv_right_compose, inv_right_compose!
export isapprox, is_point, is_vector
export exp, exp!
export identity_element, identity_element!, is_identity, inner
export inv, inv!, diff_inv, diff_inv!
export injectivity_radius
export jacobian_conjugate, jacobian_conjugate!
export jacobian_exp, jacobian_exp!
export lie_bracket, lie_bracket!, log, log!
export manifold_dimension
export norm, number_of_coordinates
export push_forward_tangent, push_forward_tangent!, pull_back_tangent, pull_back_tangent!
export rand, rand!
export switch
export sym_rem
export vee, vee!
export zero_vector, zero_vector!
end # module LieGroups
