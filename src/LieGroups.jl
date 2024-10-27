module LieGroups

using ManifoldsBase, Manifolds, LinearAlgebra

#
#
# = Compatibility (and a bit of type piracy for now)
# The following imports are necessary to use Manifolds.jl 0.10 with Lie groups
# The line is removed when the Groups are removed from possibly 0.11
import Manifolds:
    apply, apply!, compose, identity_element, is_identity, hat, hat!, vee, vee!
# Both define the following structs, so these for now lead to asking for explicit prefixes
# Manifolds: Identity, TranslationGroup
include("documentation_glossary.jl")
include("interface.jl")
include("Lie_algebra/Lie_algebra_interface.jl")
# Generic Operations
include("group_operations/addition_operation.jl")
include("group_operations/multiplication_operation.jl")

# Actions
include("group_actions/group_action_interface.jl")
include("group_actions/group_operation_action.jl")

# Lie groups
include("groups/translation_group.jl")
include("groups/general_linear_group.jl")

export LieGroup, LieAlgebra
export LieAlgebraOrthogonalBasis
#
#
# Group Operations
export AbstractGroupOperation, Identity
export AdditionGroupOperation
export AbstractMultiplicationGroupOperation
export MatrixMultiplicationGroupOperation

#
#
# Group Actions
export AbstractGroupActionType, AbstractGroupAction
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export LeftGroupOperation, RightGroupOperation
export InverseLeftGroupOperation, InverseRightGroupOperation
export GroupOperationAction

#
#
# Specific groups
export TranslationGroup, GeneralLinearGroup

export adjoint, adjoint!, apply, apply!
export base_lie_group, base_manifold
export compose, compose!
export diff_apply,
    diff_apply!,
    diff_group_apply,
    diff_group_apply!,
    diff_left_compose,
    diff_left_compose!,
    diff_right_compose,
    diff_right_compose!
export get_coordinates, get_coordinates!, get_vector, get_vector!
export hat, hat!
export inv_left_compose, inv_left_compose!, inv_right_compose, inv_right_compose!
export isapprox, is_point, is_vector
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export exp, exp!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export lie_bracket, lie_bracket!, log, log!
export norm
export switch
export vee, vee!
export zero_vector, zero_vector!
end # module LieGroups
