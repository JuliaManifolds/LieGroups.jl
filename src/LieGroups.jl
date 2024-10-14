module LieGroups

using ManifoldsBase, Manifolds, LinearAlgebra

#
#
# = Compatibility (and a bit of type piracy for now)
# The following imports shall be removed once Manifolds.jl does no longer define them (in 0.11)
import Manifolds: apply, apply!, identity_element, is_identity, compose
# Both define the following structs, so these for now lead to asking for explicit prefixes
# Manifolds: Identity, TranslationGroup
include("documentation_glossary.jl")
include("interface.jl")

# Generic Operations
include("group_operations/addition.jl")

# Actions
include("group_actions/group_action_interface.jl")
include("group_actions/group_operation_action.jl")

# Lie groups
include("groups/translation_group.jl")

export LieGroup, LieAlgebra

export AbstractGroupOperation, Identity
export AdditionGroupOperation

export AbstractGroupActionType, AbstractGroupAction
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export LeftGroupOperation, RightGroupOperation
export InverseLeftGroupOperation, InverseRightGroupOperation
export GroupOperationAction

export TranslationGroup

export adjoint, adjoint!, apply, apply!
export base_Lie_group, base_manifold
export compose, compose!
export diff_apply,
    diff_apply!,
    diff_left_compose,
    diff_left_compose!,
    diff_right_compose,
    diff_right_compose!,
    inv_left_compose,
    inv_left_compose!,
    inv_right_compose,
    inv_right_compose!
export isapprox, is_point, is_vector
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export exp, exp!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export lie_bracket, lie_bracket!, log, log!
export norm
export zero_vector
end # module LieGroups
