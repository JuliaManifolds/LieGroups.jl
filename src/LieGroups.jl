module LieGroups

using ManifoldsBase, Manifolds

# before removing them from Manifolds.jl, let's for consistency include them here but
# overwrite them for now. This makes moving away from that (in Manifolds.jl 0.11) here non-breaking.
import Manifolds: apply, apply!, identity_element
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

export apply, apply!, diff_apply, diff_apply!
export base_manifold
export adjoint,
    adjoint!,
    compose,
    compose!,
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
export exp, exp!, log, log!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export base_Lie_group, base_manifold
export lie_bracket, lie_bracket!
export inv
end # module LieGroups
