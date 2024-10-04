module LieGroups

using ManifoldsBase, Manifolds

# before removing them from Manifolds.jl – let's for consistency include them here but
# overwrite them for now. This makes moving away from that (in Manifolds.jl 0.11) here nonbreaking.
import Manifolds: apply, apply!

include("documentation_glossary.jl")
include("interface.jl")
include("group_actions/group_action_interface.jl")
include("group_actions/group_operation_action.jl")
include("additive_group_operation.jl")

include("groups/additive.jl")

export LieGroup, LieAlgebra

export AbstractGroupOperation, Identity
export AdditiveGroupOperation

export AbstractGroupActionType, AbstractGroupAction
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export LeftGroupOperation, RightGroupOperation
export InverseLeftGroupOperation, InverseRightGroupOperation
export GroupOperationAction

export AdditiveGroup

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
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export exp, exp!, log, log!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export base_Lie_group, base_manifold
export lie_bracket, lie_bracket!
export inv
end # module LieGroups
