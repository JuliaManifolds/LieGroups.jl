module LieGroups

using ManifoldsBase, Manifolds

include("documentation_glossary.jl")
include("interface.jl")
include("additive_group_operation.jl")

include("groups/additive.jl")

export LieGroup, LieAlgebra

export AbstractGroupOperation, Identity
export AdditiveGroupOperation

export AdditiveGroup

export base_manifold
export compose,
    compose!,
    compose_diff_left,
    compose_diff_left!,
    compose_diff_right,
    compose_diff_right!,
    compose_inv_left,
    compose_inv_left!,
    compose_inv_right,
    compose_inv_right!,
    conjugate,
    conjugate!
export exp, exp!, log, log!
export identity_element, identity_element!, is_identity, inv, inv!, inv_diff, inv_diff!
export Lie_bracket, Lie_bracket!
end # module LieGroups
