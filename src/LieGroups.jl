@doc raw"""
LieGroups.jl: Lie groups and Lie algebras in Julia.

The package is named after the norwegian mathematician [Sophus Lie](https://en.wikipedia.org/wiki/Sophus_Lie).

* 📚 Documentation: [manoptjl.org](https://juliamanifolds.github.io/LieGroups.jl/dev/)
* 📦 Repository: [github.com/JuliaManifolds/LieGroups.jl](https://github.com/JuliaManifolds/LieGroups.jl)
* 💬 Discussions: [github.com/JuliaManifolds/LieGroups.jl/discussions](https://github.com/JuliaManifolds/LieGroups.jl/discussions)
* 🎯 Issues: [github.com/JuliaManifolds/LieGroups.jl/issues](https://github.com/JuliaManifolds/LieGroups.jl/issues)
"""
module LieGroups

using LinearAlgebra, ManifoldsBase, Manifolds, Random

using ManifoldsBase: RealNumbers

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

# Meta Lie groups
include("groups/power_group.jl")
include("groups/product_group.jl")
include("groups/semidirect_product_group.jl")
# Lie groups
include("groups/translation_group.jl")
include("groups/general_linear_group.jl")
include("groups/heisenberg_group.jl")

# explicit method error to avoid stack overflow
for GT in [LieGroup, HeisenbergGroup]
    @eval begin
        @doc "$(_doc_log_id)"
        function ManifoldsBase.log!(G::$GT, X, e::Identity, g)
            throw(
                MethodError(
                    ManifoldsBase.log!, (typeof(G), typeof(X), typeof(e), typeof(g))
                ),
            )
        end
    end
end

export LieGroup, LieAlgebra
export PowerLieGroup, ProductLieGroup
export LeftSemidirectProductLieGroup, RightSemidirectProductLieGroup
export LieAlgebraOrthogonalBasis
export ×, ^, ⋉, ⋊
#
#
# Group Operations
export AbstractGroupOperation, Identity
export AdditionGroupOperation
export AbstractMultiplicationGroupOperation
export MatrixMultiplicationGroupOperation
export ProductGroupOperation
export LeftSemidirectProductGroupOperation, RightSemidirectProductGroupOperation

#
#
# Group Actions
export AbstractGroupActionType
export AbstractLeftGroupActionType, AbstractRightGroupActionType
export LeftGroupOperationAction, RightGroupOperationAction
export InverseLeftGroupOperationAction, InverseRightGroupOperationAction
export GroupAction, GroupOperationAction
#
#
# Specific groups
export TranslationGroup, GeneralLinearGroup, HeisenbergGroup

export adjoint, adjoint!, apply, apply!
export base_lie_group, base_manifold
export compose, compose!
export default_left_action,
    default_right_action,
    det,
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
export conjugate, conjugate!, diff_conjugate, diff_conjugate!
export exp, exp!
export identity_element, identity_element!, is_identity, inv, inv!, diff_inv, diff_inv!
export lie_bracket, lie_bracket!, log, log!
export manifold_dimension
export norm
export rand, rand!
export switch
export vee, vee!
export zero_vector, zero_vector!
end # module LieGroups
