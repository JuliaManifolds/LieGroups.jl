"""
    LieGroups.Test

This module provides tools and dummy structures to test functionality provided
within `LieGroups.jl`.

For every test function, several interactions to other functions can be activated.

The following functions are expected to be available, since their defaults just pass through to the manifold

* `is_point` both on the Lie group `G` and the Lie algebra `𝔤`
* `isapprox(G, g, h)` and `isapprox(𝔤, X, Y)`
* `copy(G, g)`
* `norm(𝔤, X)`

Most functionality only gets loaded once `Test.jl` is loaded as well, i.e.
that the functions are populated with methods by the extension.
"""
module Test
using ..LieGroups
using LieGroups: AbstractLieGroup, LieGroup, AbstractGroupOperation
using LinearAlgebra: I


struct DummyOperation <: LieGroups.AbstractGroupOperation end
struct DummySecondOperation <: LieGroups.AbstractGroupOperation end
struct DummyManifold <: LieGroups.AbstractManifold{LieGroups.ℝ} end
struct DummyActionType <: LieGroups.AbstractGroupActionType end
struct DummyLeftActionType <: LieGroups.AbstractLeftGroupActionType end
struct DummyRightActionType <: LieGroups.AbstractRightGroupActionType end
struct DummyMetric <: LieGroups.AbstractMetric end
const DummyLieGroup = LieGroup{LieGroups.ℝ, DummyOperation, DummyManifold}
DummyLieGroup() = LieGroup(DummyManifold(), DummyOperation())
LieGroups.switch(a::DummyActionType) = a
LieGroups.switch(::DummyLeftActionType) = DummyRightActionType()
LieGroups.switch(::DummyRightActionType) = DummyLeftActionType()

include("so4_edge_cases.jl")

"""
    rotate_matrix(n, k1, k2, α)

Generate a rotation matrix in ``ℝ^{n×n}`` with a rotation about ``α`` (in radians)
in the ``k_1-k_2``-plane.
"""
function rotation_matrix(n, k1, k2, α)
    R = Matrix{Float64}(I, n, n)
    R[k1, k1] = cos(α)
    R[k2, k2] = cos(α)
    R[k1, k2] = sin(α)
    R[k2, k1] = -sin(α)
    return R
end

function test_adjoint end
function test_apply end
function test_compose end
function test_conjugate end
function test_copyto end
function test_diff_apply end
function test_diff_group_apply end
function test_diff_inv end
function test_diff_left_compose end
function test_diff_right_compose end
function test_diff_conjugate end
function test_exp_log end
function test_hat_vee end
function test_identity_element end
function test_injectivity_radius end
function test_is_flat end
function test_inv_compose end
function test_inv end
function test_identity end
function test_is_identity end
function test_inner end
function test_jacobian_conjugate end
function test_jacobian_exp end
function test_lie_bracket end
function test_norm end
function test_push_pull_tangent end
function test_rand end
function test_show end

function test_lie_group end
function test_group_action end
end # module Test
