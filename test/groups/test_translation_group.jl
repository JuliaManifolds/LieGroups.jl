using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    G = TranslationGroup(3)
    # Later maybe via auto-discover?
    g1, g2, g3 = [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.1, 1.2, 3.3]
    X1, X2, X3 = [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.1, 0.2, 0.3]
    properties = Dict(
        :Name => "The Translation group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Functions => [
            adjoint,
            compose,
            conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            identity_element,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            show,
        ],
    )
    expectations = Dict(
        :repr => "TranslationGroup(3; field=ℝ)",
        :diff_inv => -X1,
        :diff_left_compose => X1,
        :diff_right_compose => X1,
        :lie_bracket => zero(X1),
    )
    test_lie_group(G, properties, expectations)

    properties2 = Dict(
        :AlgebraVectors => [X1, X2, X3],
        :Functions =>
            [apply, diff_apply, diff_group_apply, base_lie_group, base_manifold, show],
        :GroupPoints => [g1, g2, g3],
        :ManifoldPoints => [g1, g2, g3],
        :TangentVectors => [X1, X2, X3],
        :Name => "",
    )
    expectations2 = Dict(:manifold => G, :group => G, :repr => "")
    @testset "Translation group operation action" begin
        # A first group action Test
        for t in [
            RightGroupOperation(),
            LeftGroupOperation(),
            InverseLeftGroupOperation(),
            InverseRightGroupOperation(),
        ]
            A = GroupOperationAction(t, G)
            properties2[:Name] = "with $A"
            expectations2[:repr] = "GroupOperationAction($t, $G)"
            test_group_action(A, properties2, expectations2)
        end
    end
end
