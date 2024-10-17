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
    test_LieGroup(G, properties, expectations)
end
