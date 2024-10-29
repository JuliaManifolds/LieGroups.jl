using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    G = GeneralLinearGroup(2)
    g1, g2, g3 = [2.0 0.0; 0.0 1.0], [1.0 0.5; 0.5 1.0], [1.0 2.0; 3.0 4.0]
    X1, X2, X3 = [-1.0 0.0; 0.0 0.0], [0.0 0.5; 0.5 0.0], [1.0 2.0; 3.0 4.0]
    properties = Dict(
        :Name => "The general linear group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Functions => [
            compose,
            conjugate,
            diff_conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            show,
        ],
    )
    expectations = Dict(
        :repr => "GeneralLinearGroup(2; field=ℝ)", :lie_bracket => X1 * X2 - X2 * X1
    )
    test_lie_group(G, properties, expectations)
end
