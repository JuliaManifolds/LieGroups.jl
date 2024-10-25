using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    G = GeneralLinearGroup(2)
    # Later maybe via auto-discover?
    g1, g2, g3 = [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.1, 1.2, 3.3]
    X1, X2, X3 = [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.1, 0.2, 0.3]
    properties = Dict(
        :Name => "The general linear group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Functions => [show],
    )
    expectations = Dict(:repr => "GeneralLinearGroup(2; field=ℝ)")
    test_lie_group(G, properties, expectations)
end
