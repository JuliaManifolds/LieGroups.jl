using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    G = TranslationGroup(3)
    # Later maybe via auto-discover?
    p1, p2, p3 = [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.1, 1.2, 3.3]
    X1, X2, X3 = [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.1, 0.2, 0.3]
    properties = Dict(
        :Name => "The Translation group",
        :Points => [p1, p2, p3],
        :Vectors => [X1, X2, X3],
        :Functions => [compose, inv, show, is_identity],
    )
    expectations = Dict(
        :repr => "TranslationGroup(Euclidean(3; field=ℝ), AdditionGroupOperation())"
    )
    test_LieGroup(G, properties, expectations)
end
