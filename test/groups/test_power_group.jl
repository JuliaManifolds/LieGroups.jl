using LieGroups, Test, ManifoldsBase, Random

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic power Lie group" begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    pG = G^2

    properties1 = Dict(:Name => "The generic Power Manifold", :Functions => [show])
    expectations1 = Dict(
        :repr => "PowerLieGroup(LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()), 2)",
    )
    test_lie_group(pG, properties1, expectations1)

    # Explicit one to test element-wise methods
    pG2 = PowerLieGroup(TranslationGroup(2), NestedPowerRepresentation(), 2)
    g, h = [[1.0, 0.0], [0.0, 3.0]], [[0.0, 1.0], [2.0, 0.0]]
    X, Y = [[0.0, 0.1], [0.2, 0.0]], [[0.1, 0.2], [0.0, 0.3]]
    properties2 = Dict(
        :Name => "The generic nested Power Manifold",
        :Points => [g, h],
        :Vectors => [X, Y],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            compose,
            conjugate,
            diff_conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            # hat,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            # lie_bracket,
            log,
            rand,
            # vee,
        ],
    )
    expectations2 = Dict(
        :repr => "PowerLieGroup(LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()), 2)",
    )
    test_lie_group(pG2, properties2, expectations2)
end
