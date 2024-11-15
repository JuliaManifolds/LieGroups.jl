using LieGroups, Test, ManifoldsBase

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic semidirect product Lie group" begin
    M = LieGroupsTestSuite.DummyManifold()
    op1 = LieGroupsTestSuite.DummyOperation()
    G1 = LieGroup(M, op1)
    op2 = LieGroupsTestSuite.DummySecondOperation()
    G2 = LieGroup(M, op2)

    fcts = [
        # compose,
        # conjugate,
        # diff_conjugate,
        # diff_inv,
        # diff_left_compose,
        # diff_right_compose,
        # exp,
        # hat,
        # inv,
        # inv_left_compose,
        # inv_right_compose,
        # is_identity,
        # lie_bracket,
        # log,
        # rand,
        show,
        #vee,
    ]

    Gl = LeftSemidirectProductLieGroup(G1, G2, LeftGroupOperationAction())
    properties = Dict(
        :Name => "The Power Manifold",
        # :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :repr => "LeftSemidirectProductLieGroup($(G1), $(G2), LeftGroupOperationAction())"
    )
    test_lie_group(Gl, properties, expectations)

    Gr = RightSemidirectProductLieGroup(G1, G2, RightGroupOperationAction())
    properties = Dict(
        :Name => "The Power Manifold",
        # :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :repr => "RightSemidirectProductLieGroup($(G1), $(G2), RightGroupOperationAction())"
    )
    test_lie_group(Gr, properties, expectations)
end
