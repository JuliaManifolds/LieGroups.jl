using LieGroups, Test, ManifoldsBase

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic product Lie group" begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    G2 = G × G

    properties = Dict(
        :Name => "The Product Manifold",
        # :Rng => Random.MersenneTwister(),
        :Functions => [
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
            #rand,
            show,
            #vee,
        ],
    )
    expectations = Dict(
        :repr => "ProductLieGroup(LieGroupsTestSuite.DummyManifold() × LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation() × LieGroupsTestSuite.DummyOperation())",
    )
    test_lie_group(G2, properties, expectations)

    @testset "Product Operation generators" begin
        op2 = LieGroupsTestSuite.DummySecondOperation()
        O1 = op × op2
        O2 = op2 × op
        @test (O1 × op) == (op × O2)
        @test (O1 × O2) == (op × op2 × op2 × op)
    end
end
