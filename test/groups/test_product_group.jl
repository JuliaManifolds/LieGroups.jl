using LieGroups, Test, ManifoldsBase, Random, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic product Lie group" begin
    G = TranslationGroup(2) × TranslationGroup(2)
    g, h = ArrayPartition([1.0, 0.0], [0.0, 3.0]), ArrayPartition([0.0, 1.0], [2.0, 0.0])
    X, Y = ArrayPartition([0.0, 0.1], [0.2, 0.0]), ArrayPartition([0.1, 0.2], [0.0, 0.3])

    properties = Dict(
        :Name => "The Product Manifold",
        :Rng => Random.MersenneTwister(),
        :Points => [g, h],
        :Vectors => [X, Y],
        :Functions => [
            compose,
            conjugate,
            diff_conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            hat,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectations = Dict(
        :repr => "ProductLieGroup(Euclidean(2; field=ℝ) × Euclidean(2; field=ℝ), AdditionGroupOperation() × AdditionGroupOperation())",
    )
    test_lie_group(G, properties, expectations)

    @testset "Product Operation generators" begin
        op = LieGroupsTestSuite.DummyOperation()
        op2 = LieGroupsTestSuite.DummySecondOperation()
        O1 = op × op2
        O2 = op2 × op
        @test (O1 × op) == (op × O2)
        @test (O1 × O2) == (op × op2 × op2 × op)
    end
end
