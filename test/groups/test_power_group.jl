using LieGroups, Test, ManifoldsBase, Random

@testset "Generic power Lie group" begin
    M = LieGroups.Test.DummyManifold()
    op = LieGroups.Test.DummyOperation()
    G = LieGroup(M, op)
    pG = G^2

    properties1 = Dict(:Name => "The generic Power Manifold", :Functions => [show])
    expectations1 = Dict(
        :repr => "PowerLieGroup(LieGroup(LieGroups.Test.DummyManifold(), LieGroups.Test.DummyOperation()), 2)",
    )
    LieGroups.Test.test_lie_group(pG, properties1, expectations1)

    # Explicit one to test element-wise methods
    pG2 = PowerLieGroup(TranslationGroup(2), NestedPowerRepresentation(), 2)
    g, h = [[1.0, 0.0], [0.0, 3.0]], [[0.0, 1.0], [2.0, 0.0]]
    X, Y = [[0.0, 0.1], [0.2, 0.0]], [[0.1, 0.2], [0.0, 0.3]]
    @testset "convencience access" begin
        @test g[pG2, 1] === g[1]
    end

    properties2 = Dict(
        :Name => "The generic nested Power Manifold",
        :Points => [g, h],
        :Vectors => [X, Y],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            adjoint,
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
            vee,
        ],
    )
    expectations2 = Dict(
        :repr => "PowerLieGroup(LieGroup(LieGroups.Test.DummyManifold(), LieGroups.Test.DummyOperation()), 2)",
    )
    LieGroups.Test.test_lie_group(pG2, properties2, expectations2)
    @testset "Special cases with identity" begin
        e = Identity(pG)
        @test ManifoldsBase.check_size(pG, e) == nothing
        eF = Identity(AdditionGroupOperation())
        @test ManifoldsBase.check_size(pG, eF) isa DomainError
    end
end
