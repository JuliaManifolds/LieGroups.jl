using LieGroups, Test, ManifoldsBase, Random, LinearAlgebra

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
            jacobian_exp,
            lie_bracket,
            log,
            rand,
            vee,
        ],
    )
    expectations2 = Dict(
        :repr => "PowerLieGroup(LieGroup(LieGroups.Test.DummyManifold(), LieGroups.Test.DummyOperation()), 2)",
        # TranslationGroup(2)² is flat and Abelian, so the Jacobian of exp is the identity
        :jacobian_exp => Matrix{Float64}(LinearAlgebra.I, 4, 4),
    )
    LieGroups.Test.test_lie_group(pG2, properties2, expectations2)
    @testset "jacobian_exp is block-diagonal for components" begin
        B = DefaultLieAlgebraOrthogonalBasis()
        Gp = PowerLieGroup(SpecialOrthogonalGroup(3), NestedPowerRepresentation(), 2)
        𝔤p = LieAlgebra(Gp)
        Xcc = [0.3, -0.2, 0.5, 0.7, -0.4, 0.6]
        Xt = hat(𝔤p, Xcc)
        Jp = jacobian_exp(Gp, Xt)
        @test size(Jp) == (6, 6)
        # off-diagonal blocks
        @test iszero(Jp[1:3, 4:6])
        @test iszero(Jp[4:6, 1:3])
        # diagonal blocks equal the base group's own jacobian_exp
        G1 = SpecialOrthogonalGroup(3)
        @test isapprox(Jp[1:3, 1:3], jacobian_exp(G1, hat(LieAlgebra(G1), Xcc[1:3])))
        @test isapprox(Jp[4:6, 4:6], jacobian_exp(G1, hat(LieAlgebra(G1), Xcc[4:6])))
        Jp2 = copy(Jp)
        jacobian_exp!(Gp, Jp2, Xt, B)
        @test isapprox(Jp, Jp2)
    end
    @testset "Special cases with identity" begin
        e = Identity(pG)
        @test ManifoldsBase.check_size(pG, e) == nothing
        eF = Identity(AdditionGroupOperation())
        @test ManifoldsBase.check_size(pG, eF) isa DomainError
    end
end
