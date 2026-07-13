using LieGroups, Test, ManifoldsBase, Random, RecursiveArrayTools, LinearAlgebra

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
            jacobian_exp,
            lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    @test LieGroups.submanifold_components(G, Identity(G)) === (Identity{AdditionGroupOperation}(), Identity{AdditionGroupOperation}())
    expectations = Dict(
        :repr => "ProductLieGroup(Euclidean(2; field=ℝ) × Euclidean(2; field=ℝ), AdditionGroupOperation() × AdditionGroupOperation())",
        # both factors are flat and Abelian, so the Jacobian of exp is the identity
        :jacobian_exp => Matrix{Float64}(LinearAlgebra.I, 4, 4),
    )
    LieGroups.Test.test_lie_group(G, properties, expectations)
    @testset "jacobian_exp is block-diagonal over the factors" begin
        B = DefaultLieAlgebraOrthogonalBasis()
        # a product with curvature: SO(3) × SE(2) (right variant)
        Gc = SpecialOrthogonalGroup(3) × SpecialEuclideanGroup(2; variant = :right)
        𝔤c = LieAlgebra(Gc)
        Xcc = [0.3, -0.2, 0.5, 0.7, -0.4, 0.6]
        Xt = hat(𝔤c, Xcc)
        gc = exp(Gc, Xt)
        Jc = jacobian_exp(Gc, gc, Xt)
        @test size(Jc) == (6, 6)
        # off-diagonal coupling blocks vanish
        @test iszero(Jc[1:3, 4:6])
        @test iszero(Jc[4:6, 1:3])
        # diagonal blocks equal each factor's own jacobian_exp
        Gc1 = SpecialOrthogonalGroup(3)
        Xc1 = hat(LieAlgebra(Gc1), Xcc[1:3])
        @test isapprox(Jc[1:3, 1:3], jacobian_exp(Gc1, exp(Gc1, Xc1), Xc1))
        Gc2 = SpecialEuclideanGroup(2; variant = :right)
        Xc2 = hat(LieAlgebra(Gc2), Xcc[4:6])
        @test isapprox(Jc[4:6, 4:6], jacobian_exp(Gc2, exp(Gc2, Xc2), Xc2))
        # mutating matches allocating
        Jc2 = copy(Jc)
        jacobian_exp!(Gc, Jc2, gc, Xt, B)
        @test isapprox(Jc, Jc2)
    end
    @testset "A small additional size check" begin
        @test ManifoldsBase.check_size(G, Identity(G)) === nothing
        @test ManifoldsBase.check_size(G, Identity(G), X) === nothing
    end
    @testset "Product Operation generators" begin
        G = LieGroups.Test.DummyLieGroup()
        G2 = G × G
        @test ProductLieGroup(G2) === G2 #Product groups are not “wrapped twice”
        op = LieGroups.Test.DummyOperation()
        op2 = LieGroups.Test.DummySecondOperation()
        O1 = op × op2
        O2 = op2 × op
        @test (O1 × op) == (op × O2)
        @test (O1 × O2) == (op × op2 × op2 × op)
        @test O1[1] == op
        @test O1[2] == op2
        @test O1[:] == (op, op2)
        @test O1[:] == LieGroups.submanifold_components(G2, O1)
    end
    @testset "× splashes" begin
        G = LieGroups.Test.DummyLieGroup()
        G2 = G × G
        G3 = G × G × G
        H = G × G × G × G
        # one or two products in cross -> splat.
        @test G3 × G == H
        @test G × G3 == H
        @test G2 × G2 == H
        # Check constructor for more than 2 Lie groups
        @test ProductLieGroup(G, G, G, G) == H
    end
end
