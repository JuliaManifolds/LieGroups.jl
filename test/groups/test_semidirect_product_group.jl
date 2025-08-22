using LieGroups, Test, ManifoldsBase, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic semidirect product Lie group" begin
    @testset "Show" begin
        M = LieGroupsTestSuite.DummyManifold()
        op1 = LieGroupsTestSuite.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroupsTestSuite.DummySecondOperation()
        G2 = LieGroup(M, op2)
        # The rest requires tests with a concrete semi-direct one
        fcts = [show]

        Gl = LeftSemidirectProductLieGroup(G1, G2, LeftGroupOperationAction())
        properties = Dict(
            :Name => "The Left Semidirect Product Manifold",
            # :Rng => Random.MersenneTwister(),
            :Functions => fcts,
        )
        expectations = Dict(
            :repr => "LeftSemidirectProductLieGroup($(G1), $(G2), LeftGroupOperationAction())",
        )
        test_lie_group(Gl, properties, expectations)

        Gr = RightSemidirectProductLieGroup(G1, G2, RightGroupOperationAction())
        properties = Dict(
            :Name => "The Right Semidirect Product Manifold",
            # :Rng => Random.MersenneTwister(),
            :Functions => fcts,
        )
        expectations = Dict(
            :repr => "RightSemidirectProductLieGroup($(G1), $(G2), RightGroupOperationAction())",
        )
        test_lie_group(Gr, properties, expectations)
    end
    @testset "Defaults" begin
        # Let's take a bit of a funny, but nondefault semidirect product
        # Possible homomorphisms for SO(2) ⋉ SO(2) is only the trivial action.
        # This reduces the semidirect product to a direct product.
        # TODO This test can probability just be deleted, as it is not really a semidirect product?
        struct TestTrivialAction <: AbstractGroupActionType end
        function LieGroups.apply!(::GroupAction{TestTrivialAction}, k, g, h)
            return k .= h
        end

        G = SpecialOrthogonalGroup(2)
        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        X1, X2, X3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 0.2 0.0], [0.0 0.0; 0.0 0.0]

        Gl = LeftSemidirectProductLieGroup(G, G, TestTrivialAction())
        Gr = RightSemidirectProductLieGroup(G, G, TestTrivialAction())

        h1 = ArrayPartition(copy(g1), copy(g2))
        h2 = ArrayPartition(copy(g2), copy(g1))
        h3 = ArrayPartition(copy(g3), copy(g3))
        Y1 = ArrayPartition(copy(X1), copy(X2))
        Y2 = ArrayPartition(copy(X2), copy(X1))
        Y3 = ArrayPartition(copy(X3), copy(X3))
        properties = Dict(
            :Name => "The generic left semi-direct product group",
            :Points => [h1, h2, h3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, inv, show],
        )
        expectations_l = Dict(
            :repr => "LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), SpecialOrthogonalGroup(2), TestTrivialAction())",
        )
        test_lie_group(Gl, properties, expectations_l)
        properties[:Name] = "The generic left semi-direct product group"
        expectations_r = Dict(
            :repr => "RightSemidirectProductLieGroup(SpecialOrthogonalGroup(2), SpecialOrthogonalGroup(2), TestTrivialAction())",
        )
        test_lie_group(Gr, properties, expectations_r)
    end
    @testset "Interaction os semidirect and product" begin
        M = LieGroupsTestSuite.DummyManifold()
        op1 = LieGroupsTestSuite.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroupsTestSuite.DummySecondOperation()
        G2 = LieGroup(M, op2)


        Gl = LeftSemidirectProductLieGroup(G1 × G2, G2 × G1, LeftGroupOperationAction())
        @test Gl.manifold isa ProductManifold
        @test length(Gl.manifold.manifolds) == 2 # we still just have 2 manifolds, those of the semiproduct
        # And it is a product of products
        @test Gl.manifold.manifolds[1] isa ProductManifold
        @test Gl.manifold.manifolds[2] isa ProductManifold
        @test Gl.manifold.manifolds[1].manifolds[1] === G1.manifold
        @test Gl.manifold.manifolds[2].manifolds[1] === G2.manifold
    end
    @testset "Defaults with SE(2)" begin

        struct TestLeftAction <: AbstractLeftGroupActionType end
        function LieGroups.apply!(A::GroupAction{TestLeftAction}, k, g, h)
            @assert is_point(A.manifold, h)
            @assert is_point(A.group, g)
            @assert is_point(A.manifold, g * h)
            return k .= g * h
        end

        struct TestRightAction <: AbstractRightGroupActionType end
        function LieGroups.apply!(A::GroupAction{TestRightAction}, k, g, h)
            @assert is_point(A.manifold, h)
            @assert is_point(A.group, g)
            @assert is_point(A.manifold, g * h)
            return k .= inv(A.group, g) * h
        end

        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        h1 = [0.1, 0.1]
        h2 = [0.0, 1.0]
        h3 = [0.0, 0]
        Xg1, Xg2, Xg3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 0.2 0.0], [0.0 0.0; 0.0 0.0]
        Xh1, Xh2, Xh3 = [0.0, 0.0], [0.0, -0.2], [1.0, 0.0]

        p1 = ArrayPartition(copy(g1), copy(h1))
        p2 = ArrayPartition(copy(g2), copy(h2))
        p3 = ArrayPartition(copy(g3), copy(h3))
        p4 = ArrayPartition(copy(g3), copy(h1))
        Y1 = ArrayPartition(copy(Xg1), copy(Xh1))
        Y2 = ArrayPartition(copy(Xg2), copy(Xh2))
        Y3 = ArrayPartition(copy(Xg3), copy(Xh3))
        properties = Dict(
            :Name => "The generic left semi-direct product group",
            :Points => [p1, p2, p3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, is_identity, inv, compose, show],
        )

        Gll = LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2), TestLeftAction())
        expectations_l = Dict(
            :repr => "LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2; field=ℝ), TestLeftAction())",
        )
        test_lie_group(Gll, properties, expectations_l)

        Glr = LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2), TestRightAction())
        expectations_l = Dict(
            :repr => "LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2; field=ℝ), TestRightAction())",
        )
        test_lie_group(Glr, properties, expectations_l)

        p1 = ArrayPartition(copy(h1), copy(g1))
        p2 = ArrayPartition(copy(h2), copy(g2))
        p3 = ArrayPartition(copy(h3), copy(g3))
        Y1 = ArrayPartition(copy(Xh1), copy(Xg1))
        Y2 = ArrayPartition(copy(Xh2), copy(Xg2))
        Y3 = ArrayPartition(copy(Xh3), copy(Xg3))
        properties = Dict(
            :Name => "The generic right semi-direct product group",
            :Points => [p1, p2, p3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, is_identity, inv, compose, show],
        )
        Grl = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), TestLeftAction())
        expectations_r = Dict(
            :repr => "RightSemidirectProductLieGroup(TranslationGroup(2; field=ℝ), SpecialOrthogonalGroup(2), TestLeftAction())",
        )
        test_lie_group(Grl, properties, expectations_r)

        Grr = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), TestRightAction())
        expectations_r = Dict(
            :repr => "RightSemidirectProductLieGroup(TranslationGroup(2; field=ℝ), SpecialOrthogonalGroup(2), TestRightAction())",
        )
        test_lie_group(Grr, properties, expectations_r)

        #sanity checks against the matrix representation
        rp = ArrayPartition(copy(h2), copy(g2))
        lp = ArrayPartition(copy(g2), copy(h2))

        # The left action versions correspond with the matrix representation
        p_mat = [g2 h2; 0 0 1]
        inv_p_mat = inv(p_mat)
        inv_g = inv_p_mat[1:2, 1:2]
        inv_h = inv_p_mat[1:2, 3]
        @test inv(Grl, rp) ≈ ArrayPartition(inv_h, inv_g)
        @test inv(Gll, lp) ≈ ArrayPartition(inv_g, inv_h)

    end

    @testset "Heisenberg group as a semidirect product" begin
        N = TranslationGroup(2)
        G = TranslationGroup(1)

        struct TestHeisenbergLeftAction <: AbstractLeftGroupActionType end
        function LieGroups.apply!(A::GroupAction{TestHeisenbergLeftAction}, k, g, h)
            @assert is_point(A.manifold, h)
            @assert is_point(A.group, g)
            y = g[1]
            x, z = h
            return k .= [x, z + x * y]
        end

        H = RightSemidirectProductLieGroup(N, G, TestHeisenbergLeftAction())

        g1 = ArrayPartition([1.0, 2.0], [3.0])
        g2 = ArrayPartition([2.0, 3.0], [-1.0])
        g3 = ArrayPartition([0.0, 0.0], [0.0])
        X1 = ArrayPartition([0.0, 0.0], [0.1])
        X2 = ArrayPartition([0.0, 0.0], [0.2])
        X3 = ArrayPartition([0.0, 0.0], [0.0])

        properties = Dict(
            :Name => "The Heisenberg group as a semidirect product",
            :Points => [g1, g2, g3],
            :Vectors => [X1, X2, X3],
            :Functions => [identity_element, is_identity, inv, compose, show],
        )
        expectations = Dict(
            :repr => "RightSemidirectProductLieGroup(TranslationGroup(2; field=ℝ), TranslationGroup(1; field=ℝ), TestHeisenbergLeftAction())",
        )
        test_lie_group(H, properties, expectations)

        #sanity checks against the matrix representation
        # ([x z], [y])
        # ([1.0, 2.0], [3.0])
        g1m = [1.0 1.0 2.0; 0.0 1.0 3.0; 0.0 0.0 1.0]
        inv_g1m = inv(g1m)
        inv_g1 = inv(H, g1)
        @test inv_g1.x[1] ≈ inv_g1m[1, 2:3]
        @test inv_g1.x[2][1] ≈ inv_g1m[2, 3]

        g2m = [1.0 2.0 3.0; 0.0 1.0 -1.0; 0.0 0.0 1.0]
        inv_g2m = inv(g2m)
        inv_g2 = inv(H, g2)
        @test inv_g2.x[1] ≈ inv_g2m[1, 2:3]
        @test inv_g2.x[2][1] ≈ inv_g2m[2, 3]
    end

end
