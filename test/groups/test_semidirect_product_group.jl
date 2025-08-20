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

        struct TestLeftAction <: AbstractLeftGroupActionType end

        function LieGroups.apply!(A::GroupAction{TestLeftAction}, k, g, h)
            return k .= g * h
        end
        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        h1 = [0.0, 0]
        h2 = [0.1, 0.1]
        h3 = [0.0, 1.0]
        Xg1, Xg2, Xg3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 0.2 0.0], [0.0 0.0; 0.0 0.0]
        Xh1, Xh2, Xh3 = [0.0, 0.0], [0.0, -0.2], [1.0, 0.0]

        Gl = LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2), TestLeftAction())
        Gr = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), TestLeftAction())

        p1 = ArrayPartition(copy(g1), copy(h1))
        p2 = ArrayPartition(copy(g2), copy(h2))
        p3 = ArrayPartition(copy(g3), copy(h3))
        Y1 = ArrayPartition(copy(Xg1), copy(Xg1))
        Y2 = ArrayPartition(copy(Xg2), copy(Xg2))
        Y3 = ArrayPartition(copy(Xg3), copy(Xg3))
        properties = Dict(
            :Name => "The generic left semi-direct product group",
            :Points => [p1, p2, p3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, inv, compose, show],
        )
        expectations_l = Dict(
            :repr => "LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2; field=ℝ), TestLeftAction())",
        )
        test_lie_group(Gl, properties, expectations_l)

        p1 = ArrayPartition(copy(h1), copy(g1))
        p2 = ArrayPartition(copy(h2), copy(g2))
        p3 = ArrayPartition(copy(h3), copy(g3))
        Y1 = ArrayPartition(copy(Xg1), copy(Xg1))
        Y2 = ArrayPartition(copy(Xg2), copy(Xg2))
        Y3 = ArrayPartition(copy(Xg3), copy(Xg3))
        properties = Dict(
            :Name => "The generic right semi-direct product group",
            :Points => [p1, p2, p3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, inv, compose, show],
        )
        expectations_r = Dict(
            :repr => "RightSemidirectProductLieGroup(TranslationGroup(2; field=ℝ), SpecialOrthogonalGroup(2), TestLeftAction())",
        )
        test_lie_group(Gr, properties, expectations_r)
    end

    # @testset "Combining ⋊, ⋉, and ×" begin
    #     fcts = [
    #         compose,
    #         exp,
    #         get_vector,
    #         hat,
    #         identity_element,
    #         inner,
    #         inv,
    #         is_identity,
    #         lie_bracket,
    #         log,
    #         norm,
    #         rand,
    #         show,
    #         vee,
    #     ]

    #     G = TranslationGroup(2) × (TranslationGroup(2) ⋊ SpecialOrthogonalGroup(2)) × TranslationGroup(2) × (SpecialOrthogonalGroup(2) ⋉ TranslationGroup(2))
    #     ε = identity_element(G)
    #     p = rand(G)
    #     pts = [p, ε]
    #     vec = [rand(G; vector_at = p), zero_vector(G, ε)]

    #     properties = Dict(
    #         :Name => "Test combining ⋊, ⋉, and ×",
    #         :Points => pts,
    #         :Vectors => vec,
    #         :Functions => fcts,
    #     )

    #     expectations = Dict(
    #         :atol => 1.0e-14,
    #     )

    #     test_lie_group(G, properties, expectations)

    #     G2 = ProductLieGroup(
    #         TranslationGroup(2),
    #         TranslationGroup(2) ⋊ SpecialOrthogonalGroup(2),
    #         TranslationGroup(2),
    #         SpecialOrthogonalGroup(2) ⋉ TranslationGroup(2)
    #     )
    #     @test G2 == G

    #     G3 = ProductLieGroup(TranslationGroup(2))
    #     @test G3.manifold isa ProductManifold

    #     #test that the different combinations of products groups are equivalent
    #     PG = TranslationGroup(1) × SpecialEuclideanGroup(2) × TranslationGroup(2) × CircleGroup()
    #     PG1 = ProductLieGroup(
    #         TranslationGroup(1),
    #         SpecialEuclideanGroup(2) × TranslationGroup(2) × CircleGroup()
    #     )
    #     @test PG1 == PG
    #     PG2 = ProductLieGroup(
    #         TranslationGroup(1) × SpecialEuclideanGroup(2),
    #         TranslationGroup(2) × CircleGroup()
    #     )
    #     @test PG2 == PG
    #     PG3 = ProductLieGroup(
    #         TranslationGroup(1) × SpecialEuclideanGroup(2) × TranslationGroup(2),
    #         CircleGroup()
    #     )
    #     @test PG3 == PG
    #     PG4 = ProductLieGroup(PG)
    #     @test PG4 == PG
    # end
end
