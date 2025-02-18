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
        G = SpecialOrthogonalGroup(2)
        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        X1, X2, X3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 0.2 0.0], [0.0 0.0; 0.0 0.0]

        Gl = LeftSemidirectProductLieGroup(G, G, LeftGroupOperationAction())
        Gr = RightSemidirectProductLieGroup(G, G, RightGroupOperationAction())

        h1 = ArrayPartition(copy(g1), copy(g2))
        h2 = ArrayPartition(copy(g2), copy(g1))
        h3 = ArrayPartition(copy(g3), copy(g3))
        Y1 = ArrayPartition(copy(X1), copy(X2))
        Y2 = ArrayPartition(copy(X2), copy(X1))
        Y3 = ArrayPartition(copy(X3), copy(X3))
        properties = Dict(
            :Name => "The generic left semidirect product group",
            :Points => [h1, h2, h3],
            :Vectors => [Y1, Y2, Y3],
            :Functions => [identity_element, inv, show],
        )
        expectations_l = Dict(
            :repr => "LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), SpecialOrthogonalGroup(2), LeftGroupOperationAction())",
        )
        test_lie_group(Gl, properties, expectations_l)
        expectations_r = Dict(
            :repr => "RightSemidirectProductLieGroup(SpecialOrthogonalGroup(2), SpecialOrthogonalGroup(2), RightGroupOperationAction())",
        )
        test_lie_group(Gr, properties, expectations_r)
    end
end
