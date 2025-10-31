using LieGroups, ManifoldsBase, Random, Test, RecursiveArrayTools

using LieGroups: MetricLieGroup
using Manifolds

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

struct CustomTranslationMetric <: ManifoldsBase.AbstractMetric end

function ManifoldsBase.exp!(
        ::MetricLieGroup{ℝ, AdditionGroupOperation, <:Euclidean, <:TranslationGroup, CustomTranslationMetric},
        q,
        p,
        X,
    )
    return q .= p .+ 2 .* X
end

@testset "MetricLieGroup: A metric decorator for LieGroups" begin

    @testset "Pass through with a dummy metric" begin
        SO2 = SpecialOrthogonalGroup(2)
        @test metric(SO2) === DefaultMetric()
        G = MetricLieGroup(SO2, LieGroupsTestSuite.DummyMetric())
        @test metric(G) === LieGroupsTestSuite.DummyMetric()
        @test base_lie_group(G) === SO2

        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        X1, X2, X3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 2.0 0.0], [0.0 0.0; 0.0 0.0]
        properties = Dict(
            :Name => "The orthogonal group O(2)",
            :Points => [g1, g2, g3],
            :Vectors => [X1, X2, X3],
            :Rng => Random.MersenneTwister(),
            :Functions => [
                adjoint,
                compose,
                conjugate,
                diff_conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                identity_element,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                rand,
                show,
            ],
        )
        expectations = Dict(:repr => "MetricLieGroup(SpecialOrthogonalGroup(2), LieGroupsTestSuite.DummyMetric())")
        test_lie_group(G, properties, expectations)
    end
    @testset "Passthrough for index access" begin
        SE2 = SpecialEuclideanGroup(2)
        se2 = LieAlgebra(SE2)
        G = MetricLieGroup(SE2, LieGroupsTestSuite.DummyMetric())
        𝔤 = LieAlgebra(G)
        g = [1.0 0.0 2.0; 0.0 1.0 3.0; 0.0 0.0 1.0]
        X = [0.0 -0.1 0.5; 0.1 0.0 1.0; 0.0 0.0 0.0]
        @test g[SE2, :Rotation] === g[G, :Rotation]
        @test g[SE2, :Translation] === g[G, :Translation]
        @test X[se2, :Rotation] === X[𝔤, :Rotation]
        @test X[se2, :Translation] === X[𝔤, :Translation]
        gT = SpecialEuclideanMatrixPoint(g)
        XT = SpecialEuclideanMatrixTangentVector(X)
        @test gT[SE2, :Rotation] == gT[G, :Rotation]
        @test gT[SE2, :Translation] == gT[G, :Translation]
        @test XT[se2, :Rotation] == XT[𝔤, :Rotation]
        @test XT[se2, :Translation] == XT[𝔤, :Translation]
    end
    @testset "Testing that correct defaults for metric functions are used" begin
        T2 = TranslationGroup(2)
        G = MetricLieGroup(T2, CustomTranslationMetric())
        p = [0.0, 2.0]
        X = [1.0, -2.0]
        @test exp(G, p, X) == [2.0, -2.0]
    end
end
