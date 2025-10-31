using LieGroups, ManifoldsBase, Random, Test, RecursiveArrayTools

using LieGroups: MetricLieGroup

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "MetricLieGroup: A metric decorator for LieGroups" begin

    @testset "Pass through with a dummy metric" begin
        SO2 = SpecialOrthogonalGroup(2)
        G = MetricLieGroup(SO2, LieGroupsTestSuite.DummyMetric())
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
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                identity_element,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                rand,
                show,
            ],
        )
        expectations = Dict(:repr => "MetricLieGroup(SpecialOrthogonalGroup(2), LieGroupsTestSuite.DummyMetric())")
        test_lie_group(G, properties, expectations)
    end
end
