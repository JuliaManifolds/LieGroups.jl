using LieGroups, Random, Test

using ManifoldsBase: ℝ, ℂ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special linear group" begin
    @testset "Real SL(2)" begin
        G = SpecialLinearGroup(2)
        g1, g2, g3 = [2.0 0.0; 0.0 0.5], [1.25 0.5; 0.5 1.0], [1.0 1.0; 3.0 4.0]
        X1, X2, X3 = [0.0 -1.0; 0.0 0.0], [0.0 0.5; 0.5 0.0], [-1.0 0.0; 0.0 1.0]
        properties = Dict(
            :Name => "The special linear group",
            :Points => [g1, g2, g3],
            :Vectors => [X1, X2, X3],
            :Rng => Random.MersenneTwister(),
            :Functions => [
                compose,
                conjugate,
                diff_conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                # exp,
                # hat, #check whether there is a nice formula for this, probably vectorize all but diag
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                # log,
                rand,
                show,
                #vee, # see hat
            ],
        )
        expectations = Dict(:repr => "SpecialLinearGroup(2, ℝ)", :atol => 1e-10)
        test_lie_group(G, properties, expectations)
    end
    @testset "Complex SL(2)" begin
        G = SpecialLinearGroup(2, ℂ)
        g1, g2, g3 = [1.0im 0.0; 0.0 -1.0im], [1.25 0.5; 0.5 1.0], [1.0 1.0; 3.0 4.0]
        X1, X2, X3 = [0.0 -1.0; 0.0 0.0], [0.0 0.5; 0.5 0.0], [-1.0 0.0; 0.0 1.0]
        properties = Dict(
            :Name => "The complex special linear group",
            :Points => [g1, g2, g3],
            :Vectors => [X1, X2, X3],
            :Rng => Random.MersenneTwister(),
            :Functions => [
                compose,
                conjugate,
                diff_conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                # exp,
                # hat
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                # log,
                rand,
                show,
                # vee,
            ],
        )
        expectations = Dict(:repr => "SpecialLinearGroup(2, ℂ)")
        test_lie_group(G, properties, expectations)
    end
end
