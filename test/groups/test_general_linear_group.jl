using LieGroups, Random, Test

using ManifoldsBase: ℂ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "GL(2)" begin
    G = GeneralLinearGroup(2)
    g1, g2, g3 = [2.0 0.0; 0.0 1.0], [1.0 0.5; 0.5 1.0], [1.0 2.0; 3.0 4.0]
    X1, X2, X3 = [-1.0 0.0; 0.0 0.0], [0.0 0.5; 0.5 0.0], [1.0 2.0; 3.0 4.0]
    properties = Dict(
        :Name => "The general linear group",
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
            exp,
            # hat, # requires a fix in Manifolds.jl to have an ONB on invertible matrices
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            #rand, # requires a fix in Manifolds.jl to have rand on invertible matrices
            show,
            #vee, # requires a fix in Manifolds.jl to have an ONB on invertible matrices
        ],
    )
    expectations = Dict(
        :repr => "GeneralLinearGroup(2; field=ℝ)", :lie_bracket => X1 * X2 - X2 * X1
    )
    test_lie_group(G, properties, expectations)

    @test_broken is_point(G, Identity(G); error=:error)
end

@testset "GL(1, 𝔽) special cases" begin
    @testset "real" begin
        G = GeneralLinearGroup(1)
        e = Identity(G)
        p = 3.0 * ones(1, 1)
        X = 1.0 * ones(1, 1)
        @test exp(G, p, X) ≈ p * exp(X)' * exp(X - X')
        q = exp(G, p, X)
        Y = log(G, p, q)
        @test Y ≈ X
        @test exp(G, e, X) ≈ exp(X)
        @test log(G, e, exp(X)) ≈ X
        log(G, e, Identity(G)) == zeros(1, 1) # Matrix to matrix
    end
    @testset "complex" begin
        G = GeneralLinearGroup(1; field=ℂ)
        e = Identity(G)
        p = (1 + im) * ones(1, 1)
        X = (1 - im) * ones(1, 1)
        @test exp(G, p, X) ≈ p * exp(X)' * exp(X - X')
        q = exp(G, p, X)
        Y = log(G, p, q)
        @test Y ≈ X
        @test exp(G, e, X) ≈ exp(X)
        @test log(G, e, exp(X)) ≈ X
    end
end
