using LieGroups, ManifoldsBase, Random, Test, RecursiveArrayTools

using LieGroups: SpecialGalileanGroup

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

using StaticArrays
using LinearAlgebra

@testset "Special Galilean" begin
    𝔰 = sqrt(2)
    fcts = [
        compose,
        exp,
        get_vector,
        hat,
        identity_element,
        inner,
        inv,
        # is_flat,
        is_identity,
        log,
        norm,
        rand,
        show,
        vee,
        apply,
    ]

    @testset "SGal(3)" begin
        G3p = SpecialGalileanGroup(3)
        hL1 = ArrayPartition(
            ArrayPartition(
                1 / 𝔰 * [1.0 1.0 0.0; -1.0 1.0 0.0; 0.0 0.0 𝔰], [1 / 𝔰, 0.0, 0.0]
            ),
            ArrayPartition([0.0, 0.0, 0.0], [0.0]),
        )
        hL2 = ArrayPartition(
            ArrayPartition([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.1, 0.0], [0.1]),
        )
        hL3 = ArrayPartition(
            ArrayPartition([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [1.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.0, 1.0], [1.0]),
        )
        hL4 = identity_element(G3p, StaticArray)

        YL1 = ArrayPartition(
            ArrayPartition([0.0 -0.23 0.0; 0.23 0.0 0.0; 0.0 0.0 0.0], [0.0, 1.0, 0.0]),
            ArrayPartition([0.0, 0.0, 0.0], [1.0]),
        )
        YL2 = ArrayPartition(
            ArrayPartition([0.0 0.3 0.0; -0.3 0.0 0.0; 0.0 0.0 0.0], [1.0, 1.0, 0.0]),
            ArrayPartition([1.0, 0.3, 0.0], [0.1]),
        )
        YL3 = ArrayPartition(
            ArrayPartition([9.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0], [1.0, 0.0, 0.0]),
            ArrayPartition([0.0, 0.0, 0.3], [0.1]),
        )
        YL4 = zero(hL4)
        hP = [hL1, hL2, hL3, hL4]
        YP = [YL1, YL2, YL3, YL4]
        for G in [G3p], (pts, vec) in zip([hP], [YP])
            properties = Dict(
                :Name => "Special Galilean Group SGal(3)",
                :Points => pts,
                :Vectors => vec,
                :Rng => Random.MersenneTwister(),
                :Functions => fcts,
            )
            expectations = Dict(
                :atol => 1.0e-14,
                # :repr => "SpecialGalileanGroup(3)",
                # :is_flat => false
            )
            test_lie_group(G, properties, expectations)
        end
    end

    @testset "Test SGal(3) SArray" begin
        G = SpecialGalileanGroup(3)
        ε = identity_element(G, StaticArray)
        c = SA[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X = hat(LieAlgebra(G), c, typeof(ε))
        @test typeof(X) == typeof(ε)
        c1 = vee(LieAlgebra(G), X)
        @test c1 ≈ c
        @test typeof(c1) == typeof(c)
        g = exp(G, X)
        @test typeof(g) == typeof(ε)
        # pack c manually into a matrix to test against the matrix exponential
        Xmat = [
            0.0  -0.03 0.02 0.1 1.0;
            0.03  0.0 -0.01 0.2 2.0;
            -0.02  0.01 0.0  0.3 3.0;
            0.0   0.0  0.0  0.0 0.5;
            0.0   0.0  0.0  0.0 0.0
        ]
        gmat = exp(Xmat)
        @test g.x[1].x[1] ≈ gmat[1:3, 1:3]
        @test g.x[1].x[2] ≈ gmat[1:3, 4]
        @test g.x[2].x[1] ≈ gmat[1:3, 5]
        @test g.x[2].x[2][1] ≈ gmat[4, 5]
    end

    @testset "Test SGal(3) extras" begin
        G = SpecialGalileanGroup(3)
        ε = identity_element(G, StaticArray)
        c = SA[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X1 = hat(LieAlgebra(G), c, typeof(ε))
        g1 = exp(G, X1)
        inv_g1 = inv(G, g1)
        ε1 = compose(G, inv_g1, g1)

        ε = identity_element(G)
        c = [1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X2 = hat(LieAlgebra(G), c, typeof(ε))
        g2 = identity_element(G)
        exp!(G, g2, X2)
        inv_g2 = identity_element(G)
        inv!(G, inv_g2, g2)
        ε2 = identity_element(G)
        compose!(G, ε2, inv_g2, g2)

        @test isapprox(G, X1, X2)
        @test isapprox(G, g1, g2)
        @test isapprox(G, inv_g1, inv_g2)
        @test isapprox(G, ε1, ε, atol = 1.0e-14)
        @test isapprox(G, ε2, ε, atol = 1.0e-14)

    end
end
