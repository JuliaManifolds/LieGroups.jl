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
        Xⁱ = SA[1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5]
        X = hat(LieAlgebra(G), Xⁱ, typeof(ε))
        @test typeof(X) == typeof(ε)
        Xⁱ1 = vee(LieAlgebra(G), X)
        @test Xⁱ1 ≈ Xⁱ
        @test typeof(Xⁱ1) == typeof(Xⁱ)
        @test typeof(exp(G, X)) == typeof(ε)
    end

end
