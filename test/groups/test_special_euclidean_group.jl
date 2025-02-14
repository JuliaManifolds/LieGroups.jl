using LieGroups, Random, Test, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special Euclidean" begin
    𝔰 = sqrt(2)
    fcts = [
        # adjoint,
        compose,
        # conjugate,
        # diff_inv,
        # diff_left_compose,
        # diff_right_compose,
        exp,
        hat,
        identity_element,
        inv,
        # inv_left_compose,
        # inv_right_compose,
        is_identity,
        # lie_bracket,
        log,
        rand,
        show,
        vee,
    ]
    #
    # ===
    # SE(2) in 4 variants, left semidirect
    @testset "SE(2)" begin
        G2f = SpecialEuclideanGroup(2)
        G2p = SpecialEuclideanGroup(2; parameter=:type)
        g1 = 1 / 𝔰 .* [1.0 1.0 𝔰; -1.0 1.0 0.0; 0.0 0.0 𝔰]
        g2 = [0.0 -1.0 0.0; 1.0 0.0 1.0; 0.0 0.0 1.0]
        g3 = [1.0 0.0 1.0; 0.0 1.0 1.0; 0.0 0.0 1.0]
        gL1 = ArrayPartition(1 / 𝔰 * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
        gL2 = ArrayPartition([0.0 -1.0; 1.0 0.0], [0.0, 1.0])
        gL3 = ArrayPartition([1.0 0.0; 0.0 1.0], [1.0, 1.0])
        X1 = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
        X2 = [0.0 0.30 1.0; -0.30 0.0 1.0; 0.0 0.0 0.0]
        X3 = [0.0 0.1 1.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
        XL1 = ArrayPartition([0.0 -0.23; 0.23 0.0], [0.0, 1.0])
        XL2 = ArrayPartition([0.0 0.30; -0.30 0.0], [1.0, 1.0])
        XL3 = ArrayPartition([9.0 0.1; -0.1 0.0], [1.0, 0.0])
        gA = [g1, g2, g3]
        gM = SpecialEuclideanMatrixPoint.(gA)
        gP = [gL1, gL2, gL3]
        gQ = SpecialEuclideanProductPoint.(gP)
        XA = [X1, X2, X3]
        XM = SpecialEuclideanMatrixTangentVector.(XA)
        XP = [XL1, XL2, XL3]
        XQ = SpecialEuclideanProductTangentVector.(XP)
        for G in [G2f, G2p], (pts, vec) in zip([gA, gM, gP, gQ], [XA, XM, XP, XQ])
            properties = Dict(
                :Name => "The special Euclidean group ($G, $(eltype(pts)))",
                :Points => pts,
                :Vectors => vec,
                :Rng => Random.MersenneTwister(),
                :Functions => fcts,
            )
            expectations = Dict(:repr => "SpecialEuclideanGroup(2)", :atol => 1e-14)
            test_lie_group(G, properties, expectations)
        end
        #
        # Right variant – exchange product cases
        G2r = TranslationGroup(2) ⋊ SpecialOrthogonalGroup(2)
        G2rt =
            TranslationGroup(2; parameter=:type) ⋊
            SpecialOrthogonalGroup(2; parameter=:type)
        gR1 = ArrayPartition([1.0, 0.0], 1 / 𝔰 * [1.0 1.0; -1.0 1.0])
        gR2 = ArrayPartition([0.0, 1.0], [0.0 -1.0; 1.0 0.0])
        gR3 = ArrayPartition([1.0, 1.0], [1.0 0.0; 0.0 1.0])
        XR1 = ArrayPartition([0.0, 1.0], [0.0 -0.23; 0.23 0.0])
        XR2 = ArrayPartition([1.0, 1.0], [0.0 0.3; -0.3 0.0])
        XR3 = ArrayPartition([1.0, 0.0], [9.0 -0.1; 0.1 0.0])
        gS = [gR1, gR2, gR3]
        gT = SpecialEuclideanProductPoint.(gS)
        XS = [XR1, XR2, XR3]
        XT = SpecialEuclideanProductTangentVector.(XS)

        for G in [G2r, G2rt], (pts, vec) in zip([gA, gM, gS, gT], [XA, XM, XS, XT])
            properties = Dict(
                :Name => "The special Euclidean group ($G, $(eltype(pts)))",
                :Points => pts,
                :Vectors => vec,
                :Rng => Random.MersenneTwister(),
                :Functions => fcts,
            )
            expectations = Dict(
                :repr => "SpecialEuclideanGroup(2; variant=:right)", :atol => 1e-14
            )
            test_lie_group(G, properties, expectations)
        end
    end
    #
    #
    # SE(3)
    @testset "SE(3)" begin
        G3f = SpecialEuclideanGroup(3)
        G3p = SpecialEuclideanGroup(3; parameter=:type)
        h1 = 1 / 𝔰 .* [1.0 1.0 0.0 1.0; -1.0 1.0 0.0 0.0; 0.0 0.0 𝔰 0.0; 0.0 0.0 0.0 𝔰]
        h2 = [0.0 -1.0 0.0 0.0; 1.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
        h3 = [1.0 0.0 0.0 1.0; 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
        hL1 = ArrayPartition(
            1 / 𝔰 * [1.0 1.0 0.0; -1.0 1.0 0.0; 0.0 0.0 𝔰], [1 / 𝔰, 0.0, 0.0]
        )
        hL2 = ArrayPartition([0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 1.0], [0.0, 1.0, 0.0])
        hL3 = ArrayPartition([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [1.0, 1.0, 0.0])
        Y1 = [0.0 -0.23 0.0 0.0; 0.23 0.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        Y2 = [0.0 0.30 0.0 1.0; -0.30 0.0 0.0 1.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        Y3 = [0.0 0.1 0.0 1.0; -0.1 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
        YL1 = ArrayPartition([0.0 -0.23 0.0; 0.23 0.0 0.0; 0.0 0.0 0.0], [0.0, 1.0, 0.0])
        YL2 = ArrayPartition([0.0 0.30 0.0; -0.30 0.0 0.0; 0.0 0.0 0.0], [1.0, 1.0, 0.0])
        YL3 = ArrayPartition([9.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0], [1.0, 0.0, 0.0])
        hA = [h1, h2, h3]
        hM = SpecialEuclideanMatrixPoint.(hA)
        hP = [hL1, hL2, hL3]
        hQ = SpecialEuclideanProductPoint.(hP)
        YA = [Y1, Y2, Y3]
        YM = SpecialEuclideanMatrixTangentVector.(YA)
        YP = [YL1, YL2, YL3]
        YQ = SpecialEuclideanProductTangentVector.(YP)
        for G in [G3f, G3p], (pts, vec) in zip([hA, hM, hP, hQ], [YA, YM, YP, YQ])
            properties = Dict(
                :Name => "The special Euclidean group ($G, $(eltype(pts)))",
                :Points => pts,
                :Vectors => vec,
                :Rng => Random.MersenneTwister(),
                :Functions => fcts,
            )
            expectations = Dict(:repr => "SpecialEuclideanGroup(3)", :atol => 1e-14)
            test_lie_group(G, properties, expectations)
        end
    end
    #
    #
    # SE(4)

    #
    #
    # Conversions
    @testset "Conversions between representations" begin
        G2p = SpecialEuclideanGroup(2)
        g1 = 1 / 𝔰 .* [1.0 1.0 𝔰; -1.0 1.0 0.0; 0.0 0.0 𝔰]
        g2 = ArrayPartition(1 / 𝔰 * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
        g3 = SpecialEuclideanMatrixPoint(g1)
        g4 = SpecialEuclideanProductPoint(g2)
        X1 = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
        X2 = ArrayPartition([0.0 -0.23; 0.23 0.0], [0.0, 1.0])
        X3 = SpecialEuclideanMatrixTangentVector(X1)
        X4 = SpecialEuclideanProductTangentVector(X2)
        # Convert everyone to everyone! (besides 1->3 since that would have been type piracy)
        @test convert(AbstractMatrix, g3) == g1
        @test convert(AbstractMatrix, g4) == g1
        @test convert(ArrayPartition, g3) == g2
        @test convert(ArrayPartition, g4) == g2

        @test convert(SpecialEuclideanMatrixPoint, g1) == g3
        @test convert(SpecialEuclideanMatrixPoint, g2) == g3
        @test convert(SpecialEuclideanMatrixPoint, g4) == g3
        @test convert(SpecialEuclideanProductPoint, g1) == g4
        @test convert(SpecialEuclideanProductPoint, g2) == g4
        @test convert(SpecialEuclideanProductPoint, g3) == g4

        @test convert(AbstractMatrix, X3) == X1
        @test convert(AbstractMatrix, X4) == X1
        @test convert(ArrayPartition, X3) == X2
        @test convert(ArrayPartition, X4) == X2

        @test convert(SpecialEuclideanMatrixTangentVector, X1) == X3
        @test convert(SpecialEuclideanMatrixTangentVector, X2) == X3
        @test convert(SpecialEuclideanMatrixTangentVector, X4) == X3
        @test convert(SpecialEuclideanProductTangentVector, X1) == X4
        @test convert(SpecialEuclideanProductTangentVector, X2) == X4
        @test convert(SpecialEuclideanProductTangentVector, X3) == X4
    end
end
