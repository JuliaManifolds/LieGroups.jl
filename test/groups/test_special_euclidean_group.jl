using LieGroups, Random, Test, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special Euclidean" begin
    𝔰 = sqrt(2)
    fcts = [compose, exp, hat, identity_element, inv, is_identity, log, rand, show, vee]
    #
    # ===
    # SE(2) in 4 variants, left semidirect
    @testset "SE(2)" begin
        G2f = SpecialEuclideanGroup(2)
        G2l = SpecialEuclideanGroup(2; parameter=:type)
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
        for G in [G2f, G2l], (pts, vec) in zip([gA, gM, gP, gQ], [XA, XM, XP, XQ])
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
    @testset "SE(4)" begin
        G4 = SpecialEuclideanGroup(4)
        g = identity_element(G4)
        h = copy(g)
        h[1:2, 1:2] .= 1 / sqrt(2) .* [1.0 1.0; -1.0 1.0]
        h[1, 5] = 1.0

        X = log(G4, g)
        @test norm(X) == 0
        @test isapprox(G4, exp(G4, X), g)

        Y = log(G4, h)
        @test is_point(LieAlgebra(G4), Y; error=:error)
        @test isapprox(G4, exp(G4, Y), h)
    end
    #
    #
    # Conversions
    @testset "Conversions between representations and incexinv" begin
        G2l = SpecialEuclideanGroup(2)
        g1 = 1 / 𝔰 .* [1.0 1.0 𝔰; -1.0 1.0 0.0; 0.0 0.0 𝔰]
        g2 = ArrayPartition(1 / 𝔰 * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
        g3 = SpecialEuclideanMatrixPoint(g1)
        g4 = SpecialEuclideanProductPoint(g2)
        X1 = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
        X2 = ArrayPartition([0.0 -0.23; 0.23 0.0], [0.0, 1.0])
        X3 = SpecialEuclideanMatrixTangentVector(X1)
        X4 = SpecialEuclideanProductTangentVector(X2)

        # Test also right semi to array
        G2r = SpecialEuclideanGroup(2; variant=:right)
        g2r = ArrayPartition([1.0, 0.0], 1 / 𝔰 * [1.0 1.0; -1.0 1.0])
        g4r = SpecialEuclideanProductPoint(g2r)
        X2r = ArrayPartition([0.0, 1.0], [0.0 -0.23; 0.23 0.0])
        X4r = SpecialEuclideanProductTangentVector(X2r)
        @testset "Conversions" begin
            # Convert everyone to everyone! (besides 1->3 since that would have been type piracy)
            @test convert(AbstractMatrix, g3) == g1
            @test convert(AbstractMatrix, g4) == g1
            @test convert(ArrayPartition, g3) == g2
            @test convert(ArrayPartition, g4) == g2

            @test convert(SpecialEuclideanMatrixPoint, g1) == g3
            @test convert(SpecialEuclideanMatrixPoint, g2) == g3
            @test convert(SpecialEuclideanMatrixPoint, g4) == g3
            @test LieGroups.internal_value(g3) == g1
            @test convert(SpecialEuclideanProductPoint, g1) == g4
            @test convert(SpecialEuclideanProductPoint, g2) == g4
            @test convert(SpecialEuclideanProductPoint, g3) == g4
            @test LieGroups.internal_value(g4) == g2

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

            @test convert(AbstractMatrix, g4r) == g1
            @test convert(AbstractMatrix, X4r) == X1
        end
        @testset "Element Access & norm" begin
            for i in [:Rotation, :Translation]
                @test g3[G2l, i] == g4[G2l, i]
                @test g3[G2l, i] == g4r[G2r, i]
                @test X3[LieAlgebra(G2l), i] == X4[LieAlgebra(G2l), i]
                @test X4r[LieAlgebra(G2r), i] == X3[LieAlgebra(G2l), i]
            end
            @test g3[G2l, :] == g4[G2l, :]
            @test g3[G2l, :] == g4r[G2r, :][[2, 1]]
            @test X3[LieAlgebra(G2l), :] == X4[LieAlgebra(G2l), :]
            @test X4r[LieAlgebra(G2r), :] == X3[LieAlgebra(G2l), :][[2, 1]]

            @test norm(LieAlgebra(G2l), X4) == norm(X1)
        end
    end
    @testset "Zero vector special types" begin
        G2l = SpecialEuclideanGroup(2)
        𝔤l = LieAlgebra(G2l)
        X1 = zero_vector(𝔤l, ArrayPartition{Float64})
        @test X1.x[1] == zeros(2, 2)
        @test X1.x[2] == zeros(2)
        G2r = SpecialEuclideanGroup(2; variant=:right)
        𝔤r = LieAlgebra(G2r)
        X2 = zero_vector(𝔤r, ArrayPartition{Float64})
        @test X2.x[1] == zeros(2)
        @test X2.x[2] == zeros(2, 2)
    end
    @testset "Test special cases in expected failing checks and internals" begin
        G = SpecialEuclideanGroup(2)
        𝔤 = LieAlgebra(G)
        Gr = SpecialEuclideanGroup(2; variant=:right)
        g1f = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]
        eG = Identity(G)
        eGm = identity_element(G)
        @test_throws DomainError is_point(G, g1f; error=:error)
        @test is_point(G, eG; error=:error)
        @test_throws DomainError is_point(G, Identity(Gr); error=:error)

        @test is_point(Gr, Identity(Gr); error=:error)
        @test_throws DomainError is_point(Gr, eG; error=:error)
        @test_throws DomainError is_point(G, g1f; error=:error)
        # non rot
        g2f = [2.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test_throws DomainError is_point(G, g2f; error=:error)
        # 2 errors
        g3f = [2.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 2.0]
        @test_throws CompositeManifoldError is_point(G, g3f; error=:error)
        g4f = zeros(2, 2)
        @test_throws DomainError is_point(G, g4f; error=:error)

        X1 = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
        @test is_vector(G, eG, X1; error=:error)
        @test is_vector(G, eGm, X1; error=:error)
        @test is_point(𝔤, X1; error=:error)
        # not affine
        X1f = [0.0 -0.23 0.0; 0.23 0.0 1.0; 0.0 0.0 1.0]
        @test_throws DomainError is_vector(G, eG, X1f; error=:error)
        @test_throws DomainError is_vector(G, eGm, X1f; error=:error)
        @test_throws DomainError is_point(𝔤, X1f; error=:error)
        # test internal fallback as well
        X1ft = SpecialEuclideanMatrixTangentVector(X1f)
        @test ManifoldsBase.check_size(𝔤, X1ft) isa DomainError
        # not skew
        X2f = [0.0 -0.63 0.0; 0.23 0.0 1.0; 0.0 0.0 0.0]
        @test_throws DomainError is_vector(G, eG, X2f; error=:error)
        @test_throws DomainError is_vector(G, eGm, X2f; error=:error)
        @test_throws DomainError is_point(𝔤, X2f; error=:error)
        # neither
        X3f = [0.0 -0.63 0.0; 0.23 0.0 1.0; 0.0 0.0 1.0]
        @test_throws CompositeManifoldError is_vector(G, eG, X3f; error=:error)
        @test_throws CompositeManifoldError is_vector(G, eGm, X3f; error=:error)
        @test_throws CompositeManifoldError is_point(𝔤, X3f; error=:error)
        # wrong size
        X4f = zeros(2, 2)
        @test_throws DomainError is_vector(G, eG, X4f; error=:error)
        @test_throws DomainError is_vector(G, eGm, X4f; error=:error)
        @test_throws DomainError is_point(𝔤, X4f; error=:error)

        # SE2 exp with zero vector
        @test is_identity(G, exp(G, zero_vector(𝔤)))
        G3 = SpecialEuclideanGroup(3)
        @test is_identity(G3, exp(G3, zero_vector(LieAlgebra(G3))))
    end
    @testset "Identity Subcomponents and rand special case" begin
        G = SpecialEuclideanGroup(2)
        X = zeros(3, 3)
        e = Identity(G)
        𝔤 = LieAlgebra(G)
        @test is_point(𝔤, rand!(G, X; vector_at=e))
        eO = Identity(SpecialOrthogonalGroup(2))
        eT = Identity(TranslationGroup(2))
        @test e[G, :Rotation] == eO
        @test e[G, :Translation] == eT
        @test e[G, :] == (eO, eT)
        Gr = SpecialEuclideanGroup(2; variant=:right)
        er = Identity(Gr)
        @test er[Gr, :Rotation] == eO
        @test er[Gr, :Translation] == eT
        @test er[G, :] == (eO, eT)
    end
end
