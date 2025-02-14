using LieGroups, Random, Test, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special Euclidean" begin
    #
    # ===
    # SE(2) in 4 variants, left semidirect
    G2f = SpecialEuclideanGroup(2)
    G2p = SpecialEuclideanGroup(2; parameter=:type)
    g1 = 1 / sqrt(2) .* [1.0 1.0 sqrt(2); -1.0 1.0 0.0; 0.0 0.0 sqrt(2)]
    g2 = [0.0 -1.0 0.0; 1.0 0.0 1.0; 0.0 0.0 1.0]
    g3 = [1.0 0.0 1.0; 0.0 1.0 1.0; 0.0 0.0 1.0]
    gL1 = ArrayPartition(1 / sqrt(2) * [1.0 1.0; -1.0 1.0], [1.0, 0.0])
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
    G2rt = TranslationGroup(2; parameter=:type) ⋊ SpecialOrthogonalGroup(2; parameter=:type)
    gR1 = ArrayPartition([1.0, 0.0], 1 / sqrt(2) * [1.0 1.0; -1.0 1.0])
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
    #
    #
    # SE(3)

    #
    #
    # SE(4)
end
