using LieGroups, Random, Test
using LinearAlgebra: I
using ManifoldsBase
using StaticArrays

@testset "The Orthogonal Group" begin
    G = OrthogonalGroup(2)
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
            exp,
            hat,
            identity_element,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            rand,
            show,
            vee,
        ],
    )
    expectations = Dict(:repr => "OrthogonalGroup(2)")
    LieGroups.Test.test_lie_group(G, properties, expectations)

    @testset "StaticArrays specializations for O(2)" begin
        ps = SMatrix{2, 2, Float64}(I)
        Xs = SA[0.0 2.0; -2.0 0.0]
        qs = exp(G, ps, Xs)

        q_ref = SA[
            -0.9513631281258474 0.30807174236304485
            -0.30807174236304485 -0.9513631281258474
        ]
        @test qs isa SMatrix{2, 2, Float64}
        @test qs === q_ref

        qs = exp(G, Xs)
        @test qs isa SMatrix{2, 2, Float64}
        @test qs === q_ref

        qs = ManifoldsBase.exp_fused(G, ps, Xs, 1.5)
        @test qs isa SMatrix{2, 2, Float64}
        @test qs === SA[
            -0.45266185729235203 -0.8916822544789362
            0.8916822544789362 -0.45266185729235203
        ]

        @test get_coordinates(G, ps, Xs, DefaultLieAlgebraOrthogonalBasis()) isa SVector{1}
        @test get_vector(G, ps, SA[1.0], DefaultLieAlgebraOrthogonalBasis()) ===
            SA[0.0 -1.0; 1.0 0.0]
    end

    #
    #
    # O(3)
    H = OrthogonalGroup(3)
    h1 = LieGroups.Test.rotation_matrix(3, 2, 1, π / 4)
    h2 = LieGroups.Test.rotation_matrix(3, 2, 1, π / 8) * LieGroups.Test.rotation_matrix(3, 3, 1, π / 4)
    h3 = LieGroups.Test.rotation_matrix(3, 3, 1, π / 4) * LieGroups.Test.rotation_matrix(3, 3, 2, π / 8)
    Y1 = [0.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
    Y2 = [0.0 0.0 -0.2; 0.0 0.0 0.0; 0.2 0.0 0.0]
    Y3 = [0.0 0.3 0.0; -0.3 0.0 0.4; 0.0 -0.4 0.0]
    # Test only specialized functions here
    properties2 = Dict(
        :Name => "The orthogonal group O(3) – specialised funcions",
        :Points => [h1, h2, h3],
        :Vectors => [Y1, Y2, Y3],
        :Functions => [exp, hat, log, show, vee],
    )
    expectations2 = Dict(:repr => "OrthogonalGroup(3)", :atols => Dict(exp => 1.0e-15))
    LieGroups.Test.test_lie_group(H, properties2, expectations2)
    @testset "O(3) special cases" begin
        @test is_identity(H, exp(H, zeros(3, 3)))
    end

    @testset "StaticArrays specializations for O(3)" begin
        ps = SMatrix{3, 3, Float64}(I)
        Xs = SA[0.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
        qs = exp(H, ps, Xs)

        q_ref = SA[
            0.9950041652780258 0.09983341664682815 0.0
            -0.09983341664682815 0.9950041652780258 0.0
            0.0 0.0 1.0
        ]
        @test qs isa SMatrix{3, 3, Float64}
        @test qs ≈ q_ref

        qs = exp(H, Xs)
        @test qs isa SMatrix{3, 3, Float64}
        @test qs ≈ q_ref

        @test get_coordinates(H, ps, Xs, DefaultLieAlgebraOrthogonalBasis()) isa SVector{3}
        @test get_vector(H, ps, SA[1.0, 2.0, 3.0], DefaultLieAlgebraOrthogonalBasis()) ===
            SA[0.0 -3.0 2.0; 3.0 0.0 -1.0; -2.0 1.0 0.0]
    end
    #
    #
    # O(4)
    J = OrthogonalGroup(4)
    j1 = LieGroups.Test.rotation_matrix(4, 3, 1, π / 4)
    j2 = LieGroups.Test.rotation_matrix(4, 4, 1, π / 8) * LieGroups.Test.rotation_matrix(4, 3, 2, π / 4)
    j3 = LieGroups.Test.rotation_matrix(4, 4, 1, π / 4) * LieGroups.Test.rotation_matrix(4, 4, 2, π / 8)
    Z1 = [0.0 0.1 0.0 0.0; -0.1 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    Z2 = [0.0 0.0 0.2 0.0; 0.0 0.0 0.0 0.0; -0.2 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]
    Z3 = [0.0 0.1 0.0 0.3; 0.0 0.0 -0.4 0.0; 0.0 0.4 0.0 0.0; -0.3 0.0 0.0 0.0]
    # Test only specialized functions here
    properties3 = Dict(
        :Name => "The orthogonal group O(4) – specialised funcions",
        :Points => [j1, j2, j3],
        :Vectors => [Z1, Z2, Z3],
        :Functions => [exp, hat, log, show, vee],
    )
    expectations3 = Dict(:repr => "OrthogonalGroup(4)", :atols => Dict(exp => 1.0e-15))
    LieGroups.Test.test_lie_group(J, properties3, expectations3)
    @testset "𝔬(4) edge cases" begin
        𝔧 = LieAlgebra(J)
        e = Identity(J)
        X = zero_vector(𝔧)
        d = vee(𝔧, X)
        p = exp(J, X)
        for c in LieGroups.Test.𝔰𝔬4_edges_cases_explog
            @testset "$c on $J" begin
                hat!(𝔧, X, c)
                vee!(𝔧, d, X)
                @test isapprox(c, d)
                l = exp(X)
                exp!(J, p, X)
                @test l ≈ p
                p2 = exp(J, log(J, p))
                @test isapprox(J, p, p2; atol = 1.0e-6)
            end
        end
    end
end
