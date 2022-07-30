@testset "rotations" begin
    @testset "SO(2) group" begin
        θ = 45/180*π
        ϕ = 90/180*π

        alg1 = so{2}([θ])
        alg2 = so{2}([-θ])

        # closure
        @test alg1 + alg2 isa so{2}

        # identity element
        @test identity(alg1) == so{2}([0.])

        # inverse element
        @test inv(alg1) == alg2

        # associativity
        alg3 = so{2}([ϕ])
        @test (alg1 + alg2) + alg3 == alg1 + (alg2 + alg3)


        g1 = SO{2}(
            [cos(θ) -sin(θ);
             sin(θ)  cos(θ)]
        )
        g2 = SO{2}(
            [cos(-θ) -sin(-θ);
             sin(-θ)  cos(-θ)]
        )

        # closure
        @test g1 * g2 isa SO{2}

        # identity element
        @test identity(g1) == SO{2}(I(2))
        @test identity(g1) * g1 == g1
        @test g1 * identity(g1) == g1

        # inverse element
        @test inv(g1) == g2
        @test g1 * g2 == identity(g1)
        @test g2 * g1 == identity(g1)

        # associativity
        g3 = SO{2}(
            [cos(ϕ) -sin(ϕ);
             sin(ϕ)  cos(ϕ)]
        )
        @test (g1 * g2) * g3 ≈ g1 * (g2 * g3)

        # inverse composition
        @test inv(g1 * g3) == inv(g3) * inv(g1)
        @test inv(g1 * g3) == SO{2}(
            [cos(-θ-ϕ) -sin(-θ-ϕ);
             sin(-θ-ϕ)  cos(-θ-ϕ)]
        )

        # matrix representation
        @test Matrix(g1) == [cos(θ) -sin(θ);
                             sin(θ)  cos(θ)]
        
        @testset "left group action" begin
            # identity
            x = [1.5, 3.0]
            @test identity(g1) ⋉ x == x

            # compatibility
            @test g2 ⋉ (g1 ⋉ x) ≈ (g1 * g2) ⋉ x
        end

        @testset "exponential map" begin
            @test exp(alg1) ≈ g1
            @test exp(∧(alg1)) ≈ g1
            @test log(g1) ≈ alg1

            @test Matrix(∧(alg1)) == [0 -θ; θ 0]
            @test Vector(∨(alg1)) == [θ]
        end
    end

    @testset "SO(3) group" begin
        θ = 45/180*π
        ϕ = 90/180*π

        alg1 = so{3}([0., 0., θ])
        alg2 = so{3}([0., 0., -θ])

        # closure
        @test alg1 + alg2 isa so{3}

        # identity element
        @test identity(alg1) == so{3}([0., 0., 0.])

        # inverse element
        @test inv(alg1) == alg2

        # associativity
        alg3 = so{3}([ϕ, 0., θ])
        @test (alg1 + alg2) + alg3 == alg1 + (alg2 + alg3)


        g1 = SO{3}(
            [cos(θ) -sin(θ) 0;
             sin(θ)  cos(θ) 0;
                  0       0 1]
        )
        g2 = SO{3}(
            [cos(-θ) -sin(-θ) 0;
             sin(-θ)  cos(-θ) 0;
                   0        0 1]
        )

        # closure
        @test g1 * g2 isa SO{3}

        # identity element
        @test identity(g1) == SO{3}(I(3))
        @test identity(g1) * g1 == g1
        @test g1 * identity(g1) == g1

        # inverse element
        @test inv(g1) == g2
        @test g1 * g2 == identity(g1)
        @test g2 * g1 == identity(g1)

        # associativity
        g3 = SO{3}(
            [1     0        0;
             0 cos(ϕ) -sin(ϕ);
             0 sin(ϕ)  cos(ϕ)]
        )
        @test (g1 * g2) * g3 == g1 * (g2 * g3)

        # inverse composition
        @test inv(g1 * g3) == inv(g3) * inv(g1)
        @test inv(g1 * g3) == SO{3}(inv(Matrix(g3))) * SO{3}(inv(Matrix(g1)))

        # matrix representation
        @test Matrix(g1) == [cos(θ) -sin(θ) 0;
                             sin(θ)  cos(θ) 0;
                                  0       0 1]

        @testset "left group action" begin
            # identity
            x = [1., 2., 3.]
            @test identity(g1) ⋉ x == x

            # compatibility
            @test g2 ⋉ (g1 ⋉ x) ≈ (g1 * g2) ⋉ x
        end

        @testset "exponential map" begin
            @test exp(alg1) ≈ g1
            @test exp(∧(alg1)) ≈ g1
            @test log(g1) ≈ alg1

            @test Matrix(∧(alg1)) == [0 -θ 0;
                                      θ  0 0;
                                      0  0 0]
            @test Vector(∨(alg1)) == [0., 0., θ]
        end
    end
end