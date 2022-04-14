@testset "rotations" begin
    @testset "SO(2) group" begin
        # closure
        g1 = SO{2}(45/180*π)
        @test g1 isa SO{2}

        # identity element
        @test identity(IdentityRotationGroup()) == IdentityRotationGroup()
        @test identity(g1) == IdentityRotationGroup()
        @test identity(g1) * g1 == g1
        @test g1 * identity(g1) == g1

        # inverse element
        g2 = inv(g1)
        @test g2 == SO{2}(-45/180*π)
        @test g1 * g2 == identity(g1)
        @test g2 * g1 == identity(g1)

        # associativity
        g3 = SO{2}(90/180*π)
        @test (g1 * g2) * g3 == g1 * (g2 * g3)

        # inverse composition
        @test inv(g1 * g3) == inv(g3) * inv(g1)
        @test inv(g1 * g3) == SO{2}(-135/180*π)

        # matrix representation
        θ = 45/180*π
        @test Matrix(g1) == [cos(θ) -sin(θ);
                             sin(θ)  cos(θ)]
        
        @testset "left group action" begin
            # identity
            x = [1, 0]
            @test identity(g1) ⋉ x == x

            # compatibility
            @test g2 ⋉ (g1 ⋉ x) ≈ (g1 * g2) ⋉ x
        end
    end

    @testset "SO(3) group" begin
        # closure
        g1 = SO{3}(45/180*π, 3)
        @test g1 isa SO{3}

        # identity element
        @test identity(IdentityRotationGroup()) == IdentityRotationGroup()
        @test identity(g1) == IdentityRotationGroup()
        @test identity(g1) * g1 == g1
        @test g1 * identity(g1) == g1

        # inverse element
        g2 = inv(g1)
        @test g2 == SO{3}(-45/180*π, 3)
        @test g1 * g2 == identity(g1)
        @test g2 * g1 == identity(g1)

        # associativity
        g3 = SO{3}(90/180*π, 2)
        @test (g1 * g2) * g3 == g1 * (g2 * g3)

        # inverse composition
        @test inv(g1 * g3) == inv(g3) * inv(g1)
        @test inv(g1 * g3) == SO{3}(-90/180*π, 2) * SO{3}(-45/180*π, 3)

        # matrix representation
        θ = 45/180*π
        @test Matrix(g1) == [cos(θ) -sin(θ) 0;
                             sin(θ)  cos(θ) 0;
                                  0       0 1]
        
        g1 = SO{3}(45/180*π, 1) * SO{3}(30/180*π, 2)
        g2 = SO{3}(-30/180*π, 2) * SO{3}(-45/180*π, 1)
        @test g1 * g2 == IdentityRotationGroup()

        @testset "left group action" begin
            # identity
            x = [1, 0, 0]
            @test identity(g1) ⋉ x == x

            # compatibility
            @test g2 ⋉ (g1 ⋉ x) ≈ (g1 * g2) ⋉ x
        end
    end
end