@testset "rigid motions" begin
    @testset "SE(2) group" begin
        T = Float64
        δ = [1., 2]
        θ = 45/180*π
        ϕ = 90/180*π
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]

        alg1 = se{2}([δ..., θ])
        alg2 = se{2}([(-δ)..., -θ])

        # closure
        @test alg1 + alg2 isa se{2}

        # identity element
        @test identity(alg1) == se{2}([0., 0., 0.])

        # inverse element
        @test inv(alg1) == alg2

        # associativity
        alg3 = se{2}([3., -2., ϕ])
        @test (alg1 + alg2) + alg3 ≈ alg1 + (alg2 + alg3)


        g1 = SE{2}(
            [                    R δ;
            zeros(T, 1, size(R,2)) 1]
        )
        g2 = SE{2}(
            [                    R' -R'*δ;
            zeros(T, 1, size(R',2))     1]
        )

        # closure
        @test g1 * g2 isa SE{2}

        # identity element
        @test identity(g1) == SE{2}(I(3))
        @test identity(g1) * g1 == g1
        @test g1 * identity(g1) == g1

        # inverse element
        @test inv(g1) == g2
        @test g1 * g2 ≈ identity(g1)
        @test g2 * g1 ≈ identity(g1)

        # associativity
        g3 = SE{2}(
            [cos(ϕ) -sin(ϕ)  3;
             sin(ϕ)  cos(ϕ) -2;
                  0       0  1]
        )
        @test (g1 * g2) * g3 ≈ g1 * (g2 * g3)


    end
end
