using LieGroups, LinearAlgebra, Test

@testset "Multiplication Operation" begin
    @testset "Base.:* and Base.:* with the Identity" begin
        # Generic & Number
        e = Identity(MatrixMultiplicationGroupOperation())
        @test (e * e) === e
        g = 2
        @test (g * e) == g
        @test (e * g) == g
        @test (g / e) == g
        @test (e \ g) == g
        @test (e / g) == 1 / g
        @test (g \ e) == 1 / g
        @test (e / e) === e
        @test (e \ e) === e
        @test inv(e) === e
        @test det(e)
        ea = Identity(AdditionGroupOperation)
        @test ea * e === ea #Id(Add)*Id(Mul) = Id(Add) or 0*1 = 0
        @test e * ea === ea
        M = LieGroups.Test.DummyManifold()
        G = LieGroup(M, MatrixMultiplicationGroupOperation())
        # Zero-dimensional array
        g2 = fill(2.0, ())
        X2 = fill(1.0, ())
        @test diff_inv(G, g2, X2) == fill(-1.0, ())
        # Array
        g3 = [2.0 0.0; 1.0 2.0]
        @test inv(G, e) === e
        h3 = zero(g3)
        X3 = [1.0 0.0; 0.0 2.0]
        @test diff_conjugate(G, g3, g3, X3) == g3 * X3 * inv(g3)
        # inplace-multiplication with e
        h3 = zero(g3)
        mul!(h3, e, g3)
        @test h3 == g3
        h3 = zero(g3)
        mul!(h3, g3, e)
        @test h3 == g3
        mul!(h3, e, e)
        @test h3 == I
        @test mul!(e, e, e) === e
        @test one(e) === e
    end
    @testset "Abelian Multiplication edge cases" begin
        G = LieGroup(
            LieGroups.Test.DummyManifold(), AbelianMultiplicationGroupOperation()
        )
        # Edge case: Number element
        @test identity_element(G, Float64) == 1.0
        # Edge case: Mixed compose
        @test compose(G, fill(1.0), 1.0) == 1.0
        @test compose(G, 1.0, fill(1.0)) == 1.0
    end
end
