using LieGroups, Test

@testset "Addition Operation" begin
    @testset "Base.:+ and Base.:- with the Identity" begin
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
    end
end
