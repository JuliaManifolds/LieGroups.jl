using LieGroups, Test

@testset "Addition Operation" begin
    @testset "Base.:+ and Base.:- with the Identity" begin
        e = Identity(AdditionGroupOperation())
        @test (+e) === e
        @test (e + e) === e
        g = 2
        @test (g + e) == g
        @test (e + g) == g
        @test (g - e) == g
        @test (e - g) == -g
        @test (e - e) === e
    end
end
