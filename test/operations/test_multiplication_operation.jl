using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

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
        @test ea * e === e
        @test e * ea === e
        # Zero array
        g2 = fill(2.0, ())
        # Array
        g3 = [2.0 0.0; 0.0 2.0]
        M = LieGroupsTestSuite.DummyManifold()
        G = LieGroup(M, MatrixMultiplicationGroupOperation())
        @test inv(G, e) === e
    end
end
