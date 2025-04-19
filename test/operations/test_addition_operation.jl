using LieGroups, Test

s = joinpath(@__DIR__, "../LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

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
        @test (-e) === e
        G = LieGroup(LieGroupsTestSuite.DummyManifold(), AdditionGroupOperation())
        @test identity_element(G, 1.0) == 0.0
    end
end
