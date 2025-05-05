using LieGroups, Test

s = joinpath(@__DIR__, "../LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

using StaticArrays

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
        @test identity_element(G, Float64) == 0.0
        @test identity_element(G, Array{Float64,0}) == fill(0.0)
        @test identity_element(G, SArray{Tuple{},Float64}) == @SArray fill(0.0)
    end
end
