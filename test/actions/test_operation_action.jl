using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Group Operation as a Group Action" begin
    a1 = RightGroupOperation()
    a2 = LeftGroupOperation()

    @test switch(a1) == a2
    @test switch(a2) == a1

    a3 = inv(a1)
    @test inv(a3) == a1
    a4 = inv(a2)
    @test inv(a4) == a2

    @test switch(a3) == a4
    @test switch(a4) == a3
end
