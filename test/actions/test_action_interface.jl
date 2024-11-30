using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Group Action Interface" begin
    M = LieGroupsTestSuite.DummyManifold()
    T = LieGroupsTestSuite.DummyActionType()
    G = LieGroupsTestSuite.DummyLieGroup(M, LieGroupsTestSuite.DummyOperation())
    a = GroupAction(T, G, M)
    @test repr(a) === "GroupAction($T, $G, $M)"
    @test switch(T) === T
    @test switch(a) === a
end
