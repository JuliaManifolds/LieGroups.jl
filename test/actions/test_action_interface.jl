using LieGroups, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Group Action Interface" begin
    M = LieGroupsTestSuite.DummyManifold()
    T = LieGroupsTestSuite.DummyLeftActionType()
    G = LieGroupsTestSuite.DummyLieGroup(M, LieGroupsTestSuite.DummyOperation())
    a = GroupAction(G, M, T)
    @test repr(a) === "GroupAction($G, $M, $T)"
    @test switch(T) === LieGroupsTestSuite.DummyRightActionType()
    @test switch(a) === GroupAction(G, M, switch(T))
end
