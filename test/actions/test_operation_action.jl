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

    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    A1 = GroupOperationAction(a1, G)
    A2 = GroupOperationAction(a2, G)
    A3 = GroupOperationAction(a3, G)
    A4 = GroupOperationAction(a4, G)
    @test inv(A1) == A3
    @test inv(A2) == A4
    @test inv(A3) == A1
    @test inv(A4) == A2
    As = [A1, A2, A3, A4]
    for A in As
        @test base_Lie_group(A) == G
    end
end
