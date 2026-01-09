using LieGroups, Test

@testset "Group Action Interface" begin
    M = LieGroups.Test.DummyManifold()
    T = LieGroups.Test.DummyLeftActionType()
    G = LieGroups.Test.DummyLieGroup(M, LieGroups.Test.DummyOperation())
    a = GroupAction(G, M, T)
    @test repr(a) === "GroupAction($G, $M, $T)"
    @test switch(T) === LieGroups.Test.DummyRightActionType()
    @test switch(a) === GroupAction(G, M, switch(T))
end
