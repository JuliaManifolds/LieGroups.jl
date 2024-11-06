using LieGroups, Test, ManifoldsBase

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic power Lie group" begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    pG = G^3
    rs = "PowerLieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation(), 3)"

    # test_lie_group(G, properties, expectations)
end
