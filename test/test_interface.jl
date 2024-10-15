using LieGroups, Test

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    rs = "LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation())"
    @test repr(G) == rs
    𝔤 = LieAlgebra(G)
    op2 = LieGroupsTestSuite.DummySecondOperation()
    rs2 = "Lie Algebra( LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()) )"
    @test repr(𝔤) == rs2
    @test is_identity(G, Identity(op))
    @test !is_identity(G, Identity(op2))
    @test base_manifold(G) === M
    e = Identity(op)
    @test compose(G, e, e) == e
    @test compose!(G, e, e, e) === e
    @test is_point(G, e)
    @test !is_point(G, Identity(op2))
    @test_throws DomainError is_point(G, Identity(op2); error=:error)
    # Exp log Method Error fallbacks that avoid the stack overflow
    g = :none
    X = :nonetoo
    @test_throws MethodError exp!(G, g, e, X)
    @test_throws MethodError log!(G, X, e, g)
end
