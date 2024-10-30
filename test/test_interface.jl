using LieGroups, Test, ManifoldsBase

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic Lie Group Interface functions" begin
    M = LieGroupsTestSuite.DummyManifold()
    op = LieGroupsTestSuite.DummyOperation()
    G = LieGroup(M, op)
    rs = "LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation())"
    @test repr(G) == rs
    𝔤 = LieAlgebra(G)
    op2 = LieGroupsTestSuite.DummySecondOperation()
    rs2 = "LieAlgebra( LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()) )"
    @test repr(𝔤) == rs2
    @test is_identity(G, Identity(op))
    @test !is_identity(G, Identity(op2))
    @test base_manifold(G) === M
    e = Identity(op)
    @test compose(G, e, e) == e
    @test compose!(G, e, e, e) === e
    @test isapprox(G, e, Identity(op))
    @test !isapprox(G, e, Identity(op2))
    @test is_point(G, e)
    @test !is_point(G, Identity(op2))
    @test_throws DomainError is_point(G, Identity(op2); error=:error)
    @testset "Methoderrors for the non-implemented mutating variants" begin
        g = :none
        h = :alsonone
        X = :nonetoo
        begin # locally define identity element
            LieGroups.identity_element(::typeof(G)) = :id
            LieGroups.exp!(::typeof(G), h, ::typeof(e), X, t::Number=1) = :id
            @test exp(G, e, X) === :id
            # delete both methods again
            Base.delete_method(which(identity_element, (typeof(G),)))
            Base.delete_method(which(exp!, typeof.([G, h, e, X, 1])))
            #
            # same for log
            ManifoldsBase.allocate_result(::typeof(G), ::typeof(log), g) = :g
            LieGroups.log!(::typeof(G), X, ::Identity, g) = :g
            @test log(G, e, g) === :g
            # delete both methods again
            Base.delete_method(which(ManifoldsBase.allocate_result, typeof.([G, log, g])))
            Base.delete_method(which(log!, typeof.([G, X, e, g])))
        end
        # so they arae undefined here again but we checked the exp fallback
        @test_throws MethodError exp!(G, g, e, X)
        @test_throws MethodError log!(G, X, e, g)
    end
end
