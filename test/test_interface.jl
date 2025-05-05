using LieGroups, Test, ManifoldsBase, Random

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
    rs2 = "LieAlgebra(LieGroup(LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()))"
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
            LieGroups.identity_element(::typeof(G), T::Type) = :id
            LieGroups.identity_element!(::typeof(G), g) = (g[] = :id)
            ManifoldsBase.allocate_result(::typeof(G), ::typeof(LieGroups.exp), g) = :a
            LieGroups.exp!(::typeof(G), h, X) = :id
            @test exp(G, X) === :a
            #
            # same for log
            ManifoldsBase.allocate_result(::typeof(G), ::typeof(LieGroups.log), g) = :g
            LieGroups.log!(::typeof(G), X, g) = :g
            @test log(G, g) === :g
            g2 = Ref(:g)
            inv!(G, g2, e)
            @test g2[] == :id
            # delete methods again
            Base.delete_method(which(identity_element, (typeof(G),)))
            Base.delete_method(which(identity_element, (typeof(G), Type)))
            Base.delete_method(which(identity_element!, typeof.([G, g2])))
            Base.delete_method(which(LieGroups.exp!, typeof.([G, h, X])))
            Base.delete_method(which(ManifoldsBase.allocate_result, typeof.([G, log, g])))
            Base.delete_method(which(ManifoldsBase.allocate_result, typeof.([G, exp, g])))
            Base.delete_method(which(LieGroups.log!, typeof.([G, X, g])))
        end
        # both undefined again.
        @test_throws MethodError exp!(G, g, X)
        @test_throws MethodError log!(G, X, g)
    end
end
@testset "Generic Lie Algebra Interface functions" begin
    @testset "Generic get_coordinates/get_vector passthrough on 𝔤" begin
        M = ManifoldsBase.DefaultManifold(2)
        op = AdditionGroupOperation()
        G = LieGroup(M, op)
        𝔤 = LieAlgebra(G)
        @test startswith(sprint(show, "text/plain", 𝔤), "The Lie algebra of")
        B = DefaultOrthonormalBasis()
        p = [1.0, 2.0]
        q = [0.0, 0.0]
        # coordinates and vector on 𝔤 are here the same as the ones on M at 0
        # Similarly: on G they are the same even for p
        X = [1.0, 0.0]
        @test default_basis(G) == DefaultLieAlgebraOrthogonalBasis()
        @test get_coordinates(𝔤, X, B) == get_coordinates(M, q, X, B)
        @test get_coordinates(G, p, X, B) == get_coordinates(M, q, X, B)
        Y = copy(X)
        @test get_coordinates!(𝔤, Y, X, B) == get_coordinates!(M, Y, q, X, B)
        @test X == Y
        @test get_coordinates!(G, Y, p, X, B) == get_coordinates!(M, Y, q, X, B)
        @test X == Y
        c = [0.0, 1.0]
        @test get_vector(𝔤, c, B) == get_vector(M, q, c, B)
        @test get_vector(G, p, c, B) == get_vector(M, q, c, B)
        @test get_vector(𝔤, c, B; tangent_vector_type=Vector{Float64}) ==
            get_vector(M, q, c, B)
        d = copy(c)
        @test get_vector!(𝔤, d, c, B) == get_vector!(M, d, q, c, B)
        @test c == d
        @test get_vector!(G, d, p, c, B) == get_vector!(M, d, q, c, B)
        @test c == d
        @test project(G, p) == project(M, p)
        @test project(𝔤, X) == project(M, p, X)
        @test project(𝔤, X, X) == project(M, p, X)

        # B2
        B2 = DefaultLieAlgebraOrthogonalBasis()
        @test get_vector(𝔤, c, B2) == c #on Euclidean this is the same
        @test get_vector(𝔤, c, B2; tangent_vector_type=Vector{Float64}) == c
        @test get_vector(𝔤, c, B2; tangent_vector_type=Vector{Float64}) == c
        @test hat(𝔤, c) == c # and hat/vee as well
        # Real fallback test here not 100% accurate
        𝔤2 = LieAlgebra(LieGroup(ManifoldsBase.DefaultManifold(), AdditionGroupOperation()))
        @test norm(𝔤, 2) == 2.0
        # Rand cases pass through tests
        Random.seed!(42)
        @test is_point(𝔤, rand(𝔤, Vector{Float64}))
        @test all(is_point.(Ref(𝔤), rand(𝔤, 3)))
        rng = Random.MersenneTwister()
        @test is_point(𝔤, rand(rng, 𝔤, Vector{Float64}))
    end
    @testset "Defaults on a nearly empty (non-implemented) Lie group" begin
        G = LieGroup(
            LieGroupsTestSuite.DummyManifold(), LieGroupsTestSuite.DummyOperation()
        )
        # Locally define zero vector
        LieGroups.zero_vector(::typeof(LieAlgebra(G))) = :id
        e = Identity(G)
        @test log(G, e) == :id
        Base.delete_method(which(zero_vector, (typeof(LieAlgebra(G)),)))
        # verify deletion
        @test_throws ErrorException zero_vector(LieAlgebra(G))
    end
    @testset "Test (new) rand(G, T, d)" begin
        M = ManifoldsBase.DefaultManifold(2)
        op = AdditionGroupOperation()
        G = LieGroup(M, op)
        𝔤 = LieAlgebra(G)
        q = zeros(2)
        Random.seed!(42)
        p = rand(G, Vector{Float64}; vector_at=q)
        @test is_point(G, p)
        pts1 = rand(G, 2)
        @test all(is_point.(Ref(G), pts1))
        pts2 = rand(G, Vector{Float64}, 2)
        @test all(is_point.(Ref(G), pts2))
        rng = Random.MersenneTwister()
        pts3 = rand(rng, G, 2)
        @test all(is_point.(Ref(G), pts3))
        pts4 = rand(rng, G, Vector{Float64}, 2)
        @test all(is_point.(Ref(G), pts4))
        p2 = rand(rng, G, Vector{Float64}; vector_at=q)
        @test is_point(G, p2)
    end
end
