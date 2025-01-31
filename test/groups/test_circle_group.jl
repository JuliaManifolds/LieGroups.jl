using LieGroups, Random, Test

using ManifoldsBase: ℂ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "(complex) Circle Group" begin
    C1 = CircleGroup()
    z1, z2, z3 = 1, 1im , -1im
    properties = Dict(
        :Name => "The Circle group",
        :Points => [z1, z2, z3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
                # compose,
                # conjugate,
                # diff_conjugate,
                # diff_inv,
                # diff_left_compose,
                # diff_right_compose,
                # exp,
                # hat,
                # inv,
                # inv_left_compose,
                # inv_right_compose,
                # is_identity,
                # lie_bracket,
                # log,
                # rand,
                # show,
                # vee,
            ],
    )
    expectations = Dict(
        :repr => "CircleGroup()"
    )
    test_lie_group(C1, properties, expectations)

    @testset "complex" begin
        C1 = CircleGroup()
        @test manifold_dimension(C1) == 1
        e=Identity(C1)
        s = 1im
        t = 2+1im     
        @test compose(C1, e, s) == s   
        @test is_point(C1, s; error=:error)
        @test_throws DomainError is_point(C1, t; error=:error)
    end
end

@testset "real Circle Group" begin
        C2 = RealCircleGroup()
        x1, x2, x3 = 0, 1, -π
        properties = Dict(
            :Name => "The real Circle group",
            :Points => [x1, x2, x3],
            :Rng => Random.MersenneTwister(),
            :Functions => [
                # compose,
                # conjugate,
                # diff_conjugate,
                # diff_inv,
                # diff_left_compose,
                # diff_right_compose,
                # exp,
                # hat,
                # inv,
                # inv_left_compose,
                # inv_right_compose,
                # is_identity,
                # lie_bracket,
                # log,
                # rand,
                # show,
                #vee,
            ],
        )
        expectations = Dict(
        :repr => "RealCircleGroup()"
    )
    test_lie_group(C2, properties, expectations)

        @testset "real" begin
            C2 = RealCircleGroup()
            @test manifold_dimension(C2) == 1
            e=Identity(C2) 
            p = 0.5 * π
            q = 2*π
            @test compose(C2, e, p) == p                     
            @test is_point(C2, p; error=:error)
            @test_throws DomainError is_point(C2, q; error=:error)
        end
end