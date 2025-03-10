using LieGroups, Random, Test

using ManifoldsBase: ℂ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Complex Circle" begin
    C1 = CircleGroup()



    z1, z2, z3 = 1., 1.0im , -1.0im
    X1, X2, X3 = 1.0im, 0., -2.0im
    properties = Dict(
        :Name => "Array Points",
        :Points => [fill(z1), fill(z2), fill(z3)],  
        :Vectors => [fill(X1), fill(X2), fill(X3)],
        :Mutating => true,
        :Rng => Random.MersenneTwister(),
        :Functions => [
                compose,
                conjugate,
                # diff_conjugate,
                # diff_inv,
                # diff_left_compose,
                # diff_right_compose,
                exp,
                # hat,
                #
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                # lie_bracket,
                # log,
                # rand,
                show,
                # vee,
            ],
    )
    expectations = Dict(
        :repr => "CircleGroup()"
    )
    test_lie_group(C1, properties, expectations)

    #add testing for references
    #properties[:Name] = "Ref Points"
    #properties[:Mutating] = false
    #properties[:Points] = [Ref(z1), Ref(z2), Ref(z3)]

    test_lie_group(C1, properties, expectations)

    @testset "complex" begin
        C1 = CircleGroup()
        @test manifold_dimension(C1) == 1
        e=Identity(C1)
        s = 1.0im
        t = 2.0+1im     
        @test compose(C1, e, s) == s   
        @test is_point(C1, s; error=:error)
        @test_throws DomainError is_point(C1, t; error=:error)
    end
end

@testset "Real Circle" begin
        C2 = RealCircleGroup()
        x1, x2, x3 = 0., 1., -π
        properties = Dict(
            :Name => "The real circle group",
            :Points => [x1, x2, x3],
            :Rng => Random.MersenneTwister(),
            :Mutating => false,
            :Functions => [
                #compose,
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
    #test_lie_group(C2, properties, expectations)
        
end
