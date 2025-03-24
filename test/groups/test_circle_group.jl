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
        :Name => "Array Points and Array Vectors",
        :Points => [fill(z1), fill(z2), fill(z3)],  
        :Vectors => [fill(X1), fill(X2), fill(X3)],
        :Mutating => true,
        :Rng => Random.MersenneTwister(),
        :Functions => [
                compose,
                conjugate,
                diff_conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                exp,
                #hat,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                #log,
                #rand,
                show,
                #vee,
            ],
    )
    expectations = Dict(
        :repr => "CircleGroup()"
    )
    test_lie_group(C1, properties, expectations)


    properties[:Name] = "IsBit-respresentation of Points and Vectors"
    properties[:Mutating] = false
    properties[:Points] = [z1, z2, z3]
    properties[:Vectors] = [X1, X2, X3]
    test_lie_group(C1, properties, expectations)

end

@testset "Real Circle" begin
        C2 = RealCircleGroup()
        x1, x2, x3 = 0., 1., -π
        X1, X2, X3 = 1., 3., -123.
        properties = Dict(
            :Name => "Array Points and Array Vectors",
            :Points => [fill(x1), fill(x2), fill(x3)],
            :Vectors => [fill(X1), fill(X2), fill(X3)],
            :Rng => Random.MersenneTwister(),
            :Mutating => true,
            :Functions => [
                compose,
                # conjugate,
                # diff_conjugate,
                # diff_inv,
                # diff_left_compose,
                # diff_right_compose,
                # exp,
                # hat,
                #inv,
                # inv_left_compose,
                # inv_right_compose,
                is_identity,
                # lie_bracket,
                # log,
                #rand,
                # show,
                #vee,
            ],
        )
        expectations = Dict(
        :repr => "RealCircleGroup()"
    )
    #test_lie_group(C2, properties, expectations)

    properties[:Name] = "IsBit-representation of Points and Vectors"
    properties[:Mutating] = false
    properties[:Points] = [x1, x2, x3]
    properties[:Vectors] = [X1, X2, X3]
    #test_lie_group(C2, properties, expectations)
        
end
