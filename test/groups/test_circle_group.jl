using LieGroups, Random, Test

using ManifoldsBase: ℂ, ℝ
using Manifolds: Circle, Sphere

using StaticArrays

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "The Circle group" begin
    @testset "Complex Circle" begin
        C1 = CircleGroup()
        z1, z2, z3 = 1.0, 1.0im, -1.0im
        X1, X2, X3 = 1.0im, 0.0, -2.0im
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
                hat,
                inner,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                log,
                norm,
                rand,
                show,
                vee,
            ],
        )
        expectations = Dict(:repr => "CircleGroup()")
        test_lie_group(C1, properties, expectations)

        properties[:Name] = "IsBit-representation of Points and Vectors"
        properties[:Mutating] = false
        properties[:Points] = [z1, z2, z3]
        properties[:Vectors] = [X1, X2, X3]
        test_lie_group(C1, properties, expectations)
    end
    @testset "Real Circle" begin
        C2 = CircleGroup(Circle(ℝ))
        x1, x2, x3 = 0.0, 1.0, -π
        X1, X2, X3 = 1.0, 3.0, -123.0
        properties = Dict(
            :Name => "Array Points and Array Vectors",
            :Points => [fill(x1), fill(x2), fill(x3), [x1], [x2], [x3]],
            :Vectors => [fill(X1), fill(X2), fill(X3), [X1], [X2], [X3]],
            :Rng => Random.MersenneTwister(),
            :Mutating => true,
            :Functions => [
                compose,
                conjugate,
                diff_conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                exp,
                hat,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                log,
                rand,
                show,
                vee,
            ],
        )
        expectations = Dict(:repr => "CircleGroup(ℝ)")
        test_lie_group(C2, properties, expectations)

        properties[:Name] = "IsBit-representation of Points and Vectors"
        properties[:Mutating] = false
        properties[:Points] = [x1, x2, x3]
        properties[:Vectors] = [X1, X2, X3]
        test_lie_group(C2, properties, expectations)

        @testset "Edge cases" begin
            @test compose(C2, 1.0, fill(1.0)) == 2.0
            @test compose(C2, fill(1.0), 1.0) == 2.0
            @test LieGroups.sym_rem(fill(Float64(π))) == fill(Float64(-π))
            a = -π + 0.1
            b = -0.2
            @test compose(C2, a, b) == LieGroups.sym_rem(a + b)
            @test compose(C2, fill(a), fill(b)) == fill(LieGroups.sym_rem(a + b))
            @test compose(C2, fill(a), b) == LieGroups.sym_rem(a + b)
            @test compose(C2, a, fill(b)) == LieGroups.sym_rem(a + b)
        end
    end
    @testset "Planar Circle" begin
        C1 = CircleGroup(ℝ^2)
        z1, z2, z3 = [1.0, 0.0], [0.0, 1.0], [0.0, -1]
        X1, X2, X3 = [0.0, 1.0], [0.0, 0.0], [0.0, -2.0]
        properties = Dict(
            :Name => "Planar Circle",
            :Points => [z1, z2, z3],
            :Vectors => [X1, X2, X3],
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
                hat,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                log,
                rand,
                show,
                vee,
            ],
        )
        expectations = Dict(:repr => "CircleGroup(Sphere(1))")
        test_lie_group(C1, properties, expectations)
        @test identity_element(C1) == [1.0, 0.0]
    end

    @testset "StaticArrays.jl specializations" begin
        C = CircleGroup(Circle(ℝ))
        a = LieGroups.get_vector_lie(
            LieAlgebra(C), SA[2.0], DefaultLieAlgebraOrthogonalBasis(), SVector{1}
        )
        @test a === SA[2.0]
    end
end
