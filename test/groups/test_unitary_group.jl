using LieGroups, ManifoldsBase, Quaternions, Random, Test
using ManifoldsBase: ℝ, ℂ, ℍ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Unitary Group" begin
    @testset "Complex Unitary group" begin
        G = UnitaryGroup(2)
        g1 = 1 / sqrt(2) * ComplexF64[1.0 1.0; -1.0 1.0]
        g2 = ComplexF64[0.0 -1.0; 1.0 0.0]
        g3 = ComplexF64[-1.0 0.0; 0.0 1.0] #one with determinant -1
        X1 = [0.0 1.0im; 1.0im 0.0]
        X2 = ComplexF64[0.0 1.0; -1.0 0.0]
        X3 = ComplexF64[0.0 -0.5; 0.5 0.0]
        properties = Dict(
            :Name => "The unitary group",
            :Points => [g1, g2, g3],
            :Vectors => [X1, X2, X3],
            :Rng => Random.MersenneTwister(),
            :Functions => [
                adjoint,
                compose,
                conjugate,
                diff_inv,
                diff_left_compose,
                diff_right_compose,
                exp,
                # hat,
                identity_element,
                inv,
                inv_left_compose,
                inv_right_compose,
                is_identity,
                lie_bracket,
                log,
                rand,
                show,
                # vee,
            ],
        )
        expectations = Dict(:repr => "UnitaryGroup(2)", :atol => 1e-14)
        test_lie_group(G, properties, expectations)
        G2 = UnitaryGroup(2; parameter=:field)
        expectations2 = Dict(:repr => "UnitaryGroup(2; parameter=:field)", :atol => 1e-14)
        test_lie_group(G2, properties, expectations2)
    end
    @testset "Quaternion Unitary group (Numbers)" begin
        G = UnitaryGroup(1, ℍ)
        g1 = quat(1.0)
        g2 = quat(0.0, 1.0, 0.0, 0.0)
        g3 = quat(0.0, 0.0, 1.0, 0.0)
        properties = Dict(
            :Name => "The quaternion unitary group",
            :Points => [g1, g2, g3],
            :Vectors => [quat(0.0)],
            :Mutating => false,
            :Rng => Random.MersenneTwister(),
            :Functions => [compose, inv, exp, log, rand, show],
        )
        expectations = Dict(:repr => "UnitaryGroup(1, ℍ)", :atols => Dict(exp => 1e-15))
        test_lie_group(G, properties, expectations)
        G2 = UnitaryGroup(1, ℍ; parameter=:field)

        expectations2 = Dict(:repr => "UnitaryGroup(1, ℍ; parameter=:field)")
        test_lie_group(
            G2,
            Dict(:Name => "The quaternion unitary group", :Functions => [show]),
            expectations2,
        )
        @test identity_element(G, QuaternionF64) == Quaternions.quat(1.0)
    end
    @testset "Quaternion Unitary group (1x1 matrices)" begin
        G = UnitaryGroup(1, ℍ)
        g1 = fill(quat(1.0), 1, 1)
        g2 = fill(quat(0.0, 1.0, 0.0, 0.0), 1, 1)
        g3 = fill(quat(0.0, 0.0, 1.0, 0.0), 1, 1)
        X1 = fill(quat(0.0), 1, 1)
        properties = Dict(
            :Name => "The quaternion unitary group",
            :Points => [g1, g2, g3],
            :Vectors => [X1],
            :Rng => Random.MersenneTwister(),
            :Functions => [compose, inv, exp, log, rand, show],
        )
        expectations = Dict(:repr => "UnitaryGroup(1, ℍ)")
        test_lie_group(G, properties, expectations)
        G2 = UnitaryGroup(1, ℍ; parameter=:field)

        expectations2 = Dict(:repr => "UnitaryGroup(1, ℍ; parameter=:field)")
        test_lie_group(
            G2,
            Dict(:Name => "The quaternion unitary group", :Functions => [show]),
            expectations2,
        )

        @test identity_element(G) isa Matrix{QuaternionF64}
    end
end
