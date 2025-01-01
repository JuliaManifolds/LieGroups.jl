using LieGroups, Random, Test

using ManifoldsBase: ℂ

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Heisenberg group" begin
    G = HeisenbergGroup(1)
    g1, g2, g3 = [1.0 2.0 3.0; 0.0 1.0 -1.0; 0.0 0.0 1.0],
    [1.0 4.0 -3.0; 0.0 1.0 3.0; 0.0 0.0 1.0],
    [1.0 -2.0 1.0; 0.0 1.0 1.1; 0.0 0.0 1.0]

    X1, X2, X3 = [0.0 2.0 3.0; 0.0 0.0 -1.0; 0.0 0.0 0.0],
    [0.0 4.0 -3.0; 0.0 0.0 3.0; 0.0 0.0 0.0],
    [0.0 -2.0 1.0; 0.0 0.0 1.1; 0.0 0.0 0.0]

    properties = Dict(
        :Name => "Heisenberg group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            adjoint,
            compose,
            conjugate,
            diff_conjugate,
            diff_inv,
            diff_left_compose,
            diff_right_compose,
            exp,
            # hat, #???
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            rand,
            show,
            #vee, #???
        ],
    )
    expectations = Dict(:repr => "HeisenbergGroup(1)", :lie_bracket => X1 * X2 - X2 * X1)
    test_lie_group(G, properties, expectations)

    @test is_point(G, Identity(G); error=:error)
    @test_throws DomainError is_point(G, Identity(AdditionGroupOperation()); error=:error)
end
