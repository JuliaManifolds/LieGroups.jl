using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Special Unitary Group" begin
    G = SpecialUnitaryGroup(2)
    g1 = 1 / sqrt(2) * ComplexF64[1.0 1.0; -1.0 1.0]
    g2 = ComplexF64[0.0 -1.0; 1.0 0.0]
    g3 = ComplexF64[1.0 0.0; 0.0 1.0]
    X1 = [0.0 1.0im; -1.0im 0.0]
    X2 = ComplexF64[0.0 1.0; -1.0 0.0]
    X3 = ComplexF64[0.0 -0.5; 0.5 0.0]
    properties = Dict(
        :Name => "The special unitary group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => [
            # adjoint,
            # compose,
            # conjugate,
            # diff_inv,
            # diff_left_compose,
            # diff_right_compose,
            exp,
            hat,
            # identity_element,
            # inv,
            # inv_left_compose,
            # inv_right_compose,
            # is_identity,
            # lie_bracket,
            log,
            # rand,
            show,
            vee,
        ],
    )
    expectations = Dict(:repr => "SpecialUnitaryGroup(2)")
    test_lie_group(G, properties, expectations)
end
