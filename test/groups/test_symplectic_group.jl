using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "The symplectic group" begin
    G = SymplecticGroup(2)
    # for this case we need [a c; c b] such that ab-c^2=1
    g1 = [1.0 2.0; 2.0 5.0]
    g2 = [2.0 3.0; 3.0 5.0]
    g3 = [4.0 5.0; 5.0 6.5]
    # to be tangent at id we need [a b; c; -a]
    X1 = [0.0 1.0; 1.0 0.0]
    X2 = [2.0 1.0; 1.0 -2.0]
    X3 = [1.0 2.0; 3.0 -1.0]
    properties = Dict(
        :Name => "The symplectic group",
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
            #hat,
            identity_element,
            inv,
            inv_left_compose,
            inv_right_compose,
            is_identity,
            lie_bracket,
            log,
            rand,
            show,
            #vee,
        ],
    )
    expectations = Dict(
        :atols => Dict(rand => 1e-13, log => 1e-13), :repr => "SymplecticGroup(2, ℝ)"
    )
    test_lie_group(G, properties, expectations)
end
