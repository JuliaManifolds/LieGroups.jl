using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Symplectic" begin
    G = SymplecticGroup(4)
    properties = Dict(
        :Name => "The symplectic group",
        #:Points => [g1, g2, g3],
        #:Vectors => [X1, X2, X3],
        #:Rng => Random.MersenneTwister(),
        :Functions => [
            #adjoint,
            #compose,
            #conjugate,
            #diff_inv,
            #diff_left_compose,
            #diff_right_compose,
            #exp,
            #hat,
            #identity_element,
            #inv,
            #inv_left_compose,
            #inv_right_compose,
            #is_identity,
            #lie_bracket,
            #log,
            # rand,
            show,
            #vee,
        ],
    )
    expectations = Dict(:repr => "SymplecticGroup(4, ℝ)")
    test_lie_group(G, properties, expectations)
end
