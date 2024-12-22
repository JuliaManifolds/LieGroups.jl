using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite
using LieGroupsTestSuite: rotation_matrix
begin
    G = OrthogonalGroup(2)
    g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
    g2 = [0.0 -1.0; 1.0 0.0]
    g3 = [1.0 0.0; 0.0 1.0]
    X1, X2, X3 = [0.0 0.1; -0.1 0.0], [0.0 -0.2; 2.0 0.0], [0.0 0.0; 0.0 0.0]
    properties = Dict(
        :Name => "The orthogonal group O(2)",
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
    expectations = Dict(:repr => "OrthogonalGroup(2)")
    test_lie_group(G, properties, expectations)
    # O(3)
    H = OrthogonalGroup(3)
    h1 = rotation_matrix(3, 2, 1, π / 4)
    h2 = rotation_matrix(3, 2, 1, π / 8) * rotation_matrix(3, 3, 1, π / 4)
    h3 = rotation_matrix(3, 3, 1, π / 4) * rotation_matrix(3, 3, 2, π / 8)
    Y1 = [0.0 0.1 0.0; -0.1 0.0 0.0; 0.0 0.0 0.0]
    Y2 = [0.0 0.0 -0.2; 0.0 0.0 0.0; 0.2 0.0 0.0]
    Y3 = [0.0 0.3 0.0; -0.3 0.0 0.4; 0.0 -0.4 0.0]
    # Test only specialized functions here
    properties2 = Dict(
        :Name => "The orthogonal group O(3) – specialised funcions",
        :Points => [h1, h2, h3],
        :Vectors => [Y1, Y2, Y3],
        :Functions => [exp, log, show],
    )
    expectations2 = Dict(:repr => "OrthogonalGroup(3)", :atols => Dict(:exp => 1e-15))
    test_lie_group(H, properties2, expectations2)
end
