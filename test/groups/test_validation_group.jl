using LieGroups, Random, Test

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Validation Lie group" begin
    G = TranslationGroup(3)
    g1, g2, g3 = [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.1, 1.2, 3.3]
    X1, X2, X3 = [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.1, 0.2, 0.3]
    VG = ValidationLieGroup(G)
    vg1, vg2, vg3 = ValidationMPoint.([g1, g2, g3])
    vX1, vX2, vX3 = ValidationLieAlgebraTangentVector.([X1, X2, X3])
    fcts = [
        adjoint,
        compose,
        conjugate,
        diff_inv,
        diff_left_compose,
        diff_right_compose,
        exp,
        hat,
        identity_element,
        inv,
        inv_left_compose,
        inv_right_compose,
        is_identity,
        lie_bracket,
        log,
        rand,
        show,
        vee,
    ]
    properties = Dict(
        :Name => "Validation of Translation group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :repr => "ValidationLieGroup of TranslationGroup(3; field=ℝ)\n    * mode = :error\n",
        :diff_inv => -X1,
        :diff_left_compose => X1,
        :diff_right_compose => X1,
        :lie_bracket => zero(X1),
    )
    test_lie_group(VG, properties, expectations)

    properties = Dict(
        :Name => "Validation of Translation group with types points/vectors",
        :Points => [vg1, vg2, vg3],
        :Vectors => [vX1, vX2, vX3],
        :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :Name => "ValidationLieGroup of TranslationGroup(3; field=ℝ)\n    * mode = :error\n",
        :diff_left_compose => X1,
        :diff_right_compose => X1,
        :lie_bracket => zero(X1),
    )
    test_lie_group(VG, properties, expectations)
end
