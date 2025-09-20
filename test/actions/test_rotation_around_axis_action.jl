using LieGroups, Test
using Manifolds

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Rotation around axis action" begin
    M = Euclidean(3)
    G = CircleGroup(Circle(ℝ))
    axis = [sqrt(2) / 2, sqrt(2) / 2, 0.0]
    A = GroupAction(LieGroups.RotationAroundAxisAction(axis), G, M)

    types_a = [Ref(Float64)]

    types_m = [Vector{Float64}]

    @test isa(A.type, AbstractLeftGroupActionType)

    for (i, T_A, T_M) in zip(1:length(types_a), types_a, types_m)
        angles = (0.0, π / 2, 2π / 3, π / 4)
        a_pts = convert.(T_A, [angles...])
        m_pts = convert.(T_M, [[0.0, 1.0, 0.0], [-1.0, 0.0, 1.0], [1.0, 1.0, -2.0]])
        X_pts = convert.(T_M, [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

        properties = Dict(
            :AlgebraVectors => X_pts,
            :Functions => [apply, base_lie_group, base_manifold],
            :GroupPoints => a_pts,
            :ManifoldPoints => m_pts,
            :TangentVectors => m_pts,
            :Name => "",
        )
        expectations = Dict(:manifold => M, :group => G, :atol => 1.0e-15)

        test_group_action(A, properties, expectations)
    end
    # make sure angle can be in a one-element vector
    @test apply(A, [π / 2], [0.0, 1.0, 0.0]) ≈ [0.5, 0.5, sqrt(2) / 2]
end
