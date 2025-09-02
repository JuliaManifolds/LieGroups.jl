using LieGroups, Test, ManifoldsBase, RecursiveArrayTools

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Generic semidirect product Lie group" begin
    @testset "Show" begin
        M = LieGroupsTestSuite.DummyManifold()
        op1 = LieGroupsTestSuite.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroupsTestSuite.DummySecondOperation()
        G2 = LieGroup(M, op2)
        # The rest requires tests with a concrete semi-direct one
        fcts = [show]

        Gl = LeftSemidirectProductLieGroup(G1, G2, LeftGroupOperationAction())
        properties = Dict(
            :Name => "The Left Semidirect Product Manifold",
            :Functions => fcts,
        )
        expectations = Dict(
            :repr => "LeftSemidirectProductLieGroup($(G1), $(G2), LeftGroupOperationAction())",
        )
        test_lie_group(Gl, properties, expectations)

        Gr = RightSemidirectProductLieGroup(G1, G2, RightGroupOperationAction())
        properties = Dict(
            :Name => "The Right Semidirect Product Manifold",
            # :Rng => Random.MersenneTwister(),
            :Functions => fcts,
        )
        expectations = Dict(
            :repr => "RightSemidirectProductLieGroup($(G1), $(G2), RightGroupOperationAction())",
        )
        test_lie_group(Gr, properties, expectations)
    end
    @testset "Interaction os semidirect and product" begin
        # Make sure the inner product manifold stays “nested” and does not collapse to a 4-product manifold
        M = LieGroupsTestSuite.DummyManifold()
        op1 = LieGroupsTestSuite.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroupsTestSuite.DummySecondOperation()
        G2 = LieGroup(M, op2)


        Gl = LeftSemidirectProductLieGroup(G1 × G2, G2 × G1, LeftGroupOperationAction())
        @test Gl.manifold isa ProductManifold
        @test length(Gl.manifold.manifolds) == 2 # we still just have 2 manifolds, those of the semiproduct
        # And it is a product of products
        @test Gl.manifold.manifolds[1] isa ProductManifold
        @test Gl.manifold.manifolds[2] isa ProductManifold
        @test Gl.manifold.manifolds[1].manifolds[1] === G1.manifold
        @test Gl.manifold.manifolds[2].manifolds[1] === G2.manifold
    end
    @testset "8 combinations with (SO(2), ℝ²)" begin
        struct TestLeftAction <: AbstractLeftGroupActionType end
        function LieGroups.apply!(A::GroupAction{TestLeftAction}, k, g, h)
            return k .= g * h
        end
        function LieGroups.diff_apply!(A::GroupAction{TestLeftAction}, Y, g, h, X)
            return Y .= g * X
        end
        function LieGroups.diff_group_apply!(A::GroupAction{TestLeftAction}, Y, g, h, X)
            return Y .= X * h
        end

        struct TestRightAction <: AbstractRightGroupActionType end
        function LieGroups.apply!(A::GroupAction{TestRightAction}, k, g, h)
            return k .= inv(A.group, g) * h # (h' * g)'
        end
        function LieGroups.diff_apply!(A::GroupAction{TestRightAction}, Y, g, h, X)
            return Y .= inv(A.group, g) * X
        end
        function LieGroups.diff_group_apply!(A::GroupAction{TestRightAction}, Y, g, h, X)

            return Y .= diff_inv(A.group, g, X) * h
        end
        g1 = 1 / sqrt(2) * [1.0 1.0; -1.0 1.0]
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        h1 = [0.1, 0.2]
        h2 = [0.0, 1.0]
        h3 = [0.0, 0]

        X1 = [0.0 0.1; -0.1 0.0]
        Y1 = [0.0, 0.2]

        for action in (TestLeftAction(), TestRightAction()),
                action_on in (ActionActsOnLeft(), ActionActsOnRight())

            G = LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2), action; action_on)
            p1 = ArrayPartition(copy(g1), copy(h1))
            p2 = ArrayPartition(copy(g2), copy(h2))
            p3 = ArrayPartition(copy(g3), copy(h3))
            V1 = ArrayPartition(copy(X1), copy(Y1))

            properties = Dict(
                :Name => "LeftSemidirectProductLieGroup, $(supertype(typeof(action))), $(action_on)",
                :Points => [p1, p2, p3],
                :Vectors => [V1],
                :Functions => [identity_element, is_identity, inv, compose, diff_left_compose, diff_right_compose, show],
            )

            H = SpecialOrthogonalGroup(2)
            N = TranslationGroup(2)
            if (action == TestLeftAction()) != (action_on == ActionActsOnLeft())
                # Same side -> forward action -> (g,h)⁻¹ =  (g⁻¹, α(g⁻¹)(h⁻¹)
                α_g_h = apply(GroupAction(action, H, N), inv(H, g1), inv(N, h1))
            else
                # Opposite sides -> inverse action -> (g,h)⁻¹ =  (g⁻¹, α(g⁻¹)⁻¹(h⁻¹)
                α_g_h = apply(GroupAction(action, H, N), g1, inv(N, h1))
            end
            if G.op isa LeftSemidirectProductGroupOperation
                inv_p1 = ArrayPartition(inv(g1), α_g_h)
            else
                inv_p1 = ArrayPartition(α_g_h, inv(g1))
            end

            expectations = Dict(
                :atol => 1.0e-14,
                :inv => ArrayPartition(inv(g1), α_g_h),
            )
            test_lie_group(G, properties, expectations)

            G = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), action; action_on)
            q1 = ArrayPartition(copy(h1), copy(g1))
            q2 = ArrayPartition(copy(h2), copy(g2))
            q3 = ArrayPartition(copy(h3), copy(g3))
            W1 = ArrayPartition(copy(Y1), copy(X1))
            properties = Dict(
                :Name => "RightSemidirectProductLieGroup, $(supertype(typeof(action))), $(action_on)",
                :Points => [q1, q2, q3],
                :Vectors => [W1],
                :Functions => [identity_element, is_identity, inv, compose, diff_left_compose, diff_right_compose, show],
            )

            expectations = Dict(
                :atol => 1.0e-14,
                :inv => ArrayPartition(α_g_h, inv(g1)),
            )

            test_lie_group(G, properties, expectations)
        end
    end
    # sanity check with matrix representation
    # reference data for (SO(2), ℝ²) semidirect products
    R = exp(SpecialOrthogonalGroup(2), [0.0 0.1; -0.1 0.0])
    t = [0.1, 0.2]
    l1_mat = vcat([1 0 0], hcat(t, R)) # build the matrix representation for a left action semidirect product - acts on the right component
    r1_mat = hcat(vcat(R, t'), [0, 0, 1]) # build the matrix representation for a right action semidirect product - acts on the left component
    # inverse of R should be the same for both cases
    inv_R = inv(SpecialOrthogonalGroup(2), R)
    # for t components, use the matrix inverse as reference
    l_inv_t = inv(l1_mat)[2:3, 1] # extract the translation part
    r_inv_t = inv(r1_mat)[3, 1:2] # extract the translation part
    #now compare against [Left|Right]SemidirectProductLieGroup
    Gl = LeftSemidirectProductLieGroup(SpecialOrthogonalGroup(2), TranslationGroup(2), TestLeftAction(); action_on = ActionActsOnRight())
    @test isapprox(inv(Gl, ArrayPartition(R, t)), ArrayPartition(inv_R, l_inv_t), atol = 1.0e-14)
    Gr = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), TestRightAction(); action_on = ActionActsOnLeft())
    @test isapprox(inv(Gr, ArrayPartition(t, R)), ArrayPartition(r_inv_t, inv_R), atol = 1.0e-14)
end
