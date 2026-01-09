using LieGroups, Test, ManifoldsBase, RecursiveArrayTools

@testset "Generic semidirect product Lie group" begin
    @testset "Show" begin
        M = LieGroups.Test.DummyManifold()
        op1 = LieGroups.Test.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroups.Test.DummySecondOperation()
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
        LieGroups.Test.test_lie_group(Gl, properties, expectations)

        Gr = RightSemidirectProductLieGroup(G1, G2, RightGroupOperationAction())
        properties = Dict(
            :Name => "The Right Semidirect Product Manifold",
            # :Rng => Random.MersenneTwister(),
            :Functions => fcts,
        )
        expectations = Dict(
            :repr => "RightSemidirectProductLieGroup($(G1), $(G2), RightGroupOperationAction())",
        )
        LieGroups.Test.test_lie_group(Gr, properties, expectations)
    end
    @testset "Interaction os semidirect and product" begin
        # Make sure the inner product manifold stays “nested” and does not collapse to a 4-product manifold
        M = LieGroups.Test.DummyManifold()
        op1 = LieGroups.Test.DummyOperation()
        G1 = LieGroup(M, op1)
        op2 = LieGroups.Test.DummySecondOperation()
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
        g1 = exp(SpecialOrthogonalGroup(2), [0.0 0.1; -0.1 0.0])
        g2 = [0.0 -1.0; 1.0 0.0]
        g3 = [1.0 0.0; 0.0 1.0]
        h1 = [0.1, 0.2]
        h2 = [0.0, 1.0]
        h3 = [0.0, 0]

        X1 = [0.0 0.1; -0.1 0.0]
        Y1 = [0.0, 0.2]

        # sanity check inv with matrix representation
        # reference data for (SO(2), ℝ²) semidirect products
        l1_mat = vcat([1 0 0], hcat(h1, g1)) # build the matrix representation for a left action semidirect product - acts on the right component
        r1_mat = hcat(vcat(g1, h1'), [0, 0, 1]) # build the matrix representation for a right action semidirect product - acts on the left component
        # inverse of R should be the same for both cases
        inv_g1 = inv(SpecialOrthogonalGroup(2), g1)
        # for t components, use the matrix inverse as reference
        l_inv_h1 = inv(l1_mat)[2:3, 1] # extract the translation part
        r_inv_h1 = inv(r1_mat)[3, 1:2] # extract the translation part

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
                :Functions => [identity_element, is_identity, inv, compose, show],
            )

            # use the matrix inverse as reference
            if action_on == ActionActsOnLeft()
                inv_h1 = r_inv_h1 # the standard right action semidirectproduct acts on the left side
            else
                inv_h1 = l_inv_h1 # the standard left action semidirectproduct acts on the right side
            end

            # also compare the matrix inverse h-component to the applied action
            H = SpecialOrthogonalGroup(2)
            N = TranslationGroup(2)
            if (action == TestLeftAction()) != (action_on == ActionActsOnLeft())
                # left right || right left -> (g,h)⁻¹ =  (g⁻¹, α(g⁻¹)(h⁻¹)
                α_g_h = apply(GroupAction(action, H, N), inv(H, g1), inv(N, h1))
            else
                # left left || right right -> inverse action -> (g,h)⁻¹ =  (g⁻¹, α(g⁻¹)⁻¹(h⁻¹)
                α_g_h = apply(GroupAction(action, H, N), g1, inv(N, h1))
            end
            @test isapprox(α_g_h, inv_h1)

            expectations = Dict(
                :atol => 1.0e-14,
                :inv => ArrayPartition(inv_g1, inv_h1),
            )
            LieGroups.Test.test_lie_group(G, properties, expectations)

            G = RightSemidirectProductLieGroup(TranslationGroup(2), SpecialOrthogonalGroup(2), action; action_on)
            q1 = ArrayPartition(copy(h1), copy(g1))
            q2 = ArrayPartition(copy(h2), copy(g2))
            q3 = ArrayPartition(copy(h3), copy(g3))
            W1 = ArrayPartition(copy(Y1), copy(X1))
            properties = Dict(
                :Name => "RightSemidirectProductLieGroup, $(supertype(typeof(action))), $(action_on)",
                :Points => [q1, q2, q3],
                :Vectors => [W1],
                :Functions => [identity_element, is_identity, inv, compose, show],
            )

            expectations = Dict(
                :atol => 1.0e-14,
                :inv => ArrayPartition(inv_h1, inv_g1),
            )

            LieGroups.Test.test_lie_group(G, properties, expectations)
        end
    end
end
