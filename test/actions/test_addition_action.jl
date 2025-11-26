using LieGroups, Test
using Manifolds: Euclidean

@testset "Addition Action" begin
    G = TranslationGroup(3)
    M = Euclidean(3)
    a = AdditionGroupAction()
    A = GroupAction(a, G, M)
    g = [1.0, 2.0, 3.0]
    p = [4.0, 5.0, 6.0]
    gp = g + p
    @test apply(A, g, p) ≈ gp
    q = similar(p)
    apply!(A, q, g, p)
    @test q ≈ gp

    X = [1.0, 0.0, 0.0]
    Y = similar(X)
    diff_apply!(A, Y, g, p, X)
    @test Y ≈ X

    Z = similar(X)
    diff_group_apply!(A, Z, g, p, X)
    @test Z ≈ X
end
