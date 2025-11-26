using LieGroups, Test
using Manifolds: Euclidean

@testset "Multiplication Action" begin
    G = SpecialOrthogonalGroup(3)
    M = Euclidean(3)
    a = LeftMultiplicationGroupAction()
    A = GroupAction(a, G, M)

    g = [0.7620105424968724 -0.01728550575563964 0.6473338739896082; -0.3990906109981697 0.7746993466910932 0.4904769173462713; -0.5099672708485612 -0.6320934531595508 0.583430586390622]
    p = [4.0, 5.0, 6.0]
    gp = g * p
    @test apply(A, g, p) ≈ gp
    q = similar(p)
    apply!(A, q, g, p)
    @test q ≈ gp

    X = [0.0 -3.0 2.0; 3.0 0.0 -1.0; -2.0 1.0 0.0]
    Z = similar(p)
    diff_group_apply!(A, Z, Identity(G), p, X)
    @test Z ≈ X * p
    @test diff_group_apply(A, Identity(G), p, X) ≈ X * p
end
