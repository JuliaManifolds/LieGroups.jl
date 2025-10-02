using LieGroups, Test
using Manifolds: Euclidean

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Columnwise Actions" begin
    G = SpecialOrthogonalGroup(3)
    M = Euclidean(3)
    cga = ColumnwiseGroupAction(LeftMultiplicationGroupAction())
    a = GroupAction(cga, G, M)
    s = 1 / sqrt(2)
    g = [1.0 0.0 0.0; 0.0 s -s; 0.0 s s]
    p = [1.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 2.0]
    gp = g * p
    @test apply(a, g, p) ≈ gp
    q = similar(p)
    apply!(a, q, g, p)

end
