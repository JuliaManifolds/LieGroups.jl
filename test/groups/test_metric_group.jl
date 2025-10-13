using LieGroups, ManifoldsBase, Random, Test, RecursiveArrayTools

using LieGroups: MetricLieGroup

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "MetricLieGroup: A metric decorator for LieGroups" begin

end
