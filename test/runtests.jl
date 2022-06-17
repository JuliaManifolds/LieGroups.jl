using LieGroups
using LinearAlgebra
using Test

@testset "LieGroups.jl" begin
    include("rotations.jl")
    include("rigid_motions.jl")
end
