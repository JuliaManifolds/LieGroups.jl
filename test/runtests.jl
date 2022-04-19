using LieGroup
using LinearAlgebra
using Test

@testset "LieGroup.jl" begin
    include("rotations.jl")
    include("rigid_motions.jl")
end
