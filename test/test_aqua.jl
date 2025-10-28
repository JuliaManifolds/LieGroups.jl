using Aqua, LieGroups, Test, LinearAlgebra

@testset "Aqua.jl" begin
    Aqua.test_all(
        LieGroups;
        ambiguities = (
            broken = false,
            exclude = [
                getindex, # ambiguities from convenience access methods; they were manually verified as safe
            ],
        ),
    )
end
