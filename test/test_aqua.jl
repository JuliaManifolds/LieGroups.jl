using Aqua, LieGroups, Test

@testset "Aqua.jl" begin
    Aqua.test_all(LieGroups; ambiguities=(broken=false, exclude=[
        Base.:+, #ambiguities between Manifolds.Identity and LIeGroups.Identity
        Base.:-, #ambiguities between Manifolds.Identity and LIeGroups.Identity
    ]))
end
