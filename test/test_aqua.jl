using Aqua, LieGroups, Test

@testset "Aqua.jl" begin
    Aqua.test_all(LieGroups; ambiguities=(broken=false, exclude=[
        Base.:+, #temporary ambiguities between Manifolds.Identity and LieGroups.Identity
        Base.:-, #temporary ambiguities between Manifolds.Identity and LieGroups.Identity
    ]))
end
