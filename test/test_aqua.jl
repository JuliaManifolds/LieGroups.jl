using Aqua, LieGroups, Test, LinearAlgebra

@testset "Aqua.jl" begin
    Aqua.test_all(
        LieGroups;
        ambiguities = false,
        # (
        #     broken = false,
        #     exclude = [
        #         getindex, # ambiguities from convenience access methods; they were manually verified as safe
        #     ],
        # ),
    )
end

# Interims solution until we properly fix the ambiguities
@testset "Ambiguities" begin
    ms = Test.detect_ambiguities(LieGroups)
    ms = [mm for mm in ms if mm[1].name != :getindex]
    LG_LIMIT = 6
    println("Number of LieGroups.jl ambiguities: $(length(ms))")
    if length(ms) > LG_LIMIT
        for amb in ms
            println(amb)
            println()
        end
    end
    @test length(ms) <= LG_LIMIT
end
