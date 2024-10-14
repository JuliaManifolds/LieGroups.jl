using LieGroups, Test

s = joinpath(@__DIR__, "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

function include_test(path)
    @info "Testing $path"
    @time include(path)  # show basic timing, (this prints a newline at end)
end

@testset "LieGroups.jl" begin
    @testset "Lie Group Interface" begin
        include_test("test_interface.jl")
    end
    @testset "Lie Groups" begin
        include_test("groups/test_translation_group.jl")
    end
end
