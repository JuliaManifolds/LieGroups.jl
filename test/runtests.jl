using Manifolds, LieGroups, Test

function include_test(path)
    @info "Testing $path"
    @time include(path)  # show basic timing, (this will print a newline at end)
end

@testset "Lie Groups" begin
    include_test("groups/translation.jl")
end
