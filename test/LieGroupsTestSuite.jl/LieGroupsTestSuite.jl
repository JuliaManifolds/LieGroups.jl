"""
    LieGroupsTestSuite.jl

This module provides tools and dummy structures to test functionality provided
within `LieGroups.jl`
"""
module LieGroupsTestSuite
using LieGroups
using Test

"""
    test_compose(G, g, h)

Test functionality of `compose`

# Keyword arguments
* `test_inverse=true`: test that `g^{-1}g` is the identity (requires `inv`, `inv!`, and `is_identity`)
* `test_identity=true`: test that composing with the identity yields the identity (requires `identity_element`)
"""
function test_compose(G::LieGroup, g, h; test_inverse=true, test_identity=true)
    @testset "compose(G, g, h)" begin
        k1 = compose(G, g, h)
        k2 = copy(G, g)
        compose!(G, k2, g, h)
        @test isapprox(G, k1, k2)
        if test_inverse
            for g_ in [g, h]
                g_inv = inv(G, g_)
                k1 = compose(G, g_inv, g_)
                compose!(G, k2, g_inv, g_)
                @test isapprox(G, k1, k2)
                @test is_identity(G, k1)
                @test is_identity(G, k2)
            end
        end
        if test_identity
            for g_ in [g, h]
                for e in [Identity(G), identity_element(G)]
                    k1 = compose(G, g_, e)
                    compose!(G, k2, g_, e)
                    @test isapprox(G, k1, k2)
                    k1 = compose(G, e, g_)
                    compose!(G, k2, e, g_)
                    @test isapprox(G, k1, k2)
                end
            end
        end
    end
    return nothing
end

"""
    test_show(G, repr_string)

Test that show methods work as expected.
For now this (only) checks that `"\$G"` yields the `repr_string`.

requires `show` (or `repr`) to be implemented.
"""
function test_show(G::LieGroup, repr_string)
    @testset "repr(G, g, h)" begin
        @test repr(G) == repr_string
    end
    return nothing
end

"""
    test_LieGroup(G, properties, expectations)

Test the Lie group ``G`` based on a `Dict` of properties and a `Dict` of `expectations

Possible properties are

* `:Functions` is a vector of all defined functions for `G`
  Note that if `f` is in `:Functions`, and `f!` makes sense, for example for `compose`,
  it is assumed that both are defined.
* `:points` is a vector of at least three points on `G`
* `:vectors` is a vector of at least 3 elements from the Lie algebra `𝔤` og `G`
* `:name` is a name of the test. If not provided, defaults to `"\$G"`

Possible `expectations` are

* `:repr` is a sting one gets from `repr(G)`
* `:atol` a global absolute tolerance, defaults to `1e-8`
"""
function test_LieGroup(G::LieGroup, properties::Dict, expectations::Dict=Dict())
    a_tol = get(expectations, :atol, 1e-8)
    functions = get(properties, :Functions, Function[])
    points = get(properties, :Points, [])
    @assert length(points) > 2
    vectors = get(properties, :Vectors, [])
    @assert length(vectors) > 2
    test_name = get(properties, :Name, "$G")
    @testset "$(test_name)" begin
        # Call function tests based on their presence in alphabetical order
        if (compose in functions)
            ti = all(in.([inv, is_identity], Ref(functions)))
            test_compose(G, points[1], points[2]; test_inverse=ti)
        end
        if (any(in.([show, repr], Ref(functions)))) && haskey(expectations, :repr)
            test_show(G, expectations[:repr])
        end
    end
end

export test_LieGroup
end # module
