"""
    LieGroupsTestSuite.jl

This module provides tools and dummy structures to test functionality provided
within `LieGroups.jl`.

For every test function, several interactions to other functions can be activated.
The following functions are expected to be available, since their defaults just pass through to the manifold
* `is_point` both on the Lie group `G` and the Lie algebra `𝔤`
* `isapprox(G,g,h)` and `issaprox(𝔤, X, Y)`
* `copy(G, g)`
* `norm(𝔤, X)`
"""
module LieGroupsTestSuite
using LieGroups
using Test

#
#
# Dummy Types
struct DummyOperation <: AbstractGroupOperation end
struct DummySecondOperation <: AbstractGroupOperation end
struct DummyManifold <: LieGroups.AbstractManifold{LieGroups.ℝ} end
# Test functionality of single functions
#
#
# --- A

"""
    test_adjoint(G, g, X; kwargs...)

Test functionality of the `adjoint` function for a given Lie group element `g` and a Lie Algebra vector `X`

# Keyword arguments
* `expected_value=missing` provide the value expected. If none is provided, the
  default from `diff_conjugate` is used
"""
function test_adjoint(G::LieGroup, g, X; expected_value=missing)
    @testset "adjoint" begin
        v = if ismissing(expected_value)
            diff_conjugate(G, g, identity_element(G), X)
        else
            expected_value
        end
        𝔤 = LieAlgebra(G)
        Y1 = adjoint(G, g, X)
        Y2 = copy(𝔤, X)
        adjoint!(G, Y2, g, X)
        @test isapprox(𝔤, Y1, Y2)
        @test isapprox(𝔤, Y1, v)
    end
    return nothing
end
#
#
# --- C
"""
    test_compose(G, g, h; kwargs...)

Test functionality of `compose` for given Lie group elements `g`, `h`.

# Keyword arguments

* `test_inverse=true`: test that `g^{-1}g` is the identity (requires `inv`, `inv!`, and `is_identity`)
* `test_identity=true`: test that composing with the identity yields the identity (requires `identity_element`)
"""
function test_compose(G::LieGroup, g, h; test_inverse=true, test_identity=true)
    @testset "compose" begin
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
            e = Identity(G)
            k3 = copy(G, g)
            compose!(G, k3, e, e)
            @test is_identity(G, k3)
        end
    end
    return nothing
end

"""
    test_conjugate(G::LieGroup, g, h; expected_value=missing)

Test functionality of `conjugate`.

# Keyword arguments
* `expected_value=missing`: the result of the conjugate can also be provided directly,
  then neither `compose` nor `inv`  are not required.
"""
function test_conjugate(G::LieGroup, g, h; expected_value=missing)
    @testset "conjugate" begin
        v = if ismissing(expected_value)
            compose(G, g, compose(G, h, inv(G, g)))
        else
            expected_value
        end
        k1 = conjugate(G, g, h)
        k2 = copy(G, g)
        conjugate!(G, k2, g, h)
        @test isapprox(G, k1, k2)
        @test isapprox(G, k1, v)
    end
    return nothing
end

"""
    test_diff_inv(G::LieGroup, g, X; expected_value=missing)

Test functionality of `diff_inv`.

# Keyword arguments
* `expected_value=missing`: the result of the differential of the inverse, if not provided,
  only consistency between the allocating and the in-place variant is checked.
"""
function test_diff_inv(G::LieGroup, g, X; expected_value=missing)
    @testset "diff_inv" begin
        𝔤 = LieAlgebra(G)
        Y1 = diff_inv(G, g, X)
        Y2 = copy(𝔤, X)
        Y2 = diff_inv!(G, Y2, g, X)
        @test isapprox(LieAlgebra(G), Y1, Y2)
        if !ismissing(expected_value)
            @test isapprox(LieAlgebra(G), Y1, expected_value)
        end
    end
end

"""
    test_diff_left_compose(G::LieGroup, g, h, X; expected_value=missing)

Test functionality of `diff_left_compose`.

# Keyword arguments
* `expected_value=missing`: the result of the differential of the compose's left argument,
  if not provided, only consistency between the allocating and the in-place variant is checked.
"""
function test_diff_left_compose(G::LieGroup, g, h, X; expected_value=missing)
    @testset "diff_left_compose" begin
        𝔤 = LieAlgebra(G)
        Y1 = diff_left_compose(G, g, h, X)
        Y2 = copy(𝔤, X)
        Y2 = diff_left_compose!(G, Y2, g, h, X)
        @test isapprox(LieAlgebra(G), Y1, Y2)
        if !ismissing(expected_value)
            @test isapprox(LieAlgebra(G), Y1, expected_value)
        end
    end
end

"""
    test_diff_right_compose(G::LieGroup, g, h, X; expected_value=missing)

Test functionality of `diff_right_compose`.

# Keyword arguments
* `expected_value=missing`: the result of the differential of the compose's right argument,
  if not provided, only consistency between the allocating and the in-place variant is checked.
"""
function test_diff_right_compose(G::LieGroup, g, h, X; expected_value=missing)
    @testset "diff_inv" begin
        𝔤 = LieAlgebra(G)
        Y1 = diff_left_compose(G, g, h, X)
        Y2 = copy(𝔤, X)
        Y2 = diff_left_compose!(G, Y2, g, h, X)
        @test isapprox(LieAlgebra(G), Y1, Y2)
        if !ismissing(expected_value)
            @test isapprox(LieAlgebra(G), Y1, expected_value)
        end
    end
end

"""
    test_copyto(G, g)

Test that `copyto!` works also when copying over an `Identity`.

The point `g` can be any point _but_ the `identity_element`.
"""
function test_copyto(G, g)
    @testset "copyto!" begin
        k = copy(G, g)
        e = Identity(G)
        copyto!(G, k, e)
        # also check that this holds both ways
        @test isapprox(G, k, Identity(G))
        @test isapprox(G, Identity(G), k)
        # Test that copying back also works
        copyto!(G, k, g)
        @test isapprox(G, k, g)
        # copy into identity only works if provided the identity
        @test copyto!(G, e, Identity(G)) == e
        # but then also always returns e
        copyto!(G, k, e)
        @test copyto!(G, e, k) == e
        # and fails if the point to copy in is not (numerically) e
        @test_throws DomainError copyto!(G, e, g)
        return nothing
    end
end

#
#
# --- E
"""
    test_exp_log(G, g, h, X)

Test functionality of `exp` and `log` for given Lie group elements `g`, `h` and
a vector `X` from the Lie Algebra.

!!! note
    This function requires `is_point(G, g)` and `is_point(LieAlgebra(G), X)` to be implemented

# Keyword arguments

* `test_exp=true`: test the exponential map yields a point on `G`
* `test_log=true`: test the logarithmic map.
"""
function test_exp_log(G::LieGroup, g, h, X; test_exp=true, test_log=true)
    @testset "(Lie group) exp & log" begin
        𝔤 = LieAlgebra(G)
        e = Identity(G)
        if test_exp
            # Lie group exp
            k1 = exp(G, e, X)
            k2 = copy(G, g)
            exp!(G, k2, e, X)
            @test isapprox(G, k1, k2)
            @test is_point(G, k1)
            # exp
            k1 = exp(G, g, X)
            k2 = copy(G, g)
            exp!(G, k2, g, X)
            @test isapprox(G, k1, k2)
            @test is_point(G, k1)
        end
        if test_log
            # Lie group log
            Y1 = log(G, e, g)
            Y2 = copy(G, g)
            log!(G, Y2, e, g)
            @test isapprox(𝔤, Y1, Y2)
            @test is_point(𝔤, Y1)
            @test norm(𝔤, log(G, g, g)) ≈ 0
            @test norm(𝔤, log(G, h, h)) ≈ 0
            # log
            Y1 = log(G, g, h)
            Y2 = copy(G, g)
            log!(G, Y2, g, h)
            @test isapprox(𝔤, Y1, Y2)
            @test is_point(𝔤, Y1)
            # or equivalently
            @test is_vector(G, Y1)
            @test is_vector(G, Identity(G), Y1)
            @test norm(𝔤, log(G, g, g)) ≈ 0
            @test norm(𝔤, log(G, h, h)) ≈ 0
        end
        if test_exp && test_log
            # Lie group exp / log
            k1 = exp(G, e, X)
            k2 = copy(G, g)
            exp!(G, k2, e, X)
            Y1 = log(G, e, k1)
            Y2 = copy(G, g)
            log!(G, Y2, e, k2)
            @test isapprox(𝔤, Y1, Y2)
            @test isapprox(𝔤, X, Y1)
            # exp & log
            k1 = exp(G, g, X)
            k2 = copy(G, g)
            exp!(G, k2, g, X)
            Y1 = log(G, g, k1)
            Y2 = copy(G, g)
            log!(G, Y2, g, k2)
            @test isapprox(𝔤, Y1, Y2)
            @test isapprox(𝔤, X, Y1)
        end
    end
    return nothing
end

#
#
# --- S
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

# The global test function for a Lie group
#
#
"""
    test_LieGroup(G, properties, expectations)

Test the Lie group ``G`` based on a `Dict` of properties and a `Dict` of `expectations

Possible properties are

* `:Functions` is a vector of all defined functions for `G`
  Note that if `f` is in `:Functions`, and `f!` makes sense, for example for `compose`,
  it is assumed that both are defined.
* `:points` is a vector of at least three points on `G`, the first is not allowed to be the identity numerically
* `:vectors` is a vector of at least 3 elements from the Lie algebra `𝔤` og `G`
* `:Name` is a name of the test. If not provided, defaults to `"\$G"`

Possible `expectations` are

* `:repr` is a sting one gets from `repr(G)`
* `:adjoint` for the result of `conjgate` in the case where `diff_conjugate` is not implemented
* `:diff_inv` for the result of `diff_inv` with respect to the first point and the first vector.
* `:diff_left_compose` for the result of `diff_left_compose` with respect to the first two points and the first vector.
* `:diff_right_compose` for the result of `diff_right_compose` with respect to the first two points and the first vector.
* `:atol` a global absolute tolerance, defaults to `1e-8`
* `:conjugate` for the result of `conjgate in the case where `compose`, `inv` are not implemented
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
        #
        #
        # --- A
        if (adjoint in functions)
            v = get(expectations, :adjoint, missing)
            test_adjoint(G, points[1], vectors[1]; expected_value=v)
        end
        #
        #
        # --- C
        if (compose in functions)
            ti = all(in.([inv, is_identity], Ref(functions)))
            test_compose(G, points[1], points[2]; test_inverse=ti)
        end
        # since there is a default, also providing compose&inv suffices
        if (conjugate in functions) || (all(in.([compose, inv], Ref(functions))))
            v = get(expectations, :conjugate, missing)
            test_conjugate(G, points[1], points[2]; expected_value=v)
        end
        # Either `copyto` or the default with `identity_element`` available
        if any(in.([copyto!, identity_element], Ref(functions)))
            test_copyto(G, points[1])
        end
        #
        #
        # --- D
        if (diff_inv in functions)
            v = get(expectations, :diff_inv, missing)
            test_diff_inv(G, points[1], vectors[1]; expected_value=v)
        end

        if (diff_left_compose in functions)
            v = get(expectations, :diff_left_compose, missing)
            test_diff_left_compose(G, points[1], points[2], vectors[1]; expected_value=v)
        end

        if (diff_right_compose in functions)
            v = get(expectations, :diff_right_compose, missing)
            test_diff_left_compose(G, points[1], points[2], vectors[1]; expected_value=v)
        end

        #
        #
        # --- E
        if any(in.([exp, log], Ref(functions)))
            test_exp_log(
                G,
                points[1],
                points[2],
                vectors[1];
                test_exp=(exp in functions),
                test_log=(log in functions),
            )
        end
        #
        #
        # --- S
        if (any(in.([show, repr], Ref(functions)))) && haskey(expectations, :repr)
            test_show(G, expectations[:repr])
        end
    end
end

export test_LieGroup
end # module
