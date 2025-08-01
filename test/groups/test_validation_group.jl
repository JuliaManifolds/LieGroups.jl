using LieGroups, Random, Test, ManifoldsBase

s = joinpath(@__DIR__, "..", "LieGroupsTestSuite.jl")
!(s in LOAD_PATH) && (push!(LOAD_PATH, s))
using LieGroupsTestSuite

@testset "Validation Lie group" begin
    G = TranslationGroup(3)
    g1, g2, g3 = [1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [1.1, 1.2, 3.3]
    X1, X2, X3 = [0.0, 1.0, 0.0], [2.0, 0.0, 0.0], [0.1, 0.2, 0.3]
    VG = ValidationLieGroup(G)
    vg1, vg2, vg3 = ValidationMPoint.([g1, g2, g3])
    vX1, vX2, vX3 = ValidationLieAlgebraTangentVector.([X1, X2, X3])
    fcts = [
        adjoint,
        compose,
        conjugate,
        diff_conjugate,
        diff_inv,
        diff_left_compose,
        diff_right_compose,
        exp,
        hat,
        identity_element,
        inv,
        inv_left_compose,
        inv_right_compose,
        is_identity,
        lie_bracket,
        jacobian_conjugate,
        log,
        rand,
        show,
        vee,
    ]
    properties = Dict(
        :Name => "Validation of Translation group",
        :Points => [g1, g2, g3],
        :Vectors => [X1, X2, X3],
        :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :repr => "ValidationLieGroup of TranslationGroup(3; field=ℝ)\n    * mode = :error\n",
        :diff_inv => -X1,
        :diff_left_compose => X1,
        :diff_right_compose => X1,
        :lie_bracket => zero(X1),
    )
    test_lie_group(VG, properties, expectations)

    properties = Dict(
        :Name => "Validation of Translation group with types points/vectors",
        :Points => [vg1, vg2, vg3],
        :Vectors => [vX1, vX2, vX3],
        :Rng => Random.MersenneTwister(),
        :Functions => fcts,
    )
    expectations = Dict(
        :Name => "ValidationLieGroup of TranslationGroup(3; field=ℝ)\n    * mode = :error\n",
        :diff_left_compose => X1,
        :diff_right_compose => X1,
        :lie_bracket => zero(X1),
    )
    test_lie_group(VG, properties, expectations)
    @testset "Constructors and corner cases" begin
        # Constructor that does not also wrap the manifold
        @test ValidationLieGroup(G, false).lie_group isa TranslationGroup
        # the Lie algebra tangent vectors do not “double wrap”
        @test LieGroups.ValidationLieAlgebraTangentVector(vX1) === vX1
        # unwrap works also in ManifoldsBase
        ManifoldsBase.internal_value(vX1) === X1
        # non-typed default hat
        @test hat(LieAlgebra(VG), [1.0, 2.0, 3.0]).value == [1.0, 2.0, 3.0]
        e = ones(3)
        # check inplace disambiguation
        @test identity_element!(VG, e) == zeros(3)
        @test is_vector(VG, Identity(VG), X1)
        @test is_vector(VG, Identity(VG), vX1)
        # fallback on Real fails here since X=2.0 is not a tangent vector
        @test_throws DomainError norm(LieAlgebra(VG), 2.0)
        # Deactivate test
        VG2 = ValidationLieGroup(G; ignore_contexts = [:Input])
        @test norm(LieAlgebra(VG2), 2.0) isa Number
        # pass through
        @test representation_size(VG) == (3,)
        Y = zero_vector(LieAlgebra(VG))
        @test Y isa ValidationLieAlgebraTangentVector
        @test Y.value == zeros(3)
    end
    @testset "_msg defaults w/strings" begin
        @test_logs (:warn, "msg") LieGroups._msg(VG, "msg"; error = :warn)
        @test_logs (:info, "msg") LieGroups._msg(VG, "msg"; error = :info)
        @test LieGroups._msg(VG, "msg"; error = :nothing) === nothing
        @test_logs LieGroups._msg(VG, "msg"; error = :none)
        @test_throws ErrorException LieGroup._msg(VG, "msg"; error = :error)
        # same with error
        @test_logs (:warn,) LieGroups._msg(VG, DomainError("msg"); error = :warn)
        @test_logs (:info,) LieGroups._msg(VG, DomainError("msg"); error = :info)
        @test LieGroups._msg(VG, DomainError("msg"); error = :nothing) === nothing
        @test_logs LieGroups._msg(VG, DomainError("msg"); error = :none)
        @test_throws DomainError LieGroups._msg(VG, DomainError("msg"); error = :error)
    end
    @testset "_vlc responses when to exclude certain functions or contexts" begin
        VG1 = ValidationLieGroup(G; ignore_functions = Dict(exp => :All))
        VG2 = ValidationLieGroup(G; ignore_functions = Dict(exp => :Input))
        VG3 = ValidationLieGroup(G; ignore_contexts = [:Input])
        # VG1: checks disabled for all of exp
        @test !LieGroups._vLc(VG1, exp, :Input)
        @test !LieGroups._vLc(VG1, exp, :Ouput)
        # but others are not
        @test LieGroups._vLc(VG1, log, :Input)
        @test LieGroups._vLc(VG1, log, :Ouput)
        @test LieGroups._vLc(VG1, nothing, :Ouput)
        # VG2: checks disabled for input of exp, but not :Output
        @test !LieGroups._vLc(VG2, exp, :Input)
        @test LieGroups._vLc(VG2, exp, :Ouput)
        # VG3: checks disabled for all Inputs
        @test !LieGroups._vLc(VG3, exp, :Input)
        @test LieGroups._vLc(VG3, exp, :Ouput)
        @test !LieGroups._vLc(VG3, log, :Input)
        @test LieGroups._vLc(VG3, log, :Ouput)
        @test LieGroups._vLc(VG3, log, (:Ouput, :Point))
        @test !LieGroups._vLc(VG3, log, (:Ouput, :Input))
        @test LieGroups._vLc(VG3, nothing, :Ouput)
        # generic fallbacks
        @test LieGroups._vLc(:a, :b)
        @test !LieGroups._vLc((:a, :b), :b)
        @test LieGroups._vLc((:a, :c), :b)
    end
    @testset "isapprox passthrough" begin
        e = Identity(G)
        e2 = Identity(LieGroups.MatrixMultiplicationGroupOperation)
        @test !isapprox(VG, e2, e)
        @test !isapprox(VG, e, e2)
        @test isapprox(VG, e, e)
        @test_throws DomainError isapprox(VG, e2, e2)
        @test_throws DomainError isapprox(VG, e2, g1)
        @test_throws DomainError isapprox(VG, g1, e2)
    end
end
