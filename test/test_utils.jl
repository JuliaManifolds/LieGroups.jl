using LieGroups, Test

@testset "Utility functions" begin
    # `log_safe!`
    Y = zeros(2, 2)
    # Neg eigenvalues not allowed neither Hermitian nor `triu` nor else
    @test_throws DomainError LieGroups.log_safe!(Y, [-1.0 0.0; 0.0 -1.0])
    @test_throws DomainError LieGroups.log_safe!(Y, [-1.0 1.0; 0.0 -1.0])
    @test_throws DomainError LieGroups.log_safe!(Y, [-1.0 1.0; 0.5 -1.0])
    # Check `triu` case is approx the same as normal log
    LieGroups.log_safe!(Y, [1.0 0.0; 0.4 1.0]) ≈ log([1.0 0.0; 0.4 1.0])
    # Check Schur `triu` case is approx the same as normal log
    LieGroups.log_safe!(Y, [1.0 0.5; 0.4 1.0]) ≈ log([1.0 0.5; 0.4 1.0])
    # Default: Complex, use log
    Y2 = zeros(ComplexF64, 2, 2)
    LieGroups.log_safe!(Y2, [1.0 0.5; 0.4 1.0im]) ≈ log([1.0 0.5; 0.4 1.0im])
    # Special cases of the safe `usinc`
    @test LieGroups.usinc_from_cos(-1.004) == 0
    @test LieGroups.usinc_from_cos(1.004) == 1.0
end
